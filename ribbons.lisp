(in-package :cl-nurbs-tests)

(defparameter *exponent* 2)
(defparameter *resolution* 40)
(defparameter *tiny* 1.0d-5)

(defvar *alpha*)

(defun points-from-angles (angles)
  "Gives a list of points on the unit circle, divided by the given angles.
For example, for a equilateral triangle give '(120 120 120)."
  (iter (for degrees in angles)
	(for angle = (/ (* degrees pi) 180.0d0))
	(for alpha first angle then (+ alpha angle))
	(collect (list (cos alpha) (sin alpha)))))

(defun lines-from-points (points)
  "List of point pairs (last-first, first-second, second-third, ...)."
  (cons (list (car (last points)) (first points))
	(iter (for k from 0 below (1- (length points)))
	      (collect (list (elt points k) (elt points (1+ k)))))))

(defun point-line-distance (p line &optional signedp)
  (let* ((v (v- (second line) (first line)))
         (d (scalar-product (vnormalize (list (second v) (- (first v))))
                            (v- p (first line)))))
    (if signedp d (abs d))))

(defun line-point (line u)
  (affine-combine (first line) u (second line)))

(defun insidep (lines p)
  "Determines if P is inside the polygon defined by LINES."
  (every (lambda (line)
           (> (* (point-line-distance p line t)
                 (point-line-distance '(0 0) line t))
              0))
         lines))

(defun blend (d i)
  "Used by RIBBON-BLEND."
  (let ((n (length d)))
    (/ (iter (for j from 0 below n)
	     (when (/= i j)
	       (multiply (expt (elt d j) *exponent*))))
       (iter (for k from 0 below n)
	     (sum (iter (for j from 0 below n)
			(when (/= j k)
			  (multiply (expt (elt d j) *exponent*)))))))))

(defun ribbon-blend (lines p i)
  "A blend function that is 1 on the Ith line, 0 on all others.
This naturally means that we have singularities in the corners."
  (let ((d (mapcar (lambda (line) (point-line-distance p line)) lines)))
    (cond ((notany (lambda (x) (< (abs x) *tiny*)) d) (blend d i))
	  ((< (min (point-distance (first (elt lines i)) p)
		   (point-distance (second (elt lines i)) p))
	      *tiny*)
	   0.5d0)
	  ((< (point-line-distance p (elt lines i)) *tiny*) 1.0d0)
	  (t 0.0d0))))

(defun corner-blend (lines p i)
  "Gregory-style blend: 1 in the corner between the Ith and (I+1)th line,
0 in all other corners, smoothly decreasing along the Ith and (I+1)th line.
This eliminates the singularity problem in the corners."
  (let ((d (mapcar (lambda (line)
		     (let ((x (point-line-distance p line)))
		       (if (< (abs x) *tiny*) 0 x)))
		   lines))
	(n (length lines)))
    (/ (iter (for j from 0 below n)
	     (when (and (/= i j) (/= (mod (1+ i) n) j))
	       (multiply (expt (elt d j) *exponent*))))
       (iter (for k from 0 below n)
	     (sum (iter (for j from 0 below n)
			(when (and (/= k j) (/= (mod (1+ k) n) j))
			  (multiply (expt (elt d j) *exponent*)))))))))

(defun compute-alpha (lines center)
  "Computes alpha such that at the center point of the domain the weight of
the interior surface will be 1/(N+1), where N is the number of lines."
  (let ((c (mapcar (lambda (line) (point-line-distance center line)) lines))
	(n (length lines)))
    (/ (iter (for k from 0 below n)
	     (sum (iter (for j from 0 below n)
			(when (/= j k)
			  (multiply (expt (elt c j) *exponent*))))))
       (iter (for j from 0 below n)
	     (multiply (expt (elt c j) *exponent*)))
       n)))

(defun interior-blend (d i)
  "Used by INTERIOR-RIBBON-BLEND."
  (let ((n (length d)))
    (/ (if (= i n)
	   (* *alpha* (iter (for j from 0 below n)
			    (multiply (expt (elt d j) *exponent*))))
	   (iter (for j from 0 below n)
		 (when (/= i j)
		   (multiply (expt (elt d j) *exponent*)))))
       (+ (iter (for k from 0 below n)
		(sum (iter (for j from 0 below n)
			   (when (/= j k)
			     (multiply (expt (elt d j) *exponent*))))))
	  (* *alpha* (iter (for j from 0 below n)
			   (multiply (expt (elt d j) *exponent*))))))))

(defun interior-ribbon-blend (lines p i)
  "When using this blend, an extra on-off parameter should be added at the end."
  (let ((d (mapcar (lambda (line) (point-line-distance p line)) lines)))
    (cond ((notany (lambda (x) (< (abs x) *tiny*)) d) (interior-blend d i))
	  ((= i (length d))
	   (if (some (lambda (x) (< (abs x) *tiny*)) d)
	       0.0d0
	       (interior-blend d i)))
	  ((< (min (point-distance (first (elt lines i)) p)
		   (point-distance (second (elt lines i)) p))
	      *tiny*)
	   0.5d0)
	  ((< (point-line-distance p (elt lines i)) *tiny*) 1.0d0)
	  (t 0.0d0))))

;;; triangular mesh
(defun triangles (n)
  "Triangles given as vertex index triples to be used along with VERTICES."
  (iter (with result = '())
	(with inner-start = 0)
	(with outer-vert = 1)
	(for layer from 1 to *resolution*)
	(for inner-vert = inner-start)
	(for outer-start = outer-vert)
	(iter (for side from 0 below n)
	      (iter (with vert = 0)
		    (for next-vert = (if (and (= side (1- n))
					      (= vert (1- layer)))
					 outer-start
					 (1+ outer-vert)))
		    (push (list inner-vert outer-vert next-vert) result)
		    (incf outer-vert)
		    (incf vert)
		    (while (/= vert layer))
		    (for inner-next = (if (and (= side (1- n))
					       (= vert (1- layer)))
					  inner-start
					  (1+ inner-vert)))
		    (push (list inner-vert next-vert inner-next) result)
		    (setf inner-vert inner-next)))
	(setf inner-start outer-start)
	(finally (return (nreverse result)))))

(defun vertices (points)
  "Sample points for the (convex) polygon defined by POINTS."
  (let* ((lines (lines-from-points points))
	 (n (length lines))
	 (center (v* (reduce #'v+ points) (/ n)))
	 (result (list center)))
    (iter (for j from 1 to *resolution*)
	  (for coeff = (/ j *resolution*))
	  (iter (for k from 0 below n)
		(iter (for i from 0 below j)
		      (for lp = (line-point (elt lines k) (/ i j)))
		      (push (affine-combine center coeff lp) result))))
    (nreverse result)))

(defun write-blends (angles on-off filename &key (blend-function #'ribbon-blend))
  "Computes samples by a set of blend functions and writes it in a VTK file.
The ON-OFF parameter declares which blends should be turned on.
For ANGLES, see POINTS-FROM-ANGLES."
  (let* ((n (length angles))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (*alpha* (compute-alpha lines (v* (reduce #'v+ points) (/ n))))
	 (vertices (mapcar (lambda (p)
			     (cons (iter (for i from 0 below (length on-off))
					 (when (elt on-off i)
					   (sum (funcall blend-function lines p i))))
				   p))
			   (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

;;; quadrilateral mesh
(defun write-blends-quad (angles on-off filename &key (blend-function #'ribbon-blend))
  "See the documentation of WRITE-BLENDS."
  (let* ((lines (lines-from-points (points-from-angles angles)))
	 (n (length lines))
	 (res (1+ (* 2 *resolution*)))
	 (points (make-array (list res res))))
    (flet ((map-coordinates (x y)
	     (list (/ (- x *resolution*) *resolution*)
		   (/ (- y *resolution*) *resolution*))))
      (iter (for x from 0 below res)
	    (iter (for y from 0 below res)
		  (let ((p (map-coordinates x y)))
		    (setf (aref points x y)
			  (if (insidep lines p)
			      (cons (iter (for i from 0 below n)
					  (when (elt on-off i)
					    (sum (funcall blend-function lines p i))))
				    p)
			      (cons 0 p)))))))
    (write-points2-vtk points filename)))

;; (write-blends '(60 60 80 120 40) '(t t nil nil nil) "/tmp/blend.vtk")
;; (write-blends '(60 60 80 120 40) '(t t nil nil t) "/tmp/corner-blend.vtk" :blend-function #'corner-blend)
;; (write-blends '(60 60 80 120 40) '(t t nil nil nil t) "/tmp/interior-blend.vtk" :blend-function #'interior-ribbon-blend)
