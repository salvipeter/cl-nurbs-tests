(in-package :cl-nurbs-tests)

(defparameter *exponent* 2)
(defparameter *resolution* 40)
(defparameter *tiny* 1.0d-5)

(defvar *alpha*)

(defun points-from-angles (angles &optional (radius 1.0d0))
  "Gives a list of points on the unit circle, divided by the given angles.
For example, for a equilateral triangle give '(120 120 120)."
  (iter (for degrees in angles)
	(for angle = (/ (* degrees pi) 180.0d0))
	(for alpha first angle then (+ alpha angle))
	(collect (v* (list (cos alpha) (sin alpha)) radius))))

(defun lines-from-points (points)
  "List of point pairs (last-first, first-second, second-third, ...)."
  (cons (list (car (last points)) (first points))
	(iter (for k from 0 below (1- (length points)))
	      (collect (list (elt points k) (elt points (1+ k)))))))

(defun angles-from-points (points)
  (let ((sum (+ (iter (for p in points)
		      (for q in (append (last points) points))
		      (sum (point-distance p q))))))
    (cons 0 (iter (for p in points)
		  (for q in (rest points))
		  (collect (* 360.0d0(/ (point-distance p q) sum)))))))

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

(defun compute-distance (type points p &optional no-tiny-p)
  (macrolet ((tiny-lambda ((args) &body body)
	       `(lambda (,args)
		  (if no-tiny-p
		      (let ((result (progn ,@body)))
			(if (< result *tiny*) 0.0d0 result))
		      (progn ,@body)))))
    (mapcar (ecase type
	      (perpendicular (tiny-lambda (lst)
			       (perpendicular-distance points lst p 'd)))
	      (barycentric (let ((lines (lines-from-points points)))
			     (tiny-lambda (lst)
			       (barycentric-distance lines lst p 'd))))
	      (chord-based (tiny-lambda (lst)
			     (chord-based-distance points lst p 'd)))
	      (radial (tiny-lambda (lst)
			(radial-distance points lst p 'd)))
	      (line-sweep (let ((center (central-point points (lines-from-points points) t)))
			    (tiny-lambda (lst)
			      (line-sweep-distance center lst p 'd)))))
	    (iter (for i from -2 below (- (length points) 2))
		  (collect (iter (for j from 0 below 4)
				 (collect (elt points (mod (+ i j) (length points))))))))))

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

(defun ribbon-blend (d i)
  "A blend function that is 1 on the Ith line, 0 on all others.
This naturally means that we have singularities in the corners."
  (cond ((notany (lambda (x) (< (abs x) *tiny*)) d) (blend d i))
	((< (elt d i) *tiny*)
	 (if (= (length (remove-if-not (lambda (di) (< di *tiny*)) d)) 1)
	     1.0d0
	     0.5d0))
	(t 0.0d0)))

(defun corner-blend (d i)
  "Gregory-style blend: 1 in the corner between the Ith and (I+1)th line,
0 in all other corners, smoothly decreasing along the Ith and (I+1)th line.
This eliminates the singularity problem in the corners."
  (let ((n (length d)))
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
	 (center (central-point points lines t))
	 (result (list center)))
    (iter (for j from 1 to *resolution*)
	  (for coeff = (/ j *resolution*))
	  (iter (for k from 0 below n)
		(iter (for i from 0 below j)
		      (for lp = (line-point (elt lines k) (/ i j)))
		      (push (affine-combine center coeff lp) result))))
    (nreverse result)))

(defun write-blends (angles on-off filename &key (blend-function #'ribbon-blend)
		     (distance-type 'perpendicular))
  "Computes samples by a set of blend functions and writes it in a VTK file.
The ON-OFF parameter declares which blends should be turned on.
For ANGLES, see POINTS-FROM-ANGLES."
  (let* ((n (length angles))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (*alpha* (compute-alpha lines (v* (reduce #'v+ points) (/ n))))
	 (vertices (mapcar (lambda (p)
			     (let ((d (compute-distance distance-type points p t)))
			       (cons (iter (for i from 0 below (length on-off))
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))
			   (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

;;; quadrilateral mesh
(defun write-blends-quad (angles on-off filename &key (blend-function #'ribbon-blend)
			  (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (n (length lines))
	 (res (1+ (* 2 *resolution*)))
	 (result (make-array (list res res))))
    (flet ((map-coordinates (x y)
	     (list (/ (- x *resolution*) *resolution*)
		   (/ (- y *resolution*) *resolution*))))
      (iter (for x from 0 below res)
	    (iter (for y from 0 below res)
		  (let ((p (map-coordinates x y)))
		    (setf (aref result x y)
			  (if (insidep lines p)
			      (let ((d (compute-distance distance-type points p t)))
				(cons (iter (for i from 0 below n)
					    (when (elt on-off i)
					      (sum (funcall blend-function d i))))
				      p))
			      (cons 0 p)))))))
    (write-points2-vtk result filename)))

;; (write-blends '(60 60 80 120 40) '(t t nil nil nil) "/tmp/blend.vtk")
;; (write-blends '(60 60 80 120 40) '(t t nil nil t) "/tmp/corner-blend.vtk" :blend-function #'corner-blend)
;; (write-blends '(60 60 80 120 40) '(t t nil nil nil t) "/tmp/interior-blend.vtk" :blend-function #'interior-ribbon-blend)

#+nil
(write-blends '(40 20 60 100 80) '(t nil nil nil nil) "/tmp/blend.vtk"
	      :blend-function #'ribbon-blend :distance-type 'line-sweep)

(defun write-blend-all-types (angles on-off directory)
  (iter (for blend in (list #'ribbon-blend #'corner-blend))
	(for blend-name in '(ribbon corner))
	(iter (for type in '(perpendicular barycentric chord-based radial line-sweep))
	      (for name = (format nil "~(~a~)-~(~a~)-~d.vtk" blend-name type *exponent*))
	      (for filename = (make-pathname :directory directory :name name))
	      (write-blends angles on-off filename :blend-function blend :distance-type type))))

#+nil
(let ((*exponent* 2)
      (*resolution* 40))
  (write-blend-all-types '(40 20 60 100 80) '(t nil nil nil nil) #p"/tmp/"))

(defun write-vtk-polylines (lines filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 1.0~
		 ~%B-spline Surface~
		 ~%ASCII~
		 ~%DATASET POLYDATA~
		 ~%POINTS ~d float~%" (reduce #'+ (mapcar #'length lines)))
    (dolist (line lines)
      (dolist (point line)
	(format s "~{~f~^ ~}~%" point))) 
    (format s "~%LINES ~d ~d~%" (length lines)
	    (+ (length lines) (reduce #'+ (mapcar #'length lines))))
    (let ((i 0))
      (dolist (line lines)
	(let ((n (length line)))
	  (format s "~d" n)
	  (dotimes (j n)
	    (format s " ~d" i)
	    (incf i))
	  (terpri s))))))

(defun write-blends-uv-polyline (angles on-off filename &key (blend-function #'ribbon-blend)
				 (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (n (length lines))
	 (res (1+ (* 2 *resolution*))))
    (flet ((map-coordinates (x y)
	     (list (/ (- x *resolution*) *resolution*)
		   (/ (- y *resolution*) *resolution*))))
      (write-vtk-polylines
       (iter (for x from 0 below res)
	     (for line = 
		  (iter (for y from 0 below res)
			(let ((p (map-coordinates x y)))
			  (when (insidep lines p)
			    (collect
			     (let ((d (compute-distance distance-type points p t)))
			       (cons (iter (for i from 0 below n)
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))))))
	     (when line (collect line)))
       filename))))

#+nil
(write-blends-uv-polyline '(40 20 60 100 80) '(t nil nil nil nil) "/tmp/blend2.vtk"
			  :blend-function #'ribbon-blend :distance-type 'line-sweep)

(defun spider-lines (n)
  (iter (with inner-start = 0)
	(with outer-vert = 1)
	(for layer from 1 to *resolution*)
	(for inner-vert = inner-start)
	(for outer-start = outer-vert)
	(collect
	 (iter (for side from 0 below n)
	       (appending
		(iter (with vert = 0)
		      (for next-vert = (if (and (= side (1- n))
						(= vert (1- layer)))
					   outer-start
					   (1+ outer-vert)))
		      (collect inner-vert)
		      (incf outer-vert)
		      (incf vert)
		      (while (/= vert layer))
		      (for inner-next = (if (and (= side (1- n))
						 (= vert (1- layer)))
					    inner-start
					    (1+ inner-vert)))
		      (setf inner-vert inner-next)))))
	(setf inner-start outer-start)))

(defun write-blends-spider-polyline (angles on-off filename &key (blend-function #'ribbon-blend)
				     (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((n (length angles))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (*alpha* (compute-alpha lines (v* (reduce #'v+ points) (/ n))))
	 (vertices (mapcar (lambda (p)
			     (let ((d (compute-distance distance-type points p t)))
			       (cons (iter (for i from 0 below (length on-off))
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))
			   (vertices points))))
    (write-vtk-polylines
     (mapcar (lambda (lst) (mapcar (lambda (i) (elt vertices i)) lst))
	     (spider-lines n))
     filename)))

#+nil
(write-blends-spider-polyline '(40 20 60 100 80) '(t nil nil nil nil) "/tmp/blend.vtk"
			      :blend-function #'ribbon-blend :distance-type 'line-sweep)

(defun generate-colors (n)
  (assert (<= n 6) (n) "Only N <= 6 is supported.")
  (subseq '((255 0 0) (0 255 0) (0 0 255)
            (255 255 0) (255 0 255) (0 255 255))
          0 n))

(defun write-color-blend-test (angles filename r &key (blend-function #'ribbon-blend)
			       (distance-type 'perpendicular) (trim '(0.89d0 0.91d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length angles))
	   (points (points-from-angles angles))
	   (lines (lines-from-points points))
	   (colors (generate-colors n)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (compute-distance distance-type points p t)) 
			 (blends (iter (for i from 0 below n)
				       (collect (funcall blend-function d i)))))
                    (if (and trim (some (lambda (x) (< (first trim) x (second trim))) blends))
			(format s "0 0 0~%")
			(format s "~{~d~^ ~}~%"
				(mapcar #'round
					(reduce (lambda (x y) (mapcar #'+ x y))
						(mapcar #'v* colors blends))))))
                  (format s "255 255 255~%")))))))))

;;; Display only the si=... lines
(defun write-si-line-test (angles filename r &key (distance-type 'perpendicular)
			   (trim '(0.495d0 0.505d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length angles))
	   (points (points-from-angles angles))
	   (lines (lines-from-points points))
	   (colors (generate-colors n)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (acond
	       ((not (insidep lines p))
		(declare (ignore it))
		(format s "255 255 255~%"))
	       ((iter (for line in lines)
		      (for i upfrom 0)
		      (when (< (point-line-distance p line) 0.01d0)
			(leave (elt colors i))))
		(format s "~{~d~^ ~}~%" it))
	       (t
		(declare (ignore it))
		(let* ((sp (compute-parameter distance-type points p t))
		       (on (and trim (some (lambda (x)
					     (and (< (first trim) x (second trim)) x))
					   sp))))
		  (if on
		      (format s "~{~d~^ ~}~%" (elt colors (position on sp)))
		      (format s "127 127 127~%"))))))))))))

#+nil
(write-color-blend-test '(40 20 60 100 80) "/tmp/blend.ppm" 200
			:blend-function #'ribbon-blend
			:distance-type 'chord-based
			:trim '(0.495d0 0.505d0))

#+nil
(write-si-line-test '(40 20 60 100 80) "/tmp/blend.ppm" 200
			:distance-type 'line-sweep
			:trim '(0.495d0 0.505d0))

(defun write-blend-voronoi (angles filename r &key (blend-function #'ribbon-blend)
			    (distance-type 'perpendicular) (threshold 0.01d0))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length angles))
	   (points (points-from-angles angles))
	   (lines (lines-from-points points)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (compute-distance distance-type points p t)) 
			 (blends (iter (for i from 0 below n)
				       (collect (funcall blend-function d i)))))
                    (if (destructuring-bind (a b)
			    (subseq (sort blends #'>) 0 2)
			  (< (- a b) threshold))
			(format s "0 0 0~%")
			(format s "127 127 127~%")))
                  (format s "255 255 255~%")))))))))

#+nil
(write-blend-voronoi '(40 20 60 100 80) "/tmp/blend2.ppm" 200
			:blend-function #'ribbon-blend
			:distance-type 'line-sweep
			:threshold 0.03d0)
