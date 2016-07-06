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
		  (collect (* 360.0d0 (/ (point-distance p q) sum)))))))

(defun angles-from-curves (curves)
  (let ((sum (reduce #'+ (mapcar #'bezier-arc-length curves))))
    (iter (for curve in curves)
	  (collect (* 360.0d0 (/ (bezier-arc-length curve) sum))))))

(defun domain-from-curves-angular (curves)
  (labels ((angle (c1 c2)
	     (acos (scalar-product (vnormalize (bezier c1 1 1)) (vnormalize (bezier c2 0 1)))))
	   (rescale (lst)
	     (let* ((min (list (reduce #'min (mapcar #'first lst))
			       (reduce #'min (mapcar #'second lst))))
		    (max (list (reduce #'max (mapcar #'first lst))
			       (reduce #'max (mapcar #'second lst))))
		    (length (max (- (first max) (first min)) (- (second max) (second min)))))
	       (mapcar (lambda (p) (v+ (v* (v- p min) (/ 2.0d0 length)) '(-1 -1))) lst))))
    (let* ((angles (mapcar #'angle curves (append (rest curves) (list (first curves)))))
	   (angle-multiplier (/ (* (- (length curves) 2) pi) (reduce #'+ angles)))
	   (normalized-angles (mapcar (lambda (x) (* x angle-multiplier)) angles))
	   (lengths (mapcar #'bezier-arc-length curves))
	   (length-sum (reduce #'+ lengths))
	   (vertices (iter (for prev first '(0 0) then next)
			   (for length in lengths)
			   (for dir first 0 then (+ dir angle))
			   (for angle in normalized-angles)
			   (for next = (v+ prev (v* (list (cos dir) (sin dir)) length)))
			   (collect next)))
	   (difference (v* (car (last vertices)) -1)))
      (rescale
       (append (iter (for length in (butlast lengths))
		     (for accumulated first length then (+ accumulated length))
		     (for vertex in vertices)
		     (collect (v+ vertex (v* difference (/ accumulated length-sum)))))
	       '((0 0)))))))

(defun find-with-lower-limit (fn low)
  "Finds x where f(x) = 0 (with *EPSILON* tolerance).
Assumes that f(y) > 0 for y < x and f(y) < 0 for y > x."
  (let ((hi (iter (for hi first (* 2 low) then (* 2 hi))
		  (while (> (* (funcall fn low) (funcall fn hi)) 0))
		  (finally (return hi)))))
    (iter (for mid = (/ (+ low hi) 2.0d0))
	  (for fx = (funcall fn mid))
	  (while (> (abs fx) *epsilon*))
	  (if (> fx 0)
	      (setf low mid)
	      (setf hi mid))
	  (finally (return mid)))))

(defun circle-radius (lengths)
  (let ((sorted (sort (copy-list lengths) #'>)))
    (when (> (reduce #'+ (rest sorted)) (first sorted))
      (flet ((radius-error (r)
	       (- (iter (for li in lengths)
			(sum (asin (/ li (* 2 r)))))
		  pi)))
	(let ((low (/ (first sorted) 2.0d0)))
	  (if (>= (radius-error low) 0)
	      (find-with-lower-limit #'radius-error low)
	      (find-with-lower-limit
	       (lambda (r)
		 (- (asin (/ (first sorted) (* 2 r)))
		    (iter (for li in (rest sorted))
			  (sum (asin (/ li (* 2 r)))))))
	       low)))))))

(defun angles-from-lengths (lengths)
  (let ((r (circle-radius lengths)))
    (assert r (lengths) "these are not the side lengths of a cyclic polygon")
    (let* ((alpha (mapcar (lambda (li)
			    (/ (* (* 2 (asin (/ li (* 2 r)))) 180.0d0) pi))
			  lengths))
	   (max-alpha (reduce #'max alpha)))
      (substitute (- 360.0d0 (reduce #'+ (remove max-alpha alpha :count 1)))
		  max-alpha alpha :count 1))))

(defun domain-from-curves (curves &optional (type 'angular))
  (ecase type
    (regular (let ((n (length curves)))
	       (points-from-angles (cons 13 (iter (repeat (1- n)) (collect (/ 360 n)))))))
    (circular (points-from-angles (angles-from-curves curves)))
    (circular-mod (points-from-angles
		   (angles-from-lengths (mapcar #'bezier-arc-length curves))))
    (angular (domain-from-curves-angular curves))))

(defun point-line-distance (p line &optional signedp)
  (let* ((v (v- (second line) (first line)))
         (d (scalar-product (vnormalize (list (second v) (- (first v))))
                            (v- p (first line)))))
    (if signedp d (abs d))))

(defun line-point (line u)
  (affine-combine (first line) u (second line)))

(defun insidep (lines p)
  "Determines if P is inside the polygon defined by LINES."
  (let ((center (central-point (mapcar #'second lines) lines t)))
    (every (lambda (line)
	     (> (* (point-line-distance p line t)
		   (point-line-distance center line t))
		0))
	   lines)))

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

(defun compute-alpha (points alpha distance-type)
  (let* ((n (length points))
	 (lines (lines-from-points points))
	 (d (compute-parameter distance-type 'd points (central-point points lines t) t)))
    (/ (* (- 1.0d0 alpha)
	  (iter (for k from 0 below n)
		(sum (iter (for j from 0 below n)
			   (when (/= j k)
			     (multiply (expt (elt d j) *exponent*)))))))
       (* alpha
	  (iter (for j from 0 below n)
		(multiply (expt (elt d j) *exponent*)))))))

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

(defun interior-ribbon-blend (d i)
  (cond ((notany (lambda (x) (< (abs x) *tiny*)) d) (interior-blend d i))
	((= i (length d)) 0.0d0)
	((>= (elt d i) *tiny*) 0.0d0)
	((= (length (remove-if-not (lambda (di) (< di *tiny*)) d)) 2)
	 0.5d0)
	(t 1.0d0)))

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

(defun write-blends (points on-off filename &key (blend-function #'ribbon-blend)
		     (distance-type 'perpendicular))
  "Computes samples by a set of blend functions and writes it in a VTK file.
The ON-OFF parameter declares which blends should be turned on."
  (let* ((n (length points))
	 #+nil(lines (lines-from-points points))
	 (*alpha* 0 #+nil(compute-alpha lines (v* (reduce #'v+ points) (/ n)) distance-type))
	 (vertices (mapcar (lambda (p)
			     (let ((d (compute-parameter distance-type 'd points p t)))
			       (cons (iter (for i from 0 below (length on-off))
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))
			   (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

;;; quadrilateral mesh
(defun write-blends-quad (points on-off filename &key (blend-function #'ribbon-blend)
			  (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((lines (lines-from-points points))
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
			      (let ((d (compute-parameter distance-type 'd points p t)))
				(cons (iter (for i from 0 below n)
					    (when (elt on-off i)
					      (sum (funcall blend-function d i))))
				      p))
			      (cons 0 p)))))))
    (write-points2-vtk result filename)))

;; (write-blends (points-from-angles '(60 60 80 120 40)) '(t t nil nil nil) "/tmp/blend.vtk")
;; (write-blends (points-from-angles '(60 60 80 120 40)) '(t t nil nil t) "/tmp/corner-blend.vtk" :blend-function #'corner-blend)
;; (write-blends (points-from-angles '(60 60 80 120 40)) '(t t nil nil nil t) "/tmp/interior-blend.vtk" :blend-function #'interior-ribbon-blend)

#+nil
(write-blends (points-from-angles '(40 20 60 100 80)) '(t nil nil nil nil) "/tmp/blend.vtk"
	      :blend-function #'ribbon-blend :distance-type 'line-sweep)

(defun write-blend-all-types (points on-off directory)
  (iter (for blend in (list #'ribbon-blend #'corner-blend))
	(for blend-name in '(ribbon corner))
	(iter (for type in '(perpendicular barycentric chord-based radial line-sweep))
	      (for name = (format nil "~(~a~)-~(~a~)-~d.vtk" blend-name type *exponent*))
	      (for filename = (make-pathname :directory directory :name name))
	      (write-blends points on-off filename :blend-function blend :distance-type type))))

#+nil
(let ((*exponent* 2)
      (*resolution* 40))
  (write-blend-all-types (points-from-angles '(40 20 60 100 80)) '(t nil nil nil nil) #p"/tmp/"))

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

(defun write-blends-uv-polyline (points on-off filename &key (blend-function #'ribbon-blend)
				 (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((lines (lines-from-points points))
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
			     (let ((d (compute-parameter distance-type 'd points p t)))
			       (cons (iter (for i from 0 below n)
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))))))
	     (when line (collect line)))
       filename))))

#+nil
(write-blends-uv-polyline (points-from-angles '(40 20 60 100 80)) '(t nil nil nil nil)
			  "/tmp/blend2.vtk" :blend-function #'ribbon-blend
			  :distance-type 'line-sweep)

(defun spider-uv-lines (n)
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

(defun write-blends-spider-polyline (points on-off filename &key (blend-function #'ribbon-blend)
				     (distance-type 'perpendicular))
  "See the documentation of WRITE-BLENDS."
  (let* ((n (length points)) 
	 (lines (lines-from-points points))
	 (*alpha* (compute-alpha lines (v* (reduce #'v+ points) (/ n)) distance-type))
	 (vertices (mapcar (lambda (p)
			     (let ((d (compute-parameter distance-type 'd points p t)))
			       (cons (iter (for i from 0 below (length on-off))
					   (when (elt on-off i)
					     (sum (funcall blend-function d i))))
				     p)))
			   (vertices points))))
    (write-vtk-polylines
     (mapcar (lambda (lst) (mapcar (lambda (i) (elt vertices i)) lst))
	     (spider-uv-lines n))
     filename)))

#+nil
(write-blends-spider-polyline (points-from-angles '(40 20 60 100 80)) '(t nil nil nil nil)
			      "/tmp/blend.vtk" :blend-function #'ribbon-blend
			      :distance-type 'line-sweep)

(defun generate-colors (n)
  (assert (<= n 7) (n) "Only N <= 7 is supported.")
  (subseq '((255 0 0) (0 255 0) (0 0 255)
            (255 255 0) (255 0 255) (0 255 255)
            (0 0 0))
          0 n))

(defun write-color-blend-test (points filename r &key (blend-function #'ribbon-blend)
			       (distance-type 'perpendicular) (trim '(0.89d0 0.91d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
	   (lines (lines-from-points points))
	   (colors (generate-colors n)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (compute-parameter distance-type 'd points p t)) 
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
(defun write-si-line-test (points filename r &key (distance-type 'perpendicular)
			   (trim '(0.495d0 0.505d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
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
		(let* ((sp (compute-parameter distance-type 's points p t))
		       (on (and trim (some (lambda (x)
					     (and (< (first trim) x (second trim)) x))
					   sp))))
		  (if on
		      (format s "~{~d~^ ~}~%" (elt colors (position on sp)))
		      (format s "127 127 127~%"))))))))))))

#+nil
(write-color-blend-test (points-from-angles '(40 20 60 100 80)) "/tmp/blend.ppm" 200
			:blend-function #'ribbon-blend
			:distance-type 'chord-based
			:trim '(0.495d0 0.505d0))

#+nil
(write-si-line-test (points-from-angles '(40 20 60 100 80)) "/tmp/blend.ppm" 200
		    :distance-type 'line-sweep
		    :trim '(0.495d0 0.505d0))

(defun write-blend-voronoi (points filename r &key (blend-function #'ribbon-blend)
			    (distance-type 'perpendicular) (threshold 0.01d0))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
	   (lines (lines-from-points points)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (compute-parameter distance-type 'd points p t)) 
			 (blends (iter (for i from 0 below n)
				       (collect (funcall blend-function d i)))))
                    (if (destructuring-bind (a b)
			    (subseq (sort blends #'>) 0 2)
			  (< (- a b) threshold))
			(format s "0 0 0~%")
			(format s "127 127 127~%")))
                  (format s "255 255 255~%")))))))))

#+nil
(write-blend-voronoi (points-from-angles '(40 20 60 100 80)) "/tmp/blend2.ppm" 200
		     :blend-function #'ribbon-blend
		     :distance-type 'line-sweep
		     :threshold 0.03d0)

(defun floater-mean-value (points p)
  (flet ((angle (a b c)
	   (acos (scalar-product (vnormalize (v- a b))
				 (vnormalize (v- c b))))))
    (let* ((wi (iter (for vi in points)
		     (for vi-1 in (append (last points) points))
		     (for vi+1 in (append (rest points) points))
		     (for ri = (point-distance p vi))
		     (when (< (abs ri) *epsilon*)
		       (collect 1)
		       (next-iteration))
		     (for alphai = (angle vi p vi+1))
		     (for alphai-1 = (angle vi-1 p vi))
		     (collect (/ (+ (tan (/ alphai 2))
				    (tan (/ alphai-1 2)))
				 ri))))
	   (sum (reduce #'+ wi)))
      (mapcar (lambda (x) (/ x sum)) wi))))

(defun write-floater-blend-test (points filename r &key
				 (trim '(0.89d0 0.91d0)) (cornerp t))
  "TODO: How should side blends work in this case?"
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let ((wh (1+ (* 2 r)))
	  (lines (lines-from-points points))
	  (colors (generate-colors (length points))))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let ((blends (if cornerp
				    (floater-mean-value points p)
				    (let ((b (floater-mean-value points p)))
				      (mapcar (lambda (x y) (/ (+ x y) 2))
					      b (append (last b) b))))))
                    (if (and trim (some (lambda (x)
					  (< (first trim) x (second trim)))
					blends))
			(format s "0 0 0~%")
			(format s "~{~d~^ ~}~%"
				(mapcar #'round
					(reduce (lambda (x y) (mapcar #'+ x y))
						(mapcar #'v* colors blends))))))
                  (format s "255 255 255~%")))))))))

#+nil
(write-floater-blend-test (points-from-angles '(40 20 60 100 80))
			  "n-sided-paper/09-floater-corner.ppm" 400 :trim nil)
#+nil
(write-floater-blend-test (points-from-angles '(40 20 60 100 80))
			  "n-sided-paper/09-floater-corner-trim.ppm" 400)
#+nil
(write-floater-blend-test (points-from-angles '(40 20 60 100 80))
			  "n-sided-paper/09-floater-side.ppm" 400
			  :cornerp nil :trim nil)

(defun write-floater-blend-voronoi (points filename r &key
				    (threshold 0.01d0) (cornerp t))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let ((wh (1+ (* 2 r)))
	  (lines (lines-from-points points)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let ((blends (if cornerp
				    (floater-mean-value points p)
				    (let ((b (floater-mean-value points p)))
				      (mapcar (lambda (x y) (/ (+ x y) 2))
					      b (append (last b) b))))))
                    (if (destructuring-bind (a b)
			    (subseq (sort blends #'>) 0 2)
			  (< (- a b) threshold))
			(format s "0 0 0~%")
			(format s "127 127 127~%")))
                  (format s "255 255 255~%")))))))))

#+nil
(write-floater-blend-voronoi (points-from-angles '(40 20 60 100 80))
			     "n-sided-paper/09-floater-corner-voronoi.ppm" 400)

(defparameter *spider-density* 4)
(defparameter *spider-lines* 3)
(defun spider-lines (points)
  "Every polyline is computed according to *RESOLUTION*,
but the number of actual lines is also affected by two parameters:
*SPIDER-DENSITY* : draw only every ~th loop
*SPIDER-LINES*   : the number of `vertical' lines on one side."
  (let* ((lines (lines-from-points points))
	 (n (length lines))
	 (center (central-point points lines t))
	 (polylines (iter (for j from 0 to *resolution*)
			  (for coeff = (/ j *resolution*))
			  (for next = '())
			  (iter (for k from 0 below n)
				(iter (for i from 0 to *resolution*)
				      (for lp = (line-point (elt lines k) (/ i *resolution*)))
				      (push (affine-combine center coeff lp) next)))
			  (push (first next) next)
			  (collect (nreverse next)))))
    (append (iter (for line in (reverse polylines))
		  (for i upfrom 0)
		  (when (zerop (mod i *spider-density*))
		    (collect line)))
	    (iter (with gap = (floor (1+ *resolution*) (1+ *spider-lines*)))
		  (with indices = (cons 0 (iter (for j from 1 to *spider-lines*)
						   (collect (* j gap)))))
		  (for side from 0 below n)
		  (appending
		   (iter (for i in indices)
			 (for j = (+ (* side (1+ *resolution*)) i))
			 (collect (iter (for line in polylines)
					(collect (elt line j))))))))))
