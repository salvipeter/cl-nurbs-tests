(in-package :cl-nurbs-tests)

(defparameter *ribbon-multiplier* 1.0d0)
(defparameter *hermite* nil)

(defun delta (points i d)
  (if (= d 0)
      (elt points i)
      (v- (delta points (1+ i) (1- d))
	  (delta points i (1- d)))))

(defun factorial (n)
  (labels ((rec (k acc)
	     (if (<= k 1)
		 acc
		 (rec (1- k) (* k acc)))))
    (rec n 1)))

(defun combination (n k)
  (labels ((rec (n k acc)
	     (if (= k 0)
		 acc
		 (rec (1- n) (1- k) (/ (* n acc) k)))))
    (rec n k 1)))

(defun bernstein (n i u)
  (* (combination n i) (expt u i) (expt (- 1.0d0 u) (- n i))))

(defun bezier (points u &optional (derivative 0))
  (let ((n (1- (length points))))
    (v* (iter (with p = '(0 0 0))
	      (for i from 0 to (- n derivative))
	      (setf p (v+ p (v* (delta points i derivative)
				(bernstein (- n derivative) i u))))
	      (finally (return p)))
	(/ (factorial n) (factorial (- n derivative))))))

#+nil
(defun generate-patch (points twists)
  "POINTS : x4=p1 p2 p3 p4=q1 q2 q3 q4=r1 ... x3
TWISTS : x3p2 p3q2 q3r2 ... w3x2"
  (let ((n (length points)))
    (labels ((point (i) (elt points (mod i n)))
	     (twist (corner twist)
	       (v* (v- (v+ (point corner) twist)
		       (v+ (point (1- corner))
			   (point (1+ corner))))
		   9.0d0)))
      (iter (with n = (/ (length points) 3))
	    (for i from 0 below n)
	    (for j = (* i 3))
	    (collect (list (point (+ j 0))
			   (point (+ j 1))
			   (point (+ j 2))
			   (point (+ j 3)))
	      into curves)
	    (collect (list (point (- j 1))
			   (elt twists i)
			   (elt twists (mod (1+ i) n))
			   (point (+ j 4)))
	      into curves2)
	    (finally
	     (return
	       (iter (for k from 0 below n)
		     (for du = (bezier (elt curves k) 0 1))
		     (for dv = (v* (bezier (elt curves (mod (1- k) n)) 1 1) -1))
		     (collect
			 (list (elt curves k)
			       (elt curves2 k)
			       (bezier (elt curves k) 0 0)
			       (if (oddp k) dv du)
			       (if (oddp k) du dv)
			       (twist (* k 3) (elt twists k)))))))))))

#+nil
(defun uniform-angles (n)
  (let ((angle (/ 360.0d0 n)))
    (loop repeat n collect angle)))

#+nil
(defun uv-parameter (lines p)
  (let ((edge-length (* 2.0d0 (sin (/ pi (length lines))))))
    (mapcar (lambda (line)
	      (/ (point-line-distance p line)
		 edge-length))
	    lines)))

#+nil
(defun coons-evaluate (patch lines p)
  "PATCH is a list of (C1i C2i Ri Rui Rvi Ruvi ...)
For a 4-sided patch, D is (U V 1-U 1-V)"
  (let ((d (uv-parameter lines p)))
    (labels ((alpha0 (u) (if *hermite* (+ (* 2 u u u) (* -3 u u) 1) u))
	     (alpha1 (u) (if *hermite* (+ (* -2 u u u) (* 3 u u)) (- 1.0d0 u)))
	     (beta0 (u) (if *hermite* (+ (* u u u) (* -2 u u) u) 0.0d0))
	     (beta1 (u) (if *hermite* (- (* u u) (* u u u)) 0.0d0))
	     (rc (u v &optional deriv)
	       (let ((i (case deriv (u 3) (v 4) (uv 5) (t 2))))
		 (cond ((= u v 0) (elt (elt patch 0) i))
		       ((= u v 1) (elt (elt patch 2) i))
		       ((= u 0) (elt (elt patch 3) i))
		       (t (elt (elt patch 1) i)))))
	     (r (p u deriv)
	       (if deriv
		   (v* (v- (bezier (elt (elt patch p) 1) u)
			   (bezier (elt (elt patch p) 0) u))
		       3.0d0)
		   (bezier (elt (elt patch p) 0) u)))
	     (rv (u v &optional deriv)
	       (if (= u 0)
		   (r 3 v deriv)
		   (r 1 (- 1.0d0 v) deriv)))
	     (ru (u v &optional deriv)
	       (if (= v 0)
		   (r 0 u deriv)
		   (r 2 (- 1.0d0 u) deriv))))
      (let ((u (second d))
	    (v (first d)))
	(v- (v+ (v* (rv 0 v) (alpha0 u))
		(v* (rv 1 v) (alpha1 u))
		(v* (rv 0 v 'u) (beta0 u))
		(v* (rv 1 v 'u) (beta1 u))
		(v* (ru u 0) (alpha0 v))
		(v* (ru u 1) (alpha1 v))
		(v* (ru u 0 'v) (beta0 v))
		(v* (ru u 1 'v) (beta1 v)))
	    (v+ (v* (rc 0 0) (alpha0 u) (alpha0 v))
		(v* (rc 0 1) (alpha0 u) (alpha1 v))
		(v* (rc 0 0 'v) (alpha0 u) (beta0 v))
		(v* (rc 0 1 'v) (alpha0 u) (beta1 v))
		(v* (rc 1 0) (alpha1 u) (alpha0 v))
		(v* (rc 1 1) (alpha1 u) (alpha1 v))
		(v* (rc 1 0 'v) (alpha1 u) (beta0 v))
		(v* (rc 1 1 'v) (alpha1 u) (beta1 v))
		(v* (rc 0 0 'u) (beta0 u) (alpha0 v))
		(v* (rc 0 1 'u) (beta0 u) (alpha1 v))
		(v* (rc 0 0 'uv) (beta0 u) (beta0 v))
		(v* (rc 0 1 'uv) (beta0 u) (beta1 v))
		(v* (rc 1 0 'u) (beta1 u) (alpha0 v))
		(v* (rc 1 1 'u) (beta1 u) (alpha1 v))
		(v* (rc 1 0 'uv) (beta1 u) (beta0 v))
		(v* (rc 1 1 'uv) (beta1 u) (beta1 v))))))))

#+nil
(defun classic-ribbon-evaluate (patch lines p)
  (iter (with n = (length lines))
	(with d = (uv-parameter lines p))
	(with result = '(0 0 0))
	(for i from 0 below n)
	(for di = (elt d i))
	(for dj = (elt d (mod (1- i) n)))
	(for blend = (ribbon-blend lines p i))
	(for q = (bezier (elt (elt patch i) 0) dj))
	(for point = (v+ q (v* (v- (bezier (elt (elt patch i) 1) dj) q)
			       3.0d0 di *ribbon-multiplier*)))
	(setf result (v+ result (v* point blend)))
	(finally (return result))))

#+nil
(defun corner-correction (patch d i)
  (let* ((i-1 (mod (1- i) (length d)))
	 (di-1 (elt d i-1))
	 (di (elt d i)))
    (v+ (elt (elt patch i) 2)
	(v* (elt (elt patch i) (if (oddp i) 4 3)) di-1 *ribbon-multiplier*)
	(v* (elt (elt patch i) (if (oddp i) 3 4)) di *ribbon-multiplier*)
	(v* (elt (elt patch i) 5) di-1 di
	    *ribbon-multiplier* *ribbon-multiplier*))))

#+nil
(defun corner-ribbon-evaluate (patch lines p)
  (iter (with n = (length lines))
	(with d = (uv-parameter lines p))
	(with result = '(0 0 0))
	(for i from 0 below n)
	(for i-1 = (mod (1- i) n))
	(for di = (elt d i))
	(for dj = (elt d i-1))
	(for blend = (corner-blend lines p i-1))
	(for q1 = (bezier (reverse (elt (elt patch i-1) 0)) di))
	(for q2 = (bezier (elt (elt patch i) 0) dj))
	(for point = (v+ q1 (v* (v- (bezier (reverse (elt (elt patch i-1) 1)) di) q1)
				3.0d0 dj *ribbon-multiplier*)
			 q2 (v* (v- (bezier (elt (elt patch i) 1) dj) q2)
				3.0d0 di *ribbon-multiplier*)))
	(for correction = (corner-correction patch d i))
	(setf result (v+ result (v* (v- point correction) blend)))
	(finally (return result))))

#+nil
(defun double-corner-ribbon-evaluate (patch lines p)
  (iter (with n = (length lines))
	(with d = (uv-parameter lines p))
	(with result = '(0 0 0))
	(for i from 0 below n)
	(for di = (elt d i))
	(for dj = (elt d (mod (1- i) n)))
	(for blend = (+ (corner-blend lines p (mod (1- i) n))
			(corner-blend lines p i)))
	(for q = (bezier (elt (elt patch i) 0) dj))
	(for point = (v+ q (v* (v- (bezier (elt (elt patch i) 1) dj) q)
			       3.0d0 di *ribbon-multiplier*)))
	(for correction = (corner-correction patch d i))
	(for correction-blend = (corner-blend lines p (mod (1- i) n)))
	(setf result (v+ result (v- (v* point blend)
				    (v* correction correction-blend))))
	(finally (return result))))

#+nil
(defun write-patch (patch filename &optional (fn #'coons-evaluate))
  (let* ((n (length patch))
	 (angles (uniform-angles n))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (parameters (vertices angles))
	 (vertices (mapcar (lambda (p) (funcall fn patch lines p))
			   parameters)))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

#+nil
(defun write-square-patch (patch filename &optional (fn #'coons-evaluate))
  (assert (= (length patch) 4))
  (let* ((points (points-from-angles '(45 90 90 90)))
	 (lines (lines-from-points points)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "~d ~d~%" *resolution* *resolution*)
      (iter (for ui from 0 below *resolution*)
	    (for u = (* (sqrt 2.0d0) (- (/ ui (1- *resolution*)) 0.5d0)))
	    (iter (for vi from 0 below *resolution*)
		  (for v = (* (sqrt 2.0d0) (- (/ vi (1- *resolution*)) 0.5d0)))
		  (format s "~{~f~^ ~}~%"
			  (funcall fn patch lines (list u v))))))))

#+nil
(defun write-bezier-polygons (curves filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let ((n (* (length curves) 4)))
      (format s "# vtk DataFile Version 1.0~%~
               Bezier control polygons~%~
               ASCII~%~
               DATASET POLYDATA~%~%~
               POINTS ~d float~%" n)
      (iter (for p in (reduce #'append curves))
	    (format s "~{~f~^ ~}~%" p))
      (format s "~%LINES 1 ~d~%~d " (+ n 2) (1+ n))
      (iter (for i from 0 below n)
	    (format s "~d " i))
      (format s "0~%"))))

#+nil
(let ((*hermite* nil)
      (*ribbon-multiplier* 0.5d0)
      (*exponent* 2)
      (*resolution* 40))
  (write-patch
   (generate-patch
    '((0 0 0) (2 0 1) (4 0 1)
      (6 0 0) (6 2 2) (7 4 1)
      (8 6 0) (6 7 1) (3 7 1)
      (1 6 0) (0 4 1) (0 2 2))
    '((2 2 1) (4 2 1) (5 5 2) (2 5 2)))
   "/tmp/coons.vtk" #'coons-evaluate))

#+nil
(let ((*hermite* nil)
      (*ribbon-multiplier* 0.5d0)
      (*exponent* 2)
      (*resolution* 40))
  (write-patch
   (generate-patch
    '((0 0 0) (2 0 0) (4 0 0)
      (6 0 0) (6 2 0) (6 4 0)
      (6 6 0) (4 6 0) (2 6 0)
      (0 6 0) (0 4 0) (0 2 0))
    '((2 2 1) (4 2 1) (4 4 1) (2 4 1)))
   "/tmp/coons.vtk" #'coons-evaluate))


(defun generate-coordinates (lines heights)
  (iter (for line in lines)
	(for lst in heights)
	(for n = (length lst))
	(collect
	 (iter (for i from 0 below n)
	       (for x in lst)
	       (for yz = (affine-combine (first line) (/ i (1- n)) (second line)))
	       (collect (cons x yz))))))

#+nil
(generate-coordinates (lines-from-points (points-from-angles '(30 120 100)))
		      '((1 2 1) (1 3) (3 2 2 1)))

(defun ribbon-evaluate (patch i s d)
  "PATCH has two elements: (outer-bezier-points inner-bezier-points)."
  (let* ((base-point (bezier (elt (first patch) i) (elt s i)))
	 (inner-point (bezier (elt (second patch) i) (elt s i)))
	 (derivative (v* (v- inner-point base-point) 3.0d0)))
    (v+ base-point (v* derivative (elt d i) *ribbon-multiplier*))))

(defun corner-evaluate (patch i s)
  (let* ((i-1 (mod (1- i) (length s)))
	 (si-1 (elt s i-1))
	 (si (elt s i))
	 (ci-1 (bezier (elt (first patch) i-1) si-1))
	 (ci (bezier (elt (first patch) i) si))
	 (di-1 (v* (v- (bezier (elt (second patch) i-1) si-1) ci-1) 3.0d0
		   *ribbon-multiplier*))
	 (di (v* (v- (bezier (elt (second patch) i) si) ci) 3.0d0
		 *ribbon-multiplier*)))
    (v+ ci ci-1 (v* di (- 1.0d0 si-1)) (v* di-1 si))))

(defun corner-correction (patch i s)
  (let* ((i-1 (mod (1- i) (length s)))
	 (si-1 (elt s i-1))
	 (si (elt s i))
	 (previous (let ((lst (elt (first patch) i-1)))
		     (elt lst (- (length lst) 2))))
	 (corner (first (elt (first patch) i)))
	 (next (second (elt (first patch) i)))
	 (twist (second (elt (second patch) i))))
    (v+ corner
	(v* (v- previous corner) 3.0d0 (- 1.0d0 si-1) *ribbon-multiplier*)
	(v* (v- next corner) 3.0d0 si *ribbon-multiplier*)
	(v* (v- (v+ corner twist) (v+ previous next)) 9.0d0
	    (- 1.0d0 si-1) si *ribbon-multiplier* *ribbon-multiplier*))))

(defun compute-parameter (type points p &optional no-tiny-p)
  (macrolet ((tiny-lambda ((args) &body body)
	       `(lambda (,args)
		  (if no-tiny-p
		      (let ((result (progn ,@body)))
			(if (< result *tiny*) 0.0d0 result))
		      (progn ,@body)))))
    (mapcar (ecase type
	      (perpendicular (tiny-lambda (lst) (perpendicular-distance points lst p 's)))
	      (barycentric (let ((lines (lines-from-points points)))
			     (tiny-lambda (lst) (barycentric-distance lines lst p 's))))
	      (chord-based (tiny-lambda (lst) (chord-based-distance points lst p 's)))
	      (radial (tiny-lambda (lst) (radial-distance points lst p 's)))
	      (line-sweep (let ((center (central-point points (lines-from-points points) t)))
			    (tiny-lambda (lst) (line-sweep-distance center lst p 's)))))
	    (iter (for i from -2 below (- (length points) 2))
		  (collect (iter (for j from 0 below 4)
				 (collect (elt points (mod (+ i j) (length points))))))))))

(defun generate-patch (outer inner)
  (let ((n (length outer)))
    (list outer
	  (iter (for i from 0 below n)
		(for previous = (elt outer (mod (1- i) n)))
		(for next = (elt outer (mod (1+ i) n)))
		(collect (append (list (elt previous (- (length previous) 2)))
				 (elt inner i)
				 (list (elt next 1))))))))

(defun write-patch (angles type filename &key heights coords (distance-type 'perpendicular))
  (let* ((n (length angles))
	 (points (points-from-angles angles))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch
		     (generate-coordinates (lines-from-points points) (first heights))
		     (generate-coordinates (lines-from-points (points-from-angles angles 0.5d0))
					   (second heights)))))
	 (vertices (iter (for domain-point in (vertices points))
			 (for p = (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0))
					   domain-point))
			 (for d = (compute-distance distance-type points p t))
			 (for s = (compute-parameter distance-type points p t))
			 (collect
			  (iter (for i from 0 below (length points))
				(with result = '(0 0 0))
				(setf result
				      (v+ result
					  (case type
					    (ribbon (v* (ribbon-evaluate patch i s d)
							(ribbon-blend d i)))
					    (corner (v* (v- (corner-evaluate patch i s)
							    (corner-correction patch i s))
							(corner-blend d (mod (1- i) n))))
					    (hybrid (v- (v* (ribbon-evaluate patch i s d)
							    (+ (corner-blend d (mod (1- i) n))
							       (corner-blend d i)))
							(v* (corner-correction patch i s)
							    (corner-blend d (mod (1- i) n))))))))
				(finally (return result)))))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

#+nil
(let ((*ribbon-multiplier* 1.0d0))
  (write-patch '(40 20 60 100 80)
	       'corner
	       "/tmp/patch.vtk"
	       :heights
	       '(((0.0d0 0.1d0 0.1d0 0.0d0)
		  (0.0d0 0.2d0 0.3d0 0.4d0)
		  (0.4d0 0.6d0 0.6d0 0.4d0)
		  (0.4d0 0.5d0 0.6d0 0.4d0 0.2d0 0.0d0)
		  (0.0d0 0.2d0 0.1d0 0.0d0))
		 ((0.2 0.2)
		  (0.2 0.5)
		  (0.5 0.8)
		  (0.8 0.2)
		  (0.2 0.2)))
	       :distance-type 'line-sweep))

#+nil
(let ((*ribbon-multiplier* 0.5d0)
      (*resolution* 20))
  (write-patch '(50 20 50 240)
	       'ribbon
	       "/tmp/patch.vtk"
	       :heights
	       '(((0 0.2 0.2 0)
		  (0 0.4 0.4 0)
		  (0 0.2 0.2 0)
		  (0 0.4 0.4 0))
		 ((0.6 0.6)
		  (0.6 0.6)
		  (0.6 0.6)
		  (0.6 0.6)))
	       :distance-type 'line-sweep))

;;; problemak:
;;; - mi lesz a coons patch kiertekelessel?

(defun write-vtk-curves (curves filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let* ((lengths (mapcar #'length curves))
	   (n (reduce #'+ lengths)))
      (format s "# vtk DataFile Version 1.0~
		 ~%Control Curves~
		 ~%ASCII~
		 ~%DATASET POLYDATA~%~
		 ~%POINTS ~d float~%" n)
      (iter (for point in (reduce #'append curves))
	    (format s "~{~f ~}~%" point))
      (let ((nlines (- n (length curves))))
	(format s "~%LINES ~d ~d~%" nlines (* 3 nlines)))
      (iter (for k first 0 then (+ k d))
	    (for d in lengths)
	    (iter (for i from k below (+ k d -1))
		  (format s "2 ~d ~d~%" i (1+ i)))))))

(defun write-constraint-ribbons (angles filename &key (resolution 100) heights coords)
  (let* ((points (points-from-angles angles))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch
		     (generate-coordinates (lines-from-points points) (first heights))
		     (generate-coordinates (lines-from-points (points-from-angles angles 0.5d0))
					   (second heights)))))
	 (curves (iter (for curve1 in (first patch))
		       (for curve2 in (second patch))
		       (for points1 =
			    (iter (for u from 0 to 1 by (/ resolution))
				  (collect (bezier curve1 u))))
		       (for points2 =
			    (iter (for u from 0 to 1 by (/ resolution))
				  (collect (bezier curve2 u))))
		       (appending (append (list points1 points2)
					  (mapcar #'list points1 points2))))))
    (write-vtk-curves curves filename)))

(defun write-constraint-grid (angles filename &key heights coords)
  (let* ((points (points-from-angles angles))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch
		     (generate-coordinates (lines-from-points points) (first heights))
		     (generate-coordinates (lines-from-points (points-from-angles angles 0.5d0))
					   (second heights))))))
    (write-vtk-curves (append (first patch) (second patch)) filename)))

#+nil
(write-constraint-grid
 '(50 20 50 240)
 "/tmp/grid.vtk"
 :heights
 '(((0 0.2 0.2 0)
    (0 0.4 0.4 0)
    (0 0.2 0.2 0)
    (0 0.4 0.4 0))
   ((0.6 0.6)
    (0.6 0.6)
    (0.6 0.6)
    (0.6 0.6))))

#+nil
(write-constraint-ribbons
 '(90 90 90 90)
 "/tmp/grid.vtk"
 :heights
 '(((0 0.2 0.2 0)
    (0 0.4 0.4 0)
    (0 0.2 0.2 0)
    (0 0.4 0.4 0))
   ((0.6 0.6)
    (0.6 0.6)
    (0.6 0.6)
    (0.6 0.6))))
