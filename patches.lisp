(in-package :cl-nurbs-tests)

(defparameter *ribbon-multiplier* 1.0d0)
(defparameter *hermite* nil)

(defun delta (points i d)
  (if (= d 0)
      (elt points i)
      (v- (delta points (1+ i) (1- d))
	  (delta points i (1- d)))))

(defun factorial (n)
  (iter (for i from 1 to n)
	(multiply i)))

(defun combination (n k)
  (/ (factorial n) (factorial k) (factorial (- n k))))

(defun bernstein (n i u)
  (* (combination n i) (expt u i) (expt (- 1.0d0 u) (- n i))))

(defun bezier (points u &optional (derivative 0))
  (v* (iter (with p = '(0 0 0))
	    (for i from 0 to (- 3 derivative))
	    (setf p (v+ p (v* (delta points i derivative)
			      (bernstein (- 3 derivative) i u))))
	    (finally (return p)))
      (/ #.(factorial 3) (factorial (- 3 derivative)))))

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

(defun uniform-angles (n)
  (let ((angle (/ 360.0d0 n)))
    (loop repeat n collect angle)))

(defun uv-parameter (lines p)
  (let ((edge-length (* 2.0d0 (sin (/ pi (length lines))))))
    (mapcar (lambda (line)
	      (/ (point-line-distance p line)
		 edge-length))
	    lines)))

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

(defun corner-ribbon-evaluate (patch lines p)
  (iter (with n = (length lines))
	(with d = (uv-parameter lines p))
	(with result = '(0 0 0))
	(for i from 0 below n)
	(for di = (elt d i))
	(for dj = (elt d (mod (1- i) n)))
	(for blend = (corner-blend lines p i))
	(for point = 'todo)
	(setf result (v+ result (v* point blend)))
	(finally (return result))))

(defun double-corner-ribbon-evaluate (patch lines p)
  (iter (with n = (length lines))
	(with d = (uv-parameter lines p))
	(with result = '(0 0 0))
	(for i from 0 below n)
	(for di = (elt d i))
	(for dj = (elt d (mod (1- i) n)))
	(for blend = (/ (+ (corner-blend lines p (mod (1- i) n))
			   (corner-blend lines p i))
			2.0d0))
	(for q = (bezier (elt (elt patch i) 0) dj))
	(for point = (v+ q (v* (v- (bezier (elt (elt patch i) 1) dj) q)
			       3.0d0 di *ribbon-multiplier*)))
	(setf result (v+ result (v* point blend)))
	(finally (return result))))

(defun write-patch (patch filename &optional (fn #'coons-evaluate))
  (let* ((n (length patch))
	 (angles (uniform-angles n))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (parameters (vertices angles))
	 (vertices (mapcar (lambda (p) (funcall fn patch lines p))
			   parameters)))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

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
