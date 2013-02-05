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

(defun bezier-arc-length (points &key (from 0.0d0) (to 1.0d0))
  (iter (for gauss in cl-nurbs::+gaussian-quadrature+)
	(for u = (/ (+ (* (- to from) (first gauss)) from to) 2.0))
	(for dn = (vlength (bezier points u 1)))
	(sum (* dn (second gauss) (- to from) 0.5))))

(defun generate-coons-patch (outer inner)
  (assert (= (length outer) 4))
  (destructuring-bind (curves curves2)
      (generate-patch outer inner)
    (let ((n (1- (length (first outer)))))
      (labels ((twist (i)
		 (let* ((corner (first (elt curves i)))
			(twist (second (elt curves2 i)))
			(prev-curve (elt curves (mod (1- i) 4)))
			(prev (elt prev-curve (- (length prev-curve) 2)))
			(next (second (elt curves i))))
		   (v* (v- (v+ corner twist) (v+ prev next)) n n))))
	(iter (for k from 0 below 4)
	      (for du = (bezier (elt curves k) 0 1))
	      (for dv = (v* (bezier (elt curves (mod (1- k) 4)) 1 1) -1))
	      (collect
		  (list (elt curves k)
			(elt curves2 k)
			(bezier (elt curves k) 0 0)
			(if (oddp k) dv du)
			(if (oddp k) du dv)
			(twist k))))))))

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
  (let ((d (uv-parameter lines p))
	(n (1- (length (first (first patch))))))
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
		       n)
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

(defun write-coons-patch (coords filename &optional (fn #'coons-evaluate))
  (let* ((patch (generate-coons-patch (first coords) (second coords)))
	 (angles (uniform-angles 4))
	 (points (points-from-angles angles))
	 (lines (lines-from-points points))
	 (parameters (vertices points))
	 (vertices (mapcar (lambda (p) (funcall fn patch lines p))
			   parameters)))
    (write-vtk-indexed-mesh vertices (triangles 4) filename)))

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
    (v+ base-point (v* derivative (gamma (elt d i)) *ribbon-multiplier*))))

(defun corner-evaluate (patch i s &optional d)
  "The corner defined by segments I-1 and I,
thus containing point I-1 (NOT point I)."
  (let* ((i-1 (mod (1- i) (length s)))
	 (si-1 (if d (- 1.0d0 (elt d i)) (elt s i-1)))
	 (si (elt s i))
	 (ci-1 (bezier (elt (first patch) i-1) si-1))
	 (ci (bezier (elt (first patch) i) si))
	 (di-1 (v* (v- (bezier (elt (second patch) i-1) si-1) ci-1) 3.0d0))
	 (di (v* (v- (bezier (elt (second patch) i) si) ci) 3.0d0)))
    (v+ ci ci-1 (v* di (- 1.0d0 si-1)) (v* di-1 si))))

(defun corner-correction (patch i s &optional d)
  "Correction for the same corner as explained in CORNER-EVALUATE."
  (let* ((i-1 (mod (1- i) (length s)))
	 (si-1 (if d (elt d i) (- 1.0d0 (elt s i-1))))
	 (si (elt s i))
	 (previous (let ((lst (elt (first patch) i-1)))
		     (elt lst (- (length lst) 2))))
	 (corner (first (elt (first patch) i)))
	 (next (second (elt (first patch) i)))
	 (twist (second (elt (second patch) i))))
    (v+ corner
	(v* (v- previous corner) 3.0d0 (gamma si-1))
	(v* (v- next corner) 3.0d0 (gamma si))
	(v* (v- (v+ corner twist) (v+ previous next)) 9.0d0 (gamma si-1) (gamma si)))))

(defparameter *use-gamma* nil)
(defun gamma (d)
  (if *use-gamma*
      (/ d (1+ (* 2 d)))
      d))

(defun coons-ribbon-evaluate (patch i s d)
  (let ((n (length (first patch)))
	(si (if (atom s) s (elt s i)))
	(di (if (atom d) d (elt d i))))
    (flet ((linear-ribbon (i si di)
	     (let* ((base-point (bezier (elt (first patch) i) si))
		    (inner-point (bezier (elt (second patch) i) si))
		    (derivative (v* (v- inner-point base-point) 3.0d0)))
	       (v+ base-point (v* derivative (gamma di)))))
	   (correction (i si si-1)
	     (let* ((i-1 (mod (1- i) n))
		    (prev (let ((lst (elt (first patch) i-1)))
			    (elt lst (- (length lst) 2))))
		    (corner (first (elt (first patch) i)))
		    (next (second (elt (first patch) i)))
		    (twist (second (elt (second patch) i))))
	       (v+ corner
		   (v* (v- prev corner) 3.0d0 (gamma si-1))
		   (v* (v- next corner) 3.0d0 (gamma si))
		   (v* (v- (v+ corner twist) (v+ prev next)) 9.0d0 (gamma si-1) (gamma si))))))
      (let ((hs (hermite-blend-function 'point 'start si))
	    (hd (hermite-blend-function 'point 'start di)))
	(v- (v+ (v* (linear-ribbon (mod (1- i) n) (- 1.0d0 di) si) hs)  
		(v* (linear-ribbon i si di) hd)
		(v* (linear-ribbon (mod (1+ i) n) di (- 1.0d0 si)) (- 1.0d0 hs)))
	    (v+ (v* (correction i si di) hs hd)
		(v* (correction (mod (1+ i) n) di (- 1.0d0 si)) (- 1.0d0 hs) hd)))))))

(defun compute-parameter (type dir points p &optional no-tiny-p)
  (macrolet ((tiny-lambda ((args) &body body)
	       `(lambda (,args)
		  (if no-tiny-p
		      (let ((result (progn ,@body)))
			(if (< result *tiny*) 0.0d0 result))
		      (progn ,@body)))))
    (mapcar (tiny-lambda (lst) (compute-distance type points lst p dir))
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

(defun patch-evaluate (patch points type distance-type domain-point)
  (let* ((n (length points))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
	 (d (unless (and (eq type 'corner) (eq distance-type 'biquadratic))
	      (compute-parameter distance-type 'd points p t)))
	 (dc (if (and (member type '(corner hybrid)) (eq distance-type 'biquadratic))
		 (compute-parameter 'biquadratic-corner 'd points p t)
		 d))
	 (s (unless (and (eq type 'corner) (eq distance-type 'biquadratic))
	      (compute-parameter distance-type 's points p t)))
	 (sc (if (and (member type '(corner hybrid)) (eq distance-type 'biquadratic))
		 (compute-parameter 'biquadratic-corner 's points p t)
		 s))
	 (b (and (eq type 'sketches) (compute-parameter 'perpendicular 'd points p t)))
         (*use-gamma* (member type '(hybrid hybrid-coons))))
    (iter (for i from 0 below n)
	  (with result = '(0 0 0))
	  (setf result
		(v+ result
		    (case type
		      (ribbon (v* (ribbon-evaluate patch i s d)
				  (ribbon-blend d i)))
		      (ribbon-coons (v* (coons-ribbon-evaluate patch i s d)
					(ribbon-blend d i)))
		      (corner (if (eq distance-type 'biquadratic)
				  (v* (v- (corner-evaluate patch i sc dc)
					  (corner-correction patch i sc dc))
				      (corner-blend dc (mod (1- i) n)))
				  (v* (v- (corner-evaluate patch i s)
					  (corner-correction patch i s))
				      (corner-blend d (mod (1- i) n)))))
		      (hybrid (v- (v* (ribbon-evaluate patch i s d)
				      (+ (corner-blend d (mod (1- i) n))
					 (corner-blend d i)))
				  (if (eq distance-type 'biquadratic)
				      (v* (corner-correction patch i sc dc)
					  (corner-blend d (mod (1- i) n)))
				      (v* (corner-correction patch i s)
					  (corner-blend d (mod (1- i) n))))))
		      (hybrid-coons (v* (coons-ribbon-evaluate patch i s d)
					;; (+ (corner-blend d (mod (1- i) n))
					;;    (corner-blend d i))
                                        ;; (mean-side-blend points (mod (1- i) n) i p)
					(+ (wachspress-corner-blend points (mod (1- i) n) p)
					   (wachspress-corner-blend points i p))
					;; (+ (mean-corner-blend points (mod (1- i) n) p)
					;;    (mean-corner-blend points i p))
					1/2))
		      (sketches (v* (ribbon-evaluate patch i s d)
				    (ribbon-blend b i)))
		      (sketches-coons (v* (coons-ribbon-evaluate patch i s d)
					  (ribbon-blend b i))))))
	  (finally (return result)))))

(defun generate-patch-from-heights (points inner-points heights)
  (generate-patch
   (generate-coordinates (lines-from-points points) (first heights))
   (generate-coordinates (lines-from-points inner-points) (second heights))))

(defun write-patch (points type filename &key inner-points heights coords
		    (distance-type 'perpendicular) spider)
  (let* ((n (length points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (if spider
	(write-vtk-polylines
	 (iter (for line in (spider-lines points))
	       (collect (iter (for domain-point in line)
			      (collect (patch-evaluate patch points type distance-type
						       domain-point)))))
	 filename)
	(write-vtk-indexed-mesh
	 (iter (for domain-point in (vertices points))
	       (collect (patch-evaluate patch points type distance-type domain-point)))
	 (triangles n) filename))))

#+nil
(let ((*ribbon-multiplier* 1.0d0))
  (write-patch (points-from-angles '(40 20 60 100 80))
	       'corner
	       "/tmp/patch.vtk"
	       :inner-points (points-from-angles '(40 20 60 100 80) 0.5d0)
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
  (write-patch (points-from-angles '(50 20 50 240))
	       'ribbon
	       "/tmp/patch.vtk"
	       :inner-points (points-from-angles '(50 20 50 240) 0.5d0)
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

(defun write-constraint-ribbons (points filename &key (resolution 100) inner-points heights coords)
  (let* ((patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights)))
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

(defun write-constraint-coons-ribbons (points filename-base &key
				       (resolution 100) inner-points heights coords)
  (let* ((n (length points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (iter (for i from 0 below n)
	  (for s-curves =
	       (iter (for s from 0 to 1 by (/ resolution))
		     (collect
			 (iter (for d from 0 to 1 by (/ resolution))
			       (collect (coons-ribbon-evaluate patch i s d))))))
	  (for d-curves =
	       (iter (for d from 0 to 1 by (/ resolution))
		     (collect
			 (iter (for s from 0 to 1 by (/ resolution))
			       (collect (coons-ribbon-evaluate patch i s d))))))
	  (write-vtk-curves (append s-curves d-curves)
			    (format nil "~f-~d.vtk" filename-base (1+ i))))))

(defun write-constraint-grid (points filename &key inner-points heights coords)
  (let* ((patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (write-vtk-curves (append (first patch) (second patch)) filename)))

#+nil
(write-constraint-grid
 (points-from-angles '(50 20 50 240))
 "/tmp/grid.vtk"
 :inner-points (points-from-angles '(50 20 50 240) 0.5d0)
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
 (points-from-angles '(90 90 90 90))
 "/tmp/grid.vtk"
 :inner-points (points-from-angles '(90 90 90 90) 0.5d0)
 :heights
 '(((0 0.2 0.2 0)
    (0 0.4 0.4 0)
    (0 0.2 0.2 0)
    (0 0.4 0.4 0))
   ((0.6 0.6)
    (0.6 0.6)
    (0.6 0.6)
    (0.6 0.6))))

;;; Coons patch example

#+nil
(let ((*hermite* nil)
      (*resolution* 80))
  (write-coons-patch
   (generate-coons-patch
    '((0 0 0) (2 0 1) (4 0 1)
      (6 0 0) (6 2 2) (7 4 1)
      (8 6 0) (6 7 1) (3 7 1)
      (1 6 0) (0 4 1) (0 2 2))
    '((2 2 2) (4 2 2) (5 5 2) (2 5 2)))
   "n-sided-paper/coons.vtk" #'coons-evaluate))

#+nil
(write-constraint-grid
 (points-from-angles '(50 20 50 240))
 "n-sided-paper/coons-grid.vtk"
 :coords '((((0 0 0) (2 0 1) (4 0 1) (6 0 0))
	    ((6 0 0) (6 2 2) (7 4 1) (8 6 0))
	    ((8 6 0) (6 7 1) (3 7 1) (1 6 0))
	    ((1 6 0) (0 4 1) (0 2 2) (0 0 0)))
	   (((2 2 2) (4 2 2))
	    ((4 2 2) (5 5 2))
	    ((5 5 2) (2 5 2))
	    ((2 5 2) (2 2 2)))))

#+nil
(write-constraint-ribbons
 (points-from-angles '(50 20 50 240))
 "n-sided-paper/coons-ribbons.vtk"
 :coords '((((0 0 0) (2 0 1) (4 0 1) (6 0 0))
	    ((6 0 0) (6 2 2) (7 4 1) (8 6 0))
	    ((8 6 0) (6 7 1) (3 7 1) (1 6 0))
	    ((1 6 0) (0 4 1) (0 2 2) (0 0 0)))
	   (((2 2 2) (4 2 2))
	    ((4 2 2) (5 5 2))
	    ((5 5 2) (2 5 2))
	    ((2 5 2) (2 2 2))))
 :resolution 20)

(defun bezier-projection-starting-value (curve point res)
  (iter (for i from 0 below res)
	(for u = (/ i (1- res)))
	(finding u minimizing (vlength (v- point (bezier curve u))))))

(defun bezier-project-point (curve point iterations search-resolution &optional
			     (distance-tolerance 0.0) (cosine-tolerance 0.0))
  "Returns the parameter of CURVE's closest point to POINT.

The function uses the Newton-Raphson method. (NURBS Book, pp. 230-232)
SEARCH-RESOLUTION parameters are checked for a suitable initial value."
  (let ((u0 (bezier-projection-starting-value curve point search-resolution)))
    (iter (repeat iterations)
	  (for last first nil then u)
	  (for u first u0 then
	       (min (max (- u (/ (scalar-product d deviation)
				 (+ (scalar-product d2 deviation)
				    (scalar-product d d))))
			 0)
		    1))
	  (for p = (bezier curve u))
	  (for d = (bezier curve u 1))
	  (for d2 = (bezier curve u 2))
	  (for deviation = (v- p point))
	  (when (or (<= (vlength deviation) distance-tolerance)
		    (<= (/ (abs (scalar-product d deviation))
			   (* (vlength d) (vlength deviation)))
			cosine-tolerance)
		    (and last
			 (<= (vlength (v* d (- u last))) distance-tolerance)))
	    (leave (values u (vlength deviation))))
	  (finally (return (values u (vlength deviation)))))))

(defun check-patch (points type &key inner-points heights coords (distance-type 'perpendicular))
  "Only side interpolation checking for now."
  (let* ((n (length points))
	 (lines (lines-from-points points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (iter (for side from 0 below n)
	  (for curve in (first patch))
	  (for line in lines)
	  (maximize
	   (iter (for i from 0 below *resolution*)
		 (for u = (/ i (1- *resolution*)))
		 (for domain-point = (line-point line u))
		 (for p = (patch-evaluate patch points type distance-type domain-point))
		 (for v = (bezier-project-point curve p 10 *resolution*))
		 (for diff = (point-distance p (bezier curve v)))
		 (maximize diff))))))

(defun check-patch-normals (points type &key inner-points heights coords
			     (distance-type 'perpendicular) (step 0.1d-3))
  "Only side interpolation checking for now."
  (let* ((n (length points))
	 (lines (lines-from-points points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (iter (for side from 0 below n)
	  (for curve in (first patch))
	  (for inner-curve in (second patch))
	  (for line in lines)
	  (for direction = (v- (second line) (first line)))
	  (for normal = (vnormalize (list (- (second direction)) (first direction))))
	  (with center = (central-point points lines t))
	  (when (< (scalar-product (v- center (first line)) normal) 0)
	    (setf normal (v* normal -1)))
	  (maximize
	   (iter (for i from 0 below *resolution*)
		 (for u = (/ i (1- *resolution*)))
		 (for domain-point = (line-point line u))
		 (for inner-domain-point = (v+ domain-point (v* normal step)))
		 (for p = (patch-evaluate patch points type distance-type domain-point))
		 (for q = (patch-evaluate patch points type distance-type inner-domain-point))
		 (for v = (bezier-project-point curve p 10 *resolution*))
		 (for derivative = (vnormalize (bezier curve v 1)))
		 (for diff =
		      (vlength
		       (cross-product
			(vnormalize (cross-product derivative (v- q p)))
			(vnormalize (cross-product derivative
						   (v- (bezier inner-curve v)
						       (bezier curve v)))))))
		 (maximize diff))))))

(defun check-patch-tangents (points type &key inner-points heights coords
			     (distance-type 'perpendicular) (step 0.1d-3))
  "Only side interpolation checking for now."
  (let* ((n (length points))
	 (lines (lines-from-points points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (iter (for side from 0 below n)
	  (for curve in (first patch))
	  (for inner-curve in (second patch))
	  (for line in lines)
	  (for direction = (v- (second line) (first line)))
	  (for normal = (vnormalize (list (- (second direction)) (first direction))))
	  (with center = (central-point points lines t))
	  (when (< (scalar-product (v- center (first line)) normal) 0)
	    (setf normal (v* normal -1)))
	  (maximize
	   (iter (for i from 0 below *resolution*)
		 (for u = (/ i (1- *resolution*)))
		 (for domain-point = (line-point line u))
		 (for inner-domain-point = (v+ domain-point (v* normal step)))
		 (for p = (patch-evaluate patch points type distance-type domain-point))
		 (for q = (patch-evaluate patch points type distance-type inner-domain-point))
		 (for v = (bezier-project-point curve p 10 *resolution*))
		 (for diff =
		      (vlength
		       (cross-product
			(vnormalize (v- q p))
			(vnormalize (v- (bezier inner-curve v)
					(bezier curve v))))))
		 (maximize diff))))))

#+nil
(let ((*resolution* 100)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* t))
  (format t "Centralized line sweep is ~:[OFF~;ON~]~%" *centralized-line-sweep*)
  (iter (for type in '(ribbon corner hybrid sketches))
	(iter (for distance in (if (eq type 'sketches)
				   '(perpendicular)
				   '(perpendicular barycentric radial chord-based line-sweep)))
	      (format t "Maximal  point  deviation using ~a with ~a: ~f~%"
		      type distance
		      (check-patch *points* type :coords *coords* :distance-type distance))
	      (format t "Maximal tangent deviation using ~a with ~a: ~f~%"
		      type distance
		      (check-patch-tangents *points* type :coords *coords* :distance-type distance
					    :step 1.0d-4)))))

;; Centralized line sweep is ON
;; Maximal  point  deviation using RIBBON with PERPENDICULAR: 0.000000000000002737549743767154
;; Maximal tangent deviation using RIBBON with PERPENDICULAR: 0.0006094788868459588
;; Maximal  point  deviation using RIBBON with BARYCENTRIC: 0.000000000000002817430187029739
;; Maximal tangent deviation using RIBBON with BARYCENTRIC: 0.001242037224292269
;; Maximal  point  deviation using RIBBON with RADIAL: 0.0000000000000027577367830154345
;; Maximal tangent deviation using RIBBON with RADIAL: 0.0008624970038566559
;; Maximal  point  deviation using RIBBON with CHORD-BASED: 0.000000000000003575193960459936
;; Maximal tangent deviation using RIBBON with CHORD-BASED: 0.9805759386872815
;; Maximal  point  deviation using RIBBON with LINE-SWEEP: 0.0000000000000026852712547870858
;; Maximal tangent deviation using RIBBON with LINE-SWEEP: 0.002059602313486452
;; Maximal  point  deviation using CORNER with PERPENDICULAR: 0.0000000000000019389211565826713
;; Maximal tangent deviation using CORNER with PERPENDICULAR: 0.00040396101106856175
;; Maximal  point  deviation using CORNER with BARYCENTRIC: 0.00000009080625653904613
;; Maximal tangent deviation using CORNER with BARYCENTRIC: 0.00037085771283083065
;; Maximal  point  deviation using CORNER with RADIAL: 0.0000000000000045124795972965204
;; Maximal tangent deviation using CORNER with RADIAL: 0.0005765992792941358
;; Maximal  point  deviation using CORNER with CHORD-BASED: 0.000000000000003575193960459936
;; Maximal tangent deviation using CORNER with CHORD-BASED: 0.999986997783648
;; Maximal  point  deviation using CORNER with LINE-SWEEP: 0.00000000000003879258363621679
;; Maximal tangent deviation using CORNER with LINE-SWEEP: 0.0005762084466185979
;; Maximal  point  deviation using HYBRID with PERPENDICULAR: 0.749579458685565
;; Maximal tangent deviation using HYBRID with PERPENDICULAR: 0.9999779540452253
;; Maximal  point  deviation using HYBRID with BARYCENTRIC: 0.8472146626182875
;; Maximal tangent deviation using HYBRID with BARYCENTRIC: 0.5405037327175598
;; Maximal  point  deviation using HYBRID with RADIAL: 0.746717342966125
;; Maximal tangent deviation using HYBRID with RADIAL: 0.9997339367321527
;; Maximal  point  deviation using HYBRID with CHORD-BASED: 0.5724445094370448
;; Maximal tangent deviation using HYBRID with CHORD-BASED: 0.9805570318055594
;; Maximal  point  deviation using HYBRID with LINE-SWEEP: 0.000000000000003611319317287327
;; Maximal tangent deviation using HYBRID with LINE-SWEEP: 0.6858708327932557
;; Maximal  point  deviation using SKETCHES with PERPENDICULAR: 0.0000000000000026852712547870858
;; Maximal tangent deviation using SKETCHES with PERPENDICULAR: 0.0005292372078123327

;;; So HYBRID does not work, even for centralized line sweeps.

;;; As for CHORD-BASED, the derivatives seem to be OK for step=0.1
;;; (the normal vector deviates 4-5 degrees, comparable to other methods),
;;; but since it has very similar d parameter values near the sides,
;;; a small step like 1.0d-3 would lead to too similar (s,d) parameters,
;;; resulting in very large numerical computation errors.
#+nil
(let ((*resolution* 100))
  (format t "ribbon: ~f~%corner: ~f~%"
	  (check-patch-tangents *angles* 'ribbon :coords *coords* :distance-type 'chord-based
				:step 1.0d-2)
	  (check-patch-tangents *angles* 'corner :coords *coords* :distance-type 'chord-based
				:step 1.0d-2)))

(defun derivative-error-to-angle (err)
  (/ (* (asin err) 180) pi))

;;; Biquadratic tests:
#+nil
(let ((*resolution* 100)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* t)) 
  (iter (for type in '(ribbon corner hybrid))
	(iter (for distance in '(biquadratic))
	      (format t "Maximal  point  deviation using ~a with ~a: ~f~%"
		      type distance
		      (check-patch *points* type :coords *coords* :distance-type distance))
	      (format t "Maximal tangent deviation using ~a with ~a: ~f~%"
		      type distance
		      (check-patch-tangents *points* type :coords *coords* :distance-type distance
					    :step 1.0d-4)))))

;;; Maximal  point  deviation using RIBBON with BIQUADRATIC: 0.000000003326540158483576 < 10^-8
;;; Maximal tangent deviation using RIBBON with BIQUADRATIC: 0.0020216320104845657 ~ 0.1 degrees
;;; Maximal  point  deviation using CORNER with BIQUADRATIC: 0.000000000000002817430187 < 10^-14
;;; Maximal tangent deviation using CORNER with BIQUADRATIC: 0.0005757214242044044 ~ 0.03 degrees
;;; Maximal  point  deviation using HYBRID with BIQUADRATIC: 0.000000027472604973633572 < 10^-7
;;; Maximal tangent deviation using HYBRID with BIQUADRATIC: 0.0005757216288853803 ~ 0.03 degrees
