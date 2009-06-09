;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun project-to-plane (point plane)
  "PLANE is represented by a list \(origin normal-vector)."
  (let ((len (scalar-product (v- (first plane) point) (second plane))))
    (v+ point (v* (second plane) len))))

(defun bsc-ideal-control-point (curve i &key (iteration 100) (alpha 1/2)
				just-error-p)
  (flet ((knot (k) (elt (knot-vector curve) (+ k (degree curve))))
	 (cpoint (k) (elt (control-points curve) k)))
    (let* ((u (knot i))
	   (d (vnormalize (bsc-evaluate curve u :derivative 1)))
	   (n (bsc-out-direction curve u))
	   (p2 (cpoint (- i 2))) (p1 (cpoint (- i 1)))
	   (r1 (cpoint (+ i 1))) (r2 (cpoint (+ i 2)))
	   (plane (list (affine-combine p1 1/2 r1) (cross-product d n)))
	   (p2v (project-to-plane p2 plane)) (p1v (project-to-plane p1 plane))
	   (qv (project-to-plane (cpoint i) plane))
	   (r2v (project-to-plane r2 plane)) (r1v (project-to-plane r1 plane))
	   (krr (bsc-krr-better-positions curve i)))
      (labels ((from-plane (v)
		 (destructuring-bind (x y) v
		   (v+ qv (v* n x) (v* d y))))
	       (angle (v1 v2)
		 (let ((n1 (vnormalize v1))
		       (n2 (vnormalize v2)))
		   (* (acos (min (max (- (scalar-product n1 n2)) -1.0d0) 1.0d0))
		      (let ((x (cross-product n1 n2)))
			(if (> (scalar-product x (second plane)) 0) 1 -1)))))
	       (error-fn (v)
		 (let* ((q (from-plane v))
			(v-2 (v- p1v p2v)) (v-1 (v- q p1v))
			(v+1 (v- r1v q)) (v+2 (v- r2v r1v))
			(len-1 (vlength v-1)) (len+1 (vlength v+1))
			(len (interpolate len-1 1/2 len+1))
			(ang-1 (angle v-2 v-1))
			(ang   (angle v-1 v+1))
			(ang+1 (angle v+1 v+2))
			(delta (- (interpolate ang-1
					       (/ len-1 (+ len-1 len+1))
					       ang+1)
				  ang)))
		   (+ (* (- 1.0d0 alpha) (expt (* (/ delta pi) len) 2))
		      (* alpha
			 (reduce #'min
				 (mapcar (lambda (p) (vlength2 (v- p q)))
					 krr)))))))
	(if just-error-p
	    (error-fn '(0 0))
	    (from-plane
	     (downhill-simplex:minimize #'error-fn '(0 0) iteration)))))))

(defun limit-change (old new limit)
  "Moves OLD towards NEW by length LIMIT."
  (v+ old (v* (vnormalize (v- new old)) limit)))

(defun limit-deviation (origin radius old new)
  "We want to move from OLD to NEW, but don't want to go out of the circle
given by ORIGIN and RADIUS, so the vector is cut off at the intersection."
  (if (< (point-distance origin new) radius)
      new
      (let ((a (v- old origin))
	    (v (vnormalize (v- new old))))
	(v+ old (v* v (- (sqrt (max 0 (+ (expt (scalar-product a v) 2)
					 (- (vlength2 a)) (* radius radius))))
			 (scalar-product a v)))))))

(defun bsc-fair-control-point (curve i &key (iteration 100) (alpha 1/2)
			       limit dev-origin dev-limit)
  (let* ((limit (or limit
		    (/ (bsc-estimate-arc-length curve)
		       (* (length (control-points curve)) 30))))
	 (dev-limit (or dev-limit (* limit 5)))
	 (dev-origin (or dev-origin (elt (control-points curve) i))))
    (setf (elt (control-points curve) i)
	  (limit-deviation dev-origin dev-limit
			   (elt (control-points curve) i)
			   (limit-change (elt (control-points curve) i)
					 (bsc-ideal-control-point
					  curve i :iteration iteration
					  :alpha alpha)
					 limit)))))

(defun bsc-2d-to-3d (curve)
  (let* ((new-curve (copy-bspline-curve curve))
	 (cps (control-points new-curve)))
    (dotimes (i (length cps))
      (setf (elt cps i) (list (first (elt cps i)) (second (elt cps i)) 0)))
    new-curve))

(defun bsc-3d-to-2d (curve)
  (let* ((new-curve (copy-bspline-curve curve))
	 (cps (control-points new-curve)))
    (dotimes (i (length cps))
      (setf (elt cps i) (subseq (elt cps i) 0 2)))
    new-curve))

#+nil
(defun test-to-file (curve n filename &key (remembering 3) alpha limit)
  (let ((c (bsc-2d-to-3d (copy-bspline-curve curve)))
	(k (- (length (control-points curve)) 2))
	(last (iter (repeat remembering) (collect 0))))
    (flet ((select-next ()
	     (iter (for index from 2 below k)
		   (unless (member index last)
		     (for error =
			  (bsc-ideal-control-point c index :alpha alpha
						   :just-error-p t))
		     (finding index maximizing error)))))
      (iter (for i from 0 to n)
	    (for fname = (format nil "~a-~2,'0d.ps" filename i))
	    (write-ps (bsc-3d-to-2d c) fname 100 :show-knots-p t)
	    (for next = (select-next))
	    (setf last (cons next (butlast last)))
	    (bsc-fair-control-point
	     c next :alpha alpha :limit limit
	     :dev-origin (append (elt (control-points curve) next) '(0)))))
    c))

;;; (test-to-file *curve* 50 "/tmp/test" :alpha 0.7)

;;; Using the notations of
;;; S. Hahmann : Knot-Removal Surface Fairing using Search Strategies
;;; Note that there is a misprint in the definition of alpha_r, it should be
;;;            t_{j+1}   - t~_r
;;; alpha_r = ----------------- .
;;;           t~_{r+k-1} - t~_r
(defun bsc-remove-and-reinsert-knot (curve knot &key (using :lsq))
  "USING can be :LSQ or the index of the \(only) control point to change.
In the first case, it returns all the control points,
otherwise only the new coordinates of the specified one."
  (let ((j (1- knot))
	(k (1+ (degree curve)))
	(dim (bsc-dimension curve))
	(knots (knot-vector curve))
	(points (control-points curve)))
    (assert (or (eq using :lsq) (<= (- j k -2) using j)) (using)
	    "USING should be :LSQ or between ~d and ~d" (- j k -2) j)
    (assert (<= k knot (- (length knots) (max k 2) 1)) (knot)
	    "KNOT should be between ~d and ~d"
	    k (- (length knots) (max k 2) 1))
    (labels ((u (r) (elt knots r))
	     (u- (r) (if (<= r j) (u r) (u (1+ r))))
	     (alpha (r) (safe-/ (- (u (1+ j)) (u- r))
				(- (u- (+ r k -1)) (u- r)))))
      (declare (inline u u-))
      (let* ((rows (if (eq using :lsq) (- k 1) (- k 2)))
	     (matrix (make-array (list (* rows dim) (* (- k 2) dim))
				 :initial-element 0.0d0))
	     (result (make-array (* rows dim))))
	(iter (with real-row = 0)
	      (for row upfrom 0)
	      (for next upfrom (- j k -2))
	      (while (< real-row rows))
	      (when (or (eq using :lsq) (/= using next))
		(unless (= row 0)
		  (dotimes (d dim)
		    (setf (aref matrix
				(+ (* real-row dim) d) (+ (* (1- row) dim) d))
			  (- 1.0d0 (alpha next)))))
		(unless (= row (- k 2))
		  (dotimes (d dim)
		    (setf (aref matrix
				(+ (* real-row dim) d) (+ (* row dim) d))
			  (alpha next))))
		(for p = (cond ((= row 0)
				(v- (elt points next)
				    (v* (elt points (1- next))
					(- 1.0d0 (alpha next)))))
			       ((= row (- k 2))
				(v- (elt points next)
				    (v* (elt points (1+ next))
					(alpha next))))
			       (t (elt points next))))
		(dotimes (d dim)
		  (setf (elt result (+ (* real-row dim) d))
			(elt p d)))
		(incf real-row)))
	(flet ((combined (vec)
		 (concatenate 'vector
			      (list (elt points (- j k -1)))
			      (iter (for i from 0 below (/ (length vec) dim))
				    (collect
					(iter (for d from 0 below dim)
					      (collect
						  (elt vec (+ (* i dim) d))))))
			      (list (elt points (1+ j)))))
	       (compute-point (vec index)
		 (let ((v-index (- index (- j k -2))))
		   (v+ (v* (elt vec v-index)
			   (- 1.0d0 (alpha index)))
		       (v* (elt vec (1+ v-index))
			   (alpha index))))))
	  (if (eq using :lsq)
	      (let ((d (combined (lu-solver:least-squares
				  matrix (matrix:from-vector result)))))
		(concatenate 'list
			     (iter (for i from 0 below (- j k -2))
				   (collect (elt points i)))
			     (iter (for i from (- j k -2) to j)
				   (collect (compute-point d i)))
			     (iter (for i from (1+ j) below (length points))
				   (collect (elt points i)))))
	      (let ((d (combined (lu-solver:solve matrix result))))
		(compute-point d using))))))))

(defun bsc-krr-better-positions (curve index)
  "Computes better positions for the control point designated by INDEX.
Returns a list of points, containing at most \(DEGREE CURVE) elements."
  (let ((k (1+ (degree curve)))
	(n (length (knot-vector curve))))
    (iter (for knot from k to (- n (max k 2) 1))
	  (when (< (- knot k) index knot)
	    (collect (bsc-remove-and-reinsert-knot curve knot :using index))))))

#+nil
(defun test (lst)
  (let ((curve (copy-bspline-curve *curve*)))
    (dolist (i lst)
      (let ((ideal (bsc-ideal-control-point (bsc-2d-to-3d curve) i)))
	(setf (elt (control-points curve) i) (subseq ideal 0 2))))
    (write-ps curve "/tmp/proba.ps" 100)))

(defun bss-ideal-control-point (surface i j &key (iteration 100) (alpha 1/2)
				just-error-p)
  (let* ((cu (bss-construction-curve surface j :u-direction t))
	 (cv (bss-construction-curve surface i :u-direction nil))
	 (pu (bsc-ideal-control-point cu i :iteration iteration :alpha alpha
				      :just-error-p just-error-p))
	 (pv (bsc-ideal-control-point cv j :iteration iteration :alpha alpha
				      :just-error-p just-error-p)))
    (if just-error-p
	(interpolate pu 1/2 pv)
	(affine-combine pu 1/2 pv))))
