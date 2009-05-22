;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun halve-end-intervals (surface)
  (flet ((one-side (surface u-dir endp)
	   (flet ((ufirst (lst) (if u-dir (first lst) (second lst))))
	     (let* ((degree (ufirst (degrees surface)))
		    (knots (ufirst (knot-vectors surface)))
		    (u1 (if endp
			    (elt knots (- (length knots) degree 2))
			    (elt knots degree)))
		    (u2 (if endp
			    (elt knots (- (length knots) degree 1))
			    (elt knots (1+ degree))))
		    (u (interpolate u1 0.5 u2)))
	       (bss-insert-knot surface u :u-direction u-dir)))))
    (one-side (one-side (one-side (one-side surface t t) t nil) nil t) nil nil)))

(defun halve-all-intervals (surface)
  (flet ((one-side (surface u-dir)
	   (flet ((ufirst (lst) (if u-dir (first lst) (second lst))))
	     (iter (with result = surface)
		   (with knots = (ufirst (knot-vectors surface)))
		   (for k from 1 below (length knots))
		   (unless (= (elt knots (1- k)) (elt knots k))
		     (let ((u (/ (+ (elt knots (1- k)) (elt knots k)) 2)))
		       (setf result
			     (bss-insert-knot result u :u-direction u-dir))))
		   (finally (return result))))))
    (one-side (one-side surface t) nil)))

(defun halve-inner-intervals (surface)
  (flet ((one-side (surface u-dir)
	   (flet ((ufirst (lst) (if u-dir (first lst) (second lst))))
	     (iter (with result = surface)
		   (with degree = (ufirst (degrees surface)))
		   (with knots = (ufirst (knot-vectors surface)))
		   (for k from (+ degree 2)
			below (- (length knots) degree 1))
		   (unless (= (elt knots (1- k)) (elt knots k))
		     (let ((u (/ (+ (elt knots (1- k)) (elt knots k)) 2)))
		       (setf result
			     (bss-insert-knot result u :u-direction u-dir))))
		   (finally (return result))))))
    (one-side (one-side surface t) nil)))

(defun in-system (u v p)
  (let* ((2nd (if (= (second v) 0) #'third #'second))
	 (ux (first u))
	 (uy (funcall 2nd u))
	 (vx (first v))
	 (vy (funcall 2nd v))
	 (x (first p))
	 (y (funcall 2nd p))
	 (vx/vy (/ vx vy))
	 (a (/ (- x (* vx/vy y)) (- ux (* vx/vy uy)))))
    (list a (/ (- y (* a uy)) vy))))

(defun normal-curvature (surface uv direction)
  (let* ((der-00 (bss-evaluate surface uv :derivative '(1 0)))
	 (der-01 (bss-evaluate surface uv :derivative '(0 1)))
	 (der-10 (bss-evaluate surface uv :derivative '(2 0)))
	 (der-11 (bss-evaluate surface uv :derivative '(1 1)))
	 (der-12 (bss-evaluate surface uv :derivative '(0 2)))
	 (normal (vnormalize (cross-product der-00 der-01)))
	 (E (scalar-product der-00 der-00))
	 (F (scalar-product der-00 der-01))
	 (G (scalar-product der-01 der-01))
	 (L (scalar-product normal der-10))
	 (M (scalar-product normal der-11))
	 (N (scalar-product normal der-12)))
    (destructuring-bind (du dv) (in-system der-00 der-01 direction)
      (if (< du dv)
	  (let ((x (/ du dv)))
	    (/ (+ (* L x x) (* 2 M x) N)
	       (+ (* E x x) (* 2 F x) G)))
	  (let ((x (/ dv du)))
	    (/ (+ L (* 2 M x) (* N x x))
	       (+ E (* 2 F x) (* G x x))))))))

(defun control-normal (surface i j)
  (let ((net (control-net surface)))
    (vnormalize
     (cross-product (v- (aref net (1+ i) j) (aref net (1- i) j))
		    (v- (aref net i (1+ j)) (aref net i (1- j)))))))

;;; Uses the control-net-based outward direction
(defun ensure-g2-one-side-fixed-dir (surface master resolution &key u-dir endp
				     move-twists)
  "Ensure destructively G2-connectivity with MASTER at one side of SURFACE.

RESOLUTION gives the size of the LSQ-matrix used in the minimization.
The twist control points will not be moved, unless MOVE-TWISTS is T."
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
	 (usecond (lst) (if u-dir (second lst) (first lst))))
    (let* ((udegree (ufirst (degrees surface)))
	   (vdegree (usecond (degrees surface)))
	   (uknots (ufirst (knot-vectors surface)))
	   (vknots (usecond (knot-vectors surface)))
	   (net (control-net surface))
	   (n (1- (array-dimension net (if u-dir 0 1))))
	   (m (1- (array-dimension net (if u-dir 1 0))))
	   (u-low (ufirst (bss-lower-parameter surface)))
	   (u-high (ufirst (bss-upper-parameter surface)))
	   (v-low (usecond (bss-lower-parameter surface)))
	   (v-high (usecond (bss-upper-parameter surface)))
	   (u (if endp u-high u-low))
	   (length-v (- v-high v-low))
	   (index (if endp (- n 2) 2))
	   (movable (if move-twists (- m 3) (- m 5)))
	   (first-movable (if move-twists 2 3))
	   (x (make-array (list (- resolution 2) movable)))
	   (y (make-array (list (- resolution 2) 1))))
      (iter (for k from 1 below (1- resolution))
	    (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
	    (for uv = (if u-dir (list u vk) (list vk u)))
	    (for d1 = (bss-evaluate surface uv
				    :derivative (if u-dir '(1 0) '(0 1))))
	    (for d1-square = (scalar-product d1 d1))
	    (for k-master = (normal-curvature master uv d1))
	    (for k-surface = (normal-curvature surface uv d1))
	    (for normal = (bss-surface-normal surface uv))
	    (iter (for j from 0 below movable)
		  (setf (aref x (1- k) j)
			(* (bspline-basis vknots (+ j first-movable) vdegree vk)
			   (scalar-product
			      (if u-dir
				  (control-normal surface index
						  (+ j first-movable))
				  (control-normal surface (+ j first-movable)
						  index))
			      normal))))
	    (setf (aref y (1- k) 0)
		  (/ (* (- k-master k-surface) d1-square
			(- (elt uknots (+ udegree 1)) (elt uknots 2))
			(- (elt uknots (+ udegree 2)) (elt uknots 2)))
		     (* udegree (1- udegree)))))
      (let ((solution (lu-solver:least-squares x y)))
        (iter (for j from 0 below movable)
              (for i1 = (if u-dir index (+ j first-movable)))
              (for i2 = (if u-dir (+ j first-movable) index))
	      (for normal = (control-normal surface i1 i2))
              (setf (aref net i1 i2)
		    (v+ (aref net i1 i2)
			(v* normal (elt solution j))))))))
  surface)

;;; Minimizes the deviation from the original.
(defun ensure-g2-one-side-min-dev (surface master resolution &key u-dir endp
				   move-twists)
  "Ensure destructively G2-connectivity with MASTER at one side of SURFACE.

RESOLUTION gives the size of the LSQ-matrix used in the minimization.
The twist control points will not be moved, unless MOVE-TWISTS is T."
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
	 (usecond (lst) (if u-dir (second lst) (first lst))))
    (let* ((udegree (ufirst (degrees surface)))
	   (vdegree (usecond (degrees surface)))
	   (uknots (ufirst (knot-vectors surface)))
	   (vknots (usecond (knot-vectors surface)))
	   (net (control-net surface))
	   (n (1- (array-dimension net (if u-dir 0 1))))
	   (m (1- (array-dimension net (if u-dir 1 0))))
	   (u-low (ufirst (bss-lower-parameter surface)))
	   (u-high (ufirst (bss-upper-parameter surface)))
	   (v-low (usecond (bss-lower-parameter surface)))
	   (v-high (usecond (bss-upper-parameter surface)))
	   (u (if endp u-high u-low))
	   (length-v (- v-high v-low))
	   (index (if endp (- n 2) 2))
	   (movable (if move-twists (- m 3) (- m 5)))
	   (first-movable (if move-twists 2 3))
	   (x (make-array (list (* 3 (- resolution 2)) (* 3 movable))
			  :initial-element 0.0d0))
	   (y (make-array (list (* 3 (- resolution 2)) 1))))
      (iter (for k from 1 below (1- resolution))
	    (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
	    (for uv = (if u-dir (list u vk) (list vk u)))
	    (for d1 = (bss-evaluate surface uv
				    :derivative (if u-dir '(1 0) '(0 1))))
	    (for d1-square = (scalar-product d1 d1))
	    (for k-master = (normal-curvature master uv d1))
	    (for k-surface = (normal-curvature surface uv d1))
	    (for normal = (bss-surface-normal surface uv))
	    (iter (for j from 0 below movable)
		  (iter (for s from 0 below 3)
			(setf (aref x (+ (* 3 (1- k)) s) (+ (* 3 j) s))
			      (bspline-basis vknots (+ j first-movable)
					     vdegree vk))))
	    (for right-side =
		 (v* normal
		     (/ (* (- k-master k-surface) d1-square
			   (- (elt uknots (+ udegree 1)) (elt uknots 2))
			   (- (elt uknots (+ udegree 2)) (elt uknots 2)))
			(* udegree (1- udegree)))))
	    (iter (for s from 0 below 3)
		  (setf (aref y (+ (* 3 (1- k)) s) 0)
			(elt right-side s))))
      (let ((solution (lu-solver:least-squares x y)))
        (iter (for j from 0 below movable)
              (for i1 = (if u-dir index (+ j first-movable)))
              (for i2 = (if u-dir (+ j first-movable) index))
	      (for normal = (control-normal surface i1 i2))
              (setf (aref net i1 i2)
		    (v+ (aref net i1 i2)
			(list (elt solution (+ (* 3 j) 0))
			      (elt solution (+ (* 3 j) 1))
			      (elt solution (+ (* 3 j) 2)))))))))
  surface)

(defun ensure-g2-one-side (surface master resolution &key u-dir endp move-twists
			   (algorithm 'minimal-deviation))
  (ecase algorithm
    (minimal-deviation
     (ensure-g2-one-side-min-dev surface master resolution
				 :u-dir u-dir :endp endp
				 :move-twists move-twists))
    (fixed-direction
     (ensure-g2-one-side-fixed-dir surface master resolution
				   :u-dir u-dir :endp endp
				   :move-twists move-twists))))

(defun ensure-g2-continuity (surface lsurface rsurface dsurface usurface res
			     &key (algorithm 'minimal-deviation))
  (let ((result (copy-bspline-surface surface)))
    (ensure-g2-one-side result lsurface (second res)
			:u-dir t :endp nil :algorithm algorithm)
    (ensure-g2-one-side result rsurface (second res)
			:u-dir t :endp t :algorithm algorithm)
    (ensure-g2-one-side result dsurface (first res)
			:u-dir nil :endp nil :algorithm algorithm)
    (ensure-g2-one-side result usurface (first res)
			:u-dir nil :endp t :algorithm algorithm)
    result))

;;; Pl:
#+nil
(apply #'ensure-g2-continuity (halve-end-intervals (first *xnode*))
       (append (rest *xnode*) '((100 100))))

(in-package :cl-nurbs)

;;; otlet:
;;; minden kontrollponthoz generaljunk parameterertekeket, amik lefedik azt
;;; az intervallumot, ahol a kontrollpontnak hatasa van, es a fairness
;;; fuggveny csak ezeknek a valtozasat nezze.

#+nil
(defun bss-fair (surface &key (measure #'fairness) (resolution 100)
		 (target-iteration 100) (simplex-iteration 5)
		 (fairing-iteration 5) (lock-endpoints t))
  (let* ((parameters (arc-length-sampling curve from to resolution))
	 (curvatures (target-curvature curve parameters target-iteration))
	 (new-surface (copy-bspline-surface surface))
	 (net (control-net new-surface))
	 (size (array-dimensions net)))
    (dotimes (k fairing-iteration)
      (let ((low (if lock-endpoints 1 0))
	    (high-u (if lock-endpoints (1- (first size)) (first size)))
	    (high-v (if lock-endpoints (1- (second size)) (second size))))
	(iter (for i from low below high-u)
	      (iter (for j from low below high-v)
		    (flet ((fairness-fn (point)
			     (let ((old (elt points i j)))
			       (setf (elt points i j) point)
			       (prog1 (funcall measure new-surface
					       parameters curvatures)
				 (setf (elt points i j) old)))))
		      (setf (elt points i j)
			    (downhill-simplex:minimize
			     #'fairness-fn (elt points i j)
			     simplex-iteration)))))))
    new-surface))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Old implementation with two curves going into the corners ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; (defun surface-curve-derivatives (surface v &key u-dir endp type)
;;   "Returns the derivatives at the end of the curve.

;; TYPE is UP or DOWN."
;;   (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
;; 	 (usecond (lst) (if u-dir (second lst) (first lst))))
;;     (let* ((umin (ufirst (bss-lower-parameter surface)))
;; 	   (umax (ufirst (bss-upper-parameter surface)))
;; 	   (vmin (usecond (bss-lower-parameter surface)))
;; 	   (vmax (usecond (bss-upper-parameter surface)))
;; 	   (n (ufirst (array-dimensions (control-net surface))))
;; 	   (m (usecond (array-dimensions (control-net surface))))
;; 	   (uknots (ufirst (knot-vectors surface)))
;; 	   (vknots (usecond (knot-vectors surface)))
;; 	   (udegree (ufirst (degrees surface)))
;; 	   (vdegree (usecond (degrees surface)))
;; 	   (index (if endp (- n 3) 2))
;; 	   (d '(0 0 0))
;; 	   (dd '(0 0 0))
;; 	   (alpha '())
;; 	   (u (if endp umax umin))
;; 	   (du (- umax umin))
;; 	   (dv (ecase type (up (- vmax v)) (down (- vmin v)))))
;;       (iter (for i from 0 below n)
;; 	    (iter (for j from 0 below m)
;; 		  (for alpha1 =
;; 		       (+ (* (bspline-basis uknots i udegree u :derivative 1)
;; 			     (bspline-basis vknots j vdegree v)
;; 			     du)
;; 			  (* (bspline-basis uknots i udegree u)
;; 			     (bspline-basis vknots j vdegree v :derivative 1)
;; 			     dv)))
;; 		  (for alpha2 =
;; 		       (+ (* (bspline-basis uknots i udegree u :derivative 2)
;; 			     (bspline-basis vknots j vdegree v)
;; 			     du du)
;; 			  (* (bspline-basis uknots i udegree u :derivative 1)
;; 			     (bspline-basis vknots j vdegree v :derivative 1)
;; 			     2 du dv)
;; 			  (* (bspline-basis uknots i udegree u)
;; 			     (bspline-basis vknots j vdegree v :derivative 2)
;; 			     dv dv)))
;; 		  (for p = (if u-dir
;; 			       (aref (control-net surface) i j)
;; 			       (aref (control-net surface) j i)))
;; 		  (setf d (v+ d (v* p alpha1)))
;; 		  (setf dd (v+ dd (v* p alpha2)))
;; 		  (when (and (= i index) (<= 3 j (- m 4)))
;; 		    (push alpha2 alpha))))
;;       (list d dd (nreverse alpha)))))

;; (defun ensure-g2-one-side (surface master resolution &key u-dir endp)
;;   "Ensure destructively G2-connectivity with MASTER at one side of SURFACE.

;; RESOLUTION gives the size of the LSQ-matrix used in the minimization.
;; The twist control points will not be moved."
;;   (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
;; 	 (usecond (lst) (if u-dir (second lst) (first lst))))
;;     (let* ((net (control-net surface))
;; 	   (n (1- (array-dimension net (if u-dir 0 1))))
;; 	   (m (1- (array-dimension net (if u-dir 1 0))))
;; 	   (u-low (ufirst (bss-lower-parameter surface)))
;; 	   (u-high (ufirst (bss-upper-parameter surface)))
;; 	   (v-low (usecond (bss-lower-parameter surface)))
;; 	   (v-high (usecond (bss-upper-parameter surface)))
;; 	   (u (if endp u-high u-low))
;; 	   (length-v (- v-high v-low))
;; 	   (index (if endp (- n 2) 2))
;; 	   (x (make-array (list (* (- resolution 2) 2) (- m 5))))
;; 	   (y (make-array (list (* (- resolution 2) 2) 1))))
;;       (dolist (dir '(up down))
;; 	(iter (with start = (if (eq dir 'up) 0 (- resolution 2)))
;; 	      (for k from 1 below (1- resolution))
;; 	      (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
;; 	      (for uv = (if u-dir (list u vk) (list vk u)))
;; 	      (for (d1 d2 alpha) =
;; 		   (surface-curve-derivatives surface vk :u-dir u-dir
;; 					      :endp endp :type dir))
;; 	      (for target-curvature = (normal-curvature master uv d1))
;; 	      (for normal = (bss-surface-normal surface uv))
;; 	      (iter (for j from 0 to (- m 6))
;; 		    (setf (aref x (+ start k -1) j)
;; 			  (* (elt alpha j)
;; 			     (scalar-product
;; 			      (if u-dir
;; 				  (control-normal surface index (+ j 3))
;; 				  (control-normal surface (+ j 3) index))
;; 			      normal))))
;; 	      (setf (aref y (+ start k -1) 0)
;; 		    (- (* target-curvature (scalar-product d1 d1))
;; 		       (scalar-product d2 normal)))))
;;       (let ((solution (lu-solver:least-squares x y)))
;;         (iter (for j from 3 to (- m 3))
;;               (for i1 = (if u-dir index j))
;;               (for i2 = (if u-dir j index))
;; 	      (for normal = (if u-dir
;; 				(control-normal surface index j)
;; 				(control-normal surface j index)))
;;               (setf (aref net i1 i2)
;; 		    (v+ (aref net i1 i2)
;; 			(v* normal (elt solution (- j 3)))))))))
;;   surface)
