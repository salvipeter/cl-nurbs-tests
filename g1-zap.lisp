;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun zap-to-curve (surface curve)
  "Zaps the lower u isocurve of SURFACE to CURVE."
  (let ((curve (bsc-reparameterize curve
				  (second (bss-lower-parameter surface))
				  (second (bss-upper-parameter surface))))
	(surface (copy-bspline-surface surface)))
    (iter (for i upfrom 0)
	  (while (and (< i (length (knot-vector curve)))
		      (< i (length (second (knot-vectors surface))))))
	  (for curve-knot = (elt (knot-vector curve) i))
	  (for surface-knot = (elt (second (knot-vectors surface)) i))
	  (cond ((> curve-knot surface-knot)
		 (setf curve (bsc-insert-knot curve surface-knot)))
		((< curve-knot surface-knot)
		 (setf surface (bss-insert-knot surface curve-knot
						:u-direction nil)))))
    (dotimes (i (length (control-points curve)))
      (setf (aref (control-net surface) 0 i) (elt (control-points curve) i)))
    surface))

(defun zap-to-surfaces (surface lsurface rsurface dsurface usurface)
  (let ((result surface))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 lsurface (first (bss-upper-parameter lsurface))
			 :u-curve nil)))
    (setf result (bss-flip-uv result))
;;; ----
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 dsurface (second (bss-upper-parameter dsurface))
			 :u-curve t)))
    (setf result (bss-flip-uv result))
;;; ----
    (setf result (bss-reverse-parameterization result :u t :v nil))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 rsurface (first (bss-lower-parameter rsurface))
			 :u-curve nil)))
    (setf result (bss-reverse-parameterization result :u t :v nil))
;;; ----
    (setf result (bss-flip-uv result))
    (setf result (bss-reverse-parameterization result :u t :v nil))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 usurface (second (bss-lower-parameter usurface))
			 :u-curve t)))
    (setf result (bss-reverse-parameterization result :u t :v nil))
    (setf result (bss-flip-uv result))
    result))

(defun bspline-basis (knots i p u &key (derivative 0))
  (cond ((> derivative 0)
	 (* p (- (safe-/ (bspline-basis knots i (1- p) u
					:derivative (1- derivative))
			 (- (elt knots (+ i p)) (elt knots i)))
		 (safe-/ (bspline-basis knots (1+ i) (1- p) u
					:derivative (1- derivative))
			 (- (elt knots (+ i p 1)) (elt knots (1+ i)))))))
	((and (= u (elt knots (1- (length knots))))
	      (= i (position u knots :from-end t :test-not #'=)))
	 1)
	((= p 0) (if (and (<= (elt knots i) u) (< u (elt knots (1+ i)))) 1 0))
	(t (+ (* (safe-/ (- u (elt knots i))
			 (- (elt knots (+ i p)) (elt knots i)))
		 (bspline-basis knots i (1- p) u))
	      (* (safe-/ (- (elt knots (+ i p 1)) u)
			 (- (elt knots (+ i p 1)) (elt knots (1+ i))))
		 (bspline-basis knots (1+ i) (1- p) u))))))

(defun zap-surfaces (surface master &key &allow-other-keys)
  "TODO: Works only at the lower U side of SURFACE."
  (zap-to-curve surface
		(bss-get-surface-curve
		 master (first (bss-upper-parameter master))
		 :u-curve nil)))

(defun insert-some-knots (surface n &key u-dir)
  "Inserts N knots into SURFACE at arbitrary places."
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst))))
    (iter (repeat (1+ n))
	  (for result first surface then
	       (bss-insert-knot result place :u-direction u-dir))
	  (for knots = (ufirst (knot-vectors result)))
	  (for place =
	       (iter (for i from 0 below (1- (length knots)))
		     (for k1 = (elt knots i))
		     (for k2 = (elt knots (1+ i)))
		     (finding (interpolate k1 0.5d0 k2) maximizing (- k2 k1))))
	  (finally (return result)))))

(defun ensure-g0-one-side (surface master resolution &key u-dir endp)
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
	 (usecond (lst) (if u-dir (second lst) (first lst))))
    (let ((missing-knots (- (length (usecond (knot-vectors master)))
			    (length (usecond (knot-vectors surface))))))
      (when (> missing-knots 0)
	(setf surface (insert-some-knots surface missing-knots
					 :u-dir (not u-dir)))))
    (let* ((p (degrees surface))
	   (knots (knot-vectors surface))
	   (net (control-net surface))
	   (n (1- (array-dimension net (if u-dir 0 1))))
	   (m (1- (array-dimension net (if u-dir 1 0))))
	   (u-low (ufirst (bss-lower-parameter surface)))
	   (u-high (ufirst (bss-upper-parameter surface)))
	   (v-low (usecond (bss-lower-parameter surface)))
	   (v-high (usecond (bss-upper-parameter surface)))
	   (u (if endp u-high u-low))
	   (length-v (- v-high v-low))
	   (index (if endp n 0))
	   (base-u (bspline-basis (ufirst knots) index (ufirst p) u))
	   (x (make-array (list (* resolution 3) (* (1+ m) 3))
			  :initial-element 0.0d0))
	   (y (make-array (list (* resolution 3) 1))))
      (iter (for k from 0 below resolution)
	    (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
	    (for uv = (if u-dir (list u vk) (list vk u)))
	    (iter (for j from 0 to m)
		  (for base-v = (bspline-basis (usecond knots) j
					       (usecond p) vk))
		  (for c = (* base-u base-v))
		  (iter (for r from 0 below 3)
			(setf (aref x (+ (* 3 k) r) (+ (* 3 j) r)) c)))
	    (for v1 = (bss-evaluate master uv))
	    (iter (for r from 0 below 3)
		  (setf (aref y (+ (* 3 k) r) 0) (elt v1 r))))
      (let ((solution (lu-solver:least-squares x y)))
        (iter (for j from 0 to m)
              (for i1 = (if u-dir index j))
              (for i2 = (if u-dir j index))
              (setf (aref net i1 i2)
		    (list (elt solution (+ (* j 3) 0))
			  (elt solution (+ (* j 3) 1))
			  (elt solution (+ (* j 3) 2))))))))
  surface)

(defun ensure-g0-continuity (surface lsurface rsurface dsurface usurface res)
  (let ((result (copy-bspline-surface surface)))
    (ensure-g0-one-side result lsurface (second res) :u-dir t :endp nil)
    (ensure-g0-one-side result rsurface (second res) :u-dir t :endp t)
    (ensure-g0-one-side result dsurface (first res) :u-dir nil :endp nil)
    (ensure-g0-one-side result usurface (first res) :u-dir nil :endp t)
    result))

(defun ensure-g1-one-side (surface master resolution &key u-dir endp
			   move-twists)
  "Ensure destructively G1-connectivity with MASTER at one side of SURFACE.

RESOLUTION gives the size of the LSQ-matrix used in the minimization.
The twist control points will not be moved, unless MOVE-TWISTS is T."
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
	 (usecond (lst) (if u-dir (second lst) (first lst))))
    (let* ((p (degrees surface))
	   (knots (knot-vectors surface))
	   (net (control-net surface))
	   (n (1- (array-dimension net (if u-dir 0 1))))
	   (m (1- (array-dimension net (if u-dir 1 0))))
	   (u-low (ufirst (bss-lower-parameter surface)))
	   (u-high (ufirst (bss-upper-parameter surface)))
	   (v-low (usecond (bss-lower-parameter surface)))
	   (v-high (usecond (bss-upper-parameter surface)))
	   (u (if endp u-high u-low))
	   (length-v (- v-high v-low))
	   (index (if endp (1- n) 1))
	   (dbase (bspline-basis (ufirst knots) index (ufirst p) u
				 :derivative 1))
	   (movable (if move-twists (- m 1) (- m 3)))
	   (first-movable (if move-twists 1 2))
	   (x (make-array (list (* (- resolution 2) 3) (* movable 3))))
	   (y (make-array (list (* (- resolution 2) 3) 1))))
      (iter (for k from 1 below (1- resolution))
	    (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
	    (for uv = (if u-dir (list u vk) (list vk u)))
	    (for n = (bss-surface-normal master uv))
	    (for c = (iter (for j from 0 below movable)
			   (for base-v = (bspline-basis (usecond knots)
							(+ j first-movable)
							(usecond p) vk))
			   (collect (* dbase base-v))))
	    (iter (for j from 0 below movable)
		  (iter (for r from 0 below 3)
			(iter (for s from 0 below 3)
			      (setf (aref x (+ (* 3 (1- k)) r) (+ (* 3 j) s))
				    (if (= r s) (elt c j) 0)))))
	    (for v1 = (bss-evaluate surface uv
				    :derivative (if u-dir '(1 0) '(0 1))))
	    (for right-side = (v* n (scalar-product v1 n) -1.0))
	    (iter (for r from 0 below 3)
		  (setf (aref y (+ (* 3 (1- k)) r) 0) (elt right-side r))))
      (let ((solution (lu-solver:least-squares x y)))
        (iter (for j from 0 below movable)
              (for i1 = (if u-dir index (+ j first-movable)))
              (for i2 = (if u-dir (+ j first-movable) index))
              (setf (aref net i1 i2)
                    (v+ (aref net i1 i2)
                        (list (elt solution (+ (* j 3) 0))
                              (elt solution (+ (* j 3) 1))
                              (elt solution (+ (* j 3) 2)))))))))
  surface)

(defun ensure-g1-continuity (surface lsurface rsurface dsurface usurface res)
  (let ((result (copy-bspline-surface surface)))
    (ensure-g1-one-side result lsurface (second res) :u-dir t :endp nil)
    (ensure-g1-one-side result rsurface (second res) :u-dir t :endp t)
    (ensure-g1-one-side result dsurface (first res) :u-dir nil :endp nil)
    (ensure-g1-one-side result usurface (first res) :u-dir nil :endp t)
    result))
