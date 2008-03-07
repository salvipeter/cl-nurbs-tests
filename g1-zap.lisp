;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun zap-to-curve (surface curve)
  "Zaps the u=0 isocurve of SURFACE to CURVE."
  (let ((curve (bsc-reparameterize curve
				  (second (bss-lower-parameter surface))
				  (second (bss-upper-parameter surface))))
	(surface (copy-bspline-surface surface)))
    (iter (for i first 0 then (1+ i))
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
	((= p 0) (if (and (<= (elt knots i) u) (< u (elt knots (1+ i)))) 1 0))
	(t (+ (* (safe-/ (- u (elt knots i))
			 (- (elt knots (+ i p)) (elt knots i)))
		 (bspline-basis knots i (1- p) u))
	      (* (safe-/ (- (elt knots (+ i p 1)) u)
			 (- (elt knots (+ i p 1)) (elt knots (1+ i))))
		 (bspline-basis knots (1+ i) (1- p) u))))))

(defun least-squares (x y)
  "Solve (X^T X)b = X^T y, where X is an (n,m) matrix and Y is an (n) vector."
  (lu-solve (matrix-multiplication (matrix-transpose x) x)
	    (matrix->vector (matrix-multiplication (matrix-transpose x) y))))

(defun ensure-g1-one-side (surface master resolution &key u-dir endp)
  "Ensure destructively G1-connectivity with MASTER at one side of SURFACE.

RESOLUTION gives the size of the LSQ-matrix used in the minimization.
The twist control points will not be moved."
  (declare (ignore u-dir endp))		; TODO
  (let* ((p (degrees surface))
	 (knots (knot-vectors surface))
	 (net (control-net surface))
	 (m (1- (array-dimension net 1)))
	 (u (first (bss-lower-parameter surface)))
	 (base-v (second (bss-lower-parameter surface)))
	 (length-v (- (second (bss-upper-parameter surface)) base-v))
	 (dbase (bspline-basis (first knots) 1 (first p) u :derivative 1))
	 (x (make-array (list (* (- resolution 2) 3) (* (- m 3) 3))))
	 (y (make-array (list (* (- resolution 2) 3) 1))))
    (iter (for k from 1 below (1- resolution))
	  (for vk = (+ base-v (/ (* k length-v) (1- resolution))))
	  (for n = (bss-surface-normal surface (list u vk)))
	  (for v0 = (bss-evaluate master (list u vk) :derivative '(1 0)))
	  (for no = (cross-product n v0))
	  (for c = (iter (for j from 2 to (- m 2))
			 (for base-v = (bspline-basis (second knots) j
						      (second p) vk))
			 (collect (* dbase base-v))))
	  (iter (for j from 0 to (- m 4))
		(iter (for r from 0 below 3)
		      (iter (for s from 0 below 3)
			    (setf (aref x (+ (* 3 (1- k)) r) (+ (* 3 j) s))
				  (- (if (= r s) (elt c j) 0)
				     (* (elt no r) (elt c j) (elt no s)))))))
	  (for v1 = (bss-evaluate surface (list u vk) :derivative '(1 0)))
	  (for projected-v1 = (v- v1 (v* no (scalar-product v1 no))))
	  (for right-side = (v- v0 projected-v1))
	  (iter (for r from 0 below 3)
		(setf (aref y (+ (* (1- k) 3) r) 0) (elt right-side r))))
    (let ((solution (least-squares x y)))
      (iter (for j from 2 to (- m 2))
	    (setf (aref net 1 j)
		  (v+ (aref net 1 j)
		      (list (elt solution (+ (* (- j 2) 3) 0))
			    (elt solution (+ (* (- j 2) 3) 1))
			    (elt solution (+ (* (- j 2) 3) 2)))))))))

(defun ensure-g1-continuity (surface lsurface rsurface dsurface usurface res)
  (declare (ignore rsurface dsurface usurface)) ; TODO
  (let ((result (copy-bspline-surface surface)))
    (ensure-g1-one-side result lsurface (second res) :u-dir t :endp nil)
;;     (ensure-g1-one-side result rsurface (second res) :u-dir t :endp t)
;;     (ensure-g1-one-side result dsurface (first res) :u-dir nil :endp nil)
;;     (ensure-g1-one-side result usurface (first res) :u-dir nil :endp t)
    result))
