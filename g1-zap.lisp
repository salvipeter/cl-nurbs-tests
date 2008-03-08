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

(defun least-squares (x y)
  "Solve (X^T X)b = X^T y, where X is an (n,m) matrix and Y is an (n) vector."
  (lu-solve (matrix-multiplication (matrix-transpose x) x)
	    (matrix->vector (matrix-multiplication (matrix-transpose x) y))))

(defun ensure-g1-one-side (surface master resolution &key u-dir endp)
  "Ensure destructively G1-connectivity with MASTER at one side of SURFACE.

RESOLUTION gives the size of the LSQ-matrix used in the minimization.
The twist control points will not be moved."
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
	   (x (make-array (list (* (- resolution 2) 3) (* (- m 3) 3))))
	   (y (make-array (list (* (- resolution 2) 3) 1))))
      (iter (for k from 1 below (1- resolution))
	    (for vk = (+ v-low (/ (* k length-v) (1- resolution))))
	    (for uv = (if u-dir (list u vk) (list vk u)))
	    (for n = (bss-surface-normal master uv))
	    (for c = (iter (for j from 2 to (- m 2))
			   (for base-v = (bspline-basis (usecond knots) j
							(usecond p) vk))
			   (collect (* dbase base-v))))
	    (iter (for j from 0 to (- m 4))
		  (iter (for r from 0 below 3)
			(iter (for s from 0 below 3)
			      (setf (aref x (+ (* 3 (1- k)) r) (+ (* 3 j) s))
				    (if (= r s) (elt c j) 0)))))
	    (for v1 = (bss-evaluate surface uv
				    :derivative (if u-dir '(1 0) '(0 1))))
	    (for right-side = (v* n (scalar-product v1 n) -1.0))
	    (iter (for r from 0 below 3)
		  (setf (aref y (+ (* 3 (1- k)) r) 0) (elt right-side r))))
      (let ((solution (least-squares x y)))
        (iter (for j from 2 to (- m 2))
              (for i1 = (if u-dir index j))
              (for i2 = (if u-dir j index))
              (setf (aref net i1 i2)
                    (v+ (aref net i1 i2)
                        (list (elt solution (+ (* (- j 2) 3) 0))
                              (elt solution (+ (* (- j 2) 3) 1))
                              (elt solution (+ (* (- j 2) 3) 2))))))))))

(defun ensure-g1-continuity (surface lsurface rsurface dsurface usurface res)
  (let ((result (copy-bspline-surface surface)))
    (ensure-g1-one-side result lsurface (second res) :u-dir t :endp nil)
    (ensure-g1-one-side result rsurface (second res) :u-dir t :endp t)
    (ensure-g1-one-side result dsurface (first res) :u-dir nil :endp nil)
    (ensure-g1-one-side result usurface (first res) :u-dir nil :endp t)
    result))
