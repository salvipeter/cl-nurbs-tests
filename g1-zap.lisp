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

(defun ensure-g1-continuity (surface lsurface rsurface dsurface usurface)
  (let* ((cnet (control-net surface))
	 (nu (array-dimension cnet 0))
	 (nv (array-dimension cnet 1))
	 (lnet (list (control-net lsurface) t t))
	 (rnet (list (control-net rsurface) t nil))
	 (dnet (list (control-net dsurface) nil t))
	 (unet (list (control-net usurface) nil nil)))
    (flet ((get-index (array i u-dir from-end)
	     (if from-end (- (array-dimension array (if u-dir 0 1)) i 1) i))
	   (g1 (p0 p1 p)
	     (let ((n (vnormalize (v- p0 p1)))
		   (v (v- p p1)))
	       (v+ p1 (v* n (scalar-product v n))))))
      (dolist (net (list lnet rnet dnet unet))
	(iter (with u-dir = (second net))
	      (with end = (third net))
	      (for k from 1 below (if u-dir (1- nv) (1- nu)))
	      (for i = (if u-dir (get-index cnet 1 t (not end)) k))
	      (for j = (if u-dir k (get-index cnet 1 nil (not end))))
	      (for i0 = (if u-dir (get-index (first net) 0 t end) k))
	      (for j0 = (if u-dir k (get-index (first net) 0 nil end)))
	      (for i1 = (if u-dir (get-index (first net) 1 t end) k))
	      (for j1 = (if u-dir k (get-index (first net) 1 nil end)))
	      (setf (aref cnet i j)
		    (g1 (aref (first net) i0 j0)
			(aref (first net) i1 j1)
			(aref cnet i j))))))
    surface))
