;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

(defun reparametrize-curve (curve min max)
  (let* ((old-min (bsc-lower-parameter curve))
	 (old-max (bsc-upper-parameter curve))
	 (knots (knot-vector curve))
	 (n (length knots))
	 (points (control-points curve))
	 (m (length points))
	 (new-knots (make-array n))
	 (new-points (make-array m)))
    (dotimes (i n)
      (setf (elt new-knots i)
	    (+ min (* (- max min) (/ (- (elt knots i) old-min)
				     (- old-max old-min))))))
    (dotimes (i m)
      (setf (elt new-points i) (copy-list (elt points i))))
    (make-bspline-curve (degree curve) new-knots new-points)))

(defun flip-uv (surface)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((n (array-dimension net 0))
	   (m (array-dimension net 1))
	   (new-net (make-array (list m n))))
      (dotimes (i n)
	(dotimes (j m)
	  (setf (aref new-net j i) (aref net i j))))
      (make-bspline-surface (reverse degrees) (reverse knots) new-net))))

(defun reverse-parametrization (surface &key u v)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((n (array-dimension net 0))
	   (m (array-dimension net 1))
	   (k (length (first knots)))
	   (l (length (second knots)))
	   (new-knots (list (if u (make-array k) (copy-seq (first knots)))
			    (if v (make-array l) (copy-seq (second knots)))))
	   (new-net (make-array (list n m))))
      (when u
	(let ((low (elt (first knots) 0))
	      (high (elt (first knots) (1- k))))
	  (dotimes (i k)
	    (setf (elt (first new-knots) i)
		  (+ low (- high (elt (first knots) (- k i 1))))))))
      (when v
	(let ((low (elt (second knots) 0))
	      (high (elt (second knots) (1- l))))
	  (dotimes (i l)
	    (setf (elt (second new-knots) i)
		  (+ low (- high (elt (second knots) (- l i 1))))))))
      (dotimes (i n)
	(dotimes (j m)
	  (setf (aref new-net i j)
		(aref net (if u (- n i 1) i) (if v (- m j 1) j)))))
      (make-bspline-surface (copy-list degrees) new-knots new-net))))

(defun zap-to-curve (surface curve)
  "Zaps the u=0 isocurve of SURFACE to CURVE."
  (let ((curve (reparametrize-curve curve
				    (second (bss-lower-parameter surface))
				    (second (bss-upper-parameter surface))))
	(surface (copy-bspline-surface surface)))
    (loop with i = 0
	  for curve-short = (>= i (length (knot-vector curve)))
	  for surface-short = (>= i (length (second (knot-vectors surface))))
	  while (and (not curve-short) (not surface-short)) do
	  (if (or curve-short
		  (and (not surface-short)
		       (> (elt (knot-vector curve) i)
			  (elt (second (knot-vectors surface)) i))))
	      (setf curve
		    (bsc-insert-knot curve
				     (elt (second (knot-vectors surface)) i)))
	      (when (or surface-short
			(< (elt (knot-vector curve) i)
			   (elt (second (knot-vectors surface)) i)))
		(setf surface
		      (bsc-insert-knot surface
				       (elt (knot-vector curve) i)
				       :u-direction nil))))
	  (setf i (1+ i)))
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
    (setf result (flip-uv result))
;;; ----
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 dsurface (second (bss-upper-parameter dsurface))
			 :u-curve t)))
    (setf result (flip-uv result))
;;; ----
    (setf result (reverse-parametrization result :u t :v nil))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 rsurface (first (bss-lower-parameter rsurface))
			 :u-curve nil)))
    (setf result (reverse-parametrization result :u t :v nil))
;;; ----
    (setf result (flip-uv result))
    (setf result (reverse-parametrization result :u t :v nil))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 usurface (second (bss-lower-parameter usurface))
			 :u-curve t)))
    (setf result (reverse-parametrization result :u t :v nil))
    (setf result (flip-uv result))
    result))

(defun ensure-g1-continuity (surface lsurface rsurface dsurface usurface)
  (let* ((nu (array-dimension (control-net surface) 0))
	 (nv (array-dimension (control-net surface) 1))
	 (cnet (control-net surface))
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
	      (for k from 1 below (if (second net) (1- nv) (1- nu)))
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
