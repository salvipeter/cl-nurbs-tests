;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun bsc-reparameterize (curve min max)
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

(defun bsc-reverse-parameterization (curve)
  (let* ((points (control-points curve))
	 (knots (knot-vector curve))
	 (n (length points))
	 (k (length knots))
	 (new-points (make-array n))
	 (new-knots (make-array k))
	 (lower (bsc-lower-parameter curve))
	 (upper (bsc-upper-parameter curve)))
    (dotimes (i k)
      (setf (elt new-knots i) (+ lower (- upper (elt knots (- k i 1))))))
    (dotimes (i n)
      (setf (aref new-points i) (aref points (- n i 1))))
    (make-bspline-curve (degree curve) new-knots new-points)))

(defun bsc-shift-parameterization (curve shift)
  (let ((result (copy-bspline-curve curve)))
    (setf (knot-vector result)
	  (mapcar (lambda (x) (+ x shift))
		  (coerce (knot-vector result) 'list)))
    result))

(defun bss-flip-uv (surface)
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

(defun bss-reverse-parameterization (surface &key u v)
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

(defun bss-shift-parameterization (surface &key (u 0.0) (v 0.0))
  (let ((result (copy-bspline-surface surface)))
    (setf (knot-vectors result)
	  (mapcar (lambda (lst shift)
		    (mapcar (lambda (x) (+ x shift)) (coerce lst 'list)))
		  (knot-vectors result) (list u v)))
    result))

(defun bss-reparametrize (surface umin umax vmin vmax)
  (let* ((result (copy-bspline-surface surface))
	 (low (bss-lower-parameter surface))
	 (length (mapcar #'- (bss-upper-parameter surface) low)))
    (setf (knot-vectors result)
	  (mapcar (lambda (knots min new-min len scale)
		    (mapcar (lambda (x)
			      (+ new-min (/ (* (- x min) scale) len)))
			    (coerce knots 'list)))
		  (knot-vectors result) low (list umin vmin) length
		  (list (- umax umin) (- vmax vmin))))
    result))
