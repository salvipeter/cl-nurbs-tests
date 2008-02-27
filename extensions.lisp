;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun bsc-unclamp-end (curve lst)
  "Unclamps the end of the curve so that the knot vector ends with LST."
  (let* ((knots (knot-vector curve))
	 (p (degree curve))
	 (points (control-points curve))
	 (n (length points))
	 (new-knots (concatenate 'vector (subseq knots 0 n) lst))
	 (new-points (make-array n :initial-contents points)))
    (dotimes (i (1- p))
      (iter (for j from i downto 0)
	    (for alpha = (/ (- (elt new-knots n) (elt new-knots (- n j 1)))
			    (- (elt new-knots (+ (- n j) i 1))
			       (elt new-knots (- n j 1)))))
	    (setf (elt new-points (- n j 1))
		  (v* (v- (elt new-points (- n j 1))
			  (v* (elt new-points (- n j 2)) (- 1.0 alpha)))
		      (/ alpha)))))
    (make-bspline-curve p new-knots new-points)))

(defun bsc-extend-end (curve new-knot)
  (let* ((upper (bsc-upper-parameter curve))
	 (degree (degree curve))
	 (new-end (cons upper (iter (repeat degree) (collect new-knot))))
	 (unclamped (bsc-unclamp-end curve new-end))
	 (inserted (bsc-insert-knot unclamped upper)))
    (bsc-unclamp-end inserted (iter (repeat (1+ degree)) (collect new-knot)))))

(defun bsc-extend (curve percent)
  "Extends a curve at both sides by PERCENT%."
  (let* ((lower (bsc-lower-parameter curve))
	 (upper (bsc-upper-parameter curve))
	 (new-knot (+ lower (* (* percent 0.01) (- upper lower)))))
    (bsc-reverse-parameterization
     (bsc-extend-end
      (bsc-reverse-parameterization
       (bsc-shift-parameterization
	(bsc-extend-end curve (+ upper new-knot))
	(- new-knot)))
      (+ upper new-knot)))))

(defun bss-unclamp-u-end (surface lst)
  "Unclamps the end of the surface so that the U knot vector ends with LST."
  (let* ((knots (first (knot-vectors surface)))
	 (p (first (degrees surface)))
	 (net (control-net surface))
	 (n (array-dimension net 0))
	 (m (array-dimension net 1))
	 (new-knots (concatenate 'vector (subseq knots 0 n) lst))
	 (new-net (make-array (list n m))))
    (dotimes (j m)
      (dotimes (i n)
	(setf (aref new-net i j) (aref net i j))))
    (dotimes (column m)
      (dotimes (i (1- p))
	(iter (for j from i downto 0)
	      (for alpha = (/ (- (elt new-knots n) (elt new-knots (- n j 1)))
			      (- (elt new-knots (+ (- n j) i 1))
				 (elt new-knots (- n j 1)))))
	      (setf (aref new-net (- n j 1) column)
		    (v* (v- (aref new-net (- n j 1) column)
			    (v* (aref new-net (- n j 2) column) (- 1.0 alpha)))
			(/ alpha))))))
    (make-bspline-surface (list p (second (degrees surface)))
			  (list new-knots
				(copy-seq (second (knot-vectors surface))))
			  new-net)))

(defun bss-extend-u-end (surface new-knot)
  (let* ((upper (first (bss-upper-parameter surface)))
	 (degree (first (degrees surface)))
	 (new-end (cons upper (iter (repeat degree) (collect new-knot))))
	 (unclamped (bss-unclamp-u-end surface new-end))
	 (inserted (bss-insert-knot unclamped upper :u-direction t)))
    (bss-unclamp-u-end inserted (iter (repeat (1+ degree))
				      (collect new-knot)))))

(defun bss-extend-u (surface percent)
  "Extends the U domain of a surface at both sides by PERCENT%."
  (let* ((lower (first (bss-lower-parameter surface)))
	 (upper (first (bss-upper-parameter surface)))
	 (new-knot (+ lower (* (* percent 0.01) (- upper lower)))))
    (bss-reverse-parameterization
     (bss-extend-u-end
      (bss-reverse-parameterization
       (bss-shift-parameterization
	(bss-extend-u-end surface (+ upper new-knot))
	:u (- new-knot))
       :u t)
      (+ upper new-knot))
     :u t)))

(defun bss-extend (surface percent)
  "Extends the U and V domain of a surface at both sides by PERCENT%."
  (bss-flip-uv
   (bss-extend-u
    (bss-flip-uv
     (bss-extend-u surface percent))
    percent)))
