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
