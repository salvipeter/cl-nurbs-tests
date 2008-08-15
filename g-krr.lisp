;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun g-fair-krr (surface iteration &optional (g 1))
  "Warning: destructively changes surface."
  (iter (with index = (1+ g))
	(with (n m) = (array-dimensions (control-net surface)))
	(with temp = (make-array `(,n ,m)))
	(repeat iteration)
	(iter (for i from index below (- n index))
	      (iter (for j from index below (- m index))
		    (setf (aref temp i j)
			  (bss-knot-removal-reinsertion surface `(,i ,j)))))
	(iter (for i from index below (- n index))
	      (iter (for j from index below (- m index))
		    (setf (aref (control-net surface) i j) (aref temp i j)))))
  surface)

(defun g-fair-krr-additive (surface iteration &optional (g 1))
  "Warning: destructively changes surface."
  (iter (with index = (1+ g))
	(with (n m) = (array-dimensions (control-net surface)))
	(repeat iteration)
	(iter (for i from index below (- n index))
	      (iter (for j from index below (- m index))
		    (setf (aref (control-net surface) i j)
			  (bss-knot-removal-reinsertion surface `(,i ,j))))))
  surface)

;;; For using it with FAIR-AND-FIT-XNODE
(defun fair-xnode-with-krr (xnode resolution iteration max-deviation)
  (declare (ignore max-deviation))
  (sample-surface (g-fair-krr-additive (first xnode) iteration) resolution))

(defun fair-g2-halved (node &key (halving 'inner)
		       (resolution '(100 100)) (iteration 100))
  (let ((fn (ecase halving
	      (none #'identity)
	      (end #'halve-end-intervals)
	      (inner #'halve-inner-intervals)
	      (all #'halve-all-intervals))))
    (g-fair-krr-additive
     (apply #'ensure-g2-continuity (funcall fn (first node))
	    (append (rest node) (list resolution)))
     iteration 2)))
