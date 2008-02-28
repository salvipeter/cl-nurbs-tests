;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun plane-line-intersection (la lb p0 p1 p2)
  "Returns the intersection point of the line defined by the points LA and LB,
and the plane defined by the points P0, P1 and P2."
  (macrolet ((x (sexp) `(first ,sexp))
	     (y (sexp) `(second ,sexp))
	     (z (sexp) `(third ,sexp)))
    (let ((m (make-array
	      '(3 3) :initial-contents
	      `((,(- (x la) (x lb)) ,(- (x p1) (x p0)) ,(- (x p2) (x p0)))
		(,(- (y la) (y lb)) ,(- (y p1) (y p0)) ,(- (y p2) (y p0)))
		(,(- (z la) (z lb)) ,(- (z p1) (z p0)) ,(- (z p2) (z p0))))))
	  (v (make-array
	      '(3 1) :initial-contents
	      `((,(- (x la) (x p0)))
		(,(- (y la) (y p0)))
		(,(- (z la) (z p0)))))))
      (let ((result (matrix-multiplication (matrix-inverse-3x3 m) v)))
	(v+ la (v* (v- lb la) (aref result 0 0)))))))

(defun closest-point (p p-dir q q-dir)
  "Returns the closest pair of points of the two lines."
  (flet ((closest (intersector in-plane dir)
	   (plane-line-intersection
	    (first intersector) (v+ (first intersector) (second intersector))
	    (first in-plane) (v+ (first in-plane) (second in-plane))
	    (v+ (first in-plane) dir))))
    (let ((perp (cross-product p-dir q-dir)))
      (list (closest (list p p-dir) (list q q-dir) perp)
	    (closest (list q q-dir) (list p p-dir) perp)))))
