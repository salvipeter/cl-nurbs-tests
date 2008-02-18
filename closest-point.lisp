;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun matrix-multiplication (m1 m2)
  (let ((n1 (array-dimension m1 0))
	(n (array-dimension m1 1))
	(n2 (array-dimension m2 1)))
    (assert (= n (array-dimension m2 0)) (m1 m2)
	    "The two matrices are incompatible.")
    (let ((result (make-array `(,n1 ,n2))))
      (dotimes (i n1)
	(dotimes (j n2)
	  (setf (aref result i j)
		(iter (for k from 0 below n)
		      (sum (* (aref m1 i k) (aref m2 k j)))))))
      result)))

(defun matrix-inverse-3x3 (m)
  (let ((result (make-array '(3 3)))
	(det (+ (* (aref m 0 0) (- (* (aref m 2 2) (aref m 1 1))
				   (* (aref m 2 1) (aref m 1 2))))
		(- (* (aref m 1 0) (- (* (aref m 2 2) (aref m 0 1))
				      (* (aref m 2 1) (aref m 0 2)))))
		(* (aref m 2 0) (- (* (aref m 1 2) (aref m 0 1))
				   (* (aref m 1 1) (aref m 0 2)))))))
    (macrolet ((store (i j neg a1 a2 b1 b2 c1 c2 d1 d2)
		 `(setf (aref result ,i ,j)
			(/ (- (* (aref m ,a1 ,a2) (aref m ,b1 ,b2))
			      (* (aref m ,c1 ,c2) (aref m ,d1 ,d2)))
			   ,(if neg -1 1) det))))
      (store 0 0 nil 2 2 1 1 2 1 1 2)
      (store 0 1  t  2 2 0 1 2 1 0 2)
      (store 0 2 nil 1 2 0 1 1 1 0 2)
      (store 1 0  t  2 2 1 0 2 0 1 2)
      (store 1 1 nil 2 2 0 0 2 0 0 2)
      (store 1 2  t  1 2 0 0 1 0 0 2)
      (store 2 0 nil 2 1 1 0 2 0 1 1)
      (store 2 1  t  2 1 0 0 2 0 0 1)
      (store 2 2 nil 1 1 0 0 1 0 0 1))
    result))

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
