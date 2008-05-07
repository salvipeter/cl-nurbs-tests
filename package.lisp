;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs-tests
  (:nicknames :n-test)
  (:use :common-lisp :iterate :cffi :fff :cl-nurbs)
  (:import-from :cl-nurbs
		:safe-/
		:blend-function
		:sequence->double-array
		:bspline-curve-from-gcf
		:uniform-parameter-points
		:uniform-parameter-points-2d
		:bspline-surface-from-sf))
