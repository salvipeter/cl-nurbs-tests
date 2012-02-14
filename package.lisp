;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs-tests
  (:nicknames :n-test)
  (:use :common-lisp :iterate :cffi #+fff :fff :cl-nurbs)
  (:import-from :cl-nurbs
		:safe-/
		:blend-function
		#+fff :sequence->double-array
		#+fff :bspline-curve-from-gcf
		:uniform-parameter-points
		:uniform-parameter-points-2d
		#+fff :bspline-surface-from-sf))
