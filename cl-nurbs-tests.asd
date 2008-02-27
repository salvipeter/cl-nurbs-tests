;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs-tests-asd
  (:use :cl :asdf))

(in-package :cl-nurbs-tests-asd)

(asdf:defsystem :cl-nurbs-tests
  :author "Peter Salvi"
  :licence "Private"
  :description "Research test environment based on CL-NURBS."
  :depends-on (:iterate :cffi :fff :cl-nurbs :ieee-floats)
  :components ((:file "package")
	       (:file "three-curve" :depends-on ("package"))
	       (:file "three-curve-iterative" :depends-on ("three-curve"))
	       (:file "geomagic-formats" :depends-on ("package"))
	       (:file "parameterization" :depends-on ("package"))
	       (:file "g1-zap" :depends-on ("parameterization"))
	       (:file "file-preparation" :depends-on ("parameterization"))
	       (:file "extensions" :depends-on ("parameterization"))
	       (:file "closest-point" :depends-on ("package"))
	       (:file "five-surface" :depends-on ("geomagic-formats"
						  "g1-zap"
						  "closest-point"))
	       (:file "five-surface-iterative"
		      :depends-on ("geomagic-formats"))))
