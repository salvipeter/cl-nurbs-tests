;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs-tests-asd
  (:use :cl :asdf))

(in-package :cl-nurbs-tests-asd)

(asdf:defsystem :cl-nurbs-tests
  :author "Peter Salvi"
  :licence "Private"
  :description "Research test environment based on CL-NURBS."
  :depends-on (:iterate :cffi :cl-nurbs :ieee-floats :gsll :antik)
  :components ((:file "package")
	       (:file "matrix")
	       (:file "lu" :depends-on ("matrix"))
	       (:file "kobbelt" :depends-on ("lu"))
	       (:file "three-curve" :depends-on ("package"))
	       (:file "three-curve-iterative" :depends-on ("three-curve"))
	       (:file "geomagic-formats" :depends-on ("package"))
	       (:file "octree" :depends-on ("package"))
	       (:file "stl-format" :depends-on ("geomagic-formats" "octree"))
	       (:file "parameterization" :depends-on ("package"))
	       (:file "file-preparation" :depends-on ("parameterization"))
	       (:file "extensions" :depends-on ("parameterization"))
	       (:file "projection" :depends-on ("matrix" "package"))
	       (:file "closest-point" :depends-on ("matrix" "package"))
	       (:file "g1-zap" :depends-on ("parameterization" "matrix" "lu"))
	       (:file "g2" :depends-on ("g1-zap"))
	       (:file "g-krr" :depends-on ("g2"))
	       (:file "fair-surface-from-corners" :depends-on ("package"))
	       (:file "five-surface" :depends-on ("geomagic-formats"
						  "g1-zap"
						  "projection"
						  "closest-point"
						  "fair-surface-from-corners"))
	       (:file "kobbelt-fairing" :depends-on ("kobbelt" "five-surface"))
	       (:file "five-surface-iterative"
		      :depends-on ("geomagic-formats"))
	       (:file "edge-blend-fairing"
		      :depends-on ("g-krr" "kobbelt" "five-surface"))
	       (:file "n-patch" :depends-on ("g-krr" "kobbelt"))
               (:file "n-sided-paper" :depends-on ("package"))
               (:file "blends" :depends-on ("n-sided-paper"))
               (:file "kato-test" :depends-on ("blends"))
               (:file "ribbons" :depends-on ("kato-test"))
               (:file "patches" :depends-on ("ribbons"))
               (:file "barycentric-parameterization" :depends-on ("package"))))
