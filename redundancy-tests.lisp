;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

(defparameter *curve*
  (make-bspline-curve 3 '(0.0 0.0 0.0 0.0 2.84792 3.34905 3.651261 4.178889
			  6.369409 6.369409 6.369409 6.369409)
		      '((-1.0 -1.0) (0.0 1.0) (1.0 1.0) (1.3 1.0) (1.6 1.1)
			(1.9 1.0) (3.0 1.0) (4.0 0.0))))
(defvar *faired*)
(write-ps *curve* #p"results/original.ps" 100 :margin '(20 130))
(write-rbn (bsc-extrude *curve* 5) "results/original.rbn")

;;; Iterative fairing
(setf *faired* (bsc-fair *curve* :simplex-iteration 10))
(write-ps *faired* #p"results/it-faired.ps" 100 :margin '(20 130))
(write-rbn (bsc-extrude *faired* 5) "results/it-faired.rbn")

;;; Integrative fairing
(let ((points (bsc-faired-polygon *curve* 100 100 100)))
  (with-open-file (s "results/in-faired.data" :direction :output)
    (format s "~{~{~d~^ ~}~%~}" (coerce points 'list)))
  (setf *faired* (bsc-fit *curve* points 0.1)))
(write-ps *faired* #p"results/in-faired.ps" 100 :margin '(20 130))
(write-rbn (bsc-extrude *faired* 5) "results/in-faired.rbn")

;;;     _
;;; A _| |_ alaku tesztgorbe
;;;
;; (defparameter *curve*
;;   (make-bspline-curve 3 '(0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 13 13 13)
;; 		      '((-4.0 -2.0) (-3.0 -2.1) (-2.0 -1.9) (-1.0 -2.0)
;; 			(-1.1 0.0) (-0.9 1.0) (-1.0 2.0) (0.0 2.1) (1.0 1.9)
;; 			(2.0 2.0) (2.1 1.0) (1.9 0.0) (2.0 -1.0) (3.0 -1.1)
;; 			(3.5 -0.9) (4.0 -1.0))))
