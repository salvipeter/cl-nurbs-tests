;;; -*- mode: lisp; syntax: common-lisp -*-

(require :ieee-floats)
(in-package :cl-nurbs)

(defun read-double (stream)
  (let ((result 0))
    (dolist (x '(0 8 16 24 32 40 48 56))
      (setf (ldb (byte 8 x) result) (read-byte stream)))
    (ieee-floats:decode-float64 result)))

(defun write-double (stream d)
  (let ((64bit (ieee-floats:encode-float64 (coerce d 'double-float))))
    (dolist (c (mapcar (lambda (x) (ldb (byte 8 x) 64bit))
		       '(0 8 16 24 32 40 48 56)))
      (write-byte c stream))))

(defun read-unsigned-integer (stream)
  (let ((result 0))
    (dolist (x '(0 8 16 24))
      (setf (ldb (byte 8 x) result) (read-byte stream)))
    result))

(defun write-unsigned-integer (stream n)
  (dolist (c (mapcar (lambda (x) (ldb (byte 8 x) n)) '(0 8 16 24)))
    (write-byte c stream)))

(defun read-bss (filename)
  (with-open-file (s filename :element-type '(unsigned-byte 8))
    (let* ((degree-u (read-unsigned-integer s))
	   (nr-knots-u (read-unsigned-integer s))
	   (cyclic-u (read-unsigned-integer s))
	   (knots-u (iter (for i from 0 below nr-knots-u)
			  (collect (read-double s))))
	   (degree-v (read-unsigned-integer s))
	   (nr-knots-v (read-unsigned-integer s))
	   (cyclic-v (read-unsigned-integer s))
	   (knots-v (iter (for i from 0 below nr-knots-v)
			  (collect (read-double s))))
	   (nr-cpts-u (read-unsigned-integer s))
	   (nr-cpts-v (read-unsigned-integer s))
	   (net (make-array (list nr-cpts-u nr-cpts-v))))
      (declare (ignore cyclic-u cyclic-v))
      (dotimes (u nr-cpts-u)
	(dotimes (v nr-cpts-v)
	  (let* ((x (read-double s))
		 (y (read-double s))
		 (z (read-double s)))
	    (setf (aref net u v) (list x y z)))))
      (make-bspline-surface (list degree-u degree-v)
			    (list knots-u knots-v)
			    net))))

(defun write-bss-stream (stream surface)
  (with-accessors ((degrees degrees)
		     (knots knot-vectors)
		     (net control-net))
	surface
      (write-unsigned-integer stream (first degrees))
      (write-unsigned-integer stream (length (first knots)))
      (write-unsigned-integer stream 0)			; cyclic u
      (dolist (knot (coerce (first knots) 'list))
	(write-double stream knot))
      (write-unsigned-integer stream (second degrees))
      (write-unsigned-integer stream (length (second knots)))
      (write-unsigned-integer stream 0)			; cyclic v
      (dolist (knot (coerce (second knots) 'list))
	(write-double stream knot))
      (write-unsigned-integer stream (array-dimension net 0))
      (write-unsigned-integer stream (array-dimension net 1))
      (dotimes (u (array-dimension net 0))
	(dotimes (v (array-dimension net 1))
	  (dolist (coordinate (aref net u v))
	    (write-double stream coordinate))))))

(defun write-bss (surface filename)
  (with-open-file (s filename :direction :output :if-exists :supersede
		     :element-type '(unsigned-byte 8))
    (write-bss-stream s surface)))

(defun write-rdn (surface-list filename)
  (with-open-file (s filename :direction :output :if-exists :supersede
		     :element-type '(unsigned-byte 8))
    (dolist (surface surface-list)
      (write-unsigned-integer s 50)	; B-Spline Surface
      (write-unsigned-integer s 0)	; Length of INFO string
      ;; INFO string goes here
      (let ((nrknots (mapcar #'length (knot-vectors surface)))
	    (nrcpts (apply #'* (array-dimensions (control-net surface)))))
	;; Calculate the length
	(write-unsigned-integer s (+ 4 4 4 (* (first nrknots) 8)
				     4 4 4 (* (second nrknots) 8)
				     4 4 (* nrcpts 3 8))))
      (write-bss-stream s surface))))
