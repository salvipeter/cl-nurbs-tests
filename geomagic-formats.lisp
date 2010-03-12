;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

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

(defgeneric write-rdn-object (obj stream))

(defmethod write-rdn-object ((lst cons) s)
  (write-unsigned-integer s 1)
  (write-unsigned-integer s 0)		; info
  (write-unsigned-integer s 0)		; size
  (dolist (obj lst)
    (write-rdn-object obj s))
  (write-unsigned-integer s 2)
  (write-unsigned-integer s 0)		; info
  (write-unsigned-integer s 0))		; size

(defmethod write-rdn-object ((curve bspline-curve) s)
  (write-unsigned-integer s 40) ; B-Spline Curve
  (write-unsigned-integer s 0)  ; Length of INFO string
  ;; INFO string goes here
  (with-accessors ((degree degree) (knots knot-vector) (points control-points))
      curve
    (let ((nrknots (length knots))
	  (nrcpts (length points)))
      ;; Calculate the length
      (write-unsigned-integer s (+ 4 4 4 (* nrknots 8) 4 (* 3 nrcpts 8)))
      (write-unsigned-integer s degree)
      (write-unsigned-integer s nrknots)
      (write-unsigned-integer s 0)	; closedp
      (dolist (knot (coerce knots 'list))
	(write-double s knot))
      (write-unsigned-integer s nrcpts)
      (dolist (point (coerce points 'list))
	(dolist (coordinate point)
	  (write-double s coordinate))))))

(defmethod write-rdn-object ((surface bspline-surface) s)
  (write-unsigned-integer s 50)	; B-Spline Surface
  (write-unsigned-integer s 0)	; Length of INFO string
  ;; INFO string goes here
  (let ((nrknots (mapcar #'length (knot-vectors surface)))
	(nrcpts (apply #'* (array-dimensions (control-net surface)))))
    ;; Calculate the length
    (write-unsigned-integer s (+ 4 4 4 (* (first nrknots) 8)
				 4 4 4 (* (second nrknots) 8)
				 4 4 (* nrcpts 3 8))))
  (write-bss-stream s surface))

(defun write-rdn (lst filename)
  (with-open-file (s filename :direction :output :if-exists :supersede
		     :element-type '(unsigned-byte 8))
    (dolist (obj lst)
      (write-rdn-object obj s))))

(defgeneric read-rdn-object (obj stream))

(defmethod read-rdn-object ((obj (eql 'bspline-curve)) s)
  (let* ((degree (read-unsigned-integer s))
	 (k (read-unsigned-integer s))
	 (closedp (read-unsigned-integer s))
	 (knots (iter (repeat k) (collect (read-double s))))
	 (n (read-unsigned-integer s))
	 (points (iter (repeat n)
		       (collect (iter (repeat 3)
				      (collect (read-double s)))))))
    (declare (ignore closedp))
    (make-bspline-curve degree knots points)))

(defmethod read-rdn-object ((obj (eql 'bspline-surface)) s)
  (let* ((u-degree (read-unsigned-integer s))
	 (u-nknots (read-unsigned-integer s))
	 (u-closedp (read-unsigned-integer s))
	 (u-knots (iter (repeat u-nknots) (collect (read-double s))))
	 (v-degree (read-unsigned-integer s))
	 (v-nknots (read-unsigned-integer s))
	 (v-closedp (read-unsigned-integer s))
	 (v-knots (iter (repeat v-nknots) (collect (read-double s))))
	 (u-npoints (read-unsigned-integer s))
	 (v-npoints (read-unsigned-integer s))
	 (net (make-array (list u-npoints v-npoints))))
    (declare (ignore u-closedp v-closedp))
    (dotimes (u u-npoints)
      (dotimes (v v-npoints)
	(setf (aref net u v) (iter (repeat 3) (collect (read-double s))))))
    (make-bspline-surface (list u-degree v-degree)
			  (list u-knots v-knots)
			  net)))

(defun read-rdn-stream (s)
  (iter (with rdn-types = '((40 . bspline-curve) (50 . bspline-surface)))
	(for code = (handler-case (read-unsigned-integer s)
		      (end-of-file () (finish))))
	(let ((info-length (read-unsigned-integer s)))
	  (iter (repeat info-length) (read-byte s)))
	(read-unsigned-integer s)	; size
	(cond ((= code 2) (finish))	; EndGroup
	      ((= code 1)		; BeginGroup
	       (collect (read-rdn-stream s)))
	      (t (let ((type (cdr (assoc code rdn-types))))
		   (assert type (type) "Unknown type: ~a" code)
		   (collect (read-rdn-object type s)))))))

(defun read-rdn (filename)
  (with-open-file (s filename :element-type '(unsigned-byte 8))
    (read-rdn-stream s)))
