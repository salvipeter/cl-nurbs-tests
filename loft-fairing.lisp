;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

;; (defparameter *profiles* (first (read-rdn "/cygdrive/d/tmp/profiles.rdn")))

;;; Let us assume that all curves are in the parameter range [0,1].

(defun spine-fairness-deviation (points)
  "Deviation from the line segments defined by the neighboring points."
  (let ((n (length points)))
    (iter (for i from 1 below (1- n))
	  (for prev = (elt points (1- i)))
	  (for p = (elt points i))
	  (for next = (elt points (1+ i)))
	  (for dir = (vnormalize (v- next prev)))
	  (for proj = (v+ prev (v* dir (scalar-product (v- p prev) dir))))
	  (sum ((lambda (x) (* x x)) (point-distance p proj))))))

(defun spine-fairness-length (points)
  "Length."
  (iter (for p1 in points)
	(for p2 in (rest points))
	(for diff = (v- p1 p2))
	(sum (scalar-product diff diff))))

(defun spine-fairness (points)
  (+ (* 1.0d0 (spine-fairness-deviation points))
     (* 1.0d0 (spine-fairness-length points))))

(defun fair-spine-parameters (profiles &key (sections 10) (iteration 1000) (zero-start t))
  (flet ((measure (params)
	   (iter (for section from 0 below sections)
		 (for u = (/ section sections))
		 (for points =
		      (iter (for curve in profiles)
			    (for p in params)
			    (collect (bsc-evaluate curve (mod (+ u p) 1.0d0)))))
		 (sum (spine-fairness points)))))
    (mapcar (lambda (x) (mod x 1.0d0))
	    (let ((start (loop repeat (length profiles) collect
			      (if zero-start 0.0d0 (random 1.0d0)))))
	      (downhill-simplex:minimize #'measure start iteration)))))

(defclass polyline () ((points :initarg :points :accessor polyline-points)))
(defun make-polyline (points) (make-instance 'polyline :points points))
(defmethod write-rdn-object ((line polyline) s)
  (write-unsigned-integer s 30)
  (write-unsigned-integer s 0)
  (let ((points (polyline-points line)))
    (write-unsigned-integer s (+ 4 (* 3 (length points) 8)))
    (write-unsigned-integer s (length points))
    (dolist (point (coerce points 'list))
      (dolist (coordinate point)
	(write-double s coordinate)))))

(defun canonize-spine-parameters (params)
  (let ((first (first params)))
    (mapcar (lambda (x) (mod (- x first) 1.0d0)) params)))

(defun guide-polygons (profiles params &key (sections 10))
  (iter (for section from 0 below sections)
	(for u = (/ section sections))
	(collect
	 (make-polyline
	  (iter (for curve in profiles)
		(for p in (canonize-spine-parameters params))
		(collect (bsc-evaluate curve (mod (+ u p) 1.0d0))))))))
