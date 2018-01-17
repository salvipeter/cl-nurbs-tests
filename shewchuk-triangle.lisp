(in-package :cl-user)

;;; http://www.cs.cmu.edu/~quake/triangle.html

(defpackage :shewchuk-triangle
  (:use :common-lisp)
  (:export :*program-path* :mesh))

(in-package :shewchuk-triangle)

(defparameter *program-path* "/home/salvi/project/cl-nurbs/tests/shewchuk/triangle")

(defun write-poly (stream points)
  (let ((n (length points)))
    (format stream "~a 2 0 0~%" n)      ; # of vertices, dimension, # of attributes, boundary marker
    (loop for i upfrom 1 as p in points do (format stream "~a~{ ~f~}~%" i (elt points (1- i))))
    (format stream "~a 0~%" n)          ; # of segments, boundary marker
    (format stream "1 ~a 1~%" n)        ; 1st vertex (from n to 1)
    (loop for i from 1 below n do (format stream "~a ~a ~a~%" (1+ i) i (1+ i)))
    (format stream "0~%")))             ; # of holes

(defun read-node-ele (node ele)
  (let ((vertices (make-array 0 :adjustable t))
        (triangles (make-array 0 :adjustable t)))
    (with-open-file (s node)
      (let ((n (read s)))
        (adjust-array vertices n)
        (read s) (read s) (read s)      ; skip 3 values
        (dotimes (i n)
          (read s)                      ; skip vertex number
          (let* ((x (read s))
                 (y (read s)))
            (setf (aref vertices i) (list x y))))))
    (with-open-file (s ele)
      (let ((n (read s)))
        (adjust-array triangles n)
        (read s) (read s)               ; skip 2 values
        (dotimes (i n)
          (read s)                      ; skip triangle number
          (let* ((a (1- (read s)))
                 (b (1- (read s)))
                 (c (1- (read s))))
            (setf (aref triangles i) (list a b c))))))
    (list vertices triangles)))

(defun mesh (points max-triangle-area)
  "POINTS is an ordered list of 2D points defining a polygon.
Returns a list of two vectors: (VERTICES TRIANGLES),
the latter is a triples of indices, starting from 0."
  (cl-fad:with-open-temporary-file (s :template "mesh-%.poly")
    (write-poly s points)
    (finish-output s)
    (let ((name (namestring (pathname s))))
      (external-program:run *program-path* (list (format nil "-pa~fDBPIQ" max-triangle-area) name))
      (let* ((root (subseq name 0 (- (length name) 5)))
             (node (format nil "~a.node" root))
             (ele (format nil "~a.ele" root)))
        (unwind-protect (read-node-ele node ele)
          (delete-file node)
          (delete-file ele))))))
