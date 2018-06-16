(in-package :cl-user)

(defpackage :maxent
  (:use :common-lisp :cffi)
  (:export :make-maxent-coordinates
           :free-maxent-coordinates
           :with-maxent-coordinates
           :maxent-coordinates))

(in-package :maxent)

(define-foreign-library maxent-lib
  (t "/home/salvi/project/mec/libmec.so"))
(load-foreign-library 'maxent-lib)

(defcfun mec-init :pointer
  (n :int) (points :pointer))
(defcfun mec-eval :void
  (mec :pointer) (point :pointer) (coordinates :pointer))
(defcfun mec-free :void
  (mec :pointer))

(defun list->double-array (lst)
  (foreign-alloc :double :initial-contents (mapcar (lambda (x) (float x 1d0)) lst)))

(defun make-maxent-coordinates (points)
  (let ((array (list->double-array (reduce #'append points))))
    (unwind-protect
         (let ((n (length points)))
           (cons (mec-init n array) n))
      (foreign-free array))))

(defun free-maxent-coordinates (mec)
  (mec-free (car mec)))

(defun maxent-coordinates (mec p)
  (let ((p-arr (list->double-array p))
        (result (foreign-alloc :double :count (cdr mec))))
    (unwind-protect
         (progn
           (mec-eval (car mec) p-arr result)
           (loop for i from 0 below (cdr mec) 
              collect (mem-aref result :double i))) 
      (foreign-free result)
      (foreign-free p-arr))))

(defmacro with-maxent-coordinates ((var points) &body body)
  `(let ((,var (make-maxent-coordinates ,points)))
     (unwind-protect (progn ,@body)
       (free-maxent-coordinates ,var))))

#+nil
(let ((points '((0 0) (6 0) (6 6) (3.5 6) (3 3) (2.5 6) (0 6))))
  (with-maxent-coordinates (m points)
    (maxent-coordinates m '(2 2))))
;;; =>
;;; (0.426 0.170 0.036 0.051 0.140 0.066 0.111)
