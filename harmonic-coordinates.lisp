(in-package :cl-user)

(defpackage :harmonic
  (:use :common-lisp :cffi)
  (:export :make-harmonic-coordinates
           :free-harmonic-coordinates
           :with-harmonic-coordinates
           :harmonic-coordinates))

(in-package :harmonic)

(define-foreign-library harmonic-lib
  (t "/home/salvi/project/rust/harmonic/target/libharmonic.so"))
(load-foreign-library 'harmonic-lib)

(defcfun harmonic-init :pointer
  (size :unsigned-int) (points :pointer) (levels :unsigned-int) (epsilon :double))
(defcfun harmonic-free :void
  (map :pointer))
(defcfun harmonic-write-ppm :void
  (map :pointer) (filename :string))
(defcfun harmonic-eval :int
  (map :pointer) (point :pointer) (value :pointer))

(defun list->double-array (lst)
  (foreign-alloc :double :initial-contents (mapcar (lambda (x) (float x 1d0)) lst)))

(defun make-harmonic-coordinates (points levels epsilon)
  "Returns a list of harmonic maps, one for each vertex."
  (let ((n (length points))
        (arr (list->double-array (loop for p in points appending (append p '(0d0))))))
    (unwind-protect
         (loop for i from 0 below n collect
              (prog2
                  (setf (mem-aref arr :double (+ (* 3 i) 2)) 1d0)
                  (harmonic-init n arr levels epsilon)
                (setf (mem-aref arr :double (+ (* 3 i) 2)) 0d0)))
      (foreign-free arr))))

(defun free-harmonic-coordinates (maps)
  (when maps
    (unwind-protect (free-harmonic-coordinates (rest maps))
      (harmonic-free (first maps)))))

(defun harmonic-coordinates (maps p)
  (let ((p-arr (list->double-array p))
        (result (foreign-alloc :double)))
    (unwind-protect
         (mapcar (lambda (m)
                   (unless (zerop (harmonic-eval m p-arr result))
                     (mem-ref result :double)))
                 maps)
      (foreign-free result)
      (foreign-free p-arr))))

(defmacro with-harmonic-coordinates ((var points &key (levels 7) (epsilon 1d-5)) &body body)
  `(let ((,var (make-harmonic-coordinates ,points ,levels ,epsilon)))
     (unwind-protect (progn ,@body)
       (free-harmonic-coordinates ,var))))

#+nil
(let ((points '((0 0) (6 0) (6 6) (3.5 6) (3 3) (2.5 6) (0 6))))
  (with-harmonic-coordinates (m points :levels 7)
    (harmonic-coordinates m '(2 2))))
;;; =>
;;; levels=7: (0.363, 0.152, 0.019, 0.010, 0.273, 0.044, 0.120)
;;; levels=9: (0.368, 0.152, 0.019, 0.009, 0.269, 0.043, 0.122)
