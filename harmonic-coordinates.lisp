(in-package :cl-user)

(defpackage :harmonic
  (:use :common-lisp :cffi)
  (:export :make-harmonic-coordinates
           :free-harmonic-coordinates
           :with-harmonic-coordinates
           :harmonic-coordinates))

(in-package :harmonic)

(define-foreign-library libharmonic
  (t (:default "libharmonic")))

(push "/home/salvi/project/rust/harmonic/target/" cffi:*foreign-library-directories*)
;;; (load-foreign-library "libharmonic.so")
(use-foreign-library libharmonic)

(defcfun harmonic-init :pointer
  (size :unsigned-int) (points :pointer) (levels :unsigned-int) (epsilon :double))
(defcfun harmonic-free :void
  (map :pointer))
(defcfun harmonic-eval :int
  (map :pointer) (point :pointer) (value :pointer))

(defun make-harmonic-coordinates (points levels)
  "Returns a list of harmonic maps, one for each vertex."
  'TODO)

(defun free-harmonic-coordinates (maps)
  (when maps
    (unwind-protect (free-harmonic-coordinates (rest maps))
      (harmonic-free (first maps)))))

(defun harmonic-coordinates (maps p)
  (mapcar (lambda (m)
            'TODO)
          maps))

(defmacro with-harmonic-coordinates ((var points &key (levels 7)) &body body)
  `(let ((,var (make-harmonic-coordinates ,points ,levels)))
     (unwind-protect (progn ,@body)
       (free-harmonic-coordinates ,var))))

