(in-package :cl-user)

(defpackage :harmonic
  (:use :common-lisp :cffi)
  (:export :make-harmonic-coordinates
           :make-harmonic-sd-coordinates
           :free-harmonic-coordinates
           :with-harmonic-coordinates
           :with-harmonic-sd-coordinates
           :harmonic-coordinates
           :harmonic-sd-coordinates
           :*linear-constraint*))

(in-package :harmonic)

(define-foreign-library harmonic-lib
  (t "/home/salvi/project/rust/harmonic/target/libharmonic.so"))
(load-foreign-library 'harmonic-lib)

(defcfun harmonic-create :pointer
  (min :pointer) (max :pointer) (levels :unsigned-int))
(defcfun harmonic-add-point :void
  (map :pointer) (point :pointer))
(defcfun harmonic-add-line :void
  (map :pointer) (from :pointer) (to :pointer))
(defcfun harmonic-solve :void
  (map :pointer) (epsilon :double) (biharmonic :int))
(defcfun harmonic-free :void
  (map :pointer))
(defcfun harmonic-write-ppm :void
  (map :pointer) (filename :string))
(defcfun harmonic-eval :int
  (map :pointer) (point :pointer) (value :pointer))

(defvar *linear-constraint* t)

(defun list->double-array (lst)
  (foreign-alloc :double :initial-contents (mapcar (lambda (x) (float x 1d0)) lst)))

(defun bounding-box (points)
  (list (reduce (lambda (p q) (mapcar #'min p q)) points)
        (reduce (lambda (p q) (mapcar #'max p q)) points)))

(defun make-harmonic (points levels)
  (let ((bbox (mapcar #'list->double-array (bounding-box points))))
    (unwind-protect
         (harmonic-create (first bbox) (second bbox) levels)
      (mapcar #'foreign-free bbox))))

(defun make-harmonic-coordinates (points levels epsilon biharmonicp)
  "Returns a list of harmonic maps, one for each vertex."
  (let* ((n (length points))
         (bh (if biharmonicp 1 0)))
    (loop for i from 0 below n collect
         (let ((m (make-harmonic points levels)))
           (unless *linear-constraint*
             (let ((p (list->double-array (append (elt points i) '(1)))))
               (unwind-protect
                    (harmonic-add-point m p)
                 (foreign-free p))))
           (dotimes (j n)
             (let ((j+1 (mod (1+ j) n)))
               (unless (and (not *linear-constraint*) (or (= j i) (= j+1 i)))
                 (let ((from (list->double-array (append (elt points j) (if (= j i) '(1) '(0)))))
                       (to (list->double-array (append (elt points j+1) (if (= j+1 i) '(1) '(0))))))
                   (unwind-protect
                        (harmonic-add-line m from to)
                     (foreign-free from)
                     (foreign-free to))))))
           (harmonic-solve m epsilon bh)
           m))))

(defun make-harmonic-sd-coordinates (points levels epsilon biharmonicp)
  "Returns a list of harmonic maps, one for each vertex."
  (let* ((n (length points))
         (bh (if biharmonicp 1 0)))
    (loop for i from 0 below n append
         (list
          (let* ((m (make-harmonic points levels))
                 (i-2 (mod (- i 2) n))
                 (i-1 (mod (1- i) n))
                 (i+1 (mod (1+ i) n))
                 (p-2 (list->double-array (append (elt points i-2) '(0))))
                 (p-1 (list->double-array (append (elt points i-1) '(0))))
                 (p   (list->double-array (append (elt points i  ) '(1))))
                 (p+1 (list->double-array (append (elt points i+1) '(1)))))
            (unwind-protect
                 (progn
                   (harmonic-add-line m p-2 p-1)
                   (harmonic-add-line m p-1 p  )
                   (harmonic-add-line m p   p+1))
              (foreign-free p-2)
              (foreign-free p-1)
              (foreign-free p  )
              (foreign-free p+1))
            (harmonic-solve m epsilon bh)
            m)
          (let ((m (make-harmonic points levels))
                (i-1 (mod (1- i) n)))
            (dotimes (j n)
              (let ((j+1 (mod (1+ j) n)))
                (unless (and (not *linear-constraint*) (or (= j+1 i-1) (= j i)))
                  (let ((from (list->double-array
                               (append (elt points j) (if (or (= j i-1) (= j i)) '(0) '(1)))))
                        (to (list->double-array
                             (append (elt points j+1) (if (or (= j+1 i-1) (= j+1 i)) '(0) '(1))))))
                    (unwind-protect
                         (harmonic-add-line m from to)
                      (foreign-free from)
                      (foreign-free to))))))
            (harmonic-solve m epsilon bh)
            m)))))

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

(defun harmonic-sd-coordinates (maps p)
  (let ((h (harmonic-coordinates maps p)))
    (loop for i from 0 below (length h) by 2 collect
         (list (elt h i) (elt h (1+ i))))))

(defmacro with-harmonic-coordinates ((var points &key (levels 9) (epsilon 1d-5) biharmonicp)
                                     &body body)
  `(let ((,var (make-harmonic-coordinates ,points ,levels ,epsilon ,biharmonicp)))
     (unwind-protect (progn ,@body)
       (free-harmonic-coordinates ,var))))

(defmacro with-harmonic-sd-coordinates ((var points &key (levels 9) (epsilon 1d-5) biharmonicp)
                                        &body body)
  `(let ((,var (make-harmonic-sd-coordinates ,points ,levels ,epsilon ,biharmonicp)))
     (unwind-protect (progn ,@body)
       (free-harmonic-coordinates ,var))))

#+nil
(let ((points '((0 0) (6 0) (6 6) (3.5 6) (3 3) (2.5 6) (0 6))))
  (with-harmonic-coordinates (m points :levels 7)
    (harmonic-coordinates m '(2 2))))
;;; =>
;;; levels=7: (0.363, 0.152, 0.019, 0.010, 0.273, 0.044, 0.120)
;;; levels=9: (0.368, 0.152, 0.019, 0.009, 0.269, 0.043, 0.122)
