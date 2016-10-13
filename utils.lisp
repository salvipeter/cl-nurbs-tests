(in-package :cl-nurbs-tests)

(defmacro dprint (&rest vars)
  `(progn
     ,@(mapcar (lambda (var)
                 `(format t "~a: ~a~%" ',var ,var))
               vars)))

(defun bisection-search-root (fn min max iterations)
  (let ((mid (/ (+ min max) 2)))
    (if (zerop iterations)
        mid
        (let ((x (funcall fn mid)))
          (cond ((< (abs x) 1.0d-6) mid)
                ((< x 0) (bisection-search-root fn mid max (1- iterations)))
                (t (bisection-search-root fn min mid (1- iterations))))))))

(defun golden-search (f min max iterations)
  (let* ((phi (/ (1+ (sqrt 5)) 2))
         (d (* (- phi 1) (- max min)))
         (x1 (+ min d))
         (x2 (- max d))
         (f1 (funcall f x1))
         (f2 (funcall f x2)))
    (iter (with x)
          (repeat iterations)
          (if (< f1 f2)
              (setf x x1
                    min x2
                    x1 (+ x2 (* (1- phi) (- max x2)))
                    x2 x
                    f2 f1
                    f1 (funcall f x1))
              (setf x x2
                    max x1
                    x2 (- x1 (* (1- phi) (- x1 min)))
                    x1 x
                    f1 f2
                    f2 (funcall f x2)))
          (finally (return x)))))
