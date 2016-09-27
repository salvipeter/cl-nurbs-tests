(in-package :cl-nurbs-tests)

(defun generate-domain (n)
  (points-from-angles (cons (- 270 (/ 180 n)) (rest (uniform-angles n)))))

(defun tangent-error (n)
  (let* ((points (generate-domain n))
         (p (affine-combine (first points) 1/2 (second points)))
         (q (v+ p (list 0 *epsilon*)))
         (l (barycentric-coordinates points q)))
    (abs (/ (- (barycentric-d l 0) 1/2) *epsilon*))))

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

(defun find-optimal-delta (n &key (iterations 100) (min 0.0d0) (max 20.0d0))
  (flet ((f (x)
           (let ((*barycentric-dilation* x))
             (tangent-error n))))
    (golden-search #'f min max iterations)))
