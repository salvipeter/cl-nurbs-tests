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

(defun continued-fraction (x n)
  (iter (repeat n)
        (for i = (floor x))
        (collect i)
        (setf x (/ (- x i)))))

(defun eval-continued-fraction (f)
  (cond ((null f) 0)
        ((null (rest f)) (first f))
        (t (+ (first f) (/ (eval-continued-fraction (rest f)))))))

(defun best-rational-approximations (x n)
  "Returns a list of best rational approximations of X.
Each element is a pair of the approximation and its error.
N is the length of the continued fraction used in the algorithm."
  (iter (with f = (continued-fraction x n))
        (for j from 1 to n)
        (for fr = (reverse (subseq f 0 j)))
        (for p = (first fr))
        (for ok = 
             (or (oddp p)
                 (> (abs (- x (eval-continued-fraction (reverse (rest fr)))))
                    (abs (- x (eval-continued-fraction (reverse (cons (/ p 2) (rest fr)))))))))
        (appending
         (iter (for i from (+ (ceiling p 2) (if ok 0 1)) to p)
               (for y = (eval-continued-fraction (reverse (cons i (rest fr)))))
               (collect (cons y (abs (- x y))))))))
