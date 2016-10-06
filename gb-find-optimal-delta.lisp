(in-package :cl-nurbs-tests)

(defun generate-domain (n)
  (points-from-angles (cons (- 270 (/ 180 n)) (rest (uniform-angles n)))))

(defun tangent-error (n)
  (let* ((points (generate-domain n))
         (p (affine-combine (first points) 1/2 (second points)))
         (q (v+ p (list 0 *epsilon*)))
         (l (barycentric-coordinates points q)))
    (abs (/ (- (barycentric-d l 0) 1/2) *epsilon*))))

(defun find-optimal-delta (n &key (iterations 100) (min -20.0d0) (max 20.0d0))
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

(defun rational-approximations-below (x denom)
  (iter (for n upfrom 2)
        (for lst = (best-rational-approximations x n))
        (for d = (denominator (caar (last lst))))
        (while (< d denom))
        (finally (return (remove-if-not (lambda (x) (<= (denominator (car x)) denom)) lst)))))
