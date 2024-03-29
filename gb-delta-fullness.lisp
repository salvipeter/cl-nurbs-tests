(in-package :cl-nurbs-tests)

(defun fullness-height-circle (n)
  (let* ((alpha (/ (* 2 pi) n))
         (beta (- (/ pi 2) alpha))
         (OQ (cos (/ alpha 2)))
         (QC (* OQ (sin alpha)))
         (RC (* QC (/ (sin beta) (1+ (cos beta)))))
         (TC (* OQ (- 1 (cos alpha))))
         (RT (+ TC RC)))
    (- RT OQ)))

(defun fullness-height-parabola (n)
  (let* ((alpha (/ (* 2 pi) n)) 
         (OQ (cos (/ alpha 2))))
    (- (* OQ (cos alpha) 0.5))))

(defparameter *fullness-height-function* #'fullness-height-circle)

(defun fullness-height (n)
  (funcall *fullness-height-function* n))

(defun find-fullness-delta (n)
  (let* ((x (fullness-height n))
         (points (points-from-angles (cons (+ 270 (/ 180 n)) (rest (uniform-angles n)))))
         (l (barycentric-coordinates points (list 0 x)))
         (h (let ((*barycentric-dilation* 0))
              (barycentric-d l 0))))
    (flet ((li (i) (elt l (mod i n))))
      (/ (- h 0.5) (li -2) (li 1)
         (- 1 (li -2) (li -1) (li 0) (li 1))))))

(defun find-fullness-delta-triangle (n)
  (assert (= n 3))
  (let* ((x (fullness-height n))
         (points (points-from-angles (cons (+ 270 (/ 180 n)) (rest (uniform-angles n)))))
         (l (barycentric-coordinates points (list 0 x)))
         (h (let ((*barycentric-dilation* 0))
              (barycentric-d l 0))))
    (flet ((li (i) (elt l (mod i n))))
      (/ (- 0.5 h) (- 1 (li 1)) (li -1) (li 0) (li 1)))))

(defun find-point-on-05-with-h (points h)
  (bisection-search-root (lambda (x)
                           (- (barycentric-d (barycentric-coordinates points `(0 ,x)) 0) h))
                         -1 0.0 100.0))

#+nil
(defun extended-cp-ratio (n d)
  (let* ((*barycentric-dilation* (find-fullness-delta n))
         (points (points-from-angles (cons (+ 270 (/ 180 n)) (rest (uniform-angles n)))))
         (h (- 0.5 (/ d)))
         (h-cp (barycentric-d (barycentric-coordinates points '(0 0)) 0))
         (p05 (find-point-on-05-with-h points 0.5))
         (ph  (find-point-on-05-with-h points h))
         (pcp (find-point-on-05-with-h points h-cp)))
    (/ (- p05 ph) (- pcp ph))))

(defun extended-cp-ratio (n d)
  (let* ((*barycentric-dilation* (find-fullness-delta n))
         (points (points-from-angles (cons (+ 270 (/ 180 n)) (rest (uniform-angles n)))))
         (h (- 0.5 (/ d)))
         (h-cp (barycentric-d (barycentric-coordinates points '(0 0)) 0)))
    (/ (- 0.5 h) (- h-cp h))))
