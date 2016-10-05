(in-package :cl-nurbs-tests)

(defun fullness-height (n)
  "Circle version"
  (let* ((alpha (/ (* 2 pi) n))
         (beta (- (/ pi 2) alpha))
         (OQ (cos (/ alpha 2)))
         (QC (* OQ (sin alpha)))
         (RC (* QC (/ (sin beta) (1+ (cos beta)))))
         (TC (* OQ (- 1 (cos alpha))))
         (RT (+ TC RC)))
    (- RT OQ)))

(defun fullness-height (n)
  "Parabola version"
  (let* ((alpha (/ (* 2 pi) n)) 
         (OQ (cos (/ alpha 2))))
    (- (* OQ (cos alpha) 0.5))))

(defun find-fullness-delta (n)
  (let* ((x (fullness-height n))
         (points (points-from-angles (cons (+ 270 (/ 180 n)) (rest (uniform-angles n)))))
         (l (barycentric-coordinates points (list 0 x)))
         (h (let ((*barycentric-dilation* 0))
              (barycentric-d l 0))))
    (flet ((li (i) (elt l (mod i n))))
      (/ (- h 0.5) (li -2) (li 1)
         (- 1 (li -2) (li -1) (li 0) (li 1))))))
