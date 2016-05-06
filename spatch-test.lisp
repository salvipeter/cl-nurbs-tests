(in-package :cl-nurbs-tests)

;;; S-patch-like frame vs. Bezier
;;; Result: does not interpolate tangentially

(defun evaluate-part-with-bezier (sd)
  (destructuring-bind (s d) sd
    (let ((p '(0 0 0)))
      (iter (for i from 0 below 4)
            (for bi = (bernstein 3 i s))
            (iter (for j from 0 below 2)
                  (for bj = (bernstein 3 j d))
                  (setf p (v+ p (v* (elt (elt *coords* j) i) bi bj)))))
      p)))

(defun evaluate-part-with-lambda (l)
  
  (let ((li-2 (elt l 0))
        (li-1 (elt l 1))
        (li   (elt l 2))
        (li+1 (elt l 3)))
    (v+ (v* (elt (elt *coords* 0) 0) (expt li-1 3))
        (v* (elt (elt *coords* 0) 1) (* 3 (expt li-1 2) li))
        (v* (elt (elt *coords* 0) 2) (* 3 li-1 (expt li 2)))
        (v* (elt (elt *coords* 0) 3) (expt li 3))
        (v* (elt (elt *coords* 1) 0) (* 3 li-2 (expt li-1 2)))
        (v* (elt (elt *coords* 1) 1) (* 9 li-2 li-1 li))
        (v* (elt (elt *coords* 1) 2) (* 9 li-1 li li+1))
        (v* (elt (elt *coords* 1) 3) (* 3 (expt li 2) li+1)))))

(defparameter *coords*
  '(((0 0 0) (1 0 1) (2 0 2) (3 0 1))
    ((-1 1 0) (0.5 1 0.8) (1.7 1.3 1.4) (2.9 1 0.8))))

(let ((*resolution* 100)
      (*barycentric-type* 'wachspress)
      (*barycentric-normalized-p* t)
      (step 1d-4)
      (angles #+nil'(60 20 110 60 110) 5))
  (let* ((points (if (listp angles)
                     (points-from-angles angles)
                     (points-from-angles (uniform-angles angles))))
         (lines (lines-from-points points))
         (center (central-point points lines t))
         (line (third lines))
         (direction (v- (second line) (first line)))
         (normal (vnormalize (list (- (second direction)) (first direction)))))
    (when (< (scalar-product (v- center (first line)) normal) 0)
      (setf normal (v* normal -1)))    
    (iter (for i from 0 below *resolution*)
          (for u = (/ i (1- *resolution*)))
          (for domain-point = (line-point line u))
          (for inner-domain-point = (v+ domain-point (v* normal step)))
          (for p-b = (evaluate-part-with-bezier (list u 0)))
          (for q-b = (evaluate-part-with-bezier (list u step)))
          (for p-l = (evaluate-part-with-lambda
                      (barycentric-coordinates points domain-point)))
          (for q-l = (evaluate-part-with-lambda
                      (barycentric-coordinates points inner-domain-point)))
          (when (< (abs (- u 0.5)) 1.0d-2)
            (format t "B: ~a~%" (vnormalize (v- q-b p-b)))
            (format t "L: ~a~%" (vnormalize (v- q-l p-l))))
          (for derivative = (vnormalize (bezier (first *coords*) u 1)))
          (for diff =
               (safe-acos
                (scalar-product
                 (vnormalize (cross-product derivative (v- q-b p-b)))
                 (vnormalize (cross-product derivative (v- q-l p-l))))))
          (format t "~f => ~a~%" u (* (/ 180 pi) diff))
          (maximize #+nil(point-distance p-l p-b) (* (/ 180 pi) diff)))))
