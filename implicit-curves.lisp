(in-package :cl-nurbs-tests)

;;; SLICED-CONCAVE-DISTANCE-FUNCTION-TEST in concave-domain.lisp

(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (p1 '(-0.5 0))
      (p2 '(0.5 0)))
  (flet ((foo (p) (- (point-distance p p1) (point-distance p p2))))
    (sliced-concave-distance-function-test points #'foo "/tmp/pt-pt-minus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list p1 p2))))

(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (p1 '(-0.5 0))
      (p2 '(0.5 0)))
  (flet ((foo (p) (+ (point-distance p p1) (point-distance p p2))))
    (sliced-concave-distance-function-test points #'foo "/tmp/pt-pt-plus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list p1 p2))))

(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (p1 '(-0.5 0))
      (l1 '((0.5 -1) (0.5 1))))
  (flet ((foo (p) (- (point-line-distance p l1) (point-distance p p1))))
    (sliced-concave-distance-function-test points #'foo "/tmp/pt-ln-minus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list p1 l1))))

(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (p1 '(-0.5 0))
      (l1 '((0.5 -1) (0.5 1))))
  (flet ((foo (p) (+ (point-line-distance p l1) (point-distance p p1))))
    (sliced-concave-distance-function-test points #'foo "/tmp/pt-ln-plus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list p1 l1))))


(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (l1 '((-0.7 -1) (0.3 1)))
      (l2 '((0.4 -1) (-0.5 1))))
  (flet ((foo (p) (- (point-line-distance p l1) (point-line-distance p l2))))
    (sliced-concave-distance-function-test points #'foo "/tmp/ln-ln-minus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list l1 l2))))

(let ((points '((-1 -1) (1 -1) (1 1) (-1 1)))
      (l1 '((-0.7 -1) (0.3 1)))
      (l2 '((0.4 -1) (-0.5 1))))
  (flet ((foo (p) (+ (point-line-distance p l1) (point-line-distance p l2))))
    (sliced-concave-distance-function-test points #'foo "/tmp/ln-ln-plus.ps"
                                           :resolution 0.0001 :density 0.1
                                           :elements (list l1 l2))))
