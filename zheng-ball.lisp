(in-package :cl-nurbs-test)

;;; Test implementation of
;;; J.J. Zheng, A.A. Ball: Control point surfaces over non-four-sided areas (1996)

(defun boundary-index-p (l)
  "An index is a boundary index if one of its values is zero."
  (find 0 l))

(defun all-shifts (lst)
  "E.g. (A B C) => ((A B C) (B C A) (C A B))."
  (labels ((rec (a b)
             (if (null a)
                 '()
                 (cons (append a b)
                       (rec (rest a) (append b (list (first a))))))))
    (rec lst '())))

(defun all-indices (n m)
  "Parameters:
N: # of sides
M: degree"
  (let ((base (iter (for l2 from 0 to (floor m 2))
                    (appending
                     (iter (for l3 from 0 to (floor m 2))
                           (collect (append (list (- m l3) l2 l3 (- m l2))
                                            (make-list (- n 4)
                                                       :initial-element (- m (min l2 l3))))))))))
    (delete-duplicates (mapcan #'all-shifts base) :test #'equal)))

(defun zb-blend (n m l u)
  "As in Eq. (4.5). Parameters:
N: # of sides
M: degree
L: control point index (a list of N values)
U: parametric point (a list of N values)"
  )

(defun zb-blend-without-deficiency (n m l u)
  "As in Eq. (4.13). Parameters:
N: # of sides
M: degree
L: control point index (a list of N values)
U: parametric point (a list of N values)"
  (if (boundary-index-p l)
      (zb-blend n m l u)
      (let ((deficiency (1- (reduce #'+ (lambda (l)
                                          (zb-blend n m l u))
                                    (all-indices n m)))))
        (- (zb-blend n m l u)
           (/ deficiency
              (if (even m)
                  (1+ (/ (* n (- m 2) m) 4))
                  (/ (* n (1- m) (1- m)) 4)))))))


;;; Parameterization tests

;;; 5-sided parameters a la Sabin, with c5 = 1:

(defparameter *points*
  (iter (for u1 from 0 to 1 by 1/100)
        (appending
         (iter (for u3 from (- 1 u1) to 1 by 1/100)
               (cond ((zerop u1)
                      (collect (list 0 0 1 1 1))
                      (collect (list 0 1 1 1 0)))
                     ((zerop u3)
                      (collect (list 1 0 0 1 1))
                      (collect (list 1 1 0 0 1)))
                     (t (collect (list u1
                                       (/ (1- (+ u1 u3)) (* u1 u3))
                                       u3
                                       (/ (- 1 u1) u3)
                                       (/ (- 1 u3) u1)))))))))


(defparameter *points2*
  (append
   (iter (for u1 from 0 to 1 by 1/10)
         (appending
          (iter (for u2 from 0 to 1 by 1/10)
                (unless (= u1 u2 1)
                  (for u4 = (- 1 (* u1 u2)))
                  (for u5 = (/ (- 1 u2) u4))
                  (for u3 = (- 1 (* u1 u5)))
                  (collect (list u1 u2 u3 u4 u5))))))
   (iter (for u3 from 0 to 1 by 1/10)
         (collect (list 1 1 u3 0 (- 1 u3))))))

