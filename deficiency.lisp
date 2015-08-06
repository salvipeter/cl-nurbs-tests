(in-package :cl-nurbs-tests)

(defun test (n degree ptype)
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (ecase ptype
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))))
         (l (barycentric-coordinates points p)))
    (- 1
       (iter (for i from 0 below n)
             (for i-1 = (mod (1- i) n))
             (for i+1 = (mod (1+ i) n))
             (for si = (barycentric-s l i))
             (for di = (barycentric-d l i))
             (for di-1 = (barycentric-d l i-1))
             (for di+1 = (barycentric-d l i+1))
             (for alpha = (/ di-1 (+ di-1 di)))
             (for beta = (/ di+1 (+ di+1 di)))
             (for blf-sum = 0)
             (iter (for j from 0 to degree)
                   (for limit = (if (or (< j 2) (> j (- degree 2)))
                                    1
                                    (- degree 2)))
                   (iter (for k from 0 to limit)
                         (incf blf-sum
                               (* (bernstein degree j si)
                                  (bernstein degree k di)
                                  (cond ((< j 2) alpha)
                                        ((> j (- degree 2)) beta)
                                        ((= k 2)
                                         (ecase degree
                                           (4 1/2)
                                           (5 7/10)
                                           (6 4/5)))
                                        ((= k 3)
                                         (ecase degree
                                           (5 3/10)
                                           (6 1/2)))
                                        ((= k 4)
                                         (ecase degree
                                           (6 1/5)))
                                        (t 1))))))
             (sum blf-sum)))))

;;; (format t "~5f~%" (test 6 4 'center))
