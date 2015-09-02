(in-package :cl-nurbs-tests)

;;; BARYCENTRIC-COORDINATES => see n-sided-ribbon.lisp
;;; BARYCENTRIC-S, BARYCENTRIC-D => see gregory-d-test.lisp
(defun test (n degree &key (position 'center) (use-d t) (weight-type 'everywhere))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (ecase position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))))
         (l (barycentric-coordinates points p))
         (half-low (floor degree 2))
         (half-up (ceiling degree 2)))
    (- 1
       (iter (for i from 0 below n)
             (for i-1 = (mod (1- i) n))
             (for i+1 = (mod (1+ i) n))
             (for si = (barycentric-s l i))
             (for si-1 = (barycentric-s l i-1))
             (for si+1 = (barycentric-s l i+1))
             (for di = (barycentric-d l i))
             (for di-1 = (barycentric-d l i-1))
             (for di+1 = (barycentric-d l i+1))
             (for alpha =
                  (if use-d
                      (/ di-1 (+ di-1 di))
                      (/ si (+ si (- 1 si-1)))))
             (for beta =
                  (if use-d
                      (/ di+1 (+ di+1 di))
                      (/ (- 1 si) (+ (- 1 si) si+1))))
             (for blf-sum = 0)
             (iter (for row from 0 below (ceiling degree 2))
                   (iter (for col from 0 to degree)
                         (for blend = (* (bernstein degree row di)
                                         (bernstein degree col si)))
                         (for mu =
                              (ecase weight-type
                                (everywhere
                                 (cond ((and (evenp degree) (= col half-low)) 1)
                                       ((< col half-up) alpha)
                                       ((> col half-low) beta)))
                                (all-corners
                                 'todo)
                                (corners-with-halving
                                 'todo)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

;;; 
#+nil
(iter (for dp in '(t nil))
      (format t "Using ~:[S~;D~] in alpha/beta:~%" dp)
      (iter (for type in '(center edge-center-mid vertex-center-mid))
            (format t "~a:~%" type)
            (iter (for n from 3 to 8)
                  (format t "~a sides:~%" n)
                  (iter (for d from 1 to 7)
                        (format t "deg: ~a => ~a%~%"
                                d (round (* (test n d :position type :use-d dp) 100)))))))
