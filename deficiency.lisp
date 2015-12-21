(in-package :cl-nurbs-tests)

(defun deficiency (n degree &key (position 'center) (use-d t))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
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
                      (if (< (+ di-1 di) *epsilon*)
                          0.5
                          (/ di-1 (+ di-1 di)))
                      (if (< (+ si (- 1 si-1)) *epsilon*)
                          0.5
                          (/ si (+ si (- 1 si-1))))))
             (for beta =
                  (if use-d
                      (if (< (+ di+1 di) *epsilon*)
                          0.5
                          (/ di+1 (+ di+1 di)))
                      (if (< (+ (- 1 si) si+1) *epsilon*)
                          0.5
                          (/ (- 1 si) (+ (- 1 si) si+1)))))
             (for blf-sum = 0)
             (iter (for row from 0 below (ceiling degree 2))
                   (iter (for col from 0 to degree)
                         (for blend = (* (bernstein degree row di)
                                         (bernstein degree col si)))
                         (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                         ((and (< row 2) (> col (- degree 2))) beta)
                                         ((or (< col row) (> col (- degree row))) 0)
                                         ((or (= col row) (= col (- degree row))) 1/2)
                                         (t 1)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

#+nil
(iter (for dp in '(t nil))
      (format t "Using ~:[S~;D~] in alpha/beta:~%" dp)
      (iter (for type in '(center edge-center-mid vertex-center-mid))
            (format t "~a:~%" type)
            (iter (for n from 3 to 8)
                  (format t "~a sides:~%" n)
                  (iter (for d from 1 to 7)
                        (format t "deg: ~a => ~a%~%"
                                d (round (* (deficiency n d :position type :use-d dp) 100)))))))

(defun innermost-axis-weight-sum (n degree &key (position 'center))
  "Weight sum of the innermost axis control points."
  (assert (> degree 4))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (iter (for i from 0 below n)
          (for si = (barycentric-s l i))
          (for di = (barycentric-d l i))
          (for j = (1- (ceiling degree 2)))
          (sum (* (bernstein degree j di) (bernstein degree j si))))))

(defun inner-weight-sum (n degree &key (position 'center))
  "Weight sum of all inner control points."
  (assert (> degree 4))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (iter (for i from 0 below n)
          (for si = (barycentric-s l i))
          (for di = (barycentric-d l i))
          (for blf-sum = 0)
          (iter (for row from 2 below (ceiling degree 2))
                (iter (for col from row to (- degree row))
                      (for blend = (* (bernstein degree row di)
                                      (bernstein degree col si)))
                      (if (or (= col row) (= col (- degree row)))
                          (incf blf-sum (/ blend 2))
                          (incf blf-sum blend))))
          (sum blf-sum))))

#+nil
(let ((n 3)
      (*resolution* 50))
  (iter (for degree from 5 to 8)
        (iter (for p in (vertices (points-from-angles (uniform-angles n))))
              (for d = (deficiency n degree :position p :use-d t))
              (for w = (inner-weight-sum n degree :position p))
              ;(for w = (innermost-axis-weight-sum n degree :position p))
              (minimizing (+ w d) into min-wd)
              (maximizing (+ w d) into max-wd)
              (finally (format t "Degree ~a: [~5f, ~5f]~%" degree min-wd max-wd)))))



;;; AutoWachspress hack

(defun barycentric-d-autowp (l i)
  (let* ((s (barycentric-s l i))
         (d (barycentric-d l i))
         (w (elt *auto-wachspress-weights* i)))
    (* d (1+ (* (- 1 s) s (- 1 d) d
                (/ (* 4 (- *auto-wachspress-central-d* w))
                   (* (- 1 w) w w)))))))

(defun deficiency-autowp (n degree &key (position 'center) (use-d t))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (- 1
       (iter (for i from 0 below n)
             (for i-1 = (mod (1- i) n))
             (for i+1 = (mod (1+ i) n))
             (for si = (barycentric-s l i))
             (for si-1 = (barycentric-s l i-1))
             (for si+1 = (barycentric-s l i+1))
             (for di = (barycentric-d-autowp l i))
             (for di-1 = (barycentric-d-autowp l i-1))
             (for di+1 = (barycentric-d-autowp l i+1))
             (for alpha =
                  (if use-d
                      (if (< (+ di-1 di) *epsilon*)
                          0.5
                          (/ di-1 (+ di-1 di)))
                      (if (< (+ si (- 1 si-1)) *epsilon*)
                          0.5
                          (/ si (+ si (- 1 si-1))))))
             (for beta =
                  (if use-d
                      (if (< (+ di+1 di) *epsilon*)
                          0.5
                          (/ di+1 (+ di+1 di)))
                      (if (< (+ (- 1 si) si+1) *epsilon*)
                          0.5
                          (/ (- 1 si) (+ (- 1 si) si+1)))))
             (for blf-sum = 0)
             (iter (for row from 0 below (ceiling degree 2))
                   (iter (for col from 0 to degree)
                         (for blend = (* (bernstein degree row di)
                                         (bernstein degree col si)))
                         (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                         ((and (< row 2) (> col (- degree 2))) beta)
                                         ((or (< col row) (> col (- degree row))) 0)
                                         ((or (= col row) (= col (- degree row))) 1/2)
                                         (t 1)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

#+nil
(let ((*auto-wachspress-central-d* 0.6)
      (*resolution* 50)
      (n 3))
  (iter (for degree from 1 to 12)
        (for *auto-wachspress-weights* = (make-list n :initial-element (/ (- n 2) n)))
        (for d = (iter (for p in (vertices (points-from-angles (uniform-angles n))))
                       (minimizing (deficiency-autowp n degree :position p :use-d t))))
        (when (< d -1d-14)
          (warn "alpha: ~a, degree ~a: ~a" *auto-wachspress-central-d* degree d))))
