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
                                         ((< row 2) 1)
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

(defmethod compute-distance ((type (eql 'bary-autowp)) points segments p dir)
  (let* ((i (position (elt segments 2) points :test #'equal))
         (l (barycentric-coordinates points p)))
    (if (eq dir 's)
        (barycentric-s l i)
        (barycentric-d-autowp l i))))

(defparameter *gb-diagonal-weight* 1)
(defparameter *gb-side-weight* 1/2)

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
                                         ((or (= col row) (= col (- degree row))) *gb-diagonal-weight*)
                                         ((< row 2) 1)
                                         (t *gb-side-weight*)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

#+nil
(let ((degree 3)
      (*resolution* 100)
      (n 3))
  (iter (for *auto-wachspress-central-d* from 0.2 to 0.8 by 0.01)
        (for *auto-wachspress-weights* = (make-list n :initial-element (/ (- n 2) n)))
        (for d = (iter (for p in (vertices (points-from-angles (uniform-angles n))))
                       (minimizing (deficiency-autowp n degree :position p :use-d t))))
        (when (> d -1d-5)
          (format t "~a => ~a~%" *auto-wachspress-central-d* d))))

(defun bisection-search-root (fn min max iterations)
  (let ((mid (/ (+ min max) 2)))
    (if (zerop iterations)
        mid
        (let ((x (funcall fn mid)))
          (cond ((< (abs x) 1.0d-6) mid)
                ((< x 0) (bisection-search-root fn mid max (1- iterations)))
                (t (bisection-search-root fn min mid (1- iterations))))))))

(defun find-autowp-for-deficiency (n d &key (target 0.0) (iterations 100))
  (flet ((f (x)
           (let ((*auto-wachspress-central-d* x)
                 (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
             (- (deficiency-autowp n d) target))))
    (bisection-search-root #'f 0.0 1.0 iterations)))

;;; Find autowachspress settings for a given target deficiency
#+nil
(iter (for n from 3 to 8)
      (iter (for d from 3 to 10)
            (for x = (find-autowp-for-deficiency n d :target 0.1))
            (format t "n=~d,	d=~d	=>	~9,6f~%" n d x)))

;;; Check also if the patch has negative deficiency somewhere
#+nil
(let ((*resolution* 100))
  (iter (for n from 3 to 8)
        (iter (for d from 3 to 10)
              (for x = (find-autowp-for-deficiency n d :target 0.0))
              (format t "n=~d,	d=~d	=>	~9,6f~%" n d x)
              (for min =
                   (let ((*auto-wachspress-central-d* x)
                         (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
                     (iter (for p in (vertices (points-from-angles (uniform-angles n))))
                           (minimizing (deficiency-autowp n d)))))
              (when (< min -1d-5)
                (format t "Error: ~a~%" min)))))


;;; Deficiency for MP patches

(defun deficiency-midpoint (n &key (position 'center))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (flet ((hermite (x) (+ (expt (- 1 x) 3)
                           (* 3 (expt (- 1 x) 2) x))))
      (- 1
         (iter (for i from 0 below n)
               (for i-1 = (mod (1- i) n))
               (for si = (barycentric-s l i))
               (for si-1* = (- 1 (barycentric-s l i-1)))
               (for di = (barycentric-d l i))
               (for di-1 = (barycentric-d l i-1))
               (sum (/ (+ (* di (hermite si-1*) (hermite di-1))
                          (* di-1 (hermite si) (hermite di)))
                       (+ di di-1))))))))


;;; for bary-constr

(defun deficiency1 (n degree &key (position 'center) (use-d t) (parameterization 'bary-constr))
  "Might be slower, but uses COMPUTE-PARAMETER, so this can be used with BARY-CONSTR."
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (s (compute-parameter parameterization 's points p))
         (d (compute-parameter parameterization 'd points p)))
    (- 1
       (iter (for i from 0 below n)
             (for i-1 = (mod (1- i) n))
             (for i+1 = (mod (1+ i) n))
             (for si = (elt s i))
             (for si-1 = (elt s i-1))
             (for si+1 = (elt s i+1))
             (for di = (elt d i))
             (for di-1 = (elt d i-1))
             (for di+1 = (elt d i+1))
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
                                         ((< row 2) 1)
                                         (t 1)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

(defun deficiency1-squared (n degree &key (position 'center) (use-d t) (parameterization 'bary-constr))
  "Might be slower, but uses COMPUTE-PARAMETER, so this can be used with BARY-CONSTR."
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (s (compute-parameter parameterization 's points p))
         (d (compute-parameter parameterization 'd points p)))
    (flet ((sqr (x) (* x x)))
      (- 1
         (iter (for i from 0 below n)
               (for i-1 = (mod (1- i) n))
               (for i+1 = (mod (1+ i) n))
               (for si = (elt s i))
               (for si-1 = (elt s i-1))
               (for si+1 = (elt s i+1))
               (for di = (elt d i))
               (for di-1 = (elt d i-1))
               (for di+1 = (elt d i+1))
               (for alpha =
                    (if use-d
                        (if (< (+ di-1 di) *epsilon*)
                            0.5
                            (/ (sqr di-1) (+ (sqr di-1) (sqr di))))
                        (if (< (+ si (- 1 si-1)) *epsilon*)
                            0.5
                            (/ (sqr si) (+ (sqr si) (sqr (- 1 si-1)))))))
               (for beta =
                    (if use-d
                        (if (< (+ di+1 di) *epsilon*)
                            0.5
                            (/ (sqr di+1) (+ (sqr di+1) (sqr di))))
                        (if (< (+ (- 1 si) si+1) *epsilon*)
                            0.5
                            (/ (sqr (- 1 si)) (+ (sqr (- 1 si)) (sqr si+1))))))
               (for blf-sum = 0)
               (iter (for row from 0 below (ceiling degree 2))
                     (iter (for col from 0 to degree)
                           (for blend = (* (bernstein degree row di)
                                           (bernstein degree col si)))
                           (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                           ((and (< row 2) (> col (- degree 2))) beta)
                                           ((or (< col row) (> col (- degree row))) 0)
                                           ((or (= col row) (= col (- degree row))) 1/2)
                                           ((< row 2) 1)
                                           (t 1)))
                           (incf blf-sum (* mu blend))))
               (sum blf-sum))))))


(defun deficiency-squared (n degree &key (position 'center) (use-d t))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (flet ((sqr (x) (* x x)))
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
                            (/ (sqr di-1) (+ (sqr di-1) (sqr di))))
                        (if (< (+ si (- 1 si-1)) *epsilon*)
                            0.5
                            (/ (sqr si) (+ (sqr si) (sqr (- 1 si-1)))))))
               (for beta =
                    (if use-d
                        (if (< (+ di+1 di) *epsilon*)
                            0.5
                            (/ (sqr di+1) (+ (sqr di+1) (sqr di))))
                        (if (< (+ (- 1 si) si+1) *epsilon*)
                            0.5
                            (/ (sqr (- 1 si)) (+ (sqr (- 1 si)) (sqr si+1))))))
               (for blf-sum = 0)
               (iter (for row from 0 below (ceiling degree 2))
                     (iter (for col from 0 to degree)
                           (for blend = (* (bernstein degree row di)
                                           (bernstein degree col si)))
                           (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                           ((and (< row 2) (> col (- degree 2))) beta)
                                           ((or (< col row) (> col (- degree row))) 0)
                                           ((or (= col row) (= col (- degree row))) 1/2)
                                           ((< row 2) 1)
                                           (t 1)))
                           (incf blf-sum (* mu blend))))
               (sum blf-sum))))))

(defun deficiency-squared-no0 (n degree &key (position 'center) (use-d t))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (flet ((sqr (x) (* x x)))
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
                            (/ (sqr di-1) (+ (sqr di-1) (sqr di))))
                        (if (< (+ si (- 1 si-1)) *epsilon*)
                            0.5
                            (/ (sqr si) (+ (sqr si) (sqr (- 1 si-1)))))))
               (for beta =
                    (if use-d
                        (if (< (+ di+1 di) *epsilon*)
                            0.5
                            (/ (sqr di+1) (+ (sqr di+1) (sqr di))))
                        (if (< (+ (- 1 si) si+1) *epsilon*)
                            0.5
                            (/ (sqr (- 1 si)) (+ (sqr (- 1 si)) (sqr si+1))))))
               (for blf-sum = 0)
               (iter (for row from 0 below (ceiling degree 2))
                     (iter (for col from 0 to degree)
                           (for blend = (* (bernstein degree row di)
                                           (bernstein degree col si)))
                           (for mu = (cond ((< col (/ degree 2)) alpha)
                                           ((> col (/ degree 2)) beta)
                                           (t 1)))
                           (incf blf-sum (* mu blend))))
               (sum blf-sum))))))

#+nil
(let ((*barycentric-dilation* 0))
  (format t "|n|d|normal|squared|squared-no0|squared-constr|~%")
  (iter (for n from 4 to 8)
        (iter (for d from 3 to 10)
              (format t "|~d|~d|~,5f|~,5f|~,5f|~,5f|~%"
                      n d
                      (deficiency n d)
                      (deficiency-squared n d)
                      (deficiency-squared-no0 n d)
                      (deficiency1-squared n d)))))

#+nil
(let ((*barycentric-dilation* 0)
      (*resolution* 40))
  (format t "|n|d|normal|squared|squared-no0|squared-constr|~%")
  (iter (for n from 5 to 8)
        (iter (for d from 3 to 10)
              (format t "|~d|~d|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~%"
                      n d
                      (let ((x (let ((*deficiency-function* #'deficiency))
                                 (find-dilation-negative-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared))
                                 (find-dilation-negative-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency-squared n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared-no0))
                                 (find-dilation-negative-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency-squared-no0 n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency1-squared))
                                 (find-dilation-negative-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency1-squared n d))))))))

#+nil
(let ((*barycentric-dilation* 0)
      (*resolution* 40))
  (format t "|n|d|normal|squared|squared-no0|squared-constr|~%")
  (iter (for n from 5 to 8)
        (iter (for d from 3 to 10)
              (format t "|~d|~d|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~%"
                      n d
                      (let ((x (let ((*deficiency-function* #'deficiency))
                                 (find-dilation-monotone-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared))
                                 (find-dilation-monotone-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency-squared n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared-no0))
                                 (find-dilation-monotone-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency-squared-no0 n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency1-squared))
                                 (find-dilation-monotone-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency1-squared n d))))))))
