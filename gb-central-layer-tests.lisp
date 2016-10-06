(in-package :n-test)

(defun deficiency-squared-central-layer (n degree &key (position 'center) (use-d t))
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
                           #+nil(format t "~f	<-	~a,~a,~a~%" (* mu blend) i col row)
                           (incf blf-sum (* mu blend))))
               (when (evenp degree)
                 (incf blf-sum
                       (* (bernstein degree (/ degree 2) di)
                          (bernstein degree (/ degree 2) si)
                          (/ 4))))
               (sum blf-sum))))))

(defun deficiency-squared-central-layer-nondiagonal (n degree &key (position 'center) (use-d t))
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
               (iter (with layers = (ceiling degree 2))
                     (for row from 0 below layers)
                     (iter (for col from 0 to degree)
                           (for blend = (* (bernstein degree row di)
                                           (bernstein degree col si)))
                           (for mu = (cond ((< col layers) alpha)
                                           ((> col (- degree layers)) beta)
                                           (t 1)))
                           (incf blf-sum (* mu blend))))
               (when (evenp degree)
                 (incf blf-sum
                       (* (bernstein degree (/ degree 2) di)
                          (bernstein degree (/ degree 2) si)
                          (/ 4))))
               (sum blf-sum))))))

#+nil
(let ((*resolution* 2))
  (format t "Squared alpha/beta, central layer:~%|n|d|original|plus layer|negative|monotonic|~%")
  (iter (for n from 5 to 8)
        (iter (for d from 3 to 10)
              (format t "|~d|~d|~6,3f|~6,3f|~{~6,3f -> ~6,3f~}|~{~6,3f -> ~6,3f~}|~%"
                      n d
                      (let ((*barycentric-dilation* 0.0))
                        (deficiency n d))
                      (let ((*barycentric-dilation* 0.0))
                        (deficiency-squared-central-layer n d))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared-central-layer))
                                 (find-dilation-negative-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency n d))))
                      (let ((x (let ((*deficiency-function* #'deficiency-squared-central-layer))
                                 (find-dilation-monotone-boundary n d 0.0 30.0))))
                        (let ((*barycentric-dilation* x))
                          (list x (deficiency n d))))))))
