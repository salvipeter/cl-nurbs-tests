(in-package :cl-nurbs-tests)

(defparameter *barycentric-type* 'wachspress) ; or: meanvalue / harmonic
(defparameter *barycentric-normalized-p* t)
(defparameter *barycentric-squared-p* nil) ; if T then treats *BARYCENTRIC-NORMALIZED-P* as T
(defun barycentric-coordinates (points p)
  (let* ((vectors (mapcar (lambda (x) (v- p x)) points))
         (lengths (mapcar #'vlength vectors))
         (n (length points)))
    (labels ((inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (area (i)                  ; signed area = det(si,si+1)/2
               (let ((si (elt vectors i))
                     (si+1 (elt vectors (inc i))))
                 (/ (- (* (elt si 0) (elt si+1 1))
                       (* (elt si 1) (elt si+1 0)))
                    2.0d0)))
             (area-product (exceptions)
               (let ((r 1))
                 (iter (for j from 0 below n)
                       (unless (member j exceptions)
                         (setf r (* r (area j)))))
                 r))
             (bctype (x)
               (ecase *barycentric-type*
                 (harmonic (* x x))
                 (meanvalue x)
                 (wachspress 1))))
      (let* ((corner nil)
             (w (iter (for i from 0 below n)
                      (for Ai = (area-product (list i)))
                      (for Ai-1 = (area-product (list (dec i))))
                      (for Ai-1i = (area-product (list (dec i) i)))
                      (for Bi = (let ((si-1 (elt vectors (dec i)))
                                      (si+1 (elt vectors (inc i))))
                                  (/ (- (* (elt si-1 0) (elt si+1 1))
                                        (* (elt si-1 1) (elt si+1 0)))
                                     2.0d0)))
                      (for ri-1 = (bctype (elt lengths (dec i))))
                      (for ri = (bctype (elt lengths i)))
                      (for ri+1 = (bctype (elt lengths (inc i))))
                      (when (< ri *epsilon*)
                        (setf corner i))
                      (collect (+ (* ri-1 Ai-1)
                                  (* ri+1 Ai)
                                  (* -1 ri Bi Ai-1i)))))
             (wsum (reduce #'+ w)))
        (cond (corner (let ((lst (make-list n :initial-element 0)))
                        (setf (elt lst corner) 1)
                        lst))
              (*barycentric-squared-p*
               (let* ((w2 (mapcar (lambda (x) (* x x)) w))
                      (wsum2 (reduce #'+ w2)))
                 (mapcar (lambda (wi) (/ wi wsum2)) w2)))
              (*barycentric-normalized-p* (mapcar (lambda (wi) (/ wi wsum)) w))
              (t w))))))

(defun barycentric-s (l i)
  (let* ((n (length l))
         (i-1 (mod (1- i) n)))
    (safe-/ (elt l i)
            (+ (elt l i-1) (elt l i)))))

(defun barycentric-d-original (l i)
  (let* ((n (length l))
         (i-1 (mod (1- i) n)))
    (- 1 (elt l i-1) (elt l i))))

(defparameter *barycentric-d-function* #'barycentric-d-original)

(defun barycentric-d (l i)
  (funcall *barycentric-d-function* l i))

(defmethod compute-distance ((type (eql 'bary)) points segments p dir)
  (let* ((i (position (elt segments 2) points :test #'equal))
         (l (barycentric-coordinates points p)))
    (if (eq dir 's)
        (barycentric-s l i)
        (barycentric-d l i))))

(defmethod compute-distance ((type (eql 'bary-lambda)) points segments p dir)
  (let* ((i (position (elt segments 2) points :test #'equal))
         (l (barycentric-coordinates points p)))
    (if (eq dir 's)
        (error "No S in this distance function")
        (elt l i))))

(defmethod compute-distance ((type (eql 'bary-constr)) points segments p dir)
  "Constrained barycentric parameterization."
  (if (eq dir 's)
      (compute-distance 'bary points segments p 's)
      (let* ((s (compute-distance 'bary points segments p 's))
             (d (compute-distance 'bary points segments p 'd))
             (s-1 (compute-distance 'bary points (segments-prev points segments) p 's))
             (s+1 (compute-distance 'bary points (segments-next points segments) p 's))
             (lst (list d (- 1 s) (- 1 d) s)))
        (+ (* d (+ (ribbon-blend lst 0) (ribbon-blend lst 2)))
           (* s+1 (ribbon-blend lst 1))
           (* (- 1 s-1) (ribbon-blend lst 3))))))

#+nil
(let* ((points (points-from-angles (uniform-angles 7))))
  (vectorized-distance-function-test
   points '(s d s nil nil nil nil) "/tmp/proba.ps"
   :resolution 0.001d0 :density 20 :distance-type 'bary-constr :color t))

(defparameter *prior-edge-fn* 'kato)    ; kato / angle
(defun prior-distribution-coordinates (points p)
  (let ((n (length points)))
    (labels ((pi (i) (elt points (mod i n)))
             (edge-fn (i)
               (ecase *prior-edge-fn*
                 (kato (- (+ (point-distance p (pi i)) (point-distance p (pi (1- i))))
                          (point-distance (pi i) (pi (1- i)))))
                 (angle (+ (* (point-distance p (pi i)) (point-distance p (pi (1- i))))
                           (scalar-product (v- p (pi i)) (v- p (pi (1- i)))))))))
      (let* ((ds (iter (for i from 0 below n) (collect (edge-fn i))))
             (ps (iter (for i from 0 below n)
                       (collect (iter (for j from 0 below n)
                                      (for dj in ds)
                                      (unless (or (= j i) (= j (mod (1+ i) n)))
                                        (multiply dj))))))
             (sum (reduce #'+ ps)))
        (if (> sum *epsilon*)
            (mapcar (lambda (x) (/ x sum)) ps)
            (iter (for i from 0 below n)
                  (if (and (< (elt ds i) *epsilon*)
                           (< (elt ds (mod (1+ i) n)) *epsilon*))
                      (collect 1)
                      (collect 0))))))))

#+nil
(defun maximum-entropy-coordinates (points p)
  "BAD: exact Hessian is divergent, should use DFP (see below)."
  (let* ((n (length points))
         (dim (length (first points)))
         (prior (coerce (prior-distribution-coordinates points p) 'vector))
         (vectors (mapcar (lambda (v) (v- v p)) points)))
    (when (>= (count 0 prior :test #'=) (- n 2)) ; or maybe only at the vertices?
      (return-from maximum-entropy-coordinates prior))
    (assert (= dim 2) () "Matrix inversion now only works for 2 dimensions")
    (flet ((compute-gradients (l)
             (let ((gradient (make-array dim :element-type 'double-float :initial-element 0.0d0))
                   (hessian (make-array (list dim dim)
                                        :element-type 'double-float :initial-element 0.0d0))
                   (denom 0.0d0))
               (iter (for k from 0 below n)
                     (for vk in vectors)
                     (for me = (* (elt prior k) (exp (- (scalar-product l vk)))))
                     (incf denom me)
                     (iter (for d from 0 below dim)
                           (incf (elt gradient d) (* me (- (elt vk d))))
                           (iter (for d2 from 0 below dim)
                                 (incf (aref hessian d d2) (* me (elt vk d) (elt vk d2))))))
               (iter (for d from 0 below dim)
                     (setf (elt gradient d) (/ (elt gradient d) denom))
                     (iter (for d2 from 0 below dim)
                           (setf (aref hessian d d2) (/ (aref hessian d d2) denom))))
               (iter (for d from 0 below dim)
                     (iter (for d2 from 0 below dim)
                           (decf (aref hessian d d2) (* (elt gradient d) (elt gradient d2)))))
               (values gradient hessian)))
           (compute-coords (l)
             (let* ((coords (iter (for k from 0 below n)
                                  (for vk in vectors)
                                  (for me = (* (elt prior k) (exp (- (scalar-product l vk)))))
                                  (collect me)))
                    (sum (reduce #'+ coords)))
               (mapcar (lambda (x) (/ x sum)) coords))))
      (iter (with l = (make-list n :initial-element 0.0d0))
            (multiple-value-bind (g h)
                (compute-gradients l)
              (setf l (v- l (coerce
                             (matrix:to-vector
                              (matrix:multiplication
                               (matrix:inverse-2x2 h)
                               (matrix:from-vector g)))
                             'list)))
              (while (> (vlength (coerce g 'list)) *epsilon*))
              (finally (return (compute-coords l))))))))

(defun dfp-approximation (x x1 f f1 h)
  (flet ((a- (u v) (map 'vector #'- u v)))
    (matrix:m+
     h
     (matrix:m+
      (matrix:m*
       (matrix:multiplication
        (matrix:from-vector (a- x1 x))
        (matrix:transpose (matrix:from-vector (a- x1 x))))
       (/ (aref (matrix:multiplication
                 (matrix:transpose (matrix:from-vector (a- x1 x)))
                 (matrix:from-vector (a- f1 f)))
                0 0)))
      (matrix:m*
       (matrix:multiplication
        (matrix:multiplication
         h
         (matrix:from-vector (a- f1 f)))
        (matrix:transpose
         (matrix:multiplication
          h
          (matrix:from-vector (a- f1 f)))))
       (/ (- (aref (matrix:multiplication
                    (matrix:multiplication
                     (matrix:transpose (matrix:from-vector (a- f1 f)))
                     h)
                    (matrix:from-vector (a- f1 f)))
                   0 0))))))))

(defun maximum-entropy-coordinates (points p)
  "With the David-Fletcher-Powell approximation of the Hessian."
  (let* ((n (length points))
         (prior (coerce (prior-distribution-coordinates points p) 'vector))
         (vectors (mapcar (lambda (v) (v- v p)) points)))
    (let ((zeros (count 0 prior :test #'=)))
      (when (>= zeros (- n 2))
        (return-from maximum-entropy-coordinates
          (if (= zeros (1- n))
              prior
              (iter (for i from 0 below n)
                    (if (> (elt prior i) 0)
                        (let ((ip (mod (1+ i) n)))
                          (if (> (elt prior ip) 0)
                              (collect (/ (point-distance p (elt points ip))
                                          (point-distance (elt points i) (elt points ip))))
                              (let ((im (mod (1- i) n)))
                                (collect (/ (point-distance p (elt points im))
                                            (point-distance (elt points i) (elt points im)))))))
                        (collect 0.0d0)))))))
    (labels ((a- (u v) (map 'vector #'- u v))
             (a* (u v) (reduce #'+ (map 'list #'* u v)))
             (compute-gradient (l)
               (let ((gradient (make-array 2 :element-type 'double-float :initial-element 0.0d0))
                     (denom 0.0d0))
                 (iter (for k from 0 below n)
                       (for vk in vectors)
                       (for me = (* (elt prior k) (exp (- (a* l vk)))))
                       (incf denom me)
                       (iter (for d from 0 below 2)
                             (incf (elt gradient d) (* me (- (elt vk d))))))
                 (iter (for d from 0 below 2)
                       (setf (elt gradient d) (/ (elt gradient d) denom))) 
                 gradient))
             (compute-coords (l)
               (let* ((coords (iter (for k from 0 below n)
                                    (for vk in vectors)
                                    (for me = (* (elt prior k) (exp (- (a* l vk)))))
                                    (collect me)))
                      (sum (reduce #'+ coords)))
                 (mapcar (lambda (x) (/ x sum)) coords))))
      (iter (with l = (make-array n :initial-element 0.0d0))
            (with h = (make-array '(2 2) :element-type 'double-float
                                  :initial-contents '((1.0d0 0.0d0) (0.0d0 1.0d0))))
            (for prev-l = (copy-seq l))
            (for g first (compute-gradient l) then next-g)
            (setf l (a- l (matrix:to-vector (matrix:multiplication h (matrix:from-vector g)))))
            (for next-g = (compute-gradient l))
            (setf h (dfp-approximation prev-l l g next-g h))
            (while (> (vlength (coerce g 'list)) *epsilon*))
            (finally (return (compute-coords l)))))))
