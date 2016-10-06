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
