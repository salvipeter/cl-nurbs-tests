(in-package :cl-nurbs-tests)

(defun wachspress-coordinates (points p)
  (let ((lengths (mapcar (lambda (x) (point-distance p x)) points))
        (n (length points)))
    (labels ((heron% (a b c)
               (let ((s (/ (+ a b c) 2)))
                 (sqrt (* s (- s a) (- s b) (- s c)))))
             (heron (a b c)
               (let ((a (heron% a b c)))
                 (if (complexp a) 0 a)))
             (inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (area-a (i)
               (heron (elt lengths i)
                      (elt lengths (inc i))
                      (point-distance (elt points i)
                                      (elt points (inc i)))))
             (area-c (i)
               (heron (point-distance (elt points (dec i))
				      (elt points i))
		      (point-distance (elt points (dec i))
				      (elt points (inc i)))
		      (point-distance (elt points i)
				      (elt points (inc i))))))
      (let* ((A (iter (for i from 0 below n) (collect (area-a i))))
             (w (iter (for i from 0 below n)
                      (for Ci = (area-c i))
                      (collect (reduce #'* (cons Ci (iter (for j from 0 below n)
                                                          (unless (or (= j i) (= j (dec i)))
                                                            (collect (elt A j)))))))))
             (wsum (reduce #'+ w)))
        (mapcar (lambda (wi) (/ wi wsum)) w)))))

(defun wachspress-convert-point-from-regular (points reg-p)
  (let* ((reg-points (domain-from-curves points 'regular))
         (w (wachspress-coordinates reg-points reg-p)))
    (reduce #'v+ (mapcar #'v* points w))))

(defun wachspress-convert-point-to-regular (points p)
  (let* ((w (wachspress-coordinates points p))
         (reg-points (domain-from-curves points 'regular)))
    (reduce #'v+ (mapcar #'v* reg-points w))))

(defun wachspress-converted-central-point (points)
  (wachspress-convert-point-from-regular points '(0 0)))

(defmacro defconverted-distance (distance)
  "Warning: parameters are evaluated multiple times."
  `(defmethod compute-distance ((type (eql ',(intern (format nil "CONVERTED-~:@(~a~)" distance))))
				points segments p dir)
     (let* ((w (wachspress-coordinates points p))
            (reg-points (domain-from-curves points 'regular))
            (reg-p (reduce #'v+ (mapcar #'v* reg-points w)))
            (i (position (first segments) points :test #'equal))
            (n (length points))
            (reg-segments (iter (for j from 0 below 4)
                                (for i+j = (mod (+ i j) n))
                                (collect (elt reg-points i+j)))))
       (compute-distance ',distance reg-points reg-segments reg-p dir))))

(defconverted-distance bilinear)
(defconverted-distance line-sweep-mod)

;;; Interconnected
;;; --------------

#+nil
(let ((points (points-from-angles '(40 20 60 100 80))))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/inter-conv1.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-line-sweep-mod :color nil
   :other-center (wachspress-converted-central-point points))
  (vectorized-distance-function-test
   points '(nil sd nil nil nil) "/tmp/inter-conv2.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-line-sweep-mod :color nil
   :other-center (wachspress-converted-central-point points))
  (vectorized-distance-function-test
   points '(nil nil sd nil nil) "/tmp/inter-conv3.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-line-sweep-mod :color nil
   :other-center (wachspress-converted-central-point points)))

;;; Compare with

#+nil
(let ((points (points-from-angles '(40 20 60 100 80))))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/inter1.ps"
   :resolution 0.001d0 :density 6 :distance-type 'line-sweep-mod :color nil)
  (vectorized-distance-function-test
   points '(nil sd nil nil nil) "/tmp/inter2.ps"
   :resolution 0.001d0 :density 6 :distance-type 'line-sweep-mod :color nil)
  (vectorized-distance-function-test
   points '(nil nil sd nil nil) "/tmp/inter3.ps"
   :resolution 0.001d0 :density 6 :distance-type 'line-sweep-mod :color nil))

;;; And with

#+nil
(let ((points (domain-from-curves '(40 20 60 100 80) 'regular)))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/inter-reg.ps"
   :resolution 0.001d0 :density 6 :distance-type 'line-sweep-mod :color nil))

;;; Bilinear
;;; --------

#+nil
(let ((points (points-from-angles '(40 20 60 100 80))))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/bilin-conv1.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-bilinear :color nil
   :other-center (wachspress-converted-central-point points))
  (vectorized-distance-function-test
   points '(nil sd nil nil nil) "/tmp/bilin-conv2.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-bilinear :color nil
   :other-center (wachspress-converted-central-point points))
  (vectorized-distance-function-test
   points '(nil nil sd nil nil) "/tmp/bilin-conv3.ps"
   :resolution 0.001d0 :density 6 :distance-type 'converted-bilinear :color nil
   :other-center (wachspress-converted-central-point points)))

;;; Compare with

#+nil
(let ((points (points-from-angles '(40 20 60 100 80))))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/bilin1.ps"
   :resolution 0.001d0 :density 6 :distance-type 'bilinear :color nil)
  (vectorized-distance-function-test
   points '(nil sd nil nil nil) "/tmp/bilin2.ps"
   :resolution 0.001d0 :density 6 :distance-type 'bilinear :color nil)
  (vectorized-distance-function-test
   points '(nil nil sd nil nil) "/tmp/bilin3.ps"
   :resolution 0.001d0 :density 6 :distance-type 'bilinear :color nil))

;;; And with

#+nil
(let ((points (domain-from-curves '(40 20 60 100 80) 'regular)))
  (vectorized-distance-function-test
   points '(sd nil nil nil nil) "/tmp/bilin-reg.ps"
   :resolution 0.001d0 :density 6 :distance-type 'line-sweep-mod :color nil))
