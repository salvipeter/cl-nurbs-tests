(in-package :cl-nurbs-tests)

(defun efficient-mean-value (points values p)
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
                    2.0d0))))
      (let* ((w (iter (for i from 0 below n)
                      (for Ai = (area i))
                      (for Ai-1 = (area (dec i)))
                      (for Di = (scalar-product (elt vectors (inc i)) (elt vectors i)))
                      (for Di-1 = (scalar-product (elt vectors i) (elt vectors (dec i))))
                      (for ri-1 = (elt lengths (dec i)))
                      (for ri = (elt lengths i))
                      (for ri+1 = (elt lengths (inc i)))
                      (when (< (abs ri) *epsilon*)
                        (return-from efficient-mean-value (elt values i)))
                      (when (and (< (abs Ai) *epsilon*)
                                 (< Di (- *epsilon*)))
                        (return-from efficient-mean-value
                          (/ (+ (* (elt values i) ri+1)
                                (* (elt values (inc i)) ri))
                             (+ ri ri+1))))
                      (collect (+ (if (> (abs Ai-1) *epsilon*)
                                      (/ (- ri-1 (/ Di-1 ri)) Ai-1)
                                      0)
                                  (if (> (abs Ai) *epsilon*)
                                      (/ (- ri+1 (/ Di ri)) Ai)
                                      0)))))
             (wsum (reduce #'+ w)))
        (iter (for i from 0 below n)
              (sum (/ (* (elt values i) (elt w i)) wsum)))))))

(defun tangent-mean-value (points values p)
  (let* ((vectors (mapcar (lambda (x) (v- p x)) points))
         (lengths (mapcar #'vlength vectors))
         (n (length points)))
    (labels ((inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (tan/2 (i)
               (let* ((a (vnormalize (elt vectors i)))
                      (b (vnormalize (elt vectors (inc i))))
                      (cos-ab (scalar-product a b))
                      (sin-ab (- (* (first a) (second b)) (* (second a) (first b)))))
                 (/ (- 1.0 cos-ab) sin-ab))))
      (let* ((w (iter (for i from 0 below n)
                      (for ri = (elt lengths i))
                      (for tan1 = (tan/2 (dec i)))
                      (for tan2 = (tan/2 i))
                      (when (< (abs ri) *epsilon*)
                        (return-from tangent-mean-value (elt values i)))
                      (collect (/ (+ tan1 tan2) ri))))
             (wsum (reduce #'+ w)))
        (iter (for i from 0 below n)

              (sum (/ (* (elt values i) (elt w i)) wsum)))))))

;;; Kai-problem

#+nil
(let ((points '((1 0) (0 1) (0 0) (-1 0) (-1 2) (2 2) (2 0)))
      (p '(1 1/2)))
  (iter (for i from 0 below 7)
        (for values = (iter (for j from 0 below 7) (collect (if (= i j) 1 0))))
        (collect (efficient-mean-value points values p))))

;;; Tests

(defparameter *points* '((2 5) (3 3) (5 3) (6 5) (8 3) (7 0) (1 0) (0 3)))

;;; Peti version:
;;; 1. Generate mean value coordinates by the concave domain
;;; 2. Compute the corresonding point in a regular domain,
;;;    where concave edge-series are represented as one edge with inner vertices
;;; 3. Compute (s,d) in the regular domain

(defun peti-test (p)
  (let ((regular (points-from-angles (uniform-angles 6)))
        (quasi-regular '((0.5d0 0.8660254037844386d0)
                         (0.16666666666666666d0 0.8660254037844386d0)
                         (-0.16666666666666666d0 0.8660254037844386d0)
                         (-0.5d0 0.8660254037844387d0)
                         (-1.0d0 0)
                         (-0.5d0 -0.8660254037844384d0)
                         (0.5d0 -0.866025403784439d0)
                         (1.0d0 0))))
    (let ((u (efficient-mean-value *points* (mapcar #'first quasi-regular) p))
          (v (efficient-mean-value *points* (mapcar #'second quasi-regular) p)))
      (list (elt (compute-parameter 'mean-bilinear 's regular (list u v)) 1)
            (elt (compute-parameter 'mean-bilinear 'd regular (list u v)) 1)))))

;;; Tomi version:
;;; 1. Compute d-parameters the same way as Wachspress coordinates
;;; 2. Compute s-parameters by a similar construct

(defun tomi-test (p)
  (let ((values-s '(0 1/3 2/3 1 1 2/3 1/3 0))
        (values-d '(0 0 0 0 1 1 1 1)))
    (list (efficient-mean-value *points* values-s p)
          (efficient-mean-value *points* values-d p))))

;;; Simple evaluation

(defun bitmap-test (points sd-fun filename &key (size 400) (density 10))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%255~%" size size)
    (let ((lines (lines-from-points points)))
      (labels ((map-point (p) (v* p (/ 10.0d0 size)))
               (query (p)
                 (reduce #'min (mapcar (lambda (x)
                                         (multiple-value-bind (quot rem) (round x (/ density))
                                           (declare (ignore quot))
                                           (abs rem)))
                                       (funcall sd-fun p)))))      
        (iter (for y from (1- size) downto 0)
              (iter (for x from 0 below size)
                    (for p = (map-point (list x y)))
                    (for color = (min (floor (* (query p) density 255)) 255))
                    (format s "~d " color))
              (terpri s))))))

;;; (bitmap-test *points* #'peti-test "/tmp/peti.pgm")
;;; (bitmap-test *points* #'tomi-test "/tmp/tomi.pgm")

;;; Several other tests

(defun compute-sd (points values-s values-d)
  (lambda (p)
    (list (efficient-mean-value points values-s p)
          (efficient-mean-value points values-d p))))

(let ((p '((0 10) (5 5) (10 10) (10 0) (0 0)))
      (s '(0 1/2 1 1 0))
      (d '(0 0 0 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly1.pgm"))

(let ((p '((0 10) (5 7) (10 10) (10 0) (5 3) (0 0)))
      (s '(0 1/2 1 1 1/2 0))
      (d '(0 0 0 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly2.pgm"))

(let ((p '((0 7) (3 6) (10 10) (10 0) (3 4) (0 3)))
      (s '(0 0.28 1 1 0.28 0))
      (d '(0 0 0 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly3.pgm"))

(let ((p '((1 10) (3 6) (6 6) (7 10) (10 10) (10 0) (0 0) (0 10)))
      (s '(0 0.39 0.64 1 1 2/3 1/3 0))
      (d '(0 0 0 0 1 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly4.pgm"))

(let ((p '((1 10) (2 5) (3 5) (6 10) (10 10) (10 0) (0 0) (0 10)))
      (s '(0 0.43 0.51 1 1 2/3 1/3 0))
      (d '(0 0 0 0 1 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly5.pgm"))

(let ((p '((2 10) (4 6) (6 6) (8 10) (10 10) (10 0) (3 2) (0 0) (0 10)))
      (s '(0 0.41 0.59 1 1 0.68 0.44 0.32 0))
      (d '(0 0 0 0 1 1 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly6.pgm"))

(let ((p '((0 10) (2 10) (4 6) (6 6) (8 10) (10 10) (10 0) (3 2) (0 0)))
      (s '(0 1 1 1 1 0.91 0.48 0.16 0))
      (d '(0 0 1 1 1 1 1 1 1)))
  (bitmap-test p (compute-sd p s d) "/tmp/poly6b.pgm"))

;;; Tracing

(defun trace-mean-parameter (points sd-fun type parameter resolution)
  (let ((lines (lines-from-points points)))
    (labels ((query (p)
	       (abs (- (elt (funcall sd-fun p) (if (eq type 's) 0 1)) parameter)))
             (get-start (line)
	       (let ((k (/ (point-distance (first line) (second line)) resolution)))
		 (iter (repeat (1+ (floor k)))
		       (with d = (v* (v- (second line) (first line)) (/ k)))
		       (for p first (first line) then (v+ p d))
		       (collect (cons (query p) p) into tests)
		       (finally (return (first (sort tests #'< :key #'first)))))))
	     (get-fuzzy-start ()
	       (iter (with n = (length lines))
                     (for i from 0 below n)
		     (for line = (elt lines i))
                     (for (value . param) = (get-start line))
                     (finding param minimizing value)))
	     (sorted (actual)
	       (iter (for xi from -1 to 1)
		     (for x = (+ (first actual) (* xi resolution)))
		     (appending
		      (iter (for yi from -1 to 1)
			    (unless (= xi yi 0)
			      (for y = (+ (second actual) (* yi resolution)))
			      (collect (list (query (list x y)) x y))))
		      into lst)
		     (finally (return (mapcar #'rest (sort lst #'< :key #'first))))))
	     (best-one (prev actual)
	       (first (remove-if (lambda (p)
                                   (and prev
                                        (> (scalar-product (v- prev actual) (v- p actual)) 0)))
				 (sorted actual)))))
      (let ((start (get-fuzzy-start)))
	(cons start
	      (iter (for pprev first nil then prev)
		    (for prev first start then next)
		    (for next = (best-one pprev prev))
		    (while next)
		    (collect next)))))))

(defun trace-concave (points sd-fun filename &key (resolution 0.01d0) (density 10))
  "First k points constitute a concave `edge', for which we want to generate s,d lines."
  (flet ((map-point (p)
	   (list (* (+ (first p) 1.0d0) 60)
		 (* (+ (second p) 1.0d0) 60))))
    (let* ((n (length points)) 
	   (lines (lines-from-points points)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
	(format s "%!PS~%")
	(iter (for i from 0 below n)
	      (for line in lines)
	      (format s "% Segment: ~a~%" i)
	      (format s "2 setlinewidth~%~
                         newpath~%~
                         ~{~f ~}moveto~%~
                         ~{~f ~}lineto~%~
                         stroke~%~
                         1 setlinewidth~%"
		      (map-point (first line))
		      (map-point (second line)))
	      (iter (for type in '(s d))
		    (format s "% Type: ~a~%" type)
		    (iter (with d = (/ density))
			  (for parameter from d below 1 by d)
			  (format s "% Parameter: ~a~%" parameter)
			  (for trace =
                               (trace-mean-parameter points sd-fun type parameter resolution))
			  (format s "newpath~%")
			  (format s "~{~f ~}moveto~%~
                                     ~{~{~f ~}lineto~%~}"
				  (map-point (first trace))
				  (mapcar #'map-point (rest trace)))
			  (format s "stroke~%"))))
	(format s "showpage~%")))))

;;; (trace-concave *points* #'peti-test "/tmp/peti.ps")
;;; (trace-concave *points* #'tomi-test "/tmp/tomi.ps")
