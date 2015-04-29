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

(defun peti-test (points p)
  (let ((regular (points-from-angles (uniform-angles 6)))
        (quasi-regular '((0.5d0 0.8660254037844386d0)
                         (0.16666666666666666d0 0.8660254037844386d0)
                         (-0.16666666666666666d0 0.8660254037844386d0)
                         (-0.5d0 0.8660254037844387d0)
                         (-1.0d0 0)
                         (-0.5d0 -0.8660254037844384d0)
                         (0.5d0 -0.866025403784439d0)
                         (1.0d0 0))))
    (let ((u (efficient-mean-value points (mapcar #'first quasi-regular) p))
          (v (efficient-mean-value points (mapcar #'second quasi-regular) p)))
      (list (elt (compute-parameter 'mean-bilinear 's regular (list u v)) 1)
            (elt (compute-parameter 'mean-bilinear 'd regular (list u v)) 1)))))

;;; Tomi version:
;;; 1. Compute d-parameters the same way as Wachspress coordinates
;;; 2. Compute s-parameters by a similar construct

(defun tomi-test (points p)
  (let ((values-s '(0 1/3 2/3 1 1 2/3 1/3 0))
        (values-d '(0 0 0 0 1 1 1 1)))
    (list (efficient-mean-value points values-s p)
          (efficient-mean-value points values-d p))))

;;; Simple evaluation

(defun bitmap-test (points sd-fun filename &key (object-size 2.0d0) (size 400) (density 10))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%255~%" size size)
    (labels ((map-point (p)
               (v- (v* p (/ object-size size))
                   (list (/ object-size 2) (/ object-size 2))))
             (query (p)
               (reduce #'min (mapcar (lambda (x)
                                       (multiple-value-bind (quot rem) (round x (/ density))
                                         (declare (ignore quot))
                                         (abs rem)))
                                     (funcall sd-fun points p)))))
      (iter (for y from (1- size) downto 0)
            (iter (for x from 0 below size)
                  (for p = (map-point (list x y)))
                  (for color = (min (floor (* (query p) density 255)) 255))
                  (format s "~d " color))
            (terpri s)))))

;;; (bitmap-test *points* #'peti-test "/tmp/peti.pgm" :object-size 20)
;;; (bitmap-test *points* #'tomi-test "/tmp/tomi.pgm" :object-size 20)

;;; Several other tests

(defun compute-sd (values-s values-d)
  (lambda (points p)
    (list (efficient-mean-value points values-s p)
          (efficient-mean-value points values-d p))))

(let ((p '((0 10) (5 5) (10 10) (10 0) (0 0)))
      (s '(0 1/2 1 1 0))
      (d '(0 0 0 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly1.pgm"))

(let ((p '((0 10) (5 7) (10 10) (10 0) (5 3) (0 0)))
      (s '(0 1/2 1 1 1/2 0))
      (d '(0 0 0 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly2.pgm"))

(let ((p '((0 7) (3 6) (10 10) (10 0) (3 4) (0 3)))
      (s '(0 0.28 1 1 0.28 0))
      (d '(0 0 0 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly3.pgm"))

(let ((p '((1 10) (3 6) (6 6) (7 10) (10 10) (10 0) (0 0) (0 10)))
      (s '(0 0.39 0.64 1 1 2/3 1/3 0))
      (d '(0 0 0 0 1 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly4.pgm"))

(let ((p '((1 10) (2 5) (3 5) (6 10) (10 10) (10 0) (0 0) (0 10)))
      (s '(0 0.43 0.51 1 1 2/3 1/3 0))
      (d '(0 0 0 0 1 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly5.pgm"))

(let ((p '((2 10) (4 6) (6 6) (8 10) (10 10) (10 0) (3 2) (0 0) (0 10)))
      (s '(0 0.41 0.59 1 1 0.68 0.44 0.32 0))
      (d '(0 0 0 0 1 1 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly6.pgm"))

(let ((p '((0 10) (2 10) (4 6) (6 6) (8 10) (10 10) (10 0) (3 2) (0 0)))
      (s '(0 1 1 1 1 0.91 0.48 0.16 0))
      (d '(0 0 1 1 1 1 1 1 1)))
  (bitmap-test p (compute-sd s d) "/tmp/poly6b.pgm"))

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


;;; And now for something completely different.

;;; Handling concave patches by another method - proof of concept test
(defparameter *ribbons*
  (list (list (make-bspline-curve 3 '(0.0 0.0 0.0 0.0 0.2 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.8 1.0 1.0 1.0 1.0) '((1.00 2.00 0.0) (1.80 1.60 0.0) (4.29 0.67 0.0) (5.94 6.17 0.0) (8.57 4.37 0.0) (9.00 4.00 0.0) (9.00 4.00 0.0) (9.00 4.00 0.0) (9.00 4.00 0.0) (9.00 4.00 0.0) (9.00 4.00 0.0) (8.67 3.62 0.0) (7.43 1.81 0.0) (9.53 -0.30 0.0) (11.08 0.58 0.0) (12.00 1.50 0.0)))
              (make-bspline-curve 3 '(0.0 0.0 0.0 0.0 0.2 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.8 1.0 1.0 1.0 1.0) '((1.447214 2.894427 0.000000) (1.575217 2.646297 0.000000) (3.274421 1.125071 0.000000) (5.663517 7.298393 0.000000) (9.137107 5.300407 0.000000) (9.645942 4.763386 0.000000) (9.645942 4.763386 0.000000) (10.154778 4.226366 0.000000) (10.014585 3.818321 0.000000) (9.747409 3.335636 0.000000) (9.747409 3.335636 0.000000) (9.480234 2.852952 0.000000) (8.691918 1.738304 0.000000) (9.562593 1.076016 0.000000) (11.042059 2.001239 0.000000) (11.292893 2.207107 0.000000))))
        (list (make-bspline-curve 3 '(0 0 0 0 1 1 1 1)
                                  '((12 1.5 0) (11.292893 2.207107 0) (12 10 0) (13 12 0)))
              (make-bspline-curve 3 '(0 0 0 0 1 1 1 1)
                                  '((10.62 0.12 0) (9 1 0) (10 10 0) (11 12 0))))
        (list (make-bspline-curve 3 '(0 0 0 0 1 1 1 1) '((13 12 0) (11 12 0) (4 12 0) (3 12 0)))
              (make-bspline-curve 3 '(0 0 0 0 1 1 1 1) '((12 10 0) (11 10 0) (4 10 0) (3 10 0))))
        (list (make-bspline-curve 3 '(0 0 0 0 1 1 1 1)
                                  '((3 12 0) (3 10 0) (1.447214 2.894427 0) (1 2 0)))
              (make-bspline-curve 3 '(0 0 0 0 1 1 1 1)
                                  '((4 12 0) (4 10 0) (4 3 0) (2.2 1.4 0))))))

#+nil
(with-open-file (s "/tmp/ribbons" :direction :output :if-exists :supersede)
  (iter (for (out in) in *ribbons*)
        (format s "~{~{~f~^ ~}~%~}" (sample-curve out 100))
        (format s "~{~{~f~^ ~}~%~}" (sample-curve in 100))))

;;; We need our own RIBBON-EVALUATE:
(defun ribbon-evaluate (patch i s d)
  "PATCH has n elements, each consisting of a tuple: (outer-bspline inner-bspline)."
  (let* ((base-point (bsc-evaluate (first (elt patch i)) (elt s i)))
	 (inner-point (bsc-evaluate (second (elt patch i)) (elt s i)))
	 (derivative (v* (v- inner-point base-point) 3.0d0)))
    (v+ base-point (v* derivative (gamma (elt d i)) *ribbon-multiplier*))))

(defun write-concave-patch (patch type filename &key (distance-type 'perpendicular) (output 'vtk))
  "Similar to WRITE-PATCH in patches.lisp, but PATCH is a list of bspline-tuples.
OUTPUT is one of (SPIDER RIBBONS PATCH)."
  (let* ((n (length patch))
         (points (points-from-angles (uniform-angles n))))
    (ecase output
      (spider
       (write-vtk-polylines
        (iter (for line in (spider-lines points))
              (collect (iter (for domain-point in line)
                             (collect (patch-evaluate patch points type distance-type
                                                      domain-point)))))
        filename))
      (ribbons
       (write-vtk-curves
        (iter (for (curve1 curve2) in patch)
              (for points1 = (sample-curve curve1 *resolution*))
              (for points2 = (sample-curve curve2 *resolution*))
              (appending (append (list points1 points2)
                                 (mapcar #'list points1 points2))))
        filename))
      (patch
       (write-ply-indexed-mesh
        (iter (for domain-point in (vertices points))
              (collect (patch-evaluate patch points type distance-type domain-point)))
        (triangles n) filename)))))

#+nil
(let ((*resolution* 30))
  (write-concave-patch *ribbons* 'ribbon "/tmp/proba.ply" :output 'patch)
  (write-concave-patch *ribbons* 'ribbon "/tmp/proba.vtk" :output 'spider)
  (write-concave-patch *ribbons* 'ribbon "/tmp/ribbons.vtk" :output 'ribbons))


;;; Try a new blend function.
;;; This blend function ensures that ribbons vanish for d > 0.

(defun blend (d i)
  (if (> (elt d i) 1)
      0
      (let* ((n (length d))
             (im (mod (1- i) n))
             (ip (mod (1+ i) n)))
        (/ (* (expt (elt d im) *exponent*)
              (expt (elt d ip) *exponent*)
              (expt (- 1 (elt d i)) *exponent*))
           (iter (for k from 0 below n)
                 (unless (> (elt d k) 1)
                   (for km = (mod (1- k) n))
                   (for kp = (mod (1+ k) n))
                   (sum (* (expt (elt d km) *exponent*)
                           (expt (elt d kp) *exponent*)
                           (expt (- 1 (elt d k)) *exponent*)))))))))

#+nil
(let ((*resolution* 50)
      (*exponent* 3))
  (write-patch (points-from-angles '(40 20 60 100 80))
	       'ribbon
	       "/tmp/spider.vtk"
	       :inner-points (points-from-angles '(40 20 60 100 80) 0.5d0)
	       :heights
	       '(((0.0d0 0.0d0 0.0d0 0.0d0)
		  (0.0d0 0.0d0 0.0d0 0.0d0)
		  (0.0d0 0.0d0 0.0d0 0.0d0)
		  (0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0)
		  (0.0d0 0.0d0 0.0d0 0.0d0))
		 ((0.0 0.0)
		  (0.0 0.0)
		  (0.0 0.0)
		  (0.0 0.0)
		  (0.0 0.0)))
	       :distance-type 'bilinear
               :spider nil))

#+nil
(let ((*resolution* 50)
      (*exponent* 3))
  (write-patch (points-from-angles '(40 20 60 100 80))
	       'ribbon
	       "/tmp/patch.vtk"
	       :inner-points (points-from-angles '(40 20 60 100 80) 0.5d0)
	       :heights
	       '(((0.0d0 0.1d0 0.1d0 0.0d0)
		  (0.0d0 0.2d0 0.3d0 0.4d0)
		  (0.4d0 0.6d0 0.6d0 0.4d0)
		  (0.4d0 0.5d0 0.6d0 0.4d0 0.2d0 0.0d0)
		  (0.0d0 0.2d0 0.1d0 0.0d0))
		 ((0.2 0.2)
		  (0.2 0.5)
		  (0.5 0.8)
		  (0.8 0.2)
		  (0.2 0.2)))
	       :distance-type 'bilinear
               :spider nil))


;;; Mean value + Kato-style s-parameters

(defun mean-kato-sd (v1 v2)
  (lambda (points p)
    (let ((d1 (efficient-mean-value points v1 p))
          (d2 (efficient-mean-value points v2 p)))
      (list (/ d1 (+ d1 d2))))))

#+nil
(let ((points '((-0.8 -0.8) (0.4 -0.8) (0.4 0) (0 0) (0 -0.4) (-0.4 -0.4) (-0.4 0) (-0.8 0)))
      (*wachspressp* nil)
      (size 800)
      (density 20))
  (bitmap-test points (mean-kato-sd '(0 0 1 1 1 1 1 1) '(1 1 1 1 1 .5 0 0))
               "/tmp/mean-kato-0.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(0 1 1 1 1 1 1 0) '(1 0 0 1 1 1 1 1))
               "/tmp/mean-kato-1.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(0 0 1 1 1 1 1 1) '(1 1 0 0 .5 1 1 1))
               "/tmp/mean-kato-2.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(1 0 0 1 1 1 1 1) '(1 1 1 0 0 0 0 1))
               "/tmp/mean-kato-3.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(1 1 0 0 .5 1 1 1) '(1 1 1 1 1 .5 0 0))
               "/tmp/mean-kato-4.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(1 1 1 0 0 0 0 1) '(0 1 1 1 1 1 1 0))
               "/tmp/mean-kato-5.pbm" :size size :density density))

#+nil
(let ((points '((-0.8 -0.8) (0.4 -0.8) (0.4 0) (0 0) (0 -0.4) (-0.4 -0.4) (-0.4 0) (-0.8 0)))
      (*wachspressp* nil)
      (size 800)
      (density 20))
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(0 1 1 1 1 1 1 0) p)))
               "/tmp/mean-kato-d-0.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(0 0 1 1 1 1 1 1) p)))
               "/tmp/mean-kato-d-1.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 0 0 1 1 1 1 1) p)))
               "/tmp/mean-kato-d-2.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 1 0 0 .5 1 1 1) p)))
               "/tmp/mean-kato-d-3.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 1 1 0 0 0 0 1) p)))
               "/tmp/mean-kato-d-4.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 1 1 1 1 .5 0 0) p)))
               "/tmp/mean-kato-d-5.pbm" :size size :density density))

#+nil
(let ((points '((-0.8 -0.8) (0.4 -0.8) (0.4 0) (0 0) (0 -0.4) (-0.2 -0.4) (-0.4 -0.4) (-0.4 0) (-0.8 0)))
      (*wachspressp* nil)
      (size 800)
      (density 20))
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 1 0 0 .5 1 1 1 1) p)))
               "/tmp/mean-kato-d-3-tomi.pbm" :size size :density density)
  (bitmap-test points (lambda (points p) (list (efficient-mean-value points '(1 1 1 1 1 1 .5 0 0) p)))
               "/tmp/mean-kato-d-5-tomi.pbm" :size size :density density)
  (bitmap-test points (mean-kato-sd '(1 1 0 0 .5 1 1 1 1) '(1 1 1 1 1 1 .5 0 0))
               "/tmp/mean-kato-s-4-tomi.pbm" :size size :density density))


;;; Yet another idea.
;;; S,D parameters are generated by Mean Value coordinates.
;;; D is simpler: base vertices have a 0 value, all other vertices are 1.
;;; For S, we first categorize all segments as valid or invalid.
;;; If the base is P0P1, then it has the normal vector N1=(P0y - P1y, P1x - P0x).
;;; Similarly, for segment i (Pi-1Pi), the normal vector is Ni=(Pi-1y - Piy, Pix - Pi-1x)
;;; If Ni * N1 < 0, segment i is valid; otherwise invalid.
;;; Let L be the sum of all valid segments' lengths.
;;; Then the values are distributed as follows.
;;; P1 gets 1. Now 1 is the current value.
;;; From here, moving along an invalid segment does not change the current value.
;;; So if P1P2 is invalid, P2 will also get 0.
;;; Suppose that P2P3 is valid. Then P3 will get 1 - (length of P2P3) / L,
;;; and this is the new current value.
;;; When we get to P0, it will get the value 0.

(defun mean-s-values (points index)
  "INDEX is the index of the base, i.e. I for the segments with vertices I and I-1."
  (let* ((n (length points))
         (result (make-array n)))
    (labels ((p (i) (elt points (mod i n)))
             (normal (i)
               (list (- (second (p (1- i))) (second (p i)))
                     (- (first (p i)) (first (p (1- i))))))
             (validp (i) (< (scalar-product (normal index) (normal i)) 0))
             (segment-length (i) (vlength (v- (p i) (p (1- i))))))
      (let ((total-length (iter (for i from 0 below n)
                                (when (validp i)
                                  (sum (segment-length i))))))
        (iter (with current = 1)
              (for j from 1 below (1- n))
              (for i = (mod (+ index j) n))
              (when (validp i)
                (setf current (max 0 (- current (/ (segment-length i) total-length)))))
              (setf (elt result i) current))))
    (setf (elt result index) 1)
    (setf (elt result (mod (1- index) n)) 0)
    (coerce result 'list)))

(defun mean-d-values (points index)
  (iter (for i from 0 below (length points))
        (if (or (= i index) (= i (1- index)))
            (collect 0)
            (collect 1))))


;;; New test

(defparameter *concave*
  '((-10 8) (-4 3) (-10 -6) (9 -10) (8 0) (-2 -3) (6 9)))
(defparameter *convex*
  '((-10 8) (-10 3) (-10 -6) (9 -10) (8 0) (7 4) (6 9)))

(defun peti-test-concave (points p)
  (declare (ignore points))
  (let ((u (efficient-mean-value *concave* (mapcar #'first *convex*) p))
        (v (efficient-mean-value *concave* (mapcar #'second *convex*) p)))
    (list (elt (compute-parameter 'mean-bilinear 's *convex* (list u v)) 0)
          (elt (compute-parameter 'mean-bilinear 'd *convex* (list u v)) 0))))

(defun peti-test-convex (points p)
  (declare (ignore points))
  (list (elt (compute-parameter 'mean-bilinear 's *convex* p) 0)
        (elt (compute-parameter 'mean-bilinear 'd *convex* p) 0)))

#+nil
(setf *wachspressp* nil)
#+nil
(bitmap-test nil #'peti-test-concave "/tmp/peti-concave.pgm" :object-size 25)
#+nil
(bitmap-test nil #'peti-test-convex "/tmp/peti-convex.pgm" :object-size 25)


;;; Returning to the original concept: composite ribbons
;;; Implementation notes:
;;; - uses Bezier curves, so no cusp at the concave point
;;;   (so later we should use B-splines)
;;; - points: all vertices
;;; - convex-points: a list of convex polygons covering the whole concave domain
;;; - lines: all edges, the composite edge counting as one
;;; - coords: list of (OUTER INNER)
;;;   - outer: all curves as a list of Bezier control points
;;;   - inner: all inner ribbon curves as a list of Bezier control points,
;;;            but without the first and last control point
;;;            (these are extracted from OUTER)

(defun bezier-to-bspline (curve)
  (let ((n (length curve)))
    (make-bspline-curve (1- n)
                        (append (make-list n :initial-element 0) (make-list n :initial-element 1))
                        curve)))

(defparameter *points*
  '((-1 -1) (1 -1) (1 1) (0 1) (0 0) (-1 0)))

(defparameter *convex-points*
  '(((-1 -1) (1 -1) (1 0) (-1 0))
    ((1 0) (1 1) (0 1) (0 0))))

(defparameter *lines*
  '(((-1 0) (-1 -1)) ((-1 -1) (1 -1)) ((1 -1) (1 1)) ((1 1) (0 1)) ((0 1) (0 0) (-1 0))))

(defparameter *coords*
  '((((-1 0 0) (-1 -0.25 0) (-1 -0.75 0) (-1 -1 0))
     ((-1 -1 0) (-0.75 -1 0) (0.75 -1 0) (1 -1 0))
     ((1 -1 0) (1 -0.75 0) (1 0.75 0) (1 1 0))
     ((1 1 0) (0.75 1 0) (0.25 1 0) (0 1 0))
     ((0 1 0) (0 0.75 0) (0 0.25 0) (0 0 0) (0 0 0) (-0.25 0 0) (-0.75 0 0) (-1 0 0)))
    (((-0.75 -0.25 0.2) (-0.75 -0.75 0.2))
     ((-0.75 -0.75 0.2) (0.75 -0.75 0.2))
     ((0.75 -0.75 0.2) (0.75 0.75 0.2))
     ((0.75 0.75 0.2) (0.25 0.75 0.2))
     ((0.25 0.75 0.2) (0.25 0 0) (0.25 -0.25 0) (0 -0.25 0) (-0.75 -0.25 0.2)))))

(defun lines-next (lst x)
  (let ((i (position x lst :test 'equal))
        (n (length lst)))
    (elt lst (mod (1+ i) n))))

(defun lines-prev (lst x)
  (let ((i (position x lst :test 'equal))
        (n (length lst)))
    (elt lst (mod (1- i) n))))

(defun compute-concave-distance (points lines line p dir)
  (if (eq dir 's)
      (let ((dnext (compute-concave-distance points lines (lines-next lines line) p 'd))
            (dprev (compute-concave-distance points lines (lines-prev lines line) p 'd)))
        (let ((dd (+ dprev dnext)))
          (if (< dd *epsilon*)
              0.5
              (/ dprev dd))))
      (let ((values (mapcar (lambda (x) (if (member x line :test 'equal) 0 1)) points)))
        (mean-value points values p))))

(defun compute-concave-parameter (dir points lines p &optional no-tiny-p)
  (macrolet ((tiny-lambda ((args) &body body)
	       `(lambda (,args)
		  (if no-tiny-p
		      (let ((result (progn ,@body)))
			(if (< result *tiny*) 0.0d0 result))
		      (progn ,@body)))))
    (mapcar (tiny-lambda (lst) (compute-concave-distance points lines lst p dir)) lines)))

(defun concave-sp-patch-evaluate (patch points lines domain-point)
  (let* ((n (length lines))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
	 (d (compute-concave-parameter 'd points lines p t))
	 (s (compute-concave-parameter 's points lines p t))
         (*use-gamma* nil))
    (iter (for i from 0 below n)
	  (with result = '(0 0 0))
	  (setf result
		(v+ result
                    (v* (ribbon-evaluate patch i s d)
                        (ribbon-blend d i))))
	  (finally (return result)))))

(defun write-concave-sp-patch (points convex-points lines coords filename &key spider)
  (let* ((patch (generate-patch (first coords) (second coords))))
    (if spider
	(write-vtk-polylines
	 (iter (for line in (spider-lines points))
	       (collect (iter (for domain-point in line)
			      (collect (concave-sp-patch-evaluate patch points
                                                                  lines domain-point)))))
	 filename)
	(let* ((vertices (mapcar #'vertices convex-points))
               (triangles (iter (for i from 0 below (length convex-points))
                                (for j first 0 then (+ j (length (elt vertices (1- i)))))
                                (for tri = (triangles (length (elt convex-points i))))
                                (appending (mapcar (lambda (p) (mapcar (lambda (k) (+ k j)) p))
                                                   tri)))))
          (write-vtk-indexed-mesh
           (iter (for domain-point in (reduce #'append vertices))
                 (collect (concave-sp-patch-evaluate patch points lines domain-point)))
           triangles filename)))))

#+nil
(write-constraint-ribbons
 *points*
 "/tmp/ribbons4.vtk"
 :coords *coords*
 :resolution 20)

#+nil
(let ((*resolution* 50)
      (*ribbon-multiplier* 1.0d0)
      (*wachspressp* nil))
  (write-concave-sp-patch *points* *convex-points* *lines* *coords*
                          "/tmp/sp-patch4.vtk" :spider nil))

#+nil
(flet ((s-fun (i)
         (lambda (points p)
           (list (elt (compute-concave-parameter 's points *lines* p) i))))
       (d-fun (i)
         (lambda (points p)
           (list (elt (compute-concave-parameter 'd points *lines* p) i)))))
  (dotimes (i 5)
    (bitmap-test *points* (s-fun i)
                 (format nil "/tmp/~as.pgm" i)
                 :object-size 2.0d0)
    (bitmap-test *points* (d-fun i)
                 (format nil "/tmp/~ad.pgm" i)
                 :object-size 2.0d0)))
