(in-package :cl-nurbs-tests)

;;; Patches with interior

(defun interior-patch-evaluate (patch interior-patch-fn points alpha domain-point with-auxiliary-p)
  (let* ((n (length points))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
	 (d (compute-parameter 'line-sweep 'd points p t))
	 (s (compute-parameter 'line-sweep 's points p t))
	 (b (compute-parameter 'perpendicular 'd points p t)))
    (flet ((interior (patch-point)
	     (let ((closest (downhill-simplex:minimize
			     (lambda (uv)
			       (point-distance (subseq (funcall interior-patch-fn uv) 0 2)
					       (subseq patch-point 0 2)))
			     '(0 0) 100)))
	       (funcall interior-patch-fn closest))))
      (iter (for i from 0 below n)
	    (with result = '(0 0 0))
	    (setf result
		  (v+ result
		      (v* (ribbon-evaluate patch i s d)
			  (interior-ribbon-blend b i))))
	    (finally (return (let ((patch-point (patch-evaluate patch points 'ribbon
								'perpendicular domain-point)))
			       (v+ result
				   (v* (if with-auxiliary-p
					   (v* (v- (interior patch-point)
						   (v* patch-point alpha))
					       (/ (- 1.0d0 alpha)))
					   (interior patch-point))
				       (interior-ribbon-blend b n))))))))))

(defun write-interior-patch (points interior-patch-fn alpha filename
			     &key with-auxiliary-p inner-points heights coords spider)
  (let* ((n (length points))
	 (*alpha* (compute-alpha points alpha 'perpendicular))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (if spider
	(write-vtk-polylines
	 (iter (for line in (spider-lines points))
	       (collect (iter (for domain-point in line)
			      (collect (interior-patch-evaluate patch interior-patch-fn
								points alpha domain-point
								with-auxiliary-p)))))
	 filename)
	(write-vtk-indexed-mesh
	 (iter (for domain-point in (vertices points))
	       (collect (interior-patch-evaluate patch interior-patch-fn points alpha
						 domain-point with-auxiliary-p)))
	 (triangles n) filename))))

(defun scale-points (points scale)
  (mapcar (lambda (p) (v* p scale)) points))

(defun write-interior-surface-on-domain (points interior-patch-fn filename)
  (write-vtk-indexed-mesh
   (iter (for domain-point in (vertices (scale-points points 1)))
	 (collect (funcall interior-patch-fn domain-point)))
   (triangles (length points)) filename))

(defun write-interior-surface-on-circle (interior-patch-fn filename &optional (resolution 100))
  (let ((points (make-array (list resolution resolution))))
    (iter (for i from 0 below resolution)
	  (for u = (1- (* (/ i (1- resolution)) 2)))
	  (iter (for j from 0 below resolution)
		(for v = (1- (* (/ j (1- resolution)) 2)))
		(for point = (list u v))
		(when (> (vlength point) 1.0d0)
		  (setf point (vnormalize point)))
		(setf (aref points i j)
		      (funcall interior-patch-fn point))))
    (write-points2-vtk points filename)))

(defun rotational (curve)
  "CURVE is a Bezier-curve in the XZ plane, rotated around the Z axis."
  (flet ((rotate (p a)
	   (destructuring-bind (x y z) p
	     (list (- (* (cos a) x) (* (sin a) y))
		   (+ (* (sin a) x) (* (cos a) y))
		   z))))
    (lambda (xy)
      (let ((u (vlength xy))
	    (phi (* (acos (first (vnormalize xy)))
		    (if (< (second xy) 0) -1 1))))
	(rotate (bezier curve u) phi)))))

;;; Example: [uses n-sided-domain-test.lisp]
(defparameter *coords*
  (generate-planar-patch-by-points
   '((1.0d0 6.432571489356498d-16)
     (0.1304590474444296d0 0.9914536988381717d0)
     (-0.6664991597115454d0 0.7455057814019982d0)
     (-0.9998187154970686d0 -0.019040381870958342d0)
     (0.4087121522622124d0 -0.9126633424177779d0))))
(defparameter *coords*
  '((((-2 -9 0) (0 -9 1) (0 -9 1) (2 -9 0))
     ((2 -9 0) (3 -7 2) (7 1 2) (8 3 0))
     ((8 3 0) (7 5 1) (7 5 1) (6 7 0))
     ((6 7 0) (4 7 2) (-4 7 2) (-6 7 0))
     ((-6 7 0) (-7 5 1) (-7 5 1) (-8 3 0))
     ((-8 3 0) (-7 1 2) (-3 -7 2) (-2 -9 0)))
    (((-1 -7 3) (1 -7 3))
     ((1 -7 3) (6 3 3))
     ((6 3 3) (5 5 3))
     ((5 5 3) (-5 5 3))
     ((-5 5 3) (-6 3 3))
     ((-6 3 3) (-1 -7 3)))))
(defparameter *points*
  (domain-from-curves (first *coords*) 'circular-mod))

#+nil
(let ((*resolution* 40)
      (alpha 0.5d0)
      (*exponent* 2)
      (*ribbon-multiplier* 1.0d0)
      (interior (rotational '((0 0 7) (1 0 7) (2 0 7) (3 0 7)
			      (3 0 -1) (6 0 2) (8 0 1) (10 0 0)))))
  (write-interior-patch *points* interior alpha "/tmp/interior-weak.vtk"
			:coords *coords*)
  (write-interior-patch *points* interior alpha "/tmp/interior-strong.vtk"
			:coords *coords* :with-auxiliary-p t)
  (write-interior-surface-on-circle interior "/tmp/interior.vtk")
  (write-constraint-grid *points* "/tmp/grid.vtk" :coords *coords*)
  (write-constraint-ribbons *points* "/tmp/ribbons.vtk" :coords *coords* :resolution 20)
  (write-patch *points* 'ribbon "/tmp/patch.vtk" :coords *coords*))

(defun write-color-interior-blend-test (points alpha filename r &key (trim '(0.89d0 0.91d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
	   (lines (lines-from-points points))
	   (*alpha* (compute-alpha points alpha 'perpendicular))
	   (colors (cons '(255 255 255) (generate-colors n))))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (compute-parameter 'perpendicular 'd points p t)) 
			 (blends (cons (interior-ribbon-blend d n)
				       (iter (for i from 0 below n)
					     (collect (interior-ribbon-blend d i))))))
                    (if (and trim (some (lambda (x) (< (first trim) x (second trim))) blends))
			(format s "0 0 0~%")
			(format s "~{~d~^ ~}~%"
				(mapcar #'round
					(reduce (lambda (x y) (mapcar #'+ x y))
						(mapcar #'v* colors blends))))))
                  (format s "255 255 255~%")))))))))

; (defparameter *points* (points-from-angles '(40 20 60 100 80)))
#+nil
(write-color-interior-blend-test *points* 0.2d0 "/tmp/interior-blend-test.ppm" 400
				 :trim '(0.89d0 0.91d0))

(defun write-color-aux-blend-test (points point multiplier filename r &key
				   (trim '(0.89d0 0.91d0)))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
	   (lines (lines-from-points points))
	   (colors (cons '(255 255 255) (generate-colors n))))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (cons (* (point-distance p point) multiplier)
				  (compute-parameter 'perpendicular 'd points p t))) 
			 (blends (iter (for i from 0 to n)
				       (collect (ribbon-blend d i)))))
                    (if (and trim (some (lambda (x) (< (first trim) x (second trim))) blends))
			(format s "0 0 0~%")
			(format s "~{~d~^ ~}~%"
				(mapcar #'round
					(reduce (lambda (x y) (mapcar #'+ x y))
						(mapcar #'v* colors blends))))))
                  (format s "255 255 255~%")))))))))

(defun write-aux-blend-voronoi (points point multiplier filename r &key
				(threshold 0.01d0))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (n (length points))
	   (lines (lines-from-points points)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P3~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((d (cons (* (point-distance p point) multiplier)
				  (compute-parameter 'perpendicular 'd points p t))) 
			 (blends (iter (for i from 0 to n)
				       (collect (ribbon-blend d i)))))
                    (if (destructuring-bind (a b)
			    (subseq (sort blends #'>) 0 2)
			  (< (- a b) threshold))
			(format s "0 0 0~%")
			(format s "127 127 127~%")))
                  (format s "255 255 255~%")))))))))

; (defparameter *points* (points-from-angles '(40 20 60 100 80)))
#+nil
(let ((multiplier 0.5d0))
  (write-color-aux-blend-test *points* '(0.0d0 0.0d0) multiplier
			      "/tmp/aux-blend-test-black.ppm" 400 :trim '(0.89d0 0.91d0))
  (write-aux-blend-voronoi *points* '(0.0d0 0.0d0) multiplier
			   "/tmp/aux-blend-test-white.ppm" 400)
  (write-color-aux-blend-test *points* '(0.0d0 0.0d0) multiplier
			      "/tmp/aux-blend-test-clean.ppm" 400 :trim nil))
