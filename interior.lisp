(in-package :cl-nurbs-tests)

;;; Patches with interior

(defun interior-patch-evaluate (patch interior-patch-fn points domain-point)
  (let* ((n (length points))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
	 (d (compute-parameter 'line-sweep 'd points p t))
	 (s (compute-parameter 'line-sweep 's points p t))
	 (b (compute-parameter 'perpendicular 'd points p t)))
    (iter (for i from 0 below n)
	  (with result = '(0 0 0))
	  (setf result
		(v+ result
		    (v* (ribbon-evaluate patch i s d)
			(interior-ribbon-blend b i))))
	  (finally (return (v+ result
			       (v* (funcall interior-patch-fn domain-point)
				   (interior-ribbon-blend b n))))))))

(defun write-interior-patch (points interior-patch-fn alpha filename
			     &key inner-points heights coords spider)
  (let* ((n (length points))
	 (*alpha* (compute-alpha points alpha 'perpendicular))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (if spider
	(write-vtk-polylines
	 (iter (for line in (spider-lines points))
	       (collect (iter (for domain-point in line)
			      (collect (interior-patch-evaluate patch interior-patch-fn
								points domain-point)))))
	 filename)
	(write-vtk-indexed-mesh
	 (iter (for domain-point in (vertices points))
	       (collect (interior-patch-evaluate patch interior-patch-fn points domain-point)))
	 (triangles n) filename))))

(defun write-interior-surface (points interior-patch-fn filename)
  (write-vtk-indexed-mesh
   (iter (for domain-point in (vertices points))
	 (collect (funcall interior-patch-fn domain-point)))
   (triangles (length points)) filename))

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
(defparameter *points*
  (domain-from-curves (first *coords*) 'circular-mod))

#+nil
(flet ((interior (uv)
	 (append uv (list (if (< (vlength uv) 0.3)
			      1.0d0
			      (- 1.9d0 (* 3 (vlength uv))))))))
  (let ((*resolution* 40)
	(interior (rotational '((0 0 0.5) (0.1 0 0.5) (0.3 0 0) (1 0 0)))))
    (write-interior-patch *points* #'interior 0.5d0 "/tmp/proba.vtk"
			  :coords *coords*)
    (write-interior-surface *points* #'interior "/tmp/proba2.vtk")))

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

; (defparameter *points* (points-from-angles '(40 20 60 100 80)))
#+nil
(let ((multiplier 2.0d0))
  (write-color-aux-blend-test *points* '(0.0d0 0.0d0) multiplier
			      "/tmp/aux-blend-test-black.ppm" 400 :trim '(0.89d0 0.91d0))
  (write-color-aux-blend-test *points* '(0.0d0 0.0d0) multiplier
			      "/tmp/aux-blend-test-white.ppm" 400 :trim '(0.89d0 0.91d0))
  (write-color-aux-blend-test *points* '(0.0d0 0.0d0) multiplier
			      "/tmp/aux-blend-test-clean.ppm" 400 :trim '(0.89d0 0.91d0)))
