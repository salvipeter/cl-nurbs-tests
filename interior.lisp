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

(defun write-interior-patch (points interior-patch-fn filename
			     &key inner-points heights coords spider)
  (let* ((n (length points))
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

(defun write-interior-surface (interior-patch-fn filename &optional (resolution 100))
  (let ((points (make-array (list resolution resolution))))
    (iter (for i from 0 below resolution)
	  (for u = (1- (* (/ i (1- resolution)) 2)))
	  (iter (for j from 0 below resolution)
		(for v = (1- (* (/ j (1- resolution)) 2)))
		(setf (aref points i j)
		      (funcall interior-patch-fn (list u v)))))
    (write-points2-vtk points filename)))

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

(flet ((interior (uv)
	 (append uv (list (- (cos (+ (first uv) (second uv))) 0.5d0)))))
  (let ((*alpha* 100.0d0))
    (write-interior-patch *points* #'interior 
			  "/tmp/proba.vtk" :coords *coords*)
    (write-interior-surface #'interior "/tmp/proba2.vtk")))
