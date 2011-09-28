(in-package :cl-nurbs-tests)

(defun os-vertices ()
  (iter (for i from 0 below *resolution*)
	(for u = (- (* (/ i (1- *resolution*)) 2) 1))
	(appending
	 (iter (for j from 0 below *resolution*)
	       (for v = (- (* (/ j (1- *resolution*)) 2) 1))
	       (for alpha = (acos (first (vnormalize (list u v)))))
	       (when (< alpha 0) (setf alpha (- alpha)))
	       (when (> alpha (/ pi 2)) (setf alpha (- pi alpha)))
	       (when (> alpha (/ pi 4)) (setf alpha (- (/ pi 2) alpha)))
	       (collect (v* (list u v) (cos alpha)))))))

(defun os-triangles ()
  (flet ((vertex (i j) (+ (* j *resolution*) i)))
    (iter (for i from 0 below (1- *resolution*))
	  (for left = (< i (floor *resolution* 2)))
	  (appending
	   (iter (for j from 0 below (1- *resolution*))
		 (for bottom = (< j (floor *resolution* 2)))
		 (if (or (and left bottom) (not (or left bottom)))
		     (progn
		       (collect (list (vertex i j)
				      (vertex i (1+ j))
				      (vertex (1+ i) (1+ j))))
		       (collect (list (vertex i j)
				      (vertex (1+ i) (1+ j))
				      (vertex (1+ i) j))))
		     (progn
		       (collect (list (vertex i j)
				      (vertex i (1+ j))
				      (vertex (1+ i) j)))
		       (collect (list (vertex i (1+ j))
				      (vertex (1+ i) (1+ j))
				      (vertex (1+ i) j))))))))))

(defun os-compute-parameters (p)
  "Returns ((Sribbon) (Dribbon Dpoint))."
  (let* ((v (vnormalize p))
	 (alpha (acos (first v))))
    (list (list (/ (if (< (second v) 0)
		       (- (* 2 pi) alpha)
		       alpha)
		   (* 2 pi)))
	  (list (- 1.0d0 (vlength p))
		(vlength p)))))

(defun os-aux-evaluate (patch s d)
  (let ((curve-point (bezier (first (first patch)) (first s))))
    (destructuring-bind (p n scale) (second (first patch))
      (let* ((n (vnormalize n))
	     (diff (vnormalize (v- curve-point p)))
	     (dir (v- diff (v+ (v* n (scalar-product diff n))))))
	(v+ p (v* (vnormalize dir) scale (second d)))))))

(defun os-spider-lines ()
  (append
   (iter (for i from 1 to *spider-lines*)
	 (for d = (/ i *spider-lines*))
	 (collect
	  (iter (for j from 0 below *resolution*)
		(for alpha = (* (/ j (1- *resolution*)) 2.0d0 pi))
		(collect (v* (list (cos alpha) (sin alpha)) d)))))
   (iter (for j from 0 below *resolution*)
	 (for alpha = (* (/ j (1- *resolution*)) 2.0d0 pi))
	 (when (zerop (mod j *spider-density*))
	   (collect
	    (iter (for i from 0 to *spider-lines*)
		  (for d = (/ i *spider-lines*))
		  (collect (v* (list (cos alpha) (sin alpha)) d))))))))

(defun os-patch-evaluate (patch domain-point)
  (let ((p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0))
		    domain-point)))
    (destructuring-bind (s d) (os-compute-parameters p)
      (v+ (v* (ribbon-evaluate patch 0 s d) (ribbon-blend d 0))
	  (v* (os-aux-evaluate patch s d) (ribbon-blend d 1))))))

(defun os-write-patch (patch filename &key spider)
  "Patch:((RibbonCurve (AuxiliaryPoint APNormal APRadius)) (RibbonInnerCurve))."
  (if spider
      (write-vtk-polylines
       (iter (for line in (os-spider-lines))
	     (collect (iter (for domain-point in line)
			    (collect (os-patch-evaluate patch domain-point)))))
       filename)
      (write-vtk-indexed-mesh
       (iter (for domain-point in (os-vertices))
	     (collect (os-patch-evaluate patch domain-point)))
       (os-triangles) filename)))

(defun os-write-constraint-grid (patch filename)
  (write-vtk-curves (cons (first (first patch)) (second patch)) filename))

(defun os-write-constraint-ribbons (patch filename)
  "TODO: auxiliary ribbon is always displayed with a Z normal."
  (let* ((aux-point (first (second (first patch))))
	 ;; (aux-normal (vnormalize (second (second (first patch)))))
	 (aux-radius (third (second (first patch))))
	 (curves (iter (for curve1 in (first patch))
		       (for curve2 in (second patch))
		       (for points1 =
			    (iter (for u from 0 to 1 by (/ *resolution*))
				  (collect (bezier curve1 u))))
		       (for points2 =
			    (iter (for u from 0 to 1 by (/ *resolution*))
				  (collect (bezier curve2 u))))
		       (appending (append (list points1 points2)
					  (mapcar #'list points1 points2)))))
	 (aux-curves (iter (for i from 0 below *resolution*)
			   (for alpha = (* 2 pi (/ i (1- *resolution*))))
			   (for v = (v* (list (cos alpha) (sin alpha) 0)
					aux-radius))
			   (for p = (v+ aux-point v)) ;todo
			   (collect (list aux-point p) into edges)
			   (collect p into circle)
			   (finally (return (cons circle edges))))))
    (write-vtk-curves (append aux-curves curves) filename)))

(let ((patch '((((0 0 0) (3 0 0) (5 -2 0) (7 1 0) (7 4 0) (5 6 0) (3 4 0)
		 (0 4 0) (-3 4 0) (-6 6 0) (-9 4 0) (-9 0 0) (-6 -2 0) (-3 0 0)
		 (0 0 0))
		((0 2.2 1.5) (0 0 1.5) 1))
	       (((0 1 1) (2 1 1) (5 -1 1) (6 1 1) (6 4 1) (5 5 1) (2 3 1)
		 (0 3 1) (-2 3 1) (-6 5 1) (-8 3 1) (-8 1 1) (-6 -1 1) (-2 1 1)
		 (0 1 1)))))
      (*resolution* 60)
      (*ribbon-multiplier* 0.5)
      (*spider-lines* 6))
  (os-write-constraint-grid patch "/tmp/grid.vtk")
  (os-write-constraint-ribbons patch "/tmp/ribbon.vtk")
  (os-write-patch patch "/tmp/proba.vtk")
  (os-write-patch patch "/tmp/proba2.vtk" :spider t))
