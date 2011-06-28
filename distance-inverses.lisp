(in-package :cl-nurbs-tests)

(defun radial-inverse (segments points sd)
  (let* ((p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments))
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (- (second base-x)) (first base-x)))
	 (a (in-system base-x base-y (v- p0 p1)))
	 (b (in-system base-x base-y (v- p3 p2)))
	 (length (point-distance p1 p2))
	 (max (iter (for q1 in points)
		    (maximize (iter (for q2 in points)
				    (maximize (point-distance q1 q2))))))
	 (axby-aybx (- (* (first a) (second b)) (* (second a) (first b)))))
    (if (< (abs axby-aybx) *epsilon*)
	(v+ p1
	    (v* base-x length (first sd))
	    (v* base-y (* (second sd) max)))
	(let* ((q1 (v+ (list length 0) (v* b (/ (* length (second a)) axby-aybx))))
	       (q2 (list (* length (first sd)) 0))
	       (v (v* (vnormalize (v- q2 q1)) (second sd) max)))
	  (v+ p1
	      (v* base-x (first (v+ q2 v)))
	      (v* base-y (second (v+ q2 v))))))))

(defun line-sweep-inverse (segments points sd)
  (when *centralized-line-sweep*
    (error "The invese function is only implemented for standard line sweeps."))
  (let* ((p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments))
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (- (second base-x)) (first base-x)))
	 (a (in-system base-x base-y (v- p0 p1)))
	 (b (in-system base-x base-y (v- p3 p2)))
	 (length (point-distance p1 p2))
	 (max (iter (for q1 in points)
		    (maximize (iter (for q2 in points)
				    (maximize (point-distance q1 q2))))))
	 (q (list (* length (first sd)) 0))
	 (center (central-point points (lines-from-points points) t))
	 (comb (v+ (v* a (- 1 (first sd)) (- 1 (first sd)))
		   (v* (center-line a b length
				    (in-system base-x base-y (v- center p1)))
		       2 (- 1 (first sd)) (first sd))
		   (v* b (first sd) (first sd))))
	 (v (v* comb (second sd) max (/ (second comb)))))
    (v+ p1
	(v* base-x (first (v+ q v)))
	(v* base-y (second (v+ q v))))))

(defun write-d-lines (points type filename &key inner-points heights coords
		      (distance-type 'radial) (reverse-fn #'radial-inverse) (d 1/3))
  (let* ((n (length points))
	 (patch (or (and coords (generate-patch (first coords) (second coords)))
		    (generate-patch-from-heights points inner-points heights))))
    (write-vtk-polylines
     (iter (for i from 0 below n)
	   (collect
	    (iter (for j from 0 to 100)
		  (for s = (/ j 100))
		  (for domain-point =
		       (funcall reverse-fn
				(iter (for k from 0 below 4)
				      (collect (elt points (mod (+ i k) n))))
				points (list s d)))
		  (collect (patch-evaluate patch points type distance-type
					   domain-point)))))
     filename)))

(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 1.0))
  (write-d-lines (domain-from-curves (first *coords*) 'regular) 'ribbon
		 "/tmp/patch.vtk"
		 :coords *coords* :distance-type 'radial :reverse-fn #'radial-inverse))

(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 1.0))
  (write-d-lines (domain-from-curves (first *coords*) 'circular) 'ribbon
		 "/tmp/patch.vtk"
		 :coords *coords* :distance-type 'line-sweep :reverse-fn #'line-sweep-inverse))
