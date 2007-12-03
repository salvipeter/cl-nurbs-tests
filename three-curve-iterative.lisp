(in-package :cl-nurbs)

(defparameter *alpha* 0.5)
(defparameter *stability* 1000)

(defun project-to-vector (p from n)
  (scalar-product (v- p from) n))

(defun continuity-measure (curve parameters curvatures)
  (let ((squared-differences
	 (map 'list #'(lambda (x y)
			(let ((difference (- (bsc-curvature curve x) y)))
			  (* *stability* difference difference)))
	      parameters curvatures)))
    (+ (* *alpha*
	  (/ (apply #'+ (cons (first squared-differences)
			      (last squared-differences)))
	     2))
       (* (- 1 *alpha*)
	  (/ (apply #'+ (rest (butlast squared-differences)))
	     (- (length parameters) 2))))))

(defun 3it-fair (curves &key (resolution 100)
		 (target-iteration 100) (simplex-iteration 5)
		 (fairing-iteration 5))
  (when (< (length (control-points (second curves))) 4)
    (error "Too few control points in the middle curve for fairing"))
  (let* ((left (first curves))
	 (curve (second curves))
	 (right (third curves))
	 (cpts-left (control-points left))
	 (left-n (length cpts-left))
	 (cpts (control-points curve))
	 (n (length cpts))
	 (cpts-right (control-points right))
	 (start-point (elt cpts-left (- left-n 1)))
	 (end-point (elt cpts-right 0))
	 (start-tangent (vnormalize
			 (v- start-point (elt cpts-left (- left-n 2)))))
	 (end-tangent (vnormalize (v- end-point (elt cpts-right 1))))
	 (start-u (project-to-vector (elt cpts 1) start-point start-tangent))
	 (end-u (project-to-vector (elt cpts (- n 2)) end-point end-tangent))
	 (start-curvature (bsc-curvature left (bsc-upper-parameter left)))
	 (end-curvature (bsc-curvature right (bsc-lower-parameter right)))
	 (low (bsc-lower-parameter curve))
	 (high (bsc-upper-parameter curve))
	 (parameters (arc-length-sampling curve low high resolution))
	 (curvatures (target-curvature curve parameters target-iteration
				       :start-value start-curvature
				       :end-value end-curvature))
	 (new-curve (copy-bspline-curve curve))
	 (points (control-points new-curve)))
    (dotimes (i fairing-iteration)
      (let ((new-u (car (downhill-simplex:minimize
			  (lambda (x)
			    (let ((old (elt points 1)))
			      (setf (elt points 1)
				    (v+ start-point (v* start-tangent (car x))))
			      (prog1
				  (continuity-measure new-curve
						      parameters curvatures)
				(setf (elt points 1) old))))
			  (list start-u) simplex-iteration))))
	(setf start-u new-u
	      (elt points 1) (v+ start-point (v* start-tangent new-u))))
      (iter (for j from 2 below (- n 2))
	    (flet ((fairness-fn (point)
		     (let ((old (elt points j)))
		       (setf (elt points j) point)
		       (prog1
			   (continuity-measure new-curve
					       parameters curvatures)
			 (setf (elt points j) old)))))
	      (setf (elt points j)
		    (downhill-simplex:minimize
		     #'fairness-fn (elt points j) simplex-iteration))))
      (let ((new-u (car (downhill-simplex:minimize
			  (lambda (x)
			    (let ((old (elt points (- n 2))))
			      (setf (elt points (- n 2))
				    (v+ end-point (v* end-tangent (car x))))
			      (prog1
				  (continuity-measure new-curve
						      parameters curvatures)
				(setf (elt points (- n 2)) old))))
			  (list end-u) simplex-iteration))))
	(setf end-u new-u
	      (elt points (- n 2)) (v+ end-point (v* end-tangent new-u)))))
    new-curve))

(defun split-to-three (curve u1 u2)
  "Split CURVE to three subcurves at the parameter positions U1 and U2."
  (let* ((left-rest (bsc-split-curve curve u1))
	 (middle-right (bsc-split-curve (second left-rest) u2)))
    (cons (first left-rest) middle-right)))

(defun three-curve-original (curve u1 u2 extrusion filename)
  (let ((split (split-to-three curve u1 u2)))
    (write-rdn (mapcar (lambda (curve) (bsc-extrude curve extrusion))
		       split)
	       filename)
    split))

(defun three-curve-iterative-test (curve u1 u2 extrusion filename &rest keys)
  (let* ((split (split-to-three curve u1 u2))
	 (faired (apply #'3it-fair split keys))
	 (new-split (list (first split) faired (third split))))
    (write-rdn (mapcar (lambda (curve) (bsc-extrude curve extrusion))
		       new-split)
	       filename)
    new-split))
