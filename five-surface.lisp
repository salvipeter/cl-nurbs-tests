;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun xnode-resolution (xnode resolution)
  (let* ((res (if (listp resolution) resolution (list resolution resolution)))
	 (low (bss-lower-parameter (first xnode)))
	 (high (bss-upper-parameter (first xnode)))
	 (uv (mapcar (lambda (a b) (/ (+ a b) 2)) low high))
	 (cu1 (bss-get-surface-curve (second xnode) (second uv) :u-curve t))
	 (cu2 (bss-get-surface-curve (first xnode) (second uv) :u-curve t))
	 (cu3 (bss-get-surface-curve (third xnode) (second uv) :u-curve t))
	 (cv1 (bss-get-surface-curve (fourth xnode) (first uv) :u-curve nil))
	 (cv2 (bss-get-surface-curve (first xnode) (first uv) :u-curve nil))
	 (cv3 (bss-get-surface-curve (fifth xnode) (first uv) :u-curve nil))
	 (ru1 (bsc-estimate-arc-length cu1))
	 (ru2 (bsc-estimate-arc-length cu2))
	 (ru3 (bsc-estimate-arc-length cu3))
	 (rv1 (bsc-estimate-arc-length cv1))
	 (rv2 (bsc-estimate-arc-length cv2))
	 (rv3 (bsc-estimate-arc-length cv3))
	 (total (list (+ ru1 ru2 ru3) (+ rv1 rv2 rv3))))
    (mapcar (lambda (ruv)
	      (mapcar (lambda (x total-resolution total)
			(round (* total-resolution x) total))
		      ruv res total))
	    `((,ru2 ,rv2) (,ru1 ,rv2) (,ru3 ,rv2) (,ru2 ,rv1) (,ru2 ,rv3)))))

(defun sample-surface (surface n &key number-of-points-u number-of-points-v)
  "If NUMBER-OF-POINTS-(U|V) is negative, take from the end.

N is the resolution -- a list for U and V, or an integer if they are equal.
The first (last) row (column) is not included: e.g. if NUMBER-OF-POINTS-U is
negative, the the V-directional border line at Umax is left out."
  (let* ((n (if (listp n) n (list n n)))
	 (lower (bss-lower-parameter surface))
	 (upper (bss-upper-parameter surface))
	 (step (mapcar #'(lambda (x y z) (/ (- x y) (1- z))) upper lower n))
	 (from-u (cond ((null number-of-points-u) 0)
		       ((> number-of-points-u 0) 1)
		       (t (+ (first n) number-of-points-u -1))))
	 (from-v (cond ((null number-of-points-v) 0)
		       ((> number-of-points-v 0) 1)
		       (t (+ (second n) number-of-points-v -1))))
	 (to-u (cond ((null number-of-points-u) (1- (first n)))
		     ((< number-of-points-u 0) (- (first n) 2))
		     (t number-of-points-u)))
	 (to-v (cond ((null number-of-points-v) (1- (second n)))
		     ((< number-of-points-v 0) (- (second n) 2))
		     (t number-of-points-v)))
	 (points (make-array (list (- to-u from-u -1) (- to-v from-v -1)))))
    (iter (for i from from-u to to-u)
	  (iter (for j from from-v to to-v)
		(for uv = (list (+ (first lower) (* (first step) i))
				(+ (second lower) (* (second step) j))))
		(setf (aref points (- i from-u) (- j from-v))
		      (bss-evaluate surface uv))))
    points))

(defun fair-xnode-in-one-direction (xnode resolution iteration distance
				    &key (u-direction t))
  (iter (with surface = (first xnode))
	(with 1st = (if u-direction 0 1))
	(with 2nd = (if u-direction 1 0))
	(with lower = (bss-lower-parameter surface))
	(with upper = (bss-upper-parameter surface))
	(with max = (1- (nth 2nd resolution)))
	(with left = (if u-direction 1 3))
	(with right = (if u-direction 2 4))
	(for i to max)
	(for p = (nth 2nd (affine-combine lower (/ i max) upper)))
	(for lcurve = (bss-get-surface-curve (nth left xnode) p
					     :u-curve u-direction))
	(for curve = (bss-get-surface-curve surface p :u-curve u-direction))
	(for rcurve = (bss-get-surface-curve (nth right xnode) p
					     :u-curve u-direction))
	(for start-cv = (bsc-curvature lcurve (bsc-upper-parameter lcurve)))
	(for end-cv = (bsc-curvature rcurve (bsc-lower-parameter rcurve)))
	(collect (bsc-faired-polygon curve
				     (nth 1st resolution) iteration distance
				     :from (nth 1st lower)
				     :to (nth 1st upper)
				     :start-curvature start-cv
				     :end-curvature end-cv))))

#+blind-mesh-blend
(defun blend-meshes (res a b)
  "TODO: The blending of the U and V meshes should be smarter."
  (let ((result (make-array res)))
    (dotimes (i (first res))
      (dotimes (j (second res))
	(setf (aref result i j)
	      (affine-combine (elt (nth j a) i) 0.5 (elt (nth i b) j)))))
    result))

(defun blend-meshes (res a b)
  (let ((midu (/ (first res) 2))
	(midv (/ (second res) 2))
	(result (make-array res)))
    (dotimes (i (first res))
      (dotimes (j (second res))
	(let* ((ad (expt (- (abs (- midu i)) midu) 2))
	       (bd (expt (- (abs (- midv j)) midv) 2))
	       (ratio (if (= (+ ad bd) 0) 0.5 (/ ad (+ ad bd)))))
	  (setf (aref result i j)
		(affine-combine (elt (nth j a) i) ratio (elt (nth i b) j))))))
    result))

(defun fair-xnode (xnode resolution iteration max-deviation)
  (let ((mesh-u (fair-xnode-in-one-direction xnode resolution iteration
					     max-deviation :u-direction t))
	(mesh-v (fair-xnode-in-one-direction xnode resolution iteration
					     max-deviation :u-direction nil)))
    (blend-meshes resolution mesh-u mesh-v)))

(defun bss-fit-engine (degree point-group-list &key
		       number-of-control-points-u number-of-control-points-v
		       knot-vector-u knot-vector-v
		       smoothness-functional optimize-parameters)
  (assert (not (or (and number-of-control-points-u knot-vector-u)
		   (and number-of-control-points-v knot-vector-v))))
  (let* ((sf (sf-create))
	 (arrays (iter (for group in point-group-list)
		       (collect (list (coerce (car group) 'double-float)
				      (/ (length (cdr group)) 5)
				      (sequence->double-array (cdr group)))))))
    (sf-set-degree-u sf (first degree))
    (sf-set-degree-v sf (second degree))
    (sf-set-closed-u sf nil)
    (sf-set-closed-v sf nil)
    (dolist (array arrays)
      (sf-add-point-group sf (first array) (second array) (third array)))
    (when number-of-control-points-u
      (sf-set-nr-ctrl-points-u sf number-of-control-points-u))
    (when number-of-control-points-v
      (sf-set-nr-ctrl-points-v sf number-of-control-points-v))
    (sf-set-optimize-parameters sf optimize-parameters)
    (sf-set-smoothness-functional sf (or smoothness-functional :smf-none))
    (flet ((set-sf-parameter (param-fn data)
	     (when data
	       (let ((size-of-data (length data))
		     (data-array (sequence->double-array data)))
		 (push (list nil size-of-data data-array) arrays)
		 (funcall param-fn sf size-of-data data-array)))))
      (set-sf-parameter #'sf-set-knot-vector-u knot-vector-u)
      (set-sf-parameter #'sf-set-knot-vector-v knot-vector-v))
    (unwind-protect
	 (if (eql (sf-fit sf) :success)
	     (bspline-surface-from-sf sf)
	     nil)
      (progn
	(sf-destroy sf)
	(dolist (array arrays) (foreign-free (third array)))))))

(defun bss-resembling-fit (surface points tolerance &key knot-vector)
  (let ((low (bss-lower-parameter surface))
	(high (bss-upper-parameter surface)))
    (bss-fit-engine
     (degrees surface)
     (list (cons tolerance
		 (uniform-parameter-points-2d points
					      (first low) (first high)
					      (second low) (second high))))
     :number-of-control-points-u (unless knot-vector
				   (array-dimension (control-net surface) 0))
     :number-of-control-points-v (unless knot-vector
				   (array-dimension (control-net surface) 1))
     :knot-vector-u (and knot-vector (first (knot-vectors surface)))
     :knot-vector-v (and knot-vector (second (knot-vectors surface)))
     :smoothness-functional :smf-none
     :optimize-parameters nil)))

(defun create-patch (up vp uendp vendp)
  "Create a set of points at the meeting point of UP's u border and VP's v
border. The u border is at v=0 if VENDP is false; similarly for the v border."
  (let* ((nu (array-dimension up 0))
	 (nv (array-dimension vp 1))
	 (result (make-array (list nu nv))))
    (dotimes (i nu)
      (dotimes (j nv)
	(flet ((get-dvec (array udir endp)
		 (macrolet ((index (axis i)
			      `(if endp
				   (- (array-dimension array ,axis) ,(1+ i))
				   ,i)))
		   (let ((p1 (aref array
				   (if udir (index 0 0) i)
				   (if udir j (index 1 0))))
			 (p2 (aref array
				   (if udir (index 0 1) i)
				   (if udir j (index 1 1)))))
		     (list p1 (v- p1 p2))))))
	  (let* ((du (get-dvec up nil vendp))
		 (dv (get-dvec vp t uendp))
		 (closest (closest-point (first du) (second du)
					 (first dv) (second dv))))
	    (setf (aref result i j)
		  (affine-combine (first closest) 0.5 (second closest)))))))
    result))

(defun fit-extended (surface percent resolution tolerance)
  "Extends the U and V domain of a surface at both sides by PERCENT%.
It fits on RESOLUTION points with TOLERANCE."
  (flet ((generate-knot (knot deg)
	   (let* ((lower (elt knot deg))
		  (upper (elt knot (- (length knot) deg 1)))
		  (plus (* percent (- upper lower) 0.01)))
	     (concatenate 'vector
			  (iter (repeat (1+ deg)) (collect (- lower plus)))
			  (subseq knot deg (- (length knot) deg))
			  (iter (repeat (1+ deg)) (collect (+ upper plus)))))))
    (let* ((knots (knot-vectors surface))
	   (degrees (degrees surface))
	   (new-knots (mapcar #'generate-knot knots degrees))
	   (points nil))
      (iter (with lower = (bss-lower-parameter surface))
	    (with upper = (bss-upper-parameter surface))
	    (with length = (mapcar #'- upper lower))
	    (with res = (if (listp resolution)
			    resolution
			    (list resolution resolution)))
	    (for i from 0 below (first res))
	    (for u = (+ (first lower) (/ (* i (first length))
					 (1- (first res)))))
	    (iter (for j from 0 below (second res))
		  (for v = (+ (second lower) (/ (* j (second length))
						(1- (second res)))))
		  (for p = (bss-evaluate surface (list u v)))
		  (setf points (append (list u v) p points))))
      (bss-fit-engine degrees (list (cons tolerance points))
		      :knot-vector-u (first new-knots)
		      :knot-vector-v (second new-knots)
		      :smoothness-functional :smf-crv
		      :optimize-parameters nil))))

(defun guess-sample-parameters (surface point-array &optional threshold)
  "For every point in the 2D POINT-ARRAY, project the point to SURFACE to
get an approximate parameterization.

When THRESHOLD is given, the algorithm doesn't include points whose
approximation deviates more than THRESHOLD.
Returns a list of the form \(U0 V0 X0 Y0 Z0 U1 V1 X1 Y1 Z1 ...).
TODO: Bogus parameters to BSS-PROJECT-POINT."
  (let ((acc nil)
	(count 0))
    (dotimes (j (array-dimension point-array 1))
      (dotimes (i (array-dimension point-array 0))
	(let* ((p (aref point-array i j))
	       (uv (bss-project-point surface p 10 '(10 10))) ; kutykurutty
	       (q (bss-evaluate surface uv)))
	  (unless (and threshold (> (vlength (v- p q)) threshold))
	    (setf acc (append uv p acc))
	    (incf count)))))
    (let ((all (apply #'* (array-dimensions point-array))))
      (when (< count all)
	(warn "GUESS-SAMPLE-PARAMETERS: Using only ~d points of ~d."
	      count all)))
    acc))

(defun suppressed-fit-xnode (xnode faired-points resolution held-points
			     loose-tolerance tight-tolerance patch-corners)
  "TODO: Bogus parameters to FIT-EXTENDED."
  (let* ((lower (bss-lower-parameter (first xnode)))
	 (upper (bss-upper-parameter (first xnode)))
	 (domain (mapcar #'- upper lower))
	 (held-lo-u (/ (first domain) (first (first resolution))))
	 (held-hi-u (/ (* (first domain) held-points)
		       (first (first resolution))))
	 (held-lo-v (/ (second domain) (second (first resolution))))
	 (held-hi-v (/ (* (second domain) held-points)
		       (second (first resolution))))
	 (lpoints (sample-surface (second xnode) (second resolution)
				  :number-of-points-u (- held-points)))
	 (rpoints (sample-surface (third xnode) (third resolution)
				  :number-of-points-u held-points))
	 (dpoints (sample-surface (fourth xnode) (fourth resolution)
				  :number-of-points-v (- held-points)))
	 (upoints (sample-surface (fifth xnode) (fifth resolution)
				  :number-of-points-v held-points))
	 (lengths (array-dimensions (control-net (first xnode))))
	 (extended (fit-extended (first xnode)
				 20 50 tight-tolerance))) ; kutykurutty
;;     (write-points2-vtk lpoints "results/left.vtk")
;;     (write-points2-vtk rpoints "results/right.vtk")
;;     (write-points2-vtk dpoints "results/bottom.vtk")
;;     (write-points2-vtk upoints "results/top.vtk")
;;     (write-points2-vtk faired-points "results/center.vtk")
    (flet ((samples-to-parameter-points (sample dir)
	     (let ((from (case dir
			   (ld (list (- (first lower) held-hi-u)
				     (- (second lower) held-hi-v)))
			   (rd (list (+ (first upper) held-lo-u)
				     (- (second lower) held-hi-v)))
			   (lu (list (- (first lower) held-hi-u)
				     (+ (second upper) held-lo-v)))
			   (ru (list (+ (first upper) held-lo-u)
				     (+ (second upper) held-lo-v)))))
		   (to (case dir
			 (ld (list (- (first lower) held-lo-u)
				   (- (second lower) held-lo-v)))
			 (rd (list (+ (first upper) held-hi-u)
				   (- (second lower) held-lo-v)))
			 (lu (list (- (first lower) held-lo-u)
				   (+ (second upper) held-hi-v)))
			 (ru (list (+ (first upper) held-hi-u)
				   (+ (second upper) held-hi-v))))))
	       (uniform-parameter-points-2d
		sample (first from) (first to) (second from) (second to)))))
      (bss-fit-engine
       '(3 3)
       (append
	(list (cons tight-tolerance
		    (guess-sample-parameters extended lpoints tight-tolerance))
	      (cons tight-tolerance
		    (guess-sample-parameters extended rpoints tight-tolerance))
	      (cons tight-tolerance
		    (guess-sample-parameters extended dpoints tight-tolerance))
	      (cons tight-tolerance
		    (guess-sample-parameters extended upoints tight-tolerance))
	      (cons loose-tolerance (uniform-parameter-points-2d
				     faired-points
				     (first lower) (first upper)
				     (second lower) (second upper))))
	(and patch-corners
	     (let ((ldpoints (create-patch lpoints dpoints nil nil))
		   (rdpoints (create-patch rpoints dpoints t nil))
		   (lupoints (create-patch lpoints upoints nil t))
		   (rupoints (create-patch rpoints upoints t t)))
;; 	       (write-points2-vtk ldpoints "results/ld.vtk")
;; 	       (write-points2-vtk rdpoints "results/rd.vtk")
;; 	       (write-points2-vtk lupoints "results/lu.vtk")
;; 	       (write-points2-vtk rupoints "results/ru.vtk")
	       (list (cons loose-tolerance
			   (samples-to-parameter-points ldpoints 'ld))
		     (cons loose-tolerance
			   (samples-to-parameter-points rdpoints 'rd))
		     (cons loose-tolerance
			   (samples-to-parameter-points lupoints 'lu))
		     (cons loose-tolerance
			   (samples-to-parameter-points rupoints 'ru))))))
       :number-of-control-points-u (first lengths)
       :number-of-control-points-v (second lengths)
       :smoothness-functional :smf-none
       :optimize-parameters nil))))

(defun just-a-surface-projection (original faired u v)
  (let ((n (bss-surface-normal original (list u v)))
	(p (bss-evaluate original (list u v)))
	(fp (bss-evaluate faired (list u v))))
    (v+ p (v* n (scalar-product (v- fp p) n)))))

(defun just-an-evaluation (original faired u v)
  (declare (ignore original))
  (bss-evaluate faired (list u v)))

(defun grid-cut (original faired resolution tolerance fn)
  (iter (with points = (make-array resolution))
	(with low = (bss-lower-parameter original))
	(with high = (bss-upper-parameter original))
	(with len = (mapcar (lambda (x y) (- y x)) low high))
	(for i from 0 below (first resolution))
	(for u = (+ (first low) (/ (* (first len) i) (1- (first resolution)))))
	(iter (for j from 0 below (second resolution))
	      (for v = (+ (second low)
			  (/ (* (second len) j) (1- (second resolution)))))
	      (setf (aref points i j) (funcall fn original faired u v)))
	(finally (return (bss-resembling-fit original points tolerance
					     :knot-vector t)))))

(defun fair-and-fit-xnode (xnode &key (resolution 300)
			   (iteration 100) (max-deviation 1000)
			   (loose-tolerance 0.1) (tight-tolerance 0.001)
			   (number-of-held-points 5)
			   no-fairing simple-fitting patch-corners cutting)
  "The missing corner regions are created (by an educated guess) if
`PATCH-CORNERS' is set.
The key `CUTTING' is a list of symbols.
* NIL  : no cutting
* GRID : grid cut with the function `JUST-A-PROJECTION'
* EVAL : grid cut with the function `JUST-AN-EVALUATION'
* ZAP  : ensure G0 connectivity (usually used with EVAL)
* G1   : ensure G1 connectivity (usually used with ZAP)"
  (let* ((res (xnode-resolution xnode resolution))
	 (faired (if no-fairing
		     (first (sample-surface (first xnode) (first res)))
		     (fair-xnode xnode (first res) iteration max-deviation))))
    (if simple-fitting
	(bss-resembling-fit (first xnode) faired loose-tolerance
			    :knot-vector nil)
	(let ((suppressed
	       (suppressed-fit-xnode xnode faired res
				     number-of-held-points
				     loose-tolerance tight-tolerance
				     patch-corners)))
	  (assert suppressed nil "Suppressed fit failed")
	  (if cutting
	      (let ((cut-sf (cond ((member 'eval cutting)
				   (grid-cut (first xnode) suppressed
					     (first res) tight-tolerance
					     #'just-an-evaluation))
				  ((member 'grid cutting)
				   (grid-cut (first xnode) suppressed
					     (first res) tight-tolerance
					     #'just-a-surface-projection)))))
		(let ((zap-sf
		       (if (member 'zap cutting)
			   (apply #'zap-to-surfaces cut-sf (rest xnode))
			   cut-sf)))
		  (if (member 'g1 cutting)
			   (ensure-g1-continuity zap-sf
						 (second xnode) (third xnode)
						 (fourth xnode) (fifth xnode)
						 (first res))
			   zap-sf)))
	      suppressed)))))

(defun five-surface-test (xnode filename &rest keys)
  (let ((faired (apply #'fair-and-fit-xnode xnode keys)))
    (when faired
      (write-rdn (cons faired (rest xnode)) filename)
      faired)))
