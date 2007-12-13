(in-package :cl-nurbs)

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
Returns a list, the first element is an array containing the sampled values,
the second element is the parameter interval ((U1 V1) (U2 V2))."
  (let* ((n (if (listp n) n (list n n)))
	 (lower (bss-lower-parameter surface))
	 (upper (bss-upper-parameter surface))
	 (step (mapcar #'(lambda (x y z) (/ (- x y) (1- z))) upper lower n))
	 (from-u (if (or (null number-of-points-u) (> number-of-points-u 0))
		     0 (+ (first n) number-of-points-u)))
	 (from-v (if (or (null number-of-points-v) (> number-of-points-v 0))
		     0 (+ (second n) number-of-points-v)))
	 (to-u (if (or (null number-of-points-u) (< number-of-points-u 0))
		   (1- (first n)) (1- number-of-points-u)))
	 (to-v (if (or (null number-of-points-v) (< number-of-points-v 0))
		   (1- (second n)) (1- number-of-points-v)))
	 (points (make-array (list (- to-u from-u -1) (- to-v from-v -1)))))
    (iter (for i from from-u to to-u)
	  (iter (for j from from-v to to-v)
		(for uv = (list (+ (first lower) (* (first step) i))
				(+ (second lower) (* (second step) j))))
		(setf (aref points (- i from-u) (- j from-v))
		      (bss-evaluate surface uv))))
    (list points
	  (mapcar (lambda (uv)
		    (mapcar (lambda (x lower step) (+ lower (* step x)))
			    uv lower step))
		  `((,from-u ,from-v) (,to-u ,to-v))))))

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

(defun fair-xnode (xnode resolution iteration max-deviation)
  "TODO: The blending of the U and V meshes should be better."
  (let ((mesh-u (fair-xnode-in-one-direction xnode resolution iteration
					     max-deviation :u-direction t))
	(mesh-v (fair-xnode-in-one-direction xnode resolution iteration
					     max-deviation :u-direction nil))
	(result (make-array resolution)))
    (dotimes (i (first resolution))
      (dotimes (j (second resolution))
	(setf (aref result i j) (affine-combine (elt (nth j mesh-u) i)
						0.5
						(elt (nth i mesh-v) j)))))
    result))

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

(defun bss-resembling-fit (surface points tolerance)
  (let ((low (bss-lower-parameter surface))
	(high (bss-upper-parameter surface)))
    (bss-fit-engine
     (degrees surface)
     (list (cons tolerance
		 (uniform-parameter-points-2d points
					      (first low) (first high)
					      (second low) (second high))))
     :number-of-control-points-u (array-dimension (control-net surface) 0)
     :number-of-control-points-v (array-dimension (control-net surface) 1)
     ;;    :knot-vector-u (first (knot-vectors surface))
     ;;    :knot-vector-v (second (knot-vectors surface))
     :smoothness-functional :smf-none
     :optimize-parameters nil)))

(defun write-2d-gnuplot-points (points filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let ((n1 (array-dimension points 0))
	  (n2 (array-dimension points 1)))
      (dotimes (i n1)
	(dotimes (j n2)
	  (format s "~{~f~^ ~}~%" (aref points i j)))))))

(defun uniform-parameter-points-2d-inner (points &optional
					  (start-u 0.0) (end-u 1.0)
					  (start-v 0.0) (end-v 1.0))
  "Returns a list of the form (U1 V1 X1 Y1 Z1 U2 V2 X2...),
where U1..UN/V1..VN are equidistant points of the [(U1,V1), (U2,V2)] interval.
The border points are not included."
  (let ((n (array-dimension points 0))
	(m (array-dimension points 1))
	(len-u (- end-u start-u))
	(len-v (- end-v start-v)))
    (iter (for i from 1 below (1- n))
	  (with ppts)
	  (iter (for j from 1 below (1- m))
		(setf ppts (append ppts
				   (list (+ (/ (* len-u i) (1- n)) start-u)
					 (+ (/ (* len-v j) (1- m)) start-v))
				   (copy-list (aref points i j)))))
	  (finally (return ppts)))))

(defun suppressed-fit-xnode (xnode faired-points resolution held-points
			     loose-tolerance tight-tolerance)
  (let ((lpoints (sample-surface (second xnode) (second resolution)
				 :number-of-points-u (- held-points)))
	(rpoints (sample-surface (third xnode) (third resolution)
				 :number-of-points-u held-points))
	(dpoints (sample-surface (fourth xnode) (fourth resolution)
				 :number-of-points-v (- held-points)))
	(upoints (sample-surface (fifth xnode) (fifth resolution)
				 :number-of-points-v held-points))
	(lower (bss-lower-parameter (first xnode)))
	(upper (bss-upper-parameter (first xnode)))
	(lengths (array-dimensions (control-net (first xnode)))))
;;     (write-2d-gnuplot-points (first lpoints) "results/left.pts")
;;     (write-2d-gnuplot-points (first rpoints) "results/right.pts")
;;     (write-2d-gnuplot-points (first dpoints) "results/bottom.pts")
;;     (write-2d-gnuplot-points (first upoints) "results/top.pts")
;;     (write-2d-gnuplot-points faired-points "results/center.pts")
    (flet ((samples-to-parameter-points (sample)
	     (let* ((interval (second sample))
		    (from (first interval))
		    (to (second interval)))
	       (uniform-parameter-points-2d-inner (first sample)
						  (first from) (first to)
						  (second from) (second to)))))
      (bss-fit-engine
       '(3 3)
       (list (cons tight-tolerance (samples-to-parameter-points lpoints))
	     (cons tight-tolerance (samples-to-parameter-points rpoints))
 	     (cons tight-tolerance (samples-to-parameter-points dpoints))
 	     (cons tight-tolerance (samples-to-parameter-points upoints))
	     (cons loose-tolerance (uniform-parameter-points-2d
				    faired-points
				    (first lower) (first upper)
				    (second lower) (second upper))))
       :number-of-control-points-u (+ (first lengths) 2)
       :number-of-control-points-v (+ (second lengths) 2)
       :smoothness-functional :smf-none
       :optimize-parameters nil))))

(defun just-a-surface-projection (original faired u v)
  (let ((n (bss-surface-normal original (list u v)))
	(p (bss-evaluate original (list u v)))
	(fp (bss-evaluate faired (list u v))))
    (v+ p (v* n (scalar-product (v- fp p) n)))))

(defun grid-cut (original faired resolution tolerance)
  (iter (with points = (make-array resolution))
	(with low = (bss-lower-parameter original))
	(with high = (bss-upper-parameter original))
	(with len = (mapcar (lambda (x y) (- y x)) low high))
	(for i from 0 below (first resolution))
	(for u = (+ (first low) (/ (* (first len) i) (1- (first resolution)))))
	(iter (for j from 0 below (second resolution))
	      (for v = (+ (second low)
			  (/ (* (second len) j) (1- (second resolution)))))
	      (setf (aref points i j)
		    (just-a-surface-projection original faired u v)))
	(finally (return (bss-resembling-fit original points tolerance)))))

(defun fair-and-fit-xnode (xnode &key (resolution 300)
			   (iteration 100) (max-deviation 1000)
			   (loose-tolerance 0.1) (tight-tolerance 0.001)
			   (number-of-held-points 5)
			   no-fairing simple-fitting)
  (let* ((res (xnode-resolution xnode resolution))
	 (faired (if no-fairing
		     (first (sample-surface (first xnode) (first res)))
		     (fair-xnode xnode (first res) iteration max-deviation))))
    (if simple-fitting
	(bss-resembling-fit (first xnode) faired loose-tolerance)
	(let ((suppressed
	       (suppressed-fit-xnode xnode faired res
				     (+ number-of-held-points 2)
				     loose-tolerance tight-tolerance)))
	  (grid-cut (first xnode) suppressed (first res) tight-tolerance)))))

(defun five-surface-test (xnode filename &rest keys)
  (let ((faired (apply #'fair-and-fit-xnode xnode keys)))
    (when faired
      (write-rdn (cons faired (rest xnode)) filename)
      faired)))
