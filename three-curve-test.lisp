(in-package :cl-nurbs)

(defun split-to-three (curve u1 u2)
  "Split CURVE to three subcurves at the parameter positions U1 and U2."
  (let* ((left-rest (bsc-split-curve curve u1))
	 (middle-right (bsc-split-curve (second left-rest) u2)))
    (cons (first left-rest) middle-right)))

(defun sample-curve (curve n)
  (iter (with upper = (bsc-upper-parameter curve))
	(with lower = (bsc-lower-parameter curve))
	(with step  = (/ (- upper lower) (1- n)))
	(for i from 0 below n)
	(collect (bsc-evaluate curve (+ lower (* step i))))))

(defun fair-middle-curve (split resolution iteration max-deviation)
  (let* ((left-curve (first split))
	 (curve (second split))
	 (right-curve (third split))
	 (parameters (arc-length-sampling curve
					  (bsc-lower-parameter curve)
					  (bsc-upper-parameter curve)
					  resolution))
	 (start-curvature (bsc-curvature left-curve
					 (bsc-upper-parameter left-curve)))
	 (end-curvature (bsc-curvature right-curve
				       (bsc-lower-parameter right-curve)))
	 (curvatures (target-curvature curve parameters iteration
				       :start-value start-curvature
				       :end-value end-curvature))
	 (integrate-from-left (integrate curve parameters curvatures
					 max-deviation :from-right nil))
	 (integrate-from-right (integrate curve parameters curvatures
					  max-deviation :from-right t)))
    (blend-points integrate-from-left (nreverse integrate-from-right))))

(defun split-curve-resolution (split total-resolution)
  (let ((left (bsc-estimate-arc-length (first split)
				       (bsc-lower-parameter (first split))
				       (bsc-upper-parameter (first split))))
	(middle (bsc-estimate-arc-length (second split)
					 (bsc-lower-parameter (second split))
					 (bsc-upper-parameter (second split))))
	(right (bsc-estimate-arc-length (third split)
					(bsc-lower-parameter (third split))
					(bsc-upper-parameter (third split)))))
    (let ((total (+ left middle right)))
      (mapcar (lambda (x) (round (* total-resolution x) total))
	      (list left middle right)))))

(defun bsc-fit-engine (dimension degree point-group-list &key
		       number-of-control-points knot-vector
		       smoothness-functional optimize-parameters
		       start-point end-point
		       start-tangent end-tangent
		       start-curvature end-curvature)
  (assert (not (and number-of-control-points knot-vector)))
  (let* ((gcf (gcf-create dimension))
	 (arrays (iter (for group in point-group-list)
		       (collect (list (coerce (car group) 'double-float)
				      (/ (length (cdr group)) (1+ dimension))
				      (sequence->double-array (cdr group)))))))
    (gcf-set-degree gcf degree)
    (gcf-set-closed gcf nil)
    (dolist (array arrays)
      (gcf-add-point-group gcf (first array) (second array) (third array)))
    (when number-of-control-points
      (gcf-set-nr-ctrl-points gcf number-of-control-points))
    (gcf-set-optimize-parameters gcf optimize-parameters)
    (gcf-set-smoothness-functional gcf (or smoothness-functional :smf-none))
    (flet ((set-gcf-parameter (param-fn data &key with-size)
	     (when data
	       (let ((size-of-data (length data))
		     (data-array (sequence->double-array data)))
		 (push (list nil size-of-data data-array) arrays)
		 (if with-size
		     (funcall param-fn gcf size-of-data data-array)
		     (funcall param-fn gcf data-array))))))
      (set-gcf-parameter #'gcf-set-knot-vector knot-vector :with-size t)
      (set-gcf-parameter #'gcf-set-start-point start-point)
      (set-gcf-parameter #'gcf-set-end-point end-point)
      (set-gcf-parameter #'gcf-set-start-tangent start-tangent)
      (set-gcf-parameter #'gcf-set-end-tangent end-tangent)
      (set-gcf-parameter #'gcf-set-start-curvature start-curvature)
      (set-gcf-parameter #'gcf-set-end-curvature end-curvature))
    (unwind-protect
	 (if (eql (gcf-fit gcf) :success)
	     (bspline-curve-from-gcf gcf)
	     nil)
      (progn
	(gcf-destroy gcf)
	(dolist (array arrays) (foreign-free (third array)))))))

(defun bsc-resembling-fit (curve points tolerance)
  (let ((low (bsc-lower-parameter curve))
	(high (bsc-upper-parameter curve)))
    (bsc-fit-engine
     (bsc-dimension curve) (degree curve)
     (list (cons tolerance
		 (uniform-parameter-points points
					   (bsc-lower-parameter curve)
					   (bsc-upper-parameter curve))))
     :number-of-control-points (length (control-points curve))
;;      :knot-vector (knot-vector curve)
     :smoothness-functional :smf-none
     :optimize-parameters nil
     :start-point (elt points 0)
     :end-point (elt points (1- (length points)))
     :start-curvature (list (bsc-curvature curve low))
     :end-curvature (list (bsc-curvature curve high)))))

(defun suppressed-fit-middle-curve (split faired-points resolution
				    number-of-held-points
				    loose-tolerance tight-tolerance)
  (let* ((lpoints (sample-curve (first split) (first resolution)))
	 (rpoints (sample-curve (third split) (third resolution)))
	 (points (uniform-parameter-points
		  (append lpoints (coerce faired-points 'list) rpoints)
		  (bsc-lower-parameter (first split))
		  (bsc-upper-parameter (third split))))
	 (nr-lpoints (length lpoints))
	 (nr-faired (length faired-points)))
    (bsc-fit-engine
     2 3
     (list (cons tight-tolerance
		 (subseq points (* (- nr-lpoints number-of-held-points) 3)
			 (* nr-lpoints 3)))
	   (cons loose-tolerance
		 (subseq points (* nr-lpoints 3)
			 (* (+ nr-lpoints nr-faired) 3)))
	   (cons tight-tolerance
		 (subseq points (* (+ nr-lpoints nr-faired) 3)
			 (* (+ nr-lpoints nr-faired number-of-held-points)
			    3))))
     :number-of-control-points (+ (length (control-points (second split))) 2)
     :smoothness-functional :smf-none
     :optimize-parameters nil
     :start-point (elt lpoints (- nr-lpoints number-of-held-points))
     :end-point (elt rpoints (1- number-of-held-points)))))

(defun just-a-projection (original faired u)
  (let ((n (bsc-2d-normal original u))
	(p (bsc-evaluate original u))
	(fp (bsc-evaluate faired u)))
    (v+ p (v* n (scalar-product (v- fp p) n)))))

(defun grid-cut (original faired resolution tolerance)
  (iter (with low = (bsc-lower-parameter original))
	(with high = (bsc-upper-parameter original))
	(with len = (- high low))
	(for i from 0 below resolution)
	(for u = (+ low (/ (* len i) (1- resolution))))
	(collect (just-a-projection original faired u) into points)
	(finally (return (bsc-resembling-fit original points tolerance)))))

(defun fair-and-fit-middle-curve (split &key (resolution 100)
				  (iteration 100) (max-deviation 1000)
				  (loose-tolerance 0.1) (tight-tolerance 0.001)
				  (number-of-held-points 5)
				  no-fairing simple-fitting)
  (let* ((res (split-curve-resolution split resolution))
	 (faired (if no-fairing
		     (sample-curve (second split) (second res))
		     (fair-middle-curve split (second res)
					iteration max-deviation))))
    (if simple-fitting
	(bsc-resembling-fit (second split) faired loose-tolerance)
	(grid-cut (second split)
		  (suppressed-fit-middle-curve split faired res
					       number-of-held-points
					       loose-tolerance tight-tolerance)
		  (second res) tight-tolerance))))

(defun three-curve-test (curve u1 u2 extrusion filename &rest keys)
  (let* ((split (split-to-three curve u1 u2))
	 (faired (apply #'fair-and-fit-middle-curve split keys))
	 (new-split (list (first split) faired (third split))))
    (when faired
      (write-rdn (mapcar (lambda (curve) (bsc-extrude curve extrusion))
			 new-split)
		 filename)
      new-split)))

(defun three-curve-original (curve u1 u2 extrusion filename)
  (let ((split (split-to-three curve u1 u2)))
    (write-rdn (mapcar (lambda (curve) (bsc-extrude curve extrusion))
		       split)
	       filename)
    split))

(defun write-point-data (points filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "~{~{~a~^ ~}~%~}" (coerce points 'list))))
