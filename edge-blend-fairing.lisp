;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

;;; Assumptions:
;;; BLEND is a list of three surfaces: SURFACE, LEFT and RIGHT, that are
;;; connected along the u parametric direction in the order LEFT-SURFACE-RIGHT.
;;; The continuity of the connection is G1.

;;; Goal:
;;; Fair SURFACE while enhancing its connection to G2 continuity, by
;;; the process fairing-to-mesh => fastened-fit => g2.

;;; The resulting surface should be better than
;;; - fastened-fit => g2
;;; - g2 => fairing-in-place

;;; (fairing-to-mesh: salvi/kobbelt; fairing-in-place: krr)

(defun eb-resolution (blend resolution)
  (let* ((low (bss-lower-parameter (first blend)))
	 (high (bss-upper-parameter (first blend)))
	 (uv (mapcar (lambda (a b) (/ (+ a b) 2)) low high))
	 (cu1 (bss-get-surface-curve (second blend) (second uv) :u-curve t))
	 (cu2 (bss-get-surface-curve (first blend) (second uv) :u-curve t))
	 (cu3 (bss-get-surface-curve (third blend) (second uv) :u-curve t))
	 (ru1 (bsc-estimate-arc-length cu1))
	 (ru2 (bsc-estimate-arc-length cu2))
	 (ru3 (bsc-estimate-arc-length cu3))
	 (total (+ ru1 ru2 ru3)))
    (mapcar (lambda (ru)
	      (list (round (* ru (first resolution)) total)
		    (second resolution)))
	    (list ru2 ru1 ru3))))

(defun eb-fair-salvi-one-direction (blend resolution iteration distance
				    &key (u-direction t))
  "Takes curvatures from the sides, when U-DIRECTION is T."
  (iter (with surface = (first blend))
	(with 1st = (if u-direction 0 1))
	(with 2nd = (if u-direction 1 0))
	(with lower = (bss-lower-parameter surface))
	(with upper = (bss-upper-parameter surface))
	(with max = (1- (nth 2nd resolution)))
	(for i to max)
	(for p = (nth 2nd (affine-combine lower (/ i max) upper)))
	(for (start-curvature end-curvature) =
	     (when u-direction
	       (let* ((lcurve (bss-get-surface-curve
			       (second blend) p :u-curve t))
		      (rcurve (bss-get-surface-curve
			       (third blend) p :u-curve t)))
		 (list (bsc-curvature lcurve (bsc-upper-parameter lcurve))
		       (bsc-curvature rcurve (bsc-lower-parameter rcurve))))))
	(for curve = (bss-get-surface-curve surface p :u-curve u-direction))
	(collect (bsc-faired-polygon curve
				     (nth 1st resolution) iteration distance
				     :from (nth 1st lower)
				     :to (nth 1st upper)
				     :start-curvature start-curvature
				     :end-curvature end-curvature))))

(defun eb-fair-salvi (blend resolution iteration max-deviation)
  "TODO: a better blending function?"
  (let ((mesh-u (eb-fair-salvi-one-direction blend resolution iteration
					     max-deviation :u-direction t))
	(mesh-v (eb-fair-salvi-one-direction blend resolution iteration
					     max-deviation :u-direction nil)))
    (blend-meshes resolution mesh-u mesh-v)))

(defun eb-fastened-fit (blend points resolution held-points
			tight-tolerance loose-tolerance)
  "TODO: Bogus parameters to FIT-EXTENDED."
  (let ((lower (bss-lower-parameter (first blend)))
	(upper (bss-upper-parameter (first blend)))
	(lpoints (sample-surface (second blend) (second resolution)
				 :number-of-points-u (- held-points)))
	(rpoints (sample-surface (third blend) (third resolution)
				 :number-of-points-u held-points))
	(lengths (array-dimensions (control-net (first blend))))
	(extended (fit-extended (first blend)
				20 50 tight-tolerance))) ; kutykurutty
    (bss-fit-engine
     '(3 3)
     (append
      (list (cons tight-tolerance
		  (guess-sample-parameters extended lpoints tight-tolerance))
	    (cons tight-tolerance
		  (guess-sample-parameters extended rpoints tight-tolerance))
	    (cons loose-tolerance (uniform-parameter-points-2d
				   points (first lower) (first upper)
				   (second lower) (second upper)))))
     :number-of-control-points-u (first lengths)
     :knot-vector-v (second (knot-vectors (first blend)))
     :smoothness-functional :smf-none
     :optimize-parameters nil)))

(defun eb-cut (large original)
  (destructuring-bind (minu minv) (bss-lower-parameter original)
    (declare (ignore minv))
    (destructuring-bind (maxu maxv) (bss-upper-parameter original)
      (declare (ignore maxv))
      (bss-subsurface-one-direction large minu maxu t))))

(defun eb-g0 (surface left right)
  (let ((result surface))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 left (first (bss-upper-parameter left))
			 :u-curve nil)))
    (setf result (bss-reverse-parameterization result :u t :v nil))
    (setf result
	  (zap-to-curve result
			(bss-get-surface-curve
			 right (first (bss-lower-parameter right))
			 :u-curve nil)))
    (setf result (bss-reverse-parameterization result :u t :v nil))
    result))

(defun eb-g1 (surface left right res)
  (let ((result (copy-bspline-surface surface)))
    (ensure-g1-one-side result left (second res) :u-dir t :endp nil
			:move-twists t)
    (ensure-g1-one-side result right (second res) :u-dir t :endp t
			:move-twists t) 
    result))

(defun eb-g2 (surface left right res &key (algorithm 'minimal-deviation))
  (let ((result (copy-bspline-surface surface)))
    (ensure-g2-one-side result left (second res)
			:u-dir t :endp nil :algorithm algorithm :move-twists t)
    (ensure-g2-one-side result right (second res)
			:u-dir t :endp t :algorithm algorithm :move-twists t) 
    result))

(defun eb-g (surface left right res g &key (algorithm 'minimal-deviation))
  (let* ((g0 (if (member 0 g) (eb-g0 surface left right) surface))
	 (g1 (if (member 1 g) (eb-g1 g0 left right res) g0))
	 (g2 (if (member 2 g)
		 (eb-g2 g1 left right res :algorithm algorithm)
		 g1)))
    g2))

;;; User functions

(defun edge-blend-fairing (blend resolution &key
			   (fairing-method 'salvi) (iteration 100)
			   (max-deviation 0.1) (held-points 5)
			   (tight-tolerance 0.0001) (loose-tolerance 0.001))
  (let* ((res (eb-resolution blend resolution))
	 (faired (case fairing-method
		   (salvi (eb-fair-salvi blend (first res) iteration
					 max-deviation))
		   (kobbelt (bss-kobbelt-faired-points (first blend)
						       (first res) iteration))
		   (t (sample-surface (first blend) (first res))))))
    (eb-g (eb-cut (eb-fastened-fit blend faired res held-points
				   tight-tolerance loose-tolerance)
		  (first blend))
	  (second blend) (third blend) (first res) '(0 1 2))))

(defun edge-blend-krr-g2 (blend &key (halving 'inner) (resolution '(100 100))
			  (iteration 100) (algorithm 'minimal-deviation))
  "Do just a KRR inside the surface and enhance the G2 continuity."
  (let ((fn (ecase halving
	      (none #'identity)
	      (end #'halve-end-intervals)
	      (inner #'halve-inner-intervals)
	      (all #'halve-all-intervals))))
    (g-fair-krr-additive
     (apply #'eb-g2 (funcall fn (first blend))
	    (append (rest blend) (list resolution :algorithm algorithm)))
     iteration 2)))

;;; Test

#|
(defparameter *blend*
  (subseq (read-rbn "/home/salvi/project/cl-nurbs/models/xnode.rbn") 0 3))

(write-rbn (edge-blend-fairing *blend* '(1000 100)) "/tmp/fair-fit-g2.rbn")

(write-rbn (edge-blend-fairing *blend* '(1000 100) :fairing-method nil)
	   "/tmp/fit-g2.rbn")

(write-rbn (edge-blend-krr-g2 *blend*) "/tmp/g2-krr.rbn")
|#
