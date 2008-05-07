;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

#|

(defun bss-curvature-uv (surface uv &key (u-dir t))
  "Returns the U/V directional curvature at point UV."
  (let ((d1 (bss-evaluate surface uv :derivative (if u-dir '(1 0) '(0 1))))
	(d2 (bss-evaluate surface uv :derivative (if u-dir '(2 0) '(0 2)))))
    (safe-/ (vlength (cross-product d1 d2)) (expt (vlength d1) 3))))

(defun bss-arc-length-sampling (surface n &key (u-dir t))
  'TODO)

(defun bss-target-curvature (surface parameters resolution iteration &key
			     (u-dir t) curvature-start curvature-end)
  "TODO: This should be done on quasi-arc-length parametric samples."
  (let* ((result (make-array resolution :element-type 'real))
	 (low (bss-lower-parameter surface))
	 (hi (bss-upper-parameter surface))
	 (len (mapcar #'- hi low))
	 (nu (first resolution))
	 (nv (second resolution)))
    (dotimes (i nu)
      (let ((u (+ (first low) (/ (* i (first len)) (1- nu)))))
	(dotimes (j nv)
	  (let ((v (+ (second low)
		      (/ (* j (second len)) (1- nv)))))
	    (setf (aref result i j)
		  (cond ((and (= i 0) u-dir curvature-start)
			 (elt curvature-start j))
			((and (= j 0) (not u-dir) curvature-start)
			 (elt curvature-start i))
			((and (= i (1- nu)) u-dir curvature-end)
			 (elt curvature-end j))
			((and (= j (1- nv)) (not u-dir) curvature-end)
			 (elt curvature-end i))
			(t (bss-curvature-uv surface (list u v)
					     :u-dir u-dir))))))))
    (iter (repeat iteration)
	  (with tmp = (make-array resolution :element-type 'real))
	  (iter (for i from 1 below (1- nu))
		(iter (for j from 1 below (1- nv))
		      (setf (aref tmp i j)
			    (/ (+ (aref result (1- i) j)
				  (aref result (1+ i) j)
				  (aref result i (1- j))
				  (aref result i (1+ j)))
			       4))))
	  (iter (for i from 1 below (1- nu))
		(iter (for j from 1 below (1- nv))
		      (setf (aref result i j)
			    (aref tmp i j)))))
    result))

(defun bss-fair-from-corner (surface parameters targets resolution iteration
			     max-deviation &key u-endp v-endp)
  (let ((last-u (make-array resolution))
	(last-v (make-array resolution))
	(result (make-array resolution))
	(nu (first resolution))
	(nv (second resolution)))
    ;; init RESULT
    ;; init LAST-U, LAST-V
    (dotimes (i (max nu nv))
      (when (< i nu)
	)
      (when (< i nv)
	)
      (dotimes (j (1+ i))
	(when (and (< i nu) (< j nv))
	  )
	(when (and (/= i j) (< i nv) (< j nu))
	  )))
    result))

(defun blend-corner-meshes (usvs uevs usve ueve)
  "Blends the four arrays.

USVS was created from the lower U and lower V parameter corner,
UEVS was created from the higher U and lower V parameter corner, etc."
  (let ((nu (array-dimension usvs 0))
	(nv (array-dimension usvs 1))
	(result (make-array (array-dimensions usvs))))
    (dotimes (i nu)
      (let ((u (blend-function (/ i (1- nu)))))
	(dotimes (j nv)
	  (let ((v (blend-function (/ j (1- nv)))))
	    (setf (aref result i j)
		  (v+ (v* (aref usvs i j) (/ (+ (- 1.0 u) (- 1.0 v)) 4.0))
		      (v* (aref uevs i j) (/ (+ u (- 1.0 v)) 4.0))
		      (v* (aref usve i j) (/ (+ (- 1.0 u) v) 4.0))
		      (v* (aref ueve i j) (/ (+ u v) 4.0))))))))
    result))

(defun bss-fair-from-corners (surface resolution iteration max-deviation &key
			      curvature-start-u curvature-start-v
			      curvature-end-u curvature-end-v)
  "TODO: This should create a quasi-arc-length parameterization first."
  (let* ((parameters (bss-arc-length-sampling ...))
	 (target-u (bss-target-curvature surface parameters resolution
					 iteration :u-dir t
					 :curvature-start curvature-start-u
					 :curvature-end curvature-end-u))
	 (target-v (bss-target-curvature surface parameters resolution
					 iteration :u-dir nil
					 :curvature-start curvature-start-v
					 :curvature-end curvature-end-v))
	 (targets (list target-u target-v)))
    (blend-corner-meshes
     (bss-fair-from-corner surface parameters targets resolution
			   iteration max-deviation
			   :u-endp nil :v-endp nil)
     (bss-fair-from-corner surface parameters targets resolution
			   iteration max-deviation
			   :u-endp t :v-endp nil)
     (bss-fair-from-corner surface parameters targets resolution
			   iteration max-deviation
			   :u-endp nil :v-endp t)
     (bss-fair-from-corner surface parameters targets resolution
			   iteration max-deviation
			   :u-endp t :v-endp t))))

(defun fair-xnode-from-corners (xnode resolution iteration max-deviation)
  (let* ((low (bss-lower-parameter surface))
	 (hi (bss-upper-parameter surface))
	 (len (mapcar #'- hi low))
	 ;; kutykurutty: ezeket is a PARAMETERS szerint kene mintavetelezni
	 (curv-left (iter (with u = (first low))
			  (for j from 0 below (second resolution))
			  (for v = (+ (second low)
				      (/ (* j (second len))
					 (1- (second resolution)))))
			  (collect
			      (bss-curvature-uv (second xnode) (list u v)
						:u-dir t))))
	 (curv-right (iter (with u = (first hi))
			   (for j from 0 below (second resolution))
			   (for v = (+ (second low)
				       (/ (* j (second len))
					  (1- (second resolution)))))
			   (collect
			       (bss-curvature-uv (third xnode) (list u v)
						 :u-dir t))))
	 (curv-bottom (iter (with v = (second low))
			    (for i from 0 below (first resolution))
			    (for u = (+ (first low)
					(/ (* i (first len))
					   (1- (first resolution)))))
			    (collect
				(bss-curvature-uv (fourth xnode) (list u v)
						  :u-dir nil))))
	 (curve-top (iter (with v = (second hi))
			  (for i from 0 below (first resolution))
			  (for u = (+ (first low)
				      (/ (* i (first len))
					 (1- (first resolution)))))
			  (collect
			      (bss-curvature-uv (fifth xnode) (list u v)
						:u-dir nil)))))
    (bss-fair-from-corners (first xnode) resolution iteration max-deviation
			   :curvature-start-u curv-left
			   :curvature-end-u curv-right
			   :curvature-start-v curv-bottom
			   :curvature-end-v curv-top)))

|#
