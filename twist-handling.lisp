;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

;;; TODO: factor the g1- and g2-specific functions

(defun plane-plane-intersection (p1 n1 p2 n2)
  "Given two planes by an arbitrary point and a normal vector,
this function returns their intersection as the pair \(POINT DIRECTION).

The returned point always has a zero in its coordinates, at the position
where the direction vector has the largest absolute value."
  (let* ((a (vnormalize (cross-product n1 n2)))
	 (k (iter (for i upfrom 0)
		  (for ai in (mapcar #'abs a))
		  (finding i maximizing ai)))
	 (indices (remove k '(0 1 2)))
	 (i (elt indices 0))
	 (j (elt indices 1))
	 (dn (- (* (elt n1 i) (elt n2 j)) (* (elt n1 j) (elt n2 i))))
	 (pn1 (scalar-product p1 n1))
	 (pn2 (scalar-product p2 n2))
	 (qi (- (* (/ (elt n2 j) dn) pn1) (* (/ (elt n1 j) dn) pn2)))
	 (qj (- (* (/ (elt n1 i) dn) pn2) (* (/ (elt n2 i) dn) pn1)))
	 (q (list 0 0 0)))
    (setf (elt q i) qi (elt q j) qj)
    (list q a)))

(defun project-to-line (p line)
  (let ((q1 (first line))
	(v1 (second line)))
    (v+ q1 (v* v1 (scalar-product (v- p q1) v1)))))

(defun project-to-plane (p plane)
  (let ((q1 (first plane))
	(v1 (second plane)))
    (v+ p (v* v1 (scalar-product (v- q1 p) v1)))))

(defun average-plane (p1 n1 p2 n2)
  (list (affine-combine p1 1/2 p2) (vnormalize (v+ n1 n2))))

(defparameter *minimum-angle-for-valid-intersection* (/ (* 5.0d0 pi) 180.0d0))
(defun closest-point-to-plane-plane-intersection (p p1 n1 p2 n2)
  (if (< (acos (scalar-product n1 n2)) *minimum-angle-for-valid-intersection*)
      (project-to-plane p (average-plane p1 n1 p2 n2))
      (project-to-line p (plane-plane-intersection p1 n1 p2 n2))))

;;; Not used - superseded by CLOSEST-POINT-TO-PLANE-PLANE-INTERSECTION
(defun between-line-and-point (p line)
  "Given a point P and a line as an arbitrary point and a direction,
returns the middle point between P and the closest point to P on the line."
  (affine-combine p 1/2 (project-to-line p line)))

(defun greville (d knots i)
  (/ (iter (for k from 0 below d)
	   (sum (elt knots (+ i k 1))))
     d))

(defun greville-abscissae (surface i j)
  (with-accessors ((d degrees) (k knot-vectors)) surface 
    (list (greville (first d) (first k) i)
	  (greville (second d) (second k) j))))

(defun twist-control-points (surface twist uendp vendp)
  (let ((n (array-dimensions (control-net surface))))
    (list (if uendp (- (first n) twist 1) twist)
	  (if vendp (- (second n) twist 1) twist))))

(defun side-points (surface twist uendp vendp)
  "TWIST is the frame number, ie. 1 for twists, 2 for inner twists, etc.

Returns the parameters of the reference points on the sides,
first the one on the u-boundary, then the one on the v-boundary."
  (let ((side-u (if uendp
		    (first (bss-upper-parameter surface))
		    (first (bss-lower-parameter surface))))
	(side-v (if vendp
		    (second (bss-upper-parameter surface))
		    (second (bss-lower-parameter surface)))))
    (destructuring-bind (cu cv) (twist-control-points surface twist uendp vendp)
      (destructuring-bind (gu gv) (greville-abscissae surface cu cv)
	(list (list gu side-v) (list side-u gv))))))

(defun twist-position (surface twist uendp vendp)
  (destructuring-bind (cu cv) (twist-control-points surface twist uendp vendp)
    (aref (control-net surface) cu cv)))

(defun ideal-g1-twist-position (surface u-surface v-surface uendp vendp)
  "U-SURFACE and V-SURFACE are the surface connecting along an u and a v
parameter line, respectively."
  (destructuring-bind (u-uv v-uv) (side-points surface 1 uendp vendp)
    (closest-point-to-plane-plane-intersection
     (twist-position surface 1 uendp vendp)
     (bss-evaluate u-surface u-uv)
     (bss-surface-normal u-surface u-uv)
     (bss-evaluate v-surface v-uv)
     (bss-surface-normal v-surface v-uv))))

(defun smoothed-g1-twist-position/average (surface twist uendp vendp)
  "Average of the neighboring control points."
  (destructuring-bind (cu cv) (twist-control-points surface twist uendp vendp)
    (let ((net (control-net surface)))
      (v* (v+ (aref net (1- cu) cv)
	      (aref net (1+ cu) cv)
	      (aref net cu (1- cv))
	      (aref net cu (1+ cv)))
	  1/4))))

(defun smoothed-twist-position/intersection (surface twist uendp vendp)
  "A pseudo-intersection of the neighboring control points"
  (destructuring-bind (cu cv) (twist-control-points surface twist uendp vendp)
    (flet ((line (p1 p2) (list p1 (vnormalize (v- p2 p1)))))
      (let ((net (control-net surface)))
	(destructuring-bind (p p-dir)
	    (line (aref net (1- cu) cv) (aref net (1+ cu) cv))
	  (destructuring-bind (q q-dir)
	      (line (aref net cu (1- cv)) (aref net cu (1+ cv)))
	    (destructuring-bind (x1 x2) (closest-point p p-dir q q-dir)
	      (affine-combine x1 1/2 x2))))))))

(defmacro weighted-vector-sum ((dimension) &body weight-expr-pairs)
  "DIMENSION should be <= 4."
  (let ((null-vector (subseq '(0 0 0 0) 0 dimension)))
    `(v* (v+ ,@(mapcar (lambda (pair)
			 (let ((weight (gensym)))
			   `(let ((,weight ,(first pair)))
			      (if (= ,weight 0)
				  ',null-vector
				  (v* ,(second pair) ,weight)))))
		       weight-expr-pairs))
	 (/ (+ ,@(mapcar #'first weight-expr-pairs))))))

(defmacro static-weighted-vector-sum (&body weight-expr-pairs)
  "Use this only if all the weights are explicit, static numbers."
  `(v* (v+ ,@(mapcan (lambda (pair)
		       (if (= (first pair) 0)
			   nil
			   `((v* ,(second pair) ,(first pair)))))
		     weight-expr-pairs))
       (/ (+ ,@(mapcar #'first weight-expr-pairs)))))

(defun next-g1-twist-position (surface u-surface v-surface uendp vendp)
  (static-weighted-vector-sum
    (2.0d0 (ideal-g1-twist-position surface u-surface v-surface uendp vendp))
    (1.0d0 (smoothed-twist-position/intersection surface 1 uendp vendp))
    (1.0d0 (twist-position surface 1 uendp vendp))))

(defun move-g1-twists (surface lsurface rsurface dsurface usurface)
  "WARNING: Changes the surface.
Assumes standard position scheme - visualize the surface such that
the lower u,v point is in the left bottom corner, with u increasing
to the right and v upwards. Then l/r/d/u stands for left, right,
down and up."
  (with-accessors ((net control-net)) surface
    (let ((nu (array-dimension net 0))
	  (nv (array-dimension net 1)))
      (setf (aref net 1 1)
	    (next-g1-twist-position surface dsurface lsurface nil nil)
	    (aref net 1 (- nv 2))
	    (next-g1-twist-position surface usurface lsurface nil t)
	    (aref net (- nu 2) 1)
	    (next-g1-twist-position surface dsurface rsurface t nil)
	    (aref net (- nu 2) (- nv 2))
	    (next-g1-twist-position surface usurface rsurface t t))
      t)))

(defun g1-continuity-with-twists
    (surface lsurface rsurface dsurface usurface res iteration)
  (let ((result (copy-bspline-surface surface)))
    (iter (repeat iteration)
	  (move-g1-twists result lsurface rsurface dsurface usurface)
	  (ensure-g1-one-side result lsurface (second res) :u-dir t :endp nil)
	  (ensure-g1-one-side result rsurface (second res) :u-dir t :endp t)
	  (ensure-g1-one-side result dsurface (first res) :u-dir nil :endp nil)
	  (ensure-g1-one-side result usurface (first res) :u-dir nil :endp t))
    result))

(defun g2-twist-altitude (surface master uv u-dir uendp vendp)
  "How much should the twist move in the surface normal direction,
according to the MASTER surface, such that their normal curvature will be
equal in the other parametric direction.

U-DIR is T if the other direction is the u parametric direction."
  (flet ((ufirst (lst) (if u-dir (first lst) (second lst)))
	 (usecond (lst) (if u-dir (second lst) (first lst))))
    (let* ((d1 (bss-evaluate surface uv :derivative (if u-dir '(1 0) '(0 1))))
	   (k-master (normal-curvature master uv d1))
	   (k-surface (normal-curvature surface uv d1))
	   (d1-square (scalar-product d1 d1))
	   (udegree (ufirst (degrees surface)))
	   (vdegree (usecond (degrees surface)))
	   (uknots (ufirst (knot-vectors surface)))
	   (vknots (usecond (knot-vectors surface)))
	   (uknots-length (length uknots))
	   (net (control-net surface))
	   (n (1- (array-dimension net 0)))
	   (m (1- (array-dimension net 1)))
	   (base (bspline-basis vknots (if u-dir
					   (if vendp (- m 2) 2)
					   (if uendp (- n 2) 2))
				vdegree (usecond uv))))
      (flet ((uknot (ui)
	       (if (if u-dir uendp vendp)
		   (elt uknots (- uknots-length ui 1))
		 (elt uknots ui))))
	(/ (* (- k-master k-surface) d1-square
	      (- (uknot (+ udegree 1)) (uknot 2))
	      (- (uknot (+ udegree 2)) (uknot 2)))
	   (* udegree (1- udegree) base))))))

(defun ideal-g2-twist-position (surface u-surface v-surface uendp vendp)
  "U-SURFACE and V-SURFACE are the surface connecting along an u and a v
parameter line, respectively."
  (destructuring-bind (u-uv v-uv) (side-points surface 2 uendp vendp)
    (let ((twist-position (twist-position surface 2 uendp vendp)))
      (closest-point-to-plane-plane-intersection
       twist-position 
       (v+ twist-position
	   (v* (bss-surface-normal surface u-uv)
	       (g2-twist-altitude surface u-surface u-uv nil uendp vendp)))
       (bss-surface-normal surface u-uv)
       (v+ twist-position
	   (v* (bss-surface-normal surface v-uv)
	       (g2-twist-altitude surface v-surface v-uv t uendp vendp)))
       (bss-surface-normal surface v-uv)))))

(defun next-g2-twist-position (surface u-surface v-surface uendp vendp)
  (static-weighted-vector-sum
    (2.0d0 (ideal-g2-twist-position surface u-surface v-surface uendp vendp))
    (1.0d0 (smoothed-twist-position/intersection surface 2 uendp vendp))
    (1.0d0 (twist-position surface 2 uendp vendp))))

(defun move-g2-twists (surface lsurface rsurface dsurface usurface)
  "WARNING: Changes the surface.
Assumes standard position scheme - visualize the surface such that
the lower u,v point is in the left bottom corner, with u increasing
to the right and v upwards. Then l/r/d/u stands for left, right,
down and up."
  (with-accessors ((net control-net)) surface
    (let ((nu (array-dimension net 0))
	  (nv (array-dimension net 1)))
      (setf (aref net 2 2)
	    (next-g2-twist-position surface dsurface lsurface nil nil)
	    (aref net 2 (- nv 3))
	    (next-g2-twist-position surface usurface lsurface nil t)
	    (aref net (- nu 3) 2)
	    (next-g2-twist-position surface dsurface rsurface t nil)
	    (aref net (- nu 3) (- nv 3))
	    (next-g2-twist-position surface usurface rsurface t t))
      t)))

(defun g2-continuity-with-twists
    (surface lsurface rsurface dsurface usurface res iteration
     &key (algorithm 'minimal-deviation))
  (let ((result (copy-bspline-surface surface)))
    (iter (repeat iteration)
	  (move-g2-twists result lsurface rsurface dsurface usurface)
	  (ensure-g2-one-side result lsurface (second res)
			      :u-dir t :endp nil :algorithm algorithm)
	  (ensure-g2-one-side result rsurface (second res)
			      :u-dir t :endp t :algorithm algorithm)
	  (ensure-g2-one-side result dsurface (first res)
			      :u-dir nil :endp nil :algorithm algorithm)
	  (ensure-g2-one-side result usurface (first res)
			      :u-dir nil :endp t :algorithm algorithm))
    result))

#+nil
(defun test (xnode filename &key (iteration 1) (resolution '(100 100)))
  (let ((result (g2-continuity-with-twists
		 (first xnode) (second xnode) (third xnode)
		 (fourth xnode) (fifth xnode) resolution iteration)))
    (write-rbn (cons result xnode) filename)))
