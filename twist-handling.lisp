(in-package :cl-nurbs-tests)

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
	 (qi (- (* (/ (elt n2 j) dn) (scalar-product p1 n1))
		(* (/ (elt n1 j) dn) (scalar-product p2 n2))))
	 (qj (- (* (/ (elt n1 i) dn) (scalar-product p2 n2))
		(* (/ (elt n2 i) dn) (scalar-product p1 n1))))
	 (q (list 0 0 0)))
    (setf (elt q i) qi (elt q j) qj)
    (list q a)))

(defun between-line-and-point (p q1 v1)
  "Given a point P and a line as an arbitrary point and a direction,
returns the middle point between P and the closest point to P on the line."
  (let ((q (v+ q1 (v* v1 (scalar-product (v- p q1) v1)))))
    (affine-combine p 1/2 q)))

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
    (destructuring-bind (gu gv)
	(apply #'greville-abscissae surface
	       (twist-control-points surface twist uendp vendp))
      (list (list gu side-v) (list side-u gv)))))

(defun g1-twist-position (surface u-surface v-surface uendp vendp)
  "U-SURFACE and V-SURFACE are the surface connecting along an u and a v
parameter line, respectively."
  (destructuring-bind (u-uv v-uv) (side-points surface 1 uendp vendp)
    (apply #'between-line-and-point
	   (apply #'aref (control-net surface)
		  (twist-control-points surface 1 uendp vendp))
	   (plane-plane-intersection
	    (bss-evaluate u-surface u-uv)
	    (bss-surface-normal u-surface u-uv)
	    (bss-evaluate v-surface v-uv)
	    (bss-surface-normal v-surface v-uv)))))

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
	    (g1-twist-position surface dsurface lsurface nil nil)
	    (aref net 1 (- nv 2))
	    (g1-twist-position surface usurface lsurface nil t)
	    (aref net (- nu 2) 1)
	    (g1-twist-position surface dsurface rsurface t nil)
	    (aref net (- nu 2) (- nv 2))
	    (g1-twist-position surface usurface rsurface t t))
      t)))

(defun g1-continuity/twists (surface lsurface rsurface dsurface usurface res)
  (let ((result (copy-bspline-surface surface)))
    (move-g1-twists result lsurface rsurface dsurface usurface)
    (ensure-g1-one-side result lsurface (second res) :u-dir t :endp nil)
    (ensure-g1-one-side result rsurface (second res) :u-dir t :endp t)
    (ensure-g1-one-side result dsurface (first res) :u-dir nil :endp nil)
    (ensure-g1-one-side result usurface (first res) :u-dir nil :endp t)
    result))
