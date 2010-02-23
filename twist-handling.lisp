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

;;; Only a utility
(defun perturbe-g1 (surface &optional (epsilon 1.0d-4))
  (let ((result (copy-bspline-surface surface)))
    (iter (with net = (control-net result))
	  (with (n m) = (array-dimensions net))
	  (for i from 0 below n)
	  (iter (for j from 0 below m)
		(when (and (or (= i 1) (= j 1) (= i (- n 2)) (= j (- m 2)))
			   (not (or (= i 0) (= j 0) (= i (1- n)) (= j (1- m)))))
		  (setf (aref net i j)
			(mapcar (lambda (x) (+ x (- (* 2 (random epsilon)) epsilon)))
				(aref net i j))))))
    result))

;;; Not used
(defun box* (box n)
  (let ((center (affine-combine (first box) 1/2 (second box))))
    (list (v+ center (v* (v- (first box) center) n))
	  (v+ center (v* (v- (second box) center) n)))))

(defun g1-line-segment (surface u-surface v-surface uendp vendp)
  "U-SURFACE and V-SURFACE are the surface connecting along an u and a v
parameter line, respectively."
  (destructuring-bind (u-uv v-uv) (side-points surface 1 uendp vendp)
    (destructuring-bind (cu cv) (twist-control-points surface 1 uendp vendp)
      (let* ((net (control-net surface))
	     (line (plane-plane-intersection
		    (bss-evaluate u-surface u-uv)
		    (bss-surface-normal u-surface u-uv)
		    (bss-evaluate v-surface v-uv)
		    (bss-surface-normal v-surface v-uv)))
	     (twist (aref net cu cv))
	     (len (* 2 (+ (point-distance twist (aref net (1- cu) cv))
			  (point-distance twist (aref net cu (1- cv))))))
	     (center (project-to-line twist line)))
	(list (v- center (v* (second line) len))
	      (v+ center (v* (second line) len)))))))

(defun twists-vtk-g1 (xnode filename &key (res '(100 100)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 2.0~%Twists~%ASCII~%DATASET POLYDATA~%~%")
    (let* ((net (control-net (first xnode)))
	   (nm (array-dimensions net))
	   (points (reduce #'* nm)))
      (format s "POINTS ~d float~%" (+ points 16))
      (iter (for i from 0 below (first nm))
	    (iter (for j from 0 below (second nm))
		  (format s "~{~f~^ ~}~%" (aref (control-net (first xnode)) i j))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (second xnode) (second res) :u-dir t :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 1 1)))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (fourth xnode) (first res) :u-dir nil :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 1 1)))
      (format s "~{~{~f~^ ~}~%~}"
	      (g1-line-segment (first xnode) (fourth xnode) (second xnode) nil nil))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (second xnode) (second res) :u-dir t :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 1 (- (second nm) 2))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (fifth xnode) (first res) :u-dir nil :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 1 (- (second nm) 2))))
      (format s "~{~{~f~^ ~}~%~}"
	      (g1-line-segment (first xnode) (fifth xnode) (second xnode) nil t))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (third xnode) (second res) :u-dir t :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 2) (- (second nm) 2))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (fifth xnode) (first res) :u-dir nil :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 2) (- (second nm) 2))))
      (format s "~{~{~f~^ ~}~%~}"
	      (g1-line-segment (first xnode) (fifth xnode) (third xnode) t t))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (third xnode) (second res) :u-dir t :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 2) 1)))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g1-one-side temp (fourth xnode) (first res) :u-dir nil :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 2) 1)))
      (format s "~{~{~f~^ ~}~%~}"
	      (g1-line-segment (first xnode) (fourth xnode) (third xnode) t nil))
      (flet ((index (i j) (+ (* i (second nm)) j)))
	(let ((lines (+ (* 2 (reduce #'* (mapcar #'1- nm))) (reduce #'+ (mapcar #'1- nm)) 20)))
	  (format s "~%LINES ~d ~d~%" lines (* lines 3))
	  (iter (for i from 1 below (first nm))
		(iter (for j from 1 below (second nm))
		      (format s "2 ~d ~d~%" (index i j) (index i (1- j)))
		      (format s "2 ~d ~d~%" (index i j) (index (1- i) j))))
	  (iter (for i from 0 below (1- (first nm)))
		(format s "2 ~d ~d~%" (index i 0) (index (1+ i) 0)))
	  (iter (for j from 0 below (1- (second nm)))
		(format s "2 ~d ~d~%" (index 0 j) (index 0 (1+ j))))
	  (flet ((connect-to-twist (twist uendp vendp)
		   (let ((cu (if uendp
				 (list (- (first nm) 1) (- (first nm) 2))
				 '(0 1)))
			 (cv (if vendp
				 (list (- (second nm) 1) (- (second nm) 2))
				 '(0 1))))
		     (format s "2 ~d ~d~%2 ~d ~d~%"
			     (index (first cu) (second cv)) twist
			     (index (second cu) (first cv)) twist))))
	    (connect-to-twist (+ points 0) nil nil)
	    (connect-to-twist (+ points 1) nil nil)
	    (connect-to-twist (+ points 4) nil t)
	    (connect-to-twist (+ points 5) nil t)
	    (connect-to-twist (+ points 8) t t)
	    (connect-to-twist (+ points 9) t t)
	    (connect-to-twist (+ points 12) t nil)
	    (connect-to-twist (+ points 13) t nil))
	  (format s "2 ~d ~d~%" (+ points 2) (+ points 3))
	  (format s "2 ~d ~d~%" (+ points 6) (+ points 7))
	  (format s "2 ~d ~d~%" (+ points 10) (+ points 11))
	  (format s "2 ~d ~d~%" (+ points 14) (+ points 15))
	  (format s "~%CELL_DATA ~d~%SCALARS color float 1~%LOOKUP_TABLE my_colors~%"
		  lines)
	  (iter (repeat (- lines 20)) (format s "0~%"))
	  (iter (repeat 4) (format s "0.333~%0.333~%0.666~%0.666~%"))
	  (format s "1~%1~%1~%1~%")
	  (format s "LOOKUP_TABLE my_colors 4~%~
                     0 0 0 1~%~
                     1 0 0 1~%~
                     0 1 0 1~%~
                     0 0 1 1~%"))))))

(defun g2-line-segment (surface u-surface v-surface uendp vendp)
  "U-SURFACE and V-SURFACE are the surface connecting along an u and a v
parameter line, respectively."
  (destructuring-bind (u-uv v-uv) (side-points surface 2 uendp vendp)
    (destructuring-bind (cu cv) (twist-control-points surface 2 uendp vendp)
      (let* ((net (control-net surface))
	     (twist-position (twist-position surface 2 uendp vendp))
	     (line (plane-plane-intersection
		    (v+ twist-position
			(v* (bss-surface-normal surface u-uv)
			    (g2-twist-altitude surface u-surface u-uv nil uendp vendp)))
		    (bss-surface-normal surface u-uv)
		    (v+ twist-position
			(v* (bss-surface-normal surface v-uv)
			    (g2-twist-altitude surface v-surface v-uv t uendp vendp)))
		    (bss-surface-normal surface v-uv)))
	     (twist (aref net cu cv))
	     (len (* 2 (+ (point-distance twist (aref net (1- cu) cv))
			  (point-distance twist (aref net cu (1- cv))))))
	     (center (project-to-line twist line)))
	(list (v- center (v* (second line) len))
	      (v+ center (v* (second line) len)))))))

(defun twists-vtk-g2 (xnode filename &key (res '(100 100)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 2.0~%Twists~%ASCII~%DATASET POLYDATA~%~%")
    (let* ((net (control-net (first xnode)))
	   (nm (array-dimensions net))
	   (points (reduce #'* nm)))
      (format s "POINTS ~d float~%" (+ points 16))
      (iter (for i from 0 below (first nm))
	    (iter (for j from 0 below (second nm))
		  (format s "~{~f~^ ~}~%" (aref (control-net (first xnode)) i j))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (second xnode) (second res) :u-dir t :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 2 2)))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (fourth xnode) (first res) :u-dir nil :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 2 2)))
      (format s "~{~{~f~^ ~}~%~}"
	      (g2-line-segment (first xnode) (fourth xnode) (second xnode) nil nil))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (second xnode) (second res) :u-dir t :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 2 (- (second nm) 3))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (fifth xnode) (first res) :u-dir nil :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) 2 (- (second nm) 3))))
      (format s "~{~{~f~^ ~}~%~}"
	      (g2-line-segment (first xnode) (fifth xnode) (second xnode) nil t))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (third xnode) (second res) :u-dir t :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 3) (- (second nm) 3))))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (fifth xnode) (first res) :u-dir nil :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 3) (- (second nm) 3))))
      (format s "~{~{~f~^ ~}~%~}"
	      (g2-line-segment (first xnode) (fifth xnode) (third xnode) t t))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (third xnode) (second res) :u-dir t :endp t :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 3) 2)))
      (let ((temp (copy-bspline-surface (first xnode))))
	(ensure-g2-one-side temp (fourth xnode) (first res) :u-dir nil :endp nil :move-twists t)
	(format s "~{~f~^ ~}~%" (aref (control-net temp) (- (first nm) 3) 2)))
      (format s "~{~{~f~^ ~}~%~}"
	      (g2-line-segment (first xnode) (fourth xnode) (third xnode) t nil))
      (flet ((index (i j) (+ (* i (second nm)) j)))
	(let ((lines (+ (* 2 (reduce #'* (mapcar #'1- nm))) (reduce #'+ (mapcar #'1- nm)) 20)))
	  (format s "~%LINES ~d ~d~%" lines (* lines 3))
	  (iter (for i from 1 below (first nm))
		(iter (for j from 1 below (second nm))
		      (format s "2 ~d ~d~%" (index i j) (index i (1- j)))
		      (format s "2 ~d ~d~%" (index i j) (index (1- i) j))))
	  (iter (for i from 0 below (1- (first nm)))
		(format s "2 ~d ~d~%" (index i 0) (index (1+ i) 0)))
	  (iter (for j from 0 below (1- (second nm)))
		(format s "2 ~d ~d~%" (index 0 j) (index 0 (1+ j))))
	  (flet ((connect-to-twist (twist uendp vendp)
		   (let ((cu (if uendp
				 (list (- (first nm) 2) (- (first nm) 3))
				 '(1 2)))
			 (cv (if vendp
				 (list (- (second nm) 2) (- (second nm) 3))
				 '(1 2))))
		     (format s "2 ~d ~d~%2 ~d ~d~%"
			     (index (first cu) (second cv)) twist
			     (index (second cu) (first cv)) twist))))
	    (connect-to-twist (+ points 0) nil nil)
	    (connect-to-twist (+ points 1) nil nil)
	    (connect-to-twist (+ points 4) nil t)
	    (connect-to-twist (+ points 5) nil t)
	    (connect-to-twist (+ points 8) t t)
	    (connect-to-twist (+ points 9) t t)
	    (connect-to-twist (+ points 12) t nil)
	    (connect-to-twist (+ points 13) t nil))
	  (format s "2 ~d ~d~%" (+ points 2) (+ points 3))
	  (format s "2 ~d ~d~%" (+ points 6) (+ points 7))
	  (format s "2 ~d ~d~%" (+ points 10) (+ points 11))
	  (format s "2 ~d ~d~%" (+ points 14) (+ points 15))
	  (format s "~%CELL_DATA ~d~%SCALARS color float 1~%LOOKUP_TABLE my_colors~%"
		  lines)
	  (iter (repeat (- lines 20)) (format s "0~%"))
	  (iter (repeat 4) (format s "0.333~%0.333~%0.666~%0.666~%"))
	  (format s "1~%1~%1~%1~%")
	  (format s "LOOKUP_TABLE my_colors 4~%~
                     0 0 0 1~%~
                     1 0 0 1~%~
                     0 1 0 1~%~
                     0 0 1 1~%"))))))
