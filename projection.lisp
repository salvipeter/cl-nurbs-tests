;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun bsc-projection-starting-value (curve point res)
  (iter (with lower = (bsc-lower-parameter curve))
	(with upper = (bsc-upper-parameter curve))
	(with length = (- upper lower))
	(for i from 0 below res)
	(for u = (+ lower (/ (* i length) (1- res))))
	(finding u minimizing (vlength (v- point (bsc-evaluate curve u))))))

(defun bsc-project-point (curve point iterations search-resolution &optional
			  (distance-tolerance 0.0) (cosine-tolerance 0.0))
  "Returns the parameter of CURVE's closest point to POINT.

The function uses the Newton-Raphson method. (NURBS Book, pp. 230-232)
SEARCH-RESOLUTION parameters are checked for a suitable initial value."
  (let ((u0 (bsc-projection-starting-value curve point search-resolution))
	(lower (bsc-lower-parameter curve))
	(upper (bsc-upper-parameter curve)))
    (iter (repeat iterations)
	  (for last first nil then u)
	  (for u first u0 then
	       (min (max (- u (/ (scalar-product d deviation)
				 (+ (scalar-product d2 deviation)
				    (scalar-product d d))))
			 lower)
		    upper))
	  (for p = (bsc-evaluate curve u))
	  (for d = (bsc-evaluate curve u :derivative 1))
	  (for d2 = (bsc-evaluate curve u :derivative 2))
	  (for deviation = (v- p point))
	  (when (or (<= (vlength deviation) distance-tolerance)
		    (<= (/ (abs (scalar-product d deviation))
			   (* (vlength d) (vlength deviation)))
			cosine-tolerance)
		    (and last
			 (<= (vlength (v* d (- u last))) distance-tolerance)))
	    (leave u))
	  (finally (return u)))))

(defun bss-projection-starting-value (surface point res)
  (iter (with lower = (bss-lower-parameter surface))
	(with upper = (bss-upper-parameter surface))
	(with length = (mapcar #'- upper lower))
	(for i from 0 below (first res))
	(for u = (+ (first lower) (/ (* i (first length)) (1- (first res)))))
	(for (min v) =
	     (iter (for j from 0 below (second res))
		   (for v = (+ (second lower) (/ (* j (second length))
						 (1- (second res)))))
		   (for distance = (vlength
				    (v- point
					(bss-evaluate surface (list u v)))))
		   (finding (list distance v) minimizing distance)))
	(finding (list u v) minimizing min)))

(defun bss-projection-delta (surface uv r)
  (let* ((du (bss-evaluate surface uv :derivative '(1 0)))
	 (dv (bss-evaluate surface uv :derivative '(0 1)))
	 (du2 (bss-evaluate surface uv :derivative '(2 0)))
	 (duv (bss-evaluate surface uv :derivative '(1 1)))
	 (dv2 (bss-evaluate surface uv :derivative '(0 2)))
	 (k (make-array '(2 1) :initial-contents
			`((,(- (scalar-product r du)))
			  (,(- (scalar-product r dv))))))
	 (J (make-array
	     '(2 2) :initial-contents
	     `((,(+ (scalar-product du du) (scalar-product r du2))
		 ,(+ (scalar-product du dv) (scalar-product r duv)))
	       (,(+ (scalar-product du dv) (scalar-product r duv))
		 ,(+ (scalar-product dv dv) (scalar-product r dv2))))))
	 (delta (matrix-multiplication (matrix-inverse-2x2 J) k)))
    (list (aref delta 0 0) (aref delta 1 0))))

(defun bss-project-point (surface point iterations search-resolution &optional
			  (distance-tolerance 0.0) (cosine-tolerance 0.0))
  "Returns the parameters of SURFACE's closest point to POINT.

The function uses the Newton-Raphson method. (NURBS Book, pp. 232-234)
\(FIRST SEARCH-RESOLUTION) * \(SECOND SEARCH-RESOLUTION) parameters
are checked for a suitable initial value."
  (let ((uv0 (bss-projection-starting-value surface point search-resolution))
	(lower (bss-lower-parameter surface))
	(upper (bss-upper-parameter surface)))
    (iter (repeat iterations)
	  (for last first nil then uv)
	  (for uv first uv0 then
	       (mapcar (lambda (x next low upp) (min (max (+ x next) low) upp))
		       uv (bss-projection-delta surface uv deviation)
		       lower upper))
	  (for p = (bss-evaluate surface uv))
	  (for du = (bss-evaluate surface uv :derivative '(1 0)))
	  (for dv = (bss-evaluate surface uv :derivative '(0 1)))
	  (for deviation = (v- p point))
	  (when (or (<= (vlength deviation) distance-tolerance)
		    (<= (/ (abs (scalar-product du deviation))
			   (* (vlength du) (vlength deviation)))
			cosine-tolerance)
		    (<= (/ (abs (scalar-product dv deviation))
			   (* (vlength dv) (vlength deviation)))
			cosine-tolerance)
		    (and last
			 (<= (+ (vlength (v* du (- (first uv) (first last))))
				(vlength (v* dv (- (second uv) (second last)))))
			     distance-tolerance)))
	    (leave uv))
	  (finally (return uv)))))
