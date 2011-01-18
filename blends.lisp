(in-package :cl-nurbs-tests)

;; (eval-when (:compile-toplevel :load-toplevel :execute)
;;   (asdf:operate 'asdf:load-op 'trivial-shell))

(defparameter *paper* '(500 500))
(defparameter *point-size* 5)
(defparameter *font-size* 10)
(defparameter *resolution* 100)
(defparameter *curvature-comb-scale* -0.01)

(defparameter *red* '(1 0 0))
(defparameter *green* '(0 1 0))
(defparameter *blue* '(0 0 1))

(defvar *exponent* 2)
(defvar *ribbon-multiplier-start* 0.5d0)
(defvar *ribbon-multiplier-end* 0.5d0)

(defun init-font (stream)
  (format stream "/Courier findfont~%~d scalefont~%setfont~%"
	  *font-size*))

(defun title (stream str &key (color '(0 0 0)) (window 0))
  (format stream "~{~f ~}setrgbcolor~%" color)
  (let ((x0 (if (oddp window) 1 0))
	(y0 (if (> window 1) 1 0)))
    (flet ((convert (p)
	     (destructuring-bind (x y) p
	       (list (* (first *paper*) (/ (+ x0 x) 2.0d0))
		     (* (second *paper*) (/ (+ y0 y) 2.0d0))))))
      (format stream "newpath~%~{~f ~}moveto~%(~a) show~%"
	      (convert '(0.2 1)) str))))

(defun point (stream p &key (color '(0 0 0)) (window 0))
  (format stream "~{~f ~}setrgbcolor~%" color)
  (let ((x0 (if (oddp window) 1 0))
	(y0 (if (> window 1) 1 0))
	(a (/ *point-size* 2)))
    (flet ((convert (p)
	     (destructuring-bind (x y) p
	       (list (* (first *paper*) (/ (+ x0 x) 2.0d0))
		     (* (second *paper*) (/ (+ y0 y) 2.0d0))))))
      (destructuring-bind (x y) (convert p)
	(format stream "newpath~%~
                        ~f ~f moveto~%~
                        ~f ~f lineto~%~
                        ~f ~f lineto~%~
                        ~f ~f lineto~%~
                        ~f ~f lineto~%~
                        stroke~%"
		(- x a) (- y a) (+ x a) (- y a) (+ x a) (+ y a)
		(- x a) (+ y a) (- x a) (- y a))))))

(defun arc (stream p r from to &key (color '(0 0 0)) (window 0))
  (format stream "~{~f ~}setrgbcolor~%" color)
  (let ((x0 (if (oddp window) 1 0))
	(y0 (if (> window 1) 1 0)))
    (flet ((convert (p)
	     (destructuring-bind (x y) p
	       (list (* (first *paper*) (/ (+ x0 x) 2.0d0))
		     (* (second *paper*) (/ (+ y0 y) 2.0d0))))))
      (destructuring-bind (x y) (convert p)
	(format stream "newpath~%~
                        ~f ~f ~f ~f ~f arc~%~
                        stroke~%"
		x y (* (first *paper*) 0.5 r) from to)))))

(defun autozoom (points &optional (margin 0.1))
  (let ((xmin (reduce #'min (mapcar #'first points)))
	(xmax (reduce #'max (mapcar #'first points)))
	(ymin (reduce #'min (mapcar #'second points)))
	(ymax (reduce #'max (mapcar #'second points))))
    (let ((range (max (- xmax xmin) (- ymax ymin))))
      (lambda (p)
	(destructuring-bind (x y) p
	  (list (+ margin (* (- 1 (* 2 margin)) (/ (- x xmin) range)))
		(+ margin (* (- 1 (* 2 margin)) (/ (- y ymin) range)))))))))

(defun polyline (stream points &key (color '(0 0 0)) (window 0) autozoom)
  (format stream "~{~f ~}setrgbcolor~%" color)
  (let ((x0 (if (oddp window) 1 0))
	(y0 (if (> window 1) 1 0))
	(scale-fn (when autozoom (autozoom points))))
    (flet ((convert (p)
	     (destructuring-bind (x y) (if scale-fn (funcall scale-fn p) p)
	       (list (* (first *paper*) (/ (+ x0 x) 2.0d0))
		     (* (second *paper*) (/ (+ y0 y) 2.0d0))))))
      (format stream "newpath~%~{~f ~}moveto~%~{~{~f ~}lineto~%~}stroke~%"
	      (convert (first points)) (mapcar #'convert (rest points))))))

(defun planar-curve-curvature (curve-fn u)
  "CURVE-FN takes two parameters, u and derivative, eg. (CURVE-FN 0.3 1)."
  (let ((d1 (funcall curve-fn u 1))
	(d2 (funcall curve-fn u 2)))
    (safe-/ (scalar-product d1 (list (- (second d2)) (first d2)))
	    (expt (vlength d1) 3))))

(defun hermite-blend-function (type side u &key (derivative 0))
  "SIDE is one of (START END), TYPE is one of (POINT TANGENT)."
  (ecase type
    (point (ecase side
	     (start (ecase derivative
		      (0 (+ (* 2 u u u) (* -3 u u) 1))
		      (1 (- (* 6 u u) (* 6 u)))
		      (2 (- (* 12 u) 6))))
	     (end (ecase derivative
		    (0 (+ (* -2 u u u) (* 3 u u)))
		    (1 (- (* 6 u) (* 6 u u)))
		    (2 (- 6 (* 12 u)))))))
    (tangent (ecase side
	       (start (ecase derivative
			(0 (+ (* u u u) (* -2 u u) u))
			(1 (- (* 3 u u) (* 4 u) -1))
			(2 (- (* 6 u) 4))))
	       (end (ecase derivative
		      (0 (- (* u u) (* u u u)))
		      (1 (- (* 2 u) (* 3 u u)))
		      (2 (- 2 (* 6 u)))))))))

(defun hermite (p0 t0 p1 t1 u &key (derivative 0))
  (v+ (v* p0 (hermite-blend-function 'point 'start u :derivative derivative))
      (v* t0 (hermite-blend-function 'tangent 'start u :derivative derivative))
      (v* p1 (hermite-blend-function 'point 'end u :derivative derivative))
      (v* t1 (hermite-blend-function 'tangent 'end u :derivative derivative))))

(defun write-hermite (p0 t0 p1 t1 filename)
  (let (curve curvature derivative blend blend-derivative blend2 blend2-derivative)
    (iter (for i from 0 below *resolution*)
	  (for u = (/ i (1- *resolution*)))
	  (with curve-fn = (lambda (x d) (hermite p0 t0 p1 t1 x :derivative d)))
	  (for p = (funcall curve-fn u 0))
	  (collect p into tcurve)
	  (when (= (mod i 10) 5)
	    (let ((d (funcall curve-fn u 1)))
	      (collect (list p (v+ p (v* (vnormalize (list (- (second d)) (first d)))
					 (planar-curve-curvature curve-fn u)
					 *curvature-comb-scale*)))
		       into tcurvature)))
	  (collect (funcall curve-fn u 1) into tderivative)
	  (collect (list u (hermite-blend-function 'point 'start u))
		   into tblend)
	  (collect (list u (hermite-blend-function 'point 'start u :derivative 1))
		   into tblend-derivative)
	  (collect (list u (hermite-blend-function 'point 'end u))
		   into tblend2)
	  (collect (list u (hermite-blend-function 'point 'end u :derivative 1))
		   into tblend2-derivative)
	  (finally (setf curve tcurve curvature tcurvature derivative tderivative
			 blend tblend blend-derivative tblend-derivative
			 blend2 tblend2 blend2-derivative tblend2-derivative)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      (init-font s)
      ;; Curve
      (title s "Curve" :window 0)
      (point s p0 :color *red* :window 0)
      (point s (v+ p0 (v* t0 1/3)) :color *green* :window 0)
      (polyline s (list p0 (v+ p0 (v* t0 1/3))) :color *blue* :window 0)
      (point s p1 :color *red* :window 0)
      (point s (v+ p1 (v* t1 1/3)) :color *green* :window 0)
      (polyline s (list p1 (v+ p1 (v* t1 1/3))) :color *blue* :window 0)
      (polyline s curve :window 0)
      (iter (for line in curvature)
	    (polyline s line :window 0))
      ;; Derivative
      (title s "Derivative" :window 1)
      (polyline s derivative :window 1 :autozoom t)
      ;; Blend function
      (title s "Blend function" :window 2)
      (polyline s blend :window 2 :autozoom t)
      ;; (polyline s blend2 :color *blue* :window 2 :autozoom t)
      ;; Blend function derivative
      (title s "Blend function derivative" :window 3)
      (polyline s blend-derivative :window 3 :autozoom t)
      ;; (polyline s blend2-derivative :color *blue* :window 3 :autozoom t)
      ;; End of file
      (format s "showpage~%"))))

#+nil
(let ((*resolution* 200))
  (write-hermite '(0.3 0.3) '(0.1 0.7) '(0.8 0.3) '(-0.1 0.8) "/tmp/proba.ps"))

(defun tomi-blend-function (u &key (derivative 0))
  (let ((n *exponent*)
	(d0 u)
	(d1 (- 1 u)))
    (case derivative
      (0 (/ (expt d1 n) (+ (expt d0 n) (expt d1 n))))
      (1 (/ (- (+ (* n (expt d1 (1- n)) (+ (expt d0 n) (expt d1 n)))
		  (* (expt d1 n) n (- (expt d0 (1- n)) (expt d1 (1- n))))))
	    (expt (+ (expt d0 n) (expt d1 n)) 2)))
      (2 (/ (- (* (- (* n (1- n) (expt d1 (- n 2)) (+ (expt d0 n) (expt d1 n)))
		     (* (expt d1 n) n (1- n) (+ (expt d0 (- n 2)) (expt d1 (- n 2)))))
		  (expt (+ (expt d0 n) (expt d1 n)) 2))
	       (* (+ (* n (expt d1 (1- n)) (+ (expt d0 n) (expt d1 n)))
		     (* (expt d1 n) n (- (expt d0 (1- n)) (expt d1 (1- n)))))
		  (* 2 (+ (expt d0 n) (expt d1 n))
		     n (- (expt d0 (1- n)) (expt d1 (1- n))))))
	    (expt (+ (expt d0 n) (expt d1 n)) 4))))))

(defun tomi (p0 t0 p1 t1 u &key (derivative 0))
  (let ((d0 u)
	(d1 (- 1 u)))
    (case derivative
      (0 (v+ (v* (v+ p0 (v* t0 *ribbon-multiplier-start* d0))
		 (tomi-blend-function d0))
	     (v* (v+ p1 (v* t1 *ribbon-multiplier-end* d1))
		 (tomi-blend-function d1))))
      (1 (v+ (v* t0 *ribbon-multiplier-start* (tomi-blend-function d0))
	     (v* (v+ p0 (v* t0 *ribbon-multiplier-start* d0))
		 (tomi-blend-function d0 :derivative 1))
	     (v* t1 *ribbon-multiplier-end* -1 (tomi-blend-function d1))
	     (v* (v+ p1 (v* t1 *ribbon-multiplier-end* d1))
		 -1 (tomi-blend-function d1 :derivative 1))))
      (2 (v+ (v* t0 2 *ribbon-multiplier-start* (tomi-blend-function d0 :derivative 1))
	     (v* (v+ p0 (v* t0 *ribbon-multiplier-start* d0))
		 (tomi-blend-function d0 :derivative 2))
	     (v* t1 2 *ribbon-multiplier-end* (tomi-blend-function d1 :derivative 1))
	     (v* (v+ p1 (v* t1 *ribbon-multiplier-end* d1))
		 (tomi-blend-function d1 :derivative 2)))))))

(defun write-tomi (p0 t0 p1 t1 filename)
  (let (curve curvature derivative blend blend-derivative)
    (iter (for i from 0 below *resolution*)
	  (for u = (/ i (1- *resolution*)))
	  (with curve-fn = (lambda (x d) (tomi p0 t0 p1 t1 x :derivative d)))
	  (for p = (funcall curve-fn u 0))
	  (collect p into tcurve)
	  (when (= (mod i 10) 5)
	    (let ((d (funcall curve-fn u 1)))
	      (collect (list p (v+ p (v* (vnormalize (list (- (second d)) (first d)))
					 (planar-curve-curvature curve-fn u)
					 *curvature-comb-scale*)))
		       into tcurvature)))
	  (collect (funcall curve-fn u 1) into tderivative)
	  (collect (list u (tomi-blend-function u)) into tblend)
	  (collect (list u (tomi-blend-function u :derivative 1)) into tblend-derivative)
	  (finally (setf curve tcurve curvature tcurvature derivative tderivative
			 blend tblend blend-derivative tblend-derivative)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      (init-font s)
      ;; Curve
      (title s "Curve" :window 0)
      (point s p0 :color *red* :window 0)
      (point s (v+ p0 (v* t0 1/3)) :color *green* :window 0)
      (polyline s (list p0 (v+ p0 (v* t0 1/3))) :color *blue* :window 0)
      (point s p1 :color *red* :window 0)
      (point s (v+ p1 (v* t1 1/3)) :color *green* :window 0)
      (polyline s (list p1 (v+ p1 (v* t1 1/3))) :color *blue* :window 0)
      (polyline s curve :window 0)
      (iter (for line in curvature)
	    (polyline s line :window 0))
      ;; Derivative
      (title s "Derivative" :window 1)
      (polyline s derivative :window 1 :autozoom t)
      ;; Blend function
      (title s "Blend function" :window 2)
      (polyline s blend :window 2 :autozoom t)
      ;; Blend function derivative
      (title s "Blend function derivative" :window 3)
      (polyline s blend-derivative :window 3 :autozoom t)
      ;; End of file
      (format s "showpage~%"))))

#+nil
(let ((*resolution* 200))
  (write-tomi '(0.3 0.3) '(0.1 0.7) '(0.8 0.3) '(-0.1 0.8) "/tmp/proba2.ps"))

;;; Circle test

#+nil
(let ((margin 0.2)
      (width 0.1)
      (*ribbon-multiplier-start* 1.0d0)
      (*ribbon-multiplier-end* 1.0d0)
      (*resolution* 1000))
  (flet ((sqr (x) (* x x)))
    (let* ((filename "/tmp/circle.ps")
	   (p0 (list margin margin))
	   (p1 (list (- 1 margin) (- 1 margin)))
	   (t0 (list 0 (* 2 (- 1 (* 2 margin)) (1- (sqrt 2.0d0)))))
	   (t1 (list (* -2 (- 1 (* 2 margin)) (1- (sqrt 2.0d0))) 0))
	   (hermite (iter (for i from 0 below *resolution*)
			  (for u = (/ i (1- *resolution*)))
			  (collect (hermite p0 (v* t0 2) p1 (v* t1 2) u))))
	   (hermite-error (iter (for (x y) in hermite)
				(sum (sqr (- (+ (* x x) (* y y)) 1.0d0)))))
	   (tomi (iter (for i from 0 below *resolution*)
		       (for u = (/ i (1- *resolution*)))
		       (collect (tomi p0 t0 p1 t1 u))))
	   (tomi-error (iter (for (x y) in tomi)
			     (sum (sqr (- (+ (* x x) (* y y)) 1.0d0))))))
      (with-open-file (s filename :direction :output :if-exists :supersede)
	(format s "%!PS~%")
	(format s "~f setlinewidth~%" width)
	(init-font s)
	;; Hermite
	(title s "Hermite" :window 2)
	(point s p0 :color *red* :window 2)
	(point s (v+ p0 t0) :color *green* :window 2)
	(polyline s (list p0 (v+ p0 t0)) :color *blue* :window 2)
	(point s p1 :color *red* :window 2)
	(point s (v+ p1 t1) :color *green* :window 2)
	(polyline s (list p1 (v+ p1 t1)) :color *blue* :window 2)
	(arc s (list (- 1 margin) margin) (- 1 (* 2 margin)) 90 180 :color *red* :window 2)
	(polyline s hermite :window 2)
	;; Tomi
	(title s "Ribbon" :window 3)
	(point s p0 :color *red* :window 3)
	(point s (v+ p0 t0) :color *green* :window 3)
	(polyline s (list p0 (v+ p0 t0)) :color *blue* :window 3)
	(point s p1 :color *red* :window 3)
	(point s (v+ p1 t1) :color *green* :window 3)
	(polyline s (list p1 (v+ p1 t1)) :color *blue* :window 3)
	(arc s (list (- 1 margin) margin) (- 1 (* 2 margin)) 90 180 :color *red* :window 3)
	(polyline s tomi :window 3)
	;; Errors
	(title s (format nil "Error: ~f" hermite-error) :window 0)
	(title s (format nil "Error: ~f" tomi-error) :window 1)
	;; End of file
	(format s "showpage~%")))))
