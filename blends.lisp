(in-package :cl-nurbs-tests)

;; (eval-when (:compile-toplevel :load-toplevel :execute)
;;   (asdf:operate 'asdf:load-op 'trivial-shell))

(defparameter *paper* '(500 500))
(defparameter *point-size* 5)
(defparameter *font-size* 10)
(defparameter *resolution* 100)

(defparameter *red* '(1 0 0))
(defparameter *green* '(0 1 0))
(defparameter *blue* '(0 0 1))

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

(defun hermite (p0 t0 p1 t1 u)
  (let ((h00 (+ (* 2 u u u) (* -3 u u) 1))
	(h10 (+ (* u u u) (* -2 u u) u))
	(h01 (+ (* -2 u u u) (* 3 u u)))
	(h11 (- (* u u) (* u u u))))
    (v+ (v* p0 h00) (v* t0 h10) (v* p1 h01) (v* t1 h11))))

(defun hermite-derivative (p0 t0 p1 t1 u)
  (let ((h00 (- (* 6 u u) (* 6 u)))
	(h10 (+ (* 3 u u) (* -4 u) 1))
	(h01 (+ (* -6 u u) (* 6 u)))
	(h11 (- (* 2 u) (* 3 u u))))
    (v+ (v* p0 h00) (v* t0 h10) (v* p1 h01) (v* t1 h11))))

(defun hermite-blend (u)
  (list u (+ (* 2 u u u) (* -3 u u) 1)))

(defun hermite-blend-derivative (u)
  (list u (- (* 6 u u) (* 6 u))))

(defun hermite-blend2 (u)
  (list u (+ (* u u u) (* -2 u u) u)))

(defun hermite-blend2-derivative (u)
  (list u (+ (* 3 u u) (* -4 u) 1)))

(defun write-hermite (p0 t0 p1 t1 filename)
  (let (curve derivative blend blend-derivative blend2 blend2-derivative)
    (iter (for i from 0 below *resolution*)
	  (for u = (/ i (1- *resolution*)))
	  (collect (hermite p0 t0 p1 t1 u)
		   into tcurve)
	  (collect (hermite-derivative p0 t0 p1 t1 u)
		   into tderivative)
	  (collect (hermite-blend u)
		   into tblend)
	  (collect (hermite-blend-derivative u)
		   into tblend-derivative)
	  (collect (hermite-blend2 u)
		   into tblend2)
	  (collect (hermite-blend2-derivative u)
		   into tblend2-derivative)
	  (finally (setf curve tcurve derivative tderivative
			 blend tblend blend-derivative tblend-derivative
			 blend2 tblend2 blend2-derivative tblend2-derivative)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      (init-font s)
      ;; Curve
      (title s "Curve" :window 0)
      (point s p0 :color *red* :window 0)
      (point s (v+ p0 t0) :color *green* :window 0)
      (polyline s (list p0 (v+ p0 t0)) :color *blue* :window 0)
      (point s p1 :color *red* :window 0)
      (point s (v+ p1 t1) :color *green* :window 0)
      (polyline s (list p1 (v+ p1 t1)) :color *blue* :window 0)
      (polyline s curve :window 0)
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

;;; (write-hermite '(0.3 0.3) '(0.1 0.1) '(0.8 0.3) '(-0.1 0.2) "/tmp/proba.ps")

(defun tomi (p0 t0 p1 t1 u &key (r0 1.0) (r1 1.0) (exponent 2))
  (let* ((d0 u)
	 (d1 (- 1 u))
	 (alpha0 (/ (expt d1 exponent) (+ (expt d0 exponent) (expt d1 exponent))))
	 (alpha1 (/ (expt d0 exponent) (+ (expt d0 exponent) (expt d1 exponent)))))
    (v+ (v* (v+ p0 (v* t0 r0 d0)) alpha0)
	(v* (v+ p1 (v* t1 r1 d1)) alpha1))))

(defun tomi-derivative (p0 t0 p1 t1 u &key (r0 1.0) (r1 1.0) (exponent 2))
  (let ((d0 u)
	(d1 (- 1 u))
	(n exponent))
    (v+ (v* (v- (v* (v- (v* t0 r0 (expt d1 n))
			(v* (v+ p0 (v* t0 r0 d0)) n (expt d1 (1- n))))
		    (+ (expt d0 n) (expt d1 n)))
		(v* (v+ p0 (v* t0 r0 d0)) (expt d1 n)
		    (- (* n (expt d0 (1- n))) (* n (expt d1 (1- n))))))
	    (/ (+ (expt d0 (* 2 n))
		  (* (expt d0 n) (expt d1 n))
		  (expt d1 (* 2 n)))))
	(v* (v- (v* (v- (v* (v+ p1 (v* t1 r1 d1)) n (expt d0 (1- n)))
			(v* t1 r1 (expt d0 n)))
		    (+ (expt d0 n) (expt d1 n)))
		(v* (v+ p1 (v* t1 r1 d1)) (expt d0 n)
		    (- (* n (expt d0 (1- n))) (* n (expt d1 (1- n))))))
	    (/ (+ (expt d0 (* 2 n))
		  (* (expt d0 n) (expt d1 n))
		  (expt d1 (* 2 n))))))))

(defun tomi-blend (u &key (exponent 2))
  (let ((n exponent)
	(d0 u)
	(d1 (- 1 u)))
    (list u (/ (expt d1 n) (+ (expt d0 n) (expt d1 n))))))

(defun tomi-blend-derivative (u &key (exponent 2))
  (let ((n exponent)
	(d0 u)
	(d1 (- 1 u)))
    (list u
	  (/ (- (+ (* n (expt d1 (1- n)) (+ (expt d0 n) (expt d1 n)))
		   (* (expt d1 n)
		      (- (* n (expt d0 (1- n))) (* n (expt d1 (1- n)))))))
	     (+ (expt d0 (* 2 n)) (* (expt d0 n) (expt d1 n)) (expt d1 (* 2 n)))))))

(defun write-tomi (p0 t0 p1 t1 filename &key (r0 1.0) (r1 1.0) (exponent 2))
  (let (curve derivative blend blend-derivative)
    (iter (for i from 0 below *resolution*)
	  (for u = (/ i (1- *resolution*)))
	  (collect (tomi p0 t0 p1 t1 u :r0 r0 :r1 r1 :exponent exponent)
		   into tcurve)
	  (collect (tomi-derivative p0 t0 p1 t1 u :r0 r0 :r1 r1 :exponent exponent)
		   into tderivative)
	  (collect (tomi-blend u :exponent exponent)
		   into tblend)
	  (collect (tomi-blend-derivative u :exponent exponent)
		   into tblend-derivative)
	  (finally (setf curve tcurve derivative tderivative
			 blend tblend blend-derivative tblend-derivative)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      (init-font s)
      ;; Curve
      (title s "Curve" :window 0)
      (point s p0 :color *red* :window 0)
      (point s (v+ p0 t0) :color *green* :window 0)
      (polyline s (list p0 (v+ p0 t0)) :color *blue* :window 0)
      (point s p1 :color *red* :window 0)
      (point s (v+ p1 t1) :color *green* :window 0)
      (polyline s (list p1 (v+ p1 t1)) :color *blue* :window 0)
      (polyline s curve :window 0)
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

;;; (write-tomi '(0.3 0.3) '(0.1 0.1) '(0.8 0.3) '(-0.1 0.2) "/tmp/proba2.ps" :r0 0.5 :r1 0.5 :exponent 2)

;;; Circle test

(let ((margin 0.2)
      (width 0.1)
      (*resolution* 1000))
  (flet ((sqr (x) (* x x))
	 (circle-arc (u)
	   (list (+ margin (* (- 1 (* 2 margin)) (- 1 (cos (* u pi 0.25d0)))))
		 (+ margin (* (- 1 (* 2 margin)) (sin (* u pi 0.25d0)))))))
    (let* ((filename "/tmp/circle.ps")
	   (p0 (list margin margin))
	   (p1 (list (- 1 margin) (- 1 margin)))
	   (t0 (list 0 (* 2 (- 1 (* 2 margin)) (1- (sqrt 2.0d0)))))
	   (t1 (list (* -2 (- 1 (* 2 margin)) (1- (sqrt 2.0d0))) 0))
	   (circle (iter (for i from 0 below *resolution*)
			 (for u = (/ i (1- *resolution*)))
			 (collect (circle-arc u))))
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
	(polyline s circle :color *red* :window 2)
	(polyline s hermite :window 2)
	;; Tomi
	(title s "Ribbon" :window 3)
	(point s p0 :color *red* :window 3)
	(point s (v+ p0 t0) :color *green* :window 3)
	(polyline s (list p0 (v+ p0 t0)) :color *blue* :window 3)
	(point s p1 :color *red* :window 3)
	(point s (v+ p1 t1) :color *green* :window 3)
	(polyline s (list p1 (v+ p1 t1)) :color *blue* :window 3)
	(polyline s circle :color *red* :window 3)
	(polyline s tomi :window 3)
	;; Errors
	(title s (format nil "Error: ~f" hermite-error) :window 0)
	(title s (format nil "Error: ~f" tomi-error) :window 1)
	;; End of file
	(format s "showpage~%")))))
