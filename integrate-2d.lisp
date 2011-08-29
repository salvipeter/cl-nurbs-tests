(in-package :cl-nurbs-tests)

(require :cl-utilities)

;;; For second-order differential equations, there are some other possibilities:
;;; - Backward Difference Correction
;;; - Explicit Central Difference Method
;;; - Numerov's [Cowell's] Method
;;; - Implicit Central Difference Method
;;; ... but the Runge-Kutta Method works as well,
;;; see: http://mymathlib.webtrellis.net/diffeq/second_order/

(defun runge-kutta-2d (derivative-matrix step-index step start-index start)
  "Everything is a pair except START.
STEP-INDEX should be (+/-1 +/-1).
The derivative matrices are: f_u(u+h_u/2,v), f_v(u,v+h_v/2).
For derivative matrices f_u(u,v), f_v(u,v), this gives the Euler integral."
  (destructuring-bind (n m) (array-dimensions (first derivative-matrix))
    (let ((result (make-array (list n m))))
      (flet ((runge-kutta-u (i j)
	       (let ((k (- i (first step-index))))
		 (v+ (aref result k j)
		     (v* (aref (first derivative-matrix) k j)
			 (first step)))))
	     (runge-kutta-v (i j)
	       (let ((k (- j (second step-index))))
		 (v+ (aref result i k)
		     (v* (aref (second derivative-matrix) i k)
			 (second step))))))
	(iter (repeat n)
	      (for i first (first start-index) then (+ i (first step-index)))
	      (iter (repeat m)
		    (for j first (second start-index)
			 then (+ j (second step-index)))
		    (setf (aref result i j)
			  (cond ((and (= i (first start-index))
				      (= j (second start-index)))
				 start)
				((= i (first start-index))
				 (runge-kutta-v i j))
				((= j (second start-index))
				 (runge-kutta-u i j))
				(t (affine-combine
				    (runge-kutta-u i j)
				    1/2
				    (runge-kutta-v i j))))))))
      result)))

;;; Test: y^4+x^2+2xy
#+nil
(let ((test-matrices (list (make-array '(50 50))
			   (make-array '(50 50)))))
  (iter (for i from 0 below 50)
	(for x = (/ i 50))
	(for xp = (+ x 1/100))
	(iter (for j from 0 below 50)
	      (for y = (/ j 50))
	      (for yp = (+ y 1/100))
	      (setf (aref (first test-matrices) i j)
		    (list (+ (* 2 xp) (* 2 y))))
	      (setf (aref (second test-matrices) i j)
		    (list (+ (* 4 yp yp yp) (* 2 x))))))
  (with-open-file (s "/tmp/rk-test" :direction :output :if-exists :supersede)
    (let ((matrix (runge-kutta-2d test-matrices
				  '(1 -1) '(0.02d0 -0.02d0) '(0 49) '(1))))
      (iter (for i from 0 below 50)
	    (iter (for j from 0 below 50)
		  (format s "~f ~f ~f~%" (/ i 50) (/ j 50)
			  (first (aref matrix i j))))))))

#+nil
(with-open-file (s "/tmp/rk-test2" :direction :output :if-exists :supersede)
  (iter (for i from 0 below 50)
	(for x = (/ i 50))
	(iter (for j from 0 below 50)
	      (for y = (/ j 50))
	      (format s "~f ~f ~f~%" (/ i 50) (/ j 50)
		      (+ (* y y y y) (* x x) (* 2 x y))))))

(defun numerovs-method (points derivatives step-index step start-index)
  "Everything is a pair except POINTS.
STEP-INDEX should be (+/-1 +/-1).
POINTS is a sparse matrix, where only the two side rows/columns are used.
DERIVATIVES contain the second derivatives by U and V."
  (destructuring-bind (n m) (array-dimensions points)
    (let ((result (make-array (list n m))))
      (flet ((numerov-u (i j)
	       (let* ((i-1 (- i (first step-index)))
		      (i-2 (- i-1 (first step-index))))
		 (v+ (v* (aref result i-1 j) 2)
		     (v* (aref result i-2 j) -1)
		     (v* (v+ (aref (first derivatives) i-2 j)
			     (v* (aref (first derivatives) i-1 j) 10)
			     (aref (first derivatives) i j))
			 1/12 (first step) (first step)))))
	     (numerov-v (i j)
	       (let* ((j-1 (- j (second step-index)))
		      (j-2 (- j-1 (second step-index))))
		 (v+ (v* (aref result i j-1) 2)
		     (v* (aref result i j-2) -1)
		     (v* (v+ (aref (second derivatives) i j-2)
			     (v* (aref (second derivatives) i j-1) 10)
			     (aref (second derivatives) i j))
			 1/12 (second step) (second step))))))
	(iter (for ik from 0 below n)
	      (for i first (first start-index) then (+ i (first step-index)))
	      (iter (for jk from 0 below m)
		    (for j first (second start-index)
			 then (+ j (second step-index)))
		    (setf (aref result i j)
			  (cond ((or (< ik 2) (< jk 2)) (aref points i j))
				((< ik 2) (numerov-v i j))
				((< jk 2) (numerov-u i j))
				(t (affine-combine
				    (numerov-u i j) 1/2 (numerov-v i j))))))))
      result)))

(defun numerov-four-corners (points derivatives step)
  (destructuring-bind (n m) (array-dimensions points)
    (let ((m1 (numerovs-method points derivatives '(1 1) step `(0 0)))
	  (m2 (numerovs-method points derivatives '(-1 1) step `(,(1- n) 0)))
	  (m3 (numerovs-method points derivatives '(1 -1) step `(0 ,(1- m))))
	  (m4 (numerovs-method points derivatives '(-1 -1) step `(,(1- n) ,(1- m)))))
      (list m1 m2 m3 m4))))

;;; Test: y^4+x^2+2xy^2
;;; First derivative: (2x+2y^2, 4y^3+4xy)
;;; Second derivative: (2, 12y^2+4x)
#+nil
(let ((derivatives (list (make-array '(50 50)) (make-array '(50 50))))
      (points (make-array '(50 50))))
  (iter (for i from 0 below 50)
	(for x = (/ i 50))
	(iter (for j from 0 below 50)
	      (for y = (/ j 50))
	      (setf (aref points i j)
		    (list (+ (* y y y y) (* x x) (* 2 x y y))))
	      (setf (aref (first derivatives) i j)
		    (list 2))
	      (setf (aref (second derivatives) i j)
		    (list (+ (* 12 y y) (* 4 x))))))
  (let ((matrices (numerov-four-corners points derivatives '(0.02d0 0.02d0))))
    (flet ((write-matrix (m fname)
	     (with-open-file (s fname :direction :output :if-exists :supersede)
	       (iter (for i from 0 below 50)
		     (iter (for j from 0 below 50)
			   (format s "~f ~f ~f~%" (/ i 50) (/ j 50)
				   (first (aref m i j))))))))
      (write-matrix points "/tmp/numerov0")
      (write-matrix (elt matrices 0) "/tmp/numerov1")
      (write-matrix (elt matrices 1) "/tmp/numerov2")
      (write-matrix (elt matrices 2) "/tmp/numerov3")
      (write-matrix (elt matrices 3) "/tmp/numerov4"))))

;;; Surface test
(let* ((res 100)
       (derivatives (list (make-array `(,res ,res)) (make-array `(,res ,res))))
       (points (make-array `(,res ,res)))
       (surface (first (read-rbn "models/bottom.rbn")))
       (xfrom (first (bss-lower-parameter surface)))
       (xto (first (bss-upper-parameter surface)))
       (xstep (/ (- xto xfrom) (1- res)))
       (yfrom (second (bss-lower-parameter surface)))
       (yto (second (bss-upper-parameter surface)))
       (ystep (/ (- yto yfrom) (1- res))))
  (iter (for i from 0 below res)
	(for x = (interpolate xfrom (/ i (1- res)) xto))
	(iter (for j from 0 below res)
	      (for y = (interpolate yfrom (/ j (1- res)) yto))
	      (setf (aref points i j)
		    (bss-evaluate surface (list x y)))
	      (setf (aref (first derivatives) i j)
		    (bss-evaluate surface (list x y) :derivative '(2 0)))
	      (setf (aref (second derivatives) i j)
		    (bss-evaluate surface (list x y) :derivative '(0 2)))))
  (let ((matrices (numerov-four-corners points derivatives `(,xstep ,ystep))))
    (flet ((write-matrix (m fname)
	     (with-open-file (s fname :direction :output :if-exists :supersede)
	       (format s "~d ~d~%" res res)
	       (iter (for i from 0 below res)
		     (iter (for j from 0 below res)
			   (format s "~{~f~^ ~}~%" (aref m i j)))))))
      (write-matrix points "/tmp/numerov0")
      (write-matrix (elt matrices 0) "/tmp/numerov1")
      (write-matrix (elt matrices 1) "/tmp/numerov2")
      (write-matrix (elt matrices 2) "/tmp/numerov3")
      (write-matrix (elt matrices 3) "/tmp/numerov4"))))

(defun cutback-point (original new max-deviation)
  (let ((d (v- new original)))
    (flet ((sqr (x) (* x x)))
      (v+ original (v* d (- 1.0d0 (sqr (/ (vlength d) max-deviation))))))))

;;; As a Poisson PDE with Dirichlet boundary conditions
(defun integrate-by-pde (points derivatives step iteration
			 &key max-deviation)
  "POINTS contains the initial state (boundaries are preserved).
DERIVATIVES contain the second derivatives by U and V.
STEP is the difference in parameters (assumed to be the same in U and V)."
  (destructuring-bind (n m) (array-dimensions points)
    (let ((result (cl-utilities:copy-array points)))
      (iter (repeat iteration)
	    (iter (for i from 0 below n)
		  (iter (for j from 0 below m)
			(unless (or (= i 0) (= j 0) (= i (1- n)) (= j (1- m)))
			  (setf (aref result i j)
				(v+ (v* (v+ (aref result (1- i) j)
					    (aref result i (1- j))
					    (aref result (1+ i) j)
					    (aref result i (1+ j)))
					1/4)
				    (v* (v+ (aref (first derivatives) i j)
					    (aref (second derivatives) i j))
					-1 step step)))
			  (when max-deviation
			    (setf (aref result i j)
				  (cutback-point (aref points i j)
						 (aref result i j)
						 max-deviation)))))))
      result)))

(defun smooth-mesh (mesh iteration)
  "Warning: modifies mesh."
  (destructuring-bind (n m) (array-dimensions mesh)
    (iter (repeat iteration)
	  (iter (for i from 0 below n)
		(iter (for j from 0 below m)
		      (unless (or (= i 0) (= j 0) (= i (1- n)) (= j (1- m)))
			(setf (aref mesh i j)
			      (v* (v+ (aref mesh (1- i) j)
				      (aref mesh i (1- j))
				      (aref mesh (1+ i) j)
				      (aref mesh i (1+ j)))
				  1/4)))))))
  mesh)

(defun fair-by-pde (surface resolution iteration &optional max-deviation)
  (let* ((derivatives (list (make-array `(,resolution ,resolution))
			    (make-array `(,resolution ,resolution))))
	 (points (make-array `(,resolution ,resolution)))
	 (surface (bss-reparametrize surface 0 1 0 1))
	 (step (/ (1- resolution))))
    (timed-format t "Surface evaluation...~%")
    (iter (for i from 0 below resolution)
	  (for x = (/ i (1- resolution)))
	  (iter (for j from 0 below resolution)
		(for y = (/ j (1- resolution)))
		(setf (aref points i j)
		      (bss-evaluate surface (list x y)))
		(setf (aref (first derivatives) i j)
		      (bss-evaluate surface (list x y) :derivative '(2 0)))
		(setf (aref (second derivatives) i j)
		      (bss-evaluate surface (list x y) :derivative '(0 2)))))
    (timed-format t "Smoothing second derivatives...~%")
    (smooth-mesh (first derivatives) iteration)
    (smooth-mesh (second derivatives) iteration)
    (timed-format t "Gauss-Seidel iteration...~%")
    (let ((result (integrate-by-pde points derivatives step iteration
				    :max-deviation max-deviation)))
      (timed-format t "Surface fitting...~%")
      (prog1 (bss-resembling-fit surface result 0.01d0)
	(timed-format t "Done.~%")))))

;;; Example:
;;; (fair-by-pde (first (read-rbn "models/bottom.rbn")) 100 100 50.0d0)

;;; Fairing test
(let* ((res 100)
       (iteration 100)
       (max-deviation 50.0d0)
       (model "models/bottom.rbn")
       (derivatives (list (make-array `(,res ,res)) (make-array `(,res ,res))))
       (points (make-array `(,res ,res)))
       (surface (bss-reparametrize
		 (first (read-rbn model)) 0 1 0 1))
       (step (/ (1- res))))
  (timed-format t "Surface evaluation...~%")
  (iter (for i from 0 below res)
	(for x = (/ i (1- res)))
	(iter (for j from 0 below res)
	      (for y = (/ j (1- res)))
	      (setf (aref points i j)
		    (bss-evaluate surface (list x y)))
	      (setf (aref (first derivatives) i j)
		    (bss-evaluate surface (list x y) :derivative '(2 0)))
	      (setf (aref (second derivatives) i j)
		    (bss-evaluate surface (list x y) :derivative '(0 2)))))
  (timed-format t "Smoothing second derivatives...~%")
  (smooth-mesh (first derivatives) iteration)
  (smooth-mesh (second derivatives) iteration)
  (timed-format t "Gauss-Seidel iteration...~%")
  (let ((result (integrate-by-pde points derivatives step iteration
				  :max-deviation max-deviation)))
    (timed-format t "Surface fitting...~%")
    (write-rbn (bss-resembling-fit surface points 0.01d0) "/tmp/pde0.rbn")
    (write-rbn (bss-resembling-fit surface result 0.01d0) "/tmp/pde1.rbn"))
  (timed-format t "Done.~%"))