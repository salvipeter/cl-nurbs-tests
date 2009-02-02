;;; -*- mode: lisp; syntax: common-lisp -*-

;;; Usage
;;; -----
;;; First INITIALIZE (preferably with a size parameter >= the # of points).
;;; Then insert the points with INSERT-POINT, which returns indices,
;;; and also insert triangles with SET-TRIANGLE.
;;; Finally call FAIR with the desired parameters.

(in-package :cl-user)

(defpackage :kobbelt
  (:use :common-lisp :iterate :cl-nurbs)
  (:export :initialize :insert-point :set-triangle :fair))

(in-package :kobbelt)

(defstruct point coordinates neighbors closedp derivatives)
(defstruct neighbor index parameters)
(defstruct derivative u v uu uv vv)
(defclass neighbors ()
  ((size :initform 0 :accessor size)
   (points :accessor points)))

(defun neighbor-coordinates (obj neighbor)
  (point-coordinates (elt (points obj) (neighbor-index neighbor))))

(defun initialize (&optional (size 100))
  (let ((result (make-instance 'neighbors)))
    (setf (points result) (make-array size :adjustable t))
    result))

(defun insert-point (obj point)
  (with-accessors ((v points) (n size)) obj
    (when (>= n (length v))
      (adjust-array v (* 2 (length v))))
    (setf (elt v n) (make-point :coordinates point))
    (prog1 n (incf n))))

;;; Temporarily just store tuples, and merge them before parameterization.
(defun set-triangle (obj i j k)
  (with-accessors ((v points)) obj
    (flet ((set-one (i j k) (push (list j k) (point-neighbors (elt v i)))))
      (set-one i j k)
      (set-one j i k)
      (set-one k i j))
    t))

(defun merge-neighbors (obj)
  "Destructively merges the neighborhood triangles into one list per point."
  (labels ((rec (lst item)
	     (if (null lst)
		 (list item)
		 (let ((next (find item lst :test #'member)))
		   (cons item (rec (remove next lst :test #'equal)
				   (if (= (first next) item)
				       (second next)
				       (first next))))))))
    (iter (for i from 0 below (size obj))
	  (for neighborhood = (point-neighbors (elt (points obj) i)))
	  (for flattened = (reduce #'append neighborhood))
	  (for edge = (find-if (lambda (x) (= (funcall #'count x flattened) 1))
			       (remove-duplicates flattened)))
	  (setf (point-neighbors (elt (points obj) i))
		(mapcar (lambda (index) (make-neighbor :index index))
			(rec neighborhood (or edge (first flattened)))))
	  (setf (point-closedp (elt (points obj) i)) (null edge)))))

(defun parameterize (obj param-fn)
  "Calls PARAM-FN for all points in OBJ."
  (iter (for i from 0 below (size obj))
	(funcall param-fn obj i)))

(defun neighborhood-triangles (point &optional force-close-p)
  (let ((neighborhood (point-neighbors point)))
    (append (if (or force-close-p (point-closedp point))
		(list (append (last neighborhood)
			      (list (first neighborhood))))
		nil)
	    (mapcar #'list
		    (butlast neighborhood)
		    (cdr neighborhood)))))

(defun approximate-normal (obj i)
  (let* ((p (elt (points obj) i))
	 (pos (point-coordinates p)))
    (flet ((approximate (triangle)
	       (let ((q1 (neighbor-coordinates obj (first triangle)))
		     (q2 (neighbor-coordinates obj (second triangle))))
		 (vnormalize (cross-product (v- q1 pos) (v- q2 pos))))))
      (let ((normals (mapcar #'approximate (neighborhood-triangles p))))
	(vnormalize (reduce #'v+ normals))))))

(defun parameterize-projection (obj i)
  "Destructively sets the parameters of I's neighbors."
  (let* ((p (elt (points obj) i))
	 (pos (point-coordinates p))
	 (n (approximate-normal obj i))
	 (projected (mapcar (lambda (x)
			      (let ((q (neighbor-coordinates obj x)))
				(v- q (v* n (scalar-product (v- q pos) n)))))
			    (point-neighbors p)))
	 (u0 (vnormalize (v- (first projected) pos)))
	 (v0 (vnormalize (cross-product u0 n))))
    (iter (for neighbor in (point-neighbors p))
	  (for coordinates in projected)
	  (for deviation = (v- coordinates pos))
	  (setf (neighbor-parameters neighbor)
		(list (scalar-product deviation u0)
		      (scalar-product deviation v0))))))

(defun parameterize-polar (obj i)
  "Destructively sets the parameters of I's neighbors."
  (let* ((p (elt (points obj) i))
	 (pos (point-coordinates p)))
    (flet ((angle (triangle)
	     (let ((n1 (neighbor-coordinates obj (first triangle)))
		   (n2 (neighbor-coordinates obj (second triangle))))
	       (acos (scalar-product (vnormalize (v- n1 pos))
				     (vnormalize (v- n2 pos))))))
	   (deviation (x)
	     (point-distance (neighbor-coordinates obj x) pos)))
      (let* ((angles (mapcar #'angle (neighborhood-triangles p t)))
	     (weight (/ (* 2 pi) (reduce #'+ angles)))
	     (deviations (mapcar #'deviation (point-neighbors p)))
	     (max-deviation (reduce #'max deviations)))
	(iter (for neighbor in (point-neighbors p))
	      (for deviation in deviations)
	      (for alpha first 0 then (+ alpha angle))
	      (for angle in angles)
	      (setf (neighbor-parameters neighbor)
		    (list (* (/ deviation max-deviation)
			     (cos (* alpha weight)))
			  (* (/ deviation max-deviation)
			     (sin (* alpha weight))))))))))

(defun derivatives (obj)
  "Destructively sets the derivatives for every point."
  (flet ((extract-coordinates (v i)
	   (list (elt v i) (elt v (1+ i)) (elt v (+ i 2)))))
    (iter (for i from 0 below (size obj))
	  (for p = (elt (points obj) i))
	  (for pos = (point-coordinates p))
	  (for neighborhood = (point-neighbors p))
	  (for n = (length neighborhood))
	  (for x = (make-array (list (* 3 n) 15) :initial-element 0.0d0))
	  (for y = (make-array (list (* 3 n) 1)))
	  (iter (for j from 0 below n)
		(for neighbor in neighborhood)
		(for (u v) = (neighbor-parameters neighbor))
		(for coefficients = (list u v (* 1/2 u u) (* u v) (* 1/2 v v)))
		(for deviation = (v- (neighbor-coordinates obj neighbor) pos))
		(dotimes (k 3)
		  (iter (for l from 0 below 5)
			(for coefficient in coefficients)
			(setf (aref x (+ (* 3 j) k) (+ (* 3 l) k)) coefficient))
		  (setf (aref y (+ (* 3 j) k) 0) (elt deviation k))))
	  (for solution =
	       (if (>= n 5)
		   (lu-solver:least-squares x y)
		   (lu-solver:least-norm x y)))
	  (setf (point-derivatives p)
		(make-derivative :u  (extract-coordinates solution 0)
				 :v  (extract-coordinates solution 3)
				 :uu (extract-coordinates solution 6)
				 :uv (extract-coordinates solution 9)
				 :vv (extract-coordinates solution 12))))))

(defun triangle-area (p1 p2 p3)
  (let* ((a (point-distance p1 p2))
	 (b (point-distance p2 p3))
	 (c (point-distance p1 p3))
	 (s (/ (+ a b c) 2.0d0)))
    (sqrt (* s (- s a) (- s b) (- s c)))))

(defun fair (obj &key (parameterization 'projection))
  (merge-neighbors obj)
  (parameterize obj (ecase parameterization
		      (projection #'parameterize-projection)
		      (polar #'parameterize-polar)))
  'TODO)
