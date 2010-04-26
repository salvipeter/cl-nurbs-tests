;;; -*- mode: lisp; syntax: common-lisp -*-

;;; Usage
;;; -----
;;; First INITIALIZE (preferably with a size parameter >= the # of points).
;;; Then insert the points with INSERT-POINT, which returns indices,
;;; and also insert triangles with SET-TRIANGLE.
;;; Finally call FAIR with the desired parameters.

;;; Note: This implementation assumes that there are no dangling triangles,
;;; ie. each neighborhood has zero [inside] or one [edge] hole.

;;; TODO: The writing functions work only after FAIR (merge-neighbors).
;;; TODO: The weights probably go wrong, when PRESERVE-TANGENTS is NIL.
;;; TODO: Maybe the triangle-area-weight should be moved into FAIR.

(in-package :cl-user)

(defpackage :kobbelt
  (:use :common-lisp :iterate :cl-nurbs)
  (:export :initialize :insert-point :maybe-insert-point :set-triangle :fair
	   :get-point :write-vtk-mesh :write-ply-mesh))

(in-package :kobbelt)

(defstruct point coordinates neighbors border-p derivatives)
(defstruct neighbor index parameters)
(defclass neighbors ()
  ((size :initform 0 :accessor size)
   (points :accessor points)
   (weights :accessor weights)))

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

(defun maybe-insert-point (obj point)
  (with-accessors ((v points) (n size)) obj
    (or (iter (for i from 0 below n)
	      (finding i such-that (equal (point-coordinates (elt v i)) point)))
	(insert-point obj point))))

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
  (labels ((rec (lst item last-needed-p)
	     (if (null lst)
		 (if last-needed-p (list item) nil)
		 (let ((next (find item lst :test #'member)))
		   (cons item (rec (remove next lst :test #'equal)
				   (if (= (first next) item)
				       (second next)
				       (first next))
				   last-needed-p))))))
    (iter (for i from 0 below (size obj))
	  (for neighborhood = (point-neighbors (elt (points obj) i)))
	  (for flattened = (reduce #'append neighborhood))
	  (for edge = (find-if (lambda (x) (= (funcall #'count x flattened) 1))
			       (remove-duplicates flattened)))
	  (setf (point-neighbors (elt (points obj) i))
		(mapcar (lambda (index) (make-neighbor :index index))
			(rec neighborhood (or edge (first flattened)) edge)))
	  (when edge
	    (setf (point-border-p (elt (points obj) i)) 'border)
	    (iter (for neighbor in (point-neighbors (elt (points obj) i)))
		  (for npoint = (elt (points obj) (neighbor-index neighbor)))
		  (unless (point-border-p npoint)
		    (setf (point-border-p npoint) 'inner)))))))

(defun parameterize (obj param-fn)
  "Calls PARAM-FN for all points in OBJ."
  (iter (for i from 0 below (size obj))
	(funcall param-fn obj i)))

(defun neighborhood-triangles (point &optional force-close-p)
  (let ((neighborhood (point-neighbors point)))
    (append (mapcar #'list
		    (butlast neighborhood)
		    (cdr neighborhood))
	    (if (or force-close-p (not (eq (point-border-p point) 'border)))
		(list (append (last neighborhood)
			      (list (first neighborhood))))
		nil))))

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

(defun triangle-area (p1 p2 p3)
  (let* ((a (point-distance p1 p2))
	 (b (point-distance p2 p3))
	 (c (point-distance p1 p3))
	 (s (/ (+ a b c) 2.0d0)))
    (sqrt (* s (- s a) (- s b) (- s c)))))

(defun point-area (obj i)
  (let ((pos (point-coordinates (elt (points obj) i))))
    (reduce #'+ (mapcar (lambda (triangle)
			  (apply #'triangle-area pos
				 (mapcar (lambda (neighbor)
					   (neighbor-coordinates obj neighbor))
					 triangle)))
			(neighborhood-triangles (elt (points obj) i))))))

(defun compute-weights (obj)
  (let ((n (size obj)))
    (setf (weights obj) (make-array n))
    (dotimes (i n) (setf (elt (weights obj) i) (make-hash-table))))
  (iter (for i from 0 below (size obj))
	(for p = (elt (points obj) i))
	(unless (eq (point-border-p p) 'border)
	  (for area = (point-area obj i))
	  (for neighborhood = (point-neighbors p))
	  (for m = (length neighborhood))
	  (for x = (make-array (list m 5)))
	  (iter (for j from 0 below m)
		(for neighbor in neighborhood)
		(for (u v) = (neighbor-parameters neighbor))
		(setf (aref x j 0) u
		      (aref x j 1) v
		      (aref x j 2) (* 1/2 u u)
		      (aref x j 3) (* u v)
		      (aref x j 4) (* 1/2 v v)))
	  (for xt = (matrix:transpose x))
	  (for d = (matrix:multiplication
		    (lu-solver:inverse (matrix:multiplication xt x))
		    xt))
	  (for alpha-sum = (iter (for k from 0 below m) (sum (aref d 2 k))))
	  (for beta-sum = (iter (for k from 0 below m) (sum (aref d 3 k))))
	  (for gamma-sum = (iter (for k from 0 below m) (sum (aref d 4 k))))
	  (flet ((add-weight (i1 i2 w)
		   (unless (gethash i2 (elt (weights obj) i1))
		     (setf (gethash i2 (elt (weights obj) i1)) 0.0d0))
		   (incf (gethash i2 (elt (weights obj) i1)) w)))
	    (add-weight i i (* 2.0d0 area
			       (+ (* alpha-sum alpha-sum)
				  (* 2 beta-sum beta-sum)
				  (* gamma-sum gamma-sum))))
	    (iter (for n1 in neighborhood)
		  (for k1 upfrom 0)
		  (for n1-index = (neighbor-index n1))
		  (for alpha1 = (aref d 2 k1))
		  (for beta1 = (aref d 3 k1))
		  (for gamma1 = (aref d 4 k1)) 
		  (add-weight n1-index i (* -2.0d0 area
					    (+ (* alpha-sum alpha1)
					       (* 2 beta-sum beta1)
					       (* gamma-sum gamma1))))
		  (add-weight i n1-index (* -2.0d0 area
					    (+ (* alpha-sum alpha1)
					       (* 2 beta-sum beta1)
					       (* gamma-sum gamma1))))
		  (iter (for n2 in neighborhood)
			(for k2 upfrom 0)
			(for n2-index = (neighbor-index n2))
			(for weight = (* 2.0d0 area
					 (+ (* alpha1 (aref d 2 k2))
					    (* 2 beta1 (aref d 3 k2))
					    (* gamma1 (aref d 4 k2)))))
			(add-weight n1-index n2-index weight)))))))

;;; Not used at present.
(defstruct derivative u v uu uv vv)
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

;;; Not used at present.
(defun two-neighborhood (obj i)
  (let ((neighborhood (point-neighbors (elt (points obj) i))))
    (remove i (delete-duplicates
	       (append neighborhood
		       (mapcan (lambda (n)
				 (point-neighbors
				  (elt (points obj) (neighbor-index n))))
			       neighborhood)))
	    :key #'neighbor-index)))

(defun fair-more (obj iteration &key (preserve-tangents t))
  (flet ((fixed-point-p (p)
	   (let ((borderp (point-border-p p)))
	     (or (eq borderp 'border)
		 (and preserve-tangents (eq borderp 'inner))))))
    (iter (with n = (size obj))
	  (with points = (points obj))
	  (with weights = (weights obj))
	  (repeat iteration)
	  (dotimes (i n)
	    (unless (fixed-point-p (elt points i))
	      (setf (point-coordinates (elt points i))
		    (iter (with q = '(0 0 0))
			  (for (index weight) in-hashtable (elt weights i))
			  (for point = (point-coordinates (elt points index)))
			  (when (/= index i)
			    (setf q (v+ q (v* point weight))))
			  (finally
			   (let ((weight (gethash i (elt weights i))))
			     (return (v* q (/ (- weight)))))))))))))

(defun fair (obj iteration &key
	     (parameterization :projection) (preserve-tangents t))
  (merge-neighbors obj)
  (parameterize obj (ecase parameterization
		      (:projection #'parameterize-projection)
		      (:polar #'parameterize-polar)))
  (compute-weights obj)
  (fair-more obj iteration :preserve-tangents preserve-tangents))

(defun generate-triangle-list (obj)
  (iter (for i from 0 below (size obj))
	(appending (mapcan (lambda (tr)
			     (let ((indices (mapcar #'neighbor-index tr)))
			       (when (and (< i (first indices))
					  (< i (second indices)))
				 (list (cons i indices)))))
			   (neighborhood-triangles (elt (points obj) i))))))

(defun get-point (obj index)
  (point-coordinates (elt (points obj) index)))
(declaim (inline get-point))

(defun write-vtk-mesh (obj filename)
  (let ((triangles (generate-triangle-list obj))
	(points (points obj)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "# vtk DataFile Version 1.0~%~
                 Exported Mesh~%~
                 ASCII~%~
                 DATASET POLYDATA~%~
                 POINTS ~d float~%"
	      (size obj))
      (dotimes (i (size obj))
	(format s "~{~f~^ ~}~%" (point-coordinates (elt points i))))
      (let ((n (length triangles)))
	(format s "~%POLYGONS ~d ~d~%" n (* 4 n)))
      (dolist (p triangles) (format s "3~{ ~d~}~%" p)))))

(defun write-ply-mesh (obj filename)
  (let ((triangles (generate-triangle-list obj))
	(points (points obj)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "ply~%~
                 format ascii 1.0~%~
                 comment Exported Mesh~%~
                 element vertex ~d~%~
                 property float x~%~
                 property float y~%~
                 property float z~%~
                 element face ~d~%~
                 property list uchar int vertex_index~%~
                 end_header~%"
	      (size obj) (length triangles))
      (dotimes (i (size obj))
	(format s "~{~f~^ ~}~%" (point-coordinates (elt points i))))
      (dolist (p triangles) (format s "3~{ ~d~}~%" p)))))
