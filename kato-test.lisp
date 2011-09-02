(in-package :cl-nurbs-tests)

(defparameter *width* 640)
(defparameter *height* 480)
(defparameter *line-width* 2.0d0)
(defparameter *point-size* 3.0d0)
(defparameter *density* 50.0d0)
(defparameter *tolerance* 1.0d0)
(defparameter *epsilon* 1.0d-5)
(defparameter *colors*
  '((on-line . (0 0 0))
    (outside . (255 255 255))
    (nothing . (255 255 255))
    (si . (255 0 0))
    (di . (0 200 0))
    (si-1 . (0 0 255))
    (di-1 . (255 0 255))))

(defun kato-test (p1 p2 filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%255~%" *width* *height*)
    (iter (with length = (point-distance p1 p2))
	  (for y from 0 below *height*)
	  (iter (for x from 0 below *width*)
		(unless (zerop x)
		  (princ #\Space s))
		(for p = (list x y))
		(for p1d = (point-distance p p1))
		(for p2d = (point-distance p p2))
		(for distance = (+ p1d p2d (- length)))
		(cond ((and (< (abs distance) *line-width*)
			    (< p1d length) (< p2d length))
		       (princ (floor (* (abs distance) 255) *line-width*) s))
		      ((< (mod distance *density*) *line-width*)
		       (princ 0 s))
		      (t (princ 255 s))))
	  (terpri s))))

(defun on-segment-p (q1 q2 x y)
  "Test whether the point (X Y) is on the segment defined by Q1 and Q2."
  (let ((p (list x y))
	(dir (vnormalize (v- q2 q1))))
    (and (< (point-distance p (v+ q1 (v* dir (scalar-product (v- p q1) dir))))
	    *tolerance*)
	 (< (max (point-distance p q1)
		 (point-distance p q2))
	    (point-distance q1 q2)))))

(defun kato-test-2 (p0 p1 p2 p3 filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%2~%" *width* *height*)
    (let ((length1 (point-distance p0 p1))
	  (length2 (point-distance p2 p3)))
      (labels ((parameter (x y)
		 (let* ((p (list x y))
			(p0d (point-distance p p0))
			(p1d (point-distance p p1))
			(p2d (point-distance p p2))
			(p3d (point-distance p p3))
			(distance1 (+ p0d p1d (- length1)))
			(distance2 (+ p2d p3d (- length2))))
		   (/ distance1 (+ distance1 distance2))))
	       (better-points (x y base)
		 (iter (for i from -1 to 1)
		       (sum (iter (for j from -1 to 1)
				  (for current =
				       (mod (parameter (+ x i) (+ y j))
					    *density*))
				  (count (< current base)))))))
	(iter (for y from 0 below *height*)
	      (iter (for x from 0 below *width*)
		    (unless (zerop x)
		      (princ #\Space s))
		    (for current = (mod (parameter x y) *density*))
		    (cond ((or (on-segment-p p0 p1 x y)
			       (on-segment-p p1 p2 x y)
			       (on-segment-p p2 p3 x y))
			   (princ 0 s))
			  ((and (< current *line-width*)
				(<= (better-points x y current) 2))
			   (princ 1 s))
			  (t (princ 2 s))))
	      (terpri s))))))

(defun center-line (w0 w2 length center)
  "Both W0, W2 and CENTER are in local coordinates."
  (let* ((s (v* (v+ w0 w2) 1/2))
	 (alpha (/ (second center) (second s))))
    (list (- (/ (- (* 2 (first center)) length) alpha) (first s)) (second s))))

(defun line-sweep (w0 w1 w2 length p)
  (let ((a (+ (* length (- (second w2) (second w0)))
	      (* (second p) (+ (first w2) (* -2 (first w1)) (first w0)))))
	(b (+ (* length (second w0))
	      (* (first p) (- (second w0) (second w2)))
	      (* 2 (second p) (- (first w1) (first w0)))))
	(c (- (* (second p) (first w0)) (* (first p) (second w0)))))
    (second-degree-solver a b c :min 0.0d0 :max 1.0d0)))

(defun central-line-sweep (p0 p1 p2 p3 center filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%2~%" *width* *height*)
    (let* ((base-x (vnormalize (v- p2 p1)))
	   (base-y (list (- (second base-x)) (first base-x)))
	   (length (point-distance p1 p2))
	   (w0 (in-system base-x base-y (v- p0 p1)))
	   (w2 (in-system base-x base-y (v- p3 p2)))
	   (w1 (center-line w0 w2 length
			    (in-system base-x base-y (v- center p1)))))
      (labels ((parameter (x y)
		 (line-sweep w0 w1 w2 length
			     (in-system base-x base-y (v- (list x y) p1))))
	       (better-points (x y base)
		 (iter (for i from -1 to 1)
		       (sum (iter (for j from -1 to 1)
				  (for current =
				       (mod (parameter (+ x i) (+ y j))
					    *density*))
				  (count (< current base)))))))
	(iter (for y from 0 below *height*)
	      (iter (for x from 0 below *width*)
		    (unless (zerop x)
		      (princ #\Space s))
		    (for current = (mod (parameter x y) *density*))
		    (cond ((or (on-segment-p p0 p1 x y)
			       (on-segment-p p1 p2 x y)
			       (on-segment-p p2 p3 x y)
			       (< (point-distance (list x y) center)
				  *point-size*))
			   (princ 0 s)) 
			  ((and (< current *line-width*)
				(<= (better-points x y current) 2))
			   (princ 1 s))
			  (t (princ 2 s))))
	      (terpri s))))))

(defun radial-test (p0 p1 p2 p3 filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%2~%" *width* *height*)
    (let* ((base-x (vnormalize (v- p2 p1)))
	   (base-y (list (- (second base-x)) (first base-x)))
	   (a (in-system base-x base-y (v- p0 p1)))
	   (b (in-system base-x base-y (v- p3 p2)))
	   (length (point-distance p1 p2))
	   (axby-aybx (- (* (first a) (second b)) (* (second a) (first b))))
	   (laxby-div (/ (* length (first a) (second b)) axby-aybx))
	   (layby-div (/ (* length (second a) (second b)) axby-aybx)))
      (labels ((parameter (x y)
		 (destructuring-bind (x y)
		     (in-system base-x base-y (v- (list x y) p1))
		   (/ (+ x (* (/ y (- y layby-div)) (- laxby-div x))) length)))
	       (better-points (x y base)
		 (iter (for i from -1 to 1)
		       (sum (iter (for j from -1 to 1)
				  (for current =
				       (mod (parameter (+ x i) (+ y j))
					    *density*))
				  (count (< current base)))))))
	(iter (for y from 0 below *height*)
	      (iter (for x from 0 below *width*)
		    (unless (zerop x)
		      (princ #\Space s))
		    (for current = (mod (parameter x y) *density*))
		    (cond ((or (on-segment-p p0 p1 x y)
			       (on-segment-p p1 p2 x y)
			       (on-segment-p p2 p3 x y))
			   (princ 0 s))
			  ((and (< current *line-width*)
				(<= (better-points x y current) 2))
			   (princ 1 s))
			  (t (princ 2 s))))
	      (terpri s))))))

(defgeneric compute-distance (type points segments p dir))

(defmethod compute-distance ((type (eql 'perpendicular)) points segments p dir)
  (let ((p0 (elt segments 1))
	(p1 (elt segments 2)))
    (if (eq dir 'd)
	(/ (point-line-distance p (list p0 p1))
	   (iter (for q1 in points)
		 (maximize (iter (for q2 in points)
				 (maximize (point-distance q1 q2))))))
	(let ((di-1 (point-line-distance p (list (elt segments 0) p0)))
	      (di+1 (point-line-distance p (list p1 (elt segments 3)))))
	  (/ di-1 (+ di-1 di+1))))))

(defmethod compute-distance ((type (eql 'barycentric)) points segments p dir)
  (let ((lines (lines-from-points points)))
    (flet ((area (segment)
	     (let* ((a (point-distance p (first segment)))
		    (b (point-distance p (second segment)))
		    (c (apply #'point-distance segment))
		    (s (/ (+ a b c) 2)))
	       (sqrt (max 0 (* s (- s a) (- s b) (- s c)))))))
      (if (eq dir 'd)
	  (* (/ (area (subseq segments 1 3))
		(reduce #'+ (mapcar #'area lines)))
	     (/ (length lines) 2))
	  (let ((di-1 (area (subseq segments 0 2)))
		(di+1 (area (subseq segments 2))))
	    (/ di-1 (+ di-1 di+1)))))))

(defmethod compute-distance ((type (eql 'chord-based)) points segments p dir)
  (if (eq dir 'd)
      (let* ((length (point-distance (second segments) (third segments)))
	     (p1d (point-distance p (second segments)))
	     (p2d (point-distance p (third segments))))
	(/ (+ p1d p2d (- length))
	   (iter (for q1 in points)
		 (maximize (iter (for q2 in points)
				 (maximize (point-distance q1 q2)))))))
      (let* ((length1 (point-distance (first segments) (second segments)))
	     (length2 (point-distance (third segments) (fourth segments)))
	     (p0d (point-distance p (first segments)))
	     (p1d (point-distance p (second segments)))
	     (p2d (point-distance p (third segments)))
	     (p3d (point-distance p (fourth segments)))
	     (distance1 (+ p0d p1d (- length1)))
	     (distance2 (+ p2d p3d (- length2))))
	(/ distance1 (+ distance1 distance2)))))

(defun radial-intersection (segments p)
  (let* ((p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments))
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (- (second base-x)) (first base-x)))
	 (a (in-system base-x base-y (v- p0 p1)))
	 (b (in-system base-x base-y (v- p3 p2)))
	 (length (point-distance p1 p2))
	 (axby-aybx (- (* (first a) (second b)) (* (second a) (first b)))))
    (if (< (abs axby-aybx) *epsilon*)
	(let ((slength (scalar-product (v- p p1) base-x)))
	  (values (/ slength length) (v+ p1 (v* base-x slength))))
	(let ((laxby-div (/ (* length (first a) (second b)) axby-aybx))
	      (layby-div (/ (* length (second a) (second b)) axby-aybx)))
	  (destructuring-bind (x y)
	      (in-system base-x base-y (v- p p1))
	    (let ((x0 (if (< (abs (- y layby-div)) *epsilon*)
			  0.5 ; kutykurutty the point is the focus itself
			  (+ x (* (/ y (- y layby-div)) (- laxby-div x))))))
	      (values (/ x0 length) (v+ p1 (v* base-x x0)))))))))

(defmethod compute-distance ((type (eql 'radial)) points segments p dir)
  (multiple-value-bind (s q) (radial-intersection segments p)
    (if (eq dir 'd)
	(/ (point-distance p q)
	   (iter (for q1 in points)
		 (maximize (iter (for q2 in points)
				 (maximize (point-distance q1 q2))))))
	s)))

(defun segments-prev (points segments)
  (let ((i (position (first segments) points :test #'equal)))
    (cons (elt points (mod (1- i) (length points)))
	  (butlast segments))))

(defun segments-next (points segments)
  (let ((i (position (first segments) points :test #'equal)))
    (append (rest segments)
	    (list (elt points (mod (+ i 4) (length points)))))))

(defconstant +hermite-blend+ (lambda (x) (hermite-blend-function 'point 'start x)))
(defconstant +distance-blend+ (lambda (x) (blend (list x (- 1.0d0 x)) 0)))

(defmacro defmodified-distance (distance blend-fn)
  "Warning: parameters are evaluated multiple times."
  `(defmethod compute-distance ((type (eql ',(intern (format nil "~:@(~a~)-MOD" distance))))
				points segments p dir)
     (let ((si (compute-distance ',distance points segments p 's)))
       (if (eq dir 's)
	   si
	   (let ((si-1 (compute-distance ',distance points
					 (segments-prev points segments)
					 p 's)) 
		 (si+1 (compute-distance ',distance points
					 (segments-next points segments)
					 p 's)))
	     (+ (* (- 1.0d0 si-1) (funcall ,blend-fn si))
		(* si+1 (funcall ,blend-fn (- 1.0d0 si)))))))))

(defmodified-distance radial +distance-blend+)
(defmodified-distance perpendicular +distance-blend+)
(defmodified-distance barycentric +distance-blend+)
(defmodified-distance chord-based +distance-blend+)
(defmodified-distance line-sweep +distance-blend+)

(defun blend3 (x a b c)
  "Blends between three values, giving A for 0, B for 0.5, C for 1."
  (+ (* x x (+ (* 2 a) (* -4 b) (* 2 c)))
     (* x (+ (* -3 a) (* 4 b) (- c)))
     a))

(defun second-degree-solver (a b c &key min max)
  "TODO: this starts to be incomprehensible..."
  (if (< (abs a) *epsilon*)
      (- (/ c b))
      (let ((discr (- (* b b) (* 4 a c))))
	(cond ((< (abs discr) *epsilon*)
	       (/ (- b) 2 a))
	      ((< discr 0) 0)
	      (t (let* ((discr (sqrt discr))
			(u1 (/ (+ (- b) discr) (* 2 a)))
			(u2 (/ (- (- b) discr) (* 2 a)))
			(u1m (if min (max u1 min) u1))
			(u1mm (if max (min u1m max) u1m))
			(u2m (if min (max u2 min) u2))
			(u2mm (if max (min u2m max) u2m)))
		   (cond ((and min max)
			  (if (< (abs (- u1 u1mm)) (abs (- u2 u2mm)))
			      (values u1mm u2mm)
			      (values u2mm u1mm)))
			 (min
			  (cond ((or (< (abs (- u1 min)) *epsilon*)
				     (< (abs (- u2 min)) *epsilon*))
				 (values min (max u1mm u2mm)))
				((< (abs (- u1 u1mm)) (abs (- u2 u2mm)))
				 (values u1mm u2mm))
				(t (values u2mm u1mm))))
			 (max
			  (cond ((or (< (abs (- u1 max)) *epsilon*)
				     (< (abs (- u2 max)) *epsilon*))
				 (values max (min u1mm u2mm)))
				((< (abs (- u1 u1mm)) (abs (- u2 u2mm)))
				 (values u1mm u2mm))
				(t (values u2mm u1mm))))
			 (t (values u1 u2)))))))))

(defparameter *centralized-line-sweep* nil)
(defmethod compute-distance ((type (eql 'line-sweep)) points segments p dir)
  "BUG: does not work if the base segment is parallel to the Y axis."
  (let* ((c (central-point points (lines-from-points points) t))
	 (p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments)) 
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (- (second base-x)) (first base-x)))
	 (center (in-system base-x base-y (v- c p1)))
	 (p (in-system base-x base-y (v- p p1)))
	 (length (point-distance p1 p2))
	 (w0 (in-system base-x base-y (v- p0 p1)))
	 (w2 (in-system base-x base-y (v- p3 p2)))
	 (w1 (center-line w0 w2 length center))
	 (s (line-sweep w0 w1 w2 length p)))
    (if (eq dir 'd)
	(case *centralized-line-sweep*
	  ((t)
	   (/ (point-distance p (list (* s length) 0))
	      (blend3 s (vlength w0)
		      (* 2 (point-distance center (list (/ length 2) 0)))
		      (vlength w2))))
	  ((nil)
	   (/ (point-distance p (list (* s length) 0))
	      (iter (for q1 in points)
		    (maximize (iter (for q2 in points)
				    (maximize (point-distance q1 q2)))))))
	  (otherwise
	   (let* ((a (- (* (second w2) (first w0)) (* (second w0) (first w2))))
		  (b (- (+ (* (second p) (first w2)) (* (first p) (second w0)))
			(+ (* (first p) (second w2)) (* (second p) (first w0))
			   (* length (second w0)))))
		  (l (second-degree-solver a b (* (second p) length) :min 0.0d0)))
	     (* l (blend3 s 1 (/ *centralized-line-sweep*) 1)))))
	s)))

(defun central-point (points lines weightedp)
  (flet ((slength (line) (apply #'point-distance line)))
    (if weightedp
	(v* (reduce #'v+ (mapcar (lambda (p l1 l2)
				   (v* p (+ (slength l1) (slength l2))))
				 points lines (append (rest lines) lines)))
	    (/ (* 2 (reduce #'+ (mapcar #'slength lines)))))
	(v* (reduce #'v+ points) (/ (length lines))))))

(defun write-color (stream type &optional (1-alpha 0))
  (let ((alpha (- 1.0d0 1-alpha)))
    (flet ((aliasing (x) (floor (+ (* alpha x) (* 1-alpha 255)))))
      (format stream "~{~d~^ ~}" (mapcar #'aliasing (cdr (assoc type *colors*)))))))

(defun distance-function-test (points distance-type line-type filename)
  (let* ((points (mapcar (lambda (p)
			   (list (/ (+ *width* (* *width* (first p))) 2)
				 (/ (+ *height* (* *height* (second p))) 2)))
			 points))
	 (lines (lines-from-points points))
	 (center (central-point points lines t)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "P3~%~d ~d~%255~%" *width* *height*)
      (labels ((parameter (lst type x y)
		 (compute-distance distance-type points lst (list x y) type))
	       (outsidep (x y)
		 (some (lambda (line)
			 (< (* (point-line-distance (list x y) line t)
			       (point-line-distance center line t))
			     0))
			lines))
	       (better-points (lst type x y base)
		 (iter (for i from -1 to 1)
		       (sum (iter (for j from -1 to 1)
				  (for current =
				       (mod (parameter lst type (+ x i) (+ y j))
					    *density*))
				  (count (< current base))))))
	       (on-parameter-line (lst type x y)
		 (let ((current (mod (parameter lst type x y) *density*)))
		   (and (< current *line-width*)
			(<= (better-points lst type x y current) 2)
			(/ current *line-width*)))))
	(iter (for y from 0 below *height*)
	      (iter (for x from 0 below *width*)
		    (unless (zerop x)
		      (princ #\Space s))
		    (acond ((or (some (lambda (l)
					(and (on-segment-p (first l) (second l) x y)
					     (/ (point-line-distance (list x y) l)
						*tolerance*)))
				      lines)
				(let ((center-error (point-distance (list x y) center)))
				  (and (< center-error *point-size*)
				       (/ center-error *point-size*))))
			    (write-color s 'on-line it))
			   ((outsidep x y)
			    (declare (ignore it))
			    (write-color s 'outside))
			   ((on-parameter-line (subseq points 0 4)
					       (first line-type) x y) 
			    (if (eq (first line-type) 's)
				(write-color s 'si it)
				(write-color s 'di it)))
			   ((eq (first line-type) (second line-type))
			    (if (on-parameter-line (subseq points 1 5)
						   (second line-type) x y)
				(if (eq (second line-type) 's)
				    (write-color s 'si-1 it)
				    (write-color s 'di-1 it))
				(write-color s 'nothing)))
			   ((on-parameter-line (subseq points 0 4)
					       (second line-type) x y)
			    (if (eq (second line-type) 's)
				(write-color s 'si it)
				(write-color s 'di it)))
			   (t
			    (declare (ignore it))
			    (write-color s 'nothing))))
	      (terpri s))))))

(defun distance-function-complete-set (points directory)
  (iter (for type in '(perpendicular barycentric chord-based radial line-sweep))
	(iter (for lines in '((s d) (s s) (d d)))
	      (for name = (format nil "~(~a~)-~{~(~a~)~}.ppm" type lines))
	      (for filename = (make-pathname :directory directory :name name))
	      (distance-function-test points type lines filename))))

;; (kato-test '(250 240) '(390 240) "/tmp/kato.pgm")

#+nil
(let ((*density* 5.0d0)
      (*line-width* 1.0d0))
  (kato-test '(10 240) '(630 240) "/tmp/kato.pgm"))

#+nil
(let ((*density* 0.1d0)
      (*line-width* 0.01d0))
  (kato-test-2 '(180 120) '(250 240) '(390 240) '(460 120) "/tmp/kato.pgm"))

#+nil
(let ((*density* 0.1d0)
      (*line-width* 0.01d0))
  (kato-test-2 '(160 100) '(220 240) '(420 240) '(460 270) "/tmp/kato.pgm"))

#+nil
(let ((*density* 0.03d0)
      (*line-width* 0.01d0))
  (kato-test-2 '(190 100) '(220 240) '(420 240) '(500 180) "/tmp/kato.pgm"))

#+nil
(let ((*density* 0.025d0)
      (*line-width* 0.01d0))
  (central-line-sweep '(190 100) '(220 240) '(420 240) '(500 180) '(380 80)
		      "/tmp/sweep.pgm"))

#+nil
(let ((*density* 0.03d0)
      (*line-width* 0.01d0))
  (radial-test '(190 100) '(220 240) '(420 240) '(500 180) "/tmp/radial.pgm"))

#+nil
(let ((*width* 640) (*height* 480) 
      (*density* 0.05d0)
      (*line-width* 0.01d0)
      (*point-size* 3.0d0)
      (*tolerance* 1.0d0))
  (distance-function-test (points-from-angles '(40 20 60 100 80)) 'chord-based '(s d)
			  "/tmp/distance.ppm"))

#+nil
(defun trace-parameter-buggy (points i distance-type type parameter resolution)
  (let ((lines (lines-from-points points)))
    (labels ((query (p) (- (elt (compute-parameter distance-type type points p t) i) parameter))
	     (between (p1 p2)
	       (when p2
		 (affine-combine (rest p1) (/ (first p2) (- (first p2) (first p1))) (rest p2))))
	     (get-start (line)
	       (let ((k (/ (point-distance (first line) (second line)) resolution)))
		 (iter (repeat (1+ (floor k)))
		       (with d = (v* (v- (second line) (first line)) (/ k)))
		       (for p first (first line) then (v+ p d))
		       (collect (cons (abs (query p)) p) into tests)
		       (finally (return (rest (first (sort tests #'< :key #'first))))))))
	     (best-few (actual value)
	       (iter (for xi from -1 to 1)
		     (for x = (+ (first actual) (* xi resolution)))
		     (appending
		      (iter (for yi from -1 to 1)
			    (unless (= xi yi 0)
			      (for y = (+ (second actual) (* yi resolution)))
			      (collect (list (query (list x y)) x y))))
		      into lst)
		     (finally (return (sort (mapcar (lambda (x)
						      (cons (abs (first x)) (rest x)))
						    (remove-if (lambda (x)
								 (> (* (first x) value) 0))
							       lst))
					    #'< :key #'first)))))
	     (best-one (prev actual)
	       (let ((value (query actual)))
		 (between
		  (cons value actual)
		  (first (remove-if (lambda (point)
				      (let ((p (rest point)))
					(or (not (insidep lines p))
					    (and prev
						 (> (scalar-product (v- prev actual) (v- p actual))
						    0)))))
				    (best-few actual value)))))))
      (let ((start (get-start (elt lines (case type
					   (s i)
					   (d1 (mod (1- i) (length lines)))
					   (d2 (mod (1+ i) (length lines))))))))
	(cons start
	      (iter (for pprev first nil then prev)
		    (for prev first start then next)
		    (for next = (best-one pprev prev))
		    (while next)
		    (collect next)))))))

(defun trace-parameter (points i distance-type type parameter resolution)
  (let ((lines (lines-from-points points)))
    (labels ((query (p)
	       (abs (- (elt (compute-parameter distance-type type points p t) i) parameter)))
	     ;; (between (p1 v1 p2 v2) (affine-combine p1 (/ v2 (- v2 v1)) p2))
	     (get-start (line)
	       (let ((k (/ (point-distance (first line) (second line)) resolution)))
		 (iter (repeat (1+ (floor k)))
		       (with d = (v* (v- (second line) (first line)) (/ k)))
		       (for p first (first line) then (v+ p d))
		       (collect (cons (query p) p) into tests)
		       (finally (return (first (sort tests #'< :key #'first)))))))
	     (get-fuzzy-start (index)
	       (iter (for i from 0 below (length lines))
		     (for line in lines)
		     (unless (= i index)
		       (for (value . param) = (get-start line))
		       (finding param minimizing value))))
	     (best-three (actual)
	       (iter (for xi from -1 to 1)
		     (for x = (+ (first actual) (* xi resolution)))
		     (appending
		      (iter (for yi from -1 to 1)
			    (unless (= xi yi 0)
			      (for y = (+ (second actual) (* yi resolution)))
			      (collect (list (query (list x y)) x y))))
		      into lst)
		     (finally (return (subseq (mapcar #'rest (sort lst #'< :key #'first)) 0 3)))))
	     (best-one (prev actual)
	       (first (remove-if (lambda (p)
				   (or (not (insidep lines p))
				       (and prev
					    (> (scalar-product (v- prev actual) (v- p actual)) 0))))
				 (best-three actual)))))
      (let ((start (if (eq type 's) (rest (get-start (elt lines i))) (get-fuzzy-start i))))
	(cons start
	      (iter (for pprev first nil then prev)
		    (for prev first start then next)
		    (for next = (best-one pprev prev))
		    (while next)
		    (collect next)))))))

(defun trace-biquadratic (distance-type points i type parameter resolution)
  (let* ((n (length points))
	 (segments (list (elt points (mod (- i 2) n))
			 (elt points (mod (1- i) n))
			 (elt points (mod i n))
			 (elt points (mod (1+ i) n))))
	 (net (if (eq distance-type 'biquadratic)
		  (biquadratic-net points segments)
		  (biquadratic-corner-net points segments))))
    (iter (for x from 0 to 1 by resolution)
	  (collect (biquadratic-point net
				      (if (eq type 's)
					  (list parameter x)
					  (list x parameter)))))))

(defun vectorized-distance-function-test (points line-types filename &key (resolution 0.1d0)
					  (density 4) (distance-type 'perpendicular)
					  (color t))
  "LINE-TYPES is a list of symbols, each of which can be S, D or SD."
  (flet ((map-point (p)
	   (list (* (+ (first p) 1.0d0) 200)
		 (* (+ (second p) 1.0d0) 200))))
    (let* ((n (length points)) 
	   (lines (lines-from-points points))
	   (colors (or (and color (generate-colors n))
		       (iter (repeat n) (collect '(0 0 0)))))
	   (center (central-point points lines t)))
      (with-open-file (s filename :direction :output :if-exists :supersede)
	(format s "%!PS~%")
	(format s "~{~f ~}3 0 360 arc fill~%" (map-point center))
	(iter (for i from 0 below n)
	      (for color in colors)
	      (for line in lines)
	      (for line-type in line-types)
	      (format s "% Segment: ~a~%" i)
	      (format s "~{~d ~}setrgbcolor~%" color)
	      (format s "2 setlinewidth~%~
                         newpath~%~
                         ~{~f ~}moveto~%~
                         ~{~f ~}lineto~%~
                         stroke~%~
                         1 setlinewidth~%"
		      (map-point (first line))
		      (map-point (second line)))
	      (iter (for type in (case line-type (s '(s)) (d '(d)) (sd '(s d))))
		    (format s "% Type: ~a~%" type)
		    (iter (with d = (/ density))
			  (for parameter from 0 to 1 by d)
			  (format s "% Parameter: ~a~%" parameter)
			  (for trace =
			       (if (member distance-type '(biquadratic biquadratic-corner))
				   (trace-biquadratic distance-type points i type parameter resolution)
				   (trace-parameter points i distance-type type
						    parameter resolution)))
			  (format s "newpath~%")
			  (format s "~{~f ~}moveto~%~
                                     ~{~{~f ~}lineto~%~}"
				  (map-point (first trace))
				  (mapcar #'map-point (rest trace)))
			  (format s "stroke~%"))))
	(format s "showpage~%")))))

#+nil
(vectorized-distance-function-test
 (points-from-angles '(40 20 60 100 80)) '(sd nil nil nil sd) "/tmp/params.ps"
 :resolution 0.1d0 :density 4 :distance-type 'chord-based)

(defun biquadratic-net (points segments)
  "Uses Bezier curves to determine dangling corners."
  (let ((net (make-array '(3 3)))
	(n (length points)))
    (setf (aref net 0 0) (second segments)
	  (aref net 0 1) (affine-combine (second segments) 1/2 (third segments))
	  (aref net 0 2) (third segments)
	  (aref net 1 0) (affine-combine (first segments) 1/2 (second segments))
	  (aref net 1 1) (central-point points (lines-from-points points) t)
	  (aref net 1 2) (affine-combine (third segments) 1/2 (fourth segments)))
    (if (= n 3)
	(setf (aref net 2 0) (first segments)
	      (aref net 2 1) (first segments)
	      (aref net 2 2) (first segments))
	(let* ((i (position (fourth segments) points :test #'equal))
	       (curve (subseq (append points points) i (+ i n -2))))
	  (setf (aref net 2 0) (first segments)
		(aref net 2 1) (bezier curve 1/2)
		(aref net 2 2) (fourth segments))))
    net))

#+nil
(defun biquadratic-net (points segments)
  "Uses opposing corners as dangling corners."
  (let ((net (make-array '(3 3)))
	(n (length points)))
    (setf (aref net 0 0) (second segments)
	  (aref net 0 1) (affine-combine (second segments) 1/2 (third segments))
	  (aref net 0 2) (third segments)
	  (aref net 1 0) (affine-combine (first segments) 1/2 (second segments))
	  (aref net 1 1) (central-point points (lines-from-points points) t)
	  (aref net 1 2) (affine-combine (third segments) 1/2 (fourth segments)))
    (if (= n 3)
	(setf (aref net 2 0) (first segments)
	      (aref net 2 1) (first segments)
	      (aref net 2 2) (first segments))
	(let ((i (position (second segments) points :test #'equal))
	      (k (ceiling n 2)))
	  (setf (aref net 2 0) (first segments)
		(aref net 2 1) (affine-combine (elt points (mod (- i k -1) n))
					       1/2
					       (elt points (mod (+ i k) n)))
		(aref net 2 2) (fourth segments))))
    net))

(defun biquadratic-corner-net (points segments)
  "Uses Bezier curves to determine dangling corners."
  (let ((net (make-array '(3 3)))
	(n (length points)))
    (setf 
	  (aref net 0 0) (second segments)
	  (aref net 0 1) (affine-combine (second segments) 1/2 (third segments))
	  (aref net 0 2) (third segments)
	  (aref net 1 0) (affine-combine (first segments) 1/2 (second segments))
	  (aref net 1 1) (central-point points (lines-from-points points) t)
	  (aref net 2 0) (first segments))
    (if (= n 3)
	(setf (aref net 1 2) (third segments)
	      (aref net 2 1) (first segments)
	      (aref net 2 2) (affine-combine (third segments) 1/2 (fourth segments)))
	(let* ((i (position (third segments) points :test #'equal))
	       (k (ceiling n 2))
	       (curve1 (subseq (append points points) i (+ i n -2)))
	       (curve2 (subseq (append points points) (1+ i) (+ i n -1))))
	  (setf (aref net 1 2) (bezier curve1 1/2)
		(aref net 2 1) (bezier curve2 1/2)
		(aref net 2 2) (affine-combine (elt points (mod (- i 1 k) n))
					       1/2
					       (elt points (mod (+ i -1 k) n))))))
    net))

#+nil
(defun biquadratic-corner-net (points segments)
  "Uses opposing corners as dangling corners."
  (let ((net (make-array '(3 3)))
	(n (length points)))
    (setf (aref net 0 0) (second segments)
	  (aref net 0 1) (affine-combine (second segments) 1/2 (third segments))
	  (aref net 0 2) (third segments)
	  (aref net 1 0) (affine-combine (first segments) 1/2 (second segments))
	  (aref net 1 1) (central-point points (lines-from-points points) t)
	  (aref net 2 0) (first segments))
    (if (= n 3)
	(setf (aref net 1 2) (third segments)
	      (aref net 2 1) (first segments)
	      (aref net 2 2) (affine-combine (third segments) 1/2 (fourth segments)))
	(let* ((i (position (second segments) points :test #'equal))
	       (k (ceiling n 2)))
	  (setf (aref net 1 2) (affine-combine (elt points (mod (- i k) n))
					       1/2
					       (elt points (mod (+ i k -1) n)))
		(aref net 2 1) (affine-combine (elt points (mod (- i k -1) n))
					       1/2
					       (elt points (mod (+ i k) n)))
		(aref net 2 2) (affine-combine (elt points (mod (- i k) n))
					       1/2
					       (elt points (mod (+ i k) n))))))
    net))

(defun write-biquadratic-net (net filename)
  (labels ((map-point (p)
	     (list (* (+ (first p) 1.0d0) 200)
		   (* (+ (second p) 1.0d0) 200)))
	   (pp (i j) (map-point (aref net i j))))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      (dotimes (i 3)
	(dotimes (j 2)
	  (format s "newpath~%~{~f ~}moveto~%~{~f ~}lineto~%stroke~%"
		  (pp i j) (pp i (1+ j)))
	  (format s "newpath~%~{~f ~}moveto~%~{~f ~}lineto~%stroke~%"
		  (pp j i) (pp (1+ j) i))))
      (format s "showpage~%"))))

(defun biquadratic-point (net uv &optional derivative)
  (let ((u (first uv))
	(v (second uv)))
    (flet ((sqr (x) (* x x)))
      (let (u0 u1 u2 v0 v1 v2)
	(ecase derivative
	  ((nil)
	   (setf u0 (sqr (- 1 u)) u1 (* 2 u (- 1 u)) u2 (sqr u)
		 v0 (sqr (- 1 v)) v1 (* 2 v (- 1 v)) v2 (sqr v)))
	  (u
	   (setf u0 (* 2 (- u 1)) u1 (* 2 (- 1 (* 2 u))) u2 (* 2 u)
		 v0 (sqr (- 1 v)) v1 (* 2 v (- 1 v)) v2 (sqr v)))
	  (v
	   (setf u0 (sqr (- 1 u)) u1 (* 2 u (- 1 u)) u2 (sqr u)
		 v0 (* 2 (- v 1)) v1 (* 2 (- 1 (* 2 v))) v2 (* 2 v)))
	  (u2
	   (setf u0 2 u1 -4 u2 2
		 v0 (sqr (- 1 v)) v1 (* 2 v (- 1 v)) v2 (sqr v)))
	  (uv
	   (setf u0 (* 2 (- u 1)) u1 (* 2 (- 1 (* 2 u))) u2 (* 2 u)
		 v0 (* 2 (- v 1)) v1 (* 2 (- 1 (* 2 v))) v2 (* 2 v)))
	  (v2
	   (setf u0 (sqr (- 1 u)) u1 (* 2 u (- 1 u)) u2 (sqr u)
		 v0 2 v1 -4 v2 2)))
	(v+ (v* (v+ (v* (aref net 0 0) u0)
		    (v* (aref net 0 1) u1)
		    (v* (aref net 0 2) u2))
		v0)
	    (v* (v+ (v* (aref net 1 0) u0)
		    (v* (aref net 1 1) u1)
		    (v* (aref net 1 2) u2))
		v1)
	    (v* (v+ (v* (aref net 2 0) u0)
		    (v* (aref net 2 1) u1)
		    (v* (aref net 2 2) u2))
		v2))))))

;;; Newton-Raphson algorithm
;;; adapted from BSS-PROJECT-POINT & friends (projection.lisp)

(defun bq-projection-starting-value (net point res)
  (iter (with lower = '(0.0d0 0.0d0))
	(with upper = '(1.0d0 1.0d0))
	(with length = (mapcar #'- upper lower))
	(for i from 0 below (first res))
	(for u = (+ (first lower) (/ (* i (first length)) (1- (first res)))))
	(for (min v) =
	     (iter (for j from 0 below (second res))
		   (for v = (+ (second lower) (/ (* j (second length))
						 (1- (second res)))))
		   (for distance = (vlength
				    (v- point
					(biquadratic-point net (list u v)))))
		   (finding (list distance v) minimizing distance)))
	(finding (list u v) minimizing min)))

(defun bq-projection-delta (net uv r)
  (let* ((du (biquadratic-point net uv 'u))
	 (dv (biquadratic-point net uv 'v))
	 (du2 (biquadratic-point net uv 'u2))
	 (duv (biquadratic-point net uv 'uv))
	 (dv2 (biquadratic-point net uv 'v2))
	 (k (make-array '(2 1) :initial-contents
			`((,(- (scalar-product r du)))
			  (,(- (scalar-product r dv))))))
	 (J (make-array
	     '(2 2) :initial-contents
	     `((,(+ (scalar-product du du) (scalar-product r du2))
		 ,(+ (scalar-product du dv) (scalar-product r duv)))
	       (,(+ (scalar-product du dv) (scalar-product r duv))
		 ,(+ (scalar-product dv dv) (scalar-product r dv2))))))
	 (delta (matrix:multiplication (matrix:inverse-2x2 J) k)))
    (list (aref delta 0 0) (aref delta 1 0))))

(defun bq-halving-search (net point iterations &optional (distance-tolerance 0))
  (labels ((rec (lower upper iter)
	     (let* ((mid (affine-combine lower 0.5d0 upper))
		    (p (biquadratic-point net mid))
		    (distance (point-distance point p)))
	       (if (or (= iter 0) (< distance distance-tolerance))
		   (values mid distance)
		   (let ((d (v- point p))
			 (pu (biquadratic-point net mid 'u))
			 (pv (biquadratic-point net mid 'v))
			 (new-lower lower)
			 (new-upper upper))
		     (if (< (scalar-product d pu) 0)
			 (setf new-upper (list (first mid) (second new-upper)))
			 (setf new-lower (list (first mid) (second new-lower))))
		     (if (< (scalar-product d pv) 0)
			 (setf new-upper (list (first new-upper) (second mid)))
			 (setf new-lower (list (first new-lower) (second mid))))
		     (rec new-lower new-upper (1- iter)))))))
    (rec '(0 0) '(1 1) iterations)))

(defun bq-project-point (net point iterations &optional
			  (distance-tolerance 0.0) (cosine-tolerance 0.0))
  "Returns the parameters of the biquadratic defined by NET closest to POINT."
  (let ((uv0 (bq-halving-search net point iterations (* distance-tolerance 10)))
	(lower '(0.0d0 0.0d0))
	(upper '(1.0d0 2.0d0)))
    (iter (repeat iterations)
	  (for last first nil then uv)
	  (for uv first uv0 then
	       (mapcar (lambda (x next low upp) (min (max (+ x next) low) upp))
		       uv (bq-projection-delta net uv deviation)
		       lower upper))
	  (for p = (biquadratic-point net uv))
	  (for du = (biquadratic-point net uv 'u))
	  (for dv = (biquadratic-point net uv 'v))
	  (for deviation = (v- p point))
	  (when (or (<= (vlength deviation) distance-tolerance)
		    (and (<= (/ (abs (scalar-product du deviation))
				(* (vlength du) (vlength deviation)))
			     cosine-tolerance)
			 (<= (/ (abs (scalar-product dv deviation))
				(* (vlength dv) (vlength deviation)))
			     cosine-tolerance))
		    (and last
			 (<= (+ (vlength (v* du (- (first uv) (first last))))
				(vlength (v* dv (- (second uv) (second last)))))
			     distance-tolerance)))
	    (leave (values uv (vlength deviation))))
	  (finally (return (values uv (vlength deviation)))))))

(defun gsll-minimize-2d (fn start iterations deviation)
  (let ((obj (gsll:make-multi-dimensional-root-solver-fdf
	      gsll:+gnewton-mfdfsolver+
	      (list (lambda (u v)
		      (let ((f (funcall fn (list u v))))
			(values (first f) (second f))))
		    (lambda (u v)
		      (let ((fu (funcall fn (list u v) 'u))
			    (fv (funcall fn (list u v) 'v)))
			(values (first fu) (first fv)
				(second fu) (second fv))))
		    (lambda (u v)
		      (let ((f (funcall fn (list u v)))
			    (fu (funcall fn (list u v) 'u))
			    (fv (funcall fn (list u v) 'v)))
			(values (first f) (second f)
				(first fu) (first fv)
				(second fu) (second fv)))))
	      (grid:make-foreign-array 'double-float :dimensions 2
				       :initial-contents start))))
    (iter (repeat iterations)
	  (gsll:iterate obj)
	  (until (gsll:multiroot-test-residual obj deviation))
	  (finally (let ((result (gsll:solution obj)))
		     (return (list (grid:gref result 0) (grid:gref result 1))))))))

(defun bq-project-point-gsll (net point iterations deviation)
  (flet ((fn (uv &optional deriv)
	   (if deriv
	       (biquadratic-point net uv deriv)
	       (v- (biquadratic-point net uv) point))))
    (gsll-minimize-2d #'fn '(0.5d0 0.5d0) iterations deviation)))

(defmethod compute-distance ((type (eql 'biquadratic)) points segments p dir)
  (let ((sd (bq-project-point-gsll (biquadratic-net points segments) p 20 1.0d-8)))
    (if (eq dir 's)
	(first sd)
	(second sd))))

(defmethod compute-distance ((type (eql 'biquadratic-corner)) points segments p dir)
  (let ((sd (bq-project-point-gsll (biquadratic-corner-net points segments) p 20 1.0d-8)))
    (if (eq dir 's)
	(first sd)
	(second sd))))
