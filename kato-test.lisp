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
    (if (< (abs a) *epsilon*)
	(- (/ c b))
	(let ((discr (sqrt (- (* b b) (* 4 a c)))))
	  (if (complexp discr)
	      0
	      (let* ((u1 (/ (+ (- b) discr) (* 2 a)))
		     (u2 (/ (- (- b) discr) (* 2 a)))
		     (u1m (min (max u1 0.0d0) 1.0d0))
		     (u2m (min (max u2 0.0d0) 1.0d0)))
		(if (< (abs (- u1 u1m)) (abs (- u2 u2m))) u1m u2m)))))))

(defun central-line-sweep (p0 p1 p2 p3 center filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%2~%" *width* *height*)
    (let* ((base-x (vnormalize (v- p2 p1)))
	   (base-y (list (second base-x) (- (first base-x))))
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
	   (base-y (list (second base-x) (- (first base-x))))
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

(defun perpendicular-distance (points segments p type)
  (let ((p0 (elt segments 1))
	(p1 (elt segments 2)))
    (if (eq type 'd)
	(/ (point-line-distance p (list p0 p1))
	   (iter (for q1 in points)
		 (maximize (iter (for q2 in points)
				 (maximize (point-distance q1 q2))))))
	(/ (abs (scalar-product (v- p0 p) (vnormalize (v- p1 p0))))
	   (point-distance p0 p1)))))

(defun barycentric-distance (lines segments p type)
  (flet ((area (segment)
	   (let* ((a (point-distance p (first segment)))
		  (b (point-distance p (second segment)))
		  (c (apply #'point-distance segment))
		  (s (/ (+ a b c) 2)))
	     (sqrt (max 0 (* s (- s a) (- s b) (- s c)))))))
    (if (eq type 'd)
	(/ (area (subseq segments 1 3))
	   (reduce #'+ (mapcar #'area lines)))
	(let ((di-1 (area (subseq segments 0 2)))
	      (di+1 (area (subseq segments 2))))
	  (/ di-1 (+ di-1 di+1))))))

(defun chord-based-distance (points segments p type)
  (if (eq type 'd)
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

(defun radial-distance (points segments p type)
  (let* ((p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments))
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (second base-x) (- (first base-x))))
	 (a (in-system base-x base-y (v- p0 p1)))
	 (b (in-system base-x base-y (v- p3 p2)))
	 (length (point-distance p1 p2))
	 (axby-aybx (- (* (first a) (second b)) (* (second a) (first b)))))
    (if (zerop axby-aybx)
	(perpendicular-distance points segments p type)
	(let ((laxby-div (/ (* length (first a) (second b)) axby-aybx))
	      (layby-div (/ (* length (second a) (second b)) axby-aybx)))
	  (destructuring-bind (x y)
	      (in-system base-x base-y (v- p p1))
	    (let ((x0 (+ x (* (/ y (- y layby-div)) (- laxby-div x)))))
	      (if (eq type 'd)
		  (/ (point-distance p (v+ p1 (v* base-x x0)))
		     (iter (for q1 in points)
			   (maximize (iter (for q2 in points)
					   (maximize (point-distance q1 q2))))))
		  (/ x0 length))))))))

(defun blend3 (x a b c)
  "Blends between three values, giving A for 0, B for 0.5, C for 1."
  (+ (* x x (+ (* 2 a) (* -4 b) (* 2 c)))
     (* x (+ (* -3 a) (* 4 b) (- c)))
     a))

(defun line-sweep-distance (center segments p type)
  (let* ((p0 (first segments))
	 (p1 (second segments))
	 (p2 (third segments))
	 (p3 (fourth segments)) 
	 (base-x (vnormalize (v- p2 p1)))
	 (base-y (list (second base-x) (- (first base-x))))
	 (center (in-system base-x base-y (v- center p1)))
	 (p (in-system base-x base-y (v- p p1)))
	 (length (point-distance p1 p2))
	 (w0 (in-system base-x base-y (v- p0 p1)))
	 (w2 (in-system base-x base-y (v- p3 p2)))
	 (w1 (center-line w0 w2 length center))
	 (s (line-sweep w0 w1 w2 length p)))
    (if (eq type 'd)
	(/ (point-distance p (list (* s length) 0))
	   (blend3 s (vlength w0)
		   (* 2 (point-distance center (list (/ length 2) 0)))
		   (vlength w2)))
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

(defmacro acond (&rest clauses)
  (if (null clauses)
      nil
      (let ((cl1 (car clauses))
            (sym (gensym)))
        `(let ((,sym ,(car cl1)))
           (if ,sym
               (let ((it ,sym)) ,@(cdr cl1))
               (acond ,@(cdr clauses)))))))

(defun distance-function-test (angles distance-type line-type filename)
  (let* ((points (mapcar (lambda (p)
			   (list (/ (+ *width* (* *width* (first p))) 2)
				 (/ (+ *height* (* *height* (second p))) 2)))
			 (points-from-angles angles)))
	 (lines (lines-from-points points))
	 (center (central-point points lines t)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "P3~%~d ~d~%255~%" *width* *height*)
      (labels ((parameter (lst type x y)
		 (ecase distance-type
		   (perpendicular
		    (perpendicular-distance points lst (list x y) type))
		   (barycentric
		    (barycentric-distance lines lst (list x y) type))
		   (chord-based
		    (chord-based-distance points lst (list x y) type))
		   (radial
		    (radial-distance points lst (list x y) type))
		   (line-sweep
		    (line-sweep-distance center lst (list x y) type))))
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
			   ((outsidep x y) (write-color s 'outside))
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
			   (t (write-color s 'nothing))))
	      (terpri s))))))

(defun distance-function-complete-set (angles directory)
  (iter (for type in '(perpendicular barycentric chord-based radial line-sweep))
	(iter (for lines in '((s d) (s s) (d d)))
	      (for name = (format nil "~(~a~)-~{~(~a~)~}.ppm" type lines))
	      (for filename = (make-pathname :directory directory :name name))
	      (distance-function-test angles type lines filename))))

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
  (distance-function-test '(40 20 60 100 80) 'chord-based '(s d)
			  "/tmp/distance.ppm"))

(defun connect-points (points)
  (let* ((n (length points))
	 (matrix (iter (for i from 0 below n)
		       (appending
			(iter (for j from (1+ i) below n)
			      (collect
			       (list i j (point-distance (elt points i) (elt points j))))))))
	 (sorted (sort (copy-list matrix) #'< :key #'third)))
    'todo))
