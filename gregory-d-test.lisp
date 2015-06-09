(in-package :cl-nurbs-tests)

(defun corner-d-evaluate (patch i d)
  "The corner defined by segments I-1 and I,
thus containing point I-1 (NOT point I)."
  (let* ((i-1 (mod (1- i) (length d)))
	 (si-1 (- 1.0d0 (elt d i)))
	 (si (elt d i-1))
	 (ci-1 (bezier (elt (first patch) i-1) si-1))
	 (ci (bezier (elt (first patch) i) si))
	 (di-1 (v* (v- (bezier (elt (second patch) i-1) si-1) ci-1) 3.0d0))
	 (di (v* (v- (bezier (elt (second patch) i) si) ci) 3.0d0)))
    (v+ ci ci-1 (v* di (- 1.0d0 si-1)) (v* di-1 si))))

(defun corner-d-correction (patch i d)
  (let* ((i-1 (mod (1- i) (length d)))
	 (si-1 (elt d i))
	 (si (elt d i-1))
	 (previous (let ((lst (elt (first patch) i-1)))
		     (elt lst (- (length lst) 2))))
	 (corner (first (elt (first patch) i)))
	 (next (second (elt (first patch) i)))
	 (twist (second (elt (second patch) i))))
    (v+ corner
	(v* (v- previous corner) 3.0d0 (gamma si-1))
	(v* (v- next corner) 3.0d0 (gamma si))
	(v* (v- (v+ corner twist) (v+ previous next)) 9.0d0 (gamma si-1) (gamma si)))))

(defvar *quintic-baryblend-p* nil)
(defun corner-baryblend (d i)
  (let* ((n (length d))
         (i-1 (mod (1- i) n))
         (di (elt d i))
         (di-1 (elt d i-1)))
    (if *quintic-baryblend-p*
        (* (+ (bernstein 5 0 di-1)
              (bernstein 5 1 di-1)
              (bernstein 5 2 di-1))
           (+ (bernstein 5 0 di)
              (bernstein 5 1 di)
              (bernstein 5 2 di)))
        (* (+ (bernstein 3 0 di-1)
              (bernstein 3 1 di-1))
           (+ (bernstein 3 0 di)
              (bernstein 3 1 di))))))

(defun corner-sbaryblend (s i)
  (let* ((n (length s))
         (i-1 (mod (1- i) n))
         (di (- 1 (elt s i-1)))
         (di-1 (elt s i)))
    (if *quintic-baryblend-p*
        (* (+ (bernstein 5 0 di-1)
              (bernstein 5 1 di-1)
              (bernstein 5 2 di-1))
           (+ (bernstein 5 0 di)
              (bernstein 5 1 di)
              (bernstein 5 2 di)))
        (* (+ (bernstein 3 0 di-1)
              (bernstein 3 1 di-1))
           (+ (bernstein 3 0 di)
              (bernstein 3 1 di))))))

(defun corner-normalized-baryblend (d i)
  (/ (corner-baryblend d i)
     (iter (for j from 0 below (length d))
           (sum (corner-baryblend d j)))))

;;; Added (temporarily): CORNER[-D][-[NORMALIZED-][S]BARYBLEND] types
(defun patch-evaluate (patch points type distance-type domain-point)
  (let* ((n (length points))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
	 (d (unless (and (eq type 'corner) (eq distance-type 'biquadratic))
	      (compute-parameter distance-type 'd points p t)))
	 (dc (if (and (member type '(corner hybrid)) (eq distance-type 'biquadratic))
		 (compute-parameter 'biquadratic-corner 'd points p t)
		 d))
	 (s (unless (and (eq type 'corner) (eq distance-type 'biquadratic))
	      (compute-parameter distance-type 's points p t)))
	 (sc (if (and (member type '(corner hybrid)) (eq distance-type 'biquadratic))
		 (compute-parameter 'biquadratic-corner 's points p t)
		 s))
	 (b (and (eq type 'sketches) (compute-parameter 'perpendicular 'd points p t)))
         (*use-gamma* (member type '(hybrid hybrid-coons))))
    (iter (for i from 0 below n)
	  (with result = '(0 0 0))
	  (setf result
		(v+ result
		    (case type
		      (ribbon (v* (ribbon-evaluate patch i s d)
				  (ribbon-blend d i)))
		      (ribbon-coons (v* (coons-ribbon-evaluate patch i s d)
					(ribbon-blend d i)))
		      (corner (if (eq distance-type 'biquadratic)
				  (v* (v- (corner-evaluate patch i sc dc)
					  (corner-correction patch i sc dc))
				      (corner-blend dc (mod (1- i) n)))
				  (v* (v- (corner-evaluate patch i s)
					  (corner-correction patch i s))
				      (corner-blend d (mod (1- i) n)))))
                      (corner-baryblend
                       (v* (v- (corner-evaluate patch i s)
                               (corner-correction patch i s))
                           (corner-baryblend d i)))
                      (corner-normalized-baryblend
                       (v* (v- (corner-evaluate patch i s)
                               (corner-correction patch i s))
                           (corner-normalized-baryblend d i)))
                      (corner-sbaryblend
                       (v* (v- (corner-evaluate patch i s)
                               (corner-correction patch i s))
                           (corner-sbaryblend s i)))
                      (corner-d (v* (v- (corner-d-evaluate patch i d)
                                        (corner-d-correction patch i d))
                                    (corner-blend d (mod (1- i) n))))
                      (corner-d-baryblend
                       (v* (v- (corner-d-evaluate patch i d)
                               (corner-d-correction patch i d))
                           (corner-baryblend d i)))
                      (corner-d-normalized-baryblend
                       (v* (v- (corner-d-evaluate patch i d)
                               (corner-d-correction patch i d))
                           (corner-normalized-baryblend d i)))
		      (hybrid (v- (v* (ribbon-evaluate patch i s d)
				      (+ (corner-blend d (mod (1- i) n))
					 (corner-blend d i)))
				  (if (eq distance-type 'biquadratic)
				      (v* (corner-correction patch i sc dc)
					  (corner-blend d (mod (1- i) n)))
				      (v* (corner-correction patch i s)
					  (corner-blend d (mod (1- i) n))))))
		      (hybrid-coons (v* (coons-ribbon-evaluate patch i s d)
					(+ (wachspress-corner-blend points (mod (1- i) n) p)
					   (wachspress-corner-blend points i p))
					1/2))
		      (sketches (v* (ribbon-evaluate patch i s d)
				    (ribbon-blend b i)))
		      (sketches-coons (v* (coons-ribbon-evaluate patch i s d)
					  (ribbon-blend b i))))))
	  (finally (return result)))))


;;; Test:

;;; 5sided
(defparameter *points* (points-from-angles '(60 20 110 60 110)))
(defparameter *coords*
  '((((0.0d0 0.0d0 0.0d0)
      (1.0d0 0.0d0 1.0d0)
      (2.0d0 0.0d0 1.0d0)
      (2.4d0 0.0d0 0.3d0))
     ((2.4d0 0.0d0 0.3d0)
      (2.6d0 0.2d0 0.4d0)
      (2.8d0 0.4d0 0.4d0)
      (3.0d0 0.6d0 0.3d0))
     ((3.0d0 0.6d0 0.3d0)
      (3.0d0 2.0d0 1.0d0)
      (3.0d0 4.0d0 1.0d0)
      (3.0d0 6.0d0 0.0d0))
     ((3.0d0 6.0d0 0.0d0)
      (2.0d0 6.0d0 1.0d0)
      (1.0d0 6.0d0 1.0d0)
      (0.0d0 6.0d0 0.0d0))
     ((0.0d0 6.0d0 0.0d0)
      (0.0d0 4.0d0 1.0d0)
      (0.0d0 2.0d0 1.0d0)
      (0.0d0 0.0d0 0.0d0)))
    (((1.0d0 2.0d0 1.2d0)
      (2.2d0 0.7d0 1.2d0))
     ((2.2d0 0.7d0 1.2d0)
      (2.8d0 1.8d0 1.1d0))
     ((2.8d0 1.8d0 1.1d0)
      (2.0d0 4.0d0 1.2d0))
     ((2.0d0 4.0d0 1.2d0)
      (1.0d0 4.0d0 1.2d0))
     ((1.0d0 4.0d0 1.2d0)
      (1.0d0 2.0d0 1.2d0)))))


;;; 4sided
(defparameter *coords*
  '((((174.683d0 -79.2657d0 39.9866d0)
      (147.39d0 -69.055d0 34.5583d0)
      (-92.2374d0 17.9945d0 -13.5014d0)
      (-146.776d0 37.1868d0 -22.987d0))
     ((-146.776d0 37.1868d0 -22.987d0)
      (-89.7339d0 42.4423d0 -29.727d0)
      (171.56d0 64.1718d0 -60.5931d0)
      (200.335d0 67.0522d0 -63.9375d0))
     ((200.335d0 67.0522d0 -63.9375d0)
      (210.949d0 44.1498d0 -48.4583d0)
      (225.487d0 13.4434d0 -27.1227d0)
      (237.86d0 -12.5697d0 -9.64074d0))
     ((237.86d0 -12.5697d0 -9.64074d0)
      (217.464d0 -34.1532d0 6.21314d0)
      (192.862d0 -60.3934d0 26.2093d0)
      (174.683d0 -79.2657d0 39.9866d0)))
    (((164.961d0 -50.4667d0 20.407d0)
      (-39.6534d0 23.6076d0 -19.1114d0))
     ((-39.6534d0 23.6076d0 -19.1114d0)
      (181.71d0 41.4307d0 -44.9136d0))
     ((181.71d0 41.4307d0 -44.9136d0)
      (205.954d0 -8.0374d0 -10.5976d0))
     ((205.954d0 -8.0374d0 -10.5976d0)
      (164.961d0 -50.4667d0 20.407d0)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((-146.776 37.1868 -22.987)
      (200.335 67.0522 -63.9375)
      (237.86 -12.5697 -9.64074)
      (174.683 -79.2657 39.9866)))))


;;; 3sided
(defparameter *points* (points-from-angles '(5 120 150)))
(defparameter *coords*
  '((((0 0 0) (1 0 0) (2 0 0) (3 0 0))
     ((3 0 0) (2 1 0) (1 2 0) (0 3 0))
     ((0 3 0) (0 2 0) (0 1 0) (0 0 0)))
    (((1 1 1) (1.5 1 1))
     ((1.5 1 1) (0.5 2 1))
     ((0.5 2 1) (1 1 1)))))

(let ((*resolution* 100)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* nil)
      (*wachspressp* t)
      (*quintic-baryblend-p* nil)
      (*use-gamma* nil))
  (iter (for type in '(corner corner-normalized-baryblend corner-d corner-d-normalized-baryblend))
	(iter (for distance in '(mean-bilinear))
	      (format t "POS [~a / ~a]: ~f~%"
		      type distance
		      (check-patch *points* type :coords *coords* :distance-type distance))
	      (format t "TAN [~a / ~a]: ~f~%"
		      type distance
		      (/ (* 180 (check-patch-normals *points* type :coords *coords* :distance-type distance :step 1.0d-4)) pi)))))


;;; Blend deficit at a given point
#+nil
(let ((*quintic-baryblend-p* nil))
  (iter (for i from 3 to 9)
        (let* ((points (points-from-angles (uniform-angles i)))
               (n (length points))
               (lines (lines-from-points points))
               (p (affine-combine (affine-combine (first points) 0.5 (second points))
                                  0.5
                                  (central-point points lines t)))
               (d (compute-parameter 'mean-bilinear 'd points p t)))
          (format t "~d => ~a~%"
                  i (- 1 (iter (for i from 0 below n)
                               (print (corner-baryblend d i))
                               (sum (corner-baryblend d i))))))))
;;; cubic:
;; 3 => -0.6460905349794241
;; 4 => 2.220446049250313e-16
;; 5 => 0.38048000000000004
;; 6 => 0.5967078189300408
;; 7 => 0.7248765395370971
;; 8 => 0.8046874999999999
;; 9 => 0.8566614167894461
;;; quintic:
;; 3 => -0.8728852309099222
;; 4 => 2.220446049250313e-16
;; 5 => 0.4961592320000001
;; 6 => 0.7357110196616365
;; 7 => 0.8534301035344867
;; 8 => 0.9142761230468749
;; 9 => 0.9474520719011327
