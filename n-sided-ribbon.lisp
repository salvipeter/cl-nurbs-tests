(in-package :cl-nurbs-tests)

(defun wachspress-coordinates (points p)
  (let ((lengths (mapcar (lambda (x) (point-distance p x)) points))
        (n (length points)))
    (labels ((heron% (a b c)
               (let ((s (/ (+ a b c) 2)))
                 (sqrt (* s (- s a) (- s b) (- s c)))))
             (heron (a b c)
               (let ((a (heron% a b c)))
                 (if (complexp a) 0 a)))
             (inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (area-a (i)
               (heron (elt lengths i)
                      (elt lengths (inc i))
                      (point-distance (elt points i)
                                      (elt points (inc i)))))
             (area-c (i)
               (heron (point-distance (elt points (dec i))
				      (elt points i))
		      (point-distance (elt points (dec i))
				      (elt points (inc i)))
		      (point-distance (elt points i)
				      (elt points (inc i)))))
             (area-product (A i)
               (iter (for j from 0 below n)
                     (unless (or (= j i) (= j (dec i)))
                       (collect (elt A j))))))
      (let* ((A (iter (for i from 0 below n) (collect (area-a i))))
             (w (iter (for i from 0 below n)
                      (for Ci = (area-c i))
                      (collect (reduce #'* (cons Ci (area-product A i))))))
             (w (if *squared-wachspress-coordinates* (mapcar #'* w w) w))
             (wsum (reduce #'+ w)))
        (mapcar (lambda (wi) (/ wi wsum)) w)))))

;;; Test: create a (for now 2D) ribbon with n sides

(defparameter *points* '((0 0) (3 0) (4 3) (2 6) (-1 2)))
(defparameter *target* '((8 0) (15 0) (18 5) (10 10) (7 4)))
(defparameter *target-cross* '((7 4) (18 5)))

;;; Szimmetrikus pelda
;; (defparameter *points* '((0 0) (3 0) (4 2) (1.5 3) (-1 2)))
;; (defparameter *target* '((8 0) (11 0) (12 2) (9.5 3) (7 2)))
;; (defparameter *target-cross* '((7 2) (12 2)))

(defvar *multiplier*)

(defun target-ribbon (s d)
  (let ((der (affine-combine (v- (elt *target-cross* 0)
                                 (elt *target* 0))
                             s
                             (v- (elt *target-cross* 1)
                                 (elt *target* 1))))
        (p (affine-combine (elt *target* 0) s (elt *target* 1))))
    (v+ p (v* der d *multiplier*))))

(let ((*resolution* 10)
      (*squared-wachspress-coordinates* nil)
      (*multiplier* 1))
  (iter (for p in (vertices *points*))
        (for l = (wachspress-coordinates *points* p))
        #+nil(for s = (elt (compute-parameter 'bilinear 's *points* p t) 1))
        (for s = (let ((dprev (apply #'- 1 (elt l 0) (last l)))
                       (dnext (- 1 (elt l 1) (elt l 2))))
                   (/ dprev (+ dprev dnext))))
        (for d = (- 1 (elt l 0) (elt l 1)))
        (collect (v+ (v* (target-ribbon s d) (- 1 d))
                     (v* (reduce #'v+ (mapcar #'v* *target* l)) d))
          into result)
        (finally
         (with-open-file (s "/tmp/domain" :direction :output :if-exists :supersede)
           (format s "~{~{~f~^ ~}~%~}" (append *points* (list (first *points*)))))
         (with-open-file (s "/tmp/target" :direction :output :if-exists :supersede)
           (format s "~{~{~f~^ ~}~%~}" (append *target* (list (first *target*)))))
         (with-open-file (s "/tmp/result" :direction :output :if-exists :supersede)
           (format s "~{~{~f~^ ~}~%~}" result)))))

;;; plot "/tmp/domain" with lines, "/tmp/target" with lines, "/tmp/result"


;;; Testing with 3D ribbons & patches

(defparameter *coords*
  '((((-3 -25.8 0) (-1 -25.8 3) (1 -25.8 3) (3 -25.8 0))
     ((3 -25.8 0) (12 -17.2 3) (21 -8.6 3) (30 0 0))
     ((30 0 0) (25 8.6 3) (20 17.2 3) (15 25.8 0))
     ((15 25.8 0) (5 25.8 3) (-5 25.8 3) (-15 25.8 0))
     ((-15 25.8 0) (-20 17.2 3) (-25 8.6 3) (-30 0 0))
     ((-30 0 0) (-21 -8.6 3) (-12 -17.2 3) (-3 -25.8 0)))
    (((-4 -17.2 8) (4 -17.2 8))
     ((4 -17.2 8) (18 0 8))
     ((18 0 8) (10 17.2 8))
     ((10 17.2 8) (-10 17.2 8))
     ((-10 17.2 8) (-18 0 8))
     ((-18 0 8) (-4 -17.2 8)))))

(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((3 -25.8 0) (30 0 0) (15 25.8 0) (-15 25.8 0) (-30 0 0) (-3 -25.8 0)))))

(defun generate-n-sided-ribbons (outer inner)
  (generate-patch outer inner))

(defun nr-patch-vertices (ribbons)
  (let ((lst (mapcar #'first (first ribbons))))
    (append (rest lst) (list (first lst)))))

(defun n-sided-ribbon-evaluate (ribbons barycoords i)
  (let* ((n (length barycoords))
         (i-1 (mod (- i 1) n))
         (d (- 1 (elt barycoords i-1) (elt barycoords i)))
         (s (let ((dprev (- 1 (elt barycoords (mod (- i 2) n)) (elt barycoords i-1)))
                  (dnext (- 1 (elt barycoords i) (elt barycoords (mod (1+ i) n)))))
              (/ dprev (+ dprev dnext)))))
    (affine-combine (ribbon-evaluate ribbons i
                                     (make-list n :initial-element s)
                                     (make-list n :initial-element d))
                    d
                    (reduce #'v+ (mapcar #'v* (nr-patch-vertices ribbons) barycoords)))))

(defun n-sided-corner-correction (ribbons barycoords i)
  "TODO"
  (let* ((i-1 (mod (1- i) (length barycoords)))
	 (si-1 (if d (elt d i) (- 1.0d0 (elt s i-1))))
	 (si (elt s i))
	 (previous (let ((lst (elt (first ribbons) i-1)))
		     (elt lst (- (length lst) 2))))
	 (corner (first (elt (first ribbons) i)))
	 (next (second (elt (first ribbons) i)))
	 (twist (second (elt (second ribbons) i))))
    (v+ corner
	(v* (v- previous corner) 3.0d0 (gamma si-1))
	(v* (v- next corner) 3.0d0 (gamma si))
	(v* (v- (v+ corner twist) (v+ previous next)) 9.0d0 (gamma si-1) (gamma si)))))

(defun nr-patch-evaluate (ribbons points domain-point)
  (let* ((n (length points))
	 (p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
         (l (wachspress-coordinates points p))
         (d (mapcar (lambda (x y) (- 1 x y)) (cons (car (last l)) l) l))
         (*use-gamma* nil))
    (iter (for i from 0 below n)
	  (with result = '(0 0 0))
	  (setf result
		(v+ result
                    (v- (v* (n-sided-ribbon-evaluate ribbons l i)
                            (+ (corner-blend d (mod (1- i) n))
                               (corner-blend d i)))
                        (v* (n-sided-corner-correction ribbons l i)
                            (corner-blend d (mod (1- i) n))))))
	  (finally (return result)))))

(defun nr-one-ribbon-evaluate (ribbons points domain-point i)
  (let* ((p (mapcar (lambda (x) (or (and (>= (abs x) *tiny*) x) 0.0d0)) domain-point))
         (l (wachspress-coordinates points p))
         (*use-gamma* nil))
    (n-sided-ribbon-evaluate ribbons l i)))

(defun write-nr-patch (points coords filename &key spider ribbon)
  (let* ((n (length points))
	 (ribbons (generate-n-sided-ribbons (first coords) (second coords))))
    (if spider
	(write-vtk-polylines
	 (iter (for line in (spider-lines points))
	       (collect (iter (for domain-point in line)
                              (for p =
                                   (if ribbon
                                       (nr-one-ribbon-evaluate ribbons points domain-point ribbon)
                                       (nr-patch-evaluate ribbons points domain-point)))
			      (collect p))))
	 filename)
	(write-vtk-indexed-mesh
	 (iter (for domain-point in (vertices points))
               (for p =
                    (if ribbon
                        (nr-one-ribbon-evaluate ribbons points domain-point ribbon)
                        (nr-patch-evaluate ribbons points domain-point)))
	       (collect p))
	 (triangles n) filename))))

#+nil
(let ((*resolution* 30)
      (*ribbon-multiplier* 1.0d0))
  (write-constraint-ribbons *points* "/tmp/ribbons.vtk" :coords *coords* :resolution 20)
  (write-nr-patch *points* *coords* "/tmp/patch.vtk" :spider nil :ribbon 2))
