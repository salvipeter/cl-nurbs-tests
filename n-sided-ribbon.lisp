(in-package :cl-nurbs-tests)

(defparameter *barycentric-type* 'wachspress) ; or: meanvalue / harmonic
(defparameter *barycentric-normalized* t)
(defun barycentric-coordinates (points p)
  (let* ((vectors (mapcar (lambda (x) (v- p x)) points))
         (lengths (mapcar #'vlength vectors))
         (n (length points)))
    (labels ((inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (area (i)                  ; signed area = det(si,si+1)/2
               (let ((si (elt vectors i))
                     (si+1 (elt vectors (inc i))))
                 (/ (- (* (elt si 0) (elt si+1 1))
                       (* (elt si 1) (elt si+1 0)))
                    2.0d0)))
             (area-product (exceptions)
               (let ((r 1))
                 (iter (for j from 0 below n)
                       (unless (member j exceptions)
                         (setf r (* r (area j)))))
                 r))
             (bctype (x)
               (ecase *barycentric-type*
                 (harmonic (* x x))
                 (meanvalue x)
                 (wachspress 1))))
      (let* ((w (iter (for i from 0 below n)
                      (for Ai = (area-product (list i)))
                      (for Ai-1 = (area-product (list (dec i))))
                      (for Ai-1i = (area-product (list (dec i) i)))
                      (for Bi = (let ((si-1 (elt vectors (dec i)))
                                      (si+1 (elt vectors (inc i))))
                                  (/ (- (* (elt si-1 0) (elt si+1 1))
                                        (* (elt si-1 1) (elt si+1 0)))
                                     2.0d0)))
                      (for ri-1 = (bctype (elt lengths (dec i))))
                      (for ri = (bctype (elt lengths i)))
                      (for ri+1 = (bctype (elt lengths (inc i))))
                      (collect (+ (* ri-1 Ai-1)
                                  (* ri+1 Ai)
                                  (* -1 ri Bi Ai-1i)))))
             (wsum (reduce #'+ w)))
        (if *barycentric-normalized*
            (mapcar (lambda (wi) (/ wi wsum)) w)
            w)))))

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
      (*multiplier* 1))
  (iter (for p in (vertices *points*))
        (for l = (barycentric-coordinates *points* p))
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
  '((((-3 -25.8 0) (-1 -25.8 3) (0 -25.8 -1) (1 -25.8 3) (3 -25.8 0))
     ((3 -25.8 0) (12 -17.2 3) (21 -8.6 3) (30 0 0))
     ((30 0 0) (25 8.6 3) (20 17.2 3) (15 25.8 0))
     ((15 25.8 0) (5 25.8 3) (-5 25.8 3) (-15 25.8 0))
     ((-15 25.8 0) (-20 17.2 3) (-25 8.6 3) (-30 0 0))
     ((-30 0 0) (-21 -8.6 3) (-12 -17.2 3) (-3 -25.8 0)))
    (((-4 -17.2 18) (4 -17.2 18))
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
                    (* d d)
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
         (l (barycentric-coordinates points p))
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
         (l (barycentric-coordinates points p))
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

(defun write-ribbon-mesh (coords i filename)
  (let ((patch (generate-patch (first coords) (second coords)))
        (n (length (first coords)))
        (points '()))
    (iter (for ui from 0 below *resolution*)
          (for u = (/ ui (1- *resolution*)))
          (iter (for vi from 0 below *resolution*)
                (for v = (/ vi (1- *resolution*)))
                (push (ribbon-evaluate patch i
                                       (make-list n :initial-element u)
                                       (make-list n :initial-element v))
                      points)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "# vtk DataFile Version 1.0~
               ~%Bezier Surface~
               ~%ASCII~
               ~%DATASET POLYDATA~%~
               ~%POINTS ~d float~%~
               ~{~{~f~^ ~}~%~}~
               ~%POLYGONS ~d ~d~%"
              (* *resolution* *resolution*) points
              (* (1- *resolution*) (1- *resolution*) 2)
              (* (1- *resolution*) (1- *resolution*) 8))
      (dotimes (j (1- *resolution*))
        (dotimes (k (1- *resolution*))
          (format s "3 ~d ~d ~d~%"
                  (+ (* j *resolution*) k)
                  (+ (* j *resolution*) k 1)
                  (+ (* j *resolution*) k *resolution* 1))
          (format s "3 ~d ~d ~d~%"
                  (+ (* j *resolution*) k *resolution* 1)
                  (+ (* j *resolution*) k *resolution*)
                  (+ (* j *resolution*) k)))))))

#+nil
(let ((*resolution* 50)
      (*ribbon-multiplier* 1.0d0))
  (write-ribbon-mesh *coords* 2 "/tmp/ribbon.vtk"))

#+nil
(let ((*resolution* 50)
      (*ribbon-multiplier* 1.0d0))
  (write-constraint-ribbons *points* "/tmp/ribbons.vtk" :coords *coords* :resolution 20)
  (write-nr-patch *points* *coords* "/tmp/patch.vtk" :spider nil :ribbon 2))


;;; New s-parameterization test

(defparameter *use-local-d* t)
(defun compute-nr-parameter (l i dir)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (if (eq dir 's)
        (safe-/ (+ (elt l i) (elt l i+1))
                (+ (elt l i-2) (elt l i-1) (elt l i) (elt l i+1)))
        (if *use-local-d*
            (- 1 (safe-/ (+ (elt l i-1) (elt l i))
                         (+ (elt l i-2) (elt l i-1) (elt l i) (elt l i+1))))
            (- 1 (elt l i-1) (elt l i))))))

;;; See BITMAP-TEST in concave.lisp

#+nil
(let ((*barycentric-normalized* t)
      (*barycentric-type* 'meanvalue)
      (*use-local-d* nil))
  (flet ((s-fun (i)
           (lambda (points p)
             (let ((l (barycentric-coordinates points p)))
               (list (compute-nr-parameter l i 's)))))
         (d-fun (i)
           (lambda (points p)
             (let ((l (barycentric-coordinates points p)))
               (list (compute-nr-parameter l i 'd))))))
    (dotimes (i 6)
      (bitmap-test *points* (s-fun i)
                   (format nil "/tmp/side-~a-s-~:[global~;local~]-~(~a~)-~:[~;normalized~].pgm" i *use-local-d* *barycentric-type* *barycentric-normalized*)
                   :object-size 2.0d0)
      (bitmap-test *points* (d-fun i)
                   (format nil "/tmp/side-~a-d-~:[global~;local~]-~(~a~)-~:[~;normalized~].pgm" i *use-local-d* *barycentric-type* *barycentric-normalized*)
                   :object-size 2.0d0))))

;;; Note: `local' is only meaningful for d-parameters; s is always local

;;; Merge with a suitable transparent image of the polygon:
;;; for i in side*.pgm; do convert $i background.png -composite ${i%.pgm}.png; done
