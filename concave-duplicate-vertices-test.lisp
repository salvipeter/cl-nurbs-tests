(in-package :cl-nurbs-tests)

;;; As in CONCAVE.LISP
(defun bitmap-test (points sd-fun filename &key (object-size 2.0d0) (size 400) (density 10))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "P2~%~d ~d~%255~%" size size)
    (labels ((map-point (p)
               (v- (v* p (/ object-size size))
                   (list (/ object-size 2) (/ object-size 2))))
             (query (p)
               (reduce #'min (mapcar (lambda (x)
                                       (multiple-value-bind (quot rem) (round x (/ density))
                                         (declare (ignore quot))
                                         (abs rem)))
                                     (funcall sd-fun points p)))))
      (iter (for y from (1- size) downto 0)
            (iter (for x from 0 below size)
                  (for p = (map-point (list x y)))
                  (for color = (min (floor (* (query p) density 255)) 255))
                  (format s "~d " color))
            (terpri s)))))

;;; Some s/d parameter line tests

(defparameter *points*
  '((-2 0) (2 0) (2 2) (1 2) (1 1) (1 1) (-1 1) (-1 1) (-1 2) (-2 2)))

(defparameter *points2*
  '((-2 -1) (2 -1) (2 2) (1 2) (1 1) (1 1) (-1 1) (-1 1) (-1 2) (-2 2)))

#+nil
(let ((d '(1 1 1 1 1 0 0 1 1 1))
      (s '(1 1 1 1 1 1 0 0 1 1)))
  (bitmap-test *points* (compute-sd s d) "/tmp/proba.pgm"
               :object-size 5 :density 10))

#+nil
(let ((d '(1 1 1 1 1 0 0 1 1 1)))
  (bitmap-test *points2* (compute-sd d d) "/tmp/proba2.pgm"
               :object-size 5 :density 10))

#+nil
(let ((d '(1 1 1 1 1 1 0 0 1 1)))
  (bitmap-test *points2* (compute-sd d d) "/tmp/proba3.pgm"
               :object-size 5 :density 10))

#+nil
(let ((d '(0 1 1 1 1 1 1 1 1 0)))
  (bitmap-test *points2* (compute-sd d d) "/tmp/proba4.pgm"
               :object-size 5 :density 10))

(defparameter *points3*
  '((-10 -10) (-10 -10) (10 -10) (10 -10) (20 0) (0 20) (-20 0)))

#+nil
(let ((d '(1 0 0 1 1 1 1)))
  (bitmap-test *points3* (compute-sd d d) "/tmp/proba5.pgm"
               :object-size 40 :density 10))

#+nil
(let ((d '(0 0 1 1 1 1 1)))
  (bitmap-test *points3* (compute-sd d d) "/tmp/proba6.pgm"
               :object-size 40 :density 10))

#+nil
(let ((d '(0 1 1 1 1 1 0)))
  (bitmap-test *points3* (compute-sd d d) "/tmp/proba7.pgm"
               :object-size 40 :density 10))


;;; Patch tests

(defparameter *points*
  '((-1 -1) (1 -1) (1 0.3333d0) (0.3333d0 0.3333d0)
    (0.3333d0 -0.3333d0) (0.3332d0 -0.3333d0)
    (-0.3332d0 -0.3333d0) (-0.3333d0 -0.3333d0)
    (-0.3333d0 0.3333d0) (-1 0.3333d0)))

(defparameter *coords*
  '((((0 6 0) (0 5 0) (0 1 0) (0 0 0))
     ((0 0 0) (1 0 0) (8 0 0) (9 0 0))
     ((9 0 0) (9 1 0) (9 5 0) (9 6 0))
     ((9 6 0) (8 6 0) (7 6 0) (6 6 0))
     ((6 6 0) (6 5 0) (6 4 0) (6 3 0))
     ((6 3 0) (6 3 0) (5.9999 3 0) (5.9999 3 0))
     ((5.9999 3 0) (5 3 0) (4 3 0) (3.0001 3 0))
     ((3.0001 3 0) (3.0001 3 0) (3 3 0) (3 3 0))
     ((3 3 0) (3 4 0) (3 5 0) (3 6 0))
     ((3 6 0) (2 6 0) (1 6 0) (0 6 0)))
    (((1 5 1) (1 1 1))
     ((1 1 1) (8 1 1))
     ((8 1 1) (8 5 1))
     ((8 5 1) (7 5 1))
     ((7 5 1) (7 3 1))
     ((7 3 1) (6 2 1))
     ((6 2 1) (3 2 1))
     ((3 2 1) (2 3 1))
     ((2 3 1) (2 5 1))
     ((2 5 1) (1 5 1)))))

#+nil
(write-constraint-ribbons
 *points*
 "/tmp/ribbons.vtk"
 :coords *coords*
 :resolution 20)

#+nil
(let ((*resolution* 50)
      (*ribbon-multiplier* 1.0d0)
      (*wachspressp* nil))
  (write-patch *points* 'ribbon "/tmp/sp-spider.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider t)
  (write-patch *points* 'ribbon "/tmp/sp-patch.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider nil))

#+nil
(let ((*wachspressp* nil))
  (flet ((s-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 's points p) i))))
         (d-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 'd points p) i)))))
    (dotimes (i 7)
      (bitmap-test *points* (s-fun i)
                   (format nil "/tmp/~as.pgm" i)
                   :object-size 2.0d0)
      (bitmap-test *points* (d-fun i)
                   (format nil "/tmp/~ad.pgm" i)
                   :object-size 2.0d0))))


;;; Test 2

(defparameter *points*
  '((-1 -1) (-1.0d-5 -0.5) (1.0d-5 -0.5) (1 -1) (0 1)))

(defparameter *coords*
  '((((0 1 0) (-0.1 0.8 0) (-0.9 -0.8 0) (-1 -1 0))
     ((-1 -1 0) (-0.8 -0.9 0) (-0.2 -0.6 0) (-1.0d-5 -0.5 0))
     ((-1.0d-5 -0.5 0) (-1.0d-5 -0.5 0) (1.0d-5 -0.5 0) (1.0d-5 -0.5 0))
     ((1.0d-5 -0.5 0) (0.2 -0.6 0) (0.8 -0.9 0) (1 -1 0))
     ((1 -1 0) (0.9 -0.8 0) (0.1 0.8 0) (0 1 0)))
    (((0 0.6 0.1) (-0.7 -0.7 0.1))
     ((-0.7 -0.7 0.1) (-0.2 -0.4 0.1))
     ((-0.2 -0.4 0.1) (0.2 -0.4 0.1))
     ((0.2 -0.4 0.1) (0.7 -0.7 0.1))
     ((0.7 -0.7 0.1) (0 0.6 0.1)))))

#+nil
(write-constraint-ribbons
 *points*
 "/tmp/ribbons2.vtk"
 :coords *coords*
 :resolution 20)

#+nil
(let ((*resolution* 30)
      (*ribbon-multiplier* 1.0d0)
      (*wachspressp* nil))
  (write-patch *points* 'ribbon "/tmp/sp-spider2.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider t)
  (write-patch *points* 'ribbon "/tmp/sp-patch2.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider nil))

#+nil
(let ((*wachspressp* nil))
  (flet ((s-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 's points p) i))))
         (d-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 'd points p) i)))))
    (dotimes (i 3)
      (bitmap-test *points* (s-fun i)
                   (format nil "/tmp/~as.pgm" i)
                   :object-size 2.0d0)
      (bitmap-test *points* (d-fun i)
                   (format nil "/tmp/~ad.pgm" i)
                   :object-size 2.0d0))))


;;; Test 3

(defparameter *points*
  '((-1 -1) (1 -1) (1 1) (0 1) (1.0d-5 0) (-1.0d-5 0) (-1 0)))

(defparameter *coords*                  ; a konkav csucsban egy normalis van csak
  '((((-1 0 0) (-1 -0.25 0) (-1 -0.75 0) (-1 -1 0))
     ((-1 -1 0) (-0.75 -1 0) (0.75 -1 0) (1 -1 0))
     ((1 -1 0) (1 -0.75 0) (1 0.75 0) (1 1 0))
     ((1 1 0) (0.75 1 0) (0.25 1 0) (0 1 0))
     ((0 1 0) (0 0.75 0) (0 0.25 0) (0 0 0))
     ((0 0 0) (0 0 0) (0 0 0) (0 0 0))
     ((0 0 0) (-0.25 0 0) (-0.75 0 0) (-1 0 0)))
    (((-0.75 -0.25 0.2) (-0.75 -0.75 0.2))
     ((-0.75 -0.75 0.2) (0.75 -0.75 0.2))
     ((0.75 -0.75 0.2) (0.75 0.75 0.2))
     ((0.75 0.75 0.2) (0.25 0.75 0.2))
     ((0.25 0.75 0.2) (0.25 0 0))
     ((0.25 0 0) (0 -0.25 0))
     ((0 -0.25 0) (-0.75 -0.25 0.2)))))

#+nil
(write-constraint-ribbons
 *points*
 "/tmp/ribbons3.vtk"
 :coords *coords*
 :resolution 20)

#+nil
(let ((*resolution* 30)
      (*ribbon-multiplier* 1.0d0)
      (*wachspressp* nil))
  (write-patch *points* 'ribbon "/tmp/sp-spider3.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider t)
  (write-patch *points* 'ribbon "/tmp/sp-patch3.vtk"
	       :coords *coords* :distance-type 'mean-kato :spider nil))

#+nil
(let ((*wachspressp* nil))
  (flet ((s-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 's points p) i))))
         (d-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato 'd points p) i)))))
    (dotimes (i 6)
      (bitmap-test *points* (s-fun i)
                   (format nil "/tmp/~as.pgm" i)
                   :object-size 2.0d0)
      (bitmap-test *points* (d-fun i)
                   (format nil "/tmp/~ad.pgm" i)
                   :object-size 2.0d0))))


;;; Idea: compute the mean value as 1 - l_i - l_{i+1}
;;; => nothing changes :(

(defun mean-value-list (points p)
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
                    2.0d0))))
      (let* ((w (iter (for i from 0 below n)
                      (for Ai = (area i))
                      (for Ai-1 = (area (dec i)))
                      (for Di = (scalar-product (elt vectors (inc i)) (elt vectors i)))
                      (for Di-1 = (scalar-product (elt vectors i) (elt vectors (dec i))))
                      (for ri-1 = (elt lengths (dec i)))
                      (for ri = (elt lengths i))
                      (for ri+1 = (elt lengths (inc i)))
                      (collect (+ (if (> (abs Ai-1) *epsilon*)
                                      (/ (- ri-1 (/ Di-1 ri)) Ai-1)
                                      0)
                                  (if (> (abs Ai) *epsilon*)
                                      (/ (- ri+1 (/ Di ri)) Ai)
                                      0)))))
             (wsum (reduce #'+ w)))
        (mapcar (lambda (wi) (/ wi wsum)) w)))))

(defun mean-distance-direct (points segments p)
  (let* ((n (length points))
         (mv (mean-value-list points p))
         (i (position (second segments) points)))
    (- 1 (elt mv i) (elt mv (mod (1+ i) n)))))

(defmethod compute-distance ((type (eql 'mean-kato2)) points segments p dir)
  (if (eq dir 'd)
      (mean-distance-direct points segments p)
      (let ((dprev (mean-distance-direct points (segments-prev points segments) p))
            (dnext (mean-distance-direct points (segments-next points segments) p)))
        (let ((dd (+ dprev dnext)))
          (if (< dd *epsilon*)
              0.5
              (/ dprev dd))))))

#+nil
(flet ((s-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato2 's points p) i))))
         (d-fun (i)
           (lambda (points p)
             (list (elt (compute-parameter 'mean-kato2 'd points p) i)))))
    (let ((i 4))
      #+nil(bitmap-test *points* (s-fun i)
                   (format nil "/tmp/~as.pgm" i)
                   :object-size 2.0d0)
      (bitmap-test *points* (d-fun i)
                   (format nil "/tmp/~ad-c.pgm" i)
                   :object-size 2.0d0)))
