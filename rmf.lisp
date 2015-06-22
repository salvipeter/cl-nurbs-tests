(in-package :cl-nurbs-tests)

(defstruct frame position tangent binormal normal)

(defun frenet-frame (curve prev-frame u)
  (declare (ignore prev-frame))
  (let* ((p (bsc-evaluate curve u))
         (d (bsc-evaluate curve u :derivative 1))
         (dd (bsc-evaluate curve u :derivative 2))
         (n (vnormalize (cross-product d dd))))
    (make-frame :position p :tangent (vnormalize d)
                :binormal (vnormalize (cross-product n d)) :normal n)))

(defun next-frame (curve prev-frame u)
  "As described in `Computation of Rotation Minimizing Frames', Wang et al., 2008.
A limitation of this method is that it is determined by the starting frame and the curve tangents,
so an end frame cannot be supplied. As a workaround, we can add a rotation gradually,
by minimizing the total (or total squared) angular speed, see section 6.3 in the paper."
  (let* ((p (frame-position prev-frame))
         (d (frame-tangent prev-frame))
         (n (frame-normal prev-frame))
         (p-next (bsc-evaluate curve u))
         (d-next (vnormalize (bsc-evaluate curve u :derivative 1)))
         (v1 (v- p-next p))
         (c1 (vlength2 v1))
         (nL (v- n (v* v1 (safe-/ 2 c1) (scalar-product v1 n))))
         (dL (v- d (v* v1 (safe-/ 2 c1) (scalar-product v1 d))))
         (v2 (v- d-next dL))
         (c2 (vlength2 v2))
         (n-next (v- nL (v* v2 (safe-/ 2 c2) (scalar-product v2 nL)))))
    (make-frame :position p-next :tangent d-next
                :binormal (cross-product d-next n-next) :normal n-next)))

(defun collect-frames (curve frame &key (resolution 100) from-end)
  (let ((first-frame (frenet-frame curve nil (if from-end 1 0))))
    (cons first-frame
          (iter (for i from 1 below resolution)
                (for u = (/ i (1- resolution)))
                (for v = (if from-end (- 1 u) u))
                (for prev first first-frame then next)
                (for next = (funcall frame curve prev v))
                (collect next)))))

(defun combine-frames (frames1 frames2)
  "As an alternative to the rotation technique, we can also combine frames from both ends
using a Hermit blend."
  (let ((n (length frames1)))
    (iter (for f1 in frames1)
          (for f2 in (reverse frames2))
          (for i from 0 to n)
          (for u = (/ i n))
          (for hermite = (+ (* -2 u u u) (* 3 u u)))
          (for new-n = (affine-combine (frame-normal f1) hermite (frame-normal f2)))
          (collect (make-frame :position (frame-position f1) :tangent (frame-tangent f1)
                               :binormal (cross-product (frame-tangent f1) new-n) :normal new-n)))))

(defun rotation-matrix (u theta)
  (let* ((m (make-array '(3 3)))
         (c (cos theta))
         (c1 (1- c))
         (s (sin theta)))
    (destructuring-bind (x y z) u
      (flet ((sqr (x) (* x x)))
        (setf (aref m 0 0) (+ c (* (sqr x) c1))
              (aref m 1 0) (+ (* y x c1) (* z s))
              (aref m 2 0) (- (* z x c1) (* y s))
              (aref m 0 1) (- (* x y c1) (* z s))
              (aref m 1 1) (+ c (* (sqr y) c1))
              (aref m 2 1) (+ (* z y c1) (* x s))
              (aref m 0 2) (+ (* x z c1) (* y s))
              (aref m 1 2) (- (* y z c1) (* x s))
              (aref m 2 2) (+ c (* (sqr z) c1)))))
    m))

(defun rotate-frame (frame angle)
  (let* ((tangent (frame-tangent frame))
         (normal (frame-normal frame))
         (new-normal (coerce
                      (matrix:to-vector
                       (matrix:multiplication (rotation-matrix tangent angle)
                                              (matrix:from-vector (coerce normal 'vector))))
                      'list)))
    (make-frame :position (frame-position frame) :tangent tangent
                :binormal (cross-product tangent new-normal) :normal new-normal)))

(defun angle-correction (curve frames)
  (let* ((end-frame (frenet-frame curve nil 1))
         (real-end-normal (frame-normal end-frame))
         (rmf-end-normal (frame-normal (car (last frames))))
         (rotation (* (acos (scalar-product real-end-normal rmf-end-normal))
                      (signum (scalar-product (cross-product (v- rmf-end-normal real-end-normal)
                                                             real-end-normal)
                                              (frame-tangent end-frame)))))
         (resolution (length frames))
         (length (bsc-estimate-arc-length curve 0 1)))
    (iter (for i from 0 below resolution)
          (for u = (/ i (1- resolution)))
          (for s = (bsc-estimate-arc-length curve 0 u))
          (for frame in frames)
          (collect (rotate-frame frame (* rotation (/ s length)))))))

(defun print-frame (stream frame by)
  (let ((p (frame-position frame))
        (n (frame-normal frame)))
    (format stream "~{~f~^ ~}~%~{~f~^ ~}~%~{~f~^ ~}~%" p (v+ p (v* n by)) p)))

(defun show-frames (curve frame &key (by (/ (bsc-bounding-box-axis curve) 10))
                                  (filename "/tmp/rmf") (resolution 100) (type 'simple))
  "BY is the length of the spikes."
  (let* ((frames1 (collect-frames curve frame :resolution resolution))
         (frames (ecase type
                   ((simple) frames1)
                   ((hermite)
                    (combine-frames
                     frames1
                     (collect-frames curve frame :resolution resolution :from-end t)))
                   ((angle) (angle-correction curve frames1)))))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (dolist (f frames)
        (print-frame s f by)))))

(defparameter *curve*
  (make-bspline-curve 3 '(0 0 0 0 0.133053 0.232419 0.338451 0.43994 0.541278 0.667394 1 1 1 1) '((0.155279 -2.14286 0) (0.0532173 -1.80779 0) (0.190782 -1.07383 1) (1.09142 -0.69642 1) (1.81621 -0.426967 2) (2.81035 -0.467704 1) (3.07214 0.640326 1) (5.1987 -0.432608 0) (4.46645 -2.13998 1) (3.25645 -1.96703 -1))))

#+nil
(show-frames *curve* #'frenet-frame :filename "/tmp/frenet")

#+nil
(show-frames *curve* #'next-frame)

#+nil
(show-frames *curve* #'next-frame :type 'hermite :filename "/tmp/rmf-hermite")

#+nil
(show-frames *curve* #'next-frame :type 'angle :filename "/tmp/rmf-angle")

; in gnuplot: splot "/tmp/rmf" with lines
