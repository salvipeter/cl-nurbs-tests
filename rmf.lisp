(in-package :cl-nurbs-tests)

(defun frenet-frame (curve prev-frame u)
  "Returns the list (position tangent binormal normal)."
  (declare (ignore prev-frame))
  (let* ((p (bsc-evaluate curve u))
         (d (bsc-evaluate curve u :derivative 1))
         (dd (bsc-evaluate curve u :derivative 2))
         (n (vnormalize (cross-product d dd))))
    (list p (vnormalize d) (vnormalize (cross-product n d)) n)))

(defun next-frame (curve prev-frame u)
  "As described in `Computation of Rotation Minimizing Frames', Wang et al., 2008.
A limitation of this method is that it is determined by the starting frame and the curve tangents,
so an end frame cannot be supplied. As a workaround, we can add a rotation gradually,
by minimizing the total (or total squared) angular speed, see section 6.3 in the paper."
  (destructuring-bind (p d b n) prev-frame
    (declare (ignore b))
    (let* ((p-next (bsc-evaluate curve u))
           (d-next (vnormalize (bsc-evaluate curve u :derivative 1)))
           (v1 (v- p-next p))
           (c1 (vlength2 v1))
           (nL (v- n (v* v1 (safe-/ 2 c1) (scalar-product v1 n))))
           (dL (v- d (v* v1 (safe-/ 2 c1) (scalar-product v1 d))))
           (v2 (v- d-next dL))
           (c2 (vlength2 v2))
           (n-next (v- nL (v* v2 (safe-/ 2 c2) (scalar-product v2 nL)))))
      (list p-next d-next (cross-product d-next n-next) n-next))))

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
          (for new-n = (affine-combine (fourth f1) hermite (fourth f2)))
          (collect (list (first f1) (second f1) (cross-product (second f1) new-n) new-n)))))

(defun print-frame (stream frame by)
  (let ((p (first frame))
        (n (fourth frame)))
    (format stream "~{~f~^ ~}~%~{~f~^ ~}~%~{~f~^ ~}~%" p (v+ p (v* n by)) p)))

(defun show-frames (curve frame &key (by (/ (bsc-bounding-box-axis curve) 10))
                                  (filename "/tmp/rmf") (resolution 100) hermitep)
  "BY is the length of the spikes."
  (let* ((frames1 (collect-frames curve frame :resolution resolution))
         (frames2 (collect-frames curve frame :resolution resolution :from-end t))
         (frames (if hermitep (combine-frames frames1 frames2) frames1)))
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
(show-frames *curve* #'next-frame :hermitep t :filename "/tmp/rmf-hermite")

; in gnuplot: splot "/tmp/rmf" with lines
