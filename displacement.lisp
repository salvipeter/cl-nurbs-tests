(in-package :cl-nurbs-tests)

(defun bezier-surface-evaluate (points uv)
  (let ((n (1- (array-dimension points 0)))
        (m (1- (array-dimension points 1)))
        (p '(0 0 0)))
    (iter (for i from 0 to n)
          (for Bu = (bernstein n i (first uv)))
          (iter (for j from 0 to m)
                (for Bv = (bernstein m j (second uv)))
                (setf p (v+ p (v* (aref points i j) Bu Bv)))))
    p))

(defun bezier-surface-elevate-u (points)
  (let* ((n (1- (array-dimension points 0)))
         (m (1- (array-dimension points 1)))
         (result (make-array (list (+ n 2) (1+ m)))))
    (iter (for j from 0 to m)
          (setf (aref result 0 j) (aref points 0 j))
          (iter (for i from 1 to n)
                (setf (aref result i j)
                      (affine-combine (aref points i j)
                                      (/ i (1+ n))
                                      (aref points (1- i) j))))
          (setf (aref result (1+ n) j) (aref points n j)))
    result))

(defun bezier-surface-elevate-v (points)
  (let* ((n (1- (array-dimension points 0)))
         (m (1- (array-dimension points 1)))
         (result (make-array (list (1+ n) (+ m 2)))))
    (iter (for i from 0 to n)
          (setf (aref result i 0) (aref points i 0))
          (iter (for j from 1 to m)
                (setf (aref result i j)
                      (affine-combine (aref points i j)
                                      (/ j (1+ m))
                                      (aref points i (1- j)))))
          (setf (aref result i (1+ m)) (aref points i m)))
    result))

(defun bezier-surface-elevate-uv (points)
  (bezier-surface-elevate-u
   (bezier-surface-elevate-v points)))

(defun bezier-surface-sample (points)
  (let ((result (make-array (list *resolution* *resolution*))))
    (iter (for i from 0 below *resolution*)
          (for u = (/ i (1- *resolution*)))
          (iter (for j from 0 below *resolution*)
                (for v = (/ j (1- *resolution*)))
                (setf (aref result i j)
                      (bezier-surface-evaluate points (list u v)))))
    result))

(defparameter *surface*
  (make-array '(4 4) :initial-contents
              '(((0 0 0) (2 0 1) (4 0 1) (6 0 0))
                ((0 2 1) (2 2 2) (4 2 2) (6 2 1))
                ((0 4 1) (2 4 2) (4 4 2) (6 4 1))
                ((0 6 0) (2 6 1) (4 6 1) (6 6 0)))))

#+nil
(write-points2-vtk *surface* "/tmp/control.vtk")

#+nil
(write-points2-vtk (bezier-surface-elevate-uv *surface*) "/tmp/elevated.vtk")

#+nil
(let ((*resolution* 100))
  (write-points2-vtk (bezier-surface-sample *surface*) "/tmp/surface.vtk"))

;;; see BARYCENTRIC-COORDINATES and COMPUTE-NR-PARAMETER in n-sided-ribbon.lisp
(defun displacement-surface-sample (points base displacement &key only)
  (let ((result (make-array (list *resolution* *resolution*)))
        (sides (if only only (loop for i from 0 below (length points) collect i))))
    (iter (for i from 0 below *resolution*)
          (for u = (/ i (1- *resolution*)))
          (iter (for j from 0 below *resolution*)
                (for v = (/ j (1- *resolution*)))
                (setf (aref result i j)
                      (bezier-surface-evaluate base (list u v)))
                (iter (for k in sides)
                      (for (out in) = (elt displacement k))
                      (for l = (barycentric-coordinates points (list u v)))
                      (setf (aref result i j)
                            (v+ (aref result i j)
                                (v* out
                                    (bernstein 4 2 (compute-nr-parameter l k 's))
                                    (bernstein 4 0 (compute-nr-parameter l k 'd)))
                                (v* in
                                    (bernstein 4 2 (compute-nr-parameter l k 's))
                                    (bernstein 4 1 (compute-nr-parameter l k 'd))))))))
    result))

(defparameter *displacement*
  '(((0 0.1 0.3) (0.2 0.2 0.2))
    ((0 0 0) (0 0 0))
    ((0 0 1) (0 0 0))
    ((0 0 0) (0 0 0)))
  "For all sides a list of (OUTER INNER).")

#+nil
(let ((*resolution* 100)
      (*barycentric-type* 'wachspress)
      (*barycentric-normalized* t)
      (*use-local-d* nil)
      (points '((0 0) (1 0) (1 1) (0 1))))
  (write-points2-vtk (displacement-surface-sample points *surface* *displacement* :only '(2))
                     "/tmp/dsurface.vtk"))
