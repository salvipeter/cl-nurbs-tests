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

(defparameter *displacement*
  '(((0 0.1 0.3) (0.2 0.2 0.2))
    ((0 0 0) (0 0 0))
    ((0 0 1) (0 0 0))
    ((0 0 0) (0 0 0)))
  "For all sides a list of (OUTER INNER).")

(defun elevated-and-displaced (surface displacement)
  "Elevates a bicubic Bezier surface
and adds displacements to the central control points on all sides."
  (let ((s (bezier-surface-elevate-uv surface)))
    (iter (for (out in) in displacement)
          (for ((oi oj) (ii ij))
               in '(((0 2) (1 2))
                    ((2 0) (2 1))
                    ((4 2) (4 3))
                    ((2 4) (2 3))))
          (setf (aref s oi oj)
                (v+ (aref s oi oj) out))
          (setf (aref s ii ij)
                (v+ (aref s ii ij) in)))
    s))

#+nil
(write-points2-vtk
 (elevated-and-displaced *surface* *displacement*)
 "/tmp/displaced.vtk")

#+nil
(let ((*resolution* 100))
  (write-points2-vtk
   (bezier-surface-sample
    (elevated-and-displaced *surface* *displacement*))
   "/tmp/edsurface.vtk"))

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
                      (for Bs =
                           ; (bernstein 4 2 (compute-nr-parameter l k 's))
                           (* 6
                              (expt (elt l (mod (1- k) (length points))) 2)
                              (expt (elt l k) 2))
                           )
                      (setf (aref result i j)
                            (v+ (aref result i j)
                                (v* out Bs
                                    (bernstein 4 0 (compute-nr-parameter l k 'd)))
                                (v* in Bs
                                    (bernstein 4 1 (compute-nr-parameter l k 'd))))))))
    result))

#+nil
(let ((*resolution* 100)
      (*barycentric-type* 'wachspress)
      (*barycentric-normalized* t)
      (*use-local-d* nil)
      ; (points '((0 0) (1 0) (1 1) (0 1)))
      (points '((0 0.3) (0.8 0) (1 0.9) (0.1 0.1))))
  (write-points2-vtk
   (displacement-surface-sample points *surface* *displacement*)
   "/tmp/dsurface.vtk"))


;;; Blend function plots
(let ((*barycentric-normalized* t)
      (*barycentric-type* 'wachspress)
      (*use-local-d* nil)
      (*points* (points-from-angles '(20 70 30 70 30))))
  (flet ((fun (i)
           (lambda (points p)
             (let ((l (barycentric-coordinates points p)))
               (list
                ;; All blends of one side
                #-nil
                (let* ((i-1 (mod (1- i) (length points)))
                       (i+1 (mod (1+ i) (length points)))
                       (li (elt l i))
                       (li-1 (elt l i-1))
                       (s (safe-/ li (+ li li-1)))
                       (d (compute-nr-parameter l i 'd))
                       (d-1 (compute-nr-parameter l i-1 'd))
                       (d+1 (- 1 (compute-nr-parameter l i+1 'd))))
                  (+ (* (+ (bernstein 4 0 d-1)
                           (bernstein 4 1 d-1)
                           (bernstein 4 2 s)
                           (bernstein 4 3 d+1)
                           (bernstein 4 4 d+1))
                        (bernstein 4 1 d))
                     (* (+ (bernstein 4 0 d-1)
                           (bernstein 4 1 d-1)
                           (bernstein 4 2 s)
                           (bernstein 4 3 d+1)
                           (bernstein 4 4 d+1))
                        (bernstein 4 0 d))))
                ;; Sum of all coefficients
                #+nil
                (iter (for i from 0 below (length *points*))
                      (for i-1 = (mod (1- i) (length points)))
                      (for i+1 = (mod (1+ i) (length points)))
                      (for li = (elt l i))
                      (for li-1 = (elt l i-1))
                      (for s = (safe-/ li (+ li li-1)))
                      (for d = (compute-nr-parameter l i 'd))
                      (for d-1 = (compute-nr-parameter l i-1 'd))
                      (for d+1 = (- 1 (compute-nr-parameter l i+1 'd)))
                      (sum
                       (+ (* (+ (bernstein 4 1 d-1)
                                (bernstein 4 2 s))
                             (bernstein 4 1 d))
                          (* (+ (bernstein 4 0 d-1)
                                (bernstein 4 1 d-1)
                                (bernstein 4 2 s)
                                (bernstein 4 3 d+1))
                             (bernstein 4 0 d)))))
                ;; Sum of all central displacement coefficients
                #+nil
                (iter (for i from 0 below (length *points*))
                      (for li = (elt l i))
                      (for li-1 = (elt l (mod (1- i) (length points))))
                      (for li* = (safe-/ li (+ li li-1)))
                      (for li-1* = (safe-/ li-1 (+ li li-1)))
                      (sum
                       (* (bernstein 4 2 li*)
                          (+ (bernstein 4 0 (compute-nr-parameter l i 'd))
                             (bernstein 4 1 (compute-nr-parameter l i 'd))))))
                ;; Barycentric version
                #+nil
                (let* ((li (elt l i))
                       (li-1 (elt l (mod (1- i) (length points))))
                       (li* (safe-/ li (+ li li-1)))
                       (li-1* (safe-/ li-1 (+ li li-1))))
                  (* 6 (expt li* 2) (expt li-1* 2)
                     (+ (bernstein 4 0 (compute-nr-parameter l i 'd))
                        (bernstein 4 1 (compute-nr-parameter l i 'd)))))
                ;; Bernstein version
                #+nil
                (* (bernstein 4 2 (/ (1+ (first p)) 2))
                   (+ (bernstein 4 0 (/ (1+ (second p)) 2))
                      (bernstein 4 1 (/ (1+ (second p)) 2)))))))))
    (let ((i 1))
      (bitmap-test *points* (fun i) "/tmp/plot.pgm"
                   :object-size 2.0d0 :density 20))))

;;; => s_i = l_i / (l_i-1 + l_i), d_i = 1 - l_i-1 - l_i


;;; Blend function 3D test

(defun barycentric-blend (l i type)
  (let* ((n (length l))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (li (elt l i))
         (li-1 (elt l i-1))
         (s (safe-/ li (+ li li-1)))
         (d (compute-nr-parameter l i 'd))
         (d-1 (compute-nr-parameter l i-1 'd))
         (d+1 (- 1 (compute-nr-parameter l i+1 'd))))
    (declare (ignore d+1))
    (case type
      (side (* (+ (bernstein 3 0 s)
                  (bernstein 3 1 s)
                  (bernstein 3 2 s)
                  (bernstein 3 3 s))
               (+ (bernstein 3 0 d)
                  (bernstein 3 1 d))))
      (side-half (* (+ (bernstein 3 0 s)
                       (bernstein 3 1 s)
                       (bernstein 3 2 s)
                       (bernstein 3 3 s))
                    (+ (bernstein 3 0 d)
                       (bernstein 3 1 d))
                    1/2))
      (corner (* (+ (bernstein 3 0 d-1)
                    (bernstein 3 1 d-1))
                 (+ (bernstein 3 0 d)
                    (bernstein 3 1 d))))
      (deviation (- (* (+ (bernstein 3 0 s)
                          (bernstein 3 1 s)
                          (bernstein 3 2 s)
                          (bernstein 3 3 s))
                       (+ (bernstein 3 0 d)
                          (bernstein 3 1 d)))
                    (* (+ (bernstein 3 0 d-1)
                          (bernstein 3 1 d-1))
                       (+ (bernstein 3 0 d)
                          (bernstein 3 1 d))))))))

(defun write-barycentric-blends (points on-off filename &key (type 'corner))
  (let* ((n (length points))
	 (vertices
          (mapcar (lambda (p)
                    (let ((l (barycentric-coordinates points p)))
                      (cons (iter (for i from 0 below (length on-off))
                                  (when (elt on-off i)
                                    (sum (funcall #'barycentric-blend l i type))))
                            p)))
                  (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

#+nil
(let ((points (points-from-angles '(20 70 30 70 30)))
      (*resolution* 100)
      (*use-local-d* nil))
  (write-barycentric-blends points '(t nil nil nil t) "/tmp/side0.vtk")
  (write-barycentric-blends points '(t t nil nil nil) "/tmp/side1.vtk")
  (write-barycentric-blends points '(t nil nil nil nil) "/tmp/corner0.vtk")
  (write-barycentric-blends points '(nil t nil nil nil) "/tmp/corner1.vtk")
  (write-barycentric-blends points '(t t t t t) "/tmp/all.vtk"))

#+nil
(let ((points (points-from-angles '(20 70 30 70 30)))
      (*use-local-d* nil))
  (write-barycentric-blends points '(t nil nil nil nil) "/tmp/side1b.vtk"
                            :type 'side)
  (write-barycentric-blends points '(nil t nil nil nil) "/tmp/side2b.vtk"
                            :type 'side)
  (write-barycentric-blends points '(nil nil t nil nil) "/tmp/side3b.vtk"
                            :type 'side)
  (write-barycentric-blends points '(nil nil nil t nil) "/tmp/side4b.vtk"
                            :type 'side)
  (write-barycentric-blends points '(nil nil nil nil t) "/tmp/side0b.vtk"
                            :type 'side)
  (write-barycentric-blends points '(t t t t t) "/tmp/allb.vtk"
                            :type 'deviation)
  (write-barycentric-blends points '(t t t t t) "/tmp/allc.vtk"
                            :type 'side-half))
