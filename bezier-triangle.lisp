(in-package :cl-nurbs-tests)

;;; Create a quadratic Bezier triangle between two nx1-degree Bezier patches meeting at a corner
;;; Twist compatibility should be ensured first, using a user-supplied twist control point

(defun ensure-twist-compatibility (a c twist)
  "A and C are two nx1-degree Bezier patches, meeting at a concave angle at their (0,0) point.
TWIST is a twist control point given on the convex side."
  (flet ((mirror (p center) (v+ center (v- center p))))
    (macrolet ((pa (i j) `(elt (elt a ,i) ,j))
               (pc (i j) `(elt (elt c ,i) ,j)))
      (setf (pa 1 1) (mirror twist (pa 1 0))
            (pc 1 1) (mirror twist (pc 1 0))
            (pa 0 1) (mirror (pc 1 0) (pa 0 0))
            (pc 0 1) (mirror (pa 1 0) (pc 0 0))))))

(defun bezier-biquadratic-end (surface)
  "SURFACE is a nx1-degree Bezier patch.
Returns the first two control rows of its biquadratic conversion (2x3 control points)."
  (flet ((p (i j) (elt (elt surface i) j)))
    (let* ((n (1- (length surface)))
           (a00 (p 0 0))
           (a02 (p 0 1))
           (a01 (affine-combine a00 1/2 a02))
           (a10 (affine-combine (p 0 0) 1/2 (p 1 0)))
           (a12 (affine-combine (p 0 1) 1/2 (p 1 1)))
           (a11 (affine-combine a10 1/2 a12)))
      (list (list a00 a01 a02)
            (list a10 a11 a12)))))

(defun bezier-triangle-filler (a c)
  "A and C are two twist-compatible nx1-degree Bezier patches,
meeting at a concave angle at their (0,0) point.
The Bezier triangle is stored in the following format:
  ((B200 B020 B002) (B011 B101 B110))"
  (let ((a2 (bezier-biquadratic-end a))
        (c2 (bezier-biquadratic-end c)))
    (flet ((pa (i j) (elt (elt a2 i) j))
           (pc (i j) (elt (elt c2 i) j)))
      (let ((b200 (pc 0 2))
            (b020 (pa 0 0))
            (b002 (pa 0 2))
            (b011 (pa 0 1))
            (b101 (v+ (pa 0 1)
                      (v- (pa 0 1) (pa 1 1))
                      (v- (pa 0 1) (pa 0 0))
                      (v- (pa 1 0) (pa 1 1))))
            (b110 (pc 0 1)))
        (list (list b200 b020 b002)
              (list b011 b101 b110))))))

(defun quadratic-bezier-triangle-eval (tri uvw)
  "The Bezier triangle is stored in the following format:
  ((B200 B020 B002) (B011 B101 B110))"
  (let ((result '(0 0 0)))
    (iter (for i from 0 below 3)
          (setf result
                (v+ result
                    (v* (elt (first tri) i)
                        (expt (elt uvw i) 2)))))
    (iter (for i from 0 below 3)
          (setf result
                (v+ result
                    (v* (elt (second tri) i)
                        2 (reduce #'* (append (subseq uvw 0 i) (subseq uvw (1+ i))))))))
    result))

(defun sample-quadratic-bezier-triangle (tri resolution length)
  "Samples the triangle inside as well as in the negative region of the middle parameter.
The LENGTH parameter controls the maximum value of the parameters."
  (iter (for i from resolution downto 0)
        (for u = (* (/ i resolution) length))
        (for last-row first nil then row)
        (for row = (iter (for j from 0 to (- resolution i))
                         (for w = (* (/ j resolution) length))
                         (for v = (- 1 u w))
                         (collect (quadratic-bezier-triangle-eval tri (list u v w)))))
        (when last-row
          (appending (iter (for j from 0 below (- resolution i))
                           (collect (list (elt last-row j)
                                          (elt row j)
                                          (elt row (1+ j))))
                           (unless (= j 0)
                             (collect (list (elt last-row (1- j))
                                            (elt row j)
                                            (elt last-row j)))))))))

(defun bilinear-filler (a c)
  (flet ((pa (i j) (elt (elt a i) j))
         (pc (i j) (elt (elt c i) j)))
    (let ((q (v+ (pa 0 1) (v- (pa 0 1) (pa 1 1)))))
      (list (list (pa 0 0) (pa 0 1))
            (list (pc 0 1) q)))))

(defun bilinear-paralelogram (a c)
  (flet ((pa (i j) (elt (elt a i) j))
         (pc (i j) (elt (elt c i) j)))
    (let ((q (v+ (pa 0 1) (v- (pc 0 1) (pc 0 0)))))
      (list (list (pa 0 0) (pa 0 1))
            (list (pc 0 1) q)))))

(defun write-bezier-ribbon-control-points (ribbons filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "~{~{~{v~{ ~f~}~%~}~}~}" ribbons)))

(defun bezier-surface (cpts uv)
  (iter (with p = '(0 0 0))
        (with m = (1- (length cpts)))
        (with n = (1- (length (first cpts))))
        (for j from 0 to m)
        (iter (for i from 0 to n)
              (for cp = (elt (elt cpts j) i))
              (setf p (v+ p (v* cp (bernstein n i (first uv)) (bernstein m j (second uv))))))
        (finally (return p))))

(defun sample-bezier-surface (surface resolution from to)
  "FROM and TO are 2D points. Returns a list of triangles."
  (iter (with length = (v- to from))
        (for i from 0 to resolution)
        (for u = (+ (first from) (* (/ i resolution) (first length))))
        (for last-row first nil then row)
        (for row = (iter (for j from 0 to resolution)
                         (for v = (+ (second from) (* (/ j resolution) (second length))))
                         (collect (bezier-surface surface (list u v)))))
        (when last-row
          (appending (iter (for j from 0 below resolution)
                           (collect (list (elt last-row j)
                                          (elt row j)
                                          (elt row (1+ j))))
                           (collect (list (elt row (1+ j))
                                          (elt last-row (1+ j))
                                          (elt last-row j))))))))

#+nil
(let ((a '(((6 3 0) (4 0 1))
           ((4 4 1) (2 2 2))
           ((2 6 2) (1 4 2))))
      (c '(((6 3 0) (7 1 0))
           ((8 5 0) (9 2 1))
           ((10 6 1) (12 4 2))
           ((13 6 3) (14 4 3)))))
  (write-bezier-ribbon-control-points (list a c) "/tmp/points-original.obj")
  (ensure-twist-compatibility a c '(6 6 0))
  (let ((b (bezier-triangle-filler a c)))
    (write-bezier-ribbon-control-points (list a b c) "/tmp/points.obj")
    (write-stl (append (sample-bezier-surface a *resolution* '(0 0) '(2 1))
                       (sample-quadratic-bezier-triangle b *resolution* 2)
                       (sample-bezier-surface c *resolution* '(0 0) '(2 1)))
               "/tmp/surfaces.stl")))

#+nil
(let ((a '(((6 3 0) (4 0 1))
           ((4 4 1) (2 2 2))
           ((2 6 2) (1 4 2))))
      (c '(((6 3 0) (7 1 0))
           ((8 5 0) (9 2 1))
           ((10 6 1) (12 4 2))
           ((13 6 3) (14 4 3)))))
  (write-bezier-ribbon-control-points (list a c) "/tmp/points-original.obj")
  (ensure-twist-compatibility a c '(6 6 0))
  (let ((b (bilinear-filler a c)))      ; (bilinear-paralelogram a c)
    (write-bezier-ribbon-control-points (list a b c) "/tmp/points.obj")
    (write-stl (append (sample-bezier-surface a *resolution* '(0 0) '(2 1))
                       (sample-bezier-surface b *resolution* '(0 0) '(2 2))
                       (sample-bezier-surface c *resolution* '(0 0) '(2 1)))
               "/tmp/surfaces.stl")))
