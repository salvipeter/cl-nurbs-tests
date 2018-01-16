(in-package :cl-nurbs-tests)

;;; Question:
;;; Blending G2 curves can result in a non-G2 curve? (Answer: yes)

(defun parabola-blend (alpha)
  (let ((a1 '((0 0) (2 4) (6 4)))
        (a2 '((6 4) (8 4) (12 3)))
        (b '((0 0) (6 5) (12 3))))
    (flet ((curve (u)
             (let ((au (if (< u 1/2) (bezier a1 (* 2 u)) (bezier a2 (1- (* 2 u)))))
                   (bu (bezier b u)))
               (affine-combine au alpha bu))))
      (iter (for i from 0 to *resolution*)
            (for u = (/ i *resolution*))
            (collect (curve u))))))

(defun write-extruded-surface (points filename &optional (length 10))
  (let ((n (length points)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (iter (for p in points)
            (iter (for i from 0 below n)
                  (for u = (/ i (1- n)))
                  (format s "v~{ ~f~} ~f~%" p (* u length))))
      (iter (for i from 1 below n)
            (iter (for j from 1 below n)
                  (flet ((p (i j) (+ (* j n) i 1)))
                    (format s "f ~a ~a ~a ~a~%"
                            (p (1- i) (1- j)) (p (1- i) j)
                            (p i j) (p i (1- j)))))))))

#+nil
(let ((*resolution* 100))
  (write-extruded-surface (parabola-blend 0) "/tmp/a.obj")
  (write-extruded-surface (parabola-blend 1) "/tmp/b.obj")
  (write-extruded-surface (parabola-blend 1/2) "/tmp/c.obj"))
