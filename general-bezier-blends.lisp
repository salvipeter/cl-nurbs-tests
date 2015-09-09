(in-package :cl-nurbs-tests)

(defun generalized-bernstein (points p side degree col row &key (use-d t))
  (let* ((n (length points))
         (l (barycentric-coordinates points p))
         (half-low (floor degree 2))
         (i side)
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (si (barycentric-s l i))
         (si-1 (barycentric-s l i-1))
         (si+1 (barycentric-s l i+1))
         (di (barycentric-d l i))
         (di-1 (barycentric-d l i-1))
         (di+1 (barycentric-d l i+1))
         (alpha (if use-d
                    (if (< (+ di-1 di) *epsilon*)
                        0.5
                        (/ di-1 (+ di-1 di)))
                    (if (< (+ si (- 1 si-1)) *epsilon*)
                        0.5
                        (/ si (+ si (- 1 si-1))))))
         (beta (if use-d
                   (if (< (+ di+1 di) *epsilon*)
                       0.5
                       (/ di+1 (+ di+1 di)))
                   (if (< (+ (- 1 si) si+1) *epsilon*)
                       0.5
                       (/ (- 1 si) (+ (- 1 si) si+1)))))
         (blend (* (bernstein degree row di)
                   (bernstein degree col si)))
         (mu (cond ((and (evenp degree) (= col half-low)) 1)
                   ((> col half-low) beta)
                   (t alpha))))
    (* blend mu)))

;;; DEFICIENCY => see deficiency.lisp
(defun write-bernstein-blend-mesh (path n degree &key (use-d t))
  (let ((points (points-from-angles (uniform-angles n))))
    (iter (for row from 0 below (ceiling degree 2))
          (iter (for col from 0 to degree)
                (for fname = (format nil "~a/~asided-deg~a-c~a~a.obj" path n degree col row))
                (with-open-file (s fname :direction :output :if-exists :supersede)
                  (write-obj-indexed-mesh
                   (iter (for p in (vertices points))
                         (for b = (generalized-bernstein points p 0 degree col row :use-d use-d))
                         (collect (cons b p)))
                   (triangles n)
                   fname))))
    (let ((fname (format nil "~a/~asided-deg~a-center.obj" path n degree)))
      (with-open-file (s fname :direction :output :if-exists :supersede)
        (write-obj-indexed-mesh
         (iter (for p in (vertices points))
               (for def = (deficiency n degree :position p :use-d use-d))
               (collect (cons def p)))
         (triangles n)
         fname)))))

#+nil
(let ((*resolution* 30))
  (iter (for n from 3 to 8)
        (iter (for d from 1 to 7)
              (write-bernstein-blend-mesh "/tmp" n d))))

(defun write-bernstein-blend-image (path n degree &key (use-d t) (density 0.1))
  (labels ((transform (p) (list (+ (* (first p) 250) 250) (- 500 (* (second p) 250))))
           (write-poly (stream points msg)
             (format stream "~{~f ~}moveto~%~{~{~f ~}lineto~%~}closepath stroke~%"
                     (transform (first points)) (mapcar #'transform (rest points)))
             (format stream "10 790 moveto (~a, interval between lines: ~a) show showpage~%"
                     msg density)))
    (let ((fname (format nil "~a/~asided-deg~a.ps" path n degree)))
      (with-open-file (s fname :direction :output :if-exists :supersede)
        (format s "%!PS~%/Times-Roman findfont 15 scalefont setfont~%")
        (let ((points (points-from-angles (uniform-angles n))))
          (iter (for row from 0 below (ceiling degree 2))
                (iter (for col from 0 to degree)
                      (for cp-str = (format nil "Blend of control point (~a, ~a)" col row))
                      (write-ps-indexed-mesh-projection
                       (iter (for p in (vertices points))
                             (for b = (generalized-bernstein points p 0 degree col row :use-d use-d))
                             (collect (cons b p)))
                       (triangles n) s :transform #'transform :axis 0 :lines density)
                      (write-poly s points cp-str)))
          (let ((cp-str "Blend of central control point"))
            (write-ps-indexed-mesh-projection
             (iter (for p in (vertices points))
                   (for def = (deficiency n degree :position p :use-d use-d))
                   (collect (cons (if (< (abs def) *epsilon*) 0 def) p)))
             (triangles n) s :transform #'transform :axis 0 :lines density)
            (write-poly s points cp-str)))))))

#+nil
(let ((*resolution* 60))
  (iter (for n from 3 to 8)
        (iter (for d from 1 to 7)
              (write-bernstein-blend-image "/tmp" n d :density 0.05))))
