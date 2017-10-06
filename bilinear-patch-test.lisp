(in-package :cl-nurbs-tests)

(defun bezier-patch-eval (cpts u v)
  "Wrapper around the B-spline patch. Only works in [0,1]."
  (let ((cpts (if (arrayp cpts)
                  cpts
                  (make-array (list (length cpts) (length (first cpts))) :initial-contents cpts))))
    (flet ((make-knots (order)
             (append (make-list order :initial-element 0)
                     (make-list order :initial-element 1))))
      (bss-evaluate (make-bspline-surface (mapcar #'1- (array-dimensions cpts))
                                          (mapcar #'make-knots (array-dimensions cpts))
                                          cpts)
                    (list u v)))))

(defun bilinear (cpts u v)
  "Extensible bilinear Bezier patch."
  (destructuring-bind ((a b) (c d))
      cpts
    (v+ (v* a (- 1 u) (- 1 v))
        (v* b (- 1 u) v)
        (v* c u (- 1 v))
        (v* d u v))))

(defun write-tensor-patch (cpts filename &key (scale 1))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (iter (for i from 0 to *resolution*)
          (for u = (* (/ i *resolution*) scale))
          (iter (for j from 0 to *resolution*)
                (for v = (* (/ j *resolution*) scale))
                (format s "v~{ ~f~}~%"
                        (if (= scale 1)
                            (bezier-patch-eval cpts u v)
                            (bilinear cpts u v)))))
    (iter (for i from 0 below *resolution*)
          (iter (for j from 0 below *resolution*)
                (for index = (+ (* i (1+ *resolution*)) j))
                (format s "f~{ ~d~}~%"
                        (list (+ index 1)
                              (+ index 2)
                              (+ index *resolution* 3)))
                (format s "f~{ ~d~}~%"
                        (list (+ index *resolution* 3)
                              (+ index *resolution* 2)
                              (+ index 1)))))))

;;; Configuration:
;; A1 A2
;; A4 A3 CT
;; A5 CC B3 B2
;; CR B5 B4 B1
;;
;;; The 2x2 points in the corner and A4,B4 are computed automatically.
;;; There are 3 patches: [A1-5,CC], [B1-5,CC], [A5,B5,CC,CR].
;;; The last (bilinear) patch is extended to double size for testing.

(defun bilinear-test (As Bs Cs)
  (destructuring-bind ((a1 a2 a3) (b1 b2 b3) (cc ct))
      (list As Bs Cs)
    (flet ((mirror (p origin) (v+ origin (v- origin p))))
      (let* ((a4 (mirror ct a3))
             (b4 (mirror ct b3))
             (a5 (mirror b3 cc))
             (b5 (mirror a3 cc))
             (cr (mirror b4 b5)))
        (assert (< (point-distance cr (mirror a4 a5)) 1.0e-5) () "Not twist compatible")
        (let ((cpts-a `((,a1 ,a2) (,a4 ,a3) (,a5 ,cc)))
              (cpts-b `((,b1 ,b2) (,b4 ,b3) (,b5 ,cc)))
              (cpts-c `((,cc ,a5) (,b5 ,cr))))
          (write-tensor-patch cpts-a "/tmp/patch-a.obj")
          (write-tensor-patch cpts-b "/tmp/patch-b.obj")
          (write-tensor-patch cpts-c "/tmp/patch-c.obj" :scale 2))))))

#+nil
(let ((*resolution* 100))
  (bilinear-test
   '((0 3 1) (1 3 1) (1 2 0))
   '((3 0 -0.5) (3 1 -0.8) (2 1 0))
   '((1 1 0) (2 2 0.2))))
