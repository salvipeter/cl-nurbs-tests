(in-packave :cl-nurbs-tests)

(defvar *p0*)
(defvar *p1*)
(defvar *p2*)
(defvar *p3*)
(defvar *normalizep*)
(defvar *reparameterization*)
(defvar *resolution*)

(defun quadratic (x) (* x x))
(defun inverse-quadratic (x) (sqrt x))

(defun generate-points (strength)
  (let ((p0 *p0*)
        (p1 (v+ *p0* (v* (v- *p1* *p0*) strength)))
        (p2 (v+ *p3* (v* (v- *p2* *p3*) strength)))
        (p3 *p3*))
    (iter (for i from 0 to *resolution*)
          (for u = (/ i *resolution*))
          (for u1 = (funcall *reparameterization* u))
          (for u2 = (funcall *reparameterization* (- 1 u)))
          (for b0 = (bernstein 3 0 u1))
          (for b1 = (bernstein 3 1 u1))
          (for b2 = (bernstein 3 1 u2))
          (for b3 = (bernstein 3 0 u2))
          (for bsum = (+ b0 b1 b2 b3))
          (for denom = (if (and *normalizep* (> bsum *epsilon*)) (/ bsum) 1))
          (for p = (v+ (v* p0 b0 denom)
                       (v* p1 b1 denom)
                       (v* p2 b2 denom)
                       (v* p3 b3 denom)))
          (collect p))))

(let* ((*resolution* 100)
       (*normalizep* t)
       ;; (*reparameterization* #'quadratic)
       (*reparameterization* #'inverse-quadratic)
       (*p0* '(0 0))
       (*p1* `(1.5 ,(* 1.5 (sqrt 3))))
       (*p2* `(1.5 ,(* 1.5 (sqrt 3))))
       (*p3* '(3 0))
       (min 0)
       (max 1)
       (step 0.2)
       (scaling 100)
       (offset '(50 300))
       (curves1 (iter (for k from min to max by step)
                      (collect (generate-points k))))
       (*reparameterization* #'identity)
       (curves2 (iter (for k from min to max by step)
                      (collect (generate-points k)))))
  (flet ((scale (p) (v+ (v* p scaling) offset)))
    (with-open-file (s "/tmp/triangle.ps" :direction :output :if-exists :supersede)
      (format s "%!PS~%newpath~%")
      (format s "~{~f ~}moveto~%~{~f ~}lineto~%" (scale *p0*) (scale *p1*))
      (format s "~{~f ~}moveto~%~{~f ~}lineto~%" (scale *p3*) (scale *p2*))
      (format s "stroke~%")
      (iter (for curve in curves1)
            (format s "newpath~%")
            (iter (for point in curve)
                  (for command first "moveto" then "lineto")
                  (format s "~{~f ~}~a~%" (scale point) command))
            (format s "stroke~%"))
      (format s "1 0 0 setrgbcolor~%")
      (iter (for curve in curves2)
            (format s "newpath~%")
            (iter (for point in curve)
                  (for command first "moveto" then "lineto")
                  (format s "~{~f ~}~a~%" (scale point) command))
            (format s "stroke~%"))
      (format s "showpage~%"))))
