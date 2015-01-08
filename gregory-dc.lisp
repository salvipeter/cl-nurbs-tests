; Blend functions

(let* ((n 8)
       (*exponent* 2)
       (*resolution* 80)
       (points (points-from-angles (cons 5 (rest (uniform-angles n)))))
       (distance-type 'perpendicular)
       (on-off (make-array n :element-type 'boolean :initial-element nil)))
  (setf (elt on-off 0) t)
  (write-blends points (coerce on-off 'list) "/tmp/corner-blend.vtk"
		:blend-function #'corner-blend :distance-type distance-type)
  (setf (elt on-off (1- n)) t)
  (write-blends points (coerce on-off 'list) "/tmp/side-blend.vtk"
		:blend-function #'corner-blend :distance-type distance-type))


; Parameterizations

(let* ((points (points-from-angles '(5 72 72 72 72))))
  (vectorized-distance-function-test
   points '(s s s nil nil) "/tmp/radial.ps"
   :resolution 0.001d0 :density 18 :distance-type 'radial :color t)
  (vectorized-distance-function-test
   points '(nil sd nil nil nil) "/tmp/interconnected.ps"
   :resolution 0.001d0 :density 18 :distance-type 'radial-mod :color nil)
  (vectorized-distance-function-test
   points '(d s d nil nil) "/tmp/interconnected-connections.ps"
   :resolution 0.001d0 :density 18 :distance-type 'radial-mod :color t))

; Parameterization discrepancy

(defun param-discrepancy (points filename r &optional (divisor 1))
  (flet ((map-coordinates (x y) (list (/ (- x r) r) (/ (- y r) r))))
    (let* ((wh (1+ (* 2 r)))
	   (lines (lines-from-points points))
           (acc '()))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "P2~%~d ~d~%255~%" wh wh)
        (dotimes (x wh)
          (dotimes (y wh)
            (let ((p (map-coordinates x y)))
              (if (insidep lines p)
                  (let* ((di (second
                              (compute-parameter 'radial-mod 'd points p t)))
                         (si (- 1.0
                                (first
                                 (compute-parameter 'radial 's points p t))))
			 (err (abs (- si di))))
                    (push err acc)
                    (format s "~a " (round (* err 255) divisor)))
                  (format s "255 "))))
          (terpri s)))
      (reduce #'max acc))))

(let ((n 5))
  (param-discrepancy (points-from-angles (cons 5 (rest (uniform-angles n))))
                     "/tmp/discr.pgm" 100 0.5))

; Better version using VTK

(defun write-discr (points filename)
  (let* ((n (length points))
	 (*alpha* 0)
	 (vertices
          (mapcar (lambda (p)
                    (let ((d (second
                              (compute-parameter 'radial-mod 'd points p t)))
                          (s (- 1.0
                                (first
                                 (compute-parameter 'radial 's points p t)))))
                      (cons (abs (- s d)) p)))
                  (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)))

(let* ((n 5)
       (*resolution* 80)
       (points (points-from-angles (cons 5 (rest (uniform-angles n))))))
  (write-discr points "/tmp/discr.vtk"))

(defun write-discr2 (points filename)
  (let* ((n (length points))
         (acc '())
	 (*alpha* 0)
	 (vertices
          (mapcar (lambda (p)
                    (let* ((d (second
                               (compute-parameter 'radial-mod 'd points p t)))
                           (s (compute-parameter 'radial 's points p t))
                           (dis (compute-parameter 'perpendicular 'd points p t))
                           (b1 (corner-blend dis 0))
                           (b2 (corner-blend dis 1))
                           (discr (abs (- (* d (+ b1 b2))
                                          (+ (* (- 1 (first s)) b1)
                                             (* (third s) b2))))))
                      (push discr acc)
                      (cons discr p)))
                  (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)
    (reduce #'max acc)))

(let* ((n 5)
       (*resolution* 80)
       ;; (points (points-from-angles (cons 5 (rest (uniform-angles n)))))
       (points (points-from-angles '(56 30 80 74 120))))
  (write-discr2 points "/tmp/discr.vtk"))

;;; ???!!! 10^-16 a max elteres.... (ha a radial-mod a +distance-blend+ -et hasznalja)
;;; +hermite-blend+ -del is ~0.004 (5oldalu) / ~0.007 (8oldalu), tehat a max. elteres
;;; a tangens par ezrede...
;;; de ha nem regularis a poligon, akkor aranylag nagy (? pl. 0.01) elteres is lehet

(defun write-discr3 (points filename)
  (let* ((n (length points))
         (acc '())
	 (*alpha* 0)
	 (vertices
          (mapcar (lambda (p)
                    (let* ((s (second (compute-parameter 'radial 's points p t)))
                           (dis (compute-parameter 'perpendicular 'd points p t))
                           (b1 (corner-blend dis 0))
                           (b2 (corner-blend dis 1))
                           (blend (funcall +distance-blend+ s))
                           (discr (abs (- (* blend b2) (* (- 1 blend) b1)))))
                      (push discr acc)
                      (cons discr p)))
                  (vertices points))))
    (write-vtk-indexed-mesh vertices (triangles n) filename)
    (reduce #'max acc)))

(let* ((n 5)
       (*resolution* 80)
       (points (points-from-angles (cons 5 (rest (uniform-angles n))))))
  (write-discr3 points "/tmp/discr.vtk"))
