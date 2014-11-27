; Blend functions

(let* ((n 6)
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
