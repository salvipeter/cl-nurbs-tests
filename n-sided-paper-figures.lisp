;;; Parameters for Figures 1-3
(defparameter *p0* '(0.3 0.5))
(defparameter *t0* '(0.2 0.5))
(defparameter *p1* '(0.2 0.5))
(defparameter *t1* '(-1.0 1.0))


;;; Figure 1

;;; Hermite curve, showing curvature comb and point/tangent constraints
;;; (tangents are shown in 1/3 length)

(let ((*resolution* 200)
      (*curvature-comb-scale* 0.01))
  (write-hermite *p0* *t0* *p1* *t1* "n-sided-paper/01-hermite.ps"))


;;; Figure 2

;;; Ribbon curve with Hermite blend function, showing curvature comb and ribbons
;;; (ribbons are shown in 2/3 length, because the ribbon multipliers are set to 0.5)

(let ((*resolution* 200)
      (*curvature-comb-scale* 0.01)
      (*ribbon-multiplier-start* 0.5)
      (*ribbon-multiplier-end* 0.5))
  (write-tomi-hybrid *p0* *t0* *p1* *t1* "n-sided-paper/02-ribbon-hermite.ps"))


;;; Figure 3

;;; Ribbon curve with distance-based blend function, visualization same as before
;;; (ribbons are shown in 2/3 length as before)

(let ((*resolution* 200)
      (*curvature-comb-scale* 0.01)
      (*ribbon-multiplier-start* 0.5)
      (*ribbon-multiplier-end* 0.5)
      (*exponent* 2))
  (write-tomi *p0* *t0* *p1* *t1* "n-sided-paper/03-ribbon.ps"))

;;; proba
(let ((*resolution* 200)
      (*curvature-comb-scale* 0.01)
      (*ribbon-multiplier-start* 0.5)
      (*ribbon-multiplier-end* 0.5)
      (*exponent* 2)
      (filename "n-sided-paper/proba.ps")
      (type #'hermite)
      (color *green*))
  (let (curve derivative)
    (iter (for i from 0 below *resolution*)
	  (for u = (/ i (1- *resolution*)))
	  (with curve-fn = (lambda (x d) (funcall type *p0* *t0* *p1* *t1* x :derivative d)))
	  (for p = (funcall curve-fn u 0))
	  (collect p into tcurve)
	  (collect (funcall curve-fn u 1) into tderivative)
	  (finally (setf curve tcurve derivative tderivative)))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "%!PS~%")
      ;; Curve
      (polyline s curve :window 3)
      (format s "0.1 setlinewidth~%")
      (iter (for start in curve)
	    (for vec in derivative)
	    (for end = (v+ start (v* vec 1/10)))
	    (for i upfrom 0)
	    (when (zerop (mod i 20))
	      (polyline s (list start end) :window 3 :color color)))
      (format s "showpage~%"))))


;;; Figure 6

(defparameter *coords*
  '((((-3 -25.8 0) (-1 -25.8 0) (1 -25.8 0) (3 -25.8 0))
     ((3 -25.8 0) (12 -17.2 0) (21 -8.6 0) (30 0 0))
     ((30 0 0) (25 8.6 0) (20 17.2 0) (15 25.8 0))
     ((15 25.8 0) (5 25.8 0) (-5 25.8 0) (-15 25.8 0))
     ((-15 25.8 0) (-20 17.2 0) (-25 8.6 0) (-30 0 0))
     ((-30 0 0) (-21 -8.6 0) (-12 -17.2 0) (-3 -25.8 0)))
    (((-4 -17.2 0) (4 -17.2 0))
     ((4 -17.2 0) (18 0 0))
     ((18 0 0) (10 17.2 0))
     ((10 17.2 0) (-10 17.2 0))
     ((-10 17.2 0) (-18 0 0))
     ((-18 0 0) (-4 -17.2 0)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((3 -25.8 0) (30 0 0) (15 25.8 0) (-15 25.8 0) (-30 0 0) (-3 -25.8 0)))))

;;; (write-constraint-grid *points* "n-sided-paper/06x-grid.vtk" :coords *coords*)
#+nil
(write-constraint-ribbons *points* "n-sided-paper/06x-ribbons.vtk" :coords *coords*
			  :resolution 20)

#+nil
(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 2.5d0))
  (write-patch (points-from-angles (uniform-angles 6)) 'ribbon "n-sided-paper/06a-regular.vtk"
	       :coords *coords* :distance-type 'perpendicular :spider t))

#+nil
(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 2.5d0))
  (write-patch *points* 'ribbon "n-sided-paper/06b-irregular.vtk" :coords *coords*
	       :distance-type 'perpendicular :spider t))

(defparameter *coords*
  '((((0.0d0 0.0d0 0.0d0)
      (10.0d0 0.0d0 1.0d0)
      (20.0d0 0.0d0 1.0d0)
      (30.0d0 0.0d0 0.0d0))
     ((30.0d0 0.0d0 0.0d0)
      (30.0d0 2.0d0 1.0d0)
      (30.0d0 4.0d0 1.0d0)
      (30.0d0 6.0d0 0.0d0))
     ((30.0d0 6.0d0 0.0d0)
      (20.0d0 6.0d0 1.0d0)
      (10.0d0 6.0d0 1.0d0)
      (0.0d0 6.0d0 0.0d0))
     ((0.0d0 6.0d0 0.0d0)
      (0.0d0 4.0d0 1.0d0)
      (0.0d0 2.0d0 1.0d0)
      (0.0d0 0.0d0 0.0d0)))
    (((10.0d0 2.0d0 1.2d0)
      (20.0d0 2.0d0 1.2d0))
     ((20.0d0 2.0d0 1.2d0)
      (20.0d0 4.0d0 1.2d0))
     ((20.0d0 4.0d0 1.2d0)
      (10.0d0 4.0d0 1.2d0))
     ((10.0d0 4.0d0 1.2d0)
      (10.0d0 2.0d0 1.2d0)))))

(defparameter *coords*
  '((((3 0 0)
      (2 1 1)
      (1 2 1)
      (0 3 0))
     ((0 3 0)
      (-5 2 1)
      (-10 1 1)
      (-15 0 0))
     ((-15 0 0)
      (-10 -1 1)
      (-5 -2 1)
      (0 -3 0))
     ((0 -3 0)
      (1 -2 1)
      (2 -1 1)
      (3 0 0)))
    (((1 0 1.1)
      (0 1 1.1))
     ((0 1 1.1)
      (-10 0 1.1))
     ((-10 0 1.1)
      (0 -1 1.1))
     ((0 -1 1.1)
      (1 0 1.1)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((0 3) (-15 0) (0 -3) (3 0)))))

;;; Sketches rombusz:
(defparameter *coords*
  '((((174.683d0 -79.2657d0 39.9866d0)
      (147.39d0 -69.055d0 34.5583d0)
      (-92.2374d0 17.9945d0 -13.5014d0)
      (-146.776d0 37.1868d0 -22.987d0))
     ((-146.776d0 37.1868d0 -22.987d0)
      (-89.7339d0 42.4423d0 -29.727d0)
      (171.56d0 64.1718d0 -60.5931d0)
      (200.335d0 67.0522d0 -63.9375d0))
     ((200.335d0 67.0522d0 -63.9375d0)
      (210.949d0 44.1498d0 -48.4583d0)
      (225.487d0 13.4434d0 -27.1227d0)
      (237.86d0 -12.5697d0 -9.64074d0))
     ((237.86d0 -12.5697d0 -9.64074d0)
      (217.464d0 -34.1532d0 6.21314d0)
      (192.862d0 -60.3934d0 26.2093d0)
      (174.683d0 -79.2657d0 39.9866d0)))
    (((164.961d0 -50.4667d0 20.407d0)
      (-39.6534d0 23.6076d0 -19.1114d0))
     ((-39.6534d0 23.6076d0 -19.1114d0)
      (181.71d0 41.4307d0 -44.9136d0))
     ((181.71d0 41.4307d0 -44.9136d0)
      (205.954d0 -8.0374d0 -10.5976d0))
     ((205.954d0 -8.0374d0 -10.5976d0)
      (164.961d0 -50.4667d0 20.407d0)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((-146.776 37.1868 -22.987)
      (200.335 67.0522 -63.9375)
      (237.86 -12.5697 -9.64074)
      (174.683 -79.2657 39.9866)))))

;;; levagott negyzet:
(defparameter *coords*
  '((((-0.5d0 -6.0d0 0.0d0)
      (-0.16666666666666666d0 -6.0d0 0.0d0)
      (0.16666666666666666d0 -6.0d0 0.0d0)
      ( 0.5d0 -6.0d0 0.0d0))
     (( 0.5d0 -6.0d0 0.0d0)
      (1.05d0 -5.3999996d0 0.0d0)
      ( 5.5d0 -0.5d0 0.0d0)
      ( 6.0d0  0.0d0 0.0d0))
     (( 6.0d0  0.0d0 0.0d0)
      ( 5.5d0  0.5d0 0.0d0)
      ( 0.5d0  5.5d0 0.0d0)
      ( 0.0d0  6.0d0 0.0d0))
     (( 0.0d0  6.0d0 0.0d0)
      (-0.5d0  5.5d0 0.0d0)
      (-5.5d0  0.5d0 0.0d0)
      (-6.0d0  0.0d0 0.0d0))
     ((-6.0d0  0.0d0 0.0d0)
      (-5.5d0 -0.5d0 0.0d0)
      (-1.05d0 -5.3999996d0 0.0d0)
      (-0.5d0 -6.0d0 0.0d0)))
    (((-0.5d0 -5.3d0 0.0d0)
      ( 0.5d0 -5.3d0 0.0d0))
     (( 0.5d0 -5.3d0 0.0d0)
      ( 5.0d0  0.0d0 0.0d0))
     (( 5.0d0  0.0d0 0.0d0)
      ( 0.0d0  5.0d0 0.0d0))
     (( 0.0d0  5.0d0 0.0d0)
      (-5.0d0  0.0d0 0.0d0))
     ((-5.0d0  0.0d0 0.0d0)
      (-0.5d0 -5.3d0 0.0d0)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '(( 0.5d0 -6.0d0 0.0d0)
      ( 6.0d0  0.0d0 0.0d0)
      ( 0.0d0  6.0d0 0.0d0)
      (-6.0d0  0.0d0 0.0d0)
      (-0.5d0 -6.0d0 0.0d0)))))


;;; Figure 7

;;; The three blend functions using several visualization methods (slicing/gradient/uv/spider)

(let ((*exponent* 2)
      (*resolution* 40)
      (points (points-from-angles '(50 60 70 60 70 50)))
      (distance-type 'line-sweep))
  (write-blends points '(t nil nil nil nil t) "n-sided-paper/07a-side-blend.vtk"
		:blend-function #'corner-blend :distance-type distance-type)
  (write-blends points '(t nil nil nil nil nil) "n-sided-paper/07b-corner-blend.vtk"
		:blend-function #'corner-blend :distance-type distance-type)
  (write-blends points '(t nil nil nil nil nil) "n-sided-paper/07c-special-side-blend.vtk"
		:blend-function #'ribbon-blend :distance-type distance-type)
  (write-blends-uv-polyline points '(t nil nil nil nil t)
			    "n-sided-paper/07a-side-blend-uv.vtk"
			    :blend-function #'corner-blend :distance-type distance-type)
  (write-blends-uv-polyline points '(t nil nil nil nil nil)
			    "n-sided-paper/07b-corner-blend-uv.vtk"
			    :blend-function #'corner-blend :distance-type distance-type)
  (write-blends-uv-polyline points '(t nil nil nil nil nil)
			    "n-sided-paper/07c-special-side-blend-uv.vtk"
			    :blend-function #'ribbon-blend :distance-type distance-type)
  (write-blends-spider-polyline points '(t nil nil nil nil t)
				"n-sided-paper/07a-side-blend-spider.vtk"
				:blend-function #'corner-blend :distance-type distance-type)
  (write-blends-spider-polyline points '(t nil nil nil nil nil)
				"n-sided-paper/07b-corner-blend-spider.vtk"
				:blend-function #'corner-blend :distance-type distance-type)
  (write-blends-spider-polyline points '(t nil nil nil nil nil)
				"n-sided-paper/07c-special-side-blend-spider.vtk"
				:blend-function #'ribbon-blend :distance-type distance-type))


;;; Figure 8

;;; Parameterizations

(defparameter *points* (points-from-angles '(40 20 60 100 80)))

(vectorized-distance-function-test
 *points* '(nil sd nil nil nil) "n-sided-paper/08a1-perpendicular-sd.ps"
 :resolution 0.001d0 :density 6 :distance-type 'perpendicular :color nil)
(vectorized-distance-function-test
 *points* '(nil sd nil nil nil) "n-sided-paper/08b1-barycentric-sd.ps"
 :resolution 0.001d0 :density 6 :distance-type 'barycentric :color nil)
(vectorized-distance-function-test
 *points* '(nil sd nil nil nil) "n-sided-paper/08c1-radial-sd.ps"
 :resolution 0.001d0 :density 6 :distance-type 'radial :color nil)
(vectorized-distance-function-test
 *points* '(nil sd nil nil nil) "n-sided-paper/08d1-chord-sd.ps"
 :resolution 0.001d0 :density 6 :distance-type 'chord-based :color nil)
(vectorized-distance-function-test
 *points* '(nil sd nil nil nil) "n-sided-paper/08e1-sweep-sd.ps"
 :resolution 0.001d0 :density 6 :distance-type 'line-sweep :color nil)

(vectorized-distance-function-test
 *points* '(s s s s s) "n-sided-paper/08a2-perpendicular-center.ps"
 :resolution 0.001d0 :density 2 :distance-type 'perpendicular)
(vectorized-distance-function-test
 *points* '(s s s s s) "n-sided-paper/08b2-barycentric-center.ps"
 :resolution 0.001d0 :density 2 :distance-type 'barycentric)
(vectorized-distance-function-test
 *points* '(s s s s s) "n-sided-paper/08c2-radial-center.ps"
 :resolution 0.001d0 :density 2 :distance-type 'radial)
(vectorized-distance-function-test
 *points* '(s s s s s) "n-sided-paper/08d2-chord-center.ps"
 :resolution 0.001d0 :density 2 :distance-type 'chord-based)
(vectorized-distance-function-test
 *points* '(s s s s s) "n-sided-paper/08e2-sweep-center.ps"
 :resolution 0.001d0 :density 2 :distance-type 'line-sweep)

(vectorized-distance-function-test
 *points* '(s s nil nil nil) "n-sided-paper/08a3-perpendicular-ss.ps"
 :resolution 0.001d0 :density 6 :distance-type 'perpendicular :color nil)
(vectorized-distance-function-test
 *points* '(s s nil nil nil) "n-sided-paper/08b3-barycentric-ss.ps"
 :resolution 0.001d0 :density 6 :distance-type 'barycentric :color nil)
(vectorized-distance-function-test
 *points* '(s s nil nil nil) "n-sided-paper/08c3-radial-ss.ps"
 :resolution 0.001d0 :density 6 :distance-type 'radial :color nil)
(vectorized-distance-function-test
 *points* '(s s nil nil nil) "n-sided-paper/08d3-chord-ss.ps"
 :resolution 0.001d0 :density 6 :distance-type 'chord-based :color nil)
(vectorized-distance-function-test
 *points* '(s s nil nil nil) "n-sided-paper/08e3-sweep-ss.ps"
 :resolution 0.001d0 :density 6 :distance-type 'line-sweep :color nil)


;;; Figure 9

;;; Blend territories

(defparameter *points* (points-from-angles '(40 20 60 100 80)))

(write-color-blend-test *points* "n-sided-paper/09a1-corner-perpendicular.ppm" 400
			:blend-function #'corner-blend
			:distance-type 'perpendicular
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09a2-corner-radial.ppm" 400
			:blend-function #'corner-blend
			:distance-type 'radial
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09a3-corner-barycentric.ppm" 400
			:blend-function #'corner-blend
			:distance-type 'barycentric
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09a4-corner-chord.ppm" 400
			:blend-function #'corner-blend
			:distance-type 'chord-based
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09a5-corner-sweep.ppm" 400
			:blend-function #'corner-blend
			:distance-type 'line-sweep
			:trim '(0.89d0 0.91d0))

(write-color-blend-test *points* "n-sided-paper/09b1-side-perpendicular.ppm" 400
			:blend-function #'ribbon-blend
			:distance-type 'perpendicular
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09b2-side-radial.ppm" 400
			:blend-function #'ribbon-blend
			:distance-type 'radial
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09b3-side-barycentric.ppm" 400
			:blend-function #'ribbon-blend
			:distance-type 'barycentric
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09b4-side-chord.ppm" 400
			:blend-function #'ribbon-blend
			:distance-type 'chord-based
			:trim '(0.89d0 0.91d0))
(write-color-blend-test *points* "n-sided-paper/09b5-side-sweep.ppm" 400
			:blend-function #'ribbon-blend
			:distance-type 'line-sweep
			:trim '(0.89d0 0.91d0))


;;; Figure 11 - three approaches

;;; Otoldalu:
(defparameter *points* (points-from-angles '(60 20 110 60 110)))
(defparameter *coords*
  '((((0.0d0 0.0d0 0.0d0)
      (1.0d0 0.0d0 1.0d0)
      (2.0d0 0.0d0 1.0d0)
      (2.4d0 0.0d0 0.3d0))
     ((2.4d0 0.0d0 0.3d0)
      (2.6d0 0.2d0 0.4d0)
      (2.8d0 0.4d0 0.4d0)
      (3.0d0 0.6d0 0.3d0))
     ((3.0d0 0.6d0 0.3d0)
      (3.0d0 2.0d0 1.0d0)
      (3.0d0 4.0d0 1.0d0)
      (3.0d0 6.0d0 0.0d0))
     ((3.0d0 6.0d0 0.0d0)
      (2.0d0 6.0d0 1.0d0)
      (1.0d0 6.0d0 1.0d0)
      (0.0d0 6.0d0 0.0d0))
     ((0.0d0 6.0d0 0.0d0)
      (0.0d0 4.0d0 1.0d0)
      (0.0d0 2.0d0 1.0d0)
      (0.0d0 0.0d0 0.0d0)))
    (((1.0d0 2.0d0 1.2d0)
      (2.2d0 0.7d0 1.2d0))
     ((2.2d0 0.7d0 1.2d0)
      (2.8d0 1.8d0 1.1d0))
     ((2.8d0 1.8d0 1.1d0)
      (2.0d0 4.0d0 1.2d0))
     ((2.0d0 4.0d0 1.2d0)
      (1.0d0 4.0d0 1.2d0))
     ((1.0d0 4.0d0 1.2d0)
      (1.0d0 2.0d0 1.2d0)))))

#+nil
(let ((*resolution* 20)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* t))
  (write-patch *points* 'hybrid "n-sided-paper/11a-side.vtk" :coords *coords*
	       :distance-type 'line-sweep))

#+nil
(let ((*resolution* 80)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* t))
  (write-patch *points* 'corner "n-sided-paper/11b-corner.vtk" :coords *coords*
	       :distance-type 'line-sweep))

#+nil
(let ((*resolution* 80)
      (*ribbon-multiplier* 1.0d0)
      (*centralized-line-sweep* t))
  (write-patch *points* 'ribbon "n-sided-paper/11c-special-side.vtk" :coords *coords*
	       :distance-type 'line-sweep))


;;; Torzitas teszt
(defparameter *coords*
  '((((3 0 0.0d0) (2 1 0.0d0) (1 2 0.0d0) (0 3 0.0d0))
     ((0 3 0.0d0) (-5 2 0.0d0) (-10 1 0.0d0) (-15 0 0.0d0))
     ((-15 0 0.0d0) (-10 -1 0.0d0) (-5 -2 0.0d0) (0 -3 0.0d0))
     ((0 -3 0.0d0) (1 -2 0.0d0) (2 -1 0.0d0) (3 0 0.0d0)))
    (((1 0 0.0d0) (0 1 0.0d0)) ((0 1 0.0d0) (-10 0 0.0d0))
     ((-10 0 0.0d0) (0 -1 0.0d0)) ((0 -1 0.0d0) (1 0 0.0d0)))))
(defparameter *points*
  (points-from-angles
   (angles-from-points
    '((0 1/5) (-1 0) (0 -1/5) (1/5 0)))))
(defparameter *points* (domain-from-curves (first *coords*)))
#+nil
(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 1.0d0))
  (write-patch *points* 'ribbon "n-sided-paper/rombusz1.vtk" :coords *coords*
	       :distance-type 'perpendicular :spider t))

;;; Teglalap
#+nil
(domain-from-curves
 '(((0 0 0) (6 0 0) (12 0 0))
   ((12 0 0) (12 1 0) (12 2 0))
   ((12 2 0) (6 2 8) (0 2 0))
   ((0 2 0) (0 1 0) (0 0 0))))
