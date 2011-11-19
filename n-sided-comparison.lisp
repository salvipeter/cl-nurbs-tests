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

(defparameter *coords*
  '((((0.0d0 0.0d0 0.0d0)
      (0.6d0 0.0d0 1.0d0)
      (2.2d0 0.0d0 1.0d0)
      (2.4d0 0.0d0 0.3d0))
     ((2.4d0 0.0d0 0.3d0)
      (2.6d0 0.2d0 0.4d0)
      (2.8d0 0.4d0 0.4d0)
      (3.0d0 0.6d0 0.3d0))
     ((3.0d0 0.6d0 0.3d0)
      (3.0d0 1.0d0 0.6d0)
      (3.0d0 5.0d0 1.0d0)
      (3.0d0 6.0d0 0.0d0))
     ((3.0d0 6.0d0 0.0d0)
      (2.4d0 6.0d0 1.0d0)
      (0.6d0 6.0d0 1.0d0)
      (0.0d0 6.0d0 0.0d0))
     ((0.0d0 6.0d0 0.0d0)
      (0.0d0 5.0d0 1.0d0)
      (0.0d0 1.0d0 1.0d0)
      (0.0d0 0.0d0 0.0d0)))
    (((0.6d0 0.8d0 1.2d0)
      (2.4d0 0.3d0 1.2d0))
     ((2.4d0 0.3d0 1.2d0)
      (2.7d0 0.7d0 1.0d0))
     ((2.7d0 0.7d0 1.0d0)
      (2.5d0 5.0d0 1.2d0))
     ((2.5d0 5.0d0 1.2d0)
      (0.6d0 5.0d0 1.2d0))
     ((0.6d0 5.0d0 1.2d0)
      (0.6d0 0.8d0 1.2d0)))))

(write-constraint-ribbons *points* "/tmp/comp/ribbons.vtk"
			  :coords *coords* :resolution 20)
(write-constraint-coons-ribbons *points* "/tmp/comp/coons-ribbons"
				:coords *coords* :resolution 20)
(let ((*resolution* 40)
      (*ribbon-multiplier* 1.0d0))
  (write-patch *points* 'hybrid "/tmp/comp/bq-patch.vtk" :coords *coords*
	       :distance-type 'biquadratic)
  (write-patch *points* 'hybrid "/tmp/comp/bq-spider.vtk" :coords *coords*
	       :distance-type 'biquadratic :spider t)
  (write-patch *points* 'hybrid "/tmp/comp/il-bl-patch.vtk" :coords *coords*
	       :distance-type 'bilinear-mod)
  (write-patch *points* 'hybrid "/tmp/comp/il-bl-spider.vtk" :coords *coords*
	       :distance-type 'bilinear-mod :spider t)
  (write-patch *points* 'hybrid "/tmp/comp/il-cs-patch.vtk" :coords *coords*
	       :distance-type 'line-sweep-mod)
  (write-patch *points* 'hybrid "/tmp/comp/il-cs-spider.vtk" :coords *coords*
	       :distance-type 'line-sweep-mod :spider t)
  (let ((*proportion* 1/2))
    (write-patch *points* 'hybrid "/tmp/comp/tp-patch.vtk" :coords *coords*
		 :distance-type 'two-parabola)
    (write-patch *points* 'hybrid "/tmp/comp/tp-spider.vtk" :coords *coords*
		 :distance-type 'two-parabola :spider t))
  (let ((*proportion* 2/3))
    (write-patch *points* 'hybrid "/tmp/comp/tp-f-patch.vtk" :coords *coords*
		 :distance-type 'two-parabola)
    (write-patch *points* 'hybrid "/tmp/comp/tp-f-spider.vtk" :coords *coords*
		 :distance-type 'two-parabola :spider t)))

(let ((*resolution* 40)
      (*ribbon-multiplier* 1.0d0)
      (*proportion* 1/2))
  (write-patch *points* 'hybrid-coons "/tmp/comp/coons-tp-patch.vtk" :coords *coords*
		 :distance-type 'two-parabola)
  (write-patch *points* 'hybrid-coons "/tmp/comp/coons-tp-spider.vtk" :coords *coords*
	       :distance-type 'two-parabola :spider t))

(let ((points (points-from-angles '(40 20 60 100 80))))
  (labels ((filename (distance n)
	     (format nil "/tmp/comp/~(~a~)-~d.ps" distance n))
	   (write-test (type distance n lst)
	     (vectorized-distance-function-test
	      points lst (filename type n)
	      :resolution 0.001d0 :density 6
	      :distance-type distance :color nil)))
    (iter (for type in '(biquadratic biquadratic-corner bilinear-mod line-sweep-mod
			 two-parabola two-parabola-full))
	  (for *proportion* = (if (eq type 'two-parabola-full) 3/4 1/2))
	  (for distance = (if (eq type 'two-parabola-full) 'two-parabola type))
	  (write-test type distance 1 '(sd nil nil nil nil))
	  (write-test type distance 2 '(nil sd nil nil nil))
	  (write-test type distance 3 '(nil nil sd nil nil)))))
