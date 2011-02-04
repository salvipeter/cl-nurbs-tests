(defun write-spider-in-domain (points filename)
  (write-vtk-polylines
   (iter (for line in (spider-lines points))
	 (collect (iter (for domain-point in line)
			(collect (cons 0 domain-point)))))
   filename))

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

;;; Domain kepek
(let ((*resolution* 30))
  (iter (for type in '(regular circular angular))
	(write-spider-in-domain (domain-from-curves (first *coords*) type)
				(format nil "n-sided-paper/rombusz/domain-~(~a~).vtk" type))))

;;; Distance kepek
(vectorized-distance-function-test (domain-from-curves (first *coords*) 'angular) '(sd sd sd sd)
				   "n-sided-paper/rombusz/angular-perpendicular-sd.ps"
				   :resolution 0.001d0 :density 6 :distance-type 'perpendicular)

;;; Patchek
(let ((*resolution* 30)
      (*centralized-line-sweep* nil)
      (*ribbon-multiplier* 1.0d0))
  (iter (for patch-type in '(ribbon))
	(iter (for domain-type in '(regular angular circular))
	      (iter (for distance-type in '(perpendicular radial line-sweep))
		    (write-patch (domain-from-curves (first *coords*) domain-type) patch-type
				 (format nil "n-sided-paper/rombusz/~(~a~)-~(~a~)-~(~a~).vtk"
					 patch-type domain-type distance-type)
				 :coords *coords* :distance-type distance-type :spider t)))))
