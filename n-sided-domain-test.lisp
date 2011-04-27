(in-package :cl-nurbs-tests)

;;; standard tesztpatch
(defparameter *coords*
  '((((0.0d0 0.0d0 0.0d0)
      (0.8d0 0.0d0 0.0d0)
      (1.6d0 0.0d0 0.0d0)
      (2.4d0 0.0d0 0.0d0))
     ((2.4d0 0.0d0 0.0d0)
      (2.6d0 0.2d0 0.0d0)
      (2.8d0 0.4d0 0.0d0)
      (3.0d0 0.6d0 0.0d0))
     ((3.0d0 0.6d0 0.0d0)
      (3.0d0 2.4d0 0.0d0)
      (3.0d0 4.2d0 0.0d0)
      (3.0d0 6.0d0 0.0d0))
     ((3.0d0 6.0d0 0.0d0)
      (2.0d0 6.0d0 0.0d0)
      (1.0d0 6.0d0 0.0d0)
      (0.0d0 6.0d0 0.0d0))
     ((0.0d0 6.0d0 0.0d0)
      (0.0d0 4.0d0 0.0d0)
      (0.0d0 2.0d0 0.0d0)
      (0.0d0 0.0d0 0.0d0)))
    (((0.8471177944862156d0 1.4135338345864663d0 0.0d0)
      (1.9707762557077626d0 0.6356164383561644d0 0.0d0))
     ((1.9707762557077626d0 0.6356164383561644d0 0.0d0)
       (2.6098360655737705d0 1.7311475409836063d0 0.0d0))
     ((2.6098360655737705d0 1.7311475409836063d0 0.0d0)
      (2.2641509433962264d0 4.150943396226415d0 0.0d0))
     ((2.2641509433962264d0 4.150943396226415d0 0.0d0)
      (0.9354120267260579d0 4.062360801781738d0 0.0d0))
     ((0.9354120267260579d0 4.062360801781738d0 0.0d0)
      (0.8471177944862156d0 1.4135338345864663d0 0.0d0)))))

;;; uj torzitas teszt
(defparameter *coords*
  (generate-planar-patch-by-points
   '((572 783) (734 629) (626 5) (87 495) (84 593) (153 668))))

(write-constraint-ribbons nil "/tmp/ribbon.vtk" :coords *coords* :resolution 20)
(let ((*resolution* 30)
      (*centralized-line-sweep* 1.0d0)
      (*ribbon-multiplier* 1.0d0)
      (patch-type 'ribbon)
      (distance-type 'line-sweep)
      (domain-type 'angular)) ; type: regular/circular/circular-mod/angular
  (write-patch (domain-from-curves (first *coords*) domain-type) patch-type "/tmp/proba.vtk"
	       :coords *coords* :distance-type distance-type :spider t))

(let ((*centralized-line-sweep* 1.0d0)
      (distance-type 'line-sweep)
      (domain-type 'circular-mod))
  (vectorized-distance-function-test
   (domain-from-curves (first *coords*) domain-type)
   '(s s s s s s) "/tmp/proba.ps"
   :resolution 0.001d0 :density 4 :distance-type distance-type :color t))

;;; Utils

(defun line-line-intersection (l1 l2)
  (let ((l1 (mapcar (lambda (p) (subseq p 0 2)) l1))
	(l2 (mapcar (lambda (p) (subseq p 0 2)) l2)))
    (destructuring-bind (((x1 y1) (x2 y2)) ((x3 y3) (x4 y4))) (list l1 l2)
      (list (/ (- (* (- (* x1 y2) (* y1 x2)) (- x3 x4))
		  (* (- (* x3 y4) (* y3 x4)) (- x1 x2)))
	       (- (* (- x1 x2) (- y3 y4)) (* (- y1 y2) (- x3 x4))))
	    (/ (- (* (- (* x1 y2) (* y1 x2)) (- y3 y4))
		  (* (- (* x3 y4) (* y3 x4)) (- y1 y2)))
	       (- (* (- x1 x2) (- y3 y4)) (* (- y1 y2) (- x3 x4))))))))

(defun patch2d->3d (patch)
  (list (mapcar (lambda (curve) (mapcar (lambda (p) (append p (list 0))) curve)) (first patch))
	(mapcar (lambda (curve) (mapcar (lambda (p) (append p (list 0))) curve)) (second patch))))

(defun generate-planar-patch-by-points (points)
  (let ((curves (iter (for x in points)
		      (for y in (append (rest points) (list (first points))))
		      (collect (list x (affine-combine x 0.25 y) (affine-combine x 0.75 y) y)))))
    (let ((lines (iter (for c1 in (append (last curves) (butlast curves)))
		       (for c2 in (append (rest curves) (list (first curves))))
		       (collect (list (third c1) (second c2))))))
      (patch2d->3d
       (list curves
	     (iter (for l1 in (append (last lines) (butlast lines)))
		   (for l2 in lines)
		   (for l3 in (append (rest lines) (list (first lines))))
		   (collect (list (line-line-intersection l1 l2)
				  (line-line-intersection l2 l3)))))))))

