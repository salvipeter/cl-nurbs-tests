;;; N-patch evaluation (for "5sided.rbn", but should work on other models)
(in-package :cl-nurbs-tests)

(defparameter *model* (read-rbn "models/5sided.rbn"))
(defparameter *model* (append (read-rbn "models/5sided-fair-g2-inner.rbn")
			      (read-rbn "models/5sided-fair-sides.rbn")))

(defparameter *inner* (first (unified-indices *model*)))
(defparameter *inner-pairs*
  (iter (with n = (length *inner*))
	(for i from 0 below n)
	(collect (list (elt *inner* i) (elt *inner* (mod (1+ i) n))))))

(defparameter *outer*
  (iter (for i from 0 below (length *model*))
	(unless (member i *inner*)
	  (collect i))))
(defparameter *outer-pairs*
  (iter (for i1s on *outer*)
	(for i1 = (first i1s))
	(appending
	 (iter (for i2 in (rest i1s))
	       (when (unify-neighborsp (elt *model* i1) (elt *model* i2))
		 (collect (list i1 i2)))))))

(defparameter *mixed-pairs*
  (second (unified-indices *model*)))

(defun safe-acos (x)
  (cond ((> x 1) 0)
	((< x -1) pi)
	(t (acos x))))

(defun g1-deviation (p1 p2 resolution)
  (destructuring-bind (s1 s2) (unify-pair p1 p2)
    (iter (for i from 0 below resolution)
	  (for v = (/ i (1- resolution)))
	  (for uv = (list 1 v))
	  (for n1 = (bss-surface-normal s1 uv))
	  (for n2 = (bss-surface-normal s2 uv))
	  (collect (* (safe-acos (scalar-product n1 n2)) #.(/ 180.0d0 pi))))))

(defun g2-deviation (p1 p2 resolution)
  (destructuring-bind (s1 s2) (unify-pair p1 p2)
    (iter (for i from 0 below resolution)
	  (for v = (/ i (1- resolution)))
	  (for uv = (list 1 v))
	  (for d1 = (bss-evaluate s1 uv :derivative '(1 0)))
	  (for d2 = (bss-evaluate s2 uv :derivative '(1 0)))
	  (for dir = (v* (v+ (vnormalize d1) (vnormalize d2)) 1/2))
	  (for n1 = (normal-curvature s1 uv dir))
	  (for n2 = (normal-curvature s2 uv dir))
	  (collect (abs (/ (- n1 n2) n2))))))

(iter (with resolution = 100)
      (for type in '(inner outer mixed))
      (for pairs in (list *inner-pairs* *outer-pairs* *mixed-pairs*))
      (for deviations = (mapcar (lambda (pair)
				  (/ (reduce #'+
					     (g1-deviation
					      (elt *model* (first pair))
					      (elt *model* (second pair))
					      resolution))
				     resolution))
				pairs))
      (format t "[~a] Mean tangent direction deviations in degrees:~%~
                 ~{~f ~}~%(min: ~f, max: ~f, mean: ~f)~%"
	      type deviations (reduce #'min deviations)
	      (reduce #'max deviations)
	      (/ (reduce #'+ deviations) (length deviations))))

(iter (with resolution = 100)
      (for type in '(inner outer mixed))
      (for pairs in (list *inner-pairs* *outer-pairs* *mixed-pairs*))
      (for deviations = (mapcar (lambda (pair)
				  (/ (reduce #'+
					     (g2-deviation
					      (elt *model* (first pair))
					      (elt *model* (second pair))
					      resolution))
				     resolution))
				pairs))
      (format t "[~a] Mean normal curvature deviation percentages:~%~
                 ~{~f ~}~%(min: ~f, max: ~f, mean: ~f)~%"
	      type deviations (reduce #'min deviations)
	      (reduce #'max deviations)
	      (/ (reduce #'+ deviations) (length deviations))))

(defun report (lst)
  (format t "Min: ~f~%Max: ~f~%Mean: ~f~%Median: ~f~%Std. Deviation: ~f~%"
	  (reduce #'min lst) (reduce #'max lst) (statistics:mean lst)
	  (statistics:median lst) (statistics:standard-deviation lst)))

(defun evaluate-g1-xnode (xnode &optional (resolution 100))
  (format t "[LEFT SURFACE]~%")
  (report (g1-deviation (first xnode) (second xnode) resolution))
  (format t "[RIGHT SURFACE]~%")
  (report (g1-deviation (first xnode) (third xnode) resolution))
  (format t "[TOP SURFACE]~%")
  (report (g1-deviation (first xnode) (fourth xnode) resolution))
  (format t "[BOTTOM SURFACE]~%")
  (report (g1-deviation (first xnode) (fifth xnode) resolution)))

(defun evaluate-g2-xnode (xnode &optional (resolution 100))
  (format t "[LEFT SURFACE]~%")
  (report (g2-deviation (first xnode) (second xnode) resolution))
  (format t "[RIGHT SURFACE]~%")
  (report (g2-deviation (first xnode) (third xnode) resolution))
  (format t "[TOP SURFACE]~%")
  (report (g2-deviation (first xnode) (fourth xnode) resolution))
  (format t "[BOTTOM SURFACE]~%")
  (report (g2-deviation (first xnode) (fifth xnode) resolution)))
