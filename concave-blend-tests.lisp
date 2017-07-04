(in-package :cl-nurbs-tests)

(defun compute-parameter (type dir points p &optional no-tiny-p)
  "Modified version that allows NIL"
  (macrolet ((tiny-lambda ((args) &body body)
	       `(lambda (,args)
		  (if no-tiny-p
		      (let ((result (progn ,@body)))
			(if (and result (< result *tiny*)) 0.0d0 result))
		      (progn ,@body)))))
    (mapcar (tiny-lambda (lst) (compute-distance type points lst p dir))
	    (iter (for i from -2 below (- (length points) 2))
		  (collect (iter (for j from 0 below 4)
				 (collect (elt points (mod (+ i j) (length points))))))))))

(defun blend-optional (d i)
  "A variation of RIBBON-BLEND that allows infinite (NIL) distances."
  (let ((n (length d)))
    (/ (iter (for j from 0 below n)
	     (when (and (/= i j) (elt d j))
	       (multiply (expt (elt d j) *exponent*))))
       (iter (for k from 0 below n)
	     (sum (iter (for j from 0 below n)
			(when (and (/= j k) (elt d j))
			  (multiply (expt (elt d j) *exponent*)))))))))

(defmethod compute-distance ((type (eql 'concave)) points segments p dir)
  "Works only in the D direction. Returns NIL for some points."
  (let ((p0 (elt segments 1))
	(p1 (elt segments 2)))
    (assert (eq dir 'd) (dir) "This distance function is not defined in the `s' direction.")
    (let* ((d1 (point-line-distance p (list p0 p1)))
           (u (vnormalize (v- p1 p0)))
           (dl (scalar-product (v- p p0) u))
           (dr (scalar-product (v- p1 p) u)))
      (unless (or (< dl *epsilon*) (< dr *epsilon*))
        (/ d1 (min dl dr))))))

(defparameter *points* (points-from-angles '(20 100 80 120 40)))

#+nil
(write-color-blend-test
 *points*
 "/tmp/proba.ppm" 400
 :blend-function #'blend-optional
 :distance-type 'concave
 :trim '(0.89d0 0.91d0))

#+nil
(write-blends
 *points*
 '(nil t nil nil nil)
 "/tmp/proba.vtk"
 :blend-function #'blend-optional
 :distance-type 'concave)
