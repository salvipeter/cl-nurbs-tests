(in-package :cl-nurbs-tests)

;;; Problem:
;;; Given a closed, clamped bspline curve in the parametric interval [0, 1],
;;; we want to "move" the ends to a different position (x) on the curve,
;;; ie. the result should be a closed, clamped bspline curve in [0, 1], where
;;; the 0 point corresponds to the x point on the original curve.

;;; The solution works by splitting the curve in two, and then stitching them
;;; together in the other way around.

;;; (defparameter *profile* (first (first (read-rdn "/cygdrive/d/tmp/profiles.rdn"))))

(defun concatenate-curves (c1 c2)
  "There will be a knot with degree multiplicity at the joining point."
  (make-bspline-curve (degree c1)
		      (concatenate 'vector
				   (subseq (knot-vector c1)
					   0 (1- (length (knot-vector c1))))
				   (subseq (knot-vector c2) (1+ (degree c2))))
		      (concatenate 'vector
				   (control-points c1)
				   (subseq (control-points c2) 1))))

(defun append-curves (c1 c2)
  "Assumes that the parameterization is consistent
and the curves are of the same degree."
  (let ((u (elt (knot-vector c2) 0)))
    ;; BSC-REMOVE-KNOT not implemented yet
    ;; (iter (with curve = (concatenate-curves c1 c2))
    ;; 	  (repeat (degree c1))
    ;; 	  (setf curve (bsc-remove-knot curve u))
    ;; 	  (finally (return curve)))
    (concatenate-curves c1 c2)))

(defun shift-start-parameter (curve u)
  (destructuring-bind (c1 c2) (bsc-split-curve curve u)
    (bsc-shift-parameterization
     (append-curves c2 (bsc-shift-parameterization c1 1.0d0))
     (- u))))

;;; old code follows
#|
(defun insert-and-wrap (lst x)
  (flet ((less (y) (<= x y)))
    (append (list x)
	    (remove-if-not #'less lst)
	    (mapcar (lambda (y) (+ 1.0d0 y))
		    (remove-if #'less lst))
	    (list (+ 1.0d0 x)))))

(defun greville (curve i)
  (let ((d (degree curve))
	(knots (knot-vector curve)))
    (/ (iter (for k from 0 below d)
	     (sum (elt knots (+ i k 1))))
       d)))

(defun next-point (curve x)
  (iter (for i from 0 below (length (control-points curve)))
	(finding i such-that (> (greville curve i) x))))

(defun shift-start-parameter (curve shift)
  "CURVE is a closed curve in the parameter range [0, 1]."
  (assert (<= 0 shift 1))
  (with-accessors ((degree degree) (knots knot-vector) (points control-points))
      curve
    (let ((new-knots (append (iter (repeat degree) (collect shift))
			     (insert-and-wrap (subseq knots (1+ degree)
						      (- (length knots) degree 1))
					      shift)
			     (iter (repeat degree) (collect (+ shift 1)))))
	  (new-points (let ((end (bsc-evaluate curve shift))
			    (next (next-point curve shift)))
			(append (list end)
				(coerce (subseq points next (1- (length points))) 'list)
				(coerce (subseq points 1 next) 'list)
				(list end)))))
      (bsc-shift-parameterization
       (make-bspline-curve degree new-knots new-points)
       (- shift)))))
|#
