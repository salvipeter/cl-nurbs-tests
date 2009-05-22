;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun project-to-plane (point plane)
  "PLANE is represented by a list \(origin normal-vector)."
  (let ((len (scalar-product (v- (first plane) point) (second plane))))
    (v+ point (v* (second plane) len))))

(defun bsc-ideal-control-point (curve i iteration &optional (alpha 1))
  (flet ((knot (k) (elt (knot-vector curve) (+ k (degree curve))))
	 (cpoint (k) (elt (control-points curve) k)))
    (let* ((u (knot i))
	   (d (vnormalize (bsc-evaluate curve u :derivative 1)))
	   (n (bsc-out-direction curve u))
	   (p2 (cpoint (- i 2))) (p1 (cpoint (- i 1)))
	   (r1 (cpoint (+ i 1))) (r2 (cpoint (+ i 2)))
	   (plane (list (affine-combine p1 1/2 r1) (cross-product d n)))
	   (p2v (project-to-plane p2 plane)) (p1v (project-to-plane p1 plane))
	   (qv (project-to-plane (cpoint i) plane))
	   (r2v (project-to-plane r2 plane)) (r1v (project-to-plane r1 plane)))
      (labels ((from-plane (v)
		 (destructuring-bind (x y) v
		   (v+ qv (v* n x) (v* d y))))
	       (angle (v1 v2)
		 (let ((n1 (vnormalize v1))
		       (n2 (vnormalize v2)))
		   (* (acos (- (scalar-product n1 n2)))
		      (let ((x (cross-product n1 n2)))
			(if (> (scalar-product x (second plane)) 0) 1 -1)))))
	       (error-fn (v)
		 (let* ((q (from-plane v))
			(v-2 (v- p1v p2v)) (v-1 (v- q p1v))
			(v+1 (v- r1v q)) (v+2 (v- r2v r1v))
			(len-1 (- (vlength v-1) (vlength v-2)))
			(len   (- (vlength v+1) (vlength v-1)))
			(len+1 (- (vlength v+2) (vlength v+1)))
			(ang-1 (angle v-2 v-1))
			(ang   (angle v-1 v+1))
			(ang+1 (angle v+1 v+2)))
		   (+ (expt (- (interpolate ang-1 1/2 ang+1) ang) 2)
		      (* alpha
			 (+ (expt (- len len-1) 2)
			    (expt (- len+1 len) 2)))))))
	(from-plane (downhill-simplex:minimize #'error-fn '(0 0) iteration))))))

(defun limit-change (old new limit)
  "Moves OLD towards NEW by length LIMIT."
  (v+ old (v* (vnormalize (v- new old)) limit)))

;;; TODO:
;;; - egysegesiteni kell a hibafuggveny tagjait
;;; - kene vmi "ertelmes" tag is (gorbulet?)
;;; - mi alapjan valasztja ki a kovetkezo fairelendo kontrollpontot?
;;; - mennyit mozgasson egyszerre? (user parameter?)

(defun bsc-2d-to-3d (curve)
  (let* ((new-curve (copy-bspline-curve curve))
	 (cps (control-points new-curve)))
    (dotimes (i (length cps))
      (setf (elt cps i) (list (first (elt cps i)) (second (elt cps i)) 0)))
    new-curve))

(defun bsc-3d-to-2d (curve)
  (let* ((new-curve (copy-bspline-curve curve))
	 (cps (control-points new-curve)))
    (dotimes (i (length cps))
      (setf (elt cps i) (subseq (elt cps i) 0 2)))
    new-curve))

#+nil
(defun test (lst)
  (let ((curve (copy-bspline-curve *curve*)))
    (dolist (i lst)
      (let ((ideal (bsc-ideal-control-point (bsc-2d-to-3d curve) i 100)))
	(setf (elt (control-points curve) i) (subseq ideal 0 2))))
    (write-ps curve "/tmp/proba.ps" 100)))

#+nil
(defun test2 (n)
  (let ((lst (iter (with len = (length (control-points *curve*)))
		   (repeat n)
		   (collect (+ 2 (random (- len 4)))))))
    (test lst)))

#+nil
(defun bss-ideal-control-point (surface i j &key u-dir)
  )
