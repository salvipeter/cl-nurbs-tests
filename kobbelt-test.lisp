;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :kobbelt)

(defun test (surface filename &key (resolution 100) (iteration 100) projectionp)
  (let ((obj (kobbelt:initialize (* resolution resolution)))
	(lower (bss-lower-parameter surface))
	(upper (bss-upper-parameter surface)))
    (iter (for ui from 0 below resolution)
	  (for u = (interpolate (first lower)
				(/ ui (1- resolution))
				(first upper)))
	  (iter (for vi from 0 below resolution)
		(for v = (interpolate (second lower)
				      (/ vi (1- resolution))
				      (second upper)))
		(for p = (bss-evaluate surface (list u v)))
		(insert-point obj p)))
    (flet ((tr (x y) (+ (* x resolution) y)))
      (dotimes (i (1- resolution))
	(dotimes (j (1- resolution))
	  (set-triangle obj (tr i j) (tr (1+ i) j) (tr (1+ i) (1+ j)))
	  (set-triangle obj (tr i j) (tr (1+ i) (1+ j)) (tr i (1+ j))))))
    (finalize obj :parameterization (if projectionp :projection :polar))
    (fair obj iteration)
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "~d ~d~%" resolution resolution)
      (iter (for i from 0 below (size obj))
	    (for p = (elt (points obj) i))
	    (for pos = (point-coordinates p))
	    (format s "~{~f~^ ~}~%" pos)))
    obj))

#+nil
(test (first (read-rbn "models/bottom.rbn")) "/tmp/bottom-faired.pts"
      :resolution 100 :iteration 100 :projectionp t)

(defun test-fn (fn filename &key (resolution 100) (iteration 100) projectionp)
  (let ((obj (kobbelt:initialize (* resolution resolution))))
    (iter (for u from 0 below resolution)
	  (iter (for v from 0 below resolution)
		(for p = (funcall fn u v))
		(insert-point obj p)))
    (flet ((tr (x y) (+ (* x resolution) y)))
      (dotimes (i (1- resolution))
	(dotimes (j (1- resolution))
	  (set-triangle obj (tr i j) (tr (1+ i) j) (tr (1+ i) (1+ j)))
	  (set-triangle obj (tr i j) (tr (1+ i) (1+ j)) (tr i (1+ j))))))
    (fair obj iteration :parameterization (if projectionp 'projection 'polar))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "~d ~d~%" resolution resolution)
      (iter (for i from 0 below (size obj))
	    (for p = (elt (points obj) i))
	    (for pos = (point-coordinates p))
	    (format s "~{~f~^ ~}~%" pos)))
    obj))

#+nil
(test-fn (lambda (u v) (list u v (+ (* u u) v))) "/tmp/fn-faired.pts"
      :resolution 10 :iteration 100 :projectionp t)
