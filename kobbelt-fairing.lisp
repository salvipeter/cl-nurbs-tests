(in-package :cl-nurbs-tests)

(defun bss-kobbelt-faired-points (surface resolution iteration &key
				  (parameterization :projection)
				  (preserve-tangents t))
  (destructuring-bind (r1 r2) resolution
    (let ((obj (kobbelt:initialize (* r1 r2)))
	  (points (sample-surface surface resolution)))
      (dotimes (i r1)
	(dotimes (j r2)
	  (kobbelt:insert-point obj (aref points i j))))
      (flet ((p (x y) (+ (* x r2) y)))
	(dotimes (i (1- r1))
	  (dotimes (j (1- r2))
	    (kobbelt:set-triangle obj (p i j) (p (1+ i) j) (p (1+ i) (1+ j)))
	    (kobbelt:set-triangle obj (p i j) (p (1+ i) (1+ j)) (p i (1+ j))))))
      (kobbelt:fair obj iteration
		    :parameterization parameterization
		    :preserve-tangents preserve-tangents)
      (iter (with index = 0)
	    (for i from 0 below r1)
	    (iter (for j from 0 below r2)
		  (setf (aref points i j)
			(kobbelt:get-point obj index))
		  (incf index)))
      points)))
