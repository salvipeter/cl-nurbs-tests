;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun unify-take-corner (patch u-end-p v-end-p)
  (let ((net (control-net patch)))
    (aref net
	  (if u-end-p (1- (array-dimension net 0)) 0)
	  (if v-end-p (1- (array-dimension net 1)) 0))))

(defun unify-corners (patch)
  (list (unify-take-corner patch nil nil)
	(unify-take-corner patch nil  t )
	(unify-take-corner patch  t   t )
	(unify-take-corner patch  t  nil)))

(defvar *unification-tolerance* 0)

(defun unify-match-corners (patch1 patch2)
  (let ((c1 (unify-corners patch1))
	(c2 (unify-corners patch2)))
    (iter (for i from 0 below 4)
	  (for c = (elt c1 i))
	  (for j = (position c c2 :test (lambda (c p)
					  (<= (point-distance c p)
					      *unification-tolerance*))))
	  (when j (collect (list i j))))))

(defun unify-neighborsp (patch1 patch2)
  (= (length (unify-match-corners patch1 patch2)) 2))

(defun unify-npatch (npatch)
  (iter (for s in npatch)
	(for neighbors = (iter (for s2 in npatch)
			       (when (unify-neighborsp s s2)
				 (collect s2))))
	(if (= (length neighbors) 4)
	    (collect s into inside)
	    (collect (list s neighbors) into outside))
	(finally
	 (return (list (iter (for s first (first inside) then next)
			     (for neighbors =
				  (iter (for s2 in inside)
					(when (unify-neighborsp s s2)
					  (collect s2))))
			     (for pprev first nil then prev)
			     (for prev first (first inside) then next)
			     (for next first (first neighbors)
				  then (if (eq (first neighbors) pprev)
					   (second neighbors)
					   (first neighbors)))
			     (collect next)
			     (while (not (eq next (first inside)))))
		       (iter (for (s lst) in outside)
			     (collect (list s (first (remove-if-not (lambda (x)
								      (member x inside))
								    lst))))))))))

(defun unified-indices (npatch)
  "Returns a list of two lists.
The first contains the inner surface indices in a circular walking order,
the second contains all pairs of adjacent outer and inner surface indices."
  (flet ((to-index (p) (position p npatch)))
    (let ((unified (unify-npatch npatch)))
      (list (mapcar #'to-index (first unified))
	    (mapcar (lambda (x) (mapcar #'to-index x)) (second unified))))))

;; (defparameter *test* (unified-indices (read-rbn "models/5sided.rbn")))

(defun unify-pair (s1 s2)
  "Returns two surfaces that are the reparameterized versions of S1 and S2,
such that S1 is on the domain [0..1, 0..1] and S2 is on [1..2, 0..1] and
they meet at U=1."
  (labels ((u-p (corners)
	     (let ((sorted (sort corners #'<)))
	       (or (equal sorted '(1 2)) (equal sorted '(0 3)))))
	   (flip (p1 p2)
	     (let ((match (unify-match-corners p1 p2)))
	       (list (if (u-p (mapcar #'first match))
			 (bss-flip-uv p1)
			 p1)
		     (if (u-p (mapcar #'second match))
			 (bss-flip-uv p2)
			 p2))))
	   (reverse-u (p1 p2)
	     (let ((match (unify-match-corners p1 p2)))
	       (list (if (equal (sort (mapcar #'first match) #'<) '(0 1))
			 (bss-reverse-parameterization p1 :u t)
			 p1)
		     (if (equal (sort (mapcar #'second match) #'<) '(2 3))
			 (bss-reverse-parameterization p2 :u t)
			 p2))))
	   (reverse-v (p1 p2)
	     (let ((match (unify-match-corners p1 p2)))
	       (list p1
		     (if (or (equal (first match) '(3 0))
			     (equal (first match) '(2 1)))
			 p2
			 (bss-reverse-parameterization p2 :v t))))))
    (destructuring-bind (a b)
	(apply #'reverse-v (apply #'reverse-u (flip s1 s2)))
      (list (bss-reparametrize a 0.0d0 1.0d0 0.0d0 1.0d0)
	    (bss-reparametrize b 1.0d0 2.0d0 0.0d0 1.0d0)))))

(defun g2-and-fair-npatch (npatch &optional (resolution 100) (iteration 100))
  (flet ((g2-and-reverse-if-needed (p1 p2)
	   (destructuring-bind (master slave)
	       (unify-pair (elt npatch p1) (elt npatch p2))
	     (let* ((g2 (ensure-g2-one-side slave master resolution :u-dir t))
		    (low (bss-lower-parameter (elt npatch p2)))
		    (high (bss-upper-parameter (elt npatch p2)))
		    (middle (affine-combine low 1/2 high)))
	       (setf (elt npatch p2)
		     (if (< (scalar-product
			     (bss-surface-normal slave '(3/2 1/2))
			     (bss-surface-normal (elt npatch p2) middle))
			    0)
			 (bss-reverse-parameterization g2 :u t)
			 g2))))))
    (let ((unified (unified-indices npatch)))
      (iter (for (p1 p2) in (second unified))
	    (g2-and-reverse-if-needed p1 p2))
      (iter (for p1 first (car (last (first unified))) then p2)
	    (for p2 in (first unified))
	    (g2-and-reverse-if-needed p1 p2))
      (iter (for index in (first unified))
	    (setf (elt npatch index)
		  (g-fair-krr-additive (elt npatch index) iteration 2)))
      npatch)))

(defun bss-approximate-size (surface)
  (let* ((low (bss-lower-parameter surface))
	 (high (bss-upper-parameter surface))
	 (mid (affine-combine low 1/2 high)))
    (flet ((get-curve (uv u-curve-p)
	     (bss-get-surface-curve surface uv :u-curve u-curve-p)))
      (list (/ (reduce #'+
		       (mapcar (lambda (x)
				 (bsc-estimate-arc-length (get-curve x t)))
			       (mapcar #'first (list low mid high))))
	       3)
	    (/ (reduce #'+
		       (mapcar (lambda (x)
				 (bsc-estimate-arc-length (get-curve x nil)))
			       (mapcar #'second (list low mid high))))
	       3)))))

(defun proportional-resolution (lst size)
  "LST contains surfaces, SIZE is the approximate total number of points.
The result is a list of pairs, corresponding to the resolution of each patch."
  (labels ((orientation (patch1 patch2)
	     (let* ((match (unify-match-corners patch1 patch2))
		    (c1 (sort (mapcar #'first match) #'<))
		    (c2 (sort (mapcar #'second match) #'<)))
	       (list (if (or (equal c1 '(0 3)) (equal c1 '(1 2))) 'u 'v)
		     (if (or (equal c2 '(0 3)) (equal c2 '(1 2))) 'u 'v)))))
    (let* ((n (length lst))
	   (lengths (mapcar #'bss-approximate-size lst))
	   (areas (mapcar (lambda (x) (reduce #'* x)) lengths))
	   (first-one (let ((area (* (first areas) (/ size (reduce #'+ areas))))
			    (ratio (apply #'/ (first lengths))))
			(list (round (sqrt (* area ratio)))
			      (round (sqrt (/ area ratio)))))))
      (cons first-one
	    (iter (for i upfrom 1)
		  (for last-patch first (first lst) then next-patch)
		  (for next-patch in (rest lst))
		  (for last first first-one then next)
		  (for (width height) in (rest lengths))
		  (for ratio = (/ width height))
		  (for orientation = (orientation last-patch next-patch))
		  (for inherited = (if (eq (first orientation) 'u)
				       (first last)
				       (second last)))
		  (for next = (if (= i (1- n))
				  (let* ((o2 (orientation
					      (first lst) (second lst)))
					 (tmp (if (eq (first o2) 'v)
						  (first first-one)
						  (second first-one))))
				    (if (eq (second orientation) 'u)
					(list inherited tmp)
					(list tmp inherited)))
				  (list (if (eq (second orientation) 'u)
					    inherited
					    (round (* inherited ratio)))
					(if (eq (second orientation) 'v)
					    inherited
					    (round (/ inherited ratio))))))
		  (collect next))))))

(let ((order '((0 1 3) (1 2) (3 2))))
  (defun joining-side (patch1 patch2)
    "Returns (UP ENDP REVERSEP), where
UP is T if the two patches meet along the U side of PATCH1,
ENDP is T if the two patches meet at the end of the other domain of PATCH1,
REVERSEP is T if the two patches have different directions at the joining side."
    (flet ((larger (x y) (not (member y (cdr (assoc x order))))))
      (let* ((match (unify-match-corners patch1 patch2))
	     (c1 (sort (mapcar #'first match) #'<)))
	(list (or (equal c1 '(0 3)) (equal c1 '(1 2)))
	      (or (equal c1 '(1 2)) (equal c1 '(2 3)))
	      (eq (larger (first (first match)) (first (second match)))
		  (larger (second (second match)) (second (first match)))))))))

(defun npatch-center (inner-npatch)
  (let ((corners (mapcan #'unify-corners inner-npatch)))
    (iter (for corner in corners)
	  (finding corner such-that
		   (= (funcall #'count corner corners
			       :test (lambda (p q)
				       (<= (point-distance p q)
					   *unification-tolerance*)))
		      (length inner-npatch))))))

(defun kobbelt-npatch (npatch size iteration &key vtk-file ply-file
		       (parameterization :projection))
  "SIZE is the approximate number of data points, which is divided between
the surfaces proportionally to their respective (approximate) area."
  (let* ((inner (first (unify-npatch npatch)))
	 (resolutions (proportional-resolution inner size))
	 (sizes (mapcar (lambda (x) (reduce #'* x)) resolutions))
	 (obj (kobbelt:initialize (reduce #'+ sizes)))
	 (center (npatch-center inner))
	 (center-index (kobbelt:insert-point obj center)))
    (iter (with n = (length inner))
	  (with first-side = nil)
	  (for i from 0 below n)
	  (format *error-output* "Inserting patch ~d...~%" (1+ i))
	  (for last first nil then patch)
	  (for patch in inner)
	  (for (up endp revp) first '(t t t) then (joining-side patch last))
	  (for (up2 endp2 revp2) =
	       (if (/= i (1- n)) '(t t t) (joining-side patch (first inner))))
	  (for (width height) in resolutions)
	  (for (u-low v-low) = (bss-lower-parameter patch))
	  (for (u-high v-high) = (bss-upper-parameter patch))
	  (for points = (make-array (list width height)))
	  (for last-side first nil then side)
	  (iter (for ui from 0 below width)
		(for u = (interpolate u-low (/ ui (1- width)) u-high))
		(iter (for vi from 0 below height)
		      (for v = (interpolate v-low (/ vi (1- height)) v-high))
		      (cond ((and (not (= i 0))
				  (or (and (= ui 0) (not up) (not endp))
				      (and (= vi 0) up (not endp))
				      (and (= ui (1- width)) (not up) endp)
				      (and (= vi (1- height)) up endp)))
			     ;; joint side of the current and previous patch
			     (let ((k (if up ui vi))
				   (len (length last-side)))
			       (setf (aref points ui vi)
				     (if revp
					 (elt last-side (- len 1 k))
					 (elt last-side k)))))
			    ((and (= i (1- n))
				  (or (and (= ui 0) (not up2) (not endp2))
				      (and (= vi 0) up2 (not endp2))
				      (and (= ui (1- width)) (not up2) endp2)
				      (and (= vi (1- height)) up2 endp2)))
			     ;; joint side of the first and last patch
			     (let ((k (if up2 ui vi))
				   (len (length first-side)))
			       (setf (aref points ui vi)
				     (if revp2
					 (elt first-side (- len 1 k))
					 (elt first-side k)))))
			    ((and (or (= ui vi 0)
				      (and (= ui 0) (= vi (1- height)))
				      (and (= ui (1- width)) (= vi 0))
				      (and (= ui (1- width))
					   (= vi (1- height))))
				  (let ((p (bss-evaluate patch (list u v))))
				    (<= (point-distance p center)
					*unification-tolerance*)))
			     ;; center
			     (setf (aref points ui vi) center-index))
			    ;; normal case
			    (t (let ((p (bss-evaluate patch (list u v))))
				 (setf (aref points ui vi)
				       (kobbelt:insert-point obj p)))))))
	  (unless (= i (1- n))
	    (flet ((accumulate-from-patch (patch-index)
		     (iter (with (up endp revp) =
				 (joining-side patch (elt inner patch-index)))
			   (declare (ignore revp))
			   (for ui from 0 below (if up width height))
			   (with vi = (cond ((not endp) 0)
					    (up (1- height))
					    (t (1- width))))
			   (collect (if up
					(aref points ui vi)
					(aref points vi ui))))))
	      ;; saving the side that meets with the next patch
	      (for side = (accumulate-from-patch (1+ i)))
	      ;; saving the side of the first patch that meets with the last
	      (when (= i 0) (setf first-side (accumulate-from-patch (1- n))))))
	  ;; set up the triangles
	  (flet ((p (i j) (aref points i j)))
	    (dotimes (i (1- width))
	      (dotimes (j (1- height))
		(kobbelt:set-triangle obj (p i j)
				      (p (1+ i) j) (p (1+ i) (1+ j)))
		(kobbelt:set-triangle obj (p i j)
				      (p (1+ i) (1+ j)) (p i (1+ j)))))))
    (format *error-output* "Fairing...~%")
    (kobbelt:fair obj iteration
		  :parameterization parameterization
		  :preserve-tangents t)
    (when vtk-file
      (format *error-output* "Writing VTK file...~%")
      (kobbelt:write-vtk-mesh obj vtk-file))
    (when ply-file
      (format *error-output* "Writing PLY file...~%")
      (kobbelt:write-ply-mesh obj ply-file))
    obj))

;;; Test
#+nil
(defparameter *five* (read-rbn "models/5sided.rbn"))
#+nil
(kobbelt-npatch *five* 500 0 :vtk-file "/tmp/npatch-original.vtk")
#+nil
(kobbelt-npatch *five* 500 100 :vtk-file "/tmp/npatch-faired.vtk")
