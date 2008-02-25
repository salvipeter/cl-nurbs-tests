;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun shift-knot-vectors (surface &key (u 0.0) (v 0.0))
  (let ((result (copy-bspline-surface surface)))
    (setf (knot-vectors result)
	  (mapcar (lambda (lst shift)
		    (mapcar (lambda (x) (+ x shift)) (coerce lst 'list)))
		  (knot-vectors result) (list u v)))
    result))

(defun flip-reverse-to-match (surface corners)
  (dolist (flipped (list surface (flip-uv surface)))
    (dolist (reversal '((nil nil) (t nil) (nil t) (t t)))
      (let* ((reversed (reverse-parametrization flipped
						:u (first reversal)
						:v (second reversal)))
	     (net (control-net reversed))
	     (bl (aref net 0 0))
	     (br (aref net (1- (array-dimension net 0)) 0))
	     (tl (aref net 0 (1- (array-dimension net 1))))
	     (tr (aref net
		       (1- (array-dimension net 0))
		       (1- (array-dimension net 1)))))
	(when (every (lambda (actual corner)
		       (or (not corner) (equal actual corner)))
		     (list bl br tl tr) corners)
	  (return-from flip-reverse-to-match reversed)))))
  nil)

(defun flip-reverse-shift (surface bl br tl tr u v uendp vendp)
  (let ((flipped (flip-reverse-to-match surface (list bl br tl tr))))
    (when flipped
      (let ((uknot (first (knot-vectors flipped)))
	    (vknot (second (knot-vectors flipped))))
	(shift-knot-vectors flipped
			    :u (- u (if uendp
					(elt uknot (1- (length uknot)))
					(elt uknot 0)))
			    :v (- v (if vendp
					(elt vknot (1- (length vknot)))
					(elt vknot 0))))))))

(defun write-xnode (vertex s1 s2 s3 s4 filename)
  (let* ((vertex-net (control-net vertex))
	 (bottom-left (aref vertex-net 0 0))
	 (bottom-right (aref vertex-net (1- (array-dimension vertex-net 0)) 0))
	 (top-left (aref vertex-net 0 (1- (array-dimension vertex-net 1))))
	 (top-right (aref vertex-net
			  (1- (array-dimension vertex-net 0))
			  (1- (array-dimension vertex-net 1))))
	 (uknot (first (knot-vectors vertex)))
	 (vknot (second (knot-vectors vertex)))
	 (uknot-min (elt uknot 0))
	 (vknot-min (elt vknot 0))
	 (uknot-max (elt uknot (1- (length uknot))))
	 (vknot-max (elt vknot (1- (length vknot)))))
    (flet ((find-suitable (bl br tl tr u v uendp vendp)
	     (iter (for sf in (list s1 s2 s3 s4))
		   (thereis (flip-reverse-shift sf bl br tl tr
						u v uendp vendp))
		   (finally (error "The control points do not match.")))))
      (write-rbn (list vertex
		       (find-suitable nil bottom-left nil top-left
				      uknot-min vknot-min t nil)
		       (find-suitable bottom-right nil top-right nil
				      uknot-max vknot-min nil nil)
		       (find-suitable nil nil bottom-left bottom-right
				      uknot-min vknot-min nil t)
		       (find-suitable top-left top-right nil nil
				      uknot-min vknot-max nil nil))
		 filename))))
