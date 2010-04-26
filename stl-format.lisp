;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defun read-float (stream)
  (let ((result 0))
    (dolist (x '(0 8 16 24))
      (setf (ldb (byte 8 x) result) (read-byte stream)))
    (ieee-floats:decode-float32 result)))

(defun write-float (stream f)
  (let ((32bit (ieee-floats:encode-float32 (coerce f 'single-float))))
    (dolist (c (mapcar (lambda (x) (ldb (byte 8 x) 32bit)) '(0 8 16 24)))
      (write-byte c stream))))

(defun read-unsigned-short (stream)
  (let ((result 0))
    (dolist (x '(0 8))
      (setf (ldb (byte 8 x) result) (read-byte stream)))
    result))

(defun write-unsigned-short (stream n)
  (dolist (c (mapcar (lambda (x) (ldb (byte 8 x) n)) '(0 8)))
    (write-byte c stream)))

(defun read-stl-ascii (filename &optional no-reverse)
  "If NO-REVERSE is NIL, the triangles are read in reverse order.
Returns a list."
  (with-open-file (s filename)
    (assert (eq (read s) 'solid) () "Invalid header.")
    (labels ((read3d () (list (read s) (read s) (read s)))
	     (expect (lst)
	       (let ((thing (read s)))
		 (if (or (and (consp lst) (member thing lst)) (eq thing lst))
		     thing
		     (error "Expected ~:[one of ~;~]~a, got ~a."
			    (atom lst) lst thing))))
	     (read-vertex () (expect 'vertex) (read3d)))
      (let ((name (read s))
	    (result ()))
	(iter (for code = (expect '(facet endsolid)))
	      (cond ((eq code 'facet)
		     (expect 'normal)
		     (let ((normal (read3d)))
		       (declare (ignore normal))
		       (expect 'outer) (expect 'loop)
		       (let* ((v1 (read-vertex))
			      (v2 (read-vertex))
			      (v3 (read-vertex)))
			 (push (list v1 v2 v3) result))
		       (expect 'endloop)
		       (expect 'endfacet)))
		    (t (expect name) (leave))))
	(values (if no-reverse (nreverse result) result)
		name)))))

(defun read-stl-binary (filename)
  "Returns a vector."
  (with-open-file (s filename :element-type '(unsigned-byte 8))
    (let* ((header (iter (repeat 80) (collect (read-byte s))))
	   (n (read-unsigned-integer s))
	   (result (make-array n)))
      (flet ((read3d ()
	       (let* ((x (read-float s))
		      (y (read-float s))
		      (z (read-float s)))
		 (list x y z))))
	(dotimes (i n)
	  (let* ((normal (read3d))
		 (v1 (read3d))
		 (v2 (read3d))
		 (v3 (read3d))
		 (attr (read-unsigned-short s)))
	    (declare (ignore normal attr))
	    (setf (elt result i) (list v1 v2 v3)))))
      (values result
	      (coerce (mapcar #'code-char header) 'string)))))

(defun read-stl (filename)
  (let ((code (with-open-file (s filename)
		(coerce (iter (repeat 5) (collect (read-char s))) 'string))))
    (if (string-equal code "solid")
	(read-stl-ascii filename)
	(read-stl-binary filename))))

(defun write-stl-ascii (triangles filename name)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "solid ~a~%" name)
    (iter (for i from 0 below (length triangles))
	  (for (v1 v2 v3) = (elt triangles i))
	  (for normal = (vnormalize (cross-product (v- v2 v1) (v- v3 v1))))
	  (format s "facet normal~{ ~f~}~%  outer loop~%    ~
                     vertex~{ ~f~}~%    vertex~{ ~f~}~%    vertex~{ ~f~}~%  ~
                     endloop~%endfacet~%"
		  normal v1 v2 v3))
    (format s "endsolid ~a~%" name)))

(defun write-stl-binary (triangles filename header)
  (let ((header (if (<= (length header) 80) header (subseq header 0 80)))
	(n (length triangles)))
    (with-open-file (s filename :direction :output :if-exists :supersede
		       :element-type '(unsigned-byte 8))
      (iter (with len = (length header))
	    (for i from 0 below 80)
	    (write-byte (if (< i len) (char-code (elt header i)) 32) s))
      (write-unsigned-integer s n)
      (iter (for i from 0 below n)
	    (for (v1 v2 v3) = (elt triangles i))
	    (for normal = (vnormalize (cross-product (v- v2 v1) (v- v3 v1))))
	    (dolist (f normal) (write-float s f))
	    (dolist (f v1) (write-float s f))
	    (dolist (f v2) (write-float s f))
	    (dolist (f v3) (write-float s f))
	    (write-unsigned-short s 0)))))

(defun write-stl (triangles filename &key ascii name header)
  (if ascii
      (write-stl-ascii triangles filename (or name "exported"))
      (write-stl-binary triangles filename (or header "Exported by CL-NURBS."))))