;;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

;;; Symbols to export:
;;; - make-octree
;;; - octree-insert
;;; - octree-points
;;; - octree-position

(defclass octree ()
  ((bounds :initarg :bounds :reader octree-bounds)
   (regions :initarg :regions :accessor octree-regions)
   (depth :initarg :depth :reader octree-depth)
   (count :initform 0 :accessor octree-count)))

(defun make-octree (bbox depth)
  "BBOX is a list of two points that define the bounding box."
  (make-instance 'octree :bounds (list (apply #'mapcar #'min bbox)
				       (apply #'mapcar #'max bbox))
		 :depth depth :regions (list nil nil nil nil nil nil nil nil)))

(defun octree-index (octree point)
  (destructuring-bind (min max) (octree-bounds octree)
    (let ((med (affine-combine min 0.5d0 max)))
      (+ (if (< (first point) (first med)) 4 0)
	 (if (< (second point) (second med)) 2 0)
	 (if (< (third point) (third med)) 1 0)))))

(defun octree-region (octree index)
  (elt (octree-regions octree) index))

(defun (setf octree-region) (value octree index)
  (setf (elt (octree-regions octree) index)
	value))

(defun octree-subregion (octree index)
  (destructuring-bind (min max) (octree-bounds octree)
    (let ((med (affine-combine min 0.5d0 max)))
      (list (list (if (zerop (logand index 4)) (first med) (first min))
		  (if (zerop (logand index 2)) (second med) (second min))
		  (if (zerop (logand index 1)) (third med) (third min)))
	    (list (if (zerop (logand index 4)) (first max) (first med))
		  (if (zerop (logand index 2)) (second max) (second med))
		  (if (zerop (logand index 1)) (third max) (third med)))))))

(defun octree-insert (octree point)
  (let ((index (octree-index octree point)))
    (cond ((= (octree-depth octree) 1)
	   (unless (member point (octree-region octree index) :test #'equal)
	     (push point (octree-region octree index))
	     (incf (octree-count octree))))
	  ((null (octree-region octree index))
	   (let ((next (make-octree (octree-subregion octree index)
				    (1- (octree-depth octree)))))
	     (setf (octree-region octree index) next)
	     (when (octree-insert next point)
	       (incf (octree-count octree)))))
	  (t (when (octree-insert (octree-region octree index) point)
	       (incf (octree-count octree)))))))

(defun octree-points (octree)
  (if (= (octree-depth octree) 1)
      (reduce #'append (octree-regions octree))
      (mapcan #'octree-points (remove-if #'null (octree-regions octree)))))

(defun octree-position (octree point)
  "Returns the position of a point in the list returned by OCTREE-POINTS,
or NIL if the point was not in the octree.
WARNING: any point insertion may invalidate this index."
  (let* ((index (octree-index octree point))
	 (full (remove-if #'null (subseq (octree-regions octree) 0 index))))
    (if (= (octree-depth octree) 1)
	(let ((pos (position point (octree-region octree index) :test #'equal)))
	  (when pos
	    (+ (if (null full) 0 (reduce #'+ (mapcar #'length full)))
	       pos)))
	(let ((pos (octree-position (octree-region octree index) point)))
	  (when pos
	    (+ (if (null full) 0 (reduce #'+ (mapcar #'octree-count full)))
	       pos))))))
