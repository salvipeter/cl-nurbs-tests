(in-package :cl-nurbs-test)

;;; Test implementation of
;;; J.J. Zheng, A.A. Ball: Control point surfaces over non-four-sided areas (1996)


;;; Utilities

(defmacro defmemo (name args &body body)
  (cl-utilities:with-unique-names (cache rest win val)
    `(let ((,cache (make-hash-table :test #'equal)))
       (defun ,name ,args
         (let ((,rest ,(cons 'list args)))
           (multiple-value-bind (,val ,win) (gethash ,rest ,cache)
             (if ,win
                 ,val
                 (setf (gethash ,rest ,cache)
                       (progn ,@body)))))))))

(defun binomial (n k)
  (if (= k 0)
      1
      (* (binomial (1- n) (1- k))
         (/ n k))))


;;; Constants

(defparameter +side-constants+
  '(nil nil nil nil nil 1)
  "A list of c0, c1, ...")


;;; Implementation

(defun boundary-index-p (l)
  "An index is a boundary index if one of its values is zero."
  (find 0 l))

(defun all-shifts (lst)
  "E.g. (A B C) => ((A B C) (B C A) (C A B))."
  (labels ((rec (a b)
             (if (null a)
                 '()
                 (cons (append a b)
                       (rec (rest a) (append b (list (first a))))))))
    (rec lst '())))

(defmemo all-indices (n m)
  "Parameters:
N: # of sides
M: degree"
  (let ((base (iter (for l2 from 0 to (floor m 2))
                    (appending
                     (iter (for l3 from 0 to (floor m 2))
                           (collect (append (list (- m l3) l2 l3 (- m l2))
                                            (make-list (- n 4)
                                                       :initial-element (- m (min l2 l3))))))))))
    (delete-duplicates (mapcan #'all-shifts base) :test #'equal)))

(defun zb-blend (m l u)
  "As in Eq. (4.5). Parameters:
M: degree
L: control point index (a list of n values)
U: parametric point (a list of n values)"
  (if (boundary-index-p l)
      (let* ((n (length l))
             (i (position 0 l))
             (i-1 (mod (1- i) n))
             (i+1 (mod (1+ i) n))
             (k (elt l i-1)))
        (* (binomial m k)
           (expt (elt u i-1) k)
           (expt (elt u i+1) (- m k))
           (iter (for j from 0 below n)
                 (unless (member j (list i-1 i i+1))
                   (multiply (* (expt (elt u j) m)
                                (- 1 (* m (elt +side-constants+ n)
                                        (reduce #'* u)))))))))
      (let* ((min1 (reduce #'min l))
             (min2 (reduce #'min (remove min1 l :count 1))))
        (* (binomial m min1)
           (binomial m min2)
           (reduce #'* (mapcar #'expt u l))))))

(defun zb-all-blends-without-deficiencies (m u)
  "As in Eq. (4.13), called for all indices, returning a list. Parameters:
M: degree
U: parametric point (a list of n values)"
  (let* ((n (length u))
         (ls (all-indices n m))
         (blends (mapcar (lambda (l) (zb-blend m l u)) ls))
         (excess (/ (1- (reduce #'+ blends))
                    (if (evenp m)
                        (1+ (/ (* n (- m 2) m) 4))
                        (/ (* n (1- m) (1- m)) 4)))))
    (iter (for l in ls)
          (for blend in blends)
          (collect (if (boundary-index-p l)
                       blend
                       (- blend excess))))))

(defun rotate-list (i lst)
  (append (nthcdr i lst) (subseq lst 0 i)))

(defun zb-vertices (n)
  "Only works for N=5.
The center point is (phi,phi,phi,phi,phi), where phi = (sqrt(5)-1)/2.
Points on the `axes' are of the form (...,u,u,...), and we linearly
interpolate between two such values, always respecting the closer variable."
  (assert (= n 5))
  (let* ((phi (/ (1- (sqrt 5.0d0)) 2))
         (center (make-list n :initial-element phi))
	 (result (list center)))
    (iter (for j from 1 to *resolution*)
	  (for u = (* (- 1 (/ j *resolution*)) phi))
          (for u-to = (/ (1+ (sqrt (- 1 u)))))
          (for vlst = (let ((lst (iter (with max = (floor j 2))
                                       (with len = (if (or (<= j 2) (oddp j))
                                                       (- u-to u)
                                                       (* (- u-to u) (/ (- j 1) (- j 2)))))
                                       (for i from 0 below max)
                                       (collect (+ u (* (/ i max) len))))))
                        (append
                         (mapcar (lambda (x) (cons t x)) lst)
                         (if (evenp j) '() (list (cons t u-to)))
                         (nreverse (mapcar (lambda (x) (cons nil x)) lst)))))
	  (iter (for k from 0 below n)
		(iter (for (type . v) in vlst)
                      (if type
                          (let* ((u5 (- 1 (* u v)))
                                 (u1 (/ (- 1 v) u5))
                                 (u4 (- 1 (* u u1))))
                            (push (rotate-list k (list u1 u v u4 u5)) result))
                          (let* ((u4 (- 1 (* v u)))
                                 (u5 (/ (- 1 u) u4))
                                 (u3 (- 1 (* v u5))))
                            (push (rotate-list k (list v u u3 u4 u5)) result))))))
    (nreverse result)))


;;; Testing layer

(defun read-gbp (filename)
  "Reads a Generalized Bezier Patch into a hashtable.
Special items: SIDES, DEGREE and CENTER.
All control points are under (SIDE COL ROW).
Note that for even-degree patches, the center point is listed
also under (SIDE LAYER LAYER) for all sides."
  (let ((obj (make-hash-table :test 'equal)))
    (with-open-file (s filename)
      (flet ((read-point ()
               (let* ((x (read s)) (y (read s)) (z (read s)))
                 (list x y z))))
        (let* ((n (read s))
               (d (read s))
               (l (ceiling d 2))
               (cp (1+ (* n (1+ (floor d 2)) l)))
               (center (read-point))
               (i 1)
               (side 0)
               (col 0)
               (row 0))
          (setf (gethash 'sides obj) n)
          (setf (gethash 'degree obj) d)
          (setf (gethash 'center obj) center)
          (when (evenp d)
            (iter (for s from 0 below n)
                  (setf (gethash (list s l l) obj) center)))
          (loop while (< i cp) do
               (when (>= col (- d row))
                 (incf side)
                 (when (>= side n)
                   (setf side 0)
                   (incf row))
                 (setf col row))
               (let ((p (read-point)))
                 (setf (gethash (list side col row) obj) p)
                 (if (< col l)
                     (setf (gethash (list (mod (1- side) n) (- d row) col) obj) p)
                     (when (< (- d col) l)
                       (setf (gethash (list (mod (1+ side) n) row (- d col)) obj) p))))
               (incf i)
               (incf col)))))
    obj))

(defun convert-index (l)
  "Converts a Zheng-Ball type index into a (SIDE COL ROW) representation."
  (let* ((row (reduce #'min l))
         (side (position row l))
         (col (elt l (mod (1- side) (length l)))))
    (list side col row)))

(defun zb-eval (obj u)
  (let* ((n (gethash 'sides obj))
         (m (gethash 'degree obj))
         (blends (zb-all-blends-without-deficiencies m u))
         (cps (mapcar (lambda (l) (gethash (convert-index l) obj))
                      (all-indices n m))))
    (reduce #'v+ (mapcar #'v* cps blends))))


;;; Tests

#+nil
(defparameter *obj* (read-gbp "/home/salvi/project/transfinite/models/cagd86.gbp"))
#+nil
(defparameter *obj* (read-gbp "/home/salvi/Dropbox/Shares/GrafGeo/EuroGraphics/5sided.gbp"))

#+nil
(let* ((*resolution* 30)
       (points (zb-vertices 5))
       (data (mapcar (lambda (u) (zb-eval *obj* u)) points)))
  (write-obj-indexed-mesh data (triangles 5) "/tmp/proba.obj"))
