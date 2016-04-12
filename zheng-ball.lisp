(in-package :cl-nurbs-test)

;;; Test implementation of
;;; J.J. Zheng, A.A. Ball: Control point surfaces over non-four-sided areas (1996)

(defparameter +side-constants+
  '(nil nil nil nil nil 1)
  "A list of c0, c1, ...")

(defun binomial (n k)
  (if (= k 0)
      1
      (* (binomial (1- n) (1- k))
         (/ n k))))

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

(defun all-indices (n m)
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

(defun zb-blend-without-deficiency (m l u)
  "As in Eq. (4.13). Parameters:
M: degree
L: control point index (a list of n values)
U: parametric point (a list of n values)"
  (let ((n (length l)))
    (if (boundary-index-p l)
        (zb-blend m l u)
        (let ((deficiency (1- (reduce #'+
                                      (mapcar (lambda (l)
                                                (zb-blend m l u))
                                              (all-indices n m))))))
          (- (zb-blend m l u)
             (/ deficiency
                (if (evenp m)
                    (1+ (/ (* n (- m 2) m) 4))
                    (/ (* n (1- m) (1- m)) 4))))))))

(defun zb-parameters (n resolution)
  "Only works for N=5."
  (assert (= n 5))
  (iter (for u1 from 0 to 1 by (/ resolution))
        (appending
         (iter (for u3 from (- 1 u1) to 1 by (/ resolution))
               (cond ((zerop u1)
                      (collect (list 0 0 1 1 1))
                      (collect (list 0 1 1 1 0)))
                     ((zerop u3)
                      (collect (list 1 0 0 1 1))
                      (collect (list 1 1 0 0 1)))
                     (t (collect (list u1
                                       (/ (1- (+ u1 u3)) (* u1 u3))
                                       u3
                                       (/ (- 1 u1) u3)
                                       (/ (- 1 u3) u1)))))))))

(defun read-gbp (filename)
  "Reads a Generalized Bezier Patch into a hashtable.
Special items: SIDES, DEGREE and CENTER.
All control points are under (SIDE COL ROW).
Note that for even-degree patches, the center point is listed
under (SIDE LAYER LAYER) for all sides."
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
                     (setf (gethash (list side (- d row) col) obj) p)
                     (when (< (- d col) l)
                       (setf (gethash (list side row (- d col)) obj) p))))
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
  (let ((n (gethash 'sides obj))
        (m (gethash 'degree obj))
        (result '(0 0 0)))
    (iter (for l in (all-indices n m))
          (for blend = (zb-blend-without-deficiency m l u))
          (for cp = (gethash (convert-index l) obj))
          (setf result (v+ result (v* cp blend))))
    result))


;;; Tests

(defparameter *obj* (read-gbp "/home/salvi/project/transfinite/models/cagd86.gbp"))
(defparameter *points* (zb-parameters 5 30))
(defparameter *data* (mapcar (lambda (u) (zb-eval *obj* u)) *points*))
