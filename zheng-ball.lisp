(in-package :cl-nurbs-tests)

;;; Test implementation of
;;; J.J. Zheng, A.A. Ball: Control point surfaces over non-four-sided areas (1996)

(defparameter +side-constants+
  '(nil nil nil nil nil 1 2)
  "A list of c0, c1, ...")

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

(defun binomial (n k)
  (if (= k 0)
      1
      (* (binomial (1- n) (1- k))
         (/ n k))))

(defun zb-blend (m l u)
  "As in Eq. (4.5); see also Ball-Zheng (2001) for this shorter form. Parameters:
M: degree
L: control point index (a list of n values)
U: parametric point (a list of n values)"
  (let* ((min1 (reduce #'min l))
         (min2 (reduce #'min (remove min1 l :count 1))))
    (* (binomial m min1)
       (binomial m min2)
       (reduce #'* (mapcar #'expt u l))
       (if (zerop min1)
           (- 1 (* m (elt +side-constants+ (length l))
                   (reduce #'* u)))
           1))))

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
  "E.g.: 2 (A B C D E) => (C D E A B)."
  (append (nthcdr i lst) (subseq lst 0 i)))

(defgeneric zb-vertices (n))

(defmethod zb-vertices ((n (eql 5)))
  "The center point is (phi,phi,phi,phi,phi), where phi = (sqrt(5)-1)/2.
Points on the `axes' are of the form (...,u,u,...), and we linearly
interpolate between two such values, always respecting the closer variable."
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

(defun quadratic-solver (a b c min max)
  (if (< (abs a) *epsilon*)
      (if (< (abs b) *epsilon*)
          (error "0-degree equation")
          (/ (- c) b))
      (let* ((d (sqrt (- (* b b) (* 4 a c))))
             (s1 (/ (+ (- b) d) (* 2 a)))
             (s2 (/ (- (- b) d) (* 2 a)))
             (s1e (max (- min s1) (- s1 max) 0))
             (s2e (max (- min s2) (- s2 max) 0)))
        (if (< s1e s2e)
            s1
            s2))))

(defun generate-parameter6 (i x y)
  (rotate-list (- 6 i)
   (append (list x y)
           (cond ((= x 0)
                  (list 1 1 1 (- 1 y)))
                 ((= y 0)
                  (list (- 1 x) 1 1 1))
                 (t (iter (repeat 4)
                          (for u1 first x then u2)
                          (for u2 first y then r)
                          (for a = (* (- 1 (* u1 u2)) (- 1 (* 2 u1 u2))))
                          (for b = (+ (* 2 u1) (* -3 u1 u1 u2) (* u1 u2 u2)))
                          (for c = (1- (* u1 u1)))
                          (for r = (quadratic-solver a b c 0 1))
                          (collect r)))))))

(defun parameter6-on-segment-half (u)
  "Finds the parameter X such that (X U X ...) is valid."
  (iter (for x from 0.5 to (/ (sqrt 2) 2) by 1.0d-5)
        (finding x minimizing (abs (- x (elt (generate-parameter6 6 x u) 2))))))

(defmethod zb-vertices ((n (eql 6)))
  "The center point is (eta,eta,eta,eta,eta), where eta = sqrt(2)/2.
Points on the `axes' are of the form (...,u,u,...), and we linearly
interpolate between two such values, always respecting the closer variable."
  (let* ((eta (/ (sqrt 2.0d0) 2))
         (center (make-list n :initial-element eta))
	 (result (list center)))
    (iter (for j from 1 to *resolution*)
	  (for u = (* (- 1 (/ j *resolution*)) eta))
          (for u-to = (parameter6-on-segment-half u))
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
                          (push (generate-parameter6 k v u) result)
                          (push (generate-parameter6 (1+ k) u v) result)))))
    (nreverse result)))


;;; Testing layer

(defun gbp-set-center (obj p)
  (let ((n (gethash 'sides obj))
        (d (gethash 'degree obj)))
    (setf (gethash 'center obj) p)
    (when (evenp d)
      (iter (with l = (gethash 'layers obj))
            (for s from 0 below n)
            (setf (gethash (list s l l) obj) p)))))

(defun gbp-set-cp (obj side col row p)
  (let ((n (gethash 'sides obj))
        (d (gethash 'degree obj))
        (l (gethash 'layers obj)))
    (setf (gethash (list side col row) obj) p)
    (if (< col l)
        (setf (gethash (list (mod (1- side) n) (- d row) col) obj) p)
        (when (< (- d col) l)
          (setf (gethash (list (mod (1+ side) n) row (- d col)) obj) p)))))

(defun read-gbp (filename)
  "Reads a Generalized Bezier Patch into a hashtable.
Special items: SIDES, DEGREE, LAYERS and CENTER.
All control points are under (SIDE COL ROW).
Note that for even-degree patches, the center point is listed
also under (SIDE LAYERS LAYERS) for all sides."
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
          (setf (gethash 'sides obj) n
                (gethash 'degree obj) d
                (gethash 'layers obj) l)
          (gbp-set-center obj center)
          (loop while (< i cp) do
               (when (>= col (- d row))
                 (incf side)
                 (when (>= side n)
                   (setf side 0)
                   (incf row))
                 (setf col row))
               (gbp-set-cp obj side col row (read-point))
               (incf i)
               (incf col)))))
    obj))

(defun write-gbp (obj filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let* ((n (gethash 'sides obj))
           (d (gethash 'degree obj))
           (l (gethash 'layers obj))
           (cp (1+ (* n (1+ (floor d 2)) l)))
           (i 1)
           (side 0)
           (col 0)
           (row 0))
      (format s "~a ~a~%" n d)
      (format s "~{~f~^ ~}~%" (gethash 'center obj))
      (loop while (< i cp) do
           (when (>= col (- d row))
             (incf side)
             (when (>= side n)
               (setf side 0)
               (incf row))
             (setf col row))
           (format s "~{~f~^ ~}~%" (gethash (list side col row) obj))
           (incf i)
           (incf col)))))

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


;;; Degree elevation

;;; see: Ball-Zheng: Degree elevation for n-sided surfaces (2001)

(defun zb-elevate-odd (obj)
  "As in Eq. (11). Only goes from odd to even degrees.
This is exactly the same algorithm as ours for GB patches."
  (let ((result (make-hash-table :test 'equal))
        (n (gethash 'sides obj))
        (d (1+ (gethash 'degree obj)))
        (l (gethash 'layers obj)))
    (setf (gethash 'sides result) n
          (gethash 'degree result) d
          (gethash 'layers result) l)
    (let ((inner-cps (iter (for i from 0 below n)
                           (collect (gethash (list i (1- l) (1- l)) obj)))))
      (gbp-set-center result (v* (reduce #'v+ inner-cps) (/ n))))
    (iter (for i from 0 below n)
          (gbp-set-cp result i 0 0 (gethash (list i 0 0) obj))
          (iter (for k from 1 below d)
                (for p = (v* (v+ (v* (gethash (list i (1- k) 0) obj) k)
                                 (v* (gethash (list i k 0) obj) (- d k)))
                             (/ d)))
                (gbp-set-cp result i k 0 p))
          (iter (for j from 1 below l)
                (iter (for k from 1 to l)
                      (for p = (v* (v+ (v* (gethash (list i (1- k) (1- j)) obj) k j)
                                       (v* (gethash (list i (1- k) j) obj) k (- d j))
                                       (v* (gethash (list i k (1- j)) obj) (- d k) j)
                                       (v* (gethash (list i k j) obj) (- d k) (- d j)))
                                   (/ 1 d d)))
                      (gbp-set-cp result i k j p))))
    result))

(defun zb-elevate-even (obj)
  (declare (ignore obj))
  (error "TODO"))

(defun zb-elevate (obj)
  (if (oddp (gethash 'degree obj))
      (zb-elevate-odd obj)
      (zb-elevate-even obj)))


;;; Tests

#+nil
(defparameter *obj* (read-gbp "/home/salvi/project/transfinite/models/cagd86.gbp"))
#+nil
(defparameter *obj* (read-gbp "/home/salvi/Dropbox/Shares/GrafGeo/EuroGraphics/6sided-2.gbp"))

#+nil
(let* ((*resolution* 30)
       (n (gethash 'sides *obj*))
       (points (zb-vertices n))
       (data (mapcar (lambda (u) (zb-eval *obj* u)) points)))
  (write-obj-indexed-mesh data (triangles n) "/tmp/proba.obj"))
