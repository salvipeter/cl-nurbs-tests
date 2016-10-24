(in-package :cl-nurbs-tests)

(defvar *barycentric-dilation* 0)

(defun barycentric-d-1minus (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (* (- 1 (elt l i-1) (elt l i))
       (- 1 (* (elt l i-2) (elt l i+1)
               *barycentric-dilation*)))))

(defun barycentric-d-peti (l i)
  (let ((n (length l)))
    (flet ((li (k) (elt l (mod (+ i k) n))))
      (- 1
         (* (+ (li 0) (li -1))
            (1+ (* *barycentric-dilation*
                   (- 1 (li -2) (li -1) (li 0) (li 1)))))))))

(defun barycentric-d-peti-1minus (l i)
  (let ((n (length l)))
    (flet ((li (k) (elt l (mod (+ i k) n))))
      (* (- 1 (li 0) (li -1))
         (- 1 (* *barycentric-dilation*
                 (+ (li 0) (li -1))
                 (- 1 (li -2) (li -1) (li 0) (li 1))))))))

(defun barycentric-d-peti-fullness (l i)
  (let ((n (length l)))
    (flet ((li (k) (elt l (mod (+ i k) n)))
           (sqr (x) (* x x)))
      (- 1 (li 0) (li -1)
         (* *barycentric-dilation*
            (+ (li 0) (li -1) (li -2) (li 1))
            (- 1 (li -2) (li -1) (li 0) (li 1)))))))

(defun barycentric-d-pisti-all-multiplication (l i)
  (let* ((n (length l))
         (m (floor (1- n) 2)))
    (flet ((li (k) (elt l (mod (+ i k) n))))
      (- 1 (li 0) (li -1)
         (* 2
            (iter (for j from (- m) to (- 2))
                  (sum (li j)))
            (iter (for j from 1 to (1- m))
                  (sum (li j))))
         (if (oddp n)
             (- (+ (* (li 1) (li (- m)))
                   (* (li -2) (li (1- m)))))
             0)
         (if (= n 5)
             (* (li 1) (li (- 2)))
             0)
         (* *barycentric-dilation* (li -2) (li 1)
            (- 1 (li -2) (li -1) (li 0) (li 1)))))))

(defun barycentric-d-pisti-all (l i)
  "Pisti version 2 (sum + sum)."
  (let* ((n (length l))
         (m (floor (1- n) 2)))
    (flet ((li (k) (elt l (mod (+ i k) n))))
      (- 1 (li 0) (li -1)
         (* 2
            (+ (* (li 1)
                  (iter (for j from (- m) to -2)
                        (sum (li j))))
               (* (li -2)
                  (iter (for j from 2 to (1- m))
                        (sum (li j))))))
         (if (oddp n)
             (- (+ (* (li 1) (li (- m)))
                   (* (li -2) (li (1- m)))))
             0)
         (if (= n 5)
             (* (li 1) (li -2))
             0)
         (* *barycentric-dilation* (li -2) (li 1)
            (- 1 (li -2) (li -1) (li 0) (li 1)))))))

(defun barycentric-d-1term (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (- 1 (elt l i-1) (elt l i)
       (* (elt l i-2) (elt l i+1)
          *barycentric-dilation*))))

(defun barycentric-d-tomi-pisti (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (- 1 (elt l i-1) (elt l i)
       (* *barycentric-dilation*
          (elt l i-2) (elt l i+1)
          (- 1 (elt l i-2) (elt l i-1) (elt l i) (elt l i+1))))))

(defun barycentric-d-tomi2 (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (- 1 (elt l i-1) (elt l i)
       (* *barycentric-dilation*
          (+ (* (elt l i-2) (elt l i))
             (* (elt l i-1) (elt l i+1))
             (* (elt l i-2) (elt l i+1)))))))

(defun barycentric-d-tomi (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (* (- 1 (elt l i-1) (elt l i))
       (- 1 (* (+ (* (elt l i-2) (elt l i))
                  (* (elt l i+1) (elt l i-1)))
               *barycentric-dilation*)))))

(defun barycentric-d-triangle (l i)
  (let* ((n (length l))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (+ (elt l i+1)
       (* *barycentric-dilation* (- 1 (elt l i+1))
          (elt l i-1) (elt l i) (elt l i+1)))))


;;; s_i similar to h_i (not used)

;; Legyen
;; $$s_i=(1-\alpha)h_{i-1}+\alpha(1-h_{i+1}),$$
;; ahol
;; $$\alpha=(h_{i-1}+(1-h_{i+1}))/2,$$
;; Mivel az $i$-edik oldalon
;; $$h_{i-1}+h_{i+1}=1,$$
;; ezert
;; $$s_i=h_{i-1}=1-h_{i+1}$$
;; valamint az $i-1$-es oldal kozeleben $s_i$ ugy viselkedik, mint $h_{i-1}$,
;; az $i+1$-es oldalon kozeleben pedig ugy, mint $h_{i+1}$.

(defun barycentric-s-alternative (l i)
  (let* ((n (length l))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (a (barycentric-d l i-1))
         (b (- 1 (barycentric-d l i+1)))
         (alpha (/ (+ a b) 2)))
    (+ (* (- 1 alpha) a) (* alpha b))))


;;; Helper functions
;;; (basically the same as in gb-central-cp-tests, but without autowp)

(defparameter *deficiency-function* #'deficiency)

(defun deficiency-negative-p (n d)
  (iter (for p in (vertices (points-from-angles (uniform-angles n))))
        (for def = (funcall *deficiency-function* n d :position p))
        (when (< def (- *epsilon*))
          (return def))))

(defun find-dilation-negative-boundary (n d min max &key (iterations 100))
  (flet ((f (x)
           (let ((*barycentric-dilation* x))
             (- (or (deficiency-negative-p n d) 1)))))
    (bisection-search-root #'f min max iterations)))

(defun deficiency-monotone-p (n d)
  (let ((points (points-from-angles (uniform-angles n))))
    (iter (for i from 0 to *resolution*)
          (for x = (/ i *resolution*))
          (for p = (affine-combine (elt points 0) x (elt points 1)))
          (iter (for j from 0 to *resolution*)
                (for y = (/ j *resolution*))
                (for q = (v* p y))
                (for defic1 first nil then defic)
                (for defic = (funcall *deficiency-function* n d :position q))
                (when (or (< defic (- *epsilon*)) (and defic1 (> defic defic1)))
                  (return-from deficiency-monotone-p nil)))))
  t)

(defun find-dilation-monotone-boundary (n d min max &key (iterations 100))
  (flet ((f (x)
           (let ((*barycentric-dilation* x))
             (if (deficiency-monotone-p n d) -1 1))))
    (bisection-search-root #'f min max iterations)))

(defun find-dilation-for-deficiency (n d &key (min -20.0) (max 20.0) (target 0.0) (iterations 100))
  (flet ((f (x)
           (let ((*barycentric-dilation* x))
             (- target (funcall *deficiency-function* n d)))))
    (bisection-search-root #'f min max iterations)))

(defun write-bernstein-blend (path n degree &key &allow-other-keys)
  (let ((fname (format nil "~a/~asided-deg~a.obj" path n degree))
        (points (points-from-angles (uniform-angles n))))
    (write-obj-indexed-mesh
     (iter (for p in (vertices points))
           (for def = (deficiency n degree :position p))
           (when (< def (- *epsilon*))
             (warn "Negative deficiency: ~a" def))
           (collect (cons (if (< (abs def) *epsilon*) 0 def) p)))
     (triangles n) fname)))

(defun gb-find-control-point (n d i j k)
  (let* ((l (ceiling d 2))
         (cp (1+ (* n (1+ (floor d 2)) l))))
    (do ((c 1 (1+ c)) (side 0) (col 0 (1+ col)) (row 0))
        ((= c cp))
      (when (>= col (- d row))
        (incf side)
        (when (>= side n)
          (setf side 0)
          (incf row))
        (setf col row))
      (let ((side-1 (mod (1- side) n))
            (side+1 (mod (1+ side) n)))
        (when (or (and (= i side) (= j col) (= k row))
                  (and (= i side-1) (= j (- d row)) (= k col))
                  (and (= i side+1) (= j row) (= k (- d col))))
          (return-from gb-find-control-point (1+ c))))))
  0)

(defun gb-write-control-net-data (n degree data filename)
  (let* ((l (ceiling degree 2))
         (cp (1+ (* n (1+ (floor degree 2)) l))))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (let ((p (assoc 'center data)))
        (format s "v~{ ~f~}~%" (append (third p) (list (second p)))))
      (do ((i 1 (1+ i)) (side 0) (col 0 (1+ col)) (row 0))
          ((= i cp))
        (when (>= col (- degree row))
          (incf side)
          (when (>= side n)
            (setf side 0)
            (incf row))
          (setf col row))
        (let ((p (assoc (list side col row) data :test #'equal)))
          (format s "v~{ ~f~}~%" (append (third p) (list (second p))))))
      (iter (for i from 0 below n)
            (iter (for j from 0 to (floor degree 2))
                  (iter (for k from 0 below (1- l))
                        (for a = (gb-find-control-point n degree i j k))
                        (for b = (gb-find-control-point n degree i (1+ j) k))
                        (for c = (gb-find-control-point n degree i (1+ j) (1+ k)))
                        (for d = (gb-find-control-point n degree i j (1+ k)))
                        (format s "f ~a ~a ~a ~a~%" a b c d))))
      (if (evenp degree)
          (iter (for i from 0 below n)
                (for i-1 = (mod (1- i) n))
                (for a = (gb-find-control-point n degree i (1- l) (1- l)))
                (for b = (gb-find-control-point n degree i l (1- l)))
                (for c = 1)
                (for d = (gb-find-control-point n degree i-1 l (1- l)))
                (format s "f ~a ~a ~a ~a~%" a b c d))
          (iter (for i from 0 below n)
                (for a = (gb-find-control-point n degree i (1- l) (1- l)))
                (for b = (gb-find-control-point n degree i l (1- l)))
                (for c = 1)
                (format s "f ~a ~a ~a~%" a b c))))))

(defun gb-blend-maximums (n degree filename)
  (let ((data '())
        weights)
    (flet ((merge-maximums ()
             (if data
                 (setf data
                       (mapcar (lambda (a b) (if (> (second a) (second b)) a b))
                               data weights))
                 (setf data weights)))
           (add-weight (key point value)
             (let ((datum (assoc key weights :test #'equal)))
               (if datum
                   (incf (second datum) value)
                   (push (list key value point) weights)))))
      (iter (with points = (points-from-angles (uniform-angles n)))
            (with layers = (ceiling degree 2))
            (for p in (vertices points))
            (for l = (barycentric-coordinates points p))
            (setf weights '())
            (for def =
                 (iter (for i from 0 below n)
                       (for i-1 = (mod (1- i) n))
                       (for i+1 = (mod (1+ i) n))
                       (for si = (barycentric-s l i))
                       (for si-1 = (barycentric-s l i-1))
                       (for si+1 = (barycentric-s l i+1))
                       (for di = (barycentric-d l i))
                       (for di-1 = (barycentric-d l i-1))
                       (for di+1 = (barycentric-d l i+1))
                       (for alpha =
                            (if *deficiency-use-d*
                                (if (< (+ di-1 di) *epsilon*)
                                    0.5
                                    (/ di-1 (+ di-1 di)))
                                (if (< (+ si (- 1 si-1)) *epsilon*)
                                    0.5
                                    (/ si (+ si (- 1 si-1))))))
                       (for beta =
                            (if *deficiency-use-d*
                                (if (< (+ di+1 di) *epsilon*)
                                    0.5
                                    (/ di+1 (+ di+1 di)))
                                (if (< (+ (- 1 si) si+1) *epsilon*)
                                    0.5
                                    (/ (- 1 si) (+ (- 1 si) si+1)))))
                       (for blf-sum = 0)
                       (iter (for row from 0 below layers)
                             (iter (for col from 0 to degree)
                                   (for blend = (* (bernstein degree row di)
                                                   (bernstein degree col si)))
                                   (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                                   ((and (< row 2) (> col (- degree 2))) beta)
                                                   ((or (< col row) (> col (- degree row))) 0)
                                                   ((or (= col row) (= col (- degree row))) 1/2)
                                                   ((< row 2) 1)
                                                   (t 1)))
                                   (add-weight (list i col row) p (* mu blend))
                                   (when (< col layers)
                                     (add-weight (list i-1 (- degree row) col) p (* mu blend)))
                                   (when (< (- degree col) layers)
                                     (add-weight (list i+1 row (- degree col)) p (* mu blend)))
                                   (incf blf-sum (* mu blend)))) 
                       (sum blf-sum)))
            (add-weight 'center p (- 1 def))
            (merge-maximums)))
    (gb-write-control-net-data n degree data filename)))

(defun barycentric-pair-multi-max (n)
  (let ((points (points-from-angles (uniform-angles n))))
    (iter (for p in (vertices points))
          (for l = (barycentric-coordinates points p))
          (maximize (* (elt l 0) (elt l 3))))))

;;; Maximum dilation values (truncated to integers):
;; |----+----+----+----+----+----|
;; |  5 |  6 |  7 |  8 |  9 | 10 |
;; |----+----+----+----+----+----|
;; | 21 | 35 | 45 | 51 | 54 | 56 |
;; |----+----+----+----+----+----|

;;; Maximum non-negative values (truncated again):
;; |---+---+--------------+----------------|
;; | n | d | max dilation | min deficiency |
;; |---+---+--------------+----------------|
;; | 5 | 5 |            2 |             5% |
;; | 5 | 6 |            3 |             5% |
;; | 5 | 7 |            2 |            10% |
;; | 5 | 8 |            3 |             5% |
;; | 6 | 5 |            4 |            10% |
;; | 6 | 6 |            6 |             6% |
;; | 6 | 7 |            5 |             8% |
;; | 6 | 8 |            6 |             8% |
;; | 7 | 5 |            7 |             8% |
;; | 7 | 6 |            9 |            10% |
;; | 7 | 7 |            8 |            10% |
;; | 7 | 8 |           10 |             4% |
;; | 8 | 5 |           11 |             1% |
;; | 8 | 6 |           13 |             9% |
;; | 8 | 7 |           12 |             7% |
;; | 8 | 8 |           13 |            16% |
;; |---+---+--------------+----------------|

;;; Some specific values:
;; |---+---+--------+--------|
;; | n | d |    10% |    20% |
;; |---+---+--------+--------|
;; | 5 | 5 |  1.581 |  0.756 |
;; | 5 | 6 |  2.676 |  1.959 |
;; | 5 | 7 |  2.013 |  1.319 |
;; | 5 | 8 |  2.716 |  2.090 |
;; | 6 | 5 |  4.020 |  3.042 |
;; | 6 | 6 |  5.649 |  4.791 |
;; | 6 | 7 |  4.835 |  4.005 |
;; | 6 | 8 |  5.881 |  5.128 |
;; | 7 | 5 |  6.761 |  5.604 |
;; | 7 | 6 |  9.027 |  8.005 |
;; | 7 | 7 |  8.039 |  7.049 |
;; | 7 | 8 |  9.496 |  8.597 |
;; | 8 | 5 |  9.832 |  8.477 |
;; | 8 | 6 | 12.839 | 11.634 |
;; | 8 | 7 | 11.656 | 10.490 |
;; | 8 | 8 | 13.596 | 12.533 |
;; |---+---+--------+--------|
