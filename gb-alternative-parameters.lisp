(in-package :cl-nurbs-tests)

(defvar *barycentric-dilation*)

(defun barycentric-d (l i)
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n)))
    (* (- 1 (elt l i-1) (elt l i))
       (- 1 (* (elt l i-2) (elt l i+1)
               *barycentric-dilation*)))))

#+nil
(defun barycentric-d (l i)
  "With warning for negative values"
  (let* ((n (length l))
         (i-2 (mod (- i 2) n))
         (i-1 (mod (- i 1) n))
         (i+1 (mod (+ i 1) n))
         (d1 (- 1 (elt l i-1) (elt l i)))
         (d (* (max 0 d1)
               (- 1 (* (elt l i-2) (elt l i+1)
                       *barycentric-dilation*)))))
    (when (< d (- *epsilon*))
      (warn "Negative d with l = ~a, i = ~a" l i))
    d))


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

(defun barycentric-s (l i)
  (let* ((n (length l))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (a (barycentric-d l i-1))
         (b (- 1 (barycentric-d l i+1)))
         (alpha (/ (+ a b) 2)))
    (+ (* (- 1 alpha) a) (* alpha b))))


;;; Helper functions
;;; (basically the same as in gb-central-cp-tests, but without autowp)

(defun deficiency-negative-p (n d &key (use-d t))
  (iter (for p in (vertices (points-from-angles (uniform-angles n))))
        (for def = (deficiency n d :position p :use-d use-d))
        (when (< def (- *epsilon*))
          (return def))))

(defun write-bernstein-blend (path n degree &key (use-d t) &allow-other-keys)
  (let ((fname (format nil "~a/~asided-deg~a.obj" path n degree))
        (points (points-from-angles (uniform-angles n))))
    (write-obj-indexed-mesh
     (iter (for p in (vertices points))
           (for def = (deficiency n degree :position p :use-d use-d))
           (when (< def (- *epsilon*))
             (warn "Negative deficiency: ~a" def))
           (collect (cons (if (< (abs def) *epsilon*) 0 def) p)))
     (triangles n) fname)))

(defun gb-blend-maximums (n degree filename &key (use-d t))
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
                            (if use-d
                                (if (< (+ di-1 di) *epsilon*)
                                    0.5
                                    (/ di-1 (+ di-1 di)))
                                (if (< (+ si (- 1 si-1)) *epsilon*)
                                    0.5
                                    (/ si (+ si (- 1 si-1))))))
                       (for beta =
                            (if use-d
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