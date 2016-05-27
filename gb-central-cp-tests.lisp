(in-package :cl-nurbs-tests)

(flet ((p (x) (round (* 100 x))))
  (iter (for d from 4 to 10)
        (format t "d = ~a~%" d)
        (for x = (iter (for k from 2 below (1- d))
                       (for b = (bernstein d k (/ k d)))
                       (format t "  ~a => ~a%~%" k (p b))
                       (sum b)))
        (format t "===> ~a% [atlag: ~a%]~%~%" (p x) (p (/ x (- d 3))))))

(flet ((p (x) (round (* 100 x))))
  (iter (for d from 4 to 10)
        (format t "d = ~a~%" d)
        (for x = (iter (for j from 2 below (1- d))
                       (for bj = (bernstein d j (/ j d)))
                       (sum (iter (for k from 2 below (1- d))
                                  (for b = (* bj (bernstein d k (/ k d))))
                                  (format t "  ~a,~a => ~a%~%" j k (p b))
                                  (sum b)))))
        (format t "===> ~a% [atlag: ~a%]~%~%" (p x) (p (/ x (- d 3) (- d 3))))))

;;; Blend function images:
(let* ((n 6)
       (d 8)
       (b (bernstein d (/ d 2) 1/2))
       (x (* (/ 4 n) (* b b)))
       (*resolution* 60)
       (*barycentric-type* 'wachspress)
       (*auto-wachspress-central-d* 0.58 #+nil(find-autowp-for-deficiency n d :target x))
       (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
  (format t "n = ~a, d = ~a~%" n d)
  (format t "4-sided central weight: ~a~%" (* b b))
  (format t "Center deficiency: ~a~%" x)
  (format t "Autowachspress weight: ~a~%" *auto-wachspress-central-d*)
  (write-bernstein-blend-autowp "/tmp" n d :density 0.05d0))

(defun write-bernstein-blend-autowp (path n degree &key (use-d t) &allow-other-keys)
  (let ((fname (format nil "~a/~asided-deg~a.obj" path n degree))
        (points (points-from-angles (uniform-angles n))))
    (write-obj-indexed-mesh
     (iter (for p in (vertices points))
           (for def = (deficiency-autowp n degree :position p :use-d use-d))
           (when (< def (- *epsilon*))
             (warn "Negative deficiency: ~a" def))
           (collect (cons (if (< (abs def) *epsilon*) 0 def) p)))
     (triangles n) fname)))

;;; Domain parameterization:
(let* ((n 5)
       (d 6)
       (b (bernstein d (/ d 2) 1/2))
       (x (* (/ 4 n) (* b b)))
       (*resolution* 60)
       (*barycentric-type* 'wachspress)
       (*auto-wachspress-central-d* (/ (- n 2) n) #+nil(find-autowp-for-deficiency n d :target x))
       (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
  (vectorized-distance-function-test (points-from-angles (uniform-angles n))
                                     (cons 'sd (make-list (1- n)))
                                     (format nil "/tmp/~asided-deg~a-domain3.ps" n d)
                                     :density 25 :color nil :resolution 0.001d0
                                     :distance-type 'bary-autowp))

;;; Barycentric lines
(defun write-barycentric-lines (path n &key (density 0.1))
  (labels ((transform (p) (list (+ (* (first p) 250) 250) (- 500 (* (second p) 250))))
           (write-poly (stream points)
             (format stream "~{~f ~}moveto~%~{~{~f ~}lineto~%~}closepath stroke~%"
                     (transform (first points)) (mapcar #'transform (rest points)))))
    (let ((fname (format nil "~a/barycentric-~asided.ps" path n)))
      (with-open-file (s fname :direction :output :if-exists :supersede)
        (format s "%!PS~%")
        (let ((points (points-from-angles (uniform-angles n))))
          (write-ps-indexed-mesh-projection
           (iter (for p in (vertices points))
                 (for l = (elt (barycentric-coordinates points p) 0))
                 (collect (cons l p)))
           (triangles n) s :transform #'transform :axis 0 :lines density)
          (write-poly s points))))))

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

;;; Find blend maximums
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
                       (for di = (barycentric-d-autowp l i))
                       (for di-1 = (barycentric-d-autowp l i-1))
                       (for di+1 = (barycentric-d-autowp l i+1))
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
                                                   ((or (= col row) (= col (- degree row))) *gb-diagonal-weight*)
                                                   ((< row 2) 1)
                                                   (t *gb-side-weight*)))
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

;;; Find non-negative intervals
(defun deficiency-autowp-negative-p (n d &key (use-d t))
  (iter (for p in (vertices (points-from-angles (uniform-angles n))))
        (for def = (deficiency-autowp n d :position p :use-d use-d))
        (when (< def (- *epsilon*))
          (return def))))
(defun find-autowp-negative-boundary (n d min max &key (iterations 100) (use-d t))
  (flet ((f (x)
           (let ((*auto-wachspress-central-d* x)
                 (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
             (or (deficiency-autowp-negative-p n d :use-d use-d) 1))))
    (bisection-search-root #'f min max iterations)))

(let* ((n 8)
       (d 8)
       (use-d t)
       (iter 10)
       (*resolution* 40)
       (*epsilon* 1.0d-5)
       (*barycentric-type* 'wachspress)
       (target-blend 0.5)
       (*gb-side-weight* 0.9)
       (*gb-diagonal-weight* (/ *side-weight* 2))
       (*auto-wachspress-central-d* (find-autowp-for-deficiency n d :target target-blend))
       (*auto-wachspress-weights* (make-list n :initial-element (/ (- n 2) n))))
;(deficiency-autowp n d :use-d use-d)
;(deficiency-autowp-negative-p n d :use-d use-d)
;(find-autowp-negative-boundary n d 0.5 0.7 :iterations iter :use-d use-d)
;(write-bernstein-blend-image-autowp "/tmp" n d :density 0.05)
(write-bernstein-blend-autowp "/tmp" n d :density 0.05 :use-d use-d)
;(gb-blend-maximums n d "/tmp/max.obj")
  )
