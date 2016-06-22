(in-package :cl-nurbs-tests)

(defun bspline-basis (d knots i u)
  "I is the span."
  (let ((coeff (make-array (1+ d)))
        (left (make-array (1+ d)))
        (right (make-array (1+ d))))
    (setf (elt coeff 0) 1)
    (iter (for j from 1 to d)
          (setf (elt left j) (- u (elt knots (- i -1 j)))
                (elt right j) (- (elt knots (+ i j)) u))
          (iter (for saved first 0 then (* tmp (elt left (- j r))))
                (for r from 0 below j)
                (for tmp = (/ (elt coeff r)
                              (+ (elt right (1+ r)) (elt left (- j r)))))
                (setf (elt coeff r) (+ saved (* tmp (elt right (1+ r)))))
                (finally (setf (elt coeff j) saved))))
    coeff))

(defun uniform-bspline-basis (d n u)
  "N is the number of segments.
Assumes that U is in [0,1]."
  (let* ((u (max 0 (min u 1)))
         (knots (append (make-list d :initial-element 0)
                        (iter (for i from 0 to n)
                              (collect (/ i n)))
                        (make-list d :initial-element 1)))
         (upper-bound (position u knots :test #'<))
         (span (or (and upper-bound (- upper-bound 1)) (+ n d -1)))
         (basis (bspline-basis d knots span u)))
    (lambda (i)
      (if (or (< i (- span d)) (> i span))
          0
          (elt basis (+ (- i span) d))))))

(defun deficiency-bsp (n m &key (degree 3) (position 'center) (use-d t))
  (let* ((points (points-from-angles (uniform-angles n)))
         (p (case position
              (center '(0 0))
              (edge-center-mid
               (v* (v+ (first points) (second points)) 1/4))
              (vertex-center-mid
               (v* (first points) 1/2))
              (t position)))
         (l (barycentric-coordinates points p)))
    (- 1
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
             (for bsp-di = (uniform-bspline-basis degree (- m degree -1) di))
             (for bsp-si = (uniform-bspline-basis degree (- m degree -1) si))
             (for blf-sum = 0)
             (iter (for row from 0 below (ceiling m 2))
                   (iter (for col from 0 to m)
                         (for blend = (* (funcall bsp-di row) (funcall bsp-si col)))
                         (for mu = (cond ((and (< row 2) (< col 2)) alpha)
                                         ((and (< row 2) (> col (- m 2))) beta)
                                         ((or (< col row) (> col (- m row))) 0)
                                         ((or (= col row) (= col (- m row))) 1/2)
                                         ((< row 2) 1)
                                         (t 1)))
                         (incf blf-sum (* mu blend))))
             (sum blf-sum)))))

(defun write-bspline-blend (path n m &key (degree 3) (use-d t) &allow-other-keys)
  (let ((fname (format nil "~a/n~a-d~a-m~a.obj" path n degree m))
        (points (points-from-angles (uniform-angles n))))
    (write-obj-indexed-mesh
     (iter (for p in (vertices points))
           (for def = (deficiency-bsp n m :degree degree :position p :use-d use-d))
           (when (< def (- *epsilon*))
             (warn "Negative deficiency: ~a" def))
           (collect (cons (if (< (abs def) *epsilon*) 0 def) p)))
     (triangles n) fname)))
