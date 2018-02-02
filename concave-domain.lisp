(in-package :cl-nurbs-tests)

(defun concavep (points i)
  "Assumes that POINTS is given in positive order, i.e. `matter' is on the left side."
  (let* ((n (length points))
         (im (mod (1- i) n))
         (ip (mod (1+ i) n))
         (u (vnormalize (v- (elt points ip) (elt points i))))
         (v (vnormalize (v- (elt points im) (elt points i)))))
    (<= (- (* (elt u 0) (elt v 1))
           (* (elt u 1) (elt v 0)))
       0)))

(defun snip-triangle (points i)
  (append (subseq points 0 i)
          (subseq points (1+ i))))

(defun mesh-triangle (points i)
  (let* ((n (length points))
         (a (elt points (mod (1- i) n)))
         (b (elt points i))
         (c (elt points (mod (1+ i) n))))
    (iter (for i from 1 to *resolution*)
          (for u = (/ i *resolution*))
          (for p = (affine-combine b u a))
          (for q = (affine-combine b u c))
          (for last-row first (list (list b)) then (nreverse row))
          (for row = '())
          (appending
           (iter (for j from 0 below i)
                 (for tri = (list (affine-combine p (/ j i) q)
                                  (affine-combine p (/ (1+ j) i) q)
                                  (if (= j 0)
                                      (elt (elt last-row j) 0)
                                      (elt (elt last-row (1- j)) 1))))
                 (push tri row)
                 (collect tri)
                 (unless (= j 0)
                   (collect (list (affine-combine p (/ j i) q)
                                  (elt (elt last-row (1- j)) 1)
                                  (elt (elt last-row (1- j)) 0)))))))))

(defun segment-segment-intersection (a b)
  "Returns NIL if there is no intersection point."
  (destructuring-bind (q1 q2) a
    (let* ((p (first b))
           (v (v- (second b) (first b)))
           (i (if (> (abs (first v)) (abs (second v))) 0 1)))
      (macrolet ((x (var) `(elt ,var i))
                 (y (var) `(elt ,var (- 1 i))))
        (let ((denom (- (* (x v) (- (y q2) (y q1)))
                        (* (y v) (- (x q2) (x q1))))))
          (unless (< (abs denom) *epsilon*)
            (let ((u (/ (+ (* (x v) (- (y p) (y q1)))
                           (* (y v) (- (x q1) (x p))))
                        denom)))
              (unless (or (< u (- *epsilon*)) (> u (1+ *epsilon*)))
                (let* ((u1 (min (max u 0) 1))
                       (u2 (/ (+ (- (x q1) (x p))
                                 (* (- (x q2) (x q1)) u1))
                              (x v))))
                  (unless (or (< u2 (- *epsilon*)) (> u2 (1+ *epsilon*)))
                    (v+ q1 (v* (v- q2 q1) u1))))))))))))

(defun snippablep (points i)
  "Returns T if the segment from i-1 to i+1 does not intersect the polygon"
  (let* ((n (length points))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n)))
    (iter (for j from 0 below n)
          (for j+1 = (mod (1+ j) n))
          (when (and (/= j (mod (- i 2) n)) (/= j i-1) (/= j i) (/= j i+1)
                     (segment-segment-intersection
                      (list (elt points i-1) (elt points i+1))
                      (list (elt points j) (elt points j+1))))
            (return-from snippablep nil))))
  t)

(defun mesh-concave (points)
  (iter (with points = (copy-list points))
        (while (>= (length points) 3))
        (for i = (or (iter (with n = (length points))
                           (for j from 0 below n)
                           (when (and (not (concavep points j))
                                      (or (concavep points (mod (1- j) n))
                                          (concavep points (mod (1+ j) n)))
                                      (snippablep points j))
                             (return j)))
                     (iter (with n = (length points))
                           (for j from 0 below n)
                           (when (snippablep points j)
                             (return j)))))
        (appending (mesh-triangle points i))
        (setf points (snip-triangle points i))))

#+nil
(defun mesh-concave (points)
  "Just for testing convex(!) domains."
  (let ((vertices (vertices points))
        (triangles (triangles (length points))))
    (mapcar (lambda (tri)
              (mapcar (lambda (i)
                        (elt vertices i))
                      tri))
            triangles)))

#+nil
(labels ((to3d-pt (p) (append p '(0)))
         (to3d-tri (tri) (mapcar #'to3d-pt tri))
         (to3d (tris) (mapcar #'to3d-tri tris)))
  (let ((*resolution* 20)
        (points '((2 5) (3 3) (5 3) (6 5) (8 3) (7 0) (1 0) (0 3))))
    (write-stl (to3d (mesh-concave points)) "/tmp/proba.stl" :ascii t)))

(defun eval-on-concave-domain (points fn)
  (let ((mesh (mesh-concave points)))
    (mapcar (lambda (tri)
              (mapcar (lambda (p)
                        (cons (funcall fn p) p))
                      tri))
            mesh)))

(defun eval3d-on-concave-domain (points fn)
  (let ((mesh (mesh-concave points)))
    (mapcar (lambda (tri) (mapcar fn tri)) mesh)))

(defun eval-over-domain (points fn)
  (map 'vector (lambda (p) (cons (funcall fn p) p)) points))

#+nil
(let ((*resolution* 20)
      (points '((2 5) (3 3) (5 3) (6 5) (8 3) (7 0) (1 0) (0 3))))
  (flet ((fn (p)
           (mean-value points (cons 1 (make-list (1- (length points)) :initial-element 0)) p)))
    (write-stl (eval-on-concave-domain points #'fn) "/tmp/proba.stl" :ascii t)))

(defun mean-value-coordinates (points p)
  (let* ((vectors (mapcar (lambda (x) (v- p x)) points))
         (lengths (mapcar #'vlength vectors))
         (n (length points)))
    (labels ((inc (i) (mod (1+ i) n))
             (dec (i) (mod (1- i) n))
             (area (i)                  ; signed area = det(si,si+1)/2
               (let ((si (elt vectors i))
                     (si+1 (elt vectors (inc i))))
                 (/ (- (* (elt si 0) (elt si+1 1))
                       (* (elt si 1) (elt si+1 0)))
                    2.0d0))))
      (let* ((w (iter (for i from 0 below n)
                      (for Ai = (area i))
                      (for Ai-1 = (area (dec i)))
                      (for Di = (scalar-product (elt vectors (inc i)) (elt vectors i)))
                      (for Di-1 = (scalar-product (elt vectors i) (elt vectors (dec i))))
                      (for ri-1 = (elt lengths (dec i)))
                      (for ri = (elt lengths i))
                      (for ri+1 = (elt lengths (inc i)))
                      (when (< (abs ri) *epsilon*)
                        (let ((result (make-list n :initial-element 0)))
                          (setf (elt result i) 1)
                          (return-from mean-value-coordinates result)))
                      (when (and (< (abs Ai) *epsilon*)
                                 (< Di (- *epsilon*)))
                        (let ((result (make-list n :initial-element 0)))
                          (setf (elt result i) (/ ri+1 (+ ri ri+1))
                                (elt result (inc i)) (/ ri (+ ri ri+1)))
                          (return-from mean-value-coordinates result)))
                      (collect (+ (if (> (abs Ai-1) *epsilon*)
                                      (/ (- ri-1 (/ Di-1 ri)) Ai-1)
                                      0)
                                  (if (> (abs Ai) *epsilon*)
                                      (/ (- ri+1 (/ Di ri)) Ai)
                                      0)))))
             (wsum (reduce #'+ w)))
        (mapcar (lambda (wi) (/ wi wsum)) w)))))

(defun mean-kato (points i)
  (lambda (p)
    (let* ((l (mean-value-coordinates points p))
           (d (iter (with n = (length l))
                    (for i from 0 below n)
                    (for i-1 = (mod (1- i) n))
                    (collect (- 1 (elt l i-1) (elt l i))))))
      ;; (assert (every (lambda (x) (>= x -1.0d-6)) d))
      (ribbon-blend d i))))

(defun harmonic-kato (map points i)
  (lambda (p)
    (let ((l (harmonic:harmonic-coordinates map p)))
      (when (member nil l)              ; kutykurutty
        (setf l (mean-value-coordinates points p)))
      (let ((d (iter (with n = (length l))
                     (for i from 0 below n)
                     (for i-1 = (mod (1- i) n))
                     (collect (- 1 (elt l i-1) (elt l i))))))
        (ribbon-blend d i)))))

(defun mean-sweep (points i)
  (lambda (p)
    (let* ((n (length points))
           (i-1 (mod (1- i) n))
           (i+1 (mod (1+ i) n))
           (l (mean-value-coordinates points p))
           (d (iter (for i from 0 below n)
                    (for i-1 = (mod (1- i) n))
                    (collect (- 1 (elt l i-1) (elt l i))))))
      (safe-/ (elt d i-1) (+ (elt d i-1) (elt d i+1))))))

(defun mean-bernstein (points i j k)
  "Bivariate cubic Bernstein function B^3_j * B^3_k using the above parameters."
  (lambda (p)
    (let* ((n (length points))
           (i-1 (mod (1- i) n))
           (i+1 (mod (1+ i) n))
           (l (mean-value-coordinates points p))
           (d (iter (for i from 0 below n)
                    (for i-1 = (mod (1- i) n))
                    (collect (- 1 (elt l i-1) (elt l i)))))
           (si (safe-/ (elt d i-1) (+ (elt d i-1) (elt d i+1))))
           (di (elt d i)))
      (* (bernstein 3 j si) (bernstein 3 k di)))))

#+nil
(let ((*resolution* 50)
      ;; (points '((0 0) (0 7) (1 7) (1 4) (2 4) (2 7) (3 7) (3 4) (4 4) (4 7) (5 7)
      ;;           (5 0) (4 0) (4 3) (3 3) (3 0) (2 0) (2 3) (1 3) (1 0)))
      ;; (points '((0 0) (0 7) (1 7) (1 4) (2 4) (2 7) (3 7) (3 0) (2 0) (2 3) (1 3) (1 0)))
      ;; (points '((2 5) (3 3) (5 3) (6 5) (8 3) (7 0) (1 0) (0 3)))
      (points '((1 2) (1 1) (2 1) (2 0) (0 0) (0 2)))
      )
  (dotimes (j 4)
    (dotimes (k 2)
      (write-stl (eval-on-concave-domain points (mean-bernstein points 1 j k))
                 (format nil "/tmp/b~a~a.stl" j k) :ascii t))))

#+nil
(let* ((points '((0 0) (6 0) (6 6) (4 6) (4 4) (2 4) (2 6) (0 6))) ; U
       (points (scale-to-unit points)))
  (destructuring-bind (vertices triangles)
      (shewchuk-triangle:mesh points 0.0001)
    (harmonic:with-harmonic-coordinates (h points :levels 10)
      (flet ((foo (p)
               (+ (funcall (harmonic-kato h points 1) p)
                  (funcall (harmonic-kato h points 2) p))))
        (write-obj-indexed-mesh (eval-over-domain vertices #'foo)
                                triangles "/tmp/proba.obj")))))

;;; Maxima computation:
;;; B(n,k,u) := binomial(n,k)*u^k*(1-u)^(n-k);
;;; surf : sum(sum(C(i,j,coord)*B(3,i,s)*B(1,j,d),j,0,1),i,0,3)$
;;; eq : [at(surf, coord=x) = u, at(surf, coord=y) = v]$
;;; eq_s : expand(eliminate(eq, [d]))$
;;; for i:0 thru 6 do (z: coeff(eq_s, s, i), display(i,string(z)))$
;;; for i:0 thru 3 do (z: coeff(expand(at(surf, coord=x) = u), d, i), display(i,string(z)))$
;;; (lisp-szeru iras: save("filename", expr)$

(defun reverse-cubic-parameters (cpts p)
  "CPTS is ((B00 B10 B20 B30) (B01 B11 B21 B31)).
Returns parameters for which the Bezier patch given by CPTS evaluates to P."
  (macrolet ((c (i j coord)
               (if (eq coord 'x)
                   `(first (elt (elt cpts ,j) ,i))
                   `(second (elt (elt cpts ,j) ,i)))))
    (let* ((u (first p))
           (v (second p))
           (p0 (+ (* (- (c 0 1 x)) v) (* (c 0 0 x) v) (* (c 0 1 y) u) (* (- (c 0 0 y)) u) (* (- (c 0 0 x)) (c 0 1 y)) (* (c 0 0 y) (c 0 1 x))))
           (p1 (+ (* (- 3) (c 1 1 x) v) (* 3 (c 1 0 x) v) (* 3 (c 0 1 x) v) (* (- 3) (c 0 0 x) v) (* 3 (c 1 1 y) u) (* (- 3) (c 1 0 y) u) (* (- 3) (c 0 1 y) u) (* 3 (c 0 0 y) u) (* (- 3) (c 0 0 x) (c 1 1 y)) (* 3 (c 0 0 y) (c 1 1 x)) (* 3 (c 0 1 x) (c 1 0 y)) (* (- 3) (c 0 1 y) (c 1 0 x)) (* 6 (c 0 0 x) (c 0 1 y)) (* (- 6) (c 0 0 y) (c 0 1 x))))
           (p2 (+ (* (- 3) (c 2 1 x) v) (* 3 (c 2 0 x) v) (* 6 (c 1 1 x) v) (* (- 6) (c 1 0 x) v) (* (- 3) (c 0 1 x) v) (* 3 (c 0 0 x) v) (* 3 (c 2 1 y) u) (* (- 3) (c 2 0 y) u) (* (- 6) (c 1 1 y) u) (* 6 (c 1 0 y) u) (* 3 (c 0 1 y) u) (* (- 3) (c 0 0 y) u) (* (- 3) (c 0 0 x) (c 2 1 y)) (* 3 (c 0 0 y) (c 2 1 x)) (* 3 (c 0 1 x) (c 2 0 y)) (* (- 3) (c 0 1 y) (c 2 0 x)) (* (- 9) (c 1 0 x) (c 1 1 y)) (* 15 (c 0 0 x) (c 1 1 y)) (* 9 (c 1 0 y) (c 1 1 x)) (* (- 15) (c 0 0 y) (c 1 1 x)) (* (- 15) (c 0 1 x) (c 1 0 y)) (* 15 (c 0 1 y) (c 1 0 x)) (* (- 15) (c 0 0 x) (c 0 1 y)) (* 15 (c 0 0 y) (c 0 1 x))))
           (p3 (+ (* (- (c 3 1 x)) v) (* (c 3 0 x) v) (* 3 (c 2 1 x) v) (* (- 3) (c 2 0 x) v) (* (- 3) (c 1 1 x) v) (* 3 (c 1 0 x) v) (* (c 0 1 x) v) (* (- (c 0 0 x)) v) (* (c 3 1 y) u) (* (- (c 3 0 y)) u) (* (- 3) (c 2 1 y) u) (* 3 (c 2 0 y) u) (* 3 (c 1 1 y) u) (* (- 3) (c 1 0 y) u) (* (- (c 0 1 y)) u) (* (c 0 0 y) u) (* (- (c 0 0 x)) (c 3 1 y)) (* (c 0 0 y) (c 3 1 x)) (* (c 0 1 x) (c 3 0 y)) (* (- (c 0 1 y)) (c 3 0 x)) (* (- 9) (c 1 0 x) (c 2 1 y)) (* 12 (c 0 0 x) (c 2 1 y)) (* 9 (c 1 0 y) (c 2 1 x)) (* (- 12) (c 0 0 y) (c 2 1 x)) (* 9 (c 1 1 x) (c 2 0 y)) (* (- 12) (c 0 1 x) (c 2 0 y)) (* (- 9) (c 1 1 y) (c 2 0 x)) (* 12 (c 0 1 y) (c 2 0 x)) (* 36 (c 1 0 x) (c 1 1 y)) (* (- 30) (c 0 0 x) (c 1 1 y)) (* (- 36) (c 1 0 y) (c 1 1 x)) (* 30 (c 0 0 y) (c 1 1 x)) (* 30 (c 0 1 x) (c 1 0 y)) (* (- 30) (c 0 1 y) (c 1 0 x)) (* 20 (c 0 0 x) (c 0 1 y)) (* (- 20) (c 0 0 y) (c 0 1 x))))
           (p4 (+ (* (- 3) (c 1 0 x) (c 3 1 y)) (* 3 (c 0 0 x) (c 3 1 y)) (* 3 (c 1 0 y) (c 3 1 x)) (* (- 3) (c 0 0 y) (c 3 1 x)) (* 3 (c 1 1 x) (c 3 0 y)) (* (- 3) (c 0 1 x) (c 3 0 y)) (* (- 3) (c 1 1 y) (c 3 0 x)) (* 3 (c 0 1 y) (c 3 0 x)) (* (- 9) (c 2 0 x) (c 2 1 y)) (* 27 (c 1 0 x) (c 2 1 y)) (* (- 18) (c 0 0 x) (c 2 1 y)) (* 9 (c 2 0 y) (c 2 1 x)) (* (- 27) (c 1 0 y) (c 2 1 x)) (* 18 (c 0 0 y) (c 2 1 x)) (* (- 27) (c 1 1 x) (c 2 0 y)) (* 18 (c 0 1 x) (c 2 0 y)) (* 27 (c 1 1 y) (c 2 0 x)) (* (- 18) (c 0 1 y) (c 2 0 x)) (* (- 54) (c 1 0 x) (c 1 1 y)) (* 30 (c 0 0 x) (c 1 1 y)) (* 54 (c 1 0 y) (c 1 1 x)) (* (- 30) (c 0 0 y) (c 1 1 x)) (* (- 30) (c 0 1 x) (c 1 0 y)) (* 30 (c 0 1 y) (c 1 0 x)) (* (- 15) (c 0 0 x) (c 0 1 y)) (* 15 (c 0 0 y) (c 0 1 x))))
           (p5 (+ (* (- 3) (c 2 0 x) (c 3 1 y)) (* 6 (c 1 0 x) (c 3 1 y)) (* (- 3) (c 0 0 x) (c 3 1 y)) (* 3 (c 2 0 y) (c 3 1 x)) (* (- 6) (c 1 0 y) (c 3 1 x)) (* 3 (c 0 0 y) (c 3 1 x)) (* 3 (c 2 1 x) (c 3 0 y)) (* (- 6) (c 1 1 x) (c 3 0 y)) (* 3 (c 0 1 x) (c 3 0 y)) (* (- 3) (c 2 1 y) (c 3 0 x)) (* 6 (c 1 1 y) (c 3 0 x)) (* (- 3) (c 0 1 y) (c 3 0 x)) (* 18 (c 2 0 x) (c 2 1 y)) (* (- 27) (c 1 0 x) (c 2 1 y)) (* 12 (c 0 0 x) (c 2 1 y)) (* (- 18) (c 2 0 y) (c 2 1 x)) (* 27 (c 1 0 y) (c 2 1 x)) (* (- 12) (c 0 0 y) (c 2 1 x)) (* 27 (c 1 1 x) (c 2 0 y)) (* (- 12) (c 0 1 x) (c 2 0 y)) (* (- 27) (c 1 1 y) (c 2 0 x)) (* 12 (c 0 1 y) (c 2 0 x)) (* 36 (c 1 0 x) (c 1 1 y)) (* (- 15) (c 0 0 x) (c 1 1 y)) (* (- 36) (c 1 0 y) (c 1 1 x)) (* 15 (c 0 0 y) (c 1 1 x)) (* 15 (c 0 1 x) (c 1 0 y)) (* (- 15) (c 0 1 y) (c 1 0 x)) (* 6 (c 0 0 x) (c 0 1 y)) (* (- 6) (c 0 0 y) (c 0 1 x))))
           (p6 (+ (* (- (c 3 0 x)) (c 3 1 y)) (* 3 (c 2 0 x) (c 3 1 y)) (* (- 3) (c 1 0 x) (c 3 1 y)) (* (c 0 0 x) (c 3 1 y)) (* (c 3 0 y) (c 3 1 x)) (* (- 3) (c 2 0 y) (c 3 1 x)) (* 3 (c 1 0 y) (c 3 1 x)) (* (- (c 0 0 y)) (c 3 1 x)) (* (- 3) (c 2 1 x) (c 3 0 y)) (* 3 (c 1 1 x) (c 3 0 y)) (* (- (c 0 1 x)) (c 3 0 y)) (* 3 (c 2 1 y) (c 3 0 x)) (* (- 3) (c 1 1 y) (c 3 0 x)) (* (c 0 1 y) (c 3 0 x)) (* (- 9) (c 2 0 x) (c 2 1 y)) (* 9 (c 1 0 x) (c 2 1 y)) (* (- 3) (c 0 0 x) (c 2 1 y)) (* 9 (c 2 0 y) (c 2 1 x)) (* (- 9) (c 1 0 y) (c 2 1 x)) (* 3 (c 0 0 y) (c 2 1 x)) (* (- 9) (c 1 1 x) (c 2 0 y)) (* 3 (c 0 1 x) (c 2 0 y)) (* 9 (c 1 1 y) (c 2 0 x)) (* (- 3) (c 0 1 y) (c 2 0 x)) (* (- 9) (c 1 0 x) (c 1 1 y)) (* 3 (c 0 0 x) (c 1 1 y)) (* 9 (c 1 0 y) (c 1 1 x)) (* (- 3) (c 0 0 y) (c 1 1 x)) (* (- 3) (c 0 1 x) (c 1 0 y)) (* 3 (c 0 1 y) (c 1 0 x)) (* (- (c 0 0 x)) (c 0 1 y)) (* (c 0 0 y) (c 0 1 x)))))
      (let ((s (nth-degree-solver (list p6 p5 p4 p3 p2 p1 p0) :min 0 :max 1)))
        (let ((denom (+ (* -1 (c 0 0 x)) (c 0 1 x) (* 3 (c 0 0 x) s) (* -3 (c 0 1 x) s) (* -3 (c 1 0 x) s) (* 3 (c 1 1 x) s) (* -3 (c 0 0 x) (expt s 2)) (* 3 (c 0 1 x) (expt s 2)) (* 6 (c 1 0 x) (expt s 2)) (* -6 (c 1 1 x) (expt s 2)) (* -3 (c 2 0 x) (expt s 2)) (* 3 (c 2 1 x) (expt s 2)) (* (c 0 0 x) (expt s 3)) (* -1 (c 0 1 x) (expt s 3)) (* -3 (c 1 0 x) (expt s 3)) (* 3 (c 1 1 x) (expt s 3)) (* 3 (c 2 0 x) (expt s 3)) (* -3 (c 2 1 x) (expt s 3)) (* -1 (c 3 0 x) (expt s 3)) (* (c 3 1 x) (expt s 3)))))
          (if (< (abs denom) *epsilon*)
              (list s (/ (- v (+ (c 0 0 y) (* -3 (c 0 0 y) s) (* 3 (c 1 0 y) s) (* 3 (c 0 0 y) (expt s 2)) (* -6 (c 1 0 y) (expt s 2)) (* 3 (c 2 0 y) (expt s 2)) (* -1 (c 0 0 y) (expt s 3)) (* 3 (c 1 0 y) (expt s 3)) (* -3 (c 2 0 y) (expt s 3)) (* (c 3 0 y) (expt s 3))))
                         (+ (* -1 (c 0 0 y)) (c 0 1 y) (* 3 (c 0 0 y) s) (* -3 (c 0 1 y) s) (* -3 (c 1 0 y) s) (* 3 (c 1 1 y) s) (* -3 (c 0 0 y) (expt s 2)) (* 3 (c 0 1 y) (expt s 2)) (* 6 (c 1 0 y) (expt s 2)) (* -6 (c 1 1 y) (expt s 2)) (* -3 (c 2 0 y) (expt s 2)) (* 3 (c 2 1 y) (expt s 2)) (* (c 0 0 y) (expt s 3)) (* -1 (c 0 1 y) (expt s 3)) (* -3 (c 1 0 y) (expt s 3)) (* 3 (c 1 1 y) (expt s 3)) (* 3 (c 2 0 y) (expt s 3)) (* -3 (c 2 1 y) (expt s 3)) (* -1 (c 3 0 y) (expt s 3)) (* (c 3 1 y) (expt s 3)))))
              (list s (/ (- u (+ (c 0 0 x) (* -3 (c 0 0 x) s) (* 3 (c 1 0 x) s) (* 3 (c 0 0 x) (expt s 2)) (* -6 (c 1 0 x) (expt s 2)) (* 3 (c 2 0 x) (expt s 2)) (* -1 (c 0 0 x) (expt s 3)) (* 3 (c 1 0 x) (expt s 3)) (* -3 (c 2 0 x) (expt s 3)) (* (c 3 0 x) (expt s 3))))
                         denom))))))))

;;; (reverse-cubic-parameters '(((0 0) (1 0) (2 0) (3 0)) ((0 3) (1 3) (2 3) (3 3))) '(2 4))

(defun reverse-bilinear-parameters (cpts p)
  "CPTS is ((B00 B10) (B01 B11)).
Returns parameters for which the Bezier patch given by CPTS evaluates to P."
  (macrolet ((c (i j coord)
               (if (eq coord 'x)
                   `(first (elt (elt cpts ,j) ,i))
                   `(second (elt (elt cpts ,j) ,i)))))
    (let* ((u (first p))
           (v (second p))
           (p0 (+ (* (- (c 0 1 x)) v) (* (c 0 0 x) v) (* (c 0 1 y) u) (* (- (c 0 0 y)) u) (* (- (c 0 0 x)) (c 0 1 y)) (* (c 0 0 y) (c 0 1 x))))
           (p1 (+ (* (- (c 1 1 x)) v) (* (c 1 0 x) v) (* (c 0 1 x) v) (* (- (c 0 0 x)) v) (* (c 1 1 y) u) (* (- (c 1 0 y)) u) (* (- (c 0 1 y)) u) (* (c 0 0 y) u) (* (- (c 0 0 x)) (c 1 1 y)) (* (c 0 0 y) (c 1 1 x)) (* (c 0 1 x) (c 1 0 y)) (* (- (c 0 1 y)) (c 1 0 x)) (* 2 (c 0 0 x) (c 0 1 y)) (* (- 2) (c 0 0 y) (c 0 1 x))))
           (p2 (+ (* (- (c 1 0 x)) (c 1 1 y)) (* (c 0 0 x) (c 1 1 y)) (* (c 1 0 y) (c 1 1 x)) (* (- (c 0 0 y)) (c 1 1 x)) (* (- (c 0 1 x)) (c 1 0 y)) (* (c 0 1 y) (c 1 0 x)) (* (- (c 0 0 x)) (c 0 1 y)) (* (c 0 0 y) (c 0 1 x)))))
      (flet ((try (s x u)
               (when s
                 (let ((denom (+ (* (+ (elt (elt (elt cpts 1) 1) x) (- (elt (elt (elt cpts 0) 1) x)) (- (elt (elt (elt cpts 1) 0) x)) (elt (elt (elt cpts 0) 0) x)) s) (elt (elt (elt cpts 1) 0) x) (- (elt (elt (elt cpts 0) 0) x)))))
                   (when (> (abs denom) *epsilon*)
                     (list (list s (/ (+ u (- (* (elt (elt (elt cpts 0) 1) x) s)) (* (elt (elt (elt cpts 0) 0) x) s) (- (* (elt (elt (elt cpts 0) 0) x))))
                                      denom))))))))
        (multiple-value-bind (s1 s2) (second-degree-solver p2 p1 p0)
          (let ((results (append (try s1 0 u) (try s1 1 v) (try s2 0 u) (try s2 1 v))))
            (iter (for sd in results)
                  (finding sd minimizing (vlength sd)))))))))

;;; (reverse-bilinear-parameters '(((0 0) (1 0)) ((0 2) (2 3))) '(1.7 0.3))

(defun bezier-surface (cpts uv)
  (iter (with p = '(0 0 0))
        (with m = (1- (length cpts)))
        (with n = (1- (length (first cpts))))
        (for j from 0 to m)
        (iter (for i from 0 to n)
              (for cp = (elt (elt cpts j) i))
              (setf p (v+ p (v* cp (bernstein n i (first uv)) (bernstein m j (second uv))))))
        (finally (return p))))

(defun read-gbp (filename)
  (flet ((read-point (s)
           (let* ((x (read s))
                  (y (read s))
                  (z (read s)))
             (list x y z))))
    (with-open-file (s filename)
      (let* ((n (read s))
             (d (read s))
             (l (floor (1+ d) 2))
             (cp (1+ (* n (1+ (floor d 2)) l)))
             (center (read-point s)))
        (declare (ignore center))
        (iter (for i from 1 below cp)
              (for col first 0 then (1+ col))
              (with side = 0)
              (with row = 0)
              (when (>= col (- d row))
                (when (>= (incf side) n)
                  (setf side 0)
                  (incf row))
                (setf col row))
              (let ((p (read-point s)))
                (collect (list side col row p))
                (when (< col l)
                  (collect (list (mod (1- side) n) (- d row) col p)))
                (when (< (- d col) l)
                  (collect (list (mod (1+ side) n) row (- d col) p)))))))))

(defun load-ribbons (filename)
  (let ((gbp (read-gbp filename)))
    (flet ((c (i j k)
             (elt (find (list i j k) gbp :key (lambda (x) (subseq x 0 3)) :test #'equal) 3)))
      (iter (with n = (reduce #'max gbp :key #'first))
            (with d = (reduce #'max gbp :key #'second))
            (for i from 0 to n)
            (collect (iter (for k from 0 to 1)
                           (collect (iter (for j from 0 to d)
                                          (collect (c i j k))))))))))

(defun mirror-concave-corner (ribbons i &optional (alpha 1) (beta 1))
  "Creates a twist-compatible, coherent concave corner at sides I-1 and I.
Destructively modifies RIBBONS, and returns the two newly created bilinear patch,
the first one is to use with ribbon I-1, the second is with ribbon I.
ALPHA is used in the U direction of ribbon I-1, BETA in the -U direction of ribbon I."
  (flet ((mirror (p o s) (v+ o (v* (v- o p) s)))
         (p (i j k) (elt (elt (elt ribbons i) k) j))
         (setp (i j k p) (setf (elt (elt (elt ribbons i) k) j) p)))
    (let* ((n (length ribbons))
           (i-1 (mod (1- i) n))
           (corner (p i 0 0))
           (twist (p i 1 1)))
      (setp i-1 2 1 (mirror twist (p i-1 2 0) beta))
      (setp  i  1 1 (mirror twist (p i 1 0) alpha))
      (setp i-1 3 1 (mirror (p i 1 0) corner beta))
      (setp  i  0 1 (mirror (p i-1 2 0) corner alpha))
      (let ((inner-twist (mirror (p i 1 1) (p i 0 1) beta)))
        (list (list (list corner (p i 0 1))
                    (list (p i-1 3 1) inner-twist))
              (list (list (p i-1 3 1) corner)
                    (list inner-twist (p i 0 1))))))))

(defun extend-ribbons (ribbons i)
  "Destructively modifies RIBBONS."
  (flet ((p (i j k) (elt (elt (elt ribbons i) k) j))
         (setp (i j k p) (setf (elt (elt (elt ribbons i) k) j) p)))
    (let* ((n (length ribbons))
           (i-1 (mod (1- i) n)))
      (setp i-1 2 1 (v+ (p i-1 1 1) (v- (p i-1 2 0) (p i-1 1 0))))
      (setp  i  1 1 (v+ (p  i  2 1) (v- (p  i  1 0) (p  i  2 0))))
      (setp i-1 3 1 (v+ (p i-1 2 1) (v- (p i-1 3 0) (p i-1 2 0))))
      (setp  i  0 1 (v+ (p  i  1 1) (v- (p  i  0 0) (p  i  1 0)))))))

(defun write-bezier-ribbon-control-points (ribbons filename &key center)
  "Assumes that all ribbons are of the same degrees."
  (let ((n (1- (length (first (first ribbons)))))
        (m (1- (length (first ribbons)))))
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (format s "~{~{~{v~{ ~f~}~%~}~}~}" ribbons)
      (when center
        (format s "v~{ ~f~}~%" center))
      (iter (repeat (length ribbons))
            (for i first 1 then (+ i (* (1+ n) (1+ m))))
            (iter (for j from 0 below n)
                  (iter (for k from 0 below m)
                        (format s "f ~a ~a ~a ~a~%"
                                (+ i j (* k (1+ n)))
                                (+ i (1+ j) (* k (1+ n)))
                                (+ i (1+ j) (* (1+ k) (1+ n)))
                                (+ i j (* (1+ k) (1+ n)))))))
      (when center
        (format s "p ~a~%" (1+ (* (length ribbons) (1+ n) (1+ m))))))))

(defun domain-from-ribbons-angular-concave (ribbons)
  (flet ((angle (r1 r2)
           (let* ((d1 (vnormalize (bezier (first r1) 1 1)))
                  (d2 (vnormalize (bezier (first r2) 0 1)))
                  (result (acos (scalar-product d1 d2))))
             (if (< (scalar-product d1 (v- (elt (second r2) 0) (elt (first r2) 0)))
                    0)
                 result
                 (- result))))
         (flip (a) (if (< a 0) (- (+ pi a)) (- pi a))))
    (let* ((angles (mapcar #'angle ribbons (append (rest ribbons) (list (first ribbons)))))
	   (angle-multiplier (/ (* (- (length ribbons) 2 (* 2 (count 0 angles :test #'>))) pi)
                                (reduce #'+ (mapcar (lambda (x) (flip x)) angles))))
	   (normalized-angles (mapcar (lambda (x) (flip (* (flip x) angle-multiplier))) angles))
	   (lengths (mapcar #'bezier-arc-length (mapcar #'first ribbons)))
	   (length-sum (reduce #'+ lengths))
	   (vertices (iter (for prev first '(0 0) then next)
			   (for length in lengths)
			   (for dir first 0 then (+ dir angle))
			   (for angle in normalized-angles)
			   (for next = (v+ prev (v* (list (cos dir) (sin dir)) length)))
			   (collect next)))
	   (difference (v* (car (last vertices)) -1)))
      ;; (format t "~{~f~%~}~%" (mapcar (lambda (x) (round (/ (* (flip x) 180) pi))) angles))
      ;; (format t "~{~f~%~}~%" (mapcar (lambda (x) (round (/ (* (flip x) 180) pi))) normalized-angles))
      ;; (write-domain-ribbons (cons '(0 0) vertices) nil "/tmp/domain-open.ps")
      (append (iter (for length in (butlast lengths))
                    (for accumulated first length then (+ accumulated length))
                    (for vertex in vertices)
                    (collect (v+ vertex (v* difference (/ accumulated length-sum)))))
              '((0 0))))))

(defun generate-domain-ribbons (domain)
  (let ((n (length domain)))
    (flet ((p (i) (elt domain i)))
      (iter (for i from 0 below n)
            (for i-1 = (mod (1- i) n))
            (for i+1 = (mod (1+ i) n))
            (for i-2 = (mod (- i 2) n))
            (for prev = (p i-2))
            (for next = (p i+1))
            (when (concavep domain i-1)
              (setf prev (v+ (p i-1)
                             (v* (affine-combine (vnormalize (v- (p i-1) prev))
                                                 1/2
                                                 (vnormalize (v- (p i-1) (p i))))
                                 (point-distance (p i) (p i-1))))))
            (when (concavep domain i)
              (setf next (v+ (p i)
                             (v* (affine-combine (vnormalize (v- (p i) next))
                                                 1/2
                                                 (vnormalize (v- (p i) (p i-1))))
                                 (point-distance (p i-1) (p i))))))
            (collect
                (list (iter (for j from 0 to 3)
                            (collect (affine-combine (p i-1) (/ j 3) (p i))))
                      (iter (for j from 0 to 3)
                            (collect (affine-combine prev (/ j 3) next)))))))))

(defun needs-extension-p (line1 line2 p)
  "LINE1 and LINE2 are lines given with as a vector of two points (START END),
and LINE2 is in the right half-plane of LINE1,
and they look at the same direction, i.e., their absolute angle is <90 degrees.
If P is between the two lines, it returns NIL;
if it is to the left of LINE1, returns START, and if it is to the right of LINE2, returns END."
  (let ((u1 (v- (second line1) (first line1)))
        (u2 (v- (second line2) (first line2)))
        (v1 (v- p (first line1)))
        (v2 (v- p (first line2))))
    (cond ((>= (- (* (first u1) (second v1))
                  (* (second u1) (first v1)))
              0)
           'start)
          ((<= (- (* (first u2) (second v2))
                  (* (second u2) (first v2)))
              0)
           'end)
          (t nil))))

(defun ribbon-extension-test (domain-ribbon p)
  (needs-extension-p (mapcar #'first domain-ribbon)
                     (mapcar (lambda (row) (car (last row))) domain-ribbon)
                     p))

(defun write-domain-ribbons (domain ribbons filename)
  (let* ((n (length domain))
         (points-and-ribbons
          (scale-to-unit
           (append domain (reduce #'append (reduce #'append ribbons)))))
         (points (subseq points-and-ribbons 0 n))
         (ribbons (subseq points-and-ribbons n))
         (lines (lines-from-points points)))
    (flet ((map-point (p)
             (list (+ (* (+ (first p) 1.0d0) 250) 50)
                   (+ (* (+ (second p) 1.0d0) 250) 50))))
      (with-open-file (s filename :direction :output :if-exists :supersede)
        (format s "%!PS-Adobe-2.0~%")
        (format s "%%BoundingBox: 0 0 600 600~%")
        (iter (for i from 0 below n)
              (for line in lines)
              (format s "% Segment: ~a~%" i)
              (format s "2 setlinewidth~%~
                         newpath~%~
                         ~{~f ~}moveto~%~
                         ~{~f ~}lineto~%~
                         stroke~%~
                         1 setlinewidth~%"
                      (map-point (first line))
                      (map-point (second line))))
        (iter (for point in ribbons)
              (format s "~{~f ~}3 0 360 arc fill~%" (map-point point)))
        (format s "showpage~%")))))

(defun n-sided-ribbon-eval-fn (domain-ribbon ribbon)
  "Returns a function that evaluates the given ribbon for a domain point.
DOMAIN-RIBBON and RIBBON are given as ((P00 P10 P20 P30) (P01 P11 P21 P31))."
  (lambda (p)
    (let* ((extension (ribbon-extension-test domain-ribbon p))
           (h (- (point-line-distance p (first domain-ribbon) t)))
           u surf3d surf2d)
      (macrolet ((s2d (i j) `(elt (elt surf2d ,j) ,i)))
        (ecase extension
          (start (setf surf3d (list (subseq (elt ribbon 0) 0 2)
                                    (subseq (elt ribbon 1) 0 2))
                       surf2d (list (subseq (elt domain-ribbon 0) 0 2)
                                    (subseq (elt domain-ribbon 1) 0 2))
                       (s2d 1 1) (v+ (s2d 0 0) ; paralelogramize
                                     (v- (s2d 1 0) (s2d 0 0))
                                     (v- (s2d 0 1) (s2d 0 0)))
                       u (first (reverse-bilinear-parameters surf2d p))))
          (end (setf surf3d (list (subseq (elt ribbon 0) 2)
                                  (subseq (elt ribbon 1) 2))
                     surf2d (list (subseq (elt domain-ribbon 0) 2)
                                  (subseq (elt domain-ribbon 1) 2))
                     (s2d 0 1) (v+ (s2d 1 0) ; paralelogramize
                                   (v- (s2d 1 1) (s2d 1 0))
                                   (v- (s2d 0 0) (s2d 1 0)))
                     u (first (reverse-bilinear-parameters surf2d p))))
          ((nil) (setf surf3d ribbon
                       surf2d domain-ribbon
                       u (first (reverse-cubic-parameters domain-ribbon p))))))
      (let* ((q1 (bezier (first surf3d) u))
             (q2 (bezier (second surf3d) u))
             (cross (vnormalize (v- q2 q1)))
             (len3d (bezier-arc-length (first ribbon)))
             (len2d (bezier-arc-length (first domain-ribbon))))
        (v+ q1 (v* cross (/ h len2d) len3d))))))

(defun concave-patch-test (ribbons points-file ribbons-file patch-file)
  (let* ((domain (domain-from-ribbons-angular-concave ribbons))
         (domain-ribbons (generate-domain-ribbons domain)))
    (destructuring-bind (vertices triangles)
        (shewchuk-triangle:mesh domain *resolution*)
      (write-domain-ribbons domain domain-ribbons "/tmp/domain.ps")
      (write-bezier-ribbon-control-points ribbons (format nil "~a.obj" points-file))
      (dotimes (i (length ribbons))
        (format t "Writing ribbon ~a...~%" i)
        (let ((eval-fn (n-sided-ribbon-eval-fn (elt domain-ribbons i) (elt ribbons i))))
          (write-obj-indexed-mesh (map 'vector eval-fn vertices) triangles
                                  (format nil "~a-~d.obj" ribbons-file i))))
      (format t "Writing the combined patch...~%")
      (harmonic:with-harmonic-coordinates (harmonic-map domain)
        (flet ((eval-patch (p)
                 (let ((result '(0 0 0)))
                   (dotimes (i (length ribbons))
                     (let* ((eval-fn (n-sided-ribbon-eval-fn (elt domain-ribbons i) (elt ribbons i)))
                            (point (funcall eval-fn p))
                            (weight (funcall (harmonic-kato harmonic-map domain i) p)))
                       (setf result (v+ result (v* point weight)))))
                   result)))
          (write-obj-indexed-mesh (map 'vector #'eval-patch vertices) triangles
                                  (format nil "~a.obj" patch-file)))))))

(defvar *dropbox* "/home/salvi/Dropbox")

#+nil
(let* ((*resolution* 40)
       ;; (gbp (format nil "~a~a" *dropbox* "/Shares/GrafGeo/Polar/bezier-ribbon/GBConvex1.gbp"))
       ;; (gbp (format nil "~a~a" *dropbox* "/Shares/GrafGeo/Polar/bezier-ribbon/GBTest4_Cubic.gbp"))
       (gbp (format nil "~a~a" *dropbox* "/Shares/GrafGeo/Polar/bezier-ribbon/GBUTest2_Cubic.gbp"))
       (ribbons (load-ribbons gbp))
       )
  (mirror-concave-corner ribbons 2)
  (mirror-concave-corner ribbons 3)
  (concave-patch-test ribbons "/tmp/pontok" "/tmp/ribbon" "/tmp/felulet"))

(defun bezier-surface-to-bspline (surface)
  (let ((n (1- (length surface)))
        (m (1- (length (first surface)))))
    (make-bspline-surface (list n m)
                          (list (append (make-list (1+ n) :initial-element 0.0)
                                        (make-list (1+ n) :initial-element 1.0))
                                (append (make-list (1+ m) :initial-element 0.0)
                                        (make-list (1+ m) :initial-element 1.0)))
                          surface)))

(defun sample-bezier-surface (surface resolution from to)
  "FROM and TO are 2D points. Returns a list of triangles."
  (iter (with length = (v- to from))
        (for i from 0 to resolution)
        (for u = (+ (first from) (* (/ i resolution) (first length))))
        (for last-row first nil then row)
        (for row = (iter (for j from 0 to resolution)
                         (for v = (+ (second from) (* (/ j resolution) (second length))))
                         (collect (bezier-surface surface (list u v)))))
        (when last-row
          (appending (iter (for j from 0 below resolution)
                           (collect (list (elt last-row j)
                                          (elt row j)
                                          (elt row (1+ j))))
                           (collect (list (elt row (1+ j))
                                          (elt last-row (1+ j))
                                          (elt last-row j))))))))


;;; Another idea: connect an outer and an inner domain

#+nil
(let ((*resolution* 50)
      (domain '((0 0) (3 0) (3 3) (2 3) (2 2) (1 2) (1 3) (0 3)))
      (outer-cp '(((0 0 0) (3 0 0))
                  ((0 3 0) (3 3 0))))
      (inner-cp '(((1 2 1) (2 2 1))
                  ((1 2.5 0) (2 2.5 0))
                  ((1 3 0) (2 3 0)))))
  (labels ((outer-domain (p) (v* p 1/3))
           (inner-domain (p) (v- p '(1 2)))
           (foo (x alpha) (/ (* alpha x x) (+ (* alpha x x) (* (- 1 x) (- 1 x)))))
           (eval-patch (p)
             (let ((outer-p (bezier-surface outer-cp (outer-domain p)))
                   (inner-p (bezier-surface inner-cp (inner-domain p)))
                   (mean (mean-value-coordinates domain p)))
               (affine-combine outer-p
                               (foo (reduce #'+ (mapcar #'* mean '(0 0 0 1 1 1 1 0))) 1/2)
                               inner-p))))
    (write-stl (eval3d-on-concave-domain domain #'eval-patch)
               "/tmp/proba.stl" :ascii t)))

#+nil
(let ((*resolution* 50)
      (domain '((0 0) (3 0) (3 3) (1 3) (1 2) (0 2)))
      (outer-cp '(((0 0 0) (3 0 0))
                  ((0 3 0) (3 3 0))))
      (inner-cp '(((0 2 0) (0.5 2 0) (1 2 0.5))
                  ((0 2.5 0) (0.5 2.5 0) (1 2.5 0))
                  ((0 3 0) (0.5 3 0) (1 3 0)))))
  (labels ((outer-domain (p) (v* p 1/3))
           (inner-domain (p) (v- p '(0 2)))
           (foo (x alpha) (/ (* alpha x x) (+ (* alpha x x) (* (- 1 x) (- 1 x)))))
           (eval-patch (p)
             (let ((outer-p (bezier-surface outer-cp (outer-domain p)))
                   (inner-p (bezier-surface inner-cp (inner-domain p)))
                   (mean (mean-value-coordinates domain p)))
               (affine-combine outer-p
                               (foo (reduce #'+ (mapcar #'* mean '(0 0 0 1 1 1 1 0))) 1/2)
                               inner-p))))
    (write-stl (eval3d-on-concave-domain domain #'eval-patch)
               "/tmp/proba.stl" :ascii t)))

(defun sliced-concave-distance-function-test (points fn file-or-stream
                                              &key (resolution 0.001) (density 0.1) elements)
  "FN gives a value between 0 and 1 for a given point.
When FILE-OR-STREAM is a stream, no header / showpage is written."
  (let* ((n (length points)) 
         (lines (lines-from-points points)))
    (labels ((map-point (p)
               (list (+ (* (+ (first p) 1.0d0) 250) 50)
                     (+ (* (+ (second p) 1.0d0) 250) 50)))
             (writer (s)
               (iter (for i from 0 below n)
                     (for line in lines)
                     (format s "% Segment: ~a~%" i)
                     (format s "2 setlinewidth~%~
                                newpath~%~
                                ~{~f ~}moveto~%~
                                ~{~f ~}lineto~%~
                                stroke~%~
                                1 setlinewidth~%"
                             (map-point (first line))
                             (map-point (second line))))
               (iter (for element in elements)
                     (if (atom (first element))
                         (format s "% Point~%~
                                    newpath~%~
                                    ~{~f ~}3 0 360 arc closepath~%~
                                    fill~%"
                                 (map-point element))
                         (format s "% Line~%~
                                    2 setlinewidth~%~
                                    newpath~%~
                                    ~{~f ~}moveto~%~
                                    ~{~f ~}lineto~%~
                                    stroke~%~
                                    1 setlinewidth~%"
                                 (map-point (first element))
                                 (map-point (second element)))))
               (destructuring-bind (vertices triangles)
                   (shewchuk-triangle:mesh points resolution)
                 (let* ((vertices (eval-over-domain vertices fn))
                        (segments (slice-mesh vertices triangles density)))
                   (iter (for segment in segments)
                         (format s "newpath~%~
                                     ~{~f ~}moveto~%~
                                     ~{~f ~}lineto~%~
                                     stroke~%"
                                 (map-point (cdr (first segment)))
                                 (map-point (cdr (second segment)))))))))
      (if (streamp file-or-stream)
          (writer file-or-stream)
          (with-open-file (s file-or-stream :direction :output :if-exists :supersede)
            (format s "%!PS-Adobe-2.0~%")
            (format s "%%BoundingBox: 0 0 600 600~%")
            (writer s)
            (format s "showpage~%"))))))

(defun scale-to-unit (points)
  "Scale POINTS that it fits in [-1,1]x[-1x1]."
  (let* ((min (list (reduce #'min points :key #'first)
                    (reduce #'min points :key #'second)))
         (max (list (reduce #'max points :key #'first)
                    (reduce #'max points :key #'second)))
         (len (reduce #'max (v- max min))))
    (mapcar (lambda (p)
              (v- (v* (v- p min) (/ 2 len)) '(1 1)))
            points)))

#+nil
(let ((*resolution* 30)
      #+nil(points '((1 2) (1 1) (2 1) (2 0) (0 0) (0 2)))
      #+nil(points '((0 0) (0 7) (1 7) (1 4) (2 4) (2 7) (3 7) (3 0) (2 0) (2 3) (1 3) (1 0))) ; H
      (points '((0 0) (6 0) (6 6) (4 6) (4 4) (2 4) (2 6) (0 6)))) ; U
  (let ((points (scale-to-unit points)))
    (harmonic:with-harmonic-coordinates (h points)
      (flet ((fn (p)
               (let ((l (harmonic:harmonic-coordinates h p)))
                 (when (member nil l)   ; kutykurutty
                   (let ((*barycentric-type* 'meanvalue))
                     (setf l (barycentric-coordinates points p))))
                 #+nil(barycentric-d l 0)
                 (- 1 (elt l 3) (elt l 4) (elt l 5) (elt l 6)) ;; composite
                 #+nil(elt l 4))))
        (sliced-concave-distance-function-test points #'fn "/tmp/proba.ps"
                                               :resolution 0.001 :density 0.1d0)))))


;;; Grid-parameterized concave patches

(defun sigmoid-gamma (x)
  "Variations (1/5, 1/4, 1/3, 1/2, 1, 2, infty):
  x / (1 + 2 * |x|)              => 0.143, 0.167, 0.200, 0.250, 0.333, 0.400, 0.5 (w/discontinuity)
  tanh(x) / 2                    => 0.099, 0.122, 0.161, 0.231, 0.381, 0.482, 0.5 %
  x / 2 sqrt(1 + x^2)            => 0.098, 0.121, 0.158, 0.224, 0.354, 0.447, 0.5 %
  1 / (1 + e^(-4x)) - 1/2        => 0.190, 0.231, 0.291, 0.381, 0.482, 0.499, 0.5
  x / sqrt(1 + 9x^2)             => 0.171, 0.200, 0.236, 0.277, 0.316, 0.329, 0.333
  tanh(x/2) = 2 / (1 + e^-x) - 1 => 0.098, 0.124, 0.165, 0.245, 0.462, 0.762, 1 %
  x / sqrt(4 + x^2)              => 0.100, 0.124, 0.164, 0.243, 0.447, 0.707, 1 %
Those marked with `%' have a derivative other than 1 at 0."
  (if *use-gamma*
      #+nil(/ (tanh x) 2)
      #+nil(/ x 2 (sqrt (+ 1 (* x x))))
      #+nil(- (/ 1 (1+ (exp (* -4 x)))) 1/2)
      (/ x (sqrt (1+ (* 9 x x))))
      #+nil(tanh (/ x 2))
      #+nil(/ x (sqrt (+ 4 (* x x))))
      x))

(defvar *extension-degree*)
(defvar *extension-shrinking*)

(defun n-sided-ribbon-grid-eval-fn (domain-edge ribbon)
  "Returns a function that evaluates the given ribbon for a domain point.
RIBBON is given as ((P00 P10 P20 P30) (P01 P11 P21 P31))."
  (flet ((blend (x)
           (if (>= x 1)
               (1+ *extension-shrinking*)
               (+ 1 (* *extension-shrinking* (hermite-blend-function 'point 'end x))))))
    (lambda (p)
      (let* ((dir (v- (second domain-edge) (first domain-edge)))
             (s (/ (scalar-product (v- p (first domain-edge)) dir) (vlength2 dir)))
             (extension (cond ((< s 0) 'start) ((> s 1) 'end) (t nil)))
             (h (- (point-line-distance p domain-edge t))) ; signed distance from base segment (2D)
             (len2d (vlength dir))                         ; length of base segment (2D)
             (len3d (bezier-arc-length (first ribbon)))    ; length of base curve (3D)
             (surf3d (case extension
                       (start (list (subseq (first ribbon) 0 (1+ *extension-degree*))
                                    (subseq (second ribbon) 0 (1+ *extension-degree*))))
                       (end (list (subseq (first ribbon) (- 3 *extension-degree*))
                                  (subseq (second ribbon) (- 3 *extension-degree*))))
                       ((nil) ribbon)))
             (u (case extension
                  (start (* s (/ 3 *extension-degree* (blend (- s)))))
                  (end (1+ (* (1- s) (/ 3 *extension-degree* (blend (1- s))))))
                  ((nil) s)))
             (len (* (point-distance (bezier domain-edge s) p) ; signed length on sweepline (2D)
                     (if (< h 0) -1 1)))
             (q1 (bezier (first surf3d) u))
             (q2 (bezier (second surf3d) u))
             (sweep (vnormalize (v- q2 q1))))              ; normalized sweep direction (3D)
        (v+ q1 (v* sweep (sigmoid-gamma (/ len len2d)) len3d *ribbon-multiplier*))))))

(defun ribbon-force-perpendicular (ribbon)
  "Destructively changes the second control point row such that the cross-derivatives
are perpendicular to the tangent."
  (iter (with n = (1- (length (first ribbon))))
        (for i from 0 to n)
        (for u = (/ i n))
        (for d = (vnormalize (bezier (first ribbon) u 1)))
        (for p = (elt (first ribbon) i))
        (for q = (elt (second ribbon) i))
        (setf (elt (second ribbon) i)
              ;; Project Q into the plane defined by P and D
              (v+ q (v* d (scalar-product (v- p q) d))))))

(defun ribbon-uniform-length (ribbon)
  "Destructively changes the second control point row such that the cross-derivatives
are of equal length (arc length of the base divided by the degree)."
  (iter (with n = (1- (length (first ribbon))))
        (with len = (/ (bezier-arc-length (first ribbon)) n))
        (for i from 0 to n)
        (for p = (elt (first ribbon) i))
        (for q = (elt (second ribbon) i))
        (setf (elt (second ribbon) i)
              (v+ p (v* (vnormalize (v- q p)) len)))))

(defun harmonic-hermite (map points i)
  (lambda (p)
    (let ((i-1 (mod (1- i) (length points)))
          (l (harmonic:harmonic-coordinates map p)))
      (when (member nil l)              ; kutykurutty
        (setf l (mean-value-coordinates points p)))
      (/ (hermite-blend-function 'point 'start (- 1 (elt l i-1) (elt l i))) 2))))

(defun concave-grid-patch-test (ribbons points-file patch-file &optional ribbons-file)
  (let ((n (length ribbons))
        (domain (domain-from-ribbons-angular-concave ribbons)))
    (destructuring-bind (vertices triangles)
        (shewchuk-triangle:mesh domain *resolution*)
      (write-domain-ribbons domain '() "/tmp/domain.ps")
      (mapc #'ribbon-force-perpendicular ribbons)
      (mapc #'ribbon-uniform-length ribbons)
      (write-bezier-ribbon-control-points ribbons (format nil "~a.obj" points-file))
      (when ribbons-file
        (dotimes (i n)
          (format t "Writing ribbon ~a...~%" i)
          (let ((eval-fn (n-sided-ribbon-grid-eval-fn
                          (list (elt domain (mod (1- i) n)) (elt domain i)) (elt ribbons i))))
            (write-obj-indexed-mesh (map 'vector eval-fn vertices) triangles
                                    (format nil "~a-~d.obj" ribbons-file i)))))
      (format t "Writing the combined patch...~%")
      (harmonic:with-harmonic-coordinates (harmonic-map domain)
        (flet ((eval-patch (p)
                 (let ((result '(0 0 0)))
                   (dotimes (i n)
                     (let* ((eval-fn (n-sided-ribbon-grid-eval-fn
                                      (list (elt domain (mod (1- i) n)) (elt domain i))
                                      (elt ribbons i)))
                            (point (funcall eval-fn p))
                            (weight (funcall (harmonic-kato harmonic-map domain i) p)))
                       (setf result (v+ result (v* point weight)))))
                   result)))
          (write-obj-indexed-mesh (map 'vector #'eval-patch vertices) triangles
                                  (format nil "~a.obj" patch-file)))))))

#+nil
(let* ((*resolution* 2)
       (*extension-degree* 2)
       (*extension-shrinking* 0)
       (*ribbon-multiplier* 1)
       (*use-gamma* t)
       (tests "/Shares/GrafGeo/Polar/bezier-ribbon/")
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBConvex1.gbp"))         ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBTest4_Cubic.gbp"))     ; 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBUTest2_Cubic.gbp"))    ; 2 3
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "6sided.gbp"))            ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Plane.gbp")) ; 4 5
       (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Cylinder.gbp")) ; 4 5
       (ribbons (load-ribbons gbp)))
  (mirror-concave-corner ribbons 4)
  (mirror-concave-corner ribbons 5)
  (concave-grid-patch-test ribbons "/tmp/pontok" "/tmp/felulet" #+nil"/tmp/ribbon"))


;;; Concave Gregory patches

(defun parallel-parameter (domain i side p)
  "A side-parameter based on side I, parallel with one of the neighboring sides.
SIDE can be 'LEFT or 'RIGHT (referring to the previous and next sides, respectively)."
  (let* ((n (length domain))
         (p1 (elt domain i))
         (p2 (elt domain (mod (1+ i) n)))
         (dir (ecase side
                (left (v- p1 (elt domain (mod (1- i) n))))
                (right (v- p2 (elt domain (mod (+ i 2) n))))))
         (q (line-line-intersection (list p1 p2) (list p (v+ p dir)))))
    (* (/ (point-distance q p1) (point-distance p2 p1))
       (if (< (scalar-product (v- q p1) (v- p2 p1)) 0) -1 1))))

(defun concave-side-ribbon (ribbon s d)
  "RIBBON is given as ((P00 P10 P20 P30) (P01 P11 P21 P31))."
  (flet ((blend (x)
           (if (>= x 1)
               (1+ *extension-shrinking*)
               (+ 1 (* *extension-shrinking* (hermite-blend-function 'point 'end x))))))
    (let* ((extension (cond ((< s 0) 'start) ((> s 1) 'end) (t nil)))
           (surf3d (case extension
                     (start (list (subseq (first ribbon) 0 (1+ *extension-degree*))
                                  (subseq (second ribbon) 0 (1+ *extension-degree*))))
                     (end (list (subseq (first ribbon) (- 3 *extension-degree*))
                                (subseq (second ribbon) (- 3 *extension-degree*))))
                     ((nil) ribbon)))
           (u (case extension
                (start (* s (/ 3 *extension-degree* (blend (- s)))))
                (end (1+ (* (1- s) (/ 3 *extension-degree* (blend (1- s))))))
                ((nil) s)))
           (q1 (bezier (first surf3d) u))
           (q2 (bezier (second surf3d) u))
           (sweep (v* (v- q2 q1) 3)))
      (v+ q1 (v* sweep (sigmoid-gamma d) *ribbon-multiplier*)))))

(defun concave-correction-patch (ribbons i si-1 si)
  (let* ((n (length ribbons))
         (i-1 (mod (1- i) n))
         (corner (elt (first (elt ribbons i)) 0))
         (prev (elt (first (elt ribbons i-1)) 2))
         (next (elt (first (elt ribbons i)) 1))
         (di-1 (v* (v- prev corner) 3))
         (di (v* (v- next corner) 3))
         (twist-cp (elt (second (elt ribbons i-1)) 2) ; assume twist compatibility
           #+nil(v* (v+ (v* (elt (second (elt ribbons i-1)) 2) si)
                        (v* (elt (second (elt ribbons i)) 1) si-1))
                    (handler-case (/ (+ si-1 si))
                      (division-by-zero () 0))))
         (twist (v* (v- (v+ corner twist-cp) (v+ prev next)) 9)))
    (v+ corner (v* di-1 si-1 ) (v* di si) (v* twist si-1 si))))

(defun concave-corner-ribbon-parallel (domain ribbons i p)
  "RIBBONS is unfortunately indexed in a different system, where edge I is bounded by
vertices I-1 and I, but here we use vertices I and I+1..."
  (let* ((n (length domain))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (s1 (- 1 (parallel-parameter domain i-1 'right p)))
         (s2 (parallel-parameter domain i 'left p)))
    #+nil(concave-side-ribbon (elt ribbons i) (- 1 s1) s2) ; side ribbon 1
    #+nil(concave-side-ribbon (elt ribbons i+1) s2 s1)     ; side ribbon 2
    #+nil(concave-correction-patch ribbons i+1             ; correction patch
                                   (* (sigmoid-gamma s1) *ribbon-multiplier*)
                                   (* (sigmoid-gamma s2) *ribbon-multiplier*))
    (v- (v+ (concave-side-ribbon (elt ribbons i) (- 1 s1) s2)
            (concave-side-ribbon (elt ribbons i+1) s2 s1))
        (concave-correction-patch ribbons i+1
                                  (* (sigmoid-gamma s1) *ribbon-multiplier*)
                                  (* (sigmoid-gamma s2) *ribbon-multiplier*)))))

(defun concave-corner-ribbon (domain ribbons i p)
  (let* ((n (length domain))
         (i-1 (mod (1- i) n))
         (i+1 (mod (1+ i) n))
         (edge-i (list (elt domain i-1) (elt domain i)))
         (edge-i+1 (list (elt domain i) (elt domain i+1)))
         (ri (n-sided-ribbon-grid-eval-fn edge-i (elt ribbons i)))
         (ri+1 (n-sided-ribbon-grid-eval-fn edge-i+1 (elt ribbons i+1)))
         (di (point-line-distance p edge-i))
         (di+1 (point-line-distance p edge-i+1))
         (alpha (safe-/ (* di di) (+ (* di di) (* di+1 di+1)))))
    (v+ (v* (funcall ri p) (- 1 alpha))
        (v* (funcall ri+1 p) alpha))))

(defun harmonic-corner-blend (map points i)
  (lambda (p)
    (let ((l (harmonic:harmonic-coordinates map p)))
      (when (member nil l)              ; kutykurutty
        (setf l (mean-value-coordinates points p)))
      (let ((d (iter (with n = (length l))
                     (for i from 0 below n)
                     (for i-1 = (mod (1- i) n))
                     (collect (- 1 (elt l i-1) (elt l i))))))
        (corner-blend d i)))))

(defun ribbon-force-zero-twists (ribbon)
  "Destructively the twist vectors to form a paralelogram."
  (macrolet ((cp (p) `(elt (elt ribbon (second ,p)) (first ,p))))
    (flet ((set-zero-twist (p q1 q2 tw)
             (let ((corner (cp p)))
               (setf (cp tw) (v+ corner (v- (cp q1) corner) (v- (cp q2) corner))))))
      (let ((n (1- (length (first ribbon)))))
        (set-zero-twist '(0 0) '(1 0) '(0 1) '(1 1))
        (set-zero-twist `(,n 0) `(,(1- n) 0) `(,n 1) `(,(1- n) 1))))))

(defun concave-gregory-test (ribbons points-file patch-file &optional ribbons-file)
  (let ((n (length ribbons))
        (domain (domain-from-ribbons-angular-concave ribbons)))
    (destructuring-bind (vertices triangles)
        (shewchuk-triangle:mesh domain *resolution*)
      (write-domain-ribbons domain '() "/tmp/domain.ps")
      (mapc #'ribbon-force-zero-twists ribbons)
      ;; (mapc #'ribbon-force-perpendicular ribbons)
      ;; (mapc #'ribbon-uniform-length ribbons)
      (write-bezier-ribbon-control-points ribbons (format nil "~a.obj" points-file))
      (when ribbons-file
        (let ((*use-gamma* nil))
          (dotimes (i n)
            (format t "Writing ribbon ~a...~%" i)
            (flet ((eval-fn (p) (concave-corner-ribbon-parallel domain ribbons i p)))
              (write-obj-indexed-mesh (map 'vector #'eval-fn vertices) triangles
                                      (format nil "~a-~d.obj" ribbons-file i))))))
      (format t "Writing the combined patch...~%")
      (harmonic:with-harmonic-coordinates (harmonic-map domain)
        (flet ((eval-patch (p)
                 (let ((result '(0 0 0)))
                   (dotimes (i n)
                     (flet ((eval-fn (p) (concave-corner-ribbon-parallel domain ribbons i p)))
                       (let ((point (eval-fn p))
                             (weight (funcall (harmonic-corner-blend harmonic-map domain i) p)))
                         (setf result (v+ result (v* point weight))))))
                   result)))
          (write-obj-indexed-mesh (map 'vector #'eval-patch vertices) triangles
                                  (format nil "~a.obj" patch-file)))))))

#+nil
(let* ((*resolution* 2)
       (*extension-degree* 2)
       (*extension-shrinking* 0)
       (*ribbon-multiplier* 1)
       (*use-gamma* t)
       (tests "/Shares/GrafGeo/Polar/bezier-ribbon/")
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBConvex1.gbp"))            ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBTest4_Cubic.gbp"))        ; 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBUTest2_Cubic.gbp"))       ; 2 3
       (gbp (format nil "~a~a~a" *dropbox* tests "6sided.gbp"))               ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Plane.gbp"))    ; 4 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Cylinder.gbp")) ; 4 5
       (ribbons (load-ribbons gbp)))
  ;; (mirror-concave-corner ribbons 4)
  ;; (mirror-concave-corner ribbons 5)
  (concave-gregory-test ribbons "/tmp/pontok" "/tmp/felulet" "/tmp/ribbon"))


;;; Concave Generalized Bezier patch

(defun harmonic-coordinates (map points p)
  (let ((l (harmonic:harmonic-coordinates map p)))
    (when (member nil l)                ; kutykurutty
      (setf l (mean-value-coordinates points p)))
    l))

(defun concave-generalized-bernstein (map points p side degree col row &key (use-d t))
  (flet ((sqr (x) (* x x)))
    (let* ((n (length points))
           (l (harmonic-coordinates map points p))
           (i side)
           (i-1 (mod (1- i) n))
           (i+1 (mod (1+ i) n))
           (si (barycentric-s l i))
           (si-1 (barycentric-s l i-1))
           (si+1 (barycentric-s l i+1))
           (di (barycentric-d l i))
           (di-1 (barycentric-d l i-1))
           (di+1 (barycentric-d l i+1))
           (alpha (if use-d
                      (if (< (+ (sqr di-1) (sqr di)) *epsilon*)
                          0.5
                          (/ (sqr di-1) (+ (sqr di-1) (sqr di))))
                      (if (< (+ (sqr si) (sqr (- 1 si-1))) *epsilon*)
                          0.5
                          (/ (sqr si) (+ (sqr si) (sqr (- 1 si-1)))))))
           (beta (if use-d
                     (if (< (+ (sqr di+1) (sqr di)) *epsilon*)
                         0.5
                         (/ (sqr di+1) (+ (sqr di+1) (sqr di))))
                     (if (< (+ (sqr (- 1 si)) (sqr si+1)) *epsilon*)
                         0.5
                         (/ (sqr (- 1 si)) (+ (sqr (- 1 si)) (sqr si+1))))))
           (blend (* (bernstein degree row di)
                     (bernstein degree col si)))
           (mu (cond ((and (< row 2) (< col 2)) alpha)
                     ((and (< row 2) (> col (- degree 2))) beta)
                     ((or (< col row) (> col (- degree row))) 0)
                     ((or (= col row) (= col (- degree row))) 1/2)
                     (t 1))))
      (* blend mu))))

(defun concave-bezier-deficiency (map points p degree &key (use-d t))
  (- 1 (iter (for side from 0 below (length points))
             (sum (iter (for row from 0 below (ceiling degree 2))
                        (sum (iter (for col from 0 to degree)
                                   (sum (concave-generalized-bernstein
                                         map points p side degree col row :use-d use-d)))))))))

(defun write-concave-bernstein-blend-mesh (map points path degree concave
                                           &key (use-d t) (scaling 1))
  "CONCAVE is a list of concave vertex indices.
Blending functions near these vertices are computed separately."
  (destructuring-bind (vertices triangles)
      (shewchuk-triangle:mesh points *resolution*)
    (iter (with n = (length points))
          (with half = (ceiling degree 2))
          (for side from 0 below n)
          (for side-1 = (mod (1- side) n))
          (for side+1 = (mod (1+ side) n))
          (iter (for row from 0 below half)
                (iter (for col from 0 to degree)
                      (for fname = (format nil "~a-deg~a-c~a~a~a.obj" path degree side col row))
                      (flet ((foo (p)
                               (let ((b (concave-generalized-bernstein
                                         map points p side degree col row :use-d use-d))
                                     (bp (concave-generalized-bernstein
                                          map points p side-1 degree (- degree row) col
                                          :use-d use-d))
                                     (bn (concave-generalized-bernstein
                                          map points p side+1 degree row (- degree col)
                                          :use-d use-d)))
                                 (if (or (and (< col half) (member side concave))
                                         (and (>= col half) (member side+1 concave)))
                                     (* b scaling)
                                     (* (+ b bp bn) scaling)))))
                        (write-obj-indexed-mesh (eval-over-domain vertices #'foo)
                                                triangles fname)))))
    (let ((fname (format nil "~a-deg~a-center.obj" path degree)))
      (with-open-file (s fname :direction :output :if-exists :supersede)
        (flet ((foo (p)
                 (let ((d (concave-bezier-deficiency map points p degree :use-d use-d)))
                   (* (if (< (abs d) *epsilon*) 0 d) scaling))))
          (write-obj-indexed-mesh (eval-over-domain vertices #'foo) triangles fname))))))

(defun concave-generalized-bezier-eval (map points p degree ribbons center &key (use-d t))
  (declare (ignore center))
  (let ((n (length points))
        (result '(0 0 0))
        (wsum 0))
    (iter (for side from 0 below n)
          (iter (for row from 0 below (ceiling degree 2))
                (iter (for col from 0 to degree)
                      (for weight = (concave-generalized-bernstein
                                     map points p side degree col row :use-d use-d))
                      (setf wsum (+ wsum weight)
                            result (v+ result
                                       (v* (elt (elt (elt ribbons side) row) col) weight))))))
    (v* result (/ wsum))
    #+nil(v+ result (v* center (- 1 wsum)))))

(defun unify-concave-corner (ribbons i)
  "Destructively modifies RIBBONS."
  (let* ((n (length ribbons))
         (i-1 (mod (1- i) n)))
    (macrolet ((p (i j k) `(elt (elt (elt ribbons ,i) ,k) ,j)))
      (let* ((q1 (v+ (p i-1 3 1) (v- (p i-1 3 1) (p i-1 2 1))))
             (q2 (v+ (p i 0 1) (v- (p i 0 1) (p i 1 1))))
             (q (affine-combine q1 1/2 q2)))
        (setf (p i-1 3 1) q
              (p i 0 1) q)))))

(defun generalized-bezier-generate-center (ribbons)
  (let ((result '(0 0 0)))
    (dolist (r ribbons)
      (setf result (v+ result (elt (second r) 1) (elt (second r) 2))))
    (v* result (/ (* 2 (length ribbons))))))

(defun write-ribbon-surfaces (ribbons path)
  (iter (for r in ribbons)
        (for i upfrom 0)
        (for filename = (format nil "~a/ribbon-~a.stl" path i))
        (write-stl (sample-bezier-surface r *resolution* '(0 0) '(1 1)) filename)))

(defun concave-generalized-bezier-test (input-file concave output-file)
  "CONCAVE is a list of concave vertex indices."
  (let ((ribbons (load-ribbons input-file)))
    (dolist (i concave)
      (mirror-concave-corner ribbons i))
    (let ((domain (domain-from-ribbons-angular-concave ribbons))
          (center (generalized-bezier-generate-center ribbons)))
      (write-bezier-ribbon-control-points ribbons "/tmp/ribbons.obj" :center center)
      (export-cgb-ribbons ribbons center "/tmp/ribbons.cgb")
      ;; (let ((*resolution* 30))
      ;;   (write-ribbon-surfaces ribbons "/tmp"))
      (write-domain-ribbons domain '() "/tmp/domain.ps")
      (harmonic:with-harmonic-coordinates (hmap domain :levels 9)
        (destructuring-bind (vertices triangles)
            (shewchuk-triangle:mesh domain *resolution*)
          (flet ((foo (p) (concave-generalized-bezier-eval hmap domain p 3 ribbons center)))
            (write-obj-indexed-mesh (map 'vector #'foo vertices) triangles output-file)))))))

(defun concave-gb-normal-test (input-file concave &key pos-tolerance angle-tolerance step)
  "CONCAVE is a list of concave vertex indices.
Assumes that matter is always on the left side of the edges in the domain."
  (let ((ribbons (load-ribbons input-file)))
    (dolist (i concave)
      (mirror-concave-corner ribbons i))
    (let ((domain (domain-from-ribbons-angular-concave ribbons))
          (center (generalized-bezier-generate-center ribbons)))
      (harmonic:with-harmonic-coordinates (hmap domain :levels 10)
        (iter (with n = (length domain))
              (for ribbon in ribbons)
              (for i from 0 below n)
              (for edge = (list (elt domain (mod (1- i) n)) (elt domain i)))
              (for dir = (vnormalize (v- (second edge) (first edge))))
              (for cross-dir = (list (- (second dir)) (first dir)))
              (iter (for j from 0 to *resolution*)
                    (for u = (/ j *resolution*))
                    (for tangent = (bezier (first ribbon) u 1))
                    (for cross = (v- (bezier (second ribbon) u)
                                     (bezier (first ribbon) u)))
                    (for normal = (vnormalize (cross-product tangent cross)))
                    (for uv = (affine-combine (first edge) u (second edge)))
                    (for uv2 = (v+ uv (v* cross-dir step)))
                    (for p = (concave-generalized-bezier-eval hmap domain uv 3 ribbons center))
                    (when (> (point-distance p (bezier (first ribbon) u)) pos-tolerance)
                      (format t "Positional error: ~5f,	side: ~a,	u: ~5f~%~a~%vs~%~a~%"
                              (point-distance p (bezier (first ribbon) u))
                              i u p (bezier (first ribbon) u)))
                    (for p2 = (concave-generalized-bezier-eval hmap domain uv2 3 ribbons center))
                    (for approx-normal = (vnormalize (cross-product tangent (v- p2 p))))
                    (for alpha = (* (acos (scalar-product normal approx-normal)) (/ 180 pi)))
                    (when (> alpha angle-tolerance)
                      (format t "Error: ~5f	side: ~a	u: ~5f~%" alpha i u))))))))

(defun export-cgb-ribbons (ribbons center filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "~a~%~{~f~^ ~}~%~{~{~a ~a~%~{~{~{~f~^ ~}~%~}~}~}~}"
            (length ribbons) center
            (mapcar (lambda (r)
                      (list (1- (length (first r))) (length r) r))
                    ribbons))))

#+nil
(defun harmonic-coordinates (map points p)
  "Mean value (!) coordinates - for testing."
  (declare (ignore map))
  (mean-value-coordinates points p))

;;; Use mean value coordinates for this one (see above)
#+nil
(let* ((tests "/Shares/GrafGeo/Polar/bezier-ribbon/")
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBConvex1.gbp"))            ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBTest4_Cubic.gbp"))        ; 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBUTest2_Cubic.gbp"))       ; 2 3
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "6sided.gbp"))               ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Plane.gbp"))    ; 4 5
       (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Cylinder.gbp")) ; 4 5
       (*resolution* 30))
  (concave-gb-normal-test gbp '(4 5)
                          :pos-tolerance 1.0d-8 :angle-tolerance 1.0d-4 :step 0.001))

#+nil
(let ((*resolution* 30)
      #+nil(points '((1 2) (1 1) (2 1) (2 0) (0 0) (0 2)))
      #+nil(points '((0 0) (0 7) (1 7) (1 4) (2 4) (2 7) (3 7) (3 0) (2 0) (2 3) (1 3) (1 0))) ; H
      (points '((0 0) (6 0) (6 6) (4 6) (4 4) (2 4) (2 6) (0 6)))) ; U
  (let ((points (scale-to-unit points))
        (*barycentric-d-function* #'barycentric-d-1minus)
        (*barycentric-dilation* 2))
    (harmonic:with-harmonic-coordinates (h points :levels 10)
      (flet ((fn (p) (barycentric-d (harmonic-coordinates h points p) 5)))
        (sliced-concave-distance-function-test points #'fn "/tmp/proba.ps"
                                               :resolution 0.0001 :density 0.1d0)))))

#+nil
(let* ((tests "/Shares/GrafGeo/Polar/bezier-ribbon/")
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBConvex1.gbp"))            ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBTest4_Cubic.gbp"))        ; 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBUTest2_Cubic.gbp"))       ; 2 3
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "6sided.gbp"))               ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Plane.gbp"))    ; 4 5
       (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Cylinder.gbp")) ; 4 5
       (ribbons (load-ribbons gbp)))
  (mirror-concave-corner ribbons 4)
  (mirror-concave-corner ribbons 5)
  (let ((*resolution* 1)
        (domain (domain-from-ribbons-angular-concave ribbons)))
    (harmonic:with-harmonic-coordinates (harmonic-map domain :levels 10)
      (write-concave-bernstein-blend-mesh harmonic-map domain "/tmp/U" 3 '(4 5) :scaling 100))))

#+nil
(let* ((tests "/Shares/GrafGeo/Polar/bezier-ribbon/")
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBConvex1.gbp"))            ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBTest4_Cubic.gbp"))        ; 5
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "GBUTest2_Cubic.gbp"))       ; 2 3
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "6sided.gbp"))               ; -
       ;; (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Plane.gbp"))    ; 4 5
       (gbp (format nil "~a~a~a" *dropbox* tests "ConcaveTest_Cylinder.gbp")) ; 4 5
       (*resolution* 2))
  (concave-generalized-bezier-test gbp '(4 5) "/tmp/felulet.obj"))
