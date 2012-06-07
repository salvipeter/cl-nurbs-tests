(in-package :cl-nurbs-tests)

(defparameter *window* '((-1 -1) (1 1)))
(defparameter *curve*
  ;; (make-bspline-curve 3 '(0 0 0 0 0.5 1 1 1 1)
  ;;                     '((-0.5 -0.2) (-0.2 0.3) (0.1 0.4) (0.4 0) (0.6 -0.2)))
  ;; (make-bspline-curve 3 '(0 0 0 0 0.3 0.7 1 1 1 1)
  ;;                     '((0 -0.5) (-0.5 -0.5) (-0.5 0.8)
  ;;                       (0.5 0.8) (0.5 -0.5) (0 -0.5)))
  (make-bspline-curve 3 '(0 0 0 0 0.5 1 1 1 1)
                      '((-0.8 -0.3) (-0.6 0.4) (-0.3 0.5) (0.7 0.2) (0.8 -0.1)))
  "A B-spline curve with the parameter running in [0,1].")
(defparameter *resolution* 41)
(defparameter *epsilon* 1.0d-8)

(defun closest-point (p)
  (bsc-evaluate *curve* (bsc-project-point *curve* p 10 100)))

(defun closest-point-extended (p)
  (closest-point (v+ p (v* (v- (closest-point p) p) 2))))

(defun grid-point (p)
  "P is an index, e.g. (1 3)."
  (v+ (first *window*)
      (mapcar '* (v* p (/ (1- *resolution*)))
              (v* (apply #'v- *window*) -1))))

(defun closest-grid-point-low (p)
  "Finds the closest grid point indices in the lower left quarter of P."
  (let ((coord (mapcar (lambda (v l) (* v (/ l)))
                       (v- p (first *window*))
                       (v* (apply #'v- *window*) -1))))
    (mapcar #'floor (v* coord (1- *resolution*)))))

(defun closest-grid-point (p)
  (let ((q (closest-grid-point-low p)))
    (iter (for d in '((0 0) (1 0) (0 1) (1 1)))
          (for point = (grid-point (v+ q d)))
          (finding (v+ q d) minimizing (point-distance p point)))))

(let ((w (iter (for i from 0 to 3)
               (collect (/ (iter (for k from 0 to 3)
                                 (when (/= i k)
                                   (multiply (/ (- i k) 3)))))))))
  (defun uniform-cubic-barycentric-lagrange-weights (x)
    (cond ((< (abs x) *epsilon*) '(1 0 0 0))
          ((< (abs (- x 1/3)) *epsilon*) '(0 1 0 0))
          ((< (abs (- x 2/3)) *epsilon*) '(0 0 1 0))
          ((< (abs (- x 1)) *epsilon*) '(0 0 0 1))
          (t (let* ((weights (iter (for i from 0 to 3)
                                   (for wi in w)
                                   (collect (/ wi (- x (/ i 3))))))
                    (sum (reduce #'+ weights)))
               (mapcar (lambda (w) (/ w sum)) weights))))))

(defun lagrange-curve (points x)
  (let ((weights (uniform-cubic-barycentric-lagrange-weights x))
        (result '(0 0)))
    (iter (for wi in weights)
          (for p in points)
          (setf result (v+ result (v* p wi))))
    result))

(defun lagrange-surface (points p)
  "POINTS is a 3x3 array."
  (let ((wx (uniform-cubic-barycentric-lagrange-weights (first p)))
        (wy (uniform-cubic-barycentric-lagrange-weights (second p)))
        (result '(0 0 0)))
    (iter (for i from 0 to 3)
          (iter (for j from 0 to 3)
                (setf result
                      (v+ result
                          (v* (aref points i j) (elt wx i) (elt wy j))))))
    result))

(defun create-grid (density)
  "DENSITY is the sampling rate of *CURVE* for initial grid points."
  (let ((points ())
        (ghosts ())
        (change t))
    (iter (for i from 0 below density)
          (for x = (/ i (1- density)))
          (pushnew (closest-grid-point (bsc-evaluate *curve* x)) points
                   :test #'equal))
    (iter (while change)
          (setf change nil)
          (iter (for p in points)
                (for cp = (closest-point (grid-point p)))
                (for gp = (closest-grid-point-low cp))
                (iter (for dx from -1 to 2)
                      (iter (for dy from -1 to 2)
                            (for d = (list dx dy))
                            (unless (member (v+ gp d) points :test #'equal)
                              (push (v+ gp d) points)
                              (setf change t))))))
    (iter (for p in points)             ; insert `ghost points'
          (iter (for d in '((-1 0) (0 1) (1 0) (0 -1)))
                (unless (member (v+ p d) points :test #'equal)
                  (pushnew (v+ p d) ghosts :test #'equal))))
    (list points ghosts)))

(defun stencil-parameter (p base)
  (flet ((point (d) (grid-point (v+ base d))))
    (let ((bp (point '(0 0))))
      (list (+ 1/3 (/ (- (first p) (first bp))
                      (- (first (point '(1 0))) (first bp)) 3))
            (+ 1/3 (/ (- (second p) (second bp))
                      (- (second (point '(0 1))) (second bp)) 3))))))

(defun interpolation-matrix (grid)
  (let* ((n (length (first grid)))
         (m (+ n (length (second grid))))
         (result (make-array (list m n) :initial-element 0.0d0)))
    (iter (for p in (append (first grid) (second grid)))
          (for i upfrom 0)
          (for cp = (closest-point (grid-point p)))
          (for base = (closest-grid-point-low cp))
          (for param = (stencil-parameter cp base))
          (for wx = (uniform-cubic-barycentric-lagrange-weights (first param)))
          (for wy = (uniform-cubic-barycentric-lagrange-weights (second param)))
          (iter (for j from -1 to 2)
                (iter (for k from -1 to 2)
                      (for d = (list j k))
                      (for l = (position (v+ base d) (first grid) :test #'equal))
                      (setf (aref result i l)
                            (* (elt wx (1+ j)) (elt wy (1+ k)))))))
    result))

(defun derivation-matrix (grid)
  (let* ((n (length (first grid)))
         (m (+ n (length (second grid))))
         (result (make-array (list n m) :initial-element 0.0d0))
         (ghost-grid (append (first grid) (second grid))))
    (iter (for p in (first grid))
          (for i upfrom 0)
          (iter (for d in '((-1 0) (0 1) (1 0) (0 -1)))
                (for l = (position (v+ p d) ghost-grid :test #'equal))
                (setf (aref result i l) 1.0d0))
          (setf (aref result i i) -4.0d0))
    result))

(defun matrix/diagonal (matrix type)
  "TYPE is :DIAGONAL or :NO-DIAGONAL.
In the :DIAGONAL case, a square matrix is returned."
  (let* ((n (if (eq type :diagonal)
                (let ((n (reduce #'min (array-dimensions matrix))))
                  (list n n))
                (array-dimensions matrix)))
         (result (make-array n)))
    (iter (for i from 0 below (first n))
          (iter (for j from 0 below (second n))
                (setf (aref result i j)
                      (ecase type
                        (:diagonal (if (= i j) (aref matrix i j) 0.0))
                        (:no-diagonal (if (= i j) 0.0 (aref matrix i j)))))))
    result))

(defun laplace-matrix (grid)
  (let ((im (interpolation-matrix grid))
        (dm (derivation-matrix grid)))
    (matrix:m+ (matrix/diagonal dm :diagonal)
               (matrix:multiplication (matrix/diagonal dm :no-diagonal) im))))

(defun matrix-as-sequence (matrix)
  (iter (for i from 0 below (array-dimension matrix 0))
        (collect (iter (for j from 0 below (array-dimension matrix 1))
                       (collect (aref matrix i j))))))

(defun eigenvalues/vectors (matrix)
  (let ((m (grid:make-foreign-array 'double-float :dimensions 2
                                    :initial-contents (matrix-as-sequence matrix))))
    (gsll:eigenvalues-eigenvectors-nonsymm m)))

(defun grid-to-seq (grid)
  (let ((n (grid:dim0 grid)))
    (iter (for i from 0 below n)
          (collect (grid:gref grid i)))))

(defun grid-to-seq-seq (grid)
  (let ((n (grid:dim0 grid))
        (m (grid:dim1 grid)))
    (iter (for i from 0 below m)
          (collect (iter (for j from 0 below n)
                         (collect (grid:gref grid j i)))))))

;;; Tests

;; (defparameter *grid* (create-grid 100))
;; (defparameter *laplace* (laplace-matrix *grid*))
;; (multiple-value-bind (l v) (eigenvalues/vectors *laplace*)
;;   (defparameter *eigenvalues* (grid-to-seq l))
;;   (defparameter *eigenvectors* (grid-to-seq-seq v)))

(defun show-problem-setting (grid)
  (with-open-file (s "/tmp/curve" :direction :output :if-exists :supersede)
    (iter (for i from 0 below 100)
          (for u = (/ i 99))
          (format s "~{~f~^ ~}~%" (bsc-evaluate *curve* u))))
  (with-open-file (s "/tmp/grid" :direction :output :if-exists :supersede)
    (iter (for p in (first grid))
          (format s "~{~f~^ ~}~%" (grid-point p))))
  (with-open-file (s "/tmp/ghost" :direction :output :if-exists :supersede)
    (iter (for p in (second grid))
          (format s "~{~f~^ ~}~%" (grid-point p)))))

(defun show-eigenvector (grid eigen n)
  (with-open-file (pos (format nil "/tmp/positive~d" n) :direction :output :if-exists :supersede)
    (with-open-file (neg (format nil "/tmp/negative~d" n) :direction :output :if-exists :supersede)
      (iter (for i from 0 below (length (second eigen)))
            (if (>= (realpart (elt (elt (second eigen) n) i)) 0)
                (format pos "~{~f~^ ~}~%" (grid-point (elt (first grid) i)))
                (format neg "~{~f~^ ~}~%" (grid-point (elt (first grid) i)))))))
  (elt (first eigen) n))

(defun sort-eigenvalues (eigen)
  "Bubble sort by only the real part."
  (let ((values (coerce (first eigen) 'vector))
        (vectors (coerce (second eigen) 'vector)))
    (iter (for i from (length values) downto 2)
          (iter (for j from 1 below i)
                (when (< (abs (realpart (elt values j)))
                         (abs (realpart (elt values (1- j)))))
                  (rotatef (elt values (1- j)) (elt values j))
                  (rotatef (elt vectors (1- j)) (elt vectors j)))))
    (list (coerce values 'list) (coerce vectors 'list))))

(defun select-eigenvalues (lst n &optional previous (i 0) (tolerance 0.001d0))
  (unless (= n 0)
    (if (> (abs (imagpart (first lst))) *epsilon*)
        (select-eigenvalues (rest lst) n previous (1+ i))
        (let ((next (realpart (first lst))))
          (if (and previous (< (abs (- next previous)) tolerance))
              (select-eigenvalues (rest lst) n previous (1+ i))
              (cons i (select-eigenvalues (rest lst) (1- n) next (1+ i))))))))

(defun generate-gnuplot-tests ()
  (let* ((grid (create-grid 100))
         (laplace (laplace-matrix grid))
         (eigen (multiple-value-bind (l v) (eigenvalues/vectors laplace)
                  (sort-eigenvalues (list (grid-to-seq l) (grid-to-seq-seq v))))))
    (show-problem-setting grid)
    (with-open-file (s "/tmp/commands" :direction :output :if-exists :supersede)
      (format s "set terminal png~%")
      (format s "set output \"/tmp/grid.png\"~%")
      (format s "plot [-1:1] [-1:1] \"/tmp/grid\" title \"Grid\", \"/tmp/curve\" with lines title \"Curve\", \"/tmp/ghost\" title \"Ghost\"~%")
      (dolist (i (select-eigenvalues (first eigen) 6))
        (show-eigenvector grid eigen i)
        (format s "set output \"/tmp/eigen~d.png\"~%" i)
        (format s "plot [-1:1] [-1:1] \"/tmp/negative~d~:*\" title \"Negative\", \"/tmp/curve\" with lines title \"Curve\", \"/tmp/positive~d~:*\" title \"Positive\"~%" i)))))