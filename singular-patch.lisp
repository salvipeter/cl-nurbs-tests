(in-package :cl-nurbs-tests)

(defun generate-singular-patch (points width alpha)
  "WIDTH is the the ribbon width - ribbons are grown perpendicular to the sides.
ALPHA is the angle with the xy plane (in radians)."
  (let* ((lines (lines-from-points points))
         (center (central-point points lines t)))
    (list (iter (for (p q) in lines)
                (collect (list (append p '(0)) (append q '(0)))))
          (iter (for (p q) in lines)
                (for dir = (v- q p))
                (for d = (list (second dir) (- (first dir))))
                (when (< (scalar-product d (v- center p)) 0)
                  (setf d (v* d -1)))
                (for xy = (v* (vnormalize d) width (cos alpha)))
                (for z = (* width (sin alpha)))
                (for xyz = (append xy (list z)))
                (collect (list (v+ (append p '(0)) xyz)
                               (v+ (append q '(0)) xyz)))))))

#+nil
(generate-singular-patch (points-from-angles (uniform-angles 5)) 0.5 (/ (* 45 pi) 180))

(defun write-singular-patch (points type filename
                             &key (width 0.5) (alpha #.(/ (* 45 pi) 180))
                               (distance-type 'perpendicular) (output 'vtk))
  "Similar to WRITE-PATCH in patches.lisp,
but all ribbons have a fixed width and a fixed angle with the xy plane.
See the documentation of GENERATE-SINGULAR-PATCH for the setting of these values.
OUTPUT is one of (SPIDER RIBBONS PATCH)."
  (let ((n (length points))
        (patch (generate-singular-patch points width alpha)))
    (ecase output
      (spider
       (write-vtk-polylines
        (iter (for line in (spider-lines points))
              (collect (iter (for domain-point in line)
                             (collect (patch-evaluate patch points type distance-type
                                                      domain-point)))))
        filename))
      (ribbons
       (write-vtk-curves
        (iter (for curve1 in (first patch))
              (for curve2 in (second patch))
              (for points1 =
                   (iter (for u from 0 to 1 by (/ *resolution*))
                         (collect (bezier curve1 u))))
              (for points2 =
                   (iter (for u from 0 to 1 by (/ *resolution*))
                         (collect (bezier curve2 u))))
              (appending (append (list points1 points2)
                                 (mapcar #'list points1 points2))))
        filename))
      (patch
       (write-ply-indexed-mesh
        (iter (for domain-point in (vertices points))
              (collect (patch-evaluate patch points type distance-type domain-point)))
        (triangles n) filename)))))

#+nil
(let ((*resolution* 50))
  (write-singular-patch (points-from-angles (uniform-angles 5)) 'ribbon "/tmp/proba.ply"
                        :width 0.5 :alpha (/ (* 45 pi) 180) :output 'patch)
  (write-singular-patch (points-from-angles (uniform-angles 5)) 'ribbon "/tmp/proba.vtk"
                        :width 0.5 :alpha (/ (* 45 pi) 180) :output 'spider)
  (write-singular-patch (points-from-angles (uniform-angles 5)) 'ribbon "/tmp/ribbons.vtk"
                        :width 0.5 :alpha (/ (* 45 pi) 180) :output 'ribbons))
