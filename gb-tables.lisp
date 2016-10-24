(in-package :cl-nurbs-tests)

;;; Deficiency table
#+nil
(let ((*barycentric-d-function* #'barycentric-d-triangle)
      (*fullness-height-function* #'fullness-height-parabola)
      (*barycentric-dilation* 5)
      (delta-setter-function nil)
      (*deficiency-function* #'deficiency-squared-central-layer)
      (*deficiency-use-d* t))
  (dprint *barycentric-dilation* *fullness-height-function* *barycentric-dilation*
          delta-setter-function *deficiency-function* *deficiency-use-d*)
  (iter (for n from 3 to 3)
        (when delta-setter-function
          (setf *barycentric-dilation* (funcall delta-setter-function n)))
        (format t "|~a|delta=|~,3f|~%" n *barycentric-dilation*)
        (iter (for d from 3 to 10)
              (format t "|~a|~a|~,3f|~%" n d (funcall *deficiency-function* n d :position '(0.3 0.5))))))

;;; Negativity / Monotonity / Target deficiency
#+nil
(let ((*barycentric-d-function* #'barycentric-d-triangle)
      (*fullness-height-function* #'fullness-height-parabola)
      (*barycentric-dilation* 5.25)
      (delta-setter-function nil)
      (*deficiency-function* #'deficiency-squared-central-layer)
      (*deficiency-use-d* t)
      (*resolution* 50)
      (*epsilon* 1.0d-8))
  (dprint *barycentric-dilation* *fullness-height-function* *barycentric-dilation*
          delta-setter-function *deficiency-function* *deficiency-use-d* *resolution* *epsilon*)
  (iter (for n from 3 to 3)
        (for points = (points-from-angles (uniform-angles n)))
        (when delta-setter-function
          (setf *barycentric-dilation* (funcall delta-setter-function n)))
        (format t "|~a|delta=|~,3f|~%" n *barycentric-dilation*)
        (iter (for d from 3 to 10)
              (format t "Checking d = ~a...~%" d)
              (iter (for p in (vertices points))
                    (for x = (funcall *deficiency-function* n d :position p))
                    (when (< x (- *epsilon*))
                      (warn "Negative: ~f at ~{~,3f, ~,3f~}" x p))))))

;;; Boundary search
#+nil
(let ((*barycentric-d-function* #'barycentric-d-pisti-all)
      (*fullness-height-function* #'fullness-height-circle)
      (*deficiency-function* #'deficiency-squared-central-layer-nondiagonal)
      (*deficiency-use-d* t)
      (*resolution* 50)
      (*epsilon* 1.0d-5)
      (type 'monotonity)                  ; negativity / monotonity / target
      (min -10.0)
      (max 10.0)
      (target 0.0)
      (iterations 20))
  (dprint *barycentric-dilation* *fullness-height-function* *deficiency-function*
          *deficiency-use-d* *resolution* *epsilon* min max target iterations)
  (format t "Delta values for ~a:~%" type)
  (iter (for n from 5 to 10)
        (iter (for d from 3 to 10)
              (format t "|~a|~a|~,3f|~%" n d
                      (ecase type
                        (negativity
                         (find-dilation-negative-boundary n d min max :iterations iterations))
                        (monotonity
                         (find-dilation-monotone-boundary n d min max :iterations iterations))
                        (target
                         (find-dilation-for-deficiency n d :min min :max max :target target
                                                       :iterations iterations)))))))
