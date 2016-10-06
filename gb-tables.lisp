(in-package :cl-nurbs-tests)

(defparameter *delta-setter-function* nil)

(defun gb-print-info ()
  (format t "Parameterization: ~a~%" *barycentric-d-function*)
  (format t "Fullness: ~a~%" *fullness-height-function*)
  (if *delta-setter-function*
      (format t "Delta: ~a~%" *delta-setter-function*)
      (format t "Delta: ~,3f~%" *barycentric-dilation*))
  (format t "Deficiency: ~a~%" *deficiency-function*)
  (format t "Use S: ~a~%" (not *deficiency-use-d*))
  (format t "Resolution: ~a~%" *resolution*)
  (format t "Epsilon: ~f~%" *epsilon*))

;;; Deficiency table
#+nil
(let ((*barycentric-d-function* #'barycentric-d-peti)
      (*fullness-height-function* #'fullness-height-circle)
      (*barycentric-dilation* 0)
      (*delta-setter-function* #'find-optimal-delta)
      (*deficiency-function* #'deficiency-squared-central-layer)
      (*deficiency-use-d* t))
  (gb-print-info)
  (iter (for n from 5 to 10)
        (when *delta-setter-function*
          (setf *barycentric-dilation* (funcall *delta-setter-function* n)))
        (format t "|~a|delta=|~,3f|~%" n *barycentric-dilation*)
        (iter (for d from 3 to 10)
              (format t "|~a|~a|~,3f|~%" n d (funcall *deficiency-function* n d)))))

;;; Negativity
#+nil
(let ((*barycentric-d-function* #'barycentric-d-peti)
      (*fullness-height-function* #'fullness-height-circle)
      (*barycentric-dilation* 0)
      (*delta-setter-function* #'find-optimal-delta)
      (*deficiency-function* #'deficiency-squared-central-layer-nondiagonal)
      (*deficiency-use-d* t)
      (*resolution* 30)
      (*epsilon* 1.0d-8))
  (gb-print-info)
  (iter (for n from 5 to 10)
        (for points = (points-from-angles (uniform-angles n)))
        (when *delta-setter-function*
          (setf *barycentric-dilation* (funcall *delta-setter-function* n)))
        (format t "|~a|delta=|~,3f|~%" n *barycentric-dilation*)
        (iter (for d from 3 to 10)
              (format t "Checking d = ~a...~%" d)
              (iter (for p in (vertices points))
                    (for x = (funcall *deficiency-function* n d :position p))
                    (when (< x (- *epsilon*))
                      (warn "Negative: ~f at ~{~,3f, ~,3f~}" x p))))))

;;; Boundary search
#+nil
(let ((*barycentric-d-function* #'barycentric-d-peti)
      (*fullness-height-function* #'fullness-height-circle)
      (*barycentric-dilation* 0)
      (*delta-setter-function* #'find-dilation-negative-boundary)
      (*deficiency-function* #'deficiency-squared-central-layer)
      (*deficiency-use-d* t)
      (*resolution* 50)
      (*epsilon* 1.0d-5)
      (type 'target)                  ; negativity / monotonity / target
      (min -20.0)
      (max 20.0)
      (target 0.0)
      (iterations 20))
  (gb-print-info)
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
