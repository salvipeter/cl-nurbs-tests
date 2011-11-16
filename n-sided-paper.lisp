(in-package :cl-nurbs-tests)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (unless (find-package "GSLL")
    (defpackage gsll (:export make-multi-dimensional-root-solver-fdf
			      +gnewton-mfdfsolver+ iterate
			      multiroot-test-residual solution
			      polynomial-solve)))
  (unless (find-package "GRID")
    (defpackage grid (:export make-foreign-array gref))))

;;; Utils
(defmacro acond (&rest clauses)
  (if (null clauses)
      nil
      (let ((cl1 (car clauses))
            (sym (gensym)))
        `(let ((,sym ,(car cl1)))
           (if ,sym
               (let ((it ,sym)) ,@(cdr cl1))
               (acond ,@(cdr clauses)))))))

(defmacro dlet* (bindings &body body)
  (cond ((null bindings) `(progn ,@body))
	((atom (car (first bindings)))
	 `(let (,(first bindings))
	    (dlet* ,(rest bindings) ,@body)))
	(t `(destructuring-bind ,(car (first bindings))
		,(cadr (first bindings))
	      (dlet* ,(rest bindings) ,@body)))))

;;; Modules
(load (compile-file "blends"))
(load (compile-file "kato-test"))
(load (compile-file "ribbons"))
(load (compile-file "patches"))
