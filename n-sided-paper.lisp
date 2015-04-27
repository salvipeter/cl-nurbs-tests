(in-package :cl-nurbs-tests)

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
