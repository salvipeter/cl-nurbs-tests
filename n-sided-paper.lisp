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

;;; Modules
(load (compile-file "blends"))
(load (compile-file "kato-test"))
(load (compile-file "ribbons"))
(load (compile-file "patches"))
