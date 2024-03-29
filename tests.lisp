;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs-tests)

(defparameter *curve*
  (make-bspline-curve 3 '( 0 0 0 0 0.315794 0.44019 0.51789 0.585276 0.694806 1 1 1 1 ) '( (-0.10132 0.992674) (0.275082 0.412004) (0.901523 -0.583696) (2.03754 -0.979336) (2.67418 -1.1159) (3.24794 -0.9061) (4.24928 -0.418556) (4.96288 0.339319) (5.2066 1.01465) )))

(three-curve-original *curve* 0.343434 0.656566 4 "results/original.rdn")

(three-curve-test *curve* 0.343434 0.656566 4 "results/fair-szoritott-fit.rdn"
		  :resolution 100 :iteration 100 :max-deviation 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing nil :simple-fitting nil)

(three-curve-iterative-test *curve* 0.343434 0.656566 4
			    "results/fair-cont-iter.rdn"
			    :resolution 100 :target-iteration 100
			    :simplex-iteration 15 :fairing-iteration 5)

(defparameter *curve*
  (make-bspline-curve 3 '( 0 0 0 0 0.362383 0.489932 0.582052 0.721361 1 1 1 1 ) '( (0.206598 -1.35897) (0.712538 -0.884575) (1.60127 -0.659641) (2.65718 -0.513087) (3.31592 -0.444638) (4.24742 -0.699734) (4.80681 -1.44105) (4.72273 -1.96703) )))

(three-curve-original *curve* 0.30303 0.828283 4 "results/original.rdn")

(three-curve-test *curve* 0.30303 0.828283 4 "results/fair-szoritott-fit.rdn"
		  :resolution 100 :iteration 100 :max-deviation 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing nil :simple-fitting nil)

(three-curve-iterative-test *curve* 0.30303 0.828283 4
			    "results/fair-cont-iter.rdn"
			    :resolution 100 :target-iteration 100
			    :simplex-iteration 15 :fairing-iteration 5)

(defun test-series (xnode filename &rest keys)
  (let ((faired (apply #'five-surface-test xnode "results/temp.rdn"
		       :no-fairing nil :simple-fitting t
		       :patch-corners nil :cutting nil keys))
	(enlarged (apply #'five-surface-test xnode "results/temp.rdn"
			 :no-fairing nil :simple-fitting nil
			 :patch-corners t :cutting nil keys))
	(cut (apply #'five-surface-test xnode "results/temp.rdn"
		    :no-fairing nil :simple-fitting nil
		    :patch-corners t :cutting '(eval zap) keys))
	(g1 (apply #'five-surface-test xnode "results/temp.rdn"
		   :no-fairing nil :simple-fitting nil
		   :patch-corners t :cutting '(eval zap g1) keys)))
    (write-rdn (append (list (first xnode) faired enlarged cut g1)
		       (rest xnode))
	       filename)))

(defparameter *xnode* (read-rbn "models/xnode.rbn"))

(write-rdn *xnode* "results/original.rdn")

(five-surface-test *xnode* "results/fair-szoritott-fit.rdn"
		   :resolution 1000 :iteration 100 :max-deviation 100
		   :loose-tolerance 0.001 :tight-tolerance 0.0001
		   :number-of-held-points 5
		   :no-fairing nil :simple-fitting nil
		   :patch-corners t :cutting '(eval zap g1))

(test-series *xnode* "results/xnode-series.rdn"
	     :resolution 1000 :iteration 100 :max-deviation 100
	     :loose-tolerance 0.001 :tight-tolerance 0.0001
	     :number-of-held-points 5)

(five-surface-iterative-test *xnode* "results/fair-cont-iter.rdn"
			     :resolution 15 :target-iteration 100
			     :simplex-iteration 30 :fairing-iteration 20)

(defparameter *vb1* (read-rbn "models/vb1.rbn"))

(five-surface-test *vb1* "results/fair-szoritott-fit.rdn"
		   :resolution 1000 :iteration 100 :max-deviation 100
		   :loose-tolerance 0.01 :tight-tolerance 0.001
		   :number-of-held-points 5
		   :no-fairing nil :simple-fitting nil
		   :patch-corners t :cutting nil)

(test-series *vb1* "results/vb1-series.rdn"
	     :resolution 1000 :iteration 100 :max-deviation 100
	     :loose-tolerance 0.01 :tight-tolerance 0.001
	     :number-of-held-points 5)

(defparameter *vb2* (read-rbn "models/vb2.rbn"))

(five-surface-test *vb2* "results/fair-szoritott-fit.rdn"
		   :resolution 1000 :iteration 100 :max-deviation 100
		   :loose-tolerance 0.001 :tight-tolerance 0.0001
		   :number-of-held-points 5
		   :no-fairing nil :simple-fitting t
		   :patch-corners t :cutting nil)

(test-series *vb2* "results/vb2-series.rdn"
	     :resolution 1000 :iteration 100 :max-deviation 100
	     :loose-tolerance 0.001 :tight-tolerance 0.0001
	     :number-of-held-points 5)

#+emacs-lisp
(require 'ange-ftp)
#+emacs-lisp
(defun copy-result-to-dense (&optional filename)
  (interactive)
  (let* ((dfile "/home/salvi/project/cl-nurbs/results/fair-szoritott-fit.rdn")
	 (filename (or filename (read-file-name "File name: " nil dfile)))
	 (ange-ftp-ftp-program-args
	  (adjoin "-u" ange-ftp-ftp-program-args :test #'string=)))
    (copy-file filename "/ftp:dense:/Rodin-1.0.3/bss/" t)))
#+emacs-lisp
(global-set-key [f6] #'copy-result-to-dense)
