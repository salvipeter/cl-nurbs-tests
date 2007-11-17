(in-package :cl-nurbs)

(defparameter *curve*
  (make-bspline-curve 3 '(0 0 0 0 0.211829 0.317618 0.446021 0.585312 0.722782 0.826795 1 1 1 1 ) '( (0.00865101 -2.0696) (0.0753516 -1.76613) (-0.253381 -0.590111) (0.90394 -0.517662) (1.66654 -0.470963) (2.54725 -0.495634) (3.45366 -0.289979) (4.14526 -1.01476) (4.17136 -1.65527) (4.09956 -2.02564))))

(three-curve-original *curve* 0.2 0.8 4 "results/original.rdn")

(three-curve-test *curve* 0.2 0.8 4 "results/faired-001.rdn"
		  :resolution 100 :iteration 100 :max-deviation 1000
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :no-fairing nil :simple-fitting t)

(defparameter *curve*
  (make-bspline-curve 3 '( 0 0 0 0 0.315794 0.44019 0.51789 0.585276 0.694806 1 1 1 1 ) '( (-0.10132 0.992674) (0.275082 0.412004) (0.901523 -0.583696) (2.03754 -0.979336) (2.67418 -1.1159) (3.24794 -0.9061) (4.24928 -0.418556) (4.96288 0.339319) (5.2066 1.01465) )))

(three-curve-original *curve* 0.343434 0.656566 4 "results/original.rdn")

(three-curve-test *curve* 0.343434 0.656566 4 "results/csak-sima-fit.rdn"
		  :resolution 100
		  :no-fairing t :simple-fitting t)

(three-curve-test *curve* 0.343434 0.656566 4 "results/csak-szoritott-fit.rdn"
		  :resolution 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing t :simple-fitting nil)

(three-curve-test *curve* 0.343434 0.656566 4 "results/fair-szoritott-fit.rdn"
		  :resolution 100 :iteration 100 :max-deviation 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing nil :simple-fitting nil)
