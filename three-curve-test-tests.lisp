(in-package :cl-nurbs)

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

(defparameter *curve*
  (make-bspline-curve 3 '( 0 0 0 0 0.362383 0.489932 0.582052 0.721361 1 1 1 1 ) '( (0.206598 -1.35897) (0.712538 -0.884575) (1.60127 -0.659641) (2.65718 -0.513087) (3.31592 -0.444638) (4.24742 -0.699734) (4.80681 -1.44105) (4.72273 -1.96703) )))

(three-curve-original *curve* 0.30303 0.828283 4 "results/original.rdn")

(three-curve-test *curve* 0.30303 0.828283 4 "results/csak-sima-fit.rdn"
		  :resolution 100
		  :no-fairing t :simple-fitting t)

(three-curve-test *curve* 0.30303 0.828283 4 "results/csak-szoritott-fit.rdn"
		  :resolution 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing t :simple-fitting nil)

(three-curve-test *curve* 0.30303 0.828283 4 "results/fair-szoritott-fit.rdn"
		  :resolution 100 :iteration 300 :max-deviation 100
		  :loose-tolerance 0.01 :tight-tolerance 0.001
		  :number-of-held-points 5
		  :no-fairing nil :simple-fitting nil)
