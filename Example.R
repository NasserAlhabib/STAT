#=============== Import Functions ===============
source(file.choose())

#=============== Question1 ===============
Before = matrix(c(22,24,20,18,19,16,24,19,25,19,25,28,28,24,22,25,30,25,27,23), ncol = 2, byrow = TRUE)
After = matrix(c(24,25,22,20,19,17,22,18,28,18,26,28,28,26,24,27,30,27,29,24), ncol = 2, byrow = TRUE)

paired_hotelling_T2(Before, After)

#=============== Question2 ===============
Xbar = matrix(c(1.086,2.544,2.851,3.420), ncol = 1)
S = matrix(c(2.902,2.438,2.963,2.183,2.438,3.049,2.775,2.319,2.963,2.775,4.281,2.939,2.183,2.319,2.939,3.162), ncol = 4, byrow = TRUE)
C = matrix(c(1,-1, 0, 0, 0, 1,-1, 0,0, 0, 1,-1), nrow = 3, byrow = TRUE)
hotelling_T2_contrast(Xbar,S,40,C)

#=============== Question3 ===============
X1 = c(96,92.3)
X2 = c(90.25,88.25)
X3 = c(77.33,73.33)
S1 = matrix(c(1,0.5,0.5,0.33), ncol = 2, byrow = TRUE)
S2 = matrix(c(30.92,-8.42,-8.42,4.91), ncol = 2, byrow = TRUE)
S3 = matrix(c(4.33,1.83,1.83,2.33), ncol = 2, byrow = TRUE)

n = c(3,4,3)

Xbars = list(X1,X2,X3)
Ss = list(S1,S2,S3)

stat438_oneway_manova_summary(Xbars, Ss, n)

#=============== Question4 ===============
z1 = c(2,3,3,6,7,9)
y1 = c(10,5,7,19,11,18)
y2 = c(15,9,3,25,7,13)

X = cbind(z1)    
Y = cbind(y1, y2)    

stat438_reg_exam(Y = Y, X = X)

#=============== Question5 ===============
S = matrix(c(5691,600,217,600,126,24,217,24,23), ncol = 3, byrow = TRUE)
pca_exam_from_R(
  R = S,
  n = 87,
  alpha = 0.01,
  ci_method = "bonferroni_z",
  ci_for = 1:3,
  var_names = c("X1", "X2", "X3"),
  digits = 7
)

pca_exam_from_R(S, n = 87, ci_method = "bonferroni", alpha = 0.01, ci_for =  1:3)



