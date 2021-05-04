# Energy testing
N = 1000;
K = 0.3;
T = 0.2;
Imp_1 = Impurity(1, 2, 1)
Imp_2 = Impurity(2, 8, 1)
Imp_3 = Impurity(4, 4, 1)

s = System([Imp_1, Imp_2, Imp_3], T, N, K)

exact_F_I(s) - exact_neg_TS_I(s)
E_I(s)
