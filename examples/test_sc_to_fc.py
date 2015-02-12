import spherepy as sp

c = sp.zeros_coefs(5, 5)
c[4, 3] = 1.0
c[3, 0] = 1.0

p = sp.ispht(c, 804, 804)

