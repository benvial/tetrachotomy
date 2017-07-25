import sympy as sy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import time
pi=np.pi
plt.ion()
sns.set_context("poster", font_scale=1.)
sns.set_style("white")


plt.ion()
sy.init_printing()

# ---------------------------------
#          h, epsilon, mu
# ---------------------------------

j = sy.I
d, epsilon, mu, x, k = sy.symbols("d epsilon mu x k")
# epsilon2, h2 = sy.symbols("epsilon2, h2")
epsilon_sup, epsilon_sub = sy.symbols("epsilon_sup epsilon_sub")
nlay = 1
# epsilon_list = [1, epsilon, 1]
# x_list = [0, d]

epsilon_list  = list(sy.symbols('epsilon_1:%d'%(nlay+1)))
epsilon_list  = [epsilon_sub] + epsilon_list + [epsilon_sup]

h_list  = list(sy.symbols('h_1:%d'%(nlay +1)))
x_list = [0]
for i in range(nlay):
    x_list.append(x_list[i] + h_list[i])

# nlay = 2
# epsilon_list = [1, epsilon, epsilon2, 1]
# x_list = [0, h, h + h2]

a_plus  = list(sy.symbols('a^+_0:%d'%(nlay + 2)))
a_minus = list(sy.symbols('a^-_0:%d'%(nlay + 2)))
# outgoing wave conditions
a_plus[0],  a_minus[-1]= 0, 0


u = []
dudx = []
for i in range(nlay + 2):
    kr= k*sy.sqrt(epsilon_list[i])
    u.append(a_plus[i]*sy.exp(+j* kr*x) + a_minus[i]*sy.exp(-j* kr*x))
    dudx.append(sy.diff(u[i], x))

# boundary conditions
Eqs = []
for i in range(nlay + 1):
    Eq1 = u[i].replace(x, x_list[i]) - u[i+1].replace(x, x_list[i])
    Eqs.append(Eq1)
    Eq2 = dudx[i].replace(x, x_list[i]) - dudx[i+1].replace(x, x_list[i])
    Eqs.append(Eq2)


a = a_plus[1:] + a_minus[:-1]

M = sy.zeros(2*(nlay+1))
for i1 in range(2*(nlay+1)):
    for i2 in range(2*(nlay+1)):
        M[i1, i2] = Eqs[i1].coeff(a[i2])

det_system = sy.det(M)
k_sol = sy.solveset(det_system, k)
# toreplace = [(epsilon, 2), (d, 100)]

epsilon_np = [1, 2, 5, 7, 1]
h_np = [100, 150, 300]

epsilon_np = [1, 5, 1]
h_np = [100]

htot = np.sum(h_np)
toreplace = [(epsilon_list[0], epsilon_np[0]), (epsilon_list[-1], epsilon_np[-1])]
for i in range(0, nlay):
    toreplace.append( (epsilon_list[i+1], epsilon_np[i+1]))
    toreplace.append( (h_list[i], h_np[i]))

det_subs = det_system.subs(toreplace)
det_f = sy.lambdify(k, k**2/det_subs, "numpy")

# det_f = sy.lambdify(k, det_subs, "numpy")

# k_np = np.linspace(0.001, 0.1, 1111)
# det_np = det_f(k_np)
# plt.clf()
# plt.plot(k_np, np.abs(det_np/k_np**2))

kre, kim = np.meshgrid(np.linspace(0.01, 2, 111), np.linspace(-0.3, 0., 111))

kcplx = kre + 1j*kim

det_cplx = det_f(kcplx*(2*pi/htot))
# det_cplx[np.isnan(det_cplx)] = 1
yplt = np.log10(np.abs(det_cplx))


t0= time.time()
det_f(1)
t1= time.time()-t0
print('t1=', t1)


plt.figure()
fig = plt.gcf()
ax = fig.add_subplot(111)  # , aspect='equal')
cmap = sns.diverging_palette(33, 111, sep=10, n=51, as_cmap=True)
cmap = sns.cubehelix_palette(start=-1, rot=1.2, as_cmap = True)
levels = MaxNLocator(nbins=22).tick_values(yplt.min(), yplt.max() )
# levels = MaxNLocator(nbins=22).tick_values(-2, 5 )


# levels = MaxNLocator(nbins=11).tick_values(-11, 11)
plt.contourf(kre, kim, yplt, cmap=cmap, levels=levels)
plt.colorbar()
plt.xlabel(r'Re $k_n$')
plt.ylabel(r'Im $k_n$')
