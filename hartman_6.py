from einsteinpy.symbolic.riemann import RiemannCurvatureTensor, ChristoffelSymbols
from einsteinpy.symbolic import MetricTensor
from einsteinpy.symbolic.ricci import RicciTensor, RicciScalar
import sympy
from sympy import cosh, sin, cos

# See link for instructions
# https://einsteinpy.org/

sympy.init_printing()

tau, psi, theta, phi = sympy.symbols("tau psi theta phi")

l = sympy.symbols('l')

g_tau_tau = l**2
g_psi_si = l**2*cos(tau)**2
g_theta_theta = l**2*cos(tau)**2*sin(psi)**2
g_phi_phi = l**2*cos(tau)**2*sin(psi)**2*sin(theta)**2

m = sympy.diag(g_tau_tau, g_psi_si, g_theta_theta, g_phi_phi).tolist()
metric = MetricTensor(m, syms = (tau, psi, theta, phi))

chr = ChristoffelSymbols.from_metric(metric)
# print(chr.config)
# print(chr)

rm = RiemannCurvatureTensor.from_christoffels(chr)
# print(rm.config)
# print(rm)

rt = RicciTensor.from_riemann(rm)
# print(rt.config)
# print(rt)

rs = RicciScalar.from_riccitensor(rt)
# print(rs)

# print(chr[1, 1, 1]) # Since symbols were given in as tau r phi, r corresponds to 1, tau to 0 and phi to 2 so this is Gamma^r_rr
#
# print()
#
# print(chr[1, 0, 0], chr[1, 2, 2])

print(rs)

