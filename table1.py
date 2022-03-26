from sympy.printing import print_latex
from cyclo import cyclo

z = cyclo(11)
print_latex(z[0]+z[-1]) # 2cos(2pi/11)

z = cyclo(13)
print_latex(z[0]+z[-1]) # 2cos(2pi/13)

z = cyclo(17)
print_latex(z[0]+z[-1]) # 2cos(2pi/17)

z = cyclo(19)
print_latex(z[0]+z[-1]) # 2cos(2pi/19)

z = cyclo(37)
print_latex(z[0]+z[-1]) # 2cos(2pi/37)
