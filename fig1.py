import numpy as np
import matplotlib.pyplot as plt

n,g,M = 5,2,64 # fig1
# n,g,M = 7,3,64 # fig2

c = np.exp(2j*np.pi*np.arange(n)/n)
d = np.exp(2j*np.pi*np.arange(M+1)/M)

plt.plot(np.real(d), np.imag(d), 'k--', lw=1)
plt.plot(np.real(c), np.imag(c), 'ko', ms=10)

plt.text(1.05,0,'1',ha='left',va='center',fontsize=20)

for k in range(n-1):
    i = g**k%n
    s = r'$\zeta_' + str(k) + r' = \zeta'
    s+= (r'^' + str(i) + r'$') if k else r'$'
    x,y = np.real(c[i]), np.imag(c[i])
    ha = 'left' if x>0 else 'right'
    va = 'bottom' if y>0 else 'top'
    x += 0.05 if x>0 else -0.05
    plt.text(x,y,s,ha=ha,va=va,fontsize=20)

plt.axis('equal')
plt.axis('off')
plt.axis([-1.1,1.5,-1.1,1.1])
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
