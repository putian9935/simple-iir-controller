"""Assume all frequencies are normalized to sampling frequency, i.e. T_s = 1. 
"""
from numpy.polynomial import Polynomial as P

def tustin(roots):
    """Apply Tustin's method to a polynomial in its root form.

    Returns
    ------
    - ret: the resulting polynomial 
    """
    ret = P([1])
    for r in roots:
        ret *= P([-r-2, -r+2])
    return ret 

def rational2coef(num, denom):
    """Given u(z)/e(z) in rational form, find IIR coef 
    Parameters
    ---
    - num: coefficient of numerator
    - denom: coefficient of denominator, the order should be at least 2 for stable response 

    Returns
    ---
    A tuple of length 2 
    - first elem.: a NumPy array of coefficient for past errors(last elem. for latest error)
    - second elem.: a NumPy array of coefficient for past outputs
    """
    
    return (num/denom[-1], -denom[:-1]/denom[-1])

class DiscreteTransferTustin: 
    def __init__(self, inf_gain, zeroes, poles):
        self.zeroes = zeroes 
        self.poles = poles 
        self.num = P([inf_gain]) * tustin(zeroes)
        self.denom= P([1]) * tustin(poles)
        if len(zeroes) > len(poles):
            self.denom *= P([1, 1]) ** (len(zeroes) - len(poles))
        elif len(poles) > len(zeroes): 
            self.num *= P([1, 1]) ** (len(poles) - len(zeroes))
            
if __name__ == '__main__':
    zeroes, poles = [-1/1000, -1/200], [-1/5000, -1/5000, -1/2]
    inf_gain = -poles[-1]

    
    # inf_gain = 0.24
    # inf_gain = 1
    # zeroes, poles = [-1/200], [-1/5000] # pi

    # inf_gain = .01
    # zeroes, poles = [-1/20], [-1/8, 0]
    # zeroes, poles = [-1/1000, -1/30], [0, 0] # pii
    # zeroes, poles = [-1/20, -1/4], [0, -1/2.5,-1/2] # pid+cutoff 
    # zeroes, poles = [-1/2000, -1/300, -1/20], [0, 0, -1/2] # piid+cutoff 

    discrete_transfer_func = DiscreteTransferTustin(inf_gain, zeroes, poles)
    print(list(rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef)[0]))
    print(list(rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef)[1]))
    