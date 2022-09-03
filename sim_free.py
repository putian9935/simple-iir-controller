import numpy as np
from collections import deque 
from tqdm import tqdm 
import matplotlib.pyplot as plt 

def inner(la, lb):
    """Inner product of two lists 
    """
    return sum(a * b for a, b in zip(la, lb)) 
    

slope = 4.1e-4 
def plant(v):
    return slope * v 

class IIRRealSimulator:
    def __init__(self, err_coef, output_coef):
        self.err_coef = err_coef  
        self.output_coef = output_coef

    def simulate(self, noises):
        ret = []
        ret_err = []
        err = deque([0] * len(self.err_coef))
        output = deque([0] * len(self.output_coef))
        for e in noises:
            err.popleft()
            err.append(e + plant(output[-1])) 
            ret_err.append(err[-1])
            ret.append(inner(err, self.err_coef) + inner(output, self.output_coef))
            output.popleft()
            output.append(ret[-1])
            if np.isnan(ret_err[-1]):
                break 
            
        return ret_err, ret  

    def simulate_discrete(self, noises):
        ret = [0]
        ret_err = []
        err = deque([0] * len(self.err_coef))
        output = deque([0] * len(self.output_coef))
        for e in noises:
            err.popleft()
            err.append(e + plant(ret[-1])) 
            ret_err.append(err[-1])
            new_out = inner(err, self.err_coef) + inner(output, self.output_coef)
            ret.append(round(new_out, 3))
            output.popleft()
            output.append(new_out)
            if np.isnan(ret_err[-1]):
                break 
            
        return ret_err, ret  


def nm2khz(nm):
    return 3e8*(1/707e-9-1/(707e-9+nm*1e-9))/1e3

def dump_data():
    open("data", "w").write(', '.join("%.16e"%n for n in noise[:5000]))

data = np.genfromtxt("fast_b.lta", delimiter="\t")
noise = (data[:, 1] - data[0,1])

# plt.plot(noise)
# plt.xlabel('Sample index')
# plt.ylabel(r'Wavelength difference $\Delta \lambda\,/\,\mathrm{nm}$')
# plt.tight_layout()
# plt.savefig('free-running.pdf')
# plt.show()

if __name__ == '__main__':
    import scipy.signal as signal 
    freq, noise_psd = signal.welch(noise, nperseg=1024)
    # plt.loglog(freq[1:], noise_psd[1:])


    # for gain in range(0, 25, 5):
    from tustin import DiscreteTransferTustin, rational2coef

    # inf_gain = -30
    # zeroes, poles = [-1/0.02], [0] # pi
    # inf_gain = -200
    # zeroes, poles = [-1/0.3], [0] # pi

    # discrete_transfer_func = DiscreteTransferTustin(inf_gain, zeroes, poles)
    # simulator = IIRRealSimulator(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

    # print(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

    # err, out = simulator.simulate(noise)
    # freq, err_psd = signal.welch(err, nperseg=1024)
    # plt.loglog(freq[1:], err_psd[1:])
    # print(nm2khz(np.std(err)))

    # inf_gain = -1800*80
    # zeroes, poles = [-1/0.25, -1/0.05], [0, -1/0.016, -1/0.016] # pid+cutoff 

    inf_gain = -1000
    zeroes, poles = [-1/5, -1/5], [0, -1/1000, -1/2] # pii

    discrete_transfer_func = DiscreteTransferTustin(inf_gain, zeroes, poles)

    simulator = IIRRealSimulator(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

    print(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

    err, out = simulator.simulate(noise)

    np.savetxt("err", err, delimiter=",")
    w, h = signal.freqz_zpk(*signal.bilinear_zpk(zeroes, poles, inf_gain, fs=1), worN=16384)
    freq, err_psd = signal.welch(err, nperseg=1024)
    # plt.loglog(freq[1:], err_psd[1:])

    print(nm2khz(np.std(err)))
    plt.figure()
    plt.plot(err)
    plt.plot(noise)

    
    inf_gain = -1000
    zeroes, poles = [-1/5, ], [0, -1/2] # pi

    discrete_transfer_func = DiscreteTransferTustin(inf_gain, zeroes, poles)

    simulator = IIRRealSimulator(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))
    print(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

    err2, out = simulator.simulate(noise)

    w, h = signal.freqz_zpk(*signal.bilinear_zpk(zeroes, poles, inf_gain, fs=1), worN=16384)
    freq, err_psd = signal.welch(err, nperseg=1024)
    # plt.loglog(freq[1:], err_psd[1:])

    plt.plot(err2)
    plt.plot(np.array(err2) - np.array(err))
    print((np.array(err2) - np.array(err))[-10:])
    print(nm2khz(np.std(err2)))
    plt.show()



    simulator = IIRRealSimulator([-13.48164908, -1200.86011487, -2105.17165443], [0.57069787, 0.38133636])

    err, out = simulator.simulate(noise)

    w, h = signal.freqz_zpk(*signal.bilinear_zpk(zeroes, poles, inf_gain, fs=1), worN=16384)
    freq, err_psd = signal.welch(err, nperseg=1024)
    plt.loglog(freq[1:], err_psd[1:])
    print(nm2khz(np.std(err)))

    simulator = IIRRealSimulator((6.07561833e+01, -6.08612573e-01,  1.07336639e+02,  1.74027779e+02, -1.19752512e+03, -2.10938405e+03), (0.00287914, 0.02454546, 0.5705591,
      0.36614206))

    err, out = simulator.simulate(noise)

    w, h = signal.freqz_zpk(*signal.bilinear_zpk(zeroes, poles, inf_gain, fs=1), worN=16384)
    freq, err_psd = signal.welch(err, nperseg=1024)
    plt.loglog(freq[1:], err_psd[1:])
    print(nm2khz(np.std(err)))

    plt.show()
