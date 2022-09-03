from math import pi 
class IIRFirstOrderSimulator:
    def __init__(self, z, p):
        self.z = z
        self.p = p
        self.le = self.lo = 0

    def update(self, e):
        ret = ((1 + pi * self.p) * self.lo - (1 + pi * self.z) * self.le + (1 - pi * self.z) * e) / (1 - pi * self.p)
        self.lo, self.le = ret, e
        return ret


class IIRSinglePoleSimulator:
    def __init__(self, p):
        self.p, self.le, self.lo = p, 0, 0

    def update(self, e):
        ret = ((1 + pi * self.p) * self.lo + pi * (self.le + e)) / (1 - pi * self.p)
        self.lo, self.le = ret, e
        return ret


class IIRCascadeSimulator:
    def __init__(self, zs, ps, gain):
        self.gain = gain
        self.zs = sorted(zs, key=lambda _: -_)
        self.ps = sorted(ps, key=lambda _: -_)
        if len(zs) > len(ps):
            raise ValueError('There should be more poles than zeroes')
        if len(zs) and self.zs[0] > 0:
            raise ValueError('Zeroes should be in LHP')
        if len(ps) and self.ps[0] > 0:
            raise ValueError('Poles should be in LHP')

        self.simulators = []
        for z, p in zip(self.zs, self.ps):
            self.simulators.append(IIRFirstOrderSimulator(z, p))
        if len(self.ps) > len(self.zs):
            for p in self.ps[-len(self.ps)+len(self.zs):]:
                self.simulators.append(IIRSinglePoleSimulator(p))

    def update(self, e):
        for sim in self.simulators:
            e = sim.update(e)
        return e * self.gain

class IIRCascadeRealSimulator(IIRCascadeSimulator):
    def __init__(self, zs, ps, gain, plant, lower=-32768, upper=32767):
        super().__init__(zs, ps, gain)
        self.plant = plant 
        self.lo = 0
        self.lower = lower 
        self.upper = upper 
    
    def update(self, e): 
        ret = e + self.plant(self.lo)        
        new_out = super().update(ret)
        self.lo = int(min(max(new_out, self.lower), self.upper))
        return ret


if __name__ == '__main__':
    def dump_result(errors):
        import struct 
        open('from_py', 'wb').write(b''.join(struct.pack('d', e) for e in errors[:5000]))
    
    import numpy as np
    import matplotlib.pyplot as plt
    from sim_free import noise
    # the plant
    def plant(dac_num):
        return 4.1e-4  * dac_num * 6 / 65536

    zeroes, poles, gain  = [-1/5], [0, -1/3], -250 # specs for analog filter 
    iir = IIRCascadeRealSimulator(zeroes, poles, gain * 65536 / 6, plant)  # define the controller 

    errors = []
    for e in noise: 
        errors.append(iir.update(e))
        
    # plt.plot(noise, label='Free-running') 
    # plt.plot(errors, label='Locked') 

    # plt.xlabel('Sample index')
    # plt.ylabel(r'Wavelength difference $\Delta \lambda\,/\,\mathrm{nm}$')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('locked.pdf')
    # plt.show()

    import scipy.signal as signal
    freq, err_psd = signal.welch(errors, nperseg=1024)
    plt.semilogy(freq[1:], err_psd[1:], label='Locked (PI)')

    zeroes, poles, gain  = [-1/20, -1/20, -1/2.5], [0, 0, -1/3, -1/2], -800 # specs for analog filter 
    iir = IIRCascadeRealSimulator(zeroes, poles, gain * 65536 / 6, plant)  # define the controller 
    errors = []
    for e in noise: 
        errors.append(iir.update(e))

    freq, err_psd = signal.welch(errors, nperseg=1024)
    plt.semilogy(freq[1:], err_psd[1:], label='Locked (PIID)')
    freq, err_psd = signal.welch(noise, nperseg=1024)
    plt.semilogy(freq[1:], err_psd[1:], label='Free-running')

    plt.xlabel('Power spectral density$\,/\,\mathrm{nm}^2$')
    plt.ylabel(r'Normalized frequency')
    plt.legend()
    plt.tight_layout()
    plt.savefig('locked_freq.pdf')
    plt.plot()
    plt.show()



