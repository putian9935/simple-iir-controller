import numpy as np
from collections import deque 
from tqdm import tqdm 
import matplotlib.pyplot as plt 


from numba import jit

def inner(la, lb):
    """Inner product of two lists 
    """
    return sum(a * b for a, b in zip(la, lb)) 
    

@jit
def plant_servo(num):
    return 25*np.exp(((num)-500 + np.random.randn()*3)/200) + 191 + np.random.randn()*1
    # return 25*np.exp(((num)-500 + np.random.randn()*0)/200) + 191 + np.random.randn()*0
    # return 25*np.exp((70*np.log(10+20*num)-500 + np.random.randn()*3)/200) + 191 + np.random.randn()*1

class IIRRealSimulator:
    def __init__(self, err_coef, output_coef):
        self.err_coef = err_coef  
        self.output_coef = output_coef

    def simulate(self, noises):
        ret = [0]
        ret_err = []
        err = deque([0] * len(self.err_coef))
        output = deque([0] * len(self.output_coef))
        for e in noises:
            err.popleft()
            err.append(e + plant_servo(output[-1])) 
            ret_err.append(err[-1])
            new_out = inner(err, self.err_coef) + inner(output, self.output_coef)
            new_out = min(max(new_out, 0), 8000)
            ret.append(new_out)
            output.popleft()
            output.append(new_out)
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
            err.append(e + plant_servo(ret[-1])) 
            ret_err.append(err[-1])
            new_out = inner(err, self.err_coef) + inner(output, self.output_coef)
            new_out = min(max(new_out, 0), 800)
            ret.append(int(new_out))
            output.popleft()
            output.append(new_out)
            if np.isnan(ret_err[-1]):
                break 
            
        return ret_err, ret  


error = -290 * np.ones(10000)

best_num, best_denom =  [0, 0.36, -0.08, -0.44], [0, -.6, 1.6]

# nums = np.linspace(0, 8000)
# plt.plot(nums, plant_servo(nums))
# plt.show()
# exit()


def optimize_num(x):
    simulator = IIRRealSimulator(
        x, best_denom)

    err, out = simulator.simulate_discrete(error)
    ret = np.mean(np.array(err[-1000:])**2)
    print(x, ret)
    return ret

def optimize_denom(x):
    simulator = IIRRealSimulator(
        best_num, x)

    err, out = simulator.simulate_discrete(error)
    ret = np.mean(np.array(err[-1000:])**2)
    print(x, ret)
    return ret

if __name__ == '__main__':
    from tustin import DiscreteTransferTustin, rational2coef
    stds = []
    inf_gains = -10**np.linspace(0, 3, 30)
    
    zeroes, poles = [-1/500, -1/3, ], [0, 0, -1/2] 
    for inf_gain in tqdm(inf_gains):

        inf_gain = -2.5
        discrete_transfer_func = DiscreteTransferTustin(inf_gain, zeroes, poles)
  
        simulator = IIRRealSimulator(*rational2coef(discrete_transfer_func.num.coef, discrete_transfer_func.denom.coef))

        err, out = simulator.simulate_discrete(error)
        # print(np.sum(np.array(err)**2) / 20000)

        from cascade_sim import IIRCascadeRealSimulator 

        simulator = IIRCascadeRealSimulator(zeroes, poles, inf_gain, plant_servo)
        errs = []
        for e in error:
            errs.append(simulator.update(e))
        plt.plot(errs)
        plt.plot(err)
        plt.show()
        exit()

        stds.append(np.sum(np.array(err)**2) / 20000)
    plt.loglog(-inf_gains, stds)
    # plt.gca().twinx()
    # plt.plot(out, c='C1')
    plt.show()

