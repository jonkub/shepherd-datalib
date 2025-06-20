from commons import generate_h5_files, EEnvGenerator, STEP_WIDTH
from pathlib import Path
import numpy as np
import math

from scipy.special import lambertw

class PVModel:
    """
    PV-model based on a (photo-)current source, diode, shunt resistor
    and series resistor.
    """

    def get_cs(self, vs: np.ndarray, intensities: np.ndarray):
        pass


    def get_c(self, v_hrv, i_pv):
        # TODO actually verify the solution in step 2
        """
        Computes the harvesting current for the given photocurrent and harvesting voltage.

        The function was obtained in the following steps:
        1. Analyse the circuit to obtain:
        i_pv - v_d / r_sh - (v_d - v_hrv) / r_s = i_s * (exp(v_d / (n * V_t)) - 1)
        2. Use wolframalpha.com to solve for i_hrv(v_hrv, i_pv):
        https://www.wolframalpha.com/input?i=solve+k_1+*+%28+e+%5E+%28x+%2F+k_2%29+-+1%29+%3D+a+-+x+%2F+k_3+-+%28x+-+b%29+%2F+k_4+for+x%28a%2C+b%29
        -> x = (k_3 (-k_2 W_n((exp((k_3 (b + (a + k_1) k_4))/(k_2 (k_3 + k_4))) k_1 k_3 k_4)/(k_2 (k_3 + k_4))) + k_4 (a + k_1) + b) - k_2 k_4 W_n((exp((k_3 (b + (a + k_1) k_4))/(k_2 (k_3 + k_4))) k_1 k_3 k_4)/(k_2 (k_3 + k_4))))/(k_3 + k_4) and k_3 k_4!=0 and k_2 (k_3 + k_4)!=0 and k_1 k_3 k_4!=0 and n element Z 
        where: a = i_pv; b = v_hrv; k_1 = i_s; k_2 = n * v_t; k_3 = r_sh; k_4 = r_s
        """

        # TODO obtain realistic params
        # TODO move this to member variables
        r_sh = 1
        r_s = 1

        # Diode
        n = 1 # ideality factor
        v_t = 25.852e-3 # thermal voltage (27C used here)
        i_s = 1e-6 # reverse-bias saturation current

        tmp_1 = (math.exp((r_sh*(v_hrv+(i_pv+i_s)*r_s))/(n*v_t*(r_sh+r_s)))*i_s*r_sh*r_s)/(n*v_t*(r_sh+r_s))
        for k in range(10000):
            tmp_2 = lambertw(tmp_1, k=k)
            if np.isreal(tmp_2):
                break
        if not np.isreal(tmp_2):
            raise RuntimeError("")
        tmp_2 = tmp_2.real

        return (r_sh*(-n*v_t*tmp_2+r_s*(i_pv+i_s)+v_hrv)-n*v_t*r_s*tmp_2)/(r_sh+r_s)


# TODO evaluate switching from traces to surfaces (update filename and generator description accordingly)
# TODO isn't variance irrelevant if we just check rnd() > 0?
class MultivarRndWalk(EEnvGenerator):
    '''
    Generates traces according to a multivariate random walk. Each step, the voltage
    increases or decreases randomly with to the specified correlation and variance.

    Variance is used to control the volatility while correlation controls how similar
    the nodes behave.
    '''

    def __init__(self, node_count, seed, correlation, variance, voltage_min, voltage_max, current):
        super().__init__(node_count, seed)

        self.v_min = voltage_min
        self.v_max = voltage_max
        self.c = current

        self.mean = np.zeros(node_count)
        # Create a correlation matrix with off-diagonal values set to the correlation coefficient
        self.cov_matrix = self.gen_covariance_matrix(node_count, variance, correlation)

        print(self.cov_matrix)

        # Start pattern between the two bounds
        self.states = np.full(node_count, 0.5)

        self.pv_model = PVModel()

    @staticmethod
    def gen_covariance_matrix(dim, variance, correlation):
        # Create a correlation matrix with off-diagonal values set to the correlation coefficient
        correlation_matrix = np.full((dim, dim), correlation)
        np.fill_diagonal(correlation_matrix, 1)  # Diagonal should be 1 (self-correlation)

        # Convert correlation matrix to covariance matrix using variance
        covariance_matrix = correlation_matrix * variance

        return covariance_matrix

    def generate_random_pattern(self, count):
        random = self.rnd_gen.multivariate_normal(self.mean, self.cov_matrix, size=count)
        walk = self.states + random.cumsum(axis=1)
        samples = np.clip(walk, min=self.v_min, max=self.v_max)
        self.states = samples[count - 1]

        return samples

    def generate_iv_pairs(self, count):
        pattern = self.generate_random_pattern(count)
        result = [(
            self.v_min + (self.v_max - self.v_min) * pattern[::, i],
            np.full(count, self.c)
        ) for i in range(self.node_count)]
        return result

    def generate_iv_surfaces(self, count, ramp_width):
        assert count % ramp_width == 0
        ramp = np.arange(self.v_min, self.v_max, (self.v_max - self.v_min) / ramp_width)
        vs = np.tile(ramp, count // ramp_width)

        intensities = self.generate_random_pattern(count)

        result = [(
            vs,
            self.pv_model.get_cs(vs, intensities[i])
        ) for i in range(self.node_count)]

        return result

# -----
# TODO: temporary tests
# -----
import matplotlib.pyplot as plt

pv = PVModel()

for i_pv in [0, 1e-1]:
    vs = np.arange(10, 0, -5 / 100)
    cs = [pv.get_c(v_hrv=v, i_pv=i_pv) for v in vs]

    plt.plot(vs, cs)

plt.xlim(10, 0)
plt.show()

# generator = MultivarRndWalk(
#     node_count=2,
#     seed=3,
#     correlation=0.9,
#     variance=1e-3,
#     voltage_min=0.0,
#     voltage_max=5.0,
#     current=10e-3
# )
#
# count = 100
# stride = 1
# samples = generator.generate_iv_surfaces(count, 10)
#
# ts = np.arange(start=0, stop=count*10e-6, step=10e-6)
# (vs, cs) = samples[0]
#
# # plt.plot(ts[::stride], vs[::stride])
# plt.plot(ts[::stride], cs[::stride])
# plt.ylim(0, 5)
# plt.show()
