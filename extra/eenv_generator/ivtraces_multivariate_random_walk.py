from common import generate_h5_files, EEnvGenerator, STEP_WIDTH
from pathlib import Path
import numpy as np
import math

# TODO evaluate switching from traces to surfaces (update filename and generator description accordingly)
# TODO isn't variance irrelevant if we just check rnd() > 0?
class MultivarRndWalk(EEnvGenerator):
    '''
    Generates traces according to a multivariate random walk. Each step, the voltage
    increases or decreases randomly with to the specified step_size.
    '''

    def __init__(self, node_count, seed, voltage_step_size, voltage_min, voltage_max, max_current):
        super().__init__(node_count, seed)

        self.v_step = voltage_step_size
        self.v_min = voltage_min
        self.v_max = voltage_max
        self.c = max_current

        # Start pattern between the two bounds
        self.state = (voltage_max + voltage_min) / 2.0

    def generate_covariance_matrix(dim, variance, correlation):
        # Create a correlation matrix with off-diagonal values set to the correlation coefficient
        correlation_matrix = np.full((dim, dim), correlation)
        np.fill_diagonal(correlation_matrix, 1)  # Diagonal should be 1 (self-correlation)

        # Convert correlation matrix to covariance matrix using variance
        covariance_matrix = correlation_matrix * variance

        return covariance_matrix

    def generate_random_pattern(self, count):
        samples = np.zeros((count, self.node_count))

        for i in range(1, n_steps):
            rands = multivariate_normal.rvs(None, cov_matrix)
            for j in range(node_count):
                if rands[j] < 0.0:
                    samples[i, j] = samples[i - 1, j] + step_size
                else:
                    samples[i, j] = samples[i - 1, j] - step_size
                samples[i, j] = max(min(samples[i, j], bound_high), bound_low)
        return samples

    def generate_iv_pairs(self, count):
        pattern = self.generate_random_pattern(count)
        result = [(self.on_voltage * pattern[::, i], self.on_current * pattern[::, i])
                    for i in range(self.node_count)]
        return result

if __name__ == "__main__":    
    path_here = Path(__file__).parent.absolute()
    if Path("/var/shepherd/").exists():
        path_eenv = Path("/var/shepherd/content/eenv/nes_lab/")
    else:
        path_eenv = path_here / "content/eenv/nes_lab/"

    seed = 32220789340897324098232347119065234157809
    duration = 1 * 60 * 60.0

    generator = RndPeriodicWindowGenerator(node_count=10, 
                                           seed=seed,
                                           period=10e-3,
                                           duty_cycle=0.2,
                                           on_voltage=1,
                                           on_current=100e-3)

    # Create folder
    name = f"random_window_test"
    folder_path = path_eenv / name
    folder_path.mkdir(parents=True, exist_ok=False)

    generate_h5_files(folder_path, duration=duration, chunk_size=500_000, generator=generator)


def generate_cov_matrix(dim: int, variance: float, correlation: float):
    """
    Generate a covariance matrix for a multivariate normal distribution with
    strong dependence and given variance.

    Parameters:
    dim (int): Number of dimensions.
    variance (float): Desired variance for each marginal distribution.
    correlation (float): Correlation between the variables (close to 1 for strong dependence).

    Returns:
    np.ndarray: A valid covariance matrix.
    """
    # Create a correlation matrix with off-diagonal values set to the correlation coefficient
    correlation_matrix = np.full((dim, dim), correlation)
    np.fill_diagonal(correlation_matrix, 1)  # Diagonal should be 1 (self-correlation)

    # Convert correlation matrix to covariance matrix using variance
    covariance_matrix = correlation_matrix * variance

    return covariance_matrix

def multivariate_random_walk(node_count, step_size, n_steps, correlation, variance, bound_low, bound_high):
    # Generate covariance matrix
    cov_matrix = generate_cov_matrix(dim=node_count, variance=variance, correlation=correlation)

    samples = np.zeros((n_steps, node_count))
    # Start in the middle between the bounds
    samples[0] = (bound_low + bound_high) / 2.0

    for i in range(1, n_steps):
        rands = multivariate_normal.rvs(None, cov_matrix)
        for j in range(node_count):
            if rands[j] < 0.0:
                samples[i, j] = samples[i - 1, j] + step_size
            else:
                samples[i, j] = samples[i - 1, j] - step_size
            samples[i, j] = max(min(samples[i, j], bound_high), bound_low)
    return samples

# -----
# TODO: temporary tests
# -----
import matplotlib.pyplot as plt

samples = multivariate_random_walk(node_count=3,
                                   correlation=0.9,
                                   variance=0.05,
                                   bound_low=0.0,
                                   bound_high=5.0,
                                   step_size=0.01,
                                   n_steps=100_000)

plt.plot(samples)
plt.show()
