"""
Loads the stellar mass observational data.
"""

import numpy as np


class ObservationalData(object):
    """
    Holds observational data.
    """

    def __init__(
        self,
        scale_factor,
        stellar_density,
        density_error,
        scale_factor_error,
        description,
    ):
        """
        Store stuff in me!
        """
        self.scale_factor = scale_factor
        self.stellar_density = stellar_density
        self.density_error = density_error
        self.scale_factor_error = scale_factor_error
        self.description = description

        if density_error is None:
            self.old_eagle_result = True
        else:
            self.old_eagle_result = False


def read_obs_data(path="observational_data"):
    """
    Loads the observational data for the stellar mass example.
    """

    output_data = []

    h = 0.677

    # Moustakas et al 2013, Table 5
    obs7_z = np.array([0.1, 0.25, 0.35, 0.45])
    obs7_a = 1.0 / (obs7_z + 1.0)
    obs7_log10_rho_star = np.array([8.35, 8.32, 8.35, 8.31]) - np.log10(0.7 / h)
    obs7_log10_err_rho_star = np.array([0.05, 0.09, 0.06, 0.08])
    obs7_rho_star = 10 ** obs7_log10_rho_star
    obs7_rho_star_low = 10 ** (obs7_log10_rho_star - obs7_log10_err_rho_star)
    obs7_rho_star_high = 10 ** (obs7_log10_rho_star + obs7_log10_err_rho_star)

    output_data.append(
        ObservationalData(
            obs7_a,
            obs7_rho_star,
            [obs7_rho_star - obs7_rho_star_low, obs7_rho_star_high - obs7_rho_star],
            None,
            "Moustakas et al 2013, Table 5",
        )
    )

    # Muzzin et al 2013, Table 2
    obs8_z_low = np.array([0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    obs8_z_high = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
    obs8_z = 0.5 * (obs8_z_low + obs8_z_high)
    obs8_a_low = 1.0 / (obs8_z_low + 1.0)
    obs8_a_high = 1.0 / (obs8_z_high + 1.0)
    obs8_a = 1.0 / (obs8_z + 1.0)
    obs8_log10_rho_star = np.array([8.41, 8.26, 8.02, 7.79, 7.43, 7.32, 6.64])
    obs8_log10_rho_star_low = np.array([0.06, 0.03, 0.03, 0.03, 0.04, 0.09, 0.19])
    obs8_log10_rho_star_high = np.array([0.06, 0.03, 0.03, 0.05, 0.11, 0.13, 0.43])
    obs8_rho_star = 10 ** obs8_log10_rho_star
    obs8_rho_star_low = 10 ** (obs8_log10_rho_star - obs8_log10_rho_star_low)
    obs8_rho_star_high = 10 ** (obs8_log10_rho_star + obs8_log10_rho_star_high)

    output_data.append(
        ObservationalData(
            obs8_a,
            obs8_rho_star,
            [obs8_rho_star - obs8_rho_star_low, obs8_rho_star_high - obs8_rho_star],
            [obs8_a - obs8_a_low, obs8_a_high - obs8_a],
            "Muzzin et al 2013, Table 2",
        )
    )

    return output_data

