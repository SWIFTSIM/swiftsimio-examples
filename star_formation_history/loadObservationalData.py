"""
Code from Matthieu Schaller to load observational SFR data.
"""

from numpy import *

h = 0.68


class ObservationalData(object):
    """
    Holds observational data.
    """

    def __init__(self, scale_factor, sfr, error, description):
        """
        Store stuff in me!
        """
        self.scale_factor = scale_factor
        self.sfr = sfr
        self.error = error
        self.description = description

        if self.error is None:
            self.fitting_formula = True
        else:
            self.fitting_formula = False


def read_obs_data(path="observational_data"):
    """
    Reads the observational data
    """

    output = []

    # SFR Observational data from Hopkins 2004
    hcorr = log10(h) - log10(0.7)  # h^-2 for SFR, h^-3 for volume
    (z, z_up, z_down, lgrho, lgrho_up, lgrho_down) = loadtxt(
        f"{path}/sfr_hopkins2004_cor.dat", unpack=True
    )
    lgrho = lgrho + log10(0.6) + hcorr
    obs1_a = 1.0 / (1.0 + z)
    obs1_rho = 10 ** lgrho
    obs1_rho_err = array([obs1_rho - 10 ** lgrho_down, 10 ** lgrho_up - obs1_rho])

    # output.append(
    #    ObservationalData(
    #        scale_factor=obs1_a,
    #        sfr=obs1_rho,
    #        error=obs1_rho_err,
    #        description="Hopkins (2004)",
    #    )
    # )

    # SFR Observational data from Karim 2011
    (z, rho, err_up, err_down) = loadtxt(f"{path}/sfr_karim2011.dat", unpack=True)
    obs2_a = 1.0 / (1.0 + z)
    obs2_rho = rho * 0.6777 / 0.7
    obs2_rho_err = array([-err_down, err_up])

    output.append(
        ObservationalData(obs2_a, obs2_rho, obs2_rho_err, "Karim et al. (2011) [radio]")
    )

    # SFR Observational data from Bouwens 2012
    z, rhostar = loadtxt(f"{path}/bouwens_2012_sfrd_no_dust.txt", unpack=True)
    z, rhostar_dust = loadtxt(f"{path}/bouwens_2012_sfrd_dustcorr.txt", unpack=True)
    rhostar = (
        (rhostar / 1.8) * 0.6777 / 0.7
    )  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    rhostar_dust = (
        (rhostar_dust / 1.8) * 0.6777 / 0.7
    )  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    obs3_a = 1.0 / (1.0 + z)
    obs3_rho = rhostar
    obs3_rho_dust = rhostar_dust

    output.append(
        ObservationalData(
            obs3_a,
            obs3_rho,
            zeros_like(obs3_rho),
            "Bouwens et al. (2012) [UV, no dust]",
        )
    )

    # SFR Observational data from Rodighierio 2012
    z, rhostar, err_m, err_p = loadtxt(f"{path}/sfr_rodighiero2010.dat", unpack=True)
    rhostar = (
        (rhostar / 1.65) * 0.6777 / 0.7
    )  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    obs4_a = 1.0 / (1.0 + z)
    obs4_rho = rhostar
    obs4_rho_err = array([-err_m / 1.65, err_p / 1.65])

    output.append(
        ObservationalData(
            obs4_a, obs4_rho, obs4_rho_err, "Rodighiero et al. (2010) [24 $\mu$m]"
        )
    )

    # SFR Observational data from Cucciati 2012
    z, rhostar, err_m, err_p = loadtxt(f"{path}/sfr_cucciati2011.dat", unpack=True)
    rhostar = rhostar - 2.0 * log10(0.7) + 2.0 * log10(0.6777)
    obs5_a = 1.0 / (1 + z)
    obs5_rho = 10 ** rhostar / 1.65
    obs5_rho_err = 10 ** array([rhostar + (err_m), rhostar + err_p]) / 1.65
    obs5_rho_err[0] = -obs5_rho_err[0] + obs5_rho
    obs5_rho_err[1] = -obs5_rho + obs5_rho_err[1]

    output.append(
        ObservationalData(
            obs5_a, obs5_rho, obs5_rho_err, "Cucciati et al. (2012) [FUV]"
        )
    )

    ##################################################################################

    # SFR from Madau & Dickinson fitting formula (z < 10 and assuming h=0.7)
    obs6_a = logspace(log10(1.0 / (1.0 + 10.0)), 0, 100)
    obs6_z = 1.0 / obs6_a - 1.0
    obs6_rho = (
        0.015 * ((1.0 + obs6_z) ** 2.7) / (1.0 + ((1.0 + obs6_z) / 2.9) ** 5.6)
    )  # Msun / yr / Mpc^3
    obs6_rho /= 1.65  # Salpeter -> Chabrier correction

    output.append(
        ObservationalData(obs6_a, obs6_rho, None, "Madau & Dickinson (2014) [$h=0.7$]")
    )

    # Load the EAGLE NOAGN data
    eagle_data = loadtxt(f"{path}/EAGLE_NOAGN_sfr.txt")

    output.append(
        ObservationalData(
            eagle_data[:, 0], eagle_data[:, 2] / (25 ** 3), None, "EAGLE NoAGN"
        )
    )

    return output
