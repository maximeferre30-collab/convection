"""Briques physiques communes : puissance résiduelle, propriétés de l'eau, pertes en surface."""
import numpy as np

G = 9.81            # m/s²
CP = 4186.0         # J/(kg·K), eau liquide (varie de 4178 à 4216 entre 30 et 100 °C)
L_VAP = 2.257e6     # J/kg, chaleur latente à 100 °C
SIGMA = 5.670e-8    # W/(m²·K⁴)
P_ATM = 101325.0    # Pa
JOUR = 86400.0
AN = 365.25 * JOUR


def way_wigner(t, P0, T0):
    """Puissance résiduelle de Way-Wigner (W).

    t : temps depuis l'arrêt (s), T0 : durée d'irradiation (s), P0 : puissance en fonctionnement (W).
    Ordre de grandeur seulement : l'écart à ANS-5.1 croît au-delà de quelques mois (actinides, Cs/Sr).
    """
    t = np.asarray(t, dtype=float)
    return 0.0622 * P0 * (t ** -0.2 - (t + T0) ** -0.2)


def rho_eau(T):
    """Masse volumique de l'eau (kg/m³), ajustement quadratique 20–100 °C (erreur < 0,2 kg/m³)."""
    return 1002.305 - 0.13347 * T - 3.0766e-3 * T ** 2


def p_sat(T):
    """Pression de vapeur saturante (Pa), équation d'Antoine valable 1–100 °C."""
    return 133.322 * 10 ** (8.07131 - 1730.63 / (233.426 + T))


def T_sat(p):
    """Température de saturation (°C) à la pression p (Pa), Antoine 99–374 °C."""
    return 1810.94 / (8.14019 - np.log10(p / 133.322)) - 244.485


def flux_pertes(T, T_air=30.0, HR=0.5, V_air=0.0, eps=0.95):
    """Flux surfacique perdu par la surface libre (W/m²) et débit évaporé (kg/(s·m²)).

    - évaporation : corrélation ASHRAE/Carrier q = (p_w - p_a)[kPa]·(0,089 + 0,0782·V)·1000
      (établie pour des piscines à 25–30 °C : son extrapolation vers 80–100 °C est une incertitude majeure) ;
    - convection naturelle air, plaque chaude face vers le haut, turbulent : h = 1,52·ΔT^(1/3) ;
    - rayonnement vers les parois du bâtiment supposées à T_air.
    """
    T = np.asarray(T, dtype=float)
    dp_kPa = np.maximum(p_sat(T) - HR * p_sat(T_air), 0.0) / 1000.0
    q_evap = dp_kPa * (0.089 + 0.0782 * V_air) * 1000.0
    dT = np.maximum(T - T_air, 0.0)
    q_conv = 1.52 * dT ** (4.0 / 3.0)
    q_rad = eps * SIGMA * ((T + 273.15) ** 4 - (T_air + 273.15) ** 4)
    return q_evap + q_conv + q_rad, q_evap / L_VAP
