"""Niveau 1 : modèle 0D (piscine à température uniforme) après perte de refroidissement."""
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import solve_ivp

from .physique import AN, CP, JOUR, L_VAP, flux_pertes, rho_eau, way_wigner


@dataclass
class Piscine:
    longueur: float = 12.0        # m
    largeur: float = 10.0         # m
    niveau0: float = 12.5         # m, hauteur d'eau initiale
    h_racks: float = 4.5          # m, sommet des assemblages (critère de découvrement)
    v_solide: float = 70.0        # m³, volume occupé par racks + assemblages (sous h_racks)
    T_init: float = 30.0          # °C
    T_air: float = 30.0           # °C, air du hall supposé constant
    HR: float = 0.5               # humidité relative de l'air
    T_ebul: float = 100.0         # °C, saturation en surface (1 atm)
    P0: float = 1000e6            # W, puissance du combustible entreposé en fonctionnement (~1/3 de cœur 3 GWth)
    T_irr: float = 3 * AN         # s, durée d'irradiation
    t_refroid: float = 30 * JOUR  # s, temps écoulé entre l'arrêt et la perte de refroidissement
    P_init: float | None = None   # W, impose P(0) ; Way-Wigner ne donne alors que la forme de la décroissance
    pertes: bool = True           # False = adiabatique (conservatif)

    @property
    def surface(self):
        return self.longueur * self.largeur

    def puissance(self, t):
        P = way_wigner(self.t_refroid + np.asarray(t, float), self.P0, self.T_irr)
        if self.P_init is not None:
            P = P * self.P_init / way_wigner(self.t_refroid, self.P0, self.T_irr)
        return P

    def masse(self, niveau, T):
        return rho_eau(T) * (self.surface * niveau - self.v_solide)

    def niveau(self, m, T):
        return (m / rho_eau(T) + self.v_solide) / self.surface

    def avec(self, **kw):
        return replace(self, **kw)


def simuler(p: Piscine, t_max=60 * JOUR, pas=600.0):
    """Intègre le bilan m·c·dT/dt = P(t) - Φ(T), puis l'ébullition dm/dt = -P/L.

    Retourne un dict : t (s), T (°C), niveau (m), P (W), t_ebul et t_decouvr (s, None si non atteint).
    """
    A = p.surface

    def rhs1(t, y):
        T, m = y
        q, mdot = flux_pertes(T, p.T_air, p.HR) if p.pertes else (0.0, 0.0)
        return [(p.puissance(t) - q * A) / (m * CP), -mdot * A]

    def ebullition(t, y):
        return y[0] - p.T_ebul
    ebullition.terminal, ebullition.direction = True, 1

    def decouvrement(t, y):
        return p.niveau(y[1], y[0]) - p.h_racks
    decouvrement.terminal, decouvrement.direction = True, -1

    y0 = [p.T_init, p.masse(p.niveau0, p.T_init)]
    s1 = solve_ivp(rhs1, (0, t_max), y0, events=[ebullition, decouvrement], max_step=pas, rtol=1e-8)
    if s1.status < 0:  # sinon un échec passerait pour « pas d'ébullition »
        raise RuntimeError(f"échec de l'intégration 0D : {s1.message}")
    t, T, m = s1.t, s1.y[0], s1.y[1]
    t_ebul = s1.t_events[0][0] if s1.t_events[0].size else None
    t_dec = s1.t_events[1][0] if s1.t_events[1].size else None

    if t_ebul is not None and t_dec is None:
        # Phase 2 : eau à saturation, toute la puissance vaporise (conservatif : pertes convectives/radiatives négligées).
        def rhs2(t, y):
            return [-p.puissance(t) / L_VAP]

        def dec2(t, y):
            return p.niveau(y[0], p.T_ebul) - p.h_racks
        dec2.terminal, dec2.direction = True, -1

        s2 = solve_ivp(rhs2, (t_ebul, t_max), [m[-1]], events=dec2, max_step=pas * 6, rtol=1e-8)
        if s2.status < 0:
            raise RuntimeError(f"échec de l'intégration 0D (ébullition) : {s2.message}")
        t = np.concatenate([t, s2.t[1:]])
        T = np.concatenate([T, np.full(s2.t.size - 1, p.T_ebul)])
        m = np.concatenate([m, s2.y[0, 1:]])
        t_dec = s2.t_events[0][0] if s2.t_events[0].size else None

    return dict(t=t, T=T, niveau=p.niveau(m, T), P=p.puissance(t), t_ebul=t_ebul, t_decouvr=t_dec)
