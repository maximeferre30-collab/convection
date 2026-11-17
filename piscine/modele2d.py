"""Niveau 2 : convection naturelle 2D (coupe verticale), Boussinesq, racks en milieu poreux chauffant.

Formulation fonction de courant / vorticité :
    ∇²ψ = -ω,  u = ∂ψ/∂z,  w = -∂ψ/∂x
    ∂ω/∂t + u·∇ω = ν_eff ∇²ω + gβ ∂T/∂x - rot(σ u)          (σ : frein de Darcy dans les racks)
    ∂T/∂t + ∇·(u T) = α_eff ∇²T + q'''/(ρc)
ψ aux nœuds, T aux centres de mailles : les vitesses de face dérivées de ψ sont à divergence nulle
exactement, donc le bilan d'énergie du schéma volumes finis sur T est exact (vérifié par les tests).

shortcut: turbulence = viscosité effective constante (ν_eff = α_eff), pas de k-ε ; à remplacer par
Code_Saturne (k-ω SST) pour des valeurs quantitatives de vitesse et de points chauds.
"""
from dataclasses import dataclass

import numpy as np
from scipy.fft import dstn, idstn

from .modele0d import Piscine
from .physique import CP, G, P_ATM, T_sat, flux_pertes, rho_eau


@dataclass
class Config2D:
    nx: int = 48
    nz: int = 50
    rack_x: tuple = (1.0, 11.0)   # m, étendue horizontale des racks dans la coupe
    rack_z: tuple = (0.3, 4.5)    # m, du pied de rack au sommet des assemblages
    repartition: str = "uniforme" # "uniforme", "centre" ou "bord" : où sont les assemblages chauds
    facteur_chaud: float = 3.0    # densité de puissance du tiers chaud / tiers froids
    nu_eff: float = 5e-3          # m²/s, viscosité turbulente effective (ordre de grandeur)
    sigma: float = 0.2            # 1/s, frein de Darcy linéaire dans les racks
    beta: float = 5e-4            # 1/K, dilatation moyenne 30–100 °C


class Simulation2D:
    def __init__(self, p: Piscine, c: Config2D = Config2D()):
        self.p, self.c = p, c
        self.Lx, self.H = p.longueur, p.niveau0
        self.dx, self.dz = self.Lx / c.nx, self.H / c.nz
        xc = (np.arange(c.nx) + 0.5) * self.dx
        zc = (np.arange(c.nz) + 0.5) * self.dz
        self.xc, self.zc = xc, zc
        self.rhoc = rho_eau(p.T_init) * CP

        # Racks : masque, et poids de puissance selon la répartition des assemblages chauds.
        in_x = (xc >= c.rack_x[0]) & (xc <= c.rack_x[1])
        in_z = (zc >= c.rack_z[0]) & (zc <= c.rack_z[1])
        rack = in_x[:, None] & in_z[None, :]
        tiers = np.clip(((xc - c.rack_x[0]) / (c.rack_x[1] - c.rack_x[0]) * 3).astype(int), 0, 2)
        chaud = {"uniforme": None, "centre": 1, "bord": 0}[c.repartition]
        poids_x = np.where(tiers == chaud, c.facteur_chaud, 1.0) if chaud is not None else np.ones(c.nx)
        w = rack * poids_x[:, None]
        # q'''·V_cellule·profondeur somme à P : w_norm en 1/m³
        self.w_norm = w / (w.sum() * self.dx * self.dz * p.largeur)
        # σ aux nœuds intérieurs : moyenne des 4 mailles voisines
        r = rack.astype(float)
        sig = np.zeros((c.nx + 1, c.nz + 1))
        sig[1:-1, 1:-1] = c.sigma * 0.25 * (r[1:, 1:] + r[:-1, 1:] + r[1:, :-1] + r[:-1, :-1])
        self.sig = sig

        # Valeurs propres du laplacien Dirichlet (DST-I) sur les nœuds intérieurs
        kx = np.arange(1, c.nx)
        kz = np.arange(1, c.nz)
        lx = (2 * np.cos(np.pi * kx / c.nx) - 2) / self.dx ** 2
        lz = (2 * np.cos(np.pi * kz / c.nz) - 2) / self.dz ** 2
        self.lap_eig = lx[:, None] + lz[None, :]

        self.t = 0.0
        self.T = np.full((c.nx, c.nz), p.T_init)
        self.om = np.zeros((c.nx + 1, c.nz + 1))
        self.psi = np.zeros((c.nx + 1, c.nz + 1))
        # Saturation locale : la pression hydrostatique relève T_sat en profondeur (ancrée sur T_ebul en surface)
        dTsat = T_sat(P_ATM + rho_eau(p.T_init) * G * (self.H - zc)) - T_sat(P_ATM)
        self.Tsat_loc = (p.T_ebul + dTsat)[None, :] * np.ones((c.nx, 1))

    # --- briques -------------------------------------------------------------------------------
    def vitesses_faces(self):
        u = (self.psi[:, 1:] - self.psi[:, :-1]) / self.dz      # (nx+1, nz)
        w = -(self.psi[1:, :] - self.psi[:-1, :]) / self.dx     # (nx, nz+1)
        return u, w

    def flux_surface(self):
        """Flux perdu par la surface (W/m²) pour chaque colonne, évalué sur la maille du haut."""
        if not self.p.pertes:
            return np.zeros(self.c.nx)
        return flux_pertes(self.T[:, -1], self.p.T_air, self.p.HR)[0]

    def pas(self, dt_max):
        """Avance d'un pas stable, au plus dt_max (s)."""
        c, dx, dz, nu = self.c, self.dx, self.dz, self.c.nu_eff
        T, om, psi = self.T, self.om, self.psi

        # Température : volumes finis, upwind + diffusion, flux nul sauf en surface.
        u, w = self.vitesses_faces()
        dt = min(dt_max, self.dt_stable(u, w))
        Fx = np.zeros_like(u)
        ui = u[1:-1]
        Fx[1:-1] = np.where(ui > 0, ui * T[:-1], ui * T[1:]) - nu * (T[1:] - T[:-1]) / dx
        Fz = np.zeros_like(w)
        wi = w[:, 1:-1]
        Fz[:, 1:-1] = np.where(wi > 0, wi * T[:, :-1], wi * T[:, 1:]) - nu * (T[:, 1:] - T[:, :-1]) / dz
        Fz[:, -1] = self.flux_surface() / self.rhoc
        q3 = self.p.puissance(self.t) * self.w_norm / self.rhoc
        dT = -(Fx[1:] - Fx[:-1]) / dx - (Fz[:, 1:] - Fz[:, :-1]) / dz + q3

        # Vorticité aux nœuds intérieurs : upwind + diffusion + flottabilité + frein des racks.
        U = np.zeros_like(psi)
        W = np.zeros_like(psi)
        U[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * dz)
        W[1:-1, 1:-1] = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2 * dx)
        o = om[1:-1, 1:-1]
        Ui, Wi = U[1:-1, 1:-1], W[1:-1, 1:-1]
        adv = (np.maximum(Ui, 0) * (o - om[:-2, 1:-1]) + np.minimum(Ui, 0) * (om[2:, 1:-1] - o)) / dx \
            + (np.maximum(Wi, 0) * (o - om[1:-1, :-2]) + np.minimum(Wi, 0) * (om[1:-1, 2:] - o)) / dz
        lap = (om[2:, 1:-1] - 2 * o + om[:-2, 1:-1]) / dx ** 2 + (om[1:-1, 2:] - 2 * o + om[1:-1, :-2]) / dz ** 2
        dTdx = ((T[1:, :-1] - T[:-1, :-1]) + (T[1:, 1:] - T[:-1, 1:])) / (2 * dx)
        SU, SW = self.sig * U, self.sig * W
        rot_frein = (SW[2:, 1:-1] - SW[:-2, 1:-1]) / (2 * dx) - (SU[1:-1, 2:] - SU[1:-1, :-2]) / (2 * dz)
        dom = -adv + nu * lap + G * c.beta * dTdx - rot_frein

        self.T = T + dt * dT
        om[1:-1, 1:-1] = o + dt * dom
        # Poisson ∇²ψ = -ω (DST-I, Dirichlet ψ = 0 sur le contour)
        psi[1:-1, 1:-1] = idstn(dstn(-om[1:-1, 1:-1], type=1) / self.lap_eig, type=1)
        # Vorticité de paroi (Thom) : non-glissement sur fond et côtés, glissement à la surface libre
        om[0, :] = -2 * psi[1, :] / dx ** 2
        om[-1, :] = -2 * psi[-2, :] / dx ** 2
        om[:, 0] = -2 * psi[:, 1] / dz ** 2
        om[:, -1] = 0.0
        self.t += dt

    def dt_stable(self, u, w, dt_max=5.0):
        # Schéma explicite monotone : advection + diffusion + frein cumulés ≤ 1 dans chaque maille
        # (des critères CFL et diffusion pris séparément laissent apparaître un damier instable).
        taux = (np.abs(u).max() / self.dx + np.abs(w).max() / self.dz
                + 2 * self.c.nu_eff * (1 / self.dx ** 2 + 1 / self.dz ** 2) + self.c.sigma)
        return min(dt_max, 0.8 / taux)

    def energie(self):
        """Énergie sensible de la piscine (J), la coupe étant extrudée sur toute la largeur."""
        return self.rhoc * self.T.sum() * self.dx * self.dz * self.p.largeur

    def avancer(self, t_fin, sortie=600.0, champs=None, cb=None):
        """Avance jusqu'à t_fin (s). Enregistre les séries tous les `sortie` s ; `champs` = instants des
        instantanés (T, ψ) à garder. Retourne (series, instantanés)."""
        series = {k: [] for k in ("t", "T_moy", "T_max", "T_surf", "v_max", "marge_sat")}
        snaps, champs = [], sorted(champs or [])
        prochain = self.t
        while self.t < t_fin:
            if self.t >= prochain:
                u, w = self.vitesses_faces()
                series["t"].append(self.t)
                series["T_moy"].append(self.T.mean())
                series["T_max"].append(self.T.max())
                series["T_surf"].append(self.T[:, -1].mean())
                series["v_max"].append(max(np.abs(u).max(), np.abs(w).max()))
                series["marge_sat"].append((self.Tsat_loc - self.T).min())
                prochain += sortie
            while champs and self.t >= champs[0]:
                snaps.append((self.t, self.T.copy(), self.psi.copy()))
                champs.pop(0)
            self.pas(t_fin - self.t)
        return {k: np.array(v) for k, v in series.items()}, snaps
