"""Tests de vérification (pytest, ou directement : python tests/test_modeles.py)."""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from piscine.modele0d import Piscine, simuler  # noqa: E402
from piscine.modele2d import Config2D, Simulation2D  # noqa: E402
from piscine.physique import CP, JOUR, L_VAP, flux_pertes, p_sat, rho_eau, way_wigner  # noqa: E402


def test_way_wigner_ordre_de_grandeur():
    # ~0,17 % de P0 un mois après l'arrêt, pour 3 ans d'irradiation
    r = way_wigner(30 * JOUR, 1.0, 3 * 365.25 * JOUR)
    assert 1.5e-3 < r < 1.8e-3


def test_proprietes_eau():
    assert abs(rho_eau(100) - 958.35) < 0.3 and abs(rho_eau(20) - 998.2) < 0.3
    assert abs(p_sat(100) - 101325) / 101325 < 0.01


def test_0d_puissance_constante_analytique():
    # Énoncé : 1 500 m³ à 30 °C et 2 MW constants -> ~60 h avant ébullition
    p = Piscine(longueur=10, largeur=10, niveau0=15, v_solide=0, P_init=2e6,
                t_refroid=1e12, T_irr=1e12, pertes=False)  # t_refroid énorme : P quasi constante
    r = simuler(p)
    attendu = p.masse(15, 30) * CP * 70 / 2e6
    assert abs(r["t_ebul"] - attendu) / attendu < 1e-4
    assert 55 * 3600 < r["t_ebul"] < 62 * 3600
    # puis vaporisation à P/L jusqu'au sommet des racks
    m_vap = p.masse(15, 30) - p.masse(p.h_racks, 100)  # masse conservée jusqu'à l'ébullition
    assert abs((r["t_decouvr"] - r["t_ebul"]) - m_vap * L_VAP / 2e6) / (m_vap * L_VAP / 2e6) < 1e-3


def test_0d_pertes_retardent_ebullition():
    p = Piscine()
    assert simuler(p)["t_ebul"] > simuler(p.avec(pertes=False))["t_ebul"]


def test_flux_pertes_valeurs_de_reference():
    # À 100 °C, air 30 °C / HR 50 % : évaporation (101,3 - 2,1) kPa × 0,089 × 1000 ≈ 8,8 kW/m²
    # + convection 1,52·70^(4/3) ≈ 0,44 kW/m² + rayonnement ≈ 0,64 kW/m²
    q, mdot = flux_pertes(100.0)
    assert 9.5e3 < q < 10.2e3
    assert abs(mdot * L_VAP - 8.83e3) < 0.1e3
    assert flux_pertes(30.0)[0] < 300  # eau à la température de l'air : seule l'évaporation reste


def test_0d_echec_integration_signale():
    # P infinie à t = 0 : l'échec doit lever une erreur, pas renvoyer « pas d'ébullition »
    try:
        simuler(Piscine(t_refroid=0.0))
    except RuntimeError:
        return
    raise AssertionError("échec silencieux")


def test_2d_pas_de_temps_monotone():
    # Critère de positivité vérifié maille par maille (le damier de la v1 venait de sa violation)
    p = Piscine(v_solide=0)
    s = Simulation2D(p, Config2D(repartition="bord"))
    s.avancer(3 * 3600)
    u, w = s.vitesses_faces()
    uc, wc = 0.5 * np.abs(u[1:] + u[:-1]), 0.5 * np.abs(w[:, 1:] + w[:, :-1])
    nu = s.c.nu_eff
    somme = s.dt_stable(u, w) * (uc / s.dx + wc / s.dz + 2 * nu * (1 / s.dx ** 2 + 1 / s.dz ** 2))
    assert somme.max() <= 1.0
    assert np.isfinite(s.T).all() and s.T.max() - s.T.mean() < 2.0


def test_2d_conservation_energie_et_accord_0d():
    p = Piscine(v_solide=0, pertes=False)
    s = Simulation2D(p, Config2D(nx=24, nz=25, repartition="bord"))
    E0 = s.energie()
    s.avancer(2 * 3600)
    Ein = quad(p.puissance, 0, s.t)[0]
    assert abs((s.energie() - E0) / Ein - 1) < 1e-5
    r = simuler(p)
    assert abs(s.T.mean() - np.interp(s.t, r["t"], r["T"])) < 0.01
    assert s.T.max() > s.T.mean() + 0.1  # la répartition non uniforme crée un point chaud


if __name__ == "__main__":
    for nom, f in list(globals().items()):
        if nom.startswith("test_"):
            f()
            print("ok", nom)
