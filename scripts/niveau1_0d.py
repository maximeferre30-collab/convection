"""Niveau 1 : courbes du modèle 0D, étude de sensibilité, cas indicatif Fukushima Daiichi 4.

Usage : python scripts/niveau1_0d.py   (écrit figures/ et resultats/)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RACINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RACINE))
sys.stdout.reconfigure(encoding="utf-8")  # console Windows (cp1252)
from piscine.modele0d import Piscine, simuler  # noqa: E402
from piscine.physique import JOUR, flux_pertes, rho_eau  # noqa: E402

FIG, RES = RACINE / "figures", RACINE / "resultats"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)
h = lambda s: None if s is None else s / 3600
j = lambda s: None if s is None else s / JOUR
fmt = lambda v, d=1: "jamais" if v is None else f"{v:.{d}f}"

ref = Piscine()
print(f"Cas de référence : {ref.surface:.0f} m², {ref.masse(ref.niveau0, ref.T_init)/1e3:.0f} t d'eau, "
      f"P(0) = {ref.puissance(0)/1e6:.2f} MW")

# Nombre de Rayleigh (propriétés moléculaires de l'eau à ~45 °C)
nu, alpha, beta, dT = 6.0e-7, 1.5e-7, 4e-4, 10.0
print(f"Ra = gβΔTH³/(να) ≈ {9.81*beta*dT*ref.niveau0**3/(nu*alpha):.1e}  (> 1e9 : turbulent)")

# --- Figure 1 : température et niveau, adiabatique vs pertes ----------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
for p, lab in [(ref.avec(pertes=False), "adiabatique (conservatif)"), (ref, "avec pertes en surface")]:
    r = simuler(p)
    a1.plot(r["t"] / 3600, r["T"], label=f"{lab} : ébullition à {fmt(h(r['t_ebul']), 0)} h")
    a2.plot(r["t"] / JOUR, r["niveau"], label=f"{lab} : découvrement à {fmt(j(r['t_decouvr']))} j")
a1.set(xlabel="temps après perte du refroidissement (h)", ylabel="T moyenne (°C)", xlim=(0, 150))
a2.axhline(ref.h_racks, color="k", ls="--", lw=1)
a2.text(0.5, ref.h_racks + 0.2, "sommet des assemblages")
a2.set(xlabel="temps (jours)", ylabel="niveau d'eau (m)", ylim=(0, 13.5))
for a in (a1, a2):
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
fig.suptitle(f"Modèle 0D – P(0) = {ref.puissance(0)/1e6:.2f} MW, {ref.masse(ref.niveau0, 30)/1e3:.0f} t d'eau à 30 °C")
fig.tight_layout()
fig.savefig(FIG / "0d_temperature_niveau.png", dpi=150)

# --- Figure 2 : puissance résiduelle et pertes ------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
t = np.logspace(0, np.log10(3 * 365 * JOUR), 300)
a1.loglog(t / JOUR, ref.avec(t_refroid=0.0).puissance(t) / 1e6)
a1.axvline(ref.t_refroid / JOUR, color="r", ls="--", lw=1, label="perte du refroidissement (réf.)")
a1.set(xlabel="temps depuis l'arrêt (jours)", ylabel="P résiduelle (MW)", title="Way-Wigner, P0 = 1 GWth, 3 ans")
a1.legend(fontsize=8)
T = np.linspace(30, 100, 100)
q, _ = flux_pertes(T)
a2.plot(T, q * ref.surface / 1e6, label="pertes totales en surface")
a2.axhline(ref.puissance(0) / 1e6, color="r", ls="--", lw=1, label="P(0) cas de référence")
a2.set(xlabel="T eau (°C)", ylabel="MW", title=f"Pertes surface ({ref.surface:.0f} m², air 30 °C, HR 50 %)")
a2.legend(fontsize=8)
for a in (a1, a2):
    a.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(FIG / "0d_puissance_pertes.png", dpi=150)

# --- Sensibilité ------------------------------------------------------------------------------
lignes = ["parametre,valeur,P0_MW,t_ebul_h,t_decouvr_j"]
cas = [("age_combustible_j", v, ref.avec(t_refroid=v * JOUR)) for v in (7, 30, 100, 365, 1095)]
cas += [("niveau_initial_m", v, ref.avec(niveau0=v)) for v in (9.5, 11.0, 12.5)]
cas += [("T_initiale_C", v, ref.avec(T_init=v)) for v in (25, 30, 40, 50)]
cas += [("pertes_surface", v, ref.avec(pertes=v)) for v in (False, True)]
print("\nSensibilité (une variable à la fois, référence : 30 j, 12,5 m, 30 °C, pertes) :")
print(f"{'paramètre':<20}{'valeur':>8}{'P0 (MW)':>10}{'ébullition (h)':>16}{'découvrement (j)':>18}")
age = []
for nom, v, p in cas:
    r = simuler(p, t_max=400 * JOUR)
    P0 = p.puissance(0) / 1e6
    print(f"{nom:<20}{str(v):>8}{P0:>10.2f}{fmt(h(r['t_ebul'])):>16}{fmt(j(r['t_decouvr'])):>18}")
    lignes.append(f"{nom},{v},{P0:.3f},{fmt(h(r['t_ebul']))},{fmt(j(r['t_decouvr']))}")
    if nom.startswith("age"):
        age.append((v, h(r["t_ebul"]), j(r["t_decouvr"])))
(RES / "sensibilite_0d.csv").write_text("\n".join(lignes) + "\n", encoding="utf-8")

fig, a = plt.subplots(figsize=(6, 4))
a.semilogx([x[0] for x in age], [x[1] or np.nan for x in age], "o-", label="ébullition (h)")
a.semilogx([x[0] for x in age], [(x[2] or np.nan) * 24 for x in age], "s-", label="découvrement (h)")
a.set(xlabel="âge du combustible à la perte de refroidissement (jours)", ylabel="délai (h)", yscale="log")
a.grid(alpha=0.3, which="both")
a.legend()
fig.tight_layout()
fig.savefig(FIG / "0d_sensibilite_age.png", dpi=150)

# --- Cas indicatif Fukushima Daiichi 4 --------------------------------------------------------
# Valeurs couramment citées (À VÉRIFIER et sourcer dans les rapports TEPCO / NRC / AIEA avant usage) :
# piscine 12,2 × 9,9 m, ~11,5 m d'eau, P ≈ 2,26 MW le 11/03/2011, cœur déchargé ~100 j avant,
# 84 °C mesurés le 14/03 vers 4 h (≈ 61 h après le séisme).
f4 = Piscine(longueur=12.2, largeur=9.9, niveau0=11.5, P_init=2.26e6, t_refroid=100 * JOUR, T_init=27.0)
print("\nCas indicatif Fukushima Daiichi 4 (données à vérifier) :")
for p, lab in [(f4.avec(pertes=False), "adiabatique"), (f4, "avec pertes")]:
    r = simuler(p)
    T61 = np.interp(61 * 3600, r["t"], r["T"])
    print(f"  {lab:<12} T(61 h) = {T61:5.1f} °C (mesure citée : 84 °C), ébullition {fmt(h(r['t_ebul']))} h, "
          f"découvrement {fmt(j(r['t_decouvr']))} j")
