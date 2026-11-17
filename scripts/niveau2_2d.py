"""Niveau 2 : convection naturelle 2D – vérification vs 0D, points chauds, maillage, animation.

Usage : python scripts/niveau2_2d.py   (~20 min ; écrit figures/ et resultats/)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

RACINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RACINE))
sys.stdout.reconfigure(encoding="utf-8")
from piscine.modele0d import Piscine, simuler  # noqa: E402
from piscine.modele2d import Config2D, Simulation2D  # noqa: E402

FIG, RES = RACINE / "figures", RACINE / "resultats"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

# Même masse d'eau en 0D et en 2D : la coupe 2D ne représente pas le volume solide des racks.
ref = Piscine(v_solide=0.0)


def champ(ax, sim, T, psi, titre, vmin=None, vmax=None):
    c = sim.c
    im = ax.imshow(T.T, origin="lower", extent=(0, sim.Lx, 0, sim.H), cmap="inferno",
                   vmin=vmin, vmax=vmax, aspect="equal")
    xn, zn = np.linspace(0, sim.Lx, c.nx + 1), np.linspace(0, sim.H, c.nz + 1)
    ax.contour(xn, zn, psi.T, levels=14, colors="w", linewidths=0.6)  # isolignes de ψ = lignes de courant
    ax.add_patch(plt.Rectangle((c.rack_x[0], c.rack_z[0]), c.rack_x[1] - c.rack_x[0],
                               c.rack_z[1] - c.rack_z[0], fill=False, ec="c", ls="--", lw=1))
    ax.set(title=titre, xlabel="x (m)", ylabel="z (m)")
    return im


# --- 1. Vérification : bilan d'énergie 2D vs 0D, cas adiabatique ------------------------------
r0 = simuler(ref.avec(pertes=False))
s = Simulation2D(ref.avec(pertes=False))
ser_ad, _ = s.avancer(r0["t_ebul"] * 0.98, sortie=1800)
ecart = np.abs(ser_ad["T_moy"] - np.interp(ser_ad["t"], r0["t"], r0["T"])).max()
print(f"[vérif] adiabatique : écart max T_moy 2D / 0D = {ecart:.3f} K sur {ser_ad['t'][-1]/3600:.0f} h")

# --- 2. Répartition des assemblages chauds (avec pertes) --------------------------------------
r0p = simuler(ref)
t_fin = r0p["t_ebul"]
lignes = ["repartition,T_moy_fin,T_max_fin,dT_point_chaud_max,T_surf_moins_T_moy,v_max_cm_s,marge_saturation_min"]
res = {}
for rep in ("uniforme", "centre", "bord"):
    sim = Simulation2D(ref, Config2D(repartition=rep))
    ser, snaps = sim.avancer(t_fin, sortie=1800, champs=[24 * 3600])
    res[rep] = (sim, ser, snaps)
    dT = (ser["T_max"] - ser["T_moy"]).max()
    dS = (ser["T_surf"] - ser["T_moy"])[-1]
    lignes.append(f"{rep},{ser['T_moy'][-1]:.2f},{ser['T_max'][-1]:.2f},{dT:.2f},{dS:.2f},"
                  f"{100*ser['v_max'][-1]:.2f},{ser['marge_sat'].min():.2f}")
    print(f"[{rep:>8}] T_moy fin {ser['T_moy'][-1]:6.2f} °C | T_max {ser['T_max'][-1]:6.2f} °C | "
          f"surchauffe max {dT:4.2f} K | T_surf - T_moy {dS:+.2f} K | v_max {100*ser['v_max'][-1]:.1f} cm/s | "
          f"marge à saturation locale min {ser['marge_sat'].min():5.1f} K")
(RES / "repartition_2d.csv").write_text("\n".join(lignes) + "\n", encoding="utf-8")

fig, a = plt.subplots(figsize=(7, 4))
a.plot(r0p["t"] / 3600, r0p["T"], "k-", lw=2, label="0D")
a.plot(r0["t"] / 3600, r0["T"], "k:", lw=1, label="0D adiabatique")
a.plot(ser_ad["t"] / 3600, ser_ad["T_moy"], "o", ms=3, mfc="none", c="gray", label="2D adiabatique (moyenne)")
for rep, (_, ser, _) in res.items():
    l, = a.plot(ser["t"] / 3600, ser["T_moy"], "--", label=f"2D {rep} : moyenne")
    a.plot(ser["t"] / 3600, ser["T_max"], "-", lw=0.8, c=l.get_color(), label=f"2D {rep} : max")
a.set(xlabel="temps (h)", ylabel="T (°C)", xlim=(0, t_fin / 3600 * 1.02), title="Validation croisée 0D / 2D")
a.grid(alpha=0.3)
a.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(FIG / "2d_vs_0d.png", dpi=150)

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for a, (rep, (sim, _, snaps)) in zip(axs, res.items()):
    t, T, psi = snaps[0]
    im = champ(a, sim, T, psi, f"{rep} – t = {t/3600:.0f} h", vmin=T.min(), vmax=T.max())
    fig.colorbar(im, ax=a, shrink=0.8, label="T (°C)")
fig.suptitle("Champ de température et lignes de courant (racks en pointillés, tiers chaud ×3)")
fig.tight_layout()
fig.savefig(FIG / "2d_champs_repartition.png", dpi=150)

# --- 3. Convergence en maillage et sensibilité à ν_eff (6 h, répartition 'bord') -------------
# 6 h suffisent : l'écart point chaud / moyenne est établi en ~2 h ; le maillage 96 coûte 16× le 48.
t_c = 6 * 3600
lignes = ["etude,valeur,dT_point_chaud,v_max_cm_s"]
print("\n[maillage / turbulence] répartition 'bord', t = 6 h")
for etude, cfgs in [("maillage_nx", [Config2D(repartition="bord", nx=n, nz=n + n // 24) for n in (24, 48, 96)]),
                    ("nu_eff", [Config2D(repartition="bord", nu_eff=v) for v in (2.5e-3, 5e-3, 1e-2)])]:
    for c in cfgs:
        sim = Simulation2D(ref, c)
        ser, _ = sim.avancer(t_c, sortie=t_c)
        val = c.nx if etude == "maillage_nx" else c.nu_eff
        u, w = sim.vitesses_faces()
        dT, v = sim.T.max() - sim.T.mean(), 100 * max(np.abs(u).max(), np.abs(w).max())
        print(f"  {etude:<12}{val:>8} : surchauffe {dT:5.2f} K, v_max {v:5.2f} cm/s")
        lignes.append(f"{etude},{val},{dT:.3f},{v:.3f}")
(RES / "convergence_2d.csv").write_text("\n".join(lignes) + "\n", encoding="utf-8")

# --- 4. Animation de la mise en convection (première heure) ----------------------------------
sim = Simulation2D(ref, Config2D(repartition="bord"))
instants = list(np.arange(0, 3600 + 1, 60.0))
_, snaps = sim.avancer(3601, sortie=600, champs=instants)
fig, a = plt.subplots(figsize=(6, 6))


def image(k):
    a.clear()
    t, T, psi = snaps[k]
    champ(a, sim, T, psi, f"Mise en convection – t = {t/60:.0f} min", vmin=30.0, vmax=snaps[-1][1].max())


anim = FuncAnimation(fig, image, frames=len(snaps))
anim.save(FIG / "2d_mise_en_convection.gif", writer=PillowWriter(fps=8), dpi=80)
print("\nFigures écrites dans", FIG)
