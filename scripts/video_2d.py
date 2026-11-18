"""Vidéo 2D : coupe verticale de la piscine de 30 °C jusqu'à l'ébullition (~110 h), assemblages chauds au bord.

Usage : python scripts/video_2d.py   (écrit figures/2d_video_ebullition.mp4, 1920 × 1080, et un aperçu .png)

Même lecture que la coupe 3D : la couleur est l'écart à la température moyenne de la piscine (échelle
fixe), car la piscine chauffe de 70 K alors que les écarts locaux restent inférieurs à 1 K.
"""
import sys
from pathlib import Path

import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter

RACINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RACINE))
sys.stdout.reconfigure(encoding="utf-8")
from piscine.modele0d import Piscine, simuler  # noqa: E402
from piscine.modele2d import Config2D, Simulation2D  # noqa: E402

FIG = RACINE / "figures"
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

ref = Piscine(v_solide=0.0)            # même masse d'eau qu'en 2D (pas de volume solide dans la coupe)
t_ebul = simuler(ref)["t_ebul"]        # ébullition 0D avec pertes, ~108,6 h
t_fin = t_ebul + 2 * 3600
PAS = 900.0                            # une image tous les 15 min de temps physique

sim = Simulation2D(ref, Config2D(repartition="bord"))
instants = list(np.arange(0, t_fin + 1, PAS))
ser, snaps = sim.avancer(t_fin, sortie=PAS, champs=instants)
print(f"{len(snaps)} images, T_moy finale {ser['T_moy'][-1]:.1f} °C, "
      f"surchauffe max {(ser['T_max'] - ser['T_moy'])[10:].max():.2f} K")

c = sim.c
ecarts = [T - T.mean() for _, T, _ in snaps[4:]]   # la première heure (mise en route) sature l'échelle
vmax = np.percentile([e.max() for e in ecarts], 99)
vmin = -0.4 * vmax
xn, zn = np.linspace(0, sim.Lx, c.nx + 1), np.linspace(0, sim.H, c.nz + 1)
h = ser["t"] / 3600

fig = plt.figure(figsize=(16, 9), dpi=120)                 # 1920 × 1080
ax = fig.add_axes([0.04, 0.08, 0.50, 0.84])
cax = fig.add_axes([0.555, 0.15, 0.012, 0.70])
at = fig.add_axes([0.66, 0.56, 0.32, 0.34])
am = fig.add_axes([0.66, 0.10, 0.32, 0.34])

im = ax.imshow(np.zeros((c.nz, c.nx)), origin="lower", extent=(0, sim.Lx, 0, sim.H), cmap="coolwarm",
               vmin=vmin, vmax=vmax, interpolation="bicubic", aspect="equal")
fig.colorbar(im, cax=cax, label="T − T moyenne (K)")
ax.add_patch(plt.Rectangle((c.rack_x[0], c.rack_z[0]), c.rack_x[1] - c.rack_x[0], c.rack_z[1] - c.rack_z[0],
                           fill=False, ec="k", ls="--", lw=1))
tiers = (c.rack_x[1] - c.rack_x[0]) / 3
ax.add_patch(plt.Rectangle((c.rack_x[0], c.rack_z[0]), tiers, c.rack_z[1] - c.rack_z[0],
                           fill=False, ec="k", hatch="//", lw=0, alpha=0.35))
ax.text(c.rack_x[0] + tiers / 2, c.rack_z[1] + 0.25, "assemblages\nchauds", ha="center", fontsize=10)
ax.set(xlabel="x (m)", ylabel="z (m)", xlim=(0, sim.Lx), ylim=(0, sim.H))
titre = ax.set_title("", fontsize=15)
lignes = []

at.plot(h, ser["T_moy"], c="C2", label="T moyenne (2D)")
at.plot(h, ser["T_max"], c="C3", lw=0.8, label="T max (2D)")
at.axhline(100, c="gray", lw=0.8)
at.axvline(t_ebul / 3600, c="k", ls=":", lw=1, label=f"ébullition 0D ({t_ebul/3600:.0f} h)")
pt, = at.plot([], [], "o", c="k")
at.set(ylabel="T (°C)", xlim=(0, h[-1]), title="Montée en température")
at.legend(fontsize=9, loc="lower right")
am.plot(h, ser["T_max"] - ser["T_moy"], c="C3", lw=0.8, label="point chaud − moyenne (K)")
am.plot(h, ser["marge_sat"] / 20, c="C0", label="marge à l'ébullition locale / 20 (K)")
am.axhline(0, c="gray", lw=0.8)
pm = am.axvline(0, c="k", lw=1)
am.set(xlabel="temps (h)", ylabel="K", xlim=(0, h[-1]), ylim=(-0.3, 2.5), title="Point chaud et marge à l'ébullition")
am.legend(fontsize=9, loc="upper right")
for a in (at, am):
    a.grid(alpha=0.3)

ecrivain = FFMpegWriter(fps=24, codec="libx264", extra_args=["-crf", "23", "-pix_fmt", "yuv420p"])
with ecrivain.saving(fig, FIG / "2d_video_ebullition.mp4", dpi=120):
    for k, (t, T, psi) in enumerate(snaps):
        im.set_data((T - T.mean()).T)
        for l in lignes:
            l.remove()
        cs = ax.contour(xn, zn, psi.T, levels=14, colors="k", linewidths=0.5, alpha=0.6)  # lignes de courant
        lignes = [cs]
        titre.set_text(f"Coupe 2D, assemblages chauds au bord : t = {t/3600:5.1f} h, "
                       f"T moyenne = {T.mean():5.1f} °C")
        pt.set_data([t / 3600], [T.mean()])
        pm.set_xdata([t / 3600])
        ecrivain.grab_frame()
        if k == int(0.6 * len(snaps)):
            fig.savefig(FIG / "2d_video_ebullition_apercu.png", dpi=80)
print("Vidéo :", FIG / "2d_video_ebullition.mp4")
