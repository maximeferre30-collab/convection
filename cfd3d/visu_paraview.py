"""Visualisation ParaView prête à l'emploi des résultats Code_Saturne 3D.

Dans ParaView : View > Python Shell > Run Script, choisir ce fichier.
En ligne de commande : pvpython visu_paraview.py [RESULTS.case] [--png image.png | --film dossier]

Affiche une coupe verticale à travers le coin des assemblages chauds, colorée par l'écart à la
température moyenne de la piscine (ΔT = T − T_moy, échelle fixe). La piscine chauffe de 30 à 100 °C,
mais les écarts locaux ne dépassent pas 1 K : colorer par T brute donne un cube uniforme qui passe du
bleu au rouge. Les flèches montrent la vitesse de l'eau.
"""
import sys

from paraview.simple import *

CAS = (r"\\wsl.localhost\Ubuntu\home\maxim\cs\etudes\piscine\CAS_coin_nx48_ebullition"
       r"\RESU\20261118-013109\postprocessing\RESULTS.case")
args = [a for a in sys.argv[1:] if a]
if args and not args[0].startswith("--"):
    CAS = args[0]

piscine = EnSightReader(registrationName="piscine", CaseFileName=CAS)
piscine.CellArrays = ["Temperature", "Velocity"]
eau = ExtractBlock(registrationName="eau", Input=piscine, Selectors=["/Root/Fluiddomain"])  # sans les faces de bord

# écart à la moyenne volumique (mailles de même volume : moyenne simple des mailles)
ecart = PythonCalculator(registrationName="ecart", Input=eau, ArrayAssociation="Cell Data",
                         Expression="Temperature - mean(Temperature)", ArrayName="dT")
lisse = CellDatatoPointData(registrationName="lisse", Input=ecart)  # affichage lissé (mailles de 25 cm)

coupe = Slice(registrationName="coupe_y2.9", Input=lisse)
coupe.SliceType.Origin = [6.0, 2.9, 6.25]   # à travers le coin des assemblages chauds (x, y < 4 m)
coupe.SliceType.Normal = [0.0, 1.0, 0.0]
coupe.Triangulatetheslice = 0

fleches = Glyph(registrationName="vitesse", Input=coupe, GlyphType="Arrow")
fleches.OrientationArray = ["POINTS", "Velocity"]
fleches.ScaleArray = ["POINTS", "Velocity"]
fleches.ScaleFactor = 8.0                    # 10 cm/s -> flèche de 80 cm
fleches.GlyphMode = "Every Nth Point"
fleches.Stride = 3

vue = GetActiveViewOrCreate("RenderView")
vue.ViewSize = [1280, 900]
vue.OrientationAxesVisibility = 1

d = Show(coupe, vue)
ColorBy(d, ("POINTS", "dT"))
lut = GetColorTransferFunction("dT")
lut.ApplyPreset("Cool to Warm", True)
lut.AutomaticRescaleRangeMode = "Never"
lut.RescaleTransferFunction(-0.3, 0.9)
GetOpacityTransferFunction("dT").RescaleTransferFunction(-0.3, 0.9)
barre = GetScalarBar(lut, vue)
barre.Title = "T - T moyenne (K)"
barre.ComponentTitle = ""
d.SetScalarBarVisibility(vue, True)

f = Show(fleches, vue)
ColorBy(f, None)
f.AmbientColor = f.DiffuseColor = [0.1, 0.1, 0.1]

contour = Show(eau, vue)
contour.SetRepresentationType("Outline")
ColorBy(contour, None)

temps = AnnotateTimeFilter(registrationName="temps", Input=piscine)
temps.Format = "t = {time:.1f} h"
temps.Scale = 1 / 3600
Show(temps, vue).FontSize = 14
moyenne = PythonAnnotation(registrationName="T_moyenne", Input=eau, ArrayAssociation="Cell Data",
                           Expression="'T moyenne = %.1f C' % (mean(Temperature) - 273.15)")
m = Show(moyenne, vue)
m.FontSize = 14
m.WindowLocation = "Upper Right Corner"

Hide(piscine, vue)
for s in (ecart, lisse):
    Hide(s, vue)
SetActiveSource(coupe)

# regard horizontal vers +y, verticale = z (surface de l'eau en haut)
vue.CameraPosition = [6.0, -32.0, 6.25]
vue.CameraFocalPoint = [6.0, 2.9, 6.25]
vue.CameraViewUp = [0.0, 0.0, 1.0]
vue.ResetCamera(False)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = piscine.TimestepValues[-1] if piscine.TimestepValues else 0
Render(vue)

if "--png" in args:
    SaveScreenshot(args[args.index("--png") + 1], vue, ImageResolution=[1280, 900])
if "--film" in args:
    SaveAnimation(args[args.index("--film") + 1] + "/image.png", vue, ImageResolution=[1280, 900])
