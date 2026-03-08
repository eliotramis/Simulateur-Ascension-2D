import math
import numpy as np
import matplotlib.pyplot as plt
from fontTools.ttLib.tables.S_T_A_T_ import table_S_T_A_T_
from scipy.integrate import solve_ivp

# ------------------------
# 1 - PARAMETRES PHYSIQUES
# ------------------------

R_T = 6_378_000                # Rayon Terrestre (m)
G_0 = 9.81                     # "Gravité" au niveau de la mer (m/s^2)
P_0 = 101_325.0                # Pression au niveau de la mer (m/s^2)
T_0 = 288.15                   # Température au niveau de la mer (m/s^2)


# Modèle Ariane 5
M_0 = 780_000                  # masse décollage (kg)

D_EPC = 5.4                    # Diamètre de l'étage principal cryogénique (m)
S_EPC = math.pi*(D_EPC**2/4)   # Surface de référence EPC (m^2)
D_EAP = 3.05                   # Diamètre des 2 étages d'accél à poudre (m)
S_EAP = math.pi*(D_EAP**2/4)   # Surface de référence EAP (m^2)
S = S_EPC + 2*S_EAP            # Surface de référence Ariane5 (m^2)

CD = 0.5                       # Coeff de traînée en première approximation en vol subsonique, TODO : implémenter en fonction du Mach

T_vide = 15_490_000            # Thrust (poussée) dans le vide
F_T = T_vide                   # Force de thrust en première approximation égale à T_vide (N), TODO : implémenter en fonction de l'altitude (pression)

ISP = 284.0                    # Impulsion spécifique moyenne calculée en première approximation (s) TODO : implémenter en fonction EAP ou EPC et altitude (pression)


# -----------------------------
# 2 - Calcul des modèles locaux
# -----------------------------

def gravite(y: float):
    return G_0*(R_T/(R_T+y))**2

def atmosphere(y: float):

    if y > 100_000.0:
        return 0.0

    if y < 11_000.0:
        # Troposphère, modèle de T linéaire
        T = T0

# --------------------------------
# 3 - Calcul de mon vecteur d'état
# --------------------------------

def dynamique_fusee(t, state):
    x, y, vx, vy, m = state

    # conditions d'arrêt si fusée dans le sol ou masse négative
    if y<0 and y < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    if m < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    # modèles locaux
    g = gravite(y)
    rho = atmosphere(y)


if __name__ == "__main__":
    print("Hello World")
