import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, RK45
from scipy.stats import alpha

# --------------
# 1 - PARAMETRES
# --------------

# - Paramètres physiques -

R_T = 6_378_000                # Rayon Terrestre (m)
G_0 = 9.80665                  # "Gravité" au niveau de la mer (m/s^2)
R_S_air = 287.058              # Constante spécifique de l'air sec (J/(kg.K))

H_BASES = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 86000.0]) # Altitudes de chaque couche (m)
L_GRAD = np.array([-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002]) # Gradients thermiques de chaque couche (K/m)

T_BASES = np.zeros(7)
P_BASES = np.zeros(7)

T_BASES[0] = 288.15            # Température au niveau de la mer (m/s^2)
P_BASES[0] = 101_325.0         # Pression au niveau de la mer (m/s^2)

for i in range(1, 7):
    dh = H_BASES[i] - H_BASES[i-1]
    L = L_GRAD[i-1]
    T_prev = T_BASES[i-1]
    P_prev = P_BASES[i-1]

    if L != 0.0:
        T_next = T_prev + L * dh
        P_next = P_prev * (T_next / T_prev)**(-G_0 / (R_S_air * L))
    else:
        T_next = T_prev
        P_next = P_prev * np.exp(-G_0 * dh / (R_S_air * T_prev))

    T_BASES[i] = T_next
    P_BASES[i] = P_next

# - Modèle Ariane 5 -

M_0 = 780_000                  # masse décollage (kg)
M_vide = 150_000               # masse fusée sans ergols (kg)

D_EPC = 5.4                    # Diamètre de l'étage principal cryogénique (m)
S_EPC = math.pi*(D_EPC**2)/4   # Surface de référence EPC (m^2)
D_EAP = 3.05                   # Diamètre des 2 étages d'accél à poudre (m)
S_EAP = math.pi*(D_EAP**2)/4   # Surface de référence EAP (m^2)
S = S_EPC + 2*S_EAP            # Surface de référence Ariane5 (m^2)

CD = 0.5                       # Coeff de traînée en première approximation en vol subsonique, TODO : implémenter en fonction du Mach

T_vide = 15_490_000            # Thrust (poussée) cumulée dans le vide des deux EAP et du EPC (N)
F_T = T_vide                   # Force de thrust en première approximation égale à T_vide (N), TODO : implémenter en fonction de l'altitude (pression)

ISP = 284.0                    # Impulsion spécifique moyenne calculée en première approximation (s) TODO : implémenter en fonction EAP ou EPC et altitude (pression)


# - Paramètres solveur -

t_0 = 0.0                      # Temps initial pour solveur (s)
t_f = 300.0                    # Temps final pour solveur (s)
pas = 0.1                      # Pas (s)
nb_pas = int((t_f-t_0)/pas)    # Nombre de pas pour solveur

# ------------------
# 2 - CALCUL MODELES
# ------------------

def gravite(y: float):
    return G_0*((R_T/(R_T+y))**2)


def atmosphere(y):
    """
    Modèle ISA complet (7 couches).
    Convertit l'altitude géométrique y en altitude géopotentielle H,
    puis interpole la densité locale.
    """
    # Au-delà de 86 km, on considère une densité de l'air nulle
    if y >= 86000.0:
        return 0.0

    # Conversion géopotentielle
    H = (R_T * y) / (R_T + y)

    # Identification de la couche actuelle (0 à 6)
    couche = 0
    while couche < 6 and H >= H_BASES[couche + 1]:
        couche += 1

    # Extraction des constantes déjà calculées
    Hb = H_BASES[couche]
    Tb = T_BASES[couche]
    Pb = P_BASES[couche]
    L = L_GRAD[couche]

    # Hauteur géopotentielle ajustée pour prendre l'altitude de la couche comme point de réf
    dh = H - Hb

    # Résolution thermodynamique
    # Si gradient de température : modèle linéaire de la température
    if L != 0.0:
        T = Tb + L * dh
        P = Pb * (T / Tb) ** (-G_0 / (R_S_air * L))

    # Sinon, T = cte et on établit que P suit une exponentielle décroissante
    else:
        T = Tb
        P = Pb * np.exp(-G_0 * dh / (R_S_air * Tb))

    # Gaz parfaits
    rho = P / (R_S_air * T)

    return rho

# --------------------------------
# 3 - Calcul de mon vecteur d'état
# --------------------------------

def dynamique_fusee(t, state):

    x, y, vx, vy, m = state

    # Conditions d'arrêt si fusée dans le sol ou masse négative
    if y < 0 and vy < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    if m < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    # - Modèles locaux.
    # Norme du champ de gravité local
    g = gravite(y)

    # densité de l'air à l'altitude y
    rho = atmosphere(y)

    # - Calcul des forces
    # Poids
    P = m*g

    # tant qu'il y'a de l'érgols, on continue
    if m > M_vide:
        # Thrust -> constant en première approximation
        T = F_T
        # Débit massique, négatif, masse diminue
        dm_dt = -T / (ISP * G_0)

    # sinon, on coupe les moteurs
    else:
        T = 0.0
        dm_dt = 0.0

    # securité sur masse si pas trop grand
    if m < M_vide:
        m = M_vide

    # Drag
    D = (1/2)*rho*(vy**2)*CD*S

    # - PFD

    # selon x - nulle car tir vertical -> en première approximation
    dvx_dt = 0.0

    # Selon y
    dvy_dt = (T - P - D)/m

    dx_dt = vx
    dy_dt = vy

    return [dx_dt, dy_dt, dvx_dt, dvy_dt, dm_dt]



if __name__ == "__main__":

    # initialisation états initiaux
    state0 = [0.0, 0.0, 0.0, 0.0, M_0]
    # durée pour solveur
    t_ecart = (t_0, t_f)
    # tableau de mesures
    t_mesures = np.linspace(t_0, t_f, nb_pas)

    solution = solve_ivp(
        fun = dynamique_fusee,
        t_span=t_ecart,
        y0=state0,
        t_eval=t_mesures,
        method='RK45'
    )

    temps = solution.t
    altitudes = solution.y[1]
    vitesses_y = solution.y[3]

    densites = np.array([atmosphere(y) for y in altitudes])

    pressions_dyn = 0.5*densites*(vitesses_y**2)

    index_max_Q = np.argmax(pressions_dyn)

    valeur_max_Q = pressions_dyn[index_max_Q]
    altitude_max_Q = altitudes[index_max_Q]
    temps_max_Q = temps[index_max_Q]
    vitesse_max_Q = vitesses_y[index_max_Q]

    print("=== RAPPORT DE VOL ===")
    print(f"Apogée atteint : {np.max(altitudes) / 1000:.2f} km")
    print(f"Max-Q : {valeur_max_Q:.2f} Pascals")
    print(f" - Survenu à t = {temps_max_Q:.2f} s")
    print(f" - Altitude = {altitude_max_Q / 1000:.2f} km")
    print(f" - Vitesse = {vitesse_max_Q:.2f} m/s")

    plt.figure(figsize=(10, 8))

    # Graphe 1 : Profil d'altitude
    plt.subplot(3, 1, 1)
    plt.plot(temps, altitudes / 1000.0, 'b-', linewidth=2, label='Altitude (km)')
    plt.axvline(x=temps_max_Q, color='r', linestyle='--', alpha=0.5)
    plt.ylabel('Altitude [km]')
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.title('Relevé de Vol (Modèle Vertical 1D)')

    # Graphe 2 : Profil de vitesse
    plt.subplot(3, 1, 2)
    plt.plot(temps, vitesses_y, 'g-', linewidth=2, label='Vitesse (m/s)')
    plt.axvline(x=temps_max_Q, color='r', linestyle='--', alpha=1)
    plt.ylabel('Vitesse [m/s]')
    plt.grid(True, linestyle=':')
    plt.legend()

    # Graphe 3 : Pression Dynamique (Le point critique)
    plt.subplot(3, 1, 3)
    plt.plot(temps, pressions_dyn, 'r-', linewidth=2, label='Pression Dynamique $q$ (Pa)')

    # Marqueur visuel pour le Max-Q
    plt.plot(temps_max_Q, valeur_max_Q, 'ko')
    plt.annotate(f'Max-Q: {valeur_max_Q:.0f} Pa\n@ {altitude_max_Q / 1000:.1f} km',
                 xy=(temps_max_Q, valeur_max_Q),
                 xytext=(temps_max_Q + 10, valeur_max_Q * 0.7),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6))

    plt.xlabel('Temps [s]')
    plt.ylabel('Pression $q$ [Pa]')
    plt.grid(True, linestyle=':')
    plt.legend()

    plt.tight_layout()
    plt.show()



