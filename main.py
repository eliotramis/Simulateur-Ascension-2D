import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, RK45
from scipy.stats import alpha

# --------------
# 1 - PARAMETRES
# --------------

# - Paramètres physiques -

R_T = 6_378_000.0              # Rayon Terrestre (m)
G_0 = 9.80665                  # "Gravité" au niveau de la mer (m/s^2)
R_S_air = 287.058              # Constante spécifique de l'air sec (J/(kg.K))
GAMMA_air = 1.4

MACH_DATA = np.array([0.0,  0.8,  0.9,  1.0,  1.1,  1.2,  1.5,  2.0,  3.0,  5.0,  10.0, 25.0])
CD_DATA   = np.array([0.45, 0.45, 0.55, 0.75, 0.85, 0.80, 0.65, 0.50, 0.35, 0.25, 0.15, 0.10])

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

T_vide = 15_490_000            # Thrust (poussée) cumulée dans le vide des deux EAP et du EPC (N)

ISP = 284.0                    # Impulsion spécifique moyenne calculée en première approximation (s) TODO : implémenter en fonction EAP ou EPC et altitude (pression)


# - Paramètres solveur -

t_0 = 0.0                      # Temps initial pour solveur (s)
t_f = 10000.0                  # Temps final pour solveur (s)
pas = 0.1                      # Pas (s)
nb_pas = int((t_f-t_0)/pas)    # Nombre de pas pour solveur

# ------------------
# 2 - CALCUL MODELES
# ------------------

def gravite(y: float):

    if y < 0.0:
        return 0.0

    return G_0*(( R_T / (R_T + y))**2)


def atmosphere(y):

    """
    Modèle ISA complet (7 couches).
    Convertit l'altitude géométrique y en altitude géopotentielle H
    """

    # Au-delà de 86 km, on considère une densité de l'air nulle
    if y >= 86000.0:
        return 0.0, 186.87
    if y < 0.0:
        return 0.0, 288.15

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

    return rho, T

def vitesse_air(T):
    if T <= 0:
        return 0.0
    return np.sqrt(GAMMA_air * R_S_air * T)

def nombre_de_mach(T, v):
    if v <= 0.0:
        return 0.0
    return v / vitesse_air(T)

def calculer_cd(T, v):

    mach_actuel = nombre_de_mach(T, v)

    cd_actuel = np.interp(mach_actuel, MACH_DATA, CD_DATA)

    return cd_actuel

# --------------------------------
# 3 - Calcul de mon vecteur d'état
# --------------------------------

def dynamique_fusee(t, state):

    x, y, vx, vy, m = state

    # Conditions d'arrêt si fusée dans le sol ou masse négative
    if y <= 0 and vy < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    if m < 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    # tant qu'il y'a de l'érgols, on continue
    if m > M_vide:
        # Thrust -> constant en première approximation
        F_T_norme = T_vide
        # Débit massique, négatif, masse diminue
        dm_dt = - F_T_norme / (ISP * G_0)

    # sinon, on coupe les moteurs
    else:
        F_T_norme = 0.0
        dm_dt = 0.0

    # securité sur masse si plus de carburant
    if m < M_vide:
        m = M_vide

    r_vect = np.array([x, R_T + y])
    r_norme = np.linalg.norm(r_vect)
    altitude_locale = r_norme - R_T

    v_vect = np.array([vx, vy])
    v_norme = np.linalg.norm(v_vect)

    g_locale = gravite(altitude_locale)
    rho_locale, T_locale = atmosphere(altitude_locale)

    # POIDS
    poids_direction = -r_vect/r_norme
    P_vect = m * g_locale * poids_direction


    # DRAG
    D_vect = np.array([0.0, 0.0])

    C_D = calculer_cd(T_locale, v_norme)

    if v_norme > 1.0 and altitude_locale < 86_000.0:
        D_norme = 0.5 * rho_locale * (v_norme ** 2) * C_D * S
        D_vect = D_norme * (- v_vect / v_norme)


    # THRUST
    T_vect = np.array([0.0, 0.0])
    if F_T_norme > 0.0:
        if v_norme < 50.0:
            T_vect = np.array([0.0, F_T_norme])

        else:
            angle_pitch_rad = np.radians(80.0)
            T_vect = np.array([F_T_norme * math.cos(angle_pitch_rad), F_T_norme * math.sin(angle_pitch_rad)])

    PFD_vect = P_vect + T_vect + D_vect

    accel = PFD_vect/m

    return [vx, vy, accel[0], accel[1], dm_dt]



if __name__ == "__main__":

    def impact_sol(t, state):
        # index 1 du state correspond à l'altitude y
        return state[1]

    impact_sol.terminal = True  # OUI : coupe le solveur quand y = 0
    impact_sol.direction = -1  # OUI : ne s'active que si y passe par 0 en descendant

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
        events=impact_sol,
        method='RK45'

    )

    temps = solution.t
    positions = solution.y[0]
    altitudes = solution.y[1]
    vitesses = np.sqrt(solution.y[2] ** 2 + solution.y[3] ** 2)

    densites = np.array([atmosphere(y)[0] for y in altitudes])

    pressions_dyn = 0.5*densites*(vitesses**2)

    indices_montee = np.where(temps < 130.0)[0]

    # On cherche le max-Q que dans cet intervalle
    index_max_Q_montee = indices_montee[np.argmax(pressions_dyn[indices_montee])]

    valeur_max_Q = pressions_dyn[index_max_Q_montee]
    altitude_max_Q = altitudes[index_max_Q_montee]
    temps_max_Q = temps[index_max_Q_montee]
    vitesse_max_Q = vitesses[index_max_Q_montee]

    print("=== RAPPORT DE VOL ===")
    print(f"Apogée atteint : {np.max(altitudes) / 1000:.2f} km")
    print(f"Max-Q : {valeur_max_Q:.2f} Pascals")
    print(f" - Survenu à t = {temps_max_Q:.2f} s")
    print(f" - Altitude = {altitude_max_Q / 1000:.2f} km")
    print(f" - Vitesse = {vitesse_max_Q:.2f} m/s")

    plt.figure(figsize=(12, 10))

    # Graphe 1 : Profil d'altitude
    plt.subplot(4, 1, 1)
    plt.plot(temps, altitudes / 1000.0, 'b-', linewidth=2, label='Altitude (km)')
    plt.axvline(x=temps_max_Q, color='r', linestyle='--', alpha=0.5)
    plt.ylabel('Altitude [km]')
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.title('Relevé de Vol (Modèle Vertical 1D)')

    # Graphe 2 : Profil de vitesse
    plt.subplot(4, 1, 2)
    plt.plot(temps, vitesses, 'g-', linewidth=2, label='Vitesse (m/s)')
    plt.axvline(x=temps_max_Q, color='r', linestyle='--', alpha=1)
    plt.ylabel('Vitesse [m/s]')
    plt.grid(True, linestyle=':')
    plt.legend()

    # Graphe 3 : Pression Dynamique
    plt.subplot(4, 1, 3)
    plt.plot(temps, pressions_dyn, 'r-', linewidth=2, label='Pression Dynamique $q$ (Pa)')

    # Marqueur visuel pour le Max-Q
    plt.plot(temps_max_Q, valeur_max_Q, 'ko')
    plt.annotate(f'Max-Q: {valeur_max_Q:.0f} Pa\n@ {altitude_max_Q / 1000:.1f} km',
                 xy=(temps_max_Q, valeur_max_Q),
                 xytext=(temps_max_Q + 10, valeur_max_Q * 0.8))
    plt.ylim(0, valeur_max_Q * 1.2)
    plt.ylabel('Pression $q$ [Pa]')
    plt.grid(True, linestyle=':')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.plot(temps, positions / 1000.0, 'b-', linewidth=2, label='Position (km)')
    plt.axvline(x=temps_max_Q, color='r', linestyle='--', alpha=0.5)
    plt.xlabel('Temps [s]')
    plt.ylabel('Position [km]')
    plt.grid(True, linestyle=':')
    plt.legend()

    #plt.tight_layout()
    plt.show()



