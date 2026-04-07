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

# - Modèle ARIANE 5 -

M_0 = 780_000                  # masse totale décollage (kg)
M_vide = 150_000               # masse fusée sans ergols (kg)
M_PAYLOAD = 10_000.0           # Masse des satellites (kg)
M_ESC_VIDE = 4_500.0           # Masse de la structure de l'étage sup (kg)
MASSE_ORBITALE_FINALE = M_PAYLOAD + M_ESC_VIDE # Masse totale en orbite (kg)

D_EPC = 5.4                    # Diamètre de l'étage principal cryogénique (m)
S_EPC = math.pi*(D_EPC**2)/4   # Surface de référence EPC (m^2)
D_EAP = 3.05                   # Diamètre des 2 étages d'accél à poudre (m)
S_EAP = math.pi*(D_EAP**2)/4   # Surface de référence EAP (m^2)
S = S_EPC + 2*S_EAP            # Surface de référence Ariane5 (m^2)

# Profil de poussée dynamique d'un EAP (1 seul booster)
T_EAP_TEMPS = np.array([0.0, 20.0, 50.0, 90.0, 130.0, 135.0])
T_EAP_FORCE = np.array([6500000.0, 7000000.0, 4200000.0, 6000000.0, 900000.0, 0.0])
ISP_EAP = 275.0
A_E_EAP = 2.9

# Constantes de l'EPC
T_EPC_CONSTANTE = 1_350_000.0 # (N)
ISP_EPC = 432.0
A_E_EPC = 3.0

# Constantes de L'ESC
T_ESC = 64_800.0
ISP_ESC = 446.0





# ------------------
# 2 - CALCUL MODELES
# ------------------

def distance_radiale(x,y):
    return np.sqrt(x**2 + (y+R_T)**2)

def gravite_norme(r: float):

    if r < 0.0:
        return 0.0

    return G_0*(( R_T / r )**2)


def atmosphere(y):

    """
    Modèle ISA complet (7 couches).
    Convertit l'altitude géométrique y en altitude géopotentielle H
    """

    # Au-delà de 86 km, on considère une densité de l'air nulle
    if y >= 86000.0:
        return 0.0, 186.87, 0.0
    if y < 0.0:
        return 0.0, 288.15, 101325.0

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

    return rho, T, P

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
# 3 - CALCUL de mon VECTEUR d'état
# --------------------------------

def dynamique_fusee(t, state):

    x, y, vx, vy, m = state

    # masse en fonction des phases
    if t <= 135.0:
        masse_vide_actuelle = 150_000.0
    elif t <= 540.0:
        masse_vide_actuelle = 74_000.0
    else:
        masse_vide_actuelle = 14_000.0

    # cacluls parametres pour fct
    r_vect = np.array([x, R_T + y])
    r_norme = np.linalg.norm(r_vect)
    altitude_locale = r_norme - R_T

    v_vect = np.array([vx, vy])
    v_norme = np.linalg.norm(v_vect)

    r = distance_radiale(x, y)
    g_r = gravite_norme(r)
    rho_locale, T_locale, P_locale = atmosphere(altitude_locale)

    if m <= MASSE_ORBITALE_FINALE:
        m = MASSE_ORBITALE_FINALE
        F_T_norme = 0.0
        dm_dt = 0.0

    else:
        if t <= 135.0:
            # PHASE 1 : EAP + EPC
            T_EPC_vide = T_EPC_CONSTANTE
            T_EPC_reel = T_EPC_vide - (A_E_EPC * P_locale)
            dm_dt_EPC = - (T_EPC_vide / (ISP_EPC * G_0))

            T_EAP_vide_unitaire = np.interp(t, T_EAP_TEMPS, T_EAP_FORCE)
            T_EAP_reel_unitaire = T_EAP_vide_unitaire - (A_E_EAP * P_locale)
            dm_dt_EAP = - (T_EAP_vide_unitaire / (ISP_EAP * G_0))

            F_T_norme = T_EPC_reel + (2 * T_EAP_reel_unitaire)
            dm_dt = dm_dt_EPC + (2 * dm_dt_EAP)

        elif t <= 540.0:
            # PHASE 2 : EPC seul
            F_T_norme = T_EPC_CONSTANTE - (A_E_EPC * P_locale)
            dm_dt = - (T_EPC_CONSTANTE / (ISP_EPC * G_0))

        else:
            # PHASE 3 : ESC seul
            F_T_norme = T_ESC
            dm_dt = - (T_ESC / (ISP_ESC * G_0))


    # - POIDS -
    P_vect = np.array([m*(-g_r*x/r), m*(-g_r*(y+R_T)/r)])


    # - DRAG -
    D_vect = np.array([0.0, 0.0])

    C_D = calculer_cd(T_locale, v_norme)

    if v_norme > 1.0 and altitude_locale < 86_000.0:
        D_norme = 0.5 * rho_locale * (v_norme ** 2) * C_D * S
        D_vect = D_norme * (- v_vect / v_norme)


    # - THRUST -
    T_vect = np.array([0.0, 0.0])
    if F_T_norme > 0.0:
        if t < 15.0:
            theta_deg = 90.0
        elif t <= 135.0:
            # phase EAP : 90° à 30°
            theta_deg = np.interp(t, [15.0, 135.0], [90.0, 50.0])
        elif t <= 540.0:
            # phase EPC : 30° à 0°
            theta_deg = np.interp(t, [135.0, 540.0], [50.0, 15.0])
        else:
            # phase ESC : horizontale
            theta_deg = np.interp(t, [540.0, 1400.0], [15.0, 0.0])

        # conversion repère sphérique
        angle_radial = math.atan2(y + R_T, x)
        horizontale_locale = angle_radial - (math.pi / 2)

        # angle final = l'horizontale locale + consigne
        theta = horizontale_locale + math.radians(theta_deg)

        # 3. Projection du vecteur Poussée
        T_vect = np.array([F_T_norme * math.cos(theta), F_T_norme * math.sin(theta)])


    PFD_vect = P_vect + T_vect + D_vect

    accel = PFD_vect/m

    return [vx, vy, accel[0], accel[1], dm_dt]





# --------------
# 4 - RESOLUTION
# --------------

if __name__ == "__main__":

    def impact_sol(t, state):
        x = state[0]
        y = state[1]
        # Vraie altitude = distance au noyau terrestre - rayon de la Terre
        altitude_vraie = math.sqrt(x**2 + (y + R_T)**2) - R_T
        return altitude_vraie

    impact_sol.terminal = True  # oui : coupe le solveur quand y = 0
    impact_sol.direction = -1  # oui : ne s'active que si y passe par 0 en descendant

    t_f = 6500.0
    pas = 0.1

    # ===============================
    # PHASE 1 : ASCENSION (0 à 135 s)
    # ===============================
    state0 = [0.0, 0.0, 0.0, 0.0, M_0]
    t_span_1 = (0.0, 135.0)
    t_mesures_1 = np.linspace(0.0, 135.0, int(135.0 / pas))

    solution_1 = solve_ivp(fun=dynamique_fusee,
                           t_span=t_span_1,
                           y0=state0,
                           t_eval=t_mesures_1,
                           method='RK45')

    # LARGAGE DES EAP (On retire 76 tonnes)
    state_largage_1 = [solution_1.y[0][-1],
                       solution_1.y[1][-1],
                       solution_1.y[2][-1],
                       solution_1.y[3][-1],
                       solution_1.y[4][-1] - 76000.0]



    # ==========================================
    # PHASE 2 : ÉTAGE PRINCIPAL (135 à 540 s)
    # ==========================================
    t_span_2 = (135.0, 540.0)
    t_mesures_2 = np.linspace(135.0, 540.0, int((540.0 - 135.0) / pas))

    solution_2 = solve_ivp(fun=dynamique_fusee,
                           t_span=t_span_2,
                           y0=state_largage_1,
                           t_eval=t_mesures_2,
                           events=impact_sol,
                           method='RK45')

    # LARGAGE DE L'EPC (On retire 14 tonnes)
    state_largage_2 = [solution_2.y[0][-1],
                       solution_2.y[1][-1],
                       solution_2.y[2][-1],
                       solution_2.y[3][-1],
                       solution_2.y[4][-1] - 14000.0]



    # =================================================
    # PHASE 3 : ÉTAGE SUPÉRIEUR & ORBITE (540 à 4000 s)
    # =================================================
    t_span_3 = (540.0, t_f)
    t_mesures_3 = np.linspace(540.0, t_f, int((t_f - 540.0) / pas))

    solution_3 = solve_ivp(fun=dynamique_fusee,
                           t_span=t_span_3,
                           y0=state_largage_2,
                           t_eval=t_mesures_3,
                           events=impact_sol,
                           method='RK45')

    # =========================
    # CONCATENATION DES DONNÉES
    # =========================

    temps = np.concatenate((solution_1.t, solution_2.t, solution_3.t))
    positions = np.concatenate((solution_1.y[0], solution_2.y[0], solution_3.y[0]))
    y_coords = np.concatenate((solution_1.y[1], solution_2.y[1], solution_3.y[1]))
    vitesses_x = np.concatenate((solution_1.y[2], solution_2.y[2], solution_3.y[2]))
    vitesses_y = np.concatenate((solution_1.y[3], solution_2.y[3], solution_3.y[3]))

    altitudes = np.sqrt(positions ** 2 + (y_coords + R_T) ** 2) - R_T

    vitesses = np.sqrt( vitesses_x ** 2 + vitesses_y ** 2)

    densites = np.array([atmosphere(y)[0] for y in altitudes])

    pressions_dyn = 0.5*densites*(vitesses**2)

    index_montee = np.where(temps < 130.0)[0]

    # On cherche le max-Q que dans cet intervalle
    index_max_Q_montee = index_montee[np.argmax(pressions_dyn[index_montee])]

    valeur_max_Q = pressions_dyn[index_max_Q_montee]
    altitude_max_Q = altitudes[index_max_Q_montee]
    position_max_Q = positions[index_max_Q_montee]
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
    plt.title('Relevé de Vol (Modèle 2D, pitch angle non constant)')
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

    # Graphe 4 : altitude/position

    plt.subplot(4, 1, 4)

    # courbure de la Terre (un cercle complet)
    theta_terre = np.linspace(0, 2 * math.pi, 500)
    x_terre = R_T * np.cos(theta_terre)
    y_terre = R_T * np.sin(theta_terre) - R_T
    plt.plot(x_terre / 1000.0, y_terre / 1000.0, 'g-', linewidth=2, label='Terre (Rayon 6378 km)')

    # atmosphère
    x_atm = (R_T + 100000.0) * np.cos(theta_terre)
    y_atm = (R_T + 100000.0) * np.sin(theta_terre) - R_T
    plt.plot(x_atm / 1000.0, y_atm / 1000.0, 'c--', linewidth=1, alpha=0.5, label='Atmosphère (100 km)')

    # trajectoire cartésienne
    plt.plot(positions / 1000.0, y_coords / 1000.0, 'b-', linewidth=1, label='Trajectoire Fusée')

    # affichage isométrique
    plt.axis('equal')

    # dézoom (Fenêtre de 16 000 km par 16 000 km)
    plt.xlim(-4000, 4000)
    plt.ylim(-15000, 1000)  # L'origine Y=0 est en haut du globe, donc on descend jusqu'à -15000

    plt.xlabel('X [km]')
    plt.ylabel('Y [km]')
    plt.grid(True, linestyle=':')
    plt.legend(loc='upper right')

    #plt.tight_layout()
    plt.show()


