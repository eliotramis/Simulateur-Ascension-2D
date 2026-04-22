# Ariane 5 Flight Dynamics Simulator (2D)

Ce projet implémente un simulateur de trajectoire d'ascension pour le lanceur Ariane 5. Le script modélise les phases de vol, de l'allumage jusqu'à la mise en orbite, en prenant en compte les variations environnementales et les changements de configuration du lanceur.

<br>

**Fonctionnalités** :

- Gestion du largage des EAP (Boosters) et de l'EPC (Étage Principal).

- Modèle atmosphérique ISA (7 couches) et gravité variable en $\frac{1}{r^2}$ .

- Coefficient de traînée $C_d$​ variable selon le nombre de Mach.

- Calcul automatique de la pression dynamique maximale (Max-Q) pour analyse.

<br>
<br>

**Modèles Physiques** :

1. Équation de la dynamique

Le simulateur résout le principe fondamental de la dynamique dans le référentiel géocentrique supposé galiléen :
    
$$
m(t)\frac{dv}{dt} = T(P_a) + D(h,v) + m(t)g(r)
$$
    
Avec :
- $T$ la poussée, corrigée par la pression ambiante $P_a$.
- $D$ est la traînée aérodynamique.
- $g$ est le vecteur gravité au local dépendant de la distance radiale $r$. 

<br>

2. Modèle Atmosphérique (ISA)
   
Le code utilise le modèle International Standard Atmosphere pour calculer la densité de l'air $\rho$ et la pression $P$ jusqu'à 86 km. Il intègre la conversion de l'altitude géométrique $z$ en altitude géopotentielle $H$ :

$$
H = \frac{R_T \cdot z}{R_T+z}
$$

  Avec $R_T$ le rayon terrestre.

<br>

3. Déroulé du vol

Le simulateur se décompose en 3 phases successives :
- Allumage de l'EPC et des 2 EAP.
- Séparation des EAP.
- Séparation de l'EPC et allumage de l'ESC (Étage Supérieur Cryogénique).

<br>
<br>

**Hypothèses physiques** :

- Plan (x,y).
- Loi de commande d'attitude (Pitch) définie de manière empirique par interpolation linéaire sur intervalles de temps.
- Terre parfaitement sphérique, non tournante (absence des composantes de Coriolis et d'entraînement).

<br>
<br>

**Résultats** :

Voici les résultats obtenus après plusieurs campagnes d'itérations sur le pitch control et le profil de poussée des EAP :

<p align="center">
  <img width="642" height="485" alt="Capture d’écran 2026-04-21 à 21 17 01" src="https://github.com/user-attachments/assets/49960349-9f79-4da4-b973-4cff0680a0be" />
  <br><br>
  <img width="778" height="478" alt="Capture d’écran 2026-04-21 à 21 16 27" src="https://github.com/user-attachments/assets/6f951297-6948-4215-abab-72b8862ac9f0" />
</p>

<br>
<br>

**Regard Critique** :

Bien que le simulateur produise un profil de vol cohérent, l'analyse des résultats met en évidence deux biais de modélisation qui nécessiteront des itérations futures :

- Le Max-Q est atteint à une altitude d'environ 7 km, alors que la télémétrie réelle d'Ariane 5 le situe proche des 13/14 km. Cette anomalie suggère une interpolation de la poussée des EAP trop agressive dans les basses couches de l'atmosphère, ou l'absence d'une loi de throttling explicite lors du passage en supersonique.

- On observe une dégradation de l'orbite sur le long terme, malgré une vitesse finale suffisante (8.3 km/s). Cela montre qu'à l'extinction du moteur, le lanceur n'est pas tout à fait parallèle à la surface terrestre (angle de pente de vol non nul). La consigne d'inclinaison sur la dernière phase de vol mérite donc d'être un peu mieux réglée. Ou alors le solveur accumule les erreurs sur les temps longs.

