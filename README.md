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
m(t)\frac{dv}{dt} = T(P_{atm}) + D(h,v) + m(t)g(r)
$$
    
Avec :
- $T$ la poussée, corrigée par la pression ambiante $P_a$.
- $D$ est la traînée aérodynamique.
- $g$ est le vecteur gravité au local dépendant de la distance radiale $r$. 

<br>

2. Modèle Atmosphérique (ISA)
   
Le code utilise le modèle International Standard Atmosphere pour calculer la densité $\rho$ et la pression $P$ jusqu'à 86 km. Il intègre la   conversion de l'altitude géométrique $z$ en altitude géopotentielle $H$ :

$$
H = \frac{R_T*z}{R_T+z}
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
- Pitch control au doigt mouillé et interpolation linéaire sur intervalles de temps.
- Terre supposée parfaitement circulaire et immobile.

<br>
<br>

**Résultats** :

Voici les résultats obtenus après moult itérations sur le pitch control et le profil de poussée des EAP :

<p align="center">
  <img width="642" height="485" alt="Capture d’écran 2026-04-21 à 21 17 01" src="https://github.com/user-attachments/assets/49960349-9f79-4da4-b973-4cff0680a0be" />
  <br><br>
  <img width="778" height="478" alt="Capture d’écran 2026-04-21 à 21 16 27" src="https://github.com/user-attachments/assets/6f951297-6948-4215-abab-72b8862ac9f0" />
</p>

<br>
<br>

**Regard Critique** :

Au vu des courbes obtenues, je suis globalement satisfait de ce petit projet. Le max-Q à seulement 7km d'altitude me dérange pas mal, j'ai lu qu'il devrait se situer aux alentours de 12km. Je suspecte un poussée des EAP trop importante dans les premières couches de l'atmosphère. Aussi, alors que le périgée de mon satellite se situe au delà des 86.000 km, mon satellite se rapproche peu à peu de la surface (intuition confirmée quand j'augmente drastiquement le temps de simulation). 


