# Optimisation d'un Lanceur Orbital à Trois Étages

Ce dépôt contient le code source et les résultats d'un projet d'étude visant à déterminer numériquement la configuration optimale d'un lanceur orbital à trois étages.

## Contexte et Objectifs

L'objectif principal du projet est de concevoir un lanceur capable de placer une masse utile de 1000 kg sur une orbite située à 250 km d'altitude au-dessus de la surface terrestre.La vitesse cible pour atteindre cette orbite est d'environ 7754.8 m/s. 

La résolution de ce problème de modélisation repose sur deux critères majeurs :
* Minimiser la masse totale du lanceur.
* S'assurer que la charge utile atteigne l'orbite avec la vitesse cible de manière perpendiculaire, en minimisant l'erreur finale.

## Fonctionnalités Principales

* **Optimiseur SQP Personnalisé :** Un algorithme d'optimisation robuste de type SQP (Sequential Quadratic Programming) a été entièrement développé en MATLAB.Il intègre le calcul de gradients par différences finies , l'approximation des matrices hessiennes via des méthodes quasi-Newton (SR1 ou BFGS), et une étape de globalisation par recherche linéaire d'Armijo.
* **Simulateur de Trajectoire Numérique :** Un simulateur basé sur le solveur ode45 de MATLAB a été mis en place pour résoudre l'équation différentielle ordinaire régissant le vol de la fusée, de l'allumage du premier étage jusqu'au largage du troisième.
* **Itérations Étagement-Vitesse :** Implémentation d'un processus itératif complexe couplant l'optimisation de l'étagement (maximisation du ratio masse utile / masse totale) et la simulation de trajectoire pour ajuster la vitesse propulsive et les angles de poussée.

## Résultats Obtenus

L'algorithme a permis de converger vers un dimensionnement optimal respectant le cahier des charges avec une grande précision:
* **Masse totale du lanceur :** 31193.19 kg
* **Pourcentage de masse utile :** 3.20 % 
***Erreur sur l'altitude finale :** 1.443e-7 m
* **Erreur sur la vitesse finale :** 7.7087e-5 m/s 
* ]**Angles de poussée optimaux :** 39.03°, -9.47°, -26.63°, -37.12° 

## Technologies Utilisées

*MATLAB 
* Solveur ODE45 (Équations Différentielles Ordinaires)

## Structure du Dépôt

* `resolution_pb.m` : Script principal orchestrant les itérations étagement-vitesse pour déterminer le lanceur optimal.
* `test_vmax.m` : Script utilisé pour la recherche expérimentale des angles de poussée et l'optimisation de la vitesse.
* `resolution_analytique.m` : Implémentation de la méthode de Newton pour la résolution des conditions KKT du problème d'étagement.
* `SQP/` : Répertoire contenant l'interface de l'algorithme d'optimisation (`SQP.m`) et ses scripts de validation (`test_SQP.m`).
* `Simulateur_trajectoire/` : Répertoire contenant l'interface du simulateur (`simulateur_trajectoire.m`), le script de visualisation graphique (`tracer_traj.m`) et les tests de bon fonctionnement.
