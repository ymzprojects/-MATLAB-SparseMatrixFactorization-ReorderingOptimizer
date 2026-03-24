# 🚀 Sparse Matrix Factorization & Reordering Optimizer

![MATLAB](https://img.shields.io/badge/MATLAB-R2023b+-blue?logo=mathworks)
![Sparse](https://img.shields.io/badge/Matrices-Sparse_Optimization-green)
![Status](https://img.shields.io/badge/Performance-Numerical_Analysis-success)

**Sparse Matrix Factorization & Reordering Optimizer** est un moteur d'analyse numérique conçu pour optimiser la décomposition de matrices creuses de grande taille. 

Le projet automatise l'évaluation des algorithmes de réordonnancement pour minimiser le phénomène de **fill-in** (remplissage), réduisant ainsi drastiquement l'occupation mémoire et la complexité de calcul lors des factorisations **Cholesky** et **LU**.

---

## 🏗️ Architecture & Méthodologie

### 1. Algorithmes de Réordonnancement
L'outil compare plusieurs approches de permutation pour transformer la structure de sparsité originale en une forme optimisée :
* [cite_start]**Minimum Degree (AMD/COLAMD) :** Réduit le fill-in en sélectionnant les pivots de plus bas degré[cite: 3].
* [cite_start]**Symmetric Reverse Cuthill-McKee (SYMRCM) :** Réduit la largeur de bande de la matrice pour localiser les données autour de la diagonale[cite: 1, 3].
* [cite_start]**Symmetric Approximate Minimum Degree (SYMAMD) :** Optimise spécifiquement les structures symétriques pour Cholesky[cite: 1].

### 2. Moteurs Numériques
Le pipeline traite deux types de décompositions fondamentales via les scripts fournis :
* [cite_start]**Cholesky Engine (`tp_chol.m`) :** Dédié aux matrices symétriques définies positives (SPD) comme `mat3`[cite: 1].
* **LU Engine (`tp_lu.m`) :** Pour les matrices carrées générales, intégrant des permutations de colonnes externes pour la sparsité[cite: 2].

---

## 📊 Métriques d'Ingénierie

Les tests réalisés sur les jeux de données démontrent une efficacité critique :

* **Optimisation Mémoire :** Pour Cholesky, le nombre de non-zéros ($nnz$) passe de **2 814 260** à **231 220** grâce à `symamd`, soit une réduction massive du stockage[cite: 1, 2].
* [cite_start]**Gain Computationnel :** Sur LU, les opérations flottantes ($FLOPs$) chutent de **33 028** à **11 542** avec `colamd`[cite: 2, 3].
* [cite_start]**Contrôle de Précision :** Maintien d'une erreur résiduelle extrêmement faible, de l'ordre de $10^{-19}$ pour Cholesky et $10^{-17}$ pour LU[cite: 1, 2].

---

## 📦 Installation & Utilisation

```matlab
% Cloner le dépôt et ouvrir MATLAB
% Pour lancer l'analyse de Cholesky (Matrices Symétriques)
run('tp_chol.m')

% Pour lancer l'analyse LU (Matrices Générales)
run('tp_lu.m')

% Les rapports visuels et fichiers PNG sont générés automatiquement
