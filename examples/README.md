# Osinco3D - Mini-Guide Débutant

Ce guide rapide explique comment compiler et lancer une simulation **Taylor-Green Vortex** avec Osinco3D.

---

## 1. Cloner le dépôt

```bash
git clone https://github.com/jojoledemago/osinco3d.git
cd osinco3d/src/
```

---

## 2. Compiler le code

```bash
make
```

- Pour nettoyer la compilation :

```bash
make clean
make
```

- L'exécutable sera généré dans le dossier `bin/` sous le nom `osinco3d.app`.

---

## 3. Configurer la simulation Taylor-Green

1. Aller dans le dossier `bin` :

```bash
cd ../bin
```

2. Modifier le fichier `parameters.o3d` :

```fortran
typesim = 2  ! Taylor-Green Vortex
```

- Vérifiez également la taille de la grille (`nx`, `ny`, `nz`) et le pas de temps (`dt`) selon vos besoins.

---

## 4. Lancer la simulation

```bash
./osinco3d.app
```

- Les résultats seront enregistrés dans `fields.bin`.

---

## 5. Visualiser les résultats

- **2D** : Utiliser Gnuplot ou Matplotlib.
- **3D** : Utiliser ParaView pour ouvrir les fichiers `.xdmf` ou `.vtu`.

---

## 6. Astuce rapide

Pour un test rapide sans trop de calcul :

```fortran
nx = 32
ny = 32
nz = 32
dt = 0.01
```

Cela permet de voir rapidement l'évolution du vortex Taylor-Green.

---
