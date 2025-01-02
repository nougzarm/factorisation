# Factorisation de nombres entiers

Etant donné un nombre entier n (supposé positif et "assez grand"), on s'interesse aux méthodes qui permettent de trouver
deux entiers (non triviaux) a et b tels que n=ab.
La méthode naîve serait d'essayer de diviser n par i, pour i allant de 2 jusqu'à la racine de n. 
Il existe cependant des algorithmes bien plus performants pour la factorisation de nombres.

# Factorisation de Rho-Pollard 

Soit f : Z/nZ -> Z/nz une FPA, et x_0\in Z/nZ un point initial. On note la suite (x_i)_i définie par x_i = f^i(x_0).
Pour chaque k, on regarde s'il existe j<k tel que  pgcd(x_j-x_k, n) != 1,n
Si on trouve de tels k,j alors on a décomposé n.

# Factorisation via le crible quadratique

Etant données deux bornes P et A, on définit la base suivante de premiers:
    B = {p premier | p<P et jacobi(n,p)=1}
On pose ensuite l'ensemble
    S' = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
Et on prend S, le sous-ensemle de S' formé des éléments B-lisses (ie qui se décomposent dans B).

On choisit ensuite un sous-ensemble de S à |B|+1 éléments. Il se peut que le sous-ensemble
choisi ne convienne pas. Ainsi, la boucle de l'algorithme se fera sur l'ensemble des sous-ensembles 
à |B|+1 éléments de S.

Notons ce sous-ensemble : H = { (t_i)^2-n,  1 <= i <= |B|+1  }
On connait la décomposition des éléments de H et par conséquent, la décomposition des produits
possibles de ces éléments. On cherche dont un sous-ensemble F non vide de H
tel que toutes les valuations dans la base B du produit des éléments de F soient divisibles par 2.
Par conséquent, on peut extraire une racine carré de ce produit. Notons s cette racine et t le produit
des t_i tels que (t_i)^2-n est dans F. On a alors t^2 = s^2 mod n.

Pour finir, on teste si t = +- s mod n. Si ce n'est pas le cas, alors on a décomposé n et un 
de ses facteurs est pgcd(t+s, n). Sinon, on recommence avec un nouvel sous-ensemble de S

# Configuration

Nous importons la librairie GMP afin de pouvoir utiliser des nombres entiers arbitrairement grands.
Consulter le PDF pour la documentation et modifier le fichier MAKEFILE dans le cas où la librairie est installée
dans un répertoire différent.

Le MAKEFILE présent compile les fichiers pour créer un programme test permettant de tester les fonctions du fichier
source principal factorisation.c