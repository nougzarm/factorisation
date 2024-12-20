#ifndef TEST_PRIMALITE_INCLUDED
#define TEST_PRIMALITE_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

//  Structure de liste (pour stocker une base de premiers par exemple)
typedef struct{
    int* element;
    int taille;
}liste;

//  Structure permettant de représenter la décomposition d'un entier
typedef struct{
    int valeur;
    int* valuation;
    liste* base_premiers;
}decomposition;

//  Structure d'ensemble d'entiers B-lisses (où B est une base de premiers)
typedef struct{
    int cardinal;
    int* t_element;
    decomposition* element_dec;
    liste* base_premiers;
}ensemble_b_lisse;

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                         FACTORISATION DE RHO-POLLARD                                           |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Via librairie GMP
        - avec comme FPA : y -> y^2 + 1
*/
void FPA1(mpz_t* y, mpz_t* n, mpz_t* c);   // Fonction pseudo-aléatoire FPA1 : y -> y^2 + c [n]
int factorisation_rho_pollard_sm(mpz_t* y_0, mpz_t* n, mpz_t* resultat);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                   FACTORISATION VIA LE CRIBLE QUADRATIQUE                                      |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation via le crible quadratique: 
        Etant données les bornes P et A, on prendra comme base de premiers:
        B = {p premier | p<P et jacobi(n,p)=1}
*/
void base_de_premiers(int n, int P, liste* B);
/*  On définira ensuite l'ensemble
        S' = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
    Et on prend S, le sous-ensemle formé des éléments B-lisses. Remarques: 
        - S.t_element contient les indices t
        - S.element_dec contient la décomposition de t^2-n dans S.base_premiers = B  */
int valuation(int b, int p);
int decomposition_entier(int b, liste* base_premiers, decomposition* D);
void ensemble_crible_quadratique(int n, int A, liste* B, ensemble_b_lisse* S);
/*  On choisit un sous-ensemble de S à |B|+1 éléments. Il se peut que le sous-ensemble
    choisi ne convienne pas. Il nous faut ainsi une fonction permettant de 'classifier' 
    un certain nombre de ces sous-ensembles.
    Pour cela on crée une fonction choisissant une injection |B|+1 dans |S|
    Ici, on classifie ces injections suivant num_injection
    Ainsi, dans le crible quadratique on pourra parcourir ces injetions
    Remarque : ici nous considérons uniquement les injections de base qui sont les :
        x -> x + num_injection  */
void injection(int a, int b, int num_injection, liste* L);
/*  Une fois qu'on fixe un sous-ensemble noté sous_ensemble, on écris sa matrice
    qui, par définitions, contient la décomposition du i-ème éléments dans sa i-ème colonne
    La matrice résultat est de dimension |B|x|sous_ensemble| */
void matrice_de_decomposition(liste* sous_ensemble, ensemble_b_lisse* S, liste* matrice);
/*  On applique le pivot de Gauss à la matrice précédemment calculée et on extrait un élément
    non nul de son noyeau.
    Les fonctions appliquant les opérations de base nécessaires pour le pivot de pivot Gauss:
        1. Echange de lignes
        2. Addition de lignes
        3. Recherche de la prochaine ligne à déplacer
        4. L'algorithme final de pivot de Gauss
        5. Extraction d'un élément non nul du noyau (à partir d'une matrice triangulaire sup)
        6. Fonction permettant de vérifier si un élément est dans le noyau  */
void echange_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice);
void plus_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice);
int ligne_pivotable(int min, int nb_ligne, int nb_colonne, liste* matrice);
void pivot_gauss(liste* matrice, int nb_ligne, int nb_colonne);
void noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker);
int verif_noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker);

/*  Fonction finale permettant de décomposer n et de stocker son facteur dans 'resultat'.
    Voici les codes erreur possibles :
        -1 : n est premier (d'après Solovay-Strassen avec une précision k=10)
        0 : L'algorithme a réussi à décomposer n
        1 : L'algorithme n'a réussi à décomposer n avec aucun sous-ensemble  */
int crible_quadratique(int n, int P, int A, mpz_t* resultat);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                            FONCTION PRIMALITE (nécessaire pour le crible quadratique)                          |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Test de Solovay-Strassen : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(2^k)
*/
int test_solovay_strassen(int n, int k);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                              FONCTIONS PRIMAIRES                                               |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int puissance(int m, int e);    //  Renvoie m^e
int modulo(int a, int b);                       //  Renvoie a mod b
int puissance_modulo(int m, int e, int p);      //  Renvoie m^e mod p
int pgcd(int a, int b);                 //  PGCD
int partie_entiere(double n);           //  Partie entière de n
int racine_carree_entiere_mn(int n);    //  Racine carrée entière de n (algo naïf)
int jacobi(int a, int b);       //  Symbole de Jacobi
void affichage_liste(liste* l); //  Permet l'affichage d'une liste
void affichage_decomposition(decomposition* D); //  Affiche la décomposition D (de l'entier D->valeur)
void affichage_ensemble(ensemble_b_lisse* S);   //  Affiche l'ensemble S d'éléments lisses.
int indice_matrice(int ligne, int colonne, int nb_ligne, int nb_colonne);
void affichage_matrice(int nb_ligne, int nb_colonne, liste* M); //  Affiche la matrice M

#endif // TEST_PRIMALITE_INCLUDED