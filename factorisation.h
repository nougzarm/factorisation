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

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           FONCTIONS FACTORISATION                                              |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Via librairie GMP
        - avec comme FPA : y -> y^2 + 1
*/
void FPA1(mpz_t* y, mpz_t* n, mpz_t* c);   // Fonction pseudo-aléatoire FPA1 : y -> y^2 + c [n]
int factorisation_rho_pollard_sm(mpz_t* y_0, mpz_t* n, mpz_t* resultat);


/*  Factorisation via l'algorithme de Dixon: 
        - Utilisant le crible quadratique:
          Etant données les bornes P et A, on prendra comme base de premiers:
            B = {p premier | p<P et jacobi(n,p)=1}
*/
void base_de_premiers(int n, int P, liste* B);
/*  On définira ensuite l'ensemble
        S = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
    Et on cherche les entiers B-lisses de cet ensemble
    (on essaye de trouver |B|+1 tels éléments) */
int valuation(int b, int p);
int decomposition_entier(int b, liste* base_premiers, decomposition* D);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                             FONCTIONS PRIMALITE (utiles pour le crible quadratique)                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Test de Fermat : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(2^k) 
    Remarque : Nécessite l'existence de c tel que c^{n-1} != 1 [n] 
*/
int test_fermat(int n, int k);

/*  Test de Solovay-Strassen : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(2^k)
*/
int test_solovay_strassen(int n, int k);

/*  Test de Miller-Rabin : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(4^k)
*/
int test_miller_rabin(int n, int k);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                              FONCTIONS PRIMAIRES                                               |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int ordre_deux(int m);          //  La 2-valuation de m (peut être remplacé par valuation(m, 2))
int puissance(int m, int e);    //  Renvoie m^e
int modulo(int a, int b);                       //  Renvoie a mod b
int puissance_modulo(int m, int e, int p);      //  Renvoie m^e mod p
int pgcd(int a, int b);                 //  PGCD
int partie_entiere(double n);           //  Partie entière de n
int racine_carree_entiere_mn(int n);    //  Racine carrée entière de n (algo naïf)
int jacobi(int a, int b);       //  Symbole de Jacobi
void affichage_liste(liste* l); //  Permet l'affichage d'une liste

#endif // TEST_PRIMALITE_INCLUDED