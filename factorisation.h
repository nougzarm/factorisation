#ifndef TEST_PRIMALITE_INCLUDED
#define TEST_PRIMALITE_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

/*  Structure de liste, utilisée par exemple, pour:
    - Identifier un sous-ensemble (d'un ensemble S d'éléments B-lisses)
    - Matrices
    - Vecteurs du noyau des matrices considérées*/
typedef struct{
    int* element;
    int taille;
}liste;

/*  Structure de famille, utilisée pour les bases de premiers B  */
typedef struct{
    mpz_t* element;
    int taille;
}famille;

/*  Structure pour stocker la décomposition en premiers d'un entier  */
typedef struct{
    mpz_t valeur;
    int* valuation;
    famille* base_premiers;
}decomposition;

/*  Structure pour les ensembles S d'entiers B-lisses  */
typedef struct{
    int cardinal;
    mpz_t* t_element;
    decomposition* element_dec;
    famille* base_premiers;
}ensemble_b_lisse;

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                         FACTORISATION DE RHO-POLLARD                                           |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Choix du FPA : y -> y^2 + 1
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
void base_de_premiers(mpz_t* n, int P, famille* B);
/*  On définira ensuite l'ensemble
        S' = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
    Et on prend S, le sous-ensemle formé des éléments B-lisses. Remarques: 
        - S.t_element contient les indices t
        - S.element_dec contient la décomposition de t^2-n dans S.base_premiers = B  */
int valuation(mpz_t* b, mpz_t* p);
/*  Essaye de décomposer un entier b dans la base de premiers B et le stocke dans D
    - Retourne 0 si l'entier ne se décompose pas dans B (laisse D inchangé)
    - retourne 1 si l'entier se décompose (et le décompose) */
int decomposition_entier(mpz_t* b, famille* B, decomposition* D);
void ensemble_crible_quadratique(mpz_t* n, int A, famille* B, ensemble_b_lisse* S);
/*  On choisit un sous-ensemble de S à |B|+1 éléments. Il se peut que le sous-ensemble
    choisi ne convienne pas. Il nous faut ainsi une fonction permettant de 'classifier' 
    un certain nombre de ces sous-ensembles.
    Pour cela on crée une fonction choisissant une injection |B|+1 -> |S|
    Ici, on classifie ces injections suivant l'entier num_injection
    Ainsi, dans le crible quadratique on pourra parcourir ces injetions
    Remarque : ici nous considérons uniquement les injections de base qui sont les :
        x -> x + num_injection  */
void injection(int a, int b, int num_injection, liste* L);
/*  Une fois qu'on fixe un sous-ensemble noté 'sous_ensemble', on écrit sa matrice
    qui, par définitions, contient la décomposition du i-ème éléments dans sa i-ème colonne
    La matrice résultat est de dimension |B|x|sous_ensemble| */
void matrice_de_decomposition(liste* sous_ensemble, ensemble_b_lisse* S, liste* matrice);
/*  (PARTIE A OPTIMISER)
    On applique le pivot de Gauss à la matrice précédemment calculée et on extrait un élément
    non nul de son noyeau.
    Les fonctions appliquant les opérations de base nécessaires pour le pivot de Gauss:
        1. Echange de lignes : echange_lignes()
        2. Addition de lignes : plus_lignes()
        3. Recherche de la prochaine ligne à déplacer : ligne_pivotable()
        4. Premier algo permettant de se ramener à une matrice triangulaire : trigonaliser()
        5. Teste si une ligne est nulle : nullite_ligne()
        6. Teste si deux lignes peuvent être échangées sans compromettre le caractère triangulaire : teneur()
        7. Tranforme une matrice triangulaire de sorte à avoir le moins de coeff diagonaux non nuls : modif_finale()
        8. Algo final : pivot_gauss()
        9. Détermine une colonne nulle : colonne_libre()
        10. Permet/Essaye d'extraire un élément non nul du noyau : noyau()
        11. Fonction permettant de vérifier si un élément est dans le noyau : verif_noyau() */
void echange_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice);
void plus_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice);
int ligne_pivotable(int min, int nb_ligne, int nb_colonne, liste* matrice);
void trigonaliser(liste* matrice, int nb_ligne, int nb_colonne);
int nullite_ligne(liste* matrice, int nb_ligne, int nb_colonne, int ligne);
int teneur(liste* matrice, int nb_ligne, int nb_colonne, int cible, int destination);
void modif_finale(liste* matrice, int nb_ligne, int nb_colonne);
void pivot_gauss(liste* matrice, int nb_ligne, int nb_colonne);
int colonne_libre(liste* matrice, int nb_ligne, int nb_colonne);
int noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker);
int verif_noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker);

/*  Fonction finale permettant de décomposer n et de stocker son facteur dans 'resultat'.
    Voici les codes erreur possibles :
        -1 : n est premier (d'après Solovay-Strassen avec une précision k=10)
        0 : L'algorithme a réussi à décomposer n
        1 : L'algorithme n'a réussi à décomposer n avec aucun sous-ensemble
        2 : Aucun sous-ensemble possible pour les bornes A et P. Ie |S|<|B|+1
    Remarque : on décide de nouveau utiliser la librairie gmp car les produits des 
    éléments b_i B-lisses deviendront rapidement très grands  */
int crible_quadratique(mpz_t* n, int P, int A, mpz_t* resultat);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                              FONCTIONS PRIMAIRES                                               |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int modulo(int a, int b);       //  Renvoie a mod b
void affichage_liste(liste* l); 
void affichage_famille(famille* F);
void affichage_decomposition(decomposition* D); //  Affiche la décomposition D (de l'entier D->valeur)
void affichage_ensemble(ensemble_b_lisse* S);   //  Affiche l'ensemble S d'éléments lisses.
int indice_matrice(int ligne, int colonne, int nb_ligne, int nb_colonne);
void affichage_matrice(int nb_ligne, int nb_colonne, liste* M); //  Affiche la matrice M

#endif // TEST_PRIMALITE_INCLUDED