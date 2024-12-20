#include "test.h"

/* En fonction de la valeur de choix_test :
    -1. Tests de fonctions intermédiaires.
    1. Factorisation de n1 via Rho-Pollard (FPA0)
    2. Factorisation de n via crible quadratique
    3. Factorisation de n via l'algo naïf

*/

int main() {
    // CHOIX DE L'OPÉRATION --------------------------------------------------------
    int choix_test = 2;

    // CONFIGURATION - CHOIX DES PARAMETRES ----------------------------------------
    mpz_t n1;
    mpz_init_set_str(n1, "52590354472497239257283147", 10);

    /*  Exemples d'entiers à décomposer:
        n1 = 52590354472497239257283147
        n2 = 52590354464570687296135717939971 
        n = 1042387
    */

    int n = 10;
    int P = 50;
    int A = 100;

    // AFFICHAGE DE LA CONFIGURATION -----------------------------------------------
    printf("----------------------------------------------------------------- \n");
    printf("Définitions :\n");
    gmp_printf("  - n_1 = %Zd\n", n1);
    printf("  - n = %d, P = %d, A = %d\n", n, P, A);
    printf("\n");
    
    // LANCEMENT DE L'ALGORITHME ---------------------------------------------------
    test(choix_test, &n1, n, P, A);
    printf("\n----------------------------------------------------------------- \n");

    // SUPPRESSION DE LA MEMOIRE ALLOUEE
    mpz_clear(n1);
    return 1;
}


// Définition des tests
void test(int choix_test, mpz_t* n1, int n, int P, int A){
    /*  Décomposition de n1 via l'algorithme de Rho-Pollard  */
    if(choix_test ==1){
        // Choix du premier terme de la suite de Rho-Pollard
        mpz_t mpz_1;
        mpz_init_set_str(mpz_1, "1", 10);

        mpz_t resultat;
        factorisation_rho_pollard_sm(&mpz_1, n1, &resultat);
        gmp_printf("Factorisation de n_1 via Rho-Pollard; voici un facteur de n_1 : %Zd", resultat);
        mpz_clears(resultat, mpz_1);
    }

    /*  Décomposition de n via le crible quadratique  */
    else if(choix_test == 2){
        printf("Factorisation de n = %d, avec les bornes P = %d et A = %d :\n", n, P, A);
        printf("  - On a comme base de premiers B = {p premier | p<P et jacobi(n,p)=1} = ");
        liste B;
        base_de_premiers(n, P, &B);
        affichage_liste(&B); 
        printf("\n  - Puis, on extrait l'ensemble suivant des éléments B-lisses");
        printf("\n    S = {t^2-n B-lisse | sqrt(n)+1 =< t =< sqrt(n)+A} \n      = ");
        ensemble_b_lisse S;
        ensemble_crible_quadratique(n, A, &B, &S);
        affichage_ensemble(&S);
        mpz_t resultat;
        mpz_init(resultat);
        int reussite = crible_quadratique(n, P, A, &resultat);
        if(reussite == 0){
            printf("\nL'algorithme n'a pas réussi à factoriser n");
        }
        else{
            gmp_printf("\nL'algorithme a réussi à factoriser n, voici un facteur trouvé : %Zd", resultat);
        };
    }

    else if(choix_test == 3){
        int facteur = factorisation_naif(n);
        printf("Factorisation de n= via l'algo naïf; voici un facteur de n : %d", facteur);
    }
    else if(choix_test == -1){
        liste B;
        base_de_premiers(n, P, &B);
        printf("Voici la base de premiers de n=%d du crible quadratique pour la borne P=%d : B = ", n, P);
        affichage_liste(&B);
        int m=3*11*11*23*43;
        decomposition dec_m;
        int res = decomposition_entier(m, &B, &dec_m);
        printf("\nessayons de décomposer m=%d dans cette base : résultat = %d\n", m, res);
        printf("\nVoici sa décomposition : m = %d = ", m);
        affichage_decomposition(&dec_m);
        // free(dec_m.valuation);
        free(B.element);
    }

    return;
}