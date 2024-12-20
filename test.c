#include "test.h"

/* En fonction de la valeur de choix_test :
    -1. Tests de fonctions intermédiaires.
    1. Factorisation de n1 via Rho-Pollard (FPA0)
    2. Factorisation de n via le crible quadratique
*/

int main() {
    // CHOIX DE L'OPÉRATION --------------------------------------------------------
    int choix_test = 1;

    // CONFIGURATION - CHOIX DES PARAMETRES ----------------------------------------
    mpz_t n1;
    mpz_init_set_str(n1, "52590354472497239257283147", 10);

    /*  Exemples d'entiers à décomposer:
        n1 = 52590354472497239257283147
        n2 = 52590354464570687296135717939971 
        n = 1042387
    */

    int n = 1042387;
    int P = 50;
    int A = 500;

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
    if(choix_test == 1){
        // Choix du premier terme de la suite de Rho-Pollard
        mpz_t mpz_1;
        mpz_init_set_str(mpz_1, "1", 10);
        mpz_t resultat;
        mpz_init(resultat);
        factorisation_rho_pollard_sm(&mpz_1, n1, &resultat);
        gmp_printf("Factorisation de n_1 via Rho-Pollard; voici un facteur de n_1 : %Zd", resultat);
        mpz_clears(resultat, mpz_1);
    }

    /*  Décomposition de n via le crible quadratique  */
    else if(choix_test == 2){
        printf("Factorisation de n = %d, avec les bornes P = %d et A = %d :\n", n, P, A);
        mpz_t resultat;
        mpz_init(resultat);
        int reussite = crible_quadratique(n, P, A, &resultat);
        if(reussite == 1){
            printf("\nL'algorithme n'a pas réussi à factoriser n");
        }
        else if(reussite == 0){
            gmp_printf("\nL'algorithme a réussi à factoriser n, voici un facteur trouvé : %Zd", resultat);
        }
        else if(reussite == -1){
            printf("\nD'après Solovay Strassen, n est premier");
        }
        else if(reussite == 2){
            printf("\nErreur %d: Aucun sous-ensemble possible pour les bornes A et P. Ie |S|<|B|+1", reussite);
            printf("\nChoisir de nouvelles valeurs de P et/ou A");
        }
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