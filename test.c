#include "test.h"

/* En fonction de la valeur de choix_test :
    1. Factorisation de n1 via Rho-Pollard (FPA0)
    2. Factorisation de n via crible quadratique
    10. Affichage base de premiers (crible quadratique)

*/

int main() {
    // CHOIX DE L'OPÉRATION --------------------------------------------------------
    int choix_test = 10;

    // CONFIGURATION - CHOIX DES PARAMETRES ----------------------------------------
    mpz_t n1;
    mpz_init_set_str(n1, "52590354472497239257283147", 10);

    /*  Exemples d'entiers à décomposer:
        n1 = 52590354472497239257283147
        n2 = 52590354464570687296135717939971 
    */

    int n = 1042387;
    int P = 50;
    int A = 500;

    // AFFICHAGE DE LA CONFIGURATION -----------------------------------------------
    printf("----------------------------------------------------------------- \n");
    gmp_printf("Test de factorisation sur :\n  - n_1 = %Zd\n", n1);
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
    if(choix_test == 1){
        // Choix du premier terme de la suite de Rho-Pollard
        mpz_t mpz_1;
        mpz_init_set_str(mpz_1, "1", 10);

        mpz_t resultat;
        factorisation_rho_pollard_sm(&mpz_1, n1, &resultat);
        gmp_printf("Factorisation terminée, voici un facteur de n_1 : %Zd", resultat);
        mpz_clears(resultat, mpz_1);
    }
    else if(choix_test == 2){
        
    }
    else if(choix_test == 10){
        liste B;
        base_de_premiers(n, P, &B);
        printf("Voici la base de premiers de n=%d du crible quadratique pour la borne P=%d : B = ", n, P);
        affichage_liste(&B);
        free(B.element);
    }

    return;
}