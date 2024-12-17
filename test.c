#include "test.h"

/* En fonction de la valeur de choix_test :
    1. Factorisation de n1 et n2 via Rho-Pollard (FPA0)
    2. Factorisation de n via crible quadratique

*/

int main() {
    // CHOIX DE L'OPÉRATION --------------------------------------------------------

    /* Initialisation des entiers:  int mpz_init_set_str (mpz t rop, const char *str, int base) [Function]

    Initialize rop and set its value like mpz_set_str (see its documentation above for details).
    If the string is a correct base base number, the function returns 0; if an error occurs it returns
    −1. rop is initialized even if an error occurs. (I.e., you have to call mpz_clear for it.)   */

    int choix_test = 3;

    // CONFIGURATION - CHOIX DES PARAMETRES ----------------------------------------
    mpz_t n1;
    mpz_t n2;
    mpz_init_set_str(n1, "52590354472497239257283147", 10);
    mpz_init_set_str(n2, "52590354464570687296135717939971", 10);

    int n = 1042387;
    int P = 50;
    int A = 500;

    // AFFICHAGE DE LA CONFIGURATION -----------------------------------------------
    printf("----------------------------------------------------------------- \n");
    gmp_printf("Test de factorisation sur :\n  - n_1 = %Zd\n", n1);
    gmp_printf("  - n_2 = %Zd", n2); printf("\n\n");
    
    // LANCEMENT DE L'ALGORITHME ---------------------------------------------------
    test(choix_test, &n1, &n2, n, P, A);
    printf("\n----------------------------------------------------------------- \n");

    mpz_clear(n1);
    mpz_clear(n2);
    return 1;
}


// Définition des tests
void test(int choix_test, mpz_t* n1, mpz_t* n2, int n, int P, int A){
    if(choix_test == 1){
        // int resultat = factorisation_rho 
    }

    else if(choix_test == 2){
        
    }

    return;
}