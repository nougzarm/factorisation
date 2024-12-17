#include "factorisation.h"

/*  SOMMAIRE : 
    1. Fonctions maths intermédiaires
    2. Test de primalité
    3. Factorisation
        - Factorisation de Rho-Pollard 
        - Factorisation de Rho-Pollard optimisée
        - Factorisation via crible quadratique
 */
/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           FONCTIONS INTERMEDIAIRES                                             |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int ordre_deux(int m){
    int i = 0;
    while (m % 2 != 0){
        i++;
        m = m/2;
    }
    return i;
}

int puissance(int m, int e){
    if ( e==0 ){
        return 1;
    }
    else if (( e%2==0 )&( e!=0 )){
        return puissance(m*m, e/2);
    }
    else {
        return m*puissance(m*m, (e-1)/2);
    }
}

int modulo(int a, int b){
    return a%b + (b * (a%b<0));
}

int puissance_modulo(int m, int e, int p){
    int result = 1;
    for (int i = 1; i <= e; i++){
        result = modulo(result, p)*m;
    }
    return modulo(result, p);
}

int pgcd(int a, int b){
    int r0 = a, r1 = b, r2 = r0%r1;
    while ( r2 != 0 ){
        r0 = r1;
        r1 = r2;
        r2 = modulo(r0, r1);
    }
    return r1;
}

int jacobi(int a, int b){
    a = modulo(a, b);
    
    int signe = 1;
    int inter_reserve = 0;  // réserve pour le swap

    // Boucle de la 'réduction' via loi de réciprocité quadratique
    while ( ( b%a != 0 )&&( a != 2 ) ){
        if ( a%2 == 0 ){
            a = a/2;
            if ( (modulo(b, 8) == 3) || (modulo(b, 8) == 5) ){
                signe = signe*(-1);
            }
        }
        else {
            if ( ((a-1)*(b-1))%8 != 0 ){
            signe = signe*(-1);
            }
            inter_reserve = a;
            a = modulo(b, a);
            b = inter_reserve; 
        }
    }
    if ( a == 1 ){
        return signe;
    }
    if ( a == 2 ){
        int c = 1;
        if ( (modulo(b, 8) == 3) || (modulo(b, 8) == 5) ){
            c = -1;
        }
        return c*signe;
    }
    else{
        return 0;
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                             FONCTIONS PRIMALITE (utiles pour la crible quadratique)                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Test de Fermat : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(2^k) 
    Remarque : Nécessite l'existence de c tel que c^{n-1} != 1 [n] 
*/
int test_fermat(int n, int k){
    int b;
    for(int i=0; i<k; i++){
        b = 2 + (rand()%(n-4));
        if(pgcd(n, b) == 1){
            if(puissance_modulo(b, n-1, n) != 1){
                return 0;
            }
        }
        else{
            return 0;
        }
    }
    return 1;
}


/*  Test de Solovay-Strassen : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(2^k)
*/
int test_solovay_strassen(int n, int k){
    int b;
    for(int i=0; i<k; i++){
        b = 2 + (rand()%(n-4));
        if(pgcd(n, b) == 1){
            if(modulo(jacobi(b, n) - puissance_modulo(b, (n-1)/2, n), n) != 0){
                return 0;
            }
        }
        else{
            return 0;
        }
    }
    return 1;
}


/*  Test de Miller-Rabin : 
        - Si renvoie 0 alors n n'est pas premier
        - si renvoie 1 alors n est premier avec proba supérieur à 1-1/(4^k)
*/
int test_miller_rabin(int n, int k){
    int s = ordre_deux(n-1);
    int t = (n-1)/(puissance(2, s));
    int b; 
    int b_puissance_t;
    int b_puissance_finale;
    for(int j=0; j<k; j++){
        b = 2 + (rand()%(n-4));
        if(pgcd(n, b) == 1){
            b_puissance_t = puissance_modulo(b, t, n);
            if(b_puissance_t != -1 && b_puissance_t != 1){
                for(int i=1; i<s-1; i++){
                    b_puissance_finale = puissance_modulo(b_puissance_t, puissance(2, i), n);
                    if(b_puissance_finale == n-1){
                        break;
                    }
                    else if(b_puissance_finale == 1){
                        return 0;
                    }
                }
                b_puissance_finale = puissance_modulo(b_puissance_t, puissance(2, s-1), n);
                if(b_puissance_finale == n-1 || b_puissance_finale == 1){
                    return 0;
                }
            }
        }
        else{
            return 0;
        }
    }
    return 1;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                 FACTORISATION                                                  |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Via librairie GMP
        - choix du FPA FPA : y -> y^2 + c
*/
void FPA1(mpz_t* y, mpz_t* n, mpz_t* c){
    mpz_mul(*y, *y, *y);    // y -> y^2
    mpz_mod(*y, *y, *n);    // y -> y mod n
    mpz_add(*y, *y, *c);    // y -> y+c
    mpz_mod(*y, *y, *n);    // y -> y mod n
    gmp_printf ("Résultat : %Zd\n", *y);
    return;
}

int factorisation_rho_pollard_sm(mpz_t* premier, mpz_t* n, mpz_t* resultat){
    mpz_t mpz_1, mpz_c;
    mpz_init_set_str(mpz_1, "1", 10);
    mpz_init_set_str(mpz_c, "3", 10);   // La constante utilisée dans la FPA

    //  x_i (=premier) et x_2i (=second) de la suite engendrée par la FPA
    FPA1(premier, n, &mpz_c);
    mpz_t second;
    mpz_init(second);
    mpz_set(second, *premier);
    FPA1(&second, n, &mpz_c);

    //  Calcul de pgcd(x_2i-x_i, n)
    mpz_t potentiel_facteur;
    mpz_init(potentiel_facteur);
    mpz_sub(potentiel_facteur, second, *premier);   // potentiel_facteur = second-premier
    mpz_mod(potentiel_facteur, potentiel_facteur, *n);  // potentiel_facteur = potentiel_facteur [n]
    mpz_gcd(potentiel_facteur, *n, potentiel_facteur);  // potentiel_facteur = pgcd(potentiel_facteur, n);

    //  Compare le potentiel_facteur avec 1 et n
    int condition1 = mpz_cmp(potentiel_facteur, *n);
    int condition2 = mpz_cmp(potentiel_facteur, mpz_1);

    while(condition1*condition2 == 0){
        // On actualise les termes de la suite
        FPA1(premier, n, &mpz_c);
        FPA1(&second, n, &mpz_c);
        FPA1(&second, n, &mpz_c);
        // On recalcule le potentiel_facteur
        mpz_sub(potentiel_facteur, second, *premier);
        mpz_mod(potentiel_facteur, potentiel_facteur, *n);
        mpz_gcd(potentiel_facteur, *n, potentiel_facteur);
        // Actualisation des conditions
        condition1 = mpz_cmp(potentiel_facteur, *n);
        condition2 = mpz_cmp(potentiel_facteur, mpz_1);
    }
    mpz_set(*resultat, potentiel_facteur);   // Stockage du résultat
    //  Suppression de la mémoire allouée
    //  mpz_clears(second, difference, candidat_comparaison, candidat_facteur, mpz_1);
    return 1;
}

