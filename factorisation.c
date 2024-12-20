#include "factorisation.h"

/*  SOMMAIRE : 
    1. Fonctions maths intermédiaires
    2. Test de primalité
    3. Factorisation
        - Factorisation de Rho-Pollard 
        - Factorisation via crible quadratique (algorithme de Dixon)
 */
/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                         1. FONCTIONS INTERMEDIAIRES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int ordre_deux(int m){
    int i = 0;
    while (m % 2 == 0){
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

long long modulo_positif(long long a, long long b) {
    if (b <= 0) {
        printf("Erreur : le diviseur doit être strictement positif.\n");
        return -1; // Retourne une erreur si b <= 0
    }

    long long mod = a % b;
    if (mod < 0) {
        mod += b;
    }
    return mod;
}

int puissance_modulo(int m, int e, int p){
    int result = 1;
    for (int i = 1; i <= e; i++){
        result = modulo(result, p)*m;
    }
    return modulo(result, p);
}

int pgcd(int a, int b){
    if(b>a){
        int c=a;
        a=b;
        b=c; 
    }
    if(b==0){
        return a;
    }
    int r0 = a, r1 = b, r2 = r0%r1;
    while ( r2 != 0 ){
        r0 = r1;
        r1 = r2;
        r2 = modulo(r0, r1);
    }
    return r1;
}

long int long_pgcd(long int a, long int b) {
    while (b != 0) {
        long int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

int partie_entiere(double n){
    int partie_ent = (int)floor(n);
    return partie_ent;
}

int racine_carree_entiere_mn(int n){
    int racine = partie_entiere(sqrt((double)n));
    return racine;
}

int jacobi(int a, int b){
    if(b%2==0){
        return 0;
    }
    if(pgcd(a, b) != 1){
        return 0;
    }
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

void affichage_liste(liste* l){
    if(l->taille == 0){
        printf("vide\n");
    }
    else{
        printf("[");
        for(int i=0; i<l->taille-1; i++){
            printf("%d, ", l->element[i]);
        }
        printf("%d]", l->element[l->taille-1]);
    }
}

void affichage_decomposition(decomposition* D){
    for(int i=0; i<D->base_premiers->taille; i++){
        if(D->valuation[i] != 0){
            printf("%d^%d*", D->base_premiers->element[i], D->valuation[i]);
        }
    }
}

void affichage_ensemble(ensemble_b_lisse* S){
    printf("{");
    for(int i=0; i<S->cardinal; i++){
        S->element_dec[i].valeur;
        printf("%d: %d = ", S->t_element[i], S->element_dec[i].valeur);
        affichage_decomposition(&S->element_dec[i]);
        printf(", ");
    }
    printf("}");
}

int indice_matrice(int ligne, int colonne, int nb_ligne, int nb_colonne){
    return colonne*nb_ligne+ligne;
}

void affichage_matrice(int nb_ligne, int nb_colonne, liste* M){
    for(int j=0; j<nb_ligne; j++){
        for(int i=0; i<nb_colonne; i++){
            printf("%d ", M->element[indice_matrice(j, i, nb_ligne, nb_colonne)]);
        }
        printf("\n");
    }
}

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                           2. FONCTIONS PRIMALITE (utiles pour le crible quadratique)                           |   
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
    |                                               3. FACTORISATION                                                 |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Via librairie GMP
        - choix du FPA : y -> y^2 + c
*/
void FPA1(mpz_t* y, mpz_t* n, mpz_t* c){
    mpz_mul(*y, *y, *y);    // y -> y^2
    mpz_mod(*y, *y, *n);    // y -> y mod n
    mpz_add(*y, *y, *c);    // y -> y+c
    mpz_mod(*y, *y, *n);    // y -> y mod n
    return;
}

int factorisation_rho_pollard_sm(mpz_t* premier, mpz_t* n, mpz_t* resultat){
    mpz_t mpz_1, mpz_c;
    mpz_init_set_str(mpz_1, "1", 10);
    mpz_init_set_str(mpz_c, "1", 10);   // La constante utilisée dans la FPA

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
    mpz_clear(potentiel_facteur);  //  Suppression de la mémoire allouée
    return 1;
}


/*  Factorisation via l'algorithme de Dixon: 
        Utilisant le crible quadratique:
            - Etant données les bornes P et A, on prendra comme base de premiers:
                B = {p premier | p<P et jacobi(n,p)=1}
              On définira ensuite l'ensemble
                S = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
              Et on cherche les entiers B-lisses de cet ensemble
                (on s'intéresse aux sous-ensembles à |B|+1 éléments)
*/
void base_de_premiers(int n, int P, liste* B){
    B->taille = 0;
    B->element = calloc(P, sizeof(int));
    int precision_test_primalite = 10;
    B->taille++;
    B->element[B->taille-1] = 2;
    for(int i = 3; i<P; i++){
        if(jacobi(n, i) == 1){
            if(test_solovay_strassen(i, precision_test_primalite) == 1){
                B->taille++;
                B->element[B->taille-1] = i;
            }
        }
    }
    B->element = realloc(B->element, B->taille*sizeof(int));
}

// Calcule v_p(b)
int valuation(int b, int p){
    int i = 0;
    while (b % p == 0){
        i++;
        b = b/p;
    }
    return i;
}

/*  Essaye de décomposer un entier b dans la base de premiers B et le stocke dans D
        - Retourne 0 si l'entier ne se décompose pas dans B
        - retourne 1 si l'entier se décompose (et le décompose) */
int decomposition_entier(int b, liste* B, decomposition* D){
    // Potentielle decomposition
    D->valuation = calloc(B->taille, sizeof(int));
    for(int i=0; i<B->taille; i++){
        D->valuation[i] = valuation(b, B->element[i]);
    }
    // On recalcule la décomposition
    int potentiel_produit=1;
    for(int i=0; i<B->taille; i++){
        potentiel_produit = potentiel_produit*(puissance(B->element[i], D->valuation[i]));
    }
    // On teste la decomposition
    if(potentiel_produit == b){
        D->valeur = b;
        D->base_premiers = B;
        return 1;
    }
    else{
        free(D->valuation);
        return 0;
    }
}

/*  Parcourt l'ensemble S' = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A} et extrait les éléments
    B-lisses pour les stocker dans S (avec leurs décompositions et les indices t) */
void ensemble_crible_quadratique(int n, int A, liste* B, ensemble_b_lisse* S){
    S->element_dec = calloc(A, sizeof(decomposition));
    S->t_element = calloc(A, sizeof(int));
    S->base_premiers = B;   // On "relie" la base de premiers (pas nécessaire dans notre contexte)
    S->cardinal = 0;    
    int racine_n = racine_carree_entiere_mn(n);
    int candidat;
    // On extrait les éléments B-lisses en parcourant S'
    for(int t = racine_n+1; t <= racine_n+A; t++){
        candidat = puissance(t, 2)-n;
        if(decomposition_entier(candidat, B, &S->element_dec[S->cardinal]) == 1){
            S->t_element[S->cardinal] = t; // On ajoute l'indice (la décompo a été ajoutée dans la condition ci-dessus)
            S->cardinal++; 
        }
    }
    // On libère la place non occupée
    S->element_dec = realloc(S->element_dec, S->cardinal*sizeof(decomposition));
    S->t_element = realloc(S->t_element, S->cardinal*sizeof(int));
}

/*  Il se peut qu'on ait besoin de trouver plusieurs sous-ensembles à |B|+1 éléments de S
    Cette fonction permet de fixer une injection a->b, représentée par la liste L
    Les injections possibles sont indexées par num_injection et ici nous considérons uniquement
    les injections de base : x -> x + num_injection */
void injection(int a, int b, int num_injection, liste* L){
    L->taille = a;
    L->element = calloc(L->taille, sizeof(int));
    for(int i=0; i<a; i++){
        L->element[i] = i+num_injection;
    }
}

/*  Une fois qu'on a choisi notre sous-ensemble à |B|+1 éléments, on construit la matrice dont les coefficients
    sont les valuations dans la base B. Ici les éléments sont représentés par les colonnes 
      -> Remarque : la matrice est représentée par une liste dans la quelle on stocke colonne par colonne */
void matrice_de_decomposition(liste* sous_ensemble, ensemble_b_lisse* S, liste* matrice){
    int nb_ligne = (S->base_premiers)->taille;  // Le nombre total de premiers
    int nb_colonne = sous_ensemble->taille;     // Le nombre d'éléments de sous-ensemble
    matrice->taille = nb_ligne*nb_colonne;
    matrice->element = calloc(matrice->taille, sizeof(int));
    for(int i=0; i<nb_colonne; i++){
        for(int j=0; j<nb_ligne; j++){
            matrice->element[i*nb_ligne+j] = S->element_dec[sous_ensemble->element[i]].valuation[j];
            matrice->element[i*nb_ligne+j] = matrice->element[i*nb_ligne+j]%2;
        }
    }
}

/*  Fonctions appliquant les opérations de base nécessaires pour le pivot de pivot Gauss:
        1. Echange de lignes
        2. Addition de lignes
        3. Recherche de la prochaine ligne à déplacer
        4. L'algorithme final de pivot de Gauss
        5. Extraction d'un élément non nul du noyau (à partir d'une matrice triangulaire sup)
        6. Fonction permettant de vérifier si un élément est dans le noyau
 */
void echange_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice){
    int* temp = calloc(nb_colonne, sizeof(int));
    // Sauvegarde de ligne1 dans temp
    for(int i=0; i<nb_colonne; i++){
        temp[i] = matrice->element[indice_matrice(ligne1, i, nb_ligne, nb_colonne)];
    }
    // Transfert de ligne2 dans ligne1
    for(int i=0; i<nb_colonne; i++){
        matrice->element[indice_matrice(ligne1, i, nb_ligne, nb_colonne)]=matrice->element[indice_matrice(ligne2, i, nb_ligne, nb_colonne)];
    }
    // Sauvegarde de temp dans ligne2
    for(int i=0; i<nb_colonne; i++){
        matrice->element[indice_matrice(ligne2, i, nb_ligne, nb_colonne)]=temp[i];
    }
    free(temp);
}

void plus_lignes(int ligne1, int ligne2, int nb_ligne, int nb_colonne, liste* matrice){
    int ind1, ind2;
    for(int i=0; i<nb_colonne; i++){
        ind1 = indice_matrice(ligne1, i, nb_ligne, nb_colonne);
        ind2 = indice_matrice(ligne2, i, nb_ligne, nb_colonne);
        matrice->element[ind1] = modulo(matrice->element[ind1]+matrice->element[ind2], 2) ;
    }
}

int ligne_pivotable(int min, int nb_ligne, int nb_colonne, liste* matrice){
    for(int j=min; j<nb_colonne; j++){
        for(int i=min; i<nb_ligne; i++){
            if(matrice->element[indice_matrice(i, j, nb_ligne, nb_colonne)]==1){
                return i;
            }
        }
    }
    return min;
}

void pivot_gauss(liste* matrice, int nb_ligne, int nb_colonne){
    int l_pivot = ligne_pivotable(0, nb_ligne, nb_colonne, matrice);
    echange_lignes(0, l_pivot, nb_ligne, nb_colonne, matrice);
    for(int k=1; k<nb_ligne; k++){
        if(matrice->element[indice_matrice(k, 0, nb_ligne, nb_colonne)] == 1){
            plus_lignes(k, 0, nb_ligne, nb_colonne, matrice);
        }
    }
    for(int min=1; min<nb_ligne; min++){
        l_pivot = ligne_pivotable(min, nb_ligne, nb_colonne, matrice);
        echange_lignes(min, l_pivot, nb_ligne, nb_colonne, matrice);
        for(int k=min+1; k<nb_ligne; k++){
            if(matrice->element[indice_matrice(k, min, nb_ligne, nb_colonne)] == 1){
                plus_lignes(k, min, nb_ligne, nb_colonne, matrice);
            }
        }
    }
}

//  Cas particulier où nb_colonne = nb_ligne+1
void noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker){
    ker->taille = nb_colonne;
    ker->element = calloc(nb_colonne, sizeof(int));
    int a = matrice->element[indice_matrice(nb_ligne-1, nb_ligne-1, nb_ligne, nb_colonne)];
    int b = matrice->element[indice_matrice(nb_ligne-1, nb_ligne, nb_ligne, nb_colonne)];
    if((a==0 && b==0)||(a==1 && b==1)){
        ker->element[nb_colonne-1] = 1;
        ker->element[nb_colonne-2] = 1;
    }
    else if(a==1 && b==0){
        ker->element[nb_colonne-1] = 1;
        ker->element[nb_colonne-2] = 0;
    }
    else if(a==0 && b==1){
        ker->element[nb_colonne-1] = 0;
        ker->element[nb_colonne-2] = 1;
    }
    int somme_partielle;
    int c;
    for(int i=nb_ligne-2; i>=0; i--){
        c = matrice->element[indice_matrice(i, i, nb_ligne, nb_colonne)];
        somme_partielle = 0;
        for(int j=nb_colonne-1; j>=c; j--){
            somme_partielle=somme_partielle+(matrice->element[indice_matrice(i, j, nb_ligne, nb_colonne)]*ker->element[j]);
            somme_partielle=modulo(somme_partielle, 2);
        }
        if((c==1 && somme_partielle==1)||(c==0 && somme_partielle==0)){
            ker->element[i] = 1;
        }
        else if(c==1 && somme_partielle==0){
            ker->element[i] = 0;
        }
    }
}

int verif_noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker){
    int* vecteur = calloc(nb_ligne, sizeof(int));
    int somme;
    for(int i=0; i<nb_ligne; i++){
        for(int j=0; j<nb_colonne; j++){
            somme=somme+(matrice->element[indice_matrice(i, j, nb_ligne, nb_colonne)]*ker->element[j]);
        }
        vecteur[i] = modulo(somme, 2);
        somme=0;
    }
    int resultat = 1;
    for(int k=0; k<nb_ligne; k++){
        if(vecteur[k]==1){
            resultat=0;
        }
    }
    free(vecteur);
    return resultat;
}

/*  Fonction finale permettant de décomposer n et de stocker son facteur dans 'resultat'.
    Voici les codes erreur possibles :
        -1 : n est premier (d'après Solovay-Strassen avec une précision k=10)
        0 : L'algorithme a réussi à décomposer n
        1 : L'algorithme n'a réussi à décomposer n avec aucun sous-ensemble
        2 : Aucun sous-ensemble possible pour les bornes A et P. Ie |S|<|B|+1
    Remarque : on décide de nouveau utiliser la librairie gmp car les produits des 
    éléments b_i B-lisses deviendront rapidement très grands  */
int crible_quadratique(int n, int P, int A, mpz_t* resultat){
    if(test_solovay_strassen(n, 10)==1){
        return -1;  // n est premier
    }
    // Liste B contenant la base de premiers
    liste B;
    base_de_premiers(n, P, &B);
    // L'ensemble des éléments B-lisses
    ensemble_b_lisse S;
    ensemble_crible_quadratique(n, A, &B, &S);

    // Affichage de B et S
    printf("  - On a comme base de premiers B = {p premier | p<P et jacobi(n,p)=1} = ");
    base_de_premiers(n, P, &B);
    affichage_liste(&B);
    printf("\n  - Puis, on extrait l'ensemble suivant des éléments B-lisses");
    printf("\n    S = {t^2-n B-lisse | sqrt(n)+1 =< t =< sqrt(n)+A} \n      = ");
    affichage_ensemble(&S);
    printf("\n");

    /*  On vérifie si il est possible de prendre |B|+1 éléments distincts dans S
        si oui on passe à la suite */
    int nb_ss_ens = S.cardinal-(B.taille+1)+1;
    if(nb_ss_ens < 1){
        return 2; // S possède moins de |B|+1 éléments. 
    }

    // Une boucle pour chaque sous-ensemble de S à |B|+1 éléments
    liste ss_ens_S; // le sous-ensemble actuel
    liste matrice;  // la matrice des valuations correspondant au sous-ensemble actuel
    liste ker;      // un éléments non nul du noyau de la matrice actuelle

    // Les constantes qu'on comparera. On regardera si t = s [n] ou si t = -s [n]
    mpz_t t, s;
    mpz_init_set_str(t, "1", 10);
    mpz_init_set_str(s, "1", 10);
    int* valuations_s;
    mpz_t t_reduit, s_reduit, _s_reduit;
    mpz_inits(t_reduit, s_reduit, _s_reduit);

    // Utile pour les opérations dans la librairie gmp (en fin de boucle)
    mpz_t mpz_n;
    mpz_init_set_si(mpz_n, n);

    for(int i=0; i<nb_ss_ens; i++){
        printf("\nBoucle numéro %d.\n", i);
        injection(B.taille+1, S.cardinal, i, &ss_ens_S);    // On fixe le i-ème sous-ensemble

        // On calcule la matrice des valuations pour ce sous-ensemble
        matrice_de_decomposition(&ss_ens_S, &S, &matrice);

        // On trigonalise la matrice
        pivot_gauss(&matrice, B.taille, B.taille+1);
        /*  On extrait un élément non nul du noyau 
            (qui représentera les éléments du sous-ensemble à considérer dans les produits t et s) */
        liste ker;
        noyau(&matrice, B.taille, B.taille+1, &ker);

        affichage_liste(&ker);

        // Calcule du produit des b_i sélectionnés
        for(int i=0; i<B.taille+1; i++){
            if(ker.element[i] == 1){
                mpz_mul_si(t, t, S.t_element[ss_ens_S.element[i]]);
            }
        }
        printf("\n voici le produit t : %Zd", t);
        
        valuations_s = calloc(B.taille, sizeof(int));
        // Calcul des exposants de s
        for(int j=0; j<B.taille; j++){
            valuations_s[j] = 0;
            for(int i=0; i<B.taille+1; i++){
                if(ker.element[i] == 1){
                    valuations_s[j] = valuations_s[j] + S.element_dec[ss_ens_S.element[i]].valuation[j];
                }
            }
            valuations_s[j] = valuations_s[j]/2;
        }
        
        // Calcul de s
        for(int j=0; j<B.taille; j++){
            mpz_mul_si(s, s, puissance(B.element[j], valuations_s[j]));
        }
        printf("\n voici le produit s : %Zd", s);

        // Les éléments de comparaison
        mpz_neg(_s_reduit, s);
        mpz_mod(t_reduit, t, mpz_n);
        mpz_mod(s_reduit, s, mpz_n);
        mpz_mod(_s_reduit, _s_reduit, mpz_n);
        
        gmp_printf("\n t mod n : %Zd", t_reduit);
        gmp_printf("\n s mod n : %Zd", s_reduit);
        gmp_printf("\n -s mod n : %Zd", _s_reduit);

        gmp_printf("\n n : %Zd", mpz_n);
        
        /*  Si t = s [n] ou t = -s [n], alors on n'obtient pas de facteur non trivial.
            On refait tout avec un nouvel sous-ensemble  */
        if( mpz_cmp(t_reduit, s_reduit) == 0 || mpz_cmp(t_reduit, _s_reduit) == 0 ){
            printf("échoué");
            free(ss_ens_S.element);
            free(matrice.element);
            free(ker.element);
            free(valuations_s);
            mpz_set_str(t, "1", 10);
            mpz_set_str(s, "1", 10);
        }
        /*  Si t = +-s [n], on a alors un facteur non trivial qui est pgcd(t+s, n)
            On le stocke dans 'resultat' et l'algorithme se termine  */
        else{
            printf("réussi");
            mpz_t somme_t_s;
            mpz_init(somme_t_s);
            mpz_add (somme_t_s, t, s);
            mpz_gcd(*resultat, somme_t_s, mpz_n);
            return 0;
        }
    }
    mpz_clears(mpz_n, t, s, t_reduit, s_reduit, _s_reduit);
    return 1;
}