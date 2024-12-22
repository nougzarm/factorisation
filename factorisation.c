#include "factorisation.h"

/*  SOMMAIRE : 
    1. Fonctions maths intermédiaires
    2. Factorisation de Rho-Pollard 
    3. Factorisation via le crible quadratique
 */

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                         1. FONCTIONS INTERMEDIAIRES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */

int modulo(int a, int b){
    return a%b + (b * (a%b<0));
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

void affichage_famille(famille* F){
    if(F->taille == 0){
        printf("la famille est vide\n");
    }
    else{
        printf("[");
        for(int i=0; i<F->taille-1; i++){
            gmp_printf("%Zd, ", F->element[i]);
        }
        gmp_printf("%Zd]", F->element[F->taille-1]);
    }
}

void affichage_decomposition(decomposition* D){
    int max=0;
    for(int j=0; j<D->base_premiers->taille; j++){
        if(D->valuation[j] != 0){
            max = j;
        }
    }
    for(int i=0; i<max; i++){
        if(D->valuation[i] > 1){
            gmp_printf("(%Zd^%d)*", D->base_premiers->element[i], D->valuation[i]);
        }
        else if(D->valuation[i] == 1){
            gmp_printf("%Zd*", D->base_premiers->element[i]);
        }
    }
    if(D->valuation[max] > 1){
        gmp_printf("(%Zd^%d)", D->base_premiers->element[max], D->valuation[max]);
    }
    else if(D->valuation[max] == 1){
        gmp_printf("%Zd", D->base_premiers->element[max]);
    }
}

void affichage_ensemble(ensemble_b_lisse* S){
    if(S->cardinal==0){
        printf("{ }");
        return;
    }
    printf("{");
    for(int i=0; i< S->cardinal-1; i++){
        gmp_printf("(%Zd)^2-n = %Zd=", S->t_element[i], S->element_dec[i].valeur);
        affichage_decomposition(&S->element_dec[i]);
        printf(", ");
    }
    gmp_printf("(%Zd)^2-n = %Zd=", S->t_element[S->cardinal-1], S->element_dec[S->cardinal-1].valeur);
    affichage_decomposition(&S->element_dec[S->cardinal-1]);
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

void vider_base_ensemble(famille* B, ensemble_b_lisse* S){
    if(B->taille>0){
        for(int i=0; i<B->taille; i++){
            mpz_clear(B->element[i]);
        }
        free(B->element);
    }
    if(S->cardinal>0){
        for(int i=0; i<S->cardinal; i++){
            mpz_clear(S->t_element[i]);
            
            mpz_clear(S->element_dec[i].valeur);
            free(S->element_dec[i].valuation);
            // free(S->element_dec[i].base_premiers);
        }
        free(S->t_element);
        free(S->element_dec);
        /* On n'efface pas S->base_premiers car ça pointe à la même adresse que B */
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                      2. FACTORISATION DE RHO-POLLARD                                           |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation de Rho-Pollard: 
        - Choix du FPA : y -> y^2 + c
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
    mpz_set(*resultat, potentiel_facteur);  // Stockage du résultat
    /*  Suppression de la mémoire  */
    mpz_clear(potentiel_facteur);
    mpz_clear(second); mpz_clear(mpz_1); mpz_clear(mpz_c);
    return 1;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                3. FACTORISATION VIA LE CRIBLE QUADRATIQUE                                      |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
/*  Factorisation via le crible quadratique: 
        Etant données les bornes P et A, on prendra comme base de premiers:
        B = {p premier | p<P et jacobi(n,p)=1}
*/
void base_de_premiers(mpz_t* n, int P, famille* B){
    B->taille = 0;
    B->element = calloc(P, sizeof(mpz_t));
    B->taille++;
    mpz_init_set_str(B->element[B->taille-1], "2", 10);
    mpz_t mpz_indice;
    mpz_init(mpz_indice);
    for(int i = 3; i<P; i=i+2){
        mpz_set_si(mpz_indice, i);
        if(mpz_jacobi(*n, mpz_indice) == 1){
            if(mpz_probab_prime_p(mpz_indice, 10) > 0){
                B->taille++;
                mpz_init_set(B->element[B->taille-1], mpz_indice);
            }
        }
    }
    mpz_clear(mpz_indice);
    B->element = realloc(B->element, B->taille*sizeof(mpz_t));
}

int valuation(mpz_t* b, mpz_t* p){
    int i = 0;
    mpz_t b_copie;
    mpz_init_set(b_copie, *b);
    while (mpz_divisible_p(b_copie, *p) == 1){
        i++;
        mpz_cdiv_q(b_copie, b_copie, *p);
    }
    return i;
}

/*  Essaye de décomposer un entier b dans la base de premiers B et le stocke dans D
        - Retourne 0 si l'entier ne se décompose pas dans B (laisse D inchangé)
        - retourne 1 si l'entier se décompose (et le décompose) */
int decomposition_entier(mpz_t* b, famille* B, decomposition* D){
    // Potentielle decomposition
    D->valuation = calloc(B->taille, sizeof(int));
    for(int i=0; i<B->taille; i++){
        D->valuation[i] = valuation(b, &B->element[i]);
    }
    // On recalcule la décomposition
    mpz_t potentiel_produit;
    mpz_init_set_str(potentiel_produit, "1", 10);
    mpz_t facteur;
    mpz_init(facteur);
    for(int i=0; i<B->taille; i++){
        mpz_pow_ui(facteur, B->element[i], D->valuation[i]);    // facteur = p^vp
        mpz_mul(potentiel_produit, potentiel_produit, facteur); // potentiel_produit <- potentiel_produit*facteur
    }
    mpz_clear(facteur);
    // On teste la decomposition
    if(mpz_cmp(potentiel_produit, *b)==0){
        mpz_init_set(D->valeur, *b);
        D->base_premiers = B;
        mpz_clear(potentiel_produit);
        return 1;
    }
    else{
        free(D->valuation);
        mpz_clear(potentiel_produit);
        return 0;
    }
}

/*  On définira ensuite l'ensemble
        S' = {t^2-n | sqrt(n)+1 =< t =< sqrt(n)+A}
    Et on prend S, le sous-ensemle formé des éléments B-lisses. Remarques: 
        - S.t_element contient les indices t
        - S.element_dec contient la décomposition de t^2-n dans S.base_premiers = B  */
void ensemble_crible_quadratique(mpz_t* n, int A, famille* B, ensemble_b_lisse* S){
    S->element_dec = calloc(A, sizeof(decomposition));
    S->t_element = calloc(A, sizeof(mpz_t));
    S->base_premiers = B;   // On "relie" la base de premiers (pas nécessaire dans notre contexte)
    S->cardinal = 0;
    
    mpz_t racine_n;
    mpz_t candidat;
    mpz_init(racine_n);
    mpz_init(candidat);
    mpz_sqrt(racine_n, *n);
    // On extrait les éléments B-lisses en parcourant S'
    for(int t = 1; t <= A; t++){
        // candidat = (t+srt(n))^2-n
        mpz_add_ui(candidat, racine_n, t);
        mpz_pow_ui(candidat, candidat, 2);
        mpz_sub(candidat, candidat, *n);
        if(decomposition_entier(&candidat, B, &S->element_dec[S->cardinal]) == 1){
            mpz_init(S->t_element[S->cardinal]);
            mpz_add_ui(S->t_element[S->cardinal] , racine_n, t);  // On ajoute l'indice
            S->cardinal++; 
        }
    }
    /*  Suppression de la mémoire  */
    mpz_clear(racine_n); mpz_clear(candidat);
    S->element_dec = realloc(S->element_dec, S->cardinal*sizeof(decomposition));
    S->t_element = realloc(S->t_element, S->cardinal*sizeof(mpz_t));
}

/*  On choisit un sous-ensemble de S à |B|+1 éléments. Il se peut que le sous-ensemble
    choisi ne convienne pas. Il nous faut ainsi une fonction permettant de 'classifier' 
    un certain nombre de ces sous-ensembles.
    Pour cela on crée une fonction choisissant une injection |B|+1 -> |S|
    Ici, on classifie ces injections suivant l'entier num_injection
    Ainsi, dans le crible quadratique on pourra parcourir ces injetions
    Remarque : ici nous considérons uniquement les injections de base qui sont les :
        x -> x + num_injection  */
void injection(int a, int b, int num_injection, liste* L){
    L->taille = a;
    L->element = calloc(L->taille, sizeof(int));
    for(int i=0; i<a; i++){
        L->element[i] = i+num_injection;
    }
}

/*  Une fois qu'on fixe un sous-ensemble noté 'sous_ensemble', on écrit sa matrice
    qui, par définitions, contient la décomposition du i-ème éléments dans sa i-ème colonne
    La matrice résultat est de dimension |B|x|sous_ensemble| */
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

void trigonaliser(liste* matrice, int nb_ligne, int nb_colonne){
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

int nullite_ligne(liste* matrice, int nb_ligne, int nb_colonne, int ligne){
    for(int i=0; i<nb_colonne; i++){
        if(matrice->element[indice_matrice(ligne, i, nb_ligne, nb_colonne)] == 1){
            return 0;
        }
    }
    return 1;
}

int teneur(liste* matrice, int nb_ligne, int nb_colonne, int cible, int destination){
    for(int i=0; i<destination; i++){
        if(matrice->element[indice_matrice(cible, i, nb_ligne, nb_colonne)] == 1){
            return 0;
        }
    }
    return 1;
}

void modif_finale(liste* matrice, int nb_ligne, int nb_colonne){
    for(int i=nb_ligne-1; i>=0; i--){
        for(int j=0; j<i-1; j++){
            if(teneur(matrice, nb_ligne, nb_colonne, j, i)==1 && nullite_ligne(matrice, nb_ligne, nb_colonne, j)==0){
                echange_lignes(i, j, nb_ligne, nb_colonne, matrice);
                break;
            }
        }
    }
}

void pivot_gauss(liste* matrice, int nb_ligne, int nb_colonne){
    trigonaliser(matrice, nb_ligne, nb_colonne);
    modif_finale(matrice, nb_ligne, nb_colonne);
}

int colonne_libre(liste* matrice, int nb_ligne, int nb_colonne){
    int c_libre = -1;
    for(int j=nb_colonne-1; j>=0; j--){
        int est_libre = 1;
        for(int ligne=0; ligne<nb_ligne; ligne++) {
            if (matrice->element[indice_matrice(ligne, j, nb_ligne, nb_colonne)] != 0) {
                est_libre = 0;
                break;
            }
        }
        if(est_libre){
            c_libre = j;
            break;
        }
    }
    return c_libre;
}

//  Cas particulier où nb_colonne = nb_ligne+1
int noyau(liste* matrice, int nb_ligne, int nb_colonne, liste* ker){
    ker->taille = nb_colonne;
    ker->element = calloc(nb_colonne, sizeof(int));
    // Méthode 1
    int c_libre = colonne_libre(matrice, nb_ligne, nb_colonne);
    if(c_libre != -1){
        ker->element[c_libre] = 1;
        for(int ligne=nb_ligne-1; ligne>=0; ligne--){
            int somme = 0;
            for(int colonne=ligne+1; colonne<nb_colonne; colonne++){
                somme = somme + matrice->element[indice_matrice(ligne, colonne, nb_ligne, nb_colonne)]*ker->element[colonne];
            }
            somme = modulo(somme, 2);
            if(matrice->element[indice_matrice(ligne, ligne, nb_ligne, nb_colonne)] != 0){
                ker->element[ligne] = somme;
            }
            else{
                continue;
            }
        }
        return 1;   // Noyau trouvé
    }
    // Méthode 2
    else{
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
        if(verif_noyau(matrice, nb_ligne, nb_colonne, ker) == 0){
            ker->taille=0;
            free(ker->element);
            return 0;
        }
        return 1;
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
int crible_quadratique(mpz_t* n, int P, int A, mpz_t* resultat){
    // On commence par tester si n est premier
    if(mpz_probab_prime_p(*n, 10)>0){
        return -1;  // n est premier
    }
    // Calcul de la base B de premiers
    famille B;
    base_de_premiers(n, P, &B);
    int b = B.taille;

    // Calcule de l'ensemble S des éléments B-lisses
    ensemble_b_lisse S;
    ensemble_crible_quadratique(n, A, &B, &S);

/*  |----------------------------------------------------------------------------------------------------------------|
    |                                                  AFFICHAGES                                                    |   
    |----------------------------------------------------------------------------------------------------------------| */
    // Affichage de B et S
    printf("  - On a comme base de premiers B = {p premier | p<P et jacobi(n,p)=1} = ");
    base_de_premiers(n, P, &B);
    affichage_famille(&B);
    printf("\n  - Puis, on extrait l'ensemble suivant des éléments B-lisses :");
    printf("\n    S = {t^2-n B-lisse | sqrt(n)+1 =< t =< sqrt(n)+A} \n      = ");
    affichage_ensemble(&S);
    printf("\n");
/*  |----------------------------------------------------------------------------------------------------------------|   
    |----------------------------------------------------------------------------------------------------------------| */

    /*  On vérifie si il est possible de prendre |B|+1 éléments distincts dans S
        si oui on passe à la suite */
    int nb_ss_ens = S.cardinal-(b+1)+1;
    if(nb_ss_ens < 1){
        return 2; // S possède moins de |B|+1 éléments. 
    }

    // Une boucle pour chaque sous-ensemble de S à |B|+1 éléments
    liste ss_ens_S; // le sous-ensemble actuel
    liste matrice;  // la matrice des valuations correspondant au sous-ensemble actuel
    liste ker;      // un éléments non nul du noyau de 'matrice'

    // Les constantes qu'on comparera. On regardera si t = s [n] ou si t = -s [n]
    mpz_t t, s, facteur_s;
    mpz_init_set_str(t, "1", 10);
    mpz_init_set_str(s, "1", 10);
    mpz_init(facteur_s);
    int* valuations_s;
    mpz_t t_reduit, s_reduit, _s_reduit;
    mpz_inits(t_reduit, s_reduit, _s_reduit);

    for(int i=0; i<nb_ss_ens; i++){
        /* On fixe le i-ème sous-ensemble */
        injection(b+1, S.cardinal, i, &ss_ens_S); 

        /* On calcule la matrice des valuations pour ce sous-ensemble */
        matrice_de_decomposition(&ss_ens_S, &S, &matrice);

        /* On trigonalise la matrice */
        pivot_gauss(&matrice, b, b+1);

/*  |----------------------------------------------------------------------------------------------------------------|
    |                                                  AFFICHAGES                                                    |   
    |----------------------------------------------------------------------------------------------------------------| */
        printf("\n---> Boucle numéro %d :", i);
        printf("\nSous-ensemble actuel : ");
        affichage_liste(&ss_ens_S);
/*  |----------------------------------------------------------------------------------------------------------------|  
    |----------------------------------------------------------------------------------------------------------------| */

        /*  On extrait un élément non nul du noyau 
            (qui représentera les éléments du sous-ensemble à considérer dans les produits t et s) */
        if(noyau(&matrice, b, b+1, &ker)==0){
            printf("\nL'élément du noyau n'est pas bon : défaut fonction noyau()\n");
            continue;
        }
        else{
            printf("\nUn élément du noyau est trouvé avec succès: ");
        }

        /* Calcule du produit t (avec l'actuel sous-ensemble) */
        for(int i=0; i<b+1; i++){
            if(ker.element[i] == 1){
                mpz_mul(t, t, S.t_element[ss_ens_S.element[i]]);
            }
        }
        /* Calcul des exposants de s */
        valuations_s = calloc(b, sizeof(int));
        for(int j=0; j<b; j++){
            valuations_s[j] = 0;
            for(int i=0; i<b+1; i++){
                if(ker.element[i] == 1){
                    valuations_s[j] = valuations_s[j] + S.element_dec[ss_ens_S.element[i]].valuation[j];
                }
            }
            valuations_s[j] = valuations_s[j]/2;
        }
        /* Calcul du produit s */
        for(int j=0; j<b; j++){
            mpz_pow_ui(facteur_s, B.element[j], valuations_s[j]);
            mpz_mul(s, s, facteur_s);
        }

        // Les éléments de comparaison
        mpz_neg(_s_reduit, s);
        mpz_mod(t_reduit, t, *n);
        mpz_mod(s_reduit, s, *n);
        mpz_mod(_s_reduit, _s_reduit, *n);
        
/*  |----------------------------------------------------------------------------------------------------------------|
    |                                                  AFFICHAGES                                                    |   
    |----------------------------------------------------------------------------------------------------------------| */
        /* printf("\n---> Boucle numéro %d :", i);
        printf("\nSous-ensemble actuel : ");
        affichage_liste(&ss_ens_S);
        printf("\nMatrice pivoté : \n");
        affichage_matrice(B.taille, B.taille+1, &matrice); */
        affichage_liste(&ker);
        gmp_printf("\nOn a t = %Zd, s = %Zd", t, s);
        gmp_printf("\nOn a t mod n = %Zd, s mod n = %Zd, -s mod n = %Zd\n", t_reduit, s_reduit, _s_reduit);
        mpz_t diff, s_carre;
        mpz_inits(diff, s_carre);
        mpz_mul(diff, t, t);
        mpz_mul(s_carre, s, s); 
        mpz_sub(diff, diff, s_carre);
        mpz_mod(diff, diff, *n);
        gmp_printf("On teste la nullité : t^2-s^2 mod n = %Zd\n", diff);
        mpz_clears(diff, s_carre);
/*  |----------------------------------------------------------------------------------------------------------------|  
    |----------------------------------------------------------------------------------------------------------------| */

        free(ss_ens_S.element);
        free(matrice.element);
        free(ker.element);
        free(valuations_s);

        /*  Si t = s [n] ou t = -s [n], alors on n'obtient pas de facteur non trivial.
            On refait tout avec un nouvel sous-ensemble  */
        if( mpz_cmp(t_reduit, s_reduit) == 0 || mpz_cmp(t_reduit, _s_reduit) == 0 ){
            mpz_set_str(t, "1", 10);
            mpz_set_str(s, "1", 10);
        }
        /*  Si t != +-s [n], on a alors un facteur non trivial qui est pgcd(t+s, n)
            On le stocke dans 'resultat' et l'algorithme se termine  */
        else{
            mpz_t somme_t_s;
            mpz_init(somme_t_s);
            mpz_add (somme_t_s, t, s);
            mpz_gcd(*resultat, somme_t_s, *n);
            mpz_clear(somme_t_s);
            mpz_clear(t); mpz_clear(s); mpz_clear(t_reduit); mpz_clear(s_reduit); mpz_clear(_s_reduit); mpz_clear(facteur_s);
            vider_base_ensemble(&B, &S);
            return 0;
        }
    }
    mpz_clear(t); mpz_clear(s); mpz_clear(t_reduit); mpz_clear(s_reduit); mpz_clear(_s_reduit); mpz_clear(facteur_s);
    vider_base_ensemble(&B, &S);
    return 1;
}