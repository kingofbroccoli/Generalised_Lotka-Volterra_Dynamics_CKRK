/* PREPROCESSOR INCLUSION */
//#include <cmath>  // Sono già contenuti in altri include ma per fortuna sono anche dotati di include guard
//#include <ctime>  // Forse si possono cancellare
//#include <cstring>
//#include <cstdio> 
//#include <cstdlib>

#include "pietro_toolbox.h"
#include "graph_lv_toolbox.h"
#include "runge-kutta_toolbox.h"

/* PREPROCESSOR DEFINITION */

#define MAIL_NOTIFICATION
#undef MAIL_NOTIFICATION

#define MY_GNUPLOT
#undef MY_GNUPLOT

/* STRUCTURE */

/* TYPE DEFINITION */

/* PROTOTYPES */

/* MAIN BEGINNING */

int main(int argc, char **argv){
    clock_t begin = clock(); //Clock iniziale
    float time_spent;
    char *S_label = argv[1]; // Ricorda argv[0] è il nome dell'eseguibile!
    int S = atoi(S_label);
    char *sigma_label = argv[2];
    double sigma = atof(sigma_label);
    char *mu_label = argv[3];
    double mu = atof(mu_label);
    int N_estrazioni = atoi(argv[4]); 
    int N_previous_ext = atoi(argv[5]);
    int N_misure = atoi(argv[6]);
    int N_perm = 0; 

    double p_A = 1.0, p_M = 0.50;
    char interaction_kind[] = "A";
    char ia_label[] = "Antagonistic Gaussian";
    void (*ia_generator)(double*, double*, double, double, double, double, double, double) = &mixture_gaussian_extraction;
    char gr_label[] = "Trivial";
    double (*gr_generator)(double, double) = &trivial_generator; // Pointer to a function (weird syntax but somehow it has a reason)
    double gr_mu = 1.0, gr_sigma = 0.0;
    //double gr_value = 1.0; 
    //double growth_rate_interval[2] = {1.0, 1.0};
    char cc_label[] = "Trivial";
    double (*cc_generator)(double, double) = &trivial_generator; // Pointer to a function (weird syntax but somehow it has a reason)
    double cc_mu = 1.0, cc_sigma = 0.0;
    //double carrying_cap_interval[2] = {1.0, 1.0};
    //double cc_value = 1.0;
    void (*save_lv_function)(graph*, double**, double, double*, FILE*, FILE*) = &save_LotkaVolterra_nothing;

    double c = 2.0;
    bool sparsity_flag = MY_TRUE;
    double population_factor = 1.2; // Fattore estrazione popolazione iniziale
    double lambda = pow(10, -10); // Tasso di Migrazione
    double h_first = pow(10, -2);
    double h_min = pow(10, -6);
    double t_max = 1000;
    double deltat_save = 1.0;
    double eps = pow(10, -6);
    double stationary_threshold = pow(10, -4); // Soglia di stazionarietà --- //### DEVO VALUTARE MEGLIO COME CONSIDERARE L EQUILIBRIO
    double extinction_log_threshold = log(pow(10, -6)); // Soglia di estinzione
    double delta = 0.01; // Soglia per considerare due punti di equilibrio distinti

    double p_legame;
    double divergence_probability, multi_eq_probability;
    bool divergence, stepsize_trouble;
    int j, n;  
    int divergence_counter, multi_eq_counter;
    int eq_points_counter;
    int *eq_points_statistics;
    int good=0, bad=0;
    int char_lenght = 200; 
    char *name_buffer, *dir_name;
    matrix *equilibrium_points;
    graph* ecosystem;
    FILE *fp, *fp_RK, *fp_speed, *fp_eq, *fp_r, *fp_K; // Pointer to File

    printf("Stiamo eseguendo: %s\n", argv[0]);
    printf("Orario di inizio: \n");
    system("date");
    name_buffer = my_char_malloc(char_lenght); // Memoria liberata esplicitamente nel main
    dir_name = my_char_malloc(char_lenght); // Memoria liberata esplicitamente nel main
    srand48(time(0)); // Inizializzazione generatore numeri casuali, va fatto una volta sola nel main IMPORTANTE.

    // Creazione directory mu e sigma
    snprintf(name_buffer, sizeof(char) * char_lenght, "mkdir mu_%s_sigma_%s", mu_label, sigma_label); // sizeof(char)=1 so it can be omitted
    system(name_buffer); 
    // Inizializziamo la struttura generale dell'ecosistema
    ecosystem = ecosystem_initialization(S); // Memoria liberata in erase_ecosystem
    // Calcoliamo la probabilità di legame
    p_legame = c / (ecosystem->size - 1.);
    // Creazione directory
    snprintf(dir_name, char_lenght, "mu_%s_sigma_%s/S_%d_c_%.2lf", mu_label, sigma_label, ecosystem->size, c); // sizeof(char)=1 so it can be omitted
    create_bunch_of_directories(dir_name);    
    // Scriviamo su file i parametri del codice
    snprintf(name_buffer, CHAR_LENGHT, "%s/Parameters.txt", dir_name); // sizeof(char)=1 so it can be omitted
    fp = my_open_writing_file(name_buffer);
    fprintf(fp, "# PARAMETERS: \n# N_ext = %d \t N_previous_ext = %d \t N_measure = %d \t N_perm = %d\n# Graph Structure: S = %d \t c = %lf\n", N_estrazioni, N_previous_ext, N_misure, N_perm, S, c);
    fprintf(fp, "# Interactions: %s mu = %lf \t sigma = %lf\n# Growth Rates: %s mu = %lf \t sigma = %lf# Carrying Capacities: %s mu = %lf \t sigma = %lf\n", ia_label, mu, sigma, gr_label, gr_mu, gr_sigma, cc_label, cc_mu, cc_sigma);
    fprintf(fp, "# Solving log_LV with Adaptive CKRK: accuracy eps = %E \t h_first = %E \t h_min = %E \t deltat_save = %lf\n", eps, h_first, h_min, deltat_save);
    fprintf(fp, "# Stability threshold = %E or maximum time t_max = %lf\n# Immigration rate: %E \t Extinction_LogThreshold: %lf \t Popolation factor: %lf\n# Two equilibrium species are considered the same if closer than delta = %lf\n", stationary_threshold, t_max, lambda, extinction_log_threshold, population_factor, delta);
    fclose(fp);
    // I contatori di divergenze e multi_eq vanno resettati ad ogni taglia S dell'ecosistema
    divergence_counter = 0;
    multi_eq_counter = 0;
    // Estraiamo N_estrazioni volte il grafo e per ognuna di esse evolviamo N_misure volte diverse il sistema modificando ogni volta le condizioni iniziali
    for(n=(N_previous_ext+1); n<=(N_previous_ext+N_estrazioni); n++){
        // Estrazione Grafo
        graph_extraction(ecosystem, p_legame, p_A, p_M, ia_generator, mu, sigma, mu, sigma, gr_generator, gr_mu, gr_sigma, cc_generator, cc_mu, cc_sigma); // Memoria liberata in clean_up_graph
        // Salvataggio growth rate, carrying capacity e matrice di interazione nel formato desiderato
        snprintf(name_buffer, char_lenght, "%s/Extractions/Extraction_%d_Growth_Rates.txt", dir_name, n);
        fp_r = my_open_writing_file(name_buffer);
        snprintf(name_buffer, char_lenght, "%s/Extractions/Extraction_%d_Carrying_Capacities.txt", dir_name, n);
        fp_K = my_open_writing_file(name_buffer);
        snprintf(name_buffer, char_lenght, "%s/Extractions/Extraction_%d_Interaction_Matrix_Dense.txt", dir_name, n);
        fp = my_open_writing_file(name_buffer);
        save_ecosystem_to_file(ecosystem, sparsity_flag, fp_r, fp_K, fp);
        fclose(fp);
        fclose(fp_r);
        fclose(fp_K);
        // Preparativi in vista delle misure
        eq_points_counter = 0;
        equilibrium_points = my_double_matrix_malloc(1, ecosystem->size); // Memoria liberata in clean_equilibrium_points
        eq_points_statistics = my_int_calloc(1); // Memoria liberata in clean_equilibrium_points
        for(j=1; j<=N_misure; j++){
            // Resettiamo contatori e status legati alla dinamica ed estraiamo condizioni iniziali
            snprintf(name_buffer, char_lenght, "%s/Initial_Conditions/Log_Initial_Condition_Extraction_%d_Measure_%d.txt", dir_name, n, j);
            fp = my_open_writing_file(name_buffer);
            extract_and_save_random_initial_conditions(ecosystem, population_factor, fp);
            fclose(fp);
            // Integrazione numerica con RKCK_Adattivo
            snprintf(name_buffer, char_lenght, "%s/Evolutions/Lotka-Volterra_Extraction_%d_Measure_%d.txt", dir_name, n, j);
            fp_RK = my_open_writing_file(name_buffer); // Apertura file scrittura
            snprintf(name_buffer, char_lenght, "%s/Evolutions/Lotka-Volterra_Speed_Extraction_%d_Measure_%d.txt", dir_name, n, j);
            fp_speed = my_open_writing_file(name_buffer); // Apertura file scrittura
            log_RungeKutta_adaptive_driver(ecosystem, t_max, h_first, h_min, eps, stationary_threshold, j, &divergence, &stepsize_trouble, &good, &bad, lambda, deltat_save, save_lv_function, fp_RK, fp_speed);
            fflush(fp_RK);
            fclose(fp_RK);
            fclose(fp_speed);
            // Qualora sia emerso un stepsize_trouble controlliamo che sia dovuto ad una divergenza delle abundances, altrimenti solleviamo un errore.
            if(stepsize_trouble == MY_TRUE){
                // Ripartiamo dalle condizioni iniziali estratte precedentemente
                snprintf(name_buffer, CHAR_LENGHT, "%s/Initial_Conditions/Log_Initial_Condition_Extraction_%d_Measure_%d.txt", dir_name, n, j);
                fp = my_open_reading_file(name_buffer);
                load_initial_conditions(ecosystem, fp);
                fclose(fp);
                // Tentiamo integrazione numerica con RK4 con passo integrazione fisso per verificare che il stepsize_trouble sia dovuto ad una divergenza
                snprintf(name_buffer, CHAR_LENGHT, "%s/Evolutions/Lotka-Volterra_Extraction_%d_Measure_%d.txt", dir_name, n, j);
                fp_RK = my_open_writing_file(name_buffer); // Apertura file scrittura
                snprintf(name_buffer, CHAR_LENGHT, "%s/Evolutions/Lotka-Volterra_Speed_Extraction_%d_Measure_%d.txt", dir_name, n, j);
                fp_speed = my_open_writing_file(name_buffer); // Apertura file scrittura
                log_RungeKutta4_simple_driver(ecosystem, stationary_threshold, ((int)(t_max/h_first)), j, h_first, lambda,  &(divergence), save_lv_function, fp_RK, fp_speed);
                fflush(fp_RK);
                fclose(fp_RK);
                fclose(fp_speed);
                if(divergence == MY_FALSE){
                    // If we don't find any divergence we raise an ERROR.
                    printf("\n STEPSIZE ERROR: Step size too small in log_RungeKutta_adaptive_driver and no divergence identified in log_RungeKutta4_simple_driver_with_speed!\n");
                    exit(MY_TROUBLE);
                }
            }
            printf("Extraction %d Measure %d \n", n, j);
            if(divergence == MY_FALSE){ // Svolgiamo queste operazioni (delicate ma importanti) solo in assenza di divergenza
                // Classificazione dei vari status delle specie di equilibrio dopo l'evoluzione dell'ecosistema
                classify_equilibrium_populations_status(ecosystem, extinction_log_threshold);
                // Salviamo il punto di equilibrio raggiunto
                snprintf(name_buffer, char_lenght, "%s/Extractions/Extraction_%d", dir_name, n);
                add_equilibrium_point(ecosystem, equilibrium_points, &(eq_points_counter), eq_points_statistics, delta, N_misure, N_perm, population_factor, sparsity_flag, name_buffer);                            
            }
            else{
                divergence_counter++;
            }
        } // Fine loop sulle misure
        // Salvataggio statistica popolazioni di equilibrio raggiunte e pulizia memoria 
        snprintf(name_buffer, char_lenght, "%s/Equilibrium_Point_Statistics.txt", dir_name);
        fp_eq = my_open_appending_file(name_buffer);
        clean_equilibrium_points(equilibrium_points, eq_points_statistics, eq_points_counter, &(multi_eq_counter), N_misure, n, fp_eq); // Qua dentro stampiamo e puliamo
        fclose(fp_eq);
        // Pulizia interazioni grafo
        clean_up_graph(ecosystem);
    } // Fine del ciclo sulle estrazioni
    // Calcolo e scrittura probabilità divergenza
    divergence_probability = (double) divergence_counter / (N_estrazioni*N_misure);
    snprintf(name_buffer, char_lenght, "Statistics_Divergence.txt"); 
    fp = my_open_appending_file(name_buffer);
    fprintf(fp, "%lf\t%d\t%lf\t%s\t%.2f\t%d\n", sigma, ecosystem->size, divergence_probability, interaction_kind, c, N_estrazioni*N_misure);
    fclose(fp);
    // Calcolo e scrittura probabilità di avere attrattori multipli
    multi_eq_probability = (double) multi_eq_counter / N_estrazioni;
    snprintf(name_buffer, char_lenght, "Statistics_Multiple_Equilibrium.txt"); 
    fp = my_open_appending_file(name_buffer);
    fprintf(fp, "%lf\t%d\t%lf\t%s\t%.2f\t%d\n", sigma, ecosystem->size, multi_eq_probability, interaction_kind, c, N_estrazioni);
    fclose(fp);
    // Pulizia memoria
    erase_ecosystem(ecosystem); // Pulizia memoria grafo
    free(name_buffer); // Pulizia memoria stringa
    free(dir_name); // Pulizia memoria stringa
    // Calcolo tempo impiegato
    clock_t end = clock(); //Clock conclusivo
    time_spent = ((float)(end - begin)) / CLOCKS_PER_SEC; //Calcolo tempo di esecuzione
    printf("\n # Il programma ha impiegato %f secondi.\n # A breve dovrebbe arrivare la mail (se opportunamente richiesto presso i nostri uffici).\n \n", time_spent);
    printf("Orario di conclusione: \n");
    system("date");

    return MY_SUCCESS;
}

/* MAIN - THE END */

/* FUNCTIONS */