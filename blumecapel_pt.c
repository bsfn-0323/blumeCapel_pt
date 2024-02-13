/*by Stefano Bae*/
//Simulazione Blume-Capel per reticolo 2D o ER Random Graph
//Parallel Tempering
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Modificare i commenti per selezionare il reticolo o il grafo random
#define LAT  // Lattice
#define THERM 100000
//#define ERRG //ER Random Graph
#define EXPORTCONFIGS
#define STDMCSTEP 1
#define EVERY 10

#define M 8
//Modificare NNRANGE per impostare il numero massimo di vertici
#define NNRANGE 10

/*----------Variabili Globali----------*/
double inv3 = 1./3.;
int N, L;           //Numero di spins, dimensione del reticolo
int intmat = 0;     
double T, beta;     //Temperatura e temperatura inversa
double Tmin, Tmax,Bmax,Bmin;
double mu;          //Potenziale chimico
double ****pacc;     //Look up table per le probabilita' di accettazione
double p;  
 
double dmin(double x, double y){
    int t= (int)x<y;
    return (double)(t*x + (1-t)*y);
}
void exportArrayToBinary2(const char *filename, double **array, int dim1, int dim2) {
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
    }

    // Write the array data to the file
    for (int i = 0; i < dim1; i++) {
            fwrite(array[i], sizeof(double), dim2, file);
    }

    //fwrite(array, sizeof(double), dim1*dim2*dim3, file);
    fclose(file);
}
void exportArrayToBinary3(const char *filename, double ***array, int dim1, int dim2, int dim3) {
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
    }

    // Write the array data to the file
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            fwrite(array[i][j], sizeof(double), dim3, file);
        }
    }

    //fwrite(array, sizeof(double), dim1*dim2*dim3, file);
    fclose(file);
}
void exportArrayIntToBinary3(const char *filename, int ***array, int dim1, int dim2, int dim3) {
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
    }

    // Write the array data to the file
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            fwrite(array[i][j], sizeof(int), dim3, file);
        }
    }

    //fwrite(array, sizeof(double), dim1*dim2*dim3, file);
    fclose(file);
}
void exportArrayIntToBinary2(const char *filename, int **array, int dim1, int dim2) {
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
    }

    // Write the array data to the file
    for (int i = 0; i < dim1; i++) {
        fwrite(array[i], sizeof(int), dim2, file);
    }

    fclose(file);
}

void exit_failure(char *s){
    fprintf(stderr,"%s",s);
    exit(EXIT_FAILURE);
}

double drand(){
    return (double) rand()/RAND_MAX;
}

void init_rand(){
     unsigned int myrand;
    FILE *devran = fopen("/dev/random","r");
    size_t k=fread(&myrand,4,1,devran);
    fclose(devran);
    srand(myrand);
}

int new_spin(int spin){
    return spin + 3*(spin == -2) - 3*(spin == 2);
}

void init_coupling(int N,int **map,int **Jmat){
    double r;
    #ifdef LAT 
    int n,m;
    for(int i = 0;i<N;i++){
        if(intmat==0){
            //Costruisco la matrice di interazione e creo una mappa dei primi vicini
            n = i%L;
            m = floor(i/L);

            map[i][0] = (N+i-1)%N;
            map[i][1] = (N+i+1)%N;
            map[i][2] = n+L*((L+m-1)%L);
            map[i][3] = n+L*((L+m+1)%L);
            map[i][4] = -1;
            Jmat[i][(N+i+L)%N] = 1;
            Jmat[i][(N+i-1)%N] = 1;
            Jmat[i][(N+i+1)%N] = 1;
            Jmat[i][(N+i-L)%N] = 1;
        }
    }
    #endif
    #ifdef ERRG
    int count;
    int countj[N];
    for(int i = 0;i<N;i++) countj[i]=0;
    
    for(int i = 0;i<N;i++){
        if(intmat == 0){
            count = 0;
            for(int j = i+1;j<N;j++){
                //Per ogni coppia ij di spin, se p>r aggiungi il collegamento
                r = drand();
                if(p>r){
                    map[i][countj[i]+count]=j;
                    map[j][countj[j]] = i;
                    Jmat[i][j] = 1;
                    Jmat[j][i] = 1;
                    
                    count++;
                    countj[j]++;
                }
            }
            map[i][countj[i]+count]=-1;
        }
    }
    #endif
}
void init_config(int fill, int *s){
    int i;
    double r;
    if(abs(fill) > 1){
        for(i = 0;i<N;i++){
            r = drand();
            if(r<inv3) s[i] =-1;
            if((r>=inv3)&&(r<2.*inv3)) s[i] = 0;
            if(r>= 2.*inv3) s[i] = 1;
        }
    }else{
        for(i = 0;i<N;i++) s[i] = fill;
    }    

}
//Imposta le probablita' di accettazione
void init_pacc(double *betas){
    int i,j,k,t;
    for(t = 0;t<M;t++){
        for(i = 0;i<3;i++){
            for(j = 0;j<5;j++){
                for(k = 0;k<(2*NNRANGE+1);k++) pacc[t][i][j][k] = exp(-(double)betas[t]*(-(j-2)*(k-NNRANGE) + mu*(i-1)));
            }
        }   
    }
}

//Uno step singolo con condizioni elicoidali
void one_sweep_heli(int N,int **s,int **map, double *dE, int *dM,int *dRho,int rep){
    int idx;
    int count;
    int nnsum;
    int ns,os;
    double deltaE;
    int deltaRho;
    double r;

    for(int i = 0;i<N;i++){
        idx = rand()%N;
        count = 0;
        nnsum = 0;
        while(map[idx][count]!=-1){
            nnsum += s[rep][map[idx][count]];
            count++;
        }
        r = drand();
        os = s[rep][idx];
        ns = new_spin(os + (r<0.5) - (r>=0.5));
        deltaRho = ns*ns-os*os;
        deltaE = -(double)(ns-os)*nnsum + mu*deltaRho;

        if(deltaE <= 0.){
            s[rep][idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }else if(pacc[rep][deltaRho+1][(ns-os)+2][nnsum+NNRANGE] > drand()){
            s[rep][idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }

    }
}

void init_obs(int N,double *energy, int *magn, int *rho, int **s, int **map){
    
    for(int rep = 0;rep<M;rep++){
    double tmp_e = 0.;
    int tmp_rho,tmp_m;
    int nnsum =0;
    tmp_rho = tmp_m = 0;
    int count=0;
        for(int i = 0;i<N;i++){
            count = 0;
            nnsum = 0;
            while(map[i][count]!= -1){
                nnsum += s[rep][map[i][count]];
                count++;
            }
            tmp_e += -0.5*(double)s[rep][i]*nnsum + mu*s[rep][i]*s[rep][i];
            tmp_m += s[rep][i];
            tmp_rho += s[rep][i]*s[rep][i];
        }
        energy[rep] = tmp_e;
        magn[rep] = tmp_m;
        rho[rep] = tmp_rho;   
    }
}

void init_temp(double * betas){
    double db = (double)(Bmax -Bmin)/M;
    for(int i = 0;i<M;i++)
        betas[i] = Bmax - (double)db*i;
}

int main(int argc, char **argv){
    int MCS,outSize;
    //FILE *fpConfig, *fpObs,*fpJ;
    char fname[100],dirname[50],cmd[150];

    int **map,**Jmat;
    double ***output;int ***configs;

    //Impostazione parametri di simulazione
    if(argc != 7)
        exit_failure("Inserire gli argomenti nel seguente ordine: L, mu, Tmin, Tmax, p, MCS\n");
    
    L = atoi(argv[1]);
    mu = atof(argv[2]);
    Tmin = atof(argv[3]);
    Tmax = atof(argv[4]);
    Bmin = 1./Tmax;
    Bmax = 1./Tmin;
    //p = atof(argv[5]);
    MCS = atoi(argv[6]);
    outSize = (int) MCS/EVERY;
    //Setup Lattice size
    N = L*L;
    p = (double)4./N;

    int** s;
    s = (int**)malloc(sizeof(double*)*M);
    for(int i = 0;i<M;i++)
        s[i] = (int*)malloc(sizeof(double)*N);
    double betas[M];
    
    //Allocazione memoria
    Jmat = (int**)malloc(sizeof(int*)*N);
    map = (int**)malloc(sizeof(int*)*N);
    output =(double***)malloc(sizeof(double**)*M);
    configs=(int***)malloc(sizeof(int**)*M);
    pacc = (double****)malloc(M*sizeof(double***));
    for(int t = 0;t<M;t++){
        pacc[t] = (double***)malloc(3*sizeof(double**));
        for(int i = 0; i<3;i++){
            pacc[t][i] = (double**)malloc(5*sizeof(double*));
            for(int j = 0;j<5;j++){
                pacc[t][i][j] = (double*)malloc((2*NNRANGE+1)*sizeof(double));
            }
        }
    }
    for(int i = 0;i<N;i++){
        Jmat[i] = (int*)malloc(sizeof(int)*N);
        map[i] = (int*)malloc(sizeof(int)*N);
        for(int j = 0;j<N;j++){
            Jmat[i][j] = 0;
            map[i][j] = 0;
        }
    }
    for(int i = 0;i<M;i++){
        output[i] = (double**)malloc(sizeof(double*)*3);
        for(int j = 0;j<3;j++){
            output[i][j] = (double*)malloc(sizeof(double)*outSize);
        }
    }
    for(int i = 0;i<M;i++){
        configs[i] = (int**)malloc(sizeof(int*)*outSize);
        for(int j = 0;j<outSize;j++){
            configs[i][j] = (int*)malloc(sizeof(int)*N);
        }
    }
    #ifdef LAT
    sprintf(dirname,"bc_pt_L%d_Tmin%.3f_Tmax%.3f_mu%.3f",L,Tmin,Tmax,mu);
    #endif
    #ifdef ERRG
    sprintf(dirname,"bc_pt_L%d_Tmin%.3f_Tmax%.3f_mu%.3f_rg",L,Tmin,Tmax,mu);
    #endif
    sprintf(cmd,"rm -rv %s",dirname);
    system(cmd);
    sprintf(cmd,"mkdir %s",dirname);
    system(cmd);
    //Inizializzo Interazioni, Configurazioni, Accettazioni
    init_rand();
    init_coupling(N,map,Jmat);
    for(int r = 0 ;r<M;r++)
        init_config(2,s[r]);
    init_temp(betas);
    //Inizializzo probabilita' di accettazione
    init_pacc(betas);
    
    double en[M],de;
    int magn[M], rho[M], dm, drho;
    init_obs(N,en,magn,rho,s,map);
    //double pacc_t;
    int tmp_s;
    double tmp_en,tmp_m,tmp_rho;

    clock_t begin = clock();
    for(int mcStep = 0;mcStep<MCS+THERM;mcStep++){
        //Uno Step Monte Carlo e' costituito da:
        //nMetropolis+betaSwap
        for(int stdmc = 0; stdmc<STDMCSTEP; stdmc++){
            for(int rep = 0;rep<M;rep++){
                de= 0.;
                drho = dm = 0;
                one_sweep_heli(N,s,map,&de,&dm,&drho,rep);
                magn[rep] += dm;
                en[rep] += de;
                rho[rep] += drho;
            }
        }
        //rep->rep+1
        for(int rep = 0;rep<(M-1);rep++){
            de = -(betas[rep]-betas[rep+1])*(en[rep]-en[rep+1]);
            if(de<=0){
                for(int k =0;k<N;k++){
                    tmp_s = s[rep+1][k];
                    s[rep+1][k] = s[rep][k];
                    s[rep][k] = tmp_s;
                }
                tmp_en = en[rep+1];
                en[rep+1] = en[rep];
                en[rep] = tmp_en;

                tmp_m = magn[rep+1];
                magn[rep+1] = magn[rep];
                magn[rep] = tmp_m;

                tmp_rho = rho[rep+1];
                rho[rep+1] = rho[rep];
                rho[rep] = tmp_rho;
                }else if(dmin(1.,exp(-de))>drand()){
                for(int k =0;k<N;k++){
                    tmp_s = s[rep+1][k];
                    s[rep+1][k] = s[rep][k];
                    s[rep][k] = tmp_s;
                }
                tmp_en = en[rep+1];
                en[rep+1] = en[rep];
                en[rep] = tmp_en;

                tmp_m = magn[rep+1];
                magn[rep+1] = magn[rep];
                magn[rep] = tmp_m;

                tmp_rho = rho[rep+1];
                rho[rep+1] = rho[rep];
                rho[rep] = tmp_rho;
            }
        }
        if((mcStep%EVERY== 0) & (mcStep>=THERM)){
            int idx = (int)(mcStep-THERM)/EVERY;
            for(int rep = 0;rep<M;rep++){
                output[rep][0][idx] = (double)en[rep]/N;
                output[rep][1][idx] = (double)magn[rep]/N;
                output[rep][2][idx] = (double)rho[rep]/N;
                #ifdef EXPORTCONFIGS
                for(int i = 0;i<N;i++){
                    configs[rep][idx][i] = s[rep][i];
                }
                #endif
            }
        }
        
    }
    clock_t end =clock();
    double runtime = (double)(end-begin)/CLOCKS_PER_SEC;
    //printf("<e> = %f\t<m> = %f\t<rho> = %f\n",avg_en/(N*MCS),(double)avg_magn/(N*MCS),(double)avg_rho/(N*MCS));
    printf("\nRuntime = %.4f [s]\n",runtime);
    //Export Arrays
    sprintf(fname,"%s/J0.bin",dirname);
    exportArrayIntToBinary2(fname,Jmat,N,N);
    for(int i = 0;i<M;i++){
    sprintf(fname,"%s/output_B%.3f.bin",dirname,betas[i]);
    exportArrayToBinary2(fname,output[i],3,outSize);
    }
    #ifdef EXPORTCONFIGS
    for(int i = 0;i<M;i++){
        sprintf(fname,"%s/configs_B%.3f.bin",dirname,betas[i]);
        exportArrayIntToBinary2(fname,configs[i],outSize,N);
    }
    #endif
    //Free memories
    for(int t = 0;t<M;t++){
        for(int i = 0; i<3;i++){
            for(int j = 0;j<5;j++){
                free(pacc[t][i][j]); 
            }
            free(pacc[t][i]);
        }
        free(pacc[t]);
    }
    free(pacc);
    for(int i = 0;i<N;i++){
        free(Jmat[i]);
        free(map[i]);
    }
    free(Jmat);
    free(map);
    for(int i = 0;i<M;i++){
        for(int j = 0;j<3;j++){
            free(output[i][j]);
        }
        free(output[i]);
    }
    free(output);
    for(int i = 0;i<M;i++){
        free(s[i]);
    }
    free(s);
    for(int i = 0;i<M;i++){
        for(int j = 0;j<outSize;j++){
            free(configs[i][j]);
        }
        free(configs[i]);
    }
    free(configs);
    return 0;
}