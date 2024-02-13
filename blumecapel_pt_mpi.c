/*by Stefano Bae*/
//Simulazione Blume-Capel per reticolo 2D o ER Random Graph
//Parallel Tempering
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define inv3 1./3.
//Modificare i commenti per selezionare il reticolo o il grafo random
//#define LAT  // Lattice
//#define ERRG_0 //ER Random Graph
#define ERRG_1 //ER Random Graph
#define ERRG
#define THERM 100000
#define METROPOLIS_STEPS 25
#define EVERY 10

#define NNRANGE 20
//#define EXPORTCONFIGS
/*************************/
void importBinaryIntToArray2(const char *filename, int **array, int dim1, int dim2){
    FILE *file = fopen(filename, "rb");
        
    if (file == NULL) {
        fprintf(stderr, "Error opening file for reading.\n");
    }
    for (int i = 0; i < dim1; i++) {
        size_t read = fread(array[i], sizeof(int), dim2, file);
    }
    fclose(file);
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

/*************************/
void exit_failure(char *s){
    fprintf(stderr,"%s",s);
    exit(EXIT_FAILURE);
}
double dmin(double x, double y){
    int t= (int)x<y;
    return (double)(t*x + (1-t)*y);
}
int new_spin(int spin){
    return spin + 3*(spin == -2) - 3*(spin == 2);
}
/*******Random Variables********/
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
/******************************/
void init_config(int N, int *s){
    int i;
    double r;
    for(i = 0;i<N;i++){
        r = drand();
        if(r<inv3) s[i] =-1;
        if((r>=inv3)&&(r<2.*inv3)) s[i] = 0;
        if(r>= 2.*inv3) s[i] = 1;
    }

}
void init_coupling(int L,int **map,int **Jmat,double p){
    int N = L*L;
    #ifdef LAT 
    int n,m;
    for(int i = 0;i<N;i++){
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
    #endif
    #ifdef ERRG_0
    double r;
    int count;
    int countj[N];
    for(int i = 0;i<N;i++) countj[i]=0;
    
    for(int i = 0;i<N;i++){
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
    #endif
    /*A different way to generate graphs is G(N,M).
      I have N nodes and M edges. If I want a fixed connectivity 4,
      i need M = 2*N.
      Taken all the (i,j) couples, choose M of them.
    */
    #ifdef ERRG_1
    double r;
    int M = 2*N, countM = 0;
    int count[N];
    for(int k = 0;k<N;k++){
        count[k]=0;
    }
    int i,j;
    while(countM<M){
        //Extract 2 indices
        i = (int)rand()%N;
        j = (int)rand()%N;
        if((i!=j)&(Jmat[i][j]==0)){
            map[i][count[i]]=j;
            map[j][count[j]] = i;
            Jmat[i][j]=1;
            Jmat[j][i]=1;
            count[i]++;
            count[j]++;
            countM++;
        }
    }
    #endif
}
void init_pacc(double ***pacc,double beta,double mu){
    int i,j,k;
    for(i = 0;i<3;i++){
        for(j = 0;j<5;j++){
            for(k = 0;k<(2*NNRANGE+1);k++) pacc[i][j][k] = exp(-(double)beta*(-(j-2)*(k-NNRANGE) + mu*(i-1)));
        }
    }   
}
void init_obs(int N,double *energy, int *magn, int *rho, int *s, int **map,double mu){
    int nnsum =0;
    int count=0;
    (*energy) = 0.;
    (*magn) = 0;
    (*rho) = 0;
    for(int i = 0;i<N;i++){
        count = 0;
        nnsum = 0;
        while(map[i][count]!= -1){
            nnsum += s[map[i][count]];
            count++;
        }
        (*energy) += -0.5*(double)s[i]*nnsum + mu*s[i]*s[i];
        (*magn) += s[i];
        (*rho) += s[i]*s[i];
    }
    
}
/*******Metropolis*******/
void one_sweep_metropolis(int N,int *s,double ***pacc,int **map,double mu, double *dE, int *dM,int *dRho,double beta){
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
            nnsum += s[map[idx][count]];
            count++;
        }
        r = drand();
        os = s[idx];
        ns = new_spin(os + (r<0.5) - (r>=0.5));
        deltaRho = ns*ns-os*os;
        deltaE = -(double)(ns-os)*nnsum + mu*deltaRho;

        if(deltaE <= 0.){
            s[idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }
        else if(pacc[deltaRho+1][(ns-os)+2][nnsum+NNRANGE] > drand()){
            //printf("nnSum = %d, pacc = %.3f, pacc* = %.3f\n", nnsum,pacc[deltaRho+1][(ns-os)+2][nnsum+NNRANGE], exp(-beta*deltaE));
            s[idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }

    }
}

int main(int argc,char**argv){
    int rank,size;
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    //Variables
    int MCS, outSize,rep;
    int L,N;
    int magn,dm,rho,drho,magn_recv_l,rho_recv_l,magn_recv_r,rho_recv_r;
    int *s,**configs,**Jmat,**map,*s_recv_r,*s_recv_l;
    double en,den,en_recv_l,en_recv_r;
    double denl,denr;
    double ***pacc;
    double **outputs,beta,beta_recv_l,beta_recv_r,Tmin, Tmax,T,dB,mu;
    double p;
    char dirname[50],cmd[75],fname[100];
    //Input
    if(argc !=7)
        exit_failure("Inserire gli argomenti nel seguente ordine: L, mu, Tmin, Tmax, p, MCS\n");
    L = atoi(argv[1]);
    mu = atof(argv[2]);
    Tmin = atof(argv[3]);
    Tmax = atof(argv[4]);
    MCS = atoi(argv[5]);
    rep = atoi(argv[6]);
    N = L*L;
    p = (double) 4./N;
    outSize = (int)MCS/EVERY;

    //Memory allocation
    s = (int*)malloc(sizeof(int)*N);
    s_recv_l = (int*)malloc(sizeof(int)*N);
    s_recv_r = (int*)malloc(sizeof(int)*N);
    Jmat = (int**)malloc(sizeof(int*)*N);
    map = (int**)malloc(sizeof(int*)*N);
    for(int i = 0;i<N;i++){
        Jmat[i] = (int*)malloc(sizeof(int)*N);
        map[i] = (int*)malloc(sizeof(int)*N);
        for(int j = 0;j<N;j++) {
            map[i][j] = -1;
            Jmat[i][j] = 0;}
    }
    #ifdef EXPORTCONFIGS
    configs = (int**)malloc(sizeof(int*)*outSize);
    for(int i = 0;i<outSize;i++){
        configs[i] = (int*)malloc(sizeof(int)*N);
    }
    #endif
    outputs = (double**)malloc(sizeof(double*)*outSize);
    for(int i = 0;i<outSize;i++){
        outputs[i] = (double*)malloc(sizeof(double)*3);
    }
    pacc = (double***)malloc(3*sizeof(double**));
    for(int i = 0; i<3;i++){
        pacc[i] = (double**)malloc(5*sizeof(double*));
        for(int j = 0;j<5;j++){
            pacc[i][j] = (double*)malloc((2*NNRANGE+1)*sizeof(double));
        }
    }
    //Define output directory
    #ifdef LAT
    sprintf(dirname,"bc_pt_L%d_Tmin%.3f_Tmax%.3f_mu%.3f",L,Tmin,Tmax,mu);
    #endif
    #ifdef ERRG_0
    sprintf(dirname,"bc_pt_L%d_Tmin%.3f_Tmax%.3f_mu%.3f_rg_%d",L,Tmin,Tmax,mu,rep);
    #endif
    #ifdef ERRG_1
    sprintf(dirname,"bc_pt_L%d_Tmin%.3f_Tmax%.3f_mu%.3f_rg_%d",L,Tmin,Tmax,mu,rep);
    #endif

    //Initialize
    init_rand();
    //The first one generates the couplings and the mapping
    if(rank==0){
        sprintf(cmd,"rm -r %s",dirname);
        system(cmd);
        sprintf(cmd,"mkdir %s",dirname);
        system(cmd);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
        init_coupling(L,map,Jmat,p);
        #ifdef ERRG
        sprintf(fname,"%s/J0_n%d.bin",dirname,rep);
        #endif 
        #ifdef LAT
        sprintf(fname,"%s/J0.bin",dirname);
        #endif
        exportArrayIntToBinary2(fname,Jmat,N,N);
        #ifdef ERRG
        sprintf(fname,"%s/map_n%d.bin",dirname,rep);
        #endif 
        #ifdef LAT
        sprintf(fname,"%s/map.bin",dirname);
        #endif
        exportArrayIntToBinary2(fname,map,N,N);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank>0){
        //Read Jmat and map
        #ifdef ERRG
        sprintf(fname,"%s/J0_n%d.bin",dirname,rep);
        #endif 
        #ifdef LAT
        sprintf(fname,"%s/J0.bin",dirname);
        #endif
        importBinaryIntToArray2(fname,Jmat,N,N);
        #ifdef ERRG
        sprintf(fname,"%s/map_n%d.bin",dirname,rep);
        #endif 
        #ifdef LAT
        sprintf(fname,"%s/map.bin",dirname);
        #endif
        importBinaryIntToArray2(fname,map,N,N);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    init_config(N,s);
    dB = (double) (1./Tmax - 1./Tmin)/size;
    beta =(double) 1./Tmin + dB*rank;
    T = (double)1./beta;
    init_pacc(pacc,beta,mu);
    init_obs(N,&en,&magn,&rho,s,map,mu);
    double rnum, rnum_recv;
    int check_barrier = 0;
    for(int mcStep = 0;mcStep<MCS+THERM;mcStep++){
        //Metropolis
        for(int metroStep = 0;metroStep<METROPOLIS_STEPS;metroStep++){
            den = 0.;dm=0;drho=0;
            one_sweep_metropolis(N, s,pacc,map,mu, &den,&dm,&drho,beta);
            en+=den;
            magn+=dm;
            rho+=drho;
        }
        //Parallel Tempering Move
        //Receive from left and right neighbourghs
        rnum = drand();
        for(int k = 0;k<size-1;k++){
            if(rank == k){
                MPI_Send(s,N,MPI_INT,rank+1,1,MPI_COMM_WORLD);
                MPI_Send(&en,1,MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD);
                MPI_Send(&magn,1,MPI_INT,rank+1,3,MPI_COMM_WORLD);

                MPI_Send(&rho,1,MPI_INT,rank+1,4,MPI_COMM_WORLD);
                MPI_Send(&rnum,1,MPI_DOUBLE,rank+1,5,MPI_COMM_WORLD);

                MPI_Recv(s_recv_r,N,MPI_INT,rank+1,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&en_recv_r,1,MPI_DOUBLE,rank+1,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&magn_recv_r,1,MPI_INT,rank+1,33,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                MPI_Recv(&rho_recv_r,1,MPI_INT,rank+1,44,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                double paccr = dmin(1.,exp(-dB*(en - en_recv_r)));
                if(paccr > rnum){
                    for(int i = 0;i<N;i++)
                        s[i] = s_recv_r[i];
                    en = en_recv_r;
                    magn = magn_recv_r;
                    rho = rho_recv_r;
                    //if(mcStep == MCS+THERM-1)
                    //    printf("k = %d [Rank %d] en = %.3f, magn = %.3f, rho = %.3f\n",k,rank,(double)en/N, (double)magn/N,(double)rho/N);
                }
                //if(mcStep == MCS+THERM-1)
                    //printf("k = %d [Rank %d] rnum = %.3f, rnum_recv = %.3f, pacc = %.3f, cr = %d\n",k,rank,rnum, rnum_recv,paccr,paccr > rnum);
            }
            if(rank == k+1){
                MPI_Recv(s_recv_l,N,MPI_INT,rank-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&en_recv_l,1,MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&magn_recv_l,1,MPI_INT,rank-1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                MPI_Recv(&rho_recv_l,1,MPI_INT,rank-1,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                
                MPI_Recv(&rnum_recv,1,MPI_DOUBLE,rank-1,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                MPI_Send(s,N,MPI_INT,rank-1,11,MPI_COMM_WORLD);
                MPI_Send(&en,1,MPI_DOUBLE,rank-1,22,MPI_COMM_WORLD);
                MPI_Send(&magn,1,MPI_INT,rank-1,33,MPI_COMM_WORLD);

                MPI_Send(&rho,1,MPI_INT,rank-1,44,MPI_COMM_WORLD);
                double paccl = dmin(1.,exp(-dB*(en_recv_l - en)));
                if(paccl > rnum_recv){
                    for(int i = 0;i<N;i++)
                        s[i] = s_recv_l[i];
                    en = en_recv_l;
                    magn = magn_recv_l;
                    rho = rho_recv_l;
                    //if(mcStep == MCS+THERM-1)
                    //    printf("k = %d [Rank %d] en = %.3f, magn = %.3f, rho = %.3f\n",k,rank,(double)en/N, (double)magn/N,(double)rho/N);
                }
                //if(mcStep == MCS+THERM-1)
                    //printf("k = %d [Rank %d] rnum = %.3f, rnum_recv = %.3f, pacc = %.3f,cl = %d\n",k,rank,rnum, rnum_recv,paccl,paccl > rnum_recv);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        //Print
        if((mcStep%EVERY==0) & (mcStep>=THERM)){
            int idx = (int)(mcStep-THERM)/EVERY;
            outputs[idx][0] = (double)en/N;
            outputs[idx][1] = (double)magn/N;
            outputs[idx][2] = (double)rho/N;
            #ifdef EXPORTCONFIGS
            for(int i = 0;i<N;i++)
                configs[idx][i]=s[i];
            #endif
        }
        /*if((mcStep == MCS+THERM-1)){
            printf("[RANK %d] <en> = %.3f, <m> = %.3f, <r> = %.3f\n",rank,(double)en/N,(double)magn/N,(double)rho/N);
        }*/
    }
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(fname,"%s/output_B%.3f.bin",dirname,beta);
    exportArrayToBinary2(fname,outputs,outSize,3);
    #ifdef EXPORTCONFIGS
    sprintf(fname,"%s/configs_B%.3f.bin",dirname,beta);
    exportArrayIntToBinary2(fname,configs,outSize,N);
    #endif
    //Free Mem
    for(int i = 0; i<3;i++){
        for(int j = 0;j<5;j++){
            free(pacc[i][j]); 
        }
        free(pacc[i]);
    }
    free(pacc);
    for(int i = 0;i<N;i++){
        free(Jmat[i]);
        free(map[i]);
    }
    free(Jmat);
    free(map);
    for(int j = 0;j<outSize;j++){
        free(outputs[j]);
    }
    free(outputs);
    free(s);
    free(s_recv_l);
    free(s_recv_r);
    #ifdef EXPORTCONFIGS
    for(int j = 0;j<outSize;j++){
        free(configs[j]);
    }
    free(configs);
    #endif
    MPI_Finalize();
    return 0;
}