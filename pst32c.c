/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define _USE_MATH_DEFINES
#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;

params myParams;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);

//int s=stru->N;
 
 
//riguardare
//dimensione

VECTOR backbone();


VECTOR backbone(params* input,VECTOR phi,VECTOR psi){

	int n = input->N;
	int i,j=0;
	
	//mi prendo la sequenza(la dimensione) cioè N
	
	
	//distanza n-calfa,calfa-c,c-N
	
	type dist_ca_n=1.46;
	
	type dist_ca_c=1.52;
	
	type dist_c_n=1.33;
	
	//angoli std del backbone
	
	type theta1=2.028;//theta-ca-cn
	
	type theta2=2.124;//theta-cn-ca
	
	type theta3=1.940;//theta-n-ca-c
	
	//ora dovrei allocare la matrice
	
	MATRIX coords=alloc_matrix(n*3,3);//n*3 ->righe;3 colonne
	
	//setto le tre coord
	
	coords[0]=0;
	coords[1]=0;
	coords[2]=0;

	coords[3]=dist_ca_n;
	coords[4]=0;
	coords[5]=0;
	//coords[2]=0;
	
	//coord vett 1
	
	
	type v1,v2,v3;
	
	MATRIX rot;
	
	VECTOR newv;
	
	//costruire la catena
	
	for(i=0;i<n;i++){
		
		int idx=i*3;//indice base ammin.
		
		if(i>0){
			
			type v1=coords[idx-1]-coords[idx-2];//vettore calfa meno c
			
			v1=v1/abs(v1);
			
			//calcolo v1/(|v1|)
			
			rot=rotation(&v1,theta2);
			
			newv[0]=rot[3]*dist_c_n;
			newv[1]=rot[4]*dist_c_n;
			newv[2]=rot[5]*dist_c_n;

			//newv={0,dist_c_n,0}*rot;//prod matr penso non corretto,va acceduto diversamente in memoria
			
			coords[idx]=coords[idx-1]+newv[0]+newv[1]+newv[2];//effettuare questa moltiplicazione
			
			//posiziono c alfa usando phi(vettore angolo phi)
			
			v2=coords[idx]-coords[idx-1];
			
			v2=v2/abs(v2);
			
			rot=rotation(&v2,phi[i]);

			newv[0]=rot[3]*dist_ca_n;
			newv[1]=rot[4]*dist_ca_n;
			newv[2]=rot[5]*dist_ca_n;
			
			//newv={0,dist_ca_n,0}*rot;//prodotto matr.//rivedere penso non corretto va acceduto diversamente in mem
			
			coords[idx+1]=coords[idx]+newv[0]+newv[1]+newv[2];
		
		
		
		}//if
			
			
			//Posiziono c usando il vettore degli angoli psi
			
	 	
	 	coords[idx+1]-coords[idx];
			
	 	v3=v3/abs(v3);
			
	 	rot=rotation(&v3,psi[i]);

		newv[0]=rot[3]*dist_ca_c;
		newv[1]=rot[4]*dist_ca_c;
		newv[2]=rot[5]*dist_ca_c;

		//newv={0,dist_ca_c,0}*rot;//rivedere 
			
	 	coords[idx+2]=coords[idx+1]+newv[0]+newv[1]+newv[2];
	
	
	
	}//for
	
	//alla fine del for li ritorno
	
    return coords;
			
	//prima di scrivere la funzione rotation, devo allocarmi in memoria un vettore axis di 3 dimensioni
	//uso il metodo alloc matrix con una colonna
	
	VECTOR axis=alloc_matrix(3,1);//tre dimensioni: x,y,z quindi 3 righe una colonna per come è strutturato il metodo e 
	//le allocazioni in memoria
	//todo-->metterlo sopra all'inizio del metodo 



}//backbone



//funz rotation ok!C is FANDASTIC!!


 	MATRIX rotation(VECTOR axis,type theta){
 		
 		//normalizz.vett.
 		//il prod.scalare è stato realizzato all'interno della normalizzazione
 		
 		type norm=sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
 		
 		
 		VECTOR assi_norm=(VECTOR)get_block(sizeof(type),3);
 		
 		//fase di normalizzazione
 		
 		for (int i=0;i<3;i++){
 		
 			assi_norm[i]=axis[i]/norm;
 		}//for
 		
 		
 	
 	
 		//usiamo le approssimazioni per calc.i coefficenti
 		
 		type half_theta=theta/2.0;
 		
 		type cos_t=1-(half_theta*half_theta)/2 + (half_theta*half_theta*half_theta*half_theta)/24-(half_theta*half_theta*half_theta*half_theta*half_theta*half_theta)/720;//appr.cos
		
		
		type sen_t=half_theta-(half_theta*half_theta*half_theta)/6+(half_theta*half_theta*half_theta*half_theta*half_theta)/120-(half_theta*half_theta*half_theta*half_theta*half_theta*half_theta*half_theta)/5040;
		
		type a=cos_t;
		
		//definiamo b,c,d
		
		type b= -assi_norm[0]*sen_t;
		
		type c= -assi_norm[1]*sen_t;
		
		type d= -assi_norm[2]*sen_t;
		
		//allochiamo la matrice
		
		MATRIX rot_matrix=alloc_matrix(3,3);
		
		type qa=a*a;
		
		type qb=b*b;
		
		type qc=c*c;
		
		type qd=d*d;
		
		rot_matrix[0]=qa+qb-qc-qd;
		
		rot_matrix[1]=2*(b*c+a*d);
		
		rot_matrix[2]=2*(b*d-a*c);
		
		rot_matrix[3]=2*(b*c-a*d);
		
		rot_matrix[4]=qa+qc-qb-qd;
		
		rot_matrix[5]=2*(c*d+a*b);
		
		rot_matrix[6]=2*(b*d+a*c);
		
		rot_matrix[7]=2*(c*d-a*b);
		
		rot_matrix[8]=qa+qd-qb-qc;
		
		//liberiamo lo spazio allocato
		
		free_block(assi_norm);
		
	return rot_matrix;
 		
 		
 		
 		
 	
 	}//rotation

//fine rotation

void pst(params* input){ //TODO MODIFICARE FACENDO LA FUNZIONE BACKBON A SE
	//da qui inserire poi definizione funzione backbone
	
	
 	
	
	
	
	
}//fine_pst

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);//crea i vettori degli angoli
	input->psi = alloc_matrix(input->N, 1);//insomma crea le matrici monodimensionali vuote
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}
