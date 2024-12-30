

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
#define MIN(a, b) ((a) < (b) ? (a) : (b)) // dichiaro la funzione minimo

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
* 	In entrambi i casi il candidato dovr  inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente   che le matrici siano in row-major order.
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
*	Se lo si ritiene opportuno,   possibile cambiare la codifica in memoria
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
*	Se lo si ritiene opportuno,   possibile cambiare la codifica in memoria
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

float rama_energy();
//funz rotation ok!


MATRIX rotation(VECTOR axis,type theta){
 		
 		//normalizz.vett.
 		//il prod.scalare è stato realizzato all'interno della normalizzazione
 		
       // printf("axis[0]= %f\n", axis[0]);
        //printf("axis[1]= %f\n", axis[1]);
      //  printf("axis[2]= %f\n", axis[2]);

 		type norm=sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
 		
 		
 		VECTOR assi_norm=(VECTOR)get_block(sizeof(type),3);
 		
 		//fase di normalizzazione
 		
 		for (int i=0;i<3;i++){
 		
 			assi_norm[i]=axis[i]/norm;
            //printf("assi_norm [%d] = %f\n", i, assi_norm[i]);
 		}//for
 		
 		
 	
 	
 		//usiamo le approssimazioni per calc.i coefficenti
 		
 		type half_theta=theta/2.0;
        
        //printf("half_theta = %f \n", half_theta); 		

 		//type cos_t=1-(((half_theta*half_theta)/2) + ((half_theta*half_theta*half_theta*half_theta)/24)-((half_theta*half_theta*half_theta*half_theta*half_theta*half_theta)/720));//appr.cos
		
		
		//type sen_t=half_theta-(((half_theta*half_theta*half_theta)/6)+((half_theta*half_theta*half_theta*half_theta*half_theta)/120)-((half_theta*half_theta*half_theta*half_theta*half_theta*half_theta*half_theta)/5040));
		
        type cos_t = cos(half_theta);
        type sen_t = sin(half_theta);

		type a=cos_t;
		
		//definiamo b,c,d
		
		type b= -assi_norm[0]*sen_t;
		
		type c= -assi_norm[1]*sen_t;
		
		type d= -assi_norm[2]*sen_t;
		

      /*  printf("a = %f \n", a);
        printf("b = %f \n", b);
        printf("c = %f \n", c);
        printf("d = %f \n", d);
        printf("sen_t = %f\n", sen_t);
        printf("cos_t = %f\n", cos_t);
        */
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
        int i;
        /*for (i=0; i< 9; i++){
            printf("rot_matrix[%d] = %f \n", i, rot_matrix[i]);
        }*/
		
		free_block(assi_norm);
		
	return rot_matrix;
 		
 		
 		
 		
 	
 	}//rotation


VECTOR backbone(params* input,VECTOR phi,VECTOR psi){

   // MATRIX rot=alloc_matrix(3,3);
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
	
	
//	VECTOR v1,v2,v3;

	VECTOR v1 = (VECTOR)get_block(sizeof(type), 3);
	VECTOR v2 = (VECTOR)get_block(sizeof(type), 3);
	VECTOR v3 = (VECTOR)get_block(sizeof(type), 3);

	
	MATRIX rot;
	
	VECTOR newv;
	newv = alloc_matrix(3,1);
	//costruire la catena
	
	for(i=0;i<n;i++){
	    
 		int idx=i*3;//indice base ammin.
		
		if(i>0){
			
			v1[0] = coords[((idx-1)*3) + 0] - coords[((idx-2)*3) + 0]; // Componente x
            v1[1] = coords[((idx-1)*3) + 1] - coords[((idx-2)*3) + 1]; // Componente y
            v1[2] = coords[((idx-1)*3) + 2] - coords[((idx-2)*3) + 2]; // Componente z
	        


	        type norm ;
			
		    norm = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		    

            for (j = 0; j < 3; j++) {
                v1[j] /= norm;
            }


			
			//calcolo v1/(|v1|)
			
			rot=rotation(v1,theta2);
			
			newv[0]=rot[3]*dist_c_n;
			newv[1]=rot[4]*dist_c_n;
			newv[2]=rot[5]*dist_c_n;

//            printf("newv dopo v1 newv[0] = %f, newv[1]=%f, newv[2]=%f\n", newv[0], newv[1], newv[2]);
			//newv={0,dist_c_n,0}*rot;//prod matr penso non corretto,va acceduto diversamente in memoria
			
			coords[(idx * 3) + 0] = coords[((idx - 1) * 3) + 0] + newv[0];
			coords[(idx * 3) + 1] = coords[((idx - 1) * 3) + 1] + newv[1];
			coords[(idx * 3) + 2] = coords[((idx - 1) * 3) + 2] + newv[2];


			//posiziono c alfa usando phi(vettore angolo phi)
		
			v2[0] = coords[((idx)*3) + 0] - coords[((idx-1)*3) + 0]; // Componente x
            v2[1] = coords[((idx)*3) + 1] - coords[((idx-1)*3) + 1]; // Componente y
            v2[2] = coords[((idx)*3) + 2] - coords[((idx-1)*3) + 2]; // Componente z
            
           // printf("newv dopo v2 newv[0] = %f, newv[1]=%f, newv[2]=%f\n", newv[0], newv[1], newv[2]);
			
			norm = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);

            for (j = 0; j < 3; j++) {
                v2[j] /= norm;
            }
			
			rot=rotation(v2,phi[i]);

			newv[0]=rot[3]*dist_ca_n;
			newv[1]=rot[4]*dist_ca_n;
			newv[2]=rot[5]*dist_ca_n;
            
         //   printf("newv dopo v3 newv[0] = %f, newv[1]=%f, newv[2]=%f\n", newv[0], newv[1], newv[2]);
			
			//newv={0,dist_ca_n,0}*rot;//prodotto matr.//rivedere penso non corretto va acceduto diversamente in mem
			
			coords[((idx + 1) * 3) + 0] = coords[(idx * 3) + 0] + newv[0];
			coords[((idx + 1) * 3) + 1] = coords[(idx * 3) + 1] + newv[1];
			coords[((idx + 1) * 3) + 2] = coords[(idx * 3) + 2] + newv[2];

		    
		
			 	

		
		}//if
			
			
			//Posiziono c usando il vettore degli angoli psi
			
	 	v3[0] = coords[((idx+1)*3) + 0] - coords[((idx)*3) + 0]; // Componente x
        v3[1] = coords[((idx+1)*3) + 1] - coords[((idx)*3) + 1]; // Componente y
        v3[2] = coords[((idx+1)*3) + 2] - coords[((idx)*3) + 2]; // Componente z


			
		type norm = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
            for (int j = 0; j < 3; j++) {
                v3[j] /= norm;
            }

	 	rot=rotation(v3,phi[i]);

		newv[0]=rot[3]*dist_ca_c;
		newv[1]=rot[4]*dist_ca_c;
		newv[2]=rot[5]*dist_ca_c;

		//newv={0,dist_ca_c,0}*rot;//rivedere 
		
	    coords[((idx + 2) * 3) + 0] = coords[((idx + 1) * 3) + 0] + newv[0];
		coords[((idx + 2) * 3) + 1] = coords[((idx + 1) * 3) + 1] + newv[1];
		coords[((idx + 2) * 3) + 2] = coords[((idx + 1) * 3) + 2] + newv[2];

	
	
	    
	
	}//for
    /*
    for (int i=0; i<n*3; i++){
        printf("coords[%d] = %f \n", i,coords[i]);
    }   
*/

    free_block(newv);
	free_block(v1);
	free_block(v2);
	free_block(v3);

    return coords;
}//backbone



//sottofunzioni di energy

float rama_energy(params* input,VECTOR phi,VECTOR psi){

	int i;

	int n=input->N;//lunghezza

	float alfa_phi=-57.8;

	float alfa_psi=-47.0;

	float beta_phi=-119.0;

	float beta_psi=113.0;

	float energy=0;

	for(i=0;i<n;i++){

		float alfa_dist=sqrt( ( (phi[i]-alfa_phi)*(phi[i]-alfa_phi) )+( (psi[i]-alfa_psi)*(psi[i]-alfa_psi) ) );

		float beta_dist=sqrt( ( (phi[i]-beta_phi)*(phi[i]-beta_phi) )+( (psi[i]-beta_psi)*(psi[i]-beta_psi) ) );

		energy+=0.5*MIN(alfa_dist,beta_dist); 


	}//for

return energy;

}//rama


//funz. hydro energy


float hydrophobic_energy(params* input, VECTOR ca_coords){

	char* se=input->seq;

	int si,sj;

	int n=input->N;

	int i,j=0;

	int k=0;

	float energy=0;
	
	

	for(i=0;i<n-1;i++){//MESSO A N-1  
		
		for(j=i+1;j<n;j++){
		    
		    float dist = sqrt( ( (ca_coords[i*3+0]-ca_coords[j*3+0]) * (ca_coords[i*3+0]-ca_coords[j*3+0])) +
		    ((ca_coords[i*3+1]-ca_coords[j*3+1]) * (ca_coords[i*3+1]-ca_coords[j*3+1])) +
		    ((ca_coords[i*3+2]-ca_coords[j*3+2]) * (ca_coords[i*3+2]-ca_coords[j*3+2]))
		    );
		//	float dist=abs(ca_coords[i]-ca_coords[j]);
			//usando la formula di dist. euclidea consideriamo solo i valori coords_alfa[i] e coords_alfa[j]
			//e facciamo il modulo direttamente usando la formula senza fare la radice della differenza dei quadrati

			if(dist<10.0){

				si=(int)se[i]-65;

				sj=(int)se[j]-65;

				//valore ascii convertito per prendere le posizioni relative agli amminoacidi

				//65 lettera A
				if (hydrophobicity[si] != -1 && hydrophobicity[sj] != -1) {

				    energy+=(hydrophobicity[si]*hydrophobicity[sj])/dist;
				}
			}//if

		
		}//for
	
	
	}//for

	

return energy;

}//hydro



float elecrostatic_energy(params* input,VECTOR ca_coords){

	char* se=input->seq;

	int si,sj;

	int n=input->N;

	int i,j=0;

	int k=0;

	float energy=0;
	
	
	
	for(i=0;i<n-1;i++){//MESSO A N-1  
		
		for(j=i+1;j<n;j++){
		    
		    
		    float dist = sqrt( ( (ca_coords[i*3+0]-ca_coords[j*3+0]) * (ca_coords[i*3+0]-ca_coords[j*3+0])) +
		    ((ca_coords[i*3+1]-ca_coords[j*3+1]) * (ca_coords[i*3+1]-ca_coords[j*3+1])) +
		    ((ca_coords[i*3+2]-ca_coords[j*3+2]) * (ca_coords[i*3+2]-ca_coords[j*3+2]))
		    );
            
           // printf("dist = %f\n", dist);

			//float dist=abs(ca_coords[i]-ca_coords[j]);
			//usando la formula di dist. euclidea consideriamo solo i valori coords_alfa[i] e coords_alfa[j]
			//e facciamo il modulo direttamente usando la formula senza fare la radice della differenza dei quadrati
				si=(int)se[i]-65;

				sj=(int)se[j]-65;
           // printf("si = %d\n", si);
//            printf("sj = %d\n", sj);

            //if (charge[si] == -1 || charge[sj] == -1) break;

			if(dist<10.0 && charge[si]!=0 && charge[sj]!=0 && charge[si] != -1 && charge[sj] != -1){

				energy+=charge[si]*charge[sj]/(dist*4.0);
              //  printf("charge[%d] = %f", si, charge[si]);

			}//if	

				//valore ascii convertito per prendere le posizioni relative agli amminoacidi

				//65 lettera A


		
		}//for
	


	}//for


return energy;

}//electr

//portare dentro energy vettore da riutilizzare



float packing_energy(params* input,VECTOR ca_coords){

	char* se=input->seq;

	int si,sj;

	int n=input->N;

	int i,j=0;

	int k=0;

	float energy=0;
	
	

	for(i=0;i<n;i++){//MESSO A N-1  

		int density=0;
		
		for(j=0;j<n;j++){

    
    
            float dist = sqrt( ( (ca_coords[i*3+0]-ca_coords[j*3+0]) * (ca_coords[i*3+0]-ca_coords[j*3+0])) +
		    ((ca_coords[i*3+1]-ca_coords[j*3+1]) * (ca_coords[i*3+1]-ca_coords[j*3+1])) +
		    ((ca_coords[i*3+2]-ca_coords[j*3+2]) * (ca_coords[i*3+2]-ca_coords[j*3+2]))
		    );


			//float dist=abs(ca_coords[i]-ca_coords[j]);
			//usando la formula di dist. euclidea consideriamo solo i valori coords_alfa[i] e coords_alfa[j]
			//e facciamo il modulo direttamente usando la formula senza fare la radice della differenza dei quadrati
			si=(int)se[i]-65;

			sj=(int)se[j]-65;


            //if (volume[si] == -1 || volume[sj] == -1) break;

			if(i != j && dist<10.0 && volume[si] != -1 && volume[sj] != -1 ){

				density+=(volume[sj])/(dist*dist*dist);

			}//if	

				//valore ascii convertito per prendere le posizioni relative agli amminoacidi

				//65 lettera A


		
		}//for

		energy+=( (volume[si]-density) *  (volume[si]-density) );
	


	}//for


return energy;



}//pack


float energy(params* input){

	VECTOR phi= input->phi;	
	VECTOR psi= input->psi;
	MATRIX coords=backbone(input,phi,psi); 

	int n = input->N;
	VECTOR ca_coords=(VECTOR)get_block(sizeof(type),n*3);//dim n


    int i;

    for (i = 0; i < n; i++) { // Itera per tutti i residui
        
        ca_coords[i * 3 + 0] = coords[(i * 3 + 1) * 3 + 0]; // x del Cα
        ca_coords[i * 3 + 1] = coords[(i * 3 + 1) * 3 + 1]; // y del Cα
        ca_coords[i * 3 + 2] = coords[(i * 3 + 1) * 3 + 2]; // z del Cα
       
    }
    
    
  /*  for(i=0; i<n; i++){
        printf("ca_coords[%d] = %f \n", i, ca_coords[i]);
    }  */
    




	//richiamo sottofunzioni


	float rama=rama_energy(input,phi,psi);

	float hydro=hydrophobic_energy(input,ca_coords);

	float elec=elecrostatic_energy(input,ca_coords);

	float pack=packing_energy(input,ca_coords);

	//pesi per i contributi


	float omega_rama=1.0;


	float omega_hydro=0.5;



	float omega_elec=0.2;


	float omega_pack=0.3;

	float total_energy=omega_rama*rama + omega_hydro*hydro + omega_elec*elec + omega_pack*pack;

    //printf("Energie rama = %f , hydro = %f, elec = %f, pack = %f \n", rama, hydro, elec, pack);
    //printf("total energy = %f\n", total_energy);
    
	free_block(ca_coords);

	free_block(coords);

return total_energy;


}//energy


//fine rotation

void pst(params* input){ //TODO MODIFICARE FACENDO LA FUNZIONE BACKBON A SE
	//da qui inserire poi definizione funzione backbone
	
	char* s = input->seq;
 	int n = input->N;
	type T0 = input->to;
	type T = T0;
	type alpha = input->alpha;
	type k = input->k;

	VECTOR phi= input->phi;	
	VECTOR psi= input->psi;
	
	//inizializzarli random 
	int i;
	for(i=0; i<n;i++){
		phi[i]=(random()*2 * M_PI) - M_PI;
		psi[i]=(random()*2 * M_PI) - M_PI;
      //  printf("psi[%d] = %f", i,psi[i]);
     //   printf("phi[%d] = %f", i,phi[i]);
	}
   

	float e = energy(input);
	int t = 0;

	while(T>0){
		i=0;
		i=(int)(random() * (n+1));

		type y_phi = (random()*2 * M_PI) - M_PI;
		phi[i]+=y_phi;

		type y_psi = (random()*2 * M_PI) - M_PI;
		psi[i]+=y_psi;

		float newE = energy(input);

		type delta = newE-e;

		if(delta<=0){
			e = newE;
		}
		else{
			type P = exp(-delta/(k*T));
			type r = random();

			if(r<=P){
				e = newE;
                
			}
            
			else{
				phi[i]-=y_phi;
				psi[i]-=y_psi;
			}
		}

		t +=1;
		T = T0 - sqrt(alpha*t);//da controllare

    printf("energy = %f\n", e);
	}
	
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
