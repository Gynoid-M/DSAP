#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAXBLOQ 100





main(int argc, char **argv)
{
   int myrank, numprocs,numprocvalido = 0,bloqtam,bloqvalido=0,dimension, init = 0,usu = 0,destino;
   int fila,columna,i,k=0; 
   double *a;
   double *b;
   double *c;
   int * mifila; 
	
	

   MPI_Status estado; 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   
		
/*
 * Bloque de inicialización de la malla. Declaración y creación de bloques a[i,j] y b[i,j]. 
 * 
 * */		
			
		if(myrank == 0)
		{

		dimension = sqrt(numprocs);
	/*	if(dimension%2 != 0 && dimension%3 != 0) 
		{
			printf("El número de procesos no es un cuadrado perfecto\n");
			printf("Vuelve a lanzar el programa con un número de procesos que sea un cuadrado perfecto (por ejemplo 9)");
			exit(0); 
			
		}*/
			while(bloqvalido == 0)
			{
				printf("Introduce el tamaño de bloques para la matriz\n");
				scanf("%d",&bloqtam);
				if(bloqtam > 0 && bloqtam <=MAXBLOQ)
				{
					bloqvalido = 1;
				}
				else
				{
					printf("El bloque debe ser menor de 100\n");
				}	
				
			}
			
			

			for(destino = 1; destino < numprocs; destino ++)
			{
				MPI_Send(&bloqtam,1,MPI_INT,destino,8,MPI_COMM_WORLD);
				MPI_Send(&dimension,1,MPI_INT, destino,8, MPI_COMM_WORLD); 
			}
			

			usu = 1; //el usuario ya ha introducido los valores y se ha enviado todo al resto de los procesos
			
		}
		 
		if (myrank == 0 && usu == 1)
		{
			
			//inicialización de los bloques, cada uno en su memoria local. 
			//cada bloque corresponde con un proceso, por eso, se debe de calcular su fila y columna particular. 
			
			 a = malloc(bloqtam*bloqtam*sizeof(double));
			b= malloc(bloqtam*bloqtam*sizeof(double));
			c= malloc(bloqtam*bloqtam*sizeof(double));
			mifila=malloc((dimension - 1)*sizeof(int));
			inicializar_matrices(a,b, dimension, bloqtam, fila, columna, i, myrank); 
			 calcular_mifila(mifila,dimension,columna,fila,myrank);

			//printf("Bloques a y b del proceso %d inicializados\n",myrank);
			
		}
		else
		{
			
			
			MPI_Recv(&bloqtam,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado); 
			MPI_Recv(&dimension,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado); 
			 a = malloc(bloqtam*bloqtam*sizeof(double));
			b= malloc(bloqtam*bloqtam*sizeof(double));
			c= malloc(bloqtam*bloqtam*sizeof(double));
			mifila=malloc((dimension - 1)*sizeof(int));	

			fila = myrank / dimension;
			columna = myrank % dimension;

			inicializar_matrices(a,b, dimension, bloqtam, fila, columna, i, myrank); 
			calcular_mifila(mifila,dimension,columna,fila,myrank);
			//printf("Bloques a y b del proceso %d inicializados\n",myrank);
			
		}
		usu = 2; //ya se ha terminado de inicializar los bloques
	
	
	
		 //necesitamos saber los master de la primera iteracion k = 0 fila + columna
		/*
		 * Dado que hay una rotación a la derecha, los consiguientes master para las demás iteraciones, se calcularán en función del anterior - i (si las dependencias son hacia la 
		 * izquierda) o anterior + i (si las dependencias son hacia la derecha) 
		 * */
		int tarea;
		while(k < dimension) //donde k es el número de iteraciones 
		{
			if(columna == (fila + k)% dimension)
			{
				
				for(tarea = 0; tarea < dimension - 1; tarea++)
				{
					MPI_Send(a,pow(bloqtam,2),MPI_DOUBLE,mifila[tarea],8,MPI_COMM_WORLD);
				}
				mult(a,b,c,bloqtam);
			}
			else
			{
				int origen = ((fila + k)% dimension) + fila*dimension;
				printf("Yo soy %d, estoy en la iteracion %d y me tiene que llegar datos de %d\n", myrank,k,origen);
				int * atmp = malloc(bloqtam*bloqtam*sizeof(double));
				MPI_Recv(atmp,pow(bloqtam,2),MPI_DOUBLE,origen,8,MPI_COMM_WORLD,&estado);
				

				mult(atmp,b,c,bloqtam);

				free(atmp);
			}

			k++;

		}

	
	
	MPI_Finalize(); 
 } 
 
 
void mult(double a[], double b[], double *c, int m) 
{
	int i,j,k;
	for (i=0; i<m; i++)
		for (j=0; j<m; j++)
			for (k=0; k<m; k++)
				c[i*m+j]=c[i*m+j]+a[i*m+k]*b[k*m+j];
	return; 
}

void inicializar_matrices(double * a , double * b, int dimension, int bloqtam, int fila, int columna, int i, int myrank)
{
			int tam = pow(bloqtam,2);
			
			
			float pos = pow (((fila*columna)+1),2);
			
			

			int rellenob = 0;
			for(i=0; i < (tam - 1) ; i++)
			{ 
				
				a[i]=i*pos/tam;
				
				if(rellenob == bloqtam)
				{
					b[i/bloqtam] = 1;
				}
				else
					b[i] = 0;
				
				rellenob = 0;
			}
			
			
			
			
	
			
}

void calcular_mifila(int * mifila, int dimension, int columna,int fila, int myrank)
{
	int i, j = 0;
	for (i=0; i<dimension; i++)
	{
		if(i != columna)
		{
			mifila[j] = (fila*dimension) + i;
			j++;
		}
	}
}

