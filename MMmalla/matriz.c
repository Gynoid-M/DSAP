#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAXBLOQ 100





main(int argc, char **argv)
{
   int myrank, numprocs,numprocvalido = 0,bloqtam,bloqvalido=0,dimension, init = 0,usu = 0,destino;
   int fila,columna,i,k=0; 
   	
	
	

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
			
			double a[bloqtam*bloqtam];
			double b[bloqtam*bloqtam];

			inicializar_matrices(a,b, dimension, bloqtam, fila, columna, i, myrank); 
			 

			printf("Bloques a y b del proceso %d inicializados\n",myrank);
			
		}
		else
		{
			double a[bloqtam*bloqtam];
			double b[bloqtam*bloqtam];
			
			MPI_Recv(&bloqtam,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado); 
			MPI_Recv(&dimension,1,MPI_INT,0,8,MPI_COMM_WORLD,&estado); 
		
			inicializar_matrices(a,b, dimension, bloqtam, fila, columna, i, myrank); 
		
			printf("Bloques a y b del proceso %d inicializados\n",myrank);
			
		}
		
	
	
	
	
		 //necesitamos saber los master de la primera iteracion k = 0 fila + columna
		/*
		 * Dado que hay una rotación a la derecha, los consiguientes master para las demás iteraciones, se calcularán en función del anterior - i (si las dependencias son hacia la 
		 * izquierda) o anterior + i (si las dependencias son hacia la derecha) 
		 * */
	/*	while(k < dimension) //donde k es el número de iteraciones 
		{
			master = myrank/dimension + k;
			if(k > 0)
			{				
				anterior = master - 1;
				if(myrank%dimension < dimension && myrank == myrank/dimension + k || myrank== anterior - dimension) //donde i es de 0...dimension
				{
					if(myrank != myrank%dimension + i)
						MPI_Send(a,bloquetam,MPI_FLOAT,myrank - (myrank%dimension + i),1,MPI_COMM_WORLD); //enviamos a 
						
					
				}
				else
				{
					MPI_Recv(a,bloquetam,MPI_FLOAT,(myrank/dimension + k),1,MPI_COMM_WORLD,estado) ; //enviamos a
				}
			}
			else
			{
				MPI_Send(a,bloquetam,MPI_FLOAT,master - (myrank%dimension + i),1,MPI_COMM_WORLD); 
				
			}  
			                                                                                                                                              
			k++; 
		}*/
	
	
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

void inicializar_matrices(double a[] , double b[], int dimension, int bloqtam, int fila, int columna, int i, int myrank)
{
			int tam = pow(bloqtam,2);
			fila = myrank / dimension;
			columna = myrank % dimension;
			
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
			printf("%d \n ",a[0]);
			
			
			
	
			
}

