#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


typedef struct _PGMData {
    int row;
    int col;
    int max_gray;
    int **matrix;
} PGMData;

typedef struct _MinMax{
   double Min;
   double Max;
} MinMax;

int **CrearArray2D_int_consecutivo(int , int );
double **CrearArray2D_double(int , int );
void LiberarArray2D_double(int, double **);
MinMax BuscarMinMax(int **, int , int );

#define ITERMOD 500

int main(int argc, char **argv)
{
 char archivo_imagen_ori[20] = "logo.pgm";
 char archivo_imagen_aristas[20] = "logo.edge.pgm";
 char archivo_imagen_restaurada[20] = "logo.out.pgm";
 char aux[sizeof(int)];
 int k, i, j = 0, IterMax=1000,myrank, numprocs,filasproc,filasroot;
 double Tinicial, Tfinal, Ttotal;

 PGMData pgm_data, pgm_mod, pgm_aux;
 PGMData *img_data = &pgm_data;
 PGMData *img_mod = &pgm_mod; //esta es buff
 PGMData *img_aux = &pgm_aux;
 int **img_edge;
 double **img_old, **img_new;
 double norm2, cota;
    
 void readPGM(char *, PGMData *);
 void writePGM(char *, PGMData *);
 void Filtro_Laplace(PGMData *, PGMData *);
 void Escalar_Imagen(PGMData *, PGMData *);
 void Umbral_Imagen(PGMData *, PGMData *);
 double ctimer(void);

 MPI_Status estado;
    
 MPI_Request request;
 
 MPI_Init(&argc,&argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

if(myrank == 0)
{

 switch (argc)
        {
                case 3: sscanf(argv[2], "%i", &IterMax);
                case 2: strcpy(archivo_imagen_ori, argv[1]);
                        strcat(archivo_imagen_ori, ".pgm");
                        strcpy(archivo_imagen_aristas, argv[1]);
                        strcat(archivo_imagen_aristas, ".edge.pgm");
                        strcpy(archivo_imagen_restaurada, argv[1]);
                        strcat(archivo_imagen_restaurada, ".out.pgm");
                        printf("%s\n", archivo_imagen_restaurada);
                case 1: break;
                default: printf("Demasiados parametros\n");
                         return 0;
        }

 printf("\n  *************** DATOS DE LA EJECUCION ****************************\n");
 printf("  * Archivo imagen original   : %20s               *\n", archivo_imagen_ori);
 printf("  * Archivo imagen con aristas: %20s               *\n", archivo_imagen_aristas);
 printf("  * Archivo imagen restaurada : %20s               *\n",archivo_imagen_restaurada);  
 printf("  * Iteraciones = %5d                                            *\n",IterMax);
 printf("  ******************************************************************\n\n");
 printf("Leyendo imagen \"%s\" ... ", archivo_imagen_ori);

 readPGM(archivo_imagen_ori, img_data);
 printf("Dimension de la imagen: %d x %d\n", img_data->row, img_data->col);
 filasproc = img_data->row/numprocs;
 if(img_data->row%numprocs>0) filasroot = filasproc + img_data->row%numprocs;
 else filasroot = filasproc;
 
 img_old = CrearArray2D_double(filasroot + 2, img_data->col + 2);
 img_new = CrearArray2D_double(filasroot + 2, img_data->col + 2);
 img_edge = CrearArray2D_int_consecutivo(filasroot + 2, img_data->col + 2);
 Filtro_Laplace(img_data, img_mod);
 free(img_data->matrix);
 //Umbral_Imagen(img_mod, img_aux);
 writePGM(archivo_imagen_aristas, img_aux);
 free(img_aux->matrix);
 Tinicial = ctimer();
 Bcast(); //para enviar datos adicionales que necesitarán el resto de los procesos
}
else
{
 img_old = CrearArray2D_double(filasproc + 2, img_data->col + 2);
 img_new = CrearArray2D_double(filasproc + 2, img_data->col + 2);
 img_edge = CrearArray2D_int_consecutivo(filasproc + 2, img_data->col + 2);
}
/*
 * Determinar la dimensión asignada al proceso que invoca la operación. 
 */ 
 for (i=0;i<img_data->row+2;i++) {
     for (j=0;j<img_data->col+2;j++) {
        if ((i==0) || (i==img_data->row+1) || (j==0) || (j==img_data->col+1)) img_edge[i][j] = 0;
        else img_edge[i][j] = img_mod->matrix[i-1][j-1];
     }
 }

 for (i=0;i<img_data->row+2;i++) {
     for (j=0;j<img_data->col+2;j++) {
          img_old[i][j] = img_edge[i][j];
     }
 }

// for (k=0;k<IterMax;k++) {
 k=0;
 cota = sqrt((img_data->row)*(img_data->col));
 norm2 = cota +1;
 while ((norm2 > cota) && (k < IterMax)) {
    k=k+1;
    norm2=0;
    for (i=1;i<img_data->row+1;i++) {
        for (j=1;j<img_data->col+1;j++) {
            img_new[i][j] = 0.25* (img_old[i-1][j] + 
                                   img_old[i+1][j] +
                                   img_old[i][j-1] +
                                   img_old[i][j+1]   -
                                   img_edge[i][j]);
            norm2= norm2 + (img_new[i][j]-img_old[i][j])*(img_new[i][j]-img_old[i][j]);
        }
    }
    norm2 = sqrt(norm2);
    for (i=1;i<img_data->row+1;i++) {
        for (j=1;j<img_data->col+1;j++) {
            img_old[i][j] = img_new[i][j];
        }
    }
//   norm2=cota+1;
    if (k%ITERMOD == 0) printf("Iteracion %d, norma = %f\n",k, norm2);
 }

 for (i=1;i<img_data->row+1;i++) {
     for (j=1;j<img_data->col+1;j++) {
          img_mod->matrix[i-1][j-1] = img_new[i][j];
     }
 }

 Tfinal = ctimer();
 Ttotal = Tfinal-Tinicial;
 printf("\nIteracion %d, norma = %f\n",k, norm2);
 printf("Tiempo: %f segundos\n", Ttotal);
 Escalar_Imagen(img_mod, img_mod); 
 writePGM(archivo_imagen_restaurada, img_mod);

 free(img_mod->matrix);
 free(img_edge);
 LiberarArray2D_double(img_data->row + 2, img_old);
 LiberarArray2D_double(img_data->row + 2, img_new);
 
 MPI_Finalize(); 
}

#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)

void SkipComments(FILE *fp)
{
    int ch;
    char line[100];
 
    while ((ch = fgetc(fp)) != EOF && isspace(ch))
        ;
    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        SkipComments(fp);
    } else
        fseek(fp, -1, SEEK_CUR);
}

void readPGM(char *file_name, PGMData *data)
{
    FILE *pgmFile;
    char version[3];
    int i, j;
    int lo, hi;

    pgmFile = fopen(file_name, "rb");
    if (pgmFile == NULL) {
        perror("No se puede abrir el archivo para leer");
        exit(EXIT_FAILURE);
    }

    fgets(version, sizeof(version), pgmFile);
    if (strcmp(version, "P5")) {
        fprintf(stderr, "Tipo de archivo incorrecto!\n");
        exit(EXIT_FAILURE);
    }

    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->col);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->row);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->max_gray);
    fgetc(pgmFile);
 
    data->matrix = CrearArray2D_int_consecutivo(data->row, data->col);
    if (data->max_gray > 255)
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                data->matrix[i][j] = (hi << 8) + lo;
            }
    else
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = fgetc(pgmFile);
                data->matrix[i][j] = lo;
            }
 
    fclose(pgmFile);
    return;
 
}
 
void writePGM(char *filename, PGMData *data)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;
 
    pgmFile = fopen(filename, "wb");
    if (pgmFile == NULL) {
        perror("No se puede abrir el archivo para escribir");
        exit(EXIT_FAILURE);
    }
 
    fprintf(pgmFile, "P5 ");
    fprintf(pgmFile, "%d %d ", data->col, data->row);
    fprintf(pgmFile, "%d ", data->max_gray);
 
    if (data->max_gray > 255) {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                hi = HI(data->matrix[i][j]);
                lo = LO(data->matrix[i][j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }
 
        }
    } else {
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = LO(data->matrix[i][j]);
                fputc(lo, pgmFile);
            }
    }
 
    fclose(pgmFile);
}

double **CrearArray2D_double(int Filas_Cantidad, int Columnas_Cantidad)
{
    double **Array2D;
    int fila = 0, columna = 0;
//  Asignamos memoria para el array
    Array2D = (double **)malloc(Filas_Cantidad * sizeof(double *));
    if (Array2D == NULL)
    {
        perror("Problemas al asignar memoria");
        exit(EXIT_FAILURE);
    }
//     Asignamos memoria para cada una de las filas
    for (fila = 0; fila < Filas_Cantidad; fila++)
    {
        Array2D[fila] = (double *)malloc(Columnas_Cantidad * sizeof(double));
        if (Array2D[fila] == NULL)
        {
        	perror("Problemas al asignar memoria");
       	exit(EXIT_FAILURE);
        }
    }
    return Array2D;
}


void LiberarArray2D_double(int Filas_Cantidad, double **Array2D)
{
    int fila = 0, columna = 0; 
    /* Liberamos la memoria asignada a cada una de las filas */
    for ( fila = 0; fila < Filas_Cantidad; fila++ )
        free(Array2D[fila]);
    /* Liberamos la memoria asignada al array de punteros */
    free(Array2D);
    return;
}


MinMax BuscarMinMax(int **array, int fil, int col) {
   MinMax  min_max;
   int index,fila,columna,fila_next,columna_next;
   int n = fil*col; 
   if ( n%2 != 0 ){

     min_max.Min = array[0][0];
     min_max.Max = array[0][0];

     index = 1;
   }
   else{
     if ( array[0][0] < array[0][1] ){
       min_max.Min = array[0][0];
       min_max.Max = array[0][1];
      }
      else {
       min_max.Min = array[0][1];
       min_max.Max = array[0][0];
      }
       index = 2;
   }

   int big, small,i;
   for ( i = index; i < n-1; i = i+2 ){
      fila = i / col;
      columna = i % col;  
      fila_next = fila;
      columna_next = columna + 1;
      if (columna_next == col) {
         fila_next = fila +1;
         columna_next = 0;
      }
      if ( array[fila][columna] < array[fila_next][columna_next] ){ 
        small = array[fila][columna];
        big = array[fila_next][columna_next];
      }
      else{
        small = array[fila_next][columna_next];
        big = array[fila][columna];
      }
      if ( min_max.Min > small ){ 
        min_max.Min = small;
      }
      if ( min_max.Max < big ){ 
        min_max.Max = big;
      }
   }
   printf("Minimo = %f, Maximo = %f\n", min_max.Min, min_max.Max);
   return min_max;
}

void Filtro_Laplace(PGMData *img_data, PGMData *img_mod)
{
 int i,j;

 img_mod->row =  img_data->row;
 img_mod->col =  img_data->col;
 img_mod->max_gray =  img_data->max_gray;
 img_mod->matrix = CrearArray2D_int_consecutivo(img_data->row ,img_data->col);

 for (i=1; i<img_data->row -1 ;i++) {
     for (j=1; j<img_data->col -1; j++){
        img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i+1][j] +
                                img_data->matrix[i][j-1] + img_data->matrix[i][j+1] -
                               4*img_data->matrix[i][j];
     }
 }
 for (j=1; j<img_data->col-1;j++) {
     i=0;
     img_mod->matrix[i][j] = img_data->matrix[1][j] + img_data->matrix[i+1][j] +
                             img_data->matrix[i][j-1] + img_data->matrix[i][j+1] -
                           4*img_data->matrix[i][j];
     i=img_data->row-1;
     img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i][j] +
                             img_data->matrix[i][j-1] + img_data->matrix[i][j+1] -
                           4*img_data->matrix[i][j];
 }
 for (i=1; i<img_data->row-1;i++) {
     j=0;
     img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i+1][j] +
                             img_data->matrix[i][j] + img_data->matrix[i][j+1] -
                           4*img_data->matrix[i][j];
     j=img_data->col-1;
     img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i+1][j] +
                             img_data->matrix[i][j-1] + img_data->matrix[i][j] -
                           4*img_data->matrix[i][j];
 }
 i=0;
 j=0;
 img_mod->matrix[i][j] = img_data->matrix[i][j] + img_data->matrix[i+1][j] +
                         img_data->matrix[i][j] + img_data->matrix[i][j+1] -
                       4*img_data->matrix[i][j];
 i=0;
 j=img_data->col-1;
 img_mod->matrix[i][j] = img_data->matrix[i][j] + img_data->matrix[i+1][j] +
                         img_data->matrix[i][j-1] + img_data->matrix[i][j] -
                       4*img_data->matrix[i][j];
 i=img_data->row-1;
 j=0;
 img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i][j] +
                         img_data->matrix[i][j] + img_data->matrix[i][j+1] -
                       4*img_data->matrix[i][j];
 i=img_data->row-1;
 j=img_data->col-1;
 img_mod->matrix[i][j] = img_data->matrix[i-1][j] + img_data->matrix[i][j] +
                         img_data->matrix[i][j-1] + img_data->matrix[i][j] -
                       4*img_data->matrix[i][j];
 return; 
}

void Escalar_Imagen(PGMData *img_in, PGMData *img_out)
{
 MinMax BuscarMinMax(int **,int, int);
 int i,j;
 MinMax min_max_img_in;
 img_out->row =  img_in->row;
 img_out->col =  img_in->col;
 img_out->max_gray =  img_in->max_gray;
 printf("Antes de escalar: ");
 min_max_img_in = BuscarMinMax(img_in->matrix,img_in->row, img_in->col);
 for (i=0;i<img_in->row;i++) {
      for (j=0;j<img_in->col;j++) {
           img_out->matrix[i][j] = img_in->matrix[i][j] - min_max_img_in.Min;
           img_out->matrix[i][j] = img_out->matrix[i][j] * (255.0/(min_max_img_in.Max-min_max_img_in.Min));
      }
 }
 printf("Despues de escalar: ");
 min_max_img_in = BuscarMinMax(img_out->matrix,img_out->row, img_out->col);
 return;
}


void Umbral_Imagen(PGMData *img_in, PGMData *img_out)
{
 int i,j;
 img_out->row =  img_in->row;
 img_out->col =  img_in->col;
 img_out->max_gray =  img_in->max_gray;
 img_out->matrix = CrearArray2D_int_consecutivo(img_in->row ,img_in->col);
 for (i=0;i<img_in->row;i++) {
      for (j=0;j<img_in->col;j++) {
             img_out->matrix[i][j]=img_in->matrix[i][j];
             if (img_out->matrix[i][j]<0) img_out->matrix[i][j] = 0;
             if (img_out->matrix[i][j]>255) img_out->matrix[i][j] = 255;
      }
 }
 return;
}

int **CrearArray2D_int_consecutivo(int Filas, int Columnas)
{
// crea un array de 2 dimensiones en posiciones contiguas de memoria
 int *mem_matriz;
 int **matriz;
 int fila, col;
 if (Filas <=0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return;
    }
 if (Columnas <= 0)
    {
        printf("El numero de filas debe ser mayor que cero\n");
        return;
    }
 mem_matriz = malloc(Filas * Columnas * sizeof(int));
 if (mem_matriz == NULL)
        {
                printf("Insuficiente espacio de memoria\n");
                return;
        }
 matriz = malloc(Filas * sizeof(int *));
 if (matriz == NULL)
        {
                printf ("Insuficiente espacio de memoria\n");
                return;
        }
 for (fila=0; fila<Filas; fila++)
    matriz[fila] = mem_matriz + (fila*Columnas);
 return matriz;
}

#include <sys/time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

static int nclock;
double ctimer(void)
{
   struct timeval tp;
   struct timezone tzp;
   double diff;
   nclock=sysconf(_SC_CLK_TCK);
   gettimeofday(&tp, &tzp);
   diff=(double)tp.tv_sec+(double)tp.tv_usec/1.0e6;
   return diff;
}

