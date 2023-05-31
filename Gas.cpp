/*
Gas de Lennard-Jones utilizando el algoritmo de Verlet

ALGORITMO DE VERLET
Desarrollamos r_i(t+h) y v_i(t+h) en serie de Taylor. Todo son magnitudes vectoriales. Queda:

r_i(t+h) = r_i(t) + h*v_i(t) + h^2*a_i(t)/2
v_i(t+h) = v_i(t) + h*(a_i(t) + a_i(t+h))/2

Los pasos son los mismos:
1) (Solo en la primera iteración) Calcular a_i(t) a partir de la fuerza derivada del potencial y todos los r_i(t)
2) Calcular r_i(t+h) a partir de r_i(t), v_i(t) y a_i(t)
3) Calcular a_i(t+h) usando el potencial LJ
4) Calcular v_i(t+h) usando los datos anteriores


REESCALAMOS LAS VARIABLES
Tomaremos sigma=1, de manera que las distancias vienen en unidades de sigma. Asimismo, eps=kB=1, de manera que la temperatura
se medirá en unidades de 120K (1=120K). Por último, la unidad de masa es la del átomo de argón.


CONDICIONES DE CONTORNO PERIÓDICAS
Básicamente, tendrás un cuadrado 10x10, en el que identificarás cada punto con el de enfrente. Cuando un átomo atraviese una
pared, aparecerá por la contraria. Tener cuidado a la hora de calcular distancias
*/

#include <iostream>
#include <fstream> // Ficheros
#include <cmath>
#include <iomanip> // Formato en los ficheros
#include <random>

// Número de cuerpos, dimensión de la cuadrícula, y masa de las partículas
#define N 16
#define L 10
#define m 1
#define distmin 1.1

// Step temporal
#define h 1e-4
#define iter 2e6

// Pi
#define pi 3.141592

using namespace std;


void inicializar(double posiciones[][2], double velocidades[][2]);
void calcularacelentmash(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double momentosparaP);
void iteracionVerlet(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double& momentosparaP);
void actualizarposiciones(double posiciones[][2], double velocidades[][2], double& momentosparaP);
double distanciayvectordiferencia(double r1[], double r2[], double rdif[]);
void lennjon(double posiciones[][2], double aceleraciones [][2]);
double energiakin(double velocidades[][2]);
double energiapot(double posiciones[][2]);
double rVerlet(double r, double v, double a);
double vVerlet(double v, double a, double a2);

/**************************************************** FUNCIÓN PRINCIPAL ****************************************************/
int main(void) {
    double posiciones[N][2];
    double velocidades[N][2];
    double acelent[N][2];
    double acelentmash[N][2];
    double rdif[2];

    // Variables del sistema
    double energiaK;
    double energiaP;
    double momentosparaP; // Suma de los momentos de cada partícula cuando atravesiesan las paredes
    double tiempo;

    // Inicializo
    inicializar(posiciones, velocidades);

    // Calculo las aceleraciones en 0 y h
    lennjon(posiciones, acelent);
    calcularacelentmash(posiciones, velocidades, acelent, acelentmash, momentosparaP);

    // Abro los ficheros: uno para guardar, otro para el vídeo
    ofstream datos;
    ofstream datospos;
    ofstream datosvel;
    ofstream datosenergiasTyP;

    datos.open("Todo.dat");
    datospos.open("Posiciones.dat");
    datosvel.open("Corrpares.dat");
    datosenergiasTyP.open("EnergiasTyP.dat");

    // Itero en el tiempo
    momentosparaP=0;
    for(int j=0; j<=iter; j++) {
        // Para cada átomo, pego los datos en los ficheros, pero solo cada x
        if(j%1000==0) {
            tiempo=j*h;
            // Escribo el tiempo en el fichero con todo y el de velocidades
            datos << tiempo << "\n";

            for(int i=0; i<N; i++) {
                // El fichero con todo
                for(int k=0; k<2; k++) datos << setw(15) << posiciones[i][k];
                for(int k=0; k<2; k++) datos << setw(15) << velocidades[i][k];
                for(int k=0; k<2; k++) datos << setw(15) << acelent[i][k];
        
                datos << "\n";

                // El fichero de posiciones
                datospos << posiciones[i][0] << "," << posiciones[i][1] << "\n";

                // El de velocidades
                /*if(j*h>=30 && j*h<=60) {
                    datosvel << setw(15) << velocidades[i][0] << setw(15) << velocidades[i][1] << setw(15) <<
                    sqrt(velocidades[i][0]*velocidades[i][0]+velocidades[i][1]*velocidades[i][1]) << "\n";
                }*/

                // El de correlación de pares
                if(i!=0) {
                    datosvel << setw(15) << tiempo << setw(15) <<
                    distanciayvectordiferencia(posiciones[0], posiciones[i], rdif) << "\n";
                }

            }

            datos << "\n";
            datospos << "\n";

            // Escribo los datos de energía, temperatura y presión en su fichero
            // La temperatura es la energía cinética entre el número de partículas, porque k_B=1
            // La presión es momentosparaP entre el tiempo y entre 4L

            energiaK=energiakin(velocidades);
            energiaP=energiapot(posiciones);

            datosenergiasTyP << setw(15) << tiempo << setw(15) <<
            energiaK << setw(15) << energiaP << setw(15) << energiaK+energiaP << setw(15) << energiaK/N <<
            setw(15) << momentosparaP/(tiempo*4*L) << setw(15) << "\n";
        }

        // Caliento el sistema en ciertos instantes
        /*if(j*h==60 || j*h==120 || j*h==180) {
            for(int i=0; i<N; i++) for(int k=0; k<2; k++) velocidades[i][k]=1.5*velocidades[i][k];
        }*/

        // Calculo los nuevos parámetros
        iteracionVerlet(posiciones, velocidades, acelent, acelentmash, momentosparaP);
    }

    datos.close();
    datospos.close();
    datosvel.close();
    datosenergiasTyP.close();



    
    return 0;
}
/***************************************************************************************************************************/


/*Función inicializar. Genera datos aleatorios para las posiciones y direcciones iniciales, y las asigna a los
correspondientes vectores*/
void inicializar(double posiciones[][2], double velocidades[][2]) {

    random_device sem;
    default_random_engine semilla(sem());
    uniform_real_distribution<double> posicionaleat(0.0,L);
    uniform_real_distribution<double> angaleat(0.0,2*pi);

    // Distancias
    double rdif[N][N][2];

    // Si es verdadero, se agrega la partícula
    bool agregar;

    // Ángulos
    double angulos[N];

    /*Asigno valores aleatorios a la posición, entre 0 y L, y a la dirección, entre 0 y 2pi.
    Todos los módulos de la velocidad son 1*/

    // Para cada átomo
    for(int i=0; i<N; i++) {
        // Generar posiciones hasta que la nueva partícula esté a distancia >=distmin de todo el resto
        do{
            posiciones[i][0]=posicionaleat(semilla);
            posiciones[i][1]=posicionaleat(semilla);

            agregar=true;
            for(int j=0; j<i; j++) {

                // Si la distancia entre i y j es menor que la mínima, no lo agrega y lo intenta otra vez
                if(distanciayvectordiferencia(posiciones[i],posiciones[j],rdif[i][j])<distmin) {
                    agregar=false;
                    break;
                }
            }

        }while(!agregar);

        // Velocidades. Así para que si N es par, el v_CM sea 0
        if(i%2==0) {
            angulos[i]=angaleat(semilla);

            velocidades[i][0]=4*cos(angulos[i]);
            velocidades[i][1]=4*sin(angulos[i]);
        }

        else for(int k=0; k<2; k++) velocidades[i][k]=-velocidades[i-1][k];


        // Estado sólido: red cuadrada y parten del reposo
        /*posiciones[i][0]=0.5+i%4;
        posiciones[i][1]=0.5+(i-i%4)/4;

        for(int k=0; k<2; k++) velocidades[i][k]=0;*/
        
    }
    

}

/*Función calcularacelentmash. Calcula las posiciones, velocidades y aceleraciones en el instante h
a partir de las anteriores, en 0*/
void calcularacelentmash(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double momentosparaP) {

    double posicionesenh[N][2];

    // Calculo la posición en el instante h para cada átomo
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            posicionesenh[i][k]=rVerlet(posiciones[i][k], velocidades[i][k], acelent[i][k]);
        }

    }

    actualizarposiciones(posicionesenh, velocidades, momentosparaP);

    // A partir de LJ, calculo la aceleración en t+h
    lennjon(posicionesenh, acelentmash);

    return;
}

/*Función iteracionVerlet. Calcula las posiciones, velocidades y aceleraciones en el instante t+h
a partir de las anteriores, en t*/
void iteracionVerlet(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double& momentosparaP) {

    // Calculo la posición en el instante t+h para cada átomo
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            posiciones[i][k]=rVerlet(posiciones[i][k], velocidades[i][k], acelent[i][k]);
        }

    }

    actualizarposiciones(posiciones, velocidades, momentosparaP);

    // Aceleración en t+h ahora es en t
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            acelent[i][k]=acelentmash[i][k];
        }

    }

    // A partir de LJ, calculo la nueva aceleración en t+h
    lennjon(posiciones, acelentmash);

    // Calculo la velocidad en el instante t+h
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            velocidades[i][k]=vVerlet(velocidades[i][k], acelent[i][k], acelentmash[i][k]);
        }

    }

    return;
}


/*Función actualizarposiciones. Reduce las posiciones (x,y) módulo L. Además, si una partícula atraviesa
una pared, añade una cantidad 2mv_x ó y a la suma momentosparaP*/
void actualizarposiciones(double posiciones[][2], double velocidades[][2], double& momentosparaP) {

    // Para cada partícula
    for(int i=0; i<N; i++) {

        // Cada coordenada
        for(int k=0; k<2; k++) {

            // Si ha cruzado, añade la contribución correspondiente a momentosparaP y reduce
            if(posiciones[i][k]>= L || posiciones[i][k]<0) {
                momentosparaP+=2*m*abs(velocidades[i][k]);
                posiciones[i][k]+=-L*floor(posiciones[i][k]/L);
            }
        }
    }

}


/*Función distanciayvectordiferencia. Dados puntos r1, r2, calcula r1-r2 en este toro raro. Además, devuelve
la distancia entre ellos*/
double distanciayvectordiferencia(double r1[], double r2[], double rdif[]) {

    // Obtengo el cuadrado determinado por la retícula de la partícula 2 en cuyo interior se encuentra la partícula 1
    double incx=r1[0]-r2[0];
    double incxraro;

    if(r1[0]<=r2[0]) incxraro=r1[0]-r2[0]+L;
    else incxraro=r1[0]-r2[0]-L;

    double incy=r1[1]-r2[1];
    double incyraro;

    if(r1[1]<=r2[1]) incyraro=r1[1]-r2[1]+L;
    else incyraro=r1[1]-r2[1]-L;

    // Calculo las posibles distancias al cuadrado, la que sea mínima es la verdadera
    double d1=incx*incx+incy*incy;
    double d2=incxraro*incxraro+incy*incy;
    double d3=incx*incx+incyraro*incyraro;
    double d4=incxraro*incxraro+incyraro*incyraro;

    // Si la distancia usual es la menor, el vector es el usual
    if(d1<=d2 && d1<=d3 && d1<=d4) {
        rdif[0]=incx;
        rdif[1]=incy;
        return sqrt(d1);
    }

    // Si es la segunda, pues el correspondiente, y así
    else if(d2<=d1 && d2<=d3 && d2<=d4) {
        rdif[0]=incxraro;
        rdif[1]=incy;
        return sqrt(d2);
    }

    else if(d3<=d1 && d3<=d2 && d3<=d4) {
        rdif[0]=incx;
        rdif[1]=incyraro;
        return sqrt(d3);
    }

    else {
        rdif[0]=incxraro;
        rdif[1]=incyraro;
        return sqrt(d4);
    }
}

/*Función lennjon. Calcula la aceleración de todos los átomos a partir de las posiciones*/
void lennjon(double posiciones[][2], double aceleraciones [][2]) {
    double distancia[N][N];
    double rdif[N][N][2];

    // Para cada átomo
    for(int i=0; i<N; i++) {

        // Inicializo las aceleraciones a 0
        aceleraciones[i][0]=0; aceleraciones[i][1]=0;

        // Sumo la contribución del átomo j, j distinto de i
        for(int j=0; j<N; j++) {
            if(j==i) continue;
            else {
                // Si j<i, ya ha calculado el vector de j a i
                if(j<i) {
                    distancia[i][j]=distancia[j][i];

                    // k=0 es la componente x, k=1 la componente y
                    for(int k=0; k<2; k++) rdif[i][j][k]=-rdif[j][i][k];
                }

                // Si no, llamo a la función distanciayvectordiferencia que me calcula ambas cosas
                else {
                    distancia[i][j]= distanciayvectordiferencia(posiciones[i], posiciones[j], rdif[i][j]);
                }

                // Añado el correspondiente sumando
                if(distancia[i][j]>=3) continue;
                else {
                    for(int k=0; k<2; k++) {
                        aceleraciones[i][k]+=24*(2-pow(distancia[i][j],6))*rdif[i][j][k]/(m*pow(distancia[i][j],14));

                    }
                }

            }

        }


    }



    return;
}

/*Función energiakin. Calcula la energía cinética del sistema*/
double energiakin(double velocidades[][2]) {
    double energia;

    energia=0;
    for(int i=0; i<N; i++) {

        // Sumo la energía cinética de la masa i
        energia+=m*(velocidades[i][0]*velocidades[i][0]+velocidades[i][1]*velocidades[i][1])/2;
    }

    return energia;
    
}

/*Función energiapot. Calcula la energía potencial del sistema en un conjunto de posiciones*/
double energiapot(double posiciones[][2]) {
    double energia;
    double distanciaalamenos6;
    double rdif[2];

    energia=0;
    for(int i=1; i<N; i++) {

        // Sumo la interacción de la masa j<i con la masa i
        for(int j=0; j<i; j++) {
            distanciaalamenos6 = pow(distanciayvectordiferencia(posiciones[i], posiciones[j], rdif),-6);
            energia+=4*(distanciaalamenos6*distanciaalamenos6-distanciaalamenos6);
        }
        
    }

    return energia;
    
}

/*Función rVerlet. Calcula la nueva posición r a partir de los parámetros r, v y a anteriores. h es el paso temporal*/
double rVerlet(double r, double v, double a) {
    return r+h*v+h*h*a/2;
}


/*Función vVerlet. Calcula la nueva velocidad v a partir de los parámetros v, a anteriores, así como de la aceleración
del mismo paso, a2 h es el paso temporal*/
double vVerlet(double v, double a, double a2) {
    return v+h*(a+a2)/2;
}