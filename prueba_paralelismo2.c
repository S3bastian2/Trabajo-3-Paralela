#include <stdio.h>
#include <omp.h>

int binary_search(int arreglo[], int tamaño_arreglo, int valor_buscado) {
    /* Versión simplificada solicitada: recorrer secuencialmente y devolver
       el índice del primer elemento estrictamente mayor que valor_buscado. */
    for (int i = 0; i < tamaño_arreglo; ++i) {
        if (arreglo[i] > valor_buscado) {
            return i;
        }
    }
    return -1; // no existe ningún elemento mayor
}

void sequential_merge(int arregloA[], int arregloB[], int arregloC[], int tamaño_arregloA, int tamaño_arregloB, int inicio_A, int inicio_B, int inicio_C, int tiene_limite, int valor_limite) {
    int indice_arregloA = inicio_A - 1;
    int indice_arregloB = inicio_B - 1;
    int indice_arregloC = inicio_C - 1;

    while (indice_arregloA < tamaño_arregloA && indice_arregloB < tamaño_arregloB) {
        if (tiene_limite == 1) {
            if (arregloA[indice_arregloA] >= valor_limite && arregloB[indice_arregloB] >= valor_limite) {
                break;
            }
        }
        if (arregloA[indice_arregloA] <= arregloB[indice_arregloB]) {
            arregloC[indice_arregloC] = arregloA[indice_arregloA];
            indice_arregloA++;
        } else {
            arregloC[indice_arregloC] = arregloB[indice_arregloB];
            indice_arregloB++;
        }
        indice_arregloC++;
    }
    if (tiene_limite == 0) {
        while (indice_arregloA < tamaño_arregloA) {
            arregloC[indice_arregloC] = arregloA[indice_arregloA];
            indice_arregloA++;
            indice_arregloC++;
        }
        while (indice_arregloB < tamaño_arregloB) {
            arregloC[indice_arregloC] = arregloB[indice_arregloB];
            indice_arregloB++;
            indice_arregloC++;
        }
    }   
}

void crew_merge(int arregloA[], int arregloB[], int arregloC[], int tamaño_arregloA, int tamaño_arregloB, int num_procesadores) {
    struct Tripleta
    {
        int valor;
        int posicion_original;
        char origen; // Salio de A o de B.
    };
    int tamaño_arreglos = num_procesadores - 1;
    int subArregloA[tamaño_arreglos];
    int subArregloB[tamaño_arreglos];
    struct Tripleta arreglo_v[(2 * num_procesadores) - 2];
    int valor_solicitado, i, x, y, valor_posicion, valor_pivote, posicion_C, valor_limite, indice_v;

    omp_set_num_threads(num_procesadores);

    #pragma omp parallel for private(i) shared(arregloA, arregloB, subArregloA, subArregloB, tamaño_arregloA, tamaño_arregloB, num_procesadores)
    for (i = 1; i <= num_procesadores - 1; i++) {
        int indiceA = i * ((tamaño_arregloA + num_procesadores - 1) / num_procesadores);
        if (indiceA > tamaño_arregloA) {
            indiceA = tamaño_arregloA;
        }    
        subArregloA[i] = arregloA[indiceA - 1];

        int indiceB = i * ((tamaño_arregloB + num_procesadores -1) / num_procesadores);
        if (indiceB > tamaño_arregloB) {
            indiceB = tamaño_arregloB;
        }
        subArregloB[i] = arregloB[indiceB - 1];
    }

    #pragma omp parallel for private(i, valor_solicitado) shared(subArregloA, subArregloB, arreglo_v, tamaño_arreglos, num_procesadores)
    for (i = 1; i < num_procesadores - 1; i++) {
        valor_solicitado = binary_search(subArregloB, tamaño_arreglos, subArregloA[i]);
        if (valor_solicitado >= 1 && valor_solicitado <= num_procesadores - 1) {
            int posicion_v = i + valor_solicitado - 1;
            arreglo_v[posicion_v].valor = subArregloA[i];
            arreglo_v[posicion_v].posicion_original = i;
            arreglo_v[posicion_v].origen = 'A';
        } else {
            int posicion_v = i + num_procesadores - 1; 
            arreglo_v[posicion_v].valor = subArregloA[i];
            arreglo_v[posicion_v].posicion_original = i;
            arreglo_v[posicion_v].origen = 'A';
        }
    }
    #pragma omp parallel for private(i, valor_solicitado) shared(subArregloA, subArregloB, arreglo_v, tamaño_arreglos, num_procesadores)
    for (i = 1; i < num_procesadores - 1; i++) {
        valor_solicitado = binary_search(subArregloA, tamaño_arreglos, subArregloB[i]);
        if (valor_solicitado >= 1 && valor_solicitado <= num_procesadores - 1) {
            int posicion_v = i + valor_solicitado - 1;
            arreglo_v[posicion_v].valor = subArregloB[i];
            arreglo_v[posicion_v].posicion_original = i;
            arreglo_v[posicion_v].origen = 'B';
        } else {
            int posicion_v = i + num_procesadores - 1; 
            arreglo_v[posicion_v].valor = subArregloB[i];
            arreglo_v[posicion_v].posicion_original = i;
            arreglo_v[posicion_v].origen = 'B';
        }
    }

    int Q_posicionA[num_procesadores + 1];
    int Q_posicionB[num_procesadores + 1];

    Q_posicionA[1] = 1;
    Q_posicionB[1] = 1;

    #pragma omp parallel for private(i, valor_posicion, valor_pivote, valor_solicitado) shared(arreglo_v, arregloA, arregloB, Q_posicionA, Q_posicionB, tamaño_arregloA, tamaño_arregloB, num_procesadores)
    for (i = 2; i <= num_procesadores; i++) {
        indice_v = 2*i - 3;
        if (arreglo_v[indice_v].origen == 'A') {
            valor_posicion = arreglo_v[indice_v].posicion_original;
            valor_pivote = arreglo_v[indice_v].valor;
            valor_solicitado = binary_search(arregloB, tamaño_arregloB, valor_pivote);
            if (valor_solicitado == -1) {
                valor_solicitado = tamaño_arregloB + 1;            
            }
            Q_posicionA[i] = valor_posicion * ((tamaño_arregloA + num_procesadores - 1) / num_procesadores);
            Q_posicionB[i] = valor_solicitado;
        } else {
            valor_posicion = arreglo_v[indice_v].posicion_original;
            valor_pivote = arreglo_v[indice_v].valor;
            valor_solicitado = binary_search(arregloA, tamaño_arregloA, valor_pivote);
            if (valor_solicitado == -1) {
                valor_solicitado = tamaño_arregloA + 1;            
            }
            Q_posicionA[i] = valor_solicitado;
            Q_posicionB[i] = valor_posicion * ((tamaño_arregloB + num_procesadores - 1) / num_procesadores);
        }
    }

    #pragma omp parallel for private(i, x, y, posicion_C, valor_limite, indice_v) shared(Q_posicionA, Q_posicionB, arregloA, arregloB, arregloC, arreglo_v, tamaño_arregloA, tamaño_arregloB, num_procesadores)
    for (i = 1; i <= num_procesadores; i++) {
        x = Q_posicionA[i];
        y = Q_posicionB[i];
        posicion_C = x + y - 1;

        if (i <= num_procesadores - 1) {
            indice_v = (2 * i) - 1;
            valor_limite = arreglo_v[indice_v].valor;
            sequential_merge(arregloA, arregloB, arregloC, tamaño_arregloA, tamaño_arregloB, x, y, posicion_C, 1, valor_limite);
        } else {
            sequential_merge(arregloA, arregloB, arregloC, tamaño_arregloA, tamaño_arregloB, x, y, posicion_C, 0, 0);
        }
    }
}

int main() {
    int arregloA[12] = {2, 5, 8, 11, 12, 15, 18, 21, 23, 34, 45, 47};
    int arregloB[12] = {1, 3, 4, 7, 10, 13, 16, 22, 29, 30, 41, 46};
    int tamaño_arregloA = sizeof(arregloA)/sizeof(arregloA[0]);
    int tamaño_arregloB = sizeof(arregloB)/sizeof(arregloB[0]);
    int arregloC[tamaño_arregloA + tamaño_arregloB];
    int valor_a_buscar = 49;
    int resultado;
    crew_merge(arregloA, arregloB, arregloC, tamaño_arregloA, tamaño_arregloB, 3);
    printf("El arreglo C es el siguiente: \n");

    // Imprimir arregloC en formato [v1, v2, ..., vn]
    printf("[");
    for (int i = 0; i < tamaño_arregloA + tamaño_arregloB; i++) {
        if (i > 0) {
            printf(", ");
        }
        printf("%d", arregloC[i]);
    }
    printf("]\n");
    /*resultado = binary_search(arregloA, tamaño_arregloA, valor_a_buscar);
    if (resultado == -1) {
        printf("El valor %d no se encuentra en el arreglo.\n", valor_a_buscar);
        return 0;
    } else {
        printf("El valor %d se encuentra en la posicion %d del arreglo.\n", valor_a_buscar, resultado);
    return 0;
    }*/ 
}