#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define Nprocs 13
#define r 12
#define s 12

typedef struct triplet {
    int elem, id, arr_id;
} triplet;

typedef struct pair {
    int x, y;
} pair;


int rigth_insort(int *a, int n, int k){
    int low = 0, high = n - 1, mid;
    while (low <= high) {
        mid = low + (high - low) / 2;
        if (mid >= n){
            return n;
        }

        if (a[mid] == k) {
            return mid;
        } else if (a[mid] < k) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return low;
}

void seq_merge(int *a, int x, int *eofa, int *b, int y, int *eofb, int *c, int *eofc, int vlim){
    a += x;
    b += y;
    c += (x+y);

    while ((a < eofa && *a < vlim) || (b < eofb && *b < vlim)) {
        if (a == eofa || b == eofb) break;
        
        if (b == eofb || *a < *b){              // si no quedan elementos en b, o si a[i] < b[i]
            *c = *a;
            a++;
        }else if(a == eofa || *b < *a){         // si no quedan elementos en a, o si b[i] < a[i]
            *c = *b;
            b++;
        }

        c++; 
    }

    while(vlim < 0 && (b < eofb || a < eofa)){
        if(c>=eofc) break;

        if(a < eofa){
            *c = *a;
            a++;
        }else if(b < eofb){
            *c = *b;
            b++;
        }
        c++;
    }
}

int main(int argc, char **argv) {
    int N;

    N = (argc > 1) ? atoi(argv[1]) : Nprocs;

    if (N < 2) N = 2; // ensure positive sizes for aprime/bprime

    int chunk_r = (r + N - 1) / N;
    int chunk_s = (s + N - 1) / N;

    int a[r], b[s], aprime[N-1], bprime[N-1], c[r+s];
    triplet v[2*N - 2];
    pair q[N];

    for(int i = 0; i<r+s; i++) c[i] = 0;

    for(int i = 0; i<r; i++) a[i] = 2*i;
    for(int i = 0; i<s; i++) b[i] = 2*i + 1;

    omp_set_num_threads(N);

    #pragma omp parallel for
    for (int i = 0; i<N-1; i++) {
        int aidx = (i+1) * chunk_r, bidx = (i+1) * chunk_s;
        aprime[i] = a[(aidx >= r? r-1 : aidx-1)];
        bprime[i] = b[(bidx >= s? s-1 : bidx-1)];
        printf("hilo %d selecciona pivotes a'[%d] = %d, b'[%d] = %d\n", omp_get_thread_num(), i, aprime[i], i, bprime[i]);
    }

    for(int _ = 0; _<N-1;_++) printf("%d, ", aprime[_]);
    for(int _ = 0; _<N-1;_++) printf("%d, ", bprime[_]);
    printf("\n");

    // sequential build of v to avoid races
    for (int i = 0; i<N-1; i++){
        int j = rigth_insort(bprime, N-1, aprime[i]);

        if(j == N-1){
            v[i+N-1].elem = aprime[i];
            v[i+N-1].id =i;
            v[i+N-1].arr_id = 0;
        }else{
            v[i+j].elem = aprime[i];
            v[i+j].id =i;
            v[i+j].arr_id = 0;
        }

        printf("hilo %d crea tripleta v[%d] = {%d,%d,%c}\n", omp_get_thread_num(), i+j, aprime[i], i, 'A'); 
    }

    printf("\n");

    // sequential build of v from bprime
    for (int i = 0; i<N-1; i++){
        int j = rigth_insort(aprime, N-1, bprime[i]);

        if(j == N-1){
            v[i+N-1].elem = bprime[i];
            v[i+N-1].id =i;
            v[i+N-1].arr_id = 1;
        }else{
            v[i+j].elem = bprime[i];
            v[i+j].id =i;
            v[i+j].arr_id = 1;
        }
        printf("hilo %d crea tripleta v[%d] = {%d,%d,%c}\n", omp_get_thread_num(), i+j, bprime[i], i, 'B');
    }

    printf("\n");

    q[0].x = 0;
    q[0].y = 0;

    // sequential build of q to avoid races and ensure consistent indices
    for (int i = 1; i<N; i++){
        if (!v[2*i-1].arr_id) {    // from a
            int j = rigth_insort(b, s, v[2*i -1].elem);

            q[i].x = (v[2*i-1].id + 1) * chunk_r -1;
            q[i].y = j; 

        } else {
            int j = rigth_insort(a, r, v[2*i -1].elem);
            q[i].x = j;
            q[i].y = (v[2*i-1].id + 1) * chunk_s -1;
        }
        printf("Procesador %d crea dupla Q[%d] = (%d,%d)\n", omp_get_thread_num(), i, q[i].x+1, q[i].y+1);
    }

    printf("\n");
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        if (i < N-1) {
            int vlim = v[2*i+1].elem;
            seq_merge(a, q[i].x, a+r, b, q[i].y, b+s, c, c+r+s, vlim);
        } else {
            seq_merge(a, q[i].x, a+r, b, q[i].y, b+s, c, c+r+s, -1);
        }

        printf("procesador %d mezcla a partir de a[%d], b[%d], hasta hallar el limite %d\n", omp_get_thread_num(), q[i].x, q[i].y, (i<N-1?v[2*i+1].elem:-1));
    }

    printf("num processors = %d\n arreglo resultante:\nc = [",N);
    for(int _ = 0; _<r+s; _++) printf((_<r+s-1?"%d, ":"%d"), c[_]);
    printf("]");

    return 0;    
}
