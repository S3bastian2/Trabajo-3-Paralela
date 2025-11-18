#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define Nprocs 16
#define r 15
#define s 15

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

    while ((a < eofa && *a < vlim) || (b < eofb && *b < vlim) || (c < eofc && vlim < 0)) {
        if (a == eofa && b == eofb) break;
        
        if (b == eofb || *a < *b){              // si no quedan elementos en b, o si a[i] < b[i]
            *c = *a;
            a++;
        }else if(a == eofa || *b < *a){         // si no quedan elementos en a, o si b[i] < a[i]
            *c = *b;
            b++;
        }

        c++; 
    }
}

int main(int argc, char **argv) {
    int N;

    N = (argc > 1) ? atoi(argv[1]) : Nprocs;


    int a[r], b[s], aprime[N], bprime[N], c[r+s];
    triplet v[2*N - 1];
    pair q[N];

    for(int i = 0; i<r+s; i++) c[i] = 0;

    for(int i = 0; i<r; i++) a[i] = 2*i;

    for(int i = 0; i<s; i++) b[i] = 2*i + 1;

    omp_set_num_threads(N);

    #pragma omp parallel for
    for (int i = 0; i<N; i++) {
        int aidx = i*ceil(r/N), bidx = i*ceil(s/N);
        aprime[i] = a[aidx >= r? r : aidx];
        bprime[i] = bprime[bidx >= s? s : bidx];
    }

    #pragma omp parallel for
    for (int i = 0; i<N; i++){
        //int t_id = omp_get_thread_num();
        int j = rigth_insort(bprime, N, aprime[i]);


        if(j == N){                   //si es mas grande que todo bprime' se inserta al final.
            v[i+N].elem = aprime[i];
            v[i+N].id = i;
            v[i+N-1].arr_id = 0;       //arr_id 0 significa que viene de a.
        }else{
            v[i+j].elem = aprime[i];
            v[i+j].id =i;
            v[i+j].arr_id = 0;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i<N; i++){
        //int t_id = omp_get_thread_num();

        int j = rigth_insort(aprime, N, bprime[i]);

        if(j == N){
            v[i+N].elem = bprime[i];
            v[i+N].id = i;
            v[i+N].arr_id = 1;
        }else{
            v[i+j].elem = bprime[i];
            v[i+j].id =i;
            v[i+j].arr_id = 1;
        }
    }

    q[0].x = 0;
    q[0].y = 0;
    
    #pragma omp parallel for
    for (int i = 1; i<=N; i++){
        if (!v[2*i -1].arr_id) {    //solo da true si viene de a
            int j = rigth_insort(b, s, v[2*i -1].elem);

            q[i].x = (v[2*i-1].id) * ceil(r/N);
            q[i].y = j;             //si a'[k] es mayor que todo b, j será automaticamente s.

        } else if(v[2*i -1].arr_id) {
            int j = rigth_insort(a, r, v[2*i -1].elem);
            
            q[i].x = (v[2*i-1].id) * ceil(s/N);
            q[i].y = j;             //si a'[k] es mayor que todo b, j será automaticamente s.


        }
    }

    seq_merge(a, 0, a+r, b, 0, b+s, c, c+r+s, v[0].elem);

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        if (i < N - 1) {
            int vlim = v[2*i].elem;
            seq_merge(a, q[i].x, a+r, b, q[i].y, b+s, c, c+r+s, vlim);
        } else {
            seq_merge(a, q[i].x, a+r, b, q[i].y, b+s, c, c+r+s, -1);
        }
    }

    printf("num processors = %d\n",N);
    for(int _ = 0; _<r+s; _++) printf("%d ", c[_]);
    return 0;


    return 0;    
}
