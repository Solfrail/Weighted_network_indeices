/*
weighted_average_path_length(wapl.c)
-fv (int): initial vertex
-deltail filename: output each distance between 2 vertices
-w : you decide which weight you use
重み付きネットワークの平均頂点間距離を算出するプログラム
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define bufmax 500000
#define inf DBL_MAX
#define NIL -1

struct cell{int vertex; double weight; double inv_weight; double using_weight; int next;};
//vertex: 頂点番号
//weight: 頂点間の枝重み（大きい程頂点間距離が大きい）
//inv_weight: 頂点間の枝重み（大きい程頂点間距離が大きい）
//using_weight: weight か inv_weight が入る
//next: struct cell は単方向リストとして扱うため、nextには、ある頂点kと繋がるcell.vertex以外の頂点情報が入ったポインタが入る
struct heap_cell{double key; int vertex;};

//data argorhythm
void bucket_insert(struct cell *L, int *B, int u, int v, double w1, double w2, int i);
int compare_ret_large(int a, int b, int c);
int List_Search(struct cell *L, int *B, int u, int v);
double List_Search_weight(struct cell *L, int *B, int u, int v);
double List_Search_inv_weight(struct cell *L, int *B, int u, int v);
double List_Search_apl(struct cell *L, int *B, int u, int v);
int non_connect_search(int *L, int num);

//heap function
int parent(int i);
int right(int i);
int left(int i);
void insert(struct heap_cell *H, int *adr, int hsize, int a, int v);
void decrease_key(struct heap_cell *H, int *adr, int i, int a);
int delete_min(struct heap_cell *H, int *adr, int hsize);
void upheap(struct heap_cell *H, int *adr, int hsize);
void downheap(struct heap_cell *H, int *adr, int hsize);
void upheap_sort(struct heap_cell *H, int *adr, int hsize);
void downheap_sort(struct heap_cell *H, int *adr, int hsize);

int n_size_data(int *comment);
int m_size_data(int *comment);
void vertex_size(int comment, int *N);
void read_data(struct cell *L, int *bucket, int *degree, int comment, int *N);

void option(int argc, char *argv[], int *fv, int *w, int *deltail, char *fname);
int main(int argc, char *argv[]){
    static int N;//num of vertex
    static int u, v;
    static int i, j, k;
    double sum_length=0;

    int *degree;
    int *degree_dist;
    struct cell *L;//branch list
    int *bucket;//First [vertex]->next vertex list
    
    int v0; /* 始点の変数*/
    double *d; //distance
    double s, t; //tmp
    int *p, *adr; //parent, vertices involved in heap
    int hsize=0; // the number of vertices which is in heap
    double dist_tmp;

    struct heap_cell *Heap;
    
    long int connect_num;
    int *non_connect;
    int first_lap;
    
    int edge, column; //stdin line: edge, stdin column, column
    int comment = 0;

    //option
    int fv=0;
    char *fname="test.dat"; //filename
    int detail=0;
    FILE *fp;
    int weight = 1;
    
    //initialize
    option(argc, argv, &fv, &weight, &detail, fname);
    edge = n_size_data(&comment);
    column = m_size_data(&comment);
    if(column != 4){
      fprintf(stderr, "Wrong File Format\n");
      return EXIT_FAILURE;
    }

    if((L = (struct cell *)malloc(2*edge*sizeof(struct cell)))==NULL){
      fprintf(stderr, "L: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }

    for(i = 0; i < 2 * edge; i++){
      L[i].vertex = -1;
      L[i].weight = -1;
      L[i].inv_weight = -1;
      L[i].using_weight = -1;
      L[i].next = -1;
    }


    //error processes
    fprintf(stderr, "#comment: %d\n", comment);
    vertex_size(comment, &N);
    if((bucket = (int *)calloc(N, sizeof(int))) == NULL){
      fprintf(stderr, "bucket: Cannot allocate memory\n");
      return EXIT_FAILURE;
    }
    if((degree = (int *)calloc(N, sizeof(int))) == NULL){
      fprintf(stderr, "degree: Cannot allocate memory\n");
      return EXIT_FAILURE;
    }
    if((degree_dist = (int *)malloc(sizeof(int) * N))==NULL){
      fprintf(stderr, "degree_dist: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }
    if((non_connect = (int *)malloc(sizeof(int) * N))==NULL){
      fprintf(stderr, "degree_dist: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }

    if((d = (double *)malloc(sizeof(double) * N))==NULL){
      fprintf(stderr, "degree_dist: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }
    if((adr = (int *)malloc(sizeof(int) * N))==NULL){
      fprintf(stderr, "degree_dist: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }
    if((p = (int *)malloc(sizeof(int) * N))==NULL){
      fprintf(stderr, "degree_dist: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }
    if((Heap = (struct heap_cell *)malloc(sizeof(struct heap_cell) * N))==NULL){
      fprintf(stderr, "L: Cannot allocate memory.\n");
      return EXIT_FAILURE;
    }

    
    
    for(i=0;i<N;i++){
      bucket[i] = NIL;
      non_connect[i] = NIL;
    }

    fprintf(stderr, "initialize_complete\n");
    read_data(L, bucket, degree, comment, &N);

    if(weight == 1){
        for(i = 0; i < 2 * edge; ++i){
	          L[i].using_weight = L[i].weight;
	          //fprintf(stderr, "%d %lf %lf\n",i, L[i].using_weight, L[i].weight);
        }
    }else if(weight == 2){
        for(i = 0; i < 2 * edge; ++i){
	          L[i].using_weight = L[i].inv_weight;
        }
    }

    if(detail!=0){fp = fopen(argv[detail], "a");}//prepare
    
    fprintf(stderr, "#Vertex: %d\n", N);
    fprintf(stderr, "#Branch: %d\n", edge);
    if(detail != 0){
      fprintf(fp, "#Vertex: %d\n", N);
      fprintf(fp, "#Branch: %d\n", edge);
    }
    
    connect_num = N;
    
    //calculate weighted ave_path_length
    k = 0;//count non_connect_num
    first_lap = 1;
    

    for(j=0;j<N;j++){
        if(j%100==0){fprintf(stderr, "%d\n", j);}
      
        v0 = j;
        if(fv!=0 && first_lap == 1){v0 = fv;}
        else if(fv==0 || (fv!=0 && first_lap == 0) ){}
        else{
	          fprintf(stderr, "ERROR. continue impossible.\n");
	      return 0;
        }
      
        if(j>=N){break;}//break
      
        //printf("Before dijkstra completed.\n");
        //dijkstra start
        v = 0;
        while (v < N) {//initialize
	          //if (v != v0) {
            adr[v] = -1;
            d[v] = inf;
            p[v] = -1;
            Heap[v].key = inf;
            Heap[v].vertex = -1;
            //}
            v++;
        }
      
        //v0設定
        d[v0] = 0;
        p[v0] = -1;
        insert(Heap, adr, hsize, d[v0], v0);
        hsize++;
        
        while(hsize != 0){
          v = delete_min(Heap, adr, hsize);
          hsize--;
          for(u=0;u<N;u++){
            if((dist_tmp = List_Search_apl(L, bucket, u, v)) != inf){
              //fprintf(stderr, "%d, %d, %lf\n", u, v, dist_tmp);
              //fprintf(stderr, "a\n");
              if(d[u]==inf){
                d[u] = d[v] + dist_tmp;
                p[u] = v;
                insert(Heap, adr, hsize, d[u], u);
                hsize++;
              }else if(d[u] > (d[v] + dist_tmp)){
                d[u] = d[v] + dist_tmp;
                p[u] = v;
                decrease_key(Heap, adr, d[u], u);
              }
            }
          }
          u=0;
          t = inf;
          while(u<N){
            if(adr[u] <= 0){
              s = d[u];
              if(s < t){
                t = s;
                v = u;
              }
            }
            ++u;
          }
      }

      if(fv == 0 || (fv != 0 && first_lap == 0) ){
        for (i = j + 1; i < N; i++) {
          if(d[i] != inf){
            if(detail!=0){	      
              fprintf(fp, "%d %d %.9lf\n", j, i, d[i]);
            } 
            sum_length += d[i];
          }
          else{
            if(first_lap == 1 && non_connect_search(non_connect, i)==0){
              connect_num--;
              non_connect[k] = i;
              k++;
            }
          }
	      }
      }
      else if(fv != 0 && first_lap == 1){
        for(i = 0; i < N; i++){
          if(d[i] == inf){
            connect_num--;
            non_connect[k] = i;
            k++;
            j = -1;
          }
        }
      }
      first_lap = 0;
    }

    if(detail!=0){fclose(fp);}
    
    fprintf(stdout, "%lf\n", (double)sum_length/(double)(connect_num*(connect_num-1)/2));
    fprintf(stderr, "#sum_length = %lf\n", (double)sum_length/(double)(connect_num*(connect_num-1)/2));
    
    free(degree_dist);
    free(degree);
    free(L);
    free(bucket);
    free(non_connect);
    free(d);
    free(adr);
    free(p);
    free(Heap);
    return 0;
}


//option
void option(int argc, char *argv[], int *fv, int *w, int *detail, char *fname){
  int i;
  int iter=argc;
  for(i=0;i<iter;i++){
    if(strncmp("-fv", argv[i], 20)==0 && (i+1<iter)){*fv=atoi(argv[++i]);}
    else if(strncmp("-w", argv[i], 20)==0 && (i+1<iter)){*w=atoi(argv[++i]);}
    else if(strncmp("-detail", argv[i], 20)==0 && (i+1<iter)){
      //detail=argv[++i];
      //fprintf(stderr, "%s\n", fname);
      *detail=i+1;
    }
  }
}

//Reading Files
int n_size_data(int *comment){
  int n=0;
  char buf[bufmax];
  while(fgets(buf, sizeof(buf), stdin)!=NULL){
    //comment output
    if(strncmp(buf, "%", 1)==0 || strncmp(buf, "#", 1)==0){
      //fprintf(stderr, "%s", buf);
      //fprintf(stdout, "%s", buf);
      *comment += 1;
    }
    else{n++;}
  }
  fseek(stdin, 0L, SEEK_SET);
  return n;
}
int m_size_data(int *comment){
  int i, d, num=0, dim=0, n=0;
  char buf[bufmax];
  char *tok, *space=" ";
  while(fgets(buf, sizeof(buf), stdin)){
    for(i=0;i<*comment;i++){
      fgets(buf, sizeof(buf), stdin);
    }
    d = 0;
    tok = strtok(buf, space);
    while(tok!=NULL){
      d++;
      tok = strtok(NULL, space);
    }
    dim += d;
    n++;
  }
  num = n;
  dim /= num;
  fseek(stdin, 0L, SEEK_SET);
  return dim;
}
void vertex_size(int comment, int *N){
  //Road Data
  int com = comment;
  int tmp = 0;
  char buf[bufmax];
  char *tok, *space=" ";
  int u, v;
  
  while(fgets(buf, sizeof(buf), stdin)!=NULL){
    if(com == 0){
      tok = strtok(buf, space);
      u = atoi(tok);
      tok = strtok(NULL, space);
      v = atoi(tok);
      tmp = compare_ret_large(u, v, tmp);
    }else{
      com -= 1;
    }
  }
  
  *N = tmp+1; // calculate number of vertex

  fseek(stdin, 0L, SEEK_SET);
}

void read_data(struct cell *L, int *bucket, int *degree, int comment, int *N){
  int com = comment;
  char buf[bufmax];
  char *tok, *space=" ";
  int u, v;
  double weight, inv_weight;
  int i = 0;
  
  while(fgets(buf, sizeof(buf), stdin)!=NULL){
    if(com == 0){
      tok = strtok(buf, space);
      u = atoi(tok);
      tok = strtok(NULL, space);
      v = atoi(tok);
      tok = strtok(NULL, space);
      weight = atof(tok);
      tok = strtok(NULL, space);
      inv_weight = atof(tok);
      //fprintf(stderr, "%lf %lf\n", weight, inv_weight);
      bucket_insert(L, bucket, u, v, weight, inv_weight, i);
      i++;
      bucket_insert(L, bucket, v, u, weight, inv_weight, i);
      i++;
      degree[u]++;
      degree[v]++;
    }else{
      com -= 1;
    }
  }
  fseek(stdin, 0L, SEEK_SET);
}
  

//functions
void bucket_insert(struct cell *L, int *B, int u, int v, double w1, double w2, int i){
  L[i].vertex = v;
  L[i].weight = w1;
  L[i].inv_weight = w2;
  L[i].next = B[u];
  B[u] = i;
}

int compare_ret_large(int a, int b, int c){
  if(a <= c && b <= c){return c;}
  else if(a <= b && c <= b){return b;}
  else if(b <= a && c <= a){return a;}
  else{
    printf("ERROR\n");
    return 0;
  }
}

int List_Search(struct cell *L, int *B, int u, int v){
  int a;
  a = B[u];
  while(a != NIL){
    if(L[a].vertex == v) return 1;
    a = L[a].next;
  }
  return 0;
}

//you need the opreation in *1
//if cluster -> 0.0,  else if length -> inf
double List_Search_weight(struct cell *L, int *B, int u, int v){
  int a;
  a = B[u];
  while(a != NIL){
    if(L[a].vertex == v){
      //fprintf(stderr, "%lf\n", L[a].weight);
      return L[a].weight;
    }
    a = L[a].next;
  }
  return inf;  //*1
}

double List_Search_inv_weight(struct cell *L, int *B, int u, int v){
  int a;
  a = B[u];
  while(a != NIL){
    if(L[a].vertex == v){
      //fprintf(stderr, "%lf\n", L[a].weight);
      return L[a].inv_weight;
    }
    a = L[a].next;
  }
  return inf;  //*1
}

double List_Search_apl(struct cell *L, int *B, int u, int v){
  int a;
  a = B[u];
  while(a != NIL){
    if(L[a].vertex == v){
      //fprintf(stderr, "%d ", a);
      return L[a].using_weight;
    }
    a = L[a].next;
  }
  return inf;  //*1
}


int non_connect_search(int *L, int num){
  int i;
  i = 0;
  while(L[i] != NIL){
    if(num == L[i]){return 1;}
    i++;
  }
  return 0; 
}

//heap
int parent(int i){
  int ans;
  if(i%2==0){
    ans = i/2 - 1;
  }
  else{
    ans = (i-1)/2;
  }
  return ans;
}
int right(int i){
  i = 2*i+2;
  return i;
}
int left(int i){
  i = 2*i+1;
  return i;
}
void insert(struct heap_cell*H, int *adr, int hsize, int a, int v){
  hsize++;
  H[hsize-1].key = a;
  H[hsize-1].vertex = v;
  adr[v] = hsize-1;
  upheap_sort(H, adr, hsize-1);
}
void decrease_key(struct heap_cell *H, int *adr, int i, int a){
  H[a].key = i;
  upheap_sort(H, adr, adr[a]);
}
int delete_min(struct heap_cell *H, int *adr, int hsize){
  int v;
  v = H[0].vertex;
  adr[v] = -1;
  if(hsize>1){
    H[0].key=H[hsize-1].key;
    H[0].vertex=H[hsize-1].vertex;
    adr[H[0].vertex]=0;
    downheap_sort(H, adr, hsize-2);
  }
  return v;
}
void upheap_sort(struct heap_cell *H, int *adr, int hsize){
  int u, tmp;
  u = hsize;
  while((u>0)&&(H[parent(u)].key>H[u].key)){
    tmp = H[u].key;
    H[u].key = H[parent(u)].key;
    H[parent(u)].key = tmp;

    tmp = H[u].vertex;
    H[u].vertex = H[parent(u)].vertex;
    H[parent(u)].vertex = tmp;

    tmp = adr[H[u].vertex];
    adr[H[u].vertex] = adr[H[parent(u)].vertex];
    adr[H[parent(u)].vertex] = tmp;

    u = parent(u);
  }
}
void downheap_sort(struct heap_cell *H, int *adr, int hsize){
  int u, r, l, tmp;
  struct heap_cell tmps;
  u = 0;
  while(((H[u].key>H[left(u)].key)||(H[u].key>H[right(u)].key))&&((l!=u)||(r!=u))){
    if(left(u)<=hsize){
      l = left(u);
    }
    else{
      l = u;
    }
    if(right(u)<=hsize){
      r = right(u);
    }
    else{
      r = u;
    }
    if(H[l].key<H[r].key){
      tmps = H[l];
      H[l] = H[u];
      H[u] = tmps;

      tmp = adr[l];
      adr[l] = adr[u];
      adr[u] = tmp;
      u = l;
    }
    else{
      tmps = H[r];
      H[r] = H[u];
      H[u] = tmps;

      tmp = adr[r];
      adr[r] = adr[u];
      adr[u] = tmp;
      u = r;
    }
  }

}
