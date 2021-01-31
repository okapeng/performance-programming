#include <stdio.h>

int main(int argc, char **argv){
  int i,j,k;
  double d;
  float f;
  void *p;
  struct fred {
     int a;
     short s;
     double d;
     short s2;
     double d2;
  }var;
  double a[30][30];

  printf("i %d %016lx\n",sizeof(i),&i);

}
