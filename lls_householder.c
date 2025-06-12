
// Kompilieren dieses C-Quellcodes in Linux ueber die Konsoleneingabe 
// 'gcc -o lls_householder lls_householder.c -lm'


// Matrizen und Vektoren sind als 2- bzw. 1-dimensionale Arrays 
// realisiert



#include <stdio.h>
#include <math.h>
#include <float.h>



int ERROR=0;                  // zur Fehlerbehandlung
double EPSILON=DBL_EPSILON;   // 




// Ausgabe 1-dimensionales Array 
void ausgabeVek(int n, double array[n]) {
   int i;
   for(i=0; i<n; i++) {
        printf(" %12.10e ",array[i]);
   }
   printf("\n");
   printf("\n");
}


// berechne sgn(x)
int sgn(double x) {
   if(x<0) return -1;
   else return 1;
}


// loese das LGS Lu=d wobei L untere (pxn)-Dreiecksmatrix
// - die Diagonalelemente von L stehen in diag; die restl. Elemente
//   stehen in C unter der Hauptdiagonalen
// - das Ergebnis wird in die ersten p Stellen von erg geschrieben
//
void solve_lower(int p, int n, double C[p][n], double diag[p], double d[p], double erg[n]) {
     int i,k;
     double s=0;
     for(i=0; i<p; i++) {
       s=d[i];
       for(k=0; k<i; k++) s=s-(C[i][k]*erg[k]);
       erg[i]=s/diag[i];
     }
}


// loese das LGS Rv=r wobei R obere (n-p)x(n-p)-Dreiecksmatrix
// - R steht in den letzen n-p Spalten von A oberhalb der Hauptdiagonalen
// - r steht in den ersten n-p Stellen von tmp
// - das Ergebnis wird in die letzen n-p Stellen von erg geschrieben
//
void solve_upper(int m, int n, int p, double A[m][n], double erg[n], double tmp[m]) {
   int i,j;
   double s=0;
   for(i=n-p-1; i>=0; i--) {
      s=tmp[i];
      for(j=n-p-1; j>i; j--) s=s-(A[i][j+p]*erg[j+p]);
      erg[i+p]=s/A[i][i+p];
   }
}




// bestimme C^t=QR mittels Householdertransf. wobei C^t=C transponiert, C (pxn)-Matrix
// - für H_p...H_1*C^t=R ist Q=H_1...H_p (H_i Householdertransf.)
// - die Vektoren u, welche die H_i erzeugen, werden zeilenweise oberhalb der Diagonalen von C
//   gespeichert; die zugehoerigen Skalare stehen in beta
// - die Diagonalelemente von R stehen im Array diag
// - die restl. Eintraege von R stehen in C unter der Diagonalen
//
void transpose_house_decomp(int p, int n, double C[p][n], double beta[p], double diag[p]) {
   int i,j,k;
   double a, s;
   for(k=0; k<p; k++) {
      a=0;
      diag[k]=C[k][k];
      for(j=k; j<n; j++) {
        a=a+(C[k][j]*C[k][j]);
      }
      if(a<=EPSILON) {  // abbrechen, falls C numerisch nicht vollen Rang
         printf("Fehler: Martix C (numerisch) nicht vollen Zeilenrang\n");
         ERROR=1;
         break;
      } 
      a=sqrt(a);
      beta[k]=1/(a*(a+fabs(C[k][k])));
      C[k][k]=C[k][k]+sgn(C[k][k])*a;
      diag[k]=(-sgn(diag[k]))*a;
      for(i=k+1; i<p; i++) {
         s=0;
         for(j=k; j<n; j++) s=s+(C[k][j]*C[i][j]);
         s=s*beta[k];
         for(j=k; j<n; j++) C[i][j]=C[i][j]-(s*C[k][j]);
      }
   }
}



// multipliziere Produkt H_1...H_p von Householdertransf. mit Vektor erg
// (d.h. berechne H_1...H_p*erg)
// - die Vektoren u, welche die H_i erzeugen stehen in C oberhalb der Diagonalen;
//   die zugehoerigen Skalare stehen in beta
// - Ergebnis wird nach erg geschrieben
//
void mult_house_vec(int p, int n, double C[p][n], double beta[p], double erg[n]) {
   int i,j;
   double h;
   for(i=p-1; i>=0; i--) {
      h=0;
      for(j=i; j<n; j++) h=h+C[i][j]*erg[j];
      for(j=i; j<n; j++) erg[j]=erg[j]-(beta[i]*C[i][j]*h);
   }
}



// multipliziere Matrix A mit dem Produkt H_1...H_p von Householdertransf. 
// (d.h. berechne A*H_1...H_p)
// - die Vektoren u, welche die H_i erzeugen stehen in C oberhalb der Diagonalen;
//   die zugehoerigen Skalare stehen in beta
// - Ergebnis wird nach A geschrieben
//
void mult_house_mat(int m, int n, double A[m][n], int p, double C[p][n], double beta[p]) {
    int k,i,j;
    double h;
    for(k=0; k<p; k++) {
       for(i=0; i<m; i++) {
          h=0;
          for(j=k; j<n; j++) h=h+(A[i][j]*C[k][j]);
          for(j=k; j<n; j++) A[i][j]=A[i][j]-(beta[k]*h*C[k][j]);      
       }
    }
}



// berechne b-Eu wobei b,u Vektoren und E (mxp)-Matrix
// - E steht in den ersten p Spalten von A (d.h. A=[E F])
// - u steht in den ersten p Stellen von erg
// - Ergebnis wird nach tmp geschrieben
//
void comp(int m, int n, int p, double b[m], double A[m][n], double erg[n], double tmp[m]) {
     int i,j; 
     for(i=0; i<m; i++) {
        tmp[i]=0;
        for(j=0; j<p; j++) {
           tmp[i]=tmp[i]+(A[i][j]*erg[j]);
        }
     }
     for(i=0; i<m; i++) tmp[i]=b[i]-tmp[i];
}



// berechne QF=R, wobei Q=H_(n-p)....H_1 Produkt von Householdertransf., F mx(n-p)-Matrix; 
// gleichzeitig wird Q(b-Eu) berechnet wobei b,u Vektoren, E (mxp)-Matrix
// - E,F sind in A vermoege A=[E F] gespeichert
// - u steht in den ersten p Stellen von erg
// - R wird oberhalb der Diagonalen von F gespeichert
// - Ergebnis von Q(b-Eu) steht in tmp
//
void house_decomp(int m, int n, int p, double A[m][n], double erg[n], double tmp[m], double b[m]) {
     comp(m,n,p,b,A,erg,tmp);
     int i,j,k;
     double a,beta,s;
     double u[m][n-p];
     for(k=0; k<n-p; k++) {
        a=0;
        for(i=k; i<m; i++) {
            u[i][k]=A[i][p+k];
            a=a+(u[i][k]*u[i][k]);
        }
        if(a<=EPSILON) {   // abbrechen, falls Matrix 'C ueber A' numerisch nicht vollen Rang
          printf("Fehler: Matrix 'C ueber A' (numerisch) nicht vollen Spaltenrang\n");
          ERROR=1;
          break;
        }
        a=sqrt(a);
        beta=1/(a*(a+fabs(u[k][k])));
        u[k][k]=u[k][k]+(sgn(A[k][p+k])*a);
        
        A[k][p+k]=-sgn(A[k][p+k])*a;
        for(j=p+k+1; j<n; j++) {
            s=0;
            for(i=k; i<m; i++) s=s+(u[i][k]*A[i][j]);
            s=s*beta;
            for(i=k; i<m; i++) A[i][j]=A[i][j]-(s*u[i][k]);
        }
        s=0;
        for(i=k; i<m; i++) s=s+(u[i][k]*tmp[i]);
        s=s*beta;
        for(i=k; i<m; i++) tmp[i]=tmp[i]-(s*u[i][k]);
     }
}


// loese Ausgleichsproblem ||Ax-b||=min, Cx=d
// - Ergebnis x wird nach erg geschrieben, dann ausgegeben
//
void lin_Aus(int p, int n, int m, double A[m][n], double b[m], double C[p][n], double d[p]) {
     double diag[p], tmp[m], beta[p], erg[n];
     transpose_house_decomp(p,n,C,beta,diag); 
     if(ERROR==0) {                          // falls Matrix C numerisch vollen Rang hat...
         solve_lower(p,n,C,diag,d,erg);  
         mult_house_mat(m,n,A,p,C,beta);  
         house_decomp(m,n,p,A,erg,tmp,b);    
         if(ERROR==0) {                      // falls Matrix 'C ueber A' numerisch vollen Rang hat...
             solve_upper(m,n,p,A,erg,tmp);   
             mult_house_vec(p,n,C,beta,erg); 
             ausgabeVek(n,erg);   // Ausgabe Ergebnis
         }
     }   

}



// Hauptprogramm; wende Ausgleichsrechnung auf Testprobleme
// der Form ||Ax-b||=min, Cx=d an 
// - p,n,m sind Zeilen-/Spaltenzahlen der Matrizen A,C bzw.
// - der Vektoren b,d (vgl. Notation Aufabenzettel)
// - A (mxn)-Matrix
// - C (pxn)-Matrix
// - b m-stelliger Vektor
// - d p-stelliger Vektor
int main(void) {

  int p,n,m;

// Testproblem 1
          p=2;
          n=5;
          m=4;
          double b1[4]={1,1,1,1}, d1[2]={1,1};

          double C1[2][5]={{0, 0, 0, 0, 1},
                           {0, 0, 0, 1, 1}};

          double A1[4][5]={{2, 2, 3, 1, 2},
                           {0, 1, 1, 1, 2},
                           {0, 0, 1, 1, 3},
                           {1, 1, 1, 1, 2}}; 

          
          printf("Ergebnis Problem 1:  ");
          lin_Aus(p, n, m, A1, b1, C1, d1);

// Testproblem 2
          p=2;
          n=3;
          m=2;
          double b2[2]={1,-1}, d2[2]={0,3};

          double C2[2][3]={{-4, 0, 6},
                           { 0, 0, 12}};

          double A2[2][3]={{ 3,-9, 0},
                           {-1, 5,-15}}; 
   
          printf("Ergebnis Problem 2:  ");
          lin_Aus(p, n, m, A2, b2, C2, d2);
 
// Testproblem 3
          p=2;
          n=3;
          m=5;
          double b3[5]={ 0, 2, -1, 0, 0}, d3[2]={ 10, 4};

          double C3[2][3]={{13,  0, -1},
                           { 1, 22, -7}};

          double A3[5][3]={{-1, -4,  0},
                           { 6,  8,  0},
                           { 0, 10,  5},
                           { 3, -1,  3},
                           {-1, -1,  0}};   

          printf("Ergebnis Problem 3:  ");
          lin_Aus(p, n, m, A3, b3, C3, d3);

// Testproblem 4
          p=1;
          n=4;
          m=4;
          double b4[4]={ 19, 28, -23, -4}, d4[1]={ 12};

          double C4[1][4]={{3,-10, 4, -1}};

          double A4[4][4]={{ -2,  0, -2, 23},
                           { -1, 11, 10, 19},
                           {  1, 13, 14,  6},
                           {  4,-12, -8,  0}};
  
         printf("Ergebnis Problem 4:  ");
         lin_Aus(p, n, m, A4, b4, C4, d4);

// Testproblem 5
          p=1;   
          n=2;   
          m=3;   
          double b5[3]={ 1, 1, -2}, d5[1]={ 0};          

          double C5[1][2]={{ -13,  0}}; 

          double A5[3][2]={{ 24, -12},  
                           { -9,   2},
                           {  3,   8}}; 

          printf("Ergebnis Problem 5:  ");
          lin_Aus(p, n, m, A5, b5, C5, d5);


// Testproblem 6
          p=3;
          n=4;
          m=2;
          double b6[2]={ 4, -5}, d6[3]={0, 4, -3};

          double C6[3][4]={{ 1, -2, 0,  0},
                           { 0,  1, 12, 3},
                           { 4,  0, 5,  0}};

          double A6[2][4]={{ 0,  0, -3, 1},
                           { 0,  2,  3, 1}};
  
         printf("Ergebnis Problem 6:  ");
         lin_Aus(p, n, m, A6, b6, C6, d6);

     

   return 0;

}
