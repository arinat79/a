#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

double h = 0.01;

void print_matrix(double **a, int size)
{
    for (int i = 0; i< size; i++){
        for (int j = 0; j <  size; j++){
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}


double **copy_matrix(double **a, int size)
{
    double **a1 = calloc(size, sizeof(*a1));
    for (int i = 0; i < size; i++){
        a1[i] = calloc(size, sizeof(*a1[i]));
    }
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            a1[i][j] = a[i][j];
        }
    }
    return a1;
}

void scan_matrix(double **a, int size)
{
    for (int i = 0; i< size; i++){
        for (int j = 0; j <  size; j++){
            scanf("%lf", &a[i][j]);
        }
    }
}

void print_vect(double *a, int size)
{
    for (int i = 0; i < size; i++){
        printf("%f ", a[i]);
    }
    printf("\n");
}

void scan_vect(double *a, int size)
{
    for (int i = 0; i < size; i++){
        scanf("%lf", &a[i]);
    }
}



double input1 (double x, double y){
    //printf("%f ", x);
    return sin(x) - y;
}

double in1(double x, double y1, double y2)
{
    return y2 - cos(x);

}

double in2(double x, double y1, double y2)
{
    return y1 + sin(x);
}

double inp(double x)
{
    return 2;
}

double inq(double x)
{
    return -1/x;
}

double inf(double x)
{
    return -3;
}
void double_runge_kutt(FILE *out, double a, double b,double init,
                        double(*f)(double x, double y))
{
    double y = init;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        y = y + h * (f(x, y) + f(x + h, y + f(x, y) * h)) / 2;
    }
}


void square_runge_kutt(FILE *out, double a, double b, double init,
                        double(*f)(double x, double y))
{
    double y = init;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + k1 * h /2);
        double k3 = f(x + h/2, y + k2 * h/ 2);
        double k4 = f(x + h, y + h * k3);
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
}

void double_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,  double init1,  double init2,
                        double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = init1;
    double y2 = init2;
    for (double x = a; x <= b; x += h){
        fprintf(out1, "%f ", y1);
        fprintf(out2, "%f ", y2);
        double k1 = f1(x, y1, y2);
        double k21 = f2(x, y1, y2);
        double k2 = f1(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);
        double k22 = f2(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);

        y1 = y1 + h * k2;
        y2 = y2 + h * k22;
    }
}

void square_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,  double init1,  double init2,
                        double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = init1;
    double y2 = init2;
    for (double x = a; x <= b; x += h){
        fprintf(out1, "%f ", y1);
        fprintf(out2, "%f ", y2);
        double k1 = f1(x, y1, y2);
        double k21 = f2(x, y1, y2);

        double k2 = f1(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);
        double k22 = f2(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);

        double k3 = f1(x + h/2, y1 + k2 * h /2, y2 + k22 * h / 2);
        double k23 = f2(x + h/2, y1 + k2 * h /2, y2 + k22 * h / 2);

        double k4 = f1(x + h, y1 + k3 * h , y2 + k23 * h);
        double k24 = f2(x + h, y1 + k3 * h , y2 + k23 * h);

        y1 = y1 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y2 = y2 + h * (k21 + 2 * k22 + 2 * k23 + k24) / 6;
    }
}

double *tridiag_matrix(double **matrix, double *f, int n)
{
    double *solution = calloc(n, sizeof(*solution));
    double *a = calloc(n, sizeof(*a));
    double *b = calloc(n, sizeof(*b));
    double *c = calloc(n, sizeof(*c));
    double *alpha = calloc(n, sizeof(*alpha));
    double *beta = calloc(n, sizeof(*beta));

    for (int i = 1; i < n; ++i) {
        a[i] = matrix[i][i - 1];
    }

    for (int i = 0; i < n - 1; ++i) {
        b[i] = matrix[i][i + 1];
    }

    for (int i = 0; i < n; ++i) {
        c[i] = matrix[i][i];
    }

    if (fabs(c[0]) < DBL_EPSILON) {
      //  perror("Плохая матрица!\n");
        exit(1);
    }

    for (int i = 0; i < n; ++i) {
        if (fabs(c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1])) < DBL_EPSILON) {
          //  perror("Плохая матрица\n!");
            exit(1);
        }
        alpha[i] = -b[i] / (c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1]));
    }

    for (int i = 0; i < n; ++i) {
        beta[i] = (f[i] - a[i] * (i == 0 ? 0.0 : beta[i - 1])) /
                (c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1]));
    }

    for (int i = n - 1; i >= 0; --i) {
        solution[i] = beta[i] + alpha[i] * (i == n - 1 ? 0.0 : solution[i + 1]);
    }
    return solution;
}

//k[0]y(0) + k[1]y'(0) =k[2], k[3]y(1) + k[4]y'(1) =k[5]
void
boundary_problem(FILE *out, double a, double b,
                 double (*p)(double x), double (*q)(double x),
                 double (*f)(double x), double *k)
{
    int size = (b-a) / h + 1;
    double **m = calloc(size + 1, sizeof(*m));
    for (int i = 0; i < size + 1; i++){
        m[i] = calloc(size + 1, sizeof(*m[i]));
    }
    double *right = calloc(size + 1, sizeof(*right));
    m[0][0] = k[0] - k[1] / h;
    m[0][1] = k[1] / h;
    m[size][size-1] = - k[4] / h;
    m[size][size] = k[3] + k[4] / h;
    right[size] = k[5];
    right[0] = k[2];
    int j = 1;
    for (double i = a + h; i < b; i += h, j++) {
        m[j][j-1] = 1 - 0.5 * h * p(i);
        m[j][j] = q(i) * h * h - 2;
        m[j][j+1] = 1 + 0.5 * h * p(i);
        right[j] = f(i) * h * h;
    }

    for (int i = 1; i < size ; i++) {
        printf("matrix %f %f %f right %.10f\n", m[i][i-1], m[i][i], m[i][i+1], right[i]);
    }

    double *sol = tridiag_matrix(m, right, size + 1);



    for (int i = 0; i <= size; i++) {
        printf("%f\n", sol[i]);
    }

    return;
}



int main()
{
   FILE *out1 = fopen("x.txt", "w");
   for (double x = 0; x <= 10; x += h){
       fprintf(out1, "%f ", x);
   }

  /* FILE *output1 = fopen("out1.txt", "w");
   double_runge_kutt(output1, 0, 1, 10, input1);
   FILE *output2 = fopen("out2.txt", "w");
   square_runge_kutt(output2, 0, 1, 10, input1);

   FILE *output3 = fopen("out5.txt", "w");
   FILE *output4 = fopen("out4.txt", "w");
   square_runge_kutt_system(output3, output4, 0, 10, 0, 0, in1, in2);


   FILE *output6 = fopen("out6.txt", "w");
   FILE *output7 = fopen("out7.txt", "w");
   double_runge_kutt_system(output6, output7, 0, 10, 0, 0, in1, in2);
   // printf("%d Hello world!\n", file);
*/
   FILE *output_b = fopen("out_b.txt", "w");
   double *k = calloc(6, sizeof(*k));
   k[0] = 1;
   k[1] = 0;
   k[2] = 2;
   k[3] = 0.5;
   k[4] = -1;
   k[5] = 1;

   boundary_problem(output_b, 0.2, 0.5, inp, inq, inf, k);
    return 0;
}


