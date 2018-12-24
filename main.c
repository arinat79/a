#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

double h = 0.001;

void print_m(double **a, int size)
{
    for (int i = 0; i< size; i++){
        for (int j = 0; j <  size; j++){
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}


double **copy_m(double **a, int size)
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

void scan_m(double **a, int size)
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
    return sin(x) - y;
}

double ans1(double x){
    return -0.5*cos(x) + 0.5*sin(x) +10.5 * exp(-x);
}

double input2 (double x, double y){
    return 3 - y - x;
}

double ans2(double x){
    return 4 - x - 4 * exp(-x);
}

double input3 (double x, double y){
    return -y - x*x;
}

double ans3(double x){
    return -x*x*x + 2*x-2+12*exp(-x);
}


double sys11(double x, double y1, double y2)
{
    return -2*x*y1*y1 + y2*y2 -x -1;

}

double sys21(double x, double y1, double y2)
{
    return 1/(y2 * y2) - y1 - x/y1;
}


double sys12(double x, double y1, double y2)
{
    return x * y1 + y2;

}

double sys22(double x, double y1, double y2)
{
    return y1 - y2;
}


double sys13(double x, double y1, double y2)
{
    return (y1-y2) / x;

}

double sys23(double x, double y1, double y2)
{
    return (y1 + y2) / x;
}

double prec1(double x)
{
    return x * (cos(log(x)) - sin(log(x)));
}

double prec2(double x)
{
    return x * (cos(log(x)) + sin(log(x)));
}

void double_runge_kutt(FILE *out, double a, double b,double init,
                        double(*f)(double x, double y))
{
    double y = init;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        y = y + h/2 * (f(x, y) + f(x + h, y + f(x, y) * h));
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

void double_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,
                              double init1,  double init2,
                              double(*f1)(double x, double y1, double y2),
                              double(*f2)(double x, double y1, double y2))
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

void square_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,
                              double init1,  double init2,
                              double(*f1)(double x, double y1, double y2),
                              double(*f2)(double x, double y1, double y2))
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

double *tridiag(double **m, double *f, int n)
{
    double EPS = 1e-9;
    double *res = calloc(n, sizeof(*res));
    double *a = calloc(n, sizeof(*a));
    double *b = calloc(n, sizeof(*b));
    double *c = calloc(n, sizeof(*c));
    double *alpha = calloc(n, sizeof(*alpha));
    double *beta = calloc(n, sizeof(*beta));

    for (int i = 0; i < n; ++i) {
        a[i] = (i > 0 ? m[i][i - 1] : 0);
        b[i] = (i < n - 1 ? m[i][i + 1] : 0);
        c[i] = m[i][i];
    }
    if (fabs(c[0]) < EPS) {
        exit(1);
    }
    alpha[0] = -b[0] / c[0];
    beta[0] = f[0] / c[0];
    for (int i = 1; i < n; ++i) {
        if (fabs(c[i] + a[i] * alpha[i - 1]) < EPS) {
            exit(1);
        }
        alpha[i] = -b[i] / (c[i] + a[i] * alpha[i - 1]);
    }

    for (int i = 0; i < n; ++i) {
        beta[i] = (f[i] - a[i] * beta[i - 1]) /(c[i] + a[i] * alpha[i - 1]);
    }
    res[n - 1] = beta[n - 1] + alpha[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        res[i] = beta[i] + alpha[i] * res[i + 1];
    }
    return res;
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
    m[size][size-1] = - k[4] /  h;
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

    for (int i = 0; i <= size ; i++) {
        printf("m %f %f %f right %.10f\n", m[i][i-1], m[i][i], m[i][i+1], right[i]);
    }

    double *sol = tridiag(m, right, size + 1);

    for (int i = 0; i <= size; i++) {
        fprintf(out, "%f\n", sol[i]);
    }
}

double c1(double x)
{
    return 2;
}

double c2(double x)
{
    return -1/x;
}

double c3(double x)
{
    return 3;
}

double b1(double x)
{
    return 0;
}

double b2(double x)
{
    return 1;
}

double b3(double x)
{
    return 4*sin(x);
}

double sol1(double x)
{
    return sin(x) +(3-2*x)*cos(x);
}

double d1(double x)
{
    return 2;
}

double d2(double x)
{
    return 1;
}

double d3(double x)
{
    return 1;
}

double sol2(double x)
{
    return 1 - x*exp(1-x);
}


void test(FILE *o1, FILE *o2, FILE *o3, double init,
           double (*f)(double x), double (*ans)(double x))
{
    double_runge_kutt(o1, 0, 1, init, f);
    square_runge_kutt(o2, 0, 1, init, f);
    for (double x = 0; x <= 1; x += h){
       fprintf(o3, "%f ", ans(x));
   }
}

int main()
{
   FILE *outx = fopen("x.txt", "w");
   double x;
   for (x = 0.2; x <= 0.5; x += h){
       fprintf(outx, "%f ", x);
   }
   fprintf(outx, "%f ", x);

   FILE *output1 = fopen("out1.txt", "w");
   FILE *output2 = fopen("out2.txt", "w");
   FILE *output3 = fopen("out3.txt", "w");
   test(output1, output2, output3, 10, input1, ans1);

   FILE *output4 = fopen("out4.txt", "w");
   FILE *output5 = fopen("out5.txt", "w");
   FILE *output6 = fopen("out6.txt", "w");
   test(output4, output5, output6, 0, input2, ans2);

   FILE *output7 = fopen("out7.txt", "w");
   FILE *output8 = fopen("out8.txt", "w");
   FILE *output9 = fopen("out9.txt", "w");
   test(output7, output8, output9, 10, input3, ans3);

   FILE *s1 = fopen("1.txt", "w");
   FILE *s2 = fopen("2.txt", "w");
   FILE *s11 = fopen("11.txt", "w");
   FILE *s12 = fopen("12.txt", "w");
   double_runge_kutt_system(s1, s2, 0, 1, 1, 1, sys11, sys21);
   square_runge_kutt_system(s11, s12, 0, 1, 1, 1, sys11, sys21);

   FILE *s21 = fopen("21.txt", "w");
   FILE *s22 = fopen("22.txt", "w");
   FILE *s211 = fopen("211.txt", "w");
   FILE *s212 = fopen("212.txt", "w");
   double_runge_kutt_system(s21, s22, 0, 1, 0, 1, sys12, sys22);
   square_runge_kutt_system(s211, s212, 0, 1, 0, 1, sys12, sys22);

   FILE *s31 = fopen("31.txt", "w");
   FILE *s32 = fopen("32.txt", "w");
   FILE *s311 = fopen("311.txt", "w");
   FILE *s312 = fopen("312.txt", "w");
   FILE *p1 = fopen("p1.txt", "w");
   FILE *p2 = fopen("p2.txt", "w");
   double_runge_kutt_system(s31, s32, 1, 10, 1, 1, sys13, sys23);
   square_runge_kutt_system(s311, s312, 1, 10, 1, 1, sys13, sys23);
   for  (x = 1; x <= 10; x += h){
       fprintf(p1, "%f ", prec1(x));
   }
   for  (x = 1; x <= 10; x += h){
       fprintf(p2, "%f ", prec2(x));
   }

   FILE *output_c = fopen("out_c.txt", "w");
   double *k = calloc(6, sizeof(*k));
   k[0] = 1;
   k[1] = 0;
   k[2] = 2;
   k[3] = 0.5;
   k[4] = -1;
   k[5] = 1;

   boundary_problem(output_c, 0.2, 0.5, c1, c2, c3, k);

   FILE *output_b = fopen("out_b.txt", "w");
   k[0] = 1;
   k[1] = 1;
   k[2] = 2;
   k[3] = 1;
   k[4] = 1;
   k[5] = 0;
   FILE *output_b1 = fopen("out_b1.txt", "w");
   boundary_problem(output_b, 0, 1, b1, b2, b3, k);
   for  (x = 0; x <= 1; x += h){
       fprintf(output_b1, "%f ", sol1(x));
   }
   fprintf(output_b1, "%f ", sol1(x));

   FILE *output_d = fopen("out_d.txt", "w");
   k[0] = 1;
   k[1] = 0;
   k[2] = 1;
   k[3] = 1;
   k[4] = 1;
   k[5] = 0;
   FILE *output_d1 = fopen("out_d1.txt", "w");
   boundary_problem(output_d, 0, 1, d1, d2, d3, k);
   for (x = 0; x <= 1; x += h){
       fprintf(output_d1, "%f ", sol2(x));
   }
   //fprintf(output_d1, "%f ", sol2(x)); */


    return 0;
}
