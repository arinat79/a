#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double h = 0.001;

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

int main()
{
   FILE *out1 = fopen("x.txt", "w");
   for (double x = 0; x <= 10; x += h){
       fprintf(out1, "%f ", x);
   }

   FILE *output1 = fopen("out1.txt", "w");
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
    return 0;
}


