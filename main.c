#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double EPS = 1e-9;

int change_str(double **a, int x, int y)
{
    double *tmp = a[x];
    a[x] = a[y];
    a[y] = tmp;
    return -1;
}

int triangle(double **a, int size)
{
    int j = 0;
    int sign = 1;
    for (int i = 0; i < size; i++){
        // Поиск первого ненулевого элемента в i столбце начиная с i+1 столбца
        for (j = i; j < size; j++){
            if (a[j][i]){
                break;
            }
        }
        if (i != j){
            sign *= change_str(a, i, j);
        }
        for (int k = i + 1; k < size; k++){
            double k1 = a[i][i];
            double k2 = a[k][i];
            if (abs(k2) < EPS){
                continue;
            }
            for (int l = 0; l < size ; l++){
                a[k][l] -= a[i][l] * k2 / k1;
            }
        }
    }

    return sign;
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

double det(double **a1, int size)
{
    double **a = copy_matrix(a1, size);
    int sign = triangle(a, size);
    double ans = 1;
    for (int i = 0; i < size; i++){
        ans *= a[i][i];
    }
    return ans * sign;
}



void print_matrix(double **a, int size)
{
    for (int i = 0; i< size; i++){
        for (int j = 0; j <  size; j++){
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}

void scan_matrix(double **a, int size)
{
    for (int i = 0; i< size; i++){
        for (int j = 0; j <  size; j++){
            scanf("%lf", &a[i][j]);
        }
    }
}


double **inverse(double **a1, int size)
{
    double **a = copy_matrix(a1, size);
    double **res = calloc(size, sizeof(*res));
    for (int i = 0; i < size; i++){
        res[i] = calloc(res, sizeof(*res[i]));
    }
    for (int i = 0; i < size; i++){
        res[i][i]  = 1;
    }

    int j = 0;
    int sign = 1;
    for (int i = 0; i < size; i++){
        // Поиск первого ненулевого элемента в i столбце начиная с i+1 столбца
        for (j = i; j < size; j++){
            if (a[j][i]){
                break;
            }
        }
        if (i != j){
            sign *= change_str(a, i, j);
            change_str(res, i, j);
        }

        for (int k = i + 1; k < size; k++){
            double k1 = a[i][i];
            double k2 = a[k][i];
            if (fabs(k2) < EPS){
                continue;
            }
            for (int l = 0; l < size ; l++){
                a[k][l] -= a[i][l] * k2 / k1;
                res[k][l] -= res[i][l] * k2 / k1;
            }
        }
    }
    for (int i = size - 1; i > 0; i--){
        for (int j = i - 1; j >= 0; j--){
            double k1 = a[j][i] / a[i][i];
            for (int k = size - 1; k >=0; k--){
                a[j][k] -= k1 * a[i][k];
                res[j][k] -= k1 * res[i][k];
            }
        }
    }
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            res[i][j] /= a[i][i];
        }
    }
    return res;
}

double *gauss(double **a1, double *f1, int size)
{
    double **a = copy_matrix(a1, size);
    double *f = calloc(size, sizeof(*f));
    for (int i = 0; i < size; i++){
        f[i] =f1[i];
    }
    //прямой ход
    int j = 0;
    int sign = 1;
    for (int i = 0; i < size; i++){
        // Поиск первого ненулевого элемента в i столбце начиная с i+1 столбца
        for (j = i; j < size; j++){
            if (a[j][i]){
                break;
            }
        }
        if (i != j){
            sign *= change_str(a, i, j);
            double tmp = f[i];
            f[i] = f[j];
            f[j] = tmp;
        }
        for (int k = i + 1; k < size; k++){
            double k1 = a[i][i];
            double k2 = a[k][i];
            if (fabs(k2) < EPS){
                continue;
            }
            f[k] -= f[i] * k2 / k1;
            for (int l = 0; l < size ; l++){
                a[k][l] -= a[i][l] * k2 / k1;
            }
        }

    }
    //обратный ход
    for (int i = size - 1; i > 0; i--){
        for (int j = i - 1; j >= 0; j--){
            double k1 = a[j][i] / a[i][i];
            f[j] -= k1 * f[i];
            for (int k = size - 1; k >=0; k--){
                a[j][k] -= k1 * a[i][k];
            }
        }
    }
    for (int j = 0; j < size; j++){
        f[j] /= a[j][j];
    }
    return f;
}

double *gauss_main(double **a1, double *f1, int size)
{
    double **a = copy_matrix(a1, size);
    double *f = calloc(size, sizeof(*f));
    //хранит перестановки столбцов
    int *vec = calloc(size, sizeof(*vec));

    for (int i = 0; i < size; i++){
        f[i] =f1[i];
    }
    //прямой ход
    int j = 0;
    int sign = 1;
    for (int i = 0; i < size; i++){
        vec[i] = i;
    }
    for (int i = 0; i < size; i++){
        // Поиск первого ненулевого элемента в i столбце начиная с i+1 столбца
        double max = 0;
        int idx = 0;
        for (j = i; j < size; j++){
            if (fabs(a[i][j]) > max){
                max = fabs(a[i][j]);
                idx = j;
            }
        }
        vec[i] = idx;
        vec[idx] = i;
        for (int l = 0; l < size; l++){
            double tmp = a[l][i];
            a[l][i] = a[l][idx];
            a[l][idx] = tmp;
        }

        for (int k = i + 1; k < size; k++){
            double k1 = a[i][i];
            double k2 = a[k][i];
            if (abs(k2) < EPS){
                continue;
            }
            f[k] -= f[i] * k2 / k1;
            for (int l = 0; l < size ; l++){
                a[k][l] -= a[i][l] * k2 / k1;
            }
        }

    }
    //обратный ход
    for (int i = size - 1; i > 0; i--){
        for (int j = i - 1; j >= 0; j--){
            double k1 = a[j][i] / a[i][i];
            f[j] -= k1 * f[i];
            for (int k = size - 1; k >=0; k--){
                a[j][k] -= k1 * a[i][k];
            }
        }
    }
    for (int j = 0; j < size; j++){
        f[j] /= a[j][j];
    }
    double *ans = calloc(size, sizeof(*ans));
    for (int i = 0; i < size; i++){
        ans[i] = f[vec[i]];
    }
    return ans;
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

double cond_num(double **a, int n)
{
    double res1 = 0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(fabs(a[i][j]) > res1){
                res1 = fabs(a[i][j]);
            }
        }
    }
    double res2 = 0;
    double **b = inverse(a, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(fabs(a[i][j]) > res2){
                res2 = fabs(a[i][j]);
            }
        }
    }
    return res1 * res2;
}

double *relax(double **a, double *f, double w, int size, int max_iter, double solution_eps)
{
    //критерий остановки - максимальное числло итераций
    int iter_count  = 0;
    double *x = calloc(size, sizeof(*x));
    double *tmp = calloc(size, sizeof(*tmp));
    for (int k = 0; k < max_iter; k++){
        for (int i = 0; i < size; i++){
            double sub1 = 0, sub2 = 0;
            for (int j = 0; j < i; j++){
                sub1 += a[i][j] * tmp[j];
            }
            for (int j = i; j < size; j++){
                sub2 += a[i][j] * x[j];
            }
            tmp[i] = x[i] + w / a[i][i] * (f[i] - sub1 - sub2);

            //Проверяем 2 критерий остановки алгоритма - расстояние между векторами решениями
		    //на 2 последовательных шагах становится меньше заданного значения - solution_epsilon
        }

        double dist = 0;
        for (int i = 0; i < size; i++){
            dist += (x[i] -tmp[i]) * (x[i] -tmp[i]);
        }
        if (sqrt(dist) < solution_eps){
            iter_count = k;
            return tmp;
        }
        for (int i = 0; i < size; i++){
            x[i] = tmp[i];
        }
    }
    iter_count = max_iter;
    return x;
}

int main()
{
    int n;
    scanf("%d", &n);
    double **a = calloc(n, sizeof(*a));
    for (int i = 0; i < n; i++){
        a[i] = calloc(n, sizeof(*a[i]));
    }

    scan_matrix(a, n);

    double *f = calloc(n, sizeof(*f));
    scan_vect(f, n);
    double *ans = relax(a, f, 1, n, 500, 1e-5);
    print_vect(ans, n);
   // printf("%f\n", cond_num(a, n));
    return 0;
}
