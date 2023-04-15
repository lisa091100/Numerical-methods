#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long long Gauss_addition = 0;
long long Gauss_multiplication = 0;
long long Gauss_division = 0;

//Печать матрицы
void
print(int size, double **mtrx_print)
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (mtrx_print[i][j] == 0) {
            printf("%*.10g ", 16, -mtrx_print[i][j]);
            } else {
            printf("%*.10g ", 16, mtrx_print[i][j]);
}
        }
        printf("\n");
    }
}

//Выделение памяти под новую матрицу
double **
memory_creation(int rows, int cols)
{
    double **mtrx_create = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; ++i)  {
        mtrx_create[i] = calloc(cols, sizeof(double));
    }
    if (mtrx_create == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    return mtrx_create;
}

//Метод Гаусса
double *
Gauss(int size, double **mtrx, double *value_f)
{
    double **matrix_Gauss = memory_creation(size, size);
    for(int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix_Gauss[i][j] = mtrx[i][j];
        }
    }
    double *f= calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
        f[i] = value_f[i];
    }
    double EPSILON = 0.00000001;
    double *tempt = calloc(size, sizeof(double));
    if (tempt == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    int p;
    //Прямой ход метода Гаусса
    for (int t = 0; t < size - 1; ++t) {
        //Если ведущий элемент равен 0
        if ((matrix_Gauss[t][t] <= EPSILON) && (matrix_Gauss[t][t] >= -EPSILON)) {
            p = t + 1;
            //Поиск первого ненулевого элемента
            while ((matrix_Gauss[p][t] <= EPSILON) && (matrix_Gauss[p][t] >= -EPSILON)) {
                ++p;
            }
            //Меняем местами строки t и p (элементы матрицы и значения)
            for (int i = 0; i < size; ++i) {
                tempt[i] = matrix_Gauss[p][i];
            }
            double value_f = f[p];
            for (int i = 0; i < size; ++i) {
                 matrix_Gauss[p][i] = matrix_Gauss[t][i];
            }
            f[p] = f[t];
            for (int i = 0; i < size; ++i) {
                matrix_Gauss[t][i] = tempt[i];
            }
            f[t] = value_f;
        }
        for (int k = t + 1; k < size; ++k) {
            //Вычисление коэффициента
            double coefficient = matrix_Gauss[k][t] / matrix_Gauss[t][t];
            ++Gauss_division;
            for (int i = t; i < size; ++i) {
                //Вычитание из k-той строчки t-той строки, умноженной на коэффициент
                matrix_Gauss[k][i] -= coefficient * matrix_Gauss[t][i];
                ++Gauss_multiplication;
                ++Gauss_addition;
            }
            f[k] -= coefficient * f[t];
            ++Gauss_multiplication;
            ++Gauss_addition;
        }
    }
    //Обратный ход метода Гаусса
    double *res = calloc(size, sizeof(double));
    if (res == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    if ((matrix_Gauss[size - 1][size - 1] <= EPSILON) && (matrix_Gauss[size - 1][size - 1] >= -EPSILON)) {
        res[size - 1] = 0;
    } else {
        res[size - 1] = f[size - 1] / matrix_Gauss[size - 1][size - 1];
        ++Gauss_division;
    }
    for (int i = size - 2; i >= 0; --i) {
        double value = f[i];
        for (int j = size - 1; j > i; --j) {
            value -= matrix_Gauss[i][j] * res[j];
            ++Gauss_multiplication;
            ++Gauss_addition;
        }
        if (value == 0) {
            res[i] = 0;
        } else {
            res[i] = value / matrix_Gauss[i][i];
            ++Gauss_division;
        }
    }
    free(tempt);
    free(f);
     for (int i = 0; i < size; ++i) {
        if (matrix_Gauss != NULL) {
            free(matrix_Gauss[i]);
        }
    }
    free(matrix_Gauss);
    return res;
}

//Метод Гаусса с выбором главного элемента
double *
Gauss_select(int size, double **mtrx, double *value_f)
{
    double EPSILON = 0.00000001;
    double **matrix_Gauss_select = memory_creation(size, size);
    for(int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix_Gauss_select[i][j] = mtrx[i][j];
        }
    }
    double *f= calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
        f[i] = value_f[i];
    }
    double max;
    int index_max;
    double *res = calloc(size, sizeof(double));
    int *res_num = calloc(size, sizeof(double));
    if (res == NULL || res_num == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    for (int i = 0; i < size; ++i) {
        res_num[i] = i;
    }
    //Прямой ход метода Гаусса
    for (int i = 0; i < size - 1; i++) {
        //Поиск максимального по модулю элемента в строке
        max = fabs(matrix_Gauss_select[i][i]);
        index_max = i;
        for (int k = i + 1; k < size; k++) {
            if (fabs(matrix_Gauss_select[i][k]) > max) {
                max = fabs(matrix_Gauss_select[i][k]);
                index_max = k;
            }
        }
        //Меняем столбцы с номерами i и index_max местами
        if (i != index_max) {
            for (int k = 0; k < size; k++) {
                double temp = matrix_Gauss_select[k][i];
                matrix_Gauss_select[k][i] = matrix_Gauss_select[k][index_max];
                matrix_Gauss_select[k][index_max] = temp;
            }
        }
        int exchange = res_num[i];
        res_num[i] = res_num[index_max];
        res_num[index_max] = exchange;
        for (int j = i + 1; j < size; j++) {
            //Вычисление коэффициента
            double coefficient = matrix_Gauss_select[j][i]/matrix_Gauss_select[i][i];
            //Вычитание из j-той строки i-той строки, умноженной на коэффициент
            for (int k = i; k < size; k++) {
                matrix_Gauss_select[j][k] -= coefficient * matrix_Gauss_select[i][k];
            }
            f[j] -= coefficient * f[i];
        }
    }
    //Обратный ход метода Гаусса
    if ((matrix_Gauss_select[size - 1][size - 1] <= EPSILON) && (matrix_Gauss_select[size - 1][size - 1] >= -EPSILON)) {
        res[size - 1] = 0;
    } else {
        res[size - 1] = f[size - 1] / matrix_Gauss_select[size - 1][size - 1];
    }
    for (int i = size - 2; i >= 0; --i) {
        double value = f[i];
        for (int j = size - 1; j > i; --j) {
            value -= matrix_Gauss_select[i][j] * res[j];
        }
        if (value == 0) {
            res[i] = 0;
        } else {
            res[i] = value / matrix_Gauss_select[i][i];
        }
    }
    for (int i = 0; i < size; i++) {
        if (res_num[i] != i) {
            double value = res[i];
            res[i] = res[res_num[i]];
            res[res_num[i]] = value;

            int temp = res_num[i];
            res_num[i] = res_num[temp];
            res_num[temp] = temp;
        }
    }
    free(res_num);
    free(f);
     for (int i = 0; i < size; ++i) {
        if (matrix_Gauss_select != NULL) {
            free(matrix_Gauss_select[i]);
        }
    }
    free(matrix_Gauss_select);
    return res;
}

//Вычисление определителя матрицы
double
det(int size, double **mtrx)
{
    double **mtrx_det = memory_creation(size, size);
    for(int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            mtrx_det[i][j] = mtrx[i][j];
        }
    }
    double res_det = 1;
    int sign = 1;
    double EPSILON = 0.00000001;
    double *tempt = calloc(size, sizeof(double));
    if (tempt == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    int p;
    //Прямой ход метода Гаусса
    for (int t = 0; t < size - 1; ++t) {
        //Если ведущий элемент равен 0
        if ((mtrx_det[t][t] <= EPSILON) && (mtrx_det[t][t] >= -EPSILON)) {
            p = t + 1;
            //Поиск первого ненулевого элемента
            while ((mtrx_det[p][t] <= EPSILON) && (mtrx_det[p][t] >= -EPSILON)) {
                ++p;
            }
            //Меняем местами строки t и p (элементы матрицы и значения)
            for (int i = 0; i < size; ++i) {
                tempt[i] = mtrx_det[p][i];
            }
            for (int i = 0; i < size; ++i) {
                 mtrx_det[p][i] = mtrx_det[t][i];
            }
            for (int i = 0; i < size; ++i) {
                mtrx_det[t][i] = tempt[i];
            }
            sign = -sign;
        }
        for (int k = t + 1; k < size; ++k) {
            //Вычисление коэффициента
            double coefficient = mtrx_det[k][t] / mtrx_det[t][t];
            for (int i = t; i < size; ++i) {
                //Вычитание из k-той строчки t-той строки, умноженной на коэффициент
                mtrx_det[k][i] -= coefficient * mtrx_det[t][i];
            }
        }
    }
    for (int i = 0; i < size; ++i) {
        res_det *= mtrx_det[i][i];
    }
    free(tempt);
    for (int i = 0; i < size; ++i) {
        if (mtrx_det != NULL) {
            free(mtrx_det[i]);
        }
    }
    free(mtrx_det);
    return res_det * sign;
}

//Обратная матрица
double **
inv_matrix(int size, double **mtrx)
{
    //Создание единичной матрицы
    double **invmtrx = memory_creation(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (j == i) {
                invmtrx[i][j] = 1;
            } else {
                invmtrx[i][j] = 0;
            }
        }
    }
    //Копирование матрицы
    double **matrix_inv = memory_creation(size, size);
    for(int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix_inv[i][j] = mtrx[i][j];
        }
    }
    double *tempt = calloc(size, sizeof(double));
    if (tempt == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    int p;
    double EPSILON = 0.00000001;
    //Прямой ход метода Гаусса
    for (int t = 0; t < size - 1; ++t) {
        //Если ведущий элемент равен 0
        if ((matrix_inv[t][t] <= EPSILON) && (matrix_inv[t][t] >= -EPSILON)) {
            p = t + 1;
            //Поиск первого ненулевого элемента
            while ((matrix_inv[p][t] <= EPSILON) && (matrix_inv[p][t] >= -EPSILON)) {
                ++p;
            }
            //Меняем местами строки t и p (элементы обычной и единичной матрицы)
            for (int i = 0; i < size; ++i) {
                tempt[i] = matrix_inv[p][i];
            }
            for (int i = 0; i < size; ++i) {
                 matrix_inv[p][i] = matrix_inv[t][i];
            }
            for (int i = 0; i < size; ++i) {
                matrix_inv[t][i] = tempt[i];
            }
            for (int i = 0; i < size; ++i) {
                tempt[i] = invmtrx[p][i];
            }
            for (int i = 0; i < size; ++i) {
                 invmtrx[p][i] = invmtrx[t][i];
            }
            for (int i = 0; i < size; ++i) {
                invmtrx[t][i] = tempt[i];
            }
        }
        for (int k = t + 1; k < size; ++k) {
            //Вычисление коэффициента
            double coefficient = matrix_inv[k][t] / matrix_inv[t][t];
            for (int i = 0; i < size; ++i) {
                //Вычитание из k-той строчки t-той строки, умноженной на коэффициент
                matrix_inv[k][i] -= coefficient * matrix_inv[t][i];
                invmtrx[k][i] -= coefficient * invmtrx[t][i];
            }
        }
    }
    //Oбратный ход вычисления элементов обратной матрицы
    for (int i = size - 1; i >= 0; --i) {
        for (int j = i; j > 0; --j) {
            double c = matrix_inv[j - 1][i] / matrix_inv[i][i];
            for (int k = size - 1; k >= 0; k--) {
                invmtrx[j - 1][k] -= invmtrx[i][k] * c;
            }
        }
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            invmtrx[i][j] /= matrix_inv[i][i];
        }
    }
    free(tempt);
    for (int i = 0; i < size; ++i) {
        if (matrix_inv != NULL) {
            free(matrix_inv[i]);
        }
    }
    free(matrix_inv);
    return invmtrx;
}

//Матрица по примеру из приложения 2
void **
create_example(double **mtrx, double *f, double x)
{
    int n = 100, M = 6;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (i != j) {
                mtrx[i - 1][j - 1] = pow(1.001 - 2 * M * pow(10, -3), i + j) + 0.1 * (j - i);
            } else {
                mtrx[i - 1][j - 1] = pow((1.001 - 2 * M * pow(10, -3) - 1), i + j);
            }
        }
        f[i - 1] = x * exp(x / i) * cos(x / i);
    }
    return 0;
}

//Норма матрицы
double
norm(int size, double **mtrx_norm)
{
    double result_norm = 0;
    for (int i = 0; i < size; ++i) {
        double sum = 0;
        for (int j = 0; j < size; ++j) {
            sum += fabs(mtrx_norm[i][j]);
        }
        if (sum > result_norm) {
            result_norm = sum;
        }
    }
    return result_norm;
}

//Число обусловленности
double
cond(int size, double **mtrx, double **inv_mtrx)
{
    return norm(size, mtrx) * norm(size, inv_mtrx);
}

//Норма разности двух векторов
double
norm_relax(double *current_x, double *previous_x, int size)
{
    double res = 0;
    for (int i = 0; i < size; ++i) {
        res += pow((current_x[i] - previous_x[i]), 2);
    }
    return sqrt(res);
}

// Метод верхней релаксации. Вычисление итерации
double
next_iteration(double **matrix, double *value_f, double *previous, double *current, int size, double omega)
{
    double first_amount, second_amount, amount;
    for (int i = 0; i < size; ++i) {
        first_amount = 0;
        second_amount = 0;
        amount = 0;
        for (int j = 0; j < i; ++j) {
            first_amount += matrix[i][j] * current[j];
        }
        for (int j = i; j < size; ++j) {
            second_amount += matrix[i][j] * previous[j];
        }
        amount += value_f[i] - first_amount - second_amount;
        amount *= (omega / matrix[i][i]);
        current[i] = previous[i] + amount;
    }
    return norm_relax(current, previous, size);
}

int
main(int argc, char **argv)
{
    int number_of_option;//Номер задания
    int size; //Размер матрицы
    fscanf(stdin, "%d\n", &number_of_option);
    fscanf(stdin, "%d", &size);
    double **matrix = memory_creation(size, size);//Исходная матрица
    double **inverse_matrix = memory_creation(size, size); //Обратная матрица
    double *vector_value = calloc(size, sizeof(double)); //Вектор значений
    double *res = calloc(size, sizeof(double)); //Вектор решений
    if (res == NULL || vector_value == NULL) {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(1);
    }
    if (number_of_option == 1) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                scanf("%lf ", &matrix[i][j]);
            }
            scanf("%lf", &vector_value[i]);
        }
    }
    if (number_of_option == 2) {
        double x;
        scanf("%lf", &x);
        create_example(matrix, vector_value, x);
    }
    res = Gauss(size, matrix, vector_value);
    printf("\nRoots of SLAE - Gauss' method:\n"); //Решение системы методом Гаусса
    for (int i = 0; i < size; ++i) {
        printf("x%d = %lf\n", i + 1, res[i]);
    }
    res = Gauss_select(size, matrix, vector_value);
    printf("\nRoots of SLAE - Gauss' method with selection:\n"); //Решение системы методом Гаусса с выбором главного элемента
    for (int i = 0; i < size; ++i) {
        printf("x%d = %lf\n", i + 1, res[i]);
    }
    printf("\nDeterminant of matrix: %.10g\n", det(size, matrix)); //Вычисление определителя матрицы
    inverse_matrix = inv_matrix(size, matrix); // Обратная матрица
    printf("\nInverse matrix:\n");
    print(size, inverse_matrix);
    printf("\nCondition number: %lf\n", cond(size, matrix, inverse_matrix));
    printf("\nOperations needed for Gauss' method:\n");
    printf("Addition: %lld\nMultiplication: %lld\nDivision: %lld\n", Gauss_addition, Gauss_multiplication, Gauss_division);

    int converge_flag; //Флаг сходимости
    double EPS;
    printf("Please, enter the accuracy\n");
    scanf("%lf", &EPS);

    long long max_iterations;
    printf("Please, enter the max iterations number\n");
    scanf("%lld", &max_iterations);

    double omega;
    printf("Please, enter the first omega value\n");
    scanf("%lf", &omega);

    double increment;
    printf("Please, enter the omega increment\n");
    scanf("%lf", &increment);

    double *prev_vector_x = calloc(size, sizeof(double));
    double *curr_vector_x = calloc(size, sizeof(double));
    double *vector_x = calloc(size, sizeof(double));
    double *tempt = calloc(size, sizeof(double));
    double best_w;
    int min_iterations = max_iterations;
    //Сходимость для 0 < w < 2
    for (double w = omega; w < 2; w += increment) {
        for (int i = 0; i < size; ++i) {
            // Принимаем за начальное приближение нулевой вектор
            prev_vector_x[i] = 0;
        }
        long long iterations = 0;
        converge_flag = 1;
        while (next_iteration(matrix, vector_value, prev_vector_x, curr_vector_x, size, w) > EPS) {
            iterations++;
            if (iterations > max_iterations) {
                converge_flag = 0;
                break;
            }
            tempt = prev_vector_x;
            prev_vector_x = curr_vector_x;
            curr_vector_x = tempt;
        }
        if (converge_flag == 1) {
            for (int i = 0; i < size; ++i) {
                if (!isfinite(curr_vector_x[i])) {
                    converge_flag = 0;
                    break;
                }
            }
        }
        if (converge_flag == 1) {
            printf("\nOmega = %lf\nIterations = %lld\n", w, iterations);
            if (iterations < min_iterations) {
                min_iterations = iterations;
                for (int i = 0; i < size; ++i) {
                    vector_x[i] = curr_vector_x[i];
                }
                best_w = w;
            }
        } else {
            printf("\nMethod diverges, Omega = %lf \n", w);
        }
    }
    if (min_iterations != -1) {
        printf("\nMethod converges = %lf\n", best_w);
        printf("\nRoots of SLAE — successive over-relaxation method:\n");
        for (int i = 0; i < size; ++i) {
            printf("x%d = %lf\n", i + 1, vector_x[i]);
        }
    } else {
        printf("%s\n", "Successive over-relaxation method diverges");
    }
    return 0;
}

