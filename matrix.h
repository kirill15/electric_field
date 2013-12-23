#ifndef MATRIX_H
#define MATRIX_H


// Матрица в разреженном строчном формате
/* Индексы с 0! */
struct Matrix
{
    unsigned int n;   // Размерность матрицы
    unsigned int *ig; // Указатели начала строк
    unsigned int *jg; // Номера столбоцов внедиагональных элементов матрицы
    double *ggl; 	  // Элементы нижнего треуголника
    double *ggu; 	  // Элементы верхнего треуголника
    double *di;		  // Диагональные элементы матрицы

    Matrix() : ig(nullptr), jg(nullptr), ggl(nullptr), ggu(nullptr), di(nullptr) {}

    ~Matrix()
    {/*
        delete[] ig;
        delete[] jg;
        delete[] ggl;
        //if (ggu)
          //  delete[] ggu;
        delete[] di;*/
    }
};

#endif // MATRIX_H
