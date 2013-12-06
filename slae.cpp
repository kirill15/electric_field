#include "slae.h"

SLAE::SLAE()
{
}

void SLAE::multMatrixVector(Matrix &a, double *x, double *res)
{
    // Обнуляем вектор результата
    for(size_t i = 0; i < a.n; i++)
        res[i] = 0.0;

    // Умножаем диагональ
    for(size_t i = 0; i < a.n; i++)
        res[i] += a.di[i] * x[i];

    // Умножаем нижний треугольник
    for(size_t i = 1; i < a.n; i++)
    {
        int iBegin = a.ig[i],
                iEnd = a.ig[i+1];
        for(int j = iBegin; j < iEnd; j++)
            res[i] += a.ggl[j] * x[a.jg[j]];
    }

    // Умножаем верхний треугольник
    for(size_t i = 1; i < a.n; i++)
    {
        int iBegin = a.ig[i],
                iEnd = a.ig[i+1];
        double xi = x[i];
        for(int j = iBegin; j < iEnd; j++)
            res[a.jg[j]] += a.ggu[j] * xi;
    }
}

void SLAE::addVectorVector(double *a, double *b, int size, double *res)
{
    for(int i = 0; i < size; i++)
        res[i] = a[i] + b[i];
}

void SLAE::subVectorVector(double *a, double *b, int size, double *res)
{
    for(int i = 0; i < size; i++)
        res[i] = a[i] - b[i];
}

void SLAE::multVectorScalar(const double *v, double c, int size, double *res)
{
    for(int i = 0; i < size; i++)
        res[i] = v[i] * c;
}

double SLAE::scalarProduct(const double *a, const double *b, int size)
{
    double res = 0.0;
    for(int i = 0; i < size; i++)
        res += a[i] * b[i];
    return res;
}


void SLAE::factorizeLU(Matrix &a, Matrix &lu)
{
    int k = a.ig[a.n]; // Количество элементов
    double *al = new double[k],
           *au = new double[k],
           *di = new double[a.n];

    for(size_t i = 0; i < a.n; i++) // Цикл по строкам (столбцам)
    {
        int iBegin = a.ig[i],
                iEnd = a.ig[i+1];
        double sumDi = 0.0;

        for(int j = iBegin; j < iEnd; j++) // Цикл по элементам строки
        {
            int column = a.jg[j],		// Столбец текущего элемента
                jBegin = a.ig[column],	// Индекс первого элемента в столбце (сверху элемента)
                jEnd = a.ig[column+1];	// jEnd-jBegin - число элементов в этом столбце

            double sumU = 0.0,
                   sumL = 0.0;

            // Считаем скалярное произведение:
            int iLeft = iBegin,		// Начальные значения
                jUp = jBegin;		// индексов
            while(iLeft < j && jUp < jEnd)
            {
                if(a.jg[iLeft] > a.jg[jUp])	// Если номер столбца в строке элемента больше, чем номер
                    jUp++;					// строки в столбце элемента, то увеличиваем индекс
                else
                    if(a.jg[iLeft] < a.jg[jUp])
                        iLeft++;
                    else
                    {
                        sumL += al[iLeft] * au[jUp];
                        sumU += au[iLeft++] * al[jUp++];
                    }
            }

            al[j] = a.ggl[j] - sumL;
            au[j] = (a.ggu[j] - sumU) / di[column];

            sumDi += al[j] * au[j];
        }
        di[i] = a.di[i] - sumDi;
    }

    lu = Matrix(a);
    lu.ggl = al;
    lu.di = di;
    lu.ggu = au;
}


void SLAE::forwardStroke(Matrix &lu, double *b, double *y)
{
    for(size_t i = 0; i < lu.n; i++) // Цикл по строкам
    {
        double sum = b[i];
        int iBegin = lu.ig[i],
            iEnd = lu.ig[i+1];
        for(int j = iBegin; j < iEnd; j++)
            sum -= lu.ggl[j] * y[lu.jg[j]];
        y[i] = sum / lu.di[i];
    }
}

void SLAE::returnStroke(Matrix &lu, double *b, double *x)
{
    for(size_t i = 0; i < lu.n; i++)
        x[i] = b[i];
    for(int i = lu.n-1; i >= 0; i--) // Цикл по столбцам
    {
        int iBegin = lu.ig[i],
            iEnd = lu.ig[i+1];

        for(int j = iBegin; j < iEnd; j++)
            x[lu.jg[j]] -= lu.ggu[j] * x[i];
    }
}

double SLAE::iterationLU(Matrix &a, Matrix &lu, double *x, double *r, double *z, double *p)
{
    double pSqr = scalarProduct(p, p, a.n);					// (p(k-1), p(k-1))
    double alpha = scalarProduct(p, r, a.n) / pSqr;			// alpha(k) = (p(k-1), r(k-1)) / (p(k-1), p(k-1))

    static double *temp1 = new double[a.n],					// *
                  *temp2 = new double[a.n];					// *

    multVectorScalar(z, alpha, a.n, temp1);					// temp1 = alpha(k)*z(k-1)
    addVectorVector(x, temp1, a.n, x);						// x(k) = x(k-1) + alpha*z(k-1)

    multVectorScalar(p, alpha, a.n, temp1);					// temp1 = alpha(k)*p(k-1)
    subVectorVector(r, temp1, a.n, r);						// r(k) = r(k-1) - alpha*p(k-1)

    returnStroke(lu, r, temp2);								// temp2 = U^-1 * r(k)
    multMatrixVector(a, temp2, temp1);						// temp1 = A * U^-1 * r(k)
    forwardStroke(lu, temp1, temp1);						// temp1 = L^-1 * A * U^-1 * r(k)
    double betta = -scalarProduct(p, temp1, a.n) / pSqr;	// betta(k) = -(p(k-1), L^-1*A*U^-1*r(k)) / (p(k-1), p(k-1))

    multVectorScalar(z, betta, a.n, z);						// z(k) = betta(k)*z(k-1)
    addVectorVector(temp2, z, a.n, z);						// z(k) = U^-1*r + betta(k)*z(k-1)

    multVectorScalar(p, betta, a.n, temp2);					// temp2 = betta(k)*p(k-1)
    addVectorVector(temp1, temp2, a.n, p);					// p(k) = L^-1 * A * U^-1 * r(k) + betta(k)*p(k-1)

    return alpha * alpha * pSqr;							// Возвращаем произведение alpha^2 * (p(k-1), p(k-1))
}

int SLAE::solveLOS_LU(Matrix &a, double *f, double *x, double eps, int maxIter)
{
    Matrix lu;						// Факторизованная матрица
    factorizeLU(a, lu);
    double *r = new double[a.n];	// *
    multMatrixVector(a, x, r);		// r = Ax
    subVectorVector(f, r, a.n, r);	// r = f - Ax
    forwardStroke(lu, r, r);		// r = L^-1 * (f - Ax)

    double *z = new double[a.n];	// *
    returnStroke(lu, r, z);			// z = U^-1 * r

    double *p = new double[a.n];	// *
    multMatrixVector(a, z, p);		// p = Az
    forwardStroke(lu, p, p);		// p = L^-1 * Az

    iterationLU(a, lu, x, r, z, p);	// Первая итерация
    double residualSqr = scalarProduct(r, r, a.n);	// Квадрат нормы невязки (r(i), r(i))
    cout << "#1:\t\t" << residualSqr << endl;

    int i;
    for(i = 2; i <= maxIter && residualSqr > eps; i++)
    {
        residualSqr -= iterationLU(a, lu, x, r, z, p);	// Квадрат нормы невязки (r(i), r(i)) = (r(i-1), r(i-1)) - alpha^2 * (p(i-1), p(i-1))
        cout << "#" << i << ":\t\t" << residualSqr << endl;
    }

    // Освобождаем память
    delete lu.di;
    delete lu.ggl;
    delete lu.ggu;
    delete r;
    delete z;
    delete p;

    return i - 1;	// Возвращаем число итераций
}

















