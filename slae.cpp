#include "slae.h"

SLAE::SLAE()
{
}

void SLAE::refreshRk(Matrix &a, double *f, double *x, double *r)
{
    multMatrixVector(a, x, r);		// r = Ax
    subVectorVector(f, r, a.n, r);	// r = f - Ax
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

double SLAE::norm(const double *a, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++)
        sum += a[i] * a[i];
    return sqrt(sum);
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

void SLAE::factorizeLLT(Matrix &a, Matrix &l)
{
    double *al = new double[a.ig[a.n]],
           *di = new double[a.n];

    for (size_t i = 0; i < a.n; i++) // Цикл по строкам
    {
        size_t iBegin = a.ig[i]; // Индекс первого элемента в строке i
        size_t iEnd = a.ig[i + 1];

        double sumDi = 0.0;

        for (size_t j = iBegin; j < iEnd; j++) // Цикл по элементам строки
        {
            size_t column = a.jg[j]; // Столбец текущего элемента

            size_t jBegin = a.ig[column]; // Индекс первого элемента в строке j
            size_t jEnd = a.ig[column + 1];

            double sum = 0.0;

            // Скалярное произведение Lik*Ljk
            size_t ii = iBegin;
            size_t jj = jBegin;
            while (ii < j && jj < jEnd && jj < j)
            {
                if (a.jg[ii] > a.jg[jj])
                    jj++;
                else if (a.jg[ii] < a.jg[jj])
                    ii++;
                else
                    sum += al[ii++] * al[jj++];
            }

            al[j] = (a.ggl[j] - sum) / di[column];

            sumDi += al[j] * al[j];
        }

        di[i] = sqrt(a.di[i] - sumDi);
    }

    l = Matrix(a);
    l.ggl = l.ggu = al;
    l.di = di;
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


void SLAE::forwardStrokeLLT(Matrix &l, double *b, double *y)
{
    for(size_t i = 0; i < l.n; i++) // Цикл по строкам
    {
        double sum = b[i];
        int iBegin = l.ig[i],
            iEnd = l.ig[i+1];
        for(int j = iBegin; j < iEnd; j++)
            sum -= l.ggl[j] * y[l.jg[j]];
        y[i] = sum / l.di[i];
    }
}


void SLAE::returnStrokeLLT(Matrix &l, double *b, double *x)
{
    for(size_t i = 0; i < l.n; i++)
        x[i] = b[i];
    for(int i = l.n-1; i >= 0; i--) // Цикл по столбцам
    {
        x[i] /= l.di[i];

        int iBegin = l.ig[i],
            iEnd = l.ig[i+1];

        for(int j = iBegin; j < iEnd; j++)
        {
            x[l.jg[j]] -= l.ggl[j] * x[i];
        }
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


void SLAE::iterationMSGLLT(Matrix &a, Matrix &l, double *x, double *r, double *z)
{
    static double *temp1 = new double[a.n],
                  *temp2 = new double[a.n];

    // M^-1 * r(k-1)
    forwardStrokeLLT(l, r, temp1);
    returnStrokeLLT(l, temp1, temp2);

    // (M^-1 * r(k-1), r(k-1))
    double Mrr = scalarProduct(temp2, r, a.n);

    // A * z(k-1)
    multMatrixVector(a, z, temp1);

    // (A * z(k-1), z(k-1))
    double Azz = scalarProduct(temp1, z, a.n);

    double alpha = Mrr / Azz;

    // x(k) = x(k-1) + alpha*z(k-1)
    multVectorScalar(z, alpha, a.n, temp2);
    addVectorVector(x, temp2, a.n, x);

    // r(k) = r(k-1) - alpha*Az(k-1)
    multVectorScalar(temp1, alpha, a.n, temp2);
    subVectorVector(r, temp2, a.n, r);

    // M^-1 * r(k)
    forwardStrokeLLT(l, r, temp1);
    returnStrokeLLT(l, temp1, temp2);

    // (M^-1 * r(k), r(k))
    double MrrNew = scalarProduct(temp2, r, a.n);

    double beta = MrrNew / Mrr;

    // z(k) = M^-1*r(k) + beta*z(k-1)
    multVectorScalar(z, beta, a.n, temp1);
    addVectorVector(temp2, temp1, a.n, z);
}


int SLAE::solveLOS_LU(Matrix &a, double *f, double *x, double eps, int maxIter)
{
    Matrix lu;						// Факторизованная матрица
    factorizeLU(a, lu);
    //factorizeLLT(a, lu);

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

    double residualSqr0 = residualSqr;

    cout << "#1:\t\t" << residualSqr << endl;

    int i;
    for(i = 2; i <= maxIter && residualSqr > eps; i++)
    {
        residualSqr -= iterationLU(a, lu, x, r, z, p);	// Квадрат нормы невязки (r(i), r(i)) = (r(i-1), r(i-1)) - alpha^2 * (p(i-1), p(i-1))

        residualSqr = scalarProduct(r, r, a.n);	// Квадрат нормы невязки (r(i), r(i))
        residualSqr /= residualSqr0;

        cout << "#" << i << ":\t\t" << residualSqr << endl;
    }

    // Освобождаем память
    delete lu.di;
    delete lu.ggl;
    //delete lu.ggu;
    delete r;
    delete z;
    delete p;

    return i - 1;	// Возвращаем число итераций
}


int SLAE::solveMSG_LLT(Matrix &a, double *f, double *x, double eps, int maxIter)
{
    Matrix l;						// Факторизованная матрица
    factorizeLLT(a, l);

    double *r = new double[a.n];	// *
    multMatrixVector(a, x, r);		// r = Ax
    subVectorVector(f, r, a.n, r);	// r = f - Ax

    double *z = new double[a.n];	// *
    forwardStrokeLLT(l, r, z);      // *
    returnStrokeLLT(l, z, z);		// z = M^-1 * r

    iterationMSGLLT(a, l, x, r, z);	// Первая итерация

    double residual = norm(r, a.n) / norm(f, a.n);

    cout << "#1:\t\t" << residual << endl;

    int i;
    for(i = 2; i <= maxIter && residual > eps; i++)
    {
        iterationMSGLLT(a, l, x, r, z);

        residual = norm(r, a.n) / norm(f, a.n);

        cout << "#" << i << ":\t\t" << residual << endl;
    }

    // Освобождаем память
    delete l.di;
    delete l.ggl;
    delete r;
    delete z;

    return i - 1;	// Возвращаем число итераций
}














