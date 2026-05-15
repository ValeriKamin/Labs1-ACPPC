using System;
using System.Globalization;


class Program
{
    const double EPS = 1e-10;

    static void Main()
    {
        Console.OutputEncoding = System.Text.Encoding.UTF8;
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

        double[,] A = null;

        while (true)
        {
            Console.WriteLine("\nПРАКТИЧНА РОБОТА 1А");
            Console.WriteLine("");
            Console.WriteLine("1 - Ввести матрицю A");
            Console.WriteLine("2 - Показати матрицю A");
            Console.WriteLine("3 - Знайти обернену матрицю");
            Console.WriteLine("4 - Обчислити ранг матриці");
            Console.WriteLine("5 - Розв'язати СЛАР через обернену матрицю");
            Console.WriteLine("6 - Розв'язати СЛАР методом Жордана");
            Console.WriteLine("7 - Розв'язати СЛАР методом Гауса");
            Console.WriteLine("0 - Вихід");
            Console.Write("Ваш вибір: ");

            int choice = int.Parse(Console.ReadLine());

            if (choice == 0) break;

            switch (choice)
            {
                case 1:
                    A = InputMatrix();
                    break;

                case 2:
                    if (A == null) Console.WriteLine("Спочатку введіть матрицю.");
                    else PrintMatrix(A);
                    break;

                case 3:
                    if (CheckMatrix(A) && CheckSquare(A)) RunInverse(A);
                    break;

                case 4:
                    if (CheckMatrix(A)) RunRank(A);
                    break;

                case 5:
                    if (CheckMatrix(A) && CheckSquare(A)) RunByInverse(A);
                    break;

                case 6:
                    if (CheckMatrix(A) && CheckSquare(A)) RunByJordan(A);
                    break;

                case 7:
                    if (CheckMatrix(A) && CheckSquare(A)) RunByGauss(A);
                    break;

                default:
                    Console.WriteLine("Невірний вибір.");
                    break;
            }
        }
    }

    // ---------- Запуск пунктів меню ------

    static void RunInverse(double[,] A)
    {
        Console.WriteLine("\nПошук оберненої матриці методом Жорданових виключень.");

        double[,] inv = Inverse(A, true);

        if (inv == null)
        {
            Console.WriteLine("Обернена матриця не існує.");
            return;
        }

        Console.WriteLine("Обернена матриця:");
        PrintMatrix(inv);
    }

    static void RunRank(double[,] A)
    {
        Console.WriteLine("\nОбчислення рангу матриці.");
        int rank = Rank(A, true);
        Console.WriteLine($"Ранг матриці = {rank}");
    }

    static void RunByInverse(double[,] A)
    {
        double[] B = ReadVector(A.GetLength(0), "B");

        Console.WriteLine("\n1-й спосіб: X = A^(-1) * B");

        double[] X = SolveByInverse(A, B, true);

        if (X == null)
        {
            Console.WriteLine("Система не має єдиного розв'язку.");
            return;
        }

        PrintVector(X);
    }

    static void RunByJordan(double[,] A)
    {
        double[] B = ReadVector(A.GetLength(0), "B");

        Console.WriteLine("\n2-й спосіб: метод Жордана");

        double[] X = SolveByJordan(A, B, true);

        if (X == null)
        {
            Console.WriteLine("Система не має єдиного розв'язку.");
            return;
        }

        PrintVector(X);
    }

    static void RunByGauss(double[,] A)
    {
        double[] B = ReadVector(A.GetLength(0), "B");

        Console.WriteLine("\n3-й спосіб: метод Гауса");

        double[] X = SolveByGauss(A, B, true);

        if (X == null)
        {
            Console.WriteLine("Система не має єдиного розв'язку.");
            return;
        }

        PrintVector(X);
    }

    // ---------- Основні методи ----------

    static double[,] Inverse(double[,] A, bool protocol)
    {
        int n = A.GetLength(0);
        double[,] table = WithIdentity(A);

        if (protocol)
        {
            Console.WriteLine("Початкова матриця [A | E]:");
            PrintMatrix(table);
        }

        bool ok = JordanToIdentity(table, n, protocol);

        if (!ok) return null;

        return RightHalf(table);
    }

    static int Rank(double[,] A, bool protocol)
    {
        double[,] table = Copy(A);
        int rows = table.GetLength(0);
        int cols = table.GetLength(1);

        int rank = 0;
        int row = 0;

        if (protocol)
        {
            Console.WriteLine("Початкова матриця:");
            PrintMatrix(table);
        }

        for (int col = 0; col < cols && row < rows; col++)
        {
            int pivotRow = FindPivotRow(table, row, col);

            if (pivotRow == -1)
                continue;

            SwapRows(table, row, pivotRow);

            if (protocol)
            {
                Console.WriteLine($"\nКрок {rank + 1}");
                Console.WriteLine($"Ведучий елемент: A[{row + 1},{col + 1}] = {table[row, col]:F4}");
            }

            JordanPivot(table, row, col, protocol);

            row++;
            rank++;
        }

        return rank;
    }

    static double[] SolveByInverse(double[,] A, double[] B, bool protocol)
    {
        double[,] inv = Inverse(A, protocol);

        if (inv == null) return null;

        Console.WriteLine("Обчислюємо X = A^(-1) * B:");

        return Multiply(inv, B, true);
    }

    static double[] SolveByJordan(double[,] A, double[] B, bool protocol)
    {
        int n = A.GetLength(0);
        double[,] table = Augmented(A, B);

        if (protocol)
        {
            Console.WriteLine("Початкова розширена матриця [A | B]:");
            PrintMatrix(table);
        }

        bool ok = JordanToIdentity(table, n, protocol);

        if (!ok) return null;

        return LastColumn(table);
    }

    static double[] SolveByGauss(double[,] A, double[] B, bool protocol)
    {
        int n = A.GetLength(0);
        double[,] table = Augmented(A, B);

        if (protocol)
        {
            Console.WriteLine("Початкова розширена матриця [A | B]:");
            PrintMatrix(table);
            Console.WriteLine("Прямий хід:");
        }

        for (int k = 0; k < n - 1; k++)
        {
            int pivotRow = FindPivotRow(table, k, k);

            if (pivotRow == -1) return null;

            SwapRows(table, k, pivotRow);

            if (protocol)
                Console.WriteLine($"\nКрок {k + 1}, ведучий елемент = {table[k, k]:F4}");

            for (int i = k + 1; i < n; i++)
            {
                double factor = table[i, k] / table[k, k];

                if (protocol)
                    Console.WriteLine($"R{i + 1} = R{i + 1} - ({factor:F4}) * R{k + 1}");

                for (int j = k; j <= n; j++)
                    table[i, j] -= factor * table[k, j];
            }

            if (protocol)
                PrintMatrix(table);
        }

        if (Math.Abs(table[n - 1, n - 1]) < EPS)
            return null;

        double[] X = new double[n];

        Console.WriteLine("Зворотний хід:");

        for (int i = n - 1; i >= 0; i--)
        {
            double sum = table[i, n];

            for (int j = i + 1; j < n; j++)
                sum -= table[i, j] * X[j];

            X[i] = sum / table[i, i];

            if (protocol)
                Console.WriteLine($"x{i + 1} = {X[i]:F4}");
        }

        return X;
    }

    // ---------- Жорданові виключення ----------

    static bool JordanToIdentity(double[,] table, int n, bool protocol)
    {
        for (int k = 0; k < n; k++)
        {
            int pivotRow = FindPivotRow(table, k, k);

            if (pivotRow == -1)
                return false;

            SwapRows(table, k, pivotRow);

            if (protocol)
            {
                Console.WriteLine($"\nКрок {k + 1}");
                Console.WriteLine($"Ведучий елемент: A[{k + 1},{k + 1}] = {table[k, k]:F4}");
            }

            JordanPivot(table, k, k, protocol);
        }

        return true;
    }

    static void JordanPivot(double[,] table, int pivotRow, int pivotCol, bool protocol)
    {
        int rows = table.GetLength(0);
        int cols = table.GetLength(1);

        double pivot = table[pivotRow, pivotCol];

        for (int j = 0; j < cols; j++)
            table[pivotRow, j] /= pivot;

        for (int i = 0; i < rows; i++)
        {
            if (i == pivotRow) continue;

            double factor = table[i, pivotCol];

            for (int j = 0; j < cols; j++)
                table[i, j] -= factor * table[pivotRow, j];
        }

        if (protocol)
        {
            Console.WriteLine("Матриця після Жорданового виключення:");
            PrintMatrix(table);
        }
    }

    // ---------- Допоміжні методи ----------

    static double[,] InputMatrix()
    {
        Console.Write("Введіть кількість рядків y: ");
        int rows = int.Parse(Console.ReadLine());

        Console.Write("Введіть кількість стовпців x: ");
        int cols = int.Parse(Console.ReadLine());

        return ReadMatrix(rows, cols, "A");
    }

    static double[,] ReadMatrix(int rows, int cols, string name)
    {
        double[,] M = new double[rows, cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                Console.Write($"{name}[{i + 1},{j + 1}] = ");
                M[i, j] = double.Parse(Console.ReadLine());
            }

        return M;
    }

    static double[] ReadVector(int size, string name)
    {
        double[] V = new double[size];

        for (int i = 0; i < size; i++)
        {
            Console.Write($"{name}[{i + 1}] = ");
            V[i] = double.Parse(Console.ReadLine());
        }

        return V;
    }

    static double[,] Augmented(double[,] A, double[] B)
    {
        int rows = A.GetLength(0);
        int cols = A.GetLength(1);

        double[,] T = new double[rows, cols + 1];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                T[i, j] = A[i, j];

            T[i, cols] = B[i];
        }

        return T;
    }

    static double[,] WithIdentity(double[,] A)
    {
        int n = A.GetLength(0);
        double[,] T = new double[n, 2 * n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                T[i, j] = A[i, j];

            for (int j = 0; j < n; j++)
                T[i, j + n] = i == j ? 1 : 0;
        }

        return T;
    }

    static double[,] RightHalf(double[,] T)
    {
        int n = T.GetLength(0);
        double[,] R = new double[n, n];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                R[i, j] = T[i, j + n];

        return R;
    }

    static double[] LastColumn(double[,] T)
    {
        int rows = T.GetLength(0);
        int last = T.GetLength(1) - 1;
        double[] V = new double[rows];

        for (int i = 0; i < rows; i++)
            V[i] = T[i, last];

        return V;
    }

    static double[] Multiply(double[,] A, double[] B, bool protocol)
    {
        int rows = A.GetLength(0);
        int cols = A.GetLength(1);
        double[] R = new double[rows];

        for (int i = 0; i < rows; i++)
        {
            Console.Write($"x{i + 1} = ");

            for (int j = 0; j < cols; j++)
            {
                R[i] += A[i, j] * B[j];
                Console.Write($"({A[i, j]:F4} * {B[j]:F4})");

                if (j < cols - 1)
                    Console.Write(" + ");
            }

            Console.WriteLine($" = {R[i]:F4}");
        }

        return R;
    }

    static int FindPivotRow(double[,] T, int startRow, int col)
    {
        int rows = T.GetLength(0);
        int best = -1;

        for (int i = startRow; i < rows; i++)
            if (Math.Abs(T[i, col]) > EPS)
            {
                best = i;
                break;
            }

        return best;
    }

    static void SwapRows(double[,] T, int r1, int r2)
    {
        if (r1 == r2) return;

        int cols = T.GetLength(1);

        for (int j = 0; j < cols; j++)
        {
            double temp = T[r1, j];
            T[r1, j] = T[r2, j];
            T[r2, j] = temp;
        }
    }

    static double[,] Copy(double[,] A)
    {
        int rows = A.GetLength(0);
        int cols = A.GetLength(1);
        double[,] C = new double[rows, cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                C[i, j] = A[i, j];

        return C;
    }

    static bool CheckMatrix(double[,] A)
    {
        if (A == null)
        {
            Console.WriteLine("Спочатку введіть матрицю A.");
            return false;
        }

        return true;
    }

    static bool CheckSquare(double[,] A)
    {
        if (A.GetLength(0) != A.GetLength(1))
        {
            Console.WriteLine("Для цієї операції матриця має бути квадратною.");
            return false;
        }

        return true;
    }

    static void PrintMatrix(double[,] M)
    {
        int rows = M.GetLength(0);
        int cols = M.GetLength(1);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                Console.Write($"{M[i, j],10:F4}");
            Console.WriteLine();
        }

        Console.WriteLine();
    }

    static void PrintVector(double[] V)
    {
        Console.WriteLine("Розв'язок:");

        for (int i = 0; i < V.Length; i++)
            Console.WriteLine($"x{i + 1} = {V[i]:F4}");

        Console.WriteLine();
    }


}

// tyt