import java.util.Arrays;

public class Main {

    // Metode untuk menyelesaikan menggunakan metode matriks invers.
    public static double[] solveWithInverseMatrix(double[][] A, double[] B) {
        double[][] AInverse = inverseMatrix(A);
        if (AInverse == null) {
            System.out.println("Matrix is singular, cannot find inverse.");
            return null;
        }
        return multiplyMatrixVector(AInverse, B);
    }

    // Metode untuk menemukan invers dari sebuah matriks.
    public static double[][] inverseMatrix(double[][] A) {
        int n = A.length;
        double[][] identity = new double[n][n];
        double[][] AInverse = new double[n][n];

        // Initialize identity matrix
        for (int i = 0; i < n; i++) {
            identity[i][i] = 1;
        }

        // Perform Gauss-Jordan elimination
        for (int i = 0; i < n; i++) {
            // Partial pivot
            if (A[i][i] == 0) {
                return null; // Singular matrix
            }

            double pivot = A[i][i];
            for (int j = 0; j < n; j++) {
                A[i][j] /= pivot;
                identity[i][j] /= pivot;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = A[k][i];
                    for (int j = 0; j < n; j++) {
                        A[k][j] -= factor * A[i][j];
                        identity[k][j] -= factor * identity[i][j];
                    }
                }
            }
        }

        return identity;
    }

    // Metode untuk mengalikan matriks dengan vektor.
    public static double[] multiplyMatrixVector(double[][] A, double[] B) {
        int n = A.length;
        int m = B.length;
        double[] result = new double[n];

        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < m; j++) {
                sum += A[i][j] * B[j];
            }
            result[i] = sum;
        }

        return result;
    }

    // Metode untuk menyelesaikan menggunakan metode dekomposisi Gauss LU.
    public static double[] solveWithGaussLU(double[][] A, double[] B) {
        int n = A.length;
        double[][] LU = new double[n][n];
        double[] x = new double[n];

        // Decompose matrix A into LU
        decomposeLU(A, LU);

        // Solve LY = B using forward substitution
        double[] Y = forwardSubstitution(LU, B);

        // Solve UX = Y using backward substitution
        backwardSubstitution(LU, Y, x);

        return x;
    }

    // metode dari dekomposisi matriks ke LU
    public static void decomposeLU(double[][] A, double[][] LU) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
            LU[i][i] = 1;
            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += LU[i][k] * LU[k][j];
                }
                LU[i][j] = A[i][j] - sum;
            }
            for (int j = i + 1; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += LU[j][k] * LU[k][i];
                }
                LU[j][i] = (A[j][i] - sum) / LU[i][i];
            }
        }
    }

    // Metode untuk substitusi maju.
    public static double[] forwardSubstitution(double[][] LU, double[] B) {
        int n = LU.length;
        double[] Y = new double[n];

        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += LU[i][j] * Y[j];
            }
            Y[i] = (B[i] - sum) / LU[i][i];
        }

        return Y;
    }

    // Metode untuk substitusi mundur.
    public static void backwardSubstitution(double[][] LU, double[] Y, double[] x) {
        int n = LU.length;

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += LU[i][j] * x[j];
            }
            x[i] = Y[i] - sum;
        }
    }

    // Metode untuk menyelesaikan menggunakan metode dekomposisi Crout.
    public static double[] solveWithCrout(double[][] A, double[] B) {
        int n = A.length;
        double[][] LU = new double[n][n];
        double[] x = new double[n];

        // Decompose matrix A into LU
        decomposeCrout(A, LU);

        // Solve LY = B using forward substitution
        double[] Y = forwardSubstitution(LU, B);

        // Solve UX = Y using backward substitution
        backwardSubstitution(LU, Y, x);

        return x;
    }

    // Metode untuk mendekomposisi matriks menjadi LU menggunakan metode Crout.
    public static void decomposeCrout(double[][] A, double[][] LU) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
            LU[i][i] = 1;
            for (int j = i; j < n; j++) {
                double sum1 = 0;
                for (int k = 0; k < i; k++) {
                    sum1 += LU[i][k] * LU[k][j];
                }
                LU[i][j] = A[i][j] - sum1;
            }
            for (int j = i + 1; j < n; j++) {
                double sum2 = 0;
                for (int k = 0; k < i; k++) {
                    sum2 += LU[j][k] * LU[k][i];
                }
                LU[j][i] = (A[j][i] - sum2) / LU[i][i];
            }
        }
    }

    // Method to display matrix
    public static void displayMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            System.out.println(Arrays.toString(row));
        }
    }

    // Main method
    public static void main(String[] args) {
        double[][] A = {{10, 3, 5}, {10, -1, 0}, {20, 2, 7}};
        double[] B = {1, 8, -10};

        System.out.println("pakai metode Inverse Matrix :");
        double[] xInverse = solveWithInverseMatrix(A, B);
        System.out.println("Solusi: " + Arrays.toString(xInverse));

        System.out.println("\npakai metode Gauss LU Decomposition :");
        double[] xGaussLU = solveWithGaussLU(A, B);
        System.out.println("Solusi: " + Arrays.toString(xGaussLU));

        System.out.println("\npakai metode Crout Decomposition :");
        double[] xCrout = solveWithCrout(A, B);
        System.out.println("Solusi: " + Arrays.toString(xCrout));
    }
}
