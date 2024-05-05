#include <stdio.h>
#include <stdlib.h>

#define MAX_ROWS 10
#define MAX_COLS 10

float A[MAX_ROWS][MAX_COLS], B[MAX_ROWS][MAX_COLS], C[MAX_ROWS][MAX_COLS], D[MAX_ROWS][MAX_COLS];
int num_rows, num_cols;

void input_matrix(float matrix[MAX_ROWS][MAX_COLS]) {
    printf("Enter the elements of the matrix:\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            printf("Enter element [%d][%d]: ", i, j);
            scanf("%f", &matrix[i][j]);
        }
    }
}

void display_matrix(float matrix[MAX_ROWS][MAX_COLS]) {
    printf("Matrix:\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            printf("%.2f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void add_matrices() {
    printf("Adding matrices A and B:\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    display_matrix(C);
}

void subtract_matrices() {
    printf("Subtracting matrices A and B:\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    display_matrix(C);
}

void multiply_matrices() {
    printf("Multiplying matrices A and B:\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            C[i][j] = 0;
            for (int k = 0; k < num_cols; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    display_matrix(C);
}

float determinant(float matrix[MAX_ROWS][MAX_COLS], int n) {
    float det = 0;
    if (n == 1)
        return matrix[0][0];
    else if (n == 2)
        return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
    else {
        float submatrix[MAX_ROWS][MAX_COLS];
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (x % 2 == 0 ? 1 : -1) * matrix[0][x] * determinant(submatrix, n - 1);
        }
    }
    return det;
}

void cofactor(float matrix[MAX_ROWS][MAX_COLS], float temp[MAX_ROWS][MAX_COLS], int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = matrix[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

void adjoint(float matrix[MAX_ROWS][MAX_COLS], float adj[MAX_ROWS][MAX_COLS]) {
    if (num_rows == 1) {
        adj[0][0] = 1;
        return;
    }
    float temp[MAX_ROWS][MAX_COLS];
    int sign = 1;
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            cofactor(matrix, temp, i, j, num_rows);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = (sign) * (determinant(temp, num_rows - 1));
        }
    }
}

int inverse(float matrix[MAX_ROWS][MAX_COLS], float inverse_matrix[MAX_ROWS][MAX_COLS]) {
    float det = determinant(matrix, num_rows);
    if (det == 0) {
        printf("Inverse does not exist as determinant is zero.\n");
        return 0;
    }
    float adj[MAX_ROWS][MAX_COLS];
    adjoint(matrix, adj);
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            inverse_matrix[i][j] = adj[i][j] / det;
        }
    }
    return 1;
}

void LU_decomposition(float matrix[MAX_ROWS][MAX_COLS]) {
    float L[MAX_ROWS][MAX_COLS];
    float U[MAX_ROWS][MAX_COLS];
    
    // Implement LU decomposition
    for (int i = 0; i < num_rows; i++) {
        for (int k = i; k < num_cols; k++) {
            // Upper Triangular
            U[i][k] = matrix[i][k];
            for (int j = 0; j < i; j++) {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
        for (int k = i; k < num_cols; k++) {
            // Lower Triangular
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
                L[k][i] = matrix[k][i] / U[i][i];
                for (int j = 0; j < i; j++) {
                    L[k][i] -= ((L[k][j] * U[j][i]) / U[i][i]);
                }
            }
        }
    }
    
    // Print L and U matrices if needed
    printf("L matrix:\n");
    display_matrix(L);
    printf("U matrix:\n");
    display_matrix(U);
}

void solve_linear_equations(float L[MAX_ROWS][MAX_COLS], float U[MAX_ROWS][MAX_COLS], float B[MAX_ROWS]) {
    float Y[MAX_ROWS];
    float X[MAX_ROWS];

    // Solve LY = B for Y using forward substitution
    for (int i = 0; i < num_rows; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
    }

    // Solve UX = Y for X using backward substitution
    for (int i = num_rows - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < num_cols; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }

    // Display the solution
    printf("Solution:\n");
    for (int i = 0; i < num_rows; i++) {
        printf("X[%d] = %.2f\n", i, X[i]);
    }
}

int main() {
    int choice;
    
    printf("Enter the number of rows for matrices: ");
    scanf("%d", &num_rows);
    printf("Enter the number of columns for matrices: ");
    scanf("%d", &num_cols);
    
    input_matrix(A);
    input_matrix(B);
    
    while (1) {
        printf("\nChoose an option:\n");
        printf("1. Add new matrix\n");
        printf("2. Perform operations\n");
        printf("3. Find determinant\n");
        printf("4. Find inverse\n");
        printf("5. LUP decomposition\n");
        printf("6. Solve linear equations\n");
        printf("9. Exit\n");
        printf("Enter your choice: ");
        scanf("%d", &choice);
        
        switch (choice) {
            case 1:
                // Add new matrix
                break;
            case 2:
                printf("Choose operation:\n");
                printf("1. Add matrices\n");
                printf("2. Subtract matrices\n");
                printf("3. Multiply matrices\n");
                printf("Enter your choice: ");
                scanf("%d", &choice);
                switch (choice) {
                    case 1:
                        add_matrices();
                        break;
                    case 2:
                        subtract_matrices();
                        break;
                    case 3:
                        multiply_matrices();
                        break;
                    default:
                        printf("Invalid choice! Please try again.\n");
                }
                break;
            case 3:
                printf("Determinant of matrix A: %.2f\n", determinant(A, num_rows));
                break;
            case 4:
                {
                    float inverse_A[MAX_ROWS][MAX_COLS];
                    if (inverse(A, inverse_A)) {
                        printf("Inverse of matrix A:\n");
                        display_matrix(inverse_A);
                    }
                }
                break;
            case 5:
                LU_decomposition(A);
                break;
            case 6:
                {
                    float L[MAX_ROWS][MAX_COLS];
                    float U[MAX_ROWS][MAX_COLS];
                    LU_decomposition(A);
                    solve_linear_equations(L, U, B[0]);
                }
                break;
            case 9:
                printf("Exiting...\n");
                return 0;
            default:
                printf("Invalid choice! Please try again.\n");
        }
    }

    return 0;
}
