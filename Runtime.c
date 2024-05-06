#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_SIZE 10

void input_matrix(float A[MAX_SIZE][MAX_SIZE], int n) {
    printf("Enter the elements of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("Enter element [%d][%d]: ", i, j);
            scanf("%f", &A[i][j]);
        }
    }
}

void generate_random_matrix(float A[MAX_SIZE][MAX_SIZE], int n) {
    srand(time(NULL)); // Seed the random number generator
    printf("Randomly generated matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = rand() % 101; // Generate random number between 0 and 100
            //printf("%.2f\t", A[i][j]);
        }
        //printf("\n");
    }
}

void display_matrix(float A[MAX_SIZE][MAX_SIZE], int n) {
    printf("Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.2f\t", A[i][j]);
        }
        printf("\n");
    }
}

void swap(float *a, float *b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

void LUP_decomposition(float A[MAX_SIZE][MAX_SIZE], float L[MAX_SIZE][MAX_SIZE], float U[MAX_SIZE][MAX_SIZE], float P[MAX_SIZE][MAX_SIZE], int n) {
    // Initialize matrices L and P as identity matrices and U as A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                L[i][j] = 1.0;
                P[i][j] = 1.0;
            } else {
                L[i][j] = 0.0;
                P[i][j] = 0.0;
            }
            U[i][j] = A[i][j];
        }
    }

    // Start measuring runtime
    clock_t start = clock();

    // Perform Gaussian elimination with partial pivoting
    for (int k = 0; k < n - 1; k++) {
        int pivot_row = k;
        float max = abs(U[k][k]);
        
        // Find the row with maximum pivot element
        for (int i = k + 1; i < n; i++) {
            if (abs(U[i][k]) > max) {
                max = abs(U[i][k]);
                pivot_row = i;
            }
        }

        // Swap rows in U matrix
        if (pivot_row != k) {
            for (int j = 0; j < n; j++) {
                swap(&U[k][j], &U[pivot_row][j]);
            }
            // Swap rows in P matrix
            for (int j = 0; j < n; j++) {
                swap(&P[k][j], &P[pivot_row][j]);
            }
            // Track row swaps in L matrix
            for (int i = 0; i < k; i++) {
                swap(&L[k][i], &L[pivot_row][i]);
            }
        }

        for (int i = k + 1; i < n; i++) {
            float factor = U[i][k] / U[k][k];
            L[i][k] = factor;
            for (int j = k; j < n; j++) {
                U[i][j] -= factor * U[k][j];
            }
        }
    }

    // Stop measuring runtime
    clock_t end = clock();
    double runtime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Runtime for LUP decomposition: %.6f seconds\n", runtime);
}

int main() {
    int n, choice;

    printf("Enter the size of the square matrix: ");
    scanf("%d", &n);

    float A[MAX_SIZE][MAX_SIZE], L[MAX_SIZE][MAX_SIZE], U[MAX_SIZE][MAX_SIZE], P[MAX_SIZE][MAX_SIZE];

    printf("Choose an option:\n");
    printf("1. Input matrix manually\n");
    printf("2. Generate random matrix\n");
    printf("Enter your choice: ");
    scanf("%d", &choice);

    switch (choice) {
        case 1:
            input_matrix(A, n);
            break;
        case 2:
            generate_random_matrix(A, n);
            break;
        default:
            printf("Invalid choice! Exiting...\n");
            return 1;
    }

    LUP_decomposition(A, L, U, P, n);

    // printf("L matrix:\n");
    // display_matrix(L, n);
    // printf("U matrix:\n");
    // display_matrix(U, n);
    // printf("P matrix:\n");
    // display_matrix(P, n);

    return 0;
}
