// This is a function used to calculate the coefficients of the polynomial using the method of least squares
function polynomialApproximation(x, y, degree) {
    // this is the number of data points
    const n = x.length;

    // this initializes matrices to store sums of powers of x and the sum of x^i * y
    const X = [];
    const Y = [];

    // this create the matrix of sums of powers of x
    for (let i = 0; i <= degree * 2; i++) {
        X[i] = 0;
        for (let j = 0; j < n; j++) {
            X[i] += Math.pow(x[j], i);
        }
    }

    // this creates the matrix of sums of x^i * y
    for (let i = 0; i <= degree; i++) {
        Y[i] = 0;
        for (let j = 0; j < n; j++) {
            Y[i] += Math.pow(x[j], i) * y[j];
        }
    }

    // this is to initialize matrix A and vector B
    const A = [];
    const B = [];

    // thus fills matrix A with values from the sums of powers of x
    for (let i = 0; i <= degree; i++) {
        A[i] = [];
        for (let j = 0; j <= degree; j++) {
            A[i][j] = X[i + j];
        }
        B[i] = Y[i];
    }

    // Used to solve the system of linear equations using Gaussian elimination
    const coefficients = gaussianElimination(A, B);

    return coefficients;
}

// Function to solve a system of linear equations using Gaussian elimination
function gaussianElimination(A, B) {
    const n = B.length;

    // Forward elimination
    for (let i = 0; i < n; i++) {
        // Find the pivot element
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Swap the rows to put the pivot element on the diagonal
        const temp = A[i];
        A[i] = A[maxRow];
        A[maxRow] = temp;

        const tempB = B[i];
        B[i] = B[maxRow];
        B[maxRow] = tempB;

        // Make all rows below this one 0 in the current column
        for (let k = i + 1; k < n; k++) {
            const c = -A[k][i] / A[i][i];
            for (let j = i; j < n; j++) {
                if (i === j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
            B[k] += c * B[i];
        }
    }

    // Backward substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        x[i] = B[i] / A[i][i];
        for (let k = i - 1; k >= 0; k--) {
            B[k] -= A[k][i] * x[i];
        }
    }

    return x;
}

// Example usage of the polynomial approximation function
const x = [1, 2, 3, 4, 5]; // X values
const y = [1, 4, 9, 16, 25]; // Y values corresponding to X (y = x^2 as an example)

// Degree of the polynomial to approximate (e.g., 2 for quadratic)
const degree = 2;

// Get the coefficients of the polynomial
const coefficients = polynomialApproximation(x, y, degree);

console.log("Coefficients of the polynomial approximation:", coefficients);

// The result should approximate the polynomial y = x^2, so the coefficients should be close to [0, 0, 1].
