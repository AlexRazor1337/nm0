// let A = [
//     [2, 1, 1],
//     [1, -1, 0],
//     [3, -1, 2]
// ];

// let B = [2, -2, 2];

// let A = [
//     [2, 5, 4],
//     [1, 3, 2],
//     [2, 10, 9]
// ];

// let B = [30, 150, 110];

let A = [
    [1, 1, 1],
    [4, 2, 1],
    [9, 3, 1]
];

let B = [0, 1, 3];

function cpy(obj) {
    return JSON.parse(JSON.stringify(obj))
}

function det(matrix) {
    if (matrix.length == 1) return matrix[0][0];
    else {
        let nsize = matrix.length - 1;
        let sign = 1;
        let new_matrix = [...Array(nsize)].map(x => Array(nsize).fill(0))
        let det_v = 0;

        for (let c = 0; c < matrix.length; c++) {
            let m = 0;
            let n = 0;
            for (let i = 0; i < matrix.length; i++) {
                for (let j = 0; j < matrix.length; j++) {
                    if (i < nsize && j < nsize) {
                        new_matrix[i][j] = 0;
                    }

                    if (i != 0 && j != c) {
                        new_matrix[m][n] = matrix[i][j];
                        if (n < (matrix.length - 2)) {
                            n++;
                        } else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det_v = det_v + sign * (matrix[0][c] * det(new_matrix));
            sign = -1 * sign;
        }

        return det_v;
    }
}

function is_diagonally_dominant(A) {
    for (let i = 0; i < A.length; i++) {
        if (Math.abs(A[i][i]) < A[i].slice(i).reduce((a,b) => a + b)) {
            return false;
        }
    }

    return true;
}

function is_zeros_on_main_diagonal(A) {
    for (let i = 0; i < A.length; i++) {
        if (A[i][i] == 0) {
            return true;
        }
    }

    return false;
}

function is_symmetrical(A) {

}

function Cramer(A, B) { // TODO проверить сходимость, размерность, ...
    let main_det = det(A);
    let X = new Array(A.length);

    for (let i = 0; i < A.length; i++) {
        let temp = cpy(A);
        for (let j = 0; j < B.length; j++) {
            temp[j][i] = B[j];
        }

        X[i] = det(temp) / main_det;
    }

    return X;
}

function Gauss(A, B) { //TODO Переделать, добавить прове https://www.codesansar.com/numerical-methods/gauss-elimination-method-pseudocode.htm
    let X = new Array(A.length);
    for (let k = 0; k < A.length - 1; k++) {
        // Partial pivot
        let max = Math.abs(A[k][k]);
        let m = k;
        for (let i = k + 1; i < A.length; i++){ //
            if (Math.abs(A[i][k]) > max) {
                max = Math.abs(A[i][k]);
                m = i;
            }
        }

        if (m != k) {
            let b_temp = B[k];
            B[k] = B[m];
            B[m] = b_temp;
            for (let j = k; j < A.length; j++) {
                let A_temp = A[k][j];
                A[k][j] = A[m][j];
                A[m][j] = A_temp;
            }
        }
        // Forward elimination
        for (let i = k + 1; i < A.length; i++) {
            let factor = A[i][k]/A[k][k];
            for (let j = k + 1; j < A.length; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor*B[k];
        }
    }


    // Back substitution
    for (let i = A.length-1; i >= 0; i--) {
        X[i] = B[i];
        for (let j = i + 1; j < A.length; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    return X;
}



function Seidel(A, B, q = 0.001) { // TODO Добавить счетчик итераций?
    let X = new Array(A.length).fill(0);
    while (true) {
        let old_X = cpy(X);
        for (let i = 0; i < A.length; i++) {
            let S = 0;
            for (let j = 0; j < A.length; j++) {
                if (i != j) {
                    S += A[i][j] * X[j];
                }
            }
            X[i] = (1 / A[i][i]) * (B[i] - S)
        }

        let convergence = true;
        for (let i = 0; i < X.length; i++) {
            if (Math.abs(X[i] - old_X[i]) > q) {
                convergence = false;
                break;
            }
        }

        if (convergence) {
            return X;
        }
    }
}

function Gauss_Jordan(A, B) { // TODO Переделать https://www.codesansar.com/numerical-methods/gauss-jordan-method-pseudocode.htm#
    let X = new Array(A.length);
    for (let k = 0; k < A.length - 1; k++) {
        // Partial pivot
        let max = Math.abs(A[k][k]);
        let m = k;
        for (let i = k + 1; i < A.length; i++){ //
            if (Math.abs(A[i][k]) > max) {
                max = Math.abs(A[i][k]);
                m = i;
            }
        }

        if (m != k) {
            let b_temp = B[k];
            B[k] = B[m];
            B[m] = b_temp;
            for (let j = k; j < A.length; j++) {
                let A_temp = A[k][j];
                A[k][j] = A[m][j];
                A[m][j] = A_temp;
            }
        }
        // Forward elimination
        for (let i = k + 1; i < A.length; i++) {
            let factor = A[i][k]/A[k][k];
            for (let j = k + 1; j < A.length; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor*B[k];
        }
    }


    // Back substitution
    for (let i = A.length-1; i >= 0; i--) {
        X[i] = B[i];
        for (let j = i + 1; j < A.length; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    return X;
}

function Jacobi(A, B, q = 0.01) { // TODO Добавить счетчик итераций?
    let X = new Array(A.length).fill(0);
    let iter = 0;
    while (iter < 1000) {
        let old_X = cpy(X);
        for (let i = 0; i < A.length; i++) {
            let S = 0;
            for (let j = 0; j < A.length; j++) {
                if (i != j) {
                    S += A[i][j] * old_X[j];
                }
            }
            X[i] = (1 / A[i][i]) * (B[i] - S)
        }

        let convergence = true;
        for (let i = 0; i < X.length; i++) {
            if (Math.abs(X[i] - old_X[i]) > q) {
                convergence = false;
                break;
            }
        }

        if (convergence) {
            return X;
        }

        iter++;
    }

    return X;
}

function solve(A, B) {
    console.log(Cramer(cpy(A), cpy(B)));
    console.log(Gauss(cpy(A), cpy(B)));
    if (is_diagonally_dominant(A)) {
        console.log(Seidel(cpy(A), cpy(B)));
    } else {
        console.error("Can't do it with Seidel!");
    }
    console.log(Gauss_Jordan(cpy(A), cpy(B)));
    if (is_diagonally_dominant(A) && !is_zeros_on_main_diagonal(A)) {
        console.log(Jacobi(cpy(A), cpy(B)));
    } else {
        console.error("Can't do it with Jacobi, zeros on the main diagonal!");
    }
}
solve(A, B);
