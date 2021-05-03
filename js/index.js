let A = [
    [2, 1, 1],
    [1, -1, 0],
    [3, -1, 2]
];

let B = [2, -2, 2];

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

function Cramer(A, B) { // проверить сходимость, размерность, ...
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

function Gauss(A, B) {
    let X = new Array(A.length);
    let tmp = cpy(B);
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
            let b_temp = tmp[k];
            tmp[k] = tmp[m];
            tmp[m] = b_temp;
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
            tmp[i] -= factor*tmp[k];
        }
    }

    // Back substitution
    for (let i = A.length-1; i >= 0; i--) {
        X[i] = tmp[i];
        for (let j = i + 1; j < A.length; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    return X;
}


console.log(Cramer(cpy(A), cpy(B)));
console.log(Gauss(cpy(A), cpy(B)));
