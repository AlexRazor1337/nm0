
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

        if (Math.abs(A[i][i]) < A[i].reduce((a,b) => a + b) - A[i][i]) {
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
    let combined = new Array(A.length);
    for (let i = 0; i < A.length; i++) {
        combined[i] = [...A[i], B[i]];
    }

    for (let i = 0; i < A.length - 1; i++) {
        for (let j = i + 1; j < A.length; j++) {
            let ratio = combined[j][i] / combined[i][i];

            for (let k = 0; k < A.length + 1; k++) {
                combined[j][k] = combined[j][k] - ratio * combined[i][k];
            }
        }
    }

    X[A.length - 1] = combined[A.length - 1][A.length] / combined[A.length - 1][A.length - 1];

    for (let i = A.length - 2; i >= 0; i--) {
        X[i] = combined[i][A.length];
        for (let j = i + 1; j < A.length; j++) {
            X[i] = X[i] - combined[i][j] * X[j];
        }

        X[i] = X[i] / combined[i][i];
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
    let combined = new Array(A.length);
    for (let i = 0; i < A.length; i++) {
        combined[i] = [...A[i], B[i]];
    }

    for (let i = 0; i < A.length; i++) {
        for (let j = 0; j < A.length; j++) {
            if (i != j) {
                let ratio = combined[j][i] / combined[i][i];
                for (let k = 0; k < A.length + 1; k++) {
                    combined[j][k] = combined[j][k] - ratio * combined[i][k];
                }
            }
        }
    }

    X[A.length - 1] = combined[A.length - 1][A.length] / combined[A.length - 1][A.length - 1];

    for (let i = A.length - 2; i >= 0; i--) {
        X[i] = combined[i][A.length];
        for (let j = i + 1; j < A.length; j++) {
            X[i] = X[i] - combined[i][j] * X[j];
        }

        X[i] = X[i] / combined[i][i];
    }

    return X;
}

function Jacobi(A, B, q = 0.001) { // TODO Добавить счетчик итераций?
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

function check_solution(A, B, X, E = 0.001) {
    for (let i = 0; i < A.length; i++) {
        let sum = 0;
        for (let j = 0; j < A.length; j++) {
            sum += A[i][j] * X[j];
        }

        if (Math.abs(sum - B[i]) > E) return false;
    }

    return true;
}

function solve(A, B, method) {
    if (det(A) == 0) {
        alert("Матрица вырожденная!");
        return false;
    }

    let X = false;
    switch (method) {
        case 1:
            if (A.length < 5) {
                X = Cramer(cpy(A), cpy(B));
            } else {
                alert("Слишком большая матрица для этого метода!");
            }
            break;
        case 2:
            if (!is_zeros_on_main_diagonal(A)) {
                X = Gauss(cpy(A), cpy(B));
            } else {
                alert("Нули на главной диагонали!");
            }
            break;
        case 3:
            if (is_diagonally_dominant(A)) {
                X = Seidel(cpy(A), cpy(B));
            } else {
                alert("Невозможно решить Зейделем!");
            }
            break;
        case 4:
            if (!is_zeros_on_main_diagonal(A)) {
                X = Gauss_Jordan(cpy(A), cpy(B));
            } else {
                alert("Нули на главной диагонали!");
            }
            break;
        case 5:
            if (!is_zeros_on_main_diagonal(A) && is_diagonally_dominant(A)) {
                X = Jacobi(cpy(A), cpy(B));
            } else {
                alert("Невозможно использовать метод Якоби!");
            }
        default:
            break;
    }

    return X;
}



function createInput(dim, A = null, B = null) {
    if (document.getElementById('inp')) {
        document.getElementById('inp').remove()
    }

    let tbl = document.createElement('table');
    for (let i = 0; i < dim; i++) {
        let row = document.createElement('tr');
        for (let j = 0; j < dim; j++) {
            let cell = document.createElement('td');
            let inp = document.createElement('input');
            if (A) {
                inp.value = A[i][j];
            }
            inp.id = i + "-" + j
            inp.type = "number"
            let lab = document.createElement('label');
            lab.setAttribute("for", i + "-" + j);
            lab.innerText = " x" + (j+1)
            if (j != dim - 1) {
                lab.innerText += " + "
            }
            cell.appendChild(inp);
            cell.appendChild(lab);
            row.appendChild(cell);
        }

        let cell = document.createElement('td');
        let inp = document.createElement('input');
        inp.id = i + "B"
        inp.type = "number"
        if (B) {
            inp.value = B[i];
        }
        let lab = document.createElement('label');
        lab.innerText = " = "
        cell.appendChild(lab);
        cell.appendChild(inp);
        row.appendChild(cell);

        tbl.appendChild(row);
    }

    tbl.id = 'inp'
    document.getElementById('container').appendChild(tbl)
}

function solver(A, B, dim) {
    let X = solve(A, B, parseInt(document.getElementById('2').value));
    if (X) {
        let ans = document.createElement('p');
        ans.innerText = "Корни СЛАУ: "
        for (let i = 0; i < dim; i++) {
            ans.innerText += 'x' + (i + 1) + " = " + X[i];
            if (i < dim - 1) {
                ans.innerText += ',  ';
            }
        }

        document.getElementById('outter_container').appendChild(ans)
        let ch = document.createElement('p');
        ch.innerText = "Проверить решение?";
        ch.id = 'check?'
        document.getElementById('outter_container').appendChild(ch)
        let y = document.createElement('button');
        y.id = "y";
        y.onclick = function() {
            document.getElementById('container').remove()
            document.getElementById('solve').remove()
            document.getElementById('check?').remove()
            document.getElementById('y').remove()
            document.getElementById('n').remove()
            if (check_solution(A, B, X)) {
                let tr = document.createElement('p');
                tr.innerText = "Решение верное!"
                document.getElementById('outter_container').appendChild(tr)
            }
            for (let i = 0; i < dim; i++) {
                let chp = document.createElement('p');
                for (let j = 0; j < dim; j++) {
                    chp.innerText += A[i][j] + " * " + X[j];
                    if (j < dim - 1) {
                        chp.innerText += " + "
                    }

                }
                chp.innerText += " = " + B[i]
                document.getElementById('outter_container').appendChild(chp)
            }
            document.getElementById('outter_container')
        }
        y.innerText = 'Да';

        let n = document.createElement('button');
        n.id = "n"
        n.onclick = function() {
            location.reload();
        }
        n.innerText = 'Нет';

        document.getElementById('outter_container').appendChild(y);
        document.getElementById('outter_container').appendChild(n);
    }
}

window.addEventListener('DOMContentLoaded', (event) => {
    for (let i = 2; i <= 16; i++) {
        let opt = document.createElement('option');
        opt.id = "opt_" + i;
        opt.value = i;
        opt.innerText = i + "x" + i;
        document.getElementById('1').appendChild(opt)
    }
    createInput(2)
    document.getElementById('1').onchange = function() {
        createInput(document.getElementById('1').value)
    }

    document.getElementById('solve').onclick = function() {
        let A, B, dim;
        if (document.getElementById('file-input').files && document.getElementById('file-input').files[0]) {
            var reader = new FileReader();
            reader.readAsText(document.getElementById('file-input').files[0]);

            reader.onload = function(e) {
                var contents = JSON.parse(e.target.result);

                A = contents[0];
                B = contents[1];
                dim = A.length;
                createInput(dim, A, B)
                solver(A, B, dim)
            };
        } else {
            dim = parseInt(document.getElementById('1').value);
            A = new Array(dim);
            B = new Array(dim);

            for (let i = 0; i < dim; i++) {
                A[i] = new Array(A.dim);
                for (let j = 0; j < dim; j++) {
                    let val = document.getElementById(i + "-" + j).value;
                    if (val == "") return;

                    A[i][j] = parseInt(val);
                }

                let val = document.getElementById(i + "B").value;
                if (val == "") return;

                B[i] = parseInt(val);
            }
            solver(A, B, dim)
        }
    }
});
// solve(A, B);
