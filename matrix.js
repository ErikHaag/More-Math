/*
More Math library by Erik Haag version 2.0
https://github.com/ErikHaag/More-Math/
Dependencies: moreMathCore.js, rational.js
*/

class Matrix {
    //the determinate
    #determinate = null;
    //the row echelon form
    #ref = null;
    //whether the determinate method was called
    #determinateCalled = false;
    constructor(columns, indices) {
        if (indices === "identity") {
            let m = [];
            for (let i = 0; i < columns; i++) {
                let r = [];
                for (let j = 0; j < columns; j++) {
                    if (i == j) {
                        r.push(new Rational(1n));
                    } else {
                        r.push(new Rational(0n));
                    }
                }
                m.push(r);
            }
            this.indices = m;
            this.columns = columns;
            this.rows = columns;
        } else {
            let indicesLength = indices.length;
            if (indicesLength % columns == 0) {
                let m = [];
                let r = [];
                for (let i = 0; i < indicesLength; i++) {
                    if (typeof indices[i] == "number") {
                        indices[i] = new Rational(BigInt(Math.floor(indices[i])));
                    } else if (typeof indices[i] == "bigint") {
                        indices[i] = new Rational(indices[i]);
                    }
                    r.push(indices[i]);
                    if (r.length == columns) {
                        m.push(r);
                        r = [];
                    }
                }
                this.indices = m;
                this.columns = columns;
                this.rows = indicesLength / columns;
            } else {
                return new Error("The last row of the matrix is unfilled!");
            }
        }
    }
    toLatex() {
        //creates the LaTeX code for this matrix 
        let b = [];
        for (let i = 0; i < this.rows; i++) {
            let r = []
            for (let j = 0; j < this.columns; j++) {
                r.push(this.indices[i][j].toLatex());
            }
            b.push(r.join("&"));
        }
        return "\\begin{bmatrix}" + b.join("\\\\") + "\\end{bmatrix}";
    }
    swapRows(a, b) {
        [this.indices[a], this.indices[b]] = [this.indices[b], this.indices[a]];
        if (this.#determinateCalled) {
            this.#determinate.mult(-1n); 
        }
    }
    addRow(source, destination, scale = 1n) {
        for (let i = 0; i < this.columns; i++) {
            let r = this.indices[source][i].clone();
            r.mult(scale);
            this.indices[destination][i].add(r);
        }
    }
    scale(scale) {
        for (let i = 0; i < this.rows; i++) {
            this.scaleRow(i, scale);
        }
    }
    scaleRow(row, scale) {
        for (let i = 0; i < this.columns; i++) {
            this.indices[row][i].mult(scale);
        }
        if (this.#determinateCalled) {
            this.#determinate.mult(scale);
        }
    }
    transpose() {
        let t = [];
        for (let i = 0; i < this.columns; i++) {
            let r = []
            for (let j = 0; j < this.rows; j++) {
                r.push(this.indices[j][i]);
            }
            t.push(r);
        }
        this.indices = t;
        [this.columns, this.rows] = [this.rows, this.columns];
    }
    clone() {
        let m = [];
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                m.push(this.indices[i][j].clone());
            }
        }
        return new Matrix(this.columns, m.flat());
    }
    isSquare() {
        return this.columns == this.rows;
    }
    #rowReduction() {
        let h = 0;
        let k = 0;
        let neg = false;
        let m = this.clone();
        while ((h < m.rows) && (k < m.columns)) {
            let iMax = h;
            let max = new Rational(0n);
            for (let i = h; i < this.rows; i++) {
                let c = m.indices[i][k].clone();
                c.numerator = BigMathJS.abs(c.numerator);
                c.sub(max)
                if (c.numerator > 0n) {
                    max = m.indices[i][k].clone();
                    max.numerator = BigMathJS.abs(max.numerator);
                    iMax = i;
                }
            }
            if (max.numerator == 0n) {
                k += 1;
            } else {
                m.swapRows(h, iMax);
                if (h != iMax) {
                    neg = !neg;
                }
                for (let i = h + 1; i < m.rows; i++) {
                    let f = m.indices[i][k].clone();
                    f.div(m.indices[h][k]);
                    m.indices[i][k] = new Rational(0n);
                    for (let j = k + 1; j < m.columns; j++) {
                        let s = m.indices[h][j].clone();
                        s.mult(f);
                        m.indices[i][j].sub(s);
                    }
                }
                h += 1;
                k += 1;
            }
        }
        this.#ref = m;
        if (this.isSquare()) {
            let p = new Rational((neg ? -1n : 1n));
            for (let i = 0; i < m.columns; i++) {
                p.mult(m.indices[i][i]);
            }
            this.#determinate = p.clone();
        } else {
            this.#determinate = new Rational(0n);
        }
    }
    gaussianElimination() {
        if (!this.#determinateCalled) {
            this.#rowReduction();
        }
        this.indices = structuredClone(this.#ref);
        this.#determinateCalled = false;
    }
    determinate() {
        if (!this.isSquare()) {
            return new Error("Can't take determinate of non-square matrix");
        }
        if (!this.#determinateCalled) {
            this.#rowReduction();
            this.#determinateCalled = true;
        }
        return this.#determinate;
    }
    inverse() {
        if (this.isSquare()) {
            if (this.determinate().numerator != 0n) {
                let M = this.clone();
                let hMax = M.rows;
                let kMax = M.columns;
                let I = new Matrix(M.rows, "identity");
                M = M.augment(I);
                let h = 0;
                let k = 0;
                while ((h < hMax) && (k < kMax)) {
                    let iMax = h;
                    let max = new Rational(0n);
                    for (let i = h; i < M.rows; i++) {
                        let c = M.indices[i][k].clone();
                        c.numerator = BigMathJS.abs(c.numerator);
                        c.sub(max)
                        if (c.numerator > 0n) {
                            max = M.indices[i][k].clone();
                            max.numerator = BigMathJS.abs(max.numerator);
                            iMax = i;
                        }
                    }
                    if (max.numerator == 0n) {
                        k += 1;
                    } else {
                        M.swapRows(h, iMax);
                        for (let i = 0; i < M.rows; i++) {
                            if (i == h) {
                                continue;
                            }
                            let f = M.indices[i][k].clone();
                            f.div(M.indices[h][k]);
                            M.indices[i][k] = new Rational(0n);
                            for (let j = k + 1; j < M.columns; j++) {
                                let s = M.indices[h][j].clone();
                                s.mult(f);
                                M.indices[i][j].sub(s);
                            }
                        }
                        h += 1;
                        k += 1;
                    }
                }
                for (let i = 0; i < hMax; i++) {
                    M.scaleRow(i, M.indices[i][i].cloneInverse());
                }
                let m = [];
                for (let i = 0; i < hMax; i++) {
                    for (let j = kMax; j < M.columns; j++) {
                        m.push(M.indices[i][j].clone());
                    }
                }
                return new Matrix(kMax, m.flat());
            } else {
                return new Error("Matrix has a determinate of zero.")
            }
        } else {
            return new Error("Matrix isn't square.")
        }
    }
    augment(B) {
        if (this.rows == B.rows) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                let r = [];
                for (let j = 0; j < this.columns; j++) {
                    r.push(this.indices[i][j].clone());
                }
                for (let j = 0; j < B.columns; j++) {
                    r.push(B.indices[i][j].clone());
                }
                m.push(r);
            }
            this.indices = m;
            this.columns += B.columns;
            this.#determinateCalled = false;
        } else {
            return new Error("Unable to augment due too inequal rows.");
        }
    }
    add(B) {
        if (this.columns == B.columns && this.rows == B.rows) {
            for (let i = 0; i < this.rows; i++){
                for (let j = 0; j < this.columns; j++) {
                    this.indices[i][j].add(B.indices[i][j]);
                }
            }
            this.#determinateCalled = false;
        } else {
            return new Error("Matrices must be the same size");
        }
    }
    dotProduct(B) {
        //add element-wise products
        if (this.columns == 1 && B.columns == 1 && this.rows == B.rows) {
            let s = new Rational(0n);
            for (let i = 0; i < this.rows; i++) {
                let p = this.indices[i][0].clone();
                p.mult(B.indices[i][0]);
                s.add(p);
            }
            return s;
        } else {
            return new Error("Can only take dot product with vectors of equal dimension.");
        }
    }
    hadamardProduct(B) {
        //element-wise product
        if (this.columns == B.columns && this.rows == B.rows) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                let r = [];
                for (let j = 0; j < this.columns; j++) {
                    let p = this.indices[i][j].clone();
                    p.mult(B.indices[i][j]);
                    r.push(p);
                }
                m.push(r);
            }
            this.indices = m;
            this.#determinateCalled = false;
        } else {
            return new Error("Can only take Hadamard product with matrices of the same size.");
        }
    }
    kroneckerProduct(B) {
        //all the pairwise products!
        let m = [];
        for (let aI = 0; aI < this.rows; aI++) {
            for (let bI = 0; bI < B.rows; bI++) {
                let r = [];
                for (let aJ = 0; aJ < this.columns; aJ++) {
                    for (let bJ = 0; bJ < B.columns; bJ++) {                        
                        let p = this.indices[aI][aJ].clone();
                        p.mult(B.indices[bI][bJ]);
                        r.push(p);
                    }
                }
                m.push(r);
            }
        }
        this.indices = m;
        this.columns *= B.columns;
        this.rows *= B.rows;
        this.#determinateCalled = false;
    }
    outerProduct(B) {
        //multiplication table
        if (this.columns == 1 && B.columns == 1) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                let r = [];
                for (let j = 0; j < B.rows; j++) {
                    let p = this.indices[i][0].clone();
                    p.mult(B.indices[j][0]);
                    r.push(p);
                }
                m.push(r);
            }
            this.indices = m;
            this.columns = B.rows;
            this.#determinateCalled = false;
        } else {
            return new Error("Can only take outer product with vectors. (A matrix with 1 column)");
        }
    }
    product(B) {
        //standard issue product
        if (this.columns == B.rows) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                let r = [];
                for (let j = 0; j < B.columns; j++) {
                    let sum = new Rational(0n);
                    for (let k = 0; k < this.columns; k++) {
                        let p = this.indices[i][k].clone();
                        p.mult(B.indices[k][j]);
                        sum.add(p);
                    }
                    r.push(sum);
                }
                m.push(r);
            }
            this.indices = m;
            this.columns = B.columns;
            this.#determinateCalled = false;
        } else {
            return new Error("Matrices have incompatible sizes.");
        }
    }
}