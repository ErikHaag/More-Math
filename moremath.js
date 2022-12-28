/*
More Math library by Erik Haag version 1.0

Adds matrices and rationals to JavaScript, and some other things!
*/
class Rational {
    constructor(numerator, denominator = 1n) {
        if (typeof numerator == "bigint" && typeof denominator == "bigint") {
            if (0 < denominator) {
                this.numerator = numerator;
                this.denominator = denominator;
            } else {
                return new Error("Denomimator must be greater than 0.");
            }
        } else {
            return new Error("Numerator and Denominator must be bigInts.");
        }
    }
    toLatex() {
        if (this.denominator == 1) {
            return this.numerator.toString();
        } else {
            return "\\frac{" + this.numerator + "}{" + this.denominator + "}"
        }
    }
    clone() {
        this.simplify();
        return new Rational(this.numerator, this.denominator);
    }
    simplify() {
        let factor = MathJS.gcd(this.numerator, this.denominator);
        this.numerator /= factor;
        this.denominator /= factor;
        if (this.denominator < 0) {
            this.numerator *= -1n;
            this.denominator *= -1n;
        }
    }
    value() {
        return Number(this.numerator) / Number(this.denominator);
    }
    cloneInverse() {
        if (this.numerator >= 0) {
            return new Rational(this.denominator, this.numerator);
        } else {
            return new Rational(-this.denominator, -this.numerator);
        }
    }
    add(B) {
        this.numerator = this.numerator * B.denominator + this.denominator * B.numerator;
        this.denominator *= B.denominator;
        this.simplify();
    }
    sub(B) {
        this.numerator = this.numerator * B.denominator - this.denominator * B.numerator;
        this.denominator *= B.denominator;
        this.simplify();
    }
    mult(B) {
        this.numerator *= B.numerator;
        this.denominator *= B.denominator;
        this.simplify();
    }
    div(B) {
        if (B.numerator == 0) {
            return new Error("Division by zero!")
        } else {
            this.numerator *= B.denominator;
            this.denominator *= B.numerator;
            this.simplify();
        }
    }
    pow(B) {
        if (typeof B == "bigint"){
        this.numerator **= B;
        this.denominator **= B;
        } else {
            return new Error("argument must be a BigInt.")
        }
    }
}

class Matrix {
    constructor(columns, ...indices) {
        let m = [];
        let r = [];
        if (indices[0] == "indenity") {
            for (let i = 0; i < columns; i++) {
                for (let j = 0; j < columns; j++) {
                    if (i == j) {
                        r.push(new Rational(1n));
                    } else {
                        r.push(new Rational(0n));
                    }
                }
                m.push(r);
                r = [];
            }
            this.indices = m;
            this.columns = columns;
            this.rows = columns;
        } else {
            let indicesLength = indices.length;
            if (indicesLength % columns == 0) {
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
        let A = this.indices[a];
        let B = this.indices[b];
        this.indices[a] = B;
        this.indices[b] = A;
    }
    addRow(row, dest, scale = new Rational(1n)) {
        for (let i = 0; i < this.columns; i++) {
            let r = this.indices[row][i].clone();
            console.log(r);
            r.mult(scale);
            console.log(r);
            this.indices[dest][i].add(r);
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
    }
    transpose() {
        let t = [];
        for (let i = 0; i < this.columns; i++) {
            for (let j = 0; j < this.rows; j++) {
                t.push(this.indices[j][i]);
            }
        }
        return new Matrix(this.rows, ...t);
    }
    clone() {
        let m = [];
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                m.push(this.indices[i][j].clone());
            }
        }
        return new Matrix(this.columns, ...m);
    }
    isSquare() {
        return this.columns == this.rows;
    }
    gaussianElimination() {
        let h = 0;
        let k = 0;
        while ((h < this.rows) && (k < this.columns)) {
            let iMax = h;
            let max = new Rational(0n);
            for (let i = h; i < this.rows; i++) {
                let c = this.indices[i][k].clone();
                c.numerator = MathJS.abs(c.numerator);
                c.sub(max)
                if (c.numerator > 0n) {
                    max = this.indices[i][k].clone();
                    max.numerator = MathJS.abs(max.numerator);
                    iMax = i;
                }
            }
            if (max.numerator == 0n) {
                k += 1;
            } else {
                this.swapRows(h, iMax);
                for (let i = h + 1; i < this.rows; i++) {
                    let f = this.indices[i][k].clone();
                    f.div(this.indices[h][k]);
                    this.indices[i][k] = new Rational(0n);
                    for (let j = k + 1; j < this.columns; j++) {
                        let s = this.indices[h][j].clone();
                        s.mult(f);
                        this.indices[i][j].sub(s);
                    }
                }
                h += 1;
                k += 1;
            }
        }
    }
    determinate() {
        if (this.isSquare()) {
            let M = this.clone();
            let neg = false;
            let h = 0;
            let k = 0;
            while ((h < M.rows) && (k < M.columns)) {
                let iMax = h;
                let max = new Rational(0n);
                for (let i = h; i < M.rows; i++) {
                    let c = M.indices[i][k].clone();
                    c.numerator = MathJS.abs(c.numerator);
                    c.sub(max)
                    if (c.numerator > 0n) {
                        max = M.indices[i][k].clone();
                        max.numerator = MathJS.abs(max.numerator);
                        iMax = i;
                    }
                }
                if (max.numerator == 0n) {
                    k += 1;
                } else {
                    M.swapRows(h, iMax);
                    if (h != iMax) {
                        neg = !neg;
                    }
                    for (let i = h + 1; i < M.rows; i++) {
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
            let p = new Rational((neg ? -1n : 1n));
            for (let i = 0; i < M.columns; i++) {
                p.mult(M.indices[i][i]);
            }
            return p;
        } else {
            return new Error("Can't take the determinate of a non-square matrix.");
        }
    }
    inverse() {
        if (this.isSquare()) {
            if (this.determinate().numerator != 0n) {
                let M = this.clone();
                let hMax = M.rows;
                let kMax = M.columns;
                let I = new Matrix(M.rows, "indenity");
                M = M.augment(I);
                let h = 0;
                let k = 0;
                while ((h < hMax) && (k < kMax)) {
                    let iMax = h;
                    let max = new Rational(0n);
                    for (let i = h; i < M.rows; i++) {
                        let c = M.indices[i][k].clone();
                        c.numerator = MathJS.abs(c.numerator);
                        c.sub(max)
                        if (c.numerator > 0n) {
                            max = M.indices[i][k].clone();
                            max.numerator = MathJS.abs(max.numerator);
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
                return new Matrix(kMax, ...m);
            } else {
                return new Error("Matrix has a determinate of zero.")
            }
        } else {
            return new Error("Matrix isn't square.")
        }
    }
    augment(B) {
        let m = [];
        if (this.rows == B.rows) {
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    m.push(this.indices[i][j].clone());
                }
                for (let j = 0; j < B.columns; j++) {
                    m.push(B.indices[i][j].clone());
                }
            }
            return new Matrix(this.columns + B.columns, ...m);
        } else {
            return new Error("Unable to augment due too inequal rows.");
        }
    }
    dotProduct(B) {
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
        let m = [];
        if (this.columns == B.columns && this.rows == B.rows) {
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    let p = this.indices[i][j].clone();
                    p.mult(B.indices[i][j]);
                    m.push(p);
                }
            }
            return new Matrix(this.columns, ...m);
        } else {
            return new Error("Can only take Hadamard product with matrices of the same size.");
        }
    }
    kroneckerProduct(B) {
        let m = [];
        for (let aI = 0; aI < this.rows; aI++) {
            for (let bI = 0; bI < B.rows; bI++) {
                for (let aJ = 0; aJ < this.columns; aJ++) {
                    for (let bJ = 0; bJ < B.columns; bJ++) {
                        let p = this.indices[aI][aJ].clone();
                        p.mult(B.indices[bI][bJ]);
                        m.push(p);
                    }
                }
            }
        }
        return new Matrix(this.columns * B.columns, ...m)
    }
    outerProduct(B) {
        if (this.columns == 1 && B.columns == 1) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < B.rows; j++) {
                    let p = this.indices[i][0].clone();
                    p.mult(B.indices[j][0]);
                    m.push(p);
                }
            }
            return new Matrix(B.rows, ...m);
        } else {
            return new Error("Can only take outer product with vectors. (A matrix with 1 column)");
        }
    }
    product(B) {
        if (this.columns == B.rows) {
            let m = [];
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < B.columns; j++) {
                    let sum = new Rational(0n);
                    for (let k = 0; k < this.columns; k++) {
                        let p = this.indices[i][k].clone();
                        p.mult(B.indices[k][j]);
                        sum.add(p);
                    }
                    m.push(sum);
                }
            }
            return new Matrix(B.columns, ...m);
        } else {
            return new Error("Matrices have incompatible sizes.");
        }
    }
}

const MathJS = {
    gcd: function (a, b) {
        a = MathJS.abs(a);
        b = MathJS.abs(b);
        while (b > 0) {
            if (a < b) {
                let c = b;
                b = a;
                a = c;
            } else {
                a %= b;
            }
        }
        return a;
    },
    abs: function (a) {
        return (a < 0) ? -a : a;
    },
    defactor: function (a, b) {
        while (a % b == 0) {
            a /= b;
        }
        return a;
    },
    bigFactorial: function (a) {
        let F = 1n;
        for (let i = 1n; i <= a; i++) {
            F *= i;
        }
        return F;
    },
    constants: {
        pi: new Rational(22n, 7n),
        e: new Rational(2n),
        calculateConstants: function (terms = 50) {
            let piSum = new Rational(3n);
            let eSum = new Rational(0n);
            for (let i = 0; i < terms; i++) {
                let sign = (i % 2 == 0 ? 1n : -1n);
                let bigI = BigInt(i);
                let piTerm = new Rational(4n * sign, (bigI + 2n) * (bigI + 3n) * (bigI + 4n));
                let eTerm = new Rational(1n, MathJS.bigFactorial(bigI));
                piSum.add(piTerm);
                eSum.add(eTerm);
            }
            this.pi = piSum;
            this.e = eSum;
        }
    }
}

MathJS.constants.calculateConstants();