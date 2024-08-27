/*
More Math library by Erik Haag version 2.0
https://github.com/ErikHaag/More-Math/
Dependencies:
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
        return new Rational(this.numerator, this.denominator);
    }
    simplify() {
        let factor = BigMathJS.gcd(this.numerator, this.denominator);
        this.numerator /= factor;
        this.denominator /= factor;
        if (this.denominator < 0) {
            this.numerator *= -1n;
            this.denominator *= -1n;
        }
    }
    toString(decimalLength = 3n, base = 10n, decimalSeperator = ".") {
        let int = this.numerator / this.denominator;
        let frac = this.clone();
        frac.mult(new Rational(BigMathJS.sign(int)));
        let baseRational = new Rational(base);
        frac.sub(new Rational(int));
        let remainders = [];
        if (decimalLength < 0n) {
            let newRemainder = true;
            while (newRemainder) {
                remainders.push({n: frac.numerator, d: frac.denominator});
                frac.mult(baseRational);
                let fInt = frac.floor();
                frac.sub(fInt);
                for (let r of remainders) {
                    if (r.n == frac.numerator && r.d == frac.denominator) {
                        newRemainder = false;
                        break;
                    }
                }
            }
        }
    }
    floor() {
        return new Rational(this.numerator / this.denominator);
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
    compare(B) {
        let difference = this.numerator * B.denominator - B.numerator * this.denominator
        if (difference > 0) {
            return 1;
        }else if (difference < 0) {
            return -1;
        }else {
            return 0;
        }
    }
}