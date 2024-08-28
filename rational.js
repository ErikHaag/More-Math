/*
More Math library by Erik Haag version 2.0
https://github.com/ErikHaag/More-Math/
Dependencies: moreMathCore.js
*/

class Rational {
    constructor(numerator, denominator = 1n) {
        if (typeof numerator == "bigint" && typeof denominator == "bigint") {
            if (0 < denominator) {
                this.numerator = numerator;
                this.denominator = denominator;
            } else {
                return new Error("Denomimator must be greater than 0");
            }
        } else {
            return new Error("Numerator and denominator must be BigInts");
        }
    }
    toLatex() {
        if (this.denominator == 1n) {
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
    toString(hideDenominator = true) {
        if (hideDenominator && this.denominator == 1n) {
            return this.numerator.toString();
        } else {
            return this.numerator + "/" + this.denominator;
        }
    }
    toDecimal(decimalLength = 3n, base = 10n, decimalSeparator = ".") {
        //check if base is valid
        if (base < 2n || base > 36n) {
            return new Error("Invalid Base, must be between 2 and 36 (inclusive)");
        }
        //setup
        let int = this.numerator / this.denominator;
        let baseRational = new Rational(BigInt(base));
        let baseNumber = Number(base);
        let frac = this.clone();
        frac.sub(int);
        if (int < 0n) {
            frac.mult(-1n);
        }
        int = int.toString(baseNumber);
        if (decimalLength == 0n) {
            return int;
        } else {
            //long division
            let digits = 0n;
            let repeatStart = -1n;
            let remainders = [];
            let quotient = [];
            outer: while (decimalLength < 0n || digits < decimalLength) {
                remainders.push({ n: frac.numerator, d: frac.denominator });
                frac.mult(baseRational);
                let fInt = frac.clone();
                fInt.integer();
                quotient.push(fInt.numerator);
                digits++;
                frac.sub(fInt);
                for (const i in remainders) {
                    const r = remainders[i];
                    if (r.n == frac.numerator && r.d == frac.denominator) {
                        repeatStart = BigInt(i);
                        break outer;
                    }
                }
            }
            //replace integers with characters
            quotient = quotient.map(element => element.toString(baseNumber));
            if (repeatStart >= 0n) {
                //insert repeating brackets
                quotient.splice(Number(repeatStart), 0, "[");
                quotient.push("]");
            }
            //package up the string
            return int + decimalSeparator + quotient.join("");
        }
    }
    floor() {
        //bring numerator down
        this.numerator -= BigMathJS.mod(this.numerator, this.denominator);
        this.integer();
    }
    ceiling() {
        //bring numerator up
        this.numerator += BigMathJS.mod(-this.numerator, this.denominator);
        this.integer();
    }
    integer() {
        //get integer part of corresponding decimal
        this.numerator /= this.denominator;
        this.denominator = 1n;
    }
    inverse() {
        if (this.numerator == 0n) {
            return new Error("Division by zero!");
        }
        [this.numerator, this.denominator] = [this.denominator, this.numerator];
    }
    add(B) {
        if (B instanceof Rational) {
            this.numerator = this.numerator * B.denominator + this.denominator * B.numerator;
            this.denominator *= B.denominator;
            this.simplify();
        } else if (typeof B == "bigint") {
            this.numerator += this.denominator * B;
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }
    sub(B) {
        if (B instanceof Rational) {
            this.numerator = this.numerator * B.denominator - this.denominator * B.numerator;
            this.denominator *= B.denominator;
            this.simplify();
        } else if (typeof B == "bigint") {
            this.numerator -= this.denominator * B;
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }
    mult(B) {
        if (B instanceof Rational) {
            this.numerator *= B.numerator;
            this.denominator *= B.denominator;
            this.simplify();
        } else if (typeof B == "bigint") {
            if (this.denominator % B == 0n) {
                this.denominator /= B;
            } else {
                this.numerator *= B;
            }
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }
    div(B) {
        if (B instanceof Rational) {
            if (B.numerator == 0n) {
                return new Error("Division by zero!");
            }
            this.numerator *= B.denominator;
            this.denominator *= B.numerator;
            this.simplify();
        } else if (typeof B == "bigint") {
            if (B == 0n) {
                return new Error("Division by zero!");
            }
            if (this.numerator % B == 0n) {
                this.numerator /= B;
            } else {
                this.denominator *= B;
            }
        }
    }
    pow(B) {
        if (typeof B == "bigint") {
            this.numerator **= B;
            this.denominator **= B;
        } else {
            return new Error("Argument must be a BigInt");
        }
    }
    compare(B) {
        let difference = 0n;
        if (B instanceof Rational) {
            difference = this.numerator * B.denominator - B.numerator * this.denominator;
        } else if (typeof B == "bigint") {
            difference = this.numerator - B * this.denominator
        }
        if (difference > 0n) {
            return 1n;
        } else if (difference < 0n) {
            return -1n;
        } else {
            return 0n;
        }
    }
}