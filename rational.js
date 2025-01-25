/*
More Math library by Erik Haag version 2.1.0
https://github.com/ErikHaag/More-Math/
Dependencies: moreMathCore.js
*/

class Rational {
    constructor(numerator, denominator = 1n) {
        if (typeof numerator != "bigint" || typeof denominator != "bigint") {
            return new Error("Numerator and denominator must be BigInts");
        }
        if (numerator == 0n && denominator == 0n) {
            return new Error("Indeterminate form");
        }
        this.numerator = numerator;
        this.denominator = denominator;
        this.simplify();
    }

    add(B) {
        if (B instanceof Rational) {
            if (this.numerator == B.numerator && this.denominator == 0n && B.denominator == 0n) {
                //account for infinity + infinity or (-infinity) + (-infinity)
                return;
            }
            this.numerator = this.numerator * B.denominator + this.denominator * B.numerator;
            this.denominator *= B.denominator;
            let e = this.simplify();
            if (e instanceof Error) {
                return e;
            }
        } else if (typeof B == "bigint") {
            this.numerator += this.denominator * B;
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }

    ceiling() {
        if (this.denominator == 0n) return;
        //bring numerator up
        this.numerator += BigMathJS.mod(-this.numerator, this.denominator);
        this.integer();
    }

    clone() {
        return new Rational(this.numerator, this.denominator);
    }

    compare(B) {
        let difference = 0n;
        if (B instanceof Rational) {
            if (this.denominator == 0n && B.denominator == 0n) {
                //account for comparing infinities
                difference = this.numerator - B.numerator;
            } else {
                difference = this.numerator * B.denominator - B.numerator * this.denominator;
            }
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

    div(B) {
        if (B instanceof Rational) {
            this.numerator *= B.denominator;
            this.denominator *= B.numerator;
            let e = this.simplify();
            if (e instanceof Error) {
                return e;
            }
        } else if (typeof B == "bigint") {
            if (B == 0n) {
                if (this.numerator == 0n) {
                    return new Error("Indeterminate form");
                }
                this.numerator = this.numerator < 0n ? -1n : 1n;
                this.denominator = 0n;
            }
            if (this.numerator % B == 0n) {
                this.numerator /= B;
            } else {
                if (B < 0n) {
                    this.numerator *= -1n;
                    this.denominator *= -B;
                } else {
                    this.denominator *= B;
                }
            }
        }
    }

    floor() {
        if (this.denominator == 0n) return;
        //bring numerator down
        this.numerator -= BigMathJS.mod(this.numerator, this.denominator);
        this.integer();
    }

    integer() {
        if (this.denominator == 0n) return;
        //get integer part of corresponding decimal
        this.numerator /= this.denominator;
        this.denominator = 1n;
    }

    inverse() {
        [this.numerator, this.denominator] = [this.denominator, this.numerator];
    }

    mult(B) {
        if (B instanceof Rational) {
            this.numerator *= B.numerator;
            this.denominator *= B.denominator;
            let e = this.simplify();
            if (e instanceof Error) {
                return e;
            }
        } else if (typeof B == "bigint") {
            if (B == 0n && this.denominator == 0n) {
                return new Error("Indeterminate form");
            }
            if (this.denominator != 0n && this.denominator % B == 0n) {
                if (B < 0n) {
                    this.numerator *= -1;
                    this.denominator /= -B;
                } else {
                    this.denominator /= B;
                }
            } else {
                    this.numerator *= B;
            }
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }

    pow(B) {
        if (typeof B == "bigint") {
            if (B == 0n && this.denominator == 0n) {
                return new Error("Indeterminate form");
            }
            if (B < 0n) {
                [this.numerator, this.denominator] = [this.denominator, this.numerator];
                B *= -1n;
            }
            this.numerator **= B;
            this.denominator **= B;
        } else {
            return new Error("Argument must be a BigInt");
        }
    }

    simplify() {
        if (this.numerator == 0n && this.denominator == 0n) {
            return new Error("Indeterminate form");
        }
        let factor = BigMathJS.gcd(this.numerator, this.denominator);
        this.numerator /= factor;
        this.denominator /= factor;
        if (this.denominator < 0) {
            this.numerator *= -1n;
            this.denominator *= -1n;
        }
    }

    sub(B) {
        if (B instanceof Rational) {
            if (this.numerator == -B.numerator && this.denominator == 0n && B.denominator == 0n) {
                //account for infinity - (-infinity) or (-infinity) - infinity 
                return;
            }
            this.numerator = this.numerator * B.denominator - this.denominator * B.numerator;
            this.denominator *= B.denominator;
            let e = this.simplify();
            if (e instanceof Error) {
                return e;
            }
        } else if (typeof B == "bigint") {
            this.numerator -= this.denominator * B;
        } else {
            return new Error("Argument must be BigInt or Rational");
        }
    }

    toDecimal(decimalLength = 3n, base = 10n, decimalSeparator = ".") {
        if (this.denominator == 0n) {
            return (this.numerator < 0 ? "-" : "") + "Infinity";
        }
        //check if base is valid
        if (base < 2n || base > 36n) {
            return new Error("Invalid Base, must be between 2 and 36 (inclusive)");
        }
        //setup
        let negative = this.numerator < 0n;
        let int = BigMathJS.abs(this.numerator / this.denominator);
        let baseRational = new Rational(BigInt(base));
        let baseNumber = Number(base);
        let frac = this.clone();
        if (negative) {
            frac.mult(-1n);
        }
        frac.sub(int);
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
            return (negative ? "-" : "") + int + decimalSeparator + quotient.join("");
        }
    }

    toLatex() {
        if (this.denominator == 0n) {
            return (this.numerator < 0 ? "-" : "") + "\\infty";
        }
        if (this.denominator == 1n) {
            return this.numerator.toString();
        } else {
            return "\\frac{" + this.numerator + "}{" + this.denominator + "}"
        }
    }

    toString(hideDenominator = true) {
        if (this.denominator == 0) {
            return (this.numerator < 0 ? "-" : "") + "Infinity";
        } else if (hideDenominator && this.denominator == 1n) {
            return this.numerator.toString();
        } else {
            return this.numerator + "/" + this.denominator;
        }
    }
}