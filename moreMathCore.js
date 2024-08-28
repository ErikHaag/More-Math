/*
More Math library by Erik Haag version 1.2
https://github.com/ErikHaag/More-Math/
Adds matrices and rationals to JavaScript, and some other things!
*/

const BigMathJS = {
    abs: function (a) {
        return (a < 0) ? -a : a;
    },
    defactor: function (a, b) {
        while (a % b == 0) {
            a /= b;
        }
        return a;
    },
    factorial: function (a) {
        let F = 1n;
        for (let i = 1n; i <= a; i++) {
            F *= i;
        }
        return F;
    },
    gcd: function (a, b) {
        a = BigMathJS.abs(a);
        b = BigMathJS.abs(b);
        while (b > 0) {
            if (a < b) {
                [a, b] = [b, a];
            } else {
                a %= b;
            }
        }
        return a;
    },
    lcm: function(a, b) {
        return (a * b) / BigMathJS.lcm(a, b);
    },
    mod: function(a, b) {
        let r = a % b;
        return r >= 0n ? r : r + b;
    },
    sign: function(a) {
        if (a > 0n) {
            return 1;
        } else if (a < 0n) {
            return -1;
        }
        return 0n;
    }
};