namespace mod_int_define {
int MOD;
void set_mod(const int &m) {
    MOD = m;
    return;
}
int mod() {
    return MOD;
}
int ex_gcd(const int &a, const int &b, int &x, int &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }

    int d = ex_gcd(b, a % b, x, y);
    int tmp = x;
    x = y, y = tmp - a / b * y;
    return d;
}
struct modint {
private:
    int v;
protected:

public:
    modint(const int &x = 0) {
        v = (x % MOD + MOD) % MOD;
        return;
    }
    modint(const long long &x) {
        v = (x % MOD + MOD) % MOD;
        return;
    }
    modint operator = (const int &x) {
        v = (x % MOD + MOD) % MOD;
        return *this;
    }
    modint operator = (const long long &x) {
        v = (x % MOD + MOD) % MOD;
        return *this;
    }
    int val()const {
        return v;
    }
    operator int() {
        return v;
    }
    modint operator + (const modint &rhs)const {
        modint c;
        c.v = (v + rhs.v);

        if (c.v >= MOD)
            c.v -= MOD;

        return c;
    }
    modint operator - (const modint &rhs)const {
        modint c;
        c.v = (v - rhs.v);

        if (c.v < 0)
            c.v += MOD;

        return c;
    }
    modint operator * (const modint &rhs)const {
        modint c;
        c.v = 1LL * v * rhs.v % MOD;
        return c;
    }
    modint getinv()const {
        int x, y;
        int d = ex_gcd(v, MOD, x, y);

        if (d != 1)
            return -1;

        return (x + MOD) % MOD;
    }
    modint operator / (const modint &rhs)const {
        modint c = *this * rhs.getinv();
        return c;
    }
    modint operator += (const modint &rhs) {
        return *this = *this + rhs;
    }
    modint operator -= (const modint &rhs) {
        return *this = *this - rhs;
    }
    modint operator *= (const modint &rhs) {
        return *this = *this * rhs;
    }
    modint operator /= (const modint &rhs) {
        return *this = *this / rhs;
    }
    modint operator == (const modint &rhs)const {
        return v == rhs.v;
    }
    modint operator != (const modint &rhs)const {
        return v != rhs.v;
    }
    modint operator ++ () {
        v++;

        if (v >= MOD)
            v -= MOD;

        return *this;
    }
    modint operator -- () {
        v--;

        if (v < 0)
            v += MOD;

        return *this;
    }
    modint operator - ()const {
        modint c;
        c.v = -v;

        if (c.v < 0)
            c.v += MOD;

        return c;
    }
};
modint pow(const modint &aa, const long long &bb) {
    modint a = aa;
    int b = bb;
    modint res = 1;

    while (b) {
        if (b & 1)
            res = res * a;

        a = a * a, b >>= 1;
    }

    return res;
}
}
using mod_int_define::set_mod;
using mod_int_define::mod;
using mod_int_define::modint;
using mod_int_define::pow;