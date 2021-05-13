#include <cassert>
#include <vector>
#include "modint.h"
struct matrix {
private:
    int n, m;
    std::vector<std::vector<modint>>mat;
public:
    matrix(const int &_n = 0, const int &_m = 0) {
        n = _n, m = _m;
        mat.resize(n);

        for (int i = 0; i < n; i++)
            mat[i].resize(m);

        return;
    }
    void resize(const int &_n = 0, const int &_m = 0) {
        n = _n, m = _m;
        mat.resize(n);

        for (int i = 0; i < n; i++)
            mat[i].resize(m);

        return;
    }
    void clear() {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                mat[i][j] = 0;

        return;
    }
    modint det()const {
        if (n != m)
            return -1;

        std::vector<std::vector<modint>>C = mat;
        int cur = 0;
        modint res = 1;

        for (int i = 0; i < m; i++) {
            int pos = cur;

            for (int j = cur; j < n; j++)
                if (C[j][i].val() > C[pos][i].val())
                    pos = j;

            if (C[pos][i].val() == 0)
                return 0;

            if (cur != pos)
                swap(C[pos], C[cur]), res = -res;

            for (int j = cur + 1; j < n; j++)
                while (C[j][i].val()) {
                    if (C[j][i].val() < C[cur][i].val())
                        swap(C[cur], C[j]), res = -res;

                    int d = C[j][i].val() / C[cur][i].val();

                    for (int k = i; k < m; k++)
                        C[j][k] -= C[cur][k] * d;
                }

            res *= C[cur][i];
            cur++;
        }

        return res;
    }
    matrix operator + (const matrix &rhs)const {
        assert(n == rhs.n && m == rhs.m);
        matrix c(n, m);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                c.mat[i][j] = mat[i][j] + rhs.mat[i][j];

        return c;
    }
    matrix operator - (const matrix &rhs)const {
        assert(n == rhs.n && m == rhs.m);
        matrix c(n, m);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                c.mat[i][j] = mat[i][j] - rhs.mat[i][j];

        return c;
    }
    matrix operator * (const matrix &rhs)const {
        assert(m == rhs.n);
        int _n = n, _m = rhs.m;
        matrix c(_n, _m);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < rhs.m; j++)
                for (int k = 0; k < m; k++)
                    c.mat[i][j] += mat[i][k] * rhs.mat[k][j];

        return c;
    }
    matrix operator += (const matrix &rhs) {
        return *this = *this + rhs;
    }
    matrix operator -= (const matrix &rhs) {
        return *this = *this - rhs;
    }
    matrix operator *= (const matrix &rhs) {
        return *this = *this * rhs;
    }
    std::vector<modint> &operator [](const int &i) {
        return mat[i];
    }
};