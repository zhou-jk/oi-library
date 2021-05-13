#include <vector>
struct dsu {
private:
    int n;
    std::vector<int>fa, siz;
public:
    dsu(const int &_n) {
        n = _n;
        fa.resize(n);
        siz.resize(n);

        for (int i = 0; i < n; i++)
            fa[i] = i, siz[i] = 1;

        return;
    }
    int getf(const int &v) {
        if (v == fa[v])
            return v;
        else
            return fa[v] = getf(fa[v]);
    }
    int leader(const int &v) {
        return getf(v);
    }
    int merge(const int &u, const int &v) {
        int f1 = getf(u), f2 = getf(v);

        if (f1 != f2)
            fa[f2] = f1, siz[f1] += siz[f2];

        return f1;
    }
    bool same(const int &u, const int &v) {
        return getf(u) == getf(v);
    }
    int size(int v) {
        return siz[getf(v)];
    }
};