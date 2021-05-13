#include <vector>
template<typename T>
struct fenwick_tree {
private:
    int n;
    std::vector<T>C;
protected:
    T getsum(const int &x)const {
        T res = 0;

        for (int i = x + 1; i > 0; i -= (i & -i))
            res = res + C[i - 1];

        return res;
    }
public:
    fenwick_tree(const int &_n) {
        n = _n;
        C.resize(n);

        for (int i = 0; i < n; i++)
            C[i] = 0;

        return;
    }
    void add(const int &x, const T &y) {
        for (int i = x + 1; i <= n; i += (i & -i))
            C[i - 1] = C[i - 1] + y;

        return;
    }
    T sum(const int &l, const int &r)const {
        return getsum(r - 1) - getsum(l - 1);
    }
};