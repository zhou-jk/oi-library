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
#include <cstdio>
int main() {
    int n, q;
    scanf("%d %d", &n, &q);

    fenwick_tree<long long> fw(n);

    for (int i = 0; i < n; i++) {
        int a;
        scanf("%d", &a);
        fw.add(i, a);
    }

    for (int i = 0; i < q; i++) {
        int t;
        scanf("%d", &t);

        if (t == 0) {
            int p, x;
            scanf("%d %d", &p, &x);
            fw.add(p, x);
        } else {
            int l, r;
            scanf("%d %d", &l, &r);
            printf("%lld\n", fw.sum(l, r));
        }
    }
}