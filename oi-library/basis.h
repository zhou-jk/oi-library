#include <vector>
struct basis {
private:
    static const int N = 63;
    long long a[N];
    std::vector<long long>b;
    int cnt;
public:
    basis() {
        cnt = 0;

        for (int i = 0; i < N; i++)
            a[i] = 0;

        b.clear();
        return;
    }
    void clear() {
        cnt = 0;

        for (int i = 0; i < N; i++)
            a[i] = 0;

        b.clear();
        return;
    }
    std::vector<long long>basis_set()const {
        std::vector<long long>c(b.size());

        for (int i = 0; i < (int)b.size(); i++)
            c[i] = a[b[i]];

        return c;
    }
    void insert(const long long &x) {
        cnt++;
        long long t = x;

        for (int i = N - 1; i >= 0; i--)
            if ((t >> i) & 1) {
                if (a[i])
                    t ^= a[i];
                else {
                    for (int j = 0; j < i; j++)
                        if ((t >> j) & 1)
                            t ^= a[j];

                    for (int j = i + 1; j < N; j++)
                        if ((a[j] >> i) & 1)
                            a[j] ^= t;

                    a[i] = t;
                    break;
                }
            }

        b.clear();

        for (int i = 0; i < N; i++)
            if (a[i])
                b.push_back(i);

        return;
    }
    basis insert(const basis &rhs) {
        std::vector<long long>c = rhs.basis_set();

        for (int i = 0; i < (int)c.size(); i++)
            insert(c[i]);

        return *this;
    }
    int size()const {
        return b.size();
    }
    long long max_xor()const {
        long long res = 0;

        for (int i = N - 1; i >= 0; i--)
            res ^= a[i];

        return res;
    }
    long long kth_xor(const long long &x)const {
        long long k = x;

        if (size() != cnt)
            k--;

        if (k > (1LL << size()) - 1)
            return -1;

        long long res = 0;

        for (int i = 0; i < (int)b.size(); i++)
            if ((k >> i) & 1)
                res ^= a[b[i]];

        return res;
    }
    long long rank(const long long &x)const {
        long long res = 0;

        for (int i = 0; i < (int)b.size(); i++)
            if ((x >> b[i]) & 1)
                res += 1LL << i;

        return res;
    }
    basis operator + (const basis &rhs)const {
        if (size() < rhs.size()) {
            basis c = rhs;
            std::vector<long long>d = basis_set();

            for (int i = 0; i < (int)d.size(); i++)
                c.insert(d[i]);

            return c;
        } else {
            basis c = *this;
            std::vector<long long>d = rhs.basis_set();

            for (int i = 0; i < (int)d.size(); i++)
                c.insert(d[i]);

            return c;
        }
    }
    basis operator += (const basis &rhs) {
        return *this = *this + rhs;
    }
};