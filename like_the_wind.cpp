
#include <iostream>
#include<cmath>// round需要，如果嫌弃用，自己定义也行
#include<string>
#include<vector>
#include<algorithm>//shuffle需要,如果非要省时间，同时去了也行
//#include<random>
#include<ctime>
using namespace std;

const double Pi = acos(-1.0);
const int MAXN = 3e3 + 10;

class EI;
class Complex
{
    double x, y;

public:
    Complex(double xx = 0, double yy = 0) { x = xx, y = yy; }//省着写空构造了
    double get_x() const { return x; }
    double get_y() const { return y; }
    void set_x(const double& x_) { x = x_; }
    void set_y(const double& y_) { y = y_; }
    Complex operator + (const Complex& a) const { return Complex(a.get_x() + x, a.get_y() + y); }
    Complex operator - (const Complex& a) const { return Complex(x - a.get_x(), y - a.get_y()); }
    Complex operator * (const Complex& a) const { return Complex(a.get_x() * x - a.get_y() * y, a.get_x() * y + a.get_y() * x); }
    Complex interge_division(const int& p)const { return Complex(x / p, y / p); }
    bool operator == (const Complex& a) const { return a.get_x() == x && a.get_y() == y; }
    EI to_EI();

}C[MAXN], D[MAXN];
int r[MAXN];
class EI
{
    int a;
    int b;

public:
    EI() : a(0), b(0) {};
    EI(int a_, int b_) : a(a_), b(b_) {};
    EI(int a_) : a(a_), b(0) {};
    const int get_a() const { return a; }
    const int get_b() const { return b; }
    const int get_norm() const { return a * a + b * b - a * b; }
    bool operator == (const EI& alpha)const { return a == alpha.get_a() && b == alpha.get_b(); }
    bool operator != (const EI& alpha)const { return a != alpha.get_a() || b != alpha.get_b(); }
    EI operator + (const EI& alpha)const { return EI(a + alpha.a, b + alpha.b); }
    EI operator - (const EI& alpha)const { return EI(a - alpha.a, b - alpha.b); }
    EI operator * (const EI& alpha)const { return EI(a * alpha.a - b * alpha.b, a * alpha.b + b * alpha.a - b * alpha.b); }
    EI operator * (int p)const { return EI(a * p, b * p); }  //又想了，这个好像不需要了，全程没有实数，
    EI operator / (const EI& alpha)const { pair<EI, EI> results = quotient_remainder(alpha); return results.first; }
    EI operator % (const EI& alpha)const { pair<EI, EI> results = quotient_remainder(alpha); return results.second; }
    friend ostream& operator << (ostream& out, const EI& alpha) {
        string s = alpha.stringEI();
        return out << s;
    }
    Complex to_complex();

    int myRound(const int& c, const int& d) const {
        int c_ = c % d;
        if (c_ > (d >> 1)) {
            c_ -= d;
        }
        else if (c_ <= -(d >> 1)) {
            c_ += d;
        }
        return(c - c_) / d;
    }

    pair<EI, EI> quotient_remainder(const EI& beta)const {  // 移位运算符的优先级比+-还低，不加括号会错
        int kers1 = (beta.get_a() << 1) - beta.get_b();
        int kers2 = (beta.get_b() << 1) - beta.get_a();
        int Q = beta.get_norm();
        int d = Q << 1;
        int s = a * kers1 + b * kers2;
        int t = b * beta.get_a() - a * beta.get_b();
        int x0 = myRound(s, d);
        int x1 = myRound(t, d);
        EI r1 = EI(x0 + x1, x1 << 1);
        EI beta1 = *this - beta * r1;
        int s_ = s + Q;
        int t_ = t - Q;
        int y0 = myRound(s_, d);
        int y1 = myRound(t_, d);
        EI r2 = EI(y0 + y1, (y1 << 1) + 1);
        EI beta2 = *this - beta * r2;
        pair<EI, EI> result;
        if (beta1.get_norm() < beta2.get_norm())
            result = make_pair(r1, beta1);
        else if (beta1.get_norm() > beta2.get_norm())
            result = make_pair(r2, beta2);
        else if (x0 < y0)
            result = make_pair(r1, beta1);
        else
            result = make_pair(r2, beta2);
        return result;
    }

    pair<bool, EI> get_inverse(EI alpha) {//*this对alpha的逆，如果不可逆，mod alpha为最大工因子时this对应的乘积式
        if (*this == EI(1))
            return make_pair(true, *this);
        else if (*this == EI(-1))
            return make_pair(true, -1);
        else if (*this == EI(1, 1))
            return make_pair(true, EI(0, -1));
        else if (*this == EI(-1, -1))
            return make_pair(true, EI(0, 1));
        else if (*this == EI(0, -1))
            return make_pair(true, EI(1, 1));
        else if (*this == EI(0, 1))
            return make_pair(true, EI(-1, -1));
        EI t, newt, A, newA, temp, quotient;
        t = EI(0), newt = EI(1);
        A = alpha, newA = *this;
        while (newA != EI(0)) {
            pair<EI, EI> qAndR = A.quotient_remainder(newA);
            temp = newt;
            newt = t - qAndR.first * newt;
            t = temp;
            A = newA;
            newA = qAndR.second;
        }
        if (A == EI(1))
            return make_pair(true, t);
        else if (A == EI(-1))
            return make_pair(true, EI(0) - t);
        else if (A == EI(1, 1))
            return make_pair(true, EI(0, -1) * t);
        else if (A == EI(-1, -1))
            return make_pair(true, EI(0, 1) * t);
        else if (A == EI(0, -1))
            return make_pair(true, EI(1, 1) * t);
        else if (A == EI(0, 1))
            return make_pair(true, EI(-1, -1) * t);
        else
            return make_pair(false, A);
    }

    string const stringEI()const {
        string result = "";
        if (a != 0) {
            result += to_string(a);
            if (b == 1)
                result += "+w";
            else if (b == -1)
                result += "-w";
            else if (b > 0)
                result += "+" + to_string(b) + "w";
            else if (b < 0)
                result += to_string(b) + "w";
        }
        else {
            if (b == 0)
                result = "0";  // 留着空也行，输出多项式的时候再看要不要输出
            else if (b == 1)
                result = "w";
            else if (b == -1)
                result = "-w";
            else
                result = to_string(b) + "w";
        }
        return result;
    }
};

EI Complex::to_EI() {  // 因为这俩会互相调用，所以先声明后定义了，别的不存在互相，所以直接定义在结构体内了
    double a = x + y / sqrt(3);
    double b = (y * 2) / sqrt(3);//我记得偏差不会特别大，好像可以去
    int a0 = (a > 0) ? (int)(a + 0.5) : (int)(a - 0.5);
    int b0 = (b > 0) ? (int)(b + 0.5) : (int)(b - 0.5);
    return EI(a0, b0);
}

Complex EI::to_complex() {
    return Complex(a - b / 2.0, sqrt(3) * b / 2.0);
}


class Poly
{
    vector<EI> poly;//vector的==靠的是定义的结构体，即因为重载了EI的==，这里poly是可以直接判断相等的
    static int N;
public:
    bool operator == (const Poly& B)const { return poly == B.poly; }
    friend ostream& operator << (ostream& out, const Poly& A) {
        string s = A.stringPoly();
        return out << s;
    }
    vector<EI> get_poly()const { return poly; }
    void set_N(const int& N_) { N = N_; }
    void push_poly(const EI& ele) { poly.push_back(ele); }
    void push_poly(const vector<EI>& v) { poly = v; }
    void delete_front() { poly.erase(poly.begin()); }
    void shuffle() { random_shuffle(poly.begin(), poly.end()); }

    string stringPoly()const {
        string polynomials, temp;
        bool first = true;
        for (int i = 0; i < poly.size(); i++) {
            if (poly[i] == EI(0))
                continue;
            int index = poly.size() - i - 1;
            if (index == 0) {
                if (poly[index].get_a() < 0 || poly[index].get_a() == 0 && poly[index].get_b() < 0)
                    temp = "(" + poly[i].stringEI() + ")";
                else
                    temp = poly[i].stringEI();
            }
            else if (index == 1) {
                if (poly[i] == EI(1))
                    temp = "x";
                else
                    temp = "(" + poly[i].stringEI() + ")x";
            }
            else {
                if (poly[i] == EI(1))
                    temp = "x^" + to_string(index);
                else
                    temp = "(" + poly[i].stringEI() + ")x^" + to_string(index);
            }
            if (first) {
                polynomials = temp;
                first = false;
            }
            else
                polynomials += "+" + temp;
        }
        return polynomials;
    }

    void FFT(Complex* A, int type, int limit) const //不能用vector，这赋值比较乱，且超vector，只能静态数组；计算部分不能去除前导0，否则不知哪一位
    {
        for (int i = 0; i < limit; i++)
            if (i < r[i]) swap(A[i], A[r[i]]); //结合刚才的r,实现了奇偶分块，就不用再申请空间存放系数了
        for (int mid = 1; mid < limit; mid <<= 1)
        {
            Complex Wn(cos(Pi / mid), type * sin(Pi / mid));  // 单位根必须用复数实现，如果能直接由爱森斯坦数实现就好了
            for (int R = mid << 1, j = 0; j < limit; j += R)
            {
                Complex w(1, 0);
                for (int k = 0; k < mid; k++, w = w * Wn)
                {
                    Complex x = A[j + k], y = w * A[j + mid + k];
                    A[j + k] = x + y;
                    A[j + mid + k] = x - y;
                }
            }
        }
    }

    Poly multiplication_FFT(const Poly& P, const EI& p) const {//FFT在主函数调用太多就会说长度太长，计算不了，因此目前只在Poly这个类的除法以及外层test_inverse_and_multiplications中调用。如果在ETRU类中调用FFT调用几次会计算算不动了。
        Poly result, A, B;
        A.poly = poly;
        B.poly = P.poly;
        while (A.poly[0] == EI(0)) {//去0
            A.poly.erase(A.poly.begin());
        }
        while (B.poly[0] == EI(0)) {
            B.poly.erase(B.poly.begin());
        }
        int len_A = A.poly.size();
        int len_B = B.poly.size();
        int len_standard = len_A + len_B;
        int  limit = 1, l = 0;

        while (limit < len_standard - 1) {
            limit <<= 1;
            l++;
        }
        memset(r, 0, sizeof(int) * limit);//静态数组需要用到的位置清零，因为len_A一会会填充，所以只把后面的填满了就好。
        for (int i = len_A; i < limit; i++) {
            C[i] = Complex(0);
        }
        for (int i = len_B; i < limit; i++) {
            D[i] = Complex(0);
        }
        for (int i = 0; i < len_A; i++) {
            C[i] = A.poly[i].to_complex();//多项式转为虚数静态数组表达
        }
        for (int i = 0; i < len_B; i++) {
            D[i] = B.poly[i].to_complex();
        }
        for (int i = 0; i < limit; i++) {
            r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));// 使得r从中间分成两部分，两部分中心对称元素对应的二进制码互补
        }
        FFT(C, 1, limit);//转FFT
        FFT(D, 1, limit);
        for (int i = 0; i <= limit; i++) C[i] = C[i] * D[i];
        FFT(C, -1, limit);//逆FFT
        int len = len_standard - 1;
        for (int i = 0; i < len; i++) {  // 转换爱森斯坦系数并并取模
            result.poly.push_back(C[i].interge_division(limit).to_EI() % p);
        }
        if (result.poly.size() >= N) {//截断
            int len = result.poly.size();
            for (int i = N; i < len; i++) {
                result.poly[i] = result.poly[i] + result.poly[i - N];
                result.poly[i] = result.poly[i] % p;
            }
            for (int i = N; i < len; i++)
                result.poly.erase(result.poly.begin());
        }
        return result;
    }

    Poly multiplication_with_EI(const EI& b, const EI& p)const {
        Poly result;
        for (EI a : poly)
            result.poly.push_back(a * b % p);
        return result;
    }

    Poly multiplication_normal(const Poly& P, const EI& p)const {
        Poly result, A, B;
        EI temp;
        A.poly = poly;
        B.poly = P.poly;
        Poly product;
        for (int i = 0; i < B.poly.size(); i++) {
            product = A.multiplication_with_EI(B.poly[i], p);
            for (int j = 0; j < (B.poly.size() - 1 - i); j++)
                product.poly.push_back(EI(0));
            result = result.addition(product, p);
            if (result.poly.size() >= N) {
                int len = result.poly.size();
                for (int i = N; i < len; i++) {
                    result.poly[i] = result.poly[i] + result.poly[i - N];
                    result.poly[i] = result.poly[i] % p;
                }
                for (int i = N; i < len; i++)
                    result.poly.erase(result.poly.begin());
            }
        }
        return result;
    }

    Poly multiplication_convolution(const Poly& P, const EI& p)const {
        Poly result, A, B;
        EI temp;
        A.poly = poly;
        B.poly = P.poly;
        int len_A = A.poly.size(), len_B = B.poly.size();
        int i, k;
        for (k = N - 1; k >= 0; k--) {
            temp = EI(0);
            for (i = 1; i < N - k; i++) {
                if (len_A > (k + i) && len_B > (N - i))
                    temp = temp + A.poly[len_A - 1 - k - i] * B.poly[len_B - 1 - N + i];
            }
            for (i = 0; i < k + 1; i++) {
                if (len_A > k - i && len_B > i)
                    temp = temp + A.poly[len_A - 1 - k + i] * B.poly[len_B - 1 - i];
            }
            result.poly.push_back(temp % p);
        }
        return result;
    }

    Poly subtraction(const Poly& P, const EI& p)const {
        Poly result;
        Poly A, B;
        A.poly = poly;
        B.poly = P.poly;
        if (A.poly.size() < B.poly.size()) {
            int len = B.poly.size() - A.poly.size();
            for (int i = 0; i < len; i++)// 头插和反转之后尾插再反转谁的复杂度低
                A.poly.insert(A.poly.begin(), EI(0));
        }
        else if (A.poly.size() > B.poly.size()) {
            int len = A.poly.size() - B.poly.size();
            for (int i = 0; i < len; i++)
                B.poly.insert(B.poly.begin(), EI(0));
        }
        for (int i = 0; i < B.poly.size(); i++)
            result.poly.push_back((A.poly[i] - B.poly[i]) % p);
        return result;
    }

    Poly addition(const Poly& P, const EI& p)const {//循环条件在跟着改变
        Poly result;
        Poly A, B;
        A.poly = poly;
        B.poly = P.poly;
        if (A.poly.size() < B.poly.size()) {
            int len = B.poly.size() - A.poly.size();
            for (int i = 0; i < len; i++)
                A.poly.insert(A.poly.begin(), EI(0));
        }
        else if (A.poly.size() > B.poly.size()) {
            int len = A.poly.size() - B.poly.size();
            for (int i = 0; i < len; i++)
                B.poly.insert(B.poly.begin(), EI(0));
        }
        for (int i = 0; i < A.poly.size(); i++)
            result.poly.push_back((A.poly[i] + B.poly[i]) % p);
        return result;
    }

    pair<Poly, Poly> division(const Poly& P, const EI& p)const {
        Poly A, B;
        A.poly = poly;
        B.poly = P.poly;
        while (A.poly[0] == EI(0)) {
            A.poly.erase(A.poly.begin());
        }
        while (B.poly[0] == EI(0)) {
            B.poly.erase(B.poly.begin());
        }
        Poly quotients;
        int initial_len = A.poly.size() - B.poly.size() + 1;
        quotients.poly.assign(initial_len, EI(0));
        for (int i = 0; i < quotients.poly.size(); i++) {
            if (A.poly.size() >= B.poly.size()) {
                EI quotient;
                pair<bool, EI> inverse = B.poly.front().get_inverse(p);
                if (inverse.first)
                    quotient = (inverse.second * A.poly.front()) % p;
                else {
                    cout << "不互素，重新开始" << endl;
                    Poly zero;
                    zero.poly.push_back(EI(0));
                    return make_pair(zero, zero);//作为不成功的标志
                }
                int new_len = A.poly.size() - B.poly.size() + 1;  //动态
                int index = initial_len - new_len;
                quotients.poly.at(index) = quotient;   // vector 没有-1这种访问方式
                Poly temp;
                temp.poly = B.poly;
                for (int i = 0; i < new_len - 1; i++)
                    temp.poly.push_back(EI(0));
                temp = temp.multiplication_with_EI(quotients.poly[i], p);
                A = A.subtraction(temp, p);
                int len = A.poly.size();
                for (int i = 0; i < len; i++) {
                    if (A.poly.front() == EI(0)) // poly.begin()是指针拉，poly[i]才是元素值，
                        A.poly.erase(A.poly.begin());
                    else
                        break;
                }
            }
            else
                break;
        }
        return make_pair(quotients, A);
    }

    pair<bool, Poly> polyInverse(Poly P, EI p) {
        Poly A, B, Zero;
        Zero.poly.push_back(EI(0));
        B.poly = poly;   //统一了,在interface说的时候，是小的对大的的逆，但是为了除法的时候少除一次，进来先交换
        A.poly = P.poly;  //现在A是大的
        if (B.poly.size() == 1) {  //预先判断，没有也行，但是判断如果是单位数，可以少进循环一次。
            pair<bool, EI> inverse = B.poly.front().get_inverse(p);
            if (inverse.first) {  //成功返回答案，不成功zero，虽然没zero也行，但给个吃的吧
                Poly result;
                result.push_poly(inverse.second);
                return make_pair(true, result);
            }
            else {
                return make_pair(false, Zero);
            }
        }
        Poly old, neww, newer;
        old.poly.push_back(EI(0));
        neww.poly.push_back(EI(1));
        pair<Poly, Poly> quotientAndRemainder;
        while (!B.poly.empty()) {
            quotientAndRemainder = A.division(B, p);
            if (quotientAndRemainder.first == Zero && quotientAndRemainder.second == Zero)
                return make_pair(false, Zero);
            newer = old.subtraction(quotientAndRemainder.first.multiplication_convolution(neww, p), p);
            A = B, B = quotientAndRemainder.second;
            old = neww, neww = newer;
        }
        if (A.poly.size() == 1) {
            pair<bool, EI> inverse = A.poly.front().get_inverse(p);
            if (inverse.first) {
                return make_pair(true, old.multiplication_with_EI(inverse.second, p));
            }
            else {
                return make_pair(false, Zero);
            }
        }
    }

    void translation_and_output()const {//写在Poly里就不用传参了
        vector<EI>input = this->poly;
        vector<int>binary;
        int len_ = input.size();
        for (int i = 0; i < len_ - 1; i++) {
            binary.push_back(input[i].get_a());
            binary.push_back(input[i].get_b());
        }
        if (input.back().get_b() == 1) {
            binary.push_back(input.back().get_a());  //好像这么做浪费了
        }
        len_ = binary.size();
        while (len_ % 7) {
            binary.insert(binary.begin(), 0);
        }
        len_ = binary.size();
        int len = len_ / 7;
        int k = 0, ASCII, a;
        for (int i = 0; i < len; i++) {
            ASCII = 0;
            a = 1;
            for (int j = 0; j < 7; j++) {
                ASCII += binary[k] * a;
                a *= 2;
                k++;//可以不定义变量k用i*7+j,但是要用一个除法，
            }
            cout << (char)ASCII;
        }
        cout << "————————————" << endl;
    }
};

int Poly::N = 3;

Poly translation_to_Poly(string s) {//觉得不是poly的行为，所以放外面了
    int len = s.size();
    vector<int>binary;
    Poly result;
    for (int i = 0; i < len; i++) {
        int ASCII = (int)s[i];
        int init_ASCII = ASCII;
        while (ASCII != 0) {
            binary.push_back(ASCII % 2);
            ASCII /= 2;
        }
        if (init_ASCII < 64)
            binary.push_back(0);
    }
    if (len % 2 == 1) {//标记位加填充位的做法，本来直接用-1标记就可以不用两位，但是好像解密会出错
        binary.push_back(1);
    }
    else {
        binary.push_back(1);
        binary.push_back(0);
    }
    int len_ = binary.size();
    for (int i = 0; i < len_ - 1; i += 2) {
        result.push_poly(EI(binary[i], binary[i + 1]));
    }
    return result;
}

class ETRU
{
    int N, d;//
    EI p, q;
    Poly random, plaintext, Rx;
public:
    ETRU() :N(41), p(EI(2)), q(), d(4) {};
    ETRU(int N_, EI p_, EI q_, int d_) {
        N = N_;
        p = p_;
        q = q_;
        d = d_;
        Rx.push_poly(EI(1));
        for (int i = 0; i < N_ - 1; i++)
            Rx.push_poly(EI(0));
        Rx.push_poly(EI(-1));
    }
    void set_random(const Poly& r) { random = r; }
    void set_plaintext(const Poly& m) { plaintext = m; }
    Poly get_random()const { return random; }
    Poly get_plaintext()const { return plaintext; }
    Poly sample_T(bool isF = false)const {
        int num;
        num = N - 6 * d;
        Poly L;
        L.set_N(N);  //poly的N初始化到前面也行
        for (int i = 0; i < d; i++) {
            L.push_poly(EI(1));
            L.push_poly(EI(-1));
            L.push_poly(EI(0, 1));
            L.push_poly(EI(0, -1));
            L.push_poly(EI(-1, -1));
            L.push_poly(EI(1, 1));
        }
        if (isF) {
            L.push_poly(EI(1));
            num--;
        }
        for (int i = 0; i < num; i++)
            L.push_poly(EI(0));
        L.shuffle();
        return L;
    }

    vector<Poly> KeyGen() {//标记interface是字符串加解密生成了不标准的poly，还是直接给定的标准poly
        Poly f = sample_T(true);
        Poly g = sample_T();
        pair<bool, Poly> inverse1 = f.polyInverse(Rx, q);
        pair<bool, Poly> inverse2 = f.polyInverse(Rx, p);
        while (!inverse1.first || !inverse2.first) {
            f = sample_T(true);
            inverse1 = f.polyInverse(Rx, q);
            inverse2 = f.polyInverse(Rx, p);
        }
        Poly Fq = inverse1.second;
        Poly Fp = inverse2.second;
        Poly h = Fq.multiplication_convolution(g, q);
        vector<Poly> result;
        /*while (h.get_poly()[0] == EI(0)) {//只看公钥大小的话，没必要在keygen里去除先导0，测试了一下大概会消耗500ms+，但是在d很大的时候公钥长度带来的缩小不明显，所以还是注释掉了
            h.delete_front();
        }
        while (Fp.get_poly()[0] == EI(0)) {
            Fp.delete_front();
        }
        while (Fq.get_poly()[0] == EI(0)) {
            Fq.delete_front();
        }*/
        result.push_back(h);
        result.push_back(f);
        result.push_back(Fp);
        return result;
    }

    Poly Encrypt(Poly pk, bool isPoly = false) {
        if (isPoly)
            random = sample_T();//如果是进行字符串的interface，由于转写的poly不满足Lf.Lg，因此为解密必须让随机random的长度小一些，否则解密失败。而如果是标准形式的poly，则对应的random也应该是标准的Lf,Lg
        Poly e = plaintext.addition(random.multiplication_convolution(pk, q).multiplication_with_EI(p, q), q);
        return e;
    }

    Poly Decrypt(Poly e, Poly sk, Poly Fp) {
        Poly a = sk.multiplication_convolution(e, q);
        Poly b = Fp.multiplication_convolution(a, p);
        int len = b.get_poly().size();
        while (b.get_poly()[0] == EI(0)) {
            b.delete_front();
        }
        return b;
    }

    bool Verify(Poly b) {
        while (plaintext.get_poly()[0] == EI(0)) {
            plaintext.delete_front();
        }
        return b == plaintext;
    }
};

void interface_with_data() {
    srand(unsigned(time(NULL)));//随机种子
    int N, a, b, d;
    EI p, q;
    string s;
    cout << "请输入加密信息（英文数字或符号，）" << endl << "能加密的信息长度与N的选择有很大关系，一个字符占用了7位，建议字符串长度<N/10" << endl;
    getline(cin, s);
    Poly m;
    m = translation_to_Poly(s);
    cout << endl << "请输入550以下的数，600百分之百会vector太长而报错，400的话概率较小报错" << endl << "示例：400" << endl;
    cin >> N;
    cout << "请输入d，建议d<N/11" << endl;
    cin >> d;
    cout << "请输入p（示例:当想要输入爱森斯坦整数2+3w时输入：2 3）" << endl;
    cin >> a >> b;
    p = EI(a, b);
    cout << "请输入q（示例:89+0w时输入：89 0）" << endl;
    cin >> a >> b;
    q = EI(a, b);
    ETRU demo = ETRU(N, p, q, d);
    demo.set_plaintext(m);
    vector<Poly> keys = demo.KeyGen();
    Poly h, f, Fp;
    h = keys[0];
    f = keys[1];
    Fp = keys[2];
    cout << "公钥对应的多项式是：" << h << endl;
    Poly cipertext = demo.Encrypt(h);
    cout << endl<<"明文是" << m << endl;
    cout << endl<<"密文是：" << cipertext << endl << endl;
    Poly afterDecryption = demo.Decrypt(cipertext, f, Fp);
    bool success = demo.Verify(afterDecryption);
    if (success) {
        cout << endl << "有时候translataion不输出，好像说是代码太长了定位不到函数，但是正确性可以通过自己比对明文密文多项式" << endl << "解密成功，字符串是" << endl;
        afterDecryption.translation_and_output();//好像嫌太长了，有时候这句不给输出
    }
    else
        cout << "解密失败";
}

void keyGen_only() {//无交互，只生产密钥，如果不可逆，求逆失败了ETRU中会说的
    srand(unsigned(time(NULL)));//随机种子
    int N, a, b, d;
    EI p, q;
    string s;
    cout << endl << "请输入N" << endl << "示例：400" << endl;
    cin >> N;
    cout << "请输入d，建议d<N/11" << endl;
    cin >> d;
    cout << "请输入p（示例:当想要输入爱森斯坦整数2+3w时输入：2 3）" << endl;
    cin >> a >> b;
    p = EI(a, b);
    cout << "请输入q（示例:89+0w时输入：89 0）" << endl;
    cin >> a >> b;
    q = EI(a, b);
    ETRU demo = ETRU(N, p, q, d);
    clock_t start, end;
    start = clock();
    vector<Poly> keys = demo.KeyGen();
    end = clock();
    cout << "————密钥产生完毕————" << endl;
    cout << "公钥对应的多项式是：" << keys[0] << endl;
    cout << "产生时间是：" << end - start << endl;
}


void interface_with_polynomial() {
    srand(unsigned(time(NULL)));//随机种子
    int N, a, b, d;
    EI p, q;
    string s;
    cout << endl << "请输入550以下的数，否则会vector太长而报错" << endl << "示例：400" << endl;
    cin >> N;
    cout << "请输入d，建议d<N/11" << endl;
    cin >> d;
    cout << "请输入p（示例:当想要输入爱森斯坦整数2+3w时输入：2 3）" << endl;
    cin >> a >> b;
    p = EI(a, b);
    cout << "请输入q（示例:89+0w时输入：89 0）" << endl;
    cin >> a >> b;
    q = EI(a, b);
    ETRU demo = ETRU(N, p, q, d);
    Poly m = demo.sample_T();
    demo.set_plaintext(m);
    vector<Poly> keys = demo.KeyGen();
    Poly h, f, Fp;
    h = keys[0];
    f = keys[1];
    Fp = keys[2];
    cout << "公钥对应的多项式是：" << h << endl;
    Poly cipertext = demo.Encrypt(h, true);
    cout << "密文是：" << cipertext << endl;
    cout << "原始明文是：" << m << endl << endl;
    Poly afterDecryption = demo.Decrypt(cipertext, f, Fp);
    cout << "解密的密文是：" << afterDecryption << endl;
    bool success = demo.Verify(afterDecryption);
    if (success)
        cout << "解密成功";
    else
        cout << "解密失败";
}

void test_inverse_and_multiplications() {//验证多项式逆以及三种乘法
    srand(unsigned(time(NULL)));
    int N, a, b, c, d;
    EI p, q;
    cout << "请输入N(乘法截断的方式资源消耗较大，N在100以上会报错，)" << endl << "(因此这里只是对比的话建议取50以下)（示例：41）";
    cin >> N;
    cout << "请输入d（建议<N/11）";
    cin >> d;
    cout << "请输入p（示例,当想要输入爱森斯坦整数2+3w时输入：2 3）";
    cin >> a >> b;
    p = EI(a, b);
    cout << "请输入q" << endl << "q要满足远大于p，在p为2 + 3w时，建议输入q的范数在80 + ，如简单的83，89" << endl << "示例，若输入89 + 0w, 请输入：89 0）";
    cin >> a >> b;
    q = EI(a, b);
    ETRU demo = ETRU(N, p, q, d);
    Poly sample = demo.sample_T(true);
    Poly sample2 = demo.sample_T();
    Poly Rx;  //用ETRU里面的也可以，但是还得调用，于是自己重新定义了
    Rx.push_poly(EI(1));
    for (int i = 0; i < N - 1; i++)
        Rx.push_poly(EI(0));
    Rx.push_poly(EI(-1));
    Poly inverse;
    inverse = sample.polyInverse(Rx, p).second;
    cout << "——————先进行多项式求逆的验证——————" << endl << "Rp是" << Rx << endl << "模数p是：" << p;
    cout << "多项式是 " << sample << endl << "逆是  " << inverse << endl << endl;
    cout << "乘积在环上取模验证的结果是： " << endl;
    cout << inverse.multiplication_convolution(sample, p) << endl;
    cout << "——————三种多项式乘法的验证——————" << endl << "多项式1是：" << sample << endl << "多项式2是" << sample2 << endl;
    cout << "截断乘法的结果是：" << endl << sample.multiplication_normal(sample2, p) << endl;
    cout << "卷积分两类计算的结果是：" << endl << sample.multiplication_convolution(sample2, p) << endl;
    cout << "卷积分两类计算的结果是：" << endl << sample.multiplication_FFT(sample2, p) << endl;
    cout << "—————————亲爱的同学，请你看看他们仨相等码？完结撒花———————————";
}

void test_translation() {//string转多项式，与多项式转回char的测试
    cout << "请输入字符串（英文、符号与数字，长度不限，不可输出“中文”和“空格”，因为用cin采集数据，以空格作为识别结束标志）" << endl;
    string s;
    cin >> s;
    Poly result = translation_to_Poly(s);
    result.translation_and_output();
}

void test_EI_to_Complex() {//爱森斯坦数到虚数再到爱森斯坦数的测试
    EI demo = EI(1115, 887);
    cout << demo << endl;
    Complex test = demo.to_complex();
    EI t = test.to_EI();
    cout << "最初的爱森斯坦数是：" << demo << "转换回来的爱森斯坦数是：" << t << endl;
}

int main() {
    interface_with_data();
}