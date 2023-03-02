#include <complex>
#include <regex>
#include <cassert>
#include <ext/pb_ds/assoc_container.hpp>
#include <thread>
#include <random>
#include <iomanip>
#include <list>
#include <forward_list>
#include <chrono>
#include <utility>
#include <memory>
#include <valarray>
#include <cmath>
#include <queue>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <vector>
#include <bitset>
#include <numeric>
#include <functional>
#include <limits>
#include <set>
#include <unordered_set>
#include <type_traits>
#include <map>
#include <cerrno>

#define PI 22.0 / 7.0
#define none ""
#define __none []() {}
#define nl std::cout << std::endl
#define all(a) a.begin(), a.end()
#define inf std::numeric_limits<double>::infinity()
typedef unsigned long long int size;

//+---------------------------------I/O-----------------------------------+

#define in(a)          \
    ({                 \
        std::cin >> a; \
        a;             \
    })

void inm()
{
}
template <typename header, typename... trailer>
void inm(header &h, trailer &...t)
{
    in(h);
    return inm(t...);
}
// o(n)

#define out(a, separator)              \
    ({                                 \
        std::cout << (a) << separator; \
        "";                            \
    })

#define outc(__container, separator) \
    ({                               \
        for (auto &&i : __container) \
        {                            \
            out(i, separator);       \
        }                            \
    })

void outm()
{
}
template <typename header, typename... trailer>
void outm(header h, trailer... t)
{
    out(h, none);
    return outm(t...);
}

template <class T1, class T2>
std::istream &operator>>(std::istream &i, std::pair<T1, T2> &p)
{
    inm(p.first, p.second);
    return i;
}

template <class T1, class T2>
std::ostream &operator<<(std::ostream &o, const std::pair<T1, T2> &p)
{
    out(p.first, " ");
    out(p.second, none);
    return o;
}
template <class T>
std::ostream &operator<<(std::ostream &o,const std::complex<T>& c)
{
    o << c.real() << ' ';
    if (c.imag() > 0)
    {
        o << "+ " ;
        if (c.imag() != 1)
        {
            o << c.imag()<<"i";
        }
    }
    else
    {
        if (c.imag() < 0)
        {
            if(c.imag() != -1)
            {
                o << c.imag()<<"i";
            }
            else
            {
                o<<"-i";
            }
        }
    }
    return o;
}
//+---------------------------------ALIAS--------------------------------+
#define vi std::vector<int>
#define vll std::vector<long long>
#define vs std::vector<std::string>
#define pb push_back
#define pf push_front

namespace estd
{
    namespace imp
    {
        template <typename Container>
        void sort(Container &c);
        template <typename Container>
        using iterator_type = typename Container::iterator;

        template <typename Container>
        using value_type = typename Container::value_type;

        template <typename it>
        using iterator_category = typename std::iterator_traits<it>::iterator_category;

    }
    //+---------------------------------MATH----------------------------------+
    namespace math
    {
        double slope(std::pair<double, double> a, std::pair<double, double> b)
        {
            if ((a.first - b.first) == 0)
            {
                throw "the slope for vertical lines is undefined.\n";
            }
            return (a.second - b.second) / (a.first - b.first);
        }
        double fahrenheitToCelsius(double F)
        {
            return 5.0 * (F - 32.0) / 9.0;
        }
        double CelsiusToFahrenheit(double C)
        {
            return 9.0 * C / 5.0 + 32;
        }
        class RandomInt
        {
        public:
            RandomInt(int from, int to) : distro{from, to} {}
            int operator()()
            {
                return distro(en);
            }

        private:
            std::uniform_int_distribution<> distro;
            std::default_random_engine en;
        };
        // ما بقدر اعمله فنكشن اوبجكت لانه اليونيفورم دسترو هو ثابت وكل 1 , 2, 7 , 4 , 5... وكل مرة بستدعي الاو فن هيعطي اول فاليو الي هي فروم دائما

        class Line
        {
        private:
        public:
            Line(/* args */) {}
            ~Line() {}
        };
        /// @brief this function find an iterator of the smallest element in defined range.
        template <typename it, typename _Compare = std::greater<typename it::value_type>>
        it minc(it firstIt, it endIt, _Compare op = std::greater<typename it::value_type>())
        {
            it minIt = firstIt++;
            for (; firstIt != endIt; ++firstIt)
            {
                if (op(*minIt, *firstIt))
                {
                    minIt = firstIt;
                }
            }
            return minIt;
        }
        template<class T>
        T gcd(T a,T b)
        {
            !b ? throw "devision by zero\n" : nullptr;
            return a % b ? gcd(b, a % b) : b;
        }
        template <class U,class... T>
        U gcd(U a,U b ,T... c)
        {
            !b ? throw "devision by zero\n" : nullptr;
            return a % b ? gcd(b, a % b,c...) : gcd(b,c...);
        }

        long long lcm(long long a, long long b)
        {
            return (b) ? a * b / gcd(a, b) : throw "devision by zero\n";
        }

        template <class T>
        bool is_prime(T &&n)
        {
            static_assert(std::is_arithmetic<T>(), " is not an arithmetic type.");
            if (n == 2 | n == 3)
            {
                return true;
            }

            for (long long i{3}; i <= n / i; i += 2)
            {
                if (n % i == 0)
                {
                    return false;
                }
            }
            return true;
        }

        constexpr std::pair<double, double> midPoint(std::pair<double, double> a, std::pair<double, double> b)
        {
            return std::make_pair((a.first + b.first) / 2.0, (a.second + b.second) / 2.0);
        }

        constexpr double distance(std::pair<double, double> a, std::pair<double, double> b)
        {
            return std::abs(sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2)));
        }

        constexpr double circleArea(double radius)
        {
            return pow(radius, 2) * PI;
        }
        /// @brief this function solves linear equations.
        /// @param a x coefficient
        /// @param b constant
        /// @param c value after equality if exists
        /// @return value of x.
        constexpr double linear_equation(double &&a, double &&b, double &&c = 0)
        {
            if (!a)
            {
                throw " can not devision by zero.\n";
            }
            b >= 0 ? c -= b : c += -1 * b;
            return c / a;
        }

        /// @brief this function solves quadratic equations.
        /// @param a x^2 coefficient
        /// @param b x coefficient
        /// @param c constant
        /// @param d value after equality if exists
        /// @return value(s) of x if exists.
        constexpr std::pair<double, double> quadratic_equation(double a, double b, double c, double &&d = 0)
        {

            if (b * b - 4 * a * c < 0)
            {
                throw " there is no any real solution.\n";
            }
            if (!a)
            {
                throw " can not devision by zero.\n";
            }

            d >= 0 ? c -= d : c += -1 * d;
            double discriminant{std::sqrt(b * b - (4.0 * a * c))};
            auto p = std::make_pair((discriminant - b) / (2.0 * a), (-1.0 * discriminant - b) / (2.0 * a));
            return p;
        }

        constexpr double circularSectorArea(double radius, double radians)
        {
            return 0.5 * pow(radius, 2) * radians;
        }
        constexpr double radiansToDegrees(double radians)
        {
            return radians * 180.0 / PI;
        }

        constexpr double degreToRadians(double degree)
        {
            return degree * PI / 180.0;
        }

        /// @brief this function returns all primes number in [from : to]
        /// @param n end of interval
        /// @return vector of all primes number
        std::vector<int> eratosthenes_sieve(size &&from, size &&to)
        {
            std::vector<bool> is_prime(to - from + 1, true);
            vi prime;

            for (size i = 2; i <= to; ++i)
            {
                if (is_prime[i] && i >= from)
                {
                    prime.pb(i);
                    for (size j = 2; i * j <= to; ++j)
                    {
                        is_prime[j * i] = false;
                    }
                }
            }
            return prime;
        }
        class quadratic_function
        {
        private:
            double a, b, c;

        public:
            quadratic_function(double &&a, double &&b, double &&c) : a{a}, b{b}, c{c} {}
            double operator()(double x)
            {
                return a * x * x + b * x + c;
            }
            std::pair<double, double> solve()
            {
                return quadratic_equation(a, b, c);
            }
            double discriminant()
            {
                if (b * b - 4 * a * c < 0)
                {
                    throw " there is no any real solution.\n"; // change it or change complex version
                }
                if (!a)
                {
                    throw " can not devision by zero.\n";
                }

                return std::sqrt(b * b - (4.0 * a * c));
            }
            ~quadratic_function() {}
        };
    } // namespace math

    //+---------------------------------IMPLEMENTATION----------------------------------+
    namespace imp
    {
        template <class T1, class T2>
        std::istream &operator>>(std::istream &i, std::pair<T1, T2> &p)
        {
            inm(p.first, p.second);
            return i;
        }

        template <class T1, class T2>
        std::ostream &operator<<(std::ostream &o, const std::pair<T1, T2> p)
        {
            outm(p.first, p.second);
            return o;
        }

        template <typename T, typename... args>
        std::unique_ptr<T> make_unique(args &&..._Args)
        {
            return std::unique_ptr<T>{new T{std::forward<T>(_Args)...}};
        }

        template <typename T1, typename T2>
        std::pair<T1, T2> make_pair(T1 &&first, T2 &&second)
        {
            return std::pair<T1, T2>{first, second};
        }

        template <typename it>
        void sort_helper(it first, it end, std::random_access_iterator_tag)
        {
            std::sort(first, end);
        }

        template <typename it>
        void sort_helper(it first, it end, std::forward_iterator_tag)
        {
            std::vector<value_type<it>> v{first, end};
            std::sort(v.begin(), v.end());
            std::copy(v.begin(), v.end(), first);
        }

        template <typename Container>
        void sort(Container &c)
        {
            sort_helper(all(c), iterator_category<iterator_type<Container>>{}); // typename std::iterator_traits<typename Container::iterator>::iterator_category{}
        }

        template <typename Container>
        void inc_helper(Container &c, size_t size, std::random_access_iterator_tag)
        {
            value_type<Container> t;
            while (size--)
            {
                c.pb(in(t));
            }
        }

        template <typename Container>
        void inc_helper(Container &c, size_t size, std::bidirectional_iterator_tag)
        {
            value_type<Container> t;
            while (size--)
            {
                c.insert(in(t));
            }
        }

        template <typename Container>
        void inc_helper(Container &c, size_t size, std::forward_iterator_tag)
        {
            value_type<Container> t;
            while (size--)
            {
                c.pf(in(t));
            }
        }
        template <typename Container>
        void inc(Container &c, size_t size)
        {
            inc_helper(c, size, iterator_category<iterator_type<Container>>{});
        }
    } // namespace imp
    //+-----------------------------------------------------------------------ALGO-----------------------------------------------------------------------+
    namespace algo
    {
        /// @tparam _Compare function compare.
        /// @param first first element in the range.
        /// @param end last element in the range -not considered-.
        template <typename it, class _Compare = std::greater<typename it::value_type>>
        void selection_sort(it first, it end, _Compare op = std::greater<typename it::value_type>())
        {
            for (it i{first}; i != end; ++i)
            {
                std::swap(*i, *estd::math::minc(i, end, op));
            }
        }
        /// @tparam _Compare function compare.
        /// @param first first element in the range.
        /// @param end last element in the range -not considered-.
        template <typename it, class _Compare = std::greater<typename it::value_type>>
        void insertion_sort(it first, it end, _Compare op = std::greater<typename it::value_type>())
        {
            it i{first};
            ++i;
            for (; i != end; ++i)
            {
                typename it::value_type keep{*i};
                it j{i}, j1{i}; // j1: for bi-directional iterators.
                --j1;
                while ((j != first) && op(*j1, keep))
                {
                    *j = *j1;
                    --j;
                    --j1;
                }
                *j = keep;
            }
        }

        template <typename setOfPairs, class Operation>
        std::pair<std::vector<imp::iterator_type<setOfPairs>>, double> nearest_neighbor(setOfPairs &p, Operation op)
        {
            std::vector<bool> check(p.size(), true); // for checking of unvisited points
            std::vector<imp::iterator_type<setOfPairs>> pairsInOrder;

            int c = p.size() - 1;
            double sum{0.0};

            pairsInOrder.pb(p.begin());
            check[0] = false;

            op(*p.begin());
            while (c--)
            {
                imp::iterator_type<setOfPairs> temp;
                double min{std::numeric_limits<double>::max()};
                int j{0}, index{1};
                for (imp::iterator_type<setOfPairs> i{p.begin()}; i != p.end(); ++i)
                {
                    if (check[j] && math::distance(*i, *pairsInOrder.back()) < min)
                    {
                        temp = i;
                        index = j;
                        min = math::distance(*i, *pairsInOrder.back());
                    }
                    ++j;
                }
                op(*temp);
                pairsInOrder.pb(temp);
                sum += min;
                check[index] = false;
            }
            return std::make_pair(pairsInOrder, sum + math::distance(*pairsInOrder.front(), *pairsInOrder.back()));
        }
        template <typename setOfIntervals, class Operation>
        std::vector<imp::iterator_type<setOfIntervals>> scheduling(setOfIntervals &__setOfIntervals, Operation op)
        {
            std::vector<bool> check(__setOfIntervals.size(), true);
            std::vector<imp::iterator_type<setOfIntervals>> intervalsInOrder;
            imp::iterator_type<setOfIntervals> r{new imp::value_type<setOfIntervals>{
                std::numeric_limits<typename imp::value_type<setOfIntervals>::first_type>::max(),
                std::numeric_limits<typename imp::value_type<setOfIntervals>::second_type>::max()}};

            imp::iterator_type<setOfIntervals> temp{r};
            unsigned int before = intervalsInOrder.size(); // for while loop / checking after any updates done.
            int j{0}, index{j};

            for (auto i{__setOfIntervals.begin()}; i != __setOfIntervals.end(); ++i, ++j)
            {
                if (i->second < temp->second)
                {
                    temp = i;
                    index = j;
                }
            }
            check[index] = false;
            intervalsInOrder.pb(temp);
            j = 0;
            op(*temp);
            temp = r;
            while (before != intervalsInOrder.size())
            {
                before = intervalsInOrder.size();
                for (auto i{__setOfIntervals.begin()}; i != __setOfIntervals.end(); ++i, ++j)
                {

                    if (check[j] && i->first >= intervalsInOrder.back()->second) // check overlapping
                    {
                        if (i->second < temp->second) // separated for optimizing
                        {
                            temp = i;
                            index = j;
                        }
                    }
                }
                if (temp != r)
                {
                    op(*temp);
                    check[index] = false;
                    intervalsInOrder.pb(temp);
                    temp = r;
                    j = 0;
                }
            }
            return intervalsInOrder;
        }
    } // namespace algo

    //+---------------------------------OBJECT FUNCTION---------------------------------+

    template <class T>
    struct Greater
    {
        Greater() = default;
        bool
        operator()(const T &v1, const T &v2)
        {
            return v1 > v2;
        }
    };

    template <class T>
    constexpr bool is_arithmetic()
    {
        return std::is_arithmetic<T>::value;
    }

    template <class T>
    struct Natureral
    {
    public:
        static_assert(is_arithmetic<T>(), " is not an arithmetic type.");

        Natureral() = default;
        bool
        operator()(const T &v) const
        {
            return v > 0;
        }
    };
    template <class Container, class Opration>
    unsigned int count_if(Container &c, Opration o)
    {
        int _Count = 0;
        imp::iterator_type<Container> it{c.begin()};
        for (; it != c.end(); ++it)
        {
            if (o(*it))
            {
                ++_Count;
            }
        }
        return _Count;
    }

    template <typename Container>
    imp::value_type<Container> sum(Container &c)
    {
        imp::value_type<Container> v = 0;
        for (auto &&x : c)
        {
            v += x;
        }
        return v;
    }

    template <typename Container, typename Operation>
    std::vector<imp::iterator_type<Container>> find_all(Container &c, Operation op)
    {
        std::vector<imp::iterator_type<Container>> v;

        for (auto first{c.begin()}; first != c.end(); ++first)
        {
            if (op(*first))
            {
                v.push_back(first);
            }
        }
        return v;
    }
} // namespace estd

int Round(double num)
{
    return std::floor(static_cast<int>(num + 0.5)); // 5/4 style
}
int trancate(double num)
{
    return std::floor(static_cast<int>(num));
}

enum class round_style
{
    round,
    trancate
};
std::function<int(double)> R;
struct ROUND
{
    round_style r;
    ROUND(round_style R) : r(R) {}
    int operator()(double num) { return static_cast<int>(r == round_style::round ? num + 0.5 : num); }
};
struct F
{
    vi v;
    F(vi &v) : v{v} {}
    long long operator()()
    {
        return std::accumulate(all(v), 0ll);
    }
};
struct Train
{
    int key;
    bool check;
};

int main()
try
{
#ifndef ONLINE_JUDGE
    out("Fill text input :", '\n');
    system("cat > input.text");
    freopen("input.text", "r", stdin);
    auto start = std::chrono::high_resolution_clock::now();
#endif
    
    std::string x[] = {"ahmad"," ali"};
    out(x->length(),"\n");
#ifndef ONLINE_JUDGE
    nl;
    out(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count(), "ms\n");
    system("rm -f input.text");
#endif
} //(double(*)(int))

catch (const std::exception &e)
{
    std::cerr << e.what();
    nl;
}
catch (const char *e)
{
    std::cerr << e << '\n';
}
catch (...)
{
    std::cerr << "unknown exception is thrown.";
    nl;
}
/*
اللغة الالمانية(1)		12.30	14.00		
التشغيل			9.30	11.00		
أساسيات علم البيانات		14.00	15.30		
الاعجاز العلمي في القران الكريم		11.00	12.00		
هندسة البرمجيات		11.00	12.30		
*/