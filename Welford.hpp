#ifndef _F4D_WELFORD_HPP_
#define _F4D_WELFORD_HPP_

namespace cobra {

    // https://gist.github.com/alexalemi/2151722
    class Welford {

    public:
        Welford() = default;

        Welford(const Welford& other) {
            k = other.k;
            mean = other.mean;
        }

        Welford& operator=(const Welford& other) = default;

        void update(float x) {
            ++k;
            const auto new_mean = mean + (x - mean) * 1.0f / static_cast<float>(k);
            mean = new_mean;
        }

        float get_mean() const {
            return mean;
        }

        void reset() {
            k = 0;
            mean = 0.0f;
        }

    private:
        unsigned long k = 0;
        double mean = 0.0;
    };

}  // namespace cobra

#endif