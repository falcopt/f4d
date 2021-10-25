#ifndef _F4D_FLAT2DVECTOR_HPP_
#define _F4D_FLAT2DVECTOR_HPP_

#include <vector>

namespace cobra {

    template <typename T>
    class Flat2DVector {

        std::vector<T> data;

        size_t rows;
        size_t cols;

    public:
        Flat2DVector() : rows(0), cols(0) { }

        inline void resize(size_t rows_, size_t cols_) {
            rows = rows_;
            cols = cols_;
            data.resize(rows * cols, 0);
        }

        inline void at(const int i, const int j, T value) {
            data[i * cols + j] = value;
        }

        inline T at(const int i, const int j) const {
            return data[i * cols + j];
        }
    };

}  // namespace cobra

#endif