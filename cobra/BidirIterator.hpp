#ifndef CAV_BIDIRITERATOR_HPP
#define CAV_BIDIRITERATOR_HPP

#include <iterator>
#include <numeric>

namespace cav {
    template <typename Base, class OperInc, class OperDec, class OperDefer>
    struct BidirIterator {
        typedef typename std::bidirectional_iterator_tag iterator_category;
        typedef typename std::ptrdiff_t difference_type;
        typedef typename std::remove_reference_t<typename std::result_of_t<OperDefer(Base&)>> value_type;
        typedef value_type* pointer;
        typedef typename std::result_of_t<OperDefer(Base&)> reference;

    public:
        BidirIterator(Base base_) : base(base_){};

        inline reference operator*() {
            return OperDefer()(base);
        }

        inline BidirIterator operator+(difference_type i) noexcept {
            Base _base = base;
            for (; i > 0; --i) _base = OperInc()(_base);
            return BidirIterator(_base);
        }

        inline BidirIterator operator++() noexcept {
            base = OperInc()(base);
            return *this;
        }

        inline BidirIterator operator+=(difference_type i) noexcept {
            for (; i > 0; --i) base = OperInc()(base);
            return *this;
        }

        inline BidirIterator operator++(int) noexcept {
            BidirIterator curr = base;
            base = OperInc()(base);
            return curr;
        }

        inline BidirIterator operator-(difference_type i) noexcept {
            Base _base = base;
            for (; i > 0; --i) _base = OperDec()(_base);
            return BidirIterator(_base);
        }

        inline BidirIterator operator--() noexcept {
            base = OperDec()(base);
            return *this;
        }

        inline BidirIterator operator-=(difference_type i) noexcept {
            for (; i > 0; --i) base = OperDec()(base);
            return *this;
        }

        inline BidirIterator operator--(int) noexcept {
            BidirIterator curr = base;
            base = OperDec()(base);
            return curr;
        }

        inline bool operator==(const BidirIterator& other) const noexcept {
            return base == other.base;
        }
        inline bool operator!=(const BidirIterator& other) const noexcept {
            return base != other.base;
        }

    private:
        Base base;
    };
}  // namespace cav

#endif