#ifndef CAV_RANGE_HPP
#define CAV_RANGE_HPP

#include <iterator>

namespace cav {

    /**
     * @brief Just a pair of iterators for ranged for.
     *
     * @tparam IterT
     */
    template <typename IterT>
    class Range {

        static_assert(std::is_same_v<typename std::iterator_traits<IterT>::iterator_category, std::bidirectional_iterator_tag>,
                      "The iterator type need to be a bidirectional.");

    public:
        Range(IterT _first, IterT _last) : first(_first), last(_last){};

        inline auto begin() { return first; }
        inline auto begin() const { return first; }

        inline auto end() { return last; }
        inline auto end() const { return last; }

    private:
        IterT first;
        IterT last;
    };
}  // namespace cav

#endif