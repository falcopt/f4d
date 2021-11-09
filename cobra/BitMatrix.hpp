#ifndef _F4D_BITMATRIX_HPP_
#define _F4D_BITMATRIX_HPP_

#include <cmath>
#include <vector>

#include "NonCopyable.hpp"
#include "SmallFlatSet.hpp"

namespace cobra {
    template <int maxSize>
    class BitMatrix : NonCopyable<BitMatrix<maxSize>> {

        std::vector<SmallFlatSet<unsigned short, static_cast<unsigned short>(~0), maxSize>> data;

    public:
        BitMatrix(int rows) : data(rows) { }

        inline auto reset(int row) -> void {
            data[row].clear();
        }

        inline auto set(int row, int entry) -> void {
            data[row].insert(entry);
        }

        inline auto is_set(int row, int entry) -> bool {
            return static_cast<bool>(data[row].count(entry));
        }

        inline auto overwrite(int source_row, int destination_row) {
            data[destination_row] = data[source_row];
        }

        inline auto& get_set_entries_possibly_with_duplicates(int row) {
            return data[row];
        }
    };

}  // namespace cobra


#endif
