#ifndef _F4D_SMALLFLATSET_HPP_
#define _F4D_SMALLFLATSET_HPP_

#include <algorithm>

#include "functor.hpp"

namespace cobra {

    template <typename Value, Value emptyValue, int maxSize, class op = identity_functor<Value>>
    class SmallFlatSet {

        static_assert(maxSize > 0, "Needs a positive value size.");
        static_assert(std::is_integral_v<decltype(maxSize)>, "Needs a integral type for maxSize");
        static_assert(maxSize * 5 / 4 <= (1 << 16), "Choose a smalle maxSize.");
        static constexpr int next2pow(int v) {
            v--;
            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v |= v >> 16;
            return ++v;
        }
        static_assert(next2pow(maxSize * 5 / 4) * sizeof(Value) <= (1 << 16), "Maximum memory occupation of 65KB.");

    public:
        class custom_iterator {
            friend class SmallFlatSet<Value, emptyValue, maxSize, op>;

        private:
            custom_iterator(Value* _base, Value* _end) : base(_base), end(_end) {
                while (base != end && *base == emptyValue) ++base;
            };

        public:
            inline auto& operator*() {
                return *base;
            }

            inline auto operator->() {
                return base;
            }

            inline auto& operator++() {
                do {
                    ++base;
                } while (base != end && *base == emptyValue);
                return *this;
            }

            inline auto operator!=(const custom_iterator x) {
                return x.base != base;
            }

            inline auto operator==(const custom_iterator x) {
                return x.base == base;
            }

        private:
            Value* base;
            const Value* end;
        };

    public:
        SmallFlatSet() {
            for (auto& p : buffer) p = emptyValue;
        };

        inline Value& find(const Value v) {
            auto index = op()(v) & realSizem1;
            auto value = buffer[index];
            while (value != v && value != emptyValue) {
                index = (index + 1) & realSizem1;
                value = buffer[index];
            }
            return buffer[index];
        }


        inline bool insert(const Value v) {
            // at least one element must remain empty (otherwise we need to check inside
            // "find")
            // if (sz >= maxSize)
            //    return false;

            auto& candidate_place = find(v);
            if (candidate_place != emptyValue) return false;  // element already there

            candidate_place = v;
            //++sz;
            return true;
        }

        inline bool insert_or_assign(const Value v) {
            // at least one element must remain empty (otherwise we need to check inside
            // "find")
            // if (sz >= realSize)
            //    return false;

            // if not there insert it, if there assign it.
            auto& candidate_place = find(v);
            // sz += static_cast<decltype(sz)>(candidate_place == emptyValue);
            candidate_place = v;

            return true;
        }

        inline auto& operator[](Value v) {
            auto& value = find(v);
            value = v;  // no check on sz, we need to decide if we want to check or not.
            return value;
        }

        inline void clear() {
            for (auto& p : buffer) p = emptyValue;
            // sz = 0;
        }

        inline size_t count(Value v) {
            return static_cast<size_t>(find(v) != emptyValue);
        }

        // inline auto size() const { return sz; };

        // inline auto empty() const { return sz == 0; };

        static inline auto get_emptykey() {
            return emptyValue;
        }

        static inline auto get_maxsize() {
            return maxSize;
        }

        inline auto begin() {
            return custom_iterator(buffer, buffer + realSize);
        };

        inline auto end() {
            return custom_iterator(buffer + realSize, buffer + realSize);
        };

    private:
        // Why 5/4 do you ask? Clearly a well thought value, not at all the nearest
        // small fraction that transform 25 to ~32.
        constexpr static int realSize = next2pow(maxSize * 5 / 4);
        constexpr static int realSizem1 = realSize - 1;

    public:
        // unsigned sz = 0u;
        Value buffer[realSize];
        // const HashOp op;
    };


    template <typename Value, Value emptyValue, int maxSize>
    class VerySmallFlatSet {

        static_assert(maxSize > 0, "Needs a positive value size.");
        static_assert(std::is_integral_v<decltype(maxSize)>, "Needs a integral type for maxSize");
        static_assert(maxSize * 5 / 4 <= (1 << 16), "Choose a smalle maxSize.");
        static constexpr int next2pow(int v) {
            v--;
            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v |= v >> 16;
            return ++v;
        }
        static_assert(next2pow(maxSize * 5 / 4) * sizeof(Value) <= (1 << 16), "Maximum memory occupation of 65KB.");

    public:
        class custom_iterator {
            friend class VerySmallFlatSet<Value, emptyValue, maxSize>;

        private:
            custom_iterator(Value* _base, Value* _end) : base(_base), end(_end){};

        public:
            inline auto& operator*() {
                return *base;
            }

            inline auto operator->() {
                return base;
            }

            inline auto& operator++() {
                ++base;
                return *this;
            }

            inline auto operator!=(const custom_iterator x) {
                return x.base != base;
            }

            inline auto operator==(const custom_iterator x) {
                return x.base == base;
            }

        private:
            Value* base;
            const Value* end;
        };

    public:
        VerySmallFlatSet() {
            for (auto& p : buffer) p = emptyValue;
        };

        inline Value& find(const Value v) {
            int index = 0;
            while (buffer[index] != v && buffer[index] != emptyValue) {
                ++index;
                assert(index < realSize - 1);
            }
            return buffer[index];
        }


        inline bool insert(const Value v) {
            auto& candidate_place = find(v);
            if (candidate_place != emptyValue) return false;  // element already there

            candidate_place = v;
            return true;
        }

        inline bool insert_or_assign(const Value v) {
            // if not there insert it, if there assign it.
            auto& candidate_place = find(v);
            candidate_place = v;

            return true;
        }

        inline auto& operator[](Value v) {
            // auto& value = find(v);
            // value = v;  // no check on sz, we need to decide if we want to check or not.
            return find(v) = v;
        }

        inline void clear() {
            for (auto& p : buffer) p = emptyValue;
        }

        inline size_t count(Value v) {
            return static_cast<size_t>(find(v) != emptyValue);
        }

        static inline auto get_emptykey() {
            return emptyValue;
        }

        static inline auto get_maxsize() {
            return maxSize;
        }

        inline auto begin() {
            return custom_iterator(buffer, buffer + realSize);
        };

        inline auto end() {
            return custom_iterator(buffer + realSize, buffer + realSize);
        };

    private:
        // Why 5/4 do you ask? Clearly a well thought value, not at all the nearest
        // small fraction that transform 25 to ~32.
        constexpr static int realSize = next2pow(maxSize * 5 / 4);
        constexpr static int realSizem1 = realSize - 1;

    public:
        Value buffer[realSize];
    };

}  // namespace cobra

#endif