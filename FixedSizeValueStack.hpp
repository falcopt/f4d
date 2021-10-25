#ifndef _F4D_FIXEDSIZEVALUESTACK_HPP_
#define _F4D_FIXEDSIZEVALUESTACK_HPP_

#include <cassert>
#include <functional>

namespace cobra {

    template <class T>
    class FixedSizeValueStack {

    private:
        T *array = nullptr;
        int begin = 0;
        int capacity;
        std::function<T(int)> initializer = nullptr;

    public:
        explicit FixedSizeValueStack(int dimension, std::function<T(int index)> array_initializer) {
            assert(dimension > 0);
            array = new T[dimension];
            capacity = dimension;
            begin = 0;
            initializer = array_initializer;
            reset();
        }

        virtual ~FixedSizeValueStack() {
            delete[] array;
        }

        FixedSizeValueStack<T> &operator=(const FixedSizeValueStack<T> &other) {
            assert(capacity == other.capacity);

            for (int i = 0; i < capacity; i++) {
                array[i] = other.array[i];
            }
            begin = other.begin;
            initializer = other.initializer;
            return *this;
        }

        T get() {
            assert(begin < capacity);
            auto item = array[begin];
            begin++;
            return item;
        }

        void push(T item) {
            begin--;
            assert(begin >= 0);
            array[begin] = item;
        }

        void reset() {
            for (int i = 0; i < capacity; i++) {
                array[i] = initializer(i);
            }
            begin = 0;
        }

        int size() const {
            return capacity - begin;
        }

        bool is_empty() {
            return begin == capacity;
        }
    };

}  // namespace cobra

#endif