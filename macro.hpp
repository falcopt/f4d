#ifndef _F4D_MACRO_HPP_
#define _F4D_MACRO_HPP_

#define likely(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define unlikely(condition) __builtin_expect(static_cast<bool>(condition), 0)

#endif