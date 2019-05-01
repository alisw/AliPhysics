#ifndef Lifetimes_Utils_h
#define Lifetimes_Utils_h

#include <limits>

namespace Lifetimes
{

template <typename F, typename I>
F getBinCenter(I bin, F binw, F min, bool checkF = false,
               bool checkL = false)
{
  if (checkF && bin == 0)
    return -1.e32;
  else if (checkL && bin == std::numeric_limits<I>::max())
    return 1.e32;
  else
  {
    return min + (bin - 0.5) * binw;
  }
}

template <typename I, typename F>
I getBinnedValue(F val, F binw)
{
  return 1 + I(std::floor(val / binw));
}

template <typename I, typename F>
I getBinnedValue(F val, F binw, F max)
{
  return (val > max) ? std::numeric_limits<I>::max()
                     : 1 + I(std::floor(val / binw));
}

template <typename I, typename F>
I getBinnedValue(F val, F binw, F min, F max)
{
  if (val < min)
    return 0;
  else if (val > max)
    return std::numeric_limits<I>::max();
  else
    return 1 + I(std::floor((val - min) / binw));
}

template <typename I>
I flipBits(I mask, I bits, bool flag)
{
  return (mask & ~bits) | (-flag & bits);
}

} // namespace Lifetimes

#endif