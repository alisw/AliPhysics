/*  This file is part of the Vc library. {{{

    Copyright (C) 2012 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.

}}}*/

#ifndef VC_SSE_ITERATORS_H
#define VC_SSE_ITERATORS_H

namespace Vc
{
namespace Iterators
{

template<typename T> class VectorRef : public Vector<T>
{
private:
    typedef Vector<T> Base;
public:
    typedef typename Base::EntryType EntryType;
private:
    EntryType *VC_RESTRICT const m_memoryLocation;

public:
    VectorRef(EntryType *m)
        : Base(),
        m_memoryLocation(m)
    {
        Base::load(m_memoryLocation);
    }

    VectorRef(VectorRef &&a) = default;

    ~VectorRef()
    {
        Base::store(m_memoryLocation);
    }
private:
    // disable copies
    VectorRef &operator=(const VectorRef &) = delete;
    VectorRef(const VectorRef &) = delete;
};

template<typename V> class ConstRangeIterator
{
protected:
    typedef typename V::EntryType EntryType;
    size_t m_offset;

public:
    ConstRangeIterator(size_t o)
        : m_offset(o)
    {
    }

    const Vector<V> operator()(const EntryType *mem) const
    {
        return Vector<V>(mem);
    }

    bool operator!=(const ConstRangeIterator &rhs) const
    {
        return m_offset != rhs.m_offset;
    }

    ConstRangeIterator &operator++()
    {
        ++m_offset;
        return *this;
    }
    ConstRangeIterator operator++(int) {
        ConstRangeIterator ret(*this);
        ++m_offset;
        return ret;
    }
    ConstRangeIterator &operator--()
    {
        --m_offset;
        return *this;
    }
    ConstRangeIterator operator--(int) {
        ConstRangeIterator ret(*this);
        --m_offset;
        return ret;
    }
    ConstRangeIterator &operator*() {
        return *this;
    }
};

template<typename V> class RangeIterator : public ConstRangeIterator<V>
{
    typedef typename ConstRangeIterator<V>::EntryType EntryType;
public:
    RangeIterator(size_t o) : ConstRangeIterator<V>(o) {}

    VectorRef<V> operator()(EntryType *mem) const
    {
        return VectorRef<V>(mem);
    }
    RangeIterator &operator*() {
        return *this;
    }
};

template<typename V> class RangeIteration
{
private:
    size_t m_begin;
    size_t m_end;

public:
    RangeIteration(size_t b, size_t e)
        : m_begin(b),
        m_end(e)
    {
    }

    RangeIterator<V> begin() { return RangeIterator<V>(m_begin); }
    RangeIterator<V> end() { return RangeIterator<V>(m_end); }
    RangeIterator<V> begin() const { return ConstRangeIterator<V>(m_begin); }
    RangeIterator<V> end() const { return ConstRangeIterator<V>(m_end); }
};
} // namespace Iterators
} // namespace Vc

namespace std
{
    template<typename V>
    Vc::Iterators::RangeIterator<V> begin(Vc::Iterators::RangeIteration<V> &range)
    {
        return range.begin();
    }
    template<typename V>
    Vc::Iterators::RangeIterator<V> end(Vc::Iterators::RangeIteration<V> &range)
    {
        return range.end();
    }
} // namespace std

#endif // VC_SSE_ITERATORS_H
