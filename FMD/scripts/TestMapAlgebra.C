#include <AliFMDFloatMap.h>

struct SetIt : public AliFMDMap::ForOne
{
public:
  SetIt(AliFMDFloatMap* m) : fMap(m) {}
  Bool_t operator()(UShort_t d,Char_t r,UShort_t s,UShort_t t,Float_t)
  {
    UShort_t q = r == 'I' ? 0 : 1;
    Float_t  v = d * 1000 + q * 100 + s + t * 0.001;
    fMap->operator()(d, r, s, t) = v;
    return kTRUE;
  }
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,Int_t)
  {
    return kTRUE;
  }
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,UShort_t)
  {
    return kTRUE;
  }
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,Bool_t)
  {
    return kTRUE;
  }
  AliFMDFloatMap* fMap;
};

void
TestMapAlgebra()
{
  AliFMDFloatMap a(0), b(0);
  a.Reset(1);
  b.Reset(2);
  
  AliFMDFloatMap c = a + b;
  // c.Print();

  SetIt s(&b);
  b.ForEach(s);
  
  c = a * b;
  c.Print("%7.3f ");

  AliFMDFloatMap d(3, 2, 1, 3), e(2, 2, 1, 2);
  d.Reset(1);
  e.Reset(2);

  AliFMDFloatMap f = d + e;
  f.Print("%4.0f ");

}
