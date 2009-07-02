#include <STEER/AliFMDFloatMap.h>
#include <cstdio>
#include <STEER/AliLog.h>

//____________________________________________________________________
struct Printer : public AliFMDMap::ForOne
{
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
  {
    printf("FMD%d%c[%2d,%3d] = %8.3f\n", d, r, s, t, v);
    return kTRUE;
  }
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Int_t v)
  {
    printf("FMD%d%c[%2d,%3d] = %d\n", d, r, s, t, v);
    return kTRUE;
  }
};
//____________________________________________________________________
struct Tester : public AliFMDMap::ForOne
{
  Tester(AliFMDFloatMap& m) : fMap(&m) { }
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t)
  {
    Int_t    idx = fMap->CalcIndex(d, r, s, t);
    UShort_t rd, rs, rt;
    Char_t   rr;
    fMap->CalcCoords(idx, rd, rr, rs, rt);
	  
    if (d != rd || r != rr || s != rs || t != rt) 
      printf("Mismatch FMD%d%c[%2d,%3d] -> %5d -> FMD%d%c[%2d,%3d]\n",
	     d, r, s, t, idx, rd, rr, rs, rt);
    return kTRUE;
  }
  Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Int_t){return kTRUE;}
  AliFMDFloatMap* fMap;
};
//____________________________________________________________________
struct Filler : public AliFMDMap::ForOne
{
  Filler(AliFMDFloatMap& map) : fMap(&map) {}
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t)
  {
    fMap->operator()(d, r, s, t) = MakeVal(d, r, s, t);
    return kTRUE;
  }
  Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Int_t)
  {
    return kTRUE;
  }
  static Float_t MakeVal(UShort_t d, Char_t r, UShort_t s, UShort_t t)
  {
    UShort_t ir  = r == 'I' ? 0 : 1;
    Float_t  val = d * 1000 + ir * 100 + s + 0.001 * t;
    return val;
  }
  AliFMDFloatMap* fMap;
};

//____________________________________________________________________
struct Unity : public AliFMDMap::ForOne
{
  Unity(AliFMDFloatMap& map) : fMap(&map) {}
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t)
  {
    fMap->operator()(d, r, s, t) = 1;
    return kTRUE;
  }
  Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Int_t) {return kTRUE;}
  AliFMDFloatMap* fMap;
};
  
//____________________________________________________________________
void
FillMap(AliFMDFloatMap& m, Bool_t useFiller=kTRUE)
{
  if (useFiller) {
    Filler f(m);
    m.ForEach(f);
    return;
  }
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRng = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRng; q++) { 
      Char_t   r    = (q == 0 ? 'I' : 'O');
      UShort_t nSec = 1; // (q == 0 ?  20 :  40);
      UShort_t nStr = (q == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nSec; s++) { 
	for (UShort_t t = 0; t < nStr; t++) {
	  m(d, r, s, t) = Filler::MakeVal(d, r, s, t);
	}
      }
    }
  }
}
//____________________________________________________________________
void
PrintMap(AliFMDFloatMap& map)
{
  Printer p;
  map.ForEach(p);
}
//____________________________________________________________________
void
FillOne(AliFMDFloatMap& map)
{
  Unity u(map);
  map.ForEach(u);
}
  
//____________________________________________________________________
void
TestIndex(AliFMDFloatMap& map, Bool_t useTester=kTRUE)
{
  if (useTester) { 
    Tester t(map);
    map.ForEach(t);
    return;
  }
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRng = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRng; q++) { 
      Char_t   r    = (q == 0 ? 'I' : 'O');
      UShort_t nSec = 1; // (q == 0 ?  20 :  40);
      UShort_t nStr = (q == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nSec; s++) { 
	for (UShort_t t = 0; t < nStr; t++) {
	  Int_t idx = map.CalcIndex(d, r, s, t);
	  
	  UShort_t rd, rs, rt;
	  Char_t   rr;
	  map.CalcCoords(idx, rd, rr, rs, rt);
	  
	  if (d != rd || r != rr || s != rs || t != rt) 
	    printf("Mismatch FMD%d%c[%2d,%3d] -> %5d -> FMD%d%c[%2d,%3d]\n",
		   d, r, s, t, idx, rd, rr, rs, rt);
	}
      }
    }
  }
}

void 
TestFloatMap()
{
  // AliLog::SetModuleDebugLevel("FMD", 1);
  AliFMDFloatMap m1(0, 0, 0, 0);
  FillMap(m1);
  // PrintMap(m1);
  TestIndex(m1);

  AliFMDFloatMap m2;
  FillOne(m2);
  // m2 *= m1;
  // m2 /= m1;
  m2 += m2;
  // PrintMap(m2);
  
  AliFMDFloatMap m3(3, 2, 1, 512);
  TestIndex(m3, kFALSE); 
  FillMap(m3, kFALSE);
  // PrintMap(m3);
  FillOne(m2);
  m2 *= m3;
  PrintMap(m2);
}
