#include <TRandom.h>
#include <TStopwatch.h>
#include <AliFMDFloatMap.h>

Float_t 
TestValue(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
{
  return d * 1000000 + ((r=='I')+1) * 100000 + s * 1000 + t + v;
}
  
struct TestOne : public AliFMDMap::ForOne
{
  Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
  {
    // Caculate a simple number from d, r, s, t, v 
    Float_t x = TestValue(d, r, s, t, v);
    return x > -1000000000; // Always true
  }
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,Int_t) {return true;}
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,Bool_t) {return true;}
  Bool_t operator()(UShort_t,Char_t,UShort_t,UShort_t,UShort_t) {return true;}
};

void 
AccessOneByOne(AliFMDFloatMap& m)
{
  TestOne p;
  m.ForEach(p);
}

void
AccessByCoords(AliFMDFloatMap& m)
{
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nr; q++) { 
      Char_t    r = (q == 0 ? 'I' : 'O');
      UShort_t ns = (q == 0 ?  20 :  40);
      UShort_t nt = (q == 0 ? 512 : 256);
      for (UShort_t s = 0; s < ns; s++) {
	for (UShort_t t = 0; t < nt; t++) { 
	  Float_t x = TestValue(d, r, s, t, m(d,r,s,t));
	  if (x < -1000000000) break;
	}
      }
    }
  }
}
	
void
FillRandom(AliFMDFloatMap& m) 
{
  Float_t* data = m.Data();
  UInt_t   max  = 51200; // m.MaxIndex();

  for (UInt_t i = 0; i < max; i++) 
    data[i] = gRandom->Uniform(1);
}
      
void
TestMapAccess(Int_t n=1000)
{
 
  AliFMDFloatMap m(0);
  TStopwatch s;

  s.Reset();
  for (Int_t i = 0; i < n; i++) { 
    FillRandom(m);
    s.Start(kFALSE);
    AccessOneByOne(m);
    s.Stop();
  }
  s.Print();

  s.Reset();
  for (Int_t i = 0; i < n; i++) { 
    FillRandom(m);
    s.Start(kFALSE);
    AccessByCoords(m);
    s.Stop();
  }
  s.Print();
}

#ifndef __CINT__
#include <sstream>

int
main(int argc, char** argv)
{
  int n  = 1000;
  if (argc > 1) { 
    std::stringstream s(argv[1]);
    s >> n;
  }
  TestMapAccess(n);
  return 0;
}
#endif
