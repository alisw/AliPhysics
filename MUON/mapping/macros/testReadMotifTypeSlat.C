// $Id$
// $MpId: testReadMotifTypeSlat.C,v 1.1 2005/09/19 19:02:53 ivana Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMpMotifReader.h"
#include "Riostream.h"
#include "AliMpMotifType.h"
#endif

Int_t test(AliMpMotifReader& r, const char letter, Int_t from, Int_t to)
{
  char m[256];
  Int_t n = 0;
  for ( Int_t i = from; i <= to; ++i ) 
  {
    sprintf(m,"%c%d",letter,i);
    AliMpMotifType* motifType = r.BuildMotifType(m);
    if ( motifType ) 
    {
      motifType->Print("G");
      ++n;
    }
    else
    {
      cout << "Cannot read motifType " << m << endl;
    }
  }
  return n;
}

void testReadMotifTypeSlat()
{
  AliMpMotifReader r(AliMp::kStation345,AliMp::kNonBendingPlane); 
  // note that second parameter is not used for station345.

  Int_t n = 0;
  
  n += test(r,'I',1,1);
  n += test(r,'L',1,20);
  n += test(r,'O',1,19);
  n += test(r,'P',1,4);
  n += test(r,'Q',1,4);
  n += test(r,'R',1,42);
  n += test(r,'Z',1,5);
  
  cout << "==== " << n << " motifTypes successfully read in" << endl;
}  

