// $Id$
// $MpId: testSlatPads.C,v 1.1 2005/09/19 19:02:53 ivana Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpPad.h"
#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "AliMpSt345Reader.h"
#include "AliMpPlaneType.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TString.h"
#include "AliMpIntPair.h"
#include <fstream>
#else
#include "Riostream.h"
#endif

//______________________________________________________________________________
Int_t CircularTest(const AliMpSlat& slat)
{
  Int_t n = -1;
  AliMpSlatSegmentation seg(&slat);
  
  for ( Int_t i = 0; i <= seg.MaxPadIndexX(); ++i )
  {
    for ( Int_t j = 0; j <= seg.MaxPadIndexY(); ++j )
    {
      AliMpPad pad = seg.PadByIndices(AliMpIntPair(i,j),kFALSE);
      
      if ( pad.IsValid() )
      {
        ++n;
        AliMpPad xcheck = seg.PadByPosition(pad.Position()+slat.Dimensions(),kFALSE);
        if ( pad != xcheck ) 
        {
          cout << "(ix,iy)=" << i << "," << j << " ";
          pad.Print();
          cout << endl;
          xcheck.Print();
          return kFALSE;
        }             
      }
    }
  }
  return n;
}

//______________________________________________________________________________
Int_t Count(const AliMpSlat& slat)
{
  Int_t n = 0;
  AliMpSlatSegmentation seg(&slat);
  
  for ( Int_t i = 0; i <= seg.MaxPadIndexX(); ++i )
  {
    for ( Int_t j = 0; j <= seg.MaxPadIndexY(); ++j )
    {
      if ( seg.HasPad(AliMpIntPair(i,j)) ) 
      {
        ++n;
      }
    }
  }
  return n;
}

//______________________________________________________________________________
Int_t Iterate(const AliMpSlat& slat)
{
  AliMpSlatSegmentation seg(&slat);
  
  AliMpArea area(slat.Dimensions(),slat.Dimensions());
  
  AliMpVPadIterator* it = seg.CreateIterator(area);
  
  it->First();
  
  Int_t n = 0;
  
  while ( !it->IsDone() )
  {
    it->Next();
    ++n;
  }
  
  return n;
}

//______________________________________________________________________________
void testSlatPads()
{
  ifstream in("slats.list");
  char line[80];
  TObjArray slatsToTest;
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '#' ) continue;
    slatsToTest.AddLast(new TObjString(line));
  }
  
  in.close();
    
  for ( Int_t i = 0; i < slatsToTest.GetEntriesFast(); ++i )
  {
    TString slatName( ((TObjString*)slatsToTest[i])->String());
    
    AliMpSlat* bending = AliMpSt345Reader::ReadSlat(slatName.Data(),kBendingPlane);
    AliMpSlat* nonbending = AliMpSt345Reader::ReadSlat(slatName.Data(),kNonBendingPlane);

    Int_t NumberOfBendingPads = Count(*bending);
    Int_t NumberOfNonBendingPads = Count(*nonbending);
    
    Int_t xcheck_b = Iterate(*bending);   
    Int_t xcheck_nb = Iterate(*nonbending);
    
    Int_t nc_b = CircularTest(*bending);   
    Int_t nc_nb = CircularTest(*nonbending);
    
    cout << setw(10) << slatName
      << " BENDING : " << setw(5) << NumberOfBendingPads
      << " NONBENDING : " << setw(5) << NumberOfNonBendingPads
      << " CT for " << (nc_b+nc_nb) << " pads ";
    if ( nc_b>0 && nc_nb>0 ) 
    {
      cout << "OK.";
    }
    else
    {
      cout << "FAILED.";
    }
    cout << endl;
    
    if ( xcheck_b != NumberOfBendingPads )
    {
      cout << setw(20) << " Bending : HasPad and Iterator give different results !" 
      << " " << NumberOfBendingPads << " vs " << xcheck_b << endl;
    }
    if ( xcheck_nb != NumberOfNonBendingPads )
    {
      cout << setw(20) << " NonBending : HasPad and Iterator give different results !"
      << " " << NumberOfNonBendingPads << " vs " << xcheck_nb << endl;
    }
    
   }
}

int main()
{
  testSlatPads();
}
