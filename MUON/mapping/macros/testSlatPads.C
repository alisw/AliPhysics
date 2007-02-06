// $Id$
// $MpId: testSlatPads.C,v 1.1 2005/09/19 19:02:53 ivana Exp $

/// Macro to check segmentation of slats.
///
/// What is tested :
/// - first, that we can read both bending and non-bending for each slat
/// - next, for each plane (b,nb), we retrieve all the pads
///   in two different ways a) using a loop and HasPad method
///   b) and using pad iterator. We check that both lead to the same set of pads
///   and that each pad is retrieved only once ;-)
/// - finally we do a "circular test" for all pads, i.e. we get a pad p1
///   by indices and use p1's position to feed PadByPosition which gives us a
///   pad p2. p1 should be equal to p2, of course.
///
/// Usage : .L testSlatPads.C++
///         testSlatPads("slats.list"); where slats.list contains
///         the list of slats' name to be tested, e.g. 112233NR3, etc...

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpPad.h"
#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "AliMpSt345Reader.h"
#include "AliMpSlatMotifMap.h"
#include "AliMpPlaneType.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TString.h"
#include "AliMpIntPair.h"
#include "TStopwatch.h"
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
        AliMpPad xcheck = seg.PadByPosition(pad.Position(),kFALSE);
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
  Int_t n2 = seg.NofPads();
  if ( n2 != n ) 
  {
    cout << Form("Count (%d) and NofPads (%d) give different results for slat %s",
                 n,n2,slat.GetName()) << endl;
  }
  return n;
}

//______________________________________________________________________________
Int_t Iterate(const AliMpSlat& slat)
{
  AliMpSlatSegmentation seg(&slat);
  
  AliMpVPadIterator* it = seg.CreateIterator();
  
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
Int_t Contains(const TList& list, const AliMpPad& pad)
{
  TIter next(&list);
  AliMpPad* p;
  Int_t n(0);
  
  while ( ( p = (AliMpPad*)next() ) )
  {
    if ( pad == *p ) ++n;
  }
  return n;
}

//______________________________________________________________________________
void CheckMultiple(const TList& list, const char* title, const char* what)
{
  // check whether the list contains each pad only once
  
  TIter next(&list);
  AliMpPad* pad;
  while ( ( pad = (AliMpPad*)next() ) )
  {
    if ( Contains(list,*pad) != 1 )
    {
      cout << title << " " << what << " pad found more than once : " << endl;
      pad->Print();
    }
  }
      
}

//______________________________________________________________________________
void XCheck(const AliMpSlat& slat)
{
  // find out which pads are not found by iterator, as compared to HasPad method
  TList l1,l2;
  l1.SetOwner(kTRUE);
  l2.SetOwner(kTRUE);
  
  AliMpSlatSegmentation seg(&slat);
  
  for ( Int_t i = 0; i <= seg.MaxPadIndexX(); ++i )
  {
    for ( Int_t j = 0; j <= seg.MaxPadIndexY(); ++j )
    {
      AliMpPad pad = seg.PadByIndices(AliMpIntPair(i,j),kFALSE);
      if ( pad.IsValid() )
      {
        l1.Add(new AliMpPad(pad));
      }
    }
  }
  
  AliMpVPadIterator* it = seg.CreateIterator();
  
  it->First();
    
  while ( !it->IsDone() )
  {
    l2.Add(new AliMpPad(it->CurrentItem()));
    it->Next();
  }

  TIter next(&l1);
  AliMpPad* pad;
  while ( ( pad = (AliMpPad*)next() ) )
  {
    if ( Contains(l2,*pad) != 1)
    {
      cout << "The following pad is not found by iterator : " << endl;
      pad->Print();
    }
  }
  
  CheckMultiple(l2,slat.GetName(),"iterator");

  CheckMultiple(l1,slat.GetName(),"padByIndices");
}

//______________________________________________________________________________
void testSlatPads(const char* slatlist)
{
  ifstream in(slatlist);
  char line[80];
  TObjArray slatsToTest;
  slatsToTest.SetOwner(kTRUE);
  
  TStopwatch timerCount;
  TStopwatch timerIterate;
  TStopwatch timerCT;
  
  timerCount.Start(true); timerCount.Stop();
  timerIterate.Start(true); timerIterate.Stop();
  timerCT.Start(true); timerCT.Stop();
  
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '#' ) continue;
    slatsToTest.AddLast(new TObjString(line));
  }
  
  in.close();
  
  AliMpSlatMotifMap mm;
  AliMpSt345Reader reader(mm);
  
  for ( Int_t i = 0; i < slatsToTest.GetEntriesFast(); ++i )
  {
    TString slatName( ((TObjString*)slatsToTest[i])->String());
    
    AliMpSlat* bending = reader.ReadSlat(slatName.Data(),AliMp::kBendingPlane);
    AliMpSlat* nonbending = reader.ReadSlat(slatName.Data(),AliMp::kNonBendingPlane);
  
    timerCount.Start(false);
    Int_t NumberOfBendingPads = Count(*bending);
    Int_t NumberOfNonBendingPads = Count(*nonbending);
    timerCount.Stop();
    
    timerIterate.Start(false);
    Int_t xcheck_b = Iterate(*bending);   
    Int_t xcheck_nb = Iterate(*nonbending);
    timerIterate.Stop();
    
    timerCT.Start(false);
    Int_t nc_b = CircularTest(*bending);   
    Int_t nc_nb = CircularTest(*nonbending);
    timerCT.Stop();
    
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

    XCheck(*nonbending);
    XCheck(*bending);

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
  
  cout << "Count : "; 
  timerCount.Print();
  cout << "Iterate : ";
  timerIterate.Print();
  cout << "CT : ";
  timerCT.Print();
  
}
