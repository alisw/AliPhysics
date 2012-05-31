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

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSlatMotifMap.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpPad.h"
#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "AliMpSt345Reader.h"
#include "AliMpIntPair.h"
#include "slats.h"

#include <TObjArray.h>
#include <TObjString.h>
#include <TList.h>
#include <TString.h>
#include <TStopwatch.h>
#include <Riostream.h>
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
      AliMpPad pad = seg.PadByIndices(i,j,kFALSE);
      
      if ( pad.IsValid() )
      {
        ++n;
        AliMpPad xcheck 
          = seg.PadByPosition(pad.GetPositionX(),pad.GetPositionY(),kFALSE);
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
      if ( seg.HasPadByIndices(i,j) ) 
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
      AliMpPad pad = seg.PadByIndices(i,j,kFALSE);
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
void testSt345Pads()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  AliMpSt345Reader reader(dataStreams, motifMap);
  
  Int_t ok(0);
  
  for ( Int_t i = 0; i < NSLATS; ++i ) 
  {
    Bool_t slatOK(kTRUE);
    
    TString slatName( slatTypeNames[i] );
    
    AliMpSlat* bending = reader.ReadSlat(slatName.Data(),AliMp::kBendingPlane);
    AliMpSlat* nonbending = reader.ReadSlat(slatName.Data(),AliMp::kNonBendingPlane);
  
    Int_t NumberOfBendingPads(0);
    Int_t NumberOfNonBendingPads(0);
    Int_t xcheck_b(0);
    Int_t xcheck_nb(0);
    Int_t nc_b(0);
    Int_t nc_nb(0);
    
    if ( bending )
    {
      NumberOfBendingPads = Count(*bending);
      xcheck_b = Iterate(*bending);   
      nc_b = CircularTest(*bending);   
    }
    else
    {
      cout << "Could not read bending plane of slat " << slatName.Data() << endl;
    }
    
    if ( nonbending ) 
    {
      NumberOfNonBendingPads = Count(*nonbending);
      xcheck_nb = Iterate(*nonbending);
      nc_nb = CircularTest(*nonbending);
    }
    else
    {
      cout << "Could not read bending plane of slat " << slatName.Data() << endl;
    }
    
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
      slatOK = kFALSE;
      cout << "FAILED.";
    }
    cout << endl;

    if ( nonbending ) XCheck(*nonbending);
    if ( bending ) XCheck(*bending);

    if ( xcheck_b != NumberOfBendingPads )
    {
      cout << setw(20) << " Bending : HasPad and Iterator give different results !" 
      << " " << NumberOfBendingPads << " vs " << xcheck_b << endl;
      slatOK = kFALSE;
    }
    if ( xcheck_nb != NumberOfNonBendingPads )
    {
      cout << setw(20) << " NonBending : HasPad and Iterator give different results !"
      << " " << NumberOfNonBendingPads << " vs " << xcheck_nb << endl;
      slatOK = kFALSE;
    }

    if (slatOK) ++ok;
   }
  
  if ( ok == NSLATS ) 
  {
    cout << "Successfully tested " << ok << " slats" << endl;
  }
  else
  {
    cout << "Failed to read " << (NSLATS-ok) << " out of " << NSLATS << " slats" << endl;
  }
}
