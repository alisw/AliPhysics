// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSlatMotifMap.h"
#include "AliMpSt345Reader.h"
#include "AliMpSlat.h"
#include "AliMpSt345Reader.h"

#include <Riostream.h>
#include "slats.h"

#endif

void testSt345ReadSlat()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  AliMpSt345Reader r(dataStreams, motifMap);

  Int_t ok(0);
  
  for ( Int_t i = 0; i < NSLATS; ++i )
  {
    TString slat(slatTypeNames[i]);
    
    cout << "Trying to read " << slat << endl;
    AliMpSlat* b = r.ReadSlat(slat.Data(),AliMp::kBendingPlane);
    AliMpSlat* nb = r.ReadSlat(slat.Data(),AliMp::kNonBendingPlane);
    if ( !b ) cout << " Missing BENDING !" << endl;
    if ( !nb ) cout << " Missing NONBENDING !" << endl;
    if ( b && nb )
    {
      if ( b->GetSize() != nb->GetSize() ) 
      {
        cout << "NOT THE SAME NUMBER OF PCBS !" << endl;
      }
      if ( b->DX() != nb->DX() )
      {
        cout << "NOT THE SAME X-SIZE !" << endl;
      }
      cout << "Bending    : ";
      b->Print();
      cout << "NonBending : ";
      nb->Print();
      ++ok;
    }
  }
  
  if ( ok == NSLATS ) 
  {
    cout << "Successfully read " << ok << " slats" << endl;
  }
  else
  {
    cout << "Failed to read " << (NSLATS-ok) << " slats out of " << NSLATS << endl;
  }
}
