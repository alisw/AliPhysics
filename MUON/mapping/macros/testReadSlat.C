// $Id$
// $MpId: testReadSlat.C,v 1.4 2005/09/19 19:02:53 ivana Exp $

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
#include <TObjArray.h>
#include <TObjString.h>

#endif

void testReadSlat()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  AliMpSt345Reader r(dataStreams, motifMap);

  ifstream in("slats.list");
  char line[80];
  TObjArray slatsToTest;
  
  while ( in.getline(line,80) )
  {
    slatsToTest.AddLast(new TObjString(line));
  }
  
  in.close();
  
  for ( Int_t i = 0; i < slatsToTest.GetEntriesFast(); ++i )
  {
    TString slat( ((TObjString*)slatsToTest[i])->String());
    
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
      cout << "Bending : ";
      b->Print();
      cout << "NonBending : ";
      nb->Print();
    }
  }
  
  slatsToTest.SetOwner(kTRUE);
  slatsToTest.Delete();
}
