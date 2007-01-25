// $Id$
// $MpId: testReadSlat.C,v 1.4 2005/09/19 19:02:53 ivana Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "AliMpSlat.h"
#include "AliMpSt345Reader.h"
#include "TObjArray.h"
#include "TObjString.h"
#endif

void testReadSlat()
{
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
    AliMpSlat* b = AliMpSt345Reader::ReadSlat(slat.Data(),AliMp::kBendingPlane);
    AliMpSlat* nb = AliMpSt345Reader::ReadSlat(slat.Data(),AliMp::kNonBendingPlane);
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
