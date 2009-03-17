// $Id$
// $MpId: testReadPCB.C,v 1.1 2005/09/19 19:02:53 ivana Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSlatMotifMap.h"
#include "AliMpSt345Reader.h"
#include "AliMpPCB.h"
#include "AliMpMotifPosition.h"

#include <Riostream.h>

#endif

void testReadPCB()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  AliMpSt345Reader r(dataStreams, motifMap);

  const char* pcbToTest[] = { "B1", "B2", "B3+", "B3-", "N1", "N2+", "N2-", 
  "N3", "R1B", "R1N", "R2B", "R2N", "R3B", "R3N", "S2B", "S2N" };
  
  Int_t N = sizeof(pcbToTest)/sizeof(const char*);
  
  for ( Int_t i = 0; i < N; ++i )
  {
    AliMpPCB* pcb = r.ReadPCB(pcbToTest[i]);
    if (pcb)
    {
      pcb->Print();
      for ( Int_t j = 0; j < pcb->GetSize(); ++j )
      {
        AliMpMotifPosition* pos = pcb->GetMotifPosition(j);
        cout << "    " << j << " ";
        pos->Print();
      }
    }
    else
    {
      cout << "Cannot read " << pcbToTest[i] << endl;
    }
  }  
}
