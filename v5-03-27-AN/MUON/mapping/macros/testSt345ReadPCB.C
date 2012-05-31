// $Id$

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

void testSt345ReadPCB()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  AliMpSt345Reader r(dataStreams, motifMap);

  const char* pcbToTest[] = { "B1", "B2", "B3+", "B3-", "N1", "N2+", "N2-", 
  "N3", "R1B", "R1N", "R2B", "R2N", "R3B", "R3N", "S2B-", "S2B+", "S2N" };
  
  Int_t N = sizeof(pcbToTest)/sizeof(const char*);
  Int_t ok(0);
  
  for ( Int_t i = 0; i < N; ++i )
  {
    AliMpPCB* pcb = r.ReadPCB(pcbToTest[i]);
    if (pcb)
    {
      ++ok;
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
  if ( ok == N ) 
  {
    cout << "Successfully read " << ok << " PCBs" << endl;
  }
  else
  {
    cout << "Failed to read " << (N-ok) << " PCBs out of " << N << endl;
  }
}
