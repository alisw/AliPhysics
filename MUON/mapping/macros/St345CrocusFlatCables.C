///
/// Macro to help the LC2 user to go from the (local) bus patch id presented
/// by LC2 to a "human-redeable" location in the slat stations 
///
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliDAQ.h"
#include "Riostream.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpBusPatch.h"

#endif

void St345CrocusFlatCables()
{
  AliMpCDB::LoadDDLStore2();
  
  Int_t offset(AliDAQ::DdlIDOffset("MUONTRK"));
  
  while (1)
  {
    Int_t equipmentId;
    Int_t localBusPatchId;
  
    cout << Form("Equipment Id (from %d to %d, or -1 to quit) ? ",
                 AliDAQ::DdlID("MUONTRK",0),
                 AliDAQ::DdlID("MUONTRK",AliDAQ::NumberOfDdls("MUONTRK")-1));
  
    cin >> equipmentId;

    if ( equipmentId < 0 ) break;
  
    cout << "Local bus patch Id (-1 to quit) ? ";
  
    cin >> localBusPatchId;
  
    if ( localBusPatchId < 0 ) break;

    Int_t busPatchId = AliMpBusPatch::GetGlobalBusID(localBusPatchId,equipmentId-offset);

    AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
  
    cout << endl;
    
    bp->Print();
  
    cout << endl;
  }
  
}
