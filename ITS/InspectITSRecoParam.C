#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TGrid.h"
#include "TString.h"
#include "AliCDBManager.h"
#include "AliITSRecoParam.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#endif
AliITSRecoParam* InspectITSRecoParam(Int_t nrun=167713,TString selec="HighMult", Bool_t local=kFALSE ){
  // this macro retrieves one of the 3 recoparam objects stored in the OCDB:
  // according to the selection string selec:  HighMult - LowMult - Cosmic
  // if local is true, then the OCDB in $ALICE_ROOT is used
  AliCDBManager * man = AliCDBManager::Instance();
  if(local){
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  else{
    TGrid::Connect("alien://");
    man->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  }
  man->SetRun(nrun);
  AliCDBPath path("ITS","Calib","RecoParam");
  AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
  TObjArray *arr = dynamic_cast<TObjArray*>(entry->GetObject());
  if(!arr){
    cout<<"No valid TObjArray\n";
    return NULL;
  }
  Int_t elem=0;
  if(selec.Contains("Cosmic")){
    elem = 0;
  }
  else if (selec.Contains("LowMult")){
    elem = 1;
  }
  else if (selec.Contains("HighMult")){
    elem = 2;
  }
  else {
    cout<<"Invalid choice "<<selec<<endl;
    return NULL;
  }
  AliITSRecoParam* rp = (AliITSRecoParam*)arr->At(elem);
  rp->PrintParameters();
  return rp;
}
