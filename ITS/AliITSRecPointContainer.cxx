#include <TTree.h>
#include "AliLog.h"
#include "AliITSRecoParam.h"
#include "AliITSReconstructor.h"
#include "AliITSRecPointContainer.h"
#include "AliITSRecPoint.h"
#include "AliRunLoader.h"

ClassImp(AliITSRecPointContainer)

//////////////////////////////////////////////////////////////////////
// Class to store ITS RecPoints for the duration of                 //
// one event processing                                             //
// The container is cleared at each event and new RP                //
// are loaded from TTree                                            //
// Origin masera@to.infn.it  Nov. 12 2009                           //
//////////////////////////////////////////////////////////////////////

/* $Id$ */

AliITSRecPointContainer* AliITSRecPointContainer::fgInstance = 0x0;

//______________________________________________________________________
AliITSRecPointContainer::AliITSRecPointContainer(const AliITSRecoParam* krp):TObject(),
fSPDNModules(0),
fSDDNModules(0),
fSSDNModules(0),
fArray(),
fCurrentEve(-1000),
fActualSize(0),
fDet(""),
fStatusOK(kTRUE){
  // Default constructor

  if(fgkNModules != AliITSgeomTGeo::GetNModules())AliError(Form("The total number of modules is not %d, but %d",fgkNModules,AliITSgeomTGeo::GetNModules()));

  Int_t modperlay[6];
  for(Int_t i=0;i<6;i++)modperlay[i]=AliITSgeomTGeo::GetNDetectors(1+i)*AliITSgeomTGeo::GetNLadders(1+i);
  fSPDNModules=modperlay[0]+modperlay[1];
  fSDDNModules=modperlay[2]+modperlay[3];
  fSSDNModules=modperlay[4]+modperlay[5];
  //  AliInfo(Form("Total modules: %d \n SPD modules=%d , SDD modules=%d, SSD modules=%d ",fgkNModules,fSPDNModules,fSDDNModules,fSSDNModules));

  // kLimits[0:5] --> low fluw; kLimits[6,11] --> High flux
  const Int_t kLimits[12]={25,25,20,20,10,10,300,300,200,200,100,100};
  Int_t offset=0;
  if(!krp){
    AliWarning("AliITSRecoPoint is missing. Using defaults");
  }
  else {
    if(krp->GetEventSpecie() == AliRecoParam::kHighMult)offset=6;
  }
  Int_t maxval[6];
  TString values="";
  for(Int_t i=0;i<6;i++){
    maxval[i]=kLimits[i+offset];
    values+=maxval[i];
    values+=" ";
    if(i>0)modperlay[i]+=modperlay[i-1];
  }
  AliInfo(Form("Container created with sizes/layer: %s",values.Data()));
  Int_t layer=0;
  for(Int_t i=0;i<fgkNModules;i++){
    if(i>=modperlay[layer])++layer;
    fArray[i]=new TClonesArray("AliITSRecPoint",maxval[layer]);
  }
}


//______________________________________________________________________
AliITSRecPointContainer::~AliITSRecPointContainer(){
  // Destructor
  for(Int_t i=0;i<fgkNModules;i++){
    if(fArray[i]){
      fArray[i]->Delete();
      delete fArray[i];
    }
  }
}

//______________________________________________________________________
void AliITSRecPointContainer::CookEntries(){
  // From the number of entries in TTree R, the number of ITS subdetectors
  // active for the present run is inferred
  if(fActualSize == fgkNModules)fDet="ALL SPD SDD SSD ";
  if(fActualSize == fSPDNModules) fDet = "SPD ";
  if(fActualSize == fSDDNModules) fDet = "SDD ";
  if(fActualSize == fSSDNModules)fDet = "SSD ";
  if(fActualSize == (fSPDNModules+fSDDNModules)) fDet = "SPD SDD ";
  if(fActualSize == (fSPDNModules+fSSDNModules))fDet = "SPD SSD ";
  if(fActualSize == (fSDDNModules+fSSDNModules))fDet = "SDD SSD ";
  if((!fDet.Contains("SPD")) && (!fDet.Contains("SDD")) &&
     (!fDet.Contains("SSD"))){
    AliError(Form("The number of active modules %d does not correspond to any standard configuration of the detector",fActualSize));
    fStatusOK = kFALSE;
  }
}
//______________________________________________________________________
TClonesArray* AliITSRecPointContainer::FetchClusters(Int_t mod, TTree* tR){
  // retrieves Recpoints for module mod (offline mode: the event number is
  // retrieved via the AliRunLoader object)
  Int_t cureve=AliRunLoader::Instance()->GetEventNumber();
  return FetchClusters(mod,tR,cureve);
}
//______________________________________________________________________
TClonesArray* AliITSRecPointContainer::FetchClusters(Int_t mod, TTree* tR,Int_t cureve){
  // retrieves Recpoints for module mod
  // cureve is the current event number. If it is different w.r.t.
  // the event number stored in fCurrentEve, the recpoints are read from
  // the TTree. Otherwise, the RP stored in memory are used. 
  if(cureve != fCurrentEve){
    fCurrentEve = cureve;
    Reset();
    TBranch *branch = NULL;
    branch = tR->GetBranch("ITSRecPoints");
    if(!branch){
      AliError("Branch ITSRecPoints not found on ITS recpoints TTree");
      fStatusOK = kFALSE;
      return NULL;
    }

    fActualSize = branch->GetEntries();
    CookEntries();
    if(fDet.IsNull())return NULL;

    // it is assumed that the filling order of the tree is SPD, SDD, SSD
    // even if one or two subdetector are missing
    if(IsSPDActive()){
      for(Int_t i=0;i<fSPDNModules;i++){
	branch->SetAddress(&fArray[i]);
	branch->GetEvent(i);
      }
    }
    if(IsSDDActive()){
      Int_t start=0;
      if(IsSPDActive())start+=fSPDNModules;
      Int_t counter = fSPDNModules;
      for(Int_t i=start;i<start+fSDDNModules;i++){
	branch->SetAddress(&fArray[counter]);
	++counter;
	branch->GetEvent(i);
      }
    }
    if(IsSSDActive()){
      Int_t start=0;
      if(IsSPDActive())start+=fSPDNModules;
      if(IsSDDActive())start+=fSDDNModules;
      Int_t counter = fSPDNModules+fSDDNModules;
      for(Int_t i=start;i<start+fSSDNModules;i++){
	branch->SetAddress(&fArray[counter]);
	++counter;
	branch->GetEvent(i);
      }
    }
  }

  if(CheckBoundaries(mod)){
    return fArray[mod];
  }
  else {
    AliError(Form("Module %d is out of boundaries",mod));
    return NULL;
  }
  
}

//______________________________________________________________________
AliITSRecPointContainer* AliITSRecPointContainer::Instance(const AliITSRecoParam* kptr){
  // returns AliITSRecPointContainer instance (singleton)
  if(!fgInstance){
    if(!kptr){
      fgInstance =  new AliITSRecPointContainer(AliITSReconstructor::GetRecoParam());
    }
    else {
    fgInstance = new AliITSRecPointContainer(kptr);
    }
  }
  return fgInstance;
}

//______________________________________________________________________
void AliITSRecPointContainer::Reset(){
  // Resets the status of the object
  for(Int_t i=0;i<fgkNModules;i++){
    (fArray[i])->Clear();
  }
  fDet="";
}
