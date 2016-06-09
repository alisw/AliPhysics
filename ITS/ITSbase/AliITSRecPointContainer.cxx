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
fNextEvent(-1000),
fActualSize(0),
fDet(""),
fStatusOK(kTRUE){
  // Default constructor

  for(Int_t i=0;i<6;i++)fNClusters[i]=0;
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
    AliWarning("AliITSRecoParam is missing. Using defaults");
  }
  else {
    if(krp->GetEventSpecie() & AliRecoParam::kHighMult)offset=6;
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
  // The actual access to the RP TTree is done as follows:
  // If the AliRunLoader object exists, the event number is taken from it
  // If not, the data member fNextEvent is used. 
  // To set fNextEvent it is necessary to call PrepareToRead in advance.
  // if this is never done, fNextEvent will have its default negative value
  // and an error message will be delivered.
  AliRunLoader* rl = AliRunLoader::Instance();
  Int_t cureve;
  if(rl){
    cureve = rl->GetEventNumber();
  }
  else if(fNextEvent>=0){
    cureve = fNextEvent;
  }
  else {
    AliError("The RunLoader is not defined, PrepareToRead was not invoked. Revise calling sequence. Nothing done");
    return NULL;
  }
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
    Int_t modL1=AliITSgeomTGeo::GetNDetectors(1)*AliITSgeomTGeo::GetNLadders(1);
    if(IsSPDActive()){
      for(Int_t i=0;i<fSPDNModules;i++){
	branch->SetAddress(&fArray[i]);
	branch->GetEvent(i);
	if(i<modL1){
	  fNClusters[0]+=fArray[i]->GetEntries();
	}
	else {
	  fNClusters[1]+=fArray[i]->GetEntries();
	}
      }
    }
    if(IsSDDActive()){
      Int_t start=0;
      if(IsSPDActive())start+=fSPDNModules;
      Int_t modL3=AliITSgeomTGeo::GetNDetectors(3)*AliITSgeomTGeo::GetNLadders(3);
      Int_t counter = fSPDNModules;
      for(Int_t i=start;i<start+fSDDNModules;i++){
	branch->SetAddress(&fArray[counter]);
	++counter;
	branch->GetEvent(i);
	if((i-start)<modL3){
	  fNClusters[2]+=fArray[i]->GetEntries();
	}
	else {
	  fNClusters[3]+=fArray[i]->GetEntries();
	}
      }
    }
    if(IsSSDActive()){
      Int_t start=0;
      if(IsSPDActive())start+=fSPDNModules;
      if(IsSDDActive())start+=fSDDNModules;
      Int_t modL5=AliITSgeomTGeo::GetNDetectors(5)*AliITSgeomTGeo::GetNLadders(5);
      Int_t counter = fSPDNModules+fSDDNModules;
      for(Int_t i=start;i<start+fSSDNModules;i++){
	branch->SetAddress(&fArray[counter]);
	++counter;
	branch->GetEvent(i);
	if((i-start)<modL5){
	  fNClusters[4]+=fArray[i]->GetEntries();
	}
	else {
	  fNClusters[5]+=fArray[i]->GetEntries();
	}
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
UInt_t AliITSRecPointContainer::GetNClustersInLayer(Int_t lay, TTree* tR, Int_t eventN){
  // returns the number of clusters for laier lay
  // layers are numbered from 1 to 6
  if(lay<1 || lay >6){
    AliError(Form("Layer %d is out of range",lay));
    return 0;
  }
  if(eventN>=0){
    FetchClusters(0,tR,eventN);
  }
  else {
    FetchClusters(0,tR);
  }
  return fNClusters[lay-1];
}
//______________________________________________________________________
UInt_t AliITSRecPointContainer::GetNClustersInLayerFast(Int_t lay) const {
  // returns the number of clusters for laier lay
  // layers are numbered from 1 to 6
  // No checks are done on the event number: the numer of clusters 
  // for the event stored in memory is returned
  if(lay<1 || lay >6){
    AliError(Form("Layer %d is out of range",lay));
    return 0;
  }
  return fNClusters[lay-1];
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
void AliITSRecPointContainer::Destroy(){
  // deletes the singleton
  if(fgInstance){
    delete fgInstance;
    fgInstance = NULL;
  }
}

//______________________________________________________________________
void AliITSRecPointContainer::Reset(){
  // Resets the status of the object
  ClearClus(0,fgkNModules);
  fDet="";
  for(Int_t i=0;i<6;i++)fNClusters[i]=0;
}
//______________________________________________________________________
void AliITSRecPointContainer::ResetSPD(){
  // Resets only the entries in fArray concerning SPD
  // This method should be used with care only when the filling
  // of the container is not done from the RP TTree. 
  fCurrentEve = -1000;  // protection: if FetchClusters method will be used
                          // after this call, an ccess to the RP TTree will
                          // be forced
  ClearClus(0,fSPDNModules);
}

//______________________________________________________________________
void AliITSRecPointContainer::ResetSDD(){
  // Resets only the entries in fArray concerning SDD
  // This method should be used with care only when the filling
  // of the container is not done from the RP TTree. 
  fCurrentEve = -1000;  // protection: if FetchClusters method will be used
                          // after this call, an ccess to the RP TTree will
                          // be forced
  Int_t first = fSPDNModules;
  Int_t last = first + fSDDNModules; 
  ClearClus(first,last);
}

//______________________________________________________________________
void AliITSRecPointContainer::ResetSSD(){
  // Resets only the entries in fArray concerning SSD
  // This method should be used with care only when the filling
  // of the container is not done from the RP TTree. 
  fCurrentEve = -1000;  // protection: if FetchClusters method will be used
                          // after this call, an ccess to the RP TTree will
                          // be forced
  Int_t first = fSPDNModules + fSDDNModules;
  Int_t last = first + fSSDNModules; 
  ClearClus(first,last);
}

