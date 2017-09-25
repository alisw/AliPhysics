#include <TObject.h>
#include <TClonesArray.h>
#include <TH2.h>
#include <AliLog.h>
#include <AliPHOSGeometry.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliVCaloTrigger.h>
#include <AliVCluster.h>
#include <AliCaloPhoton.h>
#include <AliOADBContainer.h>

#include "AliPHOSClusterCuts.h"
#include "AliPHOSTriggerHelper.h"

// Author: Daiki Sekihata (Hiroshima University)
ClassImp(AliPHOSTriggerHelper)

//________________________________________________________________________
AliPHOSTriggerHelper::AliPHOSTriggerHelper():
  fPHOSGeo(0x0),
  fXmin(-3),
  fXmax(3),
  fZmin(-3),
  fZmax(3),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
	fTriggerInputL1(-1),//5:1PHL, 6:1PHM, 7:1PHH
	fTriggerInputL0(-1),//9:0PH0
  fIsMC(kFALSE),
  fCaloTrigger(0x0),
  fIsUserTRUBadMap(kFALSE),
  fRunNumber(-1)
{
  //Constructor
  
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

}
//________________________________________________________________________
AliPHOSTriggerHelper::AliPHOSTriggerHelper(TString trigger):
  fPHOSGeo(0x0),
  fXmin(-3),
  fXmax(3),
  fZmin(-3),
  fZmax(3),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
	fTriggerInputL1(-1),//5:1PHL, 6:1PHM, 7:1PHH
	fTriggerInputL0(-1),//9:0PH0
  fIsMC(kFALSE),
  fCaloTrigger(0x0),
  fIsUserTRUBadMap(kFALSE),
  fRunNumber(-1)
{
  //Constructor
  
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

  trigger.ToUpper();

  if(trigger.Contains("L1")){
    fTriggerInputL0 = -1;
    if(trigger.Contains("L1H"))      fTriggerInputL1 = 7;
    else if(trigger.Contains("L1M")) fTriggerInputL1 = 6;
    else if(trigger.Contains("L1L")) fTriggerInputL1 = 5;
    else                             fTriggerInputL1 = -1;
    //STU stores top-left fired channel
    fXmin = -3;
    fXmax = 0;
    fZmin = -1;
    fZmax = 2;
  }
  if(trigger.Contains("L0")){
    fTriggerInputL1 = -1;
    fTriggerInputL0 = 9;
    //TRU stores bottom-left fired channel.
    fXmin = -3;
    fXmax = 0;
    fZmin = -3;
    fZmax = 0;
  }

  if(fTriggerInputL1 > 0 && fTriggerInputL0 > 0){
    AliError("Both L1 and L0 are selected. Analyzer must select either L1 or L0. L1 has higher priority in this class.");
  }

}
//________________________________________________________________________
AliPHOSTriggerHelper::~AliPHOSTriggerHelper()
{

  for(Int_t i=0;i<6;i++){
    if(fPHOSTRUBadMap[i]){
      delete fPHOSTRUBadMap[i];
      fPHOSTRUBadMap[i] = 0x0;
    }
  }

}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsPHI7(AliVEvent *event, AliPHOSClusterCuts *cuts)
{
  fEvent    = dynamic_cast<AliVEvent*>(event);
  fESDEvent = dynamic_cast<AliESDEvent*>(event);
  fAODEvent = dynamic_cast<AliAODEvent*>(event);
  Int_t run = fEvent->GetRunNumber();

  //SetPHOSTRUBadmaps
  if(fRunNumber != run){
    fRunNumber = run;

    if(!fIsUserTRUBadMap){
      AliOADBContainer badmapContainer(Form("phosTriggerBadMap"));
      badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSTrigBadMaps.root","phosTriggerBadMap");
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(run,"phosTriggerBadMap");
      if(!maps) AliError(Form("Can not read Trigger Bad map for run %d.",run));
      else{
        AliInfo(Form("Setting PHOS Trigger bad map with name %s.",maps->GetName()));
        for(Int_t mod=0; mod<6;mod++){
          if(fPHOSTRUBadMap[mod]){
            delete fPHOSTRUBadMap[mod];
            fPHOSTRUBadMap[mod] = 0x0;
          }
          TH2I * h = (TH2I*)maps->At(mod);
          if(h) fPHOSTRUBadMap[mod] = new TH2I(*h);
        }
      }
    }
  }

  //if L1 trigger is used, L1 has higher priority
  Bool_t IsPHI7fired = kFALSE;

  if(fTriggerInputL1 > 0)
    IsPHI7fired = event->GetHeader()->GetL1TriggerInputs() & 1 << (fTriggerInputL1 - 1);//trigger input -1
  else if(fTriggerInputL0 > 0)
    IsPHI7fired = event->GetHeader()->GetL0TriggerInputs() & 1 << (fTriggerInputL0 - 1);//trigger input -1
  else{
    AliInfo("Which PHOS triggered data do you analyze? Please set trigger input (L0 or L1[H,M,L]).");
    return kFALSE;
  }

  if(!IsPHI7fired) return kFALSE;

  if(fTriggerInputL1 > 0)      AliInfo(Form("Your choice of L1 trigger input %d fired.",fTriggerInputL1));
  else if(fTriggerInputL0 > 0) AliInfo(Form("Your choice of L0 trigger input %d fired.",fTriggerInputL0));

  //TString trigClasses = event->GetFiredTriggerClasses();
  //if(245917 <= run && run <= 246994){//LHC15o
  //  if(!fIsMC && fTriggerInputL1 == 7 && !trigClasses.Contains("CINT7PHH")) return kFALSE;
  //  if(!fIsMC && fTriggerInputL1 == 6 && !trigClasses.Contains("CPER7PHM")) return kFALSE;
  //}
  //if(fTriggerInputL1 > 0) AliInfo(Form("Your choice of L1 trigger input %d matches with fired trigger classes : %s",fTriggerInputL1,trigClasses.Data()));

  if(run<209122) //Run1
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");
  else
    fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

  fCaloTrigger = (AliVCaloTrigger*)event->GetCaloTrigger("PHOS"); 
  fCaloTrigger->Reset();

  const Int_t Nfired = fCaloTrigger->GetEntries();
  AliInfo(Form("%d TRU patches fired in PHOS.",Nfired));

  TClonesArray *array = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(array){
    AliInfo("PHOSClusterArray is found!");
  }
  else{
    AliInfo("PHOSClusterArray is NOT found! return kFALSE");
    return kFALSE;
  }
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t relId[4]={};
  AliVCluster *clu1;

  Int_t L1=-999;
  if(fTriggerInputL1 > 0){
    if     (fTriggerInputL1 == 7) L1 = 0;//L1 high
    else if(fTriggerInputL1 == 6) L1 = 1;//L1 medium
    else if(fTriggerInputL1 == 5) L1 = 2;//L1 low
  }
  else if(fTriggerInputL0 > 0){
    if(fTriggerInputL0 == 9) L1 = -1;//L0
  }

  while(fCaloTrigger->Next()){
    // L1 threshold: -1-L0, 0-PHH, 1-PHM, 2-PHL

    if(fCaloTrigger->GetL1TimeSum() != L1){
      //AliInfo(Form("This patch fired as %d, which does not match with your choice %d.",fCaloTrigger->GetL1TimeSum(),L1));
      continue;
    }
    AliInfo(Form("This patch fired as %d.",fCaloTrigger->GetL1TimeSum()));
    fCaloTrigger->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

    fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelId);
    //for offline numbering, relId should be used.


    AliInfo(Form("Fired trigger channel : M%d X%d Z%d",trgrelId[0],trgrelId[2],trgrelId[3]));
    //if(!IsGoodTRUChannel("PHOS",trgrelId[0],trgrelId[2],trgrelId[3])) return kFALSE;

    Int_t multClust = array->GetEntriesFast();
    for(Int_t i=0;i<multClust;i++){
      AliCaloPhoton *ph = (AliCaloPhoton*)array->At(i);
      if(!cuts->AcceptPhoton(ph)) continue;
      if(!ph->IsTOFOK()) continue;
      clu1 = (AliVCluster*)ph->GetCluster();

      Int_t maxAbsId = FindHighestAmplitudeCellAbsId(clu1,cells);

      fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

      if(IsMatched(trgrelId,relId)) return kTRUE;

    }//end of fired trigger channel loop

  }
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsMatched(Int_t *trgrelid, Int_t *clurelid)
{
  //Returns kTRUE if cluster position coincides with 4x4 position.

  Int_t diffx = trgrelid[2] - clurelid[2];
  Int_t diffz = trgrelid[3] - clurelid[3];

  if(trgrelid[0] != clurelid[0])     return kFALSE; // different modules!//be carefull! STU kindly detects high energy hits on the border beween 2 modules.
  if(diffx < fXmin || fXmax < diffx) return kFALSE; // X-distance too large!
  if(diffz < fZmin || fZmax < diffz) return kFALSE; // Z-distance too large!
  if(fTriggerInputL1 < 0 &&  (WhichTRU(trgrelid[2],trgrelid[3]) != WhichTRU(clurelid[2],clurelid[3]))) return kFALSE;// different TRU in case of L0.

  if(!IsGoodTRUChannel("PHOS",trgrelid[0],trgrelid[2],trgrelid[3])) return kFALSE;

  AliInfo(Form("Accepted : diffx = %d , diffz = %d | fXmin = %d , fZmin = %d , fXmax = %d , fZmax = %d.",diffx,diffz,fXmin,fZmin,fXmax,fZmax));
  return kTRUE;
}
//________________________________________________________________________
Int_t AliPHOSTriggerHelper::FindHighestAmplitudeCellAbsId(AliVCluster *clu, AliVCaloCells *cells)
{
  Int_t cellAbsId=-1;
  Double_t cellamp=0;
  Int_t digMult = clu->GetNCells();
  UShort_t *dlist = clu->GetCellsAbsId();

  Double_t max=0;
  Int_t index=-1;

  for(Int_t iDigit=0;iDigit<digMult;iDigit++){
    cellAbsId = dlist[iDigit];
    cellamp = cells->GetCellAmplitude(cellAbsId) * clu->GetCellAmplitudeFraction(iDigit);

    if(cellamp > max){
      max = cellamp;
      index = cellAbsId;
    }

  }//end of cell loop

  return index;
}
//________________________________________________________________________
Int_t AliPHOSTriggerHelper::WhichTRU(Int_t cellx, Int_t cellz)
{
  Int_t tru = -1;
  if(cellx<1 || 64<cellx){
    AliError("cellx is wrong! tru=-1 will return.");
    return -1;
  }
  if(cellz<1 || 56<cellz){
    AliError("cellz is wrong! tru=-1 will return.");
    return -1;
  }

  Int_t XID = (cellx -1) / 16 + 1;
  Int_t ZID = (cellz -1) / 28;

  tru = 2*XID - ZID;
  //cout << "cellx = " << cellx << " , cellz = " << cellz << " , tru = " << tru << endl;

  return tru;
}
//________________________________________________________________________
Int_t AliPHOSTriggerHelper::WhichTRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ)
{
  //this will return TRU channel 0-111.
  Int_t ch = -1;
  if(cellx<1 || 64<cellx){
    AliError("cellx is wrong! tru=-1 will return.");
    return -1;
  }
  if(cellz<1 || 56<cellz){
    AliError("cellz is wrong! tru=-1 will return.");
    return -1;
  }

  chX = ((cellx -1)/2) %  8;
  chZ = ((cellz -1)/2) % 14;

  ch = 8*chZ + chX;
  //printf("cellx = %d , cellz = %d , chX = %d , chZ = %d , ch = %d.\n",cellx,cellz,chX,chZ,ch);

  return ch;
}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsGoodTRUChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d ",mod)) ;
      return kTRUE ;
    }
    if(!fPHOSTRUBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d",mod)) ;
      return kTRUE ;
    }
    if(fPHOSTRUBadMap[mod]->GetBinContent(ix,iz)>0)
      return kFALSE ;
    else
      return kTRUE ;
  }
  else{
    AliError(Form("Can not find bad channels for detector %s ",det)) ;
  }
  return kTRUE ;
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
