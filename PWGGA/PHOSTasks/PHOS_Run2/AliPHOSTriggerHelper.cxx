#include <vector>

#include <TObject.h>
#include <TClonesArray.h>
#include <TH2.h>
#include <TVector3.h>
#include <TMath.h>

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

using namespace std;

// Author: Daiki Sekihata (Hiroshima University)
ClassImp(AliPHOSTriggerHelper)

//________________________________________________________________________
AliPHOSTriggerHelper::AliPHOSTriggerHelper():
  fPHOSGeo(0x0),
  fXmin(-3),
  fXmax(3),
  fZmin(-3),
  fZmax(3),
  fMatchingDeltaR(0.020),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
	fTriggerInputL1(-1),//5:1PHL, 6:1PHM, 7:1PHH
	fTriggerInputL0(-1),//9:0PH0
  fIsMC(kFALSE),
  fCaloTrigger(0x0),
  fIsUserTRUBadMap(kFALSE),
  fRunNumber(-1),
  fUseDeltaRMatching(kFALSE),
  fApplyTOFCut(kFALSE),
  fDRN(-1)
{
  //Constructor
  
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

}
//________________________________________________________________________
AliPHOSTriggerHelper::AliPHOSTriggerHelper(TString trigger, Bool_t isMC):
  fPHOSGeo(0x0),
  fXmin(-3),
  fXmax(3),
  fZmin(-3),
  fZmax(3),
  fMatchingDeltaR(0.020),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
	fTriggerInputL1(-1),//5:1PHL, 6:1PHM, 7:1PHH
	fTriggerInputL0(-1),//9:0PH0
  fIsMC(kFALSE),
  fCaloTrigger(0x0),
  fIsUserTRUBadMap(kFALSE),
  fRunNumber(-1),
  fUseDeltaRMatching(kFALSE),
  fApplyTOFCut(kFALSE),
  fDRN(-1)
{
  //Constructor
   
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

  fIsMC = isMC;

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
AliPHOSTriggerHelper::AliPHOSTriggerHelper(Int_t L1triggerinput, Int_t L0triggerinput, Bool_t isMC):
  fPHOSGeo(0x0),
  fXmin(-3),
  fXmax(3),
  fZmin(-3),
  fZmax(3),
  fMatchingDeltaR(0.020),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
	fTriggerInputL1(-1),//5:1PHL, 6:1PHM, 7:1PHH
	fTriggerInputL0(-1),//9:0PH0
  fIsMC(kFALSE),
  fCaloTrigger(0x0),
  fIsUserTRUBadMap(kFALSE),
  fRunNumber(-1),
  fUseDeltaRMatching(kFALSE),
  fApplyTOFCut(kFALSE),
  fDRN(-1)
{
  //Constructor
   
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

  fIsMC = isMC;

  if(L1triggerinput > 0){
    fTriggerInputL1 = L1triggerinput;
    fTriggerInputL0 = -1;
    //STU stores top-left fired channel
    fXmin = -3;
    fXmax = 0;
    fZmin = -1;
    fZmax = 2;
  }
  else if(L0triggerinput > 0){
    fTriggerInputL1 = -1;
    fTriggerInputL0 = L0triggerinput;
    //TRU stores bottom-left fired channel.
    fXmin = -3;
    fXmax = 0;
    fZmin = -3;
    fZmax = 0;
  }

  if(fTriggerInputL1 > 0 && fTriggerInputL0 > 0){
    AliError("Both L1 and L0 are selected. Analyzer must select either L1 or L0. L1 has higher priority in this class.");
  }

  if(fTriggerInputL1 < 0 && fTriggerInputL0 < 0){
    AliError("Neither L1 nor L0 are selected. Analyzer must select either L1 or L0. L1 has higher priority in this class.");
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
Bool_t AliPHOSTriggerHelper::IsPHI7(AliVEvent *event, AliPHOSClusterCuts *cuts, Double_t Emin, Double_t ETrigger, Bool_t isCoreUsed)
{
  fEvent    = dynamic_cast<AliVEvent*>(event);
  fESDEvent = dynamic_cast<AliESDEvent*>(event);
  fAODEvent = dynamic_cast<AliAODEvent*>(event);
  Int_t run = fEvent->GetRunNumber();

  if(fIsMC && fDRN > 0){
    run = fDRN;
    AliInfo(Form("A dummy run number is set. run number = %d",run));
  }

  if(run<209122) //Run1
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");
  else
    fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

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

  fCaloTrigger = (AliVCaloTrigger*)event->GetCaloTrigger("PHOS"); 
  fCaloTrigger->Reset();

  if(fIsMC){
    //please call this line after AliVCaloTrigger and PHOSGeometry and  TRU bad maps are prepared.
      AliInfo("This is MC analysis. Always accept event.");
     return kTRUE;
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

  //AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());//only for maxabsid

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t relId[4]={};
  Float_t position[3] = {};

  Int_t L1=-999;
  if(fTriggerInputL1 > 0){
    if     (fTriggerInputL1 == 7) L1 = 0;//L1 high
    else if(fTriggerInputL1 == 6) L1 = 1;//L1 medium
    else if(fTriggerInputL1 == 5) L1 = 2;//L1 low
  }
  else if(fTriggerInputL0 > 0){
    L1 = -1;//L0
  }

  Double_t energy = 0;

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
      if(fApplyTOFCut && !ph->IsTOFOK()) continue;

      energy = ph->Energy();
      if(isCoreUsed) energy = (ph->GetMomV2())->Energy();

      //printf("cluster  E = %e GeV\n",energy);
      if(energy < Emin) continue;
      if(energy < ETrigger) continue;

      //AliVCluster *clu1 = (AliVCluster*)ph->GetCluster();//only for maxabsid
      //Int_t maxAbsId = FindHighestAmplitudeCellAbsId(clu1,cells);
      //fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);

      position[0] = ph->EMCx();
      position[1] = ph->EMCy();
      position[2] = ph->EMCz();

      relId[0] = 0; relId[1] = 0; relId[2] = 0; relId[3] = 0;
      TVector3 global1(position);
      fPHOSGeo->GlobalPos2RelId(global1,relId);

      //printf("cluster positoin M%d X%d Z%d and E = %e GeV\n",relId[0],relId[2],relId[3],energy);

      if(fUseDeltaRMatching){
        if(IsMatchedDeltaR(trgrelId,global1)) return kTRUE;
      }
      else{
        if(IsMatched(trgrelId,relId)) return kTRUE;
      }
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

  if(fTriggerInputL1 > 0){//for L1
    Int_t module = clurelid[0];//offline numbering
    if(module == 1){//M1
      //because of hardware issue, fired position on M1 is shifted by -2 in Z.
      if(diffx < fXmin || fXmax < diffx) return kFALSE; // X-distance too large!
      if(diffz < fZmin-2 || fZmax-2 < diffz) return kFALSE; // Z-distance too large!
    }  
    else{//M2,3,4
      if(diffx < fXmin || fXmax < diffx) return kFALSE; // X-distance too large!
      if(diffz < fZmin || fZmax < diffz) return kFALSE; // Z-distance too large!
    }
  }
  else if(fTriggerInputL0 > 0){//for L0
    if(diffx < fXmin || fXmax < diffx) return kFALSE; // X-distance too large!
    if(diffz < fZmin || fZmax < diffz) return kFALSE; // Z-distance too large!
  }

  //if(fTriggerInputL1 < 0 && (WhichTRU(trgrelid[2],trgrelid[3]) != WhichTRU(clurelid[2],clurelid[3]))) return kFALSE;// different TRU in case of L0.//not needed because L0 can detect very high energy photon at the border.

  if(!IsGoodTRUChannel("PHOS",trgrelid[0],trgrelid[2],trgrelid[3])) return kFALSE;

  AliInfo(Form("Accepted : diffx = %d , diffz = %d | fXmin = %d , fZmin = %d , fXmax = %d , fZmax = %d.",diffx,diffz,fXmin,fZmin,fXmax,fZmax));
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsMatchedDeltaR(Int_t *trgrelid, TVector3 position)
{
  //Returns kTRUE if cluster position (center of gravity) is close to fired TRU channel.
  Float_t x = 0;
  Float_t y = 0;
  fPHOSGeo->RelPosInModule(trgrelid,x,y); 
  TVector3 trgglobal;
  fPHOSGeo->Local2Global(trgrelid[0],x,y,trgglobal);
  Double_t eta_trigger = trgglobal.Eta();
  Double_t phi_trigger = trgglobal.Phi();
  if(phi_trigger < 0) phi_trigger += TMath::TwoPi();
 
  Int_t clurelid[4] = {}; 
  fPHOSGeo->GlobalPos2RelId(position,clurelid);
  Double_t eta_cluster = position.Eta();
  Double_t phi_cluster = position.Phi();
  if(phi_cluster < 0) phi_cluster += TMath::TwoPi();

  Double_t DeltaR = TMath::Sqrt( TMath::Power(eta_trigger - eta_cluster,2) +  TMath::Power(phi_trigger - phi_cluster,2) );
  if(DeltaR > fMatchingDeltaR) return kFALSE;

  if(trgrelid[0] != clurelid[0])     return kFALSE; // different modules!//be carefull! STU kindly detects high energy hits on the border beween 2 modules.
  if(fTriggerInputL1 < 0 &&  (WhichTRU(trgrelid[2],trgrelid[3]) != WhichTRU(clurelid[2],clurelid[3]))) return kFALSE;// different TRU in case of L0.

  if(!IsGoodTRUChannel("PHOS",trgrelid[0],trgrelid[2],trgrelid[3])) return kFALSE;

  AliInfo(Form("Accepted : DeltaR = %e | Matching criterion is less than %4.3f",DeltaR,fMatchingDeltaR));
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

  if(ix < 0 || iz < 0) return kFALSE;

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
    return kFALSE ;
  }
  return kTRUE ;
}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsOnActiveTRUChannel(AliCaloPhoton *ph)
{
  //Check if this channel belogs to the good ones
  //this function is used to check cluster hit position is on active TRU channel,
  //and to restrict hit acceptance where both FEE and TRU channels are active in trigger analysis
  //matched or non-matched do not matter here, but IsGoodTRUChannel is called for the convinience.

  Float_t position[3] = {};
  position[0] = ph->EMCx();
  position[1] = ph->EMCy();
  position[2] = ph->EMCz();

  Int_t relId[4] = {};
  TVector3 global1(position);
  fPHOSGeo->GlobalPos2RelId(global1,relId);
  
  Int_t module = relId[0];
  Int_t cellx  = relId[2];
  Int_t cellz  = relId[3];

  //cellx in [1,64]
  //cellz in [1,56]

  //convert (cellx,cellz) in FEE to (TRUchX,TRUchZ) in TRU.

  Int_t cellx00 = cellx;
  Int_t cellz00 = cellz;

  Int_t cellx01 = cellx-2;
  Int_t cellz01 = cellz;

  Int_t cellx10 = cellx;
  Int_t cellz10 = cellz-2;

  Int_t cellx11 = cellx-2;
  Int_t cellz11 = cellz-2;

  if(fTriggerInputL1 > 0){//L1 trigger analysis
    //STU stores fired position at top-left
    //L1 can detect a high energy cluster at a border of TRUs.
    if(cellx %2 == 0){
      cellx00 -= 1;
      cellx01 -= 1;
      cellx10 -= 1;
      cellx11 -= 1;
    }
    if(cellz %2 == 0){
      cellz00 += 1;
      cellz01 += 1;
      cellz10 += 1;
      cellz11 += 1;
    }
    return IsGoodTRUChannel("PHOS",module,cellx00,cellz00) | IsGoodTRUChannel("PHOS",module,cellx01,cellz01) | IsGoodTRUChannel("PHOS",module,cellx10,cellz10) | IsGoodTRUChannel("PHOS",module,cellx11,cellz11);
  }
  else if(fTriggerInputL0 >0){//L0 trigger analysis
    //TRU stores fired position at bottom-left
    //note that 2D TRU bad maps are filled in odd number bins.(1,1) (1,3) ... (55,53) (55,55), ... (63,55)
    //At maximum, 1 cluster can fire 4 TRU channels.
    //TRU can not detect a high energy cluser at a border of TRUs.
    if(cellx %2 == 0){
      cellx00 -= 1;
      cellx01 -= 1;
      cellx10 -= 1;
      cellx11 -= 1;
    }
    if(cellz %2 == 0){
      cellz00 -= 1;
      cellz01 -= 1;
      cellz10 -= 1;
      cellz11 -= 1;
    }

    if(cellx % 16 <= 2 && cellz % 28 <= 2)
      return IsGoodTRUChannel("PHOS",module,cellx00,cellz00);
    else if(cellx % 16 <= 2)
      return IsGoodTRUChannel("PHOS",module,cellx10,cellz10) | IsGoodTRUChannel("PHOS",module,cellx00,cellz00);
    else if(cellz % 28 <= 2)
      return IsGoodTRUChannel("PHOS",module,cellx00,cellz00) | IsGoodTRUChannel("PHOS",module,cellx01,cellz01);
    else
      return IsGoodTRUChannel("PHOS",module,cellx00,cellz00) | IsGoodTRUChannel("PHOS",module,cellx01,cellz01) | IsGoodTRUChannel("PHOS",module,cellx10,cellz10) | IsGoodTRUChannel("PHOS",module,cellx11,cellz11);

  }
  else{
    AliInfo("Trigger input is neither L1 nor L0. Check your configuration. return kFALSE");
    return kFALSE;
  }

}
//________________________________________________________________________
Double_t AliPHOSTriggerHelper::GetDistanceToClosestTRUChannel(AliCaloPhoton *ph)
{
  Float_t position[3] = {};
  position[0] = ph->EMCx();
  position[1] = ph->EMCy();
  position[2] = ph->EMCz();
  TVector3 global1(position); 
  Double_t eta_cluster = global1.Eta();
  Double_t phi_cluster = global1.Phi();
  if(phi_cluster < 0) phi_cluster += TMath::TwoPi();

  Int_t trgrelid[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]

  vector<Double_t> vR;

  Int_t L1=-999;
  if(fTriggerInputL1 > 0){
    if     (fTriggerInputL1 == 7) L1 = 0;//L1 high
    else if(fTriggerInputL1 == 6) L1 = 1;//L1 medium
    else if(fTriggerInputL1 == 5) L1 = 2;//L1 low
  }
  else if(fTriggerInputL0 > 0){
    L1 = -1;//L0
  }

  fCaloTrigger->Reset();

  while(fCaloTrigger->Next()){
    // L1 threshold: -1-L0, 0-PHH, 1-PHM, 2-PHL

    if(fCaloTrigger->GetL1TimeSum() != L1){
      //AliInfo(Form("This patch fired as %d, which does not match with your choice %d.",fCaloTrigger->GetL1TimeSum(),L1));
      continue;
    }
    fCaloTrigger->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

    fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelid);//for offline numbering, relId should be used.
    if(!IsGoodTRUChannel("PHOS",trgrelid[0],trgrelid[2],trgrelid[3])) continue;

    Float_t x = 0;
    Float_t y = 0;
    fPHOSGeo->RelPosInModule(trgrelid,x,y); 
    TVector3 trgglobal;
    fPHOSGeo->Local2Global(trgrelid[0],x,y,trgglobal);
    Double_t eta_trigger = trgglobal.Eta();
    Double_t phi_trigger = trgglobal.Phi();
    if(phi_trigger < 0) phi_trigger += TMath::TwoPi();

    Double_t DeltaR = TMath::Sqrt( TMath::Power(eta_trigger - eta_cluster,2) +  TMath::Power(phi_trigger - phi_cluster,2) );
    vR.push_back(DeltaR);

  }//end of trigger loop

  if(vR.size() != 0){
    Double_t min = *min_element(vR.begin(),vR.end());
    return min;
  }
  else return 999;
}
//________________________________________________________________________
//________________________________________________________________________
