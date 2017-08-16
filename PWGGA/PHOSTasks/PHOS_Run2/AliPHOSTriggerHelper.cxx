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
#include <AliPHOSTriggerHelper.h>
#include <AliCaloPhoton.h>

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
  fCaloTrigger(0x0)
{
  //Constructor
  
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

}
//________________________________________________________________________
AliPHOSTriggerHelper::AliPHOSTriggerHelper(Int_t inputL1, Int_t inputL0):
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
  fCaloTrigger(0x0)
{
  //Constructor
  
  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

  fTriggerInputL1 = inputL1;
  fTriggerInputL0 = inputL0;

  if(fTriggerInputL1 > 0 && fTriggerInputL0 > 0){
    AliError("Both L1 and L0 are selected. Analyzer must select either L1 or L0. L1 has higher priority in this class.");
  }
  AliInfo(Form("Your choice L1:%d , L0:%d",fTriggerInputL1,fTriggerInputL0));

}
//________________________________________________________________________
AliPHOSTriggerHelper::~AliPHOSTriggerHelper()
{

  for(Int_t i=0;i<6;i++){
    if(fPHOSTRUBadMap[i]) delete fPHOSTRUBadMap[i];
  }



}
//________________________________________________________________________
Bool_t AliPHOSTriggerHelper::IsPHI7(AliVEvent *event)
{

  fEvent    = dynamic_cast<AliVEvent*>(event);
  fESDEvent = dynamic_cast<AliESDEvent*>(event);
  fAODEvent = dynamic_cast<AliAODEvent*>(event);
  Int_t run = fEvent->GetRunNumber();

  //if L1 trigger is used, L1 has higher priority
  Bool_t IsPHI7fired = kFALSE;

  if(fTriggerInputL1 > 0)
    IsPHI7fired = event->GetHeader()->GetL1TriggerInputs() & 1 << fTriggerInputL1 - 1;//trigger input -1
  else if(fTriggerInputL0 > 0)
    IsPHI7fired = event->GetHeader()->GetL0TriggerInputs() & 1 << fTriggerInputL0 - 1;//trigger input -1
  else{
    AliInfo("Which PHOS triggered data do you analyze? Please set trigger input (L0 or L1).");
    return kFALSE;
  }

  if(!IsPHI7fired){
    if(fTriggerInputL1 > 0)      AliInfo(Form("Your choice of L1 trigger input %d did not fire.",fTriggerInputL1));
    else if(fTriggerInputL0 > 0) AliInfo(Form("Your choice of L0 trigger input %d did not fire.",fTriggerInputL0));
    return kFALSE;
  }

  if(fTriggerInputL1 > 0)      AliInfo(Form("Your choice of L1 trigger input %d fired.",fTriggerInputL1));
  else if(fTriggerInputL0 > 0) AliInfo(Form("Your choice of L0 trigger input %d fired.",fTriggerInputL0));

  TString trigClasses = event->GetFiredTriggerClasses();
  if(245917 <= run && run <= 246994){//LHC15o
    if(!fIsMC && fTriggerInputL1 == 7 && !trigClasses.Contains("CINT7PHH")) return kFALSE;
    if(!fIsMC && fTriggerInputL1 == 6 && !trigClasses.Contains("CPER7PHM")) return kFALSE;
  }

  if(fTriggerInputL1 > 0) AliInfo(Form("Your choice of L1 trigger input %d matches with fired trigger classes : %s",fTriggerInputL1,trigClasses.Data()));

  if(run<209122) //Run1
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");
  else
    fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

  fCaloTrigger = (AliVCaloTrigger*)event->GetCaloTrigger("PHOS"); 
  fCaloTrigger->Reset();

  const Int_t Nfired = fCaloTrigger->GetEntries();
  AliInfo(Form("%d TRU patches fired in PHOS.",Nfired));

  Int_t multClust = 0;
  TClonesArray *array = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(array){
    AliInfo("PHOSClusterArray is found!");
    multClust = array->GetEntriesFast();
  }
  else{
    AliInfo("PHOSClusterArray is NOT found! clusters from AliVEvent will be used.");
    multClust = fEvent->GetNumberOfCaloClusters();
  }
  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());

  Int_t trgrelId[4]={};
  Int_t tmod=0; //"Online" module number
  Int_t trgabsId=0; //bottom-left 4x4 edge cell absId [1-64], [1-56]
  Int_t trgmodule=0,trgcellx=0,trgcellz=0;
  Int_t maxId=0;
  Int_t relId[4]={};
  Int_t module=0,cellx=0,cellz=0;
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
      AliInfo(Form("This patch fired as %d, which does not match with your choice %d.",fCaloTrigger->GetL1TimeSum(),L1));
      continue;
    }
    AliInfo(Form("This patch fired as %d",fCaloTrigger->GetL1TimeSum()));
    fCaloTrigger->GetPosition(tmod,trgabsId);//tmod is online module numbering. i.e., tmod=1 means half module.

    //cout << "tmod = "<< tmod << " , trgabsId = " << trgabsId << endl;

    fPHOSGeo->AbsToRelNumbering(trgabsId,trgrelId);
    //for offline numbering, relId should be used.

    trgmodule = trgrelId[0];
    trgcellx  = trgrelId[2];
    trgcellz  = trgrelId[3];

    if(trgmodule == 4) continue;

    for(Int_t i=0;i<multClust;i++){
      if(array){
        AliCaloPhoton *ph = (AliCaloPhoton*)array->At(i);
        clu1 = (AliVCluster*)ph->GetCluster();
      }
      else{
        clu1 = (AliVCluster*)fEvent->GetCaloCluster(i);
        if(clu1->GetType() != AliVCluster::kPHOSNeutral
        || clu1->E() < 0.1
//        || clu1->GetNCells() < 2
          ) continue;
      }

      Int_t maxAbsId = FindHighestAmplitudeCellAbsId(clu1,cells);

      fPHOSGeo->AbsToRelNumbering(maxAbsId,relId);
      module = relId[0];
      cellx  = relId[2];
      cellz  = relId[3];

      if(module == 4) continue;

      if(IsMatched(trgrelId,relId) && IsGoodTRUChannel("PHOS",trgrelId[0],trgrelId[2],trgrelId[3])) return kTRUE;
      //if(IsMatched(trgrelId,relId)) return kTRUE;

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

  if(trgrelid[0] != clurelid[0])     return kFALSE; // different modules!
  if(diffx < fXmin || fXmax < diffx) return kFALSE; // X-distance too large!
  if(diffz < fZmin || fZmax < diffz) return kFALSE; // Z-distance too large!

  AliInfo(Form("Accepted : diffx = %d , diffz = %d | fXmin = %d , fZmin = %d , fXmax = %d , fZmax = %d.\n",diffx,diffz,fXmin,fZmin,fXmax,fZmax));

  //if(WhichTRU(trgrelid[2],trgrelid[3]) != WhichTRU(clurelid[2],clurelid[3])) return kFALSE;// different TRU.

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
