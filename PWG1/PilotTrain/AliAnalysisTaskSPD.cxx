#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMultiplicity.h"
#include "AliITSRecPoint.h"
#include "AliITSDetTypeRec.h"
#include "AliGeomManager.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliMagF.h"
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TGeoManager.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskSPD.h"
#include "AliITSsegmentationSPD.h"

ClassImp(AliAnalysisTaskSPD)
  //______________________________________________________________________________
  const Int_t AliAnalysisTaskSPD::fgkDDLModuleMap[20][12] = {
    { 4, 5, 0, 1, 80, 81, 84, 85, 88, 89, 92, 93},
    {12,13, 8, 9, 96, 97,100,101,104,105,108,109},
    {20,21,16,17,112,113,116,117,120,121,124,125},
    {28,29,24,25,128,129,132,133,136,137,140,141},
    {36,37,32,33,144,145,148,149,152,153,156,157},
    {44,45,40,41,160,161,164,165,168,169,172,173},
    {52,53,48,49,176,177,180,181,184,185,188,189},
    {60,61,56,57,192,193,196,197,200,201,204,205},
    {68,69,64,65,208,209,212,213,216,217,220,221},
    {76,77,72,73,224,225,228,229,232,233,236,237},
    { 7, 6, 3, 2, 83, 82, 87, 86, 91, 90, 95, 94},
    {15,14,11,10, 99, 98,103,102,107,106,111,110},
    {23,22,19,18,115,114,119,118,123,122,127,126},
    {31,30,27,26,131,130,135,134,139,138,143,142},
    {39,38,35,34,147,146,151,150,155,154,159,158},
    {47,46,43,42,163,162,167,166,171,170,175,174},
    {55,54,51,50,179,178,183,182,187,186,191,190},
    {63,62,59,58,195,194,199,198,203,202,207,206},
    {71,70,67,66,211,210,215,214,219,218,223,222},
    {79,78,75,74,227,226,231,230,235,234,239,238}
  };

//______________________________________________________________________________
AliAnalysisTaskSPD::AliAnalysisTaskSPD() : 
  AliAnalysisTaskSE(), 
  fOutput(0),
  fSegSPD(0)
{
  //
}
//______________________________________________________________________________
AliAnalysisTaskSPD::AliAnalysisTaskSPD(const char *name) : 
  AliAnalysisTaskSE(name), 
  fOutput(0),
  fSegSPD(0)
{
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskSPD::~AliAnalysisTaskSPD(){
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fSegSPD) delete fSegSPD;
}

//________________________________________________________________________
void AliAnalysisTaskSPD::Init()
{
  // Initialization

  if(fDebug > 1) printf("AliAnalysisTaskSPD::Init() \n");
}
//___________________________________________________________________________

void AliAnalysisTaskSPD::UserCreateOutputObjects() {
  //
  if(fDebug > 1) printf("AliAnalysisTaskSPD::UserCreateOutputObjects() \n");

  fSegSPD = new AliITSsegmentationSPD();

  OpenFile(1);
  fOutput = new TList();

  TH1I *nentries = new TH1I("events","events",1,0,1);
  fOutput->AddLast(nentries);

  TH2F *modulemap[240];

  for(Int_t imod =0; imod < 240; imod++){
    modulemap[imod] = new TH2F(Form("mod%i",imod),Form("cluster map for module %i",imod),800,-4,4,200,-1,1);
    modulemap[imod]->SetXTitle("Local Z (cm)");
    modulemap[imod]->SetYTitle("Local X (cm)");
    fOutput->AddLast(modulemap[imod]);
  }

  TH1F *hFOgood = new TH1F("hFOgood"," Events with FO and clusters ",1200,0,1200);
  hFOgood->SetXTitle("chipkey");
  fOutput->AddLast(hFOgood);

  TH1F *hFOnoisy = new TH1F("hFOnoisy"," Events with FO but no cluster",1200,0,1200);
  hFOnoisy->SetXTitle("chipkey");
  fOutput->AddLast(hFOnoisy);
  
  TH1F *hFiredChips = new TH1F("hFiredChips"," yield of fired chips",1200,0,1200);
  hFiredChips->SetXTitle("chipkey");
  fOutput->AddLast(hFiredChips);

  TH2F *hSPDphivsSPDeta= new TH2F("hSPDphivsSPDeta", "Tracklets - #varphi vs #eta",120,-3.,3,360,0.,2*TMath::Pi());
  hSPDphivsSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hSPDphivsSPDeta->GetYaxis()->SetTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphivsSPDeta);

  TH1F *hSPDphiZpos= new TH1F("hSPDphiZpos", "Tracklets - #varphi (Z>0)",360,0.,2*TMath::Pi());
  hSPDphiZpos->SetXTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphiZpos);

  TH1F *hSPDphiZneg= new TH1F("fHistSPDphivsSPDetaZneg", "Tracklets - #varphi (Z<0)",360,0.,2*TMath::Pi());
  hSPDphiZneg->SetXTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphiZneg);

}
//______________________________________________________________________________
void AliAnalysisTaskSPD::UserExec(Option_t */*option*/)
{
  if(fDebug > 1) printf("AliAnalysisTaskSPD::Exec() \n");
  

  AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!hand) {
    printf("No AliESDInputHandlerRP \n");
    return;
  }  

  AliESDEvent *ESD = hand->GetEvent();
  if(!ESD) {
    printf("No AliESDEvent \n");
    return;
  }  

  TTree * treeRP = hand->GetTreeR("ITS");
  if(!treeRP) {
    printf("No RecPoints tree \n");
    return;
  }  

  // RecPoints info

  TClonesArray statITSrec("AliITSRecPoint");
  TClonesArray *ITSCluster = &statITSrec;
  TBranch* branch=treeRP->GetBranch("ITSRecPoints");
  if(!branch) return;

  branch->SetAddress(&ITSCluster);

  Bool_t isfiredchip[20][6][10];
  for(Int_t eq=0; eq<20; eq++){
    for(Int_t hs =0; hs<6; hs++){
      for(Int_t chip=0; chip<10; chip++){
	isfiredchip[eq][hs][chip]=kFALSE;
      }
    }
  }


  for(Int_t iMod=0;iMod<240;iMod++){
    branch->GetEvent(iMod);
    Int_t nrecp = statITSrec.GetEntries();
    for(Int_t irec=0;irec<nrecp;irec++) {
      AliITSRecPoint *recp = (AliITSRecPoint*)statITSrec.At(irec);
      Int_t lay=recp->GetLayer();
      if(lay>1) continue;
      Float_t local[3]={-1,-1};
      local[1]=recp->GetDetLocalX();
      local[0]=recp->GetDetLocalZ();
      //printf("local X %f   local Z %f   in module %i \n",local[0],local[1],iMod);
      ((TH2F*)fOutput->At(iMod+1))->Fill(local[0],local[1]);
      Int_t eq = GetOnlineEqIdFromOffline(iMod);
      Int_t hs = GetOnlineHSFromOffline(iMod);
      Int_t row, col;
      fSegSPD->LocalToDet(0.5,local[0],row,col);
      Int_t chip = GetOnlineChipFromOffline(iMod,col);
      isfiredchip[eq][hs][chip]=kTRUE;
    }
  }

  //  ESDs info

  // First looking at the FO bits
  const AliMultiplicity *mult = ESD->GetMultiplicity();

  for(Int_t eq=0; eq<20; eq++){
    for(Int_t hs =0; hs<6; hs++){
      for(Int_t chip=0; chip<10; chip++){
	Int_t key = GetOfflineChipKeyFromOnline(eq,hs,chip);
	if(mult->TestFastOrFiredChips(key) && isfiredchip[eq][hs][chip]) ((TH1F*)fOutput->At(241))->Fill(key);
	if(!isfiredchip[eq][hs][chip] && mult->TestFastOrFiredChips(key)) ((TH1F*)fOutput->At(242))->Fill(key);
	if(isfiredchip[eq][hs][chip]) ((TH1F*)fOutput->At(243))->Fill(key);
      }
    }
  }

  // Then looking at tracklets 
  Bool_t eventWithVertex = kFALSE;
  const AliESDVertex* vtxESD = ESD->GetVertex();
  if(vtxESD){

    Double_t esdvtx[3];
    vtxESD->GetXYZ(esdvtx);

    for(Int_t iTracklet =0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++){

      Float_t phiTr= mult->GetPhi(iTracklet);
      Float_t etaTr =mult->GetEta(iTracklet);

      ((TH2F*)fOutput->At(244))->Fill(etaTr,phiTr);

      //layer 0
      //Float_t z = esdvtx[2] + 3.9 / TMath::Tan(2 * TMath::ATan(TMath::Exp(- etaTr)));
      Int_t isZpos = (esdvtx[2]>0 && etaTr>0) ;
      Int_t isZneg = (esdvtx[2]<0 && etaTr<0) ; 

      if(isZpos)  ((TH1F*)fOutput->At(245))->Fill(phiTr);
      if(isZneg)  ((TH1F*)fOutput->At(246))->Fill(phiTr);
    }
  }

  ((TH1I*)fOutput->At(0))->Fill(0); // monitoring plot

  PostData(1,fOutput);  
  return;
  
}
//______________________________________________________________________________
void AliAnalysisTaskSPD::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  if(fDebug > 1) printf("AliAnalysisTaskSPD:: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    //return;
  }
}
//________________________________________________________________________
UInt_t AliAnalysisTaskSPD::GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (module)
  if (eqId<20 && hs<6 && chip<10) return fgkDDLModuleMap[eqId][hs*2+chip/5];
  else return 240;
}
//________________________________________________________________________
UInt_t AliAnalysisTaskSPD::GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // online->offline (chip key: 0-1199)
  if (eqId<20 && hs<6 && chip<10) {
    UInt_t module = GetOfflineModuleFromOnline(eqId,hs,chip);
    UInt_t chipInModule = ( chip>4 ? chip-5 : chip );
    if(eqId>9) chipInModule = 4 - chipInModule;  // side C only
    return (module*5 + chipInModule);
  } else return 1200;
}
//__________________________________________________________________________
UInt_t AliAnalysisTaskSPD::GetOnlineEqIdFromOffline(UInt_t module) {
  // offline->online (eq)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return eqId;
    }
  }
  return 20; // error
}
//__________________________________________________________________________
UInt_t AliAnalysisTaskSPD::GetOnlineHSFromOffline(UInt_t module) {
  // offline->online (hs)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eqId,iModule)==(Int_t)module) return iModule/2;
    }
  }
  return 6; // error
}
//__________________________________________________________________________
UInt_t AliAnalysisTaskSPD::GetOnlineChipFromOffline(UInt_t module, UInt_t colM) {
  // offline->online (chip)
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      if (GetModuleNumber(eq,iModule)==(Int_t)module) {
        if (module<80) {
          if (eq<10) { // side A
            return (159-colM)/32 + 5*(iModule%2);
          }
          else { // side C
            return colM/32 + 5*(iModule%2);
          }
        }
        else if (module<240) {
          if (eq<10) { // side A
            return colM/32 + 5*(iModule%2);
          }
          else { // side C
            return (159-colM)/32 + 5*(iModule%2);
          }
        }
      }
    }
  }
  return 10; // error
}
//__________________________________________________________________________
Int_t AliAnalysisTaskSPD::GetModuleNumber(UInt_t iDDL, UInt_t iModule) {
  if (iDDL<20 && iModule<12) return fgkDDLModuleMap[iDDL][iModule];
  else return 240;
}


