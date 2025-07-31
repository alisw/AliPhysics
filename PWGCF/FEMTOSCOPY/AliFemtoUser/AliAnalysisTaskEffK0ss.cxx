#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliAODv0.h"
#include "AliAODpidUtil.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskEffK0ss.h"

ClassImp(AliAnalysisTaskEffK0ss)

double fV1K0s[3];


void AliAnalysisTaskEffK0ss::SetFB(int fb)
{
  fFB = fb;
}

void AliAnalysisTaskEffK0ss::SetPidMethod(PidMethod method)
{
  fPidMethod = method;
}

int AliAnalysisTaskEffK0ss::GetPidMethod()
{
  return (int)fPidMethod;
}

void AliAnalysisTaskEffK0ss::SetPidMethod(int method)
{
  switch(method){
  case 0: fPidMethod=kNSigma;
    break;
  case 1: fPidMethod=kNSigmaNoDoubleCounting;
    break;
  case 2: fPidMethod=kExclusivePID;
    break;
  case 3: fPidMethod=kExclusivePIDDiffRejection;
    break;
  }

}
void AliAnalysisTaskEffK0ss::SetNsigmaDaughters(Double_t nsigmaDaught)
{
  fNsigmaDaughters = nsigmaDaught;
}

void AliAnalysisTaskEffK0ss::SetMinCosPointingAngle(Double_t cosPA)
{
  fMinCosPointingAngle = cosPA;
}

void AliAnalysisTaskEffK0ss::SetInvMassWindow(Double_t window)
{
  fInvMassWindow = window;
}

void AliAnalysisTaskEffK0ss::SetMaxDCADaughters(Double_t dcaDaughtersMax)
{
  fMaxDcaV0Daughters = dcaDaughtersMax;
}

void AliAnalysisTaskEffK0ss::SetMinDCAToPrimVtx(Double_t dcaPVMin)
{
  fMinDCAToPrimVtx = dcaPVMin;
}

void AliAnalysisTaskEffK0ss::SetElectronRejection(Bool_t electronRejection)
{
  fElectronReject = electronRejection;
}

void AliAnalysisTaskEffK0ss::SetNsigmaElectronRejection(Double_t nsigmaErej)
{
  fNsigmaElectronRejection = nsigmaErej;
}

void AliAnalysisTaskEffK0ss::SetMinDcaPosToPrimVertex(Double_t minDCApos)
{
  fMinDcaPosToPrimVertex = minDCApos;
}

void AliAnalysisTaskEffK0ss::SetMinDcaNegToPrimVertex(Double_t minDCAneg)
{
  fMinDcaNegToPrimVertex = minDCAneg;
}

void AliAnalysisTaskEffK0ss::SetMaxCTauK0s(Double_t maxctau)
{
 fMaxCTauK0s = maxctau;
}


void AliAnalysisTaskEffK0ss::SetMinV0Radius(Double_t minV0rad)
{
  fMinV0Radius = minV0rad;
}

void AliAnalysisTaskEffK0ss::SetUseDCAcuts(Bool_t ownDCA)
{
 fUseDcaCuts = ownDCA;
}


void AliAnalysisTaskEffK0ss::SetUseDCAcutsxy(Float_t dcaxy)
{
 fDcaXYCut = dcaxy;
}

void AliAnalysisTaskEffK0ss::SetUseDCAcutsz(Float_t dcaz)
{
 fDcaZCut = dcaz;
}

//_______________________________________________________
AliAnalysisTaskEffK0ss::AliAnalysisTaskEffK0ss() :
  AliAnalysisTaskSE("name"), centrality(0), fHistoList(0), fMassInvK0sPass(0), fMassInvK0sFail(0), fMassInvK0s(0), fMassInvK0sAfterCuts(0),
fMassInvK0sPt(0), fEtaK0s(0), fPtK0s(0), fCutsK0s(0), fTruePtK0sMC(0), fRecPtK0sMC(0), fDCAtoPrimVtx(0),
fIfAliEventCuts(kTRUE), fFB(128), fPidMethod(kExclusivePIDDiffRejection), fEstEventMult(kV0M),
fpidResponse(0), fAODpidUtil(0), fEventCuts(0), fTrackPileUpRemoval(kFALSE), fV0PileUpRemoval(kFALSE),

// --- New configurable V0 and track cuts ---
fUseDcaCuts(kFALSE), 
fDcaXYCut(0.3),
fDcaZCut(1),

fMaxDcaV0Daughters(1.0),
fMinDcaPosToPrimVertex(0.1),
fMinDcaNegToPrimVertex(0.1),
fMaxCTauK0s(4.0 * 2.68),
fMinV0Radius(5.0),

// --- Existing V0 cut settings ---
fNsigmaDaughters(5.0),
fMinCosPointingAngle(0.998),
fInvMassWindow(0.015),
fMinDCAToPrimVtx(0.1),
fElectronReject(kTRUE),
fNsigmaElectronRejection(2.0),

// --- Centrality and PID ---
fCentMin(0),
fCentMax(100),
fPIDKch(0),
fPIDKeCut(0),
fPVzCut(10)

{
  for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    for(Int_t chg=0;chg<2;chg++){
      fGeneratedMCPrimaries[i][chg] = NULL;
      fMCPrimariesThatAreReconstructed[i][chg] = NULL;
      fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
      fReconstructedAfterCuts[i][chg] = NULL;
      fReconstructedNotPrimaries[i][chg] = NULL;
      fReconstructedPrimaries[i][chg] = NULL;
      fContamination[i][chg] = NULL;

      fPrimVsDCA[i][chg] = NULL;
      fSecWeakVsDCA[i][chg] = NULL;
      fSecMatVsDCA[i][chg] = NULL;
      fFakeVsDCA[i][chg] = NULL;

      fPrimVsCosPointingAngle[i][chg] = NULL;
      fSecWeakVsCosPointingAngle[i][chg] = NULL;
      fSecMatVsCosPointingAngle[i][chg] = NULL;
      fFakeVsCosPointingAngle[i][chg] = NULL;

      fMCPrimariesThatAreReconstructed4D[i][chg] = NULL;
      fGeneratedMCPrimaries4D[i][chg] = NULL;
    }
  }
  for ( Int_t i = 0; i < 11; i++) {
    if(i<4) fHistEv[i] = NULL;
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
  }

  /* init track cuts */
  //if(pidMethod!=-1) SetPidMethod(pidMethod);
  //SetFB(filterbit);



  //DefineInput(0, TChain::Class());
  //DefineOutput(0, TTree::Class());
  //DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskEffK0ss::AliAnalysisTaskEffK0ss(TString name,
                       int pidMethod,
                       int filterbit,
                       Double_t nsigmaDaught,
                       Double_t cosPA,
                       Double_t window,
                       Double_t dcaDaughtersMax,
                       Double_t dcaPVMin,
                       Bool_t electronRejection,
                       Double_t nsigmaErej,
		               Double_t minDCApos, 
		               Double_t minDCAneg,
		               Double_t maxctau, 
		               Double_t minV0rad,
		               Bool_t ownDCA, 
		               Float_t dcaxy, 
		               Float_t dcaz) : AliAnalysisTaskSE(name), centrality(0), fHistoList(0), fMassInvK0sPass(0), fMassInvK0sFail(0), fMassInvK0s(0), fMassInvK0sAfterCuts(0),
fMassInvK0sPt(0), fEtaK0s(0), fPtK0s(0), fCutsK0s(0), fTruePtK0sMC(0), fRecPtK0sMC(0), fDCAtoPrimVtx(0),
fIfAliEventCuts(kTRUE), fFB(128), fPidMethod(kExclusivePIDDiffRejection), fEstEventMult(kV0M),
fpidResponse(0), fAODpidUtil(0), fEventCuts(0), fTrackPileUpRemoval(kFALSE), fV0PileUpRemoval(kFALSE),

// --- New configurable V0 and track cuts ---
fUseDcaCuts(kFALSE), 
fDcaXYCut(0.1),
fDcaZCut(0.2),

fMaxDcaV0Daughters(1.0),
fMinDcaPosToPrimVertex(0.1),
fMinDcaNegToPrimVertex(0.1),
fMaxCTauK0s(4.0 * 2.68),
fMinV0Radius(5.0),

// --- Existing V0 cut settings ---
fNsigmaDaughters(5.0),
fMinCosPointingAngle(0.998),
fInvMassWindow(0.015),
fMinDCAToPrimVtx(0.1),
fElectronReject(kTRUE),
fNsigmaElectronRejection(2.0),

// --- Centrality and PID ---
fCentMin(0),
fCentMax(100),
fPIDKch(0),
fPIDKeCut(0),
fPVzCut(10)

{

  for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    for(Int_t chg=0;chg<2;chg++){
      fGeneratedMCPrimaries[i][chg] = NULL;
      fMCPrimariesThatAreReconstructed[i][chg] = NULL;
      fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
      fReconstructedAfterCuts[i][chg] = NULL;
      fReconstructedNotPrimaries[i][chg] = NULL;
      fReconstructedPrimaries[i][chg] = NULL;
      fContamination[i][chg] = NULL;

      fPrimVsDCA[i][chg] = NULL;
      fSecWeakVsDCA[i][chg] = NULL;
      fSecMatVsDCA[i][chg] = NULL;
      fFakeVsDCA[i][chg] = NULL;

      fPrimVsCosPointingAngle[i][chg] = NULL;
      fSecWeakVsCosPointingAngle[i][chg] = NULL;
      fSecMatVsCosPointingAngle[i][chg] = NULL;
      fFakeVsCosPointingAngle[i][chg] = NULL;

      fMCPrimariesThatAreReconstructed4D[i][chg] = NULL;
      fGeneratedMCPrimaries4D[i][chg] = NULL;
    }
  }
  for ( Int_t i = 0; i < 11; i++) {
    if(i<4) fHistEv[i] = NULL;
    fHistQA[i] = NULL;
    if(i<4) fHistQA2D[i] = NULL;
  }

  /* init track cuts */
  if(pidMethod!=-1) SetPidMethod(pidMethod);
  SetFB(filterbit);



  //DefineInput(0, TChain::Class());
  //DefineOutput(0, TTree::Class());
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskEffK0ss::~AliAnalysisTaskEffK0ss()
{
  //Destructor
  if(AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
    delete fHistoList;
}

//_______________________________________________________

void AliAnalysisTaskEffK0ss::UserCreateOutputObjects()
{

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);

  TString hname1, hname2, hname3, hname4, hname5, hname6, hname7, hname8, hname9, hname10, hname11,hname12,hname13,hname14,hname15,hname16,hname17,hname18,hname19,hname20,hname21,hname22,hname23,hname24,hname25,hname26 ;

  TString htitle1, htitle2, htitle3, htitle4,htitle6,htitle5,htitle7,htitle8,htitle9,htitle10,htitle11,htitle12,htitle13,htitle14,htitle15,htitle16,htitle17,htitle18,htitle19,htitle20,htitle21,htitle22,htitle23,htitle24,htitle25,htitle26;

  TString hname1M, hname2M, hname3M, hname4M, hname5M, hname;

  TString htitle1M, htitle2M, htitle3M, htitle4M, htitle5M, htitle;

  TString parttypename = "None";

  for(Int_t j = 0; j < PARTTYPES; j++)  {
    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";
    else if (j==4) parttypename="K0s";

    for(Int_t i = 0; i < MULTBINS; i++)  {
      hname1  = "hGeneratedMCPrimariesEffM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level eta_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries[i*PARTTYPES+j][0] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname1+="Minus";htitle1+="Minus";
      fGeneratedMCPrimaries[i*PARTTYPES+j][1] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedNoNsigmaM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) no Nsigma cut only PDG M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname2  = "hHistoReconstructedAfterCutsM"; hname2+=i; hname2+=parttypename;
      htitle2 = "Total Reconstructed tracks M"; htitle2+=i; htitle2+=parttypename;
      fReconstructedAfterCuts[i*PARTTYPES+j][0] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname2+="Minus";htitle2+="Minus";
      fReconstructedAfterCuts[i*PARTTYPES+j][1] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname4  = "hHistoReconstructedNotPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (not primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedNotPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedNotPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname4  = "hHistoReconstructedPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,100,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,100,0.,10.0);

      hname5  = "hContaminationM"; hname5+=i; hname5+=parttypename;
      htitle5 = "Contamination M"; htitle5+=i; htitle5+=parttypename;
      fContamination[i*PARTTYPES+j][0] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,50,0.,10.0); //50
      hname5+="Minus";htitle5+="Minus";
      fContamination[i*PARTTYPES+j][1] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,50,0.,10.0); //50

      fReconstructedAfterCuts[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedNotPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries[i*PARTTYPES+j][1]->Sumw2();
      fContamination[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][1]->Sumw2();


      hname6 = "hPrimVsDCAM"; hname6+=i; hname6+=parttypename;
      htitle6 = "Primaries vs DCA M"; htitle6+=i; htitle6+=parttypename;
      fPrimVsDCA[i*PARTTYPES+j][0] = new TH2F(hname6.Data(),htitle6.Data(),500,0,5.0,7,0.5,4);
      hname6+="Minus"; htitle6+="Minus";
      fPrimVsDCA[i*PARTTYPES+j][1] = new TH2F(hname6.Data(),htitle6.Data(),500,0,5.0,7,0.5,4);
      hname7 = "hSecWeakVsDCAM"; hname7+=i; hname7+=parttypename;
      htitle7 = "Sec. weak decay vs DCA M"; htitle7+=i; htitle7+=parttypename;
      fSecWeakVsDCA[i*PARTTYPES+j][0] = new TH2F(hname7.Data(),htitle7.Data(),500,0,5.0,7,0.5,4);
      hname7+="Minus"; htitle7+="Minus";
      fSecWeakVsDCA[i*PARTTYPES+j][1] = new TH2F(hname7.Data(),htitle7.Data(),500,0,5.0,7,0.5,4);
      hname8 = "hSecMatVsDCAM"; hname8+=i; hname8+=parttypename;
      htitle8 = "Sec. material vs DCA M"; htitle8+=i; htitle8+=parttypename;
      fSecMatVsDCA[i*PARTTYPES+j][0] = new TH2F(hname8.Data(),htitle8.Data(),500,0,5.0,7,0.5,4);
      hname8+="Minus"; htitle8+="Minus";
      fSecMatVsDCA[i*PARTTYPES+j][1] = new TH2F(hname8.Data(),htitle8.Data(),500,0,5.0,7,0.5,4);
      hname9 = "hFakeVsDCAM"; hname9+=i; hname9+=parttypename;
      htitle9 = "Fake vs DCA M"; htitle9+=i; htitle9+=parttypename;
      fFakeVsDCA[i*PARTTYPES+j][0] = new TH2F(hname9.Data(),htitle9.Data(),500,0,5.0,7,0.5,4);
      hname9+="Minus"; htitle9+="Minus";
      fFakeVsDCA[i*PARTTYPES+j][1] = new TH2F(hname9.Data(),htitle9.Data(),500,0,5.0,7,0.5,4);


      hname10 = "hPrimVsCosPointingAngle"; hname10+=i; hname10+=parttypename;
      htitle10 = "Primaries vs CosPointingAngle M"; htitle10+=i; htitle10+=parttypename;
      fPrimVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname10.Data(),htitle10.Data(),200,0.95,1.0,7,0.5,4);
      hname10+="Minus"; htitle10+="Minus";
      fPrimVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname10.Data(),htitle10.Data(),200,0.95,1.0,7,0.5,4);
      hname11 = "hSecWeakVsCosPointingAngleM"; hname11+=i; hname11+=parttypename;
      htitle11 = "Sec. weak decay vs CosPointingAngle M"; htitle11+=i; htitle11+=parttypename;
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname11.Data(),htitle11.Data(),200,0.95,1.0,7,0.5,4);
      hname11+="Minus"; htitle11+="Minus";
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname11.Data(),htitle11.Data(),200,0.95,1.0,7,0.5,4);
      hname12 = "hSecMatVsCosPointingAngleM"; hname12+=i; hname12+=parttypename;
      htitle12 = "Sec. material vs CosPointingAngle M"; htitle12+=i; htitle12+=parttypename;
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname12.Data(),htitle12.Data(),200,0.95,1.0,7,0.5,4);
      hname12+="Minus"; htitle12+="Minus";
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname12.Data(),htitle12.Data(),200,0.95,1.0,7,0.5,4);
      hname13 = "hFakeVsCosPointingAngleM"; hname13+=i; hname13+=parttypename;
      htitle13 = "Fake vs CosPointingAngle M"; htitle13+=i; htitle13+=parttypename;
      fFakeVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname13.Data(),htitle13.Data(),200,0.95,1.0,7,0.5,4);
      hname13+="Minus"; htitle13+="Minus";
      fFakeVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname13.Data(),htitle13.Data(),200,0.95,1.0,7,0.5,4);


      hname14 = "hPrimVsDecayRadius"; hname14+=i; hname14+=parttypename;
      htitle14 = "Primaries vs DecayRadius M"; htitle14+=i; htitle14+=parttypename;
      fPrimVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname14.Data(),htitle14.Data(),500,0,5.0,7,0.5,4);
      hname14+="Minus"; htitle14+="Minus";
      fPrimVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname14.Data(),htitle14.Data(),500,0,5.0,7,0.5,4);
      hname15 = "hSecWeakVsDecayRadiusM"; hname15+=i; hname15+=parttypename;
      htitle15 = "Sec. weak decay vs DecayRadius M"; htitle15+=i; htitle15+=parttypename;
      fSecWeakVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname15.Data(),htitle15.Data(),500,0,5.0,7,0.5,4);
      hname15+="Minus"; htitle15+="Minus";
      fSecWeakVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname15.Data(),htitle15.Data(),500,0,5.0,7,0.5,4);
      hname16 = "hSecMatVsDecayRadiusM"; hname16+=i; hname16+=parttypename;
      htitle16 = "Sec. material vs DecayRadius M"; htitle16+=i; htitle16+=parttypename;
      fSecMatVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname16.Data(),htitle16.Data(),500,0,5.0,7,0.5,4);
      hname16+="Minus"; htitle16+="Minus";
      fSecMatVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname16.Data(),htitle16.Data(),500,0,5.0,7,0.5,4);
      hname17 = "hFakeVsDecayRadiusM"; hname17+=i; hname17+=parttypename;
      htitle17 = "Fake vs DecayRadius M"; htitle17+=i; htitle17+=parttypename;
      fFakeVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fFakeVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);

      hname17 = "hAllVsDecayRadiusM"; hname17+=i; hname17+=parttypename;
      htitle17 = "All vs DecayRadius M"; htitle17+=i; htitle17+=parttypename;
      fAllVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);

      hname17 = "hAllVsCosPointingAngleM"; hname17+=i; hname17+=parttypename;
      htitle17 = "All vs CosPointingAngle M"; htitle17+=i; htitle17+=parttypename;
      fAllVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),200,0.95,1.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),200,0.95,1.0,7,0.5,4);

      hname17 = "hAllVsDCAM";hname17+=i; hname17+=parttypename;
      htitle17 = "All vs DCA M";htitle17+=i; htitle17+=parttypename;
      fAllVsDCA[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsDCA[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);



      Double_t min[]={-1.5,0,-10,0};
      Double_t max[]={1.5,10,10,2*TMath::Pi()};
      Int_t nbins[]={20,100,10,20};
      
      hname1  = "hGeneratedMCPrimariesEff4DM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level eta_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0] = new THnSparseF(hname1.Data(),htitle1.Data(),4,nbins,min,max);
      hname1+="Minus";htitle1+="Minus";
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][1] = new THnSparseF(hname1.Data(),htitle1.Data(),4,nbins,min,max);

      hname3  = "hMCPrimariesThatAreReconstructed4DM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][0] = new THnSparseF(hname3.Data(),htitle3.Data(),4,nbins,min,max);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][1] = new THnSparseF(hname3.Data(),htitle3.Data(),4,nbins,min,max);

      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][1]->Sumw2();

      fPrimVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsDCA[i*PARTTYPES+j][1]->Sumw2();

      fPrimVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();

      fPrimVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();

      fAllVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fAllVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fAllVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fAllVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fAllVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fAllVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();


    }

    hname  = "pidTPCdEdx";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum";
    fHistQAPID[0][j][0] = new TH2F(hname, htitle, 2000, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTime";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum";
    fHistQAPID[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum";
    fHistQAPID[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[2][j][1]= new TH2F(hname,htitle, 2000, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum";
    fHistQAPID[3][j][0] = new TH2F(hname,htitle, 2000, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma";
    fHistQAPID[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);

    hname  = "pidTPCdEdxFail";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum Fail";
    fHistQAPIDFail[0][j][0] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTimeFail";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum Fail";
    fHistQAPIDFail[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum Fail";
    fHistQAPIDFail[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[2][j][1]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum Fail";
    fHistQAPIDFail[3][j][0] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma Fail";
    fHistQAPIDFail[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
  }

  fHistEv[0] = new TH1F("fHistEv", "Multiplicity", 100, 0, 100);
  fHistEv[1] = new TH1F("fHistEvFB16", "Multiplicity FB16", 100, 0, 200);
  fHistEv[2] = new TH1F("fHistEvFB96", "Multiplicity FB96", 100, 0, 200);
  fHistEv[3] = new TH1F("fHistEvFB128", "Multiplicity FB128", 100, 0, 200);
  for(Int_t i = 0; i < 4; i++)
    fHistoList->Add(fHistEv[i]);

  for(Int_t i = 0; i < MULTBINS; i++)  {
    hname = "fHistEventCutsM";
    hname+= i;
    

    fHistEvCuts[i] = new TH1F(hname,Form("Event Cuts M%d",i) , 6, -0.5, 5.5);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"MultSelection");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"Centrality Cut");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"AliEventCuts");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(5,"No Vertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(6,"Pileup Generated");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(7,"z-vertex>10");
    // fHistEvCuts[i]->GetXaxis()->SetBinLabel(8,"MultCut");
    // fHistEvCuts[i]->GetXaxis()->SetBinLabel(9,"NoVertex");
    // fHistEvCuts[i]->GetXaxis()->SetBinLabel(10,"");

    fHistoList->Add(fHistEvCuts[i]);

    for(Int_t chg=0;chg<2;chg++){
      hname  = "hMisidentificationM"; hname+=i; if(chg==0) hname+="Plus"; else hname+="Minus";
      htitle = "Misidentification Fraction M"; htitle+=i; if(chg==0) htitle+="Plus"; else htitle+="Minus";
      fMisidentification[i][chg] = new TH2F(hname.Data(),htitle.Data(), 3, 0.5, 3.5, 4 , 0, 4);
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(1,"Pions, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(3,"Protons, MC");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(1,"Pions, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(3,"Protons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(4,"Other, Data");
      fHistoList->Add(fMisidentification[i][chg]);
    }
  }

  fHistQA[0] = new TH1F("fHistVtx", "Z vertex distribution", 100, -15., 15.);
  fHistQA[1] = new TH1F("fHistnTpcCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTpcClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("dcaHistDcaXY1D", "DCA XY", 210, -2.1, 2.1);
  fHistQA[4] = new TH1F("dcaHistDcaZ1D", "DCA Z", 210, -2.1, 2.1);
  fHistQA[5] = new TH1F("fHistChi2Tpc", "Chi2 TPC", 100, 0., 8.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",1000,0.,10.0);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 100, -TMath::Pi(), TMath::Pi());
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 100, -2, 2);

  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Eta");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"OutofBunch Pileup");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiEta","Eta-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());
  fHistQA2D[3] = new TH2F("fHistArmenterosPodolanski - before cuts","Armenteros-Podolanski phase space - before cuts;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fHistQA2D[4] = new TH2F("fHistArmenterosPodolanski - after cuts","Armenteros-Podolanski phase space - after cuts;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);



  fHistQAK0ss[0] = new TH2F("fHistQAK0s", "V0 Details" , 8, 0.5, 8.5,100,0,20);
  fHistQAK0ss[1] = new TH2F("fHistQAK0sMinus", "V0 Details" , 8, 0.5, 8.5,100,0,20);
  for(int i=0;i<2;i++){
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(1,"AllV0s");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(2,"DaugtersNotPrimary");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(3,"Mother Pos=Neg");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(4,"MotherFoundInMCarray");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(5,"PDG of K0s");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(6,"IsNotFrom Weak");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(7,"IsNotFrom Material");
    fHistQAK0ss[i]->GetXaxis()->SetBinLabel(8,"Is Primary");
    fHistoList->Add(fHistQAK0ss[i]);
  }

  TString originK0ss[]={"PosDaughterPDG","NegDaughterPDG","PosDaugterMotherPDG","NegDaugterMotherPDG","V0PDG"};
  for(int j=0;j<5;j++){
    fOriginK0ss[j][0] = new TH2F("fOriginK0s" + originK0ss[j],"Origin of particles "+originK0ss[j] , 4000, -4000.0, 4000.0 ,100,0,20);
    fOriginK0ss[j][1] = new TH2F("fOriginK0s" + originK0ss[j]+"Minus", "Origin of particles "+originK0ss[j]+"Minus" , 400, -4000.0, 4000.0 ,100,0,20);
    fHistoList->Add(fOriginK0ss[j][0]);
    fHistoList->Add(fOriginK0ss[j][1]);
  }

  
  for ( Int_t i = 0; i < 11; i++){
    fHistoList->Add(fHistQA[i]);
    if(i<5) fHistoList->Add(fHistQA2D[i]);
    if(i<5){
     for(Int_t j = 0 ; j<PARTTYPES; j++){
       for(int chg=0;chg<2;chg++){
         fHistoList->Add(fHistQAPID[i][j][chg]);
	 fHistoList->Add(fHistQAPIDFail[i][j][chg]);
	}
      }
    }
  }
    
  for (Int_t i = 0; i < MULTBINS*PARTTYPES; i++){
    for(Int_t chg=0;chg<2;chg++){
      fHistoList->Add(fGeneratedMCPrimaries[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed[i][chg]);
      fHistoList->Add(fGeneratedMCPrimaries4D[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed4D[i][chg]);
      if(i < MULTBINS*(PARTTYPES-2))fHistoList->Add(fMCPrimariesThatAreReconstructedNoNsigma[i][chg]);
      fHistoList->Add(fReconstructedAfterCuts[i][chg]);
      fHistoList->Add(fReconstructedNotPrimaries[i][chg]);
      fHistoList->Add(fReconstructedPrimaries[i][chg]);
      fHistoList->Add(fContamination[i][chg]);


  if(i==4){ //works only for MULTBINS == 1!!!!
	fHistoList->Add(fPrimVsDCA[i][chg]);
	fHistoList->Add(fSecWeakVsDCA[i][chg]);
	fHistoList->Add(fSecMatVsDCA[i][chg]);
	fHistoList->Add(fFakeVsDCA[i][chg]);
	fHistoList->Add(fPrimVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecWeakVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecMatVsCosPointingAngle[i][chg]);
	fHistoList->Add(fFakeVsCosPointingAngle[i][chg]);

	fHistoList->Add(fPrimVsDecayRadius[i][chg]);
	fHistoList->Add(fSecWeakVsDecayRadius[i][chg]);
	fHistoList->Add(fSecMatVsDecayRadius[i][chg]);
	fHistoList->Add(fFakeVsDecayRadius[i][chg]);
	fHistoList->Add(fAllVsDCA[i][chg]);
	fHistoList->Add(fAllVsCosPointingAngle[i][chg]);
	fHistoList->Add(fAllVsDecayRadius[i][chg]);
      }
    }

  }

  fPIDKch = new TH2F("fPIDKch","Kaon PID",500,0,10,500,0.0,1000.0);
  fPIDKch->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fPIDKch->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
  fHistoList->Add(fPIDKch);

  fPIDKeCut= new TH2F("fPIDKeCut","Kaon PID with electron Rejection",500,0,10,500,0.0,1000.0);
  fPIDKeCut->GetXaxis()->SetTitle("p_{} [GeV/c]");
  fPIDKeCut->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
  fHistoList->Add(fPIDKeCut);

  fMassInvK0sFail = new TH1D("fMassInvK0sFail","Mass Assuming K0s Hypothesis Fail", 200, 0.4, 0.6);
  fHistoList->Add(fMassInvK0sFail);

  fMassInvK0s = new TH2D("fMassInvK0s","K0s invariant mass distribution",10000, -5000, 5000., 200, 0.4, 0.6);
  fHistoList->Add(fMassInvK0s);

  fMassInvK0sAfterCuts = new TH1D("fMassInvK0sAfterCuts","K0s invariant mass distribution", 200, 0.4, 0.6);
  fHistoList->Add(fMassInvK0sAfterCuts);

  fMassInvK0sPass = new TH1D("fMassInvK0sPass","Mass Assuming K0s Hypothesis Pass", 200, 0.4, 0.6);
  fHistoList->Add(fMassInvK0sPass);

  fMassInvK0sPt = new TH2D("fMassInvK0sPt","Pt vs Inv Mass(After checking PDG) ", 200, 0.4, 0.6, 100, 0, 10.0);
  fHistoList->Add(fMassInvK0sPt);
  
  fEtaK0s = new TH1D("fEtaK0s", "|Eta| distribution of K0s", 500, -0.8, 0.8);
  fPtK0s = new TH1D("fPtK0s", "Pt distribution of K0s", 500, 0.0, 8.);
  fHistoList->Add(fEtaK0s);
  fHistoList->Add(fPtK0s);
  
  fCutsK0s = new TH1D("fCutsK0s","Cuts K0s", 20, 0.5, 20.5);
  fHistoList->Add(fCutsK0s);
  

  fTruePtK0sMC = new TH2D("fTruePtK0sMC","True pT of K0ss MC",10000, -5000, 5000.,50,0.,10.0);
  fRecPtK0sMC = new TH2D("fRecPtK0sMC","Rec pT of K0ss MC",10000, -5000, 5000.,50,0.,10.0);
  fHistoList->Add(fTruePtK0sMC);
  fHistoList->Add(fRecPtK0sMC);
  
  
  
  //********** PID ****************

  // AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  // fpidResponse = inputHandler->GetPIDResponse();
  // std::cout<<"*******"<< fpidResponse<<std::endl;



  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fAODpidUtil = aodH->GetAODpidUtil();

  if(fIfAliEventCuts){
    fEventCuts = new AliEventCuts();
  }

  PostData(1, fHistoList);
}


//_____________________________________________________________________

bool IsPionNSigmaAM(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime){
  if (mom > 0.5){
    if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 2)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCPi) < 2)
      return true;
  }
  return false;
}

bool IsPionNSigmaV0AM(float mom, float nsigmaTPCPi, float nsigmaTOFPi){
  if (TMath::Abs(nsigmaTPCPi) < 3.0) return true;
  return false;
}

bool IsPionNSigmaV0TPC5AM(float mom, float nsigmaTPCPi, float nsigmaTOFPi){
  if (TMath::Abs(nsigmaTPCPi) < 5.0) return true;
  return false;
}


bool IsProtonNSigmaV0TPC5AM(float mom, float nsigmaTPCP, float nsigmaTOFP){
  if (TMath::Abs(nsigmaTPCP) < 5.0) return true;
  return false;
}

bool IsPionNSigma3AM(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime){
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCPi) < 3)
      return true;
  }
  return false;
}

bool IsKaonNSigmaAM(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime){
  if (mom > 0.5) {
    //rejection of unwanted contamination
    if(mom>1 && TOFtime<-400)
      return false;
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 2)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 2)
      return true;
  }
  return false;
}

bool IsKaonNSigma3AM(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
    //rejection of unwanted contamination
    // if(mom>1 && TOFtime<-400)
    // return false;
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 3)
      return true;
  }
  return false;
}

bool IsProtonNSigmaAM(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {
    if(mom>1.8 && TOFtime<-300)
      return false;
    
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 2)
      return true;
	}
  else {
    if (TMath::Abs(nsigmaTPCP) < 2)
      return true;
  }
  return false;
}

bool IsProtonNSigmaV0AM(float mom, float nsigmaTPCP, float nsigmaTOFP){
  if (mom < 0.8){
    if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
  } 
  else {
    if (nsigmaTOFP < -999.) {
      if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
    } 
    else {
      if (TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0) return true;
    }
  }
  return false;
}

bool IsProtonNSigma3AM(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {   
    // if(mom>1.8 && TOFtime<-300)
    // return false;
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCP) < 3)
      return true;
  }
  return false;
}



bool AliAnalysisTaskEffK0ss::IsElectronAM(float nsigmaTPCe, float nsigmaTPCPi, float nsigmaTPCK, float nsigmaTPCP)
{
    if (TMath::Abs(nsigmaTPCe) < fNsigmaElectronRejection &&
        TMath::Abs(nsigmaTPCPi) > fNsigmaElectronRejection &&
        TMath::Abs(nsigmaTPCK) > fNsigmaElectronRejection &&
        TMath::Abs(nsigmaTPCP) > fNsigmaElectronRejection)
    {
        return true;
    }
    return false;
}

bool AliAnalysisTaskEffK0ss::IsElectronAM1(float nsigmaTPCe)
{
    if (TMath::Abs(nsigmaTPCe) < fNsigmaElectronRejection)
        return true;
    else
        return false;
}


//_______________________________________________________

void AliAnalysisTaskEffK0ss::UserExec(Option_t *)
{
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliAODEvent *fAOD = aodH->GetEvent();

  /***Get Event****/


  
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  if(aodEvent->GetNumberOfTracks() <=0)
   return;
  fHistEvCuts[0]->Fill(0);
  
  AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
  if(!MultSelection) return;


  //if (!fEventCuts->AcceptEvent(fAOD)) return;
  fHistEvCuts[0]->Fill(1);


  Double_t mult;
  if(fEstEventMult == kRefMult)
    {  
      AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
      mult = fAODheader->GetRefMultiplicity();
      //cout << "multiplicity value is" << mult << endl;
    }
  else if(fEstEventMult == kV0M)
    {
     mult = MultSelection->GetMultiplicityPercentile("V0M");
    
    }
  else if(fEstEventMult == kV0A)
    {
      mult = fEventCuts->GetCentrality(); 
      //cout << "centrality of event is " << mult << endl;
      //AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
      //mult = alicent->GetCentralityPercentile("V0A");
    }

  fHistEv[0]->Fill(mult);
  
  if ((mult < fCentMin)||(mult > fCentMax)) return;
  //cout << "centrality value is" << mult << endl;

  fHistEvCuts[0]->Fill(2);

  if(fIfAliEventCuts){
    if (!fEventCuts->AcceptEvent(aodEvent)) {
      return;
    }
  }

  fHistEvCuts[0]->Fill(3);

  // EVENT SELECTION ********************
  fHistQA[9]->Fill(1);

  //****** Multiplicity selection *********
  Int_t fcent = -999;
  //if(mult >= 0 && mult <=20)  fcent = 0;
  //else if(mult >= 20 && mult <=39) fcent = 1;
  //else if(mult >= 40 && mult <=59) fcent = 2;
  //else if(mult >= 60 && mult <=90) fcent = 3;
  //else if(mult >= 99990 && mult <=99936) fcent = 4;
  //else if(mult >= 999937 && mult <=99944) fcent = 5;
  //else if(mult >= 999945 && mult <=99957) fcent = 6;
  //else if(mult >= 999958 && mult <=99149) fcent = 6;
  //else fcent = 7;
  //if (fcent == 7) return;

  // if(mult >= 2&& mult <=150)  fcent = 0;
  // else if(mult >= 2 && mult <=19) fcent = 1;
  // else if(mult >= 20 && mult <=49) fcent = 2;
  // else if(mult >= 50 && mult <=150) fcent = 3;
  // else return;

  fcent = 0;

  //if(fcent==0)
  //else if(fcent==1)fHistEvCuts[1]->Fill(2);
  //else if(fcent==2)fHistEvCuts[2]->Fill(2);
  //else if(fcent==3)fHistEvCuts[3]->Fill(2);

  //"ESDs/pass2/AOD049/*AliAOD.root");
  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  vertex->GetPosition(fV1K0s);
  if (!vertex || vertex->GetNContributors()<=0) return;

  fHistEvCuts[0]->Fill(4);

  fHistQA[9]->Fill(2);
  // if(fcent==0)fHistEvCuts[0]->Fill(3);
  // else if(fcent==1)fHistEvCuts[1]->Fill(3);
  // else if(fcent==2)fHistEvCuts[2]->Fill(3);
  // else if(fcent==3)fHistEvCuts[3]->Fill(3);

  //********* Pile-up removal*******************
  //check this: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup
  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();

  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fOutOfBunchPlp = kFALSE;
  Bool_t fisPileUp = kFALSE;

  if(fpA2013)
    if(anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;


  //Multiple vertices with tracks
  if(fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);
  //if this (fMVPlp) is false, than rejection based on Multiple SPD vertices is used

  //out-of-bunch pile-up rejection from AnaUtil
  anaUtil->SetUseOutOfBunchPileUp(fOutOfBunchPlp);

  /*
  Int_t fMinPlpContribMV = 0; //default value: fMinPlpContribMV(5),
  Int_t fMinPlpContribSPD = 3; //default value: fMinPlpContribSPD(5),

  if(fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if(fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);
  */


  if(fisPileUp)
    if(anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;

  //pileup for LHC20e3a -> Injective Pileup over events 
  AliAODMCHeader *mcHeader = 0;
  mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
    return;
  }
  Bool_t isPileupInGeneratedEvent = kFALSE;
  isPileupInGeneratedEvent = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader,"Hijing");
  if(isPileupInGeneratedEvent) return;

  fHistEvCuts[0]->Fill(5);

  fHistQA[9]->Fill(3);
  //if(fcent==0)fHistEvCuts[0]->Fill(4);
  //else if(fcent==1)fHistEvCuts[1]->Fill(4);
  //else if(fcent==2)fHistEvCuts[2]->Fill(4);
  //else if(fcent==3)fHistEvCuts[3]->Fill(4);
//***************************************************

  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > fPVzCut) return;

  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);

  fHistEvCuts[0]->Fill(6);
  //if(fcent==0)fHistEvCuts[0]->Fill(5);
  //else if(fcent==1)fHistEvCuts[1]->Fill(5);
  //else if(fcent==2)fHistEvCuts[2]->Fill(5);
  //else if(fcent==3)fHistEvCuts[3]->Fill(5);

 //**** getting MC array ******
  TClonesArray  *arrayMC;

  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));


  //copying pid information for FB 128
  int labels[20000];
  for (int il=0; il<20000; il++) labels[il] = -1;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<aodEvent->GetNumberOfTracks();i++) {
    const AliAODTrack *aodtrack=(AliAODTrack*)aodEvent->GetTrack(i);
    if (!aodtrack->TestFilterBit(128)) {
      if(aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }



  //RECONSTRUCTED TRACKS

  TObjArray recoParticleArray[PARTTYPES];

  Int_t iTracks = (aodEvent->GetNumberOfTracks());
  //cout << "Number of tracks in an event is " << iTracks << endl;
  Double_t MultTPC = fAOD->GetNumberOfTPCTracks();
  //cout << "Number of TPC tracks in an event is " << MultTPC << endl;
  Int_t nv0s = (fAOD->GetNumberOfV0s());
  //cout << "Number of v0s in an event is " << nv0s << endl;

  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  //loop over AOD tracks

  int multFB128=0, multFB96=0, multFB16=0;

  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) // MC reconstructed loop
  {

    //get track
    //AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esdEvent),iTracks);
    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks);
    if (!track)continue;
    fHistQA[10]->Fill(2);


    if(track->TestFilterBit(128)) multFB128++;
    if(track->TestFilterBit(16)) multFB16++;
    if(track->TestFilterBit(96)) multFB96++;

    //UInt_t filterBit = (1 << (0));
    UInt_t filterBit = fFB;
    if(!track->TestFilterBit(filterBit))continue;

    // For TPC Only tracks we have to copy PID information from corresponding global tracks
    const Int_t pid_track_id = (fFB == 128)
                             ? labels[-1 - aodEvent->GetTrack(iTracks)->GetID()]
                             : iTracks;

    const auto *aodtrackpid2 = static_cast<AliAODTrack *>(aodEvent->GetTrack(pid_track_id));
    

    //Pile-up removal
    if (fTrackPileUpRemoval) {
      //method which checks if track
      //have at least 1 hit in ITS or TOF.
      bool passTrackPileUp = false;

      // does tof timing exist for our track?
      if (aodtrackpid2->GetTOFBunchCrossing() == 0) {
        passTrackPileUp = true;
      }

      // check ITS refit
      if (!(aodtrackpid2->GetStatus() & AliESDtrack::kITSrefit)) {
        continue;
      }

      // loop over 2 ITS layers and check for a hit!
      for (int i : {0, 1}) {
        if (aodtrackpid2->HasPointOnITSLayer(i)) {
          passTrackPileUp = true;
        }
      }

      if (!passTrackPileUp) {
        continue;
      }
    }

    //charge
    //  if(track->Charge() < 0 ) continue;
    Int_t charge = 0;
    if(track->Charge() > 0 ) charge=0;
    else if (track->Charge() < 0 ) charge=1;
    fHistQA[10]->Fill(3);

    if(track->Eta() < -0.8 || track->Eta() > 0.8) continue;
    fHistQA[10]->Fill(4);

    if (track->Pt() < 0.10 || track->Pt() > 10) continue; // Track pT cut
    fHistQA[10]->Fill(5);

    //pileup for LHC20e3a -> Injective Pileup over tracks 
   Bool_t isParticleFromOutOfBunchPileupCollision = kFALSE;
   isParticleFromOutOfBunchPileupCollision = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTracks,mcHeader,arrayMC);
   if(isParticleFromOutOfBunchPileupCollision) continue;

   fHistQA[10]->Fill(6);

    //single track cuts
    // if(track->Chi2perNDF() > 4.0) continue;
    // if(track->GetTPCNcls() < 70) continue;

    //DCA

    Double_t DCAXY;
    Double_t DCAZ;

    DCAXY = -TMath::Abs(track->DCA());
    DCAZ = -TMath::Abs(track->ZAtDCA());

    if(!(DCAXY==-999 || DCAZ==-999)){
	//if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
	//no DCA cut
	//if(TMath::Abs(DCAXY) > 1000.0) {continue;} //XY
	//if(TMath::Abs(DCAZ) > 1000.0) {continue;} //Z
    }
    else {
      // const AliAODVertex* vertex = (AliAODVertex*) aodEvent->GetPrimaryVertex(); (already defined above)
      float vertexX  = -999.;
      float vertexY  = -999.;
      float vertexZ  = -999.;

      if(vertex) {
      	Double32_t fCov[6];
      	vertex->GetCovarianceMatrix(fCov);
      	if(vertex->GetNContributors() > 0) {
      	  if(fCov[5] != 0) {
      	    vertexX = vertex->GetX();
      	    vertexY = vertex->GetY();
      	    vertexZ = vertex->GetZ();

      	  }
      	}
      }

      Double_t pos[3];
      track->GetXYZ(pos);

      Double_t DCAX = pos[0] - vertexX;
      Double_t DCAY = pos[1] - vertexY;
      DCAZ = pos[2] - vertexZ;
      DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
    }

    if (fUseDcaCuts) {
      if (TMath::Abs(DCAXY) > fDcaXYCut) continue;
      if (TMath::Abs(DCAZ) > fDcaZCut) continue;
  }

    AliAODTrack* aodtrackpid;

    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7)))
      aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
    else
      aodtrackpid = track;

    //Electron rejection
    double nSigmaTPCPi = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    double nSigmaTPCK = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    double nSigmaTPCP = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    double nSigmaTPCe = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);

 
    if (fElectronReject && IsElectronAM1(nSigmaTPCe)) continue;

    fHistQA[10]->Fill(7);

    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1));
    //fHistQA[2]->Fill(track->GetTPCNclsF());
    fHistQA[3]->Fill(DCAXY);
    fHistQA[4]->Fill(DCAZ);
    Float_t chi2Tpc = track->Chi2perNDF();
    fHistQA[5]->Fill(chi2Tpc);
    fHistQA[6]->Fill(track->Pt());

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();

    fHistQA[7]->Fill(ph);
    fHistQA[8]->Fill(track->Eta());
    fHistQA2D[2]->Fill(track->Eta(),ph);
    fHistQA2D[0]->Fill(tPt,DCAXY);
    fHistQA2D[1]->Fill(tPt,DCAZ);

    //PID monitors
    double nSigmaTOFPi = -1000;
    double nSigmaTOFK = -1000;
    double nSigmaTOFP = -1000;

     //

    ULong_t status = aodtrackpid->GetStatus();
    Float_t probMis;
    if (((status & AliVTrack::kTOFout) == AliVTrack::kTOFout)
	&& ((status & AliVTrack::kTIME) == AliVTrack::kTIME)) {

      probMis = fAODpidUtil->GetTOFMismatchProbability(aodtrackpid);
      if(probMis < 0.01){
	nSigmaTOFPi = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
	nSigmaTOFK = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
	nSigmaTOFP = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);
      }

    }

    float tdEdx = aodtrackpid->GetTPCsignal();
    float tTofSig = aodtrackpid->GetTOFsignal();
    double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);


    fHistQAPID[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPID[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPID[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPID[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPID[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);

    fHistQAPIDFail[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPIDFail[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPIDFail[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPIDFail[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPIDFail[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);


    bool isPionNsigma = 0;
    bool isKaonNsigma = 0;
    bool isProtonNsigma  = 0;

     if(fPidMethod==kNSigma){
    //******** With double counting *******************
      isPionNsigma = (IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      isKaonNsigma = (IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      isProtonNsigma = (IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    else if(fPidMethod==kNSigmaNoDoubleCounting){
      //******** Without double counting *******************
      double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;

      if(track->Pt()<0.5){
       nSigmaPIDPi = abs(nSigmaTPCPi);
       nSigmaPIDK  = abs(nSigmaTPCK);
       nSigmaPIDP  = abs(nSigmaTPCP);
      }
      else{
       nSigmaPIDPi = TMath::Hypot(nSigmaTPCPi,nSigmaTOFPi);
       nSigmaPIDK= TMath::Hypot(nSigmaTPCK,nSigmaTOFK);
       nSigmaPIDP= TMath::Hypot(nSigmaTPCP,nSigmaTOFP);
      }

      if(nSigmaPIDPi<nSigmaPIDK && nSigmaPIDPi<nSigmaPIDP){
       isPionNsigma = (IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      }
      else if(nSigmaPIDK<nSigmaPIDPi && nSigmaPIDK<nSigmaPIDP){
       isKaonNsigma = (IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      }
      else if(nSigmaPIDP<nSigmaPIDPi && nSigmaPIDP<nSigmaPIDK){
       isProtonNsigma = (IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      }
    }
    else if(fPidMethod==kExclusivePID){
      //******** Exclusive PID ********************
      isPionNsigma = (IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    else if(fPidMethod==kExclusivePIDDiffRejection){
      //******** Exclusive PID, different rejection  ********************
      isPionNsigma = (IsPionNSigmaAM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigma3AM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3AM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigma3AM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaAM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3AM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigma3AM(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigma3AM(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaAM(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    if (isPionNsigma){
      fHistQAPID[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    else{
      fHistQAPIDFail[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPIDFail[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPIDFail[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPIDFail[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPIDFail[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    if (isKaonNsigma){
      fHistQAPID[0][2][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
      fHistQAPID[2][2][charge]->Fill(tPt,nSigmaTOFK);
      fHistQAPID[3][2][charge]->Fill(tPt,nSigmaTPCK);
      fHistQAPID[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
    }
    else {
      fHistQAPIDFail[0][2][charge]->Fill(tPt,tdEdx);
      fHistQAPIDFail[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
      fHistQAPIDFail[2][2][charge]->Fill(tPt,nSigmaTOFK);
      fHistQAPIDFail[3][2][charge]->Fill(tPt,nSigmaTPCK);
      fHistQAPIDFail[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
    }
    if (isProtonNsigma){
      fHistQAPID[0][3][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
      fHistQAPID[2][3][charge]->Fill(tPt,nSigmaTOFP);
      fHistQAPID[3][3][charge]->Fill(tPt,nSigmaTPCP);
      fHistQAPID[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
    }
    else{
      fHistQAPIDFail[0][3][charge]->Fill(tPt,tdEdx);
      fHistQAPIDFail[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
      fHistQAPIDFail[2][3][charge]->Fill(tPt,nSigmaTOFP);
      fHistQAPIDFail[3][3][charge]->Fill(tPt,nSigmaTPCP);
      fHistQAPIDFail[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
    }
    fReconstructedAfterCuts[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());//Fills hist. for all reconstructed particles after cuts
 
   if(!arrayMC){
      continue;
    }
    //get coresponding MC particle 
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label);

   //getting no. of tracks for each particle species after all the cuts:

    //********* PID - pions ********
     if (isPionNsigma){
      if(track->Eta() < -0.8 || track->Eta() > 0.8) continue; 
      if (track->Pt() > 0.1 && track->Pt() < 10)
        fReconstructedAfterCuts[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
      if (!MCtrk) continue;
      recoParticleArray[1].Add(MCtrk);
     }
     //Fills for all identified pions found after cuts (reconstructed) - numerator for Efficiency

     //********* PID - kaons ********
     if (isKaonNsigma){
       if(track->Eta() < -0.7 || track->Eta() > 0.7) continue; 
       if (track->Pt() > 0.1 && track->Pt() < 10)
         fReconstructedAfterCuts[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[2].Add(MCtrk);
     }
     //Fills for all identified kaons found after cuts (reconstructed) - numerator for Efficiency

    //********* PID - protons ********
     /*if (isProtonNsigma){
       if(track->Eta() < -0.5 || track->Eta() > 0.5) continue;
       if (track->Pt() > 0.5 && track->Pt() < 2.5)
         fReconstructedAfterCuts[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[3].Add(MCtrk);
     }*/

     //Fills for all identified protos found after cuts (reconstructed) - numerator for Efficiency
     //******************************

    
    if (!MCtrk) continue;
    if(MCtrk->Charge()==0){std::cout<<"!!!"<<std::endl; continue;} // skip neutral tracks
    recoParticleArray[0].Add(MCtrk);


    //Fills histogram for particles that are contamination from secondaries:
    if (!MCtrk->IsPhysicalPrimary()) {
      fReconstructedNotPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }


    int PDGcode = MCtrk->GetPdgCode();

    //And secondaries for different particle species:
    if (!MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) { //secondaries in pions
      fReconstructedNotPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) {
      fReconstructedPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) { //secondaries in kaons
      fReconstructedNotPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) {
      fReconstructedPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) { //secondaries in protons
      fReconstructedNotPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) {
      fReconstructedPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
    }
    



    //step 1, TOF Matching
    UInt_t statusTOF;
    statusTOF=track->GetStatus();
    if((statusTOF&AliVTrack::kTOFout)==0 || (statusTOF&AliVTrack::kTIME)==0)
      statusTOF=0;
    if(track->Pt()<0.5) statusTOF = 1;

    //Misidentification fraction
    if(abs(PDGcode)==211)
    {
    	if(isPionNsigma)
    	  fMisidentification[fcent][charge]-> Fill(1,0.5);
    	if(isKaonNsigma)
    	  fMisidentification[fcent][charge]-> Fill(1,1.5);
    	if(isProtonNsigma)
    	  fMisidentification[fcent][charge]-> Fill(1,2.5);
    	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
    	  if(statusTOF)
    	    fMisidentification[fcent][charge]-> Fill(1,3.5);


    }
    else if(abs(PDGcode)==321)
    {
    	if(isPionNsigma)
    	  fMisidentification[fcent][charge]-> Fill(2,0.5);
    	if(isKaonNsigma)
    	  fMisidentification[fcent][charge]-> Fill(2,1.5);
    	if(isProtonNsigma)
    	  fMisidentification[fcent][charge]-> Fill(2,2.5);
    	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
    	  if(statusTOF)
    	    fMisidentification[fcent][charge]-> Fill(2,3.5);


    }
    else if(abs(PDGcode) == 2212)
    {
    	if(isPionNsigma)
    	  fMisidentification[fcent][charge]-> Fill(3,0.5);
    	if(isKaonNsigma)
    	  fMisidentification[fcent][charge]-> Fill(3,1.5);
    	if(isProtonNsigma)
    	  {
    	  fMisidentification[fcent][charge]-> Fill(3,2.5);
    	  }
    	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
    	  if(statusTOF)
    	    fMisidentification[fcent][charge]-> Fill(3,3.5);
    }


    fContamination[PARTTYPES*fcent][charge]-> Fill(PDGcode,track->Pt());
    //Contaminations: "how many pions are in the kaons sample"? etc.
    //Do not use for corrections: using those values will be dependant on i.e. Pi/K ratio in MC
    //Use misidentification fraction instead
    if(isPionNsigma){
      fContamination[PARTTYPES*fcent+1][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for pions
    }
    if(isKaonNsigma){
      fContamination[PARTTYPES*fcent+2][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for kaons
    }
    if(isProtonNsigma){
      fContamination[PARTTYPES*fcent+3][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for protons
    }

  }

  fHistEv[1]->Fill(multFB16);
  fHistEv[2]->Fill(multFB96);
  fHistEv[3]->Fill(multFB128);

  //loop over V0s
  for (Int_t i = 0; i < aodEvent->GetNumberOfV0s(); i++){

    double MassK0Short = 0.497613;
    int cutK0s = 1;
    fCutsK0s->Fill(cutK0s++);

    AliAODv0 *aodv0 = aodEvent->GetV0(i);
    if (!aodv0) continue;

    //pileup for LHC20e3a -> Injective Pileup over tracks 
    Bool_t isParticleFromOutOfBunchPileupCollision = kFALSE;
    isParticleFromOutOfBunchPileupCollision = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,mcHeader,arrayMC);
    if(isParticleFromOutOfBunchPileupCollision) continue;

    fCutsK0s->Fill(cutK0s++);
    fMassInvK0sFail->Fill(aodv0->MassK0Short()); //invariant mass distribution before cuts
    //if (aodv0->GetNDaughters() > 2) continue;
    //if (aodv0->GetNProngs() > 2) continue;
    if (aodv0->GetCharge() != 0) continue;
    //if (aodv0->ChargeProng(0) == aodv0->ChargeProng(1)) continue;

    double cospa = aodv0->CosPointingAngle(fV1K0s);
    //cout << "cos pointing angle is" << cospa << endl;
    if (aodv0->CosPointingAngle(fV1K0s) < fMinCosPointingAngle) continue;
    
    fCutsK0s->Fill(cutK0s++);

    AliAODTrack *daughterTrackPos = (AliAODTrack *)aodv0->GetDaughter(0); //getting positive daughter track
    AliAODTrack *daughterTrackNeg = (AliAODTrack *)aodv0->GetDaughter(1); //getting negative daughter track
    if (!daughterTrackPos) continue; //daughter tracks must exist
    if (!daughterTrackNeg) continue;
    //if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) continue; //and have different charge

fCutsK0s->Fill(cutK0s++);

if(aodv0->Pt() < 0.10 || aodv0->Pt() > 10) continue;
fCutsK0s->Fill(cutK0s++);

if(TMath::Abs(aodv0->Eta()) > 0.8) continue;
fCutsK0s->Fill(cutK0s++);

if(aodv0->GetOnFlyStatus() == kTRUE) continue;
fCutsK0s->Fill(cutK0s++);

if(TMath::Abs(aodv0->DcaV0Daughters()) > fMaxDcaV0Daughters) continue;  // ✅ replaced
fCutsK0s->Fill(cutK0s++);

if(TMath::Abs(aodv0->DcaPosToPrimVertex()) < fMinDcaPosToPrimVertex) continue;  // ✅ replaced
fCutsK0s->Fill(cutK0s++);

if(TMath::Abs(aodv0->DcaNegToPrimVertex()) < fMinDcaNegToPrimVertex) continue;  // ✅ replaced
fCutsK0s->Fill(cutK0s++);

// Optional: check V0 to PV DCA
// if(TMath::Abs(aodv0->DcaV0ToPrimVertex()) > 0.1) continue;
// fCutsK0s->Fill(cutK0s++);

Double_t decayL = aodv0->DecayLength(fV1K0s);
Double_t dcav0PV = TMath::Abs(aodv0->DcaV0ToPrimVertex());
Double_t cTau = decayL * (aodv0->MassK0Short() / aodv0->P());
if(cTau > fMaxCTauK0s) continue; 
fCutsK0s->Fill(cutK0s++);

Double_t radius = aodv0->RadiusV0();
if (radius < fMinV0Radius) { 
    AliDebugClass(2, "Failed fiducial volume");
    continue;
}

Double_t armpt = aodv0->PtArmV0();
Double_t alpha = aodv0->AlphaV0();
fHistQA2D[3]->Fill(alpha, armpt);
if (armpt <= 0.2*fabs(alpha)) continue;
fHistQA2D[4]->Fill(alpha, armpt);

AliAODTrack *trackpos = (AliAODTrack*)aodv0->GetDaughter(0);
AliAODTrack *trackneg = (AliAODTrack*)aodv0->GetDaughter(1);
if((!trackpos) || (!trackneg)) continue;
fCutsK0s->Fill(cutK0s++);

if(TMath::Abs(trackpos->Eta()) > 0.8) continue;
if(TMath::Abs(trackneg->Eta()) > 0.8) continue;
if(trackpos->GetTPCNcls() < 70) continue;
if(trackneg->GetTPCNcls() < 70) continue;

if(!(trackpos->GetStatus() & (AliESDtrack::kTPCrefit))) continue;
if(!(trackneg->GetStatus() & (AliESDtrack::kTPCrefit))) continue;
fCutsK0s->Fill(cutK0s++);

// PID cuts
double nSigmaTPCPiPos = fAODpidUtil->NumberOfSigmasTPC(trackpos, AliPID::kPion);
double nSigmaTPCPiNeg = fAODpidUtil->NumberOfSigmasTPC(trackneg, AliPID::kPion);

bool isPionNsigmaPos = TMath::Abs(nSigmaTPCPiPos) < fNsigmaDaughters; 
bool isPionNsigmaNeg = TMath::Abs(nSigmaTPCPiNeg) < fNsigmaDaughters; 

if (!isPionNsigmaPos || !isPionNsigmaNeg) continue;

     //bool K0s = false;

  
      //K0s
    //if(isPionNsigmaPos && isPionNsigmaNeg){
      fCutsK0s->Fill(cutK0s++);  
      // if(trackpos->Pt() < 0.10 || trackpos->Pt() > 10) continue; //pion plus
	    // fCutsK0s->Fill(cutK0s++);
      //if(trackneg->Pt() < 0.10 || trackneg->Pt() > 10) continue; //pion minus
      //fCutsK0s->Fill(cutK0s++);
      fMassInvK0sAfterCuts->Fill(aodv0->MassK0Short());
      if(aodv0->MassK0Short() < 0.48 || aodv0->MassK0Short() > 0.51) continue;

      // if(aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
        fCutsK0s->Fill(cutK0s++);
        fMassInvK0sPass->Fill(aodv0->MassK0Short());
        fPtK0s->Fill(aodv0->Pt());
        fEtaK0s->Fill(aodv0->Eta());
        fHistQAK0ss[0]->Fill(1,aodv0->Pt());
        fAllVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1K0s),aodv0->Pt());
        fAllVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	    //K0s = true;

        fReconstructedAfterCuts[PARTTYPES*fcent+4][0]->Fill(aodv0->Eta(), aodv0->Pt());
  
    //}

    //fReconstructedAfterCuts[PARTTYPES*fcent+4][0]->Fill(aodv0->Eta(), aodv0->Pt());

    //********* PID - K0s ********

    if(!arrayMC) continue;


    //get coresponding MC particles
    Int_t labelPos = TMath::Abs(trackpos->GetLabel());
    Int_t labelNeg = TMath::Abs(trackneg->GetLabel());
    AliAODMCParticle *MCtrkPos = (AliAODMCParticle*)arrayMC->At(labelPos);
    AliAODMCParticle *MCtrkNeg = (AliAODMCParticle*)arrayMC->At(labelNeg);

    Int_t motherPos = MCtrkPos->GetMother();
    Int_t motherNeg = MCtrkNeg->GetMother();

   // if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
      //********* PID - K0ss, ********
      //if(K0s){
        int charge = 0;
        fOriginK0ss[0][0]->Fill(MCtrkPos->GetPdgCode(), aodv0->Pt());
        fOriginK0ss[1][0]->Fill(MCtrkNeg->GetPdgCode(), aodv0->Pt());
      //}
    //}

    if(MCtrkPos->IsPhysicalPrimary() || MCtrkNeg->IsPhysicalPrimary()) 
      continue;

    //if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){

      fHistQAK0ss[0]->Fill(2,aodv0->Pt());
      fOriginK0ss[2][0]->Fill(((AliAODMCParticle*)arrayMC->At(motherPos))->GetPdgCode(), aodv0->Pt());
      fOriginK0ss[3][0]->Fill(((AliAODMCParticle*)arrayMC->At(motherNeg))->GetPdgCode(), aodv0->Pt());
    
  // }
     
    if(motherPos != motherNeg) continue;
    fHistQAK0ss[0]->Fill(3,aodv0->Pt());

    
    AliAODMCParticle *MCtrkMother = (AliAODMCParticle*)arrayMC->At(motherPos);
    if(!MCtrkMother) continue;
    fHistQAK0ss[0]->Fill(4,aodv0->Pt());

    int pdgMother = MCtrkMother->GetPdgCode();
    fMassInvK0s->Fill(pdgMother, aodv0->MassK0Short());


    if(MCtrkMother->IsPhysicalPrimary() && pdgMother==310) {
    //if( pdgMother==310) {
	    fReconstructedPrimaries[PARTTYPES*fcent+4][0]->Fill(aodv0->Eta(), aodv0->Pt()); // filling hHistoReconstructedPrimariesM0K0s
      fMassInvK0sPt->Fill( aodv0->MassK0Short(),aodv0->Pt());
    } 
    
    if (MCtrkMother->IsPhysicalPrimary() && pdgMother==310){
	    fPrimVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	    //if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	      fPrimVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1K0s),aodv0->Pt()); 
	      fPrimVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	    }
	 // }
      
    else if(MCtrkMother->IsSecondaryFromWeakDecay() && pdgMother==310){
	    fSecWeakVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	 //   if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	      fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1K0s),aodv0->Pt()); 
	      fSecWeakVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	    }
	 // }
     
    else if(MCtrkMother->IsSecondaryFromMaterial() && pdgMother==310){
	    fSecMatVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt());
	  //  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	      fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1K0s),aodv0->Pt()); 
	      fSecMatVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	    }
	//  }
     
    
     
    //if( aodv0->DcaV0ToPrimVertex() > fDCAtoPrimVtx) continue; // we do not longer need full DCA V0 to prim vertex sample

    //contamination from secondaries
    if (!MCtrkMother->IsPhysicalPrimary() && pdgMother==310) { //secondaries in K0ss
    	fReconstructedNotPrimaries[PARTTYPES*fcent+4][0]->Fill(aodv0->Eta(), aodv0->Pt());
    } 
    // else if(MCtrkMother->IsPhysicalPrimary() && pdgMother==310) {
	  //   fReconstructedPrimaries[PARTTYPES*fcent+4][0]->Fill(aodv0->Eta(), aodv0->Pt()); // filling hHistoReconstructedPrimariesM0K0s
    // } 
    
	  fTruePtK0sMC->Fill(pdgMother,MCtrkMother->Pt());
	  fRecPtK0sMC->Fill(pdgMother,aodv0->Pt());

	  recoParticleArray[4].Add(MCtrkMother);
	  if(pdgMother==310){
      
	    fHistQAK0ss[0]->Fill(5,aodv0->Pt());
	  } 
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromWeakDecay()) {
	    fHistQAK0ss[0]->Fill(6,aodv0->Pt());
	  } 
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromMaterial()) {
	    fHistQAK0ss[0]->Fill(7,aodv0->Pt());
	  }
	  else continue;
	  if(MCtrkMother->IsPhysicalPrimary()) {
	    fHistQAK0ss[0]->Fill(8,aodv0->Pt());
	  }

	  fOriginK0ss[4][0]->Fill(pdgMother, aodv0->Pt());
	
  } // k0s loop
      
      
   // loop over v0s end
  // MONTECARLO PARTICLES 

  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  // loop over MC stack
  for (Int_t ipart = 0; ipart < arrayMC->GetEntriesFast(); ipart++)  // MC array Truth loop
    {
    //std::cout<<"Entered MC loop"<<std::endl;

    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    //std::cout<<"particle obtained"<<std::endl;

    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode());


    //if(MCtrk->Charge() == 0) continue;
    Int_t charge=0;

    if(MCtrk->Charge() < 0) charge=1;
    else if(MCtrk->Charge() > 0) charge=0;

    if(MCtrk->Charge() == 0){	
	    if(MCtrk->GetPdgCode() == 310) charge = 0; //k0s 
    }
    if(MCtrk->Eta() < -0.8 || MCtrk->Eta() > 0.8){
      continue;
    }
	
    if(MCtrk->GetPdgCode() == 211){              // piplus 
    	if (MCtrk->Pt() < 0.10 || MCtrk->Pt() > 10){
	      continue;
	}
    }
    if(MCtrk->GetPdgCode() == 321){              // k+ we need to include k-?
    	if (MCtrk->Pt() < 0.10 || MCtrk->Pt() > 10){ 
	      continue;
	    }
    }
      if(MCtrk->GetPdgCode() == 2212){          // proton 
        if (MCtrk->Pt() < 0.10 || MCtrk->Pt() > 10){
	        continue;
	      }	
      }
      if(MCtrk->GetPdgCode() == 310){
      	if (MCtrk->Pt() < 0.10 || MCtrk->Pt() > 10){
	        continue;
	      }
      }
    
      // check physical primary 
    if(MCtrk->IsPhysicalPrimary()) // Not from weak decay!
      {

        //add pileup cut for tracks

	// Filling histograms for MC truth particles
	//if(PDGcode==211 || PDGcode==321 || PDGcode==2212) // all = charged hadrons
	  //fGeneratedMCPrimaries[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());


	Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	if(PDGcode==211 || PDGcode==321 || PDGcode==2212) // all = charged hadrons
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES][charge]->Fill(val);

      if(PDGcode==211)
      {
        fGeneratedMCPrimaries[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
        fGeneratedMCPrimaries4D[fcent*PARTTYPES+1][charge]->Fill(val);
      }
      else if(PDGcode==321)
      {
        fGeneratedMCPrimaries[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
        fGeneratedMCPrimaries4D[fcent*PARTTYPES+2][charge]->Fill(val);
      }
      else if(PDGcode==2212)
      {
        fGeneratedMCPrimaries[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
        fGeneratedMCPrimaries4D[fcent*PARTTYPES+3][charge]->Fill(val);
      }
      else if(PDGcode==310)
      {
        fGeneratedMCPrimaries[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
        fGeneratedMCPrimaries4D[fcent*PARTTYPES+4][charge]->Fill(val);
      }
      
      //Filling data from MC truth particles only for particles that were reconstruced
    	if (recoParticleArray[0].Contains(MCtrk)){ //All
	  if(PDGcode==211 || PDGcode==321 || PDGcode==2212) // all = charged hadrons
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());

    	  Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	  if(PDGcode==211 || PDGcode==321 || PDGcode==2212) // all = charged hadrons
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES][charge]->Fill(val);

    	  fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	  if(PDGcode==211)
    	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	  if(PDGcode==321)
    	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	  if(PDGcode==2212)
    	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  }


    	if (recoParticleArray[1].Contains(MCtrk)){ //Pions
    	  if(PDGcode==211){
    	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
    	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+1][charge]->Fill(val);
    	  }
    	}
    	if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
    	  if(PDGcode==321){
    	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
    	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+2][charge]->Fill(val);
    	  }
    	}
    	if (recoParticleArray[3].Contains(MCtrk)){ //Protons
    	  if(PDGcode==2212){
    	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
    	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+3][charge]->Fill(val);
    	  }
    	}
    	if (recoParticleArray[4].Contains(MCtrk)){ //K0s
    	  if(PDGcode==310){
    	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+4][0]->Fill(MCtrk->Eta(), MCtrk->Pt());
    	  }
    	}
   } // check physical primary loop end    

}  // MC array truth loop end (1638)

  PostData(1, fHistoList);

  
}
//-----------------------------------------------------------------

//void AliAnalysisTaskParticleEff::Terminate(Option_t *)
//{}
