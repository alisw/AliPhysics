class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliESDv0;
//class AliESDcascade;
class AliAODVertex;
class AliAODv0;
class AliAODcascade;

#include <Riostream.h>
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THistManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliESDcascade.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"

#include "AliAnalysisTaskStrVsMult.h"

#include "TMath.h"

ClassImp(AliAnalysisTaskStrVsMult)

AliAnalysisTaskStrVsMult::AliAnalysisTaskStrVsMult() : AliAnalysisTaskSE(),
//outputs
fHistos_eve(nullptr),
fHistos_K0S(nullptr),
fHistos_Lam(nullptr),
fHistos_ALam(nullptr),
fHistos_XiMin(nullptr),
fHistos_XiPlu(nullptr),
fHistos_OmMin(nullptr),
fHistos_OmPlu(nullptr),
//objects from the manager
fPIDResponse(0),
fTriggerMask(0),
//MC-related variables
fisMC(kFALSE),
fisMCassoc(kTRUE),
//default cuts configuration
fDefOnly(kFALSE),
fV0_Cuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
fCasc_Cuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//particle to be analysed
fParticleAnalysisStatus{true, true, true, true, true, true, true},
//variables for V0 cuts
fV0_DcaV0Daught(0),
fV0_DcaPosToPV(0),
fV0_DcaNegToPV(0),
fV0_V0CosPA(0),
fV0_V0Rad(0),
fV0_Pt(0),
fV0_yK0S(0),
fV0_yLam(0),
fV0_etaPos(0),
fV0_etaNeg(0),
fV0_InvMassK0s(0),
fV0_InvMassLam(0),
fV0_InvMassALam(0),
fV0_LeastCRaws(0),
fV0_LeastCRawsOvF(0),
fV0_NSigPosProton(0),
fV0_NSigPosPion(0),
fV0_NSigNegProton(0),
fV0_NSigNegPion(0),
fV0_DistOverTotP(0),
fV0_NegTOFBunchCrossing(0),
fV0_PosTOFBunchCrossing(0),
fV0_NegTrackStatus(0),
fV0_PosTrackStatus(0),
fV0_kinkidx(0),
//variables for Cascade analysis
fCasc_DcaCascDaught(0),
fCasc_CascCosPA(0),
fCasc_CascRad(0),
fCasc_etaPos(0),
fCasc_etaNeg(0),
fCasc_etaBac(0),
fCasc_kinkidx(0),
fCasc_NSigPosProton(0),
fCasc_NSigPosPion(0),
fCasc_NSigNegProton(0),
fCasc_NSigNegPion(0),
fCasc_NSigBacPion(0),
fCasc_NSigBacKaon(0),
fCasc_LeastCRaws(0),
fCasc_LeastCRawsOvF(0),
fCasc_InvMassLam(0),
fCasc_DcaV0Daught(0),
fCasc_V0CosPA(0),
fCasc_DcaV0ToPV(0),
fCasc_DcaBachToPV(0),
fCasc_PosTOFBunchCrossing(0),
fCasc_NegTOFBunchCrossing(0),
fCasc_BacTOFBunchCrossing(0),
fCasc_yXi(0),
fCasc_yOm(0),
fCasc_charge(0),
fCasc_Pt(0),
fCasc_DistOverTotP(0),
fCasc_InvMassXiMin(0),
fCasc_InvMassXiPlu(0),
fCasc_InvMassOmMin(0),
fCasc_InvMassOmPlu(0),
fCasc_V0Rad(0),
fCasc_DcaPosToPV(0),
fCasc_DcaNegToPV(0),
fCasc_NegTrackStatus(0),
fCasc_PosTrackStatus(0),
fCasc_BacTrackStatus(0),
fCasc_BacBarCosPA(0)
{
  //default constructor
}

AliAnalysisTaskStrVsMult::AliAnalysisTaskStrVsMult(const char *name, TString lExtraOptions) : AliAnalysisTaskSE(name),
//outputs
fHistos_eve(nullptr),
fHistos_K0S(nullptr),
fHistos_Lam(nullptr),
fHistos_ALam(nullptr),
fHistos_XiMin(nullptr),
fHistos_XiPlu(nullptr),
fHistos_OmMin(nullptr),
fHistos_OmPlu(nullptr),
//objects from the manager
fPIDResponse(0),
fTriggerMask(0),
//MC-related variables
fisMC(kFALSE),
fisMCassoc(kTRUE),
//default cuts configuration
fDefOnly(kFALSE),
fV0_Cuts{1., 0.11, 0.11, 0.97, 1., 0.5, 0.8, 70., 0.8, 5., 20., 30., -95.},
fCasc_Cuts{1., 0.99, 1., 4., 80., 0.8, 0.005, 1., 0.99, 0.1, 0.1, -95., 0.5, 0.8, 3., 3., 3., 0.2, 0.2, 0.99999},
//particle to be analysed
fParticleAnalysisStatus{true, true, true, true, true, true, true},
//variables for V0 cuts
fV0_DcaV0Daught(0),
fV0_DcaPosToPV(0),
fV0_DcaNegToPV(0),
fV0_V0CosPA(0),
fV0_V0Rad(0),
fV0_Pt(0),
fV0_yK0S(0),
fV0_yLam(0),
fV0_etaPos(0),
fV0_etaNeg(0),
fV0_InvMassK0s(0),
fV0_InvMassLam(0),
fV0_InvMassALam(0),
fV0_LeastCRaws(0),
fV0_LeastCRawsOvF(0),
fV0_NSigPosProton(0),
fV0_NSigPosPion(0),
fV0_NSigNegProton(0),
fV0_NSigNegPion(0),
fV0_DistOverTotP(0),
fV0_NegTOFBunchCrossing(0),
fV0_PosTOFBunchCrossing(0),
fV0_NegTrackStatus(0),
fV0_PosTrackStatus(0),
fV0_kinkidx(0),
//variables for Cascade analysis
fCasc_DcaCascDaught(0),
fCasc_CascCosPA(0),
fCasc_CascRad(0),
fCasc_etaPos(0),
fCasc_etaNeg(0),
fCasc_etaBac(0),
fCasc_kinkidx(0),
fCasc_NSigPosProton(0),
fCasc_NSigPosPion(0),
fCasc_NSigNegProton(0),
fCasc_NSigNegPion(0),
fCasc_NSigBacPion(0),
fCasc_NSigBacKaon(0),
fCasc_LeastCRaws(0),
fCasc_LeastCRawsOvF(0),
fCasc_InvMassLam(0),
fCasc_DcaV0Daught(0),
fCasc_V0CosPA(0),
fCasc_DcaV0ToPV(0),
fCasc_DcaBachToPV(0),
fCasc_PosTOFBunchCrossing(0),
fCasc_NegTOFBunchCrossing(0),
fCasc_BacTOFBunchCrossing(0),
fCasc_yXi(0),
fCasc_yOm(0),
fCasc_charge(0),
fCasc_Pt(0),
fCasc_DistOverTotP(0),
fCasc_InvMassXiMin(0),
fCasc_InvMassXiPlu(0),
fCasc_InvMassOmMin(0),
fCasc_InvMassOmPlu(0),
fCasc_V0Rad(0),
fCasc_DcaPosToPV(0),
fCasc_DcaNegToPV(0),
fCasc_NegTrackStatus(0),
fCasc_PosTrackStatus(0),
fCasc_BacTrackStatus(0),
fCasc_BacBarCosPA(0)
{
  //setting default cuts
  SetDefCutVals(); 
  SetDefCutVariations();
  //setting default centrality binning
  double centbins[15] = {0, 0.01, 0.1, 1., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90.};
  for (int ipart=0; ipart<knumpart; ipart++) {
    SetCentbinning(ipart, 14, centbins);
  }
  //setting default mass binning
  int massbins[4] = {90, 50, 80, 80};
  double minmass[4] = {0.41, 1.07, 1.28, 1.63};
  double maxmass[4] = {0.59, 1.17, 1.36, 1.71};
  for (int ipart=0; ipart<knumpart; ipart++) {
    SetMassbinning(ipart, massbins[ipart], minmass[ipart], maxmass[ipart]);
  }
  //setting default pt binning
  double ptbins[4][251];
  int nptbins[4] = {250, 250, 75, 75};
  double maxpt[4] = {25., 25., 15., 15.};
  for (int ipart=0; ipart<knumpart; ipart++) {
    for (int ipt=0; ipt<nptbins[ipart]+1; ipt++) {
      ptbins[ipart][ipt] = ipt*maxpt[ipart]/nptbins[ipart];
    }
    SetPtbinning(ipart, nptbins[ipart], ptbins[ipart]);
  }
  //Standard output
  DefineOutput(1, TList::Class()); // Event Histograms
  DefineOutput(2, TList::Class()); // K0S Histograms
  DefineOutput(3, TList::Class()); // Lambda Histograms
  DefineOutput(4, TList::Class()); // AntiLambda Histograms
  DefineOutput(5, TList::Class()); // XiMinus Histograms
  DefineOutput(6, TList::Class()); // XiPlus Histograms
  DefineOutput(7, TList::Class()); // OmegaMinus Histograms
  DefineOutput(8, TList::Class()); // OmegaPlus Histograms
}

AliAnalysisTaskStrVsMult::~AliAnalysisTaskStrVsMult()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    //Destroy output objects if present
//     if (fListHist) {
//         delete fListHist;
//         fListHist = 0x0;
//     }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::UserCreateOutputObjects()
{
  //histograms for event variables
  fHistos_eve = new THistManager("histos_eve");

  fHistos_eve->CreateTH1("hcent", "", 100, 0, 100, "s");  //storing #events in bins of centrality
  fHistos_eve->CreateTH1("henum", "", 1, 0, 1);  //storing total #events

  //histograms for V0 variables
  if (fParticleAnalysisStatus[kk0s]) {
    fHistos_K0S = new THistManager("histos_K0S");
    if(fisMC) fHistos_K0S->CreateTH2("h2_gen", "", fnptbins[kK0s], fptbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
    fHistos_K0S->CreateTH3("h3_ptmasscent_def", "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
  }
  if (fParticleAnalysisStatus[klam]) {
    fHistos_Lam = new THistManager("histos_Lam");
    fHistos_ALam = new THistManager("histos_ALam");
    fHistos_Lam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    fHistos_ALam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    if(fisMC) {
        fHistos_Lam->CreateTH2("h2_gen", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
        fHistos_Lam->CreateTH3("h3_FDmtxNUM_def", "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
        fHistos_Lam->CreateTH2("h2_FDmtxDEN_def", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
        fHistos_ALam->CreateTH2("h2_gen", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
        fHistos_ALam->CreateTH3("h3_FDmtxNUM_def", "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
        fHistos_ALam->CreateTH2("h2_FDmtxDEN_def", "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
    }
  }
  if (!fDefOnly && (fParticleAnalysisStatus[kk0s] || fParticleAnalysisStatus[klam])) {
    for (int icut = 0; icut<kV0cutsnum; icut++) {
      if (icut==kV0_y || icut==kV0_etaDaugh || icut==kV0_TOFBunchCrossing) continue;
      for (int ivar = 0; ivar<nvarcut_V0[icut]; ivar++) {
        if (fParticleAnalysisStatus[kk0s]) {
          if (icut!=kV0_PropLifetLam) fHistos_K0S->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
        }
        if (fParticleAnalysisStatus[klam]) {
          if (icut!=kV0_PropLifetK0s) {
            fHistos_Lam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
            fHistos_ALam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
            if(fisMC){
              fHistos_Lam->CreateTH3(Form("h3_FDmtxNUM[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
              fHistos_Lam->CreateTH2(Form("h2_FDmtxDEN[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
              fHistos_ALam->CreateTH3(Form("h3_FDmtxNUM[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnptbins[kXi], fptbinning[kXi], fncentbins[kLam], fcentbinning[kLam]);
              fHistos_ALam->CreateTH2(Form("h2_FDmtxDEN[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
            }
          }
        }
      }
    }
  }
  //histograms for Cascade variables
  if (fParticleAnalysisStatus[kxip]) {
    fHistos_XiMin = new THistManager("histos_XiMin");
    fHistos_XiPlu = new THistManager("histos_XiPlu");
    fHistos_XiMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
    fHistos_XiPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
    if(fisMC) fHistos_XiMin->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
    if(fisMC) fHistos_XiPlu->CreateTH2("h2_gen", "", fnptbins[kXi], fptbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
  }
  if (fParticleAnalysisStatus[komp]) {
    fHistos_OmMin = new THistManager("histos_OmMin");
    fHistos_OmPlu = new THistManager("histos_OmPlu");
    fHistos_OmMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
    fHistos_OmPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
    if(fisMC) fHistos_OmMin->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
    if(fisMC) fHistos_OmPlu->CreateTH2("h2_gen", "", fnptbins[kOm], fptbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
  }
  if (!fDefOnly && (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp])) {
    for (int icut=0; icut<kCasccutsnum; icut++) {
      if (icut==kCasc_y || icut==kCasc_etaDaugh || icut==kCasc_TOFBunchCrossing) continue;
      for (int ivar = 0; ivar<nvarcut_Casc[icut]; ivar++) {
        if (fParticleAnalysisStatus[kxip]) {
          if(icut!=kCasc_PropLifetOm) fHistos_XiMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
          if(icut!=kCasc_PropLifetOm) fHistos_XiPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        }
        if (fParticleAnalysisStatus[komp]) {
          if(icut!=kCasc_PropLifetXi) fHistos_OmMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
          if(icut!=kCasc_PropLifetXi) fHistos_OmPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        }
      }
    }
  }
  // PID Setup
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();

  //Output posting
  DataPosting();

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::UserExec(Option_t *)
{
  //get event from the input handler and cast it into the desired type of event
  AliVEvent *lVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!lVevent) { 
    AliWarning("ERROR: event not available \n"); 
    DataPosting(); 
    return; 
  }
  bool isESD = kFALSE;
  AliESDEvent *lESDevent = 0x0;
  AliAODEvent *lAODevent = 0x0;
  if (lVevent->InheritsFrom("AliESDEvent")) { 
    isESD = kTRUE; 
    lESDevent = dynamic_cast<AliESDEvent*>(InputEvent()); 
  } else if (lVevent->InheritsFrom("AliAODEvent")) {
    lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  }

  //get MC event
  AliMCEvent  *lMCev  = 0x0;
  if(fisMC){
    lMCev = MCEvent();
    if (!lMCev) {
      Printf("ERROR: Could not retrieve MC event in file %s\n",fInputHandler->GetTree()->GetCurrentFile()->GetName());
      DataPosting(); 
      return;
    }
  }

  //dumb histo for checking
  fHistos_eve->FillTH1("henum", 0.5);

  //get trigger information
  fTriggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  // Multiplicity Information 
  double lPercentile = -666;
  int lEvSelCode = -666;
  AliMultSelection *MultSelection = (AliMultSelection*) lVevent->FindListObject("MultSelection");
  if (!MultSelection) { 
    AliWarning("AliMultSelection object not found!"); 
    DataPosting(); 
    return; 
  } else {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    lEvSelCode = MultSelection->GetEvSelCode(); //==0 means event is good. Set by AliMultSelectionTask
  }

  //skip everything if event selection code !=0
  if(lEvSelCode != 0) { 
    DataPosting(); 
    return; 
  }

  //fill multiplicity/centrality percentile histogram
  if (((lPercentile>0. && lPercentile<10.) && fTriggerMask&AliVEvent::kCentral) || 
      ((lPercentile>30. && lPercentile<50.) && fTriggerMask&AliVEvent::kSemiCentral) || 
      (((lPercentile>10. && lPercentile<30.) || (lPercentile>50. && lPercentile<90.)) && fTriggerMask&AliVEvent::kINT7)) {
    fHistos_eve->FillTH1("hcent", lPercentile);
  } else {
    DataPosting(); 
    return; 
  }

  //get MC header and MC array
  AliAODMCHeader* header = 0x0;
  TClonesArray* MCTrackArray = 0x0;
  if(!isESD && fisMC) {
      header = static_cast<AliAODMCHeader*>(lAODevent->FindListObject(AliAODMCHeader::StdBranchName()));
      if (!header) {
        AliWarning("No header found.");
        DataPosting(); 
        return;
      } 
      MCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (MCTrackArray == NULL){
        AliWarning("No MC track array found.");
        DataPosting(); 
        return;
      }  
  }
  
  //MC truth
  if(fisMC){
    for (int i_MCtrk = 0;  i_MCtrk < lMCev->GetNumberOfTracks(); i_MCtrk++){
      AliVParticle *lPart;
      if(isESD) lPart = (AliMCParticle*) lMCev->GetTrack(i_MCtrk);
          else lPart = (AliAODMCParticle*) lMCev->GetTrack(i_MCtrk);
      bool isOOBpileup = kFALSE;
      if(isESD) isOOBpileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i_MCtrk, lMCev);
          else isOOBpileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i_MCtrk, header, MCTrackArray);
      if(!lPart || lPart->Y()<-0.5 || lPart->Y()>0.5 || !lPart->IsPhysicalPrimary() || isOOBpileup) continue;
      if(fParticleAnalysisStatus[kk0s] && lPart->PdgCode()==310   ) fHistos_K0S->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[klam] && lPart->PdgCode()==3122  ) fHistos_Lam->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[kalam]&& lPart->PdgCode()==-3122 ) fHistos_ALam->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[kxim] && lPart->PdgCode()==3312  ) fHistos_XiMin->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[kxip] && lPart->PdgCode()==-3312 ) fHistos_XiPlu->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[komm] && lPart->PdgCode()==3334  ) fHistos_OmMin->FillTH2("h2_gen",lPart->Pt(),lPercentile);
      if(fParticleAnalysisStatus[komp] && lPart->PdgCode()==-3334 ) fHistos_OmPlu->FillTH2("h2_gen",lPart->Pt(),lPercentile);
    }
  }

  //acquire best PV
  const AliVVertex *lBestPrimVtx = lVevent->GetPrimaryVertex();
  double lBestPV[3] = {-666., -666., -666.};
  lBestPrimVtx->GetXYZ(lBestPV);
  float lfBestPV[3];
  for (int i=0; i<3; i++) lfBestPV[i] = (float)lBestPV[i];

  //acquire magnetic field
  double lMagField = -666;
  lMagField = lVevent->GetMagneticField();

  //MC association
  int pdgPosDaught=0;
  int pdgNegDaught=0;
  int pdgBachelor=0;
  int pdgV0=0;
  int pdgCasc=0;
  bool assFlag[ksignednumpart]; //MC ass flags are true if associated and, by construction, always for data. They can be false only if fisMC && notassociated
  double fdmtx_ptxi = 0; //0 if no secondary lambda, value corresponding to generated-xi pT, - for lambda (from xim) and + for anti-lambda (from xip)

  // start v0 part if k=s or Lambda analysis is requested
  if (fParticleAnalysisStatus[kk0s] || fParticleAnalysisStatus[klam]) {
    int nv0s = 0;
    nv0s = lVevent->GetNumberOfV0s();
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {

      //MC association flag reset
      for(int iflag=0; iflag<kalam+1; iflag++) assFlag[iflag] = kTRUE;
      bool physprim = kTRUE;
      //feeddown lambda reset
      fdmtx_ptxi = 0;

      if (isESD) {
        // start ESD
        //get the iV0_th candidate in the event
        AliESDv0 *v0 = lESDevent->GetV0(iV0);
        if (!v0 || v0->GetOnFlyStatus()) continue;

        //get the position of the decay vertex
        double ldecayVtxV0[3];
        v0->GetXYZ(ldecayVtxV0[0], ldecayVtxV0[1], ldecayVtxV0[2]);
        fV0_V0Rad = TMath::Sqrt(ldecayVtxV0[0]*ldecayVtxV0[0]+ldecayVtxV0[1]*ldecayVtxV0[1]); //2D decay radius

        //get candidate's momentum
        double lpV0[3];
        v0->GetPxPyPz(lpV0[0], lpV0[1], lpV0[2]);
        double ltotpV0 = TMath::Sqrt(lpV0[0]*lpV0[0]+lpV0[1]*lpV0[1]+lpV0[2]*lpV0[2]); //total momentum

        //pt and rapidity
        fV0_Pt = v0->Pt();
        fV0_yK0S = v0->RapK0Short();
        fV0_yLam = v0->RapLambda();

        //retrieve daughter ESDTracks
        AliESDtrack *pTrack = lESDevent->GetTrack((UInt_t) TMath::Abs(v0->GetPindex()));
        AliESDtrack *nTrack = lESDevent->GetTrack((UInt_t) TMath::Abs(v0->GetNindex()));
        if (!pTrack || !nTrack) { 
          printf("ERROR: Could not retrieve one of the daughter tracks"); 
          continue; 
        }

        //Daughters' Eta
        fV0_etaPos = pTrack->Eta();
        fV0_etaNeg = nTrack->Eta();

        //preliminary checks
        if(((int) pTrack->GetSign()) == ((int) nTrack->GetSign())) continue; // remove like-sign V0s (if any)
        if(pTrack->GetTPCNclsF()<=0 || (int)nTrack->GetTPCNclsF()<=0) continue; //check here to avoid division by zero later

        //GetKinkIndex condition should be 0 for both if none of them is kink (disebled now)
        // fV0_kinkidx = pTrack->GetKinkIndex(0)+nTrack->GetKinkIndex(0);
        fV0_kinkidx = 0;

        //crossed raws
        double_t lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
        double_t lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
        fV0_LeastCRaws = (int) TMath::Min(lCrosRawsPos, lCrosRawsNeg);
        //crossed raws / Findable clusters
        double_t lCrosRawsOvFPos = lCrosRawsPos / ((double) (pTrack->GetTPCNclsF()));
        double_t lCrosRawsOvFNeg = lCrosRawsNeg / ((double) (nTrack->GetTPCNclsF()));
        fV0_LeastCRawsOvF = TMath::Min(lCrosRawsOvFPos, lCrosRawsOvFNeg);

        //dca info
        fV0_DcaPosToPV  = TMath::Abs(pTrack->GetD(lBestPV[0], lBestPV[1], lMagField));
        fV0_DcaNegToPV  = TMath::Abs(nTrack->GetD(lBestPV[0], lBestPV[1], lMagField));
        fV0_DcaV0Daught = v0->GetDcaV0Daughters();

        //cosPA
        fV0_V0CosPA = v0->GetV0CosineOfPointingAngle(lBestPV[0], lBestPV[1], lBestPV[2]);

        // invariant mass with different hypotheses
        v0->ChangeMassHypothesis(310);
        fV0_InvMassK0s = v0->GetEffMass();
        v0->ChangeMassHypothesis(3122);
        fV0_InvMassLam = v0->GetEffMass();
        v0->ChangeMassHypothesis(-3122);
        fV0_InvMassALam = v0->GetEffMass();

        //PID
        fV0_NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
        fV0_NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
        fV0_NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
        fV0_NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);

        //distance over total momentum
        fV0_DistOverTotP = TMath::Sqrt(TMath::Power(ldecayVtxV0[0]-lBestPV[0], 2)+TMath::Power(ldecayVtxV0[1]-lBestPV[1], 2)+TMath::Power(ldecayVtxV0[2]-lBestPV[2], 2));
        fV0_DistOverTotP /= (ltotpV0+1e-10); //avoid division by zero

        //TOF matching
        fV0_NegTOFBunchCrossing = nTrack->GetTOFBunchCrossing(lMagField);
        fV0_PosTOFBunchCrossing = pTrack->GetTOFBunchCrossing(lMagField);

        //track status: ( fV0_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
        fV0_NegTrackStatus = nTrack->GetStatus();
        fV0_PosTrackStatus = pTrack->GetStatus();

        //MC association
        if(fisMC){
          pdgPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrack->GetLabel())))-> PdgCode();
          pdgNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrack->GetLabel())))-> PdgCode();
          int labMothPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrack->GetLabel())))->GetMother();
          int labMothNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrack->GetLabel())))->GetMother();
          pdgV0 = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))-> PdgCode();
          physprim = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->IsPhysicalPrimary();
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=310)) assFlag[kk0s] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122)) assFlag[klam] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122)) assFlag[kalam] = kFALSE;
          //Feeddown matrix
          if(fParticleAnalysisStatus[klam] && !physprim){
            int labMothV0 = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
            int pdg_xi = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->PdgCode();
            double rap_xi = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Y();
            if(pdg_xi==3312 && TMath::Abs(rap_xi)<0.5) fdmtx_ptxi = -((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Pt();
              else if(pdg_xi==-3312 && TMath::Abs(rap_xi)<0.5) fdmtx_ptxi = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Pt();
          }
        }
      } else {
        //start AOD
        //get the iV0_th candidate in the event
        AliAODv0 *v0 = lAODevent->GetV0(iV0);
        if (!v0 || v0->GetOnFlyStatus()) continue;

        //get the 2D radius of the decay vertex
        fV0_V0Rad = v0->RadiusV0();

        //pt and rapidity
        fV0_Pt = v0->Pt();
        fV0_yK0S = v0->RapK0Short();
        fV0_yLam = v0->RapLambda();

        //retrieve daughter AODTracks
        AliAODTrack *pTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
        AliAODTrack *nTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);
        if (!pTrack || !nTrack) { 
          AliWarning("ERROR: Could not retrieve one of the daughter tracks"); 
          continue; 
        }

        //Daughters' Eta
        fV0_etaPos = pTrack->Eta();
        fV0_etaNeg = nTrack->Eta();

        //preliminary checks
        if(((int) pTrack->GetSign()) == ((int) nTrack->GetSign())) continue; // remove like-sign V0s (if any)

        if(pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0) continue; //check here to avoid division by zero later

        //GetKinkIndex condition should be 0 for both if none of them is kink (disabled for now)
        // AliAODVertex* pv = pTrack->GetProdVertex();
        // AliAODVertex* nv = nTrack->GetProdVertex();
        // fV0_kinkidx = (pv && (pv->GetType() & AliAODVertex::kKink)) + (nv && (nv->GetType() & AliAODVertex::kKink));
        fV0_kinkidx = 0;

        //crossed raws
        double_t lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
        double_t lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
        fV0_LeastCRaws = (int) TMath::Min(lCrosRawsPos, lCrosRawsNeg);
        //crossed raws / Findable clusters
        double_t lCrosRawsOvFPos = lCrosRawsPos/((double) (pTrack->GetTPCNclsF()));
        double_t lCrosRawsOvFNeg = lCrosRawsNeg/((double) (nTrack->GetTPCNclsF()));
        fV0_LeastCRawsOvF = TMath::Min(lCrosRawsOvFPos, lCrosRawsOvFNeg);

        //dca info
        fV0_DcaPosToPV = v0->DcaPosToPrimVertex();
        fV0_DcaNegToPV = v0->DcaNegToPrimVertex();
        fV0_DcaV0Daught = v0->DcaV0Daughters();

        //cosPA
        fV0_V0CosPA = v0->CosPointingAngle(lBestPV);

        // invariant mass with different hypotheses
        fV0_InvMassK0s = v0->MassK0Short();
        fV0_InvMassLam = v0->MassLambda();
        fV0_InvMassALam = v0->MassAntiLambda();

        //PID
        fV0_NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
        fV0_NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
        fV0_NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
        fV0_NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);

        //distance over total momentum
        fV0_DistOverTotP = v0->DecayLengthV0(lBestPV)/(v0->P()+1e-10);//avoid division by zero

        //TOF matching
        fV0_NegTOFBunchCrossing = nTrack->GetTOFBunchCrossing(lMagField);
        fV0_PosTOFBunchCrossing = pTrack->GetTOFBunchCrossing(lMagField);

        //track status: ( fV0_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
        fV0_NegTrackStatus = nTrack->GetStatus();
        fV0_PosTrackStatus = pTrack->GetStatus();

        if(fisMC){
          //MC association
          pdgPosDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrack->GetLabel())))-> GetPdgCode();
          pdgNegDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrack->GetLabel())))-> GetPdgCode();
          int labMothPosDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrack->GetLabel())))->GetMother();
          int labMothNegDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrack->GetLabel())))->GetMother();
          pdgV0 = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))-> GetPdgCode();
          physprim = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->IsPhysicalPrimary();
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122)) assFlag[klam] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122)) assFlag[kalam] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=310)) assFlag[kk0s] = kFALSE;
          //Feeddown matrix
          if(fParticleAnalysisStatus[klam] && !physprim){
            int labMothV0 = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
            int pdg_xi = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->GetPdgCode();
            double rap_xi = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Y();
            if(pdg_xi==3312 && TMath::Abs(rap_xi)<0.5) fdmtx_ptxi = -((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Pt();
              else if(pdg_xi==-3312 && TMath::Abs(rap_xi)<0.5) fdmtx_ptxi = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothV0)))->Pt();
          }
        }
      }
      //filling with default cuts
      if (fParticleAnalysisStatus[kk0s]) {
        if ( physprim && assFlag[kk0s] && ApplyCuts(kk0s)) fHistos_K0S->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassK0s, lPercentile);
      }
      if (fParticleAnalysisStatus[klam]) {
        if ( physprim && assFlag[klam] && ApplyCuts(klam)) fHistos_Lam->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassLam, lPercentile);
        if ( physprim && assFlag[kalam] && ApplyCuts(kalam)) fHistos_ALam->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassALam, lPercentile);
        if(fisMC){
          //Feeddown matrix filling
          //numerator
          if ( fdmtx_ptxi<0  && ApplyCuts(klam)) fHistos_Lam->FillTH3("h3_FDmtxNUM_def", fV0_Pt, TMath::Abs(fdmtx_ptxi), lPercentile);
          else if( fdmtx_ptxi>0 && ApplyCuts(kalam)) fHistos_ALam->FillTH3("h3_FDmtxNUM_def", fV0_Pt, TMath::Abs(fdmtx_ptxi), lPercentile);
          //denominator
          if ( assFlag[klam] && ApplyCuts(klam)) fHistos_Lam->FillTH2("h2_FDmtxDEN_def", fV0_Pt, lPercentile);
          if ( assFlag[kalam] && ApplyCuts(kalam)) fHistos_ALam->FillTH2("h2_FDmtxDEN_def", fV0_Pt, lPercentile);
        }
      }

      //filling 3D histograms
      if (!fDefOnly) FillHistCutVariations(kFALSE, lPercentile, physprim, assFlag, fdmtx_ptxi);
      
    } // end of V0 loop
  }

  // start cascade part if xi or omega analysis is requested
  if (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]) {
    int ncasc = 0;
    ncasc = lVevent->GetNumberOfCascades();
    for (int i_casc = 0; i_casc < ncasc; i_casc++) {

      //MC ass flag reset
      for(int iflag=kxip; iflag<ksignednumpart; iflag++) assFlag[iflag] = kTRUE;
      bool physprim = kTRUE;
      
      if(isESD){
        //start ESD
        //get the i_casc_th candidate in the event
        AliESDcascade *casc = lESDevent->GetCascade(i_casc);
        if (!casc) continue;

        // by default we start from the hypothesis that the cascade is a Xi-
        double quality = 0.;
        casc->ChangeMassHypothesis(quality, 3312);

        //cascade and V0 vertex positions (and calculate radii)
        double lVtxCasc[3], lVtxV0[3];
        casc->GetXYZcascade(lVtxCasc[0], lVtxCasc[1], lVtxCasc[2]);
        casc->GetXYZ(lVtxV0[0], lVtxV0[1], lVtxV0[2]); //Note that GetXYZ() is inherited from AliESDv0 and it returns the decay position of the V0 (and not of the cascade)
        fCasc_CascRad = TMath::Sqrt(lVtxCasc[0]*lVtxCasc[0]+lVtxCasc[1]*lVtxCasc[1]);
        fCasc_V0Rad = TMath::Sqrt(lVtxV0[0]*lVtxV0[0]+lVtxV0[1]*lVtxV0[1]);

        //get daughter tracks (positive, negative and bachelor)
        AliESDtrack *pTrackCasc = lESDevent->GetTrack((UInt_t) TMath::Abs(casc->GetPindex()));
        AliESDtrack *nTrackCasc = lESDevent->GetTrack((UInt_t) TMath::Abs(casc->GetNindex()));
        AliESDtrack *bTrackCasc = lESDevent->GetTrack((UInt_t) TMath::Abs(casc->GetBindex()));
        if (!pTrackCasc || !nTrackCasc || !bTrackCasc ) { 
          AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ..."); 
          continue; 
        }

        if(pTrackCasc->GetTPCNclsF()<=0 || nTrackCasc->GetTPCNclsF()<=0 || bTrackCasc->GetTPCNclsF()<=0) continue; //check here to avoid division by zero later

        //daughters' etas
        fCasc_etaPos = pTrackCasc->Eta();
        fCasc_etaNeg = nTrackCasc->Eta();
        fCasc_etaBac = bTrackCasc->Eta();

        //PID
        fCasc_NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
        fCasc_NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
        fCasc_NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
        fCasc_NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
        fCasc_NSigBacPion = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
        fCasc_NSigBacKaon = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);

        //crossed raws
        double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
        double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
        double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
        fCasc_LeastCRaws = (int) (lCrosRawsPos<lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
        //crossed raws / Findable clusters
        double lCrosRawsOvFPos = lCrosRawsPos/((double) (pTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFNeg = lCrosRawsNeg/((double) (nTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFBac = lCrosRawsBac/((double) (bTrackCasc->GetTPCNclsF()));
        fCasc_LeastCRawsOvF = lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);

        //V0 daughter mass (later to be checked against nominal)
        fCasc_InvMassLam = casc->GetEffMass(); //Note that GetEffMass() is inherited from AliESDv0 and it returns the mass of the V0 (and not of the cascade)

        //DCA info
        fCasc_DcaCascDaught = casc->GetDcaXiDaughters();
        fCasc_DcaBachToPV = TMath::Abs(bTrackCasc->GetD(lBestPV[0], lBestPV[1], lMagField));
        fCasc_DcaPosToPV = TMath::Abs(pTrackCasc->GetD(lBestPV[0], lBestPV[1], lMagField));
        fCasc_DcaNegToPV = TMath::Abs(nTrackCasc->GetD(lBestPV[0], lBestPV[1], lMagField));
        fCasc_DcaV0Daught = casc->GetDcaV0Daughters();
        fCasc_DcaV0ToPV = casc->GetD(lBestPV[0], lBestPV[1], lBestPV[2]);

        //cascade and V0 cosine of pointing angle
        fCasc_CascCosPA = casc->GetCascadeCosineOfPointingAngle(lBestPV[0], lBestPV[1], lBestPV[2]);
        fCasc_V0CosPA = casc->GetV0CosineOfPointingAngle(lBestPV[0], lBestPV[1], lBestPV[2]);

        //TOF matching
        fCasc_PosTOFBunchCrossing = pTrackCasc->GetTOFBunchCrossing(lMagField);
        fCasc_NegTOFBunchCrossing = nTrackCasc->GetTOFBunchCrossing(lMagField);
        fCasc_BacTOFBunchCrossing = bTrackCasc->GetTOFBunchCrossing(lMagField);

        //track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
        fCasc_NegTrackStatus = nTrackCasc->GetStatus();
        fCasc_PosTrackStatus = pTrackCasc->GetStatus();
        fCasc_BacTrackStatus = bTrackCasc->GetStatus();

        //candidate's rapidity (mass hypothesis dependent)
        fCasc_yXi = casc->RapXi();
        fCasc_yOm = casc->RapOmega();

        //charge
        fCasc_charge = (int) casc->Charge();

        //momentum and transverse momentum
        double lpCasc[3]; //momentum components
        casc->GetPxPyPz(lpCasc[0], lpCasc[1], lpCasc[2]);
        double ltotpCasc = TMath::Sqrt(lpCasc[0]*lpCasc[0]+lpCasc[1]*lpCasc[1]+lpCasc[2]*lpCasc[2]); //total momentum
        fCasc_Pt = casc->Pt();

        //distance over total momentum
        fCasc_DistOverTotP = TMath::Sqrt(TMath::Power(lVtxCasc[0]-lBestPV[0], 2)+TMath::Power(lVtxCasc[1]-lBestPV[1], 2)+TMath::Power(lVtxCasc[2]-lBestPV[2], 2));
        fCasc_DistOverTotP/=(ltotpCasc+1e-10); //avoid division by zero

        //candidate's invariant mass
        casc->ChangeMassHypothesis(quality, 3312); //Xi-
        fCasc_InvMassXiMin = casc->M();
        casc->ChangeMassHypothesis(quality, 3334); //Om-
        fCasc_InvMassOmMin = casc->M();
        casc->ChangeMassHypothesis(quality, -3312); //Xi+
        fCasc_InvMassXiPlu = casc->M();
        casc->ChangeMassHypothesis(quality, -3334); //Om+
        fCasc_InvMassOmPlu = casc->M();

        //calculate CosPA Bachelor-Baryon to remove "bump" structure in InvMass
        AliAnalysisTaskESDfilter esdFilter;
        if(casc->Charge()<0.){ 
          fCasc_BacBarCosPA = esdFilter.GetCosPA(pTrackCasc, bTrackCasc, lMagField, lfBestPV); 
        } else { 
          fCasc_BacBarCosPA = esdFilter.GetCosPA(nTrackCasc, bTrackCasc, lMagField, lfBestPV); 
        }
        
        //MC association
        if(fisMC){
          pdgPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))-> PdgCode();
          pdgNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))-> PdgCode();
          pdgBachelor = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))-> PdgCode();
          int labMothPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->GetMother();
          int labMothNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->GetMother();
          int labMothV0 = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
          int labMothBach = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->GetMother();
          pdgV0 = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))-> PdgCode();
          pdgCasc = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothBach)))-> PdgCode();
          physprim = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothBach)))->IsPhysicalPrimary();
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || pdgBachelor!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122 || labMothV0!=labMothBach || pdgCasc!=3312)) assFlag[kxim] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || pdgBachelor!=211 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122 || labMothV0!=labMothBach || pdgCasc!=-3312)) assFlag[kxip] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || pdgBachelor!=-321 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122 || labMothV0!=labMothBach || pdgCasc!=3334)) assFlag[komm] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || pdgBachelor!=321 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122 || labMothV0!=labMothBach || pdgCasc!=-3334)) assFlag[komp] = kFALSE;
        }

      } else {
        //start AOD part
        AliAODcascade *casc = lAODevent->GetCascade(i_casc);
        if (!casc) continue;

        //cascade and V0 2D radii
        double lVtxCasc[3];
        lVtxCasc[0] = casc->DecayVertexXiX();
        lVtxCasc[1] = casc->DecayVertexXiY();
        lVtxCasc[2] = casc->DecayVertexXiZ();
        fCasc_CascRad = TMath::Sqrt(lVtxCasc[0]*lVtxCasc[0]+lVtxCasc[1]*lVtxCasc[1]);
        fCasc_V0Rad = casc->RadiusSecVtx();

        //get daughter tracks (positive, negative and bachelor)
        AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(0));
        AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(1));
        AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDecayVertexXi()->GetDaughter(0));
        if (!pTrackCasc || !nTrackCasc || !bTrackCasc) { 
          AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n"); 
          continue; 
        }

        //preliminary check
        if(pTrackCasc->GetTPCNclsF()<=0 || nTrackCasc->GetTPCNclsF()<=0 || bTrackCasc->GetTPCNclsF()<=0) continue; //check here to avoid division by zero later
        //daughters' etas
        fCasc_etaPos = pTrackCasc->Eta();
        fCasc_etaNeg = nTrackCasc->Eta();
        fCasc_etaBac = bTrackCasc->Eta();

        //PID
        fCasc_NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
        fCasc_NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
        fCasc_NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
        fCasc_NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
        fCasc_NSigBacPion = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
        fCasc_NSigBacKaon = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);

        //crossed raws
        double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
        double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
        double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
        fCasc_LeastCRaws = (int) (lCrosRawsPos<lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
        //crossed raws / Findable clusters
        double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()));
        fCasc_LeastCRawsOvF = lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);

        //DCA info
        fCasc_DcaCascDaught = casc->DcaXiDaughters();
        fCasc_DcaBachToPV = casc->DcaBachToPrimVertex();
        fCasc_DcaPosToPV = casc->DcaPosToPrimVertex();
        fCasc_DcaNegToPV = casc->DcaNegToPrimVertex();
        fCasc_DcaV0Daught = casc->DcaV0Daughters();
        fCasc_DcaV0ToPV = casc->DcaV0ToPrimVertex();

        //cascade and V0 cosine of pointing angle
        fCasc_CascCosPA = casc->CosPointingAngleXi((const Double_t&) lBestPV[0], (const Double_t&) lBestPV[1], (const Double_t&) lBestPV[2]);
        fCasc_V0CosPA = casc->CosPointingAngle(lBestPV);

        //TOF matching
        fCasc_PosTOFBunchCrossing = pTrackCasc->GetTOFBunchCrossing(lMagField);
        fCasc_NegTOFBunchCrossing = nTrackCasc->GetTOFBunchCrossing(lMagField);
        fCasc_BacTOFBunchCrossing = bTrackCasc->GetTOFBunchCrossing(lMagField);

        //track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
        fCasc_NegTrackStatus = nTrackCasc->GetStatus();
        fCasc_PosTrackStatus = pTrackCasc->GetStatus();
        fCasc_BacTrackStatus = bTrackCasc->GetStatus();

        //candidate's rapidity (mass hypothesis dependent)
        fCasc_yXi = casc->RapXi();
        fCasc_yOm = casc->RapOmega();

        //charge
        fCasc_charge = (int) casc->ChargeXi();

        //V0 daughter mass (later to be checked against nominal)
        if (fCasc_charge<0) {
          fCasc_InvMassLam = casc->MassLambda();
        } else {
          fCasc_InvMassLam = casc->MassAntiLambda();
        }

        //transverse momentum
        fCasc_Pt = TMath::Sqrt(casc->Pt2Xi());

        //distance over total momentum
        fCasc_DistOverTotP = (TMath::Sqrt(TMath::Power(lVtxCasc[0]-lBestPV[0], 2)+TMath::Power(lVtxCasc[1]-lBestPV[1], 2)+TMath::Power(lVtxCasc[2]-lBestPV[2], 2)))/(TMath::Sqrt(casc->Ptot2Xi())+1e-10);

        //candidate's invariant mass
        fCasc_InvMassXiMin = casc->MassXi();
        fCasc_InvMassXiPlu = casc->MassXi();
        fCasc_InvMassOmMin = casc->MassOmega();
        fCasc_InvMassOmPlu = casc->MassOmega();

        //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
        fCasc_BacBarCosPA = casc->BachBaryonCosPA();

        //MC association
        if(fisMC){
          pdgPosDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))-> PdgCode();
          pdgNegDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))-> PdgCode();
          pdgBachelor = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))-> PdgCode();
          int labMothPosDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->GetMother();
          int labMothNegDaught = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->GetMother();
          int labMothV0 = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
          int labMothBach = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->GetMother();
          pdgV0 = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))-> PdgCode();
          pdgCasc = ((AliAODMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothBach)))-> PdgCode();
          physprim = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothBach)))->IsPhysicalPrimary();
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || pdgBachelor!=-211 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122 || labMothV0!=labMothBach || pdgCasc!=3312)) assFlag[kxim] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || pdgBachelor!=211 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122 || labMothV0!=labMothBach || pdgCasc!=-3312)) assFlag[kxip] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=2212 || pdgNegDaught!=-211 || pdgBachelor!=-321 || labMothPosDaught!=labMothNegDaught || pdgV0!=3122 || labMothV0!=labMothBach || pdgCasc!=3334)) assFlag[komm] = kFALSE;
          if(fisMCassoc && (pdgPosDaught!=211 || pdgNegDaught!=-2212 || pdgBachelor!=321 || labMothPosDaught!=labMothNegDaught || pdgV0!=-3122 || labMothV0!=labMothBach || pdgCasc!=-3334)) assFlag[komp] = kFALSE;
        }

      }
      //fills TH3 with default cuts
      if (fParticleAnalysisStatus[kxip]) {
        if( physprim && assFlag[kxim] && ApplyCuts(kxim)) fHistos_XiMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiMin, lPercentile);
        if( physprim && assFlag[kxip] && ApplyCuts(kxip)) fHistos_XiPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiPlu, lPercentile);
      }
      if (fParticleAnalysisStatus[komp]) {
        if( physprim && assFlag[komm] && ApplyCuts(komm)) fHistos_OmMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmMin, lPercentile);
        if( physprim && assFlag[komp] && ApplyCuts(komp)) fHistos_OmPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmPlu, lPercentile);
      }

      //filling 3D histograms
      if (!fDefOnly) FillHistCutVariations(kTRUE, lPercentile,physprim,assFlag,0);
    }
  }

  DataPosting();

}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetCutVal(bool defchange, bool iscasc, int cutnum, double cval) {
  if (!iscasc) {
    cutval_V0[cutnum] = cval;
    if (defchange) fV0_Cuts[cutnum] = cval;
  } else {
    cutval_Casc[cutnum] = cval;
    if (defchange) fCasc_Cuts[cutnum] = cval;
  }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetDefCutVals() {
  //V0 part
  if (fParticleAnalysisStatus[kk0s] || fParticleAnalysisStatus[klam]) {
    for (int iV0_cut=0; iV0_cut<kV0cutsnum; iV0_cut++) {
      SetCutVal(kFALSE, kFALSE, iV0_cut, fV0_Cuts[iV0_cut]);
    }
  }
  //Cascade part
  if (fParticleAnalysisStatus[kxip] || fParticleAnalysisStatus[komp]) {
    for (int iCasc_cut=0; iCasc_cut<kCasccutsnum; iCasc_cut++) {
      SetCutVal(kFALSE, kTRUE, iCasc_cut, fCasc_Cuts[iCasc_cut]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetCutVariation(bool iscasc, int cutnum, int nvar, double lowval, double highval) {
  if (!iscasc) {
    nvarcut_V0[cutnum] = nvar;
    varlowcut_V0[cutnum] = lowval;
    varhighcut_V0[cutnum] = highval;
  } else {
    nvarcut_Casc[cutnum] = nvar;
    varlowcut_Casc[cutnum] = lowval;
    varhighcut_Casc[cutnum] = highval;
  }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetDefCutVariations() {
  SetCutVariation(kFALSE, kV0_DcaV0Daught, 10, 0.5, 1.4);
  SetCutVariation(kFALSE, kV0_DcaPosToPV, 11, 0.10, 0.13);
  SetCutVariation(kFALSE, kV0_DcaNegToPV, 11, 0.10, 0.13);
  SetCutVariation(kFALSE, kV0_V0CosPA, 11, 0.95, 0.999);
  SetCutVariation(kFALSE, kV0_V0Rad, 11, 0.9, 1.3);
  SetCutVariation(kFALSE, kV0_LeastCRaws, 11, 60, 80);
  SetCutVariation(kFALSE, kV0_LeastCRawsOvF, 11, 0.75, 0.90);
  SetCutVariation(kFALSE, kV0_NSigPID, 6, 2, 7);
  SetCutVariation(kFALSE, kV0_PropLifetK0s, 11, 10, 40);
  SetCutVariation(kFALSE, kV0_PropLifetLam, 11, 10, 40);

  SetCutVariation(kTRUE, kCasc_DcaCascDaught, 10, 0.5, 1.4);
  SetCutVariation(kTRUE, kCasc_CascCosPA, 21, 0.95, 0.999);
  SetCutVariation(kTRUE, kCasc_CascRad, 11, 0.5, 1.5);
  SetCutVariation(kTRUE, kCasc_NSigPID, 6, 2, 7);
  SetCutVariation(kTRUE, kCasc_LeastCRaws, 11, 70, 90);
  SetCutVariation(kTRUE, kCasc_LeastCRawsOvF, 11, 0.75, 0.9);
  SetCutVariation(kTRUE, kCasc_InvMassLam, 5, 0.002, 0.006);
  SetCutVariation(kTRUE, kCasc_DcaV0Daught, 10, 0.5, 1.4);
  SetCutVariation(kTRUE, kCasc_V0CosPA, 21, 0.95, 0.999);
  SetCutVariation(kTRUE, kCasc_DcaV0ToPV, 11, 0.05, 0.15);
  SetCutVariation(kTRUE, kCasc_DcaBachToPV, 11, 0.05, 0.15);
  SetCutVariation(kTRUE, kCasc_PropLifetXi, 7, 2, 5);
  SetCutVariation(kTRUE, kCasc_PropLifetOm, 7, 2, 5);
  SetCutVariation(kTRUE, kCasc_V0Rad, 11, 1., 5.);
  SetCutVariation(kTRUE, kCasc_DcaMesToPV, 11, 0.1, 0.3);
  SetCutVariation(kTRUE, kCasc_DcaBarToPV, 11, 0.1, 0.3);
  SetCutVariation(kTRUE, kCasc_BacBarCosPA, 21, 0.9999, 0.999999);
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetDefOnly(bool isdefonly) {
  fDefOnly = isdefonly;
}

//________________________________________________________________________
bool AliAnalysisTaskStrVsMult::ApplyCuts(int part) {

  if (part<kxip) { //we are checking V0s
    // check candidate's rapidity (particle hypothesis' dependent)
    if ((part==kk0s) && TMath::Abs(fV0_yK0S)>cutval_V0[kV0_y]) return kFALSE;
    if ((part==klam || part==kalam) && TMath::Abs(fV0_yLam)>cutval_V0[kV0_y]) return kFALSE;
    // check candidate daughters' pseudo-rapidity
    if (TMath::Abs(fV0_etaPos)>cutval_V0[kV0_etaDaugh] || TMath::Abs(fV0_etaNeg)>cutval_V0[kV0_etaDaugh]) return kFALSE;
    // check candidate daughters' crossed TPC raws (note that the checked value is the lowest between the two daughter)
    if (fV0_LeastCRaws<cutval_V0[kV0_LeastCRaws]) return kFALSE;
    // check candidate daughters' crossed TPC raws over findable
    if (fV0_LeastCRawsOvF<cutval_V0[kV0_LeastCRawsOvF]) return kFALSE;
    // check candidate daughters' DCA to Primary Vertex (needs to be large because V0 decay is far from the Primary Vertex)
    if (fV0_DcaPosToPV<cutval_V0[kV0_DcaPosToPV] || fV0_DcaNegToPV<cutval_V0[kV0_DcaNegToPV]) return kFALSE;
    // check candidate daughters' DCA between them (needs to be small because they have to come from the same secondary vertex)
    if (fV0_DcaV0Daught>cutval_V0[kV0_DcaV0Daught]) return kFALSE;
    // check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
    if (fV0_V0Rad<cutval_V0[kV0_V0Rad]) return kFALSE;
    // check the cosine of the Pointing Angle (angle between candidate's momentum and vector connecting Primary and secondary vertices)
    if (fV0_V0CosPA<cutval_V0[kV0_V0CosPA]) return kFALSE;
    // check PID for all daughters (particle hypothesis' dependent)
    if ((part==kk0s) && (TMath::Abs(fV0_NSigPosPion)>cutval_V0[kV0_NSigPID] || TMath::Abs(fV0_NSigNegPion)>cutval_V0[kV0_NSigPID])) return kFALSE;
    if ((part==klam) && (TMath::Abs(fV0_NSigPosProton)>cutval_V0[kV0_NSigPID] || TMath::Abs(fV0_NSigNegPion)>cutval_V0[kV0_NSigPID])) return kFALSE;
    if ((part==kalam) && (TMath::Abs(fV0_NSigNegProton)>cutval_V0[kV0_NSigPID] || TMath::Abs(fV0_NSigPosPion)>cutval_V0[kV0_NSigPID])) return kFALSE;
    // check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
    if ((part==kk0s) && (0.497*fV0_DistOverTotP>cutval_V0[kV0_PropLifetK0s])) return kFALSE;
    if ((part>kk0s) && (1.115683*fV0_DistOverTotP>cutval_V0[kV0_PropLifetLam])) return kFALSE;
    // check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
    if (!(fV0_NegTrackStatus & AliESDtrack::kITSrefit) && !(fV0_PosTrackStatus & AliESDtrack::kITSrefit) &&
        (fV0_NegTOFBunchCrossing < cutval_V0[kV0_TOFBunchCrossing]) && (fV0_PosTOFBunchCrossing < cutval_V0[kV0_TOFBunchCrossing])) return kFALSE;
    // TPC refit, should be already verified for Offline V0s
    if (!(fV0_PosTrackStatus & AliESDtrack::kTPCrefit) || !(fV0_NegTrackStatus & AliESDtrack::kTPCrefit)) return kFALSE;
    // check that none of daughters is a kink
    if (fV0_kinkidx>0) return kFALSE;

  } else { //we are checking cascades
    // check candidate's charge
    if((part==kxim || part==komm) && fCasc_charge>0) return kFALSE;
    if((part==kxip || part==komp) && fCasc_charge<0) return kFALSE;
    // check candidate's rapidity (particle hypothesis' dependent)
    if((part==kxip || part==kxim) && TMath::Abs(fCasc_yXi)>cutval_Casc[kCasc_y]) return kFALSE;
    if((part==komp || part==komm) && TMath::Abs(fCasc_yOm)>cutval_Casc[kCasc_y]) return kFALSE;
    // check candidate daughters' pseudo-rapidity
    if(TMath::Abs(fCasc_etaPos)>cutval_Casc[kCasc_etaDaugh] || TMath::Abs(fCasc_etaNeg)>cutval_Casc[kCasc_etaDaugh] || TMath::Abs(fCasc_etaBac)>cutval_Casc[kCasc_etaDaugh]) return kFALSE;
    // check candidate daughters' crossed TPC raws (note that the checked value is the lowest among the daughters)
    if(fCasc_LeastCRaws<cutval_Casc[kCasc_LeastCRaws]) return kFALSE;
    // check candidate daughters' crossed TPC raws over findable
    if(fCasc_LeastCRawsOvF<cutval_Casc[kCasc_LeastCRawsOvF]) return kFALSE;
    // check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
    if(fCasc_CascRad<cutval_Casc[kCasc_CascRad]) return kFALSE;
    // check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
    if(fCasc_V0Rad<cutval_Casc[kCasc_V0Rad]) return kFALSE;
    // check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
    if(fCasc_CascCosPA<cutval_Casc[kCasc_CascCosPA]) return kFALSE;
    if(fCasc_V0CosPA<cutval_Casc[kCasc_V0CosPA]) return kFALSE;
    // check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
    if(fCasc_DcaBachToPV<cutval_Casc[kCasc_DcaBachToPV]) return kFALSE;
    if(fCasc_DcaV0ToPV<cutval_Casc[kCasc_DcaV0ToPV]) return kFALSE;
    // check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
    if((part==kxim || part==komm) && (fCasc_DcaPosToPV<cutval_Casc[kCasc_DcaBarToPV] || fCasc_DcaNegToPV<cutval_Casc[kCasc_DcaMesToPV])) return kFALSE; // in this case pos=p, neg=pi-, bach=pi- or K-
    if((part==kxip || part==komp) && (fCasc_DcaPosToPV<cutval_Casc[kCasc_DcaMesToPV] || fCasc_DcaNegToPV<cutval_Casc[kCasc_DcaBarToPV])) return kFALSE; // in this case pos=pi+, neg=anti-p, bach=pi+ or K+
    // check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
    if(fCasc_DcaV0Daught>cutval_Casc[kCasc_DcaV0Daught]) return kFALSE;
    // check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
    if(fCasc_DcaCascDaught>cutval_Casc[kCasc_DcaCascDaught]) return kFALSE;
    // check candidate V0 daughter's mass difference from nominal Lambda mass
    if(TMath::Abs(fCasc_InvMassLam-1.115683)>cutval_Casc[kCasc_InvMassLam]) return kFALSE;
    // check that none of daughters is a kink
    if(fCasc_kinkidx>0) return kFALSE;
    // check if at least one of candidate's daughter has a hit in the TOF or has ITSrefit flag (removes Out Of Bunch Pileup)
    if(!(fCasc_NegTrackStatus & AliESDtrack::kITSrefit) &&
       !(fCasc_PosTrackStatus & AliESDtrack::kITSrefit) &&
       !(fCasc_BacTrackStatus & AliESDtrack::kITSrefit) &&
       (fCasc_NegTOFBunchCrossing < cutval_Casc[kCasc_TOFBunchCrossing]) &&
       (fCasc_PosTOFBunchCrossing < cutval_Casc[kCasc_TOFBunchCrossing]) &&
       (fCasc_BacTOFBunchCrossing < cutval_Casc[kCasc_TOFBunchCrossing])) return kFALSE; //OOB
    // TPC refit, should be already verified for Offline V0s
    if(!(fCasc_PosTrackStatus & AliESDtrack::kTPCrefit) ||
       !(fCasc_NegTrackStatus & AliESDtrack::kTPCrefit) ||
       !(fCasc_BacTrackStatus & AliESDtrack::kTPCrefit)) return kFALSE;
    // check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
    if((part==kxim || part==kxip) && ((1.32171*fCasc_DistOverTotP)>(4.91*cutval_Casc[kCasc_PropLifetXi]))) return kFALSE;   //4.91 is the ctau of xi in cm
    if((part==komm || part==komp) && ((1.67245*fCasc_DistOverTotP)>(2.461*cutval_Casc[kCasc_PropLifetOm]))) return kFALSE;   //2.461 is the ctau of om in cm
    // check DCA bachelor-baryon. If it is too small --> bump structure in Inv Mass
    if(fCasc_BacBarCosPA>cutval_Casc[kCasc_BacBarCosPA]) return kFALSE;
    // check PID for all daughters (particle hypothesis' dependent)
    if((part==kxip) && (TMath::Abs(fCasc_NSigPosPion)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigNegProton)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacPion)>cutval_Casc[kCasc_NSigPID])) return kFALSE;
    if((part==kxim) && (TMath::Abs(fCasc_NSigNegPion)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigPosProton)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacPion)>cutval_Casc[kCasc_NSigPID])) return kFALSE;
    if((part==komp) && (TMath::Abs(fCasc_NSigPosPion)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigNegProton)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacKaon)>cutval_Casc[kCasc_NSigPID])) return kFALSE;
    if((part==komm) && (TMath::Abs(fCasc_NSigNegPion)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigPosProton)>cutval_Casc[kCasc_NSigPID] || TMath::Abs(fCasc_NSigBacKaon)>cutval_Casc[kCasc_NSigPID])) return kFALSE;
  }
  return kTRUE; //survived!
}

void AliAnalysisTaskStrVsMult::SetParticleAnalysisStatus(bool k0s, bool lambda, bool xi, bool omega) {
  fParticleAnalysisStatus[kk0s] = k0s;
  fParticleAnalysisStatus[klam] = lambda;
  fParticleAnalysisStatus[kalam] = lambda;
  fParticleAnalysisStatus[kxip] = xi;
  fParticleAnalysisStatus[kxim] = xi;
  fParticleAnalysisStatus[komp] = omega;
  fParticleAnalysisStatus[komm] = omega;
}

bool AliAnalysisTaskStrVsMult::GetParticleAnalysisStatus(int part) {
  if (part<0 || part>=ksignednumpart) {
    ::Error("AliAnalysisTaskStrVsMult::GetParticleAnalysisStatus", "Wrong particle selected: accepted values from 0 to 3");
    return false;
  }
  return fParticleAnalysisStatus[part];
}

void AliAnalysisTaskStrVsMult::SetCentbinning(int ipart, int numcentbins, double *centbins) {
  fncentbins[ipart] = numcentbins;
  for(int i=0; i<fncentbins[ipart]+1; i++) {
    fcentbinning[ipart][i] = centbins[i];
  }
}

void AliAnalysisTaskStrVsMult::SetPtbinning(int ipart, int numptbins, double *ptbins) {
  fnptbins[ipart] = numptbins;
  for(int i=0; i<fnptbins[ipart]+1; i++) {
    fptbinning[ipart][i] = ptbins[i];
  }
}

void AliAnalysisTaskStrVsMult::SetMassbinning(int ipart, int nummassbins, double valminmass, double valmaxmass) {
  fnmassbins[ipart] = nummassbins;
  for(int i=0; i<fnmassbins[ipart]+1; i++) {
    fmassbinning[ipart][i] = valminmass+i*(valmaxmass-valminmass)/fnmassbins[ipart];
  }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::DataPosting() {

  PostData(1, fHistos_eve->GetListOfHistograms());
  //Histograms for analysed particle specie
  int histnumber=1;
  if (fParticleAnalysisStatus[kk0s]) {
    histnumber++;
    PostData(histnumber, fHistos_K0S->GetListOfHistograms());
  }
  if (fParticleAnalysisStatus[klam]) {
    histnumber = histnumber+2;
    PostData(histnumber-1, fHistos_Lam->GetListOfHistograms());
    PostData(histnumber, fHistos_ALam->GetListOfHistograms());
  }
  if (fParticleAnalysisStatus[kxip]) {
    histnumber = histnumber+2;
    PostData(histnumber-1, fHistos_XiMin->GetListOfHistograms());
    PostData(histnumber, fHistos_XiPlu->GetListOfHistograms());
  }
  if (fParticleAnalysisStatus[komp]) {
    histnumber = histnumber+2;
    PostData(histnumber-1, fHistos_OmMin->GetListOfHistograms());
    PostData(histnumber, fHistos_OmPlu->GetListOfHistograms());
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::FillHistCutVariations(bool iscasc, double perc, bool phypri, bool *associFlag, double ptassxi) {

  if (!iscasc) {
    for (int i_cut=0; i_cut<kV0cutsnum; i_cut++) {
      if(i_cut==kV0_y || i_cut==kV0_etaDaugh || i_cut==kV0_TOFBunchCrossing) continue;
      for(int i_var=0; i_var<nvarcut_V0[i_cut]; i_var++){
        //K0S filling
        if (i_cut!=kV0_PropLifetLam) {
          if (fParticleAnalysisStatus[kk0s]) {
            SetCutVal(kFALSE, kFALSE, i_cut, varlowcut_V0[i_cut]+i_var*(varhighcut_V0[i_cut]-varlowcut_V0[i_cut])/(nvarcut_V0[i_cut]-1));
            if (phypri && associFlag[kk0s] && ApplyCuts(kk0s)) fHistos_K0S->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassK0s, perc);
          }
        }
        //Lam and AntiLam filling
        if (i_cut!=kV0_PropLifetK0s) {
          if (fParticleAnalysisStatus[klam]) {
            SetCutVal(kFALSE, kFALSE, i_cut, varlowcut_V0[i_cut]+i_var*(varhighcut_V0[i_cut]-varlowcut_V0[i_cut])/(nvarcut_V0[i_cut]-1));
            if (phypri && associFlag[klam] && ApplyCuts(klam)) fHistos_Lam->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassLam, perc);
            if (phypri && associFlag[kalam] && ApplyCuts(kalam)) fHistos_ALam->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassALam, perc);
            if (fisMC) {
              //Feeddown matrix filling
              //numerator
              if (ptassxi<0 && ApplyCuts(klam)) fHistos_Lam->FillTH3(Form("h3_FDmtxNUM[%d][%d]", i_cut, i_var), fV0_Pt, TMath::Abs(ptassxi), perc);
              else if (ptassxi>0 && ApplyCuts(kalam)) fHistos_ALam->FillTH3(Form("h3_FDmtxNUM[%d][%d]", i_cut, i_var), fV0_Pt, TMath::Abs(ptassxi), perc);
              //denominator
              if (associFlag[klam] && ApplyCuts(klam)) fHistos_Lam->FillTH2(Form("h2_FDmtxDEN[%d][%d]", i_cut, i_var), fV0_Pt, perc);
              if (associFlag[kalam] && ApplyCuts(kalam)) fHistos_ALam->FillTH2(Form("h2_FDmtxDEN[%d][%d]", i_cut, i_var), fV0_Pt, perc);
            }
          }
        }
      }
      SetDefCutVals(); //reset defaults
    }
  } else {
    for(int i_cut=0; i_cut<kCasccutsnum; i_cut++) {
      if(i_cut==kCasc_y || i_cut==kCasc_etaDaugh || i_cut==kCasc_TOFBunchCrossing) continue;
      for (int i_var=0; i_var<nvarcut_Casc[i_cut]; i_var++) {
        //Xi filling
        if (i_cut!=kCasc_PropLifetOm) {
          if (fParticleAnalysisStatus[kxip]) {
            SetCutVal(kFALSE, kTRUE, i_cut, varlowcut_Casc[i_cut]+i_var*(varhighcut_Casc[i_cut]-varlowcut_Casc[i_cut])/(nvarcut_Casc[i_cut]-1));
            if(phypri && associFlag[kxim] && ApplyCuts(kxim)) fHistos_XiMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiMin, perc);
            if(phypri && associFlag[kxip] && ApplyCuts(kxip)) fHistos_XiPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiPlu, perc);
          }
        }
        //Om filling
        if (i_cut!=kCasc_PropLifetXi) {
          if (fParticleAnalysisStatus[komp]) {
            SetCutVal(kFALSE, kTRUE, i_cut, varlowcut_Casc[i_cut]+i_var*(varhighcut_Casc[i_cut]-varlowcut_Casc[i_cut])/(nvarcut_Casc[i_cut]-1));
            if(phypri && associFlag[komm] && ApplyCuts(komm)) fHistos_OmMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmMin, perc);
            if(phypri && associFlag[komp] && ApplyCuts(komp)) fHistos_OmPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmPlu, perc);
          }
        }
      }
      SetDefCutVals(); //reset defaults
    }
  }
}
