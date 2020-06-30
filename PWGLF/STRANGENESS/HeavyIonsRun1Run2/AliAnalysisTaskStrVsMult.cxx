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
//default cut conf-only?
fDefOnly(kFALSE),
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
fCasc_DcaBacBar(0)
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
//default cut conf-only?
fDefOnly(kFALSE),
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
fCasc_DcaBacBar(0)
{
  //setting default cuts
  SetDefCutVals();
  if (!fDefOnly) {
    SetDefCutVariations();
  }
  //setting default centrality binning
  double centbins[16] = {0, 0.01, 0.1, 1., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
  for (int ipart=0; ipart<knumpart; ipart++) {
    SetCentbinning(ipart, 15, centbins);
  }
  //setting default mass binning
  int massbins[4] = {200, 100, 80, 80};
  double minmass[4] = {0.4, 1.07, 1.28, 1.63};
  double maxmass[4] = {0.6, 1.17, 1.36, 1.71};
  for (int ipart=0; ipart<knumpart; ipart++) {
    SetMassbinning(ipart, massbins[ipart], minmass[ipart], maxmass[ipart]);
  }
  //setting default pt binning
  double ptbins[4][600];
  int nptbins[4] = {500, 500, 75, 75};
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

  //histograms for K0S variables
  fHistos_K0S = new THistManager("histos_K0S");
  fHistos_Lam = new THistManager("histos_Lam");
  fHistos_ALam = new THistManager("histos_ALam");
  fHistos_K0S->CreateTH3("h3_ptmasscent_def", "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
  fHistos_Lam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
  fHistos_ALam->CreateTH3("h3_ptmasscent_def", "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
  if (!fDefOnly) {
    for (int icut = 0; icut<kV0cutsnum; icut++) {
      if (icut==kV0_y || icut==kV0_etaDaugh || icut==kV0_TOFBunchCrossing) continue;
      for (int ivar = 0; ivar<nvarcut_V0[icut]; ivar++) {
        if (icut!=kV0_PropLifetLam) fHistos_K0S->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kK0s], fptbinning[kK0s], fnmassbins[kK0s], fmassbinning[kK0s], fncentbins[kK0s], fcentbinning[kK0s]);
        if (icut!=kV0_PropLifetK0s) fHistos_Lam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
        if (icut!=kV0_PropLifetK0s) fHistos_ALam->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kLam], fptbinning[kLam], fnmassbins[kLam], fmassbinning[kLam], fncentbins[kLam], fcentbinning[kLam]);
      }
    }
  }
  //histograms for Xi variables
  fHistos_XiMin = new THistManager("histos_XiMin");
  fHistos_XiPlu = new THistManager("histos_XiPlu");
  fHistos_OmMin = new THistManager("histos_OmMin");
  fHistos_OmPlu = new THistManager("histos_OmPlu");
  fHistos_XiMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
  fHistos_XiPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
  fHistos_OmMin->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
  fHistos_OmPlu->CreateTH3("h3_ptmasscent_def", "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
  if (!fDefOnly) {
    for (int icut=0; icut<kCasccutsnum; icut++) {
      if (icut==kCasc_y || icut==kCasc_etaDaugh || icut==kCasc_TOFBunchCrossing) continue;
      for (int ivar = 0; ivar<nvarcut_Casc[icut]; ivar++) {
        if(icut!=kCasc_PropLifetOm) fHistos_XiMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        if(icut!=kCasc_PropLifetOm) fHistos_XiPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kXi], fptbinning[kXi], fnmassbins[kXi], fmassbinning[kXi], fncentbins[kXi], fcentbinning[kXi]);
        if(icut!=kCasc_PropLifetXi) fHistos_OmMin->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
        if(icut!=kCasc_PropLifetXi) fHistos_OmPlu->CreateTH3(Form("h3_ptmasscent[%d][%d]", icut, ivar), "", fnptbins[kOm], fptbinning[kOm], fnmassbins[kOm], fmassbinning[kOm], fncentbins[kOm], fcentbinning[kOm]);
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

  // dumb histo for checking
  fHistos_eve->FillTH1("henum", 0.5);

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
  fHistos_eve->FillTH1("hcent", lPercentile);

  //acquire best PV
  const AliVVertex *lBestPrimVtx = lVevent->GetPrimaryVertex();
  double lBestPV[3] = {-666., -666., -666.};
  lBestPrimVtx->GetXYZ(lBestPV);

  //acquire magnetic field
  double lMagField = -666;
  lMagField = lVevent->GetMagneticField();

  // start v0 part
  int nv0s = 0;
  nv0s = lVevent->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {

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

    }
    //filling with default cuts
    if (ApplyCuts(kk0s)) fHistos_K0S->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassK0s, lPercentile);
    if (ApplyCuts(klam)) fHistos_Lam->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassLam, lPercentile);
    if (ApplyCuts(kalam)) fHistos_ALam->FillTH3("h3_ptmasscent_def", fV0_Pt, fV0_InvMassALam, lPercentile);

    //filling 3D histograms
    if (!fDefOnly) FillHistCutVariations(kFALSE, lPercentile);
  } // end of V0 loop

  // start of cascades part
  //-----------------------
  int ncasc = 0;
  ncasc = lVevent->GetNumberOfCascades();

  for (int i_casc = 0; i_casc < ncasc; i_casc++) {
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
      fCasc_LeastCRawsOvF = (int) (lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac));

      //V0 daughter mass (later to be checked against nominal)
      fCasc_InvMassLam = casc->GetEffMass(); //Note that GetEffMass() is inherited from AliESDv0 and it returns the mass of the V0 (and not of the cascade)
      fCasc_InvMassLam = 1.115683;

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

      //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
      double xn, xp;
      if (casc->Charge()>0) {
        fCasc_DcaBacBar = pTrackCasc->GetDCA(bTrackCasc, lMagField, xn, xp);
      } else {
        fCasc_DcaBacBar = nTrackCasc->GetDCA(bTrackCasc, lMagField, xn, xp);
      }
      fCasc_DcaBacBar = 10.;
    } else {
      //start AOD part
      AliAODcascade *casc = lAODevent->GetCascade(i_casc);
      if (!casc) continue;

      //cascade and V0 2D radii
      fCasc_CascRad = casc->RadiusSecVtx();
      fCasc_V0Rad = casc->RadiusV0();

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
      fCasc_LeastCRawsOvF = (int) (lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac));

      //V0 daughter mass (later to be checked against nominal)
      fCasc_InvMassLam = casc->MassLambda();
      fCasc_InvMassLam = 1.115683;

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

      //transverse momentum
      fCasc_Pt = TMath::Sqrt(casc->Pt2Xi());

      //distance over total momentum
      fCasc_DistOverTotP = casc->DecayLengthXi(lBestPV[0], lBestPV[1], lBestPV[2])/(TMath::Sqrt(casc->Ptot2Xi())+1e-10);

      //candidate's invariant mass
      fCasc_InvMassXiMin = casc->MassXi();
      fCasc_InvMassXiPlu = casc->MassXi();
      fCasc_InvMassOmMin = casc->MassOmega();
      fCasc_InvMassOmPlu = casc->MassOmega();

      //calculate DCA Bachelor-Baryon to remove "bump"structure in InvMass
      fCasc_DcaBacBar = 10; //safe factor. It's not possible to calculate DCA between 2 tracks? Missing cov matrix?
      //double xn, xp;
      //if(casc->Charge()>0) fCasc_DcaBacBar = pTrackCasc->GetDCA(bTrackCasc,lMagField,xn,xp);
      //else fCasc_DcaBacBar = nTrackCasc->GetDCA(bTrackCasc,lMagField,xn,xp);
    }
      //fills TH3 with default cuts
      if(ApplyCuts(kxim)) fHistos_XiMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiMin, lPercentile);
      if(ApplyCuts(kxip)) fHistos_XiPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassXiPlu, lPercentile);
      if(ApplyCuts(komm)) fHistos_OmMin->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmMin, lPercentile);
      if(ApplyCuts(komp)) fHistos_OmPlu->FillTH3("h3_ptmasscent_def", fCasc_Pt, fCasc_InvMassOmPlu, lPercentile);

      //filling 3D histograms
      if (!fDefOnly) {
        FillHistCutVariations(kTRUE, lPercentile);
      }
  }

  DataPosting();

}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetCutVal(bool iscasc, int cutnum, double cval) {
  if (!iscasc) {
    cutval_V0[cutnum] = cval;
  } else {
    cutval_Casc[cutnum] = cval;
  }
}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::SetDefCutVals() {
  //V0 part
  SetCutVal(kFALSE, kV0_DcaV0Daught, 1.);
  SetCutVal(kFALSE, kV0_DcaPosToPV, 0.06);
  SetCutVal(kFALSE, kV0_DcaNegToPV, 0.06);
  SetCutVal(kFALSE, kV0_V0CosPA, 0.97);
  SetCutVal(kFALSE, kV0_V0Rad, 0.5);
  SetCutVal(kFALSE, kV0_y, 0.5);
  SetCutVal(kFALSE, kV0_etaDaugh, 0.8);
  SetCutVal(kFALSE, kV0_LeastCRaws, 70);
  SetCutVal(kFALSE, kV0_LeastCRawsOvF, 0.8);
  SetCutVal(kFALSE, kV0_NSigPID, 5.);
  SetCutVal(kFALSE, kV0_PropLifetK0s, 20);
  SetCutVal(kFALSE, kV0_PropLifetLam, 30);
  SetCutVal(kFALSE, kV0_TOFBunchCrossing, -95);

  //Cascade part
  SetCutVal(kTRUE, kCasc_CascRad, 0.5);
  SetCutVal(kTRUE, kCasc_V0Rad, 1.1);
  SetCutVal(kTRUE, kCasc_DcaBachToPV, 0.04);
  SetCutVal(kTRUE, kCasc_DcaV0ToPV, 0.06);
  SetCutVal(kTRUE, kCasc_DcaMesToPV, 0.04);
  SetCutVal(kTRUE, kCasc_DcaBarToPV, 0.03);
  SetCutVal(kTRUE, kCasc_DcaV0Daught, 1.5);
  SetCutVal(kTRUE, kCasc_DcaCascDaught, 1.3);
  SetCutVal(kTRUE, kCasc_CascCosPA, 0.97);
  SetCutVal(kTRUE, kCasc_V0CosPA, 0.97);
  SetCutVal(kTRUE, kCasc_InvMassLam, 0.008);
  SetCutVal(kTRUE, kCasc_NSigPID, 4);
  SetCutVal(kTRUE, kCasc_PropLifetXi, 3); //in number of ctau in this case, differently from V0s
  SetCutVal(kTRUE, kCasc_PropLifetOm, 3); //in number of ctau in this case, differently from V0s
  SetCutVal(kTRUE, kCasc_LeastCRaws, 80);
  SetCutVal(kTRUE, kCasc_LeastCRawsOvF, 0.8);
  SetCutVal(kTRUE, kCasc_TOFBunchCrossing, -95);
  SetCutVal(kTRUE, kCasc_DcaBacBar, 0.02);
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
  SetCutVariation(kFALSE, kV0_DcaV0Daught, 10, 0.5, 1.5);
  SetCutVariation(kFALSE, kV0_DcaPosToPV, 10 ,0.05, 0.08);
  SetCutVariation(kFALSE, kV0_DcaNegToPV, 10, 0.05, 0.08);
  SetCutVariation(kFALSE, kV0_V0CosPA, 10, 0.95, 0.999);
  SetCutVariation(kFALSE, kV0_V0Rad, 10, 0.3, 0.7);
  SetCutVariation(kFALSE, kV0_LeastCRaws, 5, 70, 80);
  SetCutVariation(kFALSE, kV0_LeastCRawsOvF, 10, 0.75, 0.90);
  SetCutVariation(kFALSE, kV0_NSigPID, 5, 2, 7);
  SetCutVariation(kFALSE, kV0_PropLifetK0s, 10, 10, 40);
  SetCutVariation(kFALSE, kV0_PropLifetLam, 10, 10, 40);

  SetCutVariation(kTRUE, kCasc_CascRad, 10, 0.4, 1.0);
  SetCutVariation(kTRUE, kCasc_V0Rad, 10, 1., 5.);
  SetCutVariation(kTRUE, kCasc_DcaBachToPV, 10, 0.03, 0.17);
  SetCutVariation(kTRUE, kCasc_DcaV0ToPV, 10, 0.05, 0.15);
  SetCutVariation(kTRUE, kCasc_DcaMesToPV, 10, 0.02, 0.3);
  SetCutVariation(kTRUE, kCasc_DcaBarToPV, 10, 0.02, 0.12);
  SetCutVariation(kTRUE, kCasc_DcaV0Daught, 10, 1, 2);
  SetCutVariation(kTRUE, kCasc_DcaCascDaught, 10, 0.8, 2.);
  SetCutVariation(kTRUE, kCasc_CascCosPA, 20, 0.95, 0.99);
  SetCutVariation(kTRUE, kCasc_V0CosPA, 20, 0.95, 0.99);
  SetCutVariation(kTRUE, kCasc_InvMassLam, 4, 0.006, 0.01);
  SetCutVariation(kTRUE, kCasc_NSigPID, 5, 2, 7);
  SetCutVariation(kTRUE, kCasc_PropLifetXi, 10, 2, 5);
  SetCutVariation(kTRUE, kCasc_PropLifetOm, 10, 2, 5);
  SetCutVariation(kTRUE, kCasc_LeastCRaws, 5, 70, 80);
  SetCutVariation(kTRUE, kCasc_LeastCRawsOvF, 10, 0.75, 0.9);
  SetCutVariation(kTRUE, kCasc_DcaBacBar, 10, 0., 0.05);
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

    if((part==kxip || part==komm) && fCasc_charge>0) return kFALSE;
    if((part==kxim || part==komp) && fCasc_charge<0) return kFALSE;

    // check candidate's rapidity (particle hypothesis' dependent)
    if((part==kxip || part==kxim) && TMath::Abs(fCasc_yXi)>0.5) return kFALSE;
    if((part==komp || part==komm) && TMath::Abs(fCasc_yOm)>0.5) return kFALSE;
    // check candidate daughters' pseudo-rapidity
    if(TMath::Abs(fCasc_etaPos)>0.8 || TMath::Abs(fCasc_etaNeg)>0.8 || TMath::Abs(fCasc_etaBac)>0.8) return kFALSE;
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
    if((part==komm || part==komp) && ((1.67245*fCasc_DistOverTotP)>(2.461*cutval_Casc[kCasc_PropLifetXi]))) return kFALSE;   //2.461 is the ctau of om in cm
    // check DCA bachelor-baryon. If it is too small --> bump structure in Inv Mass
    if(fCasc_DcaBacBar<cutval_Casc[kCasc_DcaBacBar]) return kFALSE;
  }
  return kTRUE; //survived!
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
  PostData(2, fHistos_K0S->GetListOfHistograms());
  PostData(3, fHistos_Lam->GetListOfHistograms());
  PostData(4, fHistos_ALam->GetListOfHistograms());
  PostData(5, fHistos_XiMin->GetListOfHistograms());
  PostData(6, fHistos_XiPlu->GetListOfHistograms());
  PostData(7, fHistos_OmMin->GetListOfHistograms());
  PostData(8, fHistos_OmPlu->GetListOfHistograms());

}

//________________________________________________________________________
void AliAnalysisTaskStrVsMult::FillHistCutVariations(bool iscasc, double perc) {

  if (!iscasc) {
    for (int i_cut=0; i_cut<kV0cutsnum; i_cut++) {
      if(i_cut==kV0_y || i_cut==kV0_etaDaugh || i_cut==kV0_TOFBunchCrossing) continue;
      SetDefCutVals(); //reset defaults
      for(int i_var=0; i_var<nvarcut_V0[i_cut]; i_var++){
        //K0S filling
        if (i_cut!=kV0_PropLifetLam) {
          SetCutVal(kFALSE, i_cut, varlowcut_V0[i_cut]+i_var*(varhighcut_V0[i_cut]-varlowcut_V0[i_cut])/nvarcut_V0[i_cut]);
          if (ApplyCuts(kk0s)) fHistos_K0S->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassK0s, perc);
        }
        //Lam and AntiLam filling
        if (i_cut!=kV0_PropLifetK0s) {
          SetCutVal(kFALSE, i_cut, varlowcut_V0[i_cut]+i_var*(varhighcut_V0[i_cut]-varlowcut_V0[i_cut])/nvarcut_V0[i_cut]);
          if(ApplyCuts(klam)) fHistos_Lam->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassLam, perc);
          if(ApplyCuts(kalam)) fHistos_ALam->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fV0_Pt, fV0_InvMassALam, perc);
        }
      }
    }
  } else{
    for(int i_cut=0; i_cut<kCasccutsnum; i_cut++) {
      if(i_cut==kCasc_y || i_cut==kCasc_etaDaugh || i_cut==kCasc_TOFBunchCrossing) continue;
      SetDefCutVals(); //reset defaults
      for (int i_var=0; i_var<nvarcut_Casc[i_cut]; i_var++) {
        //Xi filling
        if (i_cut!=kCasc_PropLifetOm) {
          SetCutVal(kTRUE, i_cut, varlowcut_Casc[i_cut]+i_var*(varhighcut_Casc[i_cut]-varlowcut_Casc[i_cut])/nvarcut_Casc[i_cut]);
          if(ApplyCuts(kxim)) fHistos_XiMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiMin, perc);
          if(ApplyCuts(kxip)) fHistos_XiPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassXiPlu, perc);
        }
        //Lam and AntiLam filling
        if (i_cut!=kCasc_PropLifetXi) {
          SetCutVal(kTRUE, i_cut, varlowcut_Casc[i_cut]+i_var*(varhighcut_Casc[i_cut]-varlowcut_Casc[i_cut])/nvarcut_Casc[i_cut]);
          if(ApplyCuts(komm)) fHistos_OmMin->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmMin, perc);
          if(ApplyCuts(komp)) fHistos_OmPlu->FillTH3(Form("h3_ptmasscent[%d][%d]", i_cut, i_var), fCasc_Pt, fCasc_InvMassOmPlu, perc);
        }
      }
    }
  }
}
