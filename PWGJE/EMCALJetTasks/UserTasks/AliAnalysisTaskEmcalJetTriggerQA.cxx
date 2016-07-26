// Jet trigger QA analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVVZERO.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAODHeader.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEmcalJetTriggerQA.h"

ClassImp(AliAnalysisTaskEmcalJetTriggerQA)

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::AliAnalysisTaskEmcalJetTriggerQA() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetTriggerQA", kTRUE),
  fDebug(kFALSE),
  fTriggerClass(""),
  fContainerFull(0),
  fContainerCharged(1),
  fMaxPatchEnergy(0),
  fMaxPatchADCEnergy(0),
  fTriggerType(-1),
  fNFastOR(16),
  fMainTrigCat(kTriggerLevel1Jet),
  fMainTrigSimple(kFALSE),
  fhNEvents(0),
  fhTriggerbit(0), 
  fHistRhovsCentFull(0),
  fHistRhovsCentCharged(0),
  fh3PtEtaPhiTracks(0),
  fh3PtEtaPhiTracksOnEmcal(0),
  fh3PtEtaPhiTracksToProp(0),
  fh3PtEtaPhiTracksProp(0),
  fh3PtEtaPhiTracksNoProp(0),
  fh2CentPtJetFull(0),
  fh2CentPtJetCharged(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fh2NJetsPtFull(0),
  fh2NJetsPtCharged(0),
  fh3PtEtaAreaJetFull(0),
  fh3PtEtaAreaJetCharged(0),
  fh2PtNConstituentsCharged(0),
  fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0),
  fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0),
  fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtzCharged(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0),
  fh3PtLeadJet2VsPatchEnergy(0),
  fh3PtLeadJet1PatchEnergyVZEROAmp(0),
  fh3PtLeadJet1RawPatchEnergyVZEROAmp(0),
  fh3PatchEnergyEtaPhiCenterJ1(0),
  fh3PatchEnergyEtaPhiCenterJ2(0),
  fh3PatchEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1(0),
  fh3PatchADCEnergyEtaPhiCenterJ2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterAll(0),
  fh3EEtaPhiCell(0),
  fh2ECellVsCent(0),
  fh2CellEnergyVsTime(0),
  fh3EClusELeadingCellVsTime(0),
  fh3JetReacCent(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::AliAnalysisTaskEmcalJetTriggerQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(kFALSE),
  fTriggerClass(""),
  fContainerFull(0),
  fContainerCharged(1),
  fMaxPatchEnergy(0),
  fMaxPatchADCEnergy(0),
  fTriggerType(-1),
  fNFastOR(16),
  fMainTrigCat(kTriggerLevel1Jet),
  fMainTrigSimple(kFALSE),
  fhNEvents(0),
  fhTriggerbit(0),
  fHistRhovsCentFull(0),
  fHistRhovsCentCharged(0),
  fh3PtEtaPhiTracks(0),
  fh3PtEtaPhiTracksOnEmcal(0),
  fh3PtEtaPhiTracksToProp(0),
  fh3PtEtaPhiTracksProp(0),
  fh3PtEtaPhiTracksNoProp(0),
  fh2CentPtJetFull(0),
  fh2CentPtJetCharged(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fh2NJetsPtFull(0),
  fh2NJetsPtCharged(0),
  fh3PtEtaAreaJetFull(0),
  fh3PtEtaAreaJetCharged(0),
  fh2PtNConstituentsCharged(0),
  fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0),
  fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0),
  fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtzCharged(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0),
  fh3PtLeadJet2VsPatchEnergy(0),
  fh3PtLeadJet1PatchEnergyVZEROAmp(0),
  fh3PtLeadJet1RawPatchEnergyVZEROAmp(0),
  fh3PatchEnergyEtaPhiCenterJ1(0),
  fh3PatchEnergyEtaPhiCenterJ2(0),
  fh3PatchEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1(0),
  fh3PatchADCEnergyEtaPhiCenterJ2(0),
  fh3PatchADCEnergyEtaPhiCenterJ1J2(0),
  fh3PatchADCEnergyEtaPhiCenterAll(0),
  fh3EEtaPhiCell(0),
  fh2ECellVsCent(0),
  fh2CellEnergyVsTime(0),
  fh3EClusELeadingCellVsTime(0),
  fh3JetReacCent(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::~AliAnalysisTaskEmcalJetTriggerQA()
{
  // Destructor.
 
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::SelectEvent() {
  //
  // Decide if event should be selected for analysis
  //

  fhNEvents->Fill(3.5);

  if(!fTriggerClass.IsNull()) {
    //Check if requested trigger was fired
    TString trigType1 = "J1";
    TString trigType2 = "J2";
    if(fTriggerClass.Contains("G")) {
      trigType1 = "G1";
      trigType2 = "G2";
    }

    TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();
    if(fTriggerClass.Contains(trigType1.Data()) && fTriggerClass.Contains(trigType2.Data())) { //if events with J1&&J2 are requested
      if(!firedTrigClass.Contains(trigType1.Data()) || !firedTrigClass.Contains(trigType2.Data()) ) //check if both are fired
        return kFALSE;
    }
    else {
      if(!firedTrigClass.Contains(fTriggerClass))
	return kFALSE;
      else if(fTriggerClass.Contains(trigType1.Data()) && firedTrigClass.Contains(trigType2.Data())) //if J2 is requested also add triggers which have J1&&J2. Reject if J1 is requested and J2 is fired
	return kFALSE;
    }
  }
  fhNEvents->Fill(1.5);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::FillTriggerPatchHistos() {

  //Fill trigger patch histos for main trigger

  AliEMCALTriggerPatchInfo *patch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);
  fMaxPatchEnergy = 0;
  fMaxPatchADCEnergy = 0;
  if(patch) {
    fMaxPatchEnergy = patch->GetPatchE();
    fMaxPatchADCEnergy = patch->GetADCAmpGeVRough();
    fh3PatchADCEnergyEtaPhiCenterAll->Fill(fMaxPatchADCEnergy,patch->GetEtaGeo(),patch->GetPhiGeo());
    if(patch->IsJetLow() && !patch->IsJetHigh()) { //main patch only fired low threshold trigger
      fh3PatchEnergyEtaPhiCenterJ2->Fill(patch->GetPatchE(),patch->GetEtaGeo(),patch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ2->Fill(fMaxPatchADCEnergy,patch->GetEtaGeo(),patch->GetPhiGeo());
    }
    else if(patch->IsJetHigh() && !patch->IsJetLow()) { //main patch only fired high threshold trigger - should never happen
      fh3PatchEnergyEtaPhiCenterJ1->Fill(patch->GetPatchE(),patch->GetEtaGeo(),patch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ1->Fill(fMaxPatchADCEnergy,patch->GetEtaGeo(),patch->GetPhiGeo());
    }
    else if(patch->IsJetHigh() && patch->IsJetLow()) { //main patch fired both triggers
      fh3PatchEnergyEtaPhiCenterJ1J2->Fill(patch->GetPatchE(),patch->GetEtaGeo(),patch->GetPhiGeo());
      fh3PatchADCEnergyEtaPhiCenterJ1J2->Fill(fMaxPatchADCEnergy,patch->GetEtaGeo(),patch->GetPhiGeo());
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhNEvents = new TH1F("fhNEvents","fhNEvents;selection;N_{evt}",5,0,5);
  fOutput->Add(fhNEvents);

  fhTriggerbit = new TProfile("fhTriggerbit","fhTriggerbit;;TriggerBit",1,0,1);
  fOutput->Add(fhTriggerbit);

  fHistRhovsCentFull = new TH2F("fHistRhovsCentFull", "fHistRhovsCentFull", 101, -1,  100, 300, 0., 300.);
  fHistRhovsCentFull->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCentFull->GetYaxis()->SetTitle("s#rho_{ch} (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCentFull);

  fHistRhovsCentCharged = new TH2F("fHistRhovsCentCharged", "fHistRhovsCentCharged", 101, -1,  100, 300, 0., 300.);
  fHistRhovsCentCharged->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCentCharged->GetYaxis()->SetTitle("#rho_{ch} (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCentCharged);
    
  Int_t fgkNCentBins = 21;
  Float_t kMinCent   = 0.;
  Float_t kMaxCent   = 105.;
  Double_t *binsCent = new Double_t[fgkNCentBins+1];
  for(Int_t i=0; i<=fgkNCentBins; i++) binsCent[i]=(Double_t)kMinCent + (kMaxCent-kMinCent)/fgkNCentBins*(Double_t)i ;
  binsCent[fgkNCentBins-1] = 100.5;
  binsCent[fgkNCentBins] = 101.5;
    
  Int_t fgkNdEPBins = 18*8;
  Float_t kMindEP   = 0.;
  Float_t kMaxdEP   = 1.*TMath::Pi()/2.;
  Double_t *binsdEP = new Double_t[fgkNdEPBins+1];
  for(Int_t i=0; i<=fgkNdEPBins; i++) binsdEP[i]=(Double_t)kMindEP + (kMaxdEP-kMindEP)/fgkNdEPBins*(Double_t)i ;

  Int_t fgkNPtBins = 200;
  Float_t kMinPt   = -50.;
  Float_t kMaxPt   = 150.;
  Double_t *binsPt = new Double_t[fgkNPtBins+1];
  for(Int_t i=0; i<=fgkNPtBins; i++) binsPt[i]=(Double_t)kMinPt + (kMaxPt-kMinPt)/fgkNPtBins*(Double_t)i ;

  Int_t fgkNPhiBins = 18*8;
  Float_t kMinPhi   = 0.;
  Float_t kMaxPhi   = 2.*TMath::Pi();
  Double_t *binsPhi = new Double_t[fgkNPhiBins+1];
  for(Int_t i=0; i<=fgkNPhiBins; i++) binsPhi[i]=(Double_t)kMinPhi + (kMaxPhi-kMinPhi)/fgkNPhiBins*(Double_t)i ;

  Int_t fgkNEtaBins = 100;
  Float_t fgkEtaMin = -1.;
  Float_t fgkEtaMax =  1.;
  Double_t *binsEta=new Double_t[fgkNEtaBins+1];
  for(Int_t i=0; i<=fgkNEtaBins; i++) binsEta[i]=(Double_t)fgkEtaMin + (fgkEtaMax-fgkEtaMin)/fgkNEtaBins*(Double_t)i ;

  Int_t fgkNAreaBins = 100;
  Float_t kMinArea   = 0.;
  Float_t kMaxArea   = 1.;
  Double_t *binsArea = new Double_t[fgkNAreaBins+1];
  for(Int_t i=0; i<=fgkNAreaBins; i++) binsArea[i]=(Double_t)kMinArea + (kMaxArea-kMinArea)/fgkNAreaBins*(Double_t)i ;

  Int_t fgkNConstBins = 100;
  Float_t kMinConst   = 0.;
  Float_t kMaxConst   = 100.;
  Double_t *binsConst = new Double_t[fgkNConstBins+1];
  for(Int_t i=0; i<=fgkNConstBins; i++) binsConst[i]=(Double_t)kMinConst + (kMaxConst-kMinConst)/fgkNConstBins*(Double_t)i ;

  Int_t fgkNMeanPtBins = 100;
  Float_t kMinMeanPt   = 0.;
  Float_t kMaxMeanPt   = 20.;
  Double_t *binsMeanPt = new Double_t[fgkNMeanPtBins+1];
  for(Int_t i=0; i<=fgkNMeanPtBins; i++) binsMeanPt[i]=(Double_t)kMinMeanPt + (kMaxMeanPt-kMinMeanPt)/fgkNMeanPtBins*(Double_t)i ;

  Int_t fgkNNEFBins = 101;
  Float_t kMinNEF   = 0.;
  Float_t kMaxNEF   = 1.01;
  Double_t *binsNEF = new Double_t[fgkNNEFBins+1];
  for(Int_t i=0; i<=fgkNNEFBins; i++) binsNEF[i]=(Double_t)kMinNEF + (kMaxNEF-kMinNEF)/fgkNNEFBins*(Double_t)i ;

  Int_t fgkNzBins = 101;
  Float_t kMinz   = 0.;
  Float_t kMaxz   = 1.01;
  Double_t *binsz = new Double_t[fgkNzBins+1];
  for(Int_t i=0; i<=fgkNzBins; i++) binsz[i]=(Double_t)kMinz + (kMaxz-kMinz)/fgkNzBins*(Double_t)i ;

  Int_t fgkNJetTypeBins = 2;
  Float_t kMinJetType   = -0.5;
  Float_t kMaxJetType   = 1.5;
  Double_t *binsJetType = new Double_t[fgkNJetTypeBins+1];
  for(Int_t i=0; i<=fgkNJetTypeBins; i++) binsJetType[i]=(Double_t)kMinJetType + (kMaxJetType-kMinJetType)/fgkNJetTypeBins*(Double_t)i ;

  Int_t fgkNTimeBins = 100;
  Float_t kMinTime   = -200.;
  Float_t kMaxTime   = 200;
  Double_t *binsTime = new Double_t[fgkNTimeBins+1];
  for(Int_t i=0; i<=fgkNTimeBins; i++) binsTime[i]=(Double_t)kMinTime + (kMaxTime-kMinTime)/fgkNTimeBins*(Double_t)i ;

  Int_t fgkNVZEROBins = 100;
  Float_t kMinVZERO   = 0.;
  Float_t kMaxVZERO   = 25000;
  Double_t *binsVZERO = new Double_t[fgkNVZEROBins+1];
  for(Int_t i=0; i<=fgkNVZEROBins; i++) binsVZERO[i]=(Double_t)kMinVZERO + (kMaxVZERO-kMinVZERO)/fgkNVZEROBins*(Double_t)i ;

  Double_t enBinEdges[3][2];
  enBinEdges[0][0] = 1.; //10 bins
  enBinEdges[0][1] = 0.1;
  enBinEdges[1][0] = 5.; //8 bins
  enBinEdges[1][1] = 0.5;
  enBinEdges[2][0] = 100.;//95 bins
  enBinEdges[2][1] = 1.;

  const Float_t enmin1 =  0;
  const Float_t enmax1 =  enBinEdges[0][0];
  const Float_t enmin2 =  enmax1 ;
  const Float_t enmax2 =  enBinEdges[1][0];
  const Float_t enmin3 =  enmax2 ;
  const Float_t enmax3 =  enBinEdges[2][0];//fgkEnMax;
  const Int_t nbin11 = (int)((enmax1-enmin1)/enBinEdges[0][1]);
  const Int_t nbin12 = (int)((enmax2-enmin2)/enBinEdges[1][1])+nbin11;
  const Int_t nbin13 = (int)((enmax3-enmin3)/enBinEdges[2][1])+nbin12;

  Int_t fgkNEnBins=nbin13;
  Double_t *binsEn=new Double_t[fgkNEnBins+1];
  for(Int_t i=0; i<=fgkNEnBins; i++) {
    if(i<=nbin11) binsEn[i]=(Double_t)enmin1 + (enmax1-enmin1)/nbin11*(Double_t)i ;
    if(i<=nbin12 && i>nbin11) binsEn[i]=(Double_t)enmin2 + (enmax2-enmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;
    if(i<=nbin13 && i>nbin12) binsEn[i]=(Double_t)enmin3 + (enmax3-enmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;
  }

  fh3PtEtaPhiTracks = new TH3F("fh3PtEtaPhiTracks","fh3PtEtaPhiTracks;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracks);

  fh3PtEtaPhiTracksOnEmcal = new TH3F("fh3PtEtaPhiTracksOnEmcal","fh3PtEtaPhiTracksOnEmcal;#it{p}_{T}^{track}_{emc};#eta_{emc};#varphi_{emc}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksOnEmcal);

  fh3PtEtaPhiTracksToProp = new TH3F("fh3PtEtaPhiTracksToProp","fh3PtEtaPhiTracksToProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksToProp);

  fh3PtEtaPhiTracksProp = new TH3F("fh3PtEtaPhiTracksProp","fh3PtEtaPhiTracksProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksProp);

  fh3PtEtaPhiTracksNoProp = new TH3F("fh3PtEtaPhiTracksNoProp","fh3PtEtaPhiTracksNoProp;#it{p}_{T}^{track}_{vtx};#eta_{vtx};#varphi_{vtx}",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiTracksNoProp);

  fh2CentPtJetFull = new TH2F("fh2CentPtJetFull","fh2CentPtJetFull;cent;#it{p}_{T}^{jet}",fgkNCentBins,binsCent,fgkNPtBins,binsPt);
  fOutput->Add(fh2CentPtJetFull);

  fh2CentPtJetCharged = new TH2F("fh2CentPtJetCharged","fh2CentPtJetCharged;cent;#it{p}_{T}^{jet}",fgkNCentBins,binsCent,fgkNPtBins,binsPt);
  fOutput->Add(fh2CentPtJetCharged);

  fh3PtEtaPhiJetFull = new TH3F("fh3PtEtaPhiJetFull","fh3PtEtaPhiJetFull;#it{p}_{T}^{jet};#eta;#varphi",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiJetFull);

  fh3PtEtaPhiJetCharged = new TH3F("fh3PtEtaPhiJetCharged","fh3PtEtaPhiJetCharged;#it{p}_{T}^{jet};#eta;#varphi",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PtEtaPhiJetCharged);

  fh2NJetsPtFull = new TH2F("fh2NJetsPtFull","fh2NJetsPtFull;N_{jets};#it{p}_{T}^{jet}",20,-0.5,19.5,fgkNPtBins,binsPt);
  fOutput->Add(fh2NJetsPtFull);

  fh2NJetsPtCharged = new TH2F("fh2NJetsPtCharged","fh2NJetsPtCharged;N_{jets};#it{p}_{T}^{jet}",20,-0.5,19.5,fgkNPtBins,binsPt);
  fOutput->Add(fh2NJetsPtCharged);

  fh3PtEtaAreaJetFull = new TH3F("fh3PtEtaAreaJetFull","fh3PtEtaAreaJetFull;#it{p}_{T}^{jet};#eta;A",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNAreaBins,binsArea);
  fOutput->Add(fh3PtEtaAreaJetFull);

  fh3PtEtaAreaJetCharged = new TH3F("fh3PtEtaAreaJetCharged","fh3PtEtaAreaJetCharged;#it{p}_{T}^{jet};#eta;A",fgkNPtBins,binsPt,fgkNEtaBins,binsEta,fgkNAreaBins,binsArea);
  fOutput->Add(fh3PtEtaAreaJetCharged);

  fh2PtNConstituentsCharged = new TH2F("fh2PtNConstituentsCharged","fh2PtNConstituentsCharged;#it{p}_{T}^{jet};N_{charged constituents}",fgkNPtBins,binsPt,fgkNConstBins,binsConst);
  fOutput->Add(fh2PtNConstituentsCharged);

  fh2PtNConstituents = new TH2F("fh2PtNConstituents","fh2PtNConstituents;#it{p}_{T}^{jet};N_{constituents}",fgkNPtBins,binsPt,fgkNConstBins,binsConst);
  fOutput->Add(fh2PtNConstituents);

  fh2PtMeanPtConstituentsCharged = new TH2F("fh2PtMeanPtConstituentsCharged","fh2PtMeanPtConstituentsCharged;#it{p}_{T}^{jet};charged #langle #it{p}_{T} #rangle",fgkNPtBins,binsPt,fgkNMeanPtBins,binsMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsCharged);

  fh2PtMeanPtConstituentsNeutral = new TH2F("fh2PtMeanPtConstituentsNeutral","fh2PtMeanPtConstituentsNeutral;#it{p}_{T}^{jet};neutral langle #it{p}_{T} #rangle",fgkNPtBins,binsPt,fgkNMeanPtBins,binsMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsNeutral);

  fh2PtNEF = new TH2F("fh2PtNEF","fh2PtNEF;#it{p}_{T}^{jet};NEF",fgkNPtBins,binsPt,fgkNNEFBins,binsNEF);
  fOutput->Add(fh2PtNEF);

  fh3NEFEtaPhi = new TH3F("fh3NEFEtaPhi","fh3NEFEtaPhi;NEF;#eta;#varphi",fgkNNEFBins,binsNEF,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3NEFEtaPhi);

  fh2NEFNConstituentsCharged = new TH2F("fh2NEFNConstituentsCharged","fh2NEFNConstituentsCharged;NEF;N_{charged constituents}",fgkNNEFBins,binsNEF,fgkNConstBins,binsConst);
  fOutput->Add(fh2NEFNConstituentsCharged);

  fh2NEFNConstituentsNeutral = new TH2F("fh2NEFNConstituentsNeutral","fh2NEFNConstituentsNeutral;NEF;N_{clusters}",fgkNNEFBins,binsNEF,fgkNConstBins,binsConst);
  fOutput->Add(fh2NEFNConstituentsNeutral);

  fh2Ptz = new TH2F("fh2Ptz","fh2Ptz;#it{p}_{T}^{jet};z=p_{t,trk}^{proj}/p_{jet}",fgkNPtBins,binsPt,fgkNzBins,binsz);
  fOutput->Add(fh2Ptz);

  fh2PtzCharged = new TH2F("fh2PtzCharged","fh2Ptz;#it{p}_{T}^{ch jet};z=p_{t,trk}^{proj}/p_{ch jet}",fgkNPtBins,binsPt,fgkNzBins,binsz);
  fOutput->Add(fh2PtzCharged);

  fh2PtLeadJet1VsLeadJet2 = new TH2F("fh2PtLeadJet1VsLeadJet2","fh2PtLeadJet1VsLeadJet2;#it{p}_{T}^{jet 1};#it{p}_{T}^{jet 2}",fgkNPtBins,binsPt,fgkNPtBins,binsPt);
  fOutput->Add(fh2PtLeadJet1VsLeadJet2);

  fh3EEtaPhiCluster = new TH3F("fh3EEtaPhiCluster","fh3EEtaPhiCluster;E_{clus};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3EEtaPhiCluster);

  fh3PtLeadJet1VsPatchEnergy = new TH3F("fh3PtLeadJet1VsPatchEnergy","fh3PtLeadJet1VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNJetTypeBins,binsJetType);
  fOutput->Add(fh3PtLeadJet1VsPatchEnergy);
  fh3PtLeadJet2VsPatchEnergy = new TH3F("fh3PtLeadJet2VsPatchEnergy","fh3PtLeadJet2VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNJetTypeBins,binsJetType);
  fOutput->Add(fh3PtLeadJet2VsPatchEnergy);

  fh3PtLeadJet1PatchEnergyVZEROAmp = new TH3F("fh3PtLeadJet1PatchEnergyVZEROAmp","fh3PtLeadJet1VsPatchEnergyVZEROAmp;#it{p}_{T}^{jet 1};Amplitude_{patch};VZERO amp",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNVZEROBins,binsVZERO);
  fOutput->Add(fh3PtLeadJet1PatchEnergyVZEROAmp);
  fh3PtLeadJet1RawPatchEnergyVZEROAmp = new TH3F("fh3PtLeadJet1RawPatchEnergyVZEROAmp","fh3PtLeadJet1RawPatchEnergyVZEROAmp;#it{p}_{T}^{jet 1};ADC Amplitude_{patch} (GeV);VZERO amp",fgkNPtBins,binsPt,fgkNPtBins,binsPt,fgkNVZEROBins,binsVZERO);
  fOutput->Add(fh3PtLeadJet1RawPatchEnergyVZEROAmp);

  fh3PatchEnergyEtaPhiCenterJ1 = new TH3F("fh3PatchEnergyEtaPhiCenterJ1","fh3PatchEnergyEtaPhiCenterJ1;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ1);

  fh3PatchEnergyEtaPhiCenterJ2 = new TH3F("fh3PatchEnergyEtaPhiCenterJ2","fh3PatchEnergyEtaPhiCenterJ2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ2);

  fh3PatchEnergyEtaPhiCenterJ1J2 = new TH3F("fh3PatchEnergyEtaPhiCenterJ1J2","fh3PatchEnergyEtaPhiCenterJ1J2;E_{patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenterJ1J2);

  fh3PatchADCEnergyEtaPhiCenterJ1 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ1","fh3PatchADCEnergyEtaPhiCenterJ1;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ1);

  fh3PatchADCEnergyEtaPhiCenterJ2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ2","fh3PatchADCEnergyEtaPhiCenterJ2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ2);

  fh3PatchADCEnergyEtaPhiCenterJ1J2 = new TH3F("fh3PatchADCEnergyEtaPhiCenterJ1J2","fh3PatchADCEnergyEtaPhiCenterJ1J2;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterJ1J2);

  fh3PatchADCEnergyEtaPhiCenterAll = new TH3F("fh3PatchADCEnergyEtaPhiCenterAll","fh3PatchADCEnergyEtaPhiCenterAll;E_{ADC,patch};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3PatchADCEnergyEtaPhiCenterAll);

  fh3EEtaPhiCell = new TH3F("fh3EEtaPhiCell","fh3EEtaPhiCell;E_{cell};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3EEtaPhiCell);

  fh2ECellVsCent = new TH2F("fh2ECellVsCent","fh2ECellVsCent;centrality;E_{cell}",101,-1,100,500,0.,5.);
  fOutput->Add(fh2ECellVsCent);

  fh2CellEnergyVsTime = new TH2F("fh2CellEnergyVsTime","fh2CellEnergyVsTime;E_{cell};time",fgkNEnBins,binsEn,fgkNTimeBins,binsTime);
  fOutput->Add(fh2CellEnergyVsTime);

  fh3EClusELeadingCellVsTime = new TH3F("fh3EClusELeadingCellVsTime","fh3EClusELeadingCellVsTime;E_{cluster};E_{leading cell};time_{leading cell}",fgkNEnBins,binsEn,fgkNEnBins,binsEn,fgkNTimeBins,binsTime);
  fOutput->Add(fh3EClusELeadingCellVsTime);
    
  fh3JetReacCent = new TH3F("fh3JetReacCent","fh3JetReacCent;E_{Jet};Centrality;dEP",fgkNEnBins,binsEn,fgkNCentBins,binsCent,fgkNdEPBins,binsdEP);
  fOutput->Add(fh3JetReacCent);
    
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
    if (h2){
      h2->Sumw2();
      continue;
    }
    TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
    if (h3){
      h3->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  delete [] binsCent;
  if(binsdEP)               delete [] binsdEP;
  if(binsEn)                delete [] binsEn;
  if(binsPt)                delete [] binsPt;
  if(binsPhi)               delete [] binsPhi;
  if(binsEta)               delete [] binsEta;
  if(binsArea)              delete [] binsArea;
  if(binsConst)             delete [] binsConst; 
  if(binsMeanPt)            delete [] binsMeanPt;  
  if(binsNEF)               delete [] binsNEF;
  if(binsz)                 delete [] binsz;
  if(binsJetType)           delete [] binsJetType;
  if(binsTime)              delete [] binsTime;
  if(binsVZERO)             delete [] binsVZERO;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::FillHistograms()
{
  // Fill histograms.

  //Tracks
  AliParticleContainer *partCont = GetParticleContainer(0);
  if (partCont) {
    partCont->ResetCurrentID();
    AliPicoTrack *track = dynamic_cast<AliPicoTrack*>(partCont->GetNextAcceptParticle());
    while(track) {
      Double_t trkphi = track->Phi()*TMath::RadToDeg();
      fh3PtEtaPhiTracks->Fill(track->Pt(),track->Eta(),track->Phi());
      //Select tracks which should be propagated
      if(track->Pt()>=0.350) {
	if (TMath::Abs(track->Eta())<=0.9 && trkphi > 10 && trkphi < 250) {
	  fh3PtEtaPhiTracksOnEmcal->Fill(track->GetTrackPtOnEMCal(),track->GetTrackEtaOnEMCal(),track->GetTrackPhiOnEMCal());
	  fh3PtEtaPhiTracksToProp->Fill(track->Pt(),track->Eta(),track->Phi());
	  if(track->GetTrackPtOnEMCal()>=0)
	    fh3PtEtaPhiTracksProp->Fill(track->Pt(),track->Eta(),track->Phi());
	  else
	    fh3PtEtaPhiTracksNoProp->Fill(track->Pt(),track->Eta(),track->Phi());
	}
      }
      track = dynamic_cast<AliPicoTrack*>(partCont->GetNextAcceptParticle());
    }
  }

  //Clusters
  AliClusterContainer  *clusCont = GetClusterContainer(0);
  if (clusCont) {
    Int_t nclusters = clusCont->GetNClusters();
    for (Int_t ic = 0; ic < nclusters; ic++) {
      AliVCluster *cluster = static_cast<AliVCluster*>(clusCont->GetCluster(ic));
      if (!cluster) {
	AliDebug(2,Form("Could not receive cluster %d", ic));
	continue;
      }
      if (!cluster->IsEMCAL()) {
	AliDebug(2,Form("%s: Cluster is not emcal",GetName()));
	continue;
      }

      TLorentzVector lp;
      cluster->GetMomentum(lp, const_cast<Double_t*>(fVertex));
      fh3EEtaPhiCluster->Fill(lp.E(),lp.Eta(),lp.Phi());
      if(fCaloCells) {
	Double_t leadCellE = GetEnergyLeadingCell(cluster);
	Double_t leadCellT = cluster->GetTOF();
	fh3EClusELeadingCellVsTime->Fill(lp.E(),leadCellE,leadCellT*1e9);
      }
    }
  }

  //Cells
  if(fCaloCells) {
    const Short_t nCells   = fCaloCells->GetNumberOfCells();

    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = fCaloCells->GetCellNumber(iCell);
      Double_t cellE = fCaloCells->GetCellAmplitude(cellId);
      Double_t cellT = fCaloCells->GetCellTime(cellId);
      TVector3 pos;
      fGeom->GetGlobal(cellId, pos);
      TLorentzVector lv(pos,cellE);
      Double_t cellEta = lv.Eta();
      Double_t cellPhi = lv.Phi();
      if(cellPhi<0.) cellPhi+=TMath::TwoPi();
      if(cellPhi>TMath::TwoPi()) cellPhi-=TMath::TwoPi();

      AliDebug(2,Form("cell energy = %f  time = %f",cellE,cellT*1e9));
      fh2CellEnergyVsTime->Fill(cellE,cellT*1e9);
      fh3EEtaPhiCell->Fill(cellE,cellEta,cellPhi);
      fh2ECellVsCent->Fill(fCent,cellE);
    }
  }

  //Jets
  Double_t ptLeadJet1 = 0.;
  Double_t ptLeadJet2 = 0.;

  fHistRhovsCentFull->Fill(fCent,GetRhoVal(fContainerFull));
  fHistRhovsCentCharged->Fill(fCent,GetRhoVal(fContainerCharged));

  TArrayI *nJetsArr = new TArrayI(fh2NJetsPtFull->GetNbinsY()+1);
  nJetsArr->Reset(0);
  nJetsArr->Set(fh2NJetsPtFull->GetNbinsY()+1);

  if (GetJetContainer(fContainerFull)) {
    const Int_t njets = GetNJets(fContainerFull);
    for (Int_t ij = 0; ij < njets; ij++) {

      AliEmcalJet* jet = GetAcceptJetFromArray(ij,fContainerFull);
      if (!jet)
	continue; //jet not selected

      Double_t jetPt = jet->Pt() - GetRhoVal(fContainerFull)*jet->Area();
      if(jetPt>ptLeadJet1) ptLeadJet1=jetPt;

      Double_t dEPJetFull = RelativeEP(jet->Phi() , fEPV0);
      fh3JetReacCent->Fill(jet->E(),fCent,dEPJetFull);
      
      fh2CentPtJetFull->Fill(fCent,jetPt);
      fh3PtEtaPhiJetFull->Fill(jetPt,jet->Eta(),jet->Phi());
      fh3PtEtaAreaJetFull->Fill(jetPt,jet->Eta(),jet->Area());

      //count jets above certain pT threshold
      Int_t ptbin = fh2NJetsPtFull->GetYaxis()->FindBin(jetPt);
      for(Int_t iptbin = ptbin; iptbin<=fh2NJetsPtFull->GetNbinsY(); iptbin++)
	nJetsArr->AddAt(nJetsArr->At(iptbin)+1,iptbin);
      
      fh2PtNConstituentsCharged->Fill(jetPt,jet->GetNumberOfTracks());
      fh2PtNConstituents->Fill(jetPt,jet->GetNumberOfConstituents());
      fh2PtNEF->Fill(jetPt,jet->NEF());
      fh3NEFEtaPhi->Fill(jet->NEF(),jet->Eta(),jet->Phi());
      fh2NEFNConstituentsCharged->Fill(jet->NEF(),jet->GetNumberOfTracks());
      fh2NEFNConstituentsNeutral->Fill(jet->NEF(),jet->GetNumberOfClusters());

      AliVParticle *vp;
      Double_t sumPtCh = 0.;
      for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
	vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
	if(!vp) continue;
	fh2Ptz->Fill(jetPt,GetZ(vp,jet));
	sumPtCh+=vp->Pt();
      }
      
      if(jet->GetNumberOfTracks()>0)
	fh2PtMeanPtConstituentsCharged->Fill(jetPt,sumPtCh/(double)(jet->GetNumberOfTracks()) );

      AliVCluster *vc = 0x0;
      Double_t sumPtNe = 0.;
      if (clusCont) {
	for(Int_t icc=0; icc<jet->GetNumberOfClusters(); icc++) {
	  vc = static_cast<AliVCluster*>(clusCont->GetCluster(icc));
	  if(!vc) continue;
	  TLorentzVector lp;
	  vc->GetMomentum(lp, const_cast<Double_t*>(fVertex));
	  sumPtNe+=lp.Pt();
	}

	if(jet->GetNumberOfClusters()>0)
	  fh2PtMeanPtConstituentsNeutral->Fill(jetPt,sumPtNe/(double)(jet->GetNumberOfClusters()) );
      }
    } //full jet loop

    for(Int_t i=1; i<=fh2NJetsPtFull->GetNbinsY(); i++) {
      Int_t nJetsInEvent = nJetsArr->At(i);
      fh2NJetsPtFull->Fill(nJetsInEvent,fh2NJetsPtFull->GetYaxis()->GetBinCenter(i));
    }
  }

  //Reset array to zero to also count charged jets
  nJetsArr->Reset(0);
  
  //Loop over charged jets
  if (GetJetContainer(fContainerCharged)) {
    const Int_t njets = GetNJets(fContainerCharged);
    for (Int_t ij = 0; ij < njets; ij++) {

      AliEmcalJet* jet = GetAcceptJetFromArray(ij,fContainerCharged);
      if (!jet)
	continue; //jet not selected

      Double_t jetPt = jet->Pt() - GetRhoVal(fContainerCharged)*jet->Area();
      if(jetPt>ptLeadJet2) ptLeadJet2=jetPt;
      fh2CentPtJetCharged->Fill(fCent,jetPt);
      fh3PtEtaPhiJetCharged->Fill(jetPt,jet->Eta(),jet->Phi());
      fh3PtEtaAreaJetCharged->Fill(jetPt,jet->Eta(),jet->Area());

      AliVParticle *vp;
      for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
	vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
	if(!vp) continue;
	fh2PtzCharged->Fill(jetPt,GetZ(vp,jet));
      }
      
      //count jets above certain pT threshold
      Int_t ptbin = fh2NJetsPtCharged->GetYaxis()->FindBin(jetPt);
      for(Int_t iptbin = ptbin; iptbin<=fh2NJetsPtCharged->GetNbinsY(); iptbin++)
	nJetsArr->AddAt(nJetsArr->At(iptbin)+1,iptbin);
      
    }//ch jet loop
    for(Int_t i=1; i<=fh2NJetsPtCharged->GetNbinsY(); i++) {
      Int_t nJetsInEvent = nJetsArr->At(i);
      fh2NJetsPtCharged->Fill(nJetsInEvent,fh2NJetsPtCharged->GetYaxis()->GetBinCenter(i));
    }
  }

  if(GetJetContainer(fContainerFull) && GetJetContainer(fContainerCharged))
    fh2PtLeadJet1VsLeadJet2->Fill(ptLeadJet1,ptLeadJet2);

  fh3PtLeadJet1VsPatchEnergy->Fill(ptLeadJet1,fMaxPatchEnergy,fTriggerType);
  fh3PtLeadJet2VsPatchEnergy->Fill(ptLeadJet2,fMaxPatchEnergy,fTriggerType);

  // Get VZERO amplitude
  Float_t VZEROAmp = InputEvent()->GetVZEROData()->GetTriggerChargeA() + InputEvent()->GetVZEROData()->GetTriggerChargeC();

  fh3PtLeadJet1PatchEnergyVZEROAmp->Fill(ptLeadJet1,fMaxPatchEnergy,VZEROAmp);
  fh3PtLeadJet1RawPatchEnergyVZEROAmp->Fill(ptLeadJet1,fMaxPatchADCEnergy,VZEROAmp);

  delete nJetsArr;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  fhTriggerbit->Fill(0.5,GetCollisionCandidates());

  //Check if event is selected (vertex & pile-up)
  if(!SelectEvent())
    return kFALSE;
  
  if(fTriggerPatchInfo) 
    FillTriggerPatchHistos();

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetZ(const AliVParticle *trk, const AliEmcalJet *jet)          const
{  
  // Get Z of constituent trk
  return GetZ(trk->Px(),trk->Py(),trk->Pz(),jet->Px(),jet->Py(),jet->Pz());
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetZ(Double_t trkPx, Double_t trkPy, Double_t trkPz, Double_t jetPx, Double_t jetPy, Double_t jetPz) const
{
  // 
  // Get the z of a constituent inside of a jet
  //
  Double_t pJetSq = jetPx*jetPx+jetPy*jetPy+jetPz*jetPz;
  if(pJetSq>0.)
    return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/pJetSq;
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f",GetName(), pJetSq)); 
    return 0;
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetTriggerQA::GetLeadingCellId(const AliVCluster *clus) const
{
  //Get energy of leading cell in cluster

  if(!fCaloCells)
    return -1;

  Double_t emax = -1.;
  Int_t iCellAbsIdMax = -1;
  Int_t nCells = clus->GetNCells();
  for(Int_t i = 0; i<nCells; i++) {
    Int_t absId = clus->GetCellAbsId(i);
    Double_t cellE = fCaloCells->GetCellAmplitude(absId);
    if(cellE>emax) {
      emax = cellE;
      iCellAbsIdMax = absId;
    }
  }
  return iCellAbsIdMax;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetEnergyLeadingCell(const AliVCluster *clus) const
{
  //Get energy of leading cell in cluster
  if(!fCaloCells)
    return -1.;

  Int_t absID = GetLeadingCellId(clus);
  if(absID>-1)
    return fCaloCells->GetCellAmplitude(absID);
  else 
    return -1.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetECross(Int_t absID) const {

  //Get Ecross = sum of energy of neighbouring cells (using uncalibrated energy)

  if(!fCaloCells)
    return -1.;

  Double_t ecross = -1.;

  Int_t absID1 = -1;
  Int_t absID2 = -1;
  Int_t absID3 = -1;
  Int_t absID4 = -1;

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  fGeom->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

  if( iphi < AliEMCALGeoParams::fgkEMCALRows-1)
    absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if( iphi > 0 )
    absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  if( ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) ) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  }
  else if( ieta == 0 && imod%2 ) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  }
  else  {
    if( ieta < AliEMCALGeoParams::fgkEMCALCols-1 )
      absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if( ieta > 0 )
      absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
  }

  Double_t ecell1 = fCaloCells->GetCellAmplitude(absID1);
  Double_t ecell2 = fCaloCells->GetCellAmplitude(absID2);
  Double_t ecell3 = fCaloCells->GetCellAmplitude(absID3);
  Double_t ecell4 = fCaloCells->GetCellAmplitude(absID4);

  ecross = ecell1+ecell2+ecell3+ecell4;

  return ecross;
}

//_________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetTriggerQA::RelativeEP(Double_t objAng, Double_t EPAng) const
{
  // function to calculate angle between object and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi = EPAng - objAng;

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() )
    dphi = dphi + 1*TMath::Pi();
  if( dphi>1*TMath::Pi())
    dphi = dphi - 1*TMath::Pi();

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  }

  return dphi;   // dphi in [0, Pi/2]
}
