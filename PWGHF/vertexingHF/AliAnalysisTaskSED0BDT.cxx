/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison with the MC truth and cut variables distributions.
//
// Authors: A.Dainese, andrea.dainese@lnl.infn.it
// Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Carmelo Di Giglio, carmelo.digiglio@ba.infn.it (like sign)
// Jeremy Wilkinson, jwilkinson@physi.uni-heidelberg.de (weighted Bayesian
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include "AliESDUtils.h"
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include "AliVertexingHFUtils.h"
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0BDT.h"
#include "AliNormalizationCounter.h"
#include "AliEventCuts.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSED0BDT);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSED0BDT::AliAnalysisTaskSED0BDT():
  AliAnalysisTaskSE(),
  fOutputMass(0),
  fOutputMassPt(0),
  fOutputMassY(0),
  fDistr(0),
  fNentries(0),
  fMCAccPrompt(0),
  fMCAccBFeed(0),
  fStepMCAcc(kTRUE),
  fCuts(0),
  fEventCuts(),
  fEnableCentralityCorrCuts(kFALSE),
  fArray(0),
  fReadMC(0),
  fCutOnDistr(0),
  fUsePid4Distr(0),
  fCounter(0),
  fNPtBins(1),
  fLsNormalization(1.),
  fFillOnlyD0D0bar(0),
  fDaughterTracks(),
  fIsSelectedCandidate(0),
  fFillVarHists(kTRUE),
  fSys(0),
  fIsRejectSDDClusters(0),
  fFillPtHist(kTRUE),
  fFillYHist(kFALSE),
  fFillImpParHist(kFALSE),
  fFillSubSampleHist(kFALSE),
  fEventCounter(0),
  fUseSelectionBit(kTRUE),
  fAODProtection(0),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fCandidateVariables(),
  fWriteProtosgnVar(kFALSE),
  fHistNtrCorrEvSel(0),
  fSelectTrueD0(kFALSE),
  fUsedMassWindow(kFALSE),
  fPIDCheck(kFALSE),
  fDrawDetSignal(kFALSE),
  fUseQuarkTagInKine(kTRUE),
  fFillSparses(0),
  fUseRejectionMethod(kFALSE),
  fRejectionFactor(0.01),
  fhStudyImpParSingleTrackSign(0),
  fhStudyImpParSingleTrackCand(0),
  fhStudyImpParSingleTrackFd(0),
  fDetSignal(0),
  fhMultVZEROTPCoutTrackCorrNoCut(0x0),
  fhMultVZEROTPCoutTrackCorr(0x0),
  fEnablePileupRejVZEROTPCout(kFALSE),
  fhMultVZEROTPCclustersCorrNoCut(0x0),
  fhMultVZEROTPCclustersCorr(0x0),
  fEnablePileupRejVZEROTPCcls(kFALSE),
  fRejectOutOfBunchPileUp(kFALSE),
  fCut4BDTptbin(0),
  fListRDHFBDT(0),
  fListBDTNames(0),
  fListBDTNtuple(0),
  fListBDTResp(0),
  fBDTSidebandSamplingFraction(0.1),
  fSampleSideband(kFALSE),
fBDTFullVarString("ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi:multiplicity"),
  fBDTClassifierVarString(""),
fSubtractTrackletsFromDau(kFALSE),
fRefMult(9.26),
fMultiplicityEstimator(kNtrk10),
fMCPrimariesEstimator(kEta10),
fListProfiles(0),
fmultiana(0),
fCounterC(0),
fDoVZER0ParamVertexCorr(1),
fYearNumber(16)
{
  /// Default constructor
    for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;
    for(Int_t i=0; i<8; i++) h3Invmass[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_19[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_1029[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_3059[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_19999[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_6099[i]=0x0;
  for(Int_t ih=0; ih<5; ih++) fHistMassPtImpParTC[ih]=0x0;
  fBDTPtCut[0]=0; fBDTPtCut[1]=1e9;
  //~ fBDTSidebandSamplingFraction=0.01;
  //~ fBDTFullVarString="ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi";
  //~ fBDTClassifierVarString="";
}

//________________________________________________________________________
AliAnalysisTaskSED0BDT::AliAnalysisTaskSED0BDT(const char *name,AliRDHFCutsD0toKpi* cuts):
  AliAnalysisTaskSE(name),
  fOutputMass(0),
  fOutputMassPt(0),
  fOutputMassY(0),
  fDistr(0),
  fNentries(0),
  fMCAccPrompt(0),
  fMCAccBFeed(0),
  fStepMCAcc(kTRUE),
  fCuts(0),
  fEventCuts(),
  fEnableCentralityCorrCuts(kFALSE),
  fArray(0),
  fReadMC(0),
  fCutOnDistr(0),
  fUsePid4Distr(0),
  fCounter(0),
  fNPtBins(1),
  fLsNormalization(1.),
  fFillOnlyD0D0bar(0),
  fDaughterTracks(),
  fIsSelectedCandidate(0),
  fFillVarHists(kTRUE),
  fSys(0),
  fIsRejectSDDClusters(0),
  fHistNtrCorrEvSel(0),
  fFillPtHist(kTRUE),
  fFillYHist(kFALSE),
  fFillImpParHist(kFALSE),
  fFillSubSampleHist(kFALSE),
  fEventCounter(0),
  fUseSelectionBit(kTRUE),
  fAODProtection(0),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fCandidateVariables(),
  fWriteProtosgnVar(kFALSE),
  fSelectTrueD0(kFALSE),
  fUsedMassWindow(kFALSE),
  fPIDCheck(kFALSE),
  fDrawDetSignal(kFALSE),
  fUseQuarkTagInKine(kTRUE),
  fFillSparses(0),
  fUseRejectionMethod(kFALSE),
  fRejectionFactor(0.01),
  fhStudyImpParSingleTrackSign(0),
  fhStudyImpParSingleTrackCand(0),
  fhStudyImpParSingleTrackFd(0),
  fDetSignal(0),
  fhMultVZEROTPCoutTrackCorrNoCut(0x0),
  fhMultVZEROTPCoutTrackCorr(0x0),
  fEnablePileupRejVZEROTPCout(kFALSE),
  fhMultVZEROTPCclustersCorrNoCut(0x0),
  fhMultVZEROTPCclustersCorr(0x0),
  fEnablePileupRejVZEROTPCcls(kFALSE),
  fRejectOutOfBunchPileUp(kFALSE),
  fCut4BDTptbin(0),
  fListRDHFBDT(0),
  fListBDTNames(0),
  fListBDTNtuple(0),
  fListBDTResp(0),
  fBDTSidebandSamplingFraction(0.1),
  fSampleSideband(kFALSE),
fBDTFullVarString("ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi:multiplicity"),
  fBDTClassifierVarString(""),
fSubtractTrackletsFromDau(kFALSE),
fRefMult(9.26),
fmultiana(0),
fMultiplicityEstimator(kNtrk10),
fMCPrimariesEstimator(kEta10),
fListProfiles(0),
fCounterC(0),
fDoVZER0ParamVertexCorr(1),
fYearNumber(16)
{
  /// Default constructor
    for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;
    for(Int_t i=0; i<8; i++) h3Invmass[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_19[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_1029[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_3059[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_19999[i]=0x0;
    for(Int_t i=0; i<8; i++) h3Invmass_6099[i]=0x0;
  fBDTPtCut[0]=0; fBDTPtCut[1]=1e9;
  //~ fBDTSidebandSamplingFraction=0.01;
  //~ fBDTFullVarString="ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi";
  //~ fBDTClassifierVarString="";

  fNPtBins=cuts->GetNPtBins();

  fCuts=cuts;
  for(Int_t ih=0; ih<5; ih++) fHistMassPtImpParTC[ih]=0x0;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TList container (distributions)
  DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TH1F container (number of events)
  DefineOutput(3,TH1F::Class());  //My private output
  // Output slot #4 writes into a TList container (cuts)
  DefineOutput(4,AliRDHFCutsD0toKpi::Class());  //My private output
  // Output slot #5 writes Normalization Counter
  DefineOutput(5,AliNormalizationCounter::Class());
  // Output slot #6 
  DefineOutput(6,TList::Class());
  DefineOutput(7,TList::Class());
    DefineOutput(8,TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSED0BDT::~AliAnalysisTaskSED0BDT()
{
    for(Int_t i=0; i<14; i++) {
      if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
    }
    delete fListProfiles;
    delete fCounterC;
  if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fOutputMassPt) {
    delete fOutputMassPt;
    fOutputMassPt = 0;
  }
  if (fOutputMassY) {
    delete fOutputMassY;
    fOutputMassY = 0;
  }
  if (fDistr) {
    delete fDistr;
    fDistr = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  for(Int_t i=0; i<5; i++){
    if(fHistMassPtImpParTC[i]) delete fHistMassPtImpParTC[i];
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
  if(fVariablesTree){
    delete fVariablesTree;
    fVariablesTree = 0;
  }
    for(Int_t i=0; i<8; i++){
      if(h3Invmass[i]) delete h3Invmass[i];
    }
    for(Int_t i=0; i<8; i++){
      if(h3Invmass_19[i]) delete h3Invmass_19[i];
    }
    for(Int_t i=0; i<8; i++){
      if(h3Invmass_1029[i]) delete h3Invmass_1029[i];
    }
    for(Int_t i=0; i<8; i++){
      if(h3Invmass_3059[i]) delete h3Invmass_3059[i];
    }
    //pxy_new
    for(Int_t i=0; i<8; i++){
      if(h3Invmass_19999[i]) delete h3Invmass_19999[i];
    }
    for(Int_t i=0; i<8; i++){
      if(h3Invmass_6099[i]) delete h3Invmass_6099[i];
    }
  if (fDetSignal) {
    delete fDetSignal;
    fDetSignal = 0;
  }
  delete fMCAccPrompt;
  delete fMCAccBFeed;
  delete fhStudyImpParSingleTrackSign;
  delete fhStudyImpParSingleTrackCand;
  delete fhStudyImpParSingleTrackFd;

  if(fhMultVZEROTPCoutTrackCorrNoCut){
    delete fhMultVZEROTPCoutTrackCorrNoCut;
    fhMultVZEROTPCoutTrackCorrNoCut = 0;
  }
  if (fhMultVZEROTPCoutTrackCorr) {
    delete fhMultVZEROTPCoutTrackCorr;
    fhMultVZEROTPCoutTrackCorr = 0;
  }

  if(fhMultVZEROTPCclustersCorrNoCut){
    delete fhMultVZEROTPCclustersCorrNoCut;
    fhMultVZEROTPCclustersCorrNoCut = 0;
  }

  if(fhMultVZEROTPCclustersCorr){
      delete fhMultVZEROTPCclustersCorr;
      fhMultVZEROTPCclustersCorr = 0;
  }
  
  if (fListRDHFBDT) {
    delete fListRDHFBDT;
    fListRDHFBDT = 0;
  }
  if (fListBDTNames) {
    delete fListBDTNames;
    fListBDTNames = 0;
  }
  if (fListBDTNtuple) {
    delete fListBDTNtuple;
    fListBDTNtuple = 0;
  }
  if (fListBDTResp) {
    delete fListBDTResp;
    fListBDTResp = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSED0BDT::Init()
{
  /// Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");


  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
  copyfCuts->SetName(nameoutput);
    fListProfiles = new TList();
    fListProfiles->SetOwner();
    TString period[14];
    Int_t nProfiles=14;
    if(fYearNumber == 10){
    period[0]="LHC10b";
    period[1]="LHC10c";
    period[2]="LHC10d";
    period[3]="LHC10e";
    nProfiles = 4;
  }else if(fYearNumber == 16){
     period[0]="LHC16d";
     period[1]="LHC16e";
     period[2]="LHC16g";
     period[3]="LHC16h1";
     period[4]="LHC16h2";
     period[5]="LHC16j";
     period[6]="LHC16k";
     period[7]="LHC16l";
     period[8]="LHC16o";
     period[9]="LHC16p";
     nProfiles = 10;
  }else if(fYearNumber == 17){
    period[0]="LHC17e";
    period[1]="LHC17f";
    period[2]="LHC17h";
    period[3]="LHC17i";
    period[4]="LHC17j";
    period[5]="LHC17k";
    period[6]="LHC17l";
    period[7]="LHC17m";
    period[8]="LHC17o";
    period[9]="LHC17r";
    nProfiles = 10;
  }else if(fYearNumber == 18){
    period[0]="LHC18b";
    period[1]="LHC18d";
    period[2]="LHC18e";
    period[3]="LHC18f";
    period[4]="LHC18g";
    period[5]="LHC18h";
    period[6]="LHC18i";
    period[7]="LHC18j";
    period[8]="LHC18k";
    period[9]="LHC18l";
    period[10]="LHC18m";
    period[11]="LHC18n";
    period[12]="LHC18o";
    period[13]="LHC18p";
    nProfiles = 14;
   }
    for(Int_t i=0; i<nProfiles; i++){
      if(fMultEstimatorAvg[i]){
        TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
        hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
        fListProfiles->Add(hprof);
      }
    }
    PostData(8,fListProfiles);
  // Post the data
  PostData(4,copyfCuts);



  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0BDT::UserCreateOutputObjects()
{

  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fOutputMassPt = new TList();
  fOutputMassPt->SetOwner();
  fOutputMassPt->SetName("listMassPt");

  fOutputMassY = new TList();
  fOutputMassY->SetOwner();
  fOutputMassY->SetName("listMassY");

  fDistr = new TList();
  fDistr->SetOwner();
  fDistr->SetName("distributionslist");

  fDetSignal = new TList();
  fDetSignal->SetOwner();
  fDetSignal->SetName("detectorsignals");

  TString nameMass=" ",nameSgn27=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassNocutsS =" ",nameMassNocutsB =" ", namedistr=" ";
  TString nameMassPt="", nameSgnPt="", nameBkgPt="", nameRflPt="";
  TString nameMassY="", nameSgnY="", nameBkgY="", nameRflY="";
  TString nameMassSubSample="";
  Int_t nbins2dPt=60; Float_t binInPt=0., binFinPt=30.;
  Int_t nbins2dY=60; Float_t binInY=-1.5, binFinY=1.5;

  for(Int_t i=0;i<fCuts->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn27="histSgn27_";
    nameSgn27+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;
    nameMassNocutsS="hMassS_";
    nameMassNocutsS+=i;
    nameMassNocutsB="hMassB_";
    nameMassNocutsB+=i;

    //histograms of cut variable distributions
    if(fFillVarHists){
      if(fReadMC){

	namedistr="hNclsD0vsptS_";
	namedistr+=i;
	TH2F *hNclsD0vsptS = new TH2F(namedistr.Data(),"N cls distrubution [S];p_{T} [GeV/c];N cls",200,0.,20.,100,0.,200.);
	namedistr="hNclsD0barvsptS_";
	namedistr+=i;
	TH2F *hNclsD0barvsptS = new TH2F(namedistr.Data(),"N cls distrubution [S];p_{T} [GeV/c];N cls",200,0.,20.,100,0.,200.);

	namedistr="hNITSpointsD0vsptS_";
	namedistr+=i;
	TH2F *hNITSpointsD0vsptS = new TH2F(namedistr.Data(),"N ITS points distrubution [S];p_{T} [GeV/c];N points",200,0.,20.,7,0.,7.);

	namedistr="hNSPDpointsD0S_";
	namedistr+=i;
	TH1I *hNSPDpointsD0S = new TH1I(namedistr.Data(),"N SPD points distrubution [S]; ;N tracks",4,0,4);
	hNSPDpointsD0S->GetXaxis()->SetBinLabel(1, "no SPD");
	hNSPDpointsD0S->GetXaxis()->SetBinLabel(2, "kOnlyFirst");
	hNSPDpointsD0S->GetXaxis()->SetBinLabel(3, "kOnlySecond");
	hNSPDpointsD0S->GetXaxis()->SetBinLabel(4, "kBoth");

	namedistr="hptD0S_";
	namedistr+=i;
	TH1F *hptD0S = new TH1F(namedistr.Data(), "p_{T} distribution [S];p_{T} [GeV/c]",200,0.,20.);
	namedistr="hptD0barS_";
	namedistr+=i;
	TH1F *hptD0barS = new TH1F(namedistr.Data(), "p_{T} distribution [S];p_{T} [GeV/c]",200,0.,20.);

	namedistr="hphiD0S_";
	namedistr+=i;
	TH1F *hphiD0S = new TH1F(namedistr.Data(), "#phi distribution [S];#phi [rad]",100,0.,2*TMath::Pi());
	namedistr="hphiD0barS_";
	namedistr+=i;
	TH1F *hphiD0barS = new TH1F(namedistr.Data(), "#phi distribution [S];#phi [rad]",100,0.,2*TMath::Pi());


	namedistr="hetaphiD0candidateS_";
	namedistr+=i;
	TH2F *hetaphiD0candidateS = new TH2F(namedistr.Data(), "D^{0} candidates #eta #phi distribution [S];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());
	namedistr="hetaphiD0barcandidateS_";
	namedistr+=i;
	TH2F *hetaphiD0barcandidateS = new TH2F(namedistr.Data(), "anti-D^{0} candidates #eta #phi distribution [S];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());

	namedistr="hetaphiD0candidatesignalregionS_";
	namedistr+=i;
	TH2F *hetaphiD0candidatesignalregionS = new TH2F(namedistr.Data(), "D^{0} candidates #eta #phi distribution [S] [mass cut];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());
	namedistr="hetaphiD0barcandidatesignalregionS_";
	namedistr+=i;
	TH2F *hetaphiD0barcandidatesignalregionS = new TH2F(namedistr.Data(), "anti-D^{0} candidates #eta #phi distribution [S] [mass cut];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());

	//  dca
	namedistr="hdcaS_";
	namedistr+=i;
	TH1F *hdcaS = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
	// impact parameter
	namedistr="hd0piS_";
	namedistr+=i;
	TH1F *hd0piS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,-0.1,0.1);

	namedistr="hd0KS_";
	namedistr+=i;
	TH1F *hd0KS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,-0.1,0.1);
	namedistr="hd0d0S_";
	namedistr+=i;
	TH1F *hd0d0S = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

	//decay lenght
	namedistr="hdeclS_";
	namedistr+=i;
	TH1F *hdeclengthS = new TH1F(namedistr.Data(), "Decay Length^{2} distribution;Decay Length^{2} [cm]",200,0,0.015);

	namedistr="hnormdeclS_";
	namedistr+=i;
	TH1F *hnormdeclengthS = new TH1F(namedistr.Data(), "Normalized Decay Length^{2} distribution;(Decay Length/Err)^{2} ",200,0,12.);

	namedistr="hdeclxyS_";
	namedistr+=i;
	TH1F* hdeclxyS=new TH1F(namedistr.Data(),"Decay Length XY distribution;Decay Length XY [cm]",200,0,0.15);
	namedistr="hnormdeclxyS_";
	namedistr+=i;
	TH1F* hnormdeclxyS=new TH1F(namedistr.Data(),"Normalized decay Length XY distribution;Decay Length XY/Err",200,0,6.);

	namedistr="hdeclxyd0d0S_";
	namedistr+=i;
	TH2F* hdeclxyd0d0S=new TH2F(namedistr.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",200,0,0.15,200,-0.001,0.001);

	namedistr="hnormdeclxyd0d0S_";
	namedistr+=i;
	TH2F* hnormdeclxyd0d0S=new TH2F(namedistr.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",200,0,6,200,-0.001,0.001);

	//  costhetapoint
	namedistr="hcosthetapointS_";
	namedistr+=i;
	TH1F *hcosthetapointS = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

	namedistr="hcosthetapointxyS_";
	namedistr+=i;
	TH1F *hcosthetapointxyS = new TH1F(namedistr.Data(), "cos#theta_{Point} XYdistribution;cos#theta_{Point}",300,0.,1.);

	TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",600,1.6248,2.2248); //range (MD0-300MeV, mD0 + 300MeV)
	tmpMS->Sumw2();

	fDistr->Add(hNclsD0vsptS);
	fDistr->Add(hNclsD0barvsptS);
	fDistr->Add(hNITSpointsD0vsptS);
	fDistr->Add(hNSPDpointsD0S);
	fDistr->Add(hptD0S);
	fDistr->Add(hphiD0S);
	fDistr->Add(hptD0barS);
	fDistr->Add(hphiD0barS);
	fDistr->Add(hetaphiD0candidateS);
	fDistr->Add(hetaphiD0candidatesignalregionS);
	fDistr->Add(hetaphiD0barcandidateS);
	fDistr->Add(hetaphiD0barcandidatesignalregionS);

	fDistr->Add(hdcaS);

	fDistr->Add(hd0piS);
	fDistr->Add(hd0KS);

	fDistr->Add(hd0d0S);

	fDistr->Add(hcosthetapointS);

	fDistr->Add(hcosthetapointxyS);

	fDistr->Add(hdeclengthS);

	fDistr->Add(hnormdeclengthS);

	fDistr->Add(hdeclxyS);

	fDistr->Add(hnormdeclxyS);

	fDistr->Add(hdeclxyd0d0S);
	fDistr->Add(hnormdeclxyd0d0S);

	fDistr->Add(tmpMS);
      }


      //Ncls, phi, pt distrubutions

      namedistr="hNclsD0vsptB_";
      namedistr+=i;
      TH2F *hNclsD0vsptB = new TH2F(namedistr.Data(),"N cls distrubution [B];p_{T} [GeV/c];N cls",200,0.,20.,100,0.,200.);
      namedistr="hNclsD0barvsptB_";
      namedistr+=i;
      TH2F *hNclsD0barvsptB = new TH2F(namedistr.Data(),"N cls distrubution [B];p_{T} [GeV/c];N cls",200,0.,20.,100,0.,200.);

      namedistr="hNITSpointsD0vsptB_";
      namedistr+=i;
      TH2F *hNITSpointsD0vsptB = new TH2F(namedistr.Data(),"N ITS points distrubution [B];p_{T} [GeV/c];N points",200,0.,20.,7,0.,7.);

      namedistr="hNSPDpointsD0B_";
      namedistr+=i;
      TH1I *hNSPDpointsD0B = new TH1I(namedistr.Data(),"N SPD points distrubution [B]; ;N tracks",4,0,4);
      hNSPDpointsD0B->GetXaxis()->SetBinLabel(1, "no SPD");
      hNSPDpointsD0B->GetXaxis()->SetBinLabel(2, "kOnlyFirst");
      hNSPDpointsD0B->GetXaxis()->SetBinLabel(3, "kOnlySecond");
      hNSPDpointsD0B->GetXaxis()->SetBinLabel(4, "kBoth");

      namedistr="hptD0B_";
      namedistr+=i;
      TH1F *hptD0B = new TH1F(namedistr.Data(), "p_{T} distribution [B];p_{T} [GeV/c]",200,0.,20.);
      namedistr="hptD0barB_";
      namedistr+=i;
      TH1F *hptD0barB = new TH1F(namedistr.Data(), "p_{T} distribution [B];p_{T} [GeV/c]",200,0.,20.);

      namedistr="hphiD0B_";
      namedistr+=i;
      TH1F *hphiD0B = new TH1F(namedistr.Data(), "#phi distribution [B];#phi [rad]",100,0.,2*TMath::Pi());
      namedistr="hphiD0barB_";
      namedistr+=i;
      TH1F *hphiD0barB = new TH1F(namedistr.Data(), "#phi distribution [B];#phi [rad]",100,0.,2*TMath::Pi());

      namedistr="hetaphiD0candidateB_";
      namedistr+=i;
      TH2F *hetaphiD0candidateB = new TH2F(namedistr.Data(), "D^{0} candidates #eta #phi distribution [B];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());
      namedistr="hetaphiD0barcandidateB_";
      namedistr+=i;
      TH2F *hetaphiD0barcandidateB = new TH2F(namedistr.Data(), "anti-D^{0} candidates #eta #phi distribution [B];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());

      namedistr="hetaphiD0candidatesignalregionB_";
      namedistr+=i;
      TH2F *hetaphiD0candidatesignalregionB = new TH2F(namedistr.Data(), "D^{0} candidates #eta #phi distribution [B] [mass cut];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());
      namedistr="hetaphiD0barcandidatesignalregionB_";
      namedistr+=i;
      TH2F *hetaphiD0barcandidatesignalregionB = new TH2F(namedistr.Data(), "anti-D^{0} candidates #eta #phi distribution [B] [mass cut];#eta;#phi [rad]",100, -1.5, 1.5, 100, 0.,2*TMath::Pi());

      //  dca
      namedistr="hdcaB_";
      namedistr+=i;
      TH1F *hdcaB = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);

      // impact parameter
      namedistr="hd0B_";
      namedistr+=i;
      TH1F *hd0B = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

      namedistr="hd0d0B_";
      namedistr+=i;
      TH1F *hd0d0B = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

      //decay lenght
      namedistr="hdeclB_";
      namedistr+=i;
      TH1F *hdeclengthB = new TH1F(namedistr.Data(), "Decay Length^{2} distribution;Decay Length^{2} [cm^{2}]",200,0,0.015);

      namedistr="hnormdeclB_";
      namedistr+=i;
      TH1F *hnormdeclengthB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution;(Decay Length/Err)^{2} ",200,0,12.);

      namedistr="hdeclxyB_";
      namedistr+=i;
      TH1F* hdeclxyB=new TH1F(namedistr.Data(),"Decay Length XY distribution;Decay Length XY [cm]",200,0,0.15);
      namedistr="hnormdeclxyB_";
      namedistr+=i;
      TH1F* hnormdeclxyB=new TH1F(namedistr.Data(),"Normalized decay Length XY distribution;Decay Length XY/Err",200,0,6.);

      namedistr="hdeclxyd0d0B_";
      namedistr+=i;
      TH2F* hdeclxyd0d0B=new TH2F(namedistr.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",200,0,0.15,200,-0.001,0.001);

      namedistr="hnormdeclxyd0d0B_";
      namedistr+=i;
      TH2F* hnormdeclxyd0d0B=new TH2F(namedistr.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",200,0,6,200,-0.001,0.001);

      //  costhetapoint
      namedistr="hcosthetapointB_";
      namedistr+=i;
      TH1F *hcosthetapointB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

      namedistr="hcosthetapointxyB_";
      namedistr+=i;
      TH1F *hcosthetapointxyB = new TH1F(namedistr.Data(), "cos#theta_{Point} XY distribution;cos#theta_{Point} XY",300,0.,1.);

      TH1F* tmpMB = new TH1F(nameMassNocutsB.Data(),"D^{0} invariant mass; M [GeV]; Entries",600,1.6248,2.2248); //range (MD0-300MeV, mD0 + 300MeV)
      tmpMB->Sumw2();


      fDistr->Add(hNclsD0vsptB);
      fDistr->Add(hNclsD0barvsptB);
      fDistr->Add(hNITSpointsD0vsptB);
      fDistr->Add(hNSPDpointsD0B);
      fDistr->Add(hptD0B);
      fDistr->Add(hphiD0B);
      fDistr->Add(hptD0barB);
      fDistr->Add(hphiD0barB);
      fDistr->Add(hetaphiD0candidateB);
      fDistr->Add(hetaphiD0candidatesignalregionB);
      fDistr->Add(hetaphiD0barcandidateB);
      fDistr->Add(hetaphiD0barcandidatesignalregionB);

      fDistr->Add(hdcaB);

      fDistr->Add(hd0B);

      fDistr->Add(hd0d0B);

      fDistr->Add(hcosthetapointB);

      fDistr->Add(hcosthetapointxyB);

      fDistr->Add(hdeclengthB);

      fDistr->Add(hnormdeclengthB);

      fDistr->Add(hdeclxyB);

      fDistr->Add(hnormdeclxyB);

      fDistr->Add(hdeclxyd0d0B);
      fDistr->Add(hnormdeclxyd0d0B);

      fDistr->Add(tmpMB);

      //histograms filled only when the secondary vertex is recalculated w/o the daughter tracks (as requested in the cut object)

      if(fCuts->GetIsPrimaryWithoutDaughters()){
	if(fReadMC){
	  namedistr="hd0vpiS_";
	  namedistr+=i;
	  TH1F *hd0vpiS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions)(vtx w/o these tracks);d0(#pi) [cm]",200,-0.1,0.1);

	  namedistr="hd0vKS_";
	  namedistr+=i;
	  TH1F *hd0vKS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons) (vtx w/o these tracks);d0(K) [cm]",200,-0.1,0.1);

	  namedistr="hd0d0vS_";
	  namedistr+=i;
	  TH1F *hd0d0vS = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

	  namedistr="hdeclvS_";
	  namedistr+=i;
	  TH1F *hdeclengthvS = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

	  namedistr="hnormdeclvS_";
	  namedistr+=i;
	  TH1F *hnormdeclengthvS = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

	  fDistr->Add(hd0vpiS);
	  fDistr->Add(hd0vKS);
	  fDistr->Add(hd0d0vS);
	  fDistr->Add(hdeclengthvS);
	  fDistr->Add(hnormdeclengthvS);

	}

	namedistr="hd0vmoresB_";
	namedistr+=i;
	TH1F *hd0vmoresB = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0vmoresB_";
	namedistr+=i;
	TH1F *hd0d0vmoresB = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.001,0.001);


	namedistr="hd0vB_";
	namedistr+=i;
	TH1F *hd0vB = new TH1F(namedistr.Data(), "Impact parameter distribution (vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0vp0B_";
	namedistr+=i;
	TH1F *hd0vp0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong + ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0vp1B_";
	namedistr+=i;
	TH1F *hd0vp1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong - ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0vB_";
	namedistr+=i;
	TH1F *hd0d0vB = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

	namedistr="hdeclvB_";
	namedistr+=i;
	TH1F *hdeclengthvB = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

	namedistr="hnormdeclvB_";
	namedistr+=i;
	TH1F *hnormdeclengthvB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

	fDistr->Add(hd0vB);
	fDistr->Add(hd0vp0B);
	fDistr->Add(hd0vp1B);
	fDistr->Add(hd0vmoresB);

	fDistr->Add(hd0d0vB);
	fDistr->Add(hd0d0vmoresB);

	fDistr->Add(hdeclengthvB);

	fDistr->Add(hnormdeclengthvB);
      }


    }
    //histograms of invariant mass distributions


    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",600,1.6248,2.2248);

      TH1F *tmpSl=(TH1F*)tmpSt->Clone();
      tmpSt->Sumw2();
      tmpSl->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",600,1.6248,2.2248);
      //TH1F *tmpRl=(TH1F*)tmpRt->Clone();
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",600,1.6248,2.2248);
      //TH1F *tmpBl=(TH1F*)tmpBt->Clone();
      tmpBt->Sumw2();
      //tmpBl->Sumw2();
      tmpRt->Sumw2();
      //tmpRl->Sumw2();

      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpBt);

    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",600,1.6248,2.2248);
    //TH1F *tmpMl=(TH1F*)tmpMt->Clone();
    tmpMt->Sumw2();
    //tmpMl->Sumw2();
    //distribution w/o cuts large range
    // TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,0.7,3.);

    fOutputMass->Add(tmpMt);

    //sub sample
    if(fFillSubSampleHist){
      nameMassSubSample="histMassvsSubSample_";
      nameMassSubSample+=i;
      TH2F* tmpMtSub = new TH2F(nameMassSubSample.Data(),"D^{0} invariant mass; M [GeV]; Sample ID; Entries",600,1.6248,2.2248,25,-0.5,24.5);
      tmpMtSub->Sumw2();
      fOutputMass->Add(tmpMtSub);
    }

    if(fSys==0){ //histograms filled only in pp to save time in PbPb
      if(fFillVarHists){

	if(fReadMC){
	  //  pT
	  namedistr="hptpiS_";
	  namedistr+=i;
	  TH1F *hptpiS = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

	  namedistr="hptKS_";
	  namedistr+=i;
	  TH1F *hptKS = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

	  //  costhetastar
	  namedistr="hcosthetastarS_";
	  namedistr+=i;
	  TH1F *hcosthetastarS = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

	  //pT no mass cut

	  namedistr="hptpiSnoMcut_";
	  namedistr+=i;
	  TH1F *hptpiSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

	  namedistr="hptKSnoMcut_";
	  namedistr+=i;
	  TH1F *hptKSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

	  fDistr->Add(hptpiS);
	  fDistr->Add(hptKS);
	  fDistr->Add(hcosthetastarS);

	  fDistr->Add(hptpiSnoMcut);
	  fDistr->Add(hptKSnoMcut);

	  //  costhetapoint vs d0 or d0d0
	  namedistr="hcosthpointd0S_";
	  namedistr+=i;
	  TH2F *hcosthpointd0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);
	  namedistr="hcosthpointd0d0S_";
	  namedistr+=i;
	  TH2F *hcosthpointd0d0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	  fDistr->Add(hcosthpointd0S);
	  fDistr->Add(hcosthpointd0d0S);

	  //to compare with AliAnalysisTaskCharmFraction
	  TH1F* tmpS27t = new TH1F(nameSgn27.Data(),"D^{0} invariant mass in M(D^{0}) +/- 27 MeV - MC; M [GeV]; Entries",600,1.6248,2.2248);
	  TH1F *tmpS27l=(TH1F*)tmpS27t->Clone();
	  tmpS27t->Sumw2();
	  tmpS27l->Sumw2();

	  fOutputMass->Add(tmpS27t);
	  fOutputMass->Add(tmpS27l);

	}

	//  pT
	namedistr="hptB_";
	namedistr+=i;
	TH1F *hptB = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

	//  costhetastar
	namedistr="hcosthetastarB_";
	namedistr+=i;
	TH1F *hcosthetastarB = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

	//pT no mass cut
	namedistr="hptB1prongnoMcut_";
	namedistr+=i;
	TH1F *hptB1pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

	namedistr="hptB2prongsnoMcut_";
	namedistr+=i;
	TH1F *hptB2pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

	fDistr->Add(hptB);
	fDistr->Add(hcosthetastarB);

	fDistr->Add(hptB1pnoMcut);
	fDistr->Add(hptB2pnoMcut);

	//impact parameter of negative/positive track
	namedistr="hd0p0B_";
	namedistr+=i;
	TH1F *hd0p0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0p1B_";
	namedistr+=i;
	TH1F *hd0p1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong -);d0 [cm]",200,-0.1,0.1);

	//impact parameter corrected for strangeness
	namedistr="hd0moresB_";
	namedistr+=i;
	TH1F *hd0moresB = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0moresB_";
	namedistr+=i;
	TH1F *hd0d0moresB = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.001,0.001);


	namedistr="hcosthetapointmoresB_";
	namedistr+=i;
	TH1F *hcosthetapointmoresB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

	//  costhetapoint vs d0 or d0d0
	namedistr="hcosthpointd0B_";
	namedistr+=i;
	TH2F *hcosthpointd0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	namedistr="hcosthpointd0d0B_";
	namedistr+=i;
	TH2F *hcosthpointd0d0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	fDistr->Add(hd0p0B);
	fDistr->Add(hd0p1B);

	fDistr->Add(hd0moresB);
	fDistr->Add(hd0d0moresB);
	fDistr->Add(hcosthetapointmoresB);


	fDistr->Add(hcosthpointd0B);


	fDistr->Add(hcosthpointd0d0B);
      }

    } //end pp histo
  }
  
  if(fFillSparses){
	if(fReadMC){
      //pt, ptB, normImpParTrk1, normImpParTrk2, decLXY, normDecLXY, iscut, ispid
      Int_t nbinsImpParStudy[8]=      {50,50,40, 40, 20,  15,  3, 4};
      Double_t limitLowImpParStudy[8]={0, 0, -5,-5., 0.,  0.,  1.,0.};
      Double_t limitUpImpParStudy[8]= {50.,50., 5, 5,  0.2, 15,  3.,4.};

      fhStudyImpParSingleTrackSign=new THnSparseF("fhStudyImpParSingleTrackSign","fhStudyImpParSingleTrackSign",8,nbinsImpParStudy,limitLowImpParStudy,limitUpImpParStudy);
      TString axTitMC[8]={"#it{p}_{T} (GeV/c)","#it{p}_{T} (GeV/c)","normalized imp par residual, trk1","normalized imp par residual, trk2","#it{L}_{xy} (cm)","norm #it{L}_{xy}","cutSel","PIDinfo"};
      for(Int_t iax=0; iax<8; iax++) fhStudyImpParSingleTrackSign->GetAxis(iax)->SetTitle(axTitMC[iax].Data());
      fOutputMass->Add(fhStudyImpParSingleTrackSign);


      fhStudyImpParSingleTrackFd=new THnSparseF("fhStudyImpParSingleTrackFd","fhStudyImpParSingleTrackFd",8,nbinsImpParStudy,limitLowImpParStudy,limitUpImpParStudy);
      for(Int_t iax=0; iax<8; iax++) fhStudyImpParSingleTrackFd->GetAxis(iax)->SetTitle(axTitMC[iax].Data());
      fOutputMass->Add(fhStudyImpParSingleTrackFd);

	  if(fStepMCAcc) CreateMCAcceptanceHistos();
    }
    //create THnSparse for impact parameter analysis in DATA sample.
    //pt, normImpParTrk1, normImpParTrk2,    decLXY, normDecLXY, massd0, cut, pid, D0D0bar
    else{
		Int_t nbinsImpParStudy[9]=      {50, 40, 40, 20, 15, 600,  3,4,2};
		Double_t limitLowImpParStudy[9]={0,  -5, -5,  0,  0,1.6248,1,0,0};
		Double_t limitUpImpParStudy[9]= {50., 5,  5, 0.2,15,2.2248,4,4,2};
		TString axTit[9]={"#it{p}_{T} (GeV/c)","normalized imp par residual, trk1","normalized imp par residual, trk2","#it{L}_{xy} (cm)","norm #it{L}_{xy}","MassD0_{k#pi} (GeV/#it{c}^{2})", "cutSel","PIDinfo","D0D0bar"};
		fhStudyImpParSingleTrackCand=new THnSparseF("fhStudyImpParSingleTrackCand","fhStudyImpParSingleTrackCand",9,nbinsImpParStudy,limitLowImpParStudy,limitUpImpParStudy);
		for(Int_t iax=0; iax<9; iax++) fhStudyImpParSingleTrackCand->GetAxis(iax)->SetTitle(axTit[iax].Data());
		fOutputMass->Add(fhStudyImpParSingleTrackCand);
	}
  }
  //for Like sign analysis
  if(fArray==1){
    namedistr="hpospair";
    TH1F* hpospair=new TH1F(namedistr.Data(),"Number of positive pairs",fCuts->GetNPtBins(),-0.5,fCuts->GetNPtBins()-0.5);
    namedistr="hnegpair";
    TH1F* hnegpair=new TH1F(namedistr.Data(),"Number of negative pairs",fCuts->GetNPtBins(),-0.5,fCuts->GetNPtBins()-0.5);
    fOutputMass->Add(hpospair);
    fOutputMass->Add(hnegpair);
  }


  // 2D Pt distributions and impact parameter histograms
  if(fFillPtHist) {

    nameMassPt="histMassPt";
    nameSgnPt="histSgnPt";
    nameBkgPt="histBkgPt";
    nameRflPt="histRflPt";

    //MC signal
    if(fReadMC){
      TH2F* tmpStPt = new TH2F(nameSgnPt.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries; Pt[GeV/c]",600,1.6248,2.2248,nbins2dPt,binInPt,binFinPt);
      TH2F *tmpSlPt=(TH2F*)tmpStPt->Clone();
      tmpStPt->Sumw2();
      tmpSlPt->Sumw2();

      //Reflection: histo filled with D0MassV1 which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH2F* tmpRtPt = new TH2F(nameRflPt.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries; Pt[GeV/c]",600,1.6248,2.2248,nbins2dPt,binInPt,binFinPt);
      TH2F* tmpBtPt = new TH2F(nameBkgPt.Data(), "Background invariant mass - MC; M [GeV]; Entries; Pt[GeV/c]",600,1.6248,2.2248,nbins2dPt,binInPt,binFinPt);
      tmpBtPt->Sumw2();
      tmpRtPt->Sumw2();

      fOutputMassPt->Add(tmpStPt);
      fOutputMassPt->Add(tmpRtPt);
      fOutputMassPt->Add(tmpBtPt);

      //       cout<<endl<<endl<<"***************************************"<<endl;
      //       cout << " created and added to the list "<< nameSgnPt.Data() <<" "<< tmpStPt <<
      // 	", "<<nameRflPt.Data() <<" " << tmpRtPt<<", "<<nameBkgPt.Data()<< " " << tmpBtPt <<endl;
      //       cout<<"***************************************"<<endl<<endl;
    }

    TH2F* tmpMtPt = new TH2F(nameMassPt.Data(),"D^{0} invariant mass; M [GeV]; Entries; Pt[GeV/c]",600,1.6248,2.2248,nbins2dPt,binInPt,binFinPt);
    tmpMtPt->Sumw2();

    fOutputMassPt->Add(tmpMtPt);
  }

  if(fFillImpParHist) CreateImpactParameterHistos();

  // 2D Y distributions

  if(fFillYHist) {
    for(Int_t i=0;i<fCuts->GetNPtBins();i++){
      nameMassY="histMassY_";
      nameMassY+=i;
      nameSgnY="histSgnY_";
      nameSgnY+=i;
      nameBkgY="histBkgY_";
      nameBkgY+=i;
      nameRflY="histRflY_";
      nameRflY+=i;
      //MC signal
      if(fReadMC){
	TH2F* tmpStY = new TH2F(nameSgnY.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries; y",600,1.6248,2.2248,nbins2dY,binInY,binFinY);
	tmpStY->Sumw2();
	//Reflection: histo filled with D0MassV1 which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
	TH2F* tmpRtY = new TH2F(nameRflY.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries; y",600,1.6248,2.2248,nbins2dY,binInY,binFinY);
	TH2F* tmpBtY = new TH2F(nameBkgY.Data(), "Background invariant mass - MC; M [GeV]; Entries; y",600,1.6248,2.2248,nbins2dY,binInY,binFinY);
	tmpBtY->Sumw2();
	tmpRtY->Sumw2();

	fOutputMassY->Add(tmpStY);
	fOutputMassY->Add(tmpRtY);
	fOutputMassY->Add(tmpBtY);
      }
      TH2F* tmpMtY = new TH2F(nameMassY.Data(),"D^{0} invariant mass; M [GeV]; Entries; y",600,1.6248,2.2248,nbins2dY,binInY,binFinY);
      tmpMtY->Sumw2();
      fOutputMassY->Add(tmpMtY);
    }
  }


  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 25,-0.5,24.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  else fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  if(fFillVarHists || fPIDCheck){
    fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
    fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
    fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
    fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
  }
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetBinLabel(18,"Phys.Sel.Rej");
  fNentries->GetXaxis()->SetBinLabel(19,"D0 failed to be filled");
  fNentries->GetXaxis()->SetBinLabel(20,"fisFilled is 0");
  fNentries->GetXaxis()->SetBinLabel(21,"fisFilled is 1");
  fNentries->GetXaxis()->SetBinLabel(22,"AOD/dAOD mismatch");
  fNentries->GetXaxis()->SetBinLabel(23,"AOD/dAOD #events ok");
  fNentries->GetXaxis()->SetBinLabel(24,"PreSelect rejection");
  if(fEnableCentralityCorrCuts)fNentries->GetXaxis()->SetBinLabel(25,"n. rejected for bad cent corr");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
  fCounter->Init();

  //
  // Output slot 7 : tree of the candidate variables after track selection - OVERRIDDEN
  //
  nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  if(fWriteVariableTree && fWriteProtosgnVar){
    AliFatal("FATAL_ERROR: Writing candidate variables both with and without --> CHOOSE ONE OF THE TWO OPTIONS!");
  }
  Int_t nVar = 15;
  TString * fCandidateVariableNames = 0x0;
  if(fWriteVariableTree){
    fCandidateVariables = new Double_t [nVar];
    fCandidateVariableNames = new TString[nVar];

    fCandidateVariableNames[0] = "massD0";
    fCandidateVariableNames[1] = "massD0bar";
    fCandidateVariableNames[2] = "pt";
    fCandidateVariableNames[3] = "dca";
    fCandidateVariableNames[4] = "costhsD0";
    fCandidateVariableNames[5] = "costhsD0bar";
    fCandidateVariableNames[6] = "ptk";
    fCandidateVariableNames[7] = "ptpi";
    fCandidateVariableNames[8] = "d0k";
    fCandidateVariableNames[9] = "d0pi";
    fCandidateVariableNames[10] = "d0xd0";
    fCandidateVariableNames[11] = "costhp";
    fCandidateVariableNames[12] = "costhpxy";
    fCandidateVariableNames[13] = "lxy";
    fCandidateVariableNames[14] = "specialcuts";
    for(Int_t ivar=0; ivar<nVar; ivar++){
      fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/d",fCandidateVariableNames[ivar].Data()));
    }
  }
  if(fWriteProtosgnVar){
    if(fReadMC) nVar = 17;
    else nVar = 16;
    fCandidateVariables = new Double_t [nVar];
    fCandidateVariableNames = new TString[nVar];

    fCandidateVariableNames[0] = "massD0";
    fCandidateVariableNames[1] = "massD0bar";
    fCandidateVariableNames[2] = "pt";
    fCandidateVariableNames[3] = "dca";
    fCandidateVariableNames[4] = "costhsD0";
    fCandidateVariableNames[5] = "costhsD0bar";
    fCandidateVariableNames[6] = "ptk";
    fCandidateVariableNames[7] = "ptpi";
    fCandidateVariableNames[8] = "d0k";
    fCandidateVariableNames[9] = "d0pi";
    fCandidateVariableNames[10] = "costhp";
    fCandidateVariableNames[11] = "costhpxy";
    fCandidateVariableNames[12] = "lxy";
    fCandidateVariableNames[13] = "specialcuts";
    fCandidateVariableNames[14] = "topomatic";
    fCandidateVariableNames[15] = "candidatetype"; //0 = D0 only; 1 = D0bar only; 2 = D0 and D0bar
    if(fReadMC)fCandidateVariableNames[16] = "ispromptc";
    for(Int_t ivar=0; ivar<nVar; ivar++){
      fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/d",fCandidateVariableNames[ivar].Data()));
    }
  }



  //
  // Output slot 8 : List for detector response histograms - OVERRIDDEN
  //
  if (fDrawDetSignal) {
    TH2F *TOFSigBefPID = new TH2F("TOFSigBefPID", "TOF signal of daughters before PID;p(daught)(GeV/c);Signal", 500, 0, 10, 5000, 0, 50e3);
    TH2F *TOFSigAftPID = new TH2F("TOFSigAftPID", "TOF signal after PID;p(daught)(GeV/c);Signal", 500, 0, 10, 5000, 0, 50e3);

    TH2F *TPCSigBefPID = new TH2F("TPCSigBefPID", "TPC dE/dx before PID;p(daught)(GeV/c);dE/dx (arb. units)", 1000, 0, 10, 1000, 0, 500);
    TH2F *TPCSigAftPID = new TH2F("TPCSigAftPID", "TPC dE/dx after PID;p(daught)(GeV/c);dE/dx (arb. units)", 1000, 0, 10, 1000, 0, 500);

    fDetSignal->Add(TOFSigBefPID);
    fDetSignal->Add(TOFSigAftPID);
    fDetSignal->Add(TPCSigBefPID);
    fDetSignal->Add(TPCSigAftPID);
  }


  fhMultVZEROTPCoutTrackCorrNoCut = new TH2F("hMultVZEROTPCoutTrackCorrNoCut", ";Tracks with kTPCout on;VZERO multiplicity", 1000, 0., 30000., 1000, 0., 40000.);
  fhMultVZEROTPCoutTrackCorr = new TH2F("hMultVZEROTPCoutTrackCorr", ";Tracks with kTPCout on;VZERO multiplicity", 1000, 0., 30000., 1000, 0., 40000.);
  fDistr->Add(fhMultVZEROTPCoutTrackCorrNoCut);
  fDistr->Add(fhMultVZEROTPCoutTrackCorr);

  fhMultVZEROTPCclustersCorrNoCut = new TH2F("fhMultVZEROTPCclustersCorrNoCut","; n^{o} TPC clusters; VZERO multiplicity", 1000,0.,10000.,6000,0.,60000.);
  fhMultVZEROTPCclustersCorr = new TH2F("fhMultVZEROTPCclustersCorr","; n^{o} TPC clusters; VZERO multiplicity", 1000,0.,10000.,6000,0.,60000.);
  fDistr->Add(fhMultVZEROTPCclustersCorrNoCut);
  fDistr->Add(fhMultVZEROTPCclustersCorr);

  if(fEnableCentralityCorrCuts){
    fEventCuts.AddQAplotsToList(fDetSignal,true);
  }
  // BDT I/O
  fListBDTNtuple = new TList(); fListBDTNtuple->SetOwner(); fListBDTNtuple->SetName("NtupleList");
  fListBDTResp = new TList(); fListBDTResp->SetOwner(); fListBDTResp->SetName("BDTResponseList");
  if(fFillSparses){
	//"ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi"
	// MC !NOT TESTED YET!
      fHistNtrCorrEvSel = new TH1F("hNtrCorrEvSel",Form("Corrected  multiplicity for selected events ; Entries"),200,0,200);
      fOutputMass->Add(fHistNtrCorrEvSel);
	if(fReadMC){
		TNtuple *NtupleD0C = new TNtuple("NtupleD0C", "MC Prompt D0", fBDTFullVarString);
		TNtuple *NtupleD0B = new TNtuple("NtupleD0B", "MC Non-prompt D0", fBDTFullVarString);
		TNtuple *NtupleRefl = new TNtuple("NtupleRefl", "MC Non-prompt D0 Reflection", fBDTFullVarString);
		
		fListBDTNtuple->Add(NtupleD0C);
		fListBDTNtuple->Add(NtupleD0B);
		fListBDTNtuple->Add(NtupleRefl);
		
	}
	else if(fSampleSideband){
		TNtuple *NtupleSB = new TNtuple("NtupleSB","Sideband",fBDTFullVarString);
		fListBDTNtuple->Add(NtupleSB);
	}
	else{
		fListRDHFBDT->SetOwner(); fListBDTNtuple->SetName("BDTList");
		for(Int_t ii=0;ii<fCut4BDTptbin->GetNPtBins();ii++){
            const Int_t NBDT = fListBDTNames->GetEntries() - 1;
            for(Int_t jj=0;jj<NBDT;jj++){
            TString BDT1Name = fListBDTNames->At(0)->GetName();
            TString BDT2Name = fListBDTNames->At(jj+1)->GetName();
            h3Invmass[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);
            h3Invmass_19[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s_19",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);
            h3Invmass_1029[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s_1029",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);
            h3Invmass_3059[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s_3059",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);
                h3Invmass_19999[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s_19999",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);
                h3Invmass_6099[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s_6099",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.6248,2.2248,80,-0.15,0.25,60,-0.05,0.25);

            //h3Invmass[jj] = new TH3F(Form("h3MassRespPt%d_%s_%s",ii,BDT1Name.Data(),BDT2Name.Data()),"Invmass",100,1.68,2.10,200,-1,1,200,-1,1);
            fListBDTResp->Add(h3Invmass[jj]);
            fListBDTResp->Add(h3Invmass_19[jj]);
            fListBDTResp->Add(h3Invmass_1029[jj]);
            fListBDTResp->Add(h3Invmass_3059[jj]);
                fListBDTResp->Add(h3Invmass_19999[jj]);
                fListBDTResp->Add(h3Invmass_6099[jj]);
        }
	}
  }
  }
    fCounterC = new AliNormalizationCounter("NormCounterCorrMult");
    fCounterC->SetStudyMultiplicity(kTRUE,1.);
    fCounterC->Init();
    fDistr->Add(fCounterC);
  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(5,fCounter);
  PostData(6,fListBDTNtuple);
  PostData(7,fListBDTResp);
  return;

}
//________________________________________________________________________
void AliAnalysisTaskSED0BDT::UserExec(Option_t */*option*/)

{
  /// Execute analysis for current event:
  /// heavy flavor candidates association to MC truth

  //cuts order
  //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
  //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
  //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
  //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
  //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
  //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
  //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(21);
      return;
    }
    fNentries->Fill(22);
  }

  TString bname;
  if(fArray==0){ //D0 candidates
    // load D0->Kpi candidates
    //cout<<"D0 candidates"<<endl;
    bname="D0toKpi";
  } else { //LikeSign candidates
    // load like sign candidates
    //cout<<"LS candidates"<<endl;
    bname="LikeSign2Prong";
  }
  TClonesArray *inputArray=0;
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else if(aod) {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }

  if(!inputArray || !aod) {
    printf("AliAnalysisTaskSED0BDT::UserExec: input branch not found!\n");
    return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

    Float_t countCorr = -99999;
    if(fmultiana){
    Int_t countTreta1=0, countTreta03=0, countTreta05=0, countTreta16=0;
    AliAODTracklets* tracklets=aod->GetTracklets();
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>-0.3 && eta<0.3) countTreta03++;
      if(eta>-0.5 && eta<0.5) countTreta05++;
      if(eta>-1.0 && eta<1.0) countTreta1++;
      if(eta>-1.6 && eta<1.6) countTreta16++;
    }
    
    Int_t vzeroMult=0, vzeroMultA=0, vzeroMultC=0;
    Int_t vzeroMultEq=0, vzeroMultAEq=0, vzeroMultCEq=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
    if(vzeroAOD) {
      vzeroMultA = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
      vzeroMultC = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
      vzeroMult = vzeroMultA + vzeroMultC;
      vzeroMultAEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aod));
      vzeroMultCEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aod));
      vzeroMultEq = vzeroMultAEq + vzeroMultCEq;
    }

    Int_t countMult = countTreta1;
    if(fMultiplicityEstimator==kNtrk03) { countMult = countTreta03; }
    else if(fMultiplicityEstimator==kNtrk05) { countMult = countTreta05; }
    else if(fMultiplicityEstimator==kNtrk10to16) { countMult = countTreta16 - countTreta1; }
    else if(fMultiplicityEstimator==kVZERO) { countMult = vzeroMult; }
    else if(fMultiplicityEstimator==kVZEROA) { countMult = vzeroMultA; }
    else if(fMultiplicityEstimator==kVZEROEq) { countMult = vzeroMultEq; }
    else if(fMultiplicityEstimator==kVZEROAEq) { countMult = vzeroMultAEq; }
    
    Double_t countTreta1corr=countTreta1;
    countCorr=countMult;
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    // In case of VZERO multiplicity, consider the zvtx correction flag
    //  fDoVZER0ParamVertexCorr: 0= none, 1= usual d2h, 2=AliESDUtils
    Bool_t isDataDrivenZvtxCorr=kTRUE;
    Bool_t isVtxOk=kFALSE;
    Int_t vzeroMultACorr=vzeroMultA, vzeroMultCCorr=vzeroMultC, vzeroMultCorr=vzeroMult;
    Int_t vzeroMultAEqCorr=vzeroMultAEq, vzeroMultCEqCorr=vzeroMultCEq, vzeroMultEqCorr=vzeroMultEq;

    if(isVtxOk){
      if( (fMultiplicityEstimator==kVZERO) || (fMultiplicityEstimator==kVZEROA) ||
      (fMultiplicityEstimator==kVZEROEq) || (fMultiplicityEstimator==kVZEROAEq) ){
        if(fDoVZER0ParamVertexCorr==0){
      // do not correct
      isDataDrivenZvtxCorr=kFALSE;
        } else if (fDoVZER0ParamVertexCorr==2){
      // use AliESDUtils correction
      Float_t zvtx = vtx1->GetZ();
      isDataDrivenZvtxCorr=kFALSE;
      vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultA,zvtx));
      vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(vzeroMultC,zvtx));
      vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr;
      vzeroMultAEqCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultAEq,zvtx));
      vzeroMultCEqCorr =static_cast<Int_t>( AliESDUtils::GetCorrV0C(vzeroMultCEq,zvtx));
      vzeroMultEqCorr = vzeroMultAEqCorr + vzeroMultCEqCorr;
      if(fMultiplicityEstimator==kVZERO) { countCorr = vzeroMultCorr; }
      else if(fMultiplicityEstimator==kVZEROA) { countCorr = vzeroMultACorr; }
      else if(fMultiplicityEstimator==kVZEROEq) { countCorr = vzeroMultEqCorr; }
      else if(fMultiplicityEstimator==kVZEROAEq) { countCorr = vzeroMultAEqCorr; }
        }
      }
    }
    // Data driven multiplicity z-vertex correction
    if(isVtxOk && isDataDrivenZvtxCorr){
      TProfile* estimatorAvg = GetEstimatorHistogram(aod);
      if(estimatorAvg){
        countTreta1corr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countTreta1,vtx1->GetZ(),fRefMult));
        // vzeroMultACorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,vzeroMultA,vtx1->GetZ(),fRefMult));
        // vzeroMultCorr= vzeroMultACorr + static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,vzeroMultC,vtx1->GetZ(),fRefMult));
        // vzeroMultAEqCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,vzeroMultAEq,vtx1->GetZ(),fRefMult));
        // vzeroMultEqCorr= vzeroMultAEqCorr + static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,vzeroMultCEq,vtx1->GetZ(),fRefMult));
        countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,vtx1->GetZ(),fRefMult));
      }
    }

    }
    fCounterC->StoreEvent(aod,fCuts,fReadMC,countCorr);
    Float_t multForCand = countCorr;
    fHistNtrCorrEvSel->Fill(countCorr,1);
    
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0BDT::UserExec: MC particles branch not found!\n");
      return;
    }

    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0BDT::UserExec: MC header branch not found!\n");
      return;
    }
  }
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());

  Int_t nTPCout=0;
  Int_t nTPCcls=aod->GetNumberOfTPCClusters();
  Float_t mTotV0=0;
  AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)aod)->GetVZEROData();
  Float_t mTotV0A=v0data->GetMTotV0A();
  Float_t mTotV0C=v0data->GetMTotV0C();
  mTotV0=mTotV0A+mTotV0C;
  Int_t ntracksEv = aod->GetNumberOfTracks();
  for(Int_t itrack=0; itrack<ntracksEv; itrack++) { // loop on tacks
    //    ... get the track
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
    if(!track) {AliFatal("Not a standard AOD");}
    if(track->GetID()<0)continue;
    if((track->GetFlags())&(AliESDtrack::kTPCout)) nTPCout++;
    else continue;
  }
  if(fhMultVZEROTPCoutTrackCorrNoCut) fhMultVZEROTPCoutTrackCorrNoCut->Fill(nTPCout,mTotV0);
  Float_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout);
  if(fEnablePileupRejVZEROTPCout){
    if(mTotV0<mV0Cut) return;
  }

  if(fhMultVZEROTPCclustersCorrNoCut) fhMultVZEROTPCclustersCorrNoCut->Fill(nTPCcls,mTotV0);
  Float_t mV0TPCclsCut=-2000.+(0.013*nTPCcls)+(1.25e-9*nTPCcls*nTPCcls);
  if(fEnablePileupRejVZEROTPCcls){ // this pile-up rejection is specific for 2018 Pb-Pb analysis
    if(fRejectOutOfBunchPileUp && mTotV0<mV0TPCclsCut) return;
    else if(!fRejectOutOfBunchPileUp && mTotV0>mV0TPCclsCut) return; //keep only out-of-bunch pile-up events
  }
  if(fhMultVZEROTPCclustersCorr) fhMultVZEROTPCclustersCorr->Fill(nTPCcls,mTotV0);

  //histogram filled with 1 for every AOD
  fNentries->Fill(0);


  fCounter->StoreEvent(aod,fCuts,fReadMC);
  //fCounter->StoreEvent(aod,fReadMC);
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
  if(fReadMC && fStepMCAcc){
    FillMCAcceptanceHistos(mcArray, mcHeader);
  }

  if(!fCuts->IsEventSelected(aod)) {
    if(fCuts->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCuts->GetWhyRejection()==2 || fCuts->GetWhyRejection()==3)) fNentries->Fill(15);
    if(fCuts->GetWhyRejection()==7) fNentries->Fill(17);
    return;
  }
  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) {
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNumberOfTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
      if(!track) AliFatal("Not a standard AOD");
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(16);
	break;
      }
    }
    if (skipEvent) return;
  }


  if(fhMultVZEROTPCoutTrackCorr)fhMultVZEROTPCoutTrackCorr->Fill(nTPCout,mTotV0);


  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Bool_t isGoodVtx=kFALSE;

  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }

  if(fEnableCentralityCorrCuts){
    fEventCuts.AcceptEvent(aod);
    if(!fEventCuts.PassedCut(AliEventCuts::kCorrelations) || !fEventCuts.PassedCut(AliEventCuts::kMultiplicity)){
      fNentries->Fill(24);
      return;
    }
  }

  fEventCounter++;

  if(fSys==1 && (fWriteProtosgnVar || fWriteVariableTree)){
    if(!fUseRejectionMethod) AliFatal("FATAL: Fill candidates variables in Pb-Pb w/o Rejection Method --> ACTIVATE IT!");
    if(fUseRejectionMethod) AliWarning("WARNING: Fill candidates variables in Pb-Pb with Rejection Method --> Check a proper rejection factor is used!");
  }

  // loop over candidates
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  // FILE *f=fopen("4display.txt","a");
  // printf("Number of D0->Kpi: %d\n",nInD0toKpi);
  Int_t nSelectedloose=0,nSelectedtight=0;

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
      if(fSubtractTrackletsFromDau){
       
      // For the D* case, subtract only the D0 daughter tracks <=== FIXME !!


      for(Int_t iDau=0; iDau<2; iDau++){
        AliAODTrack *t = NULL;

     t = (AliAODTrack*)d->GetDaughter(iDau);
        if(!t) continue;
        if(t->HasPointOnITSLayer(0) && t->HasPointOnITSLayer(1)){
          if(multForCand>0) multForCand-=1;
        }
      }
        
      }
    if(fUseSelectionBit && d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
	fNentries->Fill(2);
	continue; //skip the D0 from Dstar
      }
    if(d->GetIsFilled()==0)fNentries->Fill(19);//tmp check
    if(d->GetIsFilled()==1)fNentries->Fill(20);//tmp check
    TObjArray arrTracks(2);
    for(Int_t ipr=0;ipr<2;ipr++){
      AliAODTrack *tr=vHF->GetProng(aod,d,ipr);
      arrTracks.AddAt(tr,ipr);
    }

    if(!fCuts->PreSelect(arrTracks)){
      fNentries->Fill(23);
      continue;
    }
    if(!(vHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.
      fNentries->Fill(18); //monitor how often this fails
      continue;
    }

    if ( fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(421)) ) {
      nSelectedloose++;
      nSelectedtight++;
      if(fSys==0){
	if(fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod))fNentries->Fill(6);
      }
      Int_t ptbin=fCuts->PtBin(d->Pt());
      if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds
      fIsSelectedCandidate=fCuts->IsSelected(d,AliRDHFCuts::kAll,aod); //selected

      if(fUseRejectionMethod){
        if ((d->Pt() * 1000.) - (Int_t)(d->Pt() * 1000) > fRejectionFactor)
        continue;
      }

      if(fFillVarHists) {
        //if(!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate)) {
        fDaughterTracks.AddAt((AliAODTrack*)d->GetDaughter(0),0);
	      fDaughterTracks.AddAt((AliAODTrack*)d->GetDaughter(1),1);
	      //check daughters
	      if(!fDaughterTracks.UncheckedAt(0) || !fDaughterTracks.UncheckedAt(1)) {
	      AliDebug(1,"at least one daughter not found!");
	      fNentries->Fill(5);
	      fDaughterTracks.Clear();
	      continue;
      }
	    //}
	    FillVarHists(aod,d,mcArray,fCuts,fDistr);
    }

      if (fDrawDetSignal) {
        DrawDetSignal(d, fDetSignal);
      }

      FillMassHists(d,mcArray,mcHeader,fCuts,fOutputMass);
      FillCandVariables(aod,d,mcArray,mcHeader,fCuts);
        if(fFillSparses)  ProcessBDT(aod, d,mcArray, multForCand);
      if (fPIDCheck) {
        Int_t isSelectedPIDfill = 3;
	      if (!fReadMC || (fReadMC && fUsePid4Distr)) isSelectedPIDfill = fCuts->IsSelectedPID(d); //0 rejected,1 D0,2 Dobar, 3 both
	      if (isSelectedPIDfill == 0)fNentries->Fill(7);
        if (isSelectedPIDfill == 1)fNentries->Fill(8);
	      if (isSelectedPIDfill == 2)fNentries->Fill(9);
	      if (isSelectedPIDfill == 3)fNentries->Fill(10);
      }
    }

    fDaughterTracks.Clear();
    //if(unsetvtx) d->UnsetOwnPrimaryVtx();
  } //end for prongs
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);
  delete vHF;
  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(5,fCounter);
  PostData(6,fListBDTNtuple);
  PostData(7,fListBDTResp);

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0BDT::DrawDetSignal(AliAODRecoDecayHF2Prong *part, TList *ListDetSignal)
{
  //
  /// Function called in UserExec for drawing detector signal histograms:
  //
  fDaughterTracks.AddAt((AliAODTrack*)part->GetDaughter(0), 0);
  fDaughterTracks.AddAt((AliAODTrack*)part->GetDaughter(1), 1);

  AliESDtrack *esdtrack1 = new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(0));
  AliESDtrack *esdtrack2 = new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(1));


  //For filling detector histograms
  Int_t isSelectedPIDfill = 3;
  if (!fReadMC || (fReadMC && fUsePid4Distr)) isSelectedPIDfill = fCuts->IsSelectedPID(part); //0 rejected,1 D0,2 Dobar, 3 both


  //fill "before PID" histos with every daughter
  ((TH2F*)ListDetSignal->FindObject("TOFSigBefPID"))->Fill(esdtrack1->P(), esdtrack1->GetTOFsignal());
  ((TH2F*)ListDetSignal->FindObject("TOFSigBefPID"))->Fill(esdtrack2->P(), esdtrack2->GetTOFsignal());
  ((TH2F*)ListDetSignal->FindObject("TPCSigBefPID"))->Fill(esdtrack1->P(), esdtrack1->GetTPCsignal());
  ((TH2F*)ListDetSignal->FindObject("TPCSigBefPID"))->Fill(esdtrack2->P(), esdtrack2->GetTPCsignal());

  if (isSelectedPIDfill != 0)  { //fill "After PID" with everything that's not rejected
    ((TH2F*)ListDetSignal->FindObject("TOFSigAftPID"))->Fill(esdtrack1->P(), esdtrack1->GetTOFsignal());
    ((TH2F*)ListDetSignal->FindObject("TOFSigAftPID"))->Fill(esdtrack2->P(), esdtrack2->GetTOFsignal());
    ((TH2F*)ListDetSignal->FindObject("TPCSigAftPID"))->Fill(esdtrack1->P(), esdtrack1->GetTPCsignal());
    ((TH2F*)ListDetSignal->FindObject("TPCSigAftPID"))->Fill(esdtrack2->P(), esdtrack2->GetTPCsignal());

  }

  delete esdtrack1;
  delete esdtrack2;

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0BDT::FillVarHists(AliAODEvent* aod,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout){
  //
  /// function used in UserExec to fill variable histograms:
  //


  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t lab=-9999;
  if(fReadMC) lab=part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
  //Double_t pt = d->Pt(); //mother pt
  Int_t isSelectedPID=3;
  if(!fReadMC || (fReadMC && fUsePid4Distr)) isSelectedPID=cuts->IsSelectedPID(part); //0 rejected,1 D0,2 Dobar, 3 both
  if (!fPIDCheck) {  //if fPIDCheck, this is already filled elsewhere
    if (isSelectedPID==0)fNentries->Fill(7);
    if (isSelectedPID==1)fNentries->Fill(8);
    if (isSelectedPID==2)fNentries->Fill(9);
    if (isSelectedPID==3)fNentries->Fill(10);
    //fNentries->Fill(8+isSelectedPID);
  }

  if(fCutOnDistr && !fIsSelectedCandidate) return;
  //printf("\nif no cuts or cuts passed\n");


  //add distr here
  UInt_t pdgs[2];

  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  pdgs[0]=211;
  pdgs[1]=321;
  Double_t minvD0 = part->InvMassD0();
  pdgs[1]=211;
  pdgs[0]=321;
  Double_t minvD0bar = part->InvMassD0bar();
  //cout<<"inside mass cut"<<endl;

  Double_t invmasscut=0.03;

  TString fillthispi="",fillthisK="",fillthis="", fillthispt="", fillthisetaphi="";

  Int_t ptbin=cuts->PtBin(part->Pt());
  Double_t pt = part->Pt();

  Double_t dz1[2],dz2[2],covar1[3],covar2[3];//,d0xd0proper,errd0xd0proper;
  dz1[0]=-99; dz2[0]=-99;
  Double_t d0[2];
  Double_t decl[2]={-99,-99};
  Bool_t recalcvtx=kFALSE;



  if(fCuts->GetIsPrimaryWithoutDaughters()){
    recalcvtx=kTRUE;
    //recalculate vertex
    AliAODVertex *vtxProp=0x0;
    vtxProp=GetPrimaryVtxSkipped(aod);
    if(vtxProp) {
      part->SetOwnPrimaryVtx(vtxProp);
      //Bool_t unsetvtx=kTRUE;
      //Calculate d0 for daughter tracks with Vtx Skipped
      AliESDtrack *esdtrack1=new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(0));
      AliESDtrack *esdtrack2=new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(1));
      esdtrack1->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz1,covar1);
      esdtrack2->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz2,covar2);
      delete vtxProp; vtxProp=NULL;
      delete esdtrack1;
      delete esdtrack2;
    }

    d0[0]=dz1[0];
    d0[1]=dz2[0];

    decl[0]=part->DecayLength2();
    decl[1]=part->NormalizedDecayLength2();
    part->UnsetOwnPrimaryVtx();

  }

  Double_t cosThetaStarD0 = 99;
  Double_t cosThetaStarD0bar = 99;
  Double_t cosPointingAngle = 99;
  Double_t normalizedDecayLength2 = -1, normalizedDecayLengthxy=-1;
  Double_t decayLength2 = -1, decayLengthxy=-1;
  Double_t ptProng[2]={-99,-99};
  Double_t d0Prong[2]={-99,-99};
  Double_t etaD = 99.;
  Double_t phiD = 99.;


  //disable the PID
  if(!fUsePid4Distr) isSelectedPID=0;
  if((lab>=0 && fReadMC) || (!fReadMC && isSelectedPID)){ //signal (from MC or PID)

    //check pdg of the prongs
    AliAODTrack *prong0=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
    AliAODTrack *prong1=(AliAODTrack*)fDaughterTracks.UncheckedAt(1);

    if(!prong0 || !prong1) {
      return;
    }
    Int_t labprong[2];
    if(fReadMC){
      labprong[0]=prong0->GetLabel();
      labprong[1]=prong1->GetLabel();
    }
    AliAODMCParticle *mcprong=0;
    Int_t pdgProng[2]={0,0};

    for (Int_t iprong=0;iprong<2;iprong++){
      if(fReadMC && labprong[iprong]>=0) {
	mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	pdgProng[iprong]=mcprong->GetPdgCode();
      }
    }

    if(fSys==0){
      //no mass cut ditributions: ptbis

      fillthispi="hptpiSnoMcut_";
      fillthispi+=ptbin;

      fillthisK="hptKSnoMcut_";
      fillthisK+=ptbin;

      if ((TMath::Abs(pdgProng[0]) == 211 && TMath::Abs(pdgProng[1]) == 321)
	  || (isSelectedPID==1 || isSelectedPID==3)){
	((TH1F*)listout->FindObject(fillthispi))->Fill(prong0->Pt());
	((TH1F*)listout->FindObject(fillthisK))->Fill(prong1->Pt());
      }

      if ((TMath::Abs(pdgProng[0]) == 321 && TMath::Abs(pdgProng[1]) == 211)
	  || (isSelectedPID==2 || isSelectedPID==3)){
	((TH1F*)listout->FindObject(fillthisK))->Fill(prong0->Pt());
	((TH1F*)listout->FindObject(fillthispi))->Fill(prong1->Pt());
      }
    }

    //no mass cut ditributions: mass

    etaD = part->Eta();
    phiD = part->Phi();


    fillthis="hMassS_";
    fillthis+=ptbin;
    fillthispt="histSgnPt";

    if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)
	|| (!fReadMC && (isSelectedPID==1 || isSelectedPID==3))){//D0
      ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      if(fFillPtHist && fReadMC) ((TH2F*)fOutputMassPt->FindObject(fillthispt))->Fill(minvD0,pt);

      fillthisetaphi="hetaphiD0candidateS_";
      fillthisetaphi+=ptbin;
      ((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);

      if(TMath::Abs(minvD0-mPDG)<0.05){
	fillthisetaphi="hetaphiD0candidatesignalregionS_";
	fillthisetaphi+=ptbin;
	((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);
      }

    }
    else { //D0bar
      if(fReadMC || (!fReadMC && isSelectedPID > 1)){
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
	if(fFillPtHist && fReadMC) ((TH2F*)fOutputMassPt->FindObject(fillthispt))->Fill(minvD0bar,pt);

	fillthisetaphi="hetaphiD0barcandidateS_";
	fillthisetaphi+=ptbin;
	((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);

	if(TMath::Abs(minvD0bar-mPDG)<0.05){
	  fillthisetaphi="hetaphiD0barcandidatesignalregionS_";
	  fillthisetaphi+=ptbin;
	  ((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);
	}

      }
    }

    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){

      cosThetaStarD0 = part->CosThetaStarD0();
      cosThetaStarD0bar = part->CosThetaStarD0bar();
      cosPointingAngle = part->CosPointingAngle();
      normalizedDecayLength2 = part->NormalizedDecayLength2();
      decayLength2 = part->DecayLength2();
      decayLengthxy = part->DecayLengthXY();
      normalizedDecayLengthxy=decayLengthxy/part->DecayLengthXYError();

      ptProng[0]=prong0->Pt(); ptProng[1]=prong1->Pt();
      d0Prong[0]=part->Getd0Prong(0); d0Prong[1]=part->Getd0Prong(1);

      if(fArray==1) cout<<"LS signal: ERROR"<<endl;
      for (Int_t iprong=0; iprong<2; iprong++){
	AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(iprong);
	if (fReadMC) labprong[iprong]=prong->GetLabel();

	//cout<<"prong name = "<<prong->GetName()<<" label = "<<prong->GetLabel()<<endl;
	Int_t pdgprong=0;
	if(fReadMC && labprong[iprong]>=0) {
	  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	  pdgprong=mcprong->GetPdgCode();
	}

	Bool_t isPionHere[2]={(isSelectedPID==1 || isSelectedPID==3) ? kTRUE : kFALSE,(isSelectedPID==2 || isSelectedPID==3) ? kTRUE : kFALSE};

	if(TMath::Abs(pdgprong)==211 || isPionHere[iprong]) {
	  //cout<<"pi"<<endl;

	  if(fSys==0){
	    fillthispi="hptpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(ptProng[iprong]);
	  }

	  fillthispi="hd0piS_";
	  fillthispi+=ptbin;
	  ((TH1F*)listout->FindObject(fillthispi))->Fill(d0Prong[iprong]);
	  if(recalcvtx) {

	    fillthispi="hd0vpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(d0[iprong]);
	  }

	}

	if(TMath::Abs(pdgprong)==321 || !isPionHere[iprong]) {
	  //cout<<"kappa"<<endl;
	  if(fSys==0){
	    fillthisK="hptKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(ptProng[iprong]);
	  }


	  fillthisK="hd0KS_";
	  fillthisK+=ptbin;
	  ((TH1F*)listout->FindObject(fillthisK))->Fill(d0Prong[iprong]);
	  if (recalcvtx){
	    fillthisK="hd0vKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(d0[iprong]);
	  }
	}

	if(fSys==0){
	  fillthis="hcosthpointd0S_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[iprong]);
	}
      } //end loop on prongs

      fillthis="hdcaS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

      fillthis="hcosthetapointS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle);

      fillthis="hcosthetapointxyS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngleXY());


      fillthis="hd0d0S_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());

      fillthis="hdeclS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLength2);

      fillthis="hnormdeclS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLength2);

      fillthis="hdeclxyS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLengthxy);

      fillthis="hnormdeclxyS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy);

      fillthis="hdeclxyd0d0S_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(decayLengthxy,d0Prong[0]*d0Prong[1]);

      fillthis="hnormdeclxyd0d0S_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy,d0Prong[0]*d0Prong[1]);

      if(recalcvtx) {
	fillthis="hdeclvS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[0]);

	fillthis="hnormdeclvS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[1]);

	fillthis="hd0d0vS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
      }

      if(fSys==0){
	fillthis="hcosthetastarS_";
	fillthis+=ptbin;
	if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)) ((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	else {
	  if (fReadMC || isSelectedPID>1)((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0bar);
	  if(isSelectedPID==1 || isSelectedPID==3)((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	}

	fillthis="hcosthpointd0d0S_";
	fillthis+=ptbin;
	((TH2F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,part->Prodd0d0());
      }

      if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)){
	for(Int_t it=0; it<2; it++){
	  fillthis="hptD0S_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt());
	  fillthis="hphiD0S_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Phi());
	  Int_t nPointsITS = 0;
	  for (Int_t il=0; il<6; il++){
	    if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(il)) nPointsITS++;
	  }
	  fillthis="hNITSpointsD0vsptS_";
	  fillthis+=ptbin;
	  ((TH2F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt(),nPointsITS);
	  fillthis="hNSPDpointsD0S_";
	  fillthis+=ptbin;
	  if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && !(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1))){ //no SPD points
	    ((TH1I*)listout->FindObject(fillthis))->Fill(0);
	  }
	  if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && !(((AliAODTrack*)(fDaughterTracks.UncheckedAt(it)))->HasPointOnITSLayer(1))){ //kOnlyFirst
	    ((TH1I*)listout->FindObject(fillthis))->Fill(1);
	  }
	  if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kOnlySecond
	    ((TH1I*)listout->FindObject(fillthis))->Fill(2);
	  }
	  if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kboth
	    ((TH1I*)listout->FindObject(fillthis))->Fill(3);
	  }
	  fillthis="hNclsD0vsptS_";
	  fillthis+=ptbin;
	  Float_t mom = ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt();
	  Float_t ncls = (Float_t)((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->GetTPCNcls();
	  ((TH2F*)listout->FindObject(fillthis))->Fill(mom, ncls);
	}
      }
      else {
	if (fReadMC || isSelectedPID>1){
	  for(Int_t it=0; it<2; it++){
	    fillthis="hptD0barS_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt());
	    fillthis="hphiD0barS_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Phi());
	    fillthis="hNclsD0barvsptS_";
	    fillthis+=ptbin;
	    Float_t mom = ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt();
	    Float_t ncls = (Float_t)((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->GetTPCNcls();
	    ((TH2F*)listout->FindObject(fillthis))->Fill(mom, ncls);
	  }
	}
	if(isSelectedPID==1 || isSelectedPID==3){
	  for(Int_t it=0; it<2; it++){
	    fillthis="hptD0S_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt());
	    fillthis="hphiD0S_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Phi());
	    Int_t nPointsITS = 0;
	    for (Int_t il=0; il<6; il++){
	      if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(il)) nPointsITS++;
	    }
	    fillthis="hNITSpointsD0vsptS_";
	    fillthis+=ptbin;
	    ((TH2F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt(), nPointsITS);
	    fillthis="hNSPDpointsD0S_";
	    fillthis+=ptbin;
	    if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && !(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1))){ //no SPD points
	      ((TH1I*)listout->FindObject(fillthis))->Fill(0);
	    }
	    if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && !(((AliAODTrack*)(fDaughterTracks.UncheckedAt(it)))->HasPointOnITSLayer(1))){ //kOnlyFirst
	      ((TH1I*)listout->FindObject(fillthis))->Fill(1);
	    }
	    if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kOnlySecond
	      ((TH1I*)listout->FindObject(fillthis))->Fill(2);
	    }
	    if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kboth
	      ((TH1I*)listout->FindObject(fillthis))->Fill(3);
	    }
       	    fillthis="hNclsD0vsptS_";
	    fillthis+=ptbin;
	    Float_t mom = ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt();
	    Float_t ncls = (Float_t)((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->GetTPCNcls();
	    ((TH2F*)listout->FindObject(fillthis))->Fill(mom, ncls);
	  }
	}
      }


    } //end mass cut

  } else{ //Background or LS
    //if(!fReadMC){
    //cout<<"is background"<<endl;

    etaD = part->Eta();
    phiD = part->Phi();

    //no mass cut distributions: mass, ptbis
    fillthis="hMassB_";
    fillthis+=ptbin;
    fillthispt="histBkgPt";

    if (!fCutOnDistr || (fCutOnDistr && (fIsSelectedCandidate==1 || fIsSelectedCandidate==3))) {
      ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      if(fFillPtHist && fReadMC) ((TH2F*)fOutputMassPt->FindObject(fillthispt))->Fill(minvD0,pt);

      fillthisetaphi="hetaphiD0candidateB_";
      fillthisetaphi+=ptbin;
      ((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);

      if(TMath::Abs(minvD0-mPDG)<0.05){
	fillthisetaphi="hetaphiD0candidatesignalregionB_";
	fillthisetaphi+=ptbin;
	((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);
      }
    }
    if (!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate>1)) {
      ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
      if(fFillPtHist && fReadMC) ((TH2F*)fOutputMassPt->FindObject(fillthispt))->Fill(minvD0bar,pt);

      fillthisetaphi="hetaphiD0barcandidateB_";
      fillthisetaphi+=ptbin;
      ((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);

      if(TMath::Abs(minvD0bar-mPDG)<0.05){
	fillthisetaphi="hetaphiD0barcandidatesignalregionB_";
	fillthisetaphi+=ptbin;
	((TH2F*)listout->FindObject(fillthisetaphi))->Fill(etaD, phiD);
      }

    }
    if(fSys==0){
      fillthis="hptB1prongnoMcut_";
      fillthis+=ptbin;

      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());

      fillthis="hptB2prongsnoMcut_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());
      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());
    }


    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){
      if(fSys==0){
	ptProng[0]=((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt(); ptProng[1]=((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt();
	cosThetaStarD0 = part->CosThetaStarD0();
	cosThetaStarD0bar = part->CosThetaStarD0bar();
      }

      cosPointingAngle = part->CosPointingAngle();
      normalizedDecayLength2 = part->NormalizedDecayLength2();
      decayLength2 = part->DecayLength2();
      decayLengthxy = part->DecayLengthXY();
      normalizedDecayLengthxy=decayLengthxy/part->DecayLengthXYError();
      d0Prong[0]=part->Getd0Prong(0); d0Prong[1]=part->Getd0Prong(1);


      AliAODTrack *prongg=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
      if(!prongg) {
	if(fDebug>2) cout<<"No daughter found";
	return;
      }
      else{
	if(fArray==1){
	  if(prongg->Charge()==1) {
	    //fTotPosPairs[ptbin]++;
	    ((TH1F*)fOutputMass->FindObject("hpospair"))->Fill(ptbin);
	  } else {
	    //fTotNegPairs[ptbin]++;
	    ((TH1F*)fOutputMass->FindObject("hnegpair"))->Fill(ptbin);
	  }
	}
      }

      //fill pt and phi distrib for prongs with M cut

      if (!fCutOnDistr || (fCutOnDistr && (fIsSelectedCandidate==1 || fIsSelectedCandidate==3))){
  	for(Int_t it=0; it<2; it++){
  	  fillthis="hptD0B_";
 	  fillthis+=ptbin;
 	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt());
 	  fillthis="hphiD0B_";
 	  fillthis+=ptbin;
 	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Phi());

 	  Int_t nPointsITS = 0;
 	  for (Int_t il=0; il<6; il++){
 	    if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(il)) nPointsITS++;
 	  }
 	  fillthis="hNITSpointsD0vsptB_";
 	  fillthis+=ptbin;
	  ((TH2F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt(), nPointsITS);
	  fillthis="hNSPDpointsD0B_";
	  fillthis+=ptbin;
	  if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && !(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1))){ //no SPD points
	    ((TH1I*)listout->FindObject(fillthis))->Fill(0);
	  }
	  if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && !(((AliAODTrack*)(fDaughterTracks.UncheckedAt(it)))->HasPointOnITSLayer(1))){ //kOnlyFirst
	    ((TH1I*)listout->FindObject(fillthis))->Fill(1);
	  }
	  if(!(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0)) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kOnlySecond
	    ((TH1I*)listout->FindObject(fillthis))->Fill(2);
	  }
	  if(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(0) && ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->HasPointOnITSLayer(1)){ //kboth
	    ((TH1I*)listout->FindObject(fillthis))->Fill(3);
	  }
	  fillthis="hNclsD0vsptB_";
	  fillthis+=ptbin;
	  Float_t mom = ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt();
	  Float_t ncls = (Float_t)((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->GetTPCNcls();
	  ((TH2F*)listout->FindObject(fillthis))->Fill(mom, ncls);
	}


      }

      if (!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate>1)) {
 	for(Int_t it=0; it<2; it++){
	  fillthis="hptD0barB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt());
	  fillthis="hphiD0barB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Phi());
	  fillthis="hNclsD0barvsptB_";
	  fillthis+=ptbin;
	  Float_t mom = ((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->Pt();
	  Float_t ncls = (Float_t)((AliAODTrack*)fDaughterTracks.UncheckedAt(it))->GetTPCNcls();
	  ((TH2F*)listout->FindObject(fillthis))->Fill(mom, ncls);
	}
      }

      fillthis="hd0B_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]);
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[1]);

      if(fReadMC){
	Int_t pdgMother[2]={0,0};
	Double_t factor[2]={1,1};

	for(Int_t iprong=0;iprong<2;iprong++){
	  AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(iprong);
	  lab=prong->GetLabel();
	  if(lab>=0){
	    AliAODMCParticle* mcprong=(AliAODMCParticle*)arrMC->At(lab);
	    if(mcprong){
	      Int_t labmom=mcprong->GetMother();
	      if(labmom>=0){
		AliAODMCParticle* mcmother=(AliAODMCParticle*)arrMC->At(labmom);
		if(mcmother) pdgMother[iprong]=mcmother->GetPdgCode();
	      }
	    }
	  }

	  if(fSys==0){

	    fillthis="hd0moresB_";
	    fillthis+=ptbin;

	    if(TMath::Abs(pdgMother[iprong])==310 || TMath::Abs(pdgMother[iprong])==130 || TMath::Abs(pdgMother[iprong])==321){ //K^0_S, K^0_L, K^+-
	      if(ptProng[iprong]<=1)factor[iprong]=1./.7;
	      else factor[iprong]=1./.6;
	      fNentries->Fill(11);
	    }

	    if(TMath::Abs(pdgMother[iprong])==3122) { //Lambda
	      factor[iprong]=1./0.25;
	      fNentries->Fill(12);
	    }
	    fillthis="hd0moresB_";
	    fillthis+=ptbin;

	    ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[iprong],factor[iprong]);

	    if(recalcvtx){
	      fillthis="hd0vmoresB_";
	      fillthis+=ptbin;
	      ((TH1F*)listout->FindObject(fillthis))->Fill(d0[iprong],factor[iprong]);
	    }
	  }
	} //loop on prongs

	if(fSys==0){
	  fillthis="hd0d0moresB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0(),factor[0]*factor[1]);

	  fillthis="hcosthetapointmoresB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,factor[0]*factor[1]);

	  if(recalcvtx){
	    fillthis="hd0d0vmoresB_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1],factor[0]*factor[1]);
	  }
	}
      } //readMC

      if(fSys==0){
	//normalise pt distr to half afterwards
	fillthis="hptB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(ptProng[0]);
	((TH1F*)listout->FindObject(fillthis))->Fill(ptProng[1]);

	fillthis="hcosthetastarB_";
	fillthis+=ptbin;
	if (!fCutOnDistr || (fCutOnDistr && (fIsSelectedCandidate==1 || fIsSelectedCandidate==3)))((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	if (!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate>1))((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0bar);


	fillthis="hd0p0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]);
	fillthis="hd0p1B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[1]);

	fillthis="hcosthpointd0d0B_";
	fillthis+=ptbin;
	((TH2F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,part->Prodd0d0());

	fillthis="hcosthpointd0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[0]);
	((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[1]);


	if(recalcvtx){

	  fillthis="hd0vp0B_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	  fillthis="hd0vp1B_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);

	  fillthis="hd0vB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);

	}

      }

      fillthis="hdcaB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

      fillthis="hd0d0B_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]*d0Prong[1]);

      if(recalcvtx){
	fillthis="hd0d0vB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
      }

      fillthis="hcosthetapointB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle);

      fillthis="hcosthetapointxyB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngleXY());

      fillthis="hdeclB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLength2);

      fillthis="hnormdeclB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLength2);

      fillthis="hdeclxyB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLengthxy);

      fillthis="hnormdeclxyB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy);

      fillthis="hdeclxyd0d0B_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(decayLengthxy,d0Prong[0]*d0Prong[1]);

      fillthis="hnormdeclxyd0d0B_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy,d0Prong[0]*d0Prong[1]);


      if(recalcvtx) {

	fillthis="hdeclvB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[0]);

	fillthis="hnormdeclvB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[1]);


      }
    }//mass cut
  }//else (background)

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0BDT::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  /// function used in UserExec to fill mass histograms:
  //


  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  //cout<<"is selected = "<<fIsSelectedCandidate<<endl;

  //cout<<"check cuts = "<<endl;
  //cuts->PrintAll();
  if (!fIsSelectedCandidate){
    //cout<<" cut " << cuts << " Rejected because "<<cuts->GetWhy()<<endl;
    return;
  }


  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  //printf("SELECTED\n");
  Int_t ptbin=cuts->PtBin(part->Pt());
  Double_t pt = part->Pt();
  Double_t y = part->YD0();

  Double_t impparXY=part->ImpParXY()*10000.;
  Double_t trueImpParXY=0.;
  Double_t arrayForSparse[3]={invmassD0,pt,impparXY};
  Double_t arrayForSparseTrue[3]={invmassD0,pt,trueImpParXY};


  // AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
  // if(!prong) {
  //   AliDebug(2,"No daughter found");
  //   return;
  // }
  // else{
  // if(prong->Charge()==1) {
  //   ((TH1F*)fDistr->FindObject("hpospair"))->Fill(fCuts->GetNPtBins()+ptbin);
  //   //fTotPosPairs[ptbin]++;
  // } else {
  //   ((TH1F*)fDistr->FindObject("hnegpair"))->Fill(fCuts->GetNPtBins()+ptbin);
  //   //fTotNegPairs[ptbin]++;
  // }
  //  }

  // for(Int_t it=0;it<2;it++){

  //    //request on spd points to be addes
  //   if(/*nSPD==2 && */part->Pt() > 5. && (TMath::Abs(invmassD0-mPDG)<0.01 || TMath::Abs(invmassD0bar-mPDG)<0.01)){
  //     FILE *f=fopen("4display.txt","a");
  //     fprintf(f,"pt: %f \n Rapidity: %f \t Period Number: %x \t Run Number: %d \t BunchCrossNumb: %d \t OrbitNumber: %d\n",part->Pt(),part->Y(421),aod->GetPeriodNumber(),aod->GetRunNumber(),aod->GetBunchCrossNumber(),aod->GetOrbitNumber());
  //     fclose(f);
  //     //printf("PrimVtx NContributors: %d \n Prongs Rel Angle: %f \n \n",ncont,relangle);
  //   }
  // }

  TString fillthis="", fillthispt="", fillthismasspt="", fillthismassy="", fillthissub="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  Bool_t isPrimary=kTRUE;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //Define weights for Bayesian (if appropriate)

  Double_t weigD0=1.;
  Double_t weigD0bar=1.;
  if (fCuts->GetCombPID() && (fCuts->GetBayesianStrategy() == AliRDHFCutsD0toKpi::kBayesWeight || fCuts->GetBayesianStrategy() == AliRDHFCutsD0toKpi::kBayesWeightNoFilter)) {
    weigD0=fCuts->GetWeightsNegative()[AliPID::kKaon] * fCuts->GetWeightsPositive()[AliPID::kPion];
    weigD0bar=fCuts->GetWeightsPositive()[AliPID::kKaon] * fCuts->GetWeightsNegative()[AliPID::kPion];
    if (weigD0 < 0. || weigD0 > 1.0) {weigD0 = 0.;}
    if (weigD0bar < 0. || weigD0bar > 1.0) {weigD0bar = 0.;} //Prevents filling with weight > 1, or < 0
  }

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);

  if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2) { //D0

    arrayForSparse[0]=invmassD0; arrayForSparse[2]=impparXY;

    if(fReadMC){
      if(labD0>=0) {
	if(fArray==1) cout<<"LS signal ERROR"<<endl;

	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//	cout<<"pdg = "<<pdgD0<<endl;

	// old function	if(CheckOrigin(arrMC,partD0)==5) isPrimary=kFALSE;
        if(AliVertexingHFUtils::CheckOrigin(arrMC,partD0,fUseQuarkTagInKine)==5) isPrimary=kFALSE;
	if(!isPrimary)
	  trueImpParXY=GetTrueImpactParameter(mcHeader,arrMC,partD0)*10000.;
	arrayForSparseTrue[0]=invmassD0; arrayForSparseTrue[2]=trueImpParXY;

	if (pdgD0==421){ //D0
	  //	  cout<<"Fill S with D0"<<endl;
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,weigD0);

	  if(fFillPtHist){
	    fillthismasspt="histSgnPt";
	    ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0,pt,weigD0);
	  }
	  if(fFillImpParHist){
	    if(isPrimary) fHistMassPtImpParTC[1]->Fill(arrayForSparse,weigD0);
	    else {
	      fHistMassPtImpParTC[2]->Fill(arrayForSparse,weigD0);
	      fHistMassPtImpParTC[3]->Fill(arrayForSparseTrue,weigD0);
	    }
	  }

	  if(fFillYHist){
	    fillthismassy="histSgnY_";
	    fillthismassy+=ptbin;
	    ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0,y,weigD0);
	  }

	  if(fSys==0){
	    if(TMath::Abs(invmassD0 - mPDG) < 0.027 && fFillVarHists){
	      fillthis="histSgn27_";
	      fillthis+=ptbin;
	      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,weigD0);
	    }
	  }
	} else{ //it was a D0bar
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,weigD0);

	  if(fFillPtHist){
	    fillthismasspt="histRflPt";
	    //	    cout << " Filling "<<fillthismasspt<<" D0bar"<<endl;
	    ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0,pt,weigD0);
	  }

	  if(fFillYHist){
	    fillthismassy="histRflY_";
	    fillthismassy+=ptbin;
	    //	    cout << " Filling "<<fillthismassy<<" D0bar"<<endl;
	    ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0,y,weigD0);
	  }

	}
      } else {//background
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,weigD0);

	if(fFillPtHist){
	  fillthismasspt="histBkgPt";
	  //	  cout << " Filling "<<fillthismasspt<<" D0bar"<<endl;
	  ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0,pt,weigD0);
	}
	if(fFillImpParHist) fHistMassPtImpParTC[4]->Fill(arrayForSparse,weigD0);

	if(fFillYHist){
	  fillthismassy="histBkgY_";
	  fillthismassy+=ptbin;
	  //	  cout << " Filling "<<fillthismassy<<" D0bar"<<endl;
	  ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0,y,weigD0);
	}

      }

    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      //      cout<<"Filling "<<fillthis<<endl;

      //      printf("Fill mass with D0");
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0,weigD0);

      if(fFillSubSampleHist){
        fillthissub="histMassvsSubSample_";
        fillthissub+=ptbin;
        ((TH2F*)(listout->FindObject(fillthissub)))->Fill(invmassD0,(Double_t)(fEventCounter%24),weigD0);
      }

      if(fFillPtHist){
	fillthismasspt="histMassPt";
	//	cout<<"Filling "<<fillthismasspt<<endl;
	((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0,pt,weigD0);
      }
      if(fFillImpParHist) {
	//	cout << "Filling fHistMassPtImpParTC[0]"<<endl;
	fHistMassPtImpParTC[0]->Fill(arrayForSparse,weigD0);
      }

      if(fFillYHist){
	fillthismassy="histMassY_";
	fillthismassy+=ptbin;
	//	cout<<"Filling "<<fillthismassy<<endl;
	((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0,y,weigD0);
      }

    }

  }
  if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) { //D0bar

    arrayForSparse[0]=invmassD0bar; arrayForSparse[2]=impparXY;

    if(fReadMC){
      if(labD0>=0) {
	if(fArray==1) cout<<"LS signal ERROR"<<endl;
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//	cout<<" pdg = "<<pdgD0<<endl;

        //old function	if(CheckOrigin(arrMC,partD0)==5) isPrimary=kFALSE;
        if(AliVertexingHFUtils::CheckOrigin(arrMC,partD0,fUseQuarkTagInKine)==5) isPrimary=kFALSE;
	if(!isPrimary)
	  trueImpParXY=GetTrueImpactParameter(mcHeader,arrMC,partD0)*10000.;
	arrayForSparseTrue[0]=invmassD0bar; arrayForSparseTrue[2]=trueImpParXY;

	if (pdgD0==-421){ //D0bar
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,weigD0bar);
	  // if (TMath::Abs(invmassD0bar - mPDG) < 0.027){
	  //   fillthis="histSgn27_";
	  //   fillthis+=ptbin;
	  //   ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	  // }

	  if(fFillPtHist){
	    fillthismasspt="histSgnPt";
	    //	    cout<<" Filling "<< fillthismasspt << endl;
	    ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0bar,pt,weigD0bar);
	  }
	  if(fFillImpParHist){
	    //	    cout << " Filling impact parameter thnsparse"<<endl;
	    if(isPrimary) fHistMassPtImpParTC[1]->Fill(arrayForSparse,weigD0bar);
	    else {
	      fHistMassPtImpParTC[2]->Fill(arrayForSparse,weigD0bar);
	      fHistMassPtImpParTC[3]->Fill(arrayForSparseTrue,weigD0bar);
	    }
	  }

	  if(fFillYHist){
	    fillthismassy="histSgnY_";
	    fillthismassy+=ptbin;
	    //	    cout<<" Filling "<< fillthismassy << endl;
	    ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0bar,y,weigD0bar);
	  }

	} else{
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,weigD0bar);
	  if(fFillPtHist){
	    fillthismasspt="histRflPt";
	    //	    cout << " Filling "<<fillthismasspt<<endl;
	    ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0bar,pt,weigD0bar);
	  }
	  if(fFillYHist){
	    fillthismassy="histRflY_";
	    fillthismassy+=ptbin;
	    //	    cout << " Filling "<<fillthismassy<<endl;
	    ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0bar,y,weigD0bar);
	  }
	}
      } else {//background or LS
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar,weigD0bar);

	if(fFillPtHist){
	  fillthismasspt="histBkgPt";
	  //	  cout<<" Filling "<< fillthismasspt << endl;
	  ((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0bar,pt,weigD0bar);
	}
	if(fFillImpParHist) fHistMassPtImpParTC[4]->Fill(arrayForSparse,weigD0bar);
	if(fFillYHist){
	  fillthismassy="histBkgY_";
	  fillthismassy+=ptbin;
	  //	  cout<<" Filling "<< fillthismassy << endl;
	  ((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0bar,y,weigD0bar);
	}
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      //      printf("Fill mass with D0bar");

      ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar,weigD0bar);

      if(fFillSubSampleHist){
        fillthissub="histMassvsSubSample_";
        fillthissub+=ptbin;
        ((TH2F*)(listout->FindObject(fillthissub)))->Fill(invmassD0bar,(Double_t)(fEventCounter%24),weigD0);
      }

      if(fFillPtHist){
	fillthismasspt="histMassPt";
	//	cout<<" Filling "<< fillthismasspt << endl;
	((TH2F*)(fOutputMassPt->FindObject(fillthismasspt)))->Fill(invmassD0bar,pt,weigD0bar);
      }
      if(fFillImpParHist) fHistMassPtImpParTC[0]->Fill(arrayForSparse,weigD0bar);
      if(fFillYHist){
	fillthismassy="histMassY_";
	fillthismassy+=ptbin;
	//	cout<<" Filling "<< fillthismassy << endl;
	((TH2F*)(fOutputMassY->FindObject(fillthismassy)))->Fill(invmassD0bar,y,weigD0bar);
      }
    }
  }

  return;
}

//__________________________________________________________________________
void AliAnalysisTaskSED0BDT::FillCandVariables(AliAODEvent *aodev, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi *cuts){
  // Fill candidate variables for cut study

  if(fWriteVariableTree && fWriteProtosgnVar){
    AliFatal("FATAL_ERROR: Writing candidate variables both with and without --> CHOOSE ONE OF THE TWO OPTIONS!");
  }

  // Fill candidate variable Tree (track selection, no candidate selection)
  if( fWriteVariableTree && !part->HasBadDaughters()
   && fCuts->AreDaughtersSelected(part) && fCuts->IsSelectedPID(part) ){ //for backward compatibility
    fCandidateVariables[0] = part->InvMassD0();
    fCandidateVariables[1] = part->InvMassD0bar();
    fCandidateVariables[2] = part->Pt();
    fCandidateVariables[3] = part->GetDCA();
    Double_t ctsD0=0. ,ctsD0bar=0.; part->CosThetaStarD0(ctsD0,ctsD0bar);
    fCandidateVariables[4] = ctsD0;
    fCandidateVariables[5] = ctsD0bar;
    fCandidateVariables[6] = part->Pt2Prong(0);
    fCandidateVariables[7] = part->Pt2Prong(1);
    fCandidateVariables[8] = part->Getd0Prong(0);
    fCandidateVariables[9] = part->Getd0Prong(1);
    fCandidateVariables[10] = part->Prodd0d0();
    fCandidateVariables[11] = part->CosPointingAngle();
    fCandidateVariables[12] = part->CosPointingAngleXY();
    fCandidateVariables[13] = part->NormalizedDecayLengthXY();
    fCandidateVariables[14] = fCuts->IsSelectedSpecialCuts(part);
    fVariablesTree->Fill();
  }

  // Fill candidate variable Tree (with candidate selection)
  Double_t mPDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  Bool_t isPrimary=kTRUE;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  if (!fIsSelectedCandidate) {
    return;
  }

  if (fDebug > 2) cout << "Candidate selected" << endl;

  Float_t invmassD0 = part->InvMassD0();
  Float_t invmassD0bar = part->InvMassD0bar();
  Bool_t isD0sel = false;
  Bool_t isD0barsel = false;
  Bool_t isrealD0 = false;
  Bool_t isrealD0bar = false;
  Int_t selCand = 999.; // flag to store how the candidate is selected (0 = D0 only; 1 = D0bar only; 2 = D0 and D0bar)

  if(fWriteProtosgnVar){ //write candidate variables for proto-significance study
    if ((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) && fFillOnlyD0D0bar < 2){ //D0
      isD0sel = true;
      if (fReadMC) {
        if (labD0 >= 0) {
          if (fArray == 1)  cout << "LS signal ERROR" << endl;
          AliAODMCParticle *partD0 = (AliAODMCParticle *)arrMC->At(labD0);
          Int_t pdgD0 = partD0->GetPdgCode();
          if (AliVertexingHFUtils::CheckOrigin(arrMC, partD0, fUseQuarkTagInKine) == 5) isPrimary = kFALSE;
          if (pdgD0 == 421) { // D0
            if (fDebug > 2)  cout << "MC: D0 candidate" << endl;
            isrealD0 = true;
            if (fUsedMassWindow && (part->InvMassD0() < 1.7 || part->InvMassD0() > 2.1)) invmassD0 = 0;
          } else { // it was a D0bar
            if (fDebug > 2)  cout << "MC: D0bar candidate selected as D0 --> reflection" << endl;
          }
        } else { // background
          if (fDebug > 2)  cout << "Combinatorial background" << endl;
        }
      } else {
        if (fUsedMassWindow && (part->InvMassD0() < 1.7 || part->InvMassD0() > 2.1)) invmassD0 = 0;
      }
    }

    if(fIsSelectedCandidate > 1 && (fFillOnlyD0D0bar == 0 || fFillOnlyD0D0bar == 2)){ //D0bar
      isD0barsel=true;
      if (fReadMC) {
        if (labD0 >= 0){
          if (fArray == 1) cout << "LS signal ERROR" << endl;
          AliAODMCParticle *partD0 = (AliAODMCParticle *)arrMC->At(labD0);
          Int_t pdgD0 = partD0->GetPdgCode();
          if (AliVertexingHFUtils::CheckOrigin(arrMC, partD0,fUseQuarkTagInKine) == 5) isPrimary = kFALSE;
          if (pdgD0 == -421) { // D0bar
            if (fDebug > 2)  cout << "MC: D0bar candidate" << endl;
            isrealD0bar = true;
            if (fUsedMassWindow && (part->InvMassD0bar() < 1.7 || part->InvMassD0bar() > 2.1)) invmassD0bar = 0;
          } else {
            if (fDebug > 2)  cout << "MC: D0 candidate selected as D0bar --> reflection" << endl;
          }
        } else {
          // background or LS
          if (fDebug > 2)  cout << "Combinatorial background" << endl;
        }
      } else {
        if (fUsedMassWindow && (part->InvMassD0bar() < 1.7 || part->InvMassD0bar() > 2.1)) invmassD0bar = 0;
      }
    }

    if(!isD0sel && !isD0barsel) return;
    //assignment candidate flag: 0-->D0; 1-->D0bar; 2-->D0 and D0bar
    if(!fUsedMassWindow) AliWarning("WARNING: Mass window selection NOT used!");

    if(fReadMC && fSelectTrueD0){
      if(!isrealD0 && !isrealD0bar) return;
      if(isD0sel && isrealD0) selCand = 0;
      if(isD0barsel && isrealD0bar) selCand = 1;
    } else {
      if(isD0sel && !isD0barsel){
        selCand = 0;
      }else if(!isD0sel && isD0barsel){
        selCand = 1;
      }else if(isD0sel && isD0barsel){
        selCand = 2;
      }
    }

    fCandidateVariables[0] = invmassD0;
    fCandidateVariables[1] = invmassD0bar;
    fCandidateVariables[2] = part->Pt();
    fCandidateVariables[3] = part->GetDCA();
    Double_t ctsD0 = 0., ctsD0bar = 0.;
    part->CosThetaStarD0(ctsD0, ctsD0bar);
    fCandidateVariables[4] = ctsD0;
    fCandidateVariables[5] = ctsD0bar;
    fCandidateVariables[6] = part->Pt2Prong(0);
    fCandidateVariables[7] = part->Pt2Prong(1);
    fCandidateVariables[8] = part->Getd0Prong(0);
    fCandidateVariables[9] = part->Getd0Prong(1);
    fCandidateVariables[10] = part->CosPointingAngle();
    fCandidateVariables[11] = part->CosPointingAngleXY();
    fCandidateVariables[12] = part->NormalizedDecayLengthXY();
    fCandidateVariables[13] = fCuts->IsSelectedSpecialCuts(part);
    fCandidateVariables[14] = ComputeTopomatic(aodev, part);
    fCandidateVariables[15] = selCand;
    if(fReadMC){
      if(isPrimary) fCandidateVariables[16] = 4;
      else fCandidateVariables[16] = 5;
    }
    fVariablesTree->Fill();

  }
}

//__________________________________________________________________________
AliAODVertex* AliAnalysisTaskSED0BDT::GetPrimaryVtxSkipped(AliAODEvent *aodev){
  /// Calculate the primary vertex w/o the daughter tracks of the candidate

  Int_t skipped[2];
  Int_t nTrksToSkip=2;
  AliAODTrack *dgTrack = (AliAODTrack*)fDaughterTracks.UncheckedAt(0);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[0]=dgTrack->GetID();
  dgTrack = (AliAODTrack*)fDaughterTracks.UncheckedAt(1);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[1]=dgTrack->GetID();

  AliESDVertex *vertexESD=0x0;
  AliAODVertex *vertexAOD=0x0;
  AliVertexerTracks *vertexer = new AliVertexerTracks(aodev->GetMagneticField());

  //
  vertexer->SetSkipTracks(nTrksToSkip,skipped);
  vertexer->SetMinClusters(4);
  vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(aodev);
  if(!vertexESD) return vertexAOD;
  if(vertexESD->GetNContributors()<=0) {
    AliDebug(2,"vertexing failed");
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }

  delete vertexer; vertexer=NULL;


  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;

  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  return vertexAOD;

}


//________________________________________________________________________
void AliAnalysisTaskSED0BDT::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");


  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {
    printf("ERROR: fOutputMass not available\n");
    return;
  }
  fOutputMassPt = dynamic_cast<TList*> (GetOutputData(6));
  if ((fFillPtHist || fFillImpParHist) && !fOutputMassPt) {
    printf("ERROR: fOutputMass not available\n");
    return;
  }

  if(fFillVarHists){
    fDistr = dynamic_cast<TList*> (GetOutputData(2));
    if (!fDistr) {
      printf("ERROR: fDistr not available\n");
      return;
    }
  }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(3));

  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }
  fCuts = dynamic_cast<AliRDHFCutsD0toKpi*>(GetOutputData(4));
  if(!fCuts){
    printf("ERROR: fCuts not available\n");
    return;
  }
  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(5));
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }
  if (fDrawDetSignal) {
    fDetSignal = dynamic_cast<TList*>(GetOutputData(8));
    if (!fDetSignal) {
      printf("ERROR: fDetSignal not available\n");
      return;
    }
  }
  if(fFillYHist){
    fOutputMassY = dynamic_cast<TList*> (GetOutputData(9));
    if (fFillYHist && !fOutputMassY) {
      printf("ERROR: fOutputMassY not available\n");
      return;
    }
  }

  Int_t nptbins=fCuts->GetNPtBins();
  for(Int_t ipt=0;ipt<nptbins;ipt++){

    if(fArray==1 && fFillVarHists){
      fLsNormalization = 2.*TMath::Sqrt(((TH1F*)fOutputMass->FindObject("hpospair"))->Integral(nptbins+ipt+1,nptbins+ipt+2)*((TH1F*)fOutputMass->FindObject("hnegpair"))->Integral(nptbins+ipt+1,nptbins+ipt+2)); //after cuts


      if(fLsNormalization>1e-6) {

	TString massName="histMass_";
	massName+=ipt;
	((TH1F*)fOutputMass->FindObject(massName))->Scale((1/fLsNormalization)*((TH1F*)fOutputMass->FindObject(massName))->GetEntries());

      }


      fLsNormalization = 2.*TMath::Sqrt(((TH1F*)fOutputMass->FindObject("hpospair"))->Integral(ipt+1,ipt+2)*((TH1F*)fOutputMass->FindObject("hnegpair"))->Integral(ipt+1,ipt+2));
      //fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[4]*fTotNegPairs[4]);

      if(fLsNormalization>1e-6) {

	TString nameDistr="hdcaB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0d0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetapointB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	if(fSys==0){
	  nameDistr="hptB_";
	  nameDistr+=ipt;
	  ((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	  nameDistr="hcosthetastarB_";
	  nameDistr+=ipt;
	  ((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	  nameDistr="hcosthpointd0d0B_";
	  nameDistr+=ipt;
	  ((TH2F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH2F*)fDistr->FindObject(nameDistr))->GetEntries());
	}
      }
    }
  }
  TString cvname,cstname;

  if (fArray==0){
    cvname="D0invmass";
    cstname="cstat0";
  } else {
    cvname="LSinvmass";
    cstname="cstat1";
  }

  TCanvas *cMass=new TCanvas(cvname,cvname);
  cMass->cd();
  ((TH1F*)fOutputMass->FindObject("histMass_3"))->Draw();

  TCanvas* cStat=new TCanvas(cstname,Form("Stat%s",fArray ? "LS" : "D0"));
  cStat->cd();
  cStat->SetGridy();
  fNentries->Draw("htext0");

  // TCanvas *ccheck=new TCanvas(Form("cc%d",fArray),Form("cc%d",fArray));
  // ccheck->cd();

  return;
}


//________________________________________________________________________
void AliAnalysisTaskSED0BDT::CreateImpactParameterHistos(){
  /// Histos for impact paramter study

  Int_t nmassbins=200;
  Double_t fLowmasslimit=1.5648, fUpmasslimit=2.1648;
  Int_t fNImpParBins=400;
  Double_t fLowerImpPar=-2000., fHigherImpPar=2000.;
  Int_t nbins[3]={nmassbins,200,fNImpParBins};
  Double_t xmin[3]={fLowmasslimit,0.,fLowerImpPar};
  Double_t xmax[3]={fUpmasslimit,20.,fHigherImpPar};


  fHistMassPtImpParTC[0]=new THnSparseF("hMassPtImpParAll",
					"Mass vs. pt vs.imppar - All",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[1]=new THnSparseF("hMassPtImpParPrompt",
					"Mass vs. pt vs.imppar - promptD",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[2]=new THnSparseF("hMassPtImpParBfeed",
					"Mass vs. pt vs.imppar - DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[3]=new THnSparseF("hMassPtImpParTrueBfeed",
					"Mass vs. pt vs.true imppar -DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[4]=new THnSparseF("hMassPtImpParBkg",
				        "Mass vs. pt vs.imppar - backgr.",
					3,nbins,xmin,xmax);

  for(Int_t i=0; i<5;i++){
    fOutputMassPt->Add(fHistMassPtImpParTC[i]);
  }
}

//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSED0BDT::GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD0) const {
  /// true impact parameter calculation

  //printf(" AliAnalysisTaskSED0BDTV1::GetTrueImpactParameter() \n");

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partD0->XvYvZv(origD);
  Short_t charge=partD0->Charge();
  Double_t pXdauTrue[2],pYdauTrue[2],pZdauTrue[2];
  for(Int_t iDau=0; iDau<2; iDau++){
    pXdauTrue[iDau]=0.;
    pYdauTrue[iDau]=0.;
    pZdauTrue[iDau]=0.;
  }

  //  Int_t nDau=partD0->GetNDaughters();
  Int_t labelFirstDau = partD0->GetDaughterLabel(0);

  for(Int_t iDau=0; iDau<2; iDau++){
    Int_t ind = labelFirstDau+iDau;
    AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
    if(!part) continue;
    Int_t pdgCode = TMath::Abs( part->GetPdgCode() );
    if(!part){
      AliError("Daughter particle not found in MC array");
      return 99999.;
    }
    if(pdgCode==211 || pdgCode==321){
      pXdauTrue[iDau]=part->Px();
      pYdauTrue[iDau]=part->Py();
      pZdauTrue[iDau]=part->Pz();
    }
  }

  Double_t d0dummy[2]={0.,0.};
  AliAODRecoDecayHF aodDzeroMC(vtxTrue,origD,2,charge,pXdauTrue,pYdauTrue,pZdauTrue,d0dummy);
  return aodDzeroMC.ImpParXY();

}

//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSED0BDT::ComputeTopomatic(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part) {
  Float_t dd0max = 0.;
  for (Int_t ipr = 0; ipr < 2; ipr++) {
    Double_t diffIP, errdiffIP;
    part->Getd0MeasMinusExpProng(ipr, aod->GetMagneticField(), diffIP,
                                 errdiffIP);
    Double_t normdd0 = 0.;
    if (errdiffIP > 0.)
      normdd0 = diffIP / errdiffIP;
    if (ipr == 0)
      dd0max = normdd0;
    else if (TMath::Abs(normdd0) > TMath::Abs(dd0max))
      dd0max = normdd0;
  }
  return TMath::Abs(dd0max);
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0BDT::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {
  //  obsolete method
  /// checking whether the mother of the particles come from a charm or a bottom quark
  //
  printf(" AliAnalysisTaskSED0BDT V1::CheckOrigin() \n");

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  if(isFromB) return 5;
  else return 4;
}
//_______________________________________
void AliAnalysisTaskSED0BDT::CreateMCAcceptanceHistos(){
  /// Histos for MC Acceptance histos

  const Int_t nVarPrompt = 2;
  const Int_t nVarFD = 3;

  Int_t nbinsPrompt[nVarPrompt]={200,100};
  Int_t nbinsFD[nVarFD]={200,100,200};

  Double_t xminPrompt[nVarPrompt] = {0.,-1.};
  Double_t xmaxPrompt[nVarPrompt] = {40.,1.};

  Double_t xminFD[nVarFD] = {0.,-1.,0.};
  Double_t xmaxFD[nVarFD] = {40.,1.,40.};

  //pt, y
  fMCAccPrompt = new THnSparseF("hMCAccPrompt","kStepMCAcceptance pt vs. y - promptD",nVarPrompt,nbinsPrompt,xminPrompt,xmaxPrompt);
  fMCAccPrompt->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fMCAccPrompt->GetAxis(1)->SetTitle("y");

  //pt,y,ptB
  fMCAccBFeed = new THnSparseF("hMCAccBFeed","kStepMCAcceptance pt vs. y vs. ptB - DfromB",nVarFD,nbinsFD,xminFD,xmaxFD);
  fMCAccBFeed->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fMCAccBFeed->GetAxis(1)->SetTitle("y");
  fMCAccBFeed->GetAxis(2)->SetTitle("p_{T}^{B} (GeV/c)");

  fOutputMass->Add(fMCAccPrompt);
  fOutputMass->Add(fMCAccBFeed);
}
//___________________________________________________________________________________________________
void AliAnalysisTaskSED0BDT::FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader){
  /// Fill MC acceptance histos for cuts study
  const Int_t nProng = 2;
  Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC

  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    if (TMath::Abs(mcPart->GetPdgCode()) == 421){

      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,fUseQuarkTagInKine);//Prompt = 4, FeedDown = 5

      Int_t deca = 0;
      Bool_t isGoodDecay=kFALSE;
      Int_t labDau[4]={-1,-1,-1,-1};
      Bool_t isInAcc = kFALSE;
      Bool_t isFidAcc = kFALSE;

      deca=AliVertexingHFUtils::CheckD0Decay(arrayMC,mcPart,labDau);
      if(deca > 0) isGoodDecay=kTRUE;

      if(labDau[0]==-1){
        continue; //protection against unfilled array of labels
      }

      isFidAcc=fCuts->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y());
      isInAcc=CheckAcc(arrayMC,nProng,labDau);

      if(isGoodDecay && TMath::Abs(zMCVertex) < fCuts->GetMaxVtxZ() && isFidAcc && isInAcc) {
        //for prompt
        if(orig == 4){
          //fill histo for prompt
          Double_t arrayMCprompt[2] = {mcPart->Pt(),mcPart->Y()};
          fMCAccPrompt->Fill(arrayMCprompt);
        }
        //for FD
        else if(orig == 5){
          Double_t ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,mcPart);
          //fill histo for FD
          Double_t arrayMCFD[3] = {mcPart->Pt(),mcPart->Y(),ptB};
          fMCAccBFeed->Fill(arrayMCFD);
        }
        else
          continue;
      }
    }
  }
}
//______________________________________________________________________
Bool_t AliAnalysisTaskSED0BDT::CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) return kFALSE;
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1) return kFALSE;
  }
  return kTRUE;
}
//________________________________________
void AliAnalysisTaskSED0BDT::ProcessBDT(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part,TClonesArray *arrMC, Float_t multi){
    fDaughterTracks.AddAt((AliAODTrack*)part->GetDaughter(0), 0);
    fDaughterTracks.AddAt((AliAODTrack*)part->GetDaughter(1), 1);
    AliAODTrack *prong2=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
    AliAODTrack *prong3=(AliAODTrack*)fDaughterTracks.UncheckedAt(1);
    
    Double_t normIP[2];
    Double_t d0Prong[2];
    Double_t ptProng[2];
    
    Double_t invmassD0 = part->InvMassD0(); Double_t invmassD0bar = part->InvMassD0bar();
    Double_t ptB;
    
    Double_t cosPointingAngle;
    Double_t cosThetaStarD0 = 99;
    Double_t cosThetaStarD0bar = 99;
    
    Double_t diffIP[2], errdiffIP[2];
    part->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),diffIP[0],errdiffIP[0]);
    part->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),diffIP[1],errdiffIP[1]);
    normIP[0]=diffIP[0]/errdiffIP[0];
    normIP[1]=diffIP[1]/errdiffIP[1];
    
    AliAODVertex *secvtx  = part->GetSecondaryVtx();
    AliAODVertex *primvtx = part->GetPrimaryVtx();
    Double_t err2decaylength = secvtx->Error2DistanceXYToVertex(primvtx);
    Double_t lxy = part->AliAODRecoDecay::DecayLengthXY(primvtx);
    Bool_t isusepid = fCuts->GetIsUsePID();
    //peng
    d0Prong[0]=part->Getd0Prong(0); d0Prong[1]=part->Getd0Prong(1); //d0*d0 d0Prong[1]*d0Prong[0]
    cosPointingAngle = part->CosPointingAngle();
    ptProng[0]=part->Pt2Prong(0); ptProng[1]=part->Pt2Prong(1);
    cosThetaStarD0 = part->CosThetaStarD0();
    cosThetaStarD0bar = part->CosThetaStarD0bar();
    //  if (part->Pt() > 10) return;
    // DCA part->GetDCA();
    Float_t tmp[23];
    tmp[8]= -99;                                // Invariant Mass
    tmp[18] = -99;                                 // ptB (if accessible)
    tmp[19] = 0;                                 // PDGCode
    tmp[20] = -99;                                // Rapidity YD0
    tmp[21] = -99;                                // Azimuthal Phi
    
    tmp[0] = part->Pt();                         // ptD
    tmp[1] = normIP[0];                         // Normalized d0N-<d0N> (topo1)
    tmp[2] = normIP[1];                         // Normalized d0P-<d0P> (topo2)
    tmp[3] = lxy;                                 // Decay length
    tmp[4] = lxy/TMath::Sqrt(err2decaylength);     // Normalized decay length
    tmp[9] = d0Prong[0]*d0Prong[1];             // d0N*d0P
    tmp[10] = cosPointingAngle;                 // CosThetaPointing
    tmp[11] = part->GetDCA();                     // DCAtracks
    tmp[12] = prong2->Pt();                     // pt1
    tmp[13] = prong3->Pt();                     // pt2
    tmp[14] = part->CosPointingAngleXY();         // CosThetaPointingXY
    tmp[15] = d0Prong[0];                         // d01
    tmp[16] = d0Prong[1];                         // d02
    tmp[20] = part->YD0();
    tmp[21] = part->Phi();
    tmp[22] = multi;
    
    if(tmp[0]<fBDTPtCut[0]||tmp[0]>=fBDTPtCut[1]) return;		// Global pT cut
    // PID and Cuts
	if(isusepid)fCuts->SetUsePID(kFALSE);// if PID on, switch it off
	Int_t isCuts=fCuts->IsSelected(part,AliRDHFCuts::kAll,aod);
	if(isusepid)fCuts->SetUsePID(kTRUE);//if PID was on, switch it on
	Int_t isPid= fCuts->IsSelectedPID(part);
	tmp[5] = (Double_t)isCuts;
	tmp[6] = (Double_t)isPid;
    
    //~ std::vector<TString> vnameString; vnameString.resize(23);
    /// "ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi
    //~ vnameString.push_back("ptD"); 	vnameString.push_back("topo1"); vnameString.push_back("topo2"); 	vnameString.push_back("lxy"); 	vnameString.push_back("nlxy");
    //~ vnameString.push_back("iscut"); vnameString.push_back("ispid"); vnameString.push_back("type"); 		vnameString.push_back("mass"); 	vnameString.push_back("d0d0");
    //~ vnameString.push_back("cosp"); 	vnameString.push_back("dca"); 	vnameString.push_back("ptk"); 		vnameString.push_back("ptpi"); 	vnameString.push_back("cospxy");
    //~ vnameString.push_back("d0k"); 	vnameString.push_back("d0pi"); 	vnameString.push_back("cosstar"); 	vnameString.push_back("ptB"); 	vnameString.push_back("pdgcode");
    //~ vnameString.push_back("YD0"); 	vnameString.push_back("phi");

    if(fReadMC){  // MC THnSparse
        
        Int_t pdgDgD0toKpi[2]={321,211};
        Int_t lab=-9999;
        // MC fill this
        TNtuple *NtupleD0C = (TNtuple*)fListBDTNtuple->FindObject("NtupleD0C");
        TNtuple *NtupleD0B = (TNtuple*)fListBDTNtuple->FindObject("NtupleD0B");
        TNtuple *NtupleRefl = (TNtuple*)fListBDTNtuple->FindObject("NtupleRefl");
        lab=part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
        if(lab>=0){
            
            AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(lab);
            if(!isCuts) return;
            
            // PDGCode
            Int_t pdgcode = partD0->GetPdgCode();
            tmp[19] = (Float_t)pdgcode;
            if (pdgcode == 421) {tmp[17] = cosThetaStarD0; tmp[8] = invmassD0;}
            else if (pdgcode == -421) {tmp[17] = cosThetaStarD0bar; tmp[8] = invmassD0bar;}
            
            // Check feed down or prompt
            Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partD0,fUseQuarkTagInKine);
            if(orig==4)           tmp[7] = 4; //prompt
            else if(orig==5)       tmp[7] = 5; //feed down
            else tmp[7] = 0;
            
            // pT of Bmother if Accessible
            ptB=AliVertexingHFUtils::GetBeautyMotherPt(arrMC,partD0);
            tmp[18] = ptB;
            
            if(tmp[7]==4)         NtupleD0C->Fill(tmp);
            else if(tmp[7]==5)     NtupleD0B->Fill(tmp);
            
            // Check reflection
            if((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2){
                tmp[8] = invmassD0;
                if(pdgcode!=421) NtupleRefl->Fill(tmp); // Reflection
            }
            if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)){
                tmp[8] = invmassD0bar;
                if(pdgcode!=-421) NtupleRefl->Fill(tmp); // Reflection
            }
        }
    }
    else  // data
    {
		if(fSampleSideband)// Sideband sampling
			if(gRandom->Rndm()>fBDTSidebandSamplingFraction) return;
			
		if(!fIsSelectedCandidate) return;

		std::vector<Double_t> BDTClsVar;// BDT cls input
        if(fSys == 0)    BDTClsVar.resize(11);
        if(fSys == 1)    BDTClsVar.resize(10);
 
        if((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2){  
            tmp[7] = 1; tmp[8] = invmassD0; tmp[17] = cosThetaStarD0;
            if(tmp[8]>2.12||tmp[8]<1.65) return;
            
            // Link variables to be used as classifier
            // NOTE: for 2018 Pb-Pb the decay length lxy was not applied(tmp[4])
            if(fSys == 0){      BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[3];     BDTClsVar[3] = tmp[4];  BDTClsVar[4] = tmp[9];     BDTClsVar[5] = tmp[10];
                BDTClsVar[6] = tmp[11]; BDTClsVar[7] = tmp[14]; BDTClsVar[8] = tmp[15]; BDTClsVar[9] = tmp[16]; BDTClsVar[10] = tmp[17];}
            if(fSys == 1){      BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[3];      BDTClsVar[3] = tmp[9];     BDTClsVar[4] = tmp[10];
                BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];}
			    
            if(fSampleSideband){ // Sideband sampling
				TNtuple *NtupleSB = (TNtuple*)fListBDTNtuple->FindObject("NtupleSB");
				NtupleSB->Fill(tmp);
			}
			else{ // Data application
				Int_t thisptbin = fCut4BDTptbin->PtBin(tmp[0]);
				if(thisptbin<0) return;
				Float_t *ptbin = fCut4BDTptbin->GetPtBinLimits();
				TString ptstring = Form("_%.0f_%.0f",ptbin[thisptbin],ptbin[thisptbin+1]);
				Int_t NBDT = fListBDTNames->GetEntries();
				TString BDT1Name = fListBDTNames->At(0)->GetName();
				AliRDHFBDT *thisbdt1   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("pT_%d_%s",thisptbin,BDT1Name.Data()));
				Float_t bdt1resp = thisbdt1->GetResponse(BDTClsVar);
						
				for(Int_t ii=1;ii<NBDT;ii++){
					TString BDT2Name = fListBDTNames->At(ii)->GetName();
					AliRDHFBDT *thisbdt2   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("pT_%d_%s",thisptbin,BDT2Name.Data()));
					Float_t bdt2resp = thisbdt2->GetResponse(BDTClsVar);
                    TH3F *thish3 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_19 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_19",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_1029 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_1029",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_3059 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_3059",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_19999 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_19999",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_6099 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_6099",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    thish3->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=1&&tmp[22]<10) thish3_19->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=10&&tmp[22]<30) thish3_1029->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=30&&tmp[22]<60) thish3_3059->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=60&&tmp[22]<100) thish3_6099->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=1&&tmp[22]<10000) thish3_19999->Fill(tmp[8],bdt1resp,bdt2resp);
					// Test output info
					//~ cout<<"INFO: "<<BDT1Name.Data()<<" = "<<bdt1resp<<", "<<BDT2Name.Data()<<" = "<<bdt2resp<<endl;
					//~ cout<<"INFO: Filling TH3F "<<thish3->GetName()<<endl;
					//~ printf("INFO: %s = %.3f, %s = %.3f\n",BDT1Name.Data(),bdt1resp,BDT2Name.Data(),bdt2resp);
				}
			}
        }
        if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)){
            tmp[7] = 2; tmp[8] = invmassD0bar; tmp[17] = cosThetaStarD0bar;
            
            if(tmp[8]>2.12||tmp[8]<1.65) return;
            
            // Link variables to be used as classifier
            // NOTE: for 2018 Pb-Pb the decay length lxy was not applied(tmp[4])
            if(fSys == 0){      BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[3];     BDTClsVar[3] = tmp[4];  BDTClsVar[4] = tmp[9];     BDTClsVar[5] = tmp[10];
                BDTClsVar[6] = tmp[11]; BDTClsVar[7] = tmp[14]; BDTClsVar[8] = tmp[15]; BDTClsVar[9] = tmp[16]; BDTClsVar[10] = tmp[17];}
            if(fSys == 1){      BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[3];      BDTClsVar[3] = tmp[9];     BDTClsVar[4] = tmp[10];
                BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];}
			    
            if(fSampleSideband){ // Sideband sampling
				TNtuple *NtupleSB = (TNtuple*)fListBDTNtuple->FindObject("NtupleSB");
				NtupleSB->Fill(tmp);
			}
			else{ // Data application
				Int_t thisptbin = fCut4BDTptbin->PtBin(tmp[0]);
				if(thisptbin<0) return;
				Float_t *ptbin = fCut4BDTptbin->GetPtBinLimits();
				TString ptstring = Form("_%.0f_%.0f",ptbin[thisptbin],ptbin[thisptbin+1]);
				Int_t NBDT = fListBDTNames->GetEntries();
				TString BDT1Name = fListBDTNames->At(0)->GetName();
				AliRDHFBDT *thisbdt1   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("pT_%d_%s",thisptbin,BDT1Name.Data()));
				Float_t bdt1resp = thisbdt1->GetResponse(BDTClsVar);
						
				for(Int_t ii=1;ii<NBDT;ii++){
					TString BDT2Name = fListBDTNames->At(ii)->GetName();
					AliRDHFBDT *thisbdt2   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("pT_%d_%s",thisptbin,BDT2Name.Data()));
					Float_t bdt2resp = thisbdt2->GetResponse(BDTClsVar);
                    TH3F *thish3 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_19 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_19",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_1029 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_1029",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_3059 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_3059",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_19999 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_19999",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    TH3F *thish3_6099 = (TH3F*)fListBDTResp->FindObject(Form("h3MassRespPt%d_%s_%s_6099",thisptbin,BDT1Name.Data(),BDT2Name.Data()));
                    thish3->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=1&&tmp[22]<10) thish3_19->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=10&&tmp[22]<30) thish3_1029->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=30&&tmp[22]<60) thish3_3059->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=60&&tmp[22]<100) thish3_6099->Fill(tmp[8],bdt1resp,bdt2resp);
                    if(tmp[22]>=1&&tmp[22]<10000) thish3_19999->Fill(tmp[8],bdt1resp,bdt2resp);

					// Test output info
					//~ cout<<"INFO: "<<BDT1Name.Data()<<" = "<<bdt1resp<<", "<<BDT2Name.Data()<<" = "<<bdt2resp<<endl;
					//~ cout<<"INFO: Filling TH3F "<<thish3->GetName()<<endl;
					//~ printf("INFO: %s = %.3f, %s = %.3f\n",BDT1Name.Data(),bdt1resp,BDT2Name.Data(),bdt2resp);
				}
			}
        }
    }
    
}
//_________________________________________________________________________
TProfile* AliAnalysisTaskSED0BDT::GetEstimatorHistogram(const AliVEvent* event){
  /// Get Estimator Histogram from period event->GetRunNumber();
  ///
  /// If you select SPD tracklets in |eta|<1 you should use type == 1
  ///
    
  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;   // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  // pPb 2013: 0-LHC13b, 1-LHC13c
  // pPb 2016: 0-LHC16q: 265499->265525; 265309->265387, 1-LHC16q:265435, 2-LHC16q:265388->265427, LHC16t: 267163->267166


    if(fYearNumber==10){
    if(runNo>114930 && runNo<117223) period = 0;
    if(runNo>119158 && runNo<120830) period = 1;
    if(runNo>122373 && runNo<126438) period = 2;
    if(runNo>127711 && runNo<130851) period = 3;
    if(period<0 || period>3) return 0;
    }else if(fYearNumber==16){
    if(runNo>=252235 && runNo<=252375)period = 0;//16d
    if(runNo>=252603 && runNo<=253591)period = 1;//16e
    if(runNo>=254124 && runNo<=254332)period = 2;//16g
    if(runNo>=254378  && runNo<=255469 )period = 3;//16h_1
    if(runNo>=254418  && runNo<=254422 )period = 4;//16h_2 negative mag
    if(runNo>=256146  && runNo<=256420 )period = 5;//16j
    if(runNo>=256504  && runNo<=258537 )period = 6;//16k
    if(runNo>=258883  && runNo<=260187)period = 7;//16l
    if(runNo>=262395  && runNo<=264035 )period = 8;//16o
    if(runNo>=264076  && runNo<=264347 )period = 9;//16p
    }else if(fYearNumber==17){
    if(runNo>=270822 && runNo<=270830)period = 0;//17e
    if(runNo>=270854 && runNo<=270865)period = 1;//17f
    if(runNo>=271868 && runNo<=273103)period = 2;//17h
    if(runNo>=273591  && runNo<=274442)period = 3;//17i
    if(runNo>=274593  && runNo<=274671)period = 4;//17j
    if(runNo>=274690  && runNo<=276508)period = 5;//17k
    if(runNo>=276551  && runNo<=278216)period = 6;//17l
    if(runNo>=278914  && runNo<=280140)period = 7;//17m
    if(runNo>=280282   && runNo<=281961)period = 8;//17o
    if(runNo>=282504  && runNo<=282704)period = 9;//17r
    }else if(fYearNumber==18){
      if(runNo>=285008 && runNo<=285447)period = 0;//18b
    if(runNo>=285978 && runNo<=286350)period = 1;//18d
    if(runNo>=286380 && runNo<=286937)period = 2;//18e
    if(runNo>=287000  && runNo<=287977)period = 3;//18f
    if(runNo>=288619  && runNo<=288750)period = 4;//18g
    if(runNo>=288804  && runNo<=288806)period = 5;//18h
    if(runNo>=288861  && runNo<=288909 )period = 6;//18i
    if(runNo==288943)period = 7;//18j
    if(runNo>=289165   && runNo<=289201)period = 8;//18k
    if(runNo>=289240  && runNo<=289971)period = 9;//18l
    if(runNo>=290222  && runNo<=292839)period = 10;//18m
    if(runNo>=293357   && runNo<=293359)period = 11;//18n
    if(runNo>=293368   && runNo<=293898)period = 12;//18o
    if(runNo>=294009  && runNo<=294925)period = 13;//18p
  }
  
  return fMultEstimatorAvg[period];
}
