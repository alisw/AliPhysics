/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliAnalysisTaskSEDvsMultiplicity_BDT
// AliAnalysisTaskSE for the charmed haadrons vs. multiplcity analysis
// Authors: Renu Bala, Zaida Conesa del Valle, Francesco Prino
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDvsMultiplicity_BDT.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"
#include "AliInputEventHandler.h"
#include "AliRDHFBDT.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDvsMultiplicity_BDT);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity_BDT::AliAnalysisTaskSEDvsMultiplicity_BDT():
  AliAnalysisTaskSE(),
  fOutput(0),
  fListCuts(0),
  fOutputCounters(0),
  fListProfiles(0),
  fHistNEvents(0),
  fHistNtrEta16vsNtrEta1EvSel(0),
  fHistNtrEta05vsNtrEta1EvSel(0),
  fHistNtrEta03vsNtrEta1EvSel(0),
  fHistNtrEtaV0AvsNtrEta1EvSel(0),
  fHistNtrEtaV0MvsNtrEta1EvSel(0),
  fHistNtrEtaV0AvsV0AEqEvSel(0),
  fHistNtrEtaV0MvsV0MEqEvSel(0),
  fHistNtrCorrEta1vsNtrRawEta1EvSel(0),
  fHistMultCorrvsMultRawEvSel(0),
  fHistNtrEta16vsNtrEta1EvWithCand(0),
  fHistNtrEta05vsNtrEta1EvWithCand(0),
  fHistNtrEta03vsNtrEta1EvWithCand(0),
  fHistNtrEtaV0AvsNtrEta1EvWithCand(0),
  fHistNtrEtaV0MvsNtrEta1EvWithCand(0),
  fHistNtrEtaV0AvsV0AEqEvWithCand(0),
  fHistNtrEtaV0MvsV0MEqEvWithCand(0),
  fHistNtrCorrEta1vsNtrRawEta1EvWithCand(0),
  fHistMultCorrvsMultRawEvWithCand(0),
  fHistNtrEta16vsNtrEta1EvWithD(0),
  fHistNtrEta05vsNtrEta1EvWithD(0),
  fHistNtrEta03vsNtrEta1EvWithD(0),
  fHistNtrEtaV0AvsNtrEta1EvWithD(0),
  fHistNtrEtaV0MvsNtrEta1EvWithD(0),
  fHistNtrEtaV0AvsV0AEqEvWithD(0),
  fHistNtrEtaV0MvsV0MEqEvWithD(0),
  fHistNtrCorrEta1vsNtrRawEta1EvWithD(0),
  fHistMultCorrvsMultRawEvWithD(0),
  fHistNtrVsZvtx(0),
  fHistNtrCorrVsZvtx(0),
  fHistNtrVsNchMC(0),
  fHistNtrCorrVsNchMC(0),
  fHistNtrVsNchMCPrimary(0),
  fHistNtrCorrVsNchMCPrimary(0),
  fHistNtrVsNchMCPhysicalPrimary(0),
  fHistNtrCorrVsNchMCPhysicalPrimary(0),
  fHistGenPrimaryParticlesInelGt0(0),
  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
  fHistNtrUnCorrPSSel(0),
  fHistNtrUnCorrPSTrigSel(0),
  fHistNtrUnCorrPSTrigPileUpSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxContSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxRangeSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel(0),
  fHistNtrUnCorrEvSel(0),
  fHistNtrUnCorrEvWithCand(0),
  fHistNtrUnCorrEvWithD(0),
  fHistNtrCorrPSSel(0),
  fHistNtrCorrEvSel(0),
  fHistNtrCorrEvWithCand(0),
  fHistNtrCorrEvWithD(0),
  fPtVsMassVsMult(0),
  fPtVsMassVsMultNoPid(0),
  fPtVsMassVsMultUncorr(0),
  fPtVsMassVsMultPart(0),
  fPtVsMassVsMultAntiPart(0),
  fPtVsMassVsMultMC(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNMassBins(200),
  fRDCutsAnalysis(0),
  fCounter(0),
  fCounterC(0),
  fCounterU(0),
  fCounterCandidates(0),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fReadMC(kFALSE),
  fMCOption(0),
  fisPPbData(kFALSE),
  fUseBit(kTRUE),
  fSubtractTrackletsFromDau(kFALSE),
  fKeepCorrPlots(kFALSE),
  fAODProtection(1),
  fUseNchWeight(0),
  fHistoMCNch(0),
  fHistoMeasNch(0),
  fRefMult(9.26),
  fPdgMeson(411),
  fLctoV0(kTRUE),
  fMultiplicityEstimator(kNtrk10),
  fMCPrimariesEstimator(kEta10),
  fDoVZER0ParamVertexCorr(1),
  fYearNumber(16),
  fDaughterTracks(),
  fFillOnlyD0D0bar(0),
  fUseQuarkTagInKine(kTRUE),
  fListRDHFBDT(0),
  fListBDTNtuple(0),
  fListBDTResp(0),
  fBDTRespCut(-1),
  fBDTSidebandSamplingFraction(0.1),
  fSampleSideband(kFALSE),
  fGetRespTree(kTRUE),
  fBDTFullVarString("ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi"),
  fBDTClassifierVarString(""),
  fApplyBDT(0),
  fFillTree(0)
{
  /// Default constructor
  for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;
  fBDTPtCut[0]=0; fBDTPtCut[1]=1e9;
  fBDTSidebandCut[0]=1.792; fBDTSidebandCut[1]=1.942;
}

//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity_BDT::AliAnalysisTaskSEDvsMultiplicity_BDT(const char *name, Int_t pdgMeson,AliRDHFCuts *cuts, Bool_t switchPPb):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fListCuts(0),
  fOutputCounters(0),
  fListProfiles(0),
  fHistNEvents(0),
  fHistNtrEta16vsNtrEta1EvSel(0),
  fHistNtrEta05vsNtrEta1EvSel(0),
  fHistNtrEta03vsNtrEta1EvSel(0),
  fHistNtrEtaV0AvsNtrEta1EvSel(0),
  fHistNtrEtaV0MvsNtrEta1EvSel(0),
  fHistNtrEtaV0AvsV0AEqEvSel(0),
  fHistNtrEtaV0MvsV0MEqEvSel(0),
  fHistNtrCorrEta1vsNtrRawEta1EvSel(0),
  fHistMultCorrvsMultRawEvSel(0),
  fHistNtrEta16vsNtrEta1EvWithCand(0),
  fHistNtrEta05vsNtrEta1EvWithCand(0),
  fHistNtrEta03vsNtrEta1EvWithCand(0),
  fHistNtrEtaV0AvsNtrEta1EvWithCand(0),
  fHistNtrEtaV0MvsNtrEta1EvWithCand(0),
  fHistNtrEtaV0AvsV0AEqEvWithCand(0),
  fHistNtrEtaV0MvsV0MEqEvWithCand(0),
  fHistNtrCorrEta1vsNtrRawEta1EvWithCand(0),
  fHistMultCorrvsMultRawEvWithCand(0),
  fHistNtrEta16vsNtrEta1EvWithD(0),
  fHistNtrEta05vsNtrEta1EvWithD(0),
  fHistNtrEta03vsNtrEta1EvWithD(0),
  fHistNtrEtaV0AvsNtrEta1EvWithD(0),
  fHistNtrEtaV0MvsNtrEta1EvWithD(0),
  fHistNtrEtaV0AvsV0AEqEvWithD(0),
  fHistNtrEtaV0MvsV0MEqEvWithD(0),
  fHistNtrCorrEta1vsNtrRawEta1EvWithD(0),
  fHistMultCorrvsMultRawEvWithD(0),
  fHistNtrVsZvtx(0),
  fHistNtrCorrVsZvtx(0),
  fHistNtrVsNchMC(0),
  fHistNtrCorrVsNchMC(0),
  fHistNtrVsNchMCPrimary(0),
  fHistNtrCorrVsNchMCPrimary(0),
  fHistNtrVsNchMCPhysicalPrimary(0),
  fHistNtrCorrVsNchMCPhysicalPrimary(0),
  fHistGenPrimaryParticlesInelGt0(0),
  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
  fHistNtrUnCorrPSSel(0),
  fHistNtrUnCorrPSTrigSel(0),
  fHistNtrUnCorrPSTrigPileUpSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxContSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxRangeSel(0),
  fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel(0),
  fHistNtrUnCorrEvSel(0),
  fHistNtrUnCorrEvWithCand(0),
  fHistNtrUnCorrEvWithD(0),
  fHistNtrCorrPSSel(0),
  fHistNtrCorrEvSel(0),
  fHistNtrCorrEvWithCand(0),
  fHistNtrCorrEvWithD(0),
  fPtVsMassVsMult(0),
  fPtVsMassVsMultNoPid(0),
  fPtVsMassVsMultUncorr(0),
  fPtVsMassVsMultPart(0),
  fPtVsMassVsMultAntiPart(0),
  fPtVsMassVsMultMC(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNMassBins(200),
  fRDCutsAnalysis(cuts),
  fCounter(0),
  fCounterC(0),
  fCounterU(0),
  fCounterCandidates(0),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fReadMC(kFALSE),
  fMCOption(0),
  fisPPbData(switchPPb),
  fUseBit(kTRUE),
  fSubtractTrackletsFromDau(kFALSE),
  fKeepCorrPlots(kFALSE),
  fAODProtection(1),
  fUseNchWeight(0),
  fHistoMCNch(0),
  fHistoMeasNch(0),
  fRefMult(9.26),
  fPdgMeson(pdgMeson),
  fLctoV0(kTRUE),
  fMultiplicityEstimator(kNtrk10),
  fMCPrimariesEstimator(kEta10),
  fDoVZER0ParamVertexCorr(1),
  fYearNumber(16),
  fDaughterTracks(),
  fFillOnlyD0D0bar(0),
  fUseQuarkTagInKine(kTRUE),
  fListRDHFBDT(0),
  fListBDTNtuple(0),
  fListBDTResp(0),
  fBDTRespCut(-1),
  fBDTSidebandSamplingFraction(0.1),
  fSampleSideband(kFALSE),
  fGetRespTree(kTRUE),
  fBDTFullVarString("ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi"),
  fBDTClassifierVarString(""),
  fApplyBDT(0),
  fFillTree(0)
{
  // 
  /// Standard constructor
  //
  fBDTPtCut[0]=0; fBDTPtCut[1]=1e9;
  fBDTSidebandCut[0]=1.792; fBDTSidebandCut[1]=1.942;  
 
  for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<14; i++) fMultEstimatorAvg[i]=0;
  if(fPdgMeson==413){
    fNMassBins=200;
    SetMassLimits(0.12,0.2);
  }else if(fPdgMeson == 431) {
    Double_t MassBinSize  = 0.002;
    Int_t nInvMassBins = (Int_t)(0.7/MassBinSize+0.5);
    Double_t massDs  = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    Double_t minMass = massDs-0.5*nInvMassBins*MassBinSize;
    Double_t maxMass = massDs+0.5*nInvMassBins*MassBinSize;
    SetMassLimits(minMass,maxMass);
    SetNMassBins(nInvMassBins); 
  }else if(fPdgMeson == 4122) {
    Double_t massLc  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Int_t nInvMassBins = 1000;
    Double_t minMass = massLc-0.250;
    Double_t maxMass = massLc+0.250;
    SetMassLimits(minMass,maxMass);
    SetNMassBins(nInvMassBins);
  }
  else{
    fNMassBins=200;
    SetMassLimits(fPdgMeson,0.1);
  }
  // Default constructor
  // Otput slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes cut to private output
  DefineOutput(2,TList::Class());
  // Output slot #3 writes cut to private output
  DefineOutput(3,TList::Class()); 
  // Output slot #4 writes cut to private output
  DefineOutput(4,TList::Class()); 
  DefineOutput(5,TList::Class());
  if(!fReadMC && fApplyBDT){
     DefineOutput(6,TList::Class());
  }
}
//________________________________________________________________________
AliAnalysisTaskSEDvsMultiplicity_BDT::~AliAnalysisTaskSEDvsMultiplicity_BDT()
{
  //
  /// Destructor
  //
  delete fOutput;
  delete fHistNEvents;
  delete fListCuts;
  delete fListProfiles;
  delete fRDCutsAnalysis;
  delete fCounter;
  delete fCounterC;
  delete fCounterU;
  delete fCounterCandidates;
  for(Int_t i=0; i<14; i++) {
    if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }
  
  for(Int_t i=0; i<5; i++){
    delete fHistMassPtImpPar[i];
  }
  if(fHistoMCNch) delete fHistoMCNch;
  if(fHistoMeasNch) delete fHistoMeasNch;

  if(fListRDHFBDT) {
     delete fListRDHFBDT;
     fListRDHFBDT = 0;
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

//_________________________________________________________________
void  AliAnalysisTaskSEDvsMultiplicity_BDT::SetMassLimits(Double_t lowlimit, Double_t uplimit){
  /// set invariant mass limits
  if(uplimit>lowlimit){
    fLowmasslimit = lowlimit;
    fUpmasslimit = uplimit;
  }else{
    AliError("Wrong mass limits: upper value should be larger than lower one");
  }
}
//_________________________________________________________________
void  AliAnalysisTaskSEDvsMultiplicity_BDT::SetMassLimits(Int_t pdg, Double_t range){
  /// set invariant mass limits
  Double_t mass=TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();
  SetMassLimits(mass-range,mass+range);
}
//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::Init(){
  //
  /// Initialization
  //
  printf("AnalysisTaskSEDvsMultiplicity_0::Init() \n");

  if(fUseNchWeight && !fReadMC){ AliFatal("Nch weights can only be used in MC mode"); return; }
  if(fUseNchWeight && !fHistoMCNch){ AliFatal("Nch weights can only be used without histogram"); return; }
  if(fUseNchWeight==1 && !fHistoMeasNch) {//Nch weights
    if(fisPPbData){ AliFatal("Nch weights can only be used with MC and data histogram in pPb"); return; }
    else CreateMeasuredNchHisto();  
  }
  if(fUseNchWeight==2 && !fHistoMeasNch){ AliFatal("Ntrk weights can only be used with MC and data histogram"); return; } //for pp, pPb Ntrk weights 

  
  fListCuts=new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("CutsList");


  if(fPdgMeson==411){
    AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDplus");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==421){
    AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDzero");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==413){
    AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDStar");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==431){
    AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDs");
    fListCuts->Add(copycut);
  }else if(fPdgMeson==4122){
    if(fLctoV0){
    AliRDHFCutsLctoV0* copycut=new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsLc2pK0S");
    fListCuts->Add(copycut);
  }else{
    AliRDHFCutsLctopKpi *copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fRDCutsAnalysis)));
    copycut->SetName("LctopKpiProdCuts");
    fListCuts->Add(copycut);
  }
  }
  if(fHistoMeasNch) fListCuts->Add(fHistoMeasNch);
  if(fHistoMCNch) fListCuts->Add(fHistoMCNch);

  PostData(2,fListCuts);
  
  fListProfiles = new TList();
  fListProfiles->SetOwner();
  TString period[14];
  Int_t nProfiles=14;
  if (fisPPbData) {
    if(fYearNumber == 13) {
      period[0]="LHC13b";
      period[1]="LHC13c";
      nProfiles = 2;
    } else if(fYearNumber == 16) {
      period[0]="LHC16qt_1stBunch";
      period[1]="LHC16qt_2ndBunch";
      period[2]="LHC16qt_3rdBunch";
      period[3]="LHC16qt_4thBunch";
      nProfiles = 4;
    }
  }
  else {
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
  }
  for(Int_t i=0; i<nProfiles; i++){
    if(fMultEstimatorAvg[i]){
      TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
      hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
      fListProfiles->Add(hprof);
    }
  }

  PostData(4,fListProfiles);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  Int_t nMultBins = 200;
  Float_t firstMultBin = -0.5;
  Float_t lastMultBin = 199.5;
  Int_t nMultBinsNtrk = nMultBins;
  Float_t lastMultBinNtrk = lastMultBin;
  Int_t nMultBinsV0 = 400;
  Float_t lastMultBinV0 = 799.5;
  const char *estimatorName="tracklets";
  if(fisPPbData) {
    nMultBinsNtrk = 375;
    lastMultBinNtrk = 374.5;
    nMultBins = nMultBinsNtrk;
    lastMultBin = lastMultBinNtrk;
  }
  if(fMultiplicityEstimator==kVZERO || fMultiplicityEstimator==kVZEROA || fMultiplicityEstimator==kVZEROEq || fMultiplicityEstimator==kVZEROAEq) {
    nMultBins = nMultBinsV0;
    lastMultBin = lastMultBinV0;
    estimatorName = "vzero";
  }

  fHistNtrUnCorrPSSel = new TH1F("hNtrUnCorrPSSel",Form("Uncorrected %s multiplicity for PS selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigSel = new TH1F("hNtrUnCorrPSTrigSel",Form("Uncorrected %s multiplicity for PS + trigger name selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigPileUpSel = new TH1F("hNtrUnCorrPSTrigPileUpSel",Form("Uncorrected %s multiplicity for PS + trigger name + pileup selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigPileUpVtxSel = new TH1F("hNtrUnCorrPSTrigPileUpVtxSel",Form("Uncorrected %s multiplicity for PS + trigger name + pileup + with-vertex selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigPileUpVtxContSel = new TH1F("hNtrUnCorrPSTrigPileUpVtxContSel",Form("Uncorrected %s multiplicity for PS + trigger name + pileup + with-vertex-contrib selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigPileUpVtxRangeSel = new TH1F("hNtrUnCorrPSTrigPileUpVtxRangeSel",Form("Uncorrected %s multiplicity for PS + trigger name + pileup + with-vertex-contrib-range selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel = new TH1F("hNtrUnCorrPSTrigPileUpVtxRangeCentrSel",Form("Uncorrected %s multiplicity for PS + trigger name + pileup + with-vertex-contrib-range + centrality selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrEvSel = new TH1F("hNtrUnCorrEvSel",Form("Uncorrected %s multiplicity for selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrUnCorrEvWithCand = new TH1F("hNtrUnCorrEvWithCand",Form("Uncorrected %s multiplicity for events with D candidates; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);// Total multiplicity
  fHistNtrUnCorrEvWithD = new TH1F("hNtrUnCorrEvWithD",Form("Uncorrected %s multiplicity for events with D in mass region ; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin); //
  fHistNtrCorrPSSel = new TH1F("hNtrCorrPSSel",Form("Corrected %s multiplicity for PS selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrCorrEvSel = new TH1F("hNtrCorrEvSel",Form("Corrected %s multiplicity for selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
  fHistNtrCorrEvWithCand = new TH1F("hNtrCorrEvWithCand", Form("%s multiplicity for events with D candidates; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);// Total multiplicity
  fHistNtrCorrEvWithD = new TH1F("hNtrCorrEvWithD", Form("%s multiplicity for events with D in mass region ; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin); //

  if(fKeepCorrPlots){
    fHistNtrEta16vsNtrEta1EvSel = new TH2F("hNtrEta16vsNtrEta1EvSel","Uncorrected Eta1.6 vs Eta1.0 (events selected); Ntracklets #eta<1.0; Ntracklets #eta<1.6",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram 
    fHistNtrEta05vsNtrEta1EvSel = new TH2F("hNtrEta05vsNtrEta1EvSel","Uncorrected Eta0.5 vs Eta1.0 (events selected); Ntracklets #eta<1.0; Ntracklets #eta<0.5",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.5 vs eta 1.0 histogram 
    fHistNtrEta03vsNtrEta1EvSel = new TH2F("hNtrEta03vsNtrEta1EvSel","Uncorrected Eta0.3 vs Eta1.0 (events selected); Ntracklets #eta<1.0; Ntracklets #eta<0.3",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.3 vs eta 1.0 histogram 
    fHistNtrEtaV0AvsNtrEta1EvSel = new TH2F("hNtrEtaV0AvsNtrEta1EvSel","Uncorrected Eta-V0A vs Eta1.0 (events selected); Ntracklets #eta<1.0; Multiplicity V0A",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0A vs eta 1.0 histogram 
    fHistNtrEtaV0MvsNtrEta1EvSel = new TH2F("hNtrEtaV0MvsNtrEta1EvSel","Uncorrected Eta-V0M vs Eta1.0 (events selected); Ntracklets #eta<1.0; Multiplicity V0A+V0C",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0M vs eta 1.0 histogram
    fHistNtrEtaV0AvsV0AEqEvSel = new TH2F("hNtrEtaV0AvsV0AEqEvSel","Corrected V0A vs corrected V0A-Equalized (events selected); Vzero-A; Vzero-A Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0A - V0Aeq
    fHistNtrEtaV0MvsV0MEqEvSel = new TH2F("hNtrEtaV0MvsV0MEqEvSel","Corrected V0M vs corrected V0M-Equalized (events selected); Vzero-M; Vzero-M Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0M - V0Meq
    fHistNtrCorrEta1vsNtrRawEta1EvSel = new TH2F("hNtrCorrEta1vsNtrRawEta1EvSel","Corrected Eta1 vs Eta1.0 (events selected); Ntracklets #eta<1.0 corrected; Ntracklets #eta<1",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram
    fHistMultCorrvsMultRawEvSel = new TH2F("hMultCorrvsMultRawEvSel",Form("Corrected multiplicity vs uncorrected multiplicity (events selected); %s corrected; %s",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); // corrected vs uncorrected multiplicity

    fHistNtrEta16vsNtrEta1EvWithCand = new TH2F("hNtrEta16vsNtrEta1EvWithCand","Uncorrected Eta1.6 vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0; Ntracklets #eta<1.6",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram 
    fHistNtrEta05vsNtrEta1EvWithCand = new TH2F("hNtrEta05vsNtrEta1EvWithCand","Uncorrected Eta0.5 vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0; Ntracklets #eta<0.5",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.5 vs eta 1.0 histogram 
    fHistNtrEta03vsNtrEta1EvWithCand = new TH2F("hNtrEta03vsNtrEta1EvWithCand","Uncorrected Eta0.3 vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0; Ntracklets #eta<0.3",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.3 vs eta 1.0 histogram 
    fHistNtrEtaV0AvsNtrEta1EvWithCand = new TH2F("hNtrEtaV0AvsNtrEta1EvWithCand","Uncorrected Eta-V0A vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0; Multiplicity V0A",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0A vs eta 1.0 histogram 
    fHistNtrEtaV0MvsNtrEta1EvWithCand = new TH2F("hNtrEtaV0MvsNtrEta1EvWithCand","Uncorrected Eta-V0M vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0; Multiplicity V0A+V0C",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0M vs eta 1.0 histogram 
    fHistNtrEtaV0AvsV0AEqEvWithCand = new TH2F("hNtrEtaV0AvsV0AEqEvWithCand","Corrected V0A vs corrected V0A-Equalized (events selected with a D candidate); Vzero-A; Vzero-A Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0A - V0Aeq
    fHistNtrEtaV0MvsV0MEqEvWithCand = new TH2F("hNtrEtaV0MvsV0MEqEvWithCand","Corrected V0M vs corrected V0M-Equalized (events selected with a D candidate); Vzero-M; Vzero-M Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0M - V0Meq
    fHistNtrCorrEta1vsNtrRawEta1EvWithCand = new TH2F("hNtrCorrEta1vsNtrRawEta1EvWithCand","Corrected Eta1 vs Eta1.0 (events selected with a D candidate); Ntracklets #eta<1.0 corrected; Ntracklets #eta<1",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram
    fHistMultCorrvsMultRawEvWithCand = new TH2F("hMultCorrvsMultRawEvWithCand",Form("Corrected multiplicity vs uncorrected multiplicity (events selected) with a D candidate; %s corrected; %s",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); // corrected vs uncorrected multiplicity


    fHistNtrEta16vsNtrEta1EvWithD = new TH2F("hNtrEta16vsNtrEta1EvWithD","Uncorrected Eta1.6 vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0; Ntracklets #eta<1.6",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram 
    fHistNtrEta05vsNtrEta1EvWithD = new TH2F("hNtrEta05vsNtrEta1EvWithD","Uncorrected Eta0.5 vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0; Ntracklets #eta<0.5",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.5 vs eta 1.0 histogram 
    fHistNtrEta03vsNtrEta1EvWithD = new TH2F("hNtrEta03vsNtrEta1EvWithD","Uncorrected Eta0.3 vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0; Ntracklets #eta<0.3",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 0.3 vs eta 1.0 histogram 
    fHistNtrEtaV0AvsNtrEta1EvWithD = new TH2F("hNtrEtaV0AvsNtrEta1EvWithD","Uncorrected Eta-V0A vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0; Multiplicity V0A",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0A vs eta 1.0 histogram 
    fHistNtrEtaV0MvsNtrEta1EvWithD = new TH2F("hNtrEtaV0MvsNtrEta1EvWithD","Uncorrected Eta-V0M vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0; Multiplicity V0A+V0C",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsV0,firstMultBin,lastMultBinV0); //eta V0M vs eta 1.0 histogram 
    fHistNtrEtaV0AvsV0AEqEvWithD = new TH2F("hNtrEtaV0AvsV0AEqEvWithD","Corrected V0A vs corrected V0A-Equalized (events selected with D in mass range); Vzero-A; Vzero-A Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0A - V0Aeq
    fHistNtrEtaV0MvsV0MEqEvWithD = new TH2F("hNtrEtaV0MvsV0MEqEvWithD","Corrected V0M vs corrected V0M-Equalized (events selected with D in mass range); Vzero-M; Vzero-M Equalized",nMultBinsV0,firstMultBin,lastMultBinV0,nMultBinsV0,firstMultBin,lastMultBinV0); // comparison V0M - V0Meq
    fHistNtrCorrEta1vsNtrRawEta1EvWithD = new TH2F("hNtrCorrEta1vsNtrRawEta1EvWithD","Corrected Eta1 vs Eta1.0 (events selected with D in mass range); Ntracklets #eta<1.0 corrected; Ntracklets #eta<1",nMultBinsNtrk,firstMultBin,lastMultBinNtrk,nMultBinsNtrk,firstMultBin,lastMultBinNtrk); //eta 1.6 vs eta 1.0 histogram
    fHistMultCorrvsMultRawEvWithD = new TH2F("hMultCorrvsMultRawEvWithD",Form("Corrected multiplicity vs uncorrected multiplicity (events selected with D in mass range); %s corrected; %s",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); // corrected vs uncorrected multiplicity

  }
  fHistNtrVsZvtx = new TH2F("hNtrVsZvtx",Form("N%s vs VtxZ; VtxZ;N_{%s};",estimatorName,estimatorName),300,-15,15,nMultBins,firstMultBin,lastMultBin); //
  fHistNtrCorrVsZvtx = new TH2F("hNtrCorrVsZvtx",Form("N%s vs VtxZ; VtxZ;N_{%s};",estimatorName,estimatorName),300,-15,15,nMultBins,firstMultBin,lastMultBin); //

  fHistNtrVsNchMC = new TH2F("hNtrVsNchMC",Form("N%s vs NchMC; Nch;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
  fHistNtrCorrVsNchMC = new TH2F("hNtrCorrVsNchMC",Form("N%s vs Nch; Nch;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
  
  fHistNtrVsNchMCPrimary = new TH2F("hNtrVsNchMCPrimary",Form("N%s vs Nch (Primary); Nch (Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
  fHistNtrCorrVsNchMCPrimary = new TH2F("hNtrCorrVsNchMCPrimary",Form("N%s vs Nch (Primary); Nch(Primary) ;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //

  fHistNtrVsNchMCPhysicalPrimary = new TH2F("hNtrVsNchMCPhysicalPrimary",Form("N%s vs Nch (Physical Primary); Nch (Physical Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
  fHistNtrCorrVsNchMCPhysicalPrimary = new TH2F("hNtrCorrVsMCPhysicalPrimary",Form("N%s vs Nch (Physical Primary); Nch (Physical Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
  
  fHistGenPrimaryParticlesInelGt0 = new TH1F("hGenPrimaryParticlesInelGt0","Multiplcity of generated charged particles ; Nparticles ; Entries",nMultBins,firstMultBin,lastMultBin);

  fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary = new TH3F("fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary", "MC: Nch (Physical Primary) vs Nch (Primary) vs Nch (Generated); Nch (Generated); Nch (Primary); Nch (Physical Primary)",nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin);

  fHistNtrUnCorrPSSel->Sumw2();
  fHistNtrUnCorrPSTrigSel->Sumw2();
  fHistNtrUnCorrPSTrigPileUpSel->Sumw2();
  fHistNtrUnCorrPSTrigPileUpVtxSel->Sumw2();
  fHistNtrUnCorrPSTrigPileUpVtxContSel->Sumw2();
  fHistNtrUnCorrPSTrigPileUpVtxRangeSel->Sumw2();
  fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel->Sumw2();
  fHistNtrUnCorrEvSel->Sumw2();
  fHistNtrUnCorrEvWithCand->Sumw2();
  fHistNtrUnCorrEvWithD->Sumw2();
  fHistNtrCorrPSSel->Sumw2();
  fHistNtrCorrEvSel->Sumw2();
  fHistNtrCorrEvWithCand->Sumw2();
  fHistNtrCorrEvWithD->Sumw2();
  fHistGenPrimaryParticlesInelGt0->Sumw2();
  fOutput->Add(fHistNtrUnCorrPSSel);
  fOutput->Add(fHistNtrUnCorrPSTrigSel);
  fOutput->Add(fHistNtrUnCorrPSTrigPileUpSel);
  fOutput->Add(fHistNtrUnCorrPSTrigPileUpVtxSel);
  fOutput->Add(fHistNtrUnCorrPSTrigPileUpVtxContSel);
  fOutput->Add(fHistNtrUnCorrPSTrigPileUpVtxRangeSel);
  fOutput->Add(fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel);
  fOutput->Add(fHistNtrUnCorrEvSel);
  fOutput->Add(fHistNtrUnCorrEvWithCand);
  fOutput->Add(fHistNtrUnCorrEvWithD);
  fOutput->Add(fHistNtrCorrPSSel);
  fOutput->Add(fHistNtrCorrEvSel);
  fOutput->Add(fHistNtrCorrEvWithCand);
  fOutput->Add(fHistNtrCorrEvWithD);
  if(fKeepCorrPlots){
    fOutput->Add(fHistNtrEta16vsNtrEta1EvSel);
    fOutput->Add(fHistNtrEta05vsNtrEta1EvSel);
    fOutput->Add(fHistNtrEta03vsNtrEta1EvSel);
    fOutput->Add(fHistNtrEtaV0AvsNtrEta1EvSel);
    fOutput->Add(fHistNtrEtaV0MvsNtrEta1EvSel);
    fOutput->Add(fHistNtrEtaV0AvsV0AEqEvSel);
    fOutput->Add(fHistNtrEtaV0MvsV0MEqEvSel);
    fOutput->Add(fHistNtrCorrEta1vsNtrRawEta1EvSel);
    fOutput->Add(fHistMultCorrvsMultRawEvSel);
    fOutput->Add(fHistNtrEta16vsNtrEta1EvWithCand);
    fOutput->Add(fHistNtrEta05vsNtrEta1EvWithCand);
    fOutput->Add(fHistNtrEta03vsNtrEta1EvWithCand);
    fOutput->Add(fHistNtrEtaV0AvsNtrEta1EvWithCand);
    fOutput->Add(fHistNtrEtaV0MvsNtrEta1EvWithCand);
    fOutput->Add(fHistNtrEtaV0AvsV0AEqEvWithCand);
    fOutput->Add(fHistNtrEtaV0MvsV0MEqEvWithCand);
    fOutput->Add(fHistNtrCorrEta1vsNtrRawEta1EvWithCand);
    fOutput->Add(fHistMultCorrvsMultRawEvWithCand);
    fOutput->Add(fHistNtrEta16vsNtrEta1EvWithD);
    fOutput->Add(fHistNtrEta05vsNtrEta1EvWithD);
    fOutput->Add(fHistNtrEta03vsNtrEta1EvWithD);
    fOutput->Add(fHistNtrEtaV0AvsNtrEta1EvWithD);
    fOutput->Add(fHistNtrEtaV0MvsNtrEta1EvWithD);
    fOutput->Add(fHistNtrEtaV0AvsV0AEqEvWithD);
    fOutput->Add(fHistNtrEtaV0MvsV0MEqEvWithD);
    fOutput->Add(fHistNtrCorrEta1vsNtrRawEta1EvWithD);
    fOutput->Add(fHistMultCorrvsMultRawEvWithD);
  }
  fOutput->Add(fHistNtrVsZvtx);
  fOutput->Add(fHistNtrCorrVsZvtx);

  fOutput->Add(fHistNtrVsNchMC);
  fOutput->Add(fHistNtrCorrVsNchMC);
  fOutput->Add(fHistNtrVsNchMCPrimary);
  fOutput->Add(fHistNtrCorrVsNchMCPrimary);
  fOutput->Add(fHistNtrVsNchMCPhysicalPrimary);
  fOutput->Add(fHistNtrCorrVsNchMCPhysicalPrimary);
  fOutput->Add(fHistGenPrimaryParticlesInelGt0);
  fOutput->Add(fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary);

  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",11,-0.5,10.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents total");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents with Z vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents selected");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to phys sel");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to vertex cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to pileup");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"D after cuts (No PID)");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"D after cuts + PID)"); 
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fPtVsMassVsMult=new TH3F("hPtVsMassvsMult", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);
 
  fPtVsMassVsMultNoPid=new TH3F("hPtVsMassvsMultNoPid", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.); 

  fPtVsMassVsMultUncorr=new TH3F("hPtVsMassvsMultUncorr", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fPtVsMassVsMultPart=new TH3F("hPtVsMassvsMultPart", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fPtVsMassVsMultAntiPart=new TH3F("hPtVsMassvsMultAntiPart", "D candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fPtVsMassVsMultMC=new TH3F("hPtVsMassvsMultMC", "D true candidates: p_{t} vs mass vs tracklets multiplicity; Tracklets; Mass M [GeV/c^{2}]; p_{t} [GeV/c]",nMultBins,firstMultBin,lastMultBin,fNMassBins,fLowmasslimit,fUpmasslimit,48,0.,24.);

  fOutput->Add(fPtVsMassVsMult);
  fOutput->Add(fPtVsMassVsMultUncorr);
  fOutput->Add(fPtVsMassVsMultNoPid);
  fOutput->Add(fPtVsMassVsMultPart);
  fOutput->Add(fPtVsMassVsMultAntiPart);
  fOutput->Add(fPtVsMassVsMultMC);

  if(fDoImpPar) CreateImpactParameterHistos();
  
  //Counter for Normalization 
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
 
  fCounterC = new AliNormalizationCounter("NormCounterCorrMult");
  fCounterC->SetStudyMultiplicity(kTRUE,1.);
  fCounterC->Init();

  fCounterU = new AliNormalizationCounter("NormCounterUnCorrMult");
  fCounterU->SetStudyMultiplicity(kTRUE,1.);
  fCounterU->Init(); 
    
  fCounterCandidates = new AliNormalizationCounter("NormCounterCorrMultCandidates");
  fCounterCandidates->SetStudyMultiplicity(kTRUE,1.);
  fCounterCandidates->Init();

    
  fOutputCounters = new TList();
  fOutputCounters->SetOwner();
  fOutputCounters->SetName("OutputCounters");
  fOutputCounters->Add(fCounter);
  fOutputCounters->Add(fCounterC);
  fOutputCounters->Add(fCounterU);
  fOutputCounters->Add(fCounterCandidates);

  // BDT I/O
  fListBDTNtuple = new TList(); fListBDTNtuple->SetOwner(); fListBDTNtuple->SetName("NtupleList");
  if(!fReadMC && fApplyBDT){      
     fListRDHFBDT->SetOwner(); fListBDTNtuple->SetName("BDTList");
     fListBDTResp = new TList(); fListBDTResp->SetOwner(); fListBDTResp->SetName("BDTResponseList");
  }
  if(fReadMC){
     TNtuple *NtupleMC = new TNtuple("NtupleMC", "MC tree", fBDTFullVarString);
     fListBDTNtuple->Add(NtupleMC);
  }
  else{
    if(!fReadMC && fApplyBDT){ 
       Float_t *ptbin = fRDCutsAnalysis->GetPtBinLimits();
                Float_t BDT1(0);
                Float_t BDT2[6]={0,0,0,0,0,0};
                for(Int_t ii=0;ii<fRDCutsAnalysis->GetNPtBins();ii++){
                        TString ptstring = Form("_%.0f_%.0f",ptbin[ii],ptbin[ii+1]);
                        if(fSampleSideband){ // Deal with the sideband sampling, only need 1st step BDT
                          TNtuple *NtupleSB = new TNtuple(Form("NtupleSB%s",ptstring.Data()), "Sideband", fBDTFullVarString);
                                TTree *BDTRespTree = new TTree(Form("BDTRespTree%s",ptstring.Data()),"BDT2 Response");
                                BDTRespTree->Branch("BDT1",&BDT1,"BDT1Resp");
                                fListBDTNtuple->Add(NtupleSB);
                                fListBDTResp->Add(BDTRespTree);
                        }    
                        else{ // Deal with the real data stuff, need both steps BDT
                          if(fGetRespTree){
                                  TNtuple *NtupleData = new TNtuple(Form("NtupleData%s",ptstring.Data()), " Data", "mass");
                                  TTree *BDTRespTree = new TTree(Form("BDTRespTree%s",ptstring.Data()),"BDT2 Response");
                                  BDTRespTree->Branch("BDT1",&BDT1,"BDT1Resp");
                                  BDTRespTree->Branch("BDT2",&BDT2,"BDT2RespLL:BDT2RespL:BDT2RespML:BDT2RespMH:BDT2RespH:BDT2RespHH");
                                  fListBDTNtuple->Add(NtupleData);
                                  fListBDTResp->Add(BDTRespTree);
                              }   
                              else{
                                  TH3F *h3InvmassBDTResp[6];
                                  h3InvmassBDTResp[0] = new TH3F(Form("h3MassBDTResp%s_1",ptstring.Data()),"Invmass vs BDT1Resp vs BDT2Resp pt",100,1.68,2.10,70,-0.05,0.30,70,-0.05,0.30);     
                                  h3InvmassBDTResp[1] = (TH3F*)h3InvmassBDTResp[0]->Clone(Form("h3MassBDTResp%s_2",ptstring.Data()));
                                  h3InvmassBDTResp[2] = (TH3F*)h3InvmassBDTResp[0]->Clone(Form("h3MassBDTResp%s_3",ptstring.Data()));
                                  h3InvmassBDTResp[3] = (TH3F*)h3InvmassBDTResp[0]->Clone(Form("h3MassBDTResp%s_4",ptstring.Data()));
                                  h3InvmassBDTResp[4] = (TH3F*)h3InvmassBDTResp[0]->Clone(Form("h3MassBDTResp%s_5",ptstring.Data()));
                                  h3InvmassBDTResp[5] = (TH3F*)h3InvmassBDTResp[0]->Clone(Form("h3MassBDTResp%s_6",ptstring.Data()));
                                  fListBDTResp->Add(h3InvmassBDTResp[0]);
                                  fListBDTResp->Add(h3InvmassBDTResp[1]);
                                  fListBDTResp->Add(h3InvmassBDTResp[2]);
                                  fListBDTResp->Add(h3InvmassBDTResp[3]);
                                  fListBDTResp->Add(h3InvmassBDTResp[4]);
                                  fListBDTResp->Add(h3InvmassBDTResp[5]);
                              }   
                                  
                       }          
           }                      
      }                           
  }                           
  
  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fOutputCounters);
  PostData(4,fListProfiles);
  PostData(5,fListBDTNtuple);
  if(!fReadMC && fApplyBDT){
     PostData(6,fListBDTResp);
  }
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event:
  /// heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  
  //  AliAODTracklets* tracklets = aod->GetTracklets();
  //Int_t ntracklets = tracklets->GetNumberOfTracklets();
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return;
    }
  }

  
  TClonesArray *arrayCand = 0;
  TString arrayName="";
  UInt_t pdgDau[3];
  Int_t nDau=0;
  Int_t selbit=0;
  if(fPdgMeson==411){
    arrayName="Charm3Prong";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=211; 
    nDau=3;
    selbit=AliRDHFCuts::kDplusCuts;
  }else if(fPdgMeson==421){
    arrayName="D0toKpi";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=0;
    nDau=2;
    selbit=AliRDHFCuts::kD0toKpiCuts;
  }else if(fPdgMeson==413){
    arrayName="Dstar";
    pdgDau[0]=321; pdgDau[1]=211; pdgDau[2]=0; // Quoting here D0 daughters (D* ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kDstarCuts;
  }else if(fPdgMeson==431){
    arrayName="Charm3Prong";
    pdgDau[0]=321; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kDsCuts;
  }else if(fPdgMeson==4122){
    if(fLctoV0){
    arrayName="CascadesHF";
    pdgDau[0]=211; pdgDau[1]=211; pdgDau[2]=0; // Quoting here K0S daughters (Lc ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kLctoV0Cuts;
    }else{
    arrayName="Charm3Prong";
    pdgDau[0]=2212; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kLcCuts;
    }
  }

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
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayCand=(TClonesArray*)aodFromExt->GetList()->FindObject(arrayName.Data());
    }
  } else if(aod) {
    arrayCand=(TClonesArray*)aod->GetList()->FindObject(arrayName.Data());
  }

  if(!aod || !arrayCand) {
    printf("AliAnalysisTaskSEDvsMultiplicity_BDT::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  if(fisPPbData && fReadMC){
    Int_t runnumber = aod->GetRunNumber();
    if(aod->GetTriggerMask()==0 && 
       (runnumber>=195344 && runnumber<=195677)){
      AliDebug(3,"Event rejected because of null trigger mask");
      return;
    }
  }


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
  // Int_t countTreta1=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.);
  // Int_t countTreta03=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-0.3,0.3);
  // Int_t countTreta05=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-0.5,0.5);
  // Int_t countTreta16=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.6,1.6);
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


  fCounterU->StoreEvent(aod,fRDCutsAnalysis,fReadMC,countMult);
  fHistNEvents->Fill(0); // count event

  Double_t countTreta1corr=countTreta1;
  Double_t countCorr=countMult;
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  // In case of VZERO multiplicity, consider the zvtx correction flag
  //  fDoVZER0ParamVertexCorr: 0= none, 1= usual d2h, 2=AliESDUtils
  Bool_t isDataDrivenZvtxCorr=kTRUE;
  Bool_t isVtxOk=kFALSE;
  Int_t vzeroMultACorr=vzeroMultA, vzeroMultCCorr=vzeroMultC, vzeroMultCorr=vzeroMult;
  Int_t vzeroMultAEqCorr=vzeroMultAEq, vzeroMultCEqCorr=vzeroMultCEq, vzeroMultEqCorr=vzeroMultEq;
  if(vtx1){
    if(vtx1->GetNContributors()>0){
      fHistNEvents->Fill(1);
      isVtxOk=kTRUE;
    }
  }
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


  fCounterC->StoreEvent(aod,fRDCutsAnalysis,fReadMC,countCorr);

  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);

  if(fRDCutsAnalysis->GetWhyRejection()==5) fHistNEvents->Fill(3);
  if(fRDCutsAnalysis->GetWhyRejection()==7) fHistNEvents->Fill(4);
  if(fRDCutsAnalysis->GetWhyRejection()==6) fHistNEvents->Fill(5);
  if(fRDCutsAnalysis->GetWhyRejection()==1) fHistNEvents->Fill(6);

  Bool_t isEvPSRejected = fRDCutsAnalysis->IsEventRejectedDuePhysicsSelection();
  Bool_t isEvTrigNameRejected = fRDCutsAnalysis->IsEventRejectedDueToTrigger();
  Bool_t isEvPileUpRejected = fRDCutsAnalysis->IsEventRejectedDueToPileup();
  Bool_t isEvNoVtxRejected = fRDCutsAnalysis->IsEventRejectedDueToNotRecoVertex();
  Bool_t isEvVtxContribRejected = fRDCutsAnalysis->IsEventRejectedDueToVertexContributors();
  Bool_t isEvVtxRangeRejected= fRDCutsAnalysis->IsEventRejectedDueToZVertexOutsideFiducialRegion();
  Bool_t isEvCentralityRejected = fRDCutsAnalysis->IsEventRejectedDueToCentrality();
  if(!isEvPSRejected){
    fHistNtrUnCorrPSSel->Fill(countMult);
    fHistNtrCorrPSSel->Fill(countCorr);
    if(!isEvTrigNameRejected){
      fHistNtrUnCorrPSTrigSel->Fill(countMult);
      if(!isEvPileUpRejected){
	fHistNtrUnCorrPSTrigPileUpSel->Fill(countMult);
	if(!isEvNoVtxRejected){
	  fHistNtrUnCorrPSTrigPileUpVtxSel->Fill(countMult);
	  if(!isEvVtxContribRejected){
	    fHistNtrUnCorrPSTrigPileUpVtxContSel->Fill(countMult);
	    if(!isEvVtxRangeRejected){
	      fHistNtrUnCorrPSTrigPileUpVtxRangeSel->Fill(countMult);
	      if(!isEvCentralityRejected){
		fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel->Fill(countMult);
	      }
	    }
	  }
	}
      }
    }
  }
  
  if(!isEvSel)return;
  fHistNEvents->Fill(2);
  if(fKeepCorrPlots){
    fHistNtrEta16vsNtrEta1EvSel->Fill(countTreta1,countTreta16);
    fHistNtrEta05vsNtrEta1EvSel->Fill(countTreta1,countTreta05);
    fHistNtrEta03vsNtrEta1EvSel->Fill(countTreta1,countTreta03);
    fHistNtrEtaV0AvsNtrEta1EvSel->Fill(countTreta1,vzeroMultA);
    fHistNtrEtaV0MvsNtrEta1EvSel->Fill(countTreta1,vzeroMult);
    fHistNtrEtaV0AvsV0AEqEvSel->Fill(vzeroMultA,vzeroMultAEq);
    fHistNtrEtaV0MvsV0MEqEvSel->Fill(vzeroMult,vzeroMultEq);
    fHistNtrCorrEta1vsNtrRawEta1EvSel->Fill(countTreta1,countTreta1corr);
    fHistMultCorrvsMultRawEvSel->Fill(countMult,countCorr);
  }
  if(vtx1){
    fHistNtrVsZvtx->Fill(vtx1->GetZ(),countMult);
    fHistNtrCorrVsZvtx->Fill(vtx1->GetZ(),countCorr);
  }

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  Double_t nchWeight=1.0;

  // load MC particles
  if(fReadMC){
     
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDvsMultiplicity_BDT::UserExec: MC particles branch not found!\n");
      return;
    }  
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDvsMultiplicity_BDT::UserExec: MC header branch not found!\n");
      return;
    }
  

    // Int_t nChargedMC=AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-1.0,1.0);
    // Int_t nChargedMCPrimary=AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-1.0,1.0);
    // Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-1.0,1.0);
    //
    Int_t nChargedMCEta10=0, nChargedMCEta03=0, nChargedMCEta05=0, nChargedMCEta16=0, nChargedMCEtam37tm17=0, nChargedMCEta28t51=0;
    Int_t nChargedMCPrimaryEta10=0, nChargedMCPrimaryEta03=0, nChargedMCPrimaryEta05=0, nChargedMCPrimaryEta16=0, nChargedMCPrimaryEtam37tm17=0, nChargedMCPrimaryEta28t51=0;
    Int_t nChargedMCPhysicalPrimaryEta10=0, nChargedMCPhysicalPrimaryEta03=0, nChargedMCPhysicalPrimaryEta05=0, nChargedMCPhysicalPrimaryEta16=0, nChargedMCPhysicalPrimaryEtam37tm17=0, nChargedMCPhysicalPrimaryEta28t51=0;
    for(Int_t i=0; i<arrayMC->GetEntriesFast(); i++){
      AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(i);
      Int_t charge = part->Charge();
      Double_t eta = part->Eta();
      Bool_t isPrim = part->IsPrimary();
      Bool_t isPhysPrim = part->IsPhysicalPrimary();
      if(charge!=0) {
	if(eta>-0.3 && eta< 0.3) { 
	  nChargedMCEta03++;
	  if(isPrim) nChargedMCPrimaryEta03++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEta03++;
	}
	if(eta>-0.5 && eta< 0.5) { 
	  nChargedMCEta05++;
	  if(isPrim) nChargedMCPrimaryEta05++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEta05++;
	}
	if(eta>-1.0 && eta< 1.0) { 
	  nChargedMCEta10++;
	  if(isPrim) nChargedMCPrimaryEta10++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEta10++;
	}
	if(eta>-1.6 && eta< 1.6) { 
	  nChargedMCEta16++;
	  if(isPrim) nChargedMCPrimaryEta16++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEta16++;
	}
	if(eta>-3.7 && eta<-1.7) { 
	  nChargedMCEtam37tm17++;
	  if(isPrim) nChargedMCPrimaryEtam37tm17++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEtam37tm17++;
	}
	if(eta> 2.8 && eta< 5.1) { 
	  nChargedMCEta28t51++;
	  if(isPrim) nChargedMCPrimaryEta28t51++;
	  if(isPhysPrim) nChargedMCPhysicalPrimaryEta28t51++;
	}
      }
    }
    Int_t nChargedMC=nChargedMCEta10;
    Int_t nChargedMCPrimary=nChargedMCPrimaryEta10;
    Int_t nChargedMCPhysicalPrimary=nChargedMCPhysicalPrimaryEta10;

    // Compute the Nch weights (reference is Ntracklets within |eta|<1.0)
    if(fUseNchWeight>0){
      
      Double_t tmpweight = 1.0;
      Double_t tmpXweight=nChargedMCPhysicalPrimary; // Nch weights
      if(fUseNchWeight==2) tmpXweight=countMult;     // Ntrk weights

      if(tmpXweight<=0) tmpweight = 0.0;
      else{
	Double_t pMeas = fHistoMeasNch->GetBinContent(fHistoMeasNch->FindBin(tmpXweight));
	//printf(" pMeas=%2.2f  and histo MCNch %s \n",pMeas,fHistoMCNch);
	Double_t pMC = fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(tmpXweight));
	tmpweight = pMC>0 ? pMeas/pMC : 0.;
      }
      nchWeight *= tmpweight;
      AliDebug(2,Form("Using Nch weights, Mult=%f Weight=%f\n",tmpXweight,nchWeight));
      
    }

    // Now recompute the variables in case another MC estimator is considered
    // Int_t nChargedMCEta10 = nChargedMC;
    // Int_t nChargedMCEta16 = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-1.6,1.6);
    // Int_t nChargedMCEta05 = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-0.5,0.5);
    // Int_t nChargedMCEta03 = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-0.3,0.3);
    // Int_t nChargedMCEtam37tm17 = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-3.7,-1.7);
    // Int_t nChargedMCEta28t51 = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,2.8,5.1);
    // Int_t nChargedMCPrimaryEta10 = nChargedMCPrimary;
    // Int_t nChargedMCPrimaryEta16 = AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-1.6,1.6);
    // Int_t nChargedMCPrimaryEta05 = AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-0.5,0.5);
    // Int_t nChargedMCPrimaryEta03 = AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-0.3,0.3);
    // Int_t nChargedMCPrimaryEtam37tm17 = AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,-3.7,-1.7);
    // Int_t nChargedMCPrimaryEta28t51 = AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(arrayMC,2.8,5.1);
    // Int_t nChargedMCPhysicalPrimaryEta10 = nChargedMCPhysicalPrimary;
    // Int_t nChargedMCPhysicalPrimaryEta16 = AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-1.6,1.6);
    // Int_t nChargedMCPhysicalPrimaryEta05 = AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-0.5,0.5);
    // Int_t nChargedMCPhysicalPrimaryEta03 = AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-0.3,0.3);
    // Int_t nChargedMCPhysicalPrimaryEtam37tm17 = AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,-3.7,-1.7);
    // Int_t nChargedMCPhysicalPrimaryEta28t51 = AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(arrayMC,2.8,5.1);
    if(fMCPrimariesEstimator==kEta10to16){
      nChargedMC = nChargedMCEta16 - nChargedMCEta10;
      nChargedMCPrimary = nChargedMCPrimaryEta16 - nChargedMCPrimaryEta10;
      nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta16 - nChargedMCPhysicalPrimaryEta10;
    } else if(fMCPrimariesEstimator==kEta05){
      nChargedMC = nChargedMCEta05;
      nChargedMCPrimary = nChargedMCPrimaryEta05;
      nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta05;
    } else if(fMCPrimariesEstimator==kEta03){
      nChargedMC = nChargedMCEta03;
      nChargedMCPrimary = nChargedMCPrimaryEta03;
      nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta03;
    } else if(fMCPrimariesEstimator==kEtaVZERO){
      nChargedMC = nChargedMCEtam37tm17 + nChargedMCEta28t51;
      nChargedMCPrimary = nChargedMCPrimaryEtam37tm17 + nChargedMCPrimaryEta28t51;
      nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEtam37tm17 + nChargedMCPhysicalPrimaryEta28t51;
    } else if(fMCPrimariesEstimator==kEtaVZEROA){
      nChargedMC = nChargedMCEta28t51;
      nChargedMCPrimary = nChargedMCPrimaryEta28t51;
      nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta28t51;
    }

    // Here fill the MC correlation plots
    if(nChargedMCPhysicalPrimary>0){ // INEL>0 for |eta|<1
      fHistGenPrimaryParticlesInelGt0->Fill(nChargedMCPhysicalPrimary,nchWeight);
    }

    fHistNtrVsNchMC->Fill(nChargedMC,countMult,nchWeight);
    fHistNtrCorrVsNchMC->Fill(nChargedMC,countCorr,nchWeight);

    fHistNtrVsNchMCPrimary->Fill(nChargedMCPrimary,countMult,nchWeight);
    fHistNtrCorrVsNchMCPrimary->Fill(nChargedMCPrimary,countCorr,nchWeight);

    fHistNtrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countMult,nchWeight);
    fHistNtrCorrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countCorr,nchWeight);

    fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary->Fill(nChargedMC,nChargedMCPrimary,nChargedMCPhysicalPrimary,nchWeight);
  }
  
  Int_t nCand = arrayCand->GetEntriesFast(); 
  Int_t nSelectedNoPID=0,nSelectedPID=0,nSelectedInMassPeak=0;
  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t mDsPDG    = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Double_t mLcPDG    = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  // pdg of daughters needed for D* too
  UInt_t pdgDgDStartoD0pi[2]={421,211};

  // pdg of daughters needed for Lc2pK0S too
  UInt_t pdgDgLctopK0S[2]={2212,310};

  Double_t aveMult=0.;
  Double_t nSelCand=0.;
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayCand->UncheckedAt(iCand);
    AliAODRecoCascadeHF *dCascade = NULL;
    if(fPdgMeson==413 || (fPdgMeson==4122 && fLctoV0)) dCascade = (AliAODRecoCascadeHF*)d;
  
    fHistNEvents->Fill(7);

    if(fPdgMeson==4122){
      if(fLctoV0){
      AliAODv0 * v0part = (AliAODv0*)dCascade->Getv0();
      Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s
      if (onFlyV0){fHistNEvents->Fill(8); continue;}
      }else if(!fLctoV0 && fUseBit && !d->HasSelectionBit(selbit)) {
      fHistNEvents->Fill(8);
      continue;
    }
    }

    if(fPdgMeson!=4122 && fUseBit && !d->HasSelectionBit(selbit)) {
      fHistNEvents->Fill(8);
      continue;
    }
    
    Double_t ptCand = d->Pt();
    Double_t rapid=d->Y(fPdgMeson);
    Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
    if(!isFidAcc) continue;
   
    Int_t labD=-1;
    if(fReadMC) {
      if(fPdgMeson==413){
	labD = dCascade->MatchToMC(fPdgMeson,421,(Int_t*)pdgDgDStartoD0pi,(Int_t*)pdgDau,arrayMC);
      } else if(fPdgMeson==4122 && fLctoV0){
	labD = dCascade->MatchToMC(fPdgMeson,pdgDgLctopK0S[1],(Int_t*)pdgDgLctopK0S,(Int_t*)pdgDau,arrayMC,kTRUE);
      } else {
	labD = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
      }
      FillMCMassHistos(arrayMC,labD, countMult,nchWeight);
    }

    Int_t passAllCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
    if (fPdgMeson == 4122 && fLctoV0) passAllCuts=(((fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr));
    Int_t passTopolCuts=fRDCutsAnalysis->GetIsSelectedCuts();
    if (fPdgMeson == 4122){
      if(fLctoV0) passTopolCuts=(((fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr));
      else  passTopolCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    }
    if (fPdgMeson != 431 && passTopolCuts==0) continue;
    nSelectedNoPID++;
    fHistNEvents->Fill(9);
    if(fPdgMeson == 431 && passAllCuts==0) continue;
    if(passAllCuts){
      nSelectedPID++;
      fHistNEvents->Fill(10);
    }
    Double_t multForCand = countCorr;

    if(fSubtractTrackletsFromDau){
      if(fPdgMeson==4122 && fLctoV0) {
	// For the Lc2pK0S case, subtract only the proton daughter track <=== FIXME !!
	AliAODTrack *t = (AliAODTrack *)dCascade->GetBachelor();
	if(!t) continue;
	if(t->HasPointOnITSLayer(0) && t->HasPointOnITSLayer(1)){
	  if(multForCand>0) multForCand-=1;
	}
      } else {
	// For the D* case, subtract only the D0 daughter tracks <=== FIXME !!
	AliAODRecoDecayHF2Prong* d0fromDstar = NULL;
	if(fPdgMeson==413) d0fromDstar = (AliAODRecoDecayHF2Prong*)dCascade->Get2Prong();

	for(Int_t iDau=0; iDau<nDau; iDau++){
	  AliAODTrack *t = NULL;
	  if(fPdgMeson==413){ t = (AliAODTrack*)d0fromDstar->GetDaughter(iDau); }
	  else{ t = (AliAODTrack*)d->GetDaughter(iDau); }
	  if(!t) continue;
	  if(t->HasPointOnITSLayer(0) && t->HasPointOnITSLayer(1)){
	    if(multForCand>0) multForCand-=1;
	  }
	}
      }
    }
    Bool_t isPrimary=kTRUE;
    Double_t trueImpParXY=9999.;
    Double_t impparXY=d->ImpParXY()*10000.;
    Double_t dlen=0.1; //FIXME
    Double_t mass[2];
    if(fPdgMeson==411){
      mass[0]=d->InvMass(nDau,pdgDau);
      mass[1]=-1.;
      if(TMath::Abs(mass[0]-mDplusPDG)<0.02) nSelectedInMassPeak++; //20 MeV for now... FIXME
    }else if(fPdgMeson==421){
      UInt_t pdgdaughtersD0[2]={211,321};//pi,K 
      UInt_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      mass[0]=d->InvMass(2,pdgdaughtersD0);
      mass[1]=d->InvMass(2,pdgdaughtersD0bar);
      if(TMath::Abs(mass[0]-mD0PDG)<0.02 || TMath::Abs(mass[1]-mD0PDG)<0.02 ) nSelectedInMassPeak++; //20 MeV for now... FIXME
    }else if(fPdgMeson==413){
      // FIXME
      mass[0]=dCascade->DeltaInvMass();
      mass[1]=-1.;
      if(TMath::Abs(mass[0]-(mDstarPDG-mD0PDG))<0.0015) nSelectedInMassPeak++; //1 MeV for now... FIXME
    }else if(fPdgMeson==431){
      UInt_t pdgDsKKpi[3]={321,321,211};
      UInt_t pdgDspiKK[3]={211,321,321};
      mass[0]=d->InvMass(nDau,pdgDsKKpi);
      mass[1]=d->InvMass(nDau,pdgDspiKK);
      if(TMath::Abs(mass[0]-mDsPDG)<0.02 || TMath::Abs(mass[1]-mDsPDG)<0.02 ) nSelectedInMassPeak++; //20 MeV for now... FIXME
    }else if(fPdgMeson==4122){
      if(fLctoV0){
      mass[0]=d->InvMass(2,pdgDgLctopK0S);
      mass[1]=-1.;
      if(TMath::Abs(mass[0]-mLcPDG)<0.02) nSelectedInMassPeak++; //20 MeV for now... FIXME
     }else{
      UInt_t pdgpKpi[3]={2212,321,211};
      UInt_t pdgpiKp[3]={211,321,2212}; 
      if(passTopolCuts==3 || passTopolCuts==1)mass[0]=d->InvMass(3,pdgpKpi);
      if(passTopolCuts>=2) mass[1]=d->InvMass(3,pdgpiKp);
      if(TMath::Abs(mass[0]-mLcPDG)<0.02) nSelectedInMassPeak++; //20 MeV for now... FIXME
     }
    }

    if(fFillTree && fPdgMeson==421) ProcessBDTD0(aod,d,fRDCutsAnalysis,aod->GetMagneticField(),arrayMC,passAllCuts,multForCand);    
    if(fFillTree && fPdgMeson==431) ProcessBDTDs(d,fRDCutsAnalysis,aod->GetMagneticField(),arrayMC,passAllCuts,multForCand);

     for(Int_t iHyp=0; iHyp<2; iHyp++){
      if(mass[iHyp]<0.) continue; // for D+,D* and Lc2pK0S we have 1 mass hypothesis
      Double_t invMass=mass[iHyp];
      Double_t arrayForSparse[5]={invMass,ptCand,impparXY,dlen,multForCand};

      if(fReadMC){

 	if(fPdgMeson==413){
	  labD = dCascade->MatchToMC(fPdgMeson,421,(Int_t*)pdgDgDStartoD0pi,(Int_t*)pdgDau,arrayMC);
	}else if(fPdgMeson==4122 && fLctoV0){
	  labD = dCascade->MatchToMC(fPdgMeson,pdgDgLctopK0S[1],(Int_t*)pdgDgLctopK0S,(Int_t*)pdgDau,arrayMC,kTRUE);
  }else {
	  labD = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
	}

	Bool_t fillHisto=fDoImpPar;
	if(labD>=0){
	  AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
	  Int_t code=partD->GetPdgCode();
	  if(CheckOrigin(arrayMC,partD)==5) isPrimary=kFALSE;
	  if(code<0 && iHyp==0) fillHisto=kFALSE;
	  if(code>0 && iHyp==1) fillHisto=kFALSE;
	  if(!isPrimary){
	    if(fPdgMeson==411){
	      trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDplus(mcHeader,arrayMC,partD)*10000.;
	    }else if(fPdgMeson==421){
	      trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDzero(mcHeader,arrayMC,partD)*10000.;
	    }else if(fPdgMeson==413 || fPdgMeson==4122){
	      trueImpParXY=0.; /// FIXME
	    }
	    Double_t arrayForSparseTrue[5]={invMass,ptCand,trueImpParXY,dlen,multForCand};
	    if(fillHisto && passAllCuts){
	      fHistMassPtImpPar[2]->Fill(arrayForSparse);
	      fHistMassPtImpPar[3]->Fill(arrayForSparseTrue);
	    }
	  }else{
	    if(fillHisto && passAllCuts) fHistMassPtImpPar[1]->Fill(arrayForSparse);
	  }
	}else{
	  if(fillHisto && passAllCuts)fHistMassPtImpPar[4]->Fill(arrayForSparse);
	}

	if(TMath::Abs(labD)==fPdgMeson && fMCOption==2) continue;
	if(TMath::Abs(labD)!=fPdgMeson && fMCOption==1) continue;      

      }
      
      if(fPdgMeson==421){
	if(iHyp==0 && !(passTopolCuts&1)) continue; // candidate not passing as D0
	if(iHyp==1 && !(passTopolCuts&2)) continue; // candidate not passing as D0bar
      }

      if(fPdgMeson != 431) fPtVsMassVsMultNoPid->Fill(multForCand,invMass,ptCand);

      if(fPdgMeson==421){
	if(iHyp==0 && !(passAllCuts&1)) continue; // candidate not passing as D0
	if(iHyp==1 && !(passAllCuts&2)) continue; // candidate not passing as D0bar
      }
      if(fPdgMeson==431){
	if(iHyp==0 && !(passAllCuts&4)) continue; // candidates Ds not passing as kk(phi)pi
	if(iHyp==1 && !(passAllCuts&8)) continue; // candidates Ds not passing as pikk(phi)
      }
      if(fPdgMeson==4122 && !fLctoV0){
     if(iHyp==0 && !(passTopolCuts&1)) continue; // candidate not passing as Lc-> pKpi
     if(iHyp==1 && !(passTopolCuts&2)) continue; // candidate not passing as Lc-> piKp
      }
        
      if(passAllCuts){
	aveMult+=multForCand;
	nSelCand+=1.;
	fPtVsMassVsMult->Fill(multForCand,invMass,ptCand,nchWeight);
	fPtVsMassVsMultUncorr->Fill(countTreta1,invMass,ptCand,nchWeight);
	// Add separation between part antipart
	if(fPdgMeson==411 || fPdgMeson==413 || fPdgMeson==431 || fPdgMeson==4122){
	  if(d->GetCharge()>0) fPtVsMassVsMultPart->Fill(multForCand,invMass,ptCand,nchWeight);
	  else fPtVsMassVsMultAntiPart->Fill(multForCand,invMass,ptCand,nchWeight);
	}else if(fPdgMeson==421){
	  if(passAllCuts&1) fPtVsMassVsMultPart->Fill(multForCand,invMass,ptCand,nchWeight);
	  if(passAllCuts&2) fPtVsMassVsMultAntiPart->Fill(multForCand,invMass,ptCand,nchWeight);
	}
      	
	if(fDoImpPar){
	  fHistMassPtImpPar[0]->Fill(arrayForSparse);
	}
	
      }

    }
  }
  if(fSubtractTrackletsFromDau && nSelCand>0){
    aveMult/=nSelCand;
    fCounterCandidates->StoreEvent(aod,fRDCutsAnalysis,fReadMC,(Int_t)(aveMult+0.5001));
  }else{
    fCounterCandidates->StoreEvent(aod,fRDCutsAnalysis,fReadMC,(Int_t)countCorr);
  }

  fCounter->StoreCandidates(aod,nSelectedNoPID,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedPID,kFALSE);

  fCounterCandidates->StoreCandidates(aod,nSelectedNoPID,kTRUE);
  fCounterCandidates->StoreCandidates(aod,nSelectedPID,kFALSE);
  fHistNtrUnCorrEvSel->Fill(countMult,nchWeight);
  fHistNtrCorrEvSel->Fill(countCorr,nchWeight);
  if(nSelectedPID>0) {
    fHistNtrUnCorrEvWithCand->Fill(countMult,nchWeight);
    fHistNtrCorrEvWithCand->Fill(countCorr,nchWeight);
    if(fKeepCorrPlots){
      fHistNtrEta16vsNtrEta1EvWithCand->Fill(countTreta1,countTreta16);
      fHistNtrEta05vsNtrEta1EvWithCand->Fill(countTreta1,countTreta05);
      fHistNtrEta03vsNtrEta1EvWithCand->Fill(countTreta1,countTreta03);
      fHistNtrEtaV0AvsNtrEta1EvWithCand->Fill(countTreta1,vzeroMultA);
      fHistNtrEtaV0MvsNtrEta1EvWithCand->Fill(countTreta1,vzeroMult);
      fHistNtrEtaV0AvsV0AEqEvWithCand->Fill(vzeroMultA,vzeroMultAEq);
      fHistNtrEtaV0MvsV0MEqEvWithCand->Fill(vzeroMult,vzeroMultEq);
      fHistNtrCorrEta1vsNtrRawEta1EvWithCand->Fill(countTreta1,countTreta1corr);
      fHistMultCorrvsMultRawEvWithCand->Fill(countMult,countCorr);
    }
  }
  if(nSelectedInMassPeak>0) {
    fHistNtrUnCorrEvWithD->Fill(countMult,nchWeight);
    fHistNtrCorrEvWithD->Fill(countCorr,nchWeight);
    if(fKeepCorrPlots){
      fHistNtrEta16vsNtrEta1EvWithD->Fill(countTreta1,countTreta16);
      fHistNtrEta05vsNtrEta1EvWithD->Fill(countTreta1,countTreta05);
      fHistNtrEta03vsNtrEta1EvWithD->Fill(countTreta1,countTreta03);
      fHistNtrEtaV0AvsNtrEta1EvWithD->Fill(countTreta1,vzeroMultA);
      fHistNtrEtaV0MvsNtrEta1EvWithD->Fill(countTreta1,vzeroMult);
      fHistNtrEtaV0AvsV0AEqEvWithD->Fill(vzeroMultA,vzeroMultAEq);
      fHistNtrEtaV0MvsV0MEqEvWithD->Fill(vzeroMult,vzeroMultEq);
      fHistNtrCorrEta1vsNtrRawEta1EvWithD->Fill(countTreta1,countTreta1corr);
      fHistMultCorrvsMultRawEvWithD->Fill(countMult,countCorr);
    }
  }

  fDaughterTracks.Clear();

  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fOutputCounters);
  PostData(5,fListBDTNtuple);
  if(!fReadMC && fApplyBDT){
     PostData(6,fListBDTResp);
  }  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::CreateImpactParameterHistos(){
  /// Histos for impact paramter study
  /// mass . pt , impact parameter , decay length , multiplicity

  Int_t nbins[5]={fNMassBins,200,fNImpParBins,50,100};
  Double_t xmin[5]={fLowmasslimit,0.,fLowerImpPar,0.,0.};
  Double_t xmax[5]={fUpmasslimit,20.,fHigherImpPar,1.,100.};

  fHistMassPtImpPar[0]=new THnSparseF("hMassPtImpParAll",
				      "Mass vs. pt vs.imppar - All",
				      5,nbins,xmin,xmax);
  fHistMassPtImpPar[1]=new THnSparseF("hMassPtImpParPrompt",
				      "Mass vs. pt vs.imppar - promptD",
				      5,nbins,xmin,xmax);
  fHistMassPtImpPar[2]=new THnSparseF("hMassPtImpParBfeed",
				      "Mass vs. pt vs.imppar - DfromB",
				      5,nbins,xmin,xmax);
  fHistMassPtImpPar[3]=new THnSparseF("hMassPtImpParTrueBfeed",
				      "Mass vs. pt vs.true imppar -DfromB",
				      5,nbins,xmin,xmax);
  fHistMassPtImpPar[4]=new THnSparseF("hMassPtImpParBkg",
				      "Mass vs. pt vs.imppar - backgr.",
				      5,nbins,xmin,xmax);
  for(Int_t i=0; i<5;i++){
    fOutput->Add(fHistMassPtImpPar[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(!fHistNEvents){
    printf("ERROR: fHistNEvents not available\n");
    return;    
  }
  printf("Number of Analyzed Events = %d\n",(Int_t)fHistNEvents->GetBinContent(3));
 
  return;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDvsMultiplicity_BDT::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  /// checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  //  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      //      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isFromB) return 5;
  else return 4;
}



//____________________________________________________________________________
TProfile* AliAnalysisTaskSEDvsMultiplicity_BDT::GetEstimatorHistogram(const AliVEvent* event){
  /// Get Estimator Histogram from period event->GetRunNumber();
  ///
  /// If you select SPD tracklets in |eta|<1 you should use type == 1
  ///
    
  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;   // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  // pPb 2013: 0-LHC13b, 1-LHC13c
  // pPb 2016: 0-LHC16q: 265499->265525; 265309->265387, 1-LHC16q:265435, 2-LHC16q:265388->265427, LHC16t: 267163->267166
  if (fisPPbData) {
    if(fYearNumber==13) {
      if (runNo>195343 && runNo<195484) period = 0; //13b
      else if (runNo>195528 && runNo<195678) period = 1; //13c
      if (period < 0 || period > 1) return 0;
    }
    else if(fYearNumber==16) {
      if ((runNo>=265499 && runNo<=265525) || (runNo>=265309 && runNo<=265387)) period = 0;      //16q
      else if (runNo == 265435) period = 1;                //16q
      else if (runNo>=265388 && runNo<=265427) period = 2; //16q
      else if (runNo>=267163 && runNo<=267166) period = 3; //16t
      if (period < 0 || period > 3) return 0;
    }
  }
  else {
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
  }
  return fMultEstimatorAvg[period];
}
//__________________________________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::CreateMeasuredNchHisto(){
  /// creates historgam with measured multiplcity distribution in pp 7 TeV collisions (from Eur. Phys. J. C (2010) 68: 345–354)
  //
  // for Nch  > 70 the points were obtainedwith a double NBD distribution
  // TF1 *fit1 = new TF1("fit1","[0]*(TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])))*(TMath::Power(([2]/[1]),x))*(TMath::Power((1+([2]/[1])),-x-[1]))"); fit1->SetParameter(0,1.);// normalization constant
  // fit1->SetParameter(1,1.63); // k parameter
  // fit1->SetParameter(2,12.8); // mean multiplicity
  Double_t nchbins[82]={0.50,1.50,2.50,3.50,4.50,5.50,6.50,7.50,8.50,9.50,
			10.50,11.50,12.50,13.50,14.50,15.50,16.50,17.50,18.50,19.50,
			20.50,21.50,22.50,23.50,24.50,25.50,26.50,27.50,28.50,29.50,
			30.50,31.50,32.50,33.50,34.50,35.50,36.50,37.50,38.50,39.50,
			40.50,41.50,42.50,43.50,44.50,45.50,46.50,47.50,48.50,49.50,
			50.50,51.50,52.50,53.50,54.50,55.50,56.50,57.50,58.50,59.50,
			60.50,62.50,64.50,66.50,68.50,70.50,72.50,74.50,76.50,78.50,
			80.50,82.50,84.50,86.50,88.50,90.50,92.50,94.50,96.50,98.50, 
			100.50,102.50};
  Double_t pch[81]={0.062011,0.072943,0.070771,0.067245,0.062834,0.057383,0.051499,0.04591,0.041109,0.036954,
		    0.03359,0.030729,0.028539,0.026575,0.024653,0.0229,0.021325,0.019768,0.018561,0.017187,
		    0.01604,0.014836,0.013726,0.012576,0.011481,0.010393,0.009502,0.008776,0.008024,0.007452,
		    0.006851,0.006428,0.00594,0.005515,0.005102,0.00469,0.004162,0.003811,0.003389,0.003071,
		    0.002708,0.002422,0.002184,0.001968,0.00186,0.00165,0.001577,0.001387,0.001254,0.001118,
		    0.001037,0.000942,0.000823,0.000736,0.000654,0.000579,0.000512,0.00049,0.00045,0.000355,
		    0.000296,0.000265,0.000193,0.00016,0.000126,0.0000851, 0.0000676,0.0000537,0.0000426, 0.0000338,
		    0.0000268,0.0000213,0.0000166,0.0000133,0.0000106,0.00000837,0.00000662, 0.00000524,0.00000414, 0.00000327,
		    0.00000258};

  if(fHistoMeasNch) delete fHistoMeasNch;
  fHistoMeasNch=new TH1F("hMeaseNch","",81,nchbins);
  for(Int_t i=0; i<81; i++){
    fHistoMeasNch->SetBinContent(i+1,pch[i]);
    fHistoMeasNch->SetBinError(i+1,0.);
  }
}

//__________________________________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::FillMCMassHistos(TClonesArray *arrayMC, Int_t labD, Int_t countMult,Double_t nchWeight) 
{
  //
  /// Function to fill the true MC signal
  //
  
  if(labD>=0){
    AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
    Double_t mass = partD->M();
    Double_t pt = partD->Pt();
    fPtVsMassVsMultMC->Fill(countMult,mass,pt,nchWeight);
  }

}
//_____________________________________________________________________________________________
void AliAnalysisTaskSEDvsMultiplicity_BDT::ProcessBDTD0(AliAODEvent *aod, AliAODRecoDecayHF *partHF, AliRDHFCuts *CutsAnalysis, Double_t bfield, TClonesArray *arrayMC, Int_t passCutsValue, Int_t coutMulti)
{
   
    fDaughterTracks.AddAt((AliAODTrack*)partHF->GetDaughter(0), 0);
    fDaughterTracks.AddAt((AliAODTrack*)partHF->GetDaughter(1), 1);
    AliAODTrack *prong2 = (AliAODTrack*)fDaughterTracks.UncheckedAt(0);
    AliAODTrack *prong3 = (AliAODTrack*)fDaughterTracks.UncheckedAt(1); 

    AliAODRecoDecayHF2Prong* part = (AliAODRecoDecayHF2Prong*)partHF;
    Double_t normIP[2];
    Double_t d0Prong[2];
    Double_t ptProng[2];
    
    Double_t invmassD0 = part->InvMassD0(); Double_t invmassD0bar = part->InvMassD0bar();
    Double_t ptB;

    Double_t cosPointingAngle;
    Double_t cosThetaStarD0 = 99;
    Double_t cosThetaStarD0bar = 99;

    Double_t diffIP[2], errdiffIP[2];
    part->Getd0MeasMinusExpProng(0,bfield,diffIP[0],errdiffIP[0]);
    part->Getd0MeasMinusExpProng(1,bfield,diffIP[1],errdiffIP[1]);
    normIP[0] = diffIP[0]/errdiffIP[0];
    normIP[1] = diffIP[1]/errdiffIP[1];

    AliAODVertex *secvtx  = part->GetSecondaryVtx();
    AliAODVertex *primvtx = part->GetPrimaryVtx();
    Double_t err2decaylength = secvtx->Error2DistanceXYToVertex(primvtx);
    Double_t lxy = part->AliAODRecoDecay::DecayLengthXY(primvtx);
    Bool_t ispid = CutsAnalysis->GetIsUsePID(); 
    d0Prong[0] = part->Getd0Prong(0); 
    d0Prong[1] = part->Getd0Prong(1); //d0*d0 d0Prong[1]*d0Prong[0]
    cosPointingAngle = part->CosPointingAngle();
    ptProng[0] = part->Pt2Prong(0); ptProng[1]=part->Pt2Prong(1);
    cosThetaStarD0 = part->CosThetaStarD0();
    cosThetaStarD0bar = part->CosThetaStarD0bar(); 

    Float_t tmp[35];
    tmp[8] = -99;                                // Invariant Mass
    tmp[18] = -99;                               // ptB (if accessible)
    tmp[19] = 0;                                 // PDGCode
    tmp[20] = -99;                               // Rapidity YD0
    tmp[21] = -99;                               // Azimuthal Phi
    
    tmp[0] = part->Pt();                         // ptD
    tmp[1] = normIP[0];                          // Normalized d0N-<d0N> (topo1)
    tmp[2] = normIP[1];                          // Normalized d0P-<d0P> (topo2)
    tmp[3] = lxy;                                // Decay length
    tmp[4] = lxy/TMath::Sqrt(err2decaylength);   // Normalized decay length
    tmp[9] = d0Prong[0]*d0Prong[1];              // d0N*d0P
    tmp[10] = cosPointingAngle;                  // CosThetaPointing
    tmp[11] = part->GetDCA();                    // DCAtracks ,if (part->Pt() > 10) return; DCA part->GetDCA();
    tmp[12] = prong2->Pt();                      // pt1
    tmp[13] = prong3->Pt();                      // pt2
    tmp[14] = part->CosPointingAngleXY();        // CosThetaPointingXY
    tmp[15] = d0Prong[0];                        // d01
    tmp[16] = d0Prong[1];                        // d02
    tmp[20] = part->YD0();
    tmp[21] = part->Phi();
    tmp[34] = coutMulti;

    // PID and Cuts
    if(ispid)CutsAnalysis->SetUsePID(kFALSE);// if PID on, switch it off
    Int_t isCuts=CutsAnalysis->IsSelected(part,AliRDHFCuts::kAll,aod);
    if(ispid)CutsAnalysis->SetUsePID(kTRUE);//if PID was on, switch it on
    Int_t isPid= CutsAnalysis->IsSelectedPID(part);
    tmp[5] = (Double_t)isCuts;
    tmp[6] = (Double_t)isPid;  
 
    //PID variable (Nsigma PID for TOF and TPC , combine PID)    
    AliAODPidHF *Pid_HF = CutsAnalysis->GetPidHF();
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    AliPIDResponse *pidrespo = inputHandler->GetPIDResponse();
    //AliPIDResponse* pidrespo = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
   
    static const unsigned int knMaxProngs   = 2;  //D0
    static const unsigned int knMaxDet4Pid  = 2;  //0=>TPC; 1=>TOF
    static const unsigned int knMaxHypo4Pid = 2;
    Double_t sig[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
    Double_t PIDNsigmaVector[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
    AliPID::EParticleType parthypo[knMaxHypo4Pid] = {AliPID::kPion, AliPID::kKaon};
    Double_t sigComb[knMaxProngs][knMaxHypo4Pid]; 
 
    for(unsigned int iProng = 0; iProng < knMaxProngs; iProng++) {
        for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
            /*Double_t nSigmaTPC = -999.;
              Pid_HF->GetnSigmaTPC(prongtracks[iProng], parthypo[iPartHypo], nSigmaTPC);
              sig[iProng][0][iPartHypo] = nSigmaTPC;     
              double nSigmaTOF = -999.;
              Pid_HF->GetnSigmaTOF(prongtracks[iProng], parthypo[iPartHypo], nSigmaTOF);
              sig[iProng][1][iPartHypo] = nSigmaTOF;*/
            sig[iProng][0][iPartHypo] = -999.;
            AliAODTrack *track = (AliAODTrack*)part->GetDaughter(iProng);
            sig[iProng][0][iPartHypo] = pidrespo->NumberOfSigmasTPC(track,parthypo[iPartHypo]);
            sig[iProng][1][iPartHypo] = -999.;
            sig[iProng][1][iPartHypo] = pidrespo->NumberOfSigmasTOF(track,parthypo[iPartHypo]); 
           
            //Combine PID 
            sigComb[iProng][iPartHypo] = AliVertexingHFUtils::CombineNsigmaTPCTOF(sig[iProng][0][iPartHypo], sig[iProng][1][iPartHypo]);  
       } 
    }

    Int_t index_signalPID = 0; Int_t index_combinPID = 0;
    for(unsigned int iProng = 0; iProng <= 1; iProng++) {  
        for(int iDet = 0; iDet <=1; iDet++) {
            for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {   
                PIDNsigmaVector[iProng][iDet][iPartHypo] = sig[iProng][iDet][iPartHypo];
                tmp[22+index_signalPID] = PIDNsigmaVector[iProng][iDet][iPartHypo];
                ++index_signalPID;
            }
        }
        for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
            tmp[30+index_combinPID] = sigComb[iProng][iPartHypo];            
            ++index_combinPID;
        }
    }
 
    if(fReadMC){  // MC THnSparse
        Int_t pdgDgD0toKpi[2] = {321,211};
        Int_t lab = -9999;
        lab = part->MatchToMC(421,arrayMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
        if(lab>=0){
            AliAODMCParticle *partD0 = (AliAODMCParticle*)arrayMC->At(lab);
            // PDGCode
            Int_t pdgcode = partD0->GetPdgCode();
            tmp[19] = (Float_t)pdgcode;
            if (pdgcode == 421) {tmp[17] = cosThetaStarD0; tmp[8] = invmassD0;}
            else if (pdgcode == -421) {tmp[17] = cosThetaStarD0bar; tmp[8] = invmassD0bar;}
     
            // Check feed down or prompt
            Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partD0,fUseQuarkTagInKine);
            if(orig==4)           tmp[7] = 4; //prompt
            else if(orig==5)       tmp[7] = 5; //feed down
            else tmp[7] = 0;

            // pT of Bmother if Accessible
            ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,partD0);
            tmp[18] = ptB;
            TNtuple *NtupleMC = (TNtuple*)fListBDTNtuple->FindObject("NtupleMC");
            if(tmp[7]==4)          NtupleMC->Fill(tmp); 
            else if(tmp[7]==5)     NtupleMC->Fill(tmp);
            
            // Check reflection
            if((passCutsValue==1 || passCutsValue==3) && fFillOnlyD0D0bar<2){
                tmp[8] = invmassD0;
                if(pdgcode!= 421){
                   tmp[7] = 6; 
                   NtupleMC->Fill(tmp);} // Reflection
            }
            if (passCutsValue>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)){
                tmp[8] = invmassD0bar;
                if(pdgcode!= -421){ 
                   tmp[7] = 6; 
                   NtupleMC->Fill(tmp);} //Reflection
            }
        }
    }
    else  // data
    {
      if(!fReadMC && fApplyBDT){
        if(fSampleSideband)// Sideband sampling
           if(gRandom->Rndm()>fBDTSidebandSamplingFraction) return;
        if(!passCutsValue) return;
        
        std::vector<Double_t> BDTClsVar;// BDT cls input
        BDTClsVar.resize(10);
        // Data fill this
        Int_t thisptbin = CutsAnalysis->PtBin(tmp[0]);
        if(thisptbin<0) return;
        Float_t *ptbin = CutsAnalysis->GetPtBinLimits();
        TString ptstring = Form("_%.0f_%.0f",ptbin[thisptbin],ptbin[thisptbin+1]);
        
        if((passCutsValue==1 || passCutsValue==3) && fFillOnlyD0D0bar<2){
            tmp[7] = 1;
            tmp[8] = invmassD0;
            tmp[17] = cosThetaStarD0;
            if(tmp[8]>2.12||tmp[8]<1.65) return;
           
            // Link variables to be used as classifier
            BDTClsVar[0] = tmp[1];      BDTClsVar[1] = tmp[2];  BDTClsVar[2] = tmp[4];  BDTClsVar[3] = tmp[9]; BDTClsVar[4] = tmp[10];
            BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];
            if(fSampleSideband&&(tmp[8]<fBDTSidebandCut[0]||tmp[8]>=fBDTSidebandCut[1])){ // Sideband sampling
               TNtuple *NtupleSB = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleSB%s",ptstring.Data()));
               AliRDHFBDT *thisbdt1 = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
               TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));
               Float_t bdt1resp = thisbdt1->GetResponse(BDTClsVar);
               if(bdt1resp>1.||bdt1resp<-1.) {std::cout<<"ERROR: Response out of range, not make sense..."<<endl;return;}
               if(bdt1resp>fBDTRespCut){ // BDT response cut
                  BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                  NtupleSB->Fill(tmp);
                  BDTRespTree->Fill();
                  BDTRespTree->ResetBranchAddresses();
               }                  
            }
            else if(!fSampleSideband){ // Data application 
                AliRDHFBDT *thisbdt1   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
                AliRDHFBDT *thisbdt2ll = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_0",ptstring.Data()));
                AliRDHFBDT *thisbdt2l  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_1",ptstring.Data()));
                AliRDHFBDT *thisbdt2ml = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_2",ptstring.Data()));
                AliRDHFBDT *thisbdt2mh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_3",ptstring.Data()));
                AliRDHFBDT *thisbdt2h  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_4",ptstring.Data()));
                AliRDHFBDT *thisbdt2hh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_5",ptstring.Data()));
                Float_t bdt1resp; Float_t bdt2resp[6];
                bdt1resp = thisbdt1->GetResponse(BDTClsVar);
                bdt2resp[0] = thisbdt2ll->GetResponse(BDTClsVar); bdt2resp[1] = thisbdt2l->GetResponse(BDTClsVar); bdt2resp[2] = thisbdt2ml->GetResponse(BDTClsVar);
                bdt2resp[3] = thisbdt2mh->GetResponse(BDTClsVar); bdt2resp[4] = thisbdt2h->GetResponse(BDTClsVar); bdt2resp[5] = thisbdt2hh->GetResponse(BDTClsVar);
                // BDT Responses ready
                if(fGetRespTree){
                    TNtuple *NtupleData = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleData%s",ptstring.Data()));
                    TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));
                    if(bdt1resp>fBDTRespCut&&bdt2resp[0]>fBDTRespCut&&bdt2resp[1]>fBDTRespCut&&bdt2resp[2]>fBDTRespCut&&bdt2resp[3]>fBDTRespCut&&bdt2resp[4]>fBDTRespCut&&bdt2resp[5]>fBDTRespCut){ // BDT response cut
                        BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                        BDTRespTree->SetBranchAddress("BDT2",&bdt2resp);
                        NtupleData->Fill(tmp[8]);
                        BDTRespTree->Fill();
                        BDTRespTree->ResetBranchAddresses();
                    }
                }
                else{
                    TH3F *h3_0 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_1",ptstring.Data()));
                    TH3F *h3_1 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_2",ptstring.Data()));
                    TH3F *h3_2 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_3",ptstring.Data()));
                    TH3F *h3_3 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_4",ptstring.Data()));
                    TH3F *h3_4 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_5",ptstring.Data()));
                    TH3F *h3_5 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_6",ptstring.Data()));
                    if(bdt1resp>fBDTRespCut){
                        if(bdt2resp[0]>fBDTRespCut) h3_0->Fill(tmp[8],bdt1resp,bdt2resp[0]);
                        if(bdt2resp[1]>fBDTRespCut) h3_1->Fill(tmp[8],bdt1resp,bdt2resp[1]);
                        if(bdt2resp[2]>fBDTRespCut) h3_2->Fill(tmp[8],bdt1resp,bdt2resp[2]);
                        if(bdt2resp[3]>fBDTRespCut) h3_3->Fill(tmp[8],bdt1resp,bdt2resp[3]);
                        if(bdt2resp[4]>fBDTRespCut) h3_4->Fill(tmp[8],bdt1resp,bdt2resp[4]);
                        if(bdt2resp[5]>fBDTRespCut) h3_5->Fill(tmp[8],bdt1resp,bdt2resp[5]);
                    }
                }
            }
        }

        if (passCutsValue>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)){
            tmp[7] = 2;
            tmp[8] = invmassD0bar;
            tmp[17] = cosThetaStarD0bar;
             
            if(tmp[8]>2.12||tmp[8]<1.65) return;
            // Link variables to be used as classifier
            BDTClsVar[0] = tmp[1];      BDTClsVar[1] = tmp[2];  BDTClsVar[2] = tmp[4];  BDTClsVar[3] = tmp[9];  BDTClsVar[4] = tmp[10];
            BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];
            if(fSampleSideband&&(tmp[8]<fBDTSidebandCut[0]||tmp[8]>=fBDTSidebandCut[1])){ // Sideband sampling
                TNtuple *NtupleSB = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleSB%s",ptstring.Data()));
                AliRDHFBDT *thisbdt1 = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
                TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));

                Float_t bdt1resp = thisbdt1->GetResponse(BDTClsVar);
                if(bdt1resp>1.||bdt1resp<-1.) {std::cout<<"ERROR: Response out of range, not make sense..."<<endl;return;}
                if(bdt1resp>fBDTRespCut){ // BDT response cut  
                    BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                    NtupleSB->Fill(tmp);
                    BDTRespTree->Fill();
                    BDTRespTree->ResetBranchAddresses();
                }
            }
            else if(!fSampleSideband){ // Data application  
               AliRDHFBDT *thisbdt1 = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
               AliRDHFBDT *thisbdt2ll = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_0",ptstring.Data()));
               AliRDHFBDT *thisbdt2l  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_1",ptstring.Data()));
               AliRDHFBDT *thisbdt2ml = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_2",ptstring.Data()));
               AliRDHFBDT *thisbdt2mh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_3",ptstring.Data()));
               AliRDHFBDT *thisbdt2h  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_4",ptstring.Data()));
               AliRDHFBDT *thisbdt2hh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_5",ptstring.Data()));
               Float_t bdt1resp; Float_t bdt2resp[6];
               bdt1resp = thisbdt1->GetResponse(BDTClsVar);
               bdt2resp[0] = thisbdt2ll->GetResponse(BDTClsVar); bdt2resp[1] = thisbdt2l->GetResponse(BDTClsVar); bdt2resp[2] = thisbdt2ml->GetResponse(BDTClsVar);
               bdt2resp[3] = thisbdt2mh->GetResponse(BDTClsVar); bdt2resp[4] = thisbdt2h->GetResponse(BDTClsVar); bdt2resp[5] = thisbdt2hh->GetResponse(BDTClsVar);
               if(fGetRespTree){
                   TNtuple *NtupleData = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleData%s",ptstring.Data()));
                   TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));
                   if(bdt1resp>fBDTRespCut&&bdt2resp[0]>fBDTRespCut&&bdt2resp[1]>fBDTRespCut&&bdt2resp[2]>fBDTRespCut&&bdt2resp[3]>fBDTRespCut&&bdt2resp[4]>fBDTRespCut&&bdt2resp[5]>fBDTRespCut){ // BDT response cut
                      BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                      BDTRespTree->SetBranchAddress("BDT2",&bdt2resp);
                      NtupleData->Fill(tmp[8]);
                      BDTRespTree->Fill();
                      BDTRespTree->ResetBranchAddresses();
                   }
               }
               else{
                  TH3F *h3_0 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_1",ptstring.Data()));
                  TH3F *h3_1 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_2",ptstring.Data()));
                  TH3F *h3_2 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_3",ptstring.Data()));
                  TH3F *h3_3 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_4",ptstring.Data()));
                  TH3F *h3_4 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_5",ptstring.Data()));
                  TH3F *h3_5 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_6",ptstring.Data()));
                  if(bdt1resp>fBDTRespCut){
                     if(bdt2resp[0]>fBDTRespCut) h3_0->Fill(tmp[8],bdt1resp,bdt2resp[0]);
                     if(bdt2resp[1]>fBDTRespCut) h3_1->Fill(tmp[8],bdt1resp,bdt2resp[1]);
                     if(bdt2resp[2]>fBDTRespCut) h3_2->Fill(tmp[8],bdt1resp,bdt2resp[2]);
                     if(bdt2resp[3]>fBDTRespCut) h3_3->Fill(tmp[8],bdt1resp,bdt2resp[3]);
                     if(bdt2resp[4]>fBDTRespCut) h3_4->Fill(tmp[8],bdt1resp,bdt2resp[4]);
                     if(bdt2resp[5]>fBDTRespCut) h3_5->Fill(tmp[8],bdt1resp,bdt2resp[5]);
                  }      
              }

            }
        }
      }
    }
}
//_____________________________________________________________________________________________________________            
void AliAnalysisTaskSEDvsMultiplicity_BDT::ProcessBDTDs(AliAODRecoDecayHF *part, AliRDHFCuts *CutsAnalysis, Double_t bfield, TClonesArray *arrayMC, Int_t passCutsValue, Int_t coutMulti)
{
   AliAODPidHF *Pid_HF = CutsAnalysis->GetPidHF();

   Int_t isKKpi = passCutsValue & 1;
   Int_t ispiKK = passCutsValue & 2;
   Int_t isPhiKKpi = passCutsValue & 4;
   Int_t isPhipiKK = passCutsValue & 8;
   Int_t isK0starKKpi = passCutsValue & 16;
   Int_t isK0starpiKK = passCutsValue & 32;

  Int_t labDs = -1;
  Int_t pdgCode0 = -999;
  Int_t pdgDstoKKpi[3]={321,321,211};
  const Int_t nProng = 3;
  AliAODMCParticle *partDs = nullptr;
  Int_t orig = 0.;
  Double_t invMass_KKpi;
  Double_t invMass_piKK;
  Double_t DeltaMassPhi = 0;
  Double_t Cos_PiKPhi =-2;
  Double_t Cos_PiDs;
  Float_t tmp[34];

  if(isPhiKKpi || isPhipiKK){
     tmp[3] = part->DecayLength();
     tmp[1] = part->Pt();
     tmp[2] = part->ImpParXY()*10000.;
     tmp[4] = coutMulti;
     tmp[5] = part->DecayLengthXY();
     tmp[6] = part->NormalizedDecayLengthXY();
     tmp[7] = part->CosPointingAngle();
     tmp[8] = part->CosPointingAngleXY();
     tmp[9] = part->GetDCA();
     tmp[10] = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(part, bfield); 
  }
  
  AliAODTrack* prongtracks[3];
  Double_t PtProng[3];
  Double_t ImParProng[3];
  for(unsigned int iProng = 0; iProng < 3; iProng++) {
     prongtracks[iProng] = (AliAODTrack*)part->GetDaughter(iProng);
     PtProng[iProng] = prongtracks[iProng]->Pt();
     ImParProng[iProng] = part->Getd0Prong(iProng);
     tmp[15+iProng]=ImParProng[iProng];
     tmp[18+iProng]=PtProng[iProng];
  }


  AliAODRecoDecayHF3Prong* cand3p = (AliAODRecoDecayHF3Prong*)part;
  tmp[11] = cand3p->GetSigmaVert();

  Double_t massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  if(isPhiKKpi){
  invMass_KKpi = cand3p->InvMassDsKKpi();
  tmp[0] = invMass_KKpi;
  tmp[12] = TMath::Abs(cand3p->InvMass2Prongs(0, 1, 321, 321) - massPhi);
  Cos_PiDs = cand3p->CosPiDsLabFrameKKpi();
  tmp[13] = Cos_PiDs;
  Cos_PiKPhi = cand3p->CosPiKPhiRFrameKKpi();
  }
  if(isPhipiKK){
  invMass_piKK = cand3p->InvMassDspiKK();
  tmp[0] = invMass_piKK;
  tmp[12] = TMath::Abs(cand3p->InvMass2Prongs(1, 2, 321, 321) - massPhi);
  Cos_PiDs = cand3p->CosPiDsLabFramepiKK();
  tmp[13] = Cos_PiDs;
  Cos_PiKPhi = cand3p->CosPiKPhiRFramepiKK();
  }
  tmp[14] = Cos_PiKPhi * Cos_PiKPhi * Cos_PiKPhi;
  
  //PID variables
  static const unsigned int knMaxProngs   = 3;
  static const unsigned int knMaxDet4Pid  = 2;  //0=>TPC; 1=>TOF
  static const unsigned int knMaxHypo4Pid = 2;
  Double_t sig[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
  Double_t PIDNsigmaVector[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
  AliPID::EParticleType parthypo[knMaxHypo4Pid] = {AliPID::kPion, AliPID::kKaon};

  for(unsigned int iProng = 0; iProng < 3; iProng++) {
     for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
          Double_t nSigmaTPC = -999.;
          Pid_HF->GetnSigmaTPC(prongtracks[iProng], parthypo[iPartHypo], nSigmaTPC);
          sig[iProng][0][iPartHypo] = nSigmaTPC;

          double nSigmaTOF = -999.;
          Pid_HF->GetnSigmaTOF(prongtracks[iProng], parthypo[iPartHypo], nSigmaTOF);
          sig[iProng][1][iPartHypo] = nSigmaTOF;
      }
  }
  Int_t index = 0;
  for(unsigned int iProng = 0; iProng < 3; iProng++) {
    for(int iDet = 0; iDet <= 1; iDet++) {
       for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
            PIDNsigmaVector[iProng][iDet][iPartHypo] = sig[iProng][iDet][iPartHypo];
            tmp[21+index]=PIDNsigmaVector[iProng][iDet][iPartHypo];
            ++index;
       }
    }
  }

  if (fReadMC){
     labDs = part->MatchToMC(431,arrayMC, nProng, pdgDstoKKpi);
     if(labDs>=0){
         partDs = (AliAODMCParticle*)arrayMC->At(labDs);
         Int_t labDau0 = ((AliAODTrack*) part->GetDaughter(0))->GetLabel();
         AliAODMCParticle *p = (AliAODMCParticle *)arrayMC->UncheckedAt(TMath::Abs(labDau0));
         pdgCode0 = TMath::Abs(p->GetPdgCode());
     }
     if(partDs)
         orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partDs,kTRUE);
     if(isPhiKKpi){
        if(labDs >=0){
           if(orig == 4) tmp[33]=4; //prompt
           else if (orig == 5)  tmp[33]=5; //non-prompt
           else tmp[33] = 0;
        }
        TNtuple *NtupleDsMC = (TNtuple*)fListBDTNtuple->FindObject("NtupleMC");
        if(tmp[33]==4)          NtupleDsMC->Fill(tmp);
        else if(tmp[33]==5)     NtupleDsMC->Fill(tmp);
     }

     if(isPhipiKK){
        if(labDs >=0){
           if(orig == 4) tmp[33]=4; //prompt
           else if (orig == 5)  tmp[33]=5; //non-prompt
           else tmp[33] = 0;
        }
        TNtuple *NtupleDsMC = (TNtuple*)fListBDTNtuple->FindObject("NtupleMC");
        if(tmp[33]==4)          NtupleDsMC->Fill(tmp);
        else if(tmp[33]==5)     NtupleDsMC->Fill(tmp);
     }
  }
   else  // data
   {
     if(!fReadMC && fApplyBDT){
       if(fSampleSideband)// Sideband sampling
       if(gRandom->Rndm()>fBDTSidebandSamplingFraction) return;
       if(!passCutsValue) return;

       std::vector<Double_t> BDTClsVar;// BDT cls input
       BDTClsVar.resize(10);
       // Data fill this
       Int_t thisptbin = CutsAnalysis->PtBin(tmp[0]);
       if(thisptbin<0) return;
       Float_t *ptbin = CutsAnalysis->GetPtBinLimits();
       TString ptstring = Form("_%.0f_%.0f",ptbin[thisptbin],ptbin[thisptbin+1]);

       tmp[33] = 0;
       if(tmp[0]>2.12||tmp[0]<1.65) return;

           // Link variables to be used as classifier
           BDTClsVar[0] = tmp[3];      BDTClsVar[1] = tmp[5];  BDTClsVar[2] = tmp[6];  BDTClsVar[3] = tmp[2]; BDTClsVar[4] = tmp[9];
           BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[12]; BDTClsVar[7] = tmp[14]; BDTClsVar[8] = tmp[15]; BDTClsVar[9] = tmp[16];
           BDTClsVar[10] = tmp[17]; BDTClsVar[11] = tmp[10]; BDTClsVar[12] = tmp[21]; BDTClsVar[13] = tmp[22]; BDTClsVar[14] = tmp[23];
           BDTClsVar[15] = tmp[24]; BDTClsVar[16] = tmp[25]; BDTClsVar[17] = tmp[26]; BDTClsVar[18] = tmp[27]; BDTClsVar[19] = tmp[28];
           BDTClsVar[20] = tmp[29]; BDTClsVar[21] = tmp[30]; BDTClsVar[22] = tmp[31]; BDTClsVar[23] = tmp[32];
           if(fSampleSideband&&(tmp[0]<fBDTSidebandCut[0]||tmp[0]>=fBDTSidebandCut[1])){ // Sideband sampling
              TNtuple *NtupleSB = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleSB%s",ptstring.Data()));
              AliRDHFBDT *thisbdt1 = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
              TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));
              Float_t bdt1resp = thisbdt1->GetResponse(BDTClsVar);
              if(bdt1resp>1.||bdt1resp<-1.) {std::cout<<"ERROR: Response out of range, not make sense..."<<endl;return;}
              if(bdt1resp>fBDTRespCut){ // BDT response cut
                 BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                 NtupleSB->Fill(tmp);
                 BDTRespTree->Fill();
                 BDTRespTree->ResetBranchAddresses();
              }
           }
           else if(!fSampleSideband){ // Data application 
                AliRDHFBDT *thisbdt1   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
               AliRDHFBDT *thisbdt2ll = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_0",ptstring.Data()));
               AliRDHFBDT *thisbdt2l  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_1",ptstring.Data()));
               AliRDHFBDT *thisbdt2ml = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_2",ptstring.Data()));
               AliRDHFBDT *thisbdt2mh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_3",ptstring.Data()));
               AliRDHFBDT *thisbdt2h  = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_4",ptstring.Data()));
               AliRDHFBDT *thisbdt2hh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_5",ptstring.Data()));
               Float_t bdt1resp; Float_t bdt2resp[6];
               bdt1resp = thisbdt1->GetResponse(BDTClsVar);
               bdt2resp[0] = thisbdt2ll->GetResponse(BDTClsVar); bdt2resp[1] = thisbdt2l->GetResponse(BDTClsVar); bdt2resp[2] = thisbdt2ml->GetResponse(BDTClsVar);
               bdt2resp[3] = thisbdt2mh->GetResponse(BDTClsVar); bdt2resp[4] = thisbdt2h->GetResponse(BDTClsVar); bdt2resp[5] = thisbdt2hh->GetResponse(BDTClsVar);
               // BDT Responses ready
               if(fGetRespTree){
                   TNtuple *NtupleDsData = (TNtuple*)fListBDTNtuple->FindObject(Form("NtupleData%s",ptstring.Data()));
                   TTree *BDTRespTree = (TTree*)fListBDTResp->FindObject(Form("BDTRespTree%s",ptstring.Data()));
                   if(bdt1resp>fBDTRespCut&&bdt2resp[0]>fBDTRespCut&&bdt2resp[1]>fBDTRespCut&&bdt2resp[2]>fBDTRespCut&&bdt2resp[3]>fBDTRespCut&&bdt2resp[4]>fBDTRespCut&&bdt2resp[5]>fBDTRespCut){ // BDT response cut
                       BDTRespTree->SetBranchAddress("BDT1",&bdt1resp);
                       BDTRespTree->SetBranchAddress("BDT2",&bdt2resp);
                       NtupleDsData->Fill(tmp[8]);
                       BDTRespTree->Fill();
                       BDTRespTree->ResetBranchAddresses();
                   }
               }
               else{
                   TH3F *h3_0 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_1",ptstring.Data()));
                   TH3F *h3_1 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_2",ptstring.Data()));
                   TH3F *h3_2 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_3",ptstring.Data()));
                   TH3F *h3_3 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_4",ptstring.Data()));
                   TH3F *h3_4 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_5",ptstring.Data()));
                   TH3F *h3_5 = (TH3F*)fListBDTResp->FindObject(Form("h3MassBDTResp%s_6",ptstring.Data()));
                   if(bdt1resp>fBDTRespCut){
                       if(bdt2resp[0]>fBDTRespCut) h3_0->Fill(tmp[8],bdt1resp,bdt2resp[0]);
                       if(bdt2resp[1]>fBDTRespCut) h3_1->Fill(tmp[8],bdt1resp,bdt2resp[1]);
                       if(bdt2resp[2]>fBDTRespCut) h3_2->Fill(tmp[8],bdt1resp,bdt2resp[2]);
                       if(bdt2resp[3]>fBDTRespCut) h3_3->Fill(tmp[8],bdt1resp,bdt2resp[3]);
                       if(bdt2resp[4]>fBDTRespCut) h3_4->Fill(tmp[8],bdt1resp,bdt2resp[4]);
                       if(bdt2resp[5]>fBDTRespCut) h3_5->Fill(tmp[8],bdt1resp,bdt2resp[5]);
                   }
               }
           }
       }  
     }

}
//_________________________________________________________________________________________________
