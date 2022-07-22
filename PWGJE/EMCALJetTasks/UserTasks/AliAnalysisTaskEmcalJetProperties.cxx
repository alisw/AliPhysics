/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
// This Class is the modified version of AliAnalysisTaskJetProperties.cxx

#include <TClonesArray.h>
#include <TF1.h>  
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TList.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"

//__________________________________________________HM____________________________________//
#include "AliAnalysisTaskEmcalJetProperties.h"
#include "AliRunLoader.h"
#include "AliVVZERO.h"
#include "AliAODZDC.h"
#include "AliVZDC.h"
#include "AliMultSelection.h"

#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliParticleContainer.h"
#include "AliInputEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"

#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliRhoParameter.h"
#include "TVector3.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"
#include "AliAnalysisTaskEA.h"
//______________________________________________________________________________________//
//#include "AliUA1JetHeaderV1.h"
#include "AliAnalysisTaskEmcalJetProperties.h"

//#include "AliAnalysisTaskRhoBase.h"
//#include "AliAnalysisTaskRho.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetProperties);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetProperties::AliAnalysisTaskEmcalJetProperties() 
  : AliAnalysisTaskEmcalJet()
  ,fTrackListUE(0)
  ,matchedGenJetList(0)
  ,matchedRecJetList(0)
  ,fJets(0)
  ,fh1PtJet1D(0)

  ,fhistMultV0A(0)
  ,fhistMultV0C(0)
  ,fhistMultV0M(0)
    
  ,fh1PtLeadingJet(0)
    
  ,fh1PtLeadingJet_PtHardBinNch(0)
  ,fh1PtLeadingJet_PtHardBinFF(0)

  ,fh1PtSumInJetConeUE2(0)
  ,fh1PtSumInJetConeUE3(0)
  ,fh2NtracksLeadingJet(0)
  ,fh2NtracksLeadingJetUE2(0)
  ,fh2NtracksLeadingJetUE3(0)
  ,fProNtracksLeadingJet(0)
  ,fProNtracksLeadingJetUE2(0)
  ,fProNtracksLeadingJetUE3(0)
  ,fh2PtTrack(0)
  ,fh2PtTrackUE2(0)
  ,fh2PtTrackUE3(0)
  ,fh2FF(0)
  ,fh2FFUE2(0)
  ,fh2FFUE3(0)
  ,fh2Ksi(0)
  ,fh2KsiUE2(0)
  ,fh2KsiUE3(0)
  ,fh2DelR80pcPt(0)
  ,fh2DelR80pcPtUE2(0)
  ,fh2DelR80pcPtUE3(0)
  ,fProDelR80pcPt(0)
  ,fProDelR80pcPtUE2(0)
  ,fProDelR80pcPtUE3(0)
  ,fh3PtDelRPtSum(0)
  ,fh3PtDelRPtSumUE2(0)
  ,fh3PtDelRPtSumUE3(0)
  ,fProDelRPtDenUEN2(0)
  ,fProDelRPtDenUEN3(0)

  ,fh1DelRRCJet1(0)
  ,fh1DelRRCJet2(0)
  ,fh1nRC(0)
    

  ,hProDelRPtDenUEN(0)
  ,fh1rho(0)
  ,fh1rhoN(0)
  ,h1LeadingJetPtWUE(0)
  ,h1LeadingJetPtWoUEWoCAll(0)
  ,fProNtracksLeadingJetWUE(0)
  ,fProNtracksLeadingJetUEWoCAll(0)
  ,fProNtracksLeadingJetWoUEWoCAll(0)
  ,fProNtracksLeadingJetUEwrtSubtPtWoC(0)
  ,fProNtracksLeadingJetWoUEwrtSubtPtWoC(0)
  ,h1Cfactor(0)
  ,h1LeadingJetPtWoUEWCAll(0)
  ,fProNtracksLeadingJetUEWCAll(0)
  ,fProNtracksLeadingJetWoUEWCAll(0)
  ,fProNtracksLeadingJetUEwrtSubtPtWC(0)
  ,fProNtracksLeadingJetWoUEwrtSubtPtWC(0)

  ,h4ResNch(0)
  ,h2NchGen(0)  
  ,h2NchRec(0)

  ,h4ResNchUE(0)
  ,h2NchGenUE(0)  
  ,h2NchRecUE(0)  

  ,h4ResFF(0)
  ,h2FFGen(0)  
  ,h2FFRec(0)
  ,h4ResFFM(0)
  ,h2FFMFake(0)  
  ,h2FFMMiss(0)  
  ,h2FFMGen(0)  
  ,h2FFMRec(0)

  ,h4ResFFMUE(0)
  ,h2FFMUEFake(0)  
  ,h2FFMUEMiss(0)  
  ,h2FFMUEGen(0)  
  ,h2FFMUERec(0)

  ,h5ResPtSum(0)
  ,h5ResPtSumUE(0)

  ,h2ResJetPt(0)
  ,h2ResJetPtNch(0)


  ,fMultV0A(0.)
  ,fMultV0C(0.)
  ,fMultV0M(0.)
    

  ,fHistManager()
{
  for(Int_t ii=0; ii<21; ii++){
    fProDelRPtSum[ii]     = NULL;
    fProDelRPtSumUE2[ii]     = NULL;
    fProDelRPtSumUE3[ii]     = NULL;
    fProDelRPtDenUE2[ii]     = NULL;
    fProDelRPtDenUE3[ii]     = NULL;

    hhProDelRPtSumUE[ii]     = NULL;
  }
  // default constructor
}  

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetProperties::AliAnalysisTaskEmcalJetProperties(const char *name)
  : AliAnalysisTaskEmcalJet(name, kTRUE)

  ,fhistMultV0A(0)
  ,fhistMultV0C(0)
  ,fhistMultV0M(0)

  ,fTrackListUE(0)
  ,matchedGenJetList(0)
  ,matchedRecJetList(0)
  ,fJets(0) 
  ,fh1PtJet1D(0)
    
  ,fh1PtLeadingJet(0)

  ,fh1PtLeadingJet_PtHardBinNch(0)
  ,fh1PtLeadingJet_PtHardBinFF(0)

  ,fh1PtSumInJetConeUE2(0)
  ,fh1PtSumInJetConeUE3(0)
  ,fh2NtracksLeadingJet(0)
  ,fh2NtracksLeadingJetUE2(0)
  ,fh2NtracksLeadingJetUE3(0)
  ,fProNtracksLeadingJet(0)
  ,fProNtracksLeadingJetUE2(0)
  ,fProNtracksLeadingJetUE3(0)
  ,fh2PtTrack(0)
  ,fh2PtTrackUE2(0)
  ,fh2PtTrackUE3(0)
  ,fh2FF(0)
  ,fh2FFUE2(0)
  ,fh2FFUE3(0)
  ,fh2Ksi(0)
  ,fh2KsiUE2(0)
  ,fh2KsiUE3(0)
  ,fh2DelR80pcPt(0)
  ,fh2DelR80pcPtUE2(0)
  ,fh2DelR80pcPtUE3(0)
  ,fProDelR80pcPt(0)
  ,fProDelR80pcPtUE2(0)
  ,fProDelR80pcPtUE3(0)
  ,fh3PtDelRPtSum(0)
  ,fh3PtDelRPtSumUE2(0)
  ,fh3PtDelRPtSumUE3(0)
  ,fProDelRPtDenUEN2(0)
  ,fProDelRPtDenUEN3(0)

  ,fh1DelRRCJet1(0)
  ,fh1DelRRCJet2(0)
  ,fh1nRC(0)

  ,hProDelRPtDenUEN(0)    
  ,fh1rho(0)
  ,fh1rhoN(0)
  ,h1LeadingJetPtWUE(0)
  ,h1LeadingJetPtWoUEWoCAll(0)
  ,fProNtracksLeadingJetWUE(0)
  ,fProNtracksLeadingJetUEWoCAll(0)
  ,fProNtracksLeadingJetWoUEWoCAll(0)
  ,fProNtracksLeadingJetUEwrtSubtPtWoC(0)
  ,fProNtracksLeadingJetWoUEwrtSubtPtWoC(0)
  ,h1Cfactor(0)
  ,h1LeadingJetPtWoUEWCAll(0)
  ,fProNtracksLeadingJetUEWCAll(0)
  ,fProNtracksLeadingJetWoUEWCAll(0)
  ,fProNtracksLeadingJetUEwrtSubtPtWC(0)
  ,fProNtracksLeadingJetWoUEwrtSubtPtWC(0)

  ,h4ResNch(0)
  ,h2NchGen(0)  
  ,h2NchRec(0)  

  ,h4ResNchUE(0)
  ,h2NchGenUE(0)  
  ,h2NchRecUE(0)

  ,h4ResFF(0)
  ,h2FFGen(0)  
  ,h2FFRec(0)
  ,h4ResFFM(0)
  ,h2FFMFake(0)  
  ,h2FFMMiss(0)  
  ,h2FFMGen(0)  
  ,h2FFMRec(0) 

  ,h4ResFFMUE(0)
  ,h2FFMUEFake(0)  
  ,h2FFMUEMiss(0)  
  ,h2FFMUEGen(0)  
  ,h2FFMUERec(0)

  ,h5ResPtSum(0)
  ,h5ResPtSumUE(0)

  ,h2ResJetPt(0)
  ,h2ResJetPtNch(0)

  ,fMultV0A(0.)
  ,fMultV0C(0.)
  ,fMultV0M(0.)
    
  ,fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);
  for(Int_t ii=0; ii<21; ii++){
    fProDelRPtSum[ii]     = NULL;
    fProDelRPtSumUE2[ii]     = NULL;
    fProDelRPtSumUE3[ii]     = NULL;
    fProDelRPtDenUE2[ii]     = NULL;
    fProDelRPtDenUE3[ii]     = NULL;

    hhProDelRPtSumUE[ii]     = NULL;
  }
}//constructor

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetProperties::~AliAnalysisTaskEmcalJetProperties()
{
  if(fTrackListUE) delete fTrackListUE;
  if(matchedGenJetList) delete matchedGenJetList;
  if(matchedRecJetList) delete matchedRecJetList;
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetProperties::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateJetHistograms();

  fTrackListUE = new TList();
  fTrackListUE->SetOwner(kFALSE);
  
  matchedGenJetList = new TList();
  matchedGenJetList->SetOwner(kFALSE);
  matchedRecJetList = new TList();
  matchedRecJetList->SetOwner(kFALSE);


  h2ResJetPt = new TH2D("h2ResJetPt","h2ResJetPt", 24, 0., 120., 24, 0., 120.);
  h2ResJetPtNch = new TH2D("h2ResJetPtNch","h2ResJetPtNch", 48, 0., 120., 48, 0., 120.);

  const Int_t kNbinsNch = 30;
  const Int_t kNbinsJetPt = 48;

  //Response Object
  const Int_t fNumberOfDimensions = 4;
  const Int_t fBinNumbers[fNumberOfDimensions] = {kNbinsJetPt, kNbinsNch,kNbinsJetPt, kNbinsNch};
  const Double_t x_min_value[fNumberOfDimensions] = {0., 0., 0., 0.};
  const Double_t x_max_value[fNumberOfDimensions] = {120., 30., 120., 30.};

  h4ResNch = new THnSparseD ("h4ResNch", "Nch vs Jet pT for filling response matrix", fNumberOfDimensions, fBinNumbers, x_min_value, x_max_value);
  h2NchGen = new TH2D("h2NchGen","h2NchGen", kNbinsJetPt, 0., 120., kNbinsNch, 0., 30.);
  h2NchRec = new TH2D("h2NchRec","h2NchRec", kNbinsJetPt, 0., 120., kNbinsNch, 0., 30.);

  h4ResNchUE = new THnSparseD ("h4ResNchUE", "Nch UE vs Jet pT for filling response matrix", fNumberOfDimensions, fBinNumbers, x_min_value, x_max_value);
  h2NchGenUE = new TH2D("h2NchGenUE","h2NchGenUE", kNbinsJetPt, 0., 120., kNbinsNch, 0., 30.);
  h2NchRecUE = new TH2D("h2NchRecUE","h2NchRecUE", kNbinsJetPt, 0., 120., kNbinsNch, 0., 30.);

  const Int_t kNbinsJetPtFF = 24;
  const Int_t kNbinsFFR = 55;     //Float_t xMinFF=-0.01;    Float_t xMaxFF=1.09;
  const Int_t fBinNumbers_FF[fNumberOfDimensions] = {kNbinsJetPtFF, kNbinsFFR,kNbinsJetPtFF, kNbinsFFR};
  const Double_t x_min_value_FF[fNumberOfDimensions] = {0., 0., 0., 0.};
  const Double_t x_max_value_FF[fNumberOfDimensions] = {120., 1.1, 120., 1.1};

  h4ResFF = new THnSparseD ("h4ResFF", "FF vs Jet pT for filling response matrix", fNumberOfDimensions, fBinNumbers_FF, x_min_value_FF, x_max_value_FF);
  h4ResFFM = new THnSparseD ("h4ResFFM", "Matched FF vs Jet pT for filling response matrix", fNumberOfDimensions, fBinNumbers_FF, x_min_value_FF, x_max_value_FF);
  h4ResFFMUE = new THnSparseD ("h4ResFFMUE", "Matched FF UE vs Jet pT for filling response matrix", fNumberOfDimensions, fBinNumbers_FF, x_min_value_FF, x_max_value_FF);
  
  h2FFGen = new TH2D("h2FFGen","h2FFGen", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFRec = new TH2D("h2FFRec","h2FFRec", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMGen = new TH2D("h2FFMGen","h2FFMGen", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMRec = new TH2D("h2FFMRec","h2FFMRec", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMFake = new TH2D("h2FFMFake","h2FFMFake", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMMiss = new TH2D("h2FFMMiss","h2FFMMiss", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);

  h2FFMUEGen = new TH2D("h2FFMUEGen","h2FFMUEGen", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMUERec = new TH2D("h2FFMUERec","h2FFMUERec", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMUEFake = new TH2D("h2FFMUEFake","h2FFMUEFake", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);
  h2FFMUEMiss = new TH2D("h2FFMUEMiss","h2FFMUEMiss", kNbinsJetPtFF, 0., 120., kNbinsFFR, 0., 1.1);


  const Int_t kNbinsJetPtPtSum = 24;
  const Int_t kNbinsR = 20;
  const Int_t kNbinsPtSum = 48;

  const Int_t fBinNumbersPtSum[5] = {kNbinsJetPtPtSum, kNbinsJetPtPtSum, kNbinsR, kNbinsPtSum, kNbinsPtSum};
  const Double_t x_min_value_PtSum[5] = {0., 0., 0., 0., 0.};
  const Double_t x_max_value_PtSum[5] = {120., 120., 0.4, 120., 120.};

  h5ResPtSum = new THnSparseD ("h5ResPtSum", "PtSum vs r for filling response matrix", 5, fBinNumbersPtSum, x_min_value_PtSum, x_max_value_PtSum);
  h5ResPtSumUE = new THnSparseD ("h5ResPtSumUE", "PtSum UE vs r for filling response matrix", 5, fBinNumbersPtSum, x_min_value_PtSum, x_max_value_PtSum);

  //_____________________________________________________HM part____________________________________________________________//                               
  //name = "V0Mmult";
                                                                                                                                                             
  fhistMultV0A = new TH1D("V0Amult", "V0Amultiplicity ", 2000, 0, 2000.0);
  fhistMultV0C = new TH1D("V0Cmult", "V0Cmultiplicity ", 2000, 0, 2000.0);
  fhistMultV0M = new TH1D("V0Mmult", "V0Mmultiplicity ", 2000, 0, 2000.0);

  fOutput->Add(fhistMultV0A);
  fOutput->Add(fhistMultV0C);
  fOutput->Add(fhistMultV0M);

  //--------------------------------------------------------------------------------------------------//           
  
  fOutput->Add((TH2D*)h2ResJetPt);
  fOutput->Add((TH2D*)h2ResJetPtNch);
  fOutput->Add((THnSparse*) h4ResNch);
  fOutput->Add((TH2D*)h2NchGen);
  fOutput->Add((TH2D*)h2NchRec);

  fOutput->Add((THnSparse*) h4ResNchUE);
  fOutput->Add((TH2D*)h2NchGenUE);
  fOutput->Add((TH2D*)h2NchRecUE);

  fOutput->Add((THnSparse*) h4ResFF);
  fOutput->Add((TH2D*)h2FFGen);
  fOutput->Add((TH2D*)h2FFRec);
  fOutput->Add((THnSparse*) h4ResFFM);
  fOutput->Add((TH2D*)h2FFMGen);
  fOutput->Add((TH2D*)h2FFMRec);
  fOutput->Add((TH2D*)h2FFMFake);
  fOutput->Add((TH2D*)h2FFMMiss);
  
  fOutput->Add((THnSparse*) h4ResFFMUE);
  fOutput->Add((TH2D*)h2FFMUEGen);
  fOutput->Add((TH2D*)h2FFMUERec);
  fOutput->Add((TH2D*)h2FFMUEFake);
  fOutput->Add((TH2D*)h2FFMUEMiss);

  fOutput->Add((THnSparse*) h5ResPtSum);
  fOutput->Add((THnSparse*) h5ResPtSumUE);


  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetProperties::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.
Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);

    Int_t kNbinsPt=250;     Float_t xMinPt=0.0;      Float_t xMaxPt=250.0;
    Int_t kNbinsPtSlice=20; Float_t xMinPtSlice=0.0; Float_t xMaxPtSlice=100.0;
    Int_t kNbinsEta=40;     Float_t xMinEta=-2.0;    Float_t xMaxEta=2.0;
    Int_t kNbinsPhi=90;     Float_t xMinPhi=0.0;     Float_t xMaxPhi=TMath::TwoPi();
    Int_t kNbinsNtracks=50; Float_t xMinNtracks=0.0; Float_t xMaxNtracks=50.0;
    Int_t kNbinsFF=200;     Float_t xMinFF=-0.05;    Float_t xMaxFF=1.95;
    Int_t kNbinsKsi=80;     Float_t xMinKsi=0.;      Float_t xMaxKsi = 8.0;
    Int_t kNbinsDelR1D=50;  Float_t xMinDelR1D=0.0;  Float_t xMaxDelR1D=1.0;
    Int_t kNbinsjT=100;     Float_t xMinjT=0.0;      Float_t xMaxjT=10.0;

    Int_t kNbinsDelR=100;      Float_t xMinDelR=0.0;      Float_t xMaxDelR=1.0;
    Int_t kNbinsPtSliceJS=100; Float_t xMinPtSliceJS=0.0; Float_t xMaxPtSliceJS=250.0;

    /////// FillJetProperties 
    histname  = TString::Format("%s/fh1PtJet1D", groupname.Data());
    histtitle = TString::Format("%s;jet pt;#jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 
			   kNbinsPt,      xMinPt,      xMaxPt,"s");

    /////// FillJetShape, FillJetShapeUE2 (PC), FillJetShapeUE3 (RC)
    histname  = TString::Format("%s/fh1PtLeadingJet",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   kNbinsPt,      xMinPt,      xMaxPt,"s");

    /*****/
    histname  = TString::Format("%s/fh1PtLeadingJet_PtHardBinNch",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   48, 0.0, 120.0,"s");

    histname  = TString::Format("%s/fh1PtLeadingJet_PtHardBinFF",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   24, 0.0, 120.0,"s");

    /*****/


    histname  = TString::Format("%s/fh1PtSumInJetConeUE2", groupname.Data());
    histtitle = TString::Format("%s;p_{T}^{sum, UE}(in cone R);#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   500,      0.,      100.,"s"); // kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fh1PtSumInJetConeUE3", groupname.Data());
    histtitle = TString::Format("%s;p_{T}^{sum, UE}(in cone R);#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   500,      0.,      100.,"s"); // kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fh2NtracksLeadingJet",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   48, 0.0, 120.0,
			   30,   0.0, 30.0,"s");

    histname  = TString::Format("%s/fh2NtracksLeadingJetUE2",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});N_{ch}^{UE}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   48, 0.0, 120.0,
			   30,   0.0, 30.0,"s");

    histname  = TString::Format("%s/fh2NtracksLeadingJetUE3",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});N_{ch}^{UE}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS,
			   kNbinsNtracks,   xMinNtracks,   xMaxNtracks,"s");

    histname  = TString::Format("%s/fProNtracksLeadingJet",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetUE2",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}^{UE}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetUE3",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}^{UE}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fh2PtTrack", groupname.Data());
    histtitle = TString::Format("%s;jet pt;track pt", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   2*kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fh2PtTrackUE2", groupname.Data());
    histtitle = TString::Format("%s;jet pt;track pt", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   2*kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fh2PtTrackUE3", groupname.Data());
    histtitle = TString::Format("%s;jet pt;track pt", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   2*kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fh2FF", groupname.Data());
    histtitle = TString::Format("%s;jet pt;FF", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   24, 0., 120.,
			   55, 0., 1.1,"s");

    histname  = TString::Format("%s/fh2FFUE2", groupname.Data());
    histtitle = TString::Format("%s;jet pt;FF", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   24, 0., 120.,
			   55, 0., 1.1,"s");

    histname  = TString::Format("%s/fh2FFUE3", groupname.Data());
    histtitle = TString::Format("%s;jet pt;FF", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   kNbinsFF,      xMinFF,      xMaxFF,"s");

    histname  = TString::Format("%s/fh2Ksi", groupname.Data());
    histtitle = TString::Format("%s;jet pt;Ksi", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   kNbinsKsi,      xMinKsi,      xMaxKsi,"s");

    histname  = TString::Format("%s/fh2KsiUE2", groupname.Data());
    histtitle = TString::Format("%s;jet pt;Ksi", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   kNbinsKsi,      xMinKsi,      xMaxKsi,"s");

    histname  = TString::Format("%s/fh2KsiUE3", groupname.Data());
    histtitle = TString::Format("%s;jet pt;Ksi", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSlice, xMinPtSlice, xMaxPtSlice,
			   kNbinsKsi,      xMinKsi,      xMaxKsi,"s");

    histname  = TString::Format("%s/fh2DelR80pcPt",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R_{80}^{pT}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS,
			   kNbinsDelR,      xMinDelR,      xMaxDelR,"s");

    histname  = TString::Format("%s/fh2DelR80pcPtUE2",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R_{80}^{pT, UE}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS,
			   kNbinsDelR,      xMinDelR,      xMaxDelR,"s");

    histname  = TString::Format("%s/fh2DelR80pcPtUE3",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R_{80}^{pT, UE}", histname.Data());
    fHistManager.CreateTH2(histname, histtitle,
			   kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS,
			   kNbinsDelR,      xMinDelR,      xMaxDelR,"s");

    histname  = TString::Format("%s/fProDelR80pcPt",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<R_{80}^{pT}", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProDelR80pcPtUE2",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<R_{80}^{pT, UE}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProDelR80pcPtUE3",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<R_{80}^{pT, UE}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fh3PtDelRPtSum",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R;p_{T}^{sum}", histname.Data());
    fHistManager.CreateTH3(histname, histtitle,
			   24, 0.0, 120.0,
			   20, 0.0, 0.4,
			   48, 0.0, 120.0,"s");

    histname  = TString::Format("%s/fh3PtDelRPtSumUE2",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R;p_{T}^{sum}", histname.Data());
    fHistManager.CreateTH3(histname, histtitle,
			   24, 0.0, 120.0,
			   20, 0.0, 0.4,
			   48, 0.0, 120.0,"s");

    histname  = TString::Format("%s/fh3PtDelRPtSumUE3",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});R;p_{T}^{sum}", histname.Data());
    fHistManager.CreateTH3(histname, histtitle,
			   kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS,
			   kNbinsDelR1D,    xMinDelR1D,    xMaxDelR1D,
			   kNbinsPt,        xMinPt,        xMaxPt,"s");

    histname  = TString::Format("%s/fProDelRPtDenUEN2",groupname.Data());
    histtitle = TString::Format("%s;r;<p_{T}^{den, UE2}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

    histname  = TString::Format("%s/fProDelRPtDenUEN3",groupname.Data());
    histtitle = TString::Format("%s;r;<p_{T}^{den, UE3}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

    const int kNbinsInPt = 21;
    Float_t jetPtMinA[kNbinsInPt] = {5.,
				     10., 15., 10., 
				     20., 25., 20., 
				     30., 35., 30.,30.,
				     40., 50., 40.,
				     60., 70., 60.,50.,
				     80., 90., 80};
    Float_t jetPtMaxA[kNbinsInPt] = {10., 
				     15., 20., 20.,
				     25., 30., 30.,
				     35., 40., 40.,50., 
				     50., 60., 60.,
				     70., 80., 80.,80., 
				     90., 100., 100.};
    TString namePtBinsA[kNbinsInPt];    
    for(Int_t ii=0; ii<kNbinsInPt; ii++){
      namePtBinsA[ii] = Form("_JetPt%dto%d",(Int_t)jetPtMinA[ii],(Int_t)jetPtMaxA[ii]);
      
      histname  = TString::Format("%s/fProDelRPtSum%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{sum}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);
 
      histname  = TString::Format("%s/fProDelRPtSumUE2%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{sum}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

      histname  = TString::Format("%s/fProDelRPtSumUE3%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{sum}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

      histname  = TString::Format("%s/fProDelRPtDenUE2%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{Den,UE2}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);
      
      histname  = TString::Format("%s/fProDelRPtDenUE3%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{Den,UE3}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

      //BkgSubtracted
      histname  = TString::Format("%s/hhProDelRPtSumUE%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{sum}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

      histname  = TString::Format("%s/hhProDelRPtDenUE%s",groupname.Data(),namePtBinsA[ii].Data());
      histtitle = TString::Format("%s;R;<p_{T}^{Den}>", histname.Data());
      fHistManager.CreateTProfile(histname, histtitle,
				  kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);
      
    }      

    histname  = TString::Format("%s/fh1DelRRCJet1", groupname.Data());
    histtitle = TString::Format("%s;Distance between Jet1 and RC;", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   500,      0.,      10.,"s");

    histname  = TString::Format("%s/fh1DelRRCJet2", groupname.Data());
    histtitle = TString::Format("%s;Distance between Jet2 and RC;", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   500,      0.,      10.,"s");

    histname  = TString::Format("%s/fh1nRC", groupname.Data());
    histtitle = TString::Format("%s;N0. of random cones;", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   1,      0.,      1.,"s");

    //BkgSubtracted
    histname  = TString::Format("%s/hProDelRPtDenUEN",groupname.Data());
    histtitle = TString::Format("%s;r;<p_{T}^{den, UE}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsDelR1D,  xMinDelR1D,  xMaxDelR1D);

    histname  = TString::Format("%s/fh1rho", groupname.Data());
    histtitle = TString::Format("%s;rho;#jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 
			   400000,      0.0,      40.0,"s");

    histname  = TString::Format("%s/fh1rhoN", groupname.Data());
    histtitle = TString::Format("%s;rhoN;#jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 
			   400000,      0.0,      40.0,"s");

    histname  = TString::Format("%s/h1LeadingJetPtWUE",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/h1LeadingJetPtWoUEWoCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fProNtracksLeadingJetWUE",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetUEWoCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetWoUEWoCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetUEwrtSubtPtWoC",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetWoUEwrtSubtPtWoC",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/h1Cfactor",groupname.Data());
    histtitle = TString::Format("%s;#it{C factor};#evts", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   200,      -0.5,      1.5,"s");

    histname  = TString::Format("%s/h1LeadingJetPtWoUEWCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#leading jets", histname.Data());
    fHistManager.CreateTH1(histname, histtitle,
			   kNbinsPt,      xMinPt,      xMaxPt,"s");

    histname  = TString::Format("%s/fProNtracksLeadingJetUEWCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetWoUEWCAll",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetUEwrtSubtPtWC",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);

    histname  = TString::Format("%s/fProNtracksLeadingJetWoUEwrtSubtPtWC",groupname.Data());
    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});<N_{ch}>", histname.Data());
    fHistManager.CreateTProfile(histname, histtitle,
				kNbinsPtSliceJS, xMinPtSliceJS, xMaxPtSliceJS);
  }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetProperties::FillHistograms()
{
  //std::cout<<"----------------I am inside FillHistograms----------------------"<<std::endl;
  //std::cout<<"----------------I am inside FillHistograms----------------------"<<std::endl;
  //std::cout<<"----------------I am inside FillHistograms----------------------"<<std::endl;
  
  FillJetProperties();
  FillJetShape();
  FillJetShapeUE2(TMath::Pi()/2.);
  /////////////////////////////////////////////
  //Double_t NchSubt = BkgSubtracted(0,2);
  //cout << "NchSubt: " << NchSubt << endl;
  /////////////////////////////////////////////

  //RandomConeCalculation();
 
  //  }
  return kTRUE;
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetProperties::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetProperties::Run()
{

  AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
  fMultV0A = vzeroAOD->GetMTotV0A(); //returns total multiplicity in V0A 
  fMultV0C = vzeroAOD->GetMTotV0C(); //returns total multiplicity in V0C 
  fMultV0M  = fMultV0A + fMultV0C;   //returns total multiplicity in V0A+V0C 
                                                                                                                                                              
  fhistMultV0A->Fill(fMultV0A);
  fhistMultV0C->Fill(fMultV0C);
  fhistMultV0M->Fill(fMultV0M);

  //cout<< " V0A Multiplicity:   " << fMultV0A << "  V0C Multiplicity:   "<<fMultV0C <<  "  V0M Multiplicity:   "<<fMultV0M<< endl;
  
  if(fContainerKtRec == 0 && fContainerKtGen == 1 && fContainerAktRec == 2 && fContainerAktGen == 3)
  {
            
    //cout << "Response matrix filling" << endl;
    FillMy2DResponseMatrices();
    return kFALSE;
  }
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetProperties::Terminate(Option_t *) 
{
}

void AliAnalysisTaskEmcalJetProperties::FillJetProperties(){
  
  //cout << "I am inside FillJetProperties" << endl;
  TString groupname;
  TString histname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    //cout << "Inside FillJetProp, groupname: " << groupname << endl;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      Float_t JetPt = 0.;
      JetPt  = (Float_t)jet->Pt();

      Double_t sumPtPerp = 0;
      Int_t Ntracks_in_UE = 0;
      fTrackListUE->Clear();
      
      float fJetRadius = jetCont -> GetJetRadius();
      GetTracksTiltedwrpJetAxis(TMath::Pi()/2,fTrackListUE,jet,fJetRadius,sumPtPerp,Ntracks_in_UE);
      
      Int_t nTracksUE = fTrackListUE->GetEntries();
      //cout << "Inside FillJetProp, Jet pT: " << JetPt << " UE Ntracks: " << nTracksUE << endl;

      //cout << "FillJetProp: JetPt: " << JetPt << endl;
      histname = TString::Format("%s/fh1PtJet1D", groupname.Data());
      fHistManager.FillTH1(histname, JetPt);
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskEmcalJetProperties::FillJetShape(){
  //filling up the histograms for jet shape for leading jets only
  const int kNbinsInPt = 21;
  Float_t jetPtMinA[kNbinsInPt] = {5.,
				   10., 15., 10., 
				   20., 25., 20., 
				   30., 35., 30.,30.,
				   40., 50., 40.,
				   60., 70., 60.,50.,
				   80., 90., 80};
  Float_t jetPtMaxA[kNbinsInPt] = {10., 
				   15., 20., 20.,
				   25., 30., 30.,
				   35., 40., 40.,50., 
				   50., 60., 60.,
				   70., 80., 80.,80., 
				   90., 100., 100.};
  TString namePtBinsA[kNbinsInPt];
  for(Int_t iii=0; iii<kNbinsInPt; iii++){
    namePtBinsA[iii] = Form("_JetPt%dto%d",(Int_t)jetPtMinA[iii],(Int_t)jetPtMaxA[iii]);
  }//iii loop
  TString groupname;
  TString histname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    //cout << "Inside FillJetShape, groupname: " << groupname << endl;
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      if(count > 0)continue;//skipping all jets except leading jet
      Float_t JetEta; Float_t JetPhi; Float_t JetPt;
      JetEta = (Float_t)jet->Eta();
      JetPhi = (Float_t)jet->Phi();
      JetPt  = (Float_t)jet->Pt() ;
      
      AliEmcalJet *leadJet = jetCont->GetLeadingJet();//seleting leading jet for consistency check
      Float_t leadingJetPt  = (Float_t)leadJet->Pt();
      if(leadingJetPt - JetPt) {Printf("%s:%d Mismatch in leading jet pT: %f %f",(char*)__FILE__,__LINE__,leadingJetPt,JetPt);return;}
      
      histname  = TString::Format("%s/fh1PtLeadingJet",groupname.Data());
      fHistManager.FillTH1(histname, JetPt);

      histname  = TString::Format("%s/fh1PtLeadingJet_PtHardBinNch",groupname.Data());
      fHistManager.FillTH1(histname, JetPt);
      histname  = TString::Format("%s/fh1PtLeadingJet_PtHardBinFF",groupname.Data());
      fHistManager.FillTH1(histname, JetPt);

      Float_t PtSumA[50]          = {0.};
      Float_t delRPtSum80pc       = 0;
      //Int_t kNbinsR               = 10;
      
      Int_t nJetTracks = jet->GetNumberOfTracks();
      //cout << groupname << "**********************FillJetShape************ nJetTracks: " << nJetTracks << "   JetPt: " << JetPt << endl;     
      AliParticleContainer* tracks = jetCont->GetParticleContainer();
      
      histname = TString::Format("%s/fh2NtracksLeadingJet", groupname.Data());
      fHistManager.FillTH2(histname, JetPt,nJetTracks);
      
      histname = TString::Format("%s/fProNtracksLeadingJet", groupname.Data());
      fHistManager.FillProfile(histname, JetPt,nJetTracks);
      
      Int_t   *index     = new Int_t   [nJetTracks];//dynamic array containing index
      Float_t *delRA     = new Float_t [nJetTracks];//dynamic array containing delR
      Float_t *delEtaA   = new Float_t [nJetTracks];//dynamic array containing delEta
      Float_t *delPhiA   = new Float_t [nJetTracks];//dynamic array containing delPhi
      Float_t *trackPtA  = new Float_t [nJetTracks];//dynamic array containing pt-track
      Float_t *trackEtaA = new Float_t [nJetTracks];//dynamic array containing eta-track
      Float_t *trackPhiA = new Float_t [nJetTracks];//dynamic array containing phi-track
      for(Int_t ii=0; ii<nJetTracks; ii++){
	index[ii]     = 0;
	delRA[ii]     = 0.;
	delEtaA[ii]   = 0.;
	delPhiA[ii]   = 0.;
	trackPtA[ii]  = 0.;
	trackEtaA[ii] = 0.;
	trackPhiA[ii] = 0.;
      }//ii loop
      
      for (Int_t j =0; j< nJetTracks; j++){
	Float_t TrackEta=-99.0; Float_t TrackPt=-99.0; Float_t TrackPhi=-99.0;
	Float_t DelEta=-99.0;  Float_t DelPhi=-99.0; Float_t DelR=-99.0;
	Float_t FF=-99.0;Float_t Ksi=-99.0;

	AliVParticle* trackkine = (AliVParticle*) jet->TrackAt(j,tracks->GetArray());
	if(!trackkine)continue;
	TrackEta = trackkine->Eta();
	TrackPhi = trackkine->Phi();
	TrackPt  = trackkine->Pt();
	
	if(JetPt)FF = TrackPt/JetPt;
	if(FF)Ksi   = TMath::Log(1./FF);
	
	histname = TString::Format("%s/fh2PtTrack", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, TrackPt);
	histname = TString::Format("%s/fh2FF", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, FF);
	histname = TString::Format("%s/fh2Ksi", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, Ksi);
	
	DelEta = TMath::Abs(JetEta - TrackEta);
	DelPhi = TMath::Abs(JetPhi - TrackPhi);
	if(TMath::Abs(DelPhi)>TMath::Pi())DelPhi = TMath::Abs(DelPhi-TMath::TwoPi());
	DelR   = TMath::Sqrt(DelEta*DelEta + DelPhi*DelPhi);

	delRA[j]     = DelR;
	delEtaA[j]   = DelEta;
	delPhiA[j]   = DelPhi;
	trackPtA[j]  = TrackPt;
	trackEtaA[j] = TrackEta;
	trackPhiA[j] = TrackPhi;
	
	for(Int_t ibin=1; ibin<=50; ibin++){
	  Float_t xlow = 0.02*(ibin-1);
	  Float_t xup  = 0.02*ibin;
	  if( xlow <= DelR && DelR < xup){
	    PtSumA[ibin-1]+= TrackPt;
	  }//if loop
	}//for ibin loop
      }//track loop
      
      Float_t PtSum = 0;
      Bool_t iflagPtSum   = kFALSE;

      TMath::Sort(nJetTracks,delRA,index,0);
      for(Int_t ii=0; ii<nJetTracks; ii++){
	PtSum += trackPtA[index[ii]];

	if(!iflagPtSum){
	  if(PtSum/JetPt >= 0.8000){
	    delRPtSum80pc = delRA[index[ii]];
	    iflagPtSum = kTRUE;
	  }//if PtSum
	}//if iflag
      }//track loop 2nd
      delete [] index;
      delete [] delRA;
      delete [] delEtaA;
      delete [] delPhiA;
      delete [] trackPtA;
      delete [] trackEtaA;
      delete [] trackPhiA;

      histname  = TString::Format("%s/fh2DelR80pcPt",groupname.Data());
      fHistManager.FillTH2(histname,JetPt,delRPtSum80pc);
      histname  = TString::Format("%s/fProDelR80pcPt",groupname.Data());
      fHistManager.FillProfile(histname,JetPt,delRPtSum80pc);

      for(Int_t ibin=0; ibin<50; ibin++){
	Float_t iR = 0.02*ibin + 0.01;
	histname  = TString::Format("%s/fh3PtDelRPtSum",groupname.Data());
	fHistManager.FillTH3(histname,JetPt,iR,PtSumA[ibin]);

	for(Int_t k=0; k<kNbinsInPt; k++){
	  if(JetPt>jetPtMinA[k] && JetPt<=jetPtMaxA[k]){
	    histname  = TString::Format("%s/fProDelRPtSum%s", groupname.Data(),namePtBinsA[k].Data());
	    fHistManager.FillProfile(histname,iR,PtSumA[ibin]);
	  }//if
	}//k loop
      }//ibin loop
      count++;
    }//auto : jet
  }//while jetCont 
}//FillJetShape()
//____________________________________________________________//

void AliAnalysisTaskEmcalJetProperties::FillJetShapeUE2(const Double_t alpha){
  //filling up the histograms
  const int kNbinsInPt = 21;
  Float_t jetPtMinA[kNbinsInPt] = {5.,
				   10., 15., 10., 
				   20., 25., 20., 
				   30., 35., 30.,30.,
				   40., 50., 40.,
				   60., 70., 60.,50.,
				   80., 90., 80};
  Float_t jetPtMaxA[kNbinsInPt] = {10., 
				   15., 20., 20.,
				   25., 30., 30.,
				   35., 40., 40.,50., 
				   50., 60., 60.,
				   70., 80., 80.,80., 
				   90., 100., 100.};
  TString namePtBinsA[kNbinsInPt];
  for(Int_t iii=0; iii<kNbinsInPt; iii++){
    namePtBinsA[iii] = Form("_JetPt%dto%d",(Int_t)jetPtMinA[iii],(Int_t)jetPtMaxA[iii]);
  }
  TString groupname;
  TString histname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    //cout << "#*#****#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*groupname" << groupname << endl;
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      if(count > 0)continue;//skipping all jets except leading jet
      Float_t JetPt  = (Float_t) jet->Pt();
      
      AliEmcalJet *leadJet = jetCont->GetLeadingJet();//seleting leading jet for consistency check
      Float_t leadingJetPt  = (Float_t) leadJet->Pt();
      if(leadingJetPt - JetPt) {Printf("%s:%d Mismatch in leading jet pT: %f %f",(char*)__FILE__,__LINE__,leadingJetPt,JetPt);return;}
      
      Double_t jetMom[3];
      jet->PxPyPz(jetMom);
      
      TVector3 jet3mom(jetMom);
      Double_t etaTilted = jet3mom.Eta();
      Double_t phiTilted = TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;
      if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();
            
      Float_t PtSumA[50]          = {0.};
      Float_t delRPtSum80pc       = 0;
      Float_t Pt_density[50]          = {0.};
      Float_t annu_area[50]          = {0.};
      //Int_t kNbinsR               = 10;
      
      Double_t sumPtPerp = 0;
      Int_t Ntracks_in_UE = 0;
      fTrackListUE->Clear();
    
      float fJetRadius = jetCont -> GetJetRadius();
      GetTracksTiltedwrpJetAxis(alpha,fTrackListUE,jet,fJetRadius,sumPtPerp,Ntracks_in_UE);
      
      histname  = TString::Format("%s/fh1PtSumInJetConeUE2",groupname.Data());
      fHistManager.FillTH1(histname,sumPtPerp);

      Int_t nTracksUE = fTrackListUE->GetEntries();
      //cout << "Inside FillJetShapeUE2, Jet pT: " << JetPt << " UE Ntracks: " << nTracksUE << endl;
      histname  = TString::Format("%s/fh2NtracksLeadingJetUE2",groupname.Data());
      fHistManager.FillTH2(histname,JetPt,nTracksUE);

      histname  = TString::Format("%s/fProNtracksLeadingJetUE2",groupname.Data());
      fHistManager.FillProfile(histname,JetPt,nTracksUE);
      
      Int_t   *index     = new Int_t   [nTracksUE];//dynamic array containing index
      Float_t *delRA     = new Float_t [nTracksUE];//dynamic array containing delR
      Float_t *delEtaA   = new Float_t [nTracksUE];//dynamic array containing delEta
      Float_t *delPhiA   = new Float_t [nTracksUE];//dynamic array containing delPhi
      Float_t *trackPtA  = new Float_t [nTracksUE];//dynamic array containing pt-track
      Float_t *trackEtaA = new Float_t [nTracksUE];//dynamic array containing eta-track
      Float_t *trackPhiA = new Float_t [nTracksUE];//dynamic array containing phi-track
      for(Int_t ii=0; ii<nTracksUE; ii++){
	index[ii]     = 0;
	delRA[ii]     = 0.;
	delEtaA[ii]   = 0.;
	delPhiA[ii]   = 0.;
	trackPtA[ii]  = 0.;
	trackEtaA[ii] = 0.;
	trackPhiA[ii] = 0.;
      }//ii loop
      
      for (Int_t j =0; j< nTracksUE; j++){
	Float_t TrackEta=-99.0; Float_t TrackPt=-99.0; Float_t TrackPhi=-99.0;
	Float_t DelEta=-99.0;  Float_t DelPhi=-99.0; Float_t DelR=-99.0;
	Float_t FF=-99.0;Float_t Ksi=-99.0;

	AliVParticle* trackkine = dynamic_cast<AliVParticle*>(fTrackListUE->At(j));
	if(!trackkine)continue;
	TrackEta = trackkine->Eta();
	TrackPhi = trackkine->Phi();
	TrackPt  = trackkine->Pt();
	
	if(JetPt)FF = TrackPt/JetPt;
	if(FF)Ksi   = TMath::Log(1./FF);
	
	histname = TString::Format("%s/fh2PtTrackUE2", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, TrackPt);
	histname = TString::Format("%s/fh2FFUE2", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, FF);
	histname = TString::Format("%s/fh2KsiUE2", groupname.Data());
	fHistManager.FillTH2(histname, JetPt, Ksi);

	DelEta = TMath::Abs(etaTilted - TrackEta);
	DelPhi = TMath::Abs(phiTilted - TrackPhi);
	if(TMath::Abs(DelPhi)>TMath::Pi())DelPhi = TMath::Abs(DelPhi-TMath::TwoPi());
	DelR   = TMath::Sqrt(DelEta*DelEta + DelPhi*DelPhi);

	delRA[j]     = DelR;
	delEtaA[j]   = DelEta;
	delPhiA[j]   = DelPhi;
	trackPtA[j]  = TrackPt;
	trackEtaA[j] = TrackEta;
	trackPhiA[j] = TrackPhi;
		
	for(Int_t ibin=1; ibin<=50; ibin++){
	  Float_t xlow = 0.02*(ibin-1);
	  Float_t xup  = 0.02*ibin;
	  if( xlow <= DelR && DelR < xup){
	    PtSumA[ibin-1]+= TrackPt;
	  }//if loop
	}//for ibin loop
      }//track loop

      for(Int_t ibin = 1; ibin <=50; ibin++)
	{
	  Float_t xlow = 0.02*(ibin-1);
	  Float_t xup  = 0.02*ibin;
	  annu_area[ibin-1] = TMath::Pi()*((xup*xup)-(xlow*xlow));
	  Pt_density[ibin-1] = PtSumA[ibin-1]/annu_area[ibin-1];
	}
      //cout << "############################## Check TList before clearing" << fTrackListUE->GetEntries() << endl;
      fTrackListUE->Clear();
      //cout << "############################## Check TList after clearing" << fTrackListUE->GetEntries() << endl;
      
      Float_t PtSum = 0;
      Bool_t iflagPtSum   = kFALSE;
      TMath::Sort(nTracksUE,delRA,index,0);
      for(Int_t ii=0; ii<nTracksUE; ii++){
	PtSum += trackPtA[index[ii]];

	if(!iflagPtSum){
	  if(PtSum/JetPt >= 0.8000){
	    delRPtSum80pc = delRA[index[ii]];
	    iflagPtSum = kTRUE;
	  }//if PtSum
	}//if iflag
      }//track loop 2nd
      delete [] index;
      delete [] delRA;
      delete [] delEtaA;
      delete [] delPhiA;
      delete [] trackPtA;
      delete [] trackEtaA;
      delete [] trackPhiA;
      
      histname  = TString::Format("%s/fh2DelR80pcPtUE2",groupname.Data());
      fHistManager.FillTH2(histname,JetPt,delRPtSum80pc);
          
      histname  = TString::Format("%s/fProDelR80pcPtUE2",groupname.Data());
      fHistManager.FillProfile(histname,JetPt,delRPtSum80pc);

      for(Int_t ibin=0; ibin<50; ibin++){
	Float_t iR = 0.02*ibin + 0.01;
	histname  = TString::Format("%s/fh3PtDelRPtSumUE2",groupname.Data());
	fHistManager.FillTH3(histname,JetPt,iR,PtSumA[ibin]);
	
	//UE2
	histname  = TString::Format("%s/fProDelRPtDenUEN2", groupname.Data());
	fHistManager.FillProfile(histname,iR,Pt_density[ibin]);

	for(Int_t k=0; k<kNbinsInPt; k++){
	  if(JetPt>jetPtMinA[k] && JetPt<=jetPtMaxA[k]){
	    histname  = TString::Format("%s/fProDelRPtSumUE2%s", groupname.Data(),namePtBinsA[k].Data());
	    fHistManager.FillProfile(histname,iR,PtSumA[ibin]);

	    //UE2
	    histname  = TString::Format("%s/fProDelRPtDenUE2%s", groupname.Data(),namePtBinsA[k].Data());
	    fHistManager.FillProfile(histname,iR,Pt_density[ibin]);
	  }//if
	}//k loop
      }//ibin loop
      count++;
    }//auto : jet
  }//while jetCont
}//FillJetShapeUE2()
//___________________________________________________________

void AliAnalysisTaskEmcalJetProperties::FillJetShapeUE3(TList* outputlist, Double_t EtaTilt, Double_t PhiTilt){

  //cout << "######################### Inside FillJetShapeUE3, the value of C factor is " << fCfactor << endl;
  //filling up the histograms
  const int kNbinsInPt = 21;
  Float_t jetPtMinA[kNbinsInPt] = {5.,
				   10., 15., 10., 
				   20., 25., 20., 
				   30., 35., 30.,30.,
				   40., 50., 40.,
				   60., 70., 60.,50.,
				   80., 90., 80};
  Float_t jetPtMaxA[kNbinsInPt] = {10., 
				   15., 20., 20.,
				   25., 30., 30.,
				   35., 40., 40.,50., 
				   50., 60., 60.,
				   70., 80., 80.,80., 
				   90., 100., 100.};
  TString namePtBinsA[kNbinsInPt];
  for(Int_t iii=0; iii<kNbinsInPt; iii++){
    namePtBinsA[iii] = Form("_JetPt%dto%d",(Int_t)jetPtMinA[iii],(Int_t)jetPtMaxA[iii]);
  }
  TString groupname;
  TString histname;

  AliJetContainer* jetCont = GetJetContainer(1);
  //TIter next(&fJetCollArray);
  //while ((jetCont = static_cast<AliJetContainer*>(next()))) {
  groupname = jetCont->GetName();
  UInt_t count = 0;
  for(auto jet : jetCont->accepted()) {
    if (!jet) continue;
    if(count > 0)continue;//skipping all jets except leading jet
    Float_t JetPt  = (Float_t) jet->Pt();
    
    AliEmcalJet *leadJet = jetCont->GetLeadingJet();//seleting leading jet for consistency check
    Float_t leadingJetPt  = (Float_t) leadJet->Pt();
    if(leadingJetPt - JetPt) {Printf("%s:%d Mismatch in leading jet pT: %f %f",(char*)__FILE__,__LINE__,leadingJetPt,JetPt);return;}    
    
    Double_t etaTilted = EtaTilt;//jet3mom.Eta();
    Double_t phiTilted = PhiTilt;//TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;
    
    Float_t PtSumA[50]          = {0.};
    Float_t delRPtSum80pc       = 0;
    Float_t Pt_density[50]          = {0.};
    Float_t annu_area[50]          = {0.};
    //Int_t kNbinsR               = 10;
    
    Double_t sumPtPerp = 0;
    Int_t Ntracks_in_UE = 0;
    
    histname  = TString::Format("%s/fh1PtSumInJetConeUE3",groupname.Data());
    fHistManager.FillTH1(histname,sumPtPerp);
    
    Int_t nTracksUE = outputlist->GetEntries();
    //cout << "groupname: " << groupname << "  nTracksUE: " << nTracksUE << endl; 
    
    histname  = TString::Format("%s/fh2NtracksLeadingJetUE3",groupname.Data());
    fHistManager.FillTH2(histname,JetPt,nTracksUE);
    
    histname  = TString::Format("%s/fProNtracksLeadingJetUE3",groupname.Data());
    fHistManager.FillProfile(histname,JetPt,nTracksUE);
    
    Int_t   *index     = new Int_t   [nTracksUE];//dynamic array containing index
    Float_t *delRA     = new Float_t [nTracksUE];//dynamic array containing delR
    Float_t *delEtaA   = new Float_t [nTracksUE];//dynamic array containing delEta
    Float_t *delPhiA   = new Float_t [nTracksUE];//dynamic array containing delPhi
    Float_t *trackPtA  = new Float_t [nTracksUE];//dynamic array containing pt-track
    Float_t *trackEtaA = new Float_t [nTracksUE];//dynamic array containing eta-track
    Float_t *trackPhiA = new Float_t [nTracksUE];//dynamic array containing phi-track
    for(Int_t ii=0; ii<nTracksUE; ii++){
      index[ii]     = 0;
      delRA[ii]     = 0.;
      delEtaA[ii]   = 0.;
      delPhiA[ii]   = 0.;
      trackPtA[ii]  = 0.;
      trackEtaA[ii] = 0.;
      trackPhiA[ii] = 0.;
    }//ii loop
    
    for (Int_t j =0; j< nTracksUE; j++){
      Float_t TrackEta=-99.0; Float_t TrackPt=-99.0; Float_t TrackPhi=-99.0;
      Float_t DelEta=-99.0;  Float_t DelPhi=-99.0; Float_t DelR=-99.0;
      Float_t FF=-99.0;Float_t Ksi=-99.0;
      
      AliVParticle* trackkine = dynamic_cast<AliVParticle*>(outputlist->At(j));
      if(!trackkine)continue;
      TrackEta = trackkine->Eta();
      TrackPhi = trackkine->Phi();
      TrackPt  = trackkine->Pt();
      //cout << "TrackPt: " << TrackPt << endl;
      
      if(JetPt)FF = TrackPt/JetPt;
      if(FF)Ksi   = TMath::Log(1./FF);
      
      histname = TString::Format("%s/fh2PtTrackUE3", groupname.Data());
      fHistManager.FillTH2(histname, JetPt, TrackPt);
      histname = TString::Format("%s/fh2FFUE3", groupname.Data());
      fHistManager.FillTH2(histname, JetPt, FF);
      histname = TString::Format("%s/fh2KsiUE3", groupname.Data());
      fHistManager.FillTH2(histname, JetPt, Ksi);
      
      DelEta = TMath::Abs(etaTilted - TrackEta);
      DelPhi = TMath::Abs(phiTilted - TrackPhi);
      if(TMath::Abs(DelPhi)>TMath::Pi())DelPhi = TMath::Abs(DelPhi-TMath::TwoPi());
      DelR   = TMath::Sqrt(DelEta*DelEta + DelPhi*DelPhi);
      
      delRA[j]     = DelR;
      delEtaA[j]   = DelEta;
      delPhiA[j]   = DelPhi;
      trackPtA[j]  = TrackPt;
      trackEtaA[j] = TrackEta;
      trackPhiA[j] = TrackPhi;
      
      for(Int_t ibin=1; ibin<=50; ibin++){
	Float_t xlow = 0.02*(ibin-1);
	Float_t xup  = 0.02*ibin;
	if( xlow <= DelR && DelR < xup){
	  PtSumA[ibin-1]+= TrackPt;
	}//if loop
      }//for ibin loop
    }//track loop
    
    for(Int_t ibin = 1; ibin <=50; ibin++)
      {
	Float_t xlow = 0.02*(ibin-1);
	Float_t xup  = 0.02*ibin;
	annu_area[ibin-1] = TMath::Pi()*((xup*xup)-(xlow*xlow));
	Pt_density[ibin-1] = PtSumA[ibin-1]/annu_area[ibin-1];
      }
    
    outputlist->Clear();
    
    Float_t PtSum = 0;
    Bool_t iflagPtSum   = kFALSE;
    TMath::Sort(nTracksUE,delRA,index,0);
    for(Int_t ii=0; ii<nTracksUE; ii++){
      PtSum += trackPtA[index[ii]];
      
      if(!iflagPtSum){
	if(PtSum/JetPt >= 0.8000){
	  delRPtSum80pc = delRA[index[ii]];
	  iflagPtSum = kTRUE;
	}//if PtSum
      }//if iflag
    }//track loop 2nd
    delete [] index;
    delete [] delRA;
    delete [] delEtaA;
    delete [] delPhiA;
    delete [] trackPtA;
    delete [] trackEtaA;
    delete [] trackPhiA;
    
    histname  = TString::Format("%s/fh2DelR80pcPtUE3",groupname.Data());
    fHistManager.FillTH2(histname,JetPt,delRPtSum80pc);
    
    histname  = TString::Format("%s/fProDelR80pcPtUE3",groupname.Data());
    fHistManager.FillProfile(histname,JetPt,delRPtSum80pc);
    
    for(Int_t ibin=0; ibin<50; ibin++){
      Float_t iR = 0.02*ibin + 0.01;
      histname  = TString::Format("%s/fh3PtDelRPtSumUE3",groupname.Data());
      fHistManager.FillTH3(histname,JetPt,iR,PtSumA[ibin]);
      
      //random
      histname  = TString::Format("%s/fProDelRPtDenUEN3", groupname.Data());
      fHistManager.FillProfile(histname,iR,Pt_density[ibin]);
      
      for(Int_t k=0; k<kNbinsInPt; k++){
	if(JetPt>jetPtMinA[k] && JetPt<=jetPtMaxA[k]){
	  histname  = TString::Format("%s/fProDelRPtSumUE3%s", groupname.Data(),namePtBinsA[k].Data());
	  fHistManager.FillProfile(histname,iR,PtSumA[ibin]);
	  
	  //random
	  histname  = TString::Format("%s/fProDelRPtDenUE3%s", groupname.Data(),namePtBinsA[k].Data());
	  fHistManager.FillProfile(histname,iR,Pt_density[ibin]);
	}//if
      }//k loop
    }//ibin loop
    count++;
  }//auto : jet
  //}//while jetCont
}//FillJetShapeUE3()
//___________________________________________________________

//
void AliAnalysisTaskEmcalJetProperties::RandomConeCalculation(Bool_t CheckOverlappingWithRC)
//void AliAnalysisTaskEmcalJetProperties::RandomConeCalculation()
{
  TString groupname;
  TString histname;
  AliJetContainer* myJetCont2 = GetJetContainer(1);
  if(!myJetCont2)
    {
      //cout << "No jet container (JetCont-1) found!" << endl;
      return;
    }
  TString mygroupname2;
  mygroupname2 = myJetCont2->GetName();
  //cout << "Debugging group: " << mygroupname2 << endl;
  float fJetRadius = myJetCont2 -> GetJetRadius();
  //cout << "Jet Radius: " << fJetRadius << endl;
  
  
  Int_t Leadingjetcount = 0;
  AliEmcalJet *jet1lead = new AliEmcalJet;
  AliEmcalJet *jet2lead = new AliEmcalJet;
  Double_t jetMom1[3]; Double_t jetMom2[3];
  Double_t Jet1Eta = 0.; Double_t Jet2Eta = 0.; Double_t Jet1Phi = 0.; Double_t Jet2Phi = 0.;
  for(auto testjetakt: myJetCont2->accepted())
    {
      //cout << testjetakt->Pt() << endl;
      Leadingjetcount++;
      if(Leadingjetcount == 1)
	{
	  jet1lead = (AliEmcalJet*)testjetakt;
	  jet1lead->PxPyPz(jetMom1);
	  TVector3 jet3mom1(jetMom1);
	  Jet1Eta = jet3mom1.Eta();
	  Jet1Phi = TVector2::Phi_0_2pi(jet3mom1.Phi());
	}
      if(Leadingjetcount == 2)
	{
	  jet2lead = (AliEmcalJet*)testjetakt;
	  jet2lead->PxPyPz(jetMom2);
	  TVector3 jet3mom2(jetMom2);
	  Jet2Eta = jet3mom2.Eta();
	  Jet2Phi = TVector2::Phi_0_2pi(jet3mom2.Phi());
	}
    }
  //cout << "B4->LeadJet: " << Leadingjetcount << endl;
  if(Leadingjetcount == 0)return;
  //cout << "After->LeadJet: " << Leadingjetcount << endl;
  
  const Int_t nRC = 10; // No. of random cones
  Double_t arrRC[nRC][2];
  for(Int_t iRC = 0; iRC < nRC; iRC++)
    {
      arrRC[iRC][0] = 99.;
      arrRC[iRC][1] = 99.;
    }

  TRandom3 randomGen(0);
  //TList* outputlist;
  //outputlist->Clear();
  fTrackListUE->Clear();
  int count = 0; // count random cones
  int i = 0;
  for(i = 0; i < 1000; i++)
    {
      if(count == nRC)break;
      
      Double_t tmpEta = 0.;
      Double_t tmpPhi = 0.;
      	
      tmpEta = randomGen.Uniform(-0.5,0.5);
      tmpPhi = randomGen.Uniform(0.,TMath::TwoPi());

      Double_t DelEta1 = TMath::Abs(Jet1Eta-tmpEta);
      Double_t DelPhi1 = TMath::Abs(Jet1Phi-tmpPhi);
      if(TMath::Abs(DelPhi1)>TMath::Pi())DelPhi1 = TMath::Abs(DelPhi1-TMath::TwoPi());
      Double_t DelR1 = TMath::Sqrt((DelEta1*DelEta1)+(DelPhi1*DelPhi1));
      if(DelR1 < (2.0*fJetRadius))continue;
      Double_t DelR2 = 0.;
      if(Leadingjetcount == 2)
	{
	  Double_t DelEta2 = TMath::Abs(Jet2Eta-tmpEta);
	  Double_t DelPhi2 = TMath::Abs(Jet2Phi-tmpPhi);
	  if(TMath::Abs(DelPhi2)>TMath::Pi())DelPhi2 = TMath::Abs(DelPhi2-TMath::TwoPi());
	  DelR2 = TMath::Sqrt((DelEta2*DelEta2)+(DelPhi2*DelPhi2));
	  if(DelR2 < (2.0*fJetRadius))continue;
	}
      //--------------------------------------------------
      //check overlapping with other random cones
      //cout << "CheckOverlappingWithRC: " << CheckOverlappingWithRC << endl;
      if(CheckOverlappingWithRC)
	{
	  Int_t OverlapRC = 0;
	  for(Int_t iRC = 0; iRC < nRC; iRC++)
	    {

	      Double_t DelEtaRC = TMath::Abs(arrRC[iRC][0]-tmpEta);
	      Double_t DelPhiRC = TMath::Abs(arrRC[iRC][1]-tmpPhi);
	      if(TMath::Abs(DelPhiRC)>TMath::Pi())DelPhiRC = TMath::Abs(DelPhiRC-TMath::TwoPi());
	      Double_t DelRRC = TMath::Sqrt((DelEtaRC*DelEtaRC)+(DelPhiRC*DelPhiRC));
	      if(DelRRC < (2.0*fJetRadius))OverlapRC++;
	      if(OverlapRC)break;
	    }
	  if(OverlapRC)continue;
	}//if(!CheckOverlappingWithRC)
      //--------------------------------------------------

      AliClusterContainer* clusCont = GetClusterContainer(0);
      AliParticleContainer* partCont = 0;
      TIter next(&fParticleCollArray);
      int countPart = 0;
      int iSkip = 0;
      
      while ((partCont = static_cast<AliParticleContainer*>(next()))) 
	{
	  for(auto part : partCont->accepted()) 
	    {
	      if (!part) continue;	  
	      Double_t trackMom[3];
	      
	      part->PxPyPz(trackMom);
	      TVector3 track3mom(trackMom);
	      
	      Double_t deta = track3mom.Eta() - tmpEta;
	      Double_t dphi = TMath::Abs(track3mom.Phi() - tmpPhi);
	      if (TMath::Abs(dphi) > TMath::Pi()) dphi = TMath::Abs(TMath::TwoPi() - dphi);
	      Double_t dR = TMath::Sqrt(deta * deta + dphi * dphi);
	      
	      if(dR<=fJetRadius)
		{
		  for(int itrk = 0; itrk < jet1lead->GetNumberOfTracks(); itrk++)
		    {
		      AliVParticle* trackJ = (AliVParticle*) jet1lead->TrackAt(itrk,partCont->GetArray());
		      if(part == trackJ)iSkip++;
		    }
		  if(Leadingjetcount == 2)
		    {
		      for(int itrk2 = 0; itrk2 < jet2lead->GetNumberOfTracks(); itrk2++)
			{
			  AliVParticle* trackJ2 = (AliVParticle*) jet2lead->TrackAt(itrk2,partCont->GetArray());
			  if(part == trackJ2)iSkip++;
			}
		    }
		  if(iSkip)
		    {
		      //cout << "matches with jet particle" << endl; 
		      break;
		    }
		  //cout << " Collecting particles" << endl;
		  fTrackListUE->Add(part);
		  countPart++;
		}// if(dR<=fJetRadius)
	      if(iSkip)break;
	    }//for(auto part : partCont->accepted()) 
	  if(iSkip)break;
	}//while ((partCont = static_cast<AliParticleContainer*>(next()))) 
      
      if(iSkip)
	{
	  fTrackListUE->Clear();
	  continue;
	}
      if(countPart != 0)
	{
	  count++;
	  //cout << "count RC: " << count << endl;	  
	  if(count <= nRC)
	    {
	      arrRC[count-1][0] = tmpEta;
	      arrRC[count-1][1] = tmpPhi;
	      FillJetShapeUE3(fTrackListUE,tmpEta,tmpPhi);
	      fTrackListUE->Clear();

	      histname  = TString::Format("%s/fh1DelRRCJet1",mygroupname2.Data());
	      fHistManager.FillTH1(histname,DelR1);
	      if(Leadingjetcount == 2)
		{
		  histname  = TString::Format("%s/fh1DelRRCJet2",mygroupname2.Data());
		  fHistManager.FillTH1(histname,DelR2);
		}
	      
	      histname  = TString::Format("%s/fh1nRC",mygroupname2.Data());
	      fHistManager.FillTH1(histname,1.);

	    }//if(count <= nRC)
	}//if(countPart != 0)
    }//for(int i = 0; i < 1000; i++)
  //cout << "i: " << i << endl;
  //cout << "count: " << count << endl;
  
}//RandomConeCalculation
//___________________________________________________________


inline Bool_t AliAnalysisTaskEmcalJetProperties::IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
    {
      Int_t jet1Track = jet1->TrackAt(i);
      for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
	{
	  Int_t jet2Track = jet2->TrackAt(j);
	  if (jet1Track == jet2Track)
	    return kTRUE;
	}
    }
  return kFALSE;
}
//___________________________________________________________

void AliAnalysisTaskEmcalJetProperties::GetTracksTiltedwrpJetAxis(Float_t alpha, TList* outputlist, const AliEmcalJet* jet, Double_t radius,Double_t& sumPt,Int_t& Ntracks_in_UE){
  //This part  is inherited from AliAnalysisTaskFragmentationFunction.cxx
  // List of tracks in cone perpendicular to the jet azimuthal direction
  Double_t jetMom[3];
  jet->PxPyPz(jetMom);

  TVector3 jet3mom(jetMom);
  // Rotate phi and keep eta unchanged
  Double_t etaTilted = jet3mom.Eta();//no change in eta
  Double_t phiTilted = TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;//rotate phi by alpha
  if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();
  
  AliClusterContainer* clusCont = GetClusterContainer(0);
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    /*
    if(fContainerKtRec == 0 && fContainerKtGen == 1 && fContainerAktRec == 2 && fContainerAktGen == 3)
      {
	if(IsRec && !strcmp(partCont->GetName(),"tracks"))cout << "RecParticleContainerName: " << partCont->GetName() << endl;
	else if(!IsRec && !strcmp(partCont->GetName(),"mcparticles"))cout << "GenParticleContainerName: " << partCont->GetName() << endl;
	else
	  continue;
	  }*/
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;
      
      Double_t trackMom[3];
      
      part->PxPyPz(trackMom);
      TVector3 track3mom(trackMom);
    
      Double_t deta = track3mom.Eta() - etaTilted;
      Double_t dphi = TMath::Abs(track3mom.Phi() - phiTilted);
      if (TMath::Abs(dphi) > TMath::Pi()) dphi = TMath::Abs(TMath::TwoPi() - dphi);
      Double_t dR = TMath::Sqrt(deta * deta + dphi * dphi);
      
      if(dR<=radius){
	outputlist->Add(part);
	sumPt += part->Pt();
	Ntracks_in_UE++;
	//if(fDebug > 100) printf("GetTracksTiltewrpJetAxis itrack, trackPt %d,  %f\n",itrack,track3mom.Pt());
      }//if dR<=radius
    }//part accepted
  }//while part cont
}//initial loop
//-----------------------------------------------//

Double_t AliAnalysisTaskEmcalJetProperties::BkgSubtracted(Int_t KtJetCont, Int_t AktJetCont)
{
  Double_t Nch_subtracted = 0;  

  const int kNbinsInPt = 21;
  Float_t jetPtMinA[kNbinsInPt] = {5.,
				   10., 15., 10., 
				   20., 25., 20., 
				   30., 35., 30.,30.,
				   40., 50., 40.,
				   60., 70., 60.,50.,
				   80., 90., 80};
  Float_t jetPtMaxA[kNbinsInPt] = {10., 
				   15., 20., 20.,
				   25., 30., 30.,
				   35., 40., 40.,50., 
				   50., 60., 60.,
				   70., 80., 80.,80., 
				   90., 100., 100.};
  TString namePtBinsA[kNbinsInPt];
  for(Int_t iii=0; iii<kNbinsInPt; iii++){
    namePtBinsA[iii] = Form("_JetPt%dto%d",(Int_t)jetPtMinA[iii],(Int_t)jetPtMaxA[iii]);
  }//iii loop

  TString histname;
  TString histtitle;
  TString groupname;
  TString mygroupname;

  AliJetContainer* jetCont = GetJetContainer(AktJetCont);
  groupname = jetCont->GetName();
  
  Int_t Leadingjetcount = 0;
  AliEmcalJet *Jet1lead = NULL;
  AliEmcalJet *Jet2lead = NULL;
  
  for(auto testjetakt: jetCont->accepted())
    {
      //cout << "************************************************* Anti Jet eta " << testjetakt->Eta() << endl;
      Leadingjetcount++;
      if(Leadingjetcount == 1)
	{
	  Jet1lead = (AliEmcalJet*)testjetakt;
	  //cout << "First leading jet: " << Jet1lead->Pt() << endl;
	}
      if(Leadingjetcount == 2)
	{
	  Jet2lead = (AliEmcalJet*)testjetakt;
	  //cout << "Second leading jet: " << Jet2lead->Pt() << endl;
	}
    }
  Double_t JetEta = 0.0; 
  Double_t JetPhi = 0.0; 
  Double_t JetPt = 0.0; 
  Double_t Jet_Pt = 0.0;
  Double_t Jet_Area = 0.0;
  Double_t N_Jet_Tracks = 0;
  Double_t Leading_Jet_Pt = 0.0;
  Double_t N_Leading_Jet_Tracks = 0;
  Double_t rho = 0.0;
  Double_t rhoN = 0.0;
  Double_t tmpSumArea = 0.0;
  Double_t tmpCovArea = 0.0;
  Int_t tmpSumNo = 0;
  Int_t tmpCovNo = 0;
  Double_t Cfactor = 0.0;
  Double_t Cfactornew = 0.0;
  UInt_t count = 0;
  for(auto jet : jetCont->accepted())
    {
      if (!jet) continue;
      
      if(count == 0)
	{
	  Leading_Jet_Pt = jet->Pt();
	  N_Leading_Jet_Tracks = jet->GetNumberOfTracks();
	}
      AliParticleContainer* tracks = jetCont->GetParticleContainer();
      JetEta = jet->Eta();
      JetPhi = jet->Phi();
      JetPt  = jet->Pt() ;
      Jet_Pt = jet->Pt();
      Jet_Area = jet->Area();
      N_Jet_Tracks = jet->GetNumberOfTracks();
      Double_t pt_jet_ue_subt = 0.0;  
      pt_jet_ue_subt = Jet_Pt;
      
      Double_t pt_leading_jet_ue_subt = 0.0;  
      pt_leading_jet_ue_subt = Leading_Jet_Pt;
      
      Double_t pt_leading_jet_ue_subt_C = 0.0;  
      pt_leading_jet_ue_subt_C = Leading_Jet_Pt;
      
      Double_t n_jet_tracks_ue_subt = 0;  
      n_jet_tracks_ue_subt = N_Jet_Tracks;
      
      Double_t n_jet_tracks_ue_subt_C = 0;  
      n_jet_tracks_ue_subt_C = N_Jet_Tracks;
      
      Double_t n_leading_jet_ue = 0;  
      Double_t n_leading_jet_tracks_ue_subt = 0;  
      n_leading_jet_tracks_ue_subt = N_Leading_Jet_Tracks;
      Double_t n_leading_jet_ue_C = 0;  
      Double_t n_leading_jet_tracks_ue_subt_C = 0;  
      n_leading_jet_tracks_ue_subt_C = N_Leading_Jet_Tracks;
      
      if(count == 0)
	{
	  rho = 0.0;
	  rhoN = 0.0;
	  tmpSumArea = 0.0;
	  tmpCovArea = 0.0;
	  tmpSumNo = 0;
	  tmpCovNo = 0;
	  Cfactor = 0.0;
	  Cfactornew = 0.0;

	  Float_t PtSumA_v2[100][50] = {0.};


	  AliJetContainer* myjetCont = GetJetContainer(KtJetCont);
	  if(!myjetCont)return N_Leading_Jet_Tracks;
	  mygroupname = myjetCont->GetName();

	  if(!myjetCont->GetArrayName().IsNull())
	    {
	      fJets = myjetCont->GetArray();
	            
	      const Int_t Njets = fJets->GetEntries();
	      //cout << "No. of jets: " << Njets << endl;
	            
	      Int_t maxJetIds[]   = {-1, -1};
	      Float_t maxJetPts[] = { 0,  0};
	      if (fNExclLeadJets > 0)
		{
		  for (Int_t ij = 0; ij < Njets; ++ij)
		    {
		      AliEmcalJet *jetkt = static_cast<AliEmcalJet*>(fJets->At(ij));
		            
		      if (!jetkt)
			{
			  AliError(Form("%s: Could not receive jet %d", GetName(), ij));
			  continue;
			} 
		            
		      if(jetkt->Eta() > 0.5 || jetkt->Eta() < -0.5)continue;
		            
		      //if (!AcceptJet(jetkt))continue;
		            
		      if (jetkt->Pt() > maxJetPts[0])
			{
			  maxJetPts[1] = maxJetPts[0];
			  maxJetIds[1] = maxJetIds[0];
			  maxJetPts[0] = jetkt->Pt();
			  maxJetIds[0] = ij;
			}
		      else if (jetkt->Pt() > maxJetPts[1])
			{
			  maxJetPts[1] = jetkt->Pt();
			  maxJetIds[1] = ij;
			}
		    }
		  if (fNExclLeadJets < 2)
		    {
		      maxJetIds[1] = -1;
		      maxJetPts[1] = 0;
		    }
		}//leading jets to be excluded
	            
	      static Double_t rhovec[999];
	      static Double_t rhoNvec[999];
	      Int_t NjetAcc = 0;
	      // push all jets within selected acceptance into stack
	      for (Int_t iJets = 0; iJets < Njets; ++iJets)
		{
		  AliEmcalJet *jetkt = static_cast<AliEmcalJet*>(fJets->At(iJets));
		  if (!jetkt)
		    {
		      AliError(Form("%s: Could not receive jet %d", GetName(), iJets));
		      continue;
		    } 
		  tmpSumArea += jetkt->Area();
		  if(jetkt->Pt() > 0.15)
		    {
		      tmpCovArea += jetkt->Area();
		      //cout << "jet-kt area: " << jetkt->Area() << "  sum: " << tmpCovArea << endl;
		    }
		  if(jetkt->Eta() > 0.5 || jetkt->Eta() < -0.5)continue;
		  if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
		    {
		      //cout << "iJets: " << iJets << " Jet pT: " << jetkt->Pt() << " area: " << jetkt->Area() << " Ntracks: " << jetkt->GetNumberOfTracks() << " Nconsts: " << jetkt->GetNumberOfConstituents() << endl;
		      continue;
		    }
		  if(jetkt->Pt() > 0.15)
		    {
		      if(Leadingjetcount == 1)
			{
			  if(IsJetOverlapping(Jet1lead,jetkt))
			    {
			      continue;
			    }
			}
		      else if(Leadingjetcount >= 2)
			{
			  if(IsJetOverlapping(Jet1lead,jetkt) || IsJetOverlapping(Jet2lead,jetkt))
			    {
			      continue;
			    }
			}
		      rhovec[NjetAcc] = jetkt->Pt() / jetkt->Area();
		      rhoNvec[NjetAcc] = jetkt->GetNumberOfTracks() / jetkt->Area();
		      ++NjetAcc;    
		      
		      
		      /***************************************** dpTsum/dr starts **************************************/
		      
		      Double_t Del_eta = 0.;
		      Double_t Del_phi = 0.;
		      Double_t Del_R = 0.;
		      
		      Int_t N_Leading_Jet_v2 = 0;
		      N_Leading_Jet_v2 = jetkt->GetNumberOfTracks();
		      
		      Float_t JetEta_v2; Float_t JetPhi_v2; Float_t JetPt_v2;
		      JetEta_v2 = (Float_t)jetkt->Eta();
		      JetPhi_v2 = (Float_t)jetkt->Phi();
		      JetPt_v2  = (Float_t)jetkt->Pt();
		      //cout << "jet eta2: " << JetEta_v2 << " jet phi2: " << JetPhi_v2 << " jet pt2: " << JetPt_v2 << endl;
		      
		      Del_eta = JetEta_v2 - JetEta;
		      Del_phi = JetPhi_v2 - JetPhi;
		      if(TMath::Abs(Del_phi)>TMath::Pi())Del_phi = TMath::Abs(Del_phi-TMath::TwoPi());
		      Del_R = TMath::Sqrt(Del_eta*Del_eta + Del_phi*Del_phi);
		      
		      Int_t   *index_v2     = new Int_t   [N_Leading_Jet_v2];//dynamic array containing index
		      Float_t *delRA_v2     = new Float_t [N_Leading_Jet_v2];//dynamic array containing delR
		      Float_t *delEtaA_v2   = new Float_t [N_Leading_Jet_v2];//dynamic array containing delEta
		      Float_t *delPhiA_v2   = new Float_t [N_Leading_Jet_v2];//dynamic array containing delPhi
		      Float_t *trackPtA_v2  = new Float_t [N_Leading_Jet_v2];//dynamic array containing pt-track
		      Float_t *trackEtaA_v2 = new Float_t [N_Leading_Jet_v2];//dynamic array containing eta-track
		      Float_t *trackPhiA_v2 = new Float_t [N_Leading_Jet_v2];//dynamic array containing phi-track
		      for(Int_t ii=0; ii<N_Leading_Jet_v2; ii++){
			index_v2[ii]     = 0;
			delRA_v2[ii]     = 0.;
			delEtaA_v2[ii]   = 0.;
			delPhiA_v2[ii]   = 0.;
			trackPtA_v2[ii]  = 0.;
			trackEtaA_v2[ii] = 0.;
			trackPhiA_v2[ii] = 0.;
		      }//ii loop
		      
		      AliParticleContainer* tracks_v2 = myjetCont->GetParticleContainer();
		      //if(!tracks_v2)continue;
		      //if(!tracks_v2)cout << "Inside: Particle Container Missing." << endl;
		      
		      for (Int_t j =0; j< jetkt->GetNumberOfTracks(); j++){
			Float_t TrackEta_v2=-99.0; Float_t TrackPt_v2=-99.0; Float_t TrackPhi_v2=-99.0;
			Float_t DelEta_v2=-99.0;  Float_t DelPhi_v2=-99.0; Float_t DelR_v2=-99.0; Float_t AreaJ_v2=-99.0;
			AliVParticle* trackkine_v2 = (AliVParticle*) jetkt->TrackAt(j,tracks_v2->GetArray());

			if(!trackkine_v2)continue;
			TrackEta_v2 = trackkine_v2->Eta();
			TrackPhi_v2 = trackkine_v2->Phi();
			TrackPt_v2  = trackkine_v2->Pt();

			DelEta_v2 = TMath::Abs(JetEta_v2 - TrackEta_v2);
			DelPhi_v2 = TMath::Abs(JetPhi_v2 - TrackPhi_v2);
			if(TMath::Abs(DelPhi_v2)>TMath::Pi())DelPhi_v2 = TMath::Abs(DelPhi_v2-TMath::TwoPi());
			DelR_v2   = TMath::Sqrt(DelEta_v2*DelEta_v2 + DelPhi_v2*DelPhi_v2);			
			
			delRA_v2[j]     = DelR_v2;
			delEtaA_v2[j]   = DelEta_v2;
			delPhiA_v2[j]   = DelPhi_v2;
			trackPtA_v2[j]  = TrackPt_v2;
			trackEtaA_v2[j] = TrackEta_v2;
			trackPhiA_v2[j] = TrackPhi_v2;
			
			for(Int_t ibin=1; ibin<=50; ibin++){
			  Float_t xlow = 0.02*(ibin-1);
			  Float_t xup  = 0.02*ibin;
			  if( xlow <= DelR_v2 && DelR_v2 < xup){
			    PtSumA_v2[NjetAcc-1][ibin-1]+= TrackPt_v2;
			  }//if loop
			}//for ibin loop
		      }//track loop
		      
		      Double_t annu_area[50] = {0.};
		      Double_t Pt_density[50] = {0.};

		      for(Int_t ibin = 0; ibin < 50; ibin++)
			{
			  Float_t xlow = 0.02*(ibin);
			  Float_t xup  = 0.02*(ibin+1);
			  annu_area[ibin] = TMath::Pi()*((xup*xup)-(xlow*xlow));
			  Pt_density[ibin] = PtSumA_v2[NjetAcc-1][ibin]/annu_area[ibin];
			      
			  Float_t iR_v2 = 0.02*ibin + 0.01;
			      
			  histname  = TString::Format("%s/hProDelRPtDenUEN", groupname.Data());
			  fHistManager.FillProfile(histname,iR_v2,Pt_density[ibin]);
			      
			  for(Int_t k=0; k<kNbinsInPt; k++){
			    if(JetPt>jetPtMinA[k] && JetPt<=jetPtMaxA[k]){
			      histname  = TString::Format("%s/hhProDelRPtSumUE%s", groupname.Data(),namePtBinsA[k].Data());
			      fHistManager.FillProfile(histname,iR_v2,PtSumA_v2[NjetAcc-1][ibin]);
			      histname  = TString::Format("%s/hhProDelRPtDenUE%s", groupname.Data(),namePtBinsA[k].Data());
			      fHistManager.FillProfile(histname,iR_v2,Pt_density[ibin]);
			    }//if
			  }//k loop
			}
		    }// if(jetkt->Pt() > 0.15)
		  else
		    {
		      //cout << "Debug: Jet pt: " << jetkt->Pt() << endl;
		    }
		}//jetkt loop
	      if (NjetAcc > 0)
		{
		  //find median value
		  rho = TMath::Median(NjetAcc, rhovec);
		  rhoN = TMath::Median(NjetAcc, rhoNvec);
		  histname = TString::Format("%s/fh1rho", mygroupname.Data());
		  fHistManager.FillTH1(histname, rho);
		  histname = TString::Format("%s/fh1rhoN", mygroupname.Data());
		  fHistManager.FillTH1(histname, rhoN);
		}//if(NjetAcc > 0)
	    }//if jet container has array
	  pt_leading_jet_ue_subt = Leading_Jet_Pt - (rho * Jet_Area);
	  //cout << "Leading_Jet_Pt: " << Leading_Jet_Pt << "\t rho: " << rho << "\t Jet_Area: " << Jet_Area << "\t pt_leading_jet_ue_subt: " << pt_leading_jet_ue_subt << endl;
	  n_leading_jet_ue = rhoN * Jet_Area;
	  n_leading_jet_tracks_ue_subt = N_Leading_Jet_Tracks - (rhoN * Jet_Area);
	  //cout << "N_Leading_Jet_Tracks: " << N_Leading_Jet_Tracks << "\t rhoN: " << rhoN << "\t Jet_Area: " << Jet_Area << "\t n_leading_jet_tracks_ue_subt: " << n_leading_jet_tracks_ue_subt << endl;
	  
	  histname = TString::Format("%s/h1LeadingJetPtWUE", groupname.Data());
	  fHistManager.FillTH1(histname, Leading_Jet_Pt);
	      
	  histname = TString::Format("%s/h1LeadingJetPtWoUEWoCAll", groupname.Data());
	  fHistManager.FillTH1(histname, pt_leading_jet_ue_subt);
	  
	  histname = TString::Format("%s/fProNtracksLeadingJetWUE", groupname.Data());
	  fHistManager.FillProfile(histname, Leading_Jet_Pt, N_Leading_Jet_Tracks);
	  
	  histname = TString::Format("%s/fProNtracksLeadingJetUEWoCAll", groupname.Data());
	  fHistManager.FillProfile(histname, Leading_Jet_Pt, n_leading_jet_ue);
	  
	  histname = TString::Format("%s/fProNtracksLeadingJetWoUEWoCAll", groupname.Data());
	  fHistManager.FillProfile(histname, Leading_Jet_Pt, n_leading_jet_tracks_ue_subt);

	  histname = TString::Format("%s/fProNtracksLeadingJetUEwrtSubtPtWoC", groupname.Data());
	  fHistManager.FillProfile(histname, pt_leading_jet_ue_subt, n_leading_jet_ue);

	  histname = TString::Format("%s/fProNtracksLeadingJetWoUEwrtSubtPtWoC", groupname.Data());
	  fHistManager.FillProfile(histname, pt_leading_jet_ue_subt, n_leading_jet_tracks_ue_subt);

	  Cfactornew = (1.0*tmpCovArea)/(1.0*tmpSumArea);

	  fCfactor = Cfactornew;

	  //cout << "fCfactor: " << fCfactor << endl;
	  
	  histname = TString::Format("%s/h1Cfactor", groupname.Data());
	  fHistManager.FillTH1(histname, Cfactornew);

	  rho = rho*Cfactornew;
	  rhoN = rhoN*Cfactornew;
	  pt_leading_jet_ue_subt_C = Leading_Jet_Pt - (rho * Jet_Area);
	  //cout << "Leading_Jet_Pt: " << Leading_Jet_Pt << "\t rho: " << rho << "\t Jet_Area: " << Jet_Area << "\t pt_leading_jet_ue_subt_C: " << pt_leading_jet_ue_subt_C << endl;
	    
	  n_leading_jet_ue_C = rhoN * Jet_Area;
	  n_leading_jet_tracks_ue_subt_C = N_Leading_Jet_Tracks - (rhoN * Jet_Area);
	  //cout << "N_Leading_Jet_Tracks: " << N_Leading_Jet_Tracks << "\t rhoN: " << rhoN << "\t Jet_Area: " << Jet_Area << "\t n_leading_jet_tracks_ue_subt_C: " << n_leading_jet_tracks_ue_subt_C << endl;

	  histname = TString::Format("%s/h1LeadingJetPtWoUEWCAll", groupname.Data());
	  fHistManager.FillTH1(histname, pt_leading_jet_ue_subt_C);
	      
	  histname = TString::Format("%s/fProNtracksLeadingJetUEWCAll", groupname.Data());
	  fHistManager.FillProfile(histname, Leading_Jet_Pt, n_leading_jet_ue_C);
	      
	  histname = TString::Format("%s/fProNtracksLeadingJetWoUEWCAll", groupname.Data());
	  fHistManager.FillProfile(histname, Leading_Jet_Pt, n_leading_jet_tracks_ue_subt_C);

	  histname = TString::Format("%s/fProNtracksLeadingJetUEwrtSubtPtWC", groupname.Data());
	  fHistManager.FillProfile(histname, pt_leading_jet_ue_subt_C, n_leading_jet_ue_C);
	              
	  histname = TString::Format("%s/fProNtracksLeadingJetWoUEwrtSubtPtWC", groupname.Data());
	  fHistManager.FillProfile(histname, pt_leading_jet_ue_subt_C, n_leading_jet_tracks_ue_subt_C);
	  //cout << "count of leading jet: " << (count+1) << endl;    
	  Nch_subtracted = n_leading_jet_tracks_ue_subt_C;  
	}//if(count == 0)
      count++;  
    }// for(auto jet : jetCont->accepted())
  
  return Nch_subtracted;
}//BkgSubtracted

//This function does matching between Rec and Gen jets ----- taken from AliAnalysisTaskEmcalJettagger.cxx
//void  MatchJetsGeo(Int_t c1 = -1, Int_t c2 = -1, Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t type = 2, Bool_t bReset = kTRUE);
//void AliAnalysisTaskEmcalJetProperties::MatchJets(Int_t c1, Int_t c2)
void AliAnalysisTaskEmcalJetProperties::MatchJets(Int_t c1, Int_t c2, Double_t maxDist, TList* matchedListJet1, TList* matchedListJet2)
{
  Bool_t bReset = kTRUE;
  //if(c1<0) c1 = fContainerBase;
  //if(c2<0) c2 = fContainerTag;
  //Init();
  const Int_t nJets11 = GetNJets(c1);
  const Int_t nJets12 = GetNJets(c2);
  //cout << "Number of jets: 1-> " << nJets11 << " 2-> " << nJets12 << endl;
  if(nJets11==0 || nJets12==0) return;
  
  //AliDebugStream(1) << "Jets Base (" << GetJetContainer(c1)->GetNJets() << ", accepted " << GetJetContainer(c1)->GetNAcceptedJets() << "), jets tag(" << GetJetContainer(c2)->GetNJets() << ", accepted " << GetJetContainer(c2)->GetNAcceptedJets() << "), max distance " << maxDist << std::endl;
  
  if(bReset) {
    for(int i = 0;i<GetNJets(c1);i++){
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(GetJetFromArray(i, c1));
      if(!jet) continue;
      jet->ResetMatching();
    }
    for(int i = 0;i<GetNJets(c2);i++){
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(GetJetFromArray(i, c2));
      if(!jet) continue;
      jet->ResetMatching();
    }
    //ResetTagging(c1);
    //ResetTagging(c2);
  }
  //fMatchingDone = kFALSE;
  
  TArrayI faMatchIndex1;
  faMatchIndex1.Set(nJets12+1);
  faMatchIndex1.Reset(-1);
  
  TArrayI faMatchIndex2;
  faMatchIndex2.Set(nJets11+1);
  faMatchIndex2.Reset(-1);
  
  const Int_t nJets21 = GetNJets(c1);
  const Int_t nJets22 = GetNJets(c2);

  TArrayS iFlag((nJets11*nJets12));
  if(iFlag.GetSize()<(nJets11*nJets12)){
    iFlag.Set(nJets11*nJets12+1);
  }
  iFlag.Reset(0);
  
  // find the closest distance to the full jet
  for(int i = 0;i<nJets11;i++){
    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, c1));
    if(!jet1) continue;
    
    Float_t dist = maxDist;
    
    for(int j = 0;j <nJets12; j++){
      AliEmcalJet *jet2 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(j, c2));
      if(!jet2) continue;
      
      Double_t dR = jet1->DeltaR(jet2);
      if(dR<dist && dR<maxDist){
	faMatchIndex2[i]=j;
	dist = dR;
      }
    }//j jet loop
    if(faMatchIndex2[i]>=0) {
      iFlag[i*nJets12+faMatchIndex2[i]]+=1;//j closest to i
      //if(iDebug>10) 
      //Printf("Full Distance (%d)--(%d) %3.3f flag[%d] = %d",i,faMatchIndex2[i],dist,i*nJets12+faMatchIndex2[i],iFlag[i*nJets12+faMatchIndex2[i]]);
    }
  }//i jet loop
  
  // other way around
  for(int j = 0;j<nJets12;j++){
    AliEmcalJet *jet2 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(j, c2));
    if(!jet2)
      continue;
    
    Float_t dist = maxDist;
    for(int i = 0;i<nJets11;i++){
      AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, c1));
      if(!jet1)continue;
      
      Double_t dR = jet1->DeltaR(jet2);
      if(dR<dist && dR<maxDist){
	faMatchIndex1[j]=i;
        dist = dR;
      }   
    }
    if(faMatchIndex1[j]>=0) {
      iFlag[faMatchIndex1[j]*nJets12+j]+=2;//i closest to j
      //if(iDebug>10) 
      //Printf("Other way Distance (%d)--(%d) %3.3f flag[%d] = %d",faMatchIndex1[j],j,dist,faMatchIndex1[j]*nJets12+j,iFlag[faMatchIndex1[j]*nJets12+j]);
    }
  }
  Int_t IMatchPair = 0;
  // check for "true" correlations
  for(int i = 0;i<nJets11;i++){
    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetJetFromArray(i, c1));
    for(int j = 0;j<nJets12;j++){
      AliEmcalJet *jet2 = static_cast<AliEmcalJet*>(GetJetFromArray(j, c2));
      AliDebug(11,Form("%s: Flag[%d][%d] %d ",GetName(),i,j,iFlag[i*nJets12+j]));
      
      // we have a uniqe correlation
      if(iFlag[i*nJets12+j]==3) {
        Double_t dR = jet1->DeltaR(jet2);
	jet1->SetClosestJet(jet2,dR);
	jet2->SetClosestJet(jet1,dR);

	matchedListJet1 -> AddAt((AliEmcalJet*)jet1,IMatchPair);
	matchedListJet2 -> AddAt((AliEmcalJet*)jet2,IMatchPair);
	IMatchPair++;
      }
    }
  }

  //fMatchingDone = kTRUE;
}

//Fill 2D response matrices
void AliAnalysisTaskEmcalJetProperties::FillMy2DResponseMatrices()
{
  //cout << "-----------------------I am inside the response matrix filling function.-----------------------" << endl;

  matchedRecJetList->Clear();
  matchedGenJetList->Clear();

  MatchJets(3,2,fMaxDist,matchedGenJetList,matchedRecJetList);

  AliEmcalJet  *jet             = NULL; //jet pointer real jet
  AliEmcalJet  *jetDetMC        = NULL; //jet pointer detector level MC jet
  AliVParticle *track_PartLevel = NULL; //jet constituent
  AliVParticle *track_DetLevel  = NULL; //mc particle

  Double_t fArray_Nch[4] = {0.};
  Double_t fArray_NchUEPC[4] = {0.};
  Double_t fArray_PtSum[5] = {0.};
  Double_t fArray_PtSumUE[5] = {0.};
  Double_t fArray_FF[4] = {0.};
  Double_t fArray_FFM[4] = {0.};
  Double_t fArray_FFMUE[4] = {0.};

  Int_t mIndex = GetMLeadingJet();
  
  Int_t nMatchedJetsCheck = 0;
  nMatchedJetsCheck = matchedRecJetList->GetEntries();  

  Double_t N_Leading_Jet_Tracks = 0;
  Int_t LabelGen[100] = {0}, LabelRec[100] = {0};
  Int_t LabelUEGen[100] = {0}, LabelUERec[100] = {0};

  if(nMatchedJetsCheck)
    {
      jet = (AliEmcalJet*)matchedGenJetList->At(mIndex);
      jetDetMC = (AliEmcalJet*)matchedRecJetList->At(mIndex);

      h2ResJetPt->Fill(jetDetMC->Pt(),jet->Pt());
      h2ResJetPtNch->Fill(jetDetMC->Pt(),jet->Pt());

      
      Int_t nGenJetTracks = jet->GetNumberOfTracks();
      Int_t nRecJetTracks = jetDetMC->GetNumberOfTracks();      

      AliJetContainer *jetContGen = GetJetContainer(fContainerAktGen);
      AliJetContainer *jetContRec = GetJetContainer(fContainerAktRec);
      AliParticleContainer* tracksGen = jetContGen->GetParticleContainer();
      AliParticleContainer* tracksRec = jetContRec->GetParticleContainer();

      TList *RecUETrkList = new TList();
      TList *GenUETrkList = new TList();

      const Double_t alpha = TMath::Pi()/2.;
      float fJetRadius = 0.4;

      Double_t jetMomGen[3];
      jet->PxPyPz(jetMomGen);
      TVector3 jet3momGen(jetMomGen);
      Double_t etaTiltedGen = jet3momGen.Eta();
      Double_t phiTiltedGen = TVector2::Phi_0_2pi(jet3momGen.Phi()) + alpha;
      if(phiTiltedGen > 2*TMath::Pi()) phiTiltedGen = phiTiltedGen - 2*TMath::Pi();
      
      Double_t sumPtPerpGen = 0.;
      Int_t Ntracks_in_UE_Gen = 0;
      
      GenUETrkList->Clear();  
      IsRec = kFALSE;
      GetTracksTiltedwrpJetAxis(alpha,GenUETrkList,jet,fJetRadius,sumPtPerpGen,Ntracks_in_UE_Gen);

      Double_t jetMomRec[3];
      jetDetMC->PxPyPz(jetMomRec);
      TVector3 jet3momRec(jetMomRec);
      Double_t etaTiltedRec = jet3momRec.Eta();
      Double_t phiTiltedRec = TVector2::Phi_0_2pi(jet3momRec.Phi()) + alpha;
      if(phiTiltedRec > 2*TMath::Pi()) phiTiltedRec = phiTiltedRec - 2*TMath::Pi();
      
      Double_t sumPtPerpRec = 0.;
      Int_t Ntracks_in_UE_Rec = 0;
      
      RecUETrkList->Clear();  
      IsRec = kTRUE;
      GetTracksTiltedwrpJetAxis(alpha,RecUETrkList,jetDetMC,fJetRadius,sumPtPerpRec,Ntracks_in_UE_Rec);
      

      for(Int_t it = 0; it < Ntracks_in_UE_Gen; it++)
	{
	  AliVParticle* trackkineUEGen = dynamic_cast<AliVParticle*>(GenUETrkList->At(it));
	  LabelUEGen[it] = trackkineUEGen->GetLabel();
	}
      for(Int_t jt = 0; jt < Ntracks_in_UE_Rec; jt++)
	{
	  AliVParticle* trackkineUERec = dynamic_cast<AliVParticle*>(RecUETrkList->At(jt));
	  LabelUERec[jt] = trackkineUERec->GetLabel();
	}

      Int_t matchedUETrackIndex[100] = {0}, matchedUEMCIndex[100] = {0};
      for(Int_t iarr = 0; iarr < 100; iarr++){matchedUEMCIndex[iarr] = -99; matchedUETrackIndex[iarr] = -99;}
      Int_t nMatchedUETracks = 0;
      for(Int_t iTrkIndex = 0; iTrkIndex < Ntracks_in_UE_Rec; iTrkIndex++)
	{
	  for(Int_t iMCIndex = 0; iMCIndex < Ntracks_in_UE_Gen; iMCIndex++)
	    {
	      if(LabelUEGen[iMCIndex] == LabelUERec[iTrkIndex])
		{
		  matchedUETrackIndex[nMatchedUETracks] = iTrkIndex;
		  matchedUEMCIndex[nMatchedUETracks] = iMCIndex;
		  //cout << "Matched-> matchedUETrackIndex: " << matchedUETrackIndex[nMatchedUETracks] << "  matchedUEMCIndex[nMatchedUETracks]: " << matchedUEMCIndex[nMatchedUETracks] << "  nMatchedUETracks: " << nMatchedUETracks << endl;
		  nMatchedUETracks++;
		}
	    }
	}
      
      for(Int_t jt = 0; jt < Ntracks_in_UE_Rec; jt++)
	{
	  AliVParticle* trackkineUERec = dynamic_cast<AliVParticle*>(RecUETrkList->At(jt));
	  Int_t MatchedUERecTrack = 0;
	  for(Int_t iMIdx = 0; iMIdx < nMatchedUETracks; iMIdx++)
	    {
	      //cout << "Rec UE jt: " << jt << endl;
	      if(jt == matchedUETrackIndex[iMIdx])
		{
		  //cout << "Rec UE jt: " << jt << " matchedUETrackIndex: " << matchedUETrackIndex[iMIdx] << " matchedUEMCIndex: " << matchedUEMCIndex[iMIdx] << endl;
		  fArray_FFMUE[0] = jetDetMC->Pt();
		  fArray_FFMUE[2] = jet->Pt();
		  fArray_FFMUE[1] = GetFFUEPC(kTRUE,jetDetMC,RecUETrkList,matchedUETrackIndex[iMIdx]); 
		  fArray_FFMUE[3] = GetFFUEPC(kFALSE,jet,GenUETrkList,matchedUEMCIndex[iMIdx]);
		  MatchedUERecTrack++;
		  
		  //cout << fArray_FFMUE[0] << " " << fArray_FFMUE[1] << " " << fArray_FFMUE[2] << " " << fArray_FFMUE[3] << endl;
		  //h4
		  h4ResFFMUE->Fill(fArray_FFMUE);
		  h2FFMUERec->Fill(fArray_FFMUE[0],fArray_FFMUE[1]);
		  h2FFMUEGen->Fill(fArray_FFMUE[2],fArray_FFMUE[3]);
		  break;
		}//if(jt == matchedUETrackIndex[iMIdx])
	    }//for(Int_t iMIdx = 0; iMIdx < nMatchedUETracks; iMIdx++)
	  //cout << "Rec matchedUERecTrack: " << MatchedUERecTrack  <<  endl;
	  if(!MatchedUERecTrack)
	    {
	      //cout << "I am inside MatchedUERecTrack->Fake" << endl;
	      //h2Fake
	      h2FFMUEFake->Fill(jetDetMC->Pt(),GetFFUEPC(kTRUE,jetDetMC,RecUETrkList,jt));
	    }
	}//for(Int_t jt = 0; jt < Ntracks_in_UE_Rec; jt++)

      for(Int_t jt = 0; jt < Ntracks_in_UE_Gen; jt++)
	{
	  AliVParticle* trackkineUEGen = dynamic_cast<AliVParticle*>(GenUETrkList->At(jt));
	  Int_t MatchedUEGenTrack = 0;
	  for(Int_t iMIdx = 0; iMIdx < nMatchedUETracks; iMIdx++)
	    {
	      //cout << "Gen UE jt: " << jt << endl;
	      if(jt == matchedUEMCIndex[iMIdx])
		{
		  //cout << "Gen UE jt: " << jt << " matchedUEMCIndex: " << matchedUEMCIndex[iMIdx] << endl;
		  MatchedUEGenTrack++;
		  break;
		}//if(jt == matchedUEMCIndex[iMIdx])
	    }//for(Int_t iMIdx = 0; iMIdx < nMatchedUETracks; iMIdx++)
	  //cout << "Gen matchedUEGenTrack: " << MatchedUEGenTrack  <<  endl;
	  if(!MatchedUEGenTrack)
	    {
	      //cout << "I am inside MatchedUEGenTrack->Miss" << endl;
	      //h2Miss
	      h2FFMUEMiss->Fill(jet->Pt(),GetFFUEPC(kFALSE,jet,GenUETrkList,jt));
	    }
	}//for(Int_t jt = 0; jt < nGenJetTracks; jt++)
      
      for(Int_t it = 0; it < nGenJetTracks; it++)
	{
	  AliVParticle* trackkineGen = (AliVParticle*)jet->TrackAt(it,tracksGen->GetArray());    
	  LabelGen[it] = trackkineGen->GetLabel();
	  //cout << "Gen Label: " << LabelGen[it] << endl;
	  //cout << "gen track pt: " << trackkineGen->Pt() << endl;
	}

      for(Int_t jt = 0; jt < nRecJetTracks; jt++)
	{
	  AliVParticle* trackkineRec = (AliVParticle*)jetDetMC->TrackAt(jt,tracksRec->GetArray());    
	  LabelRec[jt] = trackkineRec->GetLabel();
	  //cout << "Rec Label: " << LabelRec[jt] << endl;
	  //cout << "rec track pt: " << trackkineRec->Pt() << endl;
	}

      Int_t matchedTrackIndex[100] = {0}, matchedMCIndex[100] = {0};
      for(Int_t iarr = 0; iarr < 100; iarr++){matchedMCIndex[iarr] = -99; matchedTrackIndex[iarr] = -99;}
      Int_t nMatchedTracks = 0;
      for(Int_t iTrkIndex = 0; iTrkIndex < nRecJetTracks; iTrkIndex++)
	{
	  for(Int_t iMCIndex = 0; iMCIndex < nGenJetTracks; iMCIndex++)
	    {
	      //cout << " iTrkIndex: " << iTrkIndex << "  iMCIndex: " << iMCIndex << endl;
	      if(LabelGen[iMCIndex] == LabelRec[iTrkIndex])
		{
		  //cout << "Matched-> iTrkIndex: " << iTrkIndex << "  iMCIndex: " << iMCIndex << endl;
		  matchedTrackIndex[nMatchedTracks] = iTrkIndex;
		  matchedMCIndex[nMatchedTracks] = iMCIndex;
		  //cout << "Matched-> matchedTrackIndex: " << matchedTrackIndex[nMatchedTracks] << "  matchedMCIndex[nMatchedTracks]: " << matchedMCIndex[nMatchedTracks] << "  nMatchedTracks: " << nMatchedTracks << endl;
		  nMatchedTracks++;

		}
	    }
	}

      for(Int_t jt = 0; jt < nRecJetTracks; jt++)
	{
	  AliVParticle* trackkineRec = (AliVParticle*)jetDetMC->TrackAt(jt,tracksRec->GetArray());
	  Int_t MatchedRecTrack = 0;
	  for(Int_t iMIdx = 0; iMIdx < nMatchedTracks; iMIdx++)
	    {
	      //cout << "Rec jt: " << jt << endl;
	      if(jt == matchedTrackIndex[iMIdx])
		{
		  //cout << "Rec jt: " << jt << " matchedTrackIndex: " << matchedTrackIndex[iMIdx] << " matchedMCIndex: " << matchedMCIndex[iMIdx] << endl;
		  fArray_FFM[0] = jetDetMC->Pt();
		  fArray_FFM[2] = jet->Pt();
		  fArray_FFM[1] = GetFF(kTRUE,jetDetMC,matchedTrackIndex[iMIdx]); 
		  fArray_FFM[3] = GetFF(kFALSE,jet,matchedMCIndex[iMIdx]);
		  MatchedRecTrack++;
		  
		  //cout << fArray_FFM[0] << " " << fArray_FFM[1] << " " << fArray_FFM[2] << " " << fArray_FFM[3] << endl;
		  //h4
		  h4ResFFM->Fill(fArray_FFM);
		  h2FFMRec->Fill(fArray_FFM[0],fArray_FFM[1]);
		  h2FFMGen->Fill(fArray_FFM[2],fArray_FFM[3]);
		  break;
		}//if(jt == matchedTrackIndex[iMIdx])
	    }//for(Int_t iMIdx = 0; iMIdx < nMatchedtracks; iMIdx++)
	  //cout << "Rec matchedRecTrack: " << MatchedRecTrack  <<  endl;
	  if(!MatchedRecTrack)
	    {
	      //cout << "I am inside MatchedRecTrack->Fake" << endl;
	      //h2Fake
	      h2FFMFake->Fill(jetDetMC->Pt(),GetFF(kTRUE,jetDetMC,jt));
	    }
	}//for(Int_t jt = 0; jt < nRecJetTracks; jt++)

      for(Int_t jt = 0; jt < nGenJetTracks; jt++)
	{
	  AliVParticle* trackkineGen = (AliVParticle*)jet->TrackAt(jt,tracksGen->GetArray());
	  Int_t MatchedGenTrack = 0;
	  for(Int_t iMIdx = 0; iMIdx < nMatchedTracks; iMIdx++)
	    {
	      //cout << "Gen jt: " << jt << endl;
	      if(jt == matchedMCIndex[iMIdx])
		{
		  //cout << "Gen jt: " << jt << " matchedMCIndex: " << matchedMCIndex[iMIdx] << endl;
		  MatchedGenTrack++;
		  break;
		}//if(jt == matchedMCIndex[iMIdx])
	    }//for(Int_t iMIdx = 0; iMIdx < nMatchedtracks; iMIdx++)
	  //cout << "Gen matchedGenTrack: " << MatchedGenTrack  <<  endl;
	  if(!MatchedGenTrack)
	    {
	      //cout << "I am inside MatchedGenTrack->Miss" << endl;
	      //h2Miss
	      h2FFMMiss->Fill(jet->Pt(),GetFF(kFALSE,jet,jt));
	    }
	}//for(Int_t jt = 0; jt < nGenJetTracks; jt++)

      Double_t PtSumARec[20] = {0.};
      Double_t PtSumAGen[20] = {0.};

      Double_t PtSumARecUE[20] = {0.};
      Double_t PtSumAGenUE[20] = {0.};

      Double_t *PtSumRecFinal;
      Double_t *PtSumGenFinal;
      Double_t *PtSumUERecFinal;
      Double_t *PtSumUEGenFinal;
      PtSumRecFinal = GetPtSum(kTRUE,jetDetMC,PtSumARec);
      PtSumGenFinal = GetPtSum(kFALSE,jet,PtSumAGen);
      PtSumUERecFinal = GetPtSumUE(kTRUE,RecUETrkList,etaTiltedRec,phiTiltedRec,PtSumARecUE);
      PtSumUEGenFinal = GetPtSumUE(kFALSE,GenUETrkList,etaTiltedGen,phiTiltedGen,PtSumAGenUE);

      for(Int_t ibin = 1; ibin <= 20; ibin++)
	{
	  Double_t iR = 0.02*(ibin-1) + 0.01;
	  //cout << "  iR: "<< iR << "  ibin: "<< ibin << " Rec->pTSum:  "<<PtSumRecFinal[ibin-1] << " Gen->pTSum:  "<<PtSumGenFinal[ibin-1] << "   UE:  Rec->pTSum:  "<<PtSumUERecFinal[ibin-1] << " Gen->pTSum:  "<<PtSumUEGenFinal[ibin-1] << endl;
	  
	  fArray_PtSum[0] = jetDetMC->Pt();
	  fArray_PtSum[1] = jet->Pt();
	  fArray_PtSum[2] = iR;
	  fArray_PtSum[3] = PtSumRecFinal[ibin-1];
	  fArray_PtSum[4] = PtSumGenFinal[ibin-1];

	  h5ResPtSum->Fill(fArray_PtSum);

	  fArray_PtSumUE[0] = jetDetMC->Pt();
	  fArray_PtSumUE[1] = jet->Pt();
	  fArray_PtSumUE[2] = iR;
	  fArray_PtSumUE[3] = PtSumUERecFinal[ibin-1];
	  fArray_PtSumUE[4] = PtSumUEGenFinal[ibin-1];

	  h5ResPtSumUE->Fill(fArray_PtSumUE);
	}
      fArray_Nch[0] = jetDetMC->Pt();
      fArray_Nch[1] = GetNch(kTRUE);
      fArray_Nch[2] = jet->Pt();
      fArray_Nch[3] = GetNch(kFALSE);

      h4ResNch->Fill(fArray_Nch);
      h2NchRec->Fill(fArray_Nch[0],fArray_Nch[1]);
      h2NchGen->Fill(fArray_Nch[2],fArray_Nch[3]);

      fArray_NchUEPC[0] = jetDetMC->Pt();
      fArray_NchUEPC[1] = GetNchUEPC(kTRUE);
      fArray_NchUEPC[2] = jet->Pt();
      fArray_NchUEPC[3] = GetNchUEPC(kFALSE);
      
      h4ResNchUE->Fill(fArray_NchUEPC);
      h2NchRecUE->Fill(fArray_NchUEPC[0],fArray_NchUEPC[1]);
      h2NchGenUE->Fill(fArray_NchUEPC[2],fArray_NchUEPC[3]);

  
      Bool_t IsRecLess = kTRUE;
      Bool_t IsEqual = kFALSE;
      Int_t nMaxTracks = 0, nMinTracks = 0;

      if(nGenJetTracks > nRecJetTracks){nMaxTracks = nGenJetTracks; nMinTracks = nRecJetTracks; IsRecLess = kTRUE;}
      if(nGenJetTracks < nRecJetTracks){nMaxTracks = nRecJetTracks; nMinTracks = nGenJetTracks; IsRecLess = kFALSE;}
      if(nGenJetTracks == nRecJetTracks){nMaxTracks = nRecJetTracks; nMinTracks = nGenJetTracks; IsEqual = kTRUE;}
      //cout << "No. of tracks checking: " << nGenJetTracks << " " << nRecJetTracks << "  " << nMaxTracks << endl;
      for(Int_t j = 0; j < nMaxTracks; j++)
	{
	  fArray_FF[0] = jetDetMC->Pt();fArray_FF[2] = jet->Pt();
	  if(IsEqual){fArray_FF[1] = GetFF(kTRUE,jetDetMC,j); fArray_FF[3] =  GetFF(kFALSE,jet,j);}
	  else if(j >= nMinTracks && IsRecLess){fArray_FF[1] = 0.; fArray_FF[3] =  GetFF(kFALSE,jet,j);}
	  else if(j >= nMinTracks && !IsRecLess){fArray_FF[3] = 0.; fArray_FF[1] = GetFF(kTRUE,jetDetMC,j);}
	  else 
	    {fArray_FF[1] = GetFF(kTRUE,jetDetMC,j); fArray_FF[3] =  GetFF(kFALSE,jet,j);}
	  //for(Int_t icheck = 0; icheck < 4; icheck++)cout << fArray_FF[icheck] << endl;

	  h4ResFF->Fill(fArray_FF);
	  h2FFRec->Fill(fArray_FF[0],fArray_FF[1]);
	  h2FFGen->Fill(fArray_FF[2],fArray_FF[3]);
	}//for(Int_t j =0; j< nMaxTracks; j++)      
    }
  
}//FillMy2DResponseMatrices


//Get PtSum of the leading jet
Double_t* AliAnalysisTaskEmcalJetProperties::GetPtSum(Bool_t Isrec, AliEmcalJet *jet, Double_t* PtSumA)
{
  IsRec = (Bool_t)Isrec;
  AliJetContainer *jetCont = NULL;
  if(IsRec)
    jetCont = GetJetContainer(fContainerAktRec);
  else
    jetCont = GetJetContainer(fContainerAktGen);

  AliParticleContainer* tracks = jetCont->GetParticleContainer();
  
  //Float_t PtSumA[50] = {0.};
  for (Int_t j = 0; j < jet->GetNumberOfTracks(); j++)
    {
      AliVParticle* trackkine = (AliVParticle*)jet->TrackAt(j,tracks->GetArray());
      Float_t DelEta=-99.0;  Float_t DelPhi=-99.0; Float_t DelR=-99.0;
      DelEta = TMath::Abs(jet->Eta() - trackkine->Eta());
      DelPhi = TMath::Abs(jet->Phi() - trackkine->Phi());
      if(TMath::Abs(DelPhi)>TMath::Pi())DelPhi = TMath::Abs(DelPhi-TMath::TwoPi());
      DelR   = TMath::Sqrt(DelEta*DelEta + DelPhi*DelPhi);
      
      for(Int_t ibin = 1; ibin <= 20; ibin++)
	{
	  Float_t xlow = 0.02*(ibin-1);
	  Float_t xup  = 0.02*ibin;
	  if( xlow <= DelR && DelR < xup)
	    {
	      //cout << "b4->ibin: " << ibin << " PtSumA: " << PtSumA[ibin-1] << endl;
	      PtSumA[ibin-1]+= trackkine->Pt();
	      //cout << "track pT: " << trackkine->Pt() << endl;
	      //cout << "after->ibin: " << ibin << " PtSumA: " << PtSumA[ibin-1] << endl;
	    }//if loop
	}//for ibin loop
    }//track loop
  
  return PtSumA;
  
}

//Get PtSumUE of the leading jet
Double_t* AliAnalysisTaskEmcalJetProperties::GetPtSumUE(Bool_t Isrec, TList *UETrkList, Double_t etaTilted, Double_t phiTilted, Double_t* PtSumA)
{
  IsRec = (Bool_t)Isrec; 
  for (Int_t j = 0; j < UETrkList->GetEntries(); j++)
    {
      AliVParticle* trackkineUE = dynamic_cast<AliVParticle*>(UETrkList->At(j));
      Float_t DelEta=-99.0;  Float_t DelPhi=-99.0; Float_t DelR=-99.0;
      DelEta = TMath::Abs(etaTilted - trackkineUE->Eta());
      DelPhi = TMath::Abs(phiTilted - trackkineUE->Phi());
      if(TMath::Abs(DelPhi)>TMath::Pi())DelPhi = TMath::Abs(DelPhi-TMath::TwoPi());
      DelR   = TMath::Sqrt(DelEta*DelEta + DelPhi*DelPhi);
      //cout << "UE" << endl;
      for(Int_t ibin = 1; ibin <= 20; ibin++)
	{
	  Float_t xlow = 0.02*(ibin-1);
	  Float_t xup  = 0.02*ibin;
	  if( xlow <= DelR && DelR < xup)
	    {
	      PtSumA[ibin-1]+= trackkineUE->Pt();
	    }//if loop
	}//for ibin loop
    }//track loop
  
  return PtSumA;
  
}

//Get FF of the leading jet
Double_t AliAnalysisTaskEmcalJetProperties::GetFF(Bool_t Isrec, AliEmcalJet *jet, Int_t iTrack)
{
  IsRec = (Bool_t)Isrec;
  AliJetContainer *jetCont = NULL;
  if(IsRec)
    jetCont = GetJetContainer(fContainerAktRec);
  else
    jetCont = GetJetContainer(fContainerAktGen);
  
  Double_t FF = 0., TrackPt = 0.;    
  AliParticleContainer* tracks = jetCont->GetParticleContainer();
  //cout << "##################### IsRec: " << IsRec << endl;
  //cout << "Hib4" << endl;
  AliVParticle* trackkine = (AliVParticle*)jet->TrackAt(iTrack,tracks->GetArray());    
  TrackPt  = trackkine->Pt();
  if(jet->Pt())FF = TrackPt/jet->Pt();
  //cout << "###################GetFF-> FF: " << FF << endl;

  return FF;
}


//Get FF UE of the leading jet
Double_t AliAnalysisTaskEmcalJetProperties::GetFFUEPC(Bool_t Isrec, AliEmcalJet *jet, TList *UETrkList, Int_t iTrack)
{
  IsRec = (Bool_t)Isrec;
  Double_t FF = 0., TrackPt = 0.;    
  AliVParticle* trackkineUE = dynamic_cast<AliVParticle*>(UETrkList->At(iTrack));
  for(Int_t i = 0; i < UETrkList->GetEntries(); i++)
    {
      AliVParticle* trackkineUEtest = dynamic_cast<AliVParticle*>(UETrkList->At(i));
      FF = trackkineUEtest->Pt()/jet->Pt();
    }
  if(jet->Pt())
    {
      FF = trackkineUE->Pt()/jet->Pt();
    }

  return FF;
}

//Get the leading matched rec. jet index
Int_t AliAnalysisTaskEmcalJetProperties::GetMLeadingJet()
{
  AliEmcalJet* mJetR = NULL;
  Int_t mIndex = 0;
  Double_t LJetPt = 0., SJetPt = 0.;
  for(Int_t iJ = 0; iJ < matchedRecJetList->GetEntries(); iJ++)
    {           
      mJetR = (AliEmcalJet*)matchedRecJetList->At(iJ);      
      SJetPt = mJetR->Pt();
      if(SJetPt > LJetPt)
	{
	  mIndex = iJ;
	  LJetPt = SJetPt;
	}//if(SJetPt > LJetPt)
    }//for(Int_t iJ = 0; iJ < matchedRecJetList->GetEntries(); iJ++)

  return mIndex;
}

//Get the No. of tracks of the leading jet
Double_t AliAnalysisTaskEmcalJetProperties::GetNch(Bool_t Isrec)
{
  IsRec = (Bool_t)Isrec;
  TList *matchedJetList;
  if(IsRec)matchedJetList = (TList*)matchedRecJetList;
  else
    matchedJetList = (TList*)matchedGenJetList;
  Double_t N_Leading_Jet_Tracks = 0;
  Int_t nMatchedJetsCheck = 0;
  nMatchedJetsCheck = matchedJetList->GetEntries();
  //cout << "Matched jets: " << nMatchedJetsCheck << endl;
  AliEmcalJet* jet = NULL;
  Int_t mIndex = 0;//matched jet index in the TList
  if(nMatchedJetsCheck)
    {
      mIndex = GetMLeadingJet();
      jet = (AliEmcalJet*)matchedJetList->At(mIndex);      
      //cout << "*************************************************** jet pT: " << jet->Pt() << endl;
      N_Leading_Jet_Tracks = jet->GetNumberOfTracks();
      //cout << "***********************************************************GetNch->N_Leading_Jet_Tracks: " << N_Leading_Jet_Tracks << endl;
    }

  return N_Leading_Jet_Tracks;
}


//Get the No. of tracks of UE (PC)
Double_t AliAnalysisTaskEmcalJetProperties::GetNchUEPC(Bool_t Isrec)//TList *matchedJetList)
{
  IsRec = (Bool_t)Isrec;
  TList *matchedJetList;
  if(IsRec)matchedJetList = (TList*)matchedRecJetList;
  else
    matchedJetList = (TList*)matchedGenJetList;
  Int_t mIndex = 0;//matched jet index in the TList
  mIndex = GetMLeadingJet();
  Double_t N_Tracks_UE = 0;
  Int_t nMatchedJetsCheck = 0;
  nMatchedJetsCheck = matchedJetList->GetEntries();
  AliEmcalJet* jet = NULL;
  float fJetRadius = 0.4;
  const Double_t alpha = TMath::Pi()/2.;
  Double_t sumPtPerp = 0;
  Int_t Ntracks_in_UE = 0;
  
  if(nMatchedJetsCheck)
    {
      jet = (AliEmcalJet*)matchedJetList->At(mIndex);
      //cout << "*************************************************** jet pT: " << jet->Pt() << endl;

      fTrackListUE->Clear();
    
      GetTracksTiltedwrpJetAxis(alpha,fTrackListUE,jet,fJetRadius,sumPtPerp,Ntracks_in_UE);
      N_Tracks_UE = fTrackListUE->GetEntries();
      fTrackListUE->Clear();
    }
  //cout << "********************************************** UE Ntracks: " << Ntracks_in_UE << endl;
  return N_Tracks_UE;
}

