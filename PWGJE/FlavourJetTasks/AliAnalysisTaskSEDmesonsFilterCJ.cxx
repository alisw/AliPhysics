// $Id$
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
//
//
// Task for selecting D mesons to be used as an input for D-Jet correlations
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// A.Grelli (Utrecht University) a.grelli@uu.nl
// Xiaoming Zhang (LBNL)  XMZhang@lbl.gov
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include "TROOT.h"
#include <TH3F.h>

#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSEDmesonsFilterCJ.h"

ClassImp(AliAnalysisTaskSEDmesonsFilterCJ)

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ() :
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fCandidateType(0),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fOutput(0),
//fOutputCandidates(0),
  fCuts(0),
  fMinMass(0.),
  fMaxMass(0.),
  fCandidateArray(0)
//fIsSelectedArray(0)
{
//
// Default constructor
//

  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
  for (Int_t i=30; i--;) fSigmaD0[i] = 0.;
}

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ(const char *name, AliRDHFCuts *cuts, ECandidateType candtype) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fCandidateType(candtype),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
//fPDGdaughters(),
  fBranchName(""),
  fOutput(0),
//fOutputCandidates(0),
  fCuts(cuts),
  fMinMass(0.),
  fMaxMass(0.),
  fCandidateArray(0)
//fIsSelectedArray(0)
{
//
// Constructor. Initialization of Inputs and Outputs
//
 
  Info("AliAnalysisTaskSEDmesonsFilterCJ","Calling Constructor");

  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
  for (Int_t i=30; i--;) fSigmaD0[i] = 0.;

  const Int_t nptbins = fCuts->GetNPtBins();
  Float_t defaultSigmaD013[13] = { 0.012, 0.012, 0.012, 0.015, 0.015, 0.018, 0.018, 0.020, 0.020, 0.030, 0.030, 0.037, 0.040 };

  switch (fCandidateType) {
  case 0 :
    fCandidateName = "D0";
    fPDGmother = 421;
    fNProngs = 2;
    fPDGdaughters[0] = 211;  // pi 
    fPDGdaughters[1] = 321;  // K
    fPDGdaughters[2] = 0;    // empty
    fPDGdaughters[3] = 0;    // empty
    fBranchName = "D0toKpi";
    break;
  case 1 :
    fCandidateName = "Dstar";
    fPDGmother = 413;
    fNProngs = 3;
    fPDGdaughters[1] = 211; // pi soft
    fPDGdaughters[0] = 421; // D0
    fPDGdaughters[2] = 211; // pi fromD0
    fPDGdaughters[3] = 321; // K from D0
    fBranchName = "Dstar";

    if (nptbins<=13) {
      for (Int_t ipt=0; ipt<nptbins;ipt++) fSigmaD0[ipt] = defaultSigmaD013[ipt];
    } else {
      AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
    }
    break;
  default :
    printf("%d not accepted!!\n",fCandidateType);
    break;
  }

  if (fCandidateType==kD0toKpi) SetMassLimits(0.15, fPDGmother);
  if (fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);

  DefineOutput(1, TList::Class());       // histos
  DefineOutput(2, AliRDHFCuts::Class()); // my cuts
}

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::~AliAnalysisTaskSEDmesonsFilterCJ()
{
//
// destructor
//

  Info("~AliAnalysisTaskSEDmesonsFilterCJ","Calling Destructor");  
 
  if (fOutput) { delete fOutput; fOutput = 0; }
  if (fCuts)   { delete fCuts;   fCuts   = 0; }
  if (fCandidateArray)  { delete fCandidateArray;  fCandidateArray  = 0; }
//if (fIsSelectedArray) { delete fIsSelectedArray; fIsSelectedArray = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Init()
{
//
// Initialization
//

  if(fDebug>1) printf("AnalysisTaskSEDmesonsForJetCorrelations::Init() \n");

  switch (fCandidateType) {
    case 0: {
              AliRDHFCutsD0toKpi* copyfCutsDzero = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
              copyfCutsDzero->SetName("AnalysisCutsDzero");
              PostData(2, copyfCutsDzero);  // Post the data
            } break;
    case 1: {
              AliRDHFCutsDStartoKpipi* copyfCutsDstar = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
              copyfCutsDstar->SetName("AnalysisCutsDStar");
              PostData(2, copyfCutsDstar); // Post the cuts
            } break;
    default: return;
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects()
{ 
//
// output 
//

  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
//OpenFile(1);
  fOutput = new TList(); fOutput->SetOwner();
  DefineHistoForAnalysis(); // define histograms
 
//fOutputCandidates= new TList();
//fOutputCandidates->SetOwner();

  fCandidateArray = new TClonesArray("AliAODRecoDecayHF",0);
  fCandidateArray->SetName(Form("fCandidateArray%s",fCandidateName.Data()));

//fOutputCandidates->Add(fCandidateArray);
//fIsSelectedArray = new TClonesArray("Int_t",0);
//fIsSelectedArray->SetName(Form("fIsSelectedArray%s",suffix.Data()));
//fOutputCandidates->Add(fIsSelectedArray);

  PostData(1, fOutput);
//PostData(3, fOutputCandidates);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserExec(Option_t *)
{
//
// user exec
//

  // add cadidate branch
  fCandidateArray->Delete();
  if (!(InputEvent()->FindListObject(Form("fCandidateArray%s",fCandidateName.Data())))) InputEvent()->AddObject(fCandidateArray);

  // Load the event
  AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
 
  TClonesArray *arrayDStartoD0pi = 0;
  if (!aodEvent && AODEvent() && IsStandardAOD()) {

    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*>(AODEvent());

    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
    }
  } else {
    arrayDStartoD0pi = (TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
  }

  if (!arrayDStartoD0pi) {
    AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
    return;
  } else {
    AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   
  }
  
  TClonesArray* mcArray = 0x0;
  if (fUseMCInfo) {
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      return;
    }
  }

  //Histograms
  TH1I* hstat = (TH1I*)fOutput->FindObject("hstat");
  TH2F* hInvMassptD = (TH2F*)fOutput->FindObject("hInvMassptD");

  TH1F* hPtPion=0x0;
  if (fCandidateType==kDstartoKpipi) hPtPion = (TH1F*)fOutput->FindObject("hPtPion");
  hstat->Fill(0);

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;

  //Event selection
  Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
//TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
  if (!iseventselected) return;
  hstat->Fill(1);

  const Int_t nD = arrayDStartoD0pi->GetEntriesFast();
  hstat->Fill(2, nD);

//Int_t icountReco = 0;  // counters for efficiencies

  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
  Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi

  //D0 from D0 bar
  Int_t pdgdaughtersD0[2] = { 211, 321 };     // pi,K 
  Int_t pdgdaughtersD0bar[2] = { 321, 211 };  // K,pi 

  Int_t iCand =0;
  Int_t isSelected = 0;
  AliAODRecoDecayHF *charmCand = 0;

  Int_t mcLabel = -9999;
  Int_t pdgMeson = 413;
  if (fCandidateType==kD0toKpi) pdgMeson = 421;

  for (Int_t icharm=0; icharm<nD; icharm++) {   //loop over D candidates
    charmCand = (AliAODRecoDecayHF*)arrayDStartoD0pi->At(icharm); // D candidates
    if (!charmCand) return;


    if (fUseMCInfo) { // Look in MC, try to simulate the z
      if (fCandidateType==kDstartoKpipi) {
	AliAODRecoCascadeHF *temp = (AliAODRecoCascadeHF*)charmCand;
	mcLabel = temp->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      }

      if (fCandidateType==kD0toKpi) mcLabel = charmCand->MatchToMC(421,mcArray,fNProngs,fPDGdaughters);

      if (mcLabel<=0) continue;	       
    }

    // region of interest + cuts
    if (!fCuts->IsInFiducialAcceptance(charmCand->Pt(),charmCand->Y(pdgMeson))) continue;    

    isSelected = fCuts->IsSelected(charmCand,AliRDHFCuts::kAll,aodEvent); //selected
    if (!isSelected) continue;

    new ((*fCandidateArray)[iCand]) AliAODRecoDecayHF(*charmCand);
//  new ((*fIsSelectedArray)[iCand]) Int_t(isSelected);
    iCand++;

    Double_t masses[2];
    Double_t ptD = charmCand->Pt();
    if (fCandidateType==kDstartoKpipi) { //D*->D0pi->Kpipi

      //softpion from D* decay
      AliAODRecoCascadeHF *temp = (AliAODRecoCascadeHF*)charmCand;
      AliAODTrack *track2 = (AliAODTrack*)temp->GetBachelor();  
      hstat->Fill(3);

      // select D* in the D0 window.
      // In the cut object window is loose to allow side bands
      Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

      // retrieve the corresponding pt bin for the candidate
      // and set the expected D0 width (x3)
//    static const Int_t n = fCuts->GetNPtBins();
      Int_t bin = fCuts->PtBin(ptD);
      if(bin<0 || bin>=fCuts->GetNPtBins()) {
	AliError(Form("Pt %.3f out of bounds",ptD));
	continue;
      }

      AliInfo(Form("Pt bin %d and sigma D0 %.4f",bin,fSigmaD0[bin]));
      if ((temp->InvMassD0()>=(mPDGD0-3.*fSigmaD0[bin])) && (temp->InvMassD0()<=(mPDGD0+3.*fSigmaD0[bin]))) {	
	masses[0] = temp->DeltaInvMass(); //D*
	masses[1] = 0.; //dummy for D*

	//D*  delta mass
	hInvMassptD->Fill(masses[0], ptD); // 2 D slice for pt bins

	// D* pt and soft pion pt for good candidates  	      	
	Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
	Double_t invmassDelta = temp->DeltaInvMass();
	if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<0.0021) hPtPion->Fill(track2->Pt());
      }
    } //Dstar specific

    if (fCandidateType==kD0toKpi) { //D0->Kpi

      //needed quantities
      masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
      masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
      hstat->Fill(3);

      // mass vs pt
      if (isSelected==1 || isSelected==3) hInvMassptD->Fill(masses[0],ptD);
      if (isSelected>=2) hInvMassptD->Fill(masses[1],ptD);
    } //D0 specific

    charmCand = 0;  
  } // end of D cand loop

//AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));

//PostData(1,fOutput);
//PostData(3,fOutputCandidates);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Terminate(Option_t*)
{
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
  
  Info("Terminate"," terminate");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t range, Int_t pdg)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits
//

  Float_t mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();

  // compute the Delta mass for the D*
  if (fCandidateType==kDstartoKpipi) mass -= TDatabasePDG::Instance()->GetParticle(421)->Mass();
//cout << "mass ---------------" << endl;

  fMinMass = mass - range;
  fMaxMass = mass + range;

  AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
  if ((fMinMass<0.) || (fMaxMass<=0.) || (fMaxMass<fMinMass)) AliFatal("Wrong mass limits!\n");

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t lowlimit, Double_t uplimit)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits
//

  if (uplimit>lowlimit) {
    fMinMass = lowlimit;
    fMaxMass = uplimit;
  } else {
    printf("Error! Lower limit larger than upper limit!\n");
    fMinMass = uplimit - uplimit*0.2;
    fMaxMass = uplimit;
  }

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar(Int_t nptbins, Float_t *width)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar
//

  if (nptbins>30) {
    AliInfo("Maximum number of bins allowed is 30!");
    return kFALSE;
  }

  if (!width) return kFALSE;
  for (Int_t ipt=0; ipt<nptbins; ipt++) fSigmaD0[ipt]=width[ipt];

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis
//

  // Statistics 
  TH1I* hstat = new TH1I("hstat","Statistics",6,-0.5,5.5);
  hstat->GetXaxis()->SetBinLabel(1, "N ev anal");
  hstat->GetXaxis()->SetBinLabel(2, "N ev sel");
  hstat->GetXaxis()->SetBinLabel(3, "N cand");
  hstat->GetXaxis()->SetBinLabel(4, "N cand sel cuts");
  hstat->GetXaxis()->SetBinLabel(5, "N jets");
  hstat->GetXaxis()->SetBinLabel(6, "N cand in jet");
/*if(fUseMCInfo) {
    hstat->GetXaxis()->SetBinLabel(7,"N D");
    hstat->GetXaxis()->SetBinLabel(8,"N D in jet");
  }*/
  hstat->SetNdivisions(1);
  fOutput->Add(hstat);

  // Invariant mass related histograms
  const Int_t nbinsmass = 200;
  TH2F *hInvMass = new TH2F("hInvMassptD", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass, 100, 0., 50.);
  hInvMass->SetStats(kTRUE);
  hInvMass->GetXaxis()->SetTitle("mass (GeV/c)");
  hInvMass->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fOutput->Add(hInvMass);
  
  if (fCandidateType==kDstartoKpipi) {
    TH1F* hPtPion = new TH1F("hPtPion", "Primary pions candidates pt", 500, 0., 10.);
    hPtPion->SetStats(kTRUE);
    hPtPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtPion->GetYaxis()->SetTitle("entries");
    fOutput->Add(hPtPion);
  }

  return kTRUE; 
}
