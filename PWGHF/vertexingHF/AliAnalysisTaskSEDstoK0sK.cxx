/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

///-------------------------------------------------------------------------
///
///                   Base class for Ds->K0s+K Analysis
///
///
///    Production of invariant mass spectra and histograms for kinematic
///                     and topological cut studies.
///      Cuts have been centralized in the class AliRDHFCutsDstoK0sK.
///
///
///    Author: J.Hamon, julien.hamon@cern.ch (IPHC)
///
///-------------------------------------------------------------------------


#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDstoK0sK.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODPidHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliNormalizationCounter.h"
#include "AliNeutralTrackParam.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDstoK0sK.h"
#include "AliVertexerTracks.h"
#include "AliVTrack.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDstoK0sK);
/// \endcond





//__________________________________________________________________________
AliAnalysisTaskSEDstoK0sK::AliAnalysisTaskSEDstoK0sK() :
AliAnalysisTaskSE()
   , fAnalysisCuts(0)
   , fCounter(0)
   , fOutputSele(0)
   , fOutputCand(0)
   , fOutputPID(0)
   , fOutputMC(0)
   , fOutputNtuple(0)
   , fHisNEvents(0)
   , fHisRapidity(0)
   , fHisRapiditySel(0)
   , fHisPidTPCKaonVsPt(0)
   , fHisPidTOFKaonVsPt(0)
   , fHisPidTPCTOFKaon(0)
   , fHisPidTPCTOFKaonSel(0)
   , fHisPidTPCKaonVsPtSel(0)
   , fHisPidTOFKaonVsPtSel(0)
   , fNPtBins(0)
   , fMassRange(0.8)
   , fMassBinSize(0.002)
   , fReadMC(kFALSE)
   , fUseSelectionBit(kFALSE)
   , fFillNtuple(kFALSE)
   , fAODProtection(0)
{
   ///
   ///  Default Constructor
   ///

   for (Int_t ihis=0; ihis<3; ihis++) {
      fHisCentrality[ihis]       = 0;
      fHisCentralityVsMult[ihis] = 0;
   }

   for (Int_t ihis=0; ihis<5; ihis++) {
      fHisInvMassDs[ihis]        = 0;
      fHisInvMassDplus[ihis]     = 0;
      fHisInvMassK0s[ihis]       = 0;
      fHisPtK0s[ihis]            = 0;
      fHisPtBachelor[ihis]       = 0;
      fHisImpParK0s[ihis]        = 0;
      fHisImpParBach[ihis]       = 0;
      fHisCTauK0s[ihis]          = 0;
      fHisCosPointingDs[ihis]    = 0;
      fHisCosPointingXYDs[ihis]  = 0;
      fHisCosThetaStarK0s[ihis]  = 0;
      fHisCosThetaStarBach[ihis] = 0;
      fHisDCAK0sBach[ihis]       = 0;
      fHisDecayLxyDs[ihis]       = 0;
      fHisNormDecayLxyDs[ihis]   = 0;
   }

   for (Int_t ihis=0; ihis<kMaxPtBins+1; ihis++) {
      fPtLimits[ihis] = 0;
   }

   for (Int_t ihis=0; ihis<kNTupleVars; ihis++) {
      fCutsMinTupleVars[ihis] = 0;
      fCutsMaxTupleVars[ihis] = 0;
   }
}





//__________________________________________________________________________
AliAnalysisTaskSEDstoK0sK::AliAnalysisTaskSEDstoK0sK(const char *name,
                                                     AliRDHFCutsDstoK0sK *cuts,
                                                     Bool_t readMC,
                                                     Bool_t fillNtuple,
                                                     Int_t nCutsTuple,
                                                     Float_t *minCutsTuple,
                                                     Float_t *maxCutsTuple) :
AliAnalysisTaskSE(name)
   , fAnalysisCuts(cuts)
   , fCounter(0)
   , fOutputSele(0)
   , fOutputCand(0)
   , fOutputPID(0)
   , fOutputMC(0)
   , fOutputNtuple(0)
   , fHisNEvents(0)
   , fHisRapidity(0)
   , fHisRapiditySel(0)
   , fHisPidTPCKaonVsPt(0)
   , fHisPidTOFKaonVsPt(0)
   , fHisPidTPCTOFKaon(0)
   , fHisPidTPCTOFKaonSel(0)
   , fHisPidTPCKaonVsPtSel(0)
   , fHisPidTOFKaonVsPtSel(0)
   , fNPtBins(0)
   , fMassRange(0.8)
   , fMassBinSize(0.002)
   , fReadMC(readMC)
   , fUseSelectionBit(kFALSE)
   , fFillNtuple(fillNtuple)
   , fAODProtection(0)
{
   ///
   ///  Standard constructor
   ///

   for (Int_t ihis=0; ihis<3; ihis++) {
      fHisCentrality[ihis]       = 0;
      fHisCentralityVsMult[ihis] = 0;
   }

   for (Int_t ihis=0; ihis<5;ihis++) {
      fHisInvMassDs[ihis]        = 0;
      fHisInvMassDplus[ihis]     = 0;
      fHisInvMassK0s[ihis]       = 0;
      fHisPtK0s[ihis]            = 0;
      fHisPtBachelor[ihis]       = 0;
      fHisImpParK0s[ihis]        = 0;
      fHisImpParBach[ihis]       = 0;
      fHisCTauK0s[ihis]          = 0;
      fHisCosPointingDs[ihis]    = 0;
      fHisCosPointingXYDs[ihis]  = 0;
      fHisCosThetaStarK0s[ihis]  = 0;
      fHisCosThetaStarBach[ihis] = 0;
      fHisDCAK0sBach[ihis]       = 0;
      fHisDecayLxyDs[ihis]       = 0;
      fHisNormDecayLxyDs[ihis]   = 0;
   }

   for (Int_t ihis=0; ihis<kMaxPtBins+1; ihis++) {
      fPtLimits[ihis] = 0;
   }

   for (Int_t ihis=0; ihis<kNTupleVars; ihis++) {
      fCutsMinTupleVars[ihis] = 0;
      fCutsMaxTupleVars[ihis] = 0;
   }



   SetPtBins(fAnalysisCuts->GetNPtBins(), fAnalysisCuts->GetPtBinLimits());
   if (fFillNtuple)
      SetCutsTupleVariables(nCutsTuple, minCutsTuple, maxCutsTuple);


   // Output slot #1 writes into a TList container (fAnalysisCuts)
   DefineOutput(1, TList::Class());

   // Output slot #2 writes into a AliNormalizationCounter container (fCounter)
   DefineOutput(2, AliNormalizationCounter::Class());

   // Output slot #3 writes into a TList container (fOutputSele)
   DefineOutput(3, TList::Class());

   if (!fFillNtuple) {
      // Output slot #4 writes into a TList container (fOutputCand)
      DefineOutput(4, TList::Class());

      // Output slot #5 writes into a TList container (fOutputPID)
      DefineOutput(5, TList::Class());

      if (fReadMC) {
         // Output slot #6 writes into a TList container (fOutputMC)
         DefineOutput(6, TList::Class());
      }
   } else {
      // Output slot #4 writes into a TNtuple container (fOutputNtuple)
      DefineOutput(4, TNtuple::Class());
   }
}





//__________________________________________________________________________
AliAnalysisTaskSEDstoK0sK::~AliAnalysisTaskSEDstoK0sK()
{
   ///
   ///  Destructor
   ///

   if (fAnalysisCuts) { delete fAnalysisCuts; fAnalysisCuts =0; }
   if (fCounter)      { delete fCounter;      fCounter      =0; }
   if (fOutputSele)   { delete fOutputSele;   fOutputSele   =0; }
   if (fOutputCand)   { delete fOutputCand;   fOutputCand   =0; }
   if (fOutputPID)    { delete fOutputPID;    fOutputPID    =0; }
   if (fOutputMC)     { delete fOutputMC;     fOutputMC     =0; }
   if (fOutputNtuple) { delete fOutputNtuple; fOutputNtuple =0; }
}





//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::Init()
{
   ///
   ///  Initialisation
   ///
   AliDebug(1, "AliAnalysisTaskSEDstoK0sK::Init()");


   //---------------------------------------------------------------------------
   // - Output slot #1:  cut objects
   //---------------------------------------------------------------------------
   AliRDHFCutsDstoK0sK *listCuts = new AliRDHFCutsDstoK0sK(*fAnalysisCuts);
   listCuts->SetName(GetOutputSlot(1)->GetContainer()->GetName());
   PostData(1, listCuts);

   return;
}





//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::UserCreateOutputObjects()
{
   ///
   ///  Create the output containers
   ///
   AliDebug(1, "AliAnalysisTaskSEDstoK0sK::UserCreateOutputObjects()");



   // - Define the Pt range and binning
   Double_t ptBinsRange[fNPtBins+1];
   for (Int_t ipt=0; ipt<fNPtBins+1; ipt++) {
      ptBinsRange[ipt] = fPtLimits[ipt];
   }



   //---------------------------------------------------------------------------
   // - Output slot #2:  counter for Normalization
   //---------------------------------------------------------------------------
   fCounter = new AliNormalizationCounter(GetOutputSlot(2)->GetContainer()->GetName());
   fCounter->Init();



   //---------------------------------------------------------------------------
   // - Output slot #3:  various histograms of selected events
   //---------------------------------------------------------------------------
   fOutputSele = new TList();
   fOutputSele->SetOwner();
   fOutputSele->SetName(GetOutputSlot(3)->GetContainer()->GetName());


   fHisNEvents = new TH1F("fHisNEvents", ";;Entries", 21, 0, 21);
   fHisNEvents->GetXaxis()->SetBinLabel(1,  "nEvents Read");
   fHisNEvents->GetXaxis()->SetBinLabel(2,  "nEvents Matched dAOD");
   fHisNEvents->GetXaxis()->SetBinLabel(3,  "nEvents Mismatched dAOD");
   fHisNEvents->GetXaxis()->SetBinLabel(4,  "Analysed events");
   fHisNEvents->GetXaxis()->SetBinLabel(5,  "Rejection: Trigger");
   fHisNEvents->GetXaxis()->SetBinLabel(6,  "Rejection: Vertex reco");
   fHisNEvents->GetXaxis()->SetBinLabel(7,  "Rejection: Pileup");
   fHisNEvents->GetXaxis()->SetBinLabel(8,  "Rejection: Centrality");
   fHisNEvents->GetXaxis()->SetBinLabel(9,  "Rejection: VtxZ fiducial acc.");
   fHisNEvents->GetXaxis()->SetBinLabel(10, "Rejection: Physics Selection");
   fHisNEvents->GetXaxis()->SetBinLabel(11, "Selected events");
   fHisNEvents->GetXaxis()->SetBinLabel(12, "Read candidates");
   fHisNEvents->GetXaxis()->SetBinLabel(13, "FillRecoCasc fails");
   fHisNEvents->GetXaxis()->SetBinLabel(14, "Offline SV reco. fails");
   fHisNEvents->GetXaxis()->SetBinLabel(15, "Cand with sec. vertex");
   fHisNEvents->GetXaxis()->SetBinLabel(16, "In fiducial acceptance");
   fHisNEvents->GetXaxis()->SetBinLabel(17, "Selected Tracks");
   fHisNEvents->GetXaxis()->SetBinLabel(18, "Selected Candidates");
   fHisNEvents->GetXaxis()->SetBinLabel(19, "Selected PID");
   fHisNEvents->GetXaxis()->SetBinLabel(20, "Selected true D^{+}_{s}");
   fHisNEvents->GetXaxis()->SetBinLabel(21, "Selected true D^{+}");
   fOutputSele->Add(fHisNEvents);


   fHisCentrality[0] = new TH1F("fHisCentrality_all", "; Centrality; Entries", 1000, 0., 100.);
   fHisCentrality[0]->SetTitle("All read events");

   fHisCentrality[1] = (TH1F*) fHisCentrality[0]->Clone("fHisCentrality_selected");
   fHisCentrality[1]->SetTitle("Selected events");

   fHisCentrality[2] = (TH1F*) fHisCentrality[0]->Clone("fHisCentrality_rejected");
   fHisCentrality[2]->SetTitle("Events rejected due to centrality");

   fOutputSele->Add(fHisCentrality[0]);
   fOutputSele->Add(fHisCentrality[1]);
   fOutputSele->Add(fHisCentrality[2]);


   fHisCentralityVsMult[0] = new TH2F("fHisCentralityVsMult_all", "; Track multiplicities; Centrality; Entries", 100, 0., 10000., 40, 0., 100.);
   fHisCentralityVsMult[0]->SetTitle("All read events");

   fHisCentralityVsMult[1] = (TH2F*) fHisCentralityVsMult[0]->Clone("fHisCentralityVsMult_selected");
   fHisCentralityVsMult[1]->SetTitle("Selected events");

   fHisCentralityVsMult[2] = (TH2F*) fHisCentralityVsMult[0]->Clone("fHisCentralityVsMult_rejected");
   fHisCentralityVsMult[2]->SetTitle("Events rejected due to centrality");

   fOutputSele->Add(fHisCentralityVsMult[0]);
   fOutputSele->Add(fHisCentralityVsMult[1]);
   fOutputSele->Add(fHisCentralityVsMult[2]);


   fHisRapidity = new TH2F("fHisRapidity_pt", "; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{y}(D^{+}_{s}); Entries", 40, 0., 20., 80, -2., 2.);
   fHisRapidity->SetTitle("Production cuts + SV reconstruction");

   fHisRapiditySel = (TH2F*) fHisRapidity->Clone("fHisRapiditySel_pt");
   fHisRapiditySel->SetTitle("Selected candidates");

   fOutputSele->Add(fHisRapidity);
   fOutputSele->Add(fHisRapiditySel);


   fHisPidTPCKaonVsPt = new TH2F("fHisPidTPCKaonVsPt", "; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c});; Entries", fNPtBins, ptBinsRange, 200, -10., 10.);
   fHisPidTPCKaonVsPt->GetYaxis()->SetTitle("n#it{#sigma}_{TPC}^{K}");
   fHisPidTPCKaonVsPt->SetTitle("Selected tracks + candidates");

   fHisPidTPCKaonVsPtSel = (TH2F*) fHisPidTPCKaonVsPt->Clone("fHisPidTPCKaonVsPtSel");
   fHisPidTPCKaonVsPtSel->SetTitle("Selected tracks + candidates + PID");

   fOutputSele->Add(fHisPidTPCKaonVsPt);
   fOutputSele->Add(fHisPidTPCKaonVsPtSel);


   fHisPidTOFKaonVsPt = (TH2F*) fHisPidTPCKaonVsPt->Clone("fHisPidTOFKaonVsPt");
   fHisPidTOFKaonVsPt->GetYaxis()->SetTitle("n#it{#sigma}_{TOF}^{K}");
   fHisPidTOFKaonVsPt->SetTitle("Selected tracks + candidates");

   fHisPidTOFKaonVsPtSel = (TH2F*) fHisPidTOFKaonVsPt->Clone("fHisPidTOFKaonVsPtSel");
   fHisPidTOFKaonVsPtSel->SetTitle("Selected tracks + candidates + PID");

   fOutputSele->Add(fHisPidTOFKaonVsPt);
   fOutputSele->Add(fHisPidTOFKaonVsPtSel);


   fHisPidTPCTOFKaon = new TH2F("fHisPidTPCTOFKaon", "; n#it{#sigma}_{TPC}^{K}; n#it{#sigma}_{TOF}^{K}; Entries", 200, -10., 10., 200, -10., 10.);
   fHisPidTPCTOFKaon->SetTitle("Selected tracks + candidates");

   fHisPidTPCTOFKaonSel = (TH2F*) fHisPidTPCTOFKaon->Clone("fHisPidTPCTOFKaonSel");
   fHisPidTPCTOFKaonSel->SetTitle("Selected tracks + candidates + PID");

   fOutputSele->Add(fHisPidTPCTOFKaon);
   fOutputSele->Add(fHisPidTPCTOFKaonSel);



   if (!fFillNtuple) {

      //---------------------------------------------------------------------------
      // - Output slot #4 and #5:  histograms at candidate/PID level
      // - Output slot #6:  histograms for MC informations
      //---------------------------------------------------------------------------
      fOutputCand = new TList();
      fOutputCand->SetOwner();
      fOutputCand->SetName(GetOutputSlot(4)->GetContainer()->GetName());

      fOutputPID = new TList();
      fOutputPID->SetOwner();
      fOutputPID->SetName(GetOutputSlot(5)->GetContainer()->GetName());

      if (fReadMC) {
         fOutputMC = new TList();
         fOutputMC->SetOwner();
         fOutputMC->SetName(GetOutputSlot(6)->GetContainer()->GetName());
      }


      // - Define the invariant mass range and binning
      Float_t massDs    = TDatabasePDG::Instance()->GetParticle(431)->Mass();
      Float_t massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();

      Int_t nBinsMass = (Int_t) (fMassRange/fMassBinSize+0.5);
      if (nBinsMass%2==1) nBinsMass++;
      Float_t minMass      = massDs - 0.5*nBinsMass*fMassBinSize;
      Float_t maxMass      = massDs + 0.5*nBinsMass*fMassBinSize;
      Float_t minMassDplus = massDplus - 0.5*nBinsMass*fMassBinSize;
      Float_t maxMassDplus = massDplus + 0.5*nBinsMass*fMassBinSize;


      TString titleHisto(""), nameHisto("");
      TString suffixType[5] = { "kCan", "kPid", "mcSig", "mcBackg", "mcDplusRefl" };
      TString titleType[5]  = { "kCandidate", "kCandidate + kPID", "MC Signal", "MC Background", "MC D^{+} reflection" };



      for (Int_t itype=0; itype<5; itype++) {

         // - Do not create MC histograms if fReadMC is not true
         if (itype>=2 && !fReadMC) break;


         nameHisto  = Form("fHisInvMassDs_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); M_{inv}(K^{0}_{S}K^{#pm}); Entries", titleType[itype].Data());
         fHisInvMassDs[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, nBinsMass, minMass, maxMass);

         nameHisto  = Form("fHisInvMassDplus_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); M_{inv}(K^{0}_{S}#pi^{#pm}); Entries", titleType[itype].Data());
         fHisInvMassDplus[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, nBinsMass, minMassDplus, maxMassDplus);

         nameHisto  = Form("fHisInvMassK0s_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); M_{inv}(#pi^{+},#pi^{-}); Entries", titleType[itype].Data());
         fHisInvMassK0s[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 200, 0.4, 0.6);

         nameHisto  = Form("fHisPtK0s_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{p}_{T}(K^{0}_{S}) (GeV/#it{c}); Entries", titleType[itype].Data());
         fHisPtK0s[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 160, 0., 40.);

         nameHisto  = Form("fHisPtBachelor_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{p}_{T}(K^{+}) (GeV/#it{c}); Entries", titleType[itype].Data());
         fHisPtBachelor[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 160, 0., 40.);

         nameHisto  = Form("fHisImpParK0s_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{d_{0}}(K^{0}_{S}) (cm); Entries", titleType[itype].Data());
         fHisImpParK0s[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 1000, -5., 5.);

         nameHisto  = Form("fHisImpParBach_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{d_{0}}(K^{+}) (cm); Entries", titleType[itype].Data());
         fHisImpParBach[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 1000, -5., 5.);

         nameHisto  = Form("fHisCTauK0s_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{c#times#tau}(K^{0}_{S}); Entries", titleType[itype].Data());
         fHisCTauK0s[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 160, 0., 40.);

         nameHisto  = Form("fHisCosPointingDs_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); cos(#it{#theta_{p}})(D_{s}^{+}); Entries", titleType[itype].Data());
         fHisCosPointingDs[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 400, 0.0, 1.);

         nameHisto  = Form("fHisCosPointingXYDs_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); cos(#it{#theta_{p}^{xy}})(D_{s}^{+}); Entries", titleType[itype].Data());
         fHisCosPointingXYDs[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 800, -1.0, 1.);

         nameHisto  = Form("fHisCosThetaStarK0s_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); cos(#it{#theta^{*}})(K_{S}^{0}); Entries", titleType[itype].Data());
         fHisCosThetaStarK0s[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 800, -1., 1.);

         nameHisto  = Form("fHisCosThetaStarBach_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); cos(#it{#theta^{*}})(K^{+}); Entries", titleType[itype].Data());
         fHisCosThetaStarBach[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 800, -1., 1.);

         nameHisto  = Form("fHisDCAK0sBach_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); DCA(K^{0}_{S},K^{+}) (cm); Entries", titleType[itype].Data());
         fHisDCAK0sBach[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 500, 0., 10.);

         nameHisto  = Form("fHisDecayLxyDs_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{L_{xy}}(D^{+}_{s}) (cm); Entries", titleType[itype].Data());
         fHisDecayLxyDs[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 400, 0., 1.0);

         nameHisto  = Form("fHisNormDecayLxyDs_%s", suffixType[itype].Data());
         titleHisto = Form("%s; #it{p}_{T}(D^{+}_{s}) (GeV/#it{c}); #it{L_{xy}}/#it{#sigma_{xy}}(D^{+}_{s}) (cm); Entries", titleType[itype].Data());
         fHisNormDecayLxyDs[itype] = new TH2F(nameHisto.Data(), titleHisto.Data(), fNPtBins, ptBinsRange, 200, 0., 40.);


         TList *theOutput = 0;
         if (itype==0)      theOutput = fOutputCand;
         else if (itype==1) theOutput = fOutputPID;
         else               theOutput = fOutputMC;

         theOutput->Add(fHisInvMassDs[itype]);
         theOutput->Add(fHisInvMassDplus[itype]);
         theOutput->Add(fHisInvMassK0s[itype]);
         theOutput->Add(fHisPtK0s[itype]);
         theOutput->Add(fHisPtBachelor[itype]);
         theOutput->Add(fHisImpParK0s[itype]);
         theOutput->Add(fHisImpParBach[itype]);
         theOutput->Add(fHisCTauK0s[itype]);
         theOutput->Add(fHisCosPointingDs[itype]);
         theOutput->Add(fHisCosPointingXYDs[itype]);
         theOutput->Add(fHisCosThetaStarK0s[itype]);
         theOutput->Add(fHisCosThetaStarBach[itype]);
         theOutput->Add(fHisDCAK0sBach[itype]);
         theOutput->Add(fHisDecayLxyDs[itype]);
         theOutput->Add(fHisNormDecayLxyDs[itype]);

      } // end for



   } else {

      //---------------------------------------------------------------------------
      // - Output slot #4:  TNtuple for candidates on data
      //---------------------------------------------------------------------------
      TString listVar("K0sInvMass");           //   0  V0:    invariant mass (pi+pi-)
      listVar.Append(":K0sPt");                //   1  V0:    Pt
      listVar.Append(":K0sRxy");               //   2  V0:    decay length XY
      listVar.Append(":K0sCTau");              //   3  V0:    pseudo proper decay length (L*m/p)
      listVar.Append(":K0sCosPA");             //   4  V0:    cosine pointing angle
      listVar.Append(":K0sd0");                //   5  V0:    impact parameter
      listVar.Append(":K0sd0DauPos");          //   6  V0:    impact parameter positive daughter
      listVar.Append(":K0sd0DauNeg");          //   7  V0:    impact parameter negative daughter
      listVar.Append(":K0sDCAdau");            //   8  V0:    DCA prong-to-prong (pi+,pi-)
      listVar.Append(":K0sPtDauPos");          //   9  V0:    Pt positive daughters
      listVar.Append(":K0sPtDauNeg");          //  10  V0:    Pt positive daughters
      listVar.Append(":BachPt");               //  11  Bach:  Pt
      listVar.Append(":Bachd0");               //  12  Bach:  impact parameter
      listVar.Append(":BachNsigmaTPCpion");    //  13  Bach:  N sigmaTPC pion
      listVar.Append(":BachNsigmaTPCkaon");    //  14  Bach:  N sigmaTPC kaon
      listVar.Append(":BachNsigmaTPCproton");  //  15  Bach:  N sigmaTPC proton
      listVar.Append(":BachNsigmaTOFpion");    //  16  Bach:  N sigmaTOF pion
      listVar.Append(":BachNsigmaTOFkaon");    //  17  Bach:  N sigmaTOF kaon
      listVar.Append(":BachNsigmaTOFproton");  //  18  Bach:  N sigmaTOF proton
      listVar.Append(":CanInvMassDs");         //  19  Ds:    invariant mass (K0sK+)
      listVar.Append(":CanInvMassDplus");      //  20  Ds:    invariant mass (K0spi+)
      listVar.Append(":CanPt");                //  21  Ds:    Pt
      listVar.Append(":CanDCAProngToProng");   //  22  Ds:    DCA prong-to-prong (K0s,K)
      listVar.Append(":CanCosThetaStarK0s");   //  23  Ds:    cosine theta* K0s
      listVar.Append(":CanCosThetaStarBach");  //  24  Ds:    cosine theta* K+
      listVar.Append(":CanCosPA");             //  25  Ds:    cosine pointing angle
      listVar.Append(":CanCosPAxy");           //  26  Ds:    cosine pointing angle XY
      listVar.Append(":CanDLengthXY");         //  27  Ds:    decay length XY
      listVar.Append(":CanNormDLengthXY");     //  28  Ds:    normalised decay length XY
      listVar.Append(":CanDLength3D");         //  29  Ds:    decay length 3D
      listVar.Append(":CanNormDLength3D");     //  30  Ds:    normalised decay length 3D
      listVar.Append(":CanSigmaVtx");          //  31  Ds:    sigma vertex
      listVar.Append(":CanNormTopoBach");      //  32  Ds:    difference between measured and expected normalised d0(K+)
      listVar.Append(":CanNormTopoK0s");       //  33  Ds:    difference between measured and expected normalised d0(K0s)
      if (fReadMC) {
         listVar.Append(":mcTruth");           //  34  MC:    MC informations: 10 (signal), 20 (D+ reflection),
                                               //                              -1 (backg w/ true K0S), -2 (backg w/ false K0S)
      }
      // listVar.Append(":CanCosThetaRFrame");    //  32  Ds:    cosine of angle between (K0s-bach) in the Ds rest frame


      fOutputNtuple = new TNtuple(GetOutputSlot(4)->GetContainer()->GetName(), "TNtuple for reconstructed candidates on real data", listVar.Data());

   }



   PostData(2, fCounter);
   PostData(3, fOutputSele);
   if (!fFillNtuple) {
      PostData(4, fOutputCand);
      PostData(5, fOutputPID);
      if (fReadMC)
         PostData(6, fOutputMC);
   } else {
      PostData(4, fOutputNtuple);
   }

   return;
}





//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::UserExec(Option_t* /*option*/)
{
   ///
   ///  Execute analysis for the current event: HF cascade candidates
   ///
   AliDebug(1, "AliAnalysisTaskSEDstoK0sK::UserExec()");


   //---------------------------------------------------------------------------
   // - Load the event and the CascadesHF branch
   //---------------------------------------------------------------------------

   AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
   fHisNEvents->Fill(0); // all events


   if(fAODProtection>=0){
      //   Protection against different number of events in the AOD and deltaAOD
      //   In case of discrepancy the event is rejected.
      Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
      if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
         // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
         fHisNEvents->Fill(2);
         return;
      }
      fHisNEvents->Fill(1);
   }


   TClonesArray *arrayCascadesHF = 0;
   if(!aod && AODEvent() && IsStandardAOD()) {
      // In case there is an AOD handler writing a standard AOD, use the AOD
      // event in memory rather than the input (ESD) event.
      aod = dynamic_cast<AliAODEvent*> (AODEvent());
      // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler *aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
      if (aodHandler->GetExtensions()) {
         AliAODExtension *ext = (AliAODExtension*) aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
         AliAODEvent *aodFromExt = ext->GetAOD();
         arrayCascadesHF = (TClonesArray*) aodFromExt->GetList()->FindObject("CascadesHF");
      }
   } else if (aod) {
      arrayCascadesHF = (TClonesArray*) aod->GetList()->FindObject("CascadesHF");
   }

   if (!aod || !arrayCascadesHF) {
      AliWarning("AliAnalysisTaskSEDstoK0sK::UserExec: CascadesHF branch not found!");
      return;
   }


   // Load MC particles
   TClonesArray   *mcArray  = 0;
   AliAODMCHeader *mcHeader = 0;

   if (fReadMC) {
      mcArray = (TClonesArray*) aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!mcArray) {
         AliWarning("AliAnalysisTaskSEDstoK0sK::UserExec: MC particles branch not found!");
         return;
      }
      mcHeader = (AliAODMCHeader*) aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if (!mcHeader) {
         AliWarning("AliAnalysisTaskSEDstoK0sK::UserExec: MC header branch not found!");
         return;
      }

      // FillMCAcceptanceHistos has to take place here
   }



   //---------------------------------------------------------------------------
   // - Event selection
   //---------------------------------------------------------------------------

   // Fix for temporary bug in ESDfilter
   // The AODs with null vertex pointer didn't pass the PhysSel
   if (!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

   fHisNEvents->Fill(3);
   fCounter->StoreEvent(aod, fAnalysisCuts, fReadMC);


   Float_t nTracks = aod->GetNumberOfTracks();
   Float_t evCentr = fAnalysisCuts->GetUseCentrality()>0 ? fAnalysisCuts->GetCentrality(aod) : 0.;

   fHisCentrality[0]->Fill(evCentr);
   fHisCentralityVsMult[0]->Fill(nTracks, evCentr);


   if ( !(fAnalysisCuts->IsEventSelected(aod)) ) {

      if (fAnalysisCuts->GetWhyRejection() == 5)  fHisNEvents->Fill(4);  // Rejection: Trigger
      if (fAnalysisCuts->GetWhyRejection() == 0)  fHisNEvents->Fill(5);  // Rejection: Vertex reco
      if (fAnalysisCuts->GetWhyRejection() == 1)  fHisNEvents->Fill(6);  // Rejection: Pileup
      if (fAnalysisCuts->GetWhyRejection() == 2)  fHisNEvents->Fill(7);  // Rejection: Centrality
      if (fAnalysisCuts->GetWhyRejection() == 6)  fHisNEvents->Fill(8);  // Rejection: VtxZ fiducial acc.
      if (fAnalysisCuts->GetWhyRejection() == 7)  fHisNEvents->Fill(9);  // Rejection: Physics selection

      if (fAnalysisCuts->GetWhyRejection() == 2) {
         fHisCentrality[2]->Fill(evCentr);
         fHisCentralityVsMult[2]->Fill(nTracks, evCentr);
      }

      return;
   }


   fHisNEvents->Fill(10);
   fHisCentrality[1]->Fill(evCentr);
   fHisCentralityVsMult[1]->Fill(nTracks, evCentr);


   Int_t nCascades = arrayCascadesHF->GetEntriesFast();
   AliDebug(1, Form("AliAnalysisTaskSEDstoK0sK::UserExec:  %d Cascades found for this event!", nCascades));




   //---------------------------------------------------------------------------
   // - Analysis: Loop over the candidates
   //---------------------------------------------------------------------------

   // vHF object is needed to call the method that refills the missing info of the candidates
   // if they have been deleted in dAOD reconstruction phase
   // in order to reduce the size of the file
   AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();



   for (Int_t iCasc=0; iCasc<nCascades; iCasc++)
   {

      //-----------------------------------------
      // - Get candidates and fill missing info
      //-----------------------------------------
      AliAODRecoCascadeHF *dCan = (AliAODRecoCascadeHF*) arrayCascadesHF->UncheckedAt(iCasc);

      if (fUseSelectionBit) {
         if ( !(dCan->CheckCascadeFlags(AliRDHFCuts::kDstoK0sCuts)) ) continue;
      }

      fHisNEvents->Fill(11);


      // Fill the data members of the candidate only if they are empty.
      if (!(vHF->FillRecoCasc(aod, dCan, kFALSE, kTRUE))) {
         fHisNEvents->Fill(12);
         continue;
      }


      // Reconstruct secondary vertex if it does not already exist (standard dAOD)
      if (dCan->GetIsFilled()==1) {
         if (!(vHF->RecoSecondaryVertexForCascades(aod, dCan))) {
            fHisNEvents->Fill(13);
            continue;
         }
      }

      fHisNEvents->Fill(14);



      //-----------------------------------------
      // - Select candidates step by step
      //-----------------------------------------
      Double_t ptCand   = dCan->Pt();
      Double_t rapidity = dCan->Y(431);

      fHisRapidity->Fill(ptCand, rapidity);


      if (!(fAnalysisCuts->IsInFiducialAcceptance(ptCand, rapidity))) continue;
      fHisNEvents->Fill(15);


      if (!(fAnalysisCuts->IsSelected(dCan, AliRDHFCuts::kTracks, aod))) continue;
      fHisNEvents->Fill(16);

      // FillTheTree(dCan, mcArray);


      if (!(fAnalysisCuts->IsSelected(dCan, AliRDHFCuts::kCandidate, aod))) continue;
      fHisNEvents->Fill(17);

      FillHistogramsVar(dCan, AliRDHFCuts::kCandidate);
      FillHistogramsPID(dCan, AliRDHFCuts::kCandidate);


      if (!(fAnalysisCuts->IsSelected(dCan, AliRDHFCuts::kPID, aod))) continue;
      fHisNEvents->Fill(18);

      fHisRapiditySel->Fill(ptCand, rapidity);

      FillHistogramsVar(dCan, AliRDHFCuts::kPID, mcArray);
      FillHistogramsPID(dCan, AliRDHFCuts::kPID, mcArray);


      FillTheTree(dCan, aod, mcArray);


   } // end loop cascades



   delete vHF;


   PostData(2, fCounter);
   PostData(3, fOutputSele);
   if (!fFillNtuple) {
      PostData(4, fOutputCand);
      PostData(5, fOutputPID);
      if (fReadMC)
         PostData(6, fOutputMC);
   } else {
      PostData(4, fOutputNtuple);
   }

   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::Terminate(Option_t* /*option*/)
{
   ///
   ///  Terminate analysis
   ///
   AliDebug(1, "AliAnalysisTaskSEDstoK0sK::Terminate()");


   // Output slot #3
   fOutputSele = dynamic_cast<TList*> (GetOutputData(3));
   if (!fOutputSele) {
      AliWarning("ERROR: fOutputSele not available");
      return;
   }


   if (!fFillNtuple) {

      // Output slot #4
      fOutputCand = dynamic_cast<TList*> (GetOutputData(4));
      if (!fOutputCand) {
         AliWarning("ERROR: fOutputCand not available");
         return;
      }

      // Output slot #5
      fOutputPID = dynamic_cast<TList*> (GetOutputData(5));
      if (!fOutputPID) {
         AliWarning("ERROR: fOutputPID not available");
         return;
      }

      if (fReadMC) {
         // Output slot #6
         fOutputMC = dynamic_cast<TList*> (GetOutputData(6));
         if (!fOutputMC) {
            AliWarning("ERROR: fOutputMC not available");
            return;
         }
      }

   } else {

      // Output slot #4
      fOutputNtuple = dynamic_cast<TNtuple*> (GetOutputData(4));
      if (!fOutputNtuple) {
         AliWarning("ERROR: fOutputNtuple not available");
         return;
      }

   }

   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::SetPtBins(Int_t nBins, Float_t* limitsPt)
{
   ///
   ///  Set the number of Pt bins and the Pt bin sizes for the analysis
   ///


   if (nBins<=kMaxPtBins) {

      fNPtBins = nBins;
      for (Int_t ipt=0; ipt<fNPtBins+1; ipt++)             fPtLimits[ipt] = limitsPt[ipt];
      for (Int_t ipt=fNPtBins+1; ipt<kMaxPtBins+1; ipt++)  fPtLimits[ipt] = 99999999.;

   } else {

      // More Pt bins than allowed by the class header
      AliWarning(Form("Maximum number of Pt bins reached: %d", kMaxPtBins));
      fNPtBins = kMaxPtBins;
      fPtLimits[0] = 2.;
      fPtLimits[1] = 4.;
      fPtLimits[2] = 6.;
      fPtLimits[3] = 8.;
      for (Int_t ipt=4; ipt<kMaxPtBins+1; ipt++)  fPtLimits[ipt] = 99999999.;

   }

   if (fDebug>1) {
      AliInfo(Form("Number of Pt bins for the analysis: %d", fNPtBins));
      for (Int_t ipt=0; ipt<fNPtBins; ipt++) {
         AliInfo(Form("  Pt bin %02d:  [%8.1f, %8.1f]", ipt, fPtLimits[ipt], fPtLimits[ipt+1]));
      }
   }

   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::SetCutsTupleVariables(Int_t nCuts, Float_t *minCuts, Float_t *maxCuts)
{
   ///
   ///  Set the minimum and maximum cut values for variables to be stored in the tuple
   ///


   if (nCuts==kNTupleVars) {

      for (Int_t ivar=0; ivar<kNTupleVars; ivar++) {
         fCutsMinTupleVars[ivar] = minCuts[ivar];
         fCutsMaxTupleVars[ivar] = maxCuts[ivar];
      }

   } else {

      // Different number of cuts than the number of variables in the tuple
      AliWarning(Form("%d cuts are set for %d variables in the tuple: default initialisation", nCuts, kNTupleVars));

      for (Int_t ivar=0; ivar<kNTupleVars; ivar++) {
         fCutsMinTupleVars[ivar] = -999.;
         fCutsMaxTupleVars[ivar] =  999.;
      }

   }


   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::FillHistogramsVar(AliAODRecoCascadeHF* dCan, AliRDHFCuts::ESelLevel selFlag, TClonesArray* mcArray)
{
   ///
   ///  Fill variable histograms depending on the selection level.
   ///  If the Ntuple has to be filled (fFillNtuple), these histograms are not filled but
   ///  only TPC/TOF control histograms (FillHistogramsPID).
   ///
   ///  In case of reading Monte Carlo (fReadMC), check if the cascade candidate is
   ///      - signal: true Ds+ -> K0S + K+
   ///      - background: wrong Ds+ -> K0S + K+
   ///      - reflection: true D+ -> K0S + pi+
   ///
   AliDebug(1, "AliAnalysisTaskSEDstoK0sK::FillHistogramsVar()");


   if (fFillNtuple) return;


   //---------------------------------------------------------------------------
   // - Monte Carlo: check the candidate at MC level
   //---------------------------------------------------------------------------
   Bool_t  okMCsignal = kFALSE;
   Bool_t  okMCbackg  = kFALSE;
   Bool_t  okMCreflec = kFALSE;

   if (selFlag==AliRDHFCuts::kPID && fReadMC) {
      if (!mcArray) return;

      if (MatchToMCDstoK0sKSignal(dCan, mcArray)>=0) {
         AliDebug(2, "This candidate matches a MC signal of Ds+ -> K0S + K+");
         okMCsignal = kTRUE;
      } else if (MatchToMCDplustoK0spiSignal(dCan, mcArray)>=0) {
         AliDebug(2, "This candidate matches a MC signal of D+ -> K0S + pi+");
         okMCreflec = kTRUE;
      } else {
         AliDebug(2, "This candidate is a MC background");
         okMCbackg = kTRUE;
      }
   }



   //---------------------------------------------------------------------------
   // - Fill the histograms depending on the selection level
   //    - if 'kCandidate' level, fill only index=0 (histograms 'kCand')
   //    - if 'kPID' level and not reading MC, fill only index=1 (histograms 'kPID')
   //    - if 'kPID' level and reading MC, fill index=1 ('kPID'), 2 ('mcSig'), 3 ('mcBackg'), 4 ('mcDplusRefl')
   //---------------------------------------------------------------------------
   Float_t  ptCand = dCan->Pt();


   for (Int_t index=0; index<5; index++) {

      if (index==0) {
         if (selFlag != AliRDHFCuts::kCandidate) continue;
      } else if (index==1) {
         if (selFlag != AliRDHFCuts::kPID) continue;
      } else if (index==2) {
         if (!okMCsignal) continue;
      } else if (index==3) {
         if (!okMCbackg) continue;
      } else if (index==4) {
         if (!okMCreflec) continue;
      }

      fHisInvMassDs[index]        ->Fill( ptCand,  dCan->InvMassDstoK0sK()              );
      fHisInvMassDplus[index]     ->Fill( ptCand,  dCan->InvMassDplustoK0spi()          );
      fHisInvMassK0s[index]       ->Fill( ptCand,  dynamic_cast<AliAODv0*>(dCan->Getv0())->MassK0Short());
      fHisPtK0s[index]            ->Fill( ptCand,  dCan->PtProng(1)                     );
      fHisPtBachelor[index]       ->Fill( ptCand,  dCan->PtProng(0)                     );
      fHisImpParK0s[index]        ->Fill( ptCand,  dCan->Getd0Prong(1)                  );
      fHisImpParBach[index]       ->Fill( ptCand,  dCan->Getd0Prong(0)                  );
      fHisCTauK0s[index]          ->Fill( ptCand,  dynamic_cast<AliAODv0*>(dCan->Getv0())->Ct(310, dCan->GetOwnPrimaryVtx()));
      fHisCosPointingDs[index]    ->Fill( ptCand,  dCan->CosPointingAngle()             );
      fHisCosPointingXYDs[index]  ->Fill( ptCand,  dCan->CosPointingAngleXY()           );
      fHisCosThetaStarK0s[index]  ->Fill( ptCand,  dCan->CosThetaStar(1, 431, 321, 310) );
      fHisCosThetaStarBach[index] ->Fill( ptCand,  dCan->CosThetaStar(0, 431, 321, 310) );
      fHisDCAK0sBach[index]       ->Fill( ptCand,  dCan->GetDCA()                       );
      fHisDecayLxyDs[index]       ->Fill( ptCand,  dCan->DecayLengthXY()                );
      fHisNormDecayLxyDs[index]   ->Fill( ptCand,  dCan->NormalizedDecayLengthXY()      );

   }


   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::FillHistogramsPID(AliAODRecoCascadeHF* dCan, AliRDHFCuts::ESelLevel selFlag, TClonesArray* mcArray)
{
   ///
   ///  Fill PID histograms (TOF and TPC) depending on the selection level.
   ///
   AliDebug(2, "AliAnalysisTaskSEDstoK0sK::FillHistogramsPID()");


   Double_t nSigmaTPC = -999;
   Double_t nSigmaTOF = -999.;
   Float_t  ptCand    = dCan->Pt();

   if (fAnalysisCuts->GetIsUsePID()) {
      AliAODPidHF *pidhf = (AliAODPidHF*) fAnalysisCuts->GetPidHF();
      if (pidhf->CheckTPCPIDStatus(dCan->GetBachelor())) pidhf->GetnSigmaTPC(dCan->GetBachelor(), AliPID::kKaon, nSigmaTPC);
      if (pidhf->CheckTOFPIDStatus(dCan->GetBachelor())) pidhf->GetnSigmaTOF(dCan->GetBachelor(), AliPID::kKaon, nSigmaTOF);
   }


   if (selFlag==AliRDHFCuts::kCandidate) {

      //---------------------------------------------------------------------------
      // - Topological selection levels
      //---------------------------------------------------------------------------
      fHisPidTPCKaonVsPt->Fill(ptCand, nSigmaTPC);
      fHisPidTOFKaonVsPt->Fill(ptCand, nSigmaTOF);
      fHisPidTPCTOFKaon->Fill(nSigmaTPC, nSigmaTOF);

   } else if (selFlag==AliRDHFCuts::kPID || selFlag==AliRDHFCuts::kAll) {

      //---------------------------------------------------------------------------
      // - Topological+PID selection levels
      //---------------------------------------------------------------------------
      fHisPidTPCKaonVsPtSel->Fill(ptCand, nSigmaTPC);
      fHisPidTOFKaonVsPtSel->Fill(ptCand, nSigmaTOF);
      fHisPidTPCTOFKaonSel->Fill(nSigmaTPC, nSigmaTOF);

      //-----------------------------------------
      // - Check the candidate at MC level
      //-----------------------------------------
      if (fReadMC && mcArray) {
         if (MatchToMCDstoK0sKSignal(dCan, mcArray)>=0) fHisNEvents->Fill(19);
         else if (MatchToMCDplustoK0spiSignal(dCan, mcArray)>=0) fHisNEvents->Fill(20);
      }

   }


   return;
}




//__________________________________________________________________________
void AliAnalysisTaskSEDstoK0sK::FillTheTree(AliAODRecoCascadeHF* dCan, AliAODEvent* aod, TClonesArray* mcArray)
{
   ///
   ///  Fill the tree for cut optimisation studies, if fFillNtuple is enabled
   ///
   AliDebug(2, "AliAnalysisTaskSEDstoK0sK::FillTheTree()");


   if (!fFillNtuple) return;

   AliAODv0 *v0 = (AliAODv0*) dCan->Getv0();
   if (!v0) return;


   //---------------------------------------------------------------------------
   // - Get PID informations if possible
   //---------------------------------------------------------------------------
   Double_t nSigmaTPCpion   = -999;
   Double_t nSigmaTPCkaon   = -999;
   Double_t nSigmaTPCproton = -999;
   Double_t nSigmaTOFpion   = -999;
   Double_t nSigmaTOFkaon   = -999;
   Double_t nSigmaTOFproton = -999;
   if (fAnalysisCuts->GetIsUsePID()) {
      if (fAnalysisCuts->GetPidHF()->CheckTPCPIDStatus(dCan->GetBachelor())) {
         fAnalysisCuts->GetPidHF()->GetnSigmaTPC(dCan->GetBachelor(), AliPID::kPion,   nSigmaTPCpion);
         fAnalysisCuts->GetPidHF()->GetnSigmaTPC(dCan->GetBachelor(), AliPID::kKaon,   nSigmaTPCkaon);
         fAnalysisCuts->GetPidHF()->GetnSigmaTPC(dCan->GetBachelor(), AliPID::kProton, nSigmaTPCproton);
      }
      if (fAnalysisCuts->GetPidHF()->CheckTOFPIDStatus(dCan->GetBachelor())) {
         fAnalysisCuts->GetPidHF()->GetnSigmaTOF(dCan->GetBachelor(), AliPID::kPion,   nSigmaTOFpion);
         fAnalysisCuts->GetPidHF()->GetnSigmaTOF(dCan->GetBachelor(), AliPID::kKaon,   nSigmaTOFkaon);
         fAnalysisCuts->GetPidHF()->GetnSigmaTOF(dCan->GetBachelor(), AliPID::kProton, nSigmaTOFproton);
      }
   }


   //---------------------------------------------------------------------------
   // - Apply cuts to the tuple variables
   //---------------------------------------------------------------------------
   Double_t diffIP[2] = {0.}, errdiffIP[2] = {1.};
   dCan->Getd0MeasMinusExpProng(0, aod->GetMagneticField(), diffIP[0], errdiffIP[0]);
   dCan->Getd0MeasMinusExpProng(1, aod->GetMagneticField(), diffIP[1], errdiffIP[1]);


   const Int_t nVarTuple = fReadMC ? kNTupleVarsMC : kNTupleVars;

   Float_t variableNtuple[nVarTuple];
   variableNtuple[ 0] = v0->MassK0Short();                      // "K0sInvMass"
   variableNtuple[ 1] = dCan->PtProng(1);                       // "K0sPt"
   variableNtuple[ 2] = v0->RadiusV0();                         // "K0sRxy"
   variableNtuple[ 3] = v0->Ct(310, dCan->GetOwnPrimaryVtx());  // "K0sCTau"
   variableNtuple[ 4] = dCan->CosV0PointingAngle();             // "K0sCosPA"
   variableNtuple[ 5] = dCan->Getd0Prong(1);                    // "K0sd0"
   variableNtuple[ 6] = v0->DcaPosToPrimVertex();               // "K0sd0DauPos"
   variableNtuple[ 7] = v0->DcaNegToPrimVertex();               // "K0sd0DauNeg"
   variableNtuple[ 8] = v0->DcaV0Daughters();                   // "K0sDCAdau"
   variableNtuple[ 9] = v0->PtProng(0);                         // "K0sPtDauPos"
   variableNtuple[10] = v0->PtProng(1);                         // "K0sPtDauNeg"
   variableNtuple[11] = dCan->PtProng(0);                       // "BachPt"
   variableNtuple[12] = dCan->Getd0Prong(0);                    // "Bachd0"
   variableNtuple[13] = nSigmaTPCpion;                          // "BachNsigmaTPCpion"
   variableNtuple[14] = nSigmaTPCkaon;                          // "BachNsigmaTPCkaon"
   variableNtuple[15] = nSigmaTPCproton;                        // "BachNsigmaTPCproton"
   variableNtuple[16] = nSigmaTOFpion;                          // "BachNsigmaTOFpion"
   variableNtuple[17] = nSigmaTOFkaon;                          // "BachNsigmaTOFkaon"
   variableNtuple[18] = nSigmaTOFproton;                        // "BachNsigmaTOFproton"
   variableNtuple[19] = dCan->InvMassDstoK0sK();                // "CanInvMassDs"
   variableNtuple[20] = dCan->InvMassDplustoK0spi();            // "CanInvMassDplus"
   variableNtuple[21] = dCan->Pt();                             // "CanPt"
   variableNtuple[22] = dCan->GetDCA();                         // "CanDCAProngToProng"
   variableNtuple[23] = dCan->CosThetaStar(1, 431, 321, 310);   // "CanCosThetaStarK0s"
   variableNtuple[24] = dCan->CosThetaStar(0, 431, 321, 310);   // "CanCosThetaStarBach"
   variableNtuple[25] = dCan->CosPointingAngle();               // "CanCosPA"
   variableNtuple[26] = dCan->CosPointingAngleXY();             // "CanCosPAxy"
   variableNtuple[27] = dCan->DecayLengthXY();                  // "CanDLengthXY"
   variableNtuple[28] = dCan->NormalizedDecayLengthXY();        // "CanNormDLengthXY"
   variableNtuple[29] = dCan->DecayLength();                    // "CanDLength3D"
   variableNtuple[30] = dCan->NormalizedDecayLength();          // "CanNormDLength3D"
   variableNtuple[31] = ComputeSigmaVert(aod, dCan);            // "CanSigmaVtx"
   variableNtuple[32] = diffIP[0]/errdiffIP[0];                 // "CanNormTopoBach"
   variableNtuple[33] = diffIP[1]/errdiffIP[1];                 // "CanNormTopoK0s"
   // variableNtuple[32] = CosThetaK0sBachRFrame(dCan);            // "CanCosThetaRFrame"


   for (Int_t ivar=0; ivar<kNTupleVars; ivar++) {
      if ( (variableNtuple[ivar]<fCutsMinTupleVars[ivar]) || (variableNtuple[ivar]>fCutsMaxTupleVars[ivar]) ) {
         AliDebug(2, "One of the tuple variable does not pass the selection cut");
         return;
      }
   }


   //---------------------------------------------------------------------------
   // - Monte Carlo: check the candidate at MC level
   //---------------------------------------------------------------------------
   if (fReadMC) {
      Bool_t  okMCsignal = (MatchToMCDstoK0sKSignal(dCan, mcArray)>=0)     ? kTRUE : kFALSE;
      Bool_t  okMCreflec = (MatchToMCDplustoK0spiSignal(dCan, mcArray)>=0) ? kTRUE : kFALSE;

      variableNtuple[34] = okMCsignal ? 10   // "CanNormDLengthXY"
                         : okMCreflec ? 20
                         : -2;

      if (!okMCsignal && !okMCreflec) {
         AliDebug(2, "The candidate is a combinatorial background, check if the V0 is a true K0S");
         Int_t pdgDgK0stoPions[2] = {211, 211};
         if (v0->MatchToMC(310, mcArray, 2, pdgDgK0stoPions)>=0) variableNtuple[29]++;
      }

   } // end if fReadMC


   //---------------------------------------------------------------------------
   // - Fill the tuple
   //---------------------------------------------------------------------------
   fOutputNtuple->Fill(variableNtuple);


   return;
}




//__________________________________________________________________________
Float_t AliAnalysisTaskSEDstoK0sK::ComputeSigmaVert(const AliAODEvent* aod, AliAODRecoCascadeHF* dCan) const
{
   ///
   ///  Compute the track dispersion around secondary vertex starting from tracks
   ///  (Source: AliAODRecoDecayHF3Prong::ComputeSigmaVert)
   ///


   // - Get bachelor and V0 from the cascade
   AliESDtrack *esdB = new AliESDtrack(dCan->GetBachelor());

   AliNeutralTrackParam *trackV0 = 0;

   const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(dCan->Getv0());
   if (trackVV0) {
      trackV0 = new AliNeutralTrackParam(trackVV0);
   }


   if (!esdB || !trackV0) {
      if (esdB)    delete esdB;
      if (trackV0) delete trackV0;
      return -999.;
   }


   // - Build ESD primary vertex
   AliVertexerTracks vertexer(aod->GetMagneticField());

   Double_t pos[3], cov[6];
   AliAODVertex *aodV = (AliAODVertex*) aod->GetPrimaryVertex();
   aodV->GetXYZ(pos);
   aodV->GetCovarianceMatrix(cov);
   Double_t chi2 = aodV->GetChi2();
   Int_t nContr  = aodV->GetNContributors();

   AliESDVertex vprim(pos, cov, chi2, nContr);
   vertexer.SetVtxStart(&vprim);


   // - Assign TObjArray of daughters
   TObjArray twoTrackArray(2);
   twoTrackArray.AddAt(esdB, 0);
   twoTrackArray.AddAt(trackV0, 1);


   // - Build ESD secondary vertex and get dispersion
   AliESDVertex *secVert = vertexer.VertexForSelectedESDTracks(&twoTrackArray, kFALSE, kTRUE, kFALSE);
   Double_t disp = secVert->GetDispersion();


   twoTrackArray.Delete();
   delete secVert;

   return disp;
}




//__________________________________________________________________________
Float_t AliAnalysisTaskSEDstoK0sK::CosThetaK0sBachRFrame(AliAODRecoCascadeHF* dCan) const
{
   ///
   ///  Computes cosine of angle between K+ and K0S in the Ds rest frame
   ///  (Source: AliAODRecoDecayHF3Prong::CosPiKPhiRFrame)
   ///


   // - Get bachelor and V0 from the cascade
   AliAODv0    *v0   = dynamic_cast<AliAODv0*> (dCan->Getv0());
   AliAODTrack *bach = dynamic_cast<AliAODTrack*> (dCan->GetBachelor());
   if (!v0 || !bach) {
      return -999;
   }


   // - Get cascade informations necessary to compute the rest frame
   Double_t eCasc  = dCan->EProng(0, 321) + dCan->EProng(1, 310);
   Double_t pxCasc = dCan->PxProng(0) + dCan->PxProng(1);
   Double_t pyCasc = dCan->PyProng(0) + dCan->PyProng(1);
   Double_t pzCasc = dCan->PzProng(0) + dCan->PzProng(1);
   Double_t bxCasc = pxCasc/eCasc;
   Double_t byCasc = pyCasc/eCasc;
   Double_t bzCasc = pzCasc/eCasc;


   // - Boost the K+ and K0s in the cascade rest frame
   TVector3 vecK0sCascFrame;
   TLorentzVector vecK0s(dCan->PxProng(1), dCan->PyProng(1), dCan->PzProng(1), dCan->EProng(1, 310));
   vecK0s.Boost(-bxCasc, -byCasc, -bzCasc);
   vecK0s.Boost(vecK0sCascFrame);
   vecK0sCascFrame = vecK0s.BoostVector();

   TVector3 vecBachCascFrame;
   TLorentzVector vecBach(dCan->PxProng(0), dCan->PyProng(0), dCan->PzProng(0), dCan->EProng(0, 321));
   vecBach.Boost(-bxCasc, -byCasc, -bzCasc);
   vecBach.Boost(vecBachCascFrame);
   vecBachCascFrame = vecBach.BoostVector();


   // - Compute the cosine of the angle between K+ and K0s
   Double_t scalProdBachK0s   = vecBachCascFrame.Dot(vecK0sCascFrame);
   Double_t normBach          = vecBachCascFrame.Mag();
   Double_t normK0s           = vecK0sCascFrame.Mag();
   Double_t cosAngleCascFrame = scalProdBachK0s/(normBach*normK0s);


   return cosAngleCascFrame;
}

