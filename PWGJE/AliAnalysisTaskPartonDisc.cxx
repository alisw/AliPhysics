/************************************************************************* 
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

////////////////////////////////////////////////////////
//                                                    //
// Analysis task for parton discrimination studies    //
//                                                    // 
// Author:                                            //
// Hermes Leon Vargas  (hleon@ikf.uni-frankfurt.de)   //
//                                                    // 
////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TList.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include <algorithm> 
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODHeader.h"
#include "AliESDtrack.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisTaskPartonDisc.h"
#include "AliAODVZERO.h"

// Analysis task for parton discrimination

ClassImp(AliAnalysisTaskPartonDisc)

Double_t *AliAnalysisTaskPartonDisc::fgContainer = 0x0; 

//________________________________________________________________________
AliAnalysisTaskPartonDisc::AliAnalysisTaskPartonDisc() 
  : AliAnalysisTaskSE(), fAOD(0), fUseAODMC(kFALSE), fPhojetMC(kFALSE), fBranchMC("jetsMC"), fBranchRec("jetsREC"), fBranchSecRec(""), fSqrts(0),  fNtX(0), fJetRadius(0.), fFlavorRadius(0.), fFilterBit(0xFF), fOutputList(0), fJetPt(0), fJetPtSec(0), fJetPtMC(0), fJetEta(0), fJetEtaSec(0), fJetPhi(0), fJetPhiSec(0), fJetEtaMC(0), fJetPhiMC(0), fPtAODMC(0), fPtAOD(0), fEtaAODMC(0), fPhiAODMC(0), fEtaAOD(0), fPhiAOD(0), fFlavor(0), fNJetsMC(0), fNJetsRD(0), fNJetsRDSeco(0), fJetsMultPtMC(0), fJetsMultPtRD(0), fNChTrRD(0), fProfNChTrRD(0), fFracQQ(0), fFracGQ(0), fFracGG(0), fFracOutGoingQQ(0), fFracOutGoingGQ(0), fFracOutGoingGG(0), fh1Xsec(0), fh1Trials(0), fMpdg(0), fProcessJetPt(0), fFlavorLead(0), fProcessLeadJetPt(0), fPDGMothLPart(0), fFlavProc(0), fAvgTrials(1), fUseAODJetInput(kFALSE), fMinTrackPtInNTX(0), fMaxTrackPtInNTX(0), fSCMRD(0), fMinpTVal(0), fZVertex(0), fh1Events(0), fUseOnlyMC(kFALSE), fCheckMCStatus(kTRUE), fEvtCount(0), fNAccJetsMC(0), fNAccJetsRD(0), fNAccJetsRDSeco(0), fEnablePrints(kFALSE), fRecJetPtInclusive(0), fMCJetPtInclusive(0), fRecJetPtLeading(0), fMCJetPtLeading(0), fSecRecJetPtInclusive(0), fSecRecJetPtLeading(0), fHasPerpCone(kTRUE), fEtaPerpCoord(0), fPhiPerpCoord(0), fPtPerpCoord(0), fJetEvent(kFALSE), fPerpCone(0), fNChTrMCPerp(0), fNChTrRecPerp(0), fSCMMCPerp(0), fSCMRecPerp(0), fIsHIevent(kFALSE), fCurrentJetMinPtNT90(0), fBckgSbsJet(0), fCurrentJetMinPtNT90Recalc(0), fNChTrCorrMCQuark(0), fNChTrCorrMCGluon(0), fNChTrCorrMCPerp(0), fIsPossibleToSubstBckg(kTRUE), fNChTrRecECorr(0), fNChTrRecPerpECorr(0), fRefMult(0), fCurrentJetCharge(0), fRefMultWOJet(0), fVZEROMult(0), fMultWOJetVZero(0), fVZero(0), fRefMultFullV0(0), fRefMultV0Corr(0), fFullV0V0Corr(0), fNTXV0MultPt(0), fNTXCBMultPt(0), fMinpTValUE(2.0), fRefMultFullV0UJ(0), fRefMultV0CorrUJ(0), fFullV0V0CorrUJ(0), fMultWOJetVZeroUJ(0), fRefMultWOJetUJ(0), fMaxpTValUE(2.0), fRefAODTrackCount(0), fRefAODTrackCountUJ(0), fTrackCountWOJet(0), fTrackCountWOJetUJ(0), fTrackCountWOJetUJMC(0), fFullV0V0CorrUJMC(0), fMinpTValMC(2.0), fIncExcR(0.0), fForceNotTR(kFALSE), fNotExtDiJEx(kFALSE), fMinTrackPtInNTXRecalc(0), fMaxTrackPtInNTXRecalc(0), fPtDistInJetConeRaw(0), fPtDistInPerpConeRaw(0), fPtInPerpCon(0), fMinTrackPtInNTXR(0), fMaxTrackPtInNTXR(0), fEventCent(0), fJetEtaAll(0), fJetEtaOnlyTPCcut(0), fNChTrRecECorrPPMult(0), fNChTrRecPerpECorrPPMult(0), fForceSkipSJ(kFALSE), fJetPtCentPbPbRaw(0), fJetPtCentPbPbCorr(0), fJetAcceptance(0.5), fIncreasingExcl(kFALSE)
{

  // Constructor

  for(Int_t a=0; a<16;a++)
    {
      fJetFlags[a]=kTRUE;
      if(a<12)
	{
	  fNChTr[a]=0;
	  fHistPtParton[a]=0;
	  fSCM[a]=0;
	  if(a<8)
	    {
	      fNChTrRDMult[a]=0;
	      fNAccJetsRDMult[a]=0;
	      fTotalJetCharge[a]=0;
	      fSCMRDMult[a]=0;
	      fNChTrRDMultMC[a]=0;
	      fSCMRDMultMC[a]=0;
	      fNChTrRDMultSE[a]=0;
	      fNAccJetsRDMultSE[a]=0;
	      fTotalJetChargeSE[a]=0;
	      fSCMRDMultSE[a]=0;
	      fNChTrRDMultOJ[a]=0;
	      fSCMRDMultOJ[a]=0;
	      fNChTrRDMultSEOJ[a]=0;
	      fSCMRDMultSEOJ[a]=0;
	      fNChTrRDMultOJMC[a]=0;
	      fSCMRDMultOJMC[a]=0;
	      fNChTrRDMultSEOJMC[a]=0;
	      fSCMRDMultSEOJMC[a]=0;
	      fNChTrRecPerpMultSEOJ[a]=0;
	    }
	  if(a<6)
	    {
	      fProcessPDG[a]=0;
	      fFragPion[a]=0;
	      fFragKaon[a]=0;
	      fFragProton[a]=0;
	      fFragChargedR4[a]=0;
	      fFragChargedR3[a]=0;
	      fFragChargedR2[a]=0;
	      fHistContainerR4[a]=0;
	      fHistContainerR3[a]=0;
	      fHistContainerR2[a]=0;
	      if(a<3)
		{
		  fJetEtaJetPt[a]=0;
		  if(a<2)
		    {
		      fFragCandidates[a]=0;
		      fMinTrackPtInNTXh[a]=0;
		      fMaxTrackPtInNTXh[a]=0; 
		    }
		}
	    }
	}
    }
}
//________________________________________________________________________
AliAnalysisTaskPartonDisc::AliAnalysisTaskPartonDisc(const char *name) 
  : AliAnalysisTaskSE(name), fAOD(0), fUseAODMC(kFALSE), fPhojetMC(kFALSE), fBranchMC("jetsMC"), fBranchRec("jetsREC"), fBranchSecRec(""), fSqrts(0),  fNtX(0), fJetRadius(0.), fFlavorRadius(0.), fFilterBit(0xFF), fOutputList(0), fJetPt(0), fJetPtSec(0), fJetPtMC(0), fJetEta(0), fJetEtaSec(0), fJetPhi(0), fJetPhiSec(0), fJetEtaMC(0), fJetPhiMC(0), fPtAODMC(0), fPtAOD(0), fEtaAODMC(0), fPhiAODMC(0), fEtaAOD(0), fPhiAOD(0), fFlavor(0), fNJetsMC(0), fNJetsRD(0), fNJetsRDSeco(0), fJetsMultPtMC(0), fJetsMultPtRD(0), fNChTrRD(0), fProfNChTrRD(0), fFracQQ(0), fFracGQ(0), fFracGG(0), fFracOutGoingQQ(0), fFracOutGoingGQ(0), fFracOutGoingGG(0), fh1Xsec(0), fh1Trials(0), fMpdg(0), fProcessJetPt(0), fFlavorLead(0), fProcessLeadJetPt(0), fPDGMothLPart(0), fFlavProc(0), fAvgTrials(1), fUseAODJetInput(kFALSE), fMinTrackPtInNTX(0), fMaxTrackPtInNTX(0), fSCMRD(0), fMinpTVal(0), fZVertex(0), fh1Events(0), fUseOnlyMC(kFALSE), fCheckMCStatus(kTRUE), fEvtCount(0), fNAccJetsMC(0), fNAccJetsRD(0), fNAccJetsRDSeco(0), fEnablePrints(kFALSE), fRecJetPtInclusive(0), fMCJetPtInclusive(0), fRecJetPtLeading(0), fMCJetPtLeading(0), fSecRecJetPtInclusive(0), fSecRecJetPtLeading(0), fHasPerpCone(kTRUE), fEtaPerpCoord(0), fPhiPerpCoord(0), fPtPerpCoord(0), fJetEvent(kFALSE), fPerpCone(0), fNChTrMCPerp(0), fNChTrRecPerp(0), fSCMMCPerp(0), fSCMRecPerp(0), fIsHIevent(kFALSE), fCurrentJetMinPtNT90(0), fBckgSbsJet(0), fCurrentJetMinPtNT90Recalc(0), fNChTrCorrMCQuark(0), fNChTrCorrMCGluon(0), fNChTrCorrMCPerp(0), fIsPossibleToSubstBckg(kTRUE), fNChTrRecECorr(0), fNChTrRecPerpECorr(0), fRefMult(0), fCurrentJetCharge(0), fRefMultWOJet(0), fVZEROMult(0), fMultWOJetVZero(0), fVZero(0), fRefMultFullV0(0), fRefMultV0Corr(0), fFullV0V0Corr(0), fNTXV0MultPt(0), fNTXCBMultPt(0), fMinpTValUE(2.0), fRefMultFullV0UJ(0), fRefMultV0CorrUJ(0), fFullV0V0CorrUJ(0), fMultWOJetVZeroUJ(0), fRefMultWOJetUJ(0), fMaxpTValUE(2.0), fRefAODTrackCount(0), fRefAODTrackCountUJ(0), fTrackCountWOJet(0), fTrackCountWOJetUJ(0), fTrackCountWOJetUJMC(0), fFullV0V0CorrUJMC(0), fMinpTValMC(2.0), fIncExcR(0.0), fForceNotTR(kFALSE), fNotExtDiJEx(kFALSE), fMinTrackPtInNTXRecalc(0), fMaxTrackPtInNTXRecalc(0), fPtDistInJetConeRaw(0), fPtDistInPerpConeRaw(0), fPtInPerpCon(0), fMinTrackPtInNTXR(0), fMaxTrackPtInNTXR(0), fEventCent(0), fJetEtaAll(0), fJetEtaOnlyTPCcut(0), fNChTrRecECorrPPMult(0), fNChTrRecPerpECorrPPMult(0), fForceSkipSJ(kFALSE), fJetPtCentPbPbRaw(0), fJetPtCentPbPbCorr(0), fJetAcceptance(0.5), fIncreasingExcl(kFALSE)
{

  // Constructor

  for(Int_t a=0; a<16;a++)
    {
      fJetFlags[a]=kTRUE;
      if(a<12)
	{
	  fNChTr[a]=0;
	  fHistPtParton[a]=0;
	  fSCM[a]=0;
	  if(a<8)
	    {
	      fNChTrRDMult[a]=0;
	      fNAccJetsRDMult[a]=0;
	      fTotalJetCharge[a]=0;
	      fSCMRDMult[a]=0;
	      fNChTrRDMultMC[a]=0;
	      fSCMRDMultMC[a]=0;
	      fNChTrRDMultSE[a]=0;
	      fNAccJetsRDMultSE[a]=0;
	      fTotalJetChargeSE[a]=0;
	      fSCMRDMultSE[a]=0;
	      fNChTrRDMultOJ[a]=0;
	      fSCMRDMultOJ[a]=0;
	      fNChTrRDMultSEOJ[a]=0;
	      fSCMRDMultSEOJ[a]=0;
	      fNChTrRDMultOJMC[a]=0;
	      fSCMRDMultOJMC[a]=0;
	      fNChTrRDMultSEOJMC[a]=0;
	      fSCMRDMultSEOJMC[a]=0;
	      fNChTrRecPerpMultSEOJ[a]=0;		
	    }
	  if(a<6)
	    {
	      fProcessPDG[a]=0;
	      fFragPion[a]=0;
	      fFragKaon[a]=0;
	      fFragProton[a]=0;
	      fFragChargedR4[a]=0;
	      fFragChargedR3[a]=0;
	      fFragChargedR2[a]=0;
	      fHistContainerR4[a]=0;
	      fHistContainerR3[a]=0;
	      fHistContainerR2[a]=0;
	      if(a<3)
		{
		  fJetEtaJetPt[a]=0;
		  if(a<2)
		    {
		      fFragCandidates[a]=0;
		      fMinTrackPtInNTXh[a]=0;
		      fMaxTrackPtInNTXh[a]=0; 
		    }
		}
	    }
	}
    }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::UserNotify()
{
  //
  // read the cross sections
  // and number of trials from pyxsec.root
  // from AliAnalysisTaskJetSpectrum2
  //

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;
  Int_t nevents = 0;

  fAvgTrials = 1;
  if(tree)
    {
      TFile *curfile = tree->GetCurrentFile();
      if (!curfile) 
	{
	  Error("Notify","No current file");
	  return kFALSE;
	}
      if(!fh1Xsec||!fh1Trials)
	{
	  Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
	  return kFALSE;
	}
      AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
      fh1Xsec->Fill("<#sigma>",xsection);
      // construct a poor man average trials (per event!?)
      Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
      if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
      // number of events read out to create the AOD
      NumberOfReadEventsAOD(curfile->GetName(),nevents);
      fh1Events->Fill("#sum{nevents}",nevents); //  filled once per file
    }  
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisTaskPartonDisc::UserCreateOutputObjects()
{
  // Create histograms
  // Called once  

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fJetPt = new TH1F("fJetPt", "p_{T} distribution of reco jets", 60, 0., 300.);
  fJetPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fJetPt->GetYaxis()->SetTitle("dN/dp_{T} (c/GeV)");
  fJetPt->Sumw2();
  fOutputList->Add(fJetPt);

  fJetPtSec = new TH1F("fJetPtSec", "p_{T} distribution of reco jets, seconday rec branch", 60, 0., 300.);
  fJetPtSec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fJetPtSec->GetYaxis()->SetTitle("dN/dp_{T} (c/GeV)");
  fJetPtSec->Sumw2();
  fOutputList->Add(fJetPtSec);

  fJetPtMC = new TH1F("fJetPtMC", "p_{T} distribution of mc jets", 60, 0., 300.);
  fJetPtMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fJetPtMC->GetYaxis()->SetTitle("dN/dp_{T} (c/GeV)");
  fJetPtMC->Sumw2();
  fOutputList->Add(fJetPtMC);

  fJetEta = new TH2F("fJetEta", "Eta distribution of reconstructed jets", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fJetEta->GetXaxis()->SetTitle("#eta");
  fJetEta->GetYaxis()->SetTitle("#eta");
  fJetEta->Sumw2();
  fOutputList->Add(fJetEta);

  fJetEtaSec = new TH2F("fJetEtaSec", "Eta distribution of reconstructed jets, secondary branch", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fJetEtaSec->GetXaxis()->SetTitle("#eta");
  fJetEtaSec->GetYaxis()->SetTitle("#eta");
  fJetEtaSec->Sumw2();
  fOutputList->Add(fJetEtaSec);

  fJetPhi = new TH2F("fJetPhi", "Phi distribution of reconstructed jets", 50, 0., TMath::TwoPi(), 50, 0., TMath::TwoPi());
  fJetPhi->GetXaxis()->SetTitle("#phi");
  fJetPhi->GetYaxis()->SetTitle("#phi");
  fJetPhi->Sumw2();
  fOutputList->Add(fJetPhi);

  fJetPhiSec = new TH2F("fJetPhiSec", "Phi distribution of reconstructed jets, secondary branch", 50, 0., TMath::TwoPi(), 50, 0., TMath::TwoPi());
  fJetPhiSec->GetXaxis()->SetTitle("#phi");
  fJetPhiSec->GetYaxis()->SetTitle("#phi");
  fJetPhiSec->Sumw2();
  fOutputList->Add(fJetPhiSec);

  fJetEtaMC = new TH2F("fJetEtaMC", "Eta distribution of MC jets", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fJetEtaMC->GetXaxis()->SetTitle("#eta");
  fJetEtaMC->GetYaxis()->SetTitle("#eta");
  fJetEtaMC->Sumw2();
  fOutputList->Add(fJetEtaMC);

  fJetPhiMC = new TH2F("fJetPhiMC", "Phi distribution of MC jets", 50, 0., TMath::TwoPi(), 50, 0., TMath::TwoPi());
  fJetPhiMC->GetXaxis()->SetTitle("#phi");
  fJetPhiMC->GetYaxis()->SetTitle("#phi");
  fJetPhiMC->Sumw2();
  fOutputList->Add(fJetPhiMC);

  fPtAODMC = new TH2F("fPtAODMC", "p_{T} distribution of mc tracks", 200, 0., 100., 200, 0., 100.);
  fPtAODMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fPtAODMC->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fPtAODMC->Sumw2();
  fOutputList->Add(fPtAODMC);

  fPtAOD = new TH2F("fPtAOD", "p_{T} distribution of aod tracks", 200, 0., 100., 200, 0., 100.);
  fPtAOD->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fPtAOD->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fPtAOD->Sumw2();
  fOutputList->Add(fPtAOD);

  fEtaAODMC = new TH2F("fEtaAODMC", "Eta distribution of MC aod tracks", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fEtaAODMC->GetXaxis()->SetTitle("#eta");
  fEtaAODMC->GetYaxis()->SetTitle("#eta");
  fEtaAODMC->Sumw2();
  fOutputList->Add(fEtaAODMC);

  fPhiAODMC = new TH2F("fPhiAODMC", "Phi distribution of MC aod tracks", 50, 0., TMath::TwoPi(), 50, 0., TMath::TwoPi());
  fPhiAODMC->GetXaxis()->SetTitle("#phi");
  fPhiAODMC->GetYaxis()->SetTitle("#phi");
  fPhiAODMC->Sumw2();
  fOutputList->Add(fPhiAODMC);

  fEtaAOD = new TH2F("fEtaAOD", "Eta distribution of aod tracks", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fEtaAOD->GetXaxis()->SetTitle("#eta");
  fEtaAOD->GetYaxis()->SetTitle("#eta");
  fEtaAOD->Sumw2();
  fOutputList->Add(fEtaAOD);

  fPhiAOD = new TH2F("fPhiAOD", "Phi distribution of aod tracks", 50, 0.,TMath::TwoPi(), 50, 0.,TMath::TwoPi());
  fPhiAOD->GetXaxis()->SetTitle("#phi");
  fPhiAOD->GetYaxis()->SetTitle("#phi");
  fPhiAOD->Sumw2();
  fOutputList->Add(fPhiAOD);

  fFlavor = new TH2F("fFlavor", "Flavor distribution of jets", 27, -5.5, 21.5, 60, 0., 300.);
  fFlavor->GetXaxis()->SetTitle("PDG code");
  fFlavor->GetYaxis()->SetTitle("p_{T}^{JET}");
  fFlavor->Sumw2();
  fOutputList->Add(fFlavor);

  fNJetsMC = new TH2F("fNJetsMC", "Number of simulated jets per event, as a function of p_T", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNJetsMC->GetXaxis()->SetTitle("Number of jets");
  fNJetsMC->GetYaxis()->SetTitle("Number of jets");
  fNJetsMC->Sumw2();
  fOutputList->Add(fNJetsMC);

  fNJetsRD = new TH2F("fNJetsRD", "Number of jets per event in real data", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNJetsRD->GetXaxis()->SetTitle("Number of jets");
  fNJetsRD->GetYaxis()->SetTitle("Number of jets");
  fNJetsRD->Sumw2();
  fOutputList->Add(fNJetsRD);

  fNJetsRDSeco = new TH2F("fNJetsRDSeco", "Number of jets per event in real data, secondary branch", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNJetsRDSeco->GetXaxis()->SetTitle("Number of jets");
  fNJetsRDSeco->GetYaxis()->SetTitle("Number of jets");
  fNJetsRDSeco->Sumw2();
  fOutputList->Add(fNJetsRDSeco);

  fJetsMultPtMC = new TH2F("fJetsMultPtMC", "Jet multiplicity associated to jet pT, as a function of p_T MC", 60, 0., 300., 7, -0.5, 6.5);
  fJetsMultPtMC->GetXaxis()->SetTitle("p_{T}^{JET}");
  fJetsMultPtMC->GetYaxis()->SetTitle("Number of jets");
  fJetsMultPtMC->Sumw2();
  fOutputList->Add(fJetsMultPtMC);

  fJetsMultPtRD = new TH2F("fJetsMultPtRD", "Jet multiplicity associated to jet pT, as a function of p_T Reco Data",60, 0., 300., 7, -0.5, 6.5);
  fJetsMultPtRD->GetXaxis()->SetTitle("p_{T}^{JET}");
  fJetsMultPtRD->GetYaxis()->SetTitle("Number of jets");
  fJetsMultPtRD->Sumw2();
  fOutputList->Add(fJetsMultPtRD);

  fNChTrRD = new TH2F("fNChTrRD","Number of tracks to recover transverse energy, jet_{p_{T}}",101,-0.5,100.5, 60, 0., 300.);
  fNChTrRD->GetXaxis()->SetTitle("NTracks");
  fNChTrRD->GetYaxis()->SetTitle("p_{T}^{JET}");
  fNChTrRD->Sumw2();
  fOutputList->Add(fNChTrRD);

  fProfNChTrRD = new TProfile("fProfNChTrRD","Number of tracks to recover transverse energy, jet_{p_{T}}", 50, 0, 250);
  fProfNChTrRD->GetXaxis()->SetTitle("p_{T}^{JET}");
  fProfNChTrRD->GetYaxis()->SetTitle("NTracks");
  fProfNChTrRD->Sumw2();
  fOutputList->Add(fProfNChTrRD);

  fFracQQ = new TH1F("fFracQQ","QQ",1000,0.,0.5);
  fFracQQ->GetXaxis()->SetTitle("x_{T}");
  fFracQQ->GetYaxis()->SetTitle("Entries");
  fFracQQ->Sumw2();
  fOutputList->Add(fFracQQ);

  fFracGQ = new TH1F("fFracGQ","GQ",1000,0.,0.5);
  fFracGQ->GetXaxis()->SetTitle("x_{T}");
  fFracGQ->GetYaxis()->SetTitle("Entries");
  fFracGQ->Sumw2();
  fOutputList->Add(fFracGQ);

  fFracGG = new TH1F("fFracGG","GG",1000,0.,0.5);
  fFracGG->GetXaxis()->SetTitle("x_{T}");
  fFracGG->GetYaxis()->SetTitle("Entries");
  fFracGG->Sumw2();
  fOutputList->Add(fFracGG);

  fFracOutGoingQQ = new TH1F("fFracOutGoingQQ","QQ",1000,0.,0.5);
  fFracOutGoingQQ->GetXaxis()->SetTitle("x_{T}");
  fFracOutGoingQQ->GetYaxis()->SetTitle("Entries");
  fFracOutGoingQQ->Sumw2();
  fOutputList->Add(fFracOutGoingQQ);

  fFracOutGoingGQ = new TH1F("fFracOutGoingGQ","GQ",1000,0.,0.5);
  fFracOutGoingGQ->GetXaxis()->SetTitle("x_{T}");
  fFracOutGoingGQ->GetYaxis()->SetTitle("Entries");
  fFracOutGoingGQ->Sumw2();
  fOutputList->Add(fFracOutGoingGQ);

  fFracOutGoingGG = new TH1F("fFracOutGoingGG","GG",1000,0.,0.5);
  fFracOutGoingGG->GetXaxis()->SetTitle("x_{T}");
  fFracOutGoingGG->GetYaxis()->SetTitle("Entries");
  fFracOutGoingGG->Sumw2();
  fOutputList->Add(fFracOutGoingGG);

  fh1Xsec =  new TProfile("h1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->SetXTitle("<#sigma>");
  fh1Xsec->Sumw2();
  fOutputList->Add(fh1Xsec);  

  fh1Trials = new TH1F("h1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->SetXTitle("#sum{ntrials}");
  fh1Trials->Sumw2();
  fOutputList->Add(fh1Trials);
           
  fProcessJetPt = new TH2F("fProcessJetPt","Pythia process number, jet_{p_{T}}",70, 0.5, 70.5, 60, 0., 300.);
  fProcessJetPt->GetXaxis()->SetTitle("Pythia process");
  fProcessJetPt->GetYaxis()->SetTitle("p_{T}^{JET}");
  fProcessJetPt->Sumw2();
  fOutputList->Add(fProcessJetPt);

  fFlavorLead = new TH2F("fFlavorLead", "Flavor distribution of leading jets", 27, -5.5, 21.5, 60, 0., 300.);
  fFlavorLead->GetXaxis()->SetTitle("PDG code");
  fFlavorLead->GetYaxis()->SetTitle("p_{T}^{JET}");
  fFlavorLead->Sumw2();
  fOutputList->Add(fFlavorLead);

  fProcessLeadJetPt = new TH2F("fProcessLeadJetPt","Pythia process number for leading, jet_{p_{T}}",70, 0.5, 70.5, 60, 0., 300.);
  fProcessLeadJetPt->GetXaxis()->SetTitle("Pythia process");
  fProcessLeadJetPt->GetYaxis()->SetTitle("p_{T}^{JET}");
  fProcessLeadJetPt->Sumw2();
  fOutputList->Add(fProcessLeadJetPt);
 
  fPDGMothLPart = new TH3F("fPDGMothLPart","Mother of leading parton, leading parton, jet p_{T}", 27, -5.5, 21.5, 27, -5.5, 21.5, 60, 0., 300.);
  fPDGMothLPart->GetXaxis()->SetTitle("Mother of leading parton");
  fPDGMothLPart->GetYaxis()->SetTitle("Leading parton");
  fPDGMothLPart->GetZaxis()->SetTitle("p_{T}^{JET}");
  fPDGMothLPart->Sumw2();
  fOutputList->Add(fPDGMothLPart);
  
  fFlavProc = new TH2F("fFlavProc","Flavor, Flavor status code", 27, -5.5, 21.5, 101, -0.5, 100.5);
  fFlavProc->GetXaxis()->SetTitle("Jet flavor");
  fFlavProc->GetYaxis()->SetTitle("Parton status code");
  fFlavProc->Sumw2();
  fOutputList->Add(fFlavProc);

  fSCMRD = new TH2F("fSCMRD","Second Central Moment, jet_{p_{T}}",200,0.,0.2, 60, 0., 300.);
  fSCMRD->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
  fSCMRD->GetYaxis()->SetTitle("p_{T}^{JET}");
  fSCMRD->Sumw2();
  fOutputList->Add(fSCMRD);
 
  fZVertex = new TH2F("fZVertex","Z vertex position, Z vertex position}",40,-20.,20., 40, -20., 20.);
  fZVertex->GetXaxis()->SetTitle("Vertex z");
  fZVertex->GetYaxis()->SetTitle("Vertex z");
  fZVertex->Sumw2();
  fOutputList->Add(fZVertex);

  fh1Events = new TH1F("fh1Events","nevents from PWG4_JetTasksOutput.root",1,0,1);
  fh1Events->SetXTitle("#sum{nevents}");
  fh1Events->Sumw2();
  fOutputList->Add(fh1Events);

  fNAccJetsMC = new TH2F("fNAccJetsMC", "Number accepted simulated jets per event", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNAccJetsMC->GetXaxis()->SetTitle("Number of jets");
  fNAccJetsMC->GetYaxis()->SetTitle("Number of jets");
  fNAccJetsMC->Sumw2();
  fOutputList->Add(fNAccJetsMC);

  fNAccJetsRD = new TH2F("fNAccJetsRD", "Number of accepted jets per event in real data", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNAccJetsRD->GetXaxis()->SetTitle("Number of jets");
  fNAccJetsRD->GetYaxis()->SetTitle("Number of jets");
  fNAccJetsRD->Sumw2();
  fOutputList->Add(fNAccJetsRD);

  fNAccJetsRDSeco = new TH2F("fNAccJetsRDSeco", "Number of accepted jets per event in real data, secondary branch", 101, -0.5, 100.5, 101, -0.5, 100.5);
  fNAccJetsRDSeco->GetXaxis()->SetTitle("Number of jets");
  fNAccJetsRDSeco->GetYaxis()->SetTitle("Number of jets");
  fNAccJetsRDSeco->Sumw2();
  fOutputList->Add(fNAccJetsRDSeco);

  fRecJetPtInclusive = new TH1F("fRecJetPtInclusive", "p_{T} distribution of inclusive reco jets", 60, 0., 300.);
  fRecJetPtInclusive->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fRecJetPtInclusive->GetYaxis()->SetTitle("d#sigma (mb)");
  fRecJetPtInclusive->Sumw2();
  fOutputList->Add(fRecJetPtInclusive);

  fMCJetPtInclusive = new TH1F("fMCJetPtInclusive", "p_{T} distribution of inclusive MC jets", 60, 0., 300.);
  fMCJetPtInclusive->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fMCJetPtInclusive->GetYaxis()->SetTitle("d#sigma (mb)");
  fMCJetPtInclusive->Sumw2();
  fOutputList->Add(fMCJetPtInclusive);

  fRecJetPtLeading = new TH1F("fRecJetPtLeading", "p_{T} distribution of leading reco jets", 60, 0., 300.);
  fRecJetPtLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fRecJetPtLeading->GetYaxis()->SetTitle("d#sigma (mb)");
  fRecJetPtLeading->Sumw2();
  fOutputList->Add(fRecJetPtLeading);

  fMCJetPtLeading = new TH1F("fMCJetPtLeading", "p_{T} distribution of leading MC jets", 60, 0., 300.);
  fMCJetPtLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fMCJetPtLeading->GetYaxis()->SetTitle("d#sigma (mb)");
  fMCJetPtLeading->Sumw2();
  fOutputList->Add(fMCJetPtLeading);

  fSecRecJetPtInclusive = new TH1F("fSecRecJetPtInclusive", "p_{T} distribution of inclusive reco jets (2nd branch)", 60, 0., 300.);
  fSecRecJetPtInclusive->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fSecRecJetPtInclusive->GetYaxis()->SetTitle("d#sigma (mb)");
  fSecRecJetPtInclusive->Sumw2();
  fOutputList->Add(fSecRecJetPtInclusive);

  fSecRecJetPtLeading = new TH1F("fSecRecJetPtLeading", "p_{T} distribution of leading reco jets (2nd branch)", 60, 0., 300.);
  fSecRecJetPtLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fSecRecJetPtLeading->GetYaxis()->SetTitle("d#sigma (mb)");
  fSecRecJetPtLeading->Sumw2();
  fOutputList->Add(fSecRecJetPtLeading);

  fNChTrMCPerp = new TH2F("fNChTrMCPerp","Number of tracks to recover transverse energy of perp. cone, jet_{p_{T} MC}",101,-0.5,100.5, 60, 0., 300.);
  fNChTrMCPerp->GetXaxis()->SetTitle("NTracks Perp");
  fNChTrMCPerp->GetYaxis()->SetTitle("p_{T}^{MC JET}");
  fNChTrMCPerp->Sumw2();
  fOutputList->Add(fNChTrMCPerp);

  fNChTrRecPerp = new TH2F("fNChTrRecPerp","Number of tracks to recover transverse energy of perp. cone, jet_{p_{T} Rec}",101,-0.5,100.5, 60, 0., 300.);
  fNChTrRecPerp->GetXaxis()->SetTitle("NTracks Perp");
  fNChTrRecPerp->GetYaxis()->SetTitle("p_{T}^{RECO JET}");
  fNChTrRecPerp->Sumw2();
  fOutputList->Add(fNChTrRecPerp);

  fSCMMCPerp = new TH2F("fSCMMCPerp","Second Central Moment of perp. cone, jet_{p_{T} MC}",200,0.,0.2, 60, 0., 300.);
  fSCMMCPerp->GetXaxis()->SetTitle("<#delta R_{c}^{2}> Perp");
  fSCMMCPerp->GetYaxis()->SetTitle("p_{T}^{JET}");
  fSCMMCPerp->Sumw2();
  fOutputList->Add(fSCMMCPerp);

  fSCMRecPerp = new TH2F("fSCMRecPerp","Second Central Moment of perp. cone, jet_{p_{T} Reco}",200,0.,0.2, 60, 0., 300.);
  fSCMRecPerp->GetXaxis()->SetTitle("<#delta R_{c}^{2}> Perp");
  fSCMRecPerp->GetYaxis()->SetTitle("p_{T}^{JET}");
  fSCMRecPerp->Sumw2();
  fOutputList->Add(fSCMRecPerp);

  fNChTrCorrMCQuark = new TH2F("fNChTrCorrMCQuark","Number of tracks to recover corrected transverse energy, MC quarks",101,-0.5,100.5, 60, 0., 300.);
  fNChTrCorrMCQuark->GetXaxis()->SetTitle("NTracks");
  fNChTrCorrMCQuark->GetYaxis()->SetTitle("p_{T}^{MC Corr. JET}");
  fNChTrCorrMCQuark->Sumw2();
  fOutputList->Add(fNChTrCorrMCQuark);

  fNChTrCorrMCGluon = new TH2F("fNChTrCorrMCGluon","Number of tracks to recover corrected transverse energy, MC gluons",101,-0.5,100.5, 60, 0., 300.);
  fNChTrCorrMCGluon->GetXaxis()->SetTitle("NTracks");
  fNChTrCorrMCGluon->GetYaxis()->SetTitle("p_{T}^{MC Corr. JET}");
  fNChTrCorrMCGluon->Sumw2();
  fOutputList->Add(fNChTrCorrMCGluon);

  fNChTrCorrMCPerp = new TH2F("fNChTrCorrMCPerp","Number of tracks to recover perp. cone. after corrected jet pT",101,-0.5,100.5, 60, 0., 300.);
  fNChTrCorrMCPerp->GetXaxis()->SetTitle("NTracks");
  fNChTrCorrMCPerp->GetYaxis()->SetTitle("p_{T}^{MC Corr. JET}");
  fNChTrCorrMCPerp->Sumw2();
  fOutputList->Add(fNChTrCorrMCPerp);

  // 9 selection bins: (nuevo)
  // 1st. Proton collisions                  fill 1  Bin1 [0.5,1.5)
  // 2nd. PbPb collisions, Bin  0-10         fill 2  Bin2 [1.5,2.5)
  // 3rd. PbPb collisions, Bin 10-20         fill 3  Bin3 [2.5,3.5)
  // 4rd. PbPb collisions, Bin 20-30         fill 4  Bin4 [3.5,4.5)
  // 5th. PbPb collisions, Bin 30-40         fill 5  Bin5 [4.5,5.5)
  // 6th. PbPb collisions, Bin 40-50         fill 6  Bin6 [5.5,6.5)
  // 7th. PbPb collisions, Bin 50-60         fill 7  Bin7 [6.5,7.5)
  // 8th. PbPb collisions, Bin 60-70         fill 8  Bin8 [7.5,8.5)
  // 9th. PbPb collisions, Bin 70-80         fill 9  Bin9 [8.5,9.5)
  // 10th. PbPb collisions, Bin 80-100.1    fill 10  Bin10 [9.5,10.5)

  fNChTrRecECorr = new TH3F("fNChTrRecECorr","NTX in ener. corr. jet , corr. jet pT, centrality",101,-0.5,100.5, 60, 0., 300.,10,0.5,10.5);
  fNChTrRecECorr->GetXaxis()->SetTitle("NTracks");
  fNChTrRecECorr->GetYaxis()->SetTitle("p_{T}^{JET}");
  fNChTrRecECorr->GetZaxis()->SetTitle("Selection Bin");
  fNChTrRecECorr->Sumw2();
  fOutputList->Add(fNChTrRecECorr);

  fNChTrRecPerpECorr = new TH3F("fNChTrRecPerpECorr","Tracks above min in perp.cone , corr. jet pT, centrality",101,-0.5,100.5, 60, 0., 300.,10,0.5,10.5);
  fNChTrRecPerpECorr->GetXaxis()->SetTitle("NTracks");
  fNChTrRecPerpECorr->GetYaxis()->SetTitle("p_{T}^{JET}");
  fNChTrRecPerpECorr->GetZaxis()->SetTitle("Selection Bin");
  fNChTrRecPerpECorr->Sumw2();
  fOutputList->Add(fNChTrRecPerpECorr);
  
  fRefMult = new TH1F("fRefMult", "Reference multiplicity in the AOD", 301, -0.5, 300.5);
  fRefMult->GetXaxis()->SetTitle("Reference multiplicity");
  fRefMult->Sumw2();
  fOutputList->Add(fRefMult);

  fRefMultWOJet = new TH2F("fRefMultWOJet", "Reference multiplicity in the AOD, multiplicity without jets", 301, -0.5, 300.5, 301, -0.5, 300.5);
  fRefMultWOJet->GetXaxis()->SetTitle("Reference multiplicity");
  fRefMultWOJet->GetYaxis()->SetTitle("Multiplicity without jets");
  fRefMultWOJet->Sumw2();
  fOutputList->Add(fRefMultWOJet);

  fVZEROMult = new TH2F("fVZEROMult", "Multiplicity V0A and V0C", 501, -0.5, 500.5, 501, -0.5, 500.5);
  fVZEROMult->GetXaxis()->SetTitle("Multiplicity V0A");
  fVZEROMult->GetYaxis()->SetTitle("Multiplicity V0C");
  fVZEROMult->Sumw2();
  fOutputList->Add(fVZEROMult);

  fMultWOJetVZero = new TH2F("fMultWOJetVZero", "Multiplicity without jets and VZERO mult.",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fMultWOJetVZero->GetXaxis()->SetTitle("Multiplicity without jets TPC");
  fMultWOJetVZero->GetYaxis()->SetTitle("Multiplicity full V0");
  fMultWOJetVZero->Sumw2();
  fOutputList->Add(fMultWOJetVZero);

  fRefMultFullV0 = new TH2F("fRefMultFullV0", "Reference multiplicity in the AOD, multiplicity from full V0",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fRefMultFullV0->GetXaxis()->SetTitle("Reference multiplicity in AOD");
  fRefMultFullV0->GetYaxis()->SetTitle("Multiplicity full V0");
  fRefMultFullV0->Sumw2();
  fOutputList->Add(fRefMultFullV0);

  fRefMultV0Corr = new TH2F("fRefMultV0Corr", "Reference multiplicity in the AOD, multiplicity from corrected V0",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fRefMultV0Corr->GetXaxis()->SetTitle("Reference multiplicity in AOD");
  fRefMultV0Corr->GetYaxis()->SetTitle("Multiplicity V0 no jets");
  fRefMultV0Corr->Sumw2();
  fOutputList->Add(fRefMultV0Corr);

  fFullV0V0Corr = new TH2F("fFullV0V0Corr", "Multiplicity from full V0, multiplicity from corrected V0",1001, -0.5, 1000.5, 1001, -0.5, 1000.5);
  fFullV0V0Corr->GetXaxis()->SetTitle("Multiplicity from full V0");
  fFullV0V0Corr->GetYaxis()->SetTitle("Multiplicity V0 no jets");
  fFullV0V0Corr->Sumw2();
  fOutputList->Add(fFullV0V0Corr);

  fNTXV0MultPt = new TH3F("fNTXV0MultPt", "NTX, Multiplicity from corrected V0, jet pT",101,-0.5,100.5, 1001, -0.5, 1000.5, 60, 0., 300.);
  fNTXV0MultPt->GetXaxis()->SetTitle("NTracks");
  fNTXV0MultPt->GetYaxis()->SetTitle("Multiplicity V0 no jets");
  fNTXV0MultPt->GetZaxis()->SetTitle("p_{T}^{JET}");
  fNTXV0MultPt->Sumw2();
  fOutputList->Add(fNTXV0MultPt);

  fNTXCBMultPt = new TH3F("fNTXCBMultPt", "NTX, Multiplicity from corrected Central Barrel, jet pT",101,-0.5,100.5, 301, -0.5, 300.5, 60, 0., 300.);
  fNTXCBMultPt->GetXaxis()->SetTitle("NTracks");
  fNTXCBMultPt->GetYaxis()->SetTitle("Multiplicity corrected Central Barrel");
  fNTXCBMultPt->GetZaxis()->SetTitle("p_{T}^{JET}");
  fNTXCBMultPt->Sumw2();
  fOutputList->Add(fNTXCBMultPt);

  fRefMultFullV0UJ = new TH2F("fRefMultFullV0UJ", "Reference multiplicity in the AOD, multiplicity from full V0, 1 jet event",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fRefMultFullV0UJ->GetXaxis()->SetTitle("Reference multiplicity in AOD");
  fRefMultFullV0UJ->GetYaxis()->SetTitle("Multiplicity full V0");
  fRefMultFullV0UJ->Sumw2();
  fOutputList->Add(fRefMultFullV0UJ);

  fRefMultV0CorrUJ = new TH2F("fRefMultV0CorrUJ", "Reference multiplicity in the AOD, multiplicity from corrected V0, 1 jet event",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fRefMultV0CorrUJ->GetXaxis()->SetTitle("Reference multiplicity in AOD");
  fRefMultV0CorrUJ->GetYaxis()->SetTitle("Multiplicity V0 no jets");
  fRefMultV0CorrUJ->Sumw2();
  fOutputList->Add(fRefMultV0CorrUJ);

  fFullV0V0CorrUJ = new TH2F("fFullV0V0CorrUJ", "Multiplicity from full V0, multiplicity from corrected V0, 1 jet event",1001, -0.5, 1000.5, 1001, -0.5, 1000.5);
  fFullV0V0CorrUJ->GetXaxis()->SetTitle("Multiplicity from full V0");
  fFullV0V0CorrUJ->GetYaxis()->SetTitle("Multiplicity V0 no jets");
  fFullV0V0CorrUJ->Sumw2();
  fOutputList->Add(fFullV0V0CorrUJ);

  fMultWOJetVZeroUJ = new TH2F("fMultWOJetVZeroUJ", "Multiplicity without jets and VZERO mult., 1 jet event",301, -0.5, 300.5, 1001, -0.5, 1000.5);
  fMultWOJetVZeroUJ->GetXaxis()->SetTitle("Multiplicity without jets TPC");
  fMultWOJetVZeroUJ->GetYaxis()->SetTitle("Multiplicity full V0");
  fMultWOJetVZeroUJ->Sumw2();
  fOutputList->Add(fMultWOJetVZeroUJ);

  fRefMultWOJetUJ = new TH2F("fRefMultWOJetUJ", "Reference multiplicity in the AOD, multiplicity without jets, 1 jet event", 301, -0.5, 300.5, 301, -0.5, 300.5);
  fRefMultWOJetUJ->GetXaxis()->SetTitle("Reference multiplicity");
  fRefMultWOJetUJ->GetYaxis()->SetTitle("Multiplicity without jets");
  fRefMultWOJetUJ->Sumw2();
  fOutputList->Add(fRefMultWOJetUJ);

  fRefAODTrackCount = new TH2F("fRefAODTrackCount", "Reference multiplicity in the AOD, my own referece mult.", 301, -0.5, 300.5, 301, -0.5, 300.5);
  fRefAODTrackCount->GetXaxis()->SetTitle("AOD Reference multiplicity");
  fRefAODTrackCount->GetYaxis()->SetTitle("My Reference multiplicity");
  fRefAODTrackCount->Sumw2();
  fOutputList->Add(fRefAODTrackCount);

  fRefAODTrackCountUJ = new TH2F("fRefAODTrackCountUJ", "Reference multiplicity in the AOD, my own referece mult., 1 jet event", 301, -0.5, 300.5, 301, -0.5, 300.5);
  fRefAODTrackCountUJ->GetXaxis()->SetTitle("AOD Reference multiplicity");
  fRefAODTrackCountUJ->GetYaxis()->SetTitle("My Reference multiplicity");
  fRefAODTrackCountUJ->Sumw2();
  fOutputList->Add(fRefAODTrackCountUJ);

  fTrackCountWOJet = new TH2F("fTrackCountWOJet", "My own total referece mult., soft mult", 151, -0.5, 150.5, 151, -0.5, 150.5);
  fTrackCountWOJet->GetXaxis()->SetTitle("Total TPC multiplicity");
  fTrackCountWOJet->GetYaxis()->SetTitle("Soft TPC multiplicity");
  fTrackCountWOJet->Sumw2();
  fOutputList->Add(fTrackCountWOJet);

  fTrackCountWOJetUJ = new TH2F("fTrackCountWOJetUJ", "My own total referece mult., soft mult, 1 jet", 151, -0.5, 150.5, 151, -0.5, 150.5);
  fTrackCountWOJetUJ->GetXaxis()->SetTitle("Total TPC multiplicity");
  fTrackCountWOJetUJ->GetYaxis()->SetTitle("Soft TPC multiplicity");
  fTrackCountWOJetUJ->Sumw2();
  fOutputList->Add(fTrackCountWOJetUJ);

  fTrackCountWOJetUJMC = new TH2F("fTrackCountWOJetUJMC", "My own total referece mult., soft mult, 1 jet, MC!", 151, -0.5, 150.5, 151, -0.5, 150.5);
  fTrackCountWOJetUJMC->GetXaxis()->SetTitle("Total TPC (eta) multiplicity");
  fTrackCountWOJetUJMC->GetYaxis()->SetTitle("Soft TPC (eta) multiplicity");
  fTrackCountWOJetUJMC->Sumw2();
  fOutputList->Add(fTrackCountWOJetUJMC);

  fFullV0V0CorrUJMC = new TH2F("fFullV0V0CorrUJMC", "Multiplicity from full V0, multiplicity from corrected V0, 1 jet event, MC!",1001, -0.5, 1000.5, 1001, -0.5, 1000.5);
  fFullV0V0CorrUJMC->GetXaxis()->SetTitle("Multiplicity from full V0 (acceptance)");
  fFullV0V0CorrUJMC->GetYaxis()->SetTitle("Multiplicity V0(acceptance) no jets");
  fFullV0V0CorrUJMC->Sumw2();
  fOutputList->Add(fFullV0V0CorrUJMC);

  fMinTrackPtInNTXRecalc = new TH3F("fMinTrackPtInNTXRecalc", "Minimum track pT for the jets after pT correction, raw jet pT", 200, 0., 100., 60, 0., 300.,10,0.5,10.5);
  fMinTrackPtInNTXRecalc->GetXaxis()->SetTitle("p_{T}^{TRACK} (GeV/c)");
  fMinTrackPtInNTXRecalc->GetYaxis()->SetTitle("p_{T}^{JET} (GeV/c)");
  fMinTrackPtInNTXRecalc->GetZaxis()->SetTitle("Selection Bin");  // 9 selections bins as fNChTrRecECorr
  fMinTrackPtInNTXRecalc->Sumw2();
  fOutputList->Add(fMinTrackPtInNTXRecalc);

  fMaxTrackPtInNTXRecalc = new TH2F("fMaxTrackPtInNTXRecalc", "Maximum track pT for the jets after pT correction, raw jet pT", 200, 0., 100., 60, 0., 300.);
  fMaxTrackPtInNTXRecalc->GetXaxis()->SetTitle("p_{T}^{TRACK} (GeV/c)");
  fMaxTrackPtInNTXRecalc->GetYaxis()->SetTitle("p_{T}^{JET} (GeV/c)");
  fMaxTrackPtInNTXRecalc->Sumw2();
  fOutputList->Add(fMaxTrackPtInNTXRecalc);

  fPtDistInJetConeRaw = new TH3F("fPtDistInJetConeRaw","pT of tracks in cone, raw jet pT bin, centrality", 200, 0., 100., 8, 0.5, 8.5, 10, 0.5, 10.5);
  fPtDistInJetConeRaw->GetXaxis()->SetTitle("p_{T}^{TRACK} (GeV/c)");
  fPtDistInJetConeRaw->GetYaxis()->SetTitle("p_{T}^{JET} Bin");
  fPtDistInJetConeRaw->GetZaxis()->SetTitle("Centrality Bin");
  fPtDistInJetConeRaw->Sumw2();
  fOutputList->Add(fPtDistInJetConeRaw);

  fPtDistInPerpConeRaw = new TH3F("fPtDistInPerpConeRaw","pT of tracks in cone, raw jet pT bin, centrality", 200, 0., 100., 8, 0.5, 8.5, 10, 0.5, 10.5);
  fPtDistInPerpConeRaw->GetXaxis()->SetTitle("p_{T}^{TRACK} (GeV/c)");
  fPtDistInPerpConeRaw->GetYaxis()->SetTitle("p_{T}^{JET} Bin");
  fPtDistInPerpConeRaw->GetZaxis()->SetTitle("Centrality Bin");
  fPtDistInPerpConeRaw->Sumw2();
  fOutputList->Add(fPtDistInPerpConeRaw);

  fPtInPerpCon = new TH3F("fPtInPerpCon","Summed pT of perpendicular cone, raw jet pT bin, centrality", 200, 0., 100., 8, 0.5, 8.5, 10, 0.5, 10.5);
  fPtInPerpCon->GetXaxis()->SetTitle("p_{T}^{PERP.CONE} (GeV/c)");
  fPtInPerpCon->GetYaxis()->SetTitle("p_{T}^{JET} Bin");
  fPtInPerpCon->GetZaxis()->SetTitle("Centrality Bin");
  fPtInPerpCon->Sumw2();
  fOutputList->Add(fPtInPerpCon);

  fJetEtaAll = new TH1F("fJetEtaAll", "Eta distribution of reconstructed jets, no cuts", 50, -1.5, 1.5);
  fJetEtaAll->GetXaxis()->SetTitle("#eta");
  fJetEtaAll->GetYaxis()->SetTitle("entries");
  fJetEtaAll->Sumw2();
  fOutputList->Add(fJetEtaAll);

  fJetEtaOnlyTPCcut = new TH1F("fJetEtaOnlyTPCcut", "Eta distribution of reconstructed jets, only tpc acceptance cut", 50, -1.5, 1.5);
  fJetEtaOnlyTPCcut->GetXaxis()->SetTitle("#eta");
  fJetEtaOnlyTPCcut->GetYaxis()->SetTitle("entries");
  fJetEtaOnlyTPCcut->Sumw2();
  fOutputList->Add(fJetEtaOnlyTPCcut);

  // 9 multiplicity bins
  // 1st.     <5    TPC tracks       fill 1  Bin1 [0.5,1.5)
  // 2nd. >= 5  <10 TPC tracks       fill 2  Bin2 [1.5,2.5)
  // 3rd. >= 10 <15 TPC tracks       fill 3  Bin3 [2.5,3.5)
  // 4rd. >= 15 <20 TPC tracks       fill 4  Bin4 [3.5,4.5)
  // 5th. >= 20 <30 TPC tracks       fill 5  Bin5 [4.5,5.5)
  // 6th. >= 30 <40 TPC tracks       fill 6  Bin6 [5.5,6.5)
  // 7th. >= 40 <50 TPC tracks       fill 7  Bin7 [6.5,7.5)
  // 8th.    >50    TPC tracks       fill 8  Bin8 [7.5,8.5)

  fNChTrRecECorrPPMult = new TH3F("fNChTrRecECorrPPMult","NTX in ener. corr. jet , corr. jet pT, pp mult.",101,-0.5,100.5, 60, 0., 300.,8,0.5,8.5);
  fNChTrRecECorrPPMult->GetXaxis()->SetTitle("NTracks_Corrected");
  fNChTrRecECorrPPMult->GetYaxis()->SetTitle("p_{T}^{JET}");
  fNChTrRecECorrPPMult->GetZaxis()->SetTitle("Multiplicity Bin");
  fNChTrRecECorrPPMult->Sumw2();
  fOutputList->Add(fNChTrRecECorrPPMult);

  fNChTrRecPerpECorrPPMult = new TH3F("fNChTrRecPerpECorrPPMult","Tracks above min in perp.cone , corr. jet pT, centrality",101,-0.5,100.5, 60, 0., 300.,8,0.5,8.5);
  fNChTrRecPerpECorrPPMult->GetXaxis()->SetTitle("NTracks_Corrected");
  fNChTrRecPerpECorrPPMult->GetYaxis()->SetTitle("p_{T}^{JET}");
  fNChTrRecPerpECorrPPMult->GetZaxis()->SetTitle("Multiplicity Bin");
  fNChTrRecPerpECorrPPMult->Sumw2();
  fOutputList->Add(fNChTrRecPerpECorrPPMult);

  fJetPtCentPbPbRaw = new TH2F("fJetPtCentPbPbRaw", "raw p_{T} distribution of reco jets", 60, 0., 300.,10,0.5,10.5);
  fJetPtCentPbPbRaw->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fJetPtCentPbPbRaw->GetYaxis()->SetTitle("Selection Bin");
  fJetPtCentPbPbRaw->Sumw2();
  fOutputList->Add(fJetPtCentPbPbRaw);

  fJetPtCentPbPbCorr = new TH2F("fJetPtCentPbPbCorr", "Corrected p_{T} distribution of reco jets", 60, 0., 300.,10,0.5,10.5);
  fJetPtCentPbPbCorr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fJetPtCentPbPbCorr->GetYaxis()->SetTitle("Selection Bin");
  fJetPtCentPbPbCorr->Sumw2();
  fOutputList->Add(fJetPtCentPbPbCorr);

  for(Int_t ipt=0;ipt<12;ipt++)
    {
      fNChTr[ipt] = new TH2F(Form("fNChTr[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}",101,-0.5,100.5, 60, 0., 300.);
      fNChTr[ipt]->GetXaxis()->SetTitle("NTracks");
      fNChTr[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
      fNChTr[ipt]->Sumw2();
      fOutputList->Add(fNChTr[ipt]);

      fHistPtParton[ipt] = new TH1F(Form("fHistPtParton[%i]",ipt),"pT distribution of jets",50,0.,250.);
      fHistPtParton[ipt]->GetXaxis()->SetTitle("p_{T}^{JET}");
      fHistPtParton[ipt]->GetYaxis()->SetTitle("Entries");
      fHistPtParton[ipt]->Sumw2();
      fOutputList->Add(fHistPtParton[ipt]);

      fSCM[ipt] = new TH2F(Form("fSCM[%i]",ipt),"Second Central Moment, jet_{p_{T}}",200,0.,0.2, 60, 0., 300.);
      fSCM[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
      fSCM[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
      fSCM[ipt]->Sumw2();
      fOutputList->Add(fSCM[ipt]);

      if(ipt<8) 
	{ 
	  fNChTrRDMult[ipt] = new TH2F(Form("fNChTrRDMult[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMult[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMult[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMult[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMult[ipt]);
	  
	  fNAccJetsRDMult[ipt] = new TH1F(Form("fNAccJetsRDMult[%i]",ipt),"Number of accepted jets per event in real data", 101, -0.5, 100.5);
	  fNAccJetsRDMult[ipt]->GetXaxis()->SetTitle("Number of jets");
	  fNAccJetsRDMult[ipt]->Sumw2();
	  fOutputList->Add(fNAccJetsRDMult[ipt]);
	  
	  fTotalJetCharge[ipt] = new TH1F(Form("fTotalJetCharge[%i]",ipt),"Charge in the jet", 41, -20.5, 20.5);
	  fTotalJetCharge[ipt]->GetXaxis()->SetTitle("Charge in jet");
	  fTotalJetCharge[ipt]->Sumw2();
	  fOutputList->Add(fTotalJetCharge[ipt]);
	  
	  fSCMRDMult[ipt] = new TH2F(Form("fSCMRDMult[%i]",ipt),"Second Central Moment, jet_{p_{T}}",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMult[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMult[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMult[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMult[ipt]);
	  
	  fNChTrRDMultMC[ipt] = new TH2F(Form("fNChTrRDMultMC[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultMC[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultMC[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultMC[ipt]);
	  
	  fSCMRDMultMC[ipt] = new TH2F(Form("fSCMRDMultMC[%i]",ipt),"Second Central Moment, jet_{p_{T}}",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultMC[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultMC[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultMC[ipt]);

	  //Second multiplicity estimator, removing jets and an area
	  fNChTrRDMultSE[ipt] = new TH2F(Form("fNChTrRDMultSE[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultSE[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultSE[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultSE[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultSE[ipt]);
	  
	  fNAccJetsRDMultSE[ipt] = new TH1F(Form("fNAccJetsRDMultSE[%i]",ipt),"Number of accepted jets per event in real data", 101, -0.5, 100.5);
	  fNAccJetsRDMultSE[ipt]->GetXaxis()->SetTitle("Number of jets");
	  fNAccJetsRDMultSE[ipt]->Sumw2();
	  fOutputList->Add(fNAccJetsRDMultSE[ipt]);
	  
	  fTotalJetChargeSE[ipt] = new TH1F(Form("fTotalJetChargeSE[%i]",ipt),"Charge in the jet", 41, -20.5, 20.5);
	  fTotalJetChargeSE[ipt]->GetXaxis()->SetTitle("Charge in jet");
	  fTotalJetChargeSE[ipt]->Sumw2();
	  fOutputList->Add(fTotalJetChargeSE[ipt]);
	  
	  fSCMRDMultSE[ipt] = new TH2F(Form("fSCMRDMultSE[%i]",ipt),"Second Central Moment, jet_{p_{T}}",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultSE[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultSE[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultSE[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultSE[ipt]);

          fNChTrRDMultOJ[ipt] = new TH2F(Form("fNChTrRDMultOJ[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}, 1 jet}",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultOJ[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultOJ[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultOJ[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultOJ[ipt]);

	  fSCMRDMultOJ[ipt] = new TH2F(Form("fSCMRDMultOJ[%i]",ipt),"Second Central Moment, jet_{p_{T}}, 1 jet",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultOJ[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultOJ[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultOJ[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultOJ[ipt]);

	  fNChTrRDMultSEOJ[ipt] = new TH2F(Form("fNChTrRDMultSEOJ[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}, 1 jet",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultSEOJ[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultSEOJ[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultSEOJ[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultSEOJ[ipt]);

	  fSCMRDMultSEOJ[ipt] = new TH2F(Form("fSCMRDMultSEOJ[%i]",ipt),"Second Central Moment, jet_{p_{T}}, 1 jet",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultSEOJ[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultSEOJ[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultSEOJ[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultSEOJ[ipt]);

          fNChTrRDMultOJMC[ipt] = new TH2F(Form("fNChTrRDMultOJMC[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}, 1 jet, MC}",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultOJMC[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultOJMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultOJMC[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultOJMC[ipt]);

	  fSCMRDMultOJMC[ipt] = new TH2F(Form("fSCMRDMultOJMC[%i]",ipt),"Second Central Moment, jet_{p_{T}}, 1 jet, MC",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultOJMC[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultOJMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultOJMC[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultOJMC[ipt]);

	  fNChTrRDMultSEOJMC[ipt] = new TH2F(Form("fNChTrRDMultSEOJMC[%i]",ipt),"Number of tracks to recover transverse energy, jet_{p_{T}}, 1 jet, MC",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRDMultSEOJMC[ipt]->GetXaxis()->SetTitle("NTracks");
	  fNChTrRDMultSEOJMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRDMultSEOJMC[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRDMultSEOJMC[ipt]);

	  fSCMRDMultSEOJMC[ipt] = new TH2F(Form("fSCMRDMultSEOJMC[%i]",ipt),"Second Central Moment, jet_{p_{T}}, 1 jet, MC",200,0.,0.2, 60, 0., 300.);
	  fSCMRDMultSEOJMC[ipt]->GetXaxis()->SetTitle("<#delta R_{c}^{2}>");
	  fSCMRDMultSEOJMC[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fSCMRDMultSEOJMC[ipt]->Sumw2();
	  fOutputList->Add(fSCMRDMultSEOJMC[ipt]);

	  fNChTrRecPerpMultSEOJ[ipt] = new TH2F(Form("fNChTrRecPerpMultSEOJ[%i]",ipt),"Number of tracks above the min pT used in NTX_Raw, jet_{p_{T}}, 1 jet",101,-0.5,100.5, 60, 0., 300.);
	  fNChTrRecPerpMultSEOJ[ipt]->GetXaxis()->SetTitle("NTracks_{Exc.}");
	  fNChTrRecPerpMultSEOJ[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
	  fNChTrRecPerpMultSEOJ[ipt]->Sumw2();
	  fOutputList->Add(fNChTrRecPerpMultSEOJ[ipt]);
	} // end if <8

      if(ipt<6)  // only entries for reconstructed || simulated jets
	{ 
	  fProcessPDG[ipt] = new TH2I(Form("fProcessPDG[%i]",ipt),"Pythia process and jet p_T", 60, 0., 300., 70, 0.5, 70.5);
	  fProcessPDG[ipt]->GetXaxis()->SetTitle("p_{T}^{JET}");
	  fProcessPDG[ipt]->GetYaxis()->SetTitle("Pythia process code");
	  fProcessPDG[ipt]->Sumw2();
	  fOutputList->Add(fProcessPDG[ipt]);

	  fFragPion[ipt] = new TH2F(Form("fFragPion[%i]",ipt),"Fragmentation in pions",35,0.,8.,50,0.,250.);
	  fFragPion[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragPion[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragPion[ipt]->Sumw2();
	  fOutputList->Add(fFragPion[ipt]);
	  
	  fFragKaon[ipt] = new TH2F(Form("fFragKaon[%i]",ipt),"Fragmentation in kaons",35,0.,8.,50,0.,250.);
	  fFragKaon[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragKaon[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragKaon[ipt]->Sumw2();
	  fOutputList->Add(fFragKaon[ipt]);

	  fFragProton[ipt] = new TH2F(Form("fFragProton[%i]",ipt),"Fragmentation in protons",35,0.,8.,50,0.,250.);
	  fFragProton[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragProton[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragProton[ipt]->Sumw2();
	  fOutputList->Add(fFragProton[ipt]);

	  fFragChargedR4[ipt] = new TH2F(Form("fFragChargedR4[%i]",ipt),"Fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fFragChargedR4[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragChargedR4[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragChargedR4[ipt]->Sumw2();
	  fOutputList->Add(fFragChargedR4[ipt]);
	  
	  fFragChargedR3[ipt] = new TH2F(Form("fFragChargedR3[%i]",ipt),"Fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fFragChargedR3[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragChargedR3[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragChargedR3[ipt]->Sumw2();
	  fOutputList->Add(fFragChargedR3[ipt]);

	  fFragChargedR2[ipt] = new TH2F(Form("fFragChargedR2[%i]",ipt),"Fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fFragChargedR2[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fFragChargedR2[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fFragChargedR2[ipt]->Sumw2();
	  fOutputList->Add(fFragChargedR2[ipt]);

	  // do not add the temporary containers
	  fHistContainerR4[ipt] = new TH2F(Form("fHistContainerR4[%i]",ipt),"Temporary fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fHistContainerR4[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fHistContainerR4[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fHistContainerR4[ipt]->Sumw2();
	  
	  fHistContainerR3[ipt] = new TH2F(Form("fHistContainerR3[%i]",ipt),"Temporary fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fHistContainerR3[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fHistContainerR3[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fHistContainerR3[ipt]->Sumw2();
	  
	  fHistContainerR2[ipt] = new TH2F(Form("fHistContainerR2[%i]",ipt),"Temporary fragmentation in charged particles",35,0.,8.,50,0.,250.);
	  fHistContainerR2[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
	  fHistContainerR2[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
	  fHistContainerR2[ipt]->Sumw2();

	  if(ipt<3)
	    {
	      fJetEtaJetPt[ipt] = new TH1F(Form("fJetEtaJetPt[%i]",ipt), "Eta distribution of reconstructed jets, all cut, with pT upper boundary", 50, -1.5, 1.5);
	      fJetEtaJetPt[ipt]->GetXaxis()->SetTitle("#eta");
	      fJetEtaJetPt[ipt]->GetYaxis()->SetTitle("entries");
	      fJetEtaJetPt[ipt]->Sumw2();
	      fOutputList->Add(fJetEtaJetPt[ipt]);

	      if(ipt<2)
		{
		  fFragCandidates[ipt] = new TH2F(Form("fFragCandidates[%i]",ipt),"Parton identified candidates",35,0.,8.,50,0.,250.);
		  fFragCandidates[ipt]->GetXaxis()->SetTitle("#xi=ln[Jet_{E_{T}}/H_{p_{T}}]");
		  fFragCandidates[ipt]->GetYaxis()->SetTitle("Jet_{E_{T}}");
		  fFragCandidates[ipt]->Sumw2();
		  fOutputList->Add(fFragCandidates[ipt]);
		  
		  fMinTrackPtInNTXh[ipt] = new TH3F(Form("fMinTrackPtInNTXh[%i]",ipt), "Minimum track pT for the jets", 200, 0., 100., 60, 0., 300.,10,0.5,10.5);
		  fMinTrackPtInNTXh[ipt]->GetXaxis()->SetTitle("p_{T}^{TRACK}");
		  fMinTrackPtInNTXh[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
		  fMinTrackPtInNTXh[ipt]->GetZaxis()->SetTitle("Selection Bin"); //9 selection bins as fNChTrRecECorr
		  fMinTrackPtInNTXh[ipt]->Sumw2();
		  fOutputList->Add(fMinTrackPtInNTXh[ipt]);
		  
		  fMaxTrackPtInNTXh[ipt] = new TH2F(Form("fMaxTrackPtInNTXh[%i]",ipt), "Maximum track pT for the jets", 200, 0., 100., 60, 0., 300.);
		  fMaxTrackPtInNTXh[ipt]->GetXaxis()->SetTitle("p_{T}^{TRACK}");
		  fMaxTrackPtInNTXh[ipt]->GetYaxis()->SetTitle("p_{T}^{JET}");
		  fMaxTrackPtInNTXh[ipt]->Sumw2();
		  fOutputList->Add(fMaxTrackPtInNTXh[ipt]);  
		} // index < 2 
	    } // index < 3 
	} // index < 6 
    } // index < 12 

  fPerpCone = new AliAODJet();
  fBckgSbsJet = new Double_t[3];

  PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPartonDisc::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  //  fAOD
  if(fUseAODJetInput)
    {    
      fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
      if(!fAOD)
	{
	  Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODJetInput);
	  return;
	}
      // fetch the header
    }
  else
    {
      //  assume that the AOD is in the general output...
      fAOD  = AODEvent();
      if(!fAOD)
	{
	  Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
	  return;
	}
    }

  // fin de test para fAOD

  if(!fInputEvent)
    {
      Error("UserExec","No event found!");
      return;
    }
   
  AliAODHandler *aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!aodHandler)
    {
      AliError("No AOD Handler!");
      return;
    }

  fEventCent=900; //dummy val for debugging

  if(fIsHIevent)
    {
      AliAODHeader *aodHeader = fAOD->GetHeader();
      fEventCent = aodHeader->GetCentrality();
    }

  // Jet eta exclusion
  if(fIncreasingExcl)
    fJetAcceptance = 0.5 - fIncExcR; // if the increase is 0.1 -> only jets within |eta|<0.4 

  // First test of reference multiplicity
  Int_t refMultiplicity = fAOD->GetHeader()->GetRefMultiplicity();
  fRefMult->Fill(refMultiplicity);

  // Multiplicity from V0 (V0A+V0C)
  fVZero = fAOD->GetVZEROData();
  Float_t multV0A = 0.0;
  Float_t multV0C = 0.0;
  Float_t multFullV0 = 0.0;
  if(fVZero)
    {
      multV0A = fVZero->GetMTotV0A();
      multV0C = fVZero->GetMTotV0C();
      multFullV0 = multV0A+multV0C; 
    }
  fVZEROMult->Fill(multV0A,multV0C);

  fEvtCount++;
  Double_t jfr = fJetRadius;   // radius used during jet finding
  Int_t ntx = fNtX;   // NTX value
  const Int_t maxJetNum=6; // maximum number of generated jets to process
  AliAODJet genJets[6];  // containers for the
  AliAODJet recJets[6];  // correlation gen-reco
  Int_t nGenJets=0; 
  Int_t nRecJets=0;
  Int_t genJetsFlavor[6]={0};    // flavor of the generated jets
  Int_t evtype = 0; //pythia event type
  // Variables para la variable de estructura
  Double_t deltaPhiPt = 0.0;
  Double_t deltaEtaPt = 0.0;
  Double_t deltaPhiSqPt = 0.0;
  Double_t deltaEtaSqPt = 0.0;
  Double_t totalTrackPt = 0.0; 
  Double_t firstMomDeltPhi = 0.0;
  Double_t firstMomDeltEta = 0.0;
  Double_t secondMomDeltPhi = 0.0;
  Double_t secondMomDeltEta = 0.0;
  Double_t secondCentralPhi = 0.0;
  Double_t secondCentralEta = 0.0;
  Double_t secondCentralR = 0.0; 

  // Variables para la variable de estructura
  // del cono perpendicular
  Double_t deltaPhiPtPerp = 0.0;
  Double_t deltaEtaPtPerp = 0.0;
  Double_t deltaPhiSqPtPerp = 0.0;
  Double_t deltaEtaSqPtPerp = 0.0;
  Double_t totalTrackPtPerp = 0.0;
  Double_t firstMomDeltPhiPerp = 0.0;
  Double_t firstMomDeltEtaPerp = 0.0;
  Double_t secondMomDeltPhiPerp = 0.0;
  Double_t secondMomDeltEtaPerp = 0.0;
  Double_t secondCentralPhiPerp = 0.0;
  Double_t secondCentralEtaPerp = 0.0;
  Double_t secondCentralRPerp = 0.0;  

  Double_t perpendicularPt;
  Float_t px,py,pz,en; // jet 4-vector à la UA1
  Float_t pTbs, etabs, phibs; // energy corrected jet properties

  // Process the MC info from the AOD
  if(fUseAODMC)
    {
      // Get the MC array
      TClonesArray *mcarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!mcarray)
	{
	  AliError("ERROR:No MC info in the AOD input");
	  return;
	} 

      AliMCEvent* mcEvent = MCEvent();
      if(mcEvent)
	{
	  if(!fPhojetMC) // if it is pythia
	    evtype = GetMCEventType(mcEvent);
	  // From pythia 6.2 manual pp 414
	  // QCD Hard Processes
	  // 11 f_{i}+f_{j} -> f_{i}+f_{j} com77, ben84, eic84, chi90 
	  // 12 f_{i}+barf_{i} -> f_{k}+barf_{k}
	  // 13 f_{i}+barf_{i} -> g+g
	  // 28 f_{i}+g -> f_{i}+g
	  // 53 g+g -> f_{k}+barf_{k}
	  // 68 g+g -> g+g
	  if(fPhojetMC) // if it is phojet
	    evtype = GetPhojetEventType(mcEvent); 
	}
      if(!mcEvent) // if a pure AOD event
	{
	  AliDebug(2,Form("%s:%d No MCEvent \n",(char*)__FILE__,__LINE__));  
	  AliDebug(2,Form("Trying to get the MC header \n"));  
	  AliAODMCHeader *genEvH = static_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));
	  if(!genEvH)
	    {
	      AliDebug(2,Form(" %s:%d No Pythia header!",(char*)__FILE__,__LINE__));  
	      evtype = 0;
	    }
	  if(genEvH)
	    evtype = genEvH->GetEventType();
	}   
      // Get the branch with the MC jets
      TClonesArray *aodMCJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchMC.Data()));
      if(!aodMCJets)
	{
	  AliDebug(2,Form("%s:%d no MC Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchMC.Data()));  
	  return;
	}
      AliDebug(2,Form("There are %d MC jets in this event\n", aodMCJets->GetEntries())); 
      Int_t mcjets =  aodMCJets->GetEntries();
      fNJetsMC->Fill(mcjets,mcjets); // number of jets FROM the branch, filled per event, this gives the event normalization...     
      HasOverlapedCones(aodMCJets); // Procedure for tagging usable jets
                                    // Up 16 jets are flagged

      // Loop over AODMC physical primary charged particles
      // for the complete event
      Int_t tracksMC = mcarray->GetEntriesFast();	  
      AliDebug(2,Form("There are %i tracks in the mcarray",tracksMC));  
      Double_t aodMCTrackEta = 0.0; 
      perpendicularPt = 0.0;
      px=0.0;
      py=0.0;
      pz=0.0;
      en=0.0;
      pTbs=0.0;
      etabs=0.0;
      phibs=0.0;
      fBckgSbsJet[0]=0.0;
      fBckgSbsJet[1]=0.0;
      fBckgSbsJet[2]=0.0;
      Int_t softRefMcNoJets = 0;
      Int_t myTotalMultiplicityMc = 0;
      Int_t v0LikeTotalMcMult = 0;
      for(Int_t aodMCTrack = 0; aodMCTrack < tracksMC; aodMCTrack++ )
	{
	  AliAODMCParticle *mctrackf = (AliAODMCParticle*) mcarray->At(aodMCTrack);
	  if(!mctrackf) continue;
	  if(!mctrackf->IsPhysicalPrimary()) continue;
	  if(mctrackf->Charge()==0||mctrackf->Charge()==-99) continue;
	  //Lo del V0, voy a contar particulas primarias cargadas
	  if(mctrackf->Pt()>fMinpTValMC) // cut off en MC
	    {
	      //V0A
	      if(((mctrackf->Eta())>(2.8))&&((mctrackf->Eta())<(5.1)))
		v0LikeTotalMcMult++;
	      //V0C
	      if(((mctrackf->Eta())>(-3.7))&&((mctrackf->Eta())<(-1.7)))
		v0LikeTotalMcMult++;
	    }
	  //Fin de lo del V0
	  aodMCTrackEta = TMath::Abs(mctrackf->Eta());
	  if(aodMCTrackEta>0.9) continue;
	  fPtAODMC->Fill(mctrackf->Pt(),mctrackf->Pt());
	  fEtaAODMC->Fill(mctrackf->Eta(),mctrackf->Eta());
	  fPhiAODMC->Fill(mctrackf->Phi(),mctrackf->Phi());
	  if(fJetEvent) // if has an accepted jet, calculate the perpendicular cone
	    {
	      if(HasPerpendicularCone()) // If there is a perpendicular cone available
		{
		  if(mctrackf->Pt()>fMinpTVal)
		    {
		    if(GetDeltaR(fEtaPerpCoord,fPhiPerpCoord,mctrackf->Eta(),mctrackf->Phi())<fJetRadius)
		      perpendicularPt = perpendicularPt + mctrackf->Pt();
		    }
		}
	    } // end IF jet event
	  if(mctrackf->Pt()>fMinpTValMC) // cut off en MC
	    {
	      myTotalMultiplicityMc++; // total multiplicity TPC like
	      if(mctrackf->Pt()<fMinpTValUE) continue; // pT cut  fMinpTValUE
	      if(mctrackf->Pt()>fMaxpTValUE) continue; // pT cut  fMaxpTValUE
	      if(!IsTrackInsideExcludedArea(mctrackf->Eta(), mctrackf->Phi(), aodMCJets)) 
		softRefMcNoJets++;
	    }
	} // end loop over particles

      Int_t correctedV0LikeMult= v0LikeTotalMcMult-GetV0LikeExcludedMultMC(aodMCJets,mcarray);

      //estimadores
      if(mcjets==1) // correlation for only monojet events
	{
	  fFullV0V0CorrUJMC->Fill(v0LikeTotalMcMult,correctedV0LikeMult);
	  fTrackCountWOJetUJMC->Fill(myTotalMultiplicityMc,softRefMcNoJets); 
	}

      if(fJetEvent) // if has an accepted jet, calculate the perpendicular cone
	{
	  if(HasPerpendicularCone()) // If there is a perpendicular cone available
	    {
	      px = perpendicularPt*TMath::Cos(fPhiPerpCoord);
	      py = perpendicularPt*TMath::Sin(fPhiPerpCoord);
	      pz = perpendicularPt/TMath::Tan(2.0*TMath::ATan(TMath::Exp(-fEtaPerpCoord)));
	      en = TMath::Sqrt(px*px + py*py + pz*pz);
	      fPerpCone->SetPxPyPzE(px, py, pz, en);
	    }
	  if(!HasPerpendicularCone())
	    AliDebug(2,"No perpendicular cone!!!");  
	}


      fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
 
      Int_t flavor = 0;     // flavor of the jet      
      Int_t nTracksPerc;    // ntx for the original jet
      Int_t nTracksPercBckgSubst;  // ntx for the energy corrected jet
      Double_t jetPt=0;
      Int_t pdgOfMCt;
      Float_t trackxi;
      Double_t jetXt;
      Double_t jetPts[7]={0};  // to store the pt of the jets
      Int_t mcJetCounter=0;    // counter of MC jets
      Int_t nTracksAboveThresholdPerp=0;  // n tracks of the perpendicular cone
      Int_t nTrUpThrPerpBckSubs=0;  // n tracks of the perpendicular cone, after the minimum pT recalculation
      fIsPossibleToSubstBckg = kTRUE; // Initialize before the loop
      if(fJetEvent) // si tiene jets validos
	{
	  if(!HasPerpendicularCone()) // pero no encontro un cono perpendicular libre
	    fIsPossibleToSubstBckg = kFALSE; // if not perpendicular cone, set to kFALSE, so no perpendicular calculations available
	}
      // Loop to fill a pT spectrum of the mc jets
      Int_t imcj=0; // index for montecarlo jets to correlate
      for (Int_t indxmc = 0; indxmc < mcjets; indxmc++) 
	{
	  AliAODJet *mcjet = dynamic_cast<AliAODJet*>(aodMCJets->At(indxmc));
	  if (!mcjet) 
	    {
	      AliDebug(2,Form("ERROR: Could not receive jet %d\n", indxmc));  
	      continue;
	    }
	  
	  ///////////////////////////////////////////////////////////////////////////////
	  ///// Part for Chritians plot of inclusive and leading jets comp at 2.76 TeV //
	  if(!IsInsideAcceptance(mcjet))  // old condition
	    continue;
	  if(indxmc==0) // leading jet
	    fMCJetPtLeading->Fill(mcjet->Pt());
	  fMCJetPtInclusive->Fill(mcjet->Pt()); // all
	  ///// End of Christians Plot MC
	  ///////////////////////////////////////////////////////////////////////////////
	  
	  if(indxmc>15)
	    continue;

	  if(!fJetFlags[indxmc]) // If the jet is flaged kFALSE, not usable
	    continue;

	  //Initialize variables for this jet
	  //adiciones para la variable de estructura
	  nTracksPerc = 0;
	  nTracksPercBckgSubst = 0;
	  fMinTrackPtInNTX=200.0;  //Initialize for each jet, overflown
	  fMaxTrackPtInNTX=200.0;  //Initialize for each jet, overflown
	  fMinTrackPtInNTXR=200.0;  //Initialize for each jet, overflown
	  fMaxTrackPtInNTXR=200.0;  //Initialize for each jet, overflown
	  deltaPhiPt = 0.0;
	  deltaEtaPt = 0.0;
	  deltaPhiSqPt = 0.0;
	  deltaEtaSqPt = 0.0;
	  totalTrackPt = 0.0;
	  firstMomDeltPhi = 0.0;
	  firstMomDeltEta = 0.0;
	  secondMomDeltPhi = 0.0;
	  secondMomDeltEta = 0.0;
	  secondCentralPhi = 0.0;
	  secondCentralEta = 0.0;
	  secondCentralR = 0.0; 
	  
	  if(imcj<maxJetNum)
	    genJets[imcj]= *mcjet;
	  if(mcJetCounter<maxJetNum)
	    jetPts[mcJetCounter]=mcjet->Pt();
	  mcJetCounter++;  // number of jets in the acceptance
	  jetPt = mcjet->Pt();
	  flavor =  GetJetFlavour(mcjet,tracksMC,mcarray);
	  if(imcj<maxJetNum)
	    genJetsFlavor[imcj] = flavor;
	  fJetPtMC->Fill(mcjet->Pt());
	  fJetEtaMC->Fill(mcjet->Eta(),mcjet->Eta());
	  fJetPhiMC->Fill(mcjet->Phi(),mcjet->Phi());
	  fFlavor->Fill(flavor,jetPt);
	  AliDebug(4,Form("Sabor del jet con pt=%f es :%d \n",jetPt,flavor)); 
	  nTracksPerc = GetNumberOfMcChargedTracks(ntx,mcjet,tracksMC,mcarray,jfr); // este fija el min track pT, si es posible substraer el bckg
	  if(fIsPossibleToSubstBckg&&!IsEqualRel(fCurrentJetMinPtNT90, 7000.)) //calculating only if there is a perpendicular cone available //IsEqualRel(jetpT, 0.0) //fCurrentJetMinPtNT90!=7000.
	    {                                                     //and only if the method worked
	      AliDebug(4,Form("For this jet and I have a perpendicular cone available")); 
	      // Aqui lo que debo contar es el numero de tracks arriba del min pT del jet correspondiente
	      // que es fCurrentJetMinPtNT90
	      nTracksAboveThresholdPerp = GetNMcChargedTracksAboveThreshold(fPerpCone,tracksMC,mcarray,jfr);
	    }
	  // Corrected jet (pT)
	  if(fIsPossibleToSubstBckg)   // for the current jet 
	    {                                
	      pTbs= mcjet->Pt()-fPerpCone->Pt();
	      etabs= mcjet->Eta();
	      phibs= mcjet->Phi();
	      fBckgSbsJet[0]=pTbs; //pT
	      fBckgSbsJet[1]=etabs; //eta
	      fBckgSbsJet[2]=phibs; //phi
	      // Now re-calculate nt90 for the energy corrected jet
	      nTracksPercBckgSubst = GetRecalcNTXMc(ntx,mcjet,tracksMC,mcarray,jfr);
	      // Now re-calculate the perpendicular cone NT90 background
	      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.)) //calculating for the current jet, if the previos method worked //!IsEqualRel(fCurrentJetMinPtNT90, 7000.) //fCurrentJetMinPtNT90Recalc!=7000.
		{
		  // Aqui lo que debo contar es el numero de tracks arriba del min pT del jet correspondiente
		  // despues de la correccion de energia del jet
		  nTrUpThrPerpBckSubs = GetRecalcNMcChTrUpThr(fPerpCone,tracksMC,mcarray,jfr);
		}
	    }

	  //check cross sections incoming partons
	  jetXt= 2*jetPt/fSqrts;
	  if(evtype==11||evtype==12||evtype==13) //QQ
	    fFracQQ->Fill(jetXt);
	  if(evtype==28) //GQ
	    fFracGQ->Fill(jetXt);
	  if(evtype==53||evtype==68) //GG
	    fFracGG->Fill(jetXt);

	  //check cross sections outgoing partons
	  if(evtype==11||evtype==12||evtype==53) //QQ
	    fFracOutGoingQQ->Fill(jetXt);
	  if(evtype==28) //GQ
	    fFracOutGoingGQ->Fill(jetXt);
	  if(evtype==13||evtype==68) //GG
	    fFracOutGoingGG->Fill(jetXt);

	  fProcessJetPt->Fill(evtype,jetPt);  // pythia process, filled for each jet in acceptance

	  //Fill jet flavor as a function of pT and the pythia process but only leading jet
	  if(imcj==0)  //leading jet
	    {
	      fFlavorLead->Fill(flavor,jetPt);
	      fProcessLeadJetPt->Fill(evtype,jetPt); 
	    }
	  AliDebug(4,Form("Before the check of comparison")); 
	  // To check tracks related to this MC jet
	  // RefTracks check
	  Bool_t rTrkFlagMC = kFALSE;
	  Int_t trkinmcjet = mcjet->GetRefTracks()->GetEntriesFast();
	  if(trkinmcjet!=0&&!fForceNotTR)
	    rTrkFlagMC = kTRUE;
	  AliDebug(4,Form("Number of tracks in RefTracks MC jet:%i \n",trkinmcjet));  
	  if(rTrkFlagMC)  // If there are tracks refs available
	    {
	      AliDebug(4,Form("Checking composition in MC with track refs")); 
	      for(Int_t aodMCT = 0; aodMCT < trkinmcjet; aodMCT++ )
		{
		  pdgOfMCt=0;
		  trackxi=0;
		  AliAODMCParticle *mctrack = (AliAODMCParticle*) mcjet->GetRefTracks()->At(aodMCT);
		  if(!mctrack) continue;
		  if(!mctrack->IsPhysicalPrimary()) continue;
		  if(mctrack->Charge()==0||mctrack->Charge()==-99) continue;
		  if(mctrack->Pt()<fMinpTVal) continue; // MC no cut in the case of track reference, should be in, NO, cut anyhow to be safe
		  deltaPhiPt += DeltaPhiMC(mcjet, mctrack)*mctrack->Pt();
		  deltaEtaPt += DeltaEtaMC(mcjet, mctrack)*mctrack->Pt();
		  deltaPhiSqPt += DeltaPhiSqMC(mcjet, mctrack)*mctrack->Pt();
		  deltaEtaSqPt += DeltaEtaSqMC(mcjet, mctrack)*mctrack->Pt();
		  totalTrackPt += mctrack->Pt();

		  pdgOfMCt=abs(mctrack->GetPdgCode());
		  if(!IsEqualRel(mctrack->Pt(), 0.0)) //!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.) // mctrack->Pt()!=0
		    trackxi= log(jetPt/mctrack->Pt());
		  switch(abs(flavor))
		    {
		    case 1:
		      if(pdgOfMCt==321)
			fFragKaon[0]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[0]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[0]->Fill(trackxi,jetPt);
		      break;
		    case 2:
		      if(pdgOfMCt==321)
			fFragKaon[1]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[1]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[1]->Fill(trackxi,jetPt);
		      break;
		    case 3:
		      if(pdgOfMCt==321)
			fFragKaon[2]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[2]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[2]->Fill(trackxi,jetPt);
		      break;
		    case 4:
		      if(pdgOfMCt==321)
			fFragKaon[3]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[3]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[3]->Fill(trackxi,jetPt);
		      break;
		    case 5:
		      if(pdgOfMCt==321)
			fFragKaon[4]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[4]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[4]->Fill(trackxi,jetPt);
		      break;
		    case 21:
		      if(pdgOfMCt==321)
			fFragKaon[5]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[5]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[5]->Fill(trackxi,jetPt);
		      break;	      
		    default:
		      break;
		    } // end switch flavor
		}// end loop over AODMC particles of trackrefs
	      if(!IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) //totalTrackPt!=0.0
		{
		  firstMomDeltPhi = deltaPhiPt/totalTrackPt;
		  firstMomDeltEta = deltaEtaPt/totalTrackPt;
		  secondMomDeltPhi = deltaPhiSqPt/totalTrackPt;
		  secondMomDeltEta = deltaEtaSqPt/totalTrackPt;
		  secondCentralPhi = secondMomDeltPhi - firstMomDeltPhi*firstMomDeltPhi;
		  secondCentralEta = secondMomDeltEta - firstMomDeltEta*firstMomDeltEta;
		  secondCentralR = secondCentralPhi + secondCentralEta;
		} // end if totalTrackPt!=0.0
	      if(IsEqualRel(totalTrackPt, 0.0))  //totalTrackPt==0.0
		secondCentralR = 10.0; //overflow
	    }// end version with ref tracks (flag check)
	    
	  if(!rTrkFlagMC)  // No ref tracks available
	    {
	      AliDebug(4,Form("Checking composition in MC without track refs")); 
	      for(Int_t aodMCT = 0; aodMCT < tracksMC; aodMCT++ )
		{
		  pdgOfMCt=0;
		  trackxi=0;
		  AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(aodMCT);
		  if(!mctrack) continue;
		  if(!mctrack->IsPhysicalPrimary()) continue;
		  if(mctrack->Charge()==0||mctrack->Charge()==-99) continue;
		  if(!IsMCTrackInsideThisJet(mctrack, mcjet, jfr)) continue;
		  if(mctrack->Pt()<fMinpTVal) continue; // MC: HERE PT CUT, NO TRACK REF
		  deltaPhiPt += DeltaPhiMC(mcjet, mctrack)*mctrack->Pt();
		  deltaEtaPt += DeltaEtaMC(mcjet, mctrack)*mctrack->Pt();
		  deltaPhiSqPt += DeltaPhiSqMC(mcjet, mctrack)*mctrack->Pt();
		  deltaEtaSqPt += DeltaEtaSqMC(mcjet, mctrack)*mctrack->Pt();
		  totalTrackPt += mctrack->Pt();

		  pdgOfMCt=abs(mctrack->GetPdgCode());
		  if(!IsEqualRel(mctrack->Pt(), 0.0)) //!IsEqualRel(mctrack->Pt(), 0.0) // mctrack->Pt()!=0
		    trackxi= log(jetPt/mctrack->Pt());
		  switch(flavor)
		    {
		    case 1:
		      if(pdgOfMCt==321)
			fFragKaon[0]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[0]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[0]->Fill(trackxi,jetPt);
		      break;
		    case 2:
		      if(pdgOfMCt==321)
			fFragKaon[1]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[1]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[1]->Fill(trackxi,jetPt);
		      break;
		    case 3:
		      if(pdgOfMCt==321)
			fFragKaon[2]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[2]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[2]->Fill(trackxi,jetPt);
		      break;
		    case 4:
		      if(pdgOfMCt==321)
			fFragKaon[3]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[3]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[3]->Fill(trackxi,jetPt);
		      break;
		    case 5:
		      if(pdgOfMCt==321)
			fFragKaon[4]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[4]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
			fFragProton[4]->Fill(trackxi,jetPt);
		      break;
		    case 21:
		      if(pdgOfMCt==321)
			fFragKaon[5]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==211)
			fFragPion[5]->Fill(trackxi,jetPt);
		      if(pdgOfMCt==2212)
		    fFragProton[5]->Fill(trackxi,jetPt);
		      break;	      
		    default:
		      break;
		    } // end switch flavor
		}// end loop over AODMC particles
	      if(!IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPt!=0.0
		{
		  firstMomDeltPhi = deltaPhiPt/totalTrackPt;
		  firstMomDeltEta = deltaEtaPt/totalTrackPt;
		  secondMomDeltPhi = deltaPhiSqPt/totalTrackPt;
		  secondMomDeltEta = deltaEtaSqPt/totalTrackPt;
		  secondCentralPhi = secondMomDeltPhi - firstMomDeltPhi*firstMomDeltPhi;
		  secondCentralEta = secondMomDeltEta - firstMomDeltEta*firstMomDeltEta;
		  secondCentralR = secondCentralPhi + secondCentralEta;
		} // end if totalTrackPt!=0.0
	      if(IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) //totalTrackPt==0.0
		secondCentralR = 10.0; //overflow
	    } //End old version (no ref tracks)

	  if(fIsPossibleToSubstBckg)
	    {
	      // To make sure, re-initialize
	      deltaPhiPtPerp = 0.0;
	      deltaEtaPtPerp = 0.0;
	      deltaPhiSqPtPerp = 0.0;
	      deltaEtaSqPtPerp = 0.0;
	      totalTrackPtPerp = 0.0;
	      firstMomDeltPhiPerp = 0.0;
	      firstMomDeltEtaPerp = 0.0;
	      secondMomDeltPhiPerp = 0.0;
	      secondMomDeltEtaPerp = 0.0;
	      secondCentralPhiPerp = 0.0;
	      secondCentralEtaPerp = 0.0;
	      secondCentralRPerp = 0.0;

	      AliDebug(4,Form("Checking SCM in MC for the perpendicular cone")); 
	      for(Int_t aodMCperp = 0; aodMCperp < tracksMC; aodMCperp++ )
		{
		  AliAODMCParticle *mctrackperp = (AliAODMCParticle*) mcarray->At(aodMCperp);
		  if(!mctrackperp) continue;
		  if(!mctrackperp->IsPhysicalPrimary()) continue;
		  if(mctrackperp->Charge()==0||mctrackperp->Charge()==-99) continue;
		  if(!IsMCTrackInsideThisJet(mctrackperp, fPerpCone, jfr)) continue;
		  if(mctrackperp->Pt()<fMinpTVal) continue; // MC: HERE PT CUT   
		  deltaPhiPtPerp += DeltaPhiMC(fPerpCone, mctrackperp)*mctrackperp->Pt();
		  deltaEtaPtPerp += DeltaEtaMC(fPerpCone, mctrackperp)*mctrackperp->Pt();
		  deltaPhiSqPtPerp += DeltaPhiSqMC(fPerpCone, mctrackperp)*mctrackperp->Pt();
		  deltaEtaSqPtPerp += DeltaEtaSqMC(fPerpCone, mctrackperp)*mctrackperp->Pt();
		  totalTrackPtPerp += mctrackperp->Pt();
		}// end loop over AODMC particles
	      if(!IsEqualRel(totalTrackPtPerp, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPtPerp!=0.0
		{
		  firstMomDeltPhiPerp = deltaPhiPtPerp/totalTrackPtPerp;
		  firstMomDeltEtaPerp = deltaEtaPtPerp/totalTrackPtPerp;
		  secondMomDeltPhiPerp = deltaPhiSqPtPerp/totalTrackPtPerp;
		  secondMomDeltEtaPerp = deltaEtaSqPtPerp/totalTrackPtPerp;
		  secondCentralPhiPerp = secondMomDeltPhiPerp - firstMomDeltPhiPerp*firstMomDeltPhiPerp;
		  secondCentralEtaPerp = secondMomDeltEtaPerp - firstMomDeltEtaPerp*firstMomDeltEtaPerp;
		  secondCentralRPerp = secondCentralPhiPerp + secondCentralEtaPerp;
		} // end if totalTrackPt!=0.0
	      if(IsEqualRel(totalTrackPtPerp, 0.0)) //!IsEqualRel(totalTrackPtPerp, 0.0) //totalTrackPtPerp==0.0
		secondCentralRPerp = 10.0; //overflow
	    }
	  ///// end of adding the SCM for the perpendicular cone

	  if(mcjets==1) // if only one jet in the whole event, and inside acceptance
	    {
	      // reference multiplicity stuff in pp, also filled in PbPb, but does not matter.
	      // set to: V0 like corrected multiplicity: correctedV0LikeMult
	      if(correctedV0LikeMult<25)
		{
		  fNChTrRDMultOJMC[0]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[0]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=25&&correctedV0LikeMult<50)
		{
		  fNChTrRDMultOJMC[1]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[1]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=50&&correctedV0LikeMult<90)
		{
		  fNChTrRDMultOJMC[2]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[2]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=90&&correctedV0LikeMult<120)
		{
		  fNChTrRDMultOJMC[3]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[3]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=120&&correctedV0LikeMult<150)
		{
		  fNChTrRDMultOJMC[4]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[4]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=150&&correctedV0LikeMult<200)
		{
		  fNChTrRDMultOJMC[5]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[5]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=200&&correctedV0LikeMult<300)
		{
		  fNChTrRDMultOJMC[6]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[6]->Fill(secondCentralR,jetPt);
		}
	      if(correctedV0LikeMult>=300)
		{
		  fNChTrRDMultOJMC[7]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultOJMC[7]->Fill(secondCentralR,jetPt);
		}
	      //Results for inclusive jets
	      // 2nd. Reference: set to: TPC tracks minus jet, minus dijet area
	      if(softRefMcNoJets<5)
		{
		  fNChTrRDMultSEOJMC[0]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[0]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=5&&softRefMcNoJets<10)
		{
		  fNChTrRDMultSEOJMC[1]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[1]->Fill(secondCentralR,jetPt); 
		}
	      if(softRefMcNoJets>=10&&softRefMcNoJets<15)
		{
		  fNChTrRDMultSEOJMC[2]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[2]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=15&&softRefMcNoJets<20)
		{
		  fNChTrRDMultSEOJMC[3]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[3]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=20&&softRefMcNoJets<30)
		{
		  fNChTrRDMultSEOJMC[4]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[4]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=30&&softRefMcNoJets<40)
		{
		  fNChTrRDMultSEOJMC[5]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[5]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=40&&softRefMcNoJets<50)
		{
		  fNChTrRDMultSEOJMC[6]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[6]->Fill(secondCentralR,jetPt);
		}
	      if(softRefMcNoJets>=50)
		{	    
		  fNChTrRDMultSEOJMC[7]->Fill(nTracksPerc,jetPt);
		  fSCMRDMultSEOJMC[7]->Fill(secondCentralR,jetPt);
		}
	    }
	  //End results for inclusive jets,starts parton by parton
	  
	  switch(abs(flavor))
	    {
	    case 1:
	      fNChTr[0]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[0]->Fill(jetPt,evtype);
	      fHistPtParton[0]->Fill(jetPt);
	      fSCM[0]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCQuark->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.)) //!IsEqualRel(mctrack->Pt(), 0.0) // fCurrentJetMinPtNT90!=7000.
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.)) // !IsEqualRel(fCurrentJetMinPtNT90, 7000.) // fCurrentJetMinPtNT90Recalc!=7000.
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;
	    case 2:
	      fNChTr[1]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[1]->Fill(jetPt,evtype);
	      fHistPtParton[1]->Fill(jetPt);
	      fSCM[1]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCQuark->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;
	    case 3:
	      fNChTr[2]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[2]->Fill(jetPt,evtype);
	      fHistPtParton[2]->Fill(jetPt);
	      fSCM[2]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCQuark->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;
	    case 4:
	      fNChTr[3]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[3]->Fill(jetPt,evtype);
	      fHistPtParton[3]->Fill(jetPt);
	      fSCM[3]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCQuark->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;
	    case 5:
	      fNChTr[4]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[4]->Fill(jetPt,evtype);
	      fHistPtParton[4]->Fill(jetPt);
	      fSCM[4]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCQuark->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;
	    case 21:
	      fNChTr[5]->Fill(nTracksPerc,jetPt);
	      fProcessPDG[5]->Fill(jetPt,evtype);
	      fHistPtParton[5]->Fill(jetPt);
	      fSCM[5]->Fill(secondCentralR,jetPt);
	      fMinTrackPtInNTXh[0]->Fill(fMinTrackPtInNTX,jetPt,1); // 0 for pp MC
	      fMaxTrackPtInNTXh[0]->Fill(fMaxTrackPtInNTX,jetPt); // 0 for MC
	      if(fIsPossibleToSubstBckg)
		{
		  fNChTrCorrMCGluon->Fill(nTracksPercBckgSubst,pTbs); 
		  fSCMMCPerp->Fill(secondCentralR,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrMCPerp->Fill(nTracksAboveThresholdPerp,jetPt);
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    fNChTrCorrMCPerp->Fill(nTrUpThrPerpBckSubs,pTbs);
		}
	      break;	      
	    default:
	      break;
	    }
	  AliDebug(4,Form("Sabor del jet numero:%d es: %d y se necesitaron %d tracks \n",indxmc,flavor,nTracksPerc)); 
	  imcj++;
	} // MC jets for cycle
      nGenJets=imcj;
      for(Int_t u=0 ; u<mcJetCounter  ;u++)
	{
	  if(u<7)
	    fJetsMultPtMC->Fill(jetPts[u],mcJetCounter);
	}
      // if(fEnablePrints)
      //  {
      //   if(mcJetCounter>=3)
      //    printf("%i Jets inside acceptance at event number:%i \n",mcJetCounter,fEvtCount-1);
      //  }
      fNAccJetsMC->Fill(mcJetCounter,mcJetCounter);
    } // end if MC info in AOD
  
  if(!fUseOnlyMC) 
    {  
      // Primero que todo, debe de ir la seleccion de eventos reconstruidos:
      // 1. Que tenga un vertice reconstruido dentro de 10 cm.
      // Vertex info for reconstructed events
      AliAODVertex *pvx = fAOD->GetPrimaryVertex();
      if(!pvx)
	{
	  AliError("No primary vertex!");
	  return;
	}
      if(TMath::Abs(pvx->GetZ())>10.) // if the event vertex is larger than 10 cm, reject
	return;
      fZVertex->Fill(pvx->GetZ(),pvx->GetZ()); // vertex, provide number of accepted events as entries for reco jets

      ///////////////////////////////////////
      // SECONDARY RECO BRANCH STUFF       // 
      // Get the secondary branch with the reconstructed jets 
      if(fBranchSecRec!="")
	{
	  AliDebug(4,Form("fBranchSecRec was not default \n")); 
	  TClonesArray *aodSecRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchSecRec.Data()));
	  if(!aodSecRecJets)
	    {
	      AliError(Form("%s:%d no reconstructed Secondary Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchSecRec.Data())); 
	      return;  //stop the analysis
	    }
	  AliDebug(4,Form("There are %d reconstructed jets from the secondary branch in this event \n", aodSecRecJets->GetEntries())); 
	  Int_t recojetsSEC =  aodSecRecJets->GetEntries();
	  fNJetsRDSeco->Fill(recojetsSEC,recojetsSEC);  // number of jets in the secondary branch

	  HasOverlapedCones(aodSecRecJets); // Procedure for tagging usable jets
	                                    // Up 16 jets are flagged
	  
	  AliDebug(4,"Antes de realizar el loop jets reconstruidos del segundo branch \n"); 
	  Int_t secondjetacccounter = 0;
	  for (Int_t IDXS = 0; IDXS < recojetsSEC; IDXS++) 
	    {
	      AliDebug(4,Form("Number of current jet:%i \n",IDXS));
	      AliAODJet *rjetsec = dynamic_cast<AliAODJet*>(aodSecRecJets->At(IDXS));
	      if (!rjetsec) 
		{
		  AliDebug(2,Form("ERROR: Could not receive jet %d from the second branch\n", IDXS)); 
		  continue;
		}
	      
	      ///////////////////////////////////////////////////////////////////////////////
	      ///// Part for Chritians plot of inclusive and leading jets comp at 2.76 TeV //
	      if(!IsInsideAcceptance(rjetsec))  // old condition
		continue;
	      if(IDXS==0) // leading jet
		fSecRecJetPtLeading->Fill(rjetsec->Pt());
	      fSecRecJetPtInclusive->Fill(rjetsec->Pt()); // all
	      ///// End of Christians Plot reco 2nd branch
	      ///////////////////////////////////////////////////////////////////////////////

	      if(IDXS>15)
		continue;

	      if(!fJetFlags[IDXS]) // If the jet is flaged kFALSE, not usable
		continue;

	      fJetPtSec->Fill(rjetsec->Pt());
	      fJetEtaSec->Fill(rjetsec->Eta(),rjetsec->Eta());
	      fJetPhiSec->Fill(rjetsec->Phi(),rjetsec->Phi());
	      secondjetacccounter++;
	    }
	  fNAccJetsRDSeco->Fill(secondjetacccounter,secondjetacccounter);
	}
      // END OF SECONDARY BRANCH STUFF     //
      ///////////////////////////////////////

      // Get the branch with the reconstructed jets
      TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRec.Data()));
      if(!aodRecJets)
	{
	  AliError(Form("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data())); 
	  return;
	}
      
      AliDebug(4,Form("There are %d reconstructed jets in this event\n", aodRecJets->GetEntries())); 
      Int_t recojets =  aodRecJets->GetEntries();
      fNJetsRD->Fill(recojets,recojets); // numero de jets directamente del branch

      HasOverlapedCones(aodRecJets); // Procedure for tagging usable jets
                                    // Up 16 jets are flagged

      AliDebug(4,"Antes de realizar el loop sobre AOD tracks \n"); 
      // Loop over AOD tracks
      Int_t tracksAOD = fAOD->GetNumberOfTracks(); 
      AliDebug(4,Form("Numero de tracks en el AOD:%d \n",tracksAOD));
      Double_t aodtracketa = 0.0;
      perpendicularPt = 0.0;
      px=0.0;
      py=0.0;
      pz=0.0;
      en=0.0;
      pTbs=0.0;
      etabs=0.0;
      phibs=0.0;
      fBckgSbsJet[0]=0.0;
      fBckgSbsJet[1]=0.0;
      fBckgSbsJet[2]=0.0;
      Int_t refNJMult = 0;
      Int_t myTotalMultRef = 0; 
      Int_t myTotalSoftMultRef = 0; 
      for(Int_t aodT = 0; aodT < tracksAOD; aodT++ )
	{
	  AliAODTrack *aodtrack = fAOD->GetTrack(aodT);
	  if(!aodtrack) continue;
	  aodtracketa = TMath::Abs(aodtrack->Eta());
	  if(aodtracketa>0.9) continue;
	  if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
	  if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries, maybe is redundant with the previous selection...
	  fEtaAOD->Fill(aodtrack->Eta(),aodtrack->Eta());
	  fPhiAOD->Fill(aodtrack->Phi(),aodtrack->Phi());
	  fPtAOD->Fill(aodtrack->Pt(),aodtrack->Pt());
	  if(fJetEvent) // if has an accepted jet, calculate the perpendicular cone
	    {
	      if(HasPerpendicularCone()) // If there is a perpendicular cone available
		{
		  if(aodtrack->Pt()>fMinpTVal)
		    {
		      if(GetDeltaR(fEtaPerpCoord,fPhiPerpCoord,aodtrack->Eta(),aodtrack->Phi())<fJetRadius)
			perpendicularPt = perpendicularPt + aodtrack->Pt();
		    }
		}
	    } // end if jet event
	  //Total TPC multiplicity of primaries
	  myTotalMultRef++;
	  if(aodtrack->Pt()<fMinpTValUE) continue; // pT cut  fMinpTValUE
	  if(aodtrack->Pt()>fMaxpTValUE) continue; // pT cut  fMaxpTValUE
	  myTotalSoftMultRef++;  
	  if(!IsTrackInsideExcludedArea(aodtrack->Eta(), aodtrack->Phi(), aodRecJets)) 
	    refNJMult++;
	} // end track loop over the event...

      fRefMultWOJet->Fill(refMultiplicity,refNJMult); 
      fMultWOJetVZero->Fill(refNJMult,multFullV0);
      Double_t v0CorrMult = multFullV0 - GetV0ExcludedMultiplicity(aodRecJets);
      fRefMultFullV0->Fill(refMultiplicity,multFullV0);
      fRefMultV0Corr->Fill(refMultiplicity,v0CorrMult);
      fFullV0V0Corr->Fill(multFullV0,v0CorrMult);
      fRefAODTrackCount->Fill(refMultiplicity,myTotalMultRef);
      fTrackCountWOJet->Fill(myTotalMultRef,refNJMult);

      if(recojets==1) // correlation for only monojet events
	{
	  fRefMultFullV0UJ->Fill(refMultiplicity,multFullV0);
	  fRefMultV0CorrUJ->Fill(refMultiplicity,v0CorrMult);
	  fFullV0V0CorrUJ->Fill(multFullV0,v0CorrMult);
	  fMultWOJetVZeroUJ->Fill(refNJMult,multFullV0);
	  fRefMultWOJetUJ->Fill(refMultiplicity,refNJMult);
	  fRefAODTrackCountUJ->Fill(refMultiplicity,myTotalMultRef);
	  fTrackCountWOJetUJ->Fill(myTotalMultRef,refNJMult); 
	}

      if(fJetEvent) // if has an accepted jet, calculate the perpendicular cone
	{
	  if(HasPerpendicularCone()) // If there is a perpendicular cone available
	    {
	      px = perpendicularPt*TMath::Cos(fPhiPerpCoord);
	      py = perpendicularPt*TMath::Sin(fPhiPerpCoord);
	      pz = perpendicularPt/TMath::Tan(2.0*TMath::ATan(TMath::Exp(-fEtaPerpCoord)));
	      en = TMath::Sqrt(px*px + py*py + pz*pz);
	      fPerpCone->SetPxPyPzE(px, py, pz, en);
	    }
	}
            
      // Loop to fill a pT spectrum of the reco jets
      Int_t irecj=0; // index for reconstructed jets to correlate
      Int_t nrectracks[6]={0};
      Double_t ptrecjet[6]={0};
      Double_t scmr[6]={0};
      Double_t aodtrackxi=0;
      Int_t ntxreco;
      Int_t nTRecAboveThresholdPerp=0; 
      Int_t ntxrecoRecalc;
      Int_t nTRecAboveThresholdPerpRecalc=0; 
	    
      for(Int_t i=0; i<6; i++) // Reset per event
	{
	  fHistContainerR4[i]->Reset();
	  fHistContainerR3[i]->Reset();
	  fHistContainerR2[i]->Reset();
	}
      
      Double_t jetPtsR[7]={0};  // to store the pt of the jets
      Int_t rJetCounter=0;    // counter of accepted reco jets 
      fIsPossibleToSubstBckg = kTRUE; // Initialize before the loop
      if(fJetEvent) // si tiene jets validos
	{
	  if(!HasPerpendicularCone()) // pero no encontro un cono perpendicular libre
	    fIsPossibleToSubstBckg = kFALSE; // if not perpendicular cone, set to kFALSE, so no perpendicular calculations available
	}
      
      AliDebug(4,"Antes de realizar el loop jets reconstruidos \n"); 
      for (Int_t indxrec = 0; indxrec < recojets; indxrec++) 
	{
	  AliDebug(4,Form("Number of current jet:%i \n",indxrec));
	  ntxreco = 0;
	  ntxrecoRecalc = 0;
	  fMinTrackPtInNTX=200.0;  //Initialize for each jet, overflown
	  fMaxTrackPtInNTX=200.0;  //Initialize for each jet, overflown
	  fMinTrackPtInNTXR=200.0;  //Initialize for each jet, overflown
	  fMaxTrackPtInNTXR=200.0;  //Initialize for each jet, overflown
	  deltaPhiPt = 0.0;
	  deltaEtaPt = 0.0;
	  deltaPhiSqPt = 0.0;
	  deltaEtaSqPt = 0.0;
	  totalTrackPt = 0.0;
	  firstMomDeltPhi = 0.0;
	  firstMomDeltEta = 0.0;
	  secondMomDeltPhi = 0.0;
	  secondMomDeltEta = 0.0;
	  secondCentralPhi = 0.0;
	  secondCentralEta = 0.0;
	  secondCentralR = 0.0;
	  
	  AliAODJet *rjet = dynamic_cast<AliAODJet*>(aodRecJets->At(indxrec));
	  if (!rjet) 
	    {
	      AliDebug(2,Form("ERROR: Could not receive jet %d\n", indxrec)); 
	      continue;
	    }
	  fJetEtaAll->Fill(rjet->Eta());// all jets

	  ///////////////////////////////////////////////////////////////////////////////
	  ///// Part for Chritians plot of inclusive and leading jets comp at 2.76 TeV //
	  if(!IsInsideAcceptance(rjet))  // old condition
	    continue;
	  if(indxrec==0) // leading jet
	    fRecJetPtLeading->Fill(rjet->Pt());
	  fRecJetPtInclusive->Fill(rjet->Pt()); // all
	  fJetEtaOnlyTPCcut->Fill(rjet->Eta());// only eta acceptance cut for TPC
	  ///// End of Christians Plot reco
	  ///////////////////////////////////////////////////////////////////////////////

	  if(indxrec>15)
	    continue;

	  if(!fJetFlags[indxrec]) // If the jet is flaged kFALSE, not usable
	    continue;
	  
	  AliDebug(4,Form("Jet #%i is in the acceptance \n",indxrec));
	  if(rJetCounter<7)
	    jetPtsR[rJetCounter]=rjet->Pt();
	  rJetCounter++;
	  fJetPt->Fill(rjet->Pt());
	  fJetEta->Fill(rjet->Eta(),rjet->Eta());
	  fJetPhi->Fill(rjet->Phi(),rjet->Phi());

	  if(rjet->Pt()>10.)
	    fJetEtaJetPt[0]->Fill(rjet->Eta());
	  if(rjet->Pt()>30.)
	    fJetEtaJetPt[1]->Fill(rjet->Eta());
	  if(rjet->Pt()>50.)
	    fJetEtaJetPt[2]->Fill(rjet->Eta());
	  
	  // Reco RefTracks check
	  Bool_t rTrkFlagRec = kFALSE;
	  Int_t trkinrecjet = rjet->GetRefTracks()->GetEntriesFast();
	  if(trkinrecjet!=0&&!fForceNotTR)
	    rTrkFlagRec = kTRUE;
	  AliDebug(4,Form("Number of tracks in RefTracks reco jet:%i \n",trkinrecjet));
	  if(rTrkFlagRec)
	    {
	      // Check the properties of the tracks in this jet with track refs
	      AliDebug(4,Form("Checking composition in Reco jets with track refs")); 
	      for(Int_t aodT = 0; aodT <trkinrecjet; aodT++ )
		{
		  aodtrackxi=0;
		  AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(rjet->GetRefTracks()->At(aodT));
		  if(!aodtrack)
		    {
		      AliError("Error, no AOD Track!");
		      continue;
		    }
		  if(!aodtrack->TestFilterBit(fFilterBit))
		    {
		      //		      printf("Rejecting track from track refs due to wrong filterbit! \n");
		      continue; //track filter selection
		    }
		  if(!aodtrack->IsPrimaryCandidate())
		    {
		      //		      printf("Rejecting track from track refs due to no primary candidate! \n");
		      continue; // only primaries
		    }
		  deltaPhiPt += DeltaPhiTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaEtaPt += DeltaEtaTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaPhiSqPt += DeltaPhiSqTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaEtaSqPt += DeltaEtaSqTrack(rjet, aodtrack)*aodtrack->Pt();
		  totalTrackPt += aodtrack->Pt();
		  
		  if(!IsEqualRel(aodtrack->Pt(), 0.0)) //!IsEqualRel(totalTrackPtPerp, 0.0) //aodtrack->Pt()!=0
		    aodtrackxi= log(rjet->Pt()/aodtrack->Pt());
		  if(irecj<maxJetNum)
		    {
		      fHistContainerR4[irecj]->Fill(aodtrackxi,rjet->Pt());
		      if(!IsTrackInsideThisJet(aodtrack, rjet, 0.3)) continue;
		      fHistContainerR3[irecj]->Fill(aodtrackxi,rjet->Pt());
		      if(!IsTrackInsideThisJet(aodtrack, rjet, 0.2)) continue;
		      fHistContainerR2[irecj]->Fill(aodtrackxi,rjet->Pt());
		    }	    
		} //end loop over track references
	      if(!IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0)  // totalTrackPt!=0.0
		{
		  firstMomDeltPhi = deltaPhiPt/totalTrackPt;
		  firstMomDeltEta = deltaEtaPt/totalTrackPt;
		  secondMomDeltPhi = deltaPhiSqPt/totalTrackPt;
		  secondMomDeltEta = deltaEtaSqPt/totalTrackPt;
		  secondCentralPhi = secondMomDeltPhi - firstMomDeltPhi*firstMomDeltPhi;
		  secondCentralEta = secondMomDeltEta - firstMomDeltEta*firstMomDeltEta;
		  secondCentralR = secondCentralPhi + secondCentralEta;
		} // end if totalTrackPt!=0.0
	      if(IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPt==0.0
		secondCentralR = 10.0; //overflow value
	    } // end if there are track references
	  
	  if(!rTrkFlagRec)
	    {
	      // Check properties of the tracks in this jet without track refs
	      AliDebug(4,Form("Checking composition in Reco jets without track refs")); 
	      for(Int_t aodT = 0; aodT < tracksAOD; aodT++ )
		{
		  AliAODTrack *aodtrack = fAOD->GetTrack(aodT);
		  if(!aodtrack) continue;
		  if(!IsTrackInsideThisJet(aodtrack, rjet, jfr)) continue;
		  if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
		  if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries, perhaps redundant with the previous
		  if(aodtrack->Pt()<fMinpTVal) continue; //DATA: PT CUT
		  deltaPhiPt += DeltaPhiTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaEtaPt += DeltaEtaTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaPhiSqPt += DeltaPhiSqTrack(rjet, aodtrack)*aodtrack->Pt();
		  deltaEtaSqPt += DeltaEtaSqTrack(rjet, aodtrack)*aodtrack->Pt();
		  totalTrackPt += aodtrack->Pt();
		  
		  if(!IsEqualRel(aodtrack->Pt(), 0.0)) //!IsEqualRel(totalTrackPt, 0.0) //aodtrack->Pt()!=0
		    aodtrackxi= log(rjet->Pt()/aodtrack->Pt());
		  if(irecj<maxJetNum)
		    {
		      fHistContainerR4[irecj]->Fill(aodtrackxi,rjet->Pt());
		      if(!IsTrackInsideThisJet(aodtrack, rjet, 0.3)) continue;
		      fHistContainerR3[irecj]->Fill(aodtrackxi,rjet->Pt());
		      if(!IsTrackInsideThisJet(aodtrack, rjet, 0.2)) continue;
		      fHistContainerR2[irecj]->Fill(aodtrackxi,rjet->Pt());
		    }
		} // end loop over tracks
	      if(!IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) //totalTrackPt!=0.0
		{
		  firstMomDeltPhi = deltaPhiPt/totalTrackPt;
		  firstMomDeltEta = deltaEtaPt/totalTrackPt;
		  secondMomDeltPhi = deltaPhiSqPt/totalTrackPt;
		  secondMomDeltEta = deltaEtaSqPt/totalTrackPt;
		  secondCentralPhi = secondMomDeltPhi - firstMomDeltPhi*firstMomDeltPhi;
		  secondCentralEta = secondMomDeltEta - firstMomDeltEta*firstMomDeltEta;
		  secondCentralR = secondCentralPhi + secondCentralEta;
		} // end if totalTrackPt!=0.0 
	      if(IsEqualRel(totalTrackPt, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPt==0.0
		secondCentralR = 10.0; //overflow value
	    } // end of no track references
	  //Esto es lo anterior, toma al jet como es, y calcula NT90	  
	  ntxreco=GetNumberOfChargedTracks(ntx,rjet, tracksAOD, fAOD, jfr); // this call fixes the minimum pT track
	  //Y aqui calcula cuantos tracks se necesitan arriba del threshold establecido en la linea anterior
	  //esto debe ser para cada jet. Lo unico que se calcula una sola vez es el cono perpendicular  
	  if(fIsPossibleToSubstBckg&&!IsEqualRel(fCurrentJetMinPtNT90, 7000.)) //and only if the method worked
	    nTRecAboveThresholdPerp = GetNRecChargedTracksAboveThreshold(fPerpCone,tracksAOD, fAOD,jfr); //here one changes NTX

	  // correct the jet pT
	  if(fIsPossibleToSubstBckg) // If there is a perpendicular cone available, substract backg and fill the new jet pT
	    {
	      pTbs= rjet->Pt()-fPerpCone->Pt();
	      etabs= rjet->Eta();
	      phibs= rjet->Phi();
	      fBckgSbsJet[0]=pTbs; //pT
	      fBckgSbsJet[1]=etabs; //eta
	      fBckgSbsJet[2]=phibs; //phi
	      // Now re-calculate nt90 for the energy corrected jet
	      ntxrecoRecalc = GetRecalcNTXRec(ntx,rjet, tracksAOD, fAOD, jfr); //This call saves the new min pT
	      // Now re-calculate the perpendicular cone NT90 background
	      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.)) //calculating for the current jet, if the previous method worked
		nTRecAboveThresholdPerpRecalc = GetRecalcNRecChTrUpThr(fPerpCone,tracksAOD, fAOD,jfr);
	    }

	  // SCM perpendicular cone
	  if(fIsPossibleToSubstBckg)
	    {
	      // To make sure, re-initialize
	      deltaPhiPtPerp = 0.0;
	      deltaEtaPtPerp = 0.0;
	      deltaPhiSqPtPerp = 0.0;
	      deltaEtaSqPtPerp = 0.0;
	      totalTrackPtPerp = 0.0;
	      firstMomDeltPhiPerp = 0.0;
	      firstMomDeltEtaPerp = 0.0;
	      secondMomDeltPhiPerp = 0.0;
	      secondMomDeltEtaPerp = 0.0;
	      secondCentralPhiPerp = 0.0;
	      secondCentralEtaPerp = 0.0;
	      secondCentralRPerp = 0.0;
	      AliDebug(4,Form("Checking SCM of perpendicular cone in Reco jets")); 
	      for(Int_t aodTperp = 0; aodTperp < tracksAOD; aodTperp++ )
		{ //fPerpCone
		  AliAODTrack *aodtrackperprec = fAOD->GetTrack(aodTperp);
		  if(!aodtrackperprec) continue;
		  if(!IsTrackInsideThisJet(aodtrackperprec, fPerpCone, jfr)) continue;
		  if(!aodtrackperprec->TestFilterBit(fFilterBit)) continue; //track filter selection
		  if(!aodtrackperprec->IsPrimaryCandidate()) continue; // only primaries, perhaps redundant with the previous
		  if(aodtrackperprec->Pt()<fMinpTVal) continue; //DATA: PT CUT
		  deltaPhiPtPerp += DeltaPhiTrack(fPerpCone, aodtrackperprec)*aodtrackperprec->Pt();
		  deltaEtaPtPerp += DeltaEtaTrack(fPerpCone, aodtrackperprec)*aodtrackperprec->Pt();
		  deltaPhiSqPtPerp += DeltaPhiSqTrack(fPerpCone, aodtrackperprec)*aodtrackperprec->Pt();
		  deltaEtaSqPtPerp += DeltaEtaSqTrack(fPerpCone, aodtrackperprec)*aodtrackperprec->Pt();
		  totalTrackPtPerp += aodtrackperprec->Pt();
		} // end loop over tracks
	      if(!IsEqualRel(totalTrackPtPerp, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPtPerp!=0.0
		{
		  firstMomDeltPhiPerp = deltaPhiPtPerp/totalTrackPtPerp;
		  firstMomDeltEtaPerp = deltaEtaPtPerp/totalTrackPtPerp;
		  secondMomDeltPhiPerp = deltaPhiSqPtPerp/totalTrackPtPerp;
		  secondMomDeltEtaPerp = deltaEtaSqPtPerp/totalTrackPtPerp;
		  secondCentralPhiPerp = secondMomDeltPhiPerp - firstMomDeltPhiPerp*firstMomDeltPhiPerp;
		  secondCentralEtaPerp = secondMomDeltEtaPerp - firstMomDeltEtaPerp*firstMomDeltEtaPerp;
		  secondCentralRPerp = secondCentralPhiPerp + secondCentralEtaPerp;
		} // end if totalTrackPt!=0.0
	      if(IsEqualRel(totalTrackPtPerp, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // totalTrackPtPerp==0.0
		secondCentralRPerp = 10.0; //overflow
	    } 

	  ///// end of adding the SCM for the perpendicular cone

	  if(irecj<maxJetNum)
	    {
	      recJets[irecj]= *rjet;
	      nrectracks[irecj] = ntxreco;
	      ptrecjet[irecj] = rjet->Pt();
	      scmr[irecj] = secondCentralR;
	      AliDebug(4,Form("Para el jet reco num: %d se necesitaron %d tracks \n",irecj,nrectracks[irecj])); 
	    }
	  AliDebug(4,"Before filling the histograms for this jet \n");
	  fNChTrRD->Fill(ntxreco,rjet->Pt());
	  fSCMRD->Fill(secondCentralR,rjet->Pt());
	  fProfNChTrRD->Fill(rjet->Pt(),ntxreco);

	  fNTXV0MultPt->Fill(ntxreco,v0CorrMult,rjet->Pt());
	  fNTXCBMultPt->Fill(ntxreco,refNJMult,rjet->Pt());

	  //refNJMult
	  // reference multiplicity stuff in pp, also filled in PbPb, but does not matter.
	  // set to: V0 corrected multiplicity: v0CorrMult
	  if(v0CorrMult<25)
	    {
	      fNChTrRDMult[0]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[0]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[0]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[0]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[0]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=25&&v0CorrMult<50)
	    {
	      fNChTrRDMult[1]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[1]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[1]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[1]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[1]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=50&&v0CorrMult<90)
	    {
	      fNChTrRDMult[2]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[2]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[2]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[2]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[2]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=90&&v0CorrMult<120)
	    {
	      fNChTrRDMult[3]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[3]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[3]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[3]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[3]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=120&&v0CorrMult<150)
	    {
	      fNChTrRDMult[4]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[4]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[4]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[4]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[4]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=150&&v0CorrMult<200)
	    {
	      fNChTrRDMult[5]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[5]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[5]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[5]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[5]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=200&&v0CorrMult<300)
	    {
	      fNChTrRDMult[6]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[6]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[6]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[6]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[6]->Fill(secondCentralR,rjet->Pt());
		}
	    }
	  if(v0CorrMult>=300)
	    {
	      fNChTrRDMult[7]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMult[7]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetCharge[7]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultOJ[7]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultOJ[7]->Fill(secondCentralR,rjet->Pt());
		}
	    }

	  // 2nd. Reference: set to: TPC tracks minus jet, minus dijet area
	  if(refNJMult<5) //&&refNJMult>1
	    {
	      fNChTrRDMultSE[0]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[0]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[0]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[0]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[0]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[0]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,1); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,1); //filling mult bin
		    }
		}
	    }
	  if(refNJMult>=5&&refNJMult<10)
	    {
	      fNChTrRDMultSE[1]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[1]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[1]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[1]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[1]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[1]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,2); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,2); //filling mult bin
		    }
		}
	    }
	  if(refNJMult>=10&&refNJMult<15)
	    {
	      fNChTrRDMultSE[2]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[2]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[2]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[2]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[2]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[2]->Fill(nTRecAboveThresholdPerp,rjet->Pt()); 
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,3); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,3); //filling mult bin
		    } 
		}
	    }
	  if(refNJMult>=15&&refNJMult<20)
	    {
	      fNChTrRDMultSE[3]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[3]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[3]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[3]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[3]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[3]->Fill(nTRecAboveThresholdPerp,rjet->Pt()); 
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,4); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,4); //filling mult bin
		    } 
		}
	    }
	  if(refNJMult>=20&&refNJMult<30)
	    {
	      fNChTrRDMultSE[4]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[4]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[4]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[4]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[4]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[4]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,5); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,5); //filling mult bin
		    }
		}
	    }
	  if(refNJMult>=30&&refNJMult<40)
	    {
	      fNChTrRDMultSE[5]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[5]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[5]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[5]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[5]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[5]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,6); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,6); //filling mult bin
		    }
		}
	    }
	  if(refNJMult>=40&&refNJMult<50)
	    {
	      fNChTrRDMultSE[6]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[6]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[6]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[6]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[6]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[6]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,7); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,7); //filling mult bin
		    }
		}
	    }
	  if(refNJMult>=50)
	    {
	      fNChTrRDMultSE[7]->Fill(ntxreco,rjet->Pt());
	      fSCMRDMultSE[7]->Fill(secondCentralR,rjet->Pt());
	      fTotalJetChargeSE[7]->Fill(fCurrentJetCharge);
	      if(recojets==1) // if only one jet in the whole event, and inside acceptance
		{
		  fNChTrRDMultSEOJ[7]->Fill(ntxreco,rjet->Pt());
		  fSCMRDMultSEOJ[7]->Fill(secondCentralR,rjet->Pt());
		  if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		    fNChTrRecPerpMultSEOJ[7]->Fill(nTRecAboveThresholdPerp,rjet->Pt());  
		  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
		    {
		      fNChTrRecECorrPPMult->Fill(ntxrecoRecalc,pTbs,8); //filling mult bin
		      fNChTrRecPerpECorrPPMult->Fill(nTRecAboveThresholdPerpRecalc,pTbs,8); //filling mult bin
		    }
		}
	    }

	  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
	    {
	      if(!IsEqualRel(fCurrentJetMinPtNT90, 7000.))
		fNChTrRecPerp->Fill(nTRecAboveThresholdPerp,rjet->Pt());  // these are my previous histos
	      fSCMRecPerp->Fill(secondCentralRPerp,rjet->Pt());  // this are my previous histos
	      if(!fIsHIevent) //if is a proton proton event
		{
		  fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,1); //filling proton bin
		  if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
		    {
		      fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,1); //filling proton bin
		      if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),1,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 1);
			}
		      if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),2,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 2);
			}
		      if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),3,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 3);
			}
		      if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),4,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 4);
			}
		      if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),5,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 5);
			}
		      if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),6,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 6);
			}
		      if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),7,1);
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 7);
			}
		      if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			{
			  fPtInPerpCon->Fill(fPerpCone->Pt(),8,1); 
			  FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 1, 8);
			}
		    }
		}
	      if(fIsHIevent) //if is a PbPb event     
		{
		  if(fEventCent>=0&&fEventCent<10.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,2); //filling first centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),2); 
		      fJetPtCentPbPbCorr->Fill(pTbs,2); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,2); //filling first centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),2); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),2); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,2);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 2, 8);
			    }
			}
		    }
		  if(fEventCent>=10&&fEventCent<20.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,3); //filling second centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),3); 
		      fJetPtCentPbPbCorr->Fill(pTbs,3); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,3); //filling second centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),3); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),3); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,3);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 3, 8);
			    }
			}
		    }
		  if(fEventCent>=20&&fEventCent<30.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,4); //filling third centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),4); 
		      fJetPtCentPbPbCorr->Fill(pTbs,4); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,4); //filling third centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),4); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),4); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,4);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 4, 8);
			    }
			}
		    }
		  if(fEventCent>=30&&fEventCent<40.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,5); //filling fourth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),5); 
		      fJetPtCentPbPbCorr->Fill(pTbs,5); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,5); //filling fourth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),5); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),5); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,5);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 5, 8);
			    }
			}
		    }
		  if(fEventCent>=40&&fEventCent<50.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,6); //filling fourth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),6); 
		      fJetPtCentPbPbCorr->Fill(pTbs,6); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,6); //filling fourth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),6); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),6); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,6);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 6, 8);
			    }
			}
		    }
		  if(fEventCent>=50&&fEventCent<60.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,7); //filling fourth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),7); 
		      fJetPtCentPbPbCorr->Fill(pTbs,7); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,7); //filling fourth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),7); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),7); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,7);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 7, 8);
			    }
			}
		    }
		  if(fEventCent>=60&&fEventCent<70.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,8); //filling fourth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),8); 
		      fJetPtCentPbPbCorr->Fill(pTbs,8); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,8); //filling fourth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),8);
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),8);  
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,8);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 8, 8);
			    }
			}
		    }
		  if(fEventCent>=70&&fEventCent<80.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,9); //filling fourth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),9); 
		      fJetPtCentPbPbCorr->Fill(pTbs,9); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,9); //filling fourth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),9); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),9); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,9);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 9, 8);
			    }
			}
		    }
		  if(fEventCent>=80&&fEventCent<100.)
		    {
		      fNChTrRecECorr->Fill(ntxrecoRecalc,pTbs,10); //filling sixth centrality bin
		      fJetPtCentPbPbRaw->Fill(rjet->Pt(),10); 
		      fJetPtCentPbPbCorr->Fill(pTbs,10); 
		      if(!IsEqualRel(fCurrentJetMinPtNT90Recalc, 7000.))
			{
			  fNChTrRecPerpECorr->Fill(nTRecAboveThresholdPerpRecalc,pTbs,10); //filling sixth centrality bin
			  fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),10); 
			  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),10); 
			  if((rjet->Pt()>10.)&&(rjet->Pt()<20.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),1,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 1);
			    }
			  if((rjet->Pt()>20.)&&(rjet->Pt()<30.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),2,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 2);
			    }
			  if((rjet->Pt()>30.)&&(rjet->Pt()<40.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),3,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 3);
			    }
			  if((rjet->Pt()>40.)&&(rjet->Pt()<50.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),4,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 4);
			    }
			  if((rjet->Pt()>50.)&&(rjet->Pt()<60.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),5,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 5);
			    }
			  if((rjet->Pt()>60.)&&(rjet->Pt()<70.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),6,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 6);
			    }
			  if((rjet->Pt()>70.)&&(rjet->Pt()<80.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),7,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 7);
			    }
			  if((rjet->Pt()>80.)&&(rjet->Pt()<90.))
			    {
			      fPtInPerpCon->Fill(fPerpCone->Pt(),8,10);
			      FillPerpConeHisto(fPtDistInPerpConeRaw, tracksAOD, fAOD, 10, 8);
			    }
			}
		    }
		}
	    }

	  AliDebug(4,"Before filling the histograms for this jet of min and max pT in NTX \n"); 
	  AliDebug(4,Form("Min:%f, Max:%f \n",fMinTrackPtInNTX,fMaxTrackPtInNTX));
	  fMinTrackPtInNTXh[1]->Fill(fMinTrackPtInNTX,rjet->Pt(),1); // for rec pp or inclusive PbPb
	  fJetPtCentPbPbRaw->Fill(rjet->Pt(),1); // for rec pp or inclusive PbPb
	  fJetPtCentPbPbCorr->Fill(pTbs,1); // for rec pp or inclusive PbPb
	  fMaxTrackPtInNTXh[1]->Fill(fMaxTrackPtInNTX,rjet->Pt());
	  if(fIsPossibleToSubstBckg) // if it was possible to calculate a perpendicular cone
	    {
	      fMinTrackPtInNTXRecalc->Fill(fMinTrackPtInNTXR,rjet->Pt(),1); // for rec pp or inclusive PbPb
	      fMaxTrackPtInNTXRecalc->Fill(fMaxTrackPtInNTXR,rjet->Pt());
	    }
	  AliDebug(4,"After filling the histograms for this jet of min and max pT in NTX \n");
	  for(Int_t v=0 ; v<rJetCounter; v++)
	    {
	      if(v<7)
		fJetsMultPtRD->Fill(jetPtsR[v],rJetCounter);
	    }
	  irecj++;
	} //reco jets
      nRecJets=irecj;
      fNAccJetsRD->Fill(rJetCounter,rJetCounter);
      // reference multiplicity stuff in pp, also filled in PbPb, but does not matter.
      if(v0CorrMult<25)
	fNAccJetsRDMult[0]->Fill(rJetCounter);
      if(v0CorrMult>=25&&v0CorrMult<50)
	fNAccJetsRDMult[1]->Fill(rJetCounter);
      if(v0CorrMult>=50&&v0CorrMult<90)
	fNAccJetsRDMult[2]->Fill(rJetCounter);
      if(v0CorrMult>=90&&v0CorrMult<120)
	fNAccJetsRDMult[3]->Fill(rJetCounter);
      if(v0CorrMult>=120&&v0CorrMult<150)
	fNAccJetsRDMult[4]->Fill(rJetCounter);
      if(v0CorrMult>=150&&v0CorrMult<200)
	fNAccJetsRDMult[5]->Fill(rJetCounter);
      if(v0CorrMult>=200&&v0CorrMult<300)
	fNAccJetsRDMult[6]->Fill(rJetCounter);
      if(v0CorrMult>=300)
	fNAccJetsRDMult[7]->Fill(rJetCounter);

      // reference multiplicity from tracks in TPC minus jet related tracks
      if(refNJMult<5) // &&refNJMult>1
	fNAccJetsRDMultSE[0]->Fill(rJetCounter);
      if(refNJMult>=5&&refNJMult<10)
	fNAccJetsRDMultSE[1]->Fill(rJetCounter);
      if(refNJMult>=10&&refNJMult<15)
	fNAccJetsRDMultSE[2]->Fill(rJetCounter);
      if(refNJMult>=15&&refNJMult<20)
	fNAccJetsRDMultSE[3]->Fill(rJetCounter);
      if(refNJMult>=20&&refNJMult<30)
	fNAccJetsRDMultSE[4]->Fill(rJetCounter);
      if(refNJMult>=30&&refNJMult<40)
	fNAccJetsRDMultSE[5]->Fill(rJetCounter);
      if(refNJMult>=40&&refNJMult<50)
	fNAccJetsRDMultSE[6]->Fill(rJetCounter);
      if(refNJMult>=50)
	fNAccJetsRDMultSE[7]->Fill(rJetCounter);

      AliDebug(4,"Checking if this data contains MC, in order to relate the jets \n");
      if(fUseAODMC)
	{
	  AliDebug(4,"Relating MC jets with reconstructed jets"); 
	  // Relate the jets
	  Int_t fCJDebug = 0; //debug level for getclosest jet
	  nGenJets = TMath::Min(nGenJets,maxJetNum);  
	  nRecJets = TMath::Min(nRecJets,maxJetNum);
	  Int_t iGenIndex[maxJetNum];    // Index of the generated jet for i-th rec -1 if none
	  Int_t iRecIndex[maxJetNum];    // Index of the reco jet for i-th gen -1 if none
	  for(int i = 0;i<maxJetNum;++i)
	    {
	      iGenIndex[i] = iRecIndex[i] = -1;
	    }
	  
	  AliAnalysisHelperJetTasks::GetClosestJets(genJets,nGenJets,recJets,nRecJets,
						    iGenIndex,iRecIndex,fCJDebug);
	  if(fCJDebug > 10)
	    AliDebug(4,Form("%s:%d",(char*)__FILE__,__LINE__)); 
	  if(fCJDebug > 3)
	    {
	      for(int i = 0;i<maxJetNum;++i)
		{
		  if(iGenIndex[i]>=0)
		    {
		      AliDebug(4,Form("iGenFound: %d -> %d",i,iGenIndex[i])); 
		      //  para el i-esimo jet reconstruido corresponde el jet iGenIndex[i]
		      AliDebug(4,Form("El jet reconstruido numero %d tiene sabor %d",i, genJetsFlavor[iGenIndex[i]])); 
		    }
		  if(iRecIndex[i]>=0)
		    {
		      AliDebug(4,Form("iRecFound: %d -> %d",i,iRecIndex[i]));  
		      //  para el i-esimo jet generado corresponde el jet iRecIndex[i]
		    }
		}
	    }
	  AliDebug(4,"Helper part finished"); 
	  
	  // Llenar los histogramas para los jets reconstruidos
	  
	  Int_t crf = 0;      // current reco jet flavor
	  Int_t crt = 0;      // current reco jet tracks
	  Double_t crpt = 0;  // current jet pt
	  Double_t cscm = 0.0;// current second central moment
	  
	  for(Int_t ixr=0; ixr<maxJetNum; ixr++)
	    {
	      AliDebug(4,Form("Processing jet number:%i",ixr)); 
	      if(iGenIndex[ixr]>=0)
		{
		  crf = genJetsFlavor[iGenIndex[ixr]]; //para el reco jet con indice ixr
		  crt = nrectracks[ixr];  // se necesitaron este numero de tracks
		  crpt = ptrecjet[ixr];
		  cscm = scmr[ixr];
		  
		  //Fill candidates histos
		  if(crt>=7) //gluon candidate
		    fFragCandidates[0]->Add(fHistContainerR4[ixr]);	      
		  if(crt<=4) //quark candidate
		    fFragCandidates[1]->Add(fHistContainerR4[ixr]);
		  
		  switch(abs(crf))
		    {
		    case 1:
		      fNChTr[6]->Fill(crt,crpt);
		      fHistPtParton[6]->Fill(crpt);
		      fSCM[6]->Fill(cscm,crpt);
		      fFragChargedR4[0]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[0]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[0]->Add(fHistContainerR2[ixr]);
		      break;
		    case 2:
		      fNChTr[7]->Fill(crt,crpt);
		      fHistPtParton[7]->Fill(crpt);
		      fSCM[7]->Fill(cscm,crpt);
		      fFragChargedR4[1]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[1]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[1]->Add(fHistContainerR2[ixr]);
		      break;
		    case 3:
		      fNChTr[8]->Fill(crt,crpt);
		      fHistPtParton[8]->Fill(crpt);
		      fSCM[8]->Fill(cscm,crpt);
		      fFragChargedR4[2]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[2]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[2]->Add(fHistContainerR2[ixr]);
		      break;
		    case 4:
		      fNChTr[9]->Fill(crt,crpt);
		      fHistPtParton[9]->Fill(crpt);
		      fSCM[9]->Fill(cscm,crpt);
		      fFragChargedR4[3]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[3]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[3]->Add(fHistContainerR2[ixr]);
		      break;
		    case 5:
		      fNChTr[10]->Fill(crt,crpt);
		      fHistPtParton[10]->Fill(crpt);
		      fSCM[10]->Fill(cscm,crpt);
		      fFragChargedR4[4]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[4]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[4]->Add(fHistContainerR2[ixr]);
		      break;
		    case 21:
		      fNChTr[11]->Fill(crt,crpt);
		      fHistPtParton[11]->Fill(crpt);
		      fSCM[11]->Fill(cscm,crpt);
		      fFragChargedR4[5]->Add(fHistContainerR4[ixr]);
		      fFragChargedR3[5]->Add(fHistContainerR3[ixr]);
		      fFragChargedR2[5]->Add(fHistContainerR2[ixr]);
		      break;	      
		    default:
		      break;
		    }  // end switch
		  AliDebug(4,Form("Sabor del reco jet con pt:%f y numero:%d es: %d y se necesitaron %d tracks \n",crpt,ixr,crf,crt));  
		} // end index condition
	    }  // end for cycle correlation gen-reco flavor
	} // end if MC info
    } // end of only MC info in the AOD, no reco jets
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskPartonDisc::Terminate(const Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetMCEventType(AliMCEvent *mcEvent) 
{
  //
  // Get the event type from the pythia headers
  //

  Int_t processNumber = 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);

  if(!pythiaGenHeader)
    {
      AliDebug(4,Form(" %s:%d No Pythia header!",(char*)__FILE__,__LINE__));  
      return 0;
    }

  processNumber = pythiaGenHeader->ProcessType();

  return processNumber;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetPhojetEventType(AliMCEvent *mcEvent) 
{
  //
  // Get the event type from the phojet header
  //

  Int_t processNumber = 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenDPMjetEventHeader* phojetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);

  if(!phojetGenHeader)
    {
      AliDebug(4,Form(" %s:%d No Phojet header!",(char*)__FILE__,__LINE__));  
      return 0;
    }

  processNumber = phojetGenHeader->ProcessType();

  return processNumber;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::IsInsideAcceptance(AliAODJet *jet)
{
  //
  // Check if the jet is inside abs(eta)<=fJetAcceptance
  //

  Double_t jeteta = jet->Eta();
  if(TMath::Abs(jeteta)<=fJetAcceptance)
    return kTRUE;
  else 
    return kFALSE;

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetJetFlavour(AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray)
{
  //
  // Get the jet flavour, using the definition:
  // The flavour will be that of the parton with the highest energy
  // within an angular distance <= 0.3
  // This method also keeps track of the mother of the parton giving the flavor
  //

  Int_t flavour = 0;
  Int_t pdgCode;
  Int_t aPDGcode;
  Double_t currentEnergy = 0;
  Double_t maxEnergy = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t tracketa = 0;
  Double_t trackphi = 0;
  Double_t flavrad = fFlavorRadius;   //radius used for flavor
  Int_t indexM = 0;
  UInt_t status=0; //status code of leading parton  
  fMpdg=0; // Initialize before each selection
  AliAODMCParticle *moftrack=0;
      for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
	  AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(iTracks);
	  if(!mctrack) continue;
	  pdgCode = mctrack->GetPdgCode();
	  aPDGcode = abs(pdgCode);  
	  if(aPDGcode==1||aPDGcode==2||aPDGcode==3||aPDGcode==4||aPDGcode==5||aPDGcode==6||aPDGcode==9||aPDGcode==21)
	    {
	      tracketa = mctrack->Eta();
	      trackphi = mctrack->Phi();
	      if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=flavrad)  
		{
		  currentEnergy = mctrack->E();
		  if(currentEnergy>maxEnergy)
		    {
		      maxEnergy = currentEnergy;
		      flavour = pdgCode;
		      indexM = mctrack->GetMother();
		      if(fCheckMCStatus) // default is true
			status =  mctrack->GetStatus();
		      //		      fFlavProc->Fill(flavour,status);
		      //testing
		      if(indexM<=0)
			fMpdg = 0; // Unknown
		      else
			{
			  moftrack = (AliAODMCParticle*) mcarray->At(indexM);
			  fMpdg = moftrack->GetPdgCode();
			}
		    } 
		}
	    }
	}

      // //Restriction to pythia string fragmentation
      // if(!fPhojetMC)  //if pythia
      // 	{
      // 	  if(fMpdg!=0) // if not from the string
      // 	    flavour=0;
      // 	}

      //  PDG code of the mother of the leading parton, leading parton, jet pT
      fPDGMothLPart->Fill(fMpdg,flavour,jet->Pt());
      fFlavProc->Fill(flavour,status);

  return flavour;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::GetDeltaR(Double_t eta1, Double_t phi1,Double_t eta2, Double_t phi2)
{
  //
  // Return R between the two jets or particles
  //
  
  Double_t deltaphi = TMath::Abs(phi2-phi1);
  if (deltaphi > TMath::Pi()) 
    deltaphi = 2.0 * TMath::Pi() - deltaphi;


  Double_t deltaR = sqrt((eta2-eta1)*(eta2-eta1)+deltaphi*deltaphi);
  return deltaR;

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetNumberOfMcChargedTracks(Int_t percentage,AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr)
{
  //
  // Calculate the number of charged particles necessary to
  // add the given percentage of the jet energy (transverse energy)
  // for the MC case

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;
  Int_t arraysize = 1000;
  Bool_t rfTrkFlag  = kFALSE;
  fCurrentJetCharge=0;

  AliDebug(4,Form("Eta of the jet:%f ",jeteta));  
  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) // jetpT==0
    return 0;

  // RefTracks
  Int_t trkinjet = jet->GetRefTracks()->GetEntriesFast();
  if(trkinjet!=0)
    rfTrkFlag = kTRUE;
  AliDebug(4,Form("Number of tracks in this mc jet by RefTracks:%i \n",trkinjet));  

  AllocateStaticContainer(arraysize);
  InitializeStaticContainer(arraysize);

  if(!rfTrkFlag)  // if not track ref, check track by track
    {
      AliDebug(4,Form("Empty Track Refs (mc)!"));  
      for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
	  AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(iTracks);
	  if(!mctrack) continue;
	  tracketa = mctrack->Eta();
	  trackphi = mctrack->Phi();
	  if(mctrack->Pt()<fMinpTVal) continue; // pT cut, not using track refs	      	  
	  if(mctrack->IsPhysicalPrimary())
	    {
	      if(mctrack->Charge()!=0&&mctrack->Charge()!=-99)
		{
		  if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
		    {	       
		      currentNumber++;
		      fCurrentJetCharge=fCurrentJetCharge+mctrack->Charge(); //add the charge of this track		 
		      fgContainer[currentNumber-1] = mctrack->Pt();  // save the current pt in the container
		    } // end if inside jet
		} // end charged
	    } // end physical primary
	}   // end for tracks
    } // end rfTrkFlag
  
  if(rfTrkFlag)  // if track ref, use them
    {
      AliDebug(4,Form("Using Track Refs (mc)!")); 
      for(Int_t ixt=0; ixt<trkinjet;ixt++)
	{
	  AliAODMCParticle *vtrack = dynamic_cast<AliAODMCParticle*>(jet->GetRefTracks()->At(ixt));
	  if(!vtrack) continue;
	  if(vtrack->Charge()!=0&&vtrack->Charge()!=-99)
	    {
	      currentNumber++;
	      fCurrentJetCharge=fCurrentJetCharge+vtrack->Charge(); //add the charge of this track		 
	      fgContainer[currentNumber-1] = vtrack->Pt();  // save the current pt in the container
	    } 
	} // end trk in jet
    } // end if trk ref

  // sort the contents of the container
  SortArray(fgContainer,arraysize);
  // loop over the contents and count how many tracks are necessary to recover the percentage of the energy
  numberOfChargedTracks = TracksForPercentage(fgContainer, arraysize, percentage, jetpT);
  AliDebug(4,Form("Number of tracks was:%i, returning",numberOfChargedTracks));  
  return numberOfChargedTracks;

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetNumberOfChargedTracks(Int_t percentage,AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr)
{
  //
  // Calculate the number of charged particles necessary to
  // add the given percentage of the jet energy (transverse energy)
  // for the AOD track case

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;
  Int_t arraysize = 1000;
  Bool_t rfTrkFlag  = kFALSE;
  fCurrentJetCharge=0;

  if(IsEqualRel(jetpT, 0.0)) // IsEqualRel(jetpT, 0.0) // jetpT==0
    return 0;

  // RefTracks
  Int_t trkinjet = jet->GetRefTracks()->GetEntriesFast();
  if(trkinjet!=0&&!fForceNotTR)
    rfTrkFlag = kTRUE;
  AliDebug(4,Form("Number of tracks in this reco jet by RefTracks:%i \n",trkinjet));  

  AllocateStaticContainer(arraysize);
  InitializeStaticContainer(arraysize);
  if(!rfTrkFlag)  // if not track ref, check track by track
    {
      AliDebug(4,Form("Empty Track Refs (reco)!"));  
      for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
	  AliAODTrack *aodtrack = aode->GetTrack(iTracks);
	  if(!aodtrack) continue;
	  tracketa = aodtrack->Eta();
	  trackphi = aodtrack->Phi();
	  if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
	    {
	      if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
	      if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries
	      if(aodtrack->Pt()<fMinpTVal) continue; // pT cut, not using track refs	      
	      currentNumber++;
	      fCurrentJetCharge=fCurrentJetCharge+aodtrack->Charge();
	      fgContainer[currentNumber-1] = aodtrack->Pt();  // save the current pt in the container

	      ////////start centrality dependent pT spec ////////
	      if(!fIsHIevent) //if is a proton proton event
		{
		  if((jet->Pt()>10.)&&(jet->Pt()<20.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,1);
		  if((jet->Pt()>20.)&&(jet->Pt()<30.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,1);
		  if((jet->Pt()>30.)&&(jet->Pt()<40.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,1);
		  if((jet->Pt()>40.)&&(jet->Pt()<50.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,1);
		  if((jet->Pt()>50.)&&(jet->Pt()<60.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,1);
		  if((jet->Pt()>60.)&&(jet->Pt()<70.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,1);
		  if((jet->Pt()>70.)&&(jet->Pt()<80.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,1);
		  if((jet->Pt()>80.)&&(jet->Pt()<90.))
		    fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,1);
		} // end of pp event
	      if(fIsHIevent) //if is a PbPb event     
		{
		  if(fEventCent>=0&&fEventCent<10.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,2);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,2);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,2);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,2);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,2);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,2);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,2);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,2);
		    } // end of 0-10
		  if(fEventCent>=10&&fEventCent<20.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,3);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,3);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,3);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,3);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,3);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,3);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,3);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,3);
		    } // end of 10-20
		  if(fEventCent>=20&&fEventCent<30.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,4);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,4);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,4);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,4);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,4);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,4);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,4);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,4);
		    } // end of 20-30
		  if(fEventCent>=30&&fEventCent<40.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,5);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,5);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,5);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,5);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,5);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,5);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,5);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,5);
		    } // end of 30-40
		  if(fEventCent>=40&&fEventCent<50.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,6);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,6);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,6);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,6);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,6);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,6);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,6);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,6);
		    }  // end of 40-50
		  if(fEventCent>=50&&fEventCent<60.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,7);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,7);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,7);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,7);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,7);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,7);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,7);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,7);
		    }  // end of 50-60
		  if(fEventCent>=60&&fEventCent<70.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,8);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,8);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,8);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,8);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,8);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,8);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,8);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,8);
		    }  // end of 60-70
		  if(fEventCent>=70&&fEventCent<80.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,9);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,9);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,9);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,9);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,9);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,9);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,9);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,9);
		    }  // end of 70-80
		  if(fEventCent>=80&&fEventCent<100.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),1,10);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),2,10);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),3,10);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),4,10);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),5,10);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),6,10);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),7,10);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(aodtrack->Pt(),8,10);
		    }  // end of 80-100
		}  //end of PbPb event
	      ////////end centrality dependent pT spec/////
	    } // end if inside jet	
	} // end for tracks
    } // end of no track ref
  
  if(rfTrkFlag)  // if track ref, use them
    {
      AliDebug(4,Form("Using Track Refs (reco)!")); 
      for(Int_t ixt=0; ixt<trkinjet;ixt++)
	{
	  AliVParticle *vtrack = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(ixt));
	  if(!vtrack) continue;
	  // No further checks, all cuts should be in in the track refs
	  if(vtrack->Charge()!=0&&vtrack->Charge()!=-99)
	    {
	      currentNumber++;
	      fCurrentJetCharge=fCurrentJetCharge+vtrack->Charge();
	      fgContainer[currentNumber-1] = vtrack->Pt();  // save the current pt in the container

	      ////////start centrality dependent pT spec ////////
	      if(!fIsHIevent) //if is a proton proton event
		{
		  if((jet->Pt()>10.)&&(jet->Pt()<20.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,1);
		  if((jet->Pt()>20.)&&(jet->Pt()<30.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,1);
		  if((jet->Pt()>30.)&&(jet->Pt()<40.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,1);
		  if((jet->Pt()>40.)&&(jet->Pt()<50.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,1);
		  if((jet->Pt()>50.)&&(jet->Pt()<60.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,1);
		  if((jet->Pt()>60.)&&(jet->Pt()<70.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,1);
		  if((jet->Pt()>70.)&&(jet->Pt()<80.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,1);
		  if((jet->Pt()>80.)&&(jet->Pt()<90.))
		    fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,1);
		} // end of pp event
	      if(fIsHIevent) //if is a PbPb event     
		{
		  if(fEventCent>=0&&fEventCent<10.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,2);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,2);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,2);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,2);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,2);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,2);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,2);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,2);
		    } // end of 0-10
		  if(fEventCent>=10&&fEventCent<20.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,3);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,3);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,3);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,3);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,3);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,3);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,3);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,3);
		    } // end of 10-20
		  if(fEventCent>=20&&fEventCent<30.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,4);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,4);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,4);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,4);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,4);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,4);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,4);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,4);
		    } // end of 20-30
		  if(fEventCent>=30&&fEventCent<40.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,5);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,5);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,5);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,5);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,5);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,5);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,5);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,5);
		    } // end of 30-40
		  if(fEventCent>=40&&fEventCent<50.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,6);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,6);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,6);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,6);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,6);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,6);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,6);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,6);
		    }  // end of 40-50
		  if(fEventCent>=50&&fEventCent<60.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,7);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,7);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,7);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,7);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,7);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,7);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,7);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,7);
		    }  // end of 50-60
		  if(fEventCent>=60&&fEventCent<70.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,8);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,8);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,8);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,8);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,8);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,8);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,8);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,8);
		    }  // end of 60-70
		  if(fEventCent>=70&&fEventCent<80.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,9);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,9);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,9);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,9);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,9);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,9);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,9);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,9);
		    }  // end of 70-80
		  if(fEventCent>=80&&fEventCent<100.)
		    {
		      if((jet->Pt()>10.)&&(jet->Pt()<20.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),1,10);
		      if((jet->Pt()>20.)&&(jet->Pt()<30.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),2,10);
		      if((jet->Pt()>30.)&&(jet->Pt()<40.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),3,10);
		      if((jet->Pt()>40.)&&(jet->Pt()<50.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),4,10);
		      if((jet->Pt()>50.)&&(jet->Pt()<60.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),5,10);
		      if((jet->Pt()>60.)&&(jet->Pt()<70.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),6,10);
		      if((jet->Pt()>70.)&&(jet->Pt()<80.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),7,10);
		      if((jet->Pt()>80.)&&(jet->Pt()<90.))
			fPtDistInJetConeRaw->Fill(vtrack->Pt(),8,10);
		    }  // end of 80-100
		}  //end of PbPb event
	      ////////end centrality dependent pT spec/////
	    } 
	} // end trk in jet
    } // end if trk ref
  

  // sort the contents of the container
  SortArray(fgContainer,arraysize);
  // loop over the contents and count how many tracks are necessary to recover the percetage of the energy
  numberOfChargedTracks = TracksForPercentage(fgContainer, arraysize, percentage, jetpT);
  return numberOfChargedTracks;
}
//________________________________________________________________________
void AliAnalysisTaskPartonDisc::AllocateStaticContainer(Int_t size)
{
  //
  // Allocate the static container with the given dimensions
  //

  if(fgContainer) 
    return;
  fgContainer = new Double_t[size];

}
//________________________________________________________________________
void AliAnalysisTaskPartonDisc::InitializeStaticContainer(Int_t size)
{
  //
  // Initialize the static container with the given dimensions
  //
    
  memset(fgContainer,0,size*sizeof(Double_t));

}
//________________________________________________________________________
void AliAnalysisTaskPartonDisc::SortArray(Double_t *pointer, Int_t arraySize)
{
  //
  // Sort the contents of the array
  // From lower to higher value
  //

  std::sort(pointer,pointer+arraySize);

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::TracksForPercentage(Double_t *array, Int_t arraysize, Int_t percentage, Double_t jetenergy)
{
  //
  // Loop over the contents and count how many tracks are necessary to recover 
  // the given percentage of the energy.
  // If all tracks did not sum the required fraction, it returns 500
  //

  AliDebug(4,Form("Calculating the number of tracks for a jet with energy:%f \n",jetenergy));
  Double_t ptsum = 0;
  Double_t threshold = jetenergy*percentage/100;
  Int_t tracknummer=0;
  fCurrentJetMinPtNT90 = 7000.; //dummy value for debugging
  for(Int_t inverse=arraysize; inverse>0; inverse--)
    {
      ptsum= ptsum + array[inverse-1];
      if(inverse==arraysize) //if the highest value
	fMaxTrackPtInNTX=array[inverse-1]; //saving the highest pT track value
      tracknummer++;
      fMinTrackPtInNTX=array[inverse-1]; // this is the current lowest pT track used
      if(fIsPossibleToSubstBckg) //Store the value if it was possible to find a perpendicular cone
	fCurrentJetMinPtNT90=array[inverse-1]; 
      if(ptsum>=threshold)  // the threshold was reached
	break;
      if((inverse==1)&&(ptsum<threshold)) //if it was not possible to reach the threshold
	{
	  tracknummer = 500; //dummy value for debugging
	  fCurrentJetMinPtNT90 = 7000.; //dummy value for debugging
	}
    }

  AliDebug(4,"Done calculating the number of tracks, returning to main code \n");
  return tracknummer;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::IsMCTrackInsideThisJet(AliAODMCParticle *MCParticle, AliAODJet *Jet, Double_t jr)
{
  //
  // Return kTrue if the mc track is inside the area covered by the cone of the jet
  //

  Double_t etapart = MCParticle->Eta();
  Double_t phipart = MCParticle->Phi();
  Double_t etajet = Jet->Eta();
  Double_t phijet = Jet->Phi(); 
  Double_t deltaeta = etajet-etapart;
  Double_t deltaphi = TMath::Abs(phijet-phipart);
  if (deltaphi > TMath::Pi()) 
    deltaphi = 2.0 * TMath::Pi() - deltaphi;

  Double_t deltar = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
  if(deltar<=jr)
    return kTRUE;
  else 
    return kFALSE;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::IsTrackInsideThisJet(AliAODTrack *aodT, AliAODJet *Jet, Double_t jr)
{
  //
  // Return kTrue if the track is inside the area covered by the cone of the jet
  //

  Double_t etapart = aodT->Eta();
  Double_t phipart = aodT->Phi();
  Double_t etajet = Jet->Eta();
  Double_t phijet = Jet->Phi(); 
  Double_t deltaeta = etajet-etapart;
  Double_t deltaphi = TMath::Abs(phijet-phipart);
  if (deltaphi > TMath::Pi()) 
    deltaphi = 2.0 * TMath::Pi() - deltaphi;

  Double_t deltar = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
  if(deltar<=jr)
    return kTRUE;
  else 
    return kFALSE;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::VertexInJet(AliAODVertex *pvtx, AliAODVertex *vtx, AliAODJet *jet, Double_t jr)
{
  //
  // Return kTRUE if the cone covers the vector 
  // from the primary vertex to this vertex
  //

  Double_t pvx =  pvtx->GetX(); //primary vertex x
  Double_t pvy =  pvtx->GetY(); //primary vertex y
  Double_t pvz =  pvtx->GetZ(); //primary vertex z

  Double_t vx =  vtx->GetX(); // vertex x
  Double_t vy =  vtx->GetY(); // vertex y
  Double_t vz =  vtx->GetZ(); // vertex z

  if(IsEqualRel(pvx, vx) && IsEqualRel(pvy, vy) && IsEqualRel(pvz, vz)) //!IsEqualRel(totalTrackPt, 0.0) // pvx==vx && pvy==vy && pvz==vz
    return kFALSE;

  Double_t thetaval = GetThetaAngle(vx-pvx,vy-pvy,vz-pvz);
  Double_t etaval = GetEtaValue(thetaval);
  Double_t phival = GetPhiAngle(vx-pvx,vy-pvy);

  Double_t etajet = jet->Eta();
  Double_t phijet = jet->Phi();

  if(GetDeltaR(etajet,phijet,etaval,phival)<=jr)
    return kTRUE;

  return kFALSE;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::GetEtaValue(Double_t theta) const
{
  //
  // Get the eta value
  // 

  Double_t eta = -TMath::Log(theta/2);
  return eta;

}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::GetThetaAngle(Double_t xval, Double_t yval, Double_t zval)
{
  //
  // Get the theta angle related to these coordinates
  // 

  if(IsEqualRel(zval, 0.0)) //!IsEqualRel(totalTrackPt, 0.0) // zval==0
    return TMath::PiOver2();

  Double_t theta = 0;
  Double_t erre = TMath::Sqrt(xval*xval+yval*yval+zval*zval);

  if(zval>0)
    theta = TMath::ACos(zval/erre);  

  if(zval<0)
    theta = TMath::Pi() - TMath::ACos(TMath::Abs(zval)/erre);

  if(IsEqualRel(theta, 0.0)) // theta==0
    AliError("ERROR in GetThetaAngle!");
 
  return theta;

}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::GetPhiAngle(Double_t xval, Double_t yval)
{
  //
  // Get the phi angle related to these coordinates
  // 

  if(IsEqualRel(xval, 0.0)) //IsEqualRel(zval, 0.0) // xval==0
    {
      if(yval>0)
	return TMath::PiOver2();
      if(yval<0)
	return (3/2*TMath::Pi());
    }
  Double_t phi = 0;

  if(xval>0)
    {
      if(yval>0)
	phi = TMath::ATan(yval/xval);
      if(yval<0)
	phi = 2*TMath::Pi()- TMath::ATan(TMath::Abs(yval)/xval); 
    }
  if(xval<0)
    {
      if(yval>0)
	phi = TMath::Pi() - TMath::ATan(yval/TMath::Abs(xval));
      if(yval<0)
	phi = TMath::Pi() + TMath::ATan(TMath::Abs(yval/xval));
    }

  return phi;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaPhiMC(AliAODJet *jet, AliAODMCParticle *particle)
{
  //
  // Get delta-phi MC jet-track
  // 

  Double_t deltaphi = jet->Phi()-particle->Phi();
  if (deltaphi > TMath::Pi()) 
    deltaphi = 2.0 * TMath::Pi() - deltaphi;
  if (deltaphi < -TMath::Pi()) 
    deltaphi =  -deltaphi- 2.0 * TMath::Pi();

  return deltaphi;

}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaEtaMC(AliAODJet *jet, AliAODMCParticle *particle)
{
  //
  // Get delta-eta MC jet-track
  // 

  Double_t deltaetaMC = jet->Eta() - particle->Eta();
  return deltaetaMC;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaPhiSqMC(AliAODJet *jet, AliAODMCParticle *particle)
{
  //
  // Get delta-phi^2 MC jet-track
  // 

  Double_t deltaphi = DeltaPhiMC(jet,particle);
  Double_t deltaphiSqMC = deltaphi*deltaphi;
  return deltaphiSqMC;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaEtaSqMC(AliAODJet *jet, AliAODMCParticle *particle)
{
  //
  // Get delta-eta^2 MC jet-track
  // 

  Double_t deltaeta = DeltaEtaMC(jet,particle);
  Double_t deltaetaSqMC = deltaeta*deltaeta;
  return deltaetaSqMC;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaPhiTrack(AliAODJet *jet, AliAODTrack *track)
{
  //
  // Get delta-phi Track jet-track
  // 

  // Double_t deltaphiTrack = jet->Phi() - track->Phi();
  // return deltaphiTrack;

  Double_t deltaphi = jet->Phi() - track->Phi();
  if (deltaphi > TMath::Pi()) 
    deltaphi = 2.0 * TMath::Pi() - deltaphi;
  if (deltaphi < -TMath::Pi()) 
    deltaphi =  -deltaphi- 2.0 * TMath::Pi();

  return deltaphi;


}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaEtaTrack(AliAODJet *jet, AliAODTrack *track)
{
  //
  // Get delta-eta Track jet-track
  // 

  Double_t deltaetaTrack = jet->Eta() - track->Eta();
  return deltaetaTrack;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaPhiSqTrack(AliAODJet *jet, AliAODTrack *track)
{
  //
  // Get delta-phi^2 Track jet-track
  // 

  Double_t deltaphi = DeltaPhiTrack(jet,track);
  Double_t deltaphiSqTrack = deltaphi*deltaphi;
  return deltaphiSqTrack;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::DeltaEtaSqTrack(AliAODJet *jet, AliAODTrack *track)
{
  //
  // Get delta-eta^2 Track jet-track
  // 

  Double_t deltaeta = DeltaEtaTrack(jet,track);
  Double_t deltaetaSqTrack = deltaeta*deltaeta;
  return deltaetaSqTrack;
}
//_________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::NumberOfReadEventsAOD(const char* currFile, Int_t &fNEvents)
{
  //
  // get the number of events read out to create this AOD
  // from the njets distribution from UA1 
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // code from AliAnalysisTaskJetSpectrum2 

  TString file(currFile);  
  fNEvents = 1;

  if(file.Contains("root_archive.zip#"))
    {
      Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
      Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
      file.Replace(pos+1,20,"");
    }
  else 
    {
      // not an archive take the basename....
      file.ReplaceAll(gSystem->BaseName(file.Data()),"");
    }
  Printf("%s",file.Data());  

  TFile *fnev = TFile::Open(Form("%s%s",file.Data(),"PWG4_JetTasksOutput.root")); 
  if(!fnev)
    {
      return kFALSE;
    } // no PWG4_JetTasksOutput.root
  else 
    {
      TList *list; 
      gDirectory->GetObject("PWG4_jethist_aodmc_ua104/jethist_aodmc_ua104;1",list);
      if(!list)
	{
	  fnev->Close();
	  return kFALSE;
	}
      fNEvents = ((TH1*)list->FindObject("NJetsH"))->GetEntries();
      fnev->Close();
    }
  return kTRUE;
}
//___________________________________________________________________________
void AliAnalysisTaskPartonDisc::HasOverlapedCones(TClonesArray *JetArray)
{
  //
  // Check if there are jet cones that overlap on the current event
  // for UA1 and SISCone, based on the cone axis and the radius.
  // There can be maximum 7.85 jet in the phi acceptance, with no overlap in
  // the phi direction (2pi/0.8).
  // Plus there can be up to two jets in contained in the eta acceptance,
  // if they are centered in -0.5 and 0.5 respectively, per phi interval
  // In total there can be up to 15.7 jets
  // Check up to 16 jets, inside the acceptance,
  // set the flags for up to those 8 jets.
  //

  // Now also check if there is a perpendicular area to the leading jet
  // that does not contain a jet, if so, set a flag to kTRUE
  // Possibility to remove single track jets if there are track references

  fJetEvent=kFALSE;
  fHasPerpCone=kTRUE;
  fEtaPerpCoord=0.0;
  fPhiPerpCoord=0.0;
  fPtPerpCoord=0.0;

  ResetJetFlags(); // reset the flags

  Int_t njets = JetArray->GetEntries();
  Int_t maxNjets = 16;
  Double_t etaCoordinates[16]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  Double_t phiCoordinates[16]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  for (Int_t idxj = 0; idxj < njets; idxj++) 
    {
      if(idxj>15)
	continue;  //maximum 16 jets

      AliAODJet *currentJet = dynamic_cast<AliAODJet*>(JetArray->At(idxj));
      if (!currentJet) 
	{
	  AliDebug(2,Form("ERROR: Could not receive jet %d \n", idxj)); 
	  continue;
	}
      // Store the axis coordinates of all sixteen jets
      etaCoordinates[idxj]=currentJet->Eta();
      phiCoordinates[idxj]=currentJet->Phi();
    }  //end of first cycle

  // Ahora checar que no tengan overlap
  Double_t deltaeta = 0.0;
  Double_t deltaphi = 0.0;
  Double_t deltar = 0.0;
  Int_t currentIdxJet = 0;
  Int_t trkinjet = 0;

  for (Int_t inxoj = 0; inxoj < maxNjets; inxoj++) 
    {
      // only check up to the number of found jets in the event
      if(inxoj==njets)
	break;

      // Optional RefTracks check for single track jets
      if(fForceSkipSJ)
	{
	  AliAODJet *currentJet = dynamic_cast<AliAODJet*>(JetArray->At(inxoj));
	  if (!currentJet) 
	    {
	      AliDebug(2,Form("ERROR: Could not receive jet %d \n", inxoj)); 
	      continue;
	    }
	  trkinjet = currentJet->GetRefTracks()->GetEntriesFast();
	  if(trkinjet==1) // si tiene 1 solo track se marca como malo
	    {
	      fJetFlags[inxoj]=kFALSE;
	      continue;
	    }
	} // end of skip of single track jets condition

      // First check if the current jet has its axis inside acceptance
      if(!(TMath::Abs(etaCoordinates[inxoj])<=fJetAcceptance)) //
	{
	  fJetFlags[inxoj]=kFALSE;
	  continue;
	}
      currentIdxJet = inxoj;
      //Check this jet with the rest of the jets
      for (Int_t idx2 = 0; idx2 < maxNjets; idx2++) 
	{
	  if(idx2==njets) // just check the number of found jets
	    break;

	  if(currentIdxJet==idx2) // if the same jet
	    continue;
	  if(!fJetFlags[idx2]) // if the other jet is already not usable
	    continue;
	  if(!(TMath::Abs(etaCoordinates[idx2])<=fJetAcceptance)) // if the jet is outside acceptance
	    continue;
	  deltaeta = etaCoordinates[currentIdxJet]-etaCoordinates[idx2];
	  deltaphi = phiCoordinates[currentIdxJet]-phiCoordinates[idx2];
	  if (deltaphi > TMath::Pi()) 
	    deltaphi = 2.0 * TMath::Pi() - deltaphi;
	  deltar = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
	  if(deltar<=2.0*fJetRadius) // if the distance between jet axis is less than 2r with any jet, mark as not usable
	    {
	      fJetFlags[currentIdxJet]=kFALSE;
	      fJetFlags[idx2]=kFALSE; // both jets are not usable
	      // if(fEnablePrints)
	      // 	{
	      // 	  printf("Rejecting jet #:%i because it had overlap with jet #:%i \n",currentIdxJet, idx2);
	      // 	  printf("Eta1:%f, Phi1:%f \n", etaCoordinates[currentIdxJet],phiCoordinates[currentIdxJet]);
	      // 	  printf("Eta2:%f, Phi2:%f \n", etaCoordinates[idx2],phiCoordinates[idx2]);
	      // 	  printf("Current event number:%i \n",fEvtCount-1);
	      // 	}
	    } // end of actual flagging
	} // end of checking the current jet with the rest
    } // end of loop over the jets of the branch

  // Check if there is at least one accepted jet, so it makes sense to calculate the perpendicular cone
  Int_t accJets = 0;
  for (Int_t checkjets = 0; checkjets < maxNjets; checkjets++) 
    {
      // only check up to the number of found jets in the event
      if(checkjets==njets)
	break;
      if(fJetFlags[checkjets]) // find the accepted leading jet in acceptance
	accJets++;
    }
  if(accJets>0)
    fJetEvent=kTRUE;

  // Check for the leading jet on the event
  for (Int_t searchlead = 0; searchlead < maxNjets; searchlead++) 
    {
      // only check up to the number of found jets in the event
      if(searchlead==njets)
	break;
      if(fJetFlags[searchlead]) // find the accepted leading jet in acceptance
	{
	  // Phi + pi/2
	  fPhiPerpCoord = phiCoordinates[searchlead] + 0.5*TMath::Pi();
	  if(fPhiPerpCoord>2.0*TMath::Pi())
	    fPhiPerpCoord = fPhiPerpCoord - 2.0*TMath::Pi();
	  // Same eta
	  fEtaPerpCoord = etaCoordinates[searchlead];
	  // Now check if this cone overlaps with any found jet
	  for (Int_t jets = 0; jets < maxNjets; jets++) 
	    {
	      // only check up to the number of found jets in the event
	      if(jets==njets)
		break;
	      // now check that this jet is not the same as the leading
	      if(jets==searchlead)
		continue;
	      
	      deltaphi = phiCoordinates[jets]-fPhiPerpCoord;
	      if (deltaphi > TMath::Pi()) 
		deltaphi = 2.0 * TMath::Pi() - deltaphi;
	      if(deltaphi<=2.0*fJetRadius) // if the distance between cone axis is less than 2r with any jet, mark as not usable
		{
		  fHasPerpCone=kFALSE;
		}
		// }
	    } // loop over accepted jets
	  break; // done doing stuff with the leading
	} // if for the first accepted jet (leading)
    }
}
//_______________________________________________________________________
void AliAnalysisTaskPartonDisc::ResetJetFlags()
{
  //
  // Reset the flags used for tagging jets from the branches
  // Use before calling HasOverlapedCones(jetbranch)
  //

  for(Int_t a=0; a<16;a++)
    {
      fJetFlags[a]=kTRUE;
    }
  
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetNMcChargedTracksAboveThreshold(AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr)
{  
  // Calculate the number of charged particles above threshold
  // inside this jet for the MC case
  // the threshold is fCurrentJetMinPtNT90

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(zval, 0.0) //jetpT==0
    return 0;    

  if(IsEqualRel(fCurrentJetMinPtNT90, 7000.)) // fCurrentJetMinPtNT90==7000.
    return 1000; //dummy val for debugging

  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
    {
      AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(iTracks);
      if(!mctrack) continue;
      tracketa = mctrack->Eta();
      trackphi = mctrack->Phi();
      if(mctrack->Pt()<fCurrentJetMinPtNT90) continue; 	      	  
      if(mctrack->IsPhysicalPrimary())
	{
	  if(mctrack->Charge()!=0&&mctrack->Charge()!=-99)
	    {
	      if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
		{	       
		  currentNumber++;		 
		} // end if inside jet
	    } // end charged
	} // end physical primary
    }   // end for tracks  
  
  numberOfChargedTracks = currentNumber;
  AliDebug(4,Form("Number of tracks above threshold MC was:%i, returning",numberOfChargedTracks));  
  return numberOfChargedTracks;

}
//________________________________________________________________________
  Int_t AliAnalysisTaskPartonDisc::GetRecalcNTXMc(Int_t percentage, AliAODJet *originaljet, Int_t ntracks, TClonesArray *mcarray, Double_t jr)
{
  //
  // Calculate the number of charged particles necessary to
  // add the given percentage of the jet energy (transverse energy)
  // after background substraction, for the MC case

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jetpT = fBckgSbsJet[0]; //pT
  Double_t jeteta = fBckgSbsJet[1]; //eta 
  Double_t jetphi = fBckgSbsJet[2]; //phi
  Double_t tracketa = 0;
  Double_t trackphi = 0;
  Int_t arraysize = 1000;
  Bool_t rfTrkFlag  = kFALSE;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) //jetpT==0
    return 0;

  // RefTracks
  Int_t trkinjet = originaljet->GetRefTracks()->GetEntriesFast();
  if(trkinjet!=0)
    rfTrkFlag = kTRUE;
  AliDebug(4,Form("Number of tracks in this mc jet by RefTracks:%i \n",trkinjet));  

  AllocateStaticContainer(arraysize);
  InitializeStaticContainer(arraysize);

  if(!rfTrkFlag)  // if not track ref, check track by track
    {
      AliDebug(4,Form("Empty Track Refs (mc)!"));  
      for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
	  AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(iTracks);
	  if(!mctrack) continue;
	  tracketa = mctrack->Eta();
	  trackphi = mctrack->Phi();
	  if(mctrack->Pt()<fMinpTVal) continue; // pT cut, not using track refs	      	  
	  if(mctrack->IsPhysicalPrimary())
	    {
	      if(mctrack->Charge()!=0&&mctrack->Charge()!=-99)
		{
		  if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
		    {	       
		      currentNumber++;		 
		      fgContainer[currentNumber-1] = mctrack->Pt();  // save the current pt in the container
		    } // end if inside jet
		} // end charged
	    } // end physical primary
	}   // end for tracks
    } // end rfTrkFlag
  
  if(rfTrkFlag)  // if track ref, use them
    {
      AliDebug(4,Form("Using Track Refs (mc)!")); 
      for(Int_t ixt=0; ixt<trkinjet;ixt++)
	{
	  AliAODMCParticle *vtrack = dynamic_cast<AliAODMCParticle*>(originaljet->GetRefTracks()->At(ixt));
	  if(!vtrack) continue;
	  if(vtrack->Charge()!=0&&vtrack->Charge()!=-99)
	    {
	      currentNumber++;
	      fgContainer[currentNumber-1] = vtrack->Pt();  // save the current pt in the container
	    } 
	} // end trk in jet
    } // end if trk ref

  // sort the contents of the container
  SortArray(fgContainer,arraysize);
  // loop over the contents and count how many tracks are necessary to recover the percentage of the energy
  numberOfChargedTracks = TracksForPercentageRecalc(fgContainer, arraysize, percentage, jetpT); //este es el que tengo que modificar...
  AliDebug(4,Form("Number of tracks was:%i, returning",numberOfChargedTracks));  
  return numberOfChargedTracks;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::TracksForPercentageRecalc(Double_t *array, Int_t arraysize, Int_t percentage, Double_t jetenergy)
{
  //
  // Loop over the contents and count how many tracks are necessary to recover 
  // the given percentage of the energy.
  // If all tracks did not sum the required fraction, it returns 500
  // this saves the minimum pT used, during the nt90 recalculation
  // after jet energy correction

  AliDebug(4,Form("Re-calculating the number of tracks for a jet with corrected pT:%f \n",jetenergy));
  Double_t ptsum = 0;
  Double_t threshold = jetenergy*percentage/100;
  Int_t tracknummer=0;
  fCurrentJetMinPtNT90Recalc=7000.; //dummy value for debugging
  for(Int_t inverse=arraysize; inverse>0; inverse--)
    {
      ptsum= ptsum + array[inverse-1];
      if(inverse==arraysize) //if the highest value
	fMaxTrackPtInNTXR=array[inverse-1]; //saving the highest pT track value in the recalculation
      tracknummer++;
      fMinTrackPtInNTXR=array[inverse-1]; // this is the current lowest pT track used in the recalculation
      if(fIsPossibleToSubstBckg) //Store the value for the current jet, during recalculation
	fCurrentJetMinPtNT90Recalc=array[inverse-1]; 
      if(ptsum>=threshold)  // the threshold was reached
	break;
      if((inverse==1)&&(ptsum<threshold)) //if it was not possible to reach the threshold
	{
	  tracknummer = 500;
	  fCurrentJetMinPtNT90Recalc=7000.; //dummy values for debugging
	}
    }

  AliDebug(4,"Done re-calculating the number of tracks, returning to main code \n");
  return tracknummer;

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetRecalcNMcChTrUpThr(AliAODJet *jet, Int_t ntracks, TClonesArray *mcarray, Double_t jr)
{  
  // Calculate the number of charged particles above threshold
  // inside this jet for the MC case
  // the threshold is fCurrentJetMinPtNT90Recalc

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) //jetpT==0
    return 0;    

  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
    {
      AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(iTracks);
      if(!mctrack) continue;
      tracketa = mctrack->Eta();
      trackphi = mctrack->Phi();
      if(mctrack->Pt()<fCurrentJetMinPtNT90Recalc) continue; 	      	  
      if(mctrack->IsPhysicalPrimary())
	{
	  if(mctrack->Charge()!=0&&mctrack->Charge()!=-99)
	    {
	      if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
		{	       
		  currentNumber++;		 
		} // end if inside jet
	    } // end charged
	} // end physical primary
    }   // end for tracks  
  
  numberOfChargedTracks = currentNumber;
  AliDebug(4,Form("Recalculated number of tracks above threshold MC was:%i, returning",numberOfChargedTracks));  
  return numberOfChargedTracks;

}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetNRecChargedTracksAboveThreshold(AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr)
{
  //
  // Calculate the number of charged particles 
  // above the threshold set by the NTX calculation
  // for the AOD track case, the threshold is fCurrentJetMinPtNT90
  // the fCurrentJetMinPtNT90 was set when calling the NTX method before

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) // jetpT==0
      return 0;

  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
    {
      AliAODTrack *aodtrack = aode->GetTrack(iTracks);
      if(!aodtrack) continue;
      tracketa = aodtrack->Eta();
      trackphi = aodtrack->Phi();
      if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
	{
	  if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
	  if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries
	  if(aodtrack->Pt()<fCurrentJetMinPtNT90) continue;	      
	  currentNumber++;		 
	} // end if inside jet	
    } // end for tracks
  
  numberOfChargedTracks = currentNumber;
  return numberOfChargedTracks;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetRecalcNTXRec(Int_t percentage,AliAODJet *originaljet, Int_t ntracks, AliAODEvent *aode, Double_t jr)
{
  //
  // Calculate the number of charged particles necessary to
  // add the given percentage of the jet energy (transverse energy)
  // after pT recalculation, for the AOD track case

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jetpT = fBckgSbsJet[0]; //pT
  Double_t jeteta = fBckgSbsJet[1]; //eta 
  Double_t jetphi = fBckgSbsJet[2]; //phi
  Double_t tracketa = 0;
  Double_t trackphi = 0;
  Int_t arraysize = 1000;
  Bool_t rfTrkFlag  = kFALSE;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) //jetpT==0
    return 0;

  // RefTracks
  Int_t trkinjet = originaljet->GetRefTracks()->GetEntriesFast();
  if(trkinjet!=0&&!fForceNotTR)
    rfTrkFlag = kTRUE;
  AliDebug(4,Form("Number of tracks in this reco jet by RefTracks:%i \n",trkinjet));  

  AllocateStaticContainer(arraysize);
  InitializeStaticContainer(arraysize);
  if(!rfTrkFlag)  // if not track ref, check track by track
    {
      AliDebug(4,Form("Empty Track Refs (reco)!"));  
      for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
	  AliAODTrack *aodtrack = aode->GetTrack(iTracks);
	  if(!aodtrack) continue;
	  tracketa = aodtrack->Eta();
	  trackphi = aodtrack->Phi();
	  if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
	    {
	      if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
	      if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries
	      if(aodtrack->Pt()<fMinpTVal) continue; // pT cut, not using track refs	      
	      currentNumber++;		 
	      fgContainer[currentNumber-1] = aodtrack->Pt();  // save the current pt in the container
	    } // end if inside jet	
	} // end for tracks
    } // end of no track ref
  
  if(rfTrkFlag)  // if track ref, use them
    {
      AliDebug(4,Form("Using Track Refs (reco)!")); 
      for(Int_t ixt=0; ixt<trkinjet;ixt++)
	{
	  AliVParticle *vtrack = dynamic_cast<AliVParticle*>(originaljet->GetRefTracks()->At(ixt));
	  if(!vtrack) continue;
	  // No further checks, all cuts should be in in the track refs
	  if(vtrack->Charge()!=0&&vtrack->Charge()!=-99)
	    {
	      currentNumber++;
	      fgContainer[currentNumber-1] = vtrack->Pt();  // save the current pt in the container
	    } 
	} // end trk in jet
    } // end if trk ref
  

  // sort the contents of the container
  SortArray(fgContainer,arraysize);
  // loop over the contents and count how many tracks are necessary to recover the percetage of the energy
  numberOfChargedTracks = TracksForPercentageRecalc(fgContainer, arraysize, percentage, jetpT);
  return numberOfChargedTracks;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetRecalcNRecChTrUpThr(AliAODJet *jet, Int_t ntracks, AliAODEvent *aode, Double_t jr)
{
  //
  // Calculate the number of charged particles 
  // above the threshold set by the NTX calculation
  // for the AOD track case, the threshold is fCurrentJetMinPtNT90Recalc
  // the fCurrentJetMinPtNT90Recalc was set when calling the NTX method before

  Int_t numberOfChargedTracks = 0;
  Int_t currentNumber = 0;
  Double_t jeteta = jet->Eta();
  Double_t jetphi = jet->Phi();
  Double_t jetpT = jet->Pt();
  Double_t tracketa = 0;
  Double_t trackphi = 0;

  if(IsEqualRel(jetpT, 0.0)) //IsEqualRel(jetpT, 0.0) // jetpT==0
      return 0;

  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
    {
      AliAODTrack *aodtrack = aode->GetTrack(iTracks);
      if(!aodtrack) continue;
      tracketa = aodtrack->Eta();
      trackphi = aodtrack->Phi();
      if(GetDeltaR(jeteta, jetphi, tracketa, trackphi)<=jr)
	{
	  if(!aodtrack->TestFilterBit(fFilterBit)) continue; //track filter selection
	  if(!aodtrack->IsPrimaryCandidate()) continue; // only primaries
	  if(aodtrack->Pt()<fCurrentJetMinPtNT90Recalc) continue;	      
	  currentNumber++;		 
	} // end if inside jet	
    } // end for tracks
  
  numberOfChargedTracks = currentNumber;
  return numberOfChargedTracks;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::IsTrackInsideExcludedArea(Double_t tracketa, Double_t trackphi, TClonesArray *recojetsA)
{
  //
  // Check if these track coordinates are in an excluded area
  //

  //Primero: checar si esta dentro de cualquier jet, no importa si esta dentro de la
  //aceptancia o no(en base al eje), aun podria haber tracks del jet dentro del
  //evento.
  Int_t numbRecJ =  recojetsA->GetEntries();
  Double_t etaCurrJet = 0.0;
  Double_t phiCurrJet = 0.0;
  Double_t deltaeta = 0.0;
  Double_t deltaphi = 0.0;
  Double_t deltar = 0.0;
  Double_t extendedRadius = fJetRadius+fIncExcR; 
  Double_t extendedRadiusDiJet = 0.0;
  Double_t phiCurrJetExArL = 0.0; // preliminary left boundary of exluded phi of current jet
  Double_t phiCurrJetExArR = 0.0; // preliminary right boundary of exluded phi of current jet
  Double_t leftBoundary = 0.0;
  Double_t rightBoundary = 0.0;

  if(fNotExtDiJEx)
    extendedRadiusDiJet = fJetRadius;  // old behaviour
  if(!fNotExtDiJEx)
    extendedRadiusDiJet = fJetRadius+fIncExcR;  // new behaviour

  for (Int_t indxrec = 0; indxrec < numbRecJ; indxrec++) 
    {
      AliDebug(4,Form("Number of current jet:%i \n",indxrec));
      AliAODJet *rjet = dynamic_cast<AliAODJet*>(recojetsA->At(indxrec));
      if (!rjet) 
	{
	  AliDebug(2,Form("ERROR: Could not receive jet %d\n", indxrec)); 
	  continue;
	}
      etaCurrJet = rjet->Eta();
      phiCurrJet = rjet->Phi();
      deltaeta = etaCurrJet-tracketa;
      deltaphi = phiCurrJet-trackphi;
      if (deltaphi > TMath::Pi()) 
	deltaphi = 2.0 * TMath::Pi() - deltaphi;
      deltar = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
      if(deltar<=extendedRadius) // if the track is inside the jet, reject (extended radius)
	{
	  //	  printf("Excluding track for being in the extended jet area. Eta:%f, Phi:%f, deltar=%f \n",tracketa,trackphi,deltar);
	  return kTRUE;
	}
      // Now check if it is in the expected dijet area for the two hardest jets
      if(indxrec==0||indxrec==1) //two hardest jets in the event
	{
	  //left boundary
	  phiCurrJetExArL= phiCurrJet + TMath::Pi() - extendedRadiusDiJet;
	  if(phiCurrJetExArL>TMath::TwoPi())
	    phiCurrJetExArL=phiCurrJetExArL-TMath::TwoPi();
	  //right boundary
	  phiCurrJetExArR= phiCurrJet + TMath::Pi() + extendedRadiusDiJet;
	  if(phiCurrJetExArR>TMath::TwoPi())
	    phiCurrJetExArR=phiCurrJetExArR-TMath::TwoPi();
	  // //Assign left and right boundary
	  leftBoundary=phiCurrJetExArL;
	  rightBoundary=phiCurrJetExArR;
	  // now check if inside the excluded area
	  if(trackphi>=leftBoundary&&trackphi<=rightBoundary)
	    {
	      // printf("Excluding track for being in di-jet excluded area!!! \n");
	      // printf("phi of track:%f, left bound:%f, right bound:%f \n",trackphi,leftBoundary,rightBoundary);
	      return kTRUE;
	    }
	}
    }
  // Si sobrevive todos los tests, regresar kFALSE
  return kFALSE;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPartonDisc::GetV0ExcludedMultiplicity(TClonesArray *recojets)
{
  //
  // Calculate the V0 MIP multiplicity that should be excluded due
  // to possible jet activity
  //

  //Check if there is V0 info
  if(!fVZero)
    {
      AliDebug(2,"ERROR: There is not VZERO info in the AOD"); 
      //      printf("No hay objeto fVZero \n");
      return 0.0;
    }

  Int_t numbRecJ =  recojets->GetEntries();
  Double_t phiCurrJet = 0.0;
  //Check if there are jets
  if(numbRecJ==0)
    {
      AliDebug(2,"ERROR: There are not jets in the event"); 
      return 0.0;
    }

  //Flags for V0A sectors
  Bool_t removeS0cells = kFALSE;
  Bool_t removeS1cells = kFALSE;
  Bool_t removeS2cells = kFALSE;
  Bool_t removeS3cells = kFALSE;
  Bool_t removeS4cells = kFALSE;
  Bool_t removeS5cells = kFALSE;
  Bool_t removeS6cells = kFALSE;
  Bool_t removeS7cells = kFALSE;
  //Flags for V0C sectors
  Bool_t removeS0cellsV0C = kFALSE;
  Bool_t removeS1cellsV0C = kFALSE;
  Bool_t removeS2cellsV0C = kFALSE;
  Bool_t removeS3cellsV0C = kFALSE;
  Bool_t removeS4cellsV0C = kFALSE;
  Bool_t removeS5cellsV0C = kFALSE;
  Bool_t removeS6cellsV0C = kFALSE;
  Bool_t removeS7cellsV0C = kFALSE;
  //Excedent multiplicity in V0A, V0C and total
  Double_t excedentV0A = 0.0;
  Double_t excedentV0C = 0.0;
  Double_t excedentV0Total = 0.0;
  Double_t phiValLow = 0.0;
  Double_t phiValUp  = 0.0;
  Double_t extendedR = fJetRadius + fIncExcR;

  for (Int_t indxrec = 0; indxrec < numbRecJ; indxrec++) 
    {
      // Now check the two hardest jets
      if(indxrec==0||indxrec==1) //two hardest jets in the event
	{
	  AliAODJet *rjet = dynamic_cast<AliAODJet*>(recojets->At(indxrec));
	  if (!rjet) 
	    {
	      AliDebug(2,Form("ERROR: Could not receive jet %d\n", indxrec)); 
	      continue;
	    }
	  phiCurrJet = rjet->Phi();

	  //Get the excluded phi boundaries
	  Double_t b1 = phiCurrJet  + TMath::Pi() - extendedR;
	  Double_t b2  = phiCurrJet  + TMath::Pi() + extendedR;
	  if(b1>TMath::TwoPi())
	    b1 = b1 - TMath::TwoPi();
	  if(b2>TMath::TwoPi())
	    b2 = b2 - TMath::TwoPi();

	  phiValLow = b1;
	  phiValUp = b2; 

	  if(phiValLow>TMath::TwoPi()||phiValUp>TMath::TwoPi()||phiValLow<0.||phiValUp<0.)
	    return 500.0;

	  // First V0A (EN EVE LA NUMERACION EMPIEZA EN 1...)
	  // Cells in sector S0 ( 0    to  pi/4): 32,40,48,56   
	  // Cells in sector S1 ( pi/4 to  pi/2): 33,41,49,57   
	  // Cells in sector S2 ( pi/2 to 3pi/4): 34,42,50,58   
	  // Cells in sector S3 (3pi/4 to  pi  ): 35,43,51,59
	  // Cells in sector S4 ( pi   to 5pi/4): 36,44,52,60
	  // Cells in sector S5 (5pi/4 to 3pi/2): 37,45,53,61
	  // Cells in sector S6 (3pi/2 to 7pi/4): 38,46,54,62 
	  // Cells in sector S7 (7pi/4 to 2pi  ): 39,47,55,63
	  
	  //Check if falls in sector S0
	  if(phiValLow>=(0.)&&phiValLow<=(1./4.*TMath::Pi()))
	    removeS0cells = kTRUE;
	  if(phiValUp>=(0.)&&phiValUp<=(1./4.*TMath::Pi()))
	    removeS0cells = kTRUE;
	  //Check if falls in sector S1
	  if(phiValLow>=(1./4.*TMath::Pi())&&phiValLow<=(1./2.*TMath::Pi()))
	    removeS1cells = kTRUE;
	  if(phiValUp>=(1./4.*TMath::Pi())&&phiValUp<=(1./2.*TMath::Pi()))
	    removeS1cells = kTRUE;
	  //Check if falls in sector S2
	  if(phiValLow>=(1./2.*TMath::Pi())&&phiValLow<=(3./4.*TMath::Pi()))
	    removeS2cells = kTRUE;
	  if(phiValUp>=(1./2.*TMath::Pi())&&phiValUp<=(3./4.*TMath::Pi()))
	    removeS2cells = kTRUE;
	  //Check if falls in sector S3
	  if(phiValLow>=(3./4.*TMath::Pi())&&phiValLow<=(TMath::Pi()))
	    removeS3cells = kTRUE;
	  if(phiValUp>=(3./4.*TMath::Pi())&&phiValUp<=(TMath::Pi()))
	    removeS3cells = kTRUE;
	  //Check if falls in sector S4
	  if(phiValLow>=(TMath::Pi())&&phiValLow<=(5./4.*TMath::Pi()))
	    removeS4cells = kTRUE;
	  if(phiValUp>=(TMath::Pi())&&phiValUp<=(5./4.*TMath::Pi()))
	    removeS4cells = kTRUE;
	  //Check if falls in sector S5
	  if(phiValLow>=(5./4.*TMath::Pi())&&phiValLow<=(3./2.*TMath::Pi()))
	    removeS5cells = kTRUE;
	  if(phiValUp>=(5./4.*TMath::Pi())&&phiValUp<=(3./2.*TMath::Pi()))
	    removeS5cells = kTRUE;
	  //Check if falls in sector S6
	  if(phiValLow>=(3./2.*TMath::Pi())&&phiValLow<=(7./4.*TMath::Pi()))
	    removeS6cells = kTRUE;
	  if(phiValUp>=(3./2.*TMath::Pi())&&phiValUp<=(7./4.*TMath::Pi()))
	    removeS6cells = kTRUE;
	  //Check if falls in sector S7
	  if(phiValLow>=(7./4.*TMath::Pi())&&phiValLow<=(TMath::TwoPi()))
	    removeS7cells = kTRUE;
	  if(phiValUp>=(7./4.*TMath::Pi())&&phiValUp<=(TMath::TwoPi()))
	    removeS7cells = kTRUE;
	  
	  /////////////////////////////////////////////////////////////////////////////////
	  
	  // Now V0C (EN EVE LA NUMERACION EMPIEZA EN 1...)
	  // Cells in sector S0 ( pi/2 to 3pi/4): 0,8,16,24 
	  // Cells in sector S1 (3pi/4 to  pi  ): 1,9,17,25
	  // Cells in sector S2 ( pi   to 5pi/4): 2,10,18,26
	  // Cells in sector S3 (5pi/4 to 3pi/2): 3,11,19,27
	  // Cells in sector S4 (3pi/2 to 7pi/4): 4,12,20,28
	  // Cells in sector S5 (7pi/4 to 2pi  ): 5,13,21,29
	  // Cells in sector S6 ( 0    to  pi/4): 6,14,22,30  
	  // Cells in sector S7 ( pi/4 to  pi/2): 7,15,23,31  
	  
	  //Check if falls in sector S0
	  if(phiValLow>=(1./2.*TMath::Pi())&&phiValLow<=(3./4.*TMath::Pi()))
	    removeS0cellsV0C = kTRUE;
	  if(phiValUp>=(1./2.*TMath::Pi())&&phiValUp<=(3./4.*TMath::Pi()))
	    removeS0cellsV0C = kTRUE;
	  //Check if falls in sector S1
	  if(phiValLow>=(3./4.*TMath::Pi())&&phiValLow<=(TMath::Pi()))
	    removeS1cellsV0C = kTRUE;
	  if(phiValUp>=(3./4.*TMath::Pi())&&phiValUp<=(TMath::Pi()))
	    removeS1cellsV0C = kTRUE;
	  //Check if falls in sector S2
	  if(phiValLow>=(TMath::Pi())&&phiValLow<=(5./4.*TMath::Pi()))
	    removeS2cellsV0C = kTRUE;
	  if(phiValUp>=(TMath::Pi())&&phiValUp<=(5./4.*TMath::Pi()))
	    removeS2cellsV0C = kTRUE;
	  //Check if falls in sector S3
	  if(phiValLow>=(5./4.*TMath::Pi())&&phiValLow<=(3./2.*TMath::Pi()))
	    removeS3cellsV0C = kTRUE;
	  if(phiValUp>=(5./4.*TMath::Pi())&&phiValUp<=(3./2.*TMath::Pi()))
	    removeS3cellsV0C = kTRUE;
	  //Check if falls in sector S4
	  if(phiValLow>=(3./2.*TMath::Pi())&&phiValLow<=(7./4.*TMath::Pi()))
	    removeS4cellsV0C = kTRUE;
	  if(phiValUp>=(3./2.*TMath::Pi())&&phiValUp<=(7./4.*TMath::Pi()))
	    removeS4cellsV0C = kTRUE;
	  //Check if falls in sector S5
	  if(phiValLow>=(7./4.*TMath::Pi())&&phiValLow<=(TMath::TwoPi()))
	    removeS5cellsV0C = kTRUE;
	  if(phiValUp>=(7./4.*TMath::Pi())&&phiValUp<=(TMath::TwoPi()))
	    removeS5cellsV0C = kTRUE;
	  //Check if falls in sector S6
	  if(phiValLow>=(0.)&&phiValLow<=(1./4.*TMath::Pi()))
	    removeS6cellsV0C = kTRUE;
	  if(phiValUp>=(0.)&&phiValUp<=(1./4.*TMath::Pi()))
	    removeS6cellsV0C = kTRUE;
	  //Check if falls in sector S7
	  if(phiValLow>=(1./4.*TMath::Pi())&&phiValLow<=(1./2.*TMath::Pi()))
	    removeS7cellsV0C = kTRUE;
	  if(phiValUp>=(1./4.*TMath::Pi())&&phiValUp<=(1./2.*TMath::Pi()))
	    removeS7cellsV0C = kTRUE;

	  //	  printf("phi del jet:%f, philow:%f, phiup:%f \n",phiCurrJet,phiValLow,phiValUp);
	} // end if leading jets (2)
    } // end jet loop

  // printf("_________V0A____________\n");
  // printf("Status sector S0:%i \n",removeS0cells);
  // printf("Status sector S1:%i \n",removeS1cells);
  // printf("Status sector S2:%i \n",removeS2cells);
  // printf("Status sector S3:%i \n",removeS3cells);
  // printf("Status sector S4:%i \n",removeS4cells);
  // printf("Status sector S5:%i \n",removeS5cells);
  // printf("Status sector S6:%i \n",removeS6cells);
  // printf("Status sector S7:%i \n",removeS7cells);
  // printf("_______________________\n");
  	  
  if(removeS0cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(32)+fVZero->GetMultiplicity(40)+fVZero->GetMultiplicity(48)+fVZero->GetMultiplicity(56);
  if(removeS1cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(33)+fVZero->GetMultiplicity(41)+fVZero->GetMultiplicity(49)+fVZero->GetMultiplicity(57);
  if(removeS2cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(34)+fVZero->GetMultiplicity(42)+fVZero->GetMultiplicity(50)+fVZero->GetMultiplicity(58);
  if(removeS3cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(35)+fVZero->GetMultiplicity(43)+fVZero->GetMultiplicity(51)+fVZero->GetMultiplicity(59);
  if(removeS4cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(36)+fVZero->GetMultiplicity(44)+fVZero->GetMultiplicity(52)+fVZero->GetMultiplicity(60);
  if(removeS5cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(37)+fVZero->GetMultiplicity(45)+fVZero->GetMultiplicity(53)+fVZero->GetMultiplicity(61);
  if(removeS6cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(38)+fVZero->GetMultiplicity(46)+fVZero->GetMultiplicity(54)+fVZero->GetMultiplicity(62);
  if(removeS7cells)
    excedentV0A = excedentV0A + fVZero->GetMultiplicity(39)+fVZero->GetMultiplicity(47)+fVZero->GetMultiplicity(55)+fVZero->GetMultiplicity(63);
  
  // printf("________V0C____________\n");
  // printf("Status sector S0:%i \n",removeS0cellsV0C);
  // printf("Status sector S1:%i \n",removeS1cellsV0C);
  // printf("Status sector S2:%i \n",removeS2cellsV0C);
  // printf("Status sector S3:%i \n",removeS3cellsV0C);
  // printf("Status sector S4:%i \n",removeS4cellsV0C);
  // printf("Status sector S5:%i \n",removeS5cellsV0C);
  // printf("Status sector S6:%i \n",removeS6cellsV0C);
  // printf("Status sector S7:%i \n",removeS7cellsV0C);
  // printf("_______________________\n");
  
  if(removeS0cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(0)+fVZero->GetMultiplicity(8)+fVZero->GetMultiplicity(16)+fVZero->GetMultiplicity(24);
  if(removeS1cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(1)+fVZero->GetMultiplicity(9)+fVZero->GetMultiplicity(17)+fVZero->GetMultiplicity(25);
  if(removeS2cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(2)+fVZero->GetMultiplicity(10)+fVZero->GetMultiplicity(18)+fVZero->GetMultiplicity(26);
  if(removeS3cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(3)+fVZero->GetMultiplicity(11)+fVZero->GetMultiplicity(19)+fVZero->GetMultiplicity(27);
  if(removeS4cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(4)+fVZero->GetMultiplicity(12)+fVZero->GetMultiplicity(20)+fVZero->GetMultiplicity(28);
  if(removeS5cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(5)+fVZero->GetMultiplicity(13)+fVZero->GetMultiplicity(21)+fVZero->GetMultiplicity(29);
  if(removeS6cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(6)+fVZero->GetMultiplicity(14)+fVZero->GetMultiplicity(22)+fVZero->GetMultiplicity(30);
  if(removeS7cellsV0C)
    excedentV0C = excedentV0C + fVZero->GetMultiplicity(7)+fVZero->GetMultiplicity(15)+fVZero->GetMultiplicity(23)+fVZero->GetMultiplicity(31);
  
  excedentV0Total = excedentV0A+excedentV0C;
	  
  //  printf("La multiplicidad total en V0 es:%f, la multiplicidad excedente en V0A es:%f, y en V0C es:%f, la multiplicidad corregida es:%f \n",fVZero->GetMTotV0A()+fVZero->GetMTotV0C(),excedentV0A,excedentV0C,fVZero->GetMTotV0A()+fVZero->GetMTotV0C()-excedentV0Total);

  return excedentV0Total;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPartonDisc::GetV0LikeExcludedMultMC(TClonesArray *mcjets, TClonesArray *mcparticles)
{
  //
  // Calculate the V0 like multiplicity that should be excluded due
  // to possible jet activity in MC events
  //

  Int_t nMcJets =  mcjets->GetEntries();
  Int_t tracksMC = mcparticles->GetEntriesFast();
  Double_t phiCurrJet = 0.0;
  //Check if there are jets
  if(nMcJets==0)
    {
      AliDebug(2,"ERROR: There are no MC jets in the event"); 
      return 0;
    }
  Int_t excludedMultV0Like = 0;
  Double_t trackphi=0.0;
  Double_t extendedR = fJetRadius + fIncExcR;
  for (Int_t indxMC = 0; indxMC < nMcJets; indxMC++) 
    {
      // Now check the two hardest jets
      if(indxMC==0||indxMC==1) //two hardest jets in the event
	{
	  AliAODJet *mcJet = dynamic_cast<AliAODJet*>(mcjets->At(indxMC));
	  if (!mcJet) 
	    {
	      AliDebug(2,Form("ERROR: Could not receive jet %d\n", indxMC)); 
	      continue;
	    }
	  phiCurrJet = mcJet->Phi();
	  //Get the excluded phi boundaries
	  Double_t b1 = phiCurrJet  + TMath::Pi() - extendedR;
	  Double_t b2  = phiCurrJet  + TMath::Pi() + extendedR;
	  if(b1>TMath::TwoPi())
	    b1 = b1 - TMath::TwoPi();
	  if(b2>TMath::TwoPi())
	    b2 = b2 - TMath::TwoPi();
	  // now check the charged tracks in the V0 acceptance
	  for(Int_t aodMCTrack = 0; aodMCTrack < tracksMC; aodMCTrack++ )
	    {
	      AliAODMCParticle *mctrackf = (AliAODMCParticle*) mcparticles->At(aodMCTrack);
	      if(!mctrackf) continue;
	      if(!mctrackf->IsPhysicalPrimary()) continue;
	      if(mctrackf->Charge()==0||mctrackf->Charge()==-99) continue;
	      if(mctrackf->Pt()<fMinpTValMC) continue; // cut off en MC	    
	      trackphi = mctrackf->Phi();
	      if(trackphi>=b1&&trackphi<=b2)
		{
		  //V0A
		  if(((mctrackf->Eta())>(2.8))&&((mctrackf->Eta())<(5.1)))
		    excludedMultV0Like++;
		  //V0C
		  if(((mctrackf->Eta())>(-3.7))&&((mctrackf->Eta())<(-1.7)))
		    excludedMultV0Like++;
		}
	    }
	} // end of 2 hardest jets
    } // end jet loop
  return excludedMultV0Like;
}
//________________________________________________________________________
void AliAnalysisTaskPartonDisc::FillPerpConeHisto(TH3F *currenthisto, Int_t ntracks, AliAODEvent *aode, Int_t CentralityBin, Int_t pTBin)
{

  // Fills the histrogram of the pT distribution in the perpendicular cone

  Double_t aodtracketaC = 0.;

  for(Int_t aodT = 0; aodT < ntracks; aodT++ )
    {
      AliAODTrack *aodtrackC = aode->GetTrack(aodT);
      if(!aodtrackC) continue;
      aodtracketaC = TMath::Abs(aodtrackC->Eta());
      if(aodtracketaC>0.9) continue;
      if(!aodtrackC->TestFilterBit(fFilterBit)) continue; //track filter selection
      if(!aodtrackC->IsPrimaryCandidate()) continue; // only primaries, maybe is redundant with the previous selection...
      if(fJetEvent) // if has an accepted jet, calculate the perpendicular cone
	{
	  if(HasPerpendicularCone()) // If there is a perpendicular cone available
	    {
	      if(aodtrackC->Pt()>fMinpTVal)
		{
		  if(GetDeltaR(fEtaPerpCoord,fPhiPerpCoord,aodtrackC->Eta(),aodtrackC->Phi())<fJetRadius)
		    {
		      currenthisto->Fill(aodtrackC->Pt(),pTBin,CentralityBin);
		    }
		}
	    }
	} // end if jet event
    }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPartonDisc::IsEqualRel(Double_t vA, Double_t vB)
{
  // Comparison of Double_t values

  Double_t epsVal = 0.000001;
  return TMath::Abs(vA-vB) <= epsVal*TMath::Abs(vA);

}
