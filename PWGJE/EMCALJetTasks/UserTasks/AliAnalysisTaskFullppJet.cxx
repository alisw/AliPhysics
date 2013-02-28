// **************************************
// Full pp jet task - ESD input only
// Extract the jet spectrum and all the 
// systematic uncertainties
// -R. Ma, Mar 2011
// **************************************

#include <TCanvas.h>
#include <TChain.h>
#include <TFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TTree.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskFullppJet.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliVParticle.h"
#include "AliInputEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "TGeoManager.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TGeoGlobalMagField.h"

#include "AliAODJet.h"
#include "AliFJWrapper.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

ClassImp(AliAnalysisTaskFullppJet)

const Float_t kRadius[3] = {0.4,0.2,0.3};
const Int_t kNRadii = 3;
const Double_t kPI = TMath::Pi();
const Double_t kdRCut[3] = {0.25,0.1,0.15};

//________________________________________________________________________
AliAnalysisTaskFullppJet::AliAnalysisTaskFullppJet() : 
  AliAnalysisTaskSE("default"), 
  fVerbosity(0), fEDSFileCounter(0), fNTracksPerChunk(0), fRejectPileup(kFALSE), fRejectExoticTrigger(kTRUE),
  fAnaType(0), fPeriod("lhc11a"), fESD(0), fAOD(0), fMC(0), fStack(0), fTrackArray(0x0), fClusterArray(0x0), fMcPartArray(0x0),
  fIsMC(kFALSE), fPhySelForMC(kFALSE), fChargedMC(kFALSE), fXsecScale(0.),
  fCentrality(99), fZVtxMax(10), 
  fTriggerType(-1), fCheckTriggerMask(kTRUE), fIsTPCOnlyVtx(kFALSE),
  fIsExoticEvent3GeV(kFALSE), fIsExoticEvent5GeV(kFALSE), fIsEventTriggerBit(kFALSE), fOfflineTrigger(kFALSE), fTriggerMask(0x0),
  fGeom(0x0), fRecoUtil(0x0),
  fEsdTrackCuts(0x0), fHybridTrackCuts1(0x0), fHybridTrackCuts2(0x0),fTrackCutsType(0), fKinCutType(0), fTrkEtaMax(0.9),
  fdEdxMin(73), fdEdxMax(90), fEoverPMin(0.8), fEoverPMax(1.2),
  fMatchType(0), fRejectExoticCluster(kTRUE), fRemoveBadChannel(kFALSE), fUseGoodSM(kFALSE),
  fStudySubEInHC(kFALSE), fStudyMcOverSubE(kFALSE),
  fElectronRejection(kFALSE), fHadronicCorrection(kFALSE), fFractionHC(0.), fHCLowerPtCutMIP(1e4), fClusterEResolution(0x0),
  fJetNEFMin(0.01), fJetNEFMax(0.98),
  fSpotGoodJet(kFALSE), fFindChargedOnlyJet(kFALSE), fFindNeutralOnlyJet(kFALSE),
  fCheckTrkEffCorr(kFALSE), fTrkEffCorrCutZ(0.3), fRandomGen(0x0), fRunUE(0), fCheckTPCOnlyVtx(0), fRunSecondaries(0),
  fSysJetTrigEff(kFALSE), fVaryJetTrigEff(0), 
  fSysTrkPtRes(kFALSE), fVaryTrkPtRes(0), fSysTrkEff(kFALSE), fVaryTrkEff(0.), 
  fSysTrkClsMth(kFALSE), fCutdEta(0.05), fCutdPhi(0.05),
  fSysNonLinearity(kFALSE), fNonLinear(0x0), fSysClusterEScale(kFALSE), fVaryClusterEScale(0), fSysClusterERes(kFALSE), fVaryClusterERes(0),
  fSysClusterizer(0),
  fNonStdBranch(""), fNonStdFile(""), fAlgorithm("aKt"), fRadius("0.4"), fRecombinationScheme(5),
  fConstrainChInEMCal(kFALSE), fRejectNK(kFALSE), fRejectWD(kFALSE), fSmearMC(kFALSE), fTrkPtResData(0x0),
  fOutputList(0x0), fSaveQAHistos(kTRUE), fhJetEventCount(0x0), fhJetEventStat(0x0), fhEventStatTPCVtx(0x0), fhChunkQA(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  // Output slot #0 id reserved by the base class for AOD
  for(Int_t i=0; i<2; i++)
    {
      fTrkPtMin[i] = 0.15;
      fTrkPtMax[i] = 200;
      fClsEtMin[i] = 0.15;
      fClsEtMax[i] = 200;

      fVertexGenZ[i]       = 0x0;
      fEventZ[i]           = 0x0;
      fhNTrials[i]         = 0x0;
      fhNMatchedTrack[i]   = 0x0;
      fhClsE[i]            = 0x0;
      for(Int_t k=0; k<2; k++)
	{
	  fhSysClusterE[i][k]     = 0x0;
	  fhSysNCellVsClsE[i][k]  = 0x0;
	}
    }

  for(Int_t r=0; r<kNRadii; r++)
    {
      for(Int_t j=0; j<5; j++)
	{
	  fHCOverSubE[r][j]       = 0x0;
	  fHCOverSubEFrac[r][j]   = 0x0;
	  for(Int_t k=0; k<2; k++)
	    {
	      fhSubClsEVsJetPt[k][r][j] = 0x0;
	      fhHCTrkPtClean[k][r][j]   = 0x0;
	      fhHCTrkPtAmbig[k][r][j]   = 0x0;
	    }

	  if(j<4) fhSubEVsTrkPt[r][j] = 0x0; 

	  if(j<3)
	    {
	      fJetCount[j][r]          = 0x0;
	      fhNeutralPtInJet[j][r]   = 0x0;
	      fhTrigNeuPtInJet[j][r]   = 0x0;
	      fhChargedPtInJet[j][r]   = 0x0;
	      fhLeadNePtInJet[j][r]    = 0x0;
	      fhLeadChPtInJet[j][r]    = 0x0;
	      fJetEnergyFraction[j][r] = 0x0;
	      fJetNPartFraction[j][r]  = 0x0;
	    }
	  if(j<2)
	    {
	      fRelTrkCon[j][r]         = 0x0;
	      fhFcrossVsZleading[j][r] = 0x0;
	      fhChLeadZVsJetPt[j][r]   = 0x0;
	      fhJetPtWithTrkThres[j][r]= 0x0;
	      fhJetPtWithClsThres[j][r]= 0x0;
	      for(Int_t k=0; k<2; k++)
		{
		  fhJetPtVsLowPtCons[j][r][k] = 0x0;
		}
	      fhJetPtInExoticEvent[j][r] = 0x0;
	      fhJetInTPCOnlyVtx[j][r]  = 0x0;
	      fhCorrTrkEffPtBin[j][r] = 0x0;
	      for(Int_t i=0; i<kNBins; i++)
		{
		  fhCorrTrkEffSample[j][r][i] = 0x0;
		}
	    }

	  if(r==0 && j<3)
	    {
	      for(Int_t k=0; k<2; k++)
		{
		  for(Int_t l=0; l<2; l++)
		    {
		      fhUEJetPtNorm[j][k][l] = 0x0;
		      fhUEJetPtVsSumPt[j][k][l] = 0x0;
		      fhUEJetPtVsConsPt[j][k][l] = 0x0;
		    }
		}
	    }
	}
      for(Int_t i=0; i<2; i++)
	{
	  fhNKFracVsJetPt[i][r]     = 0x0;
	  fhWeakFracVsJetPt[i][r]   = 0x0;
	  fhJetResponseNK[i][r]     = 0x0;
	  fhJetResponseWP[i][r]     = 0x0;
	  fhJetResolutionNK[i][r]   = 0x0;
	  fhJetResolutionWP[i][r]   = 0x0;
	  fhJetResponseNKSM[i][r]   = 0x0;
	  fhJetResponseWPSM[i][r]   = 0x0;
	  fhJetResolutionNKSM[i][r] = 0x0;
	  fhJetResolutionWPSM[i][r] = 0x0;
	}
    }

  for(Int_t s=0; s<3; s++)
    {
      for(Int_t a=0; a<2; a++)
	{
	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      fDetJetFinder[s][a][r] = 0x0;
	      fJetTCA[s][a][r] = 0x0;
	      if(s==0 && a==0)
		{
		  fTrueJetFinder[r] = 0x0;
		  fMcTruthAntikt[r] = 0x0;
		}
	    }
	}
    }

  for(Int_t j=0; j<3; j++)
    {
      fTrkEffFunc[j] = 0x0;
      fhSecondaryResponse[j] = 0x0;
    }
  for(Int_t i=0; i<10; i++)
    { 
      fTriggerCurve[i] = 0x0;
      fTriggerEfficiency[i] = 0x0; 
    }

}

//________________________________________________________________________
AliAnalysisTaskFullppJet::AliAnalysisTaskFullppJet(const char *name) : 
  AliAnalysisTaskSE(name), 
  fVerbosity(0), fEDSFileCounter(0), fNTracksPerChunk(0), fRejectPileup(kFALSE), fRejectExoticTrigger(kTRUE), 
  fAnaType(0), fPeriod("lhc11a"), fESD(0), fAOD(0), fMC(0), fStack(0), fTrackArray(0x0), fClusterArray(0x0), fMcPartArray(0x0),
  fIsMC(kFALSE), fPhySelForMC(kFALSE), fChargedMC(kFALSE), fXsecScale(0.),
  fCentrality(99), fZVtxMax(10), 
  fTriggerType(-1), fCheckTriggerMask(kTRUE), fIsTPCOnlyVtx(kFALSE),
  fIsExoticEvent3GeV(kFALSE), fIsExoticEvent5GeV(kFALSE), fIsEventTriggerBit(kFALSE), fOfflineTrigger(kFALSE), fTriggerMask(0x0),
  fGeom(0x0), fRecoUtil(0x0),
  fEsdTrackCuts(0x0), fHybridTrackCuts1(0x0), fHybridTrackCuts2(0x0), fTrackCutsType(0), fKinCutType(0), fTrkEtaMax(0.9),
  fdEdxMin(73), fdEdxMax(90), fEoverPMin(0.8), fEoverPMax(1.2),
  fMatchType(0), fRejectExoticCluster(kTRUE), fRemoveBadChannel(kFALSE), fUseGoodSM(kFALSE),
  fStudySubEInHC(kFALSE), fStudyMcOverSubE(kFALSE),
  fElectronRejection(kFALSE), fHadronicCorrection(kFALSE), fFractionHC(0.), fHCLowerPtCutMIP(1e4), fClusterEResolution(0x0),
  fJetNEFMin(0.01), fJetNEFMax(0.98),
  fSpotGoodJet(kFALSE), fFindChargedOnlyJet(kFALSE), fFindNeutralOnlyJet(kFALSE),
  fCheckTrkEffCorr(kFALSE), fTrkEffCorrCutZ(0.3), fRandomGen(0x0), fRunUE(0), fCheckTPCOnlyVtx(0), fRunSecondaries(0),
  fSysJetTrigEff(kFALSE), fVaryJetTrigEff(0), 
  fSysTrkPtRes(kFALSE), fVaryTrkPtRes(0), fSysTrkEff(kFALSE), fVaryTrkEff(0.), 
  fSysTrkClsMth(kFALSE), fCutdEta(0.05), fCutdPhi(0.05),
  fSysNonLinearity(kFALSE), fNonLinear(0x0), fSysClusterEScale(kFALSE), fVaryClusterEScale(0), fSysClusterERes(kFALSE), fVaryClusterERes(0),
  fSysClusterizer(0), 
  fNonStdBranch(""), fNonStdFile(""), fAlgorithm("aKt"), fRadius("0.4"), fRecombinationScheme(5),
  fConstrainChInEMCal(kFALSE), fRejectNK(kFALSE), fRejectWD(kFALSE), fSmearMC(kFALSE), fTrkPtResData(0x0),
  fOutputList(0x0), fSaveQAHistos(kTRUE), fhJetEventCount(0x0), fhJetEventStat(0x0), fhEventStatTPCVtx(0x0), fhChunkQA(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  // Output slot #0 id reserved by the base class for AOD
  for(Int_t i=0; i<2; i++)
    {
      fTrkPtMin[i] = 0.15;
      fTrkPtMax[i] = 200;
      fClsEtMin[i] = 0.15;
      fClsEtMax[i] = 200;

      fVertexGenZ[i]       = 0x0;
      fEventZ[i]           = 0x0;
      fhNTrials[i]         = 0x0;
      fhNMatchedTrack[i]   = 0x0;
      fhClsE[i]            = 0x0;
      for(Int_t k=0; k<2; k++)
	{
	  fhSysClusterE[i][k]     = 0x0;
	  fhSysNCellVsClsE[i][k]  = 0x0;
	}
    }

  for(Int_t r=0; r<kNRadii; r++)
    {
      for(Int_t j=0; j<5; j++)
	{
	  fHCOverSubE[r][j]       = 0x0;
	  fHCOverSubEFrac[r][j]   = 0x0;
	  for(Int_t k=0; k<2; k++)
	    {
	      fhSubClsEVsJetPt[k][r][j] = 0x0;
	      fhHCTrkPtClean[k][r][j]   = 0x0;
	      fhHCTrkPtAmbig[k][r][j]   = 0x0;
	    }

	  if(j<4) fhSubEVsTrkPt[r][j] = 0x0; 

	  if(j<3)
	    {
	      fJetCount[j][r]          = 0x0;
	      fhNeutralPtInJet[j][r]   = 0x0;
	      fhTrigNeuPtInJet[j][r]   = 0x0;
	      fhChargedPtInJet[j][r]   = 0x0;
	      fhLeadNePtInJet[j][r]    = 0x0;
	      fhLeadChPtInJet[j][r]    = 0x0;
	      fJetEnergyFraction[j][r] = 0x0;
	      fJetNPartFraction[j][r]  = 0x0;
	    }
	  if(j<2)
	    {
	      fRelTrkCon[j][r]         = 0x0;
	      fhFcrossVsZleading[j][r] = 0x0;
	      fhChLeadZVsJetPt[j][r]   = 0x0;
	      fhJetPtWithTrkThres[j][r]= 0x0;
	      fhJetPtWithClsThres[j][r]= 0x0;
	      for(Int_t k=0; k<2; k++)
		{
		  fhJetPtVsLowPtCons[j][r][k] = 0x0;
		}
	      fhJetPtInExoticEvent[j][r] = 0x0;
	      fhJetInTPCOnlyVtx[j][r]  = 0x0;
	      fhCorrTrkEffPtBin[j][r] = 0x0;
	      for(Int_t i=0; i<kNBins; i++)
		{
		  fhCorrTrkEffSample[j][r][i] = 0x0;
		}
	    }

	  if(r==0 && j<3)
	    {
	      for(Int_t k=0; k<2; k++)
		{
		  for(Int_t l=0; l<2; l++)
		    {
		      fhUEJetPtNorm[j][k][l] = 0x0;
		      fhUEJetPtVsSumPt[j][k][l] = 0x0;
		      fhUEJetPtVsConsPt[j][k][l] = 0x0;
		    }
		}
	    }
	}
      for(Int_t i=0; i<2; i++)
	{
	  fhNKFracVsJetPt[i][r]     = 0x0;
	  fhWeakFracVsJetPt[i][r]   = 0x0;
	  fhJetResponseNK[i][r]     = 0x0;
	  fhJetResponseWP[i][r]     = 0x0;
	  fhJetResolutionNK[i][r]   = 0x0;
	  fhJetResolutionWP[i][r]   = 0x0;
	  fhJetResponseNKSM[i][r]   = 0x0;
	  fhJetResponseWPSM[i][r]   = 0x0;
	  fhJetResolutionNKSM[i][r] = 0x0;
	  fhJetResolutionWPSM[i][r] = 0x0;
	}
    }

  for(Int_t s=0; s<3; s++)
    {
      for(Int_t a=0; a<2; a++)
	{
	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      fDetJetFinder[s][a][r] = 0x0;
	      fJetTCA[s][a][r] = 0x0;
	      if(s==0 && a==0)
		{
		  fTrueJetFinder[r] = 0x0;
		  fMcTruthAntikt[r] = 0x0;
		}
	    }
	}
    }

  for(Int_t j=0; j<3; j++)
    {
      fTrkEffFunc[j] = 0x0;
      fhSecondaryResponse[j] = 0x0;
    }
  for(Int_t i=0; i<10; i++)
    { 
      fTriggerCurve[i] = 0x0;
      fTriggerEfficiency[i] = 0x0; 
    }
}

//________________________________________________________________________
AliAnalysisTaskFullppJet::~AliAnalysisTaskFullppJet()
{
  //Destructor

  if(fEsdTrackCuts) delete fEsdTrackCuts;
  if(fHybridTrackCuts1) delete fHybridTrackCuts1;
  if(fHybridTrackCuts2) delete fHybridTrackCuts2;
  if(fOutputList)
    { fOutputList->Delete(); delete fOutputList;}
  for(Int_t s=0; s<3; s++)
    {
      for(Int_t a=0; a<2; a++)
	{
	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      if(fDetJetFinder[s][a][r]) delete fDetJetFinder[s][a][r];
	      if(fJetTCA[s][a][r]) { fJetTCA[s][a][r]->Delete(); delete fJetTCA[s][a][r]; }
	      if(s==0 && a==0)
		{
		  if(fTrueJetFinder[r]) delete fTrueJetFinder[r];
		  if(fMcTruthAntikt[r]) { fMcTruthAntikt[r]->Delete(); delete fMcTruthAntikt[r]; }
		}
	    }
	}
    }
  if(fRandomGen) delete fRandomGen;
  for(Int_t i=0; i<3; i++)
    {
      if(fTrkEffFunc[i]) delete fTrkEffFunc[i];
      if(fhSecondaryResponse[i]) delete fhSecondaryResponse[i];
    }
  for(Int_t i=0; i<10; i++)
    { 
      if(fTriggerEfficiency[i]) delete fTriggerEfficiency[i];
      if(fTriggerCurve[i]) delete fTriggerCurve[i];
    }
  for(Int_t r=0; r<kNRadii; r++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  if(fhCorrTrkEffPtBin[j][r]) delete fhCorrTrkEffPtBin[j][r];
	  for(Int_t i=0; i<kNBins; i++)
	    {
	      if(fhCorrTrkEffSample[j][r][i]) delete fhCorrTrkEffSample[j][r][i];
	    }
	}
    }
  if(fTrackArray)   { fTrackArray->Delete(); delete fTrackArray; }
  if(fClusterArray) { fClusterArray->Delete(); delete fClusterArray; }
  if(fMcPartArray)  { fMcPartArray->Reset(); delete fMcPartArray; }

  if(fRecoUtil) delete fRecoUtil;
  if(fClusterEResolution) delete fClusterEResolution;
  if(fNonLinear) delete fNonLinear;
  if(fTrkPtResData) delete fTrkPtResData;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::Notify()
{
  //
  // Fill the number of tracks per chunk
  // 

  fhChunkQA->SetBinContent(fEDSFileCounter,fNTracksPerChunk);
  fNTracksPerChunk = 0;
  fEDSFileCounter++;
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  //

  if(fRunUE) fFindChargedOnlyJet = kTRUE;
  
  const Int_t nTrkPtBins = 100;
  const Float_t lowTrkPtBin=0, upTrkPtBin=100.;

  const Int_t nbins = 220;
  Double_t xbins[221];
  for(Int_t i=0; i<nbins+1; i++)
    xbins[i] = i;

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fhJetEventStat = new TH1F("fhJetEventStat","Event statistics for jet analysis",12,0,12);
  fhJetEventStat->GetXaxis()->SetBinLabel(1,"ALL");
  fhJetEventStat->GetXaxis()->SetBinLabel(2,"MB");
  fhJetEventStat->GetXaxis()->SetBinLabel(3,"MB+vtx+10cm");
  fhJetEventStat->GetXaxis()->SetBinLabel(4,"EMC");
  fhJetEventStat->GetXaxis()->SetBinLabel(5,"EMC+vtx+10cm");
  fhJetEventStat->GetXaxis()->SetBinLabel(6,"MB+vtx");
  fhJetEventStat->GetXaxis()->SetBinLabel(7,"EMC+vtx");
  fhJetEventStat->GetXaxis()->SetBinLabel(8,"TriggerBit");
  fhJetEventStat->GetXaxis()->SetBinLabel(9,"LED");
  fhJetEventStat->GetXaxis()->SetBinLabel(10,"MB+TVtx+10cm");
  fhJetEventStat->GetXaxis()->SetBinLabel(11,"EMC+TVtx+10cm");
  fhJetEventStat->GetXaxis()->SetBinLabel(12,"ALL-Pileup");
  fOutputList->Add(fhJetEventStat);

  fhJetEventCount = new TH1F("fhJetEventCount","Event statistics for jet analysis",12,0,12); 
  fhJetEventCount->GetXaxis()->SetBinLabel(1,"ALL");
  fhJetEventCount->GetXaxis()->SetBinLabel(2,"MB");
  fhJetEventCount->GetXaxis()->SetBinLabel(3,"MB-pileup");
  fhJetEventCount->GetXaxis()->SetBinLabel(4,"MB+Vtx");
  fhJetEventCount->GetXaxis()->SetBinLabel(5,"MB+Vtx+10cm");
  fhJetEventCount->GetXaxis()->SetBinLabel(6,"EMC");
  fhJetEventCount->GetXaxis()->SetBinLabel(7,"EMC-pileup");
  fhJetEventCount->GetXaxis()->SetBinLabel(8,"EMC+Vtx");
  fhJetEventCount->GetXaxis()->SetBinLabel(9,"EMC+Vtx+10cm");
  fhJetEventCount->GetXaxis()->SetBinLabel(10,"Good EMC");
  fOutputList->Add(fhJetEventCount);

  if(fCheckTPCOnlyVtx)
    {
      fhEventStatTPCVtx = new TH1F("fhEventStatTPCVtx","Event statistics for TPC only vertex",9,0,9);
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(1,"FastOnly");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(2,"FastOnly+PVtx");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(3,"FastOnly+TVtx");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(4,"MB");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(5,"MB+PVtx");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(6,"MB+TVtx");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(7,"EMC");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(8,"EMC+PVtx");
      fhEventStatTPCVtx->GetXaxis()->SetBinLabel(9,"EMC+TVtx");
      fOutputList->Add(fhEventStatTPCVtx);
    }

  fhChunkQA = new TH1F("fhChunkQA","# of hybrid tracks per chunk",200,0,200);
  fOutputList->Add(fhChunkQA);

  const Int_t dim1 = 3;
  Int_t nBins1[dim1]     = {200,50,110};
  Double_t lowBin1[dim1] = {0,0,0,};
  Double_t upBin1[dim1]  = {200,0.5,1.1};

  const Int_t dim2 = 3;
  Int_t nBins2[dim2]     = {200,50,50};
  Double_t lowBin2[dim2] = {0,0,0,};
  Double_t upBin2[dim2]  = {200,0.5,50};

  const char* triggerName[3] = {"MB","EMC","MC"};
  const char* triggerTitle[3] = {"MB","EMCal-trigger","MC true"};
  const char* fraction[5] = {"MIP","30","50","70","100"};
  const char* exotic[2] = {"3GeV","5GeV"};
  const char* vertexType[2] = {"All MB","MB with vertex"};
  const char *vertexName[2] = {"All","Vertex"};
  const char *clusterizerName[2] = {"before","after"};
  const char *UEName[2] = {"charged","charged_neutral"};
  const char *UETitle[2] = {"charged","charged+neutral"};
  const char *UEEventName[2] = {"LeadingJet","Back-To-Back"};

  if(fIsMC)
    { 
      for(Int_t i=0; i<2; i++)
	{
	  fhNTrials[i] = new TH1F(Form("MC_%s_fhNTrials",triggerName[i]),Form("MC-%s: # of trials",triggerName[i]),1,0,1);
	  fOutputList->Add(fhNTrials[i]);
	  
	  fVertexGenZ[i] = new TH1F(Form("%s_fVertexGenZ",vertexName[i]),Form("Distribution of vertex z (%s); z (cm)",vertexType[i]),60,-30,30);
	  fOutputList->Add(fVertexGenZ[i]);
	}
    }

  for(Int_t i=0; i<3; i++)
    {
      if(!fIsMC && i==2) continue;

      if(fSaveQAHistos)
	{
	  if(i<2)
	    {	      
	      fEventZ[i] = new TH1F(Form("%s_fEventZ",triggerName[i]),Form("%s: Distribution of vertex z; z (cm)",triggerTitle[i]),60,-30,30);
	      fOutputList->Add(fEventZ[i]);
	    }

	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      if(!fRadius.Contains("0.4") && r==0) continue;
	      if(!fRadius.Contains("0.2") && r==1) continue;
	      if(!fRadius.Contains("0.3") && r==2) continue;

	      fJetCount[i][r] = new TH1F(Form("%s_fJetCount_%1.1f",triggerName[i],kRadius[r]),Form("%s: jet p_{T} in EMCal (R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins);
	      fOutputList->Add(fJetCount[i][r]);
	      
	      fhNeutralPtInJet[i][r] = new TH2F(Form("%s_fhNeutralPtInJet_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of neutral constituents vs jet p_{T} in EMCal(R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); p_{T} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,nTrkPtBins*10,lowTrkPtBin,upTrkPtBin);
	      fOutputList->Add(fhNeutralPtInJet[i][r]);

	      fhTrigNeuPtInJet[i][r] = new TH2F(Form("%s_fhTrigNeuPtInJet_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of triggered neutral constituents vs jet p_{T} in EMCal(R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); p_{T} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,nTrkPtBins*10,lowTrkPtBin,upTrkPtBin);
	      fOutputList->Add(fhTrigNeuPtInJet[i][r]);
	      
	      fhChargedPtInJet[i][r] = new TH2F(Form("%s_fhChargedPtInJet_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of charged constituents vs jet p_{T} in EMCal (R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); p_{T} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,nTrkPtBins*10,lowTrkPtBin,upTrkPtBin);
	      fOutputList->Add(fhChargedPtInJet[i][r]);
	      
	      fhLeadNePtInJet[i][r] = new TH2F(Form("%s_fhLeadNePtInJet_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of leading neutral constituent vs jet p_{T} in EMCal(R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); p_{T}^{ne} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,100,0,100);
	      fOutputList->Add(fhLeadNePtInJet[i][r]);
	      
	      fhLeadChPtInJet[i][r] = new TH2F(Form("%s_fhLeadChPtInJet_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of leading charged constituent vs jet p_{T} in EMCal(R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); p_{T}^{ch} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,100,0,100);
	      fOutputList->Add(fhLeadChPtInJet[i][r]);
	      
	      fhJetPtVsZ[i][r] = new TH3F(Form("%s_fhJetPtVsZ_%1.1f",triggerName[i],kRadius[r]),Form("%s: jet p_{T} vs Z_{h} vs constituent type in EMCal (R=%1.1f);p_{T}^{jet} (GeV/c); Z_{h};constituent type",triggerTitle[i],kRadius[r]),200,0,200,110,-0.05,1.05,4,0,4);
	      fOutputList->Add(fhJetPtVsZ[i][r]);
	      
	      fJetEnergyFraction[i][r] = new THnSparseF(Form("%s_fJetEnergyFraction_%1.1f",triggerName[i],kRadius[r]),Form("%s: Jet p_{T} vs radius vs energy fraction in EMCal (R=%1.1f,Z<0.98); p_{T};R;Fraction",triggerName[i],kRadius[r]),dim1,nBins1,lowBin1,upBin1);
	      fOutputList->Add(fJetEnergyFraction[i][r]);
	      
	      fJetNPartFraction[i][r] = new THnSparseF(Form("%s_fJetNPartFraction_%1.1f",triggerName[i],kRadius[r]),Form("%s: Jet p_{T} vs radius vs NPart in EMCal (R=%1.1f,Z<0.98); p_{T};R;NPart",triggerName[i],kRadius[r]),dim2,nBins2,lowBin2,upBin2);
	      fOutputList->Add(fJetNPartFraction[i][r]);
	      
	      if(i<2)
		{
		  fhJetPtInExoticEvent[i][r] = new TH1F(Form("EMC_fhJetPtInExoticEvent_%1.1f_%s",kRadius[r],exotic[i]),Form("EMC: jet p_{T} in events with exotic cluster > %s (R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c)",exotic[i],kRadius[r]),nbins,xbins);
		  fOutputList->Add(fhJetPtInExoticEvent[i][r]);
		  
		  fRelTrkCon[i][r] = new TH3F(Form("%s_fRelTrkCon_%1.1f",triggerName[i],kRadius[r]),Form("%s: jet p_{T} vs (sum p_{T}^{ch})/p_{T}^{jet} vs track class in EMCal (R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); (sum p_{T}^{ch})/p_{T}^{jet}; track class",triggerTitle[i],kRadius[r]),200,0,200,110,-0.05,1.05,3,0,3);
		  fOutputList->Add(fRelTrkCon[i][r]);
		  
		  fhFcrossVsZleading[i][r] = new TH3F(Form("%s_fhFcrossVsZleading_%1.1f",triggerName[i],kRadius[r]),Form("%s: jet p_{T} vs F_{cross} vs Z_{leading}^{ne} in EMCal (R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c);F_{cross};Z_{leading}^{ne}",triggerTitle[i],kRadius[r]),200,0,200,55,-0.05,1.05,55,-0.05,1.05);
		  fOutputList->Add(fhFcrossVsZleading[i][r]);
		  
		  fhChLeadZVsJetPt[i][r]   = new TH2F(Form("%s_fhChLeadZVsJetPt_%1.1f",triggerName[i],kRadius[r]),Form("%s: Z of leading charged constituent vs jet p_{T} in EMCal(R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c); Z_{leading}^{ch}",triggerTitle[i],kRadius[r]),nbins,xbins,100,0,1);
		  fOutputList->Add(fhChLeadZVsJetPt[i][r]);
		  
		  fhJetPtWithTrkThres[i][r]   = new TH2F(Form("%s_fhJetPtWithTrkThres_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of jets containing tracks above certain threshold (15,25,40 GeV/c) (EMCal, R=%1.1f,Z<0.98);Threshold type;p_{T}^{jet} (GeV/c)",triggerTitle[i],kRadius[r]),3,0,3,nbins,xbins);
		  fOutputList->Add(fhJetPtWithTrkThres[i][r]);
		  
		  fhJetPtWithClsThres[i][r]   = new TH2F(Form("%s_fhJetPtWithClsThres_%1.1f",triggerName[i],kRadius[r]),Form("%s: p_{T} of jets containing clusters above certain threshold (15,25,40 GeV) (EMCal, R=%1.1f,Z<0.98);Threshold type;p_{T}^{jet} (GeV/c)",triggerTitle[i],kRadius[r]),3,0,3,nbins,xbins);
		  fOutputList->Add(fhJetPtWithClsThres[i][r]);
		  
		  fhJetPtVsLowPtCons[i][r][0]   = new TH2F(Form("%s_fhJetPtVsLowPtCons_150-300MeV_%1.1f",triggerName[i],kRadius[r]),Form("%s: energy carried by constituents in 150-300MeV (EMCal, R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c);p_{T,em}^{low} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,100,0,1);
		  fOutputList->Add(fhJetPtVsLowPtCons[i][r][0]);

		  fhJetPtVsLowPtCons[i][r][1]   = new TH2F(Form("%s_fhJetPtVsLowPtCons_300-500MeV_%1.1f",triggerName[i],kRadius[r]),Form("%s: energy carried by constituents in 300-500MeV (EMCal, R=%1.1f,Z<0.98);p_{T}^{jet} (GeV/c);p_{T,em}^{low} (GeV/c)",triggerTitle[i],kRadius[r]),nbins,xbins,100,0,1);
		  fOutputList->Add(fhJetPtVsLowPtCons[i][r][1]);

		  if(fCheckTPCOnlyVtx)
		    {
		      fhJetInTPCOnlyVtx[i][r] = new TH3F(Form("%s_fhJetInTPCOnlyVtx_%1.1f",triggerName[i],kRadius[r]),Form("%s: jet pt in events with only TPC vertex (Full, R=%1.1f);p_{T}^{jet} (GeV/c);#phi;#eta",triggerTitle[i],kRadius[r]),20,0,100,36,0,360,20,-1,1);
		      fOutputList->Add(fhJetInTPCOnlyVtx[i][r]);
		    }

		  if(fStudySubEInHC)
		    {
		      for(Int_t k=0; k<5; k++)
			{
			  fhSubClsEVsJetPt[i][r][k] = new TH2F(Form("%s_fhSubClsEVsJetPt_%s_%1.1f",triggerName[i],fraction[k],kRadius[r]),Form("%s: relative %s%% subtracted cluster Et vs jet pt (R=%1.1f);jet p_{T} (GeV/c);E_{t}^{sub}/p_{T,jet}",triggerTitle[i],fraction[k], kRadius[r]),nbins,xbins,50,0,0.5);
			  fOutputList->Add(fhSubClsEVsJetPt[i][r][k]);
		      
			  fhHCTrkPtClean[i][r][k] = new TH2F(Form("%s_fhHCTrkPtClean_%s_%1.1f",triggerName[i],fraction[k],kRadius[r]),Form("%s: sum of track p_{T} that are cleanly subtracted  vs jet pt (%s%%, R=%1.1f);jet p_{T} (GeV/c);#sum(p_{T,trk}^{clean})/p_{T,jet}",triggerTitle[i],fraction[k], kRadius[r]),nbins,xbins,100,-0.005,0.995);
			  fOutputList->Add(fhHCTrkPtClean[i][r][k]);
			  
			  fhHCTrkPtAmbig[i][r][k] = new TH2F(Form("%s_fhHCTrkPtAmbig_%s_%1.1f",triggerName[i],fraction[k],kRadius[r]),Form("%s: sum of track p_{T} that are ambiguously subtracted vs jet pt (%s%%, R=%1.1f);jet p_{T} (GeV/c);#sum(p_{T,trk}^{ambig})/p_{T,jet}",triggerTitle[i],fraction[k], kRadius[r]),nbins,xbins,100,-0.005,0.995);
			  fOutputList->Add(fhHCTrkPtAmbig[i][r][k]);
			}
		    }
		}
	    }
	  if(i<2)
	    {
	      fhNMatchedTrack[i] = new TH1F(Form("%s_fhNMatchedTrack",triggerName[i]),Form("%s: # of matched tracks per cluster; N_{mth}",triggerTitle[i]),5,0,5);
	      fOutputList->Add(fhNMatchedTrack[i]);
	      
	      for(Int_t j=0; j<4; j++) 
		{
		  fhSubEVsTrkPt[i][j] = new TH2F(Form("%s_fhSubEVsTrkPt_%s",triggerName[i],fraction[j+1]),Form("%s: fractional subtracted energy (%s%% HC);#sum(p_{ch,T}^{mth}) (GeV/c);E_{sub}/#sum(P_{ch}^{mth})",triggerTitle[i],fraction[j+1]),50,0,50,110,0,1.1);
		  fOutputList->Add(fhSubEVsTrkPt[i][j]);
		} 
	    }
	}

      if(fRunUE)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      for(Int_t l=0; l<2; l++)
		{
		  fhUEJetPtNorm[i][k][l] = new TH1F(Form("%s_fhUEJetPtNorm_%s_%s",triggerName[i],UEName[k],UEEventName[l]),Form("%s: leading jet p_{T} in TPC (%s in %s event);p_{T,jet}^{ch} (GeV/c)",triggerTitle[i],UETitle[k], UEEventName[l]),nbins,xbins);
		  fOutputList->Add(fhUEJetPtNorm[i][k][l]);
	  
		  fhUEJetPtVsSumPt[i][k][l] = new TH2F(Form("%s_fhUEJetPtVsSumPt_%s_%s",triggerName[i],UEName[k],UEEventName[l]),Form("%s: leading jet p_{T} vs underlying event contribution (R=0.4,%s in %s event);p_{T,jet}^{ch} (GeV/c);p_{T,UE}^{ch} (GeV/c)",triggerTitle[i],UETitle[k], UEEventName[l]),nbins,xbins,500,0,50);
		  fOutputList->Add(fhUEJetPtVsSumPt[i][k][l]);
	  
		  fhUEJetPtVsConsPt[i][k][l] = new TH2F(Form("%s_fhUEJetPtVsConsPt_%s_%s",triggerName[i],UEName[k],UEEventName[l]),Form("%s: leading jet p_{T} vs constituent pt in UE (R=0.4,%s in %s event);p_{T,jet}^{ch} (GeV/c);p_{T,cons}^{ch} (GeV/c)",triggerTitle[i],UETitle[k], UEEventName[l]),nbins,xbins,500,0,50);
		  fOutputList->Add(fhUEJetPtVsConsPt[i][k][l]);
		}
	    }
	}

      if(fSysJetTrigEff)
	{
	  fhClsE[i] = new TH1F(Form("%s_fhClsE",triggerName[i]),Form("%s: cluster E;E (GeV)",triggerTitle[i]),1000,0,100);
	  fOutputList->Add(fhClsE[i]);
	}

      if(fSysClusterizer)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      fhSysClusterE[i][k] = new TH1F(Form("%s_fhSysClusterE_%sHC",triggerName[i],clusterizerName[k]),Form("%s: cluster E %s hadronic correction;E (GeV)",triggerTitle[i],clusterizerName[k]),100,0,100);
	      fOutputList->Add(fhSysClusterE[i][k]);

	      fhSysNCellVsClsE[i][k] = new TH2F(Form("%s_fhSysNCellVsClsE_%sHC",triggerName[i],clusterizerName[k]),Form("%s: NCell vs cluster E %s hadronic correction;E (GeV);NCell",triggerTitle[i],clusterizerName[k]),100,0,100,50,0,50);
	      fOutputList->Add(fhSysNCellVsClsE[i][k]);
	    }
	}
    }


  if(fIsMC)
    {
      if(fStudyMcOverSubE)
	{
	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      if(!fRadius.Contains("0.4") && r==0) continue;
	      if(!fRadius.Contains("0.2") && r==1) continue;
	      if(!fRadius.Contains("0.3") && r==2) continue;
	      for(Int_t i=0; i<5; i++)
		{
		  fHCOverSubE[r][i] = new TH2F(Form("%s_HC_over_sub_e_%s_%1.1f",triggerName[2],fraction[i],kRadius[r]),Form("%s: oversubtracted neutral Et by %s%% HC (R=%1.1f);particle jet p_{T} (GeV/c);#DeltaE_{t} (GeV)",triggerName[2],fraction[i],kRadius[r]),nbins,xbins,200,-49.75,50.25);
		  fOutputList->Add(fHCOverSubE[r][i]);
		  fHCOverSubEFrac[r][i] = new TH2F(Form("%s_HC_over_sub_e_frac_%s_%1.1f",triggerName[2],fraction[i],kRadius[r]),Form("%s: relative oversubtracted neutral Et fraction by %s%% HC (R=%1.1f);jet p_{T} (GeV/c);#DeltaE_{t}/p_{T}^{jet}",triggerName[2],fraction[i],kRadius[r]),nbins,xbins,200,-0.995,1.005);
		  fOutputList->Add(fHCOverSubEFrac[r][i]);
		}
	    }
	}
      if(fRunSecondaries && fAnaType==0)
	{
	  for(Int_t i=0; i<2; i++)
	    {
	      for(Int_t r=0; r<3; r++)
		{
		  fhNKFracVsJetPt[i][r]  = new TH2F(Form("%s_fhNKFracVsJetPt_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: energy fraction carried by n,k^{0}_{L} vs jet p_{T} (R=%1.1f,|#eta|<%1.1f);p_{T,jet} (GeV/c);fraction",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,200,0,1);
		  fOutputList->Add(fhNKFracVsJetPt[i][r]);

		  fhWeakFracVsJetPt[i][r]  = new TH2F(Form("%s_fhWeakFracVsJetPt_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: energy fraction carried by k^{0}_{S},hyperon vs jet p_{T} (R=%1.1f,|#eta|<%1.1f);p_{T,jet} (GeV/c);fraction",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,200,0,1);
		  fOutputList->Add(fhWeakFracVsJetPt[i][r]);
		  
		  fhJetResponseNK[i][r]  = new TH2F(Form("%s_fhJetResponseNK_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to missing n and k^{0}_{L} (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c)",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,nbins,xbins);
		  fOutputList->Add(fhJetResponseNK[i][r]);
		  
		  fhJetResponseWP[i][r]  = new TH2F(Form("%s_fhJetResponseWP_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to k^{0}_{S}, #Lambda and hyperon (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c)",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,nbins,xbins);
		  fOutputList->Add(fhJetResponseWP[i][r]);

		  fhJetResolutionNK[i][r] = new TH2F(Form("%s_fhJetResolutionNK_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to missing n and k^{0}_{L}: (p_{T,jet}^{rec}-p_{T,jet}^{gen})/p_{T,jet}^{rec} vs p_{T,jet}^{rec} (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{rec} (GeV/c);#Deltap_{T}/p_{T}",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,200,-0.995,1.005);
		  fOutputList->Add(fhJetResolutionNK[i][r]);
		  
		  fhJetResolutionWP[i][r] = new TH2F(Form("%s_fhJetResolutionWP_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to k^{0}_{S}, #Lambda and hyperon: (p_{T,jet}^{rec}-p_{T,jet}^{gen})/p_{T,jet}^{rec} vs p_{T,jet}^{rec} (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{rec} (GeV/c);#Deltap_{T}/p_{T}",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,200,-0.995,1.005);
		  fOutputList->Add(fhJetResolutionWP[i][r]);
		  
		  fhJetResponseNKSM[i][r]  = new TH2F(Form("%s_fhJetResponseNKSM_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to missing n and k^{0}_{L} via matching (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c)",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,nbins,xbins);
		  fOutputList->Add(fhJetResponseNKSM[i][r]);
		  
		  fhJetResponseWPSM[i][r]  = new TH2F(Form("%s_fhJetResponseWPSM_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet response due to k^{0}_{S}, #Lambda and hyperon via matching(R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{gen} (GeV/c);p_{T,jet}^{rec} (GeV/c)",triggerName[2],kRadius[r],i*0.5+0.5),nbins,xbins,nbins,xbins);
		  fOutputList->Add(fhJetResponseWPSM[i][r]);
		  
		  fhJetResolutionNKSM[i][r] = new TH3F(Form("%s_fhJetResolutionNKSM_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet resolution due to missing n and k^{0}_{L} via matching (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{rec} (GeV/c);#Deltap_{T,jet}/p_{T,jet}^{rec};dR",triggerName[2],kRadius[r],i*0.5+0.5),220,0,220,200,-0.995,1.005,100,0,1);
		  fOutputList->Add(fhJetResolutionNKSM[i][r]);
		  
		  fhJetResolutionWPSM[i][r] = new TH3F(Form("%s_fhJetResolutionWPSM_%1.1f_EtaCut%1.1f",triggerName[2],kRadius[r],i*0.5+0.5),Form("%s: jet resolution due tok^{0}_{S}, #Lambda and hyperon via matching (R=%1.1f,|#eta|<%1.1f);p_{T,jet}^{rec} (GeV/c);#Deltap_{T,jet}/p_{T,jet}^{rec};dR",triggerName[2],kRadius[r],i*0.5+0.5),220,0,220,200,-0.995,1.005,100,0,1);
		  fOutputList->Add(fhJetResolutionWPSM[i][r]);
		}
	    }
	}
    }

  printf("\n=======================================\n");
  printf("===== Jet task configuration ==========\n");

  if(fNonStdBranch.Length()!=0)
    {      
      const char* species[3] = {"in","ch","ne"};
      const char* algorithm[2] = {"akt","kt"};
      const char* radii[kNRadii] = {"04","02","03"};
      for(Int_t s=0; s<3; s++)
	{
	  if(!fFindChargedOnlyJet && s==1) continue;
	  if(!fFindNeutralOnlyJet && s==2) continue;
	  for(Int_t a=0; a<2; a++)
	    {
	      if(!fAlgorithm.Contains("aKt") && a==0) continue;
	      if(!fAlgorithm.Contains("kt")  && a==1) continue;
	      for(Int_t r=0; r<kNRadii; r++)
		{
		  if(!fRadius.Contains("0.4") && r==0) continue;
		  if(!fRadius.Contains("0.2") && r==1) continue;
		  if(!fRadius.Contains("0.3") && r==2) continue;
		  if(fAnaType==0)
		    {
		      fJetTCA[s][a][r] = new TClonesArray("AliAODJet",0);
		      fJetTCA[s][a][r]->SetName(Form("Jet_%s_%s_%s_%s",species[s],algorithm[a],radii[r],fNonStdBranch.Data()));
		      AddAODBranch("TClonesArray",&fJetTCA[s][a][r],fNonStdFile.Data());
		      printf("Add branch: Jet_%s_%s_%s_%s\n",species[s],algorithm[a],radii[r],fNonStdBranch.Data());
		    }

		  fDetJetFinder[s][a][r] = new AliFJWrapper(Form("DetJetFinder_%s_%s_%s_%s",species[s],algorithm[a],radii[r],fNonStdBranch.Data()),Form("DetJetFinder_%s_%s_%s_%s",species[s],algorithm[a],radii[r],fNonStdBranch.Data()));
		  if(a==0) fDetJetFinder[s][a][r]->SetAlgorithm(fastjet::antikt_algorithm);
		  if(a==1) fDetJetFinder[s][a][r]->SetAlgorithm(fastjet::kt_algorithm);
		  if(fRecombinationScheme==0)     fDetJetFinder[s][a][r]->SetRecombScheme(fastjet::E_scheme);
		  fDetJetFinder[s][a][r]->SetR(kRadius[r]);
		  fDetJetFinder[s][a][r]->SetMaxRap(0.9);
		  fDetJetFinder[s][a][r]->Clear();

		  if(s==0 && a==0)
		    {
		      if(fIsMC && fNonStdBranch.Contains("Baseline",TString::kIgnoreCase))
			{
			  if(fAnaType==0)
			    {
			      fMcTruthAntikt[r] = new TClonesArray("AliAODJet",0);
			      fMcTruthAntikt[r]->SetName(Form("Jet_mc_truth_in_akt_%s_%s",radii[r],fNonStdBranch.Data()));
			      AddAODBranch("TClonesArray",&fMcTruthAntikt[r],fNonStdFile.Data());
			      printf("Add branch: Jet_mc_truth_in_akt_%s_%s\n",radii[r],fNonStdBranch.Data());
			    }

			  fTrueJetFinder[r] = new AliFJWrapper(Form("TrueJetFinder_%s_%s_%s_%s",species[s],algorithm[a],radii[r],fNonStdBranch.Data()),Form("TrueJetFinder_%s_%s_%s_%s",species[s],algorithm[a],radii[r],fNonStdBranch.Data()));
			  fTrueJetFinder[r]->SetAlgorithm(fastjet::antikt_algorithm);
			  fTrueJetFinder[r]->SetR(kRadius[r]);
			  fTrueJetFinder[r]->SetMaxRap(0.9);
			  if(fRecombinationScheme==0) fTrueJetFinder[r]->SetRecombScheme(fastjet::E_scheme);
			  fTrueJetFinder[r]->Clear();
			}
		    }
		}
	    }
	}
    }

  fRandomGen = new TRandom3(0);
  if(fCheckTrkEffCorr && fAnaType==0)
    {
      TFile f ("/project/projectdirs/alice/marr/Analysis/2.76Jet/CorrectionFunctions/TrkEffFit.root","read");
      for(Int_t i=0; i<3; i++)
	{
	  fTrkEffFunc[i] = new TF1(*((TF1*)f.Get(Form("Trk_eff_fit_%d",i+1))));
	}
      f.Close();

      if(fPeriod.CompareTo("lhc11a",TString::kIgnoreCase)==0 || fPeriod.CompareTo("lhc11aJet",TString::kIgnoreCase)==0)
	{
	  TFile f1 ("/project/projectdirs/alice/marr/Analysis/2.76Jet/CorrectionFunctions/TrkEffSampling.root","read");
	  for(Int_t j=0; j<2; j++)
	    {
	      Int_t tmp = j;
	      if(fPeriod.CompareTo("lhc11aJet",TString::kIgnoreCase)==0) tmp = 2;
	      for(Int_t r=0; r<2; r++)
		{
		  fhCorrTrkEffPtBin[j][r] = new TH1F(*((TH1F*)f1.Get(Form("%s_%s_NTrackPerPtBin_%1.1f",fPeriod.Data(),triggerName[tmp],kRadius[r]))));
		  for(Int_t i=0; i<kNBins; i++)
		    {
		      fhCorrTrkEffSample[j][r][i] = new TH1F(*((TH1F*)f1.Get(Form("%s_%s_ChTrkPt_Bin%d_%1.1f",fPeriod.Data(),triggerName[tmp],i+1,kRadius[r]))));
		    }
		}
	    }
	  f1.Close();
	}
    }

  if(fRunSecondaries && fAnaType==0)
    {
      const char *secondaryName[3] = {"k0S","lamda","hyperon"};
      TFile f2 ("/project/projectdirs/alice/marr/Analysis/2.76Jet/CorrectionFunctions/SecondaryResponse.root","read");
      for(Int_t j=0; j<3; j++)
	{
	  fhSecondaryResponse[j] = new TH2F(*((TH2F*)f2.Get(Form("DetectorReponse_%s",secondaryName[j]))));
	}
    }

  if(fCheckTriggerMask)
    {
      char *name = "TriggerCurve.root";
      if(fAnaType==0) name = "/project/projectdirs/alice/marr/Analysis/2.76Jet/CorrectionFunctions/TriggerCurve.root";
      TFile f (name,"read");
      fTriggerMask = new TH2F(*((TH2F*)f.Get("lhc11a_TriggerMask")));
      if(fOfflineTrigger)
	{
	  for(Int_t i=0; i<10; i++)
	    {
	      fTriggerEfficiency[i] = new TF1(*((TF1*)f.Get(Form("lhc11a_TriggerEfficiency_SM%d_fit",i))));
	      fTriggerCurve[i] = new TH1D(*((TH1D*)f.Get(Form("lhc11a_TriggerCurve_SM%d",i))));
	    }   
	}
      f.Close();
    }

  fClusterEResolution = new TF1("fClusterEResolution","sqrt([0]^2+[1]^2*x+([2]*x)^2)*0.01");
  fClusterEResolution->SetParameters(4.35,9.07,1.63);

  fGeom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  if(!fGeom)
    {
      AliError("No EMCal geometry is available!");
      return;
    }
  fRecoUtil = new AliEMCALRecoUtils();
  if(fRejectExoticCluster || fRejectExoticTrigger)
    fRecoUtil->SwitchOnRejectExoticCluster();
  else
    {
      fRecoUtil->SwitchOffRejectExoticCluster();
      fRecoUtil->SwitchOffRejectExoticCell();
    }
  if(fRemoveBadChannel)
    {
      fRecoUtil->SwitchOnBadChannelsRemoval();
      // Remove the problematic SM4 region due to wrong event sequence
      for(Int_t ieta=0; ieta<36; ieta++)
	for(Int_t iphi=0; iphi<8; iphi++)
	  fRecoUtil->SetEMCALChannelStatus(4,ieta,iphi);
    }
  else   
    fRecoUtil->SwitchOffBadChannelsRemoval();

  fRecoUtil->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
  if(fSysNonLinearity) 
    {
      fNonLinear = new TF1("TB_oldBest","([0])*(1./(1.+[1]*exp(-x/[2]))*1./(1.+[3]*exp((x-[4])/[5])))*1/([6])",0.1,110.);
      fNonLinear->SetParameters(0.99078, 0.161499, 0.655166, 0.134101, 163.282, 23.6904, 0.978);
    }
  if(fSysTrkClsMth)
    {
      fRecoUtil->SwitchOnCutEtaPhiSeparate();
      fRecoUtil->SetCutEta(fCutdEta);
      fRecoUtil->SetCutPhi(fCutdPhi);
    }

  fTrackArray   = new TObjArray();
  fTrackArray->SetOwner(1);
  fClusterArray = new TObjArray();
  fClusterArray->SetOwner(1);
  fMcPartArray = new TArrayI();


  //error calculation in THnSparse
  Int_t nObj = fOutputList->GetEntries();
  for(Int_t i=0; i<nObj; i++)
    {
      TObject *obj = (TObject*) fOutputList->At(i);
      if (obj->IsA()->InheritsFrom( "THnSparse" ))
      	{
      	  THnSparseF *hn = (THnSparseF*)obj;
      	  hn->Sumw2();
      	}
    }


  PrintConfig();
  BookHistos();
  PostData(1, fOutputList);
}



//________________________________________________________________________
void AliAnalysisTaskFullppJet::BookHistos()
{
  // Book histograms.

  if(fVerbosity>10) printf("[i] Booking histograms \n");
  return;
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::UserExec(Option_t *) 
{
  //
  // Main loop, called for each event.
  //

  // Get event pointers, check for signs of life
  Double_t vtxTrueZ = -100;
  if(fIsMC)
    {
      fMC = MCEvent();
      if (!fMC) 
	{
	  printf("ERROR: Could not retrieve MC event");
	  return;
	}
      fStack = fMC->Stack();
      TParticle *particle = fStack->Particle(0);
      if(particle) vtxTrueZ = particle->Vz(); // True vertex z in MC events
      if(fVerbosity>10) printf("Generated vertex coordinate: (x,y,z) = (%4.2f, %4.2f, %4.2f)\n", particle->Vx(), particle->Vy(), particle->Vz());
    }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) 
    {
      fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    }

  if (!fESD && !fAOD) 
    {
      AliError("Neither fESD nor fAOD available");
      return;
    }

  fhJetEventStat->Fill(0.5);
  fhJetEventCount->Fill(0.5);
  if(fIsMC)
    {
      GetMCInfo();
      if(fVertexGenZ[0]) fVertexGenZ[0]->Fill(vtxTrueZ);
    }

  // Centrality, vertex, other event variables...
  if(fESD)
    {
      UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      if (trigger==0)  return; 
      if(fCheckTPCOnlyVtx) CheckTPCOnlyVtx(trigger);
      if(!fIsMC)
	{
	  if (trigger & AliVEvent::kFastOnly) return;  // Reject fast trigger cluster
	  if(fPeriod.CompareTo("lhc11a",TString::kIgnoreCase)==0)
	    {
	      if (trigger & AliVEvent::kMB)       fTriggerType = 0;
	      else if(trigger & AliVEvent::kEMC1) fTriggerType = 1;
	      else fTriggerType = -1;
	    }
	  else if (fPeriod.CompareTo("lhc11c",TString::kIgnoreCase)==0 || fPeriod.CompareTo("lhc11d",TString::kIgnoreCase)==0)
	    {
	      if (trigger & AliVEvent::kINT7)     fTriggerType = 0;
	      else if(trigger & AliVEvent::kEMC7) fTriggerType = 1;
	      else fTriggerType = -1;
	    }
	  else
	    {
	      return;
	    }
	}
      else
	{
	  if(!fPhySelForMC) fTriggerType = 0;
	  else if (trigger & AliVEvent::kAnyINT) fTriggerType = 0;
	  else fTriggerType = -1;

	  if(fOfflineTrigger)
	    {
	      RunOfflineTrigger();
	      if(fIsEventTriggerBit) fTriggerType = 1;
	    }
	}

      if(fTriggerType==-1)
	{
	  if(fVerbosity>10) printf("Error: worng trigger type %s\n",(fESD->GetFiredTriggerClasses()).Data());
	  return;
	}
    }

  fhJetEventCount->Fill(1.5+fTriggerType*4);
  if(fRejectPileup && fESD->IsPileupFromSPD() ) return; // reject pileup
  fhJetEventStat->Fill(11.5);
  fhJetEventCount->Fill(2.5+fTriggerType*4);

 
  fIsTPCOnlyVtx = 0;
  // Reject LED events
  if (IsLEDEvent()) 
    {
      fhJetEventStat->Fill(8.5);
      return;
    }

  fhJetEventStat->Fill(1.5+fTriggerType*2);

  // Check if primary vertex exists
  if (!HasPrimaryVertex()) return;
  fhJetEventStat->Fill(5.5+fTriggerType);
  fhJetEventCount->Fill(3.5+fTriggerType*4);

  const AliESDVertex* vtx = fESD->GetPrimaryVertex();
  Double_t zVertex    = vtx->GetZ();
  if(fEventZ[fTriggerType]) fEventZ[fTriggerType]->Fill(zVertex);
  if(fVertexGenZ[1]) fVertexGenZ[1]->Fill(vtxTrueZ);

  // Check if |Z_vtx|<10cm
  if( TMath::Abs(zVertex) > fZVtxMax ) return;
  fhJetEventCount->Fill(4.5+fTriggerType*4);
  
  // Check if only TPC vertex exists primitive
  fIsTPCOnlyVtx = IsTPCOnlyVtx();

  if(!fIsMC && fTriggerType==1) 
    {
      // Check if event has valid trigger bit
      CheckEventTriggerBit();
      if(fIsEventTriggerBit) 
	{
	  fhJetEventStat->Fill(7.5);
	  fhJetEventCount->Fill(9.5);
	}
      else return;
    }
  fhJetEventStat->Fill(2.5+fTriggerType*2);
  if(fIsTPCOnlyVtx) fhJetEventStat->Fill(9.5+fTriggerType);

  // Check if event contains exotic clusters
  CheckExoticEvent();

  // Write jet tree into AOD output in local analysis mode
  if(fAnaType==0) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
   
  // Clean up arrays
  for(Int_t s=0; s<3; s++)
    {
      for(Int_t a=0; a<2; a++)
	{
	  for(Int_t r=0; r<kNRadii; r++)
	    {
	      if(fJetTCA[s][a][r])       fJetTCA[s][a][r]->Delete();
	      if(fDetJetFinder[s][a][r]) fDetJetFinder[s][a][r]->Clear();

	      if(s==0 && a==0)
		{
		  if(fMcTruthAntikt[r]) fMcTruthAntikt[r]->Delete();
		  if(fTrueJetFinder[r]) fTrueJetFinder[r]->Clear();
		}
	    }
	}
    }

  if(fVerbosity>5) printf("# of jets after clear: %d\n",fJetTCA[0][0][0]->GetEntries());

  if(fTrackArray) fTrackArray->Delete();
  if(fClusterArray) fClusterArray->Delete();
  fMcPartArray->Reset();

  if (fESD) 
     {
       // get the tracks and fill the input vector for the jet finders
       GetESDTrax();
       
       // get EMCal clusters and fill the input vector for the jet finders
       GetESDEMCalClusters();

       for(Int_t s=0; s<3; s++)
	 {
	   for(Int_t a=0; a<2; a++)
	     {
	       for(Int_t r=0; r<kNRadii; r++)
		 {
		   //Detector jets
		   if(fDetJetFinder[s][a][r]) 
		     {
		       FindDetJets(s,a,r);
		       if(fJetTCA[s][a][r]) FillAODJets(fJetTCA[s][a][r], fDetJetFinder[s][a][r], 0);
		       if(s==0 && a==0 && fNonStdBranch.Contains("Baseline",TString::kIgnoreCase))
			 {
			   if(fSaveQAHistos)    AnalyzeJets(fDetJetFinder[s][a][r],fTriggerType, r); // run analysis on jets
			 }
		     }
		 }
	     }
	 }

       if(fRunUE && fDetJetFinder[1][0][0]) RunAnalyzeUE(fDetJetFinder[1][0][0],fTriggerType,0); // Run analysis on underlying event
     }

  if(fIsMC)
    {
      for(Int_t r=0; r<kNRadii; r++)
	if(fTrueJetFinder[r]) ProcessMC(r); // find particle level jets
    }

  // Fast Jet calls END --------------------------------------------------------
  
  PostData(1, fOutputList);
  return;
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::FindDetJets(const Int_t s, const Int_t a, const Int_t r)
{
  // 
  // Jet finding is happening here
  //

  Bool_t isCh = kTRUE;
  Bool_t isNe = kTRUE;
  if(s==1) isNe = kFALSE;
  if(s==2) isCh = kFALSE;

  if(isCh)
    {
      // Feed in charged tracks
      Int_t countTracks = 0;
      Int_t ntracks = fTrackArray->GetEntriesFast();
      for(Int_t it=0; it<ntracks; it++)
	{
	  AliESDtrack *t = (AliESDtrack*) fTrackArray->At(it);
	  if(!t) continue;
	  countTracks++;
	  Double_t e = t->P();
	  Double_t pt = t->Pt();
	  if(s==1 && fRunUE && pt<1) continue;
	  if(fRecombinationScheme==0) e = TMath::Sqrt(t->P()*t->P()+0.139*0.139);
	  if(fSysTrkPtRes) pt += GetSmearedTrackPt(t);
	  Double_t phi = t->Phi();
	  Double_t eta = t->Eta();
	  Double_t px = pt*TMath::Cos(phi);
	  Double_t py = pt*TMath::Sin(phi);
	  if(fConstrainChInEMCal && ( TMath::Abs(eta)>0.7 || phi>kPI || phi<TMath::DegToRad()*80) ) continue;
	  fDetJetFinder[s][a][r]->AddInputVector(px,py,t->Pz(), e, it+1);
	  if(fVerbosity>10) printf("%d: m_{ch}=%f\n",it+1,t->P()*t->P()-t->Px()*t->Px()-t->Py()*t->Py()-t->Pz()*t->Pz());
	}
      if(fVerbosity>5) printf("[i] # of tracks filled: %d\n",countTracks);
    }
  if(isNe)
    {
      // Feed in EMCal clusters
      Double_t vertex[3] = {0, 0, 0};
      fESD->GetVertex()->GetXYZ(vertex) ;
      TLorentzVector gamma;
      
      Int_t countClusters = 0;
      Int_t nclusters = fClusterArray->GetEntriesFast();
      for (Int_t ic = 0; ic < nclusters; ic++)
	{
	  AliESDCaloCluster * cl = (AliESDCaloCluster *)fClusterArray->At(ic);
	  if(!cl) continue;
	  cl->GetMomentum(gamma, vertex);
	  countClusters++;
	  Double_t e = gamma.P();
	  if(fRecombinationScheme==0) e = TMath::Sqrt(gamma.P()*gamma.P()+0.139*0.139);
	  fDetJetFinder[s][a][r]->AddInputVector(gamma.Px(), gamma.Py(), gamma.Pz(), e, (ic+1)*(-1));
	}
      if(fVerbosity>5) printf("[i] # of clusters filled: %d\n",countClusters);
    }
  // Run jet finding
  fDetJetFinder[s][a][r]->Run();
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::FillAODJets(TClonesArray *fJetArray, AliFJWrapper *jetFinder, const Bool_t isTruth)
{
  //
  // Fill the found jets into AOD output
  // Only consider jets pointing to EMCal and with pt above 1GeV/c
  //

  Int_t radiusIndex = 0;
  TString arrayName = fJetArray->GetName();
  if(arrayName.Contains("_02_"))  radiusIndex = 1;
  if(arrayName.Contains("_03_"))  radiusIndex = 2;
  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();
  if(fVerbosity>5 && radiusIndex==0) printf("[i] # of jets in %s : %d\n",fJetArray->GetName(),(Int_t)jetsIncl.size());
  AliAODJet *jet = 0;
  Int_t jetCount = 0;
  for(UInt_t ij=0; ij<jetsIncl.size(); ij++)
    {
      if(fVerbosity>10) printf("fastjet: eta=%f, phi=%f\n",jetsIncl[ij].eta(),jetsIncl[ij].phi()*TMath::RadToDeg());
      if(!IsGoodJet(jetsIncl[ij],0)) continue;

      AliAODJet tmpJet (jetsIncl[ij].px(), jetsIncl[ij].py(), jetsIncl[ij].pz(), jetsIncl[ij].E());
      jet = new ((*fJetArray)[jetCount]) AliAODJet(tmpJet);
      jetCount++;
      if(fVerbosity>10 && radiusIndex==0) printf("AOD jet: ij=%d, pt=%f, eta=%f, phi=%f\n",ij, jet->Pt(), jet->Eta(),jet->Phi()*TMath::RadToDeg());

      jet->GetRefTracks()->Clear();
      std::vector<fastjet::PseudoJet> constituents = jetFinder->GetJetConstituents(ij);
      Double_t totE=0, totPt=0, totChPt=0, leadChPt=0, neE=0, totNePt=0, leadNePt=0, leadPt=0;
      Double_t leadTrkType=0, nChPart = 0, nNePart = 0;
      Bool_t isHighPtTrigger = kFALSE, isTriggering = kFALSE, isHyperTrack = kFALSE;
      Int_t leadIndex = -1;
      for(UInt_t ic=0; ic<constituents.size(); ic++)
	{
	  if(fVerbosity>10 && radiusIndex==0) printf("ic=%d: user_index=%d, E=%f\n",ic,constituents[ic].user_index(),constituents[ic].E());
	  if(constituents[ic].user_index()<0)
	    {
	      nNePart ++;
	      totNePt += constituents[ic].perp();
	      if(constituents[ic].perp()>leadNePt)
		leadNePt = constituents[ic].perp();

	      neE += constituents[ic].E();
	      if(constituents[ic].perp()>fClsEtMax[fTriggerType])
		isHighPtTrigger = kTRUE;
              if(!isTruth)
                {
                  AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(constituents[ic].user_index()*(-1)-1);
                  if(cluster->Chi2()>0.5) isTriggering = kTRUE;
                }
	    }
	  else
	    {
	      totChPt += constituents[ic].perp();
	      nChPart ++;
	      if(constituents[ic].perp()>leadChPt)
		{
		  leadChPt = constituents[ic].perp();
		  if(!isTruth)
		    {
		      AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[ic].user_index()-1);
		      leadTrkType = track->GetTRDQuality();
		      if(track->GetIntegratedLength()>500) isHyperTrack = kTRUE;
		    }    
		}
	      if(constituents[ic].perp()>fTrkPtMax[fTriggerType])
		isHighPtTrigger = kTRUE;
	    }
          TParticle part;
          part.SetMomentum(constituents[ic].px(),constituents[ic].py(),constituents[ic].pz(),constituents[ic].E());
          jet->AddTrack(&part); //The references are not usable, this line is aimed to count the number of contituents
	  totE  += constituents[ic].E();
	  totPt += constituents[ic].perp();
	  if(constituents[ic].perp()>leadPt)
	    {
	      leadPt = constituents[ic].perp();
	      leadIndex = ic;
	    } 
	}
      if(fIsEventTriggerBit)   jet->SetTrigger(AliAODJet::kTRDTriggered);
      if(isHighPtTrigger)  jet->SetTrigger(AliAODJet::kHighTrackPtTriggered);
      if(fTriggerType==1)  jet->SetTrigger(AliAODJet::kEMCALTriggered);
      if(GetLeadingZ(ij,jetFinder) > 0.98 )  jet->SetTrigger(AliAnalysisTaskFullppJet::kHighZ);
      if(isTriggering) jet->SetTrigger(AliAnalysisTaskFullppJet::kTrigger);
      if(isHyperTrack) jet->SetTrigger(AliAnalysisTaskFullppJet::kSuspicious);
      if(fIsTPCOnlyVtx) jet->SetTrigger(AliAnalysisTaskFullppJet::kTPCOnlyVtx);
      if(constituents[leadIndex].user_index()>0) jet->SetTrigger(AliAnalysisTaskFullppJet::kLeadCh);

      if(jetsIncl[ij].E()>0)  jet->SetNEF(neE/jetsIncl[ij].E());
      jet->SetPtLeading(leadPt);
      jet->SetPtSubtracted(leadTrkType,0);

      Double_t  effAErrCh = 0, effAErrNe = leadTrkType;
      Double_t chBgEnergy = 10, neBgEnergy = nNePart;
      if(!isTruth)
	{
	  effAErrCh = GetMeasuredJetPtResolution(ij,jetFinder);
	  if(fIsMC && constituents[leadIndex].user_index()>0)
	    {
	      AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[leadIndex].user_index()-1);
	      Double_t pt = track->Pt();
	      Int_t ipart = track->GetLabel();
	      if(ipart>-1 && ipart<fMC->GetNumberOfTracks())
		{
		  AliVParticle* vParticle = fMC->GetTrack(ipart);
		  if(vParticle)
		    {
		      Double_t truePt = vParticle->Pt();
		      chBgEnergy = (truePt-pt)/truePt;
		    }
		}
	    }
	}

      jet->SetEffArea(leadPt, jetFinder->GetJetArea(ij),effAErrCh,effAErrNe);
      jet->SetBgEnergy(chBgEnergy,neBgEnergy);
      if(fVerbosity>10 && isTruth) cout<<jet->ErrorEffectiveAreaCharged()<<"  "<<jet->ErrorEffectiveAreaNeutral()<<endl;
      if(fVerbosity>10) cout<<"jet pt="<<jetsIncl[ij].perp()<<", nef="<<jet->GetNEF()<<", trk eff corr="<<chBgEnergy<<endl;

      if(fVerbosity>5)
	printf("# of ref tracks: %d = %d, and nef=%f\n",jet->GetRefTracks()->GetEntries(), (Int_t)constituents.size(), jet->GetNEF());

      // For catch good high-E jets
      if(fSpotGoodJet && !isTruth)
	{
	  if(jetsIncl[ij].perp()>100)  
	    {
	      printf("\n\n--- HAHAHA: High pt jets ---\n");
	      printf("File: %s, event = %d, pt=%f\n",CurrentFileName(),(Int_t)Entry(),jetsIncl[ij].perp());
	      printf("%s , pt < %f\n", fJetArray->GetName(), fTrkPtMax[1]);
	      printf("Jet: E=%f, eta=%f, phi=%f, # of constituents=%d, nef=%f\n",jetsIncl[ij].E(),jetsIncl[ij].eta(), jetsIncl[ij].phi()*TMath::RadToDeg(), (Int_t)constituents.size(),jet->GetNEF());
	      for(UInt_t ic=0; ic<constituents.size(); ic++)
		{
		  if(constituents[ic].user_index()<0)
		    {
		      AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(constituents[ic].user_index()*(-1)-1);
		      printf("id = %d, cluster with pt = %f, ncell = %d\n",ic,constituents[ic].perp(), cluster->GetNCells());
		    }
		  else
		    {
		      AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[ic].user_index()-1);
		      printf("id = %d, track with pt = %f, track class = %d, originalPt = %f\n",ic,constituents[ic].perp(),(Int_t)track->GetTRDQuality(), track->GetIntegratedLength());

		    }
		}
	      printf("==============================\n\n");
	    }
	}
      // End of catching good high-E jets
    }
  if(fVerbosity>5) printf("%s has %d jets\n",fJetArray->GetName(), fJetArray->GetEntries());
}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsGoodJet(fastjet::PseudoJet jet, Double_t rad)
{
  //
  // Check if the jet pt and direction fulfill the requirement
  //

  if(jet.perp()<1) return kFALSE;
  if(TMath::Abs(jet.eta())>(0.7-rad)) return kFALSE;
  if(jet.phi() < (80*TMath::DegToRad()+rad) || jet.phi() > (180*TMath::DegToRad()-rad) ) return kFALSE;
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::ProcessMC(const Int_t r)
{
  //
  // Main function for jet finding in MC
  //

  Int_t npart = fMC->GetNumberOfTracks();
  fMcPartArray->Set(npart);
  Int_t countPart = 0;
  for(Int_t ipart=0; ipart<npart; ipart++)
    {
      AliVParticle* vParticle = fMC->GetTrack(ipart);
      if(!IsGoodMcPartilce(vParticle,ipart)) continue;

      Int_t pdgCode = vParticle->PdgCode();
      if( fRejectNK && (pdgCode==130 || pdgCode==2112) ) continue;

      if( fRejectWD && (pdgCode==310 || pdgCode==3112 || pdgCode==3122 || pdgCode==3222 || pdgCode==3312 || pdgCode==3322 || pdgCode==3334)) continue;

      if( fChargedMC && vParticle->Charge()==0 ) continue;

      Int_t index = 1;
      if(vParticle->Charge()==0) { index=-1; }     
      fMcPartArray->AddAt(ipart, countPart);
      countPart++;

      fTrueJetFinder[r]->AddInputVector(vParticle->Px(), vParticle->Py(), vParticle->Pz(), vParticle->P(), (ipart+1)*index);
      if(fVerbosity>10) 
	printf("Input particle: ipart=%d, pdg=%d, species=%s, charge=%d, E=%4.3f\n",(ipart+1)*index,pdgCode,vParticle->GetName(), vParticle->Charge(),vParticle->P());
    }
  fMcPartArray->Set(countPart);
  fTrueJetFinder[r]->Run();

  if(fMcTruthAntikt[r]) FillAODJets(fMcTruthAntikt[r], fTrueJetFinder[r], 1);
  if(fSaveQAHistos)     AnalyzeJets(fTrueJetFinder[r], 2, r);
  if(fRunUE && r==0)    RunAnalyzeUE(fTrueJetFinder[r], 2, 1);

  // Run analysis on secondary particles
  if(fRunSecondaries && fAnaType==0) 
    {
      for(Int_t i=0; i<2; i++)
	{
	  AnalyzeSecondaryContribution(fTrueJetFinder[r],r,i);
	  AnalyzeSecondaryContributionViaMatching(fTrueJetFinder[r],r,0,i);
	  AnalyzeSecondaryContributionViaMatching(fTrueJetFinder[r],r,1,i);
	}
    }
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::AnalyzeJets(AliFJWrapper *jetFinder, const Int_t type,  const Int_t r)
{
  //
  // Fill all the QA plots for jets
  // Especailly all the constituents plot must be filled here since the information
  // is not available in the output AOD

  Double_t vertex[3] = {0, 0, 0};
  fESD->GetVertex()->GetXYZ(vertex) ;
  TLorentzVector gamma;

  const Int_t nBins = fJetEnergyFraction[type][r]->GetAxis(1)->GetNbins();
  Float_t radiusCut[nBins];
  Float_t eFraction[nBins];
  Int_t   nPart[nBins];
  for(Int_t i=0; i<nBins; i++)
    {
      radiusCut[i] = (fJetEnergyFraction[type][r]->GetAxis(1)->GetBinWidth(1)) * (i+1);
      eFraction[i] = 0.;
      nPart[nBins] = 0;
    }
  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();

  for(UInt_t ij=0; ij<jetsIncl.size(); ij++)
    {
      if(type<2 && fCheckTPCOnlyVtx && fIsTPCOnlyVtx) 
      	{
      	  fhJetInTPCOnlyVtx[type][r]->Fill(jetsIncl[ij].perp(),jetsIncl[ij].phi()*TMath::RadToDeg(),jetsIncl[ij].eta());
      	}
      if(!IsGoodJet(jetsIncl[ij],kRadius[r])) continue; // Fidiual cut
      Float_t jetEta = jetsIncl[ij].eta();
      Float_t jetPhi = jetsIncl[ij].phi();
      Float_t jetE   = jetsIncl[ij].E();
      Float_t jetPt  = jetsIncl[ij].perp();

      std::vector<fastjet::PseudoJet> constituents = jetFinder->GetJetConstituents(ij);
      Double_t neE=0, leadingZ = 0, maxPt = 0;
      Int_t constituentType = -1, leadingIndex = 0; 
      for(UInt_t ic=0; ic<constituents.size(); ic++)
	{
	  if(constituents[ic].perp()>maxPt)
	    {
	      maxPt = constituents[ic].perp();
	      leadingIndex = constituents[ic].user_index();
	    }

	  if(constituents[ic].user_index()<0)
	    {
	      neE += constituents[ic].E();
	      constituentType = 3;
	    }
	  else
	    {
	      if(type==2) constituentType = 0;
	      else
		{
		  AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[ic].user_index()-1);
		  constituentType = (Int_t)track->GetTRDQuality();
		}
	    }
	  Double_t cz = GetZ(constituents[ic].px(),constituents[ic].py(),constituents[ic].pz(),jetsIncl[ij].px(),jetsIncl[ij].py(),jetsIncl[ij].pz());
	  fhJetPtVsZ[type][r]->Fill(jetPt,cz, constituentType);
	  if(cz>leadingZ) leadingZ = cz;
	}

      if(type<2 && leadingIndex<0)
	{
	  AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(leadingIndex*(-1)-1);
	  fhFcrossVsZleading[type][r]->Fill(jetPt,GetExoticEnergyFraction(cluster),leadingZ);
	}

      if(leadingZ > 0.98) continue;  // Z cut

      fJetCount[type][r]->Fill(jetPt);
      if(type==1 && fIsExoticEvent3GeV) fhJetPtInExoticEvent[0][r]->Fill(jetPt);
      if(type==1 && fIsExoticEvent5GeV) fhJetPtInExoticEvent[1][r]->Fill(jetPt);

      Double_t totTrkPt[3] = {0.,0.,0.};
      Double_t chPt = 0;
      for(Int_t i=0; i<nBins; i++) { eFraction[i] = 0.; nPart[i] = 0;}
      Double_t mcSubE=0.;
      Double_t subClsE[5] = {0,0,0,0,0};
      Double_t subTrkPtClean[5] = {0,0,0,0,0};
      Double_t subTrkPtAmbig[5] = {0,0,0,0,0};
      Double_t fraction[5] = {270,0.3,0.5,0.7,1};
      Double_t leadChPt=0., leadNePt=0.;
      Int_t leadChIndex = -1;
      for(UInt_t ic=0; ic<constituents.size(); ic++)
	{
	  Float_t partEta = constituents[ic].eta();
	  Float_t partPhi = constituents[ic].phi();
	  Float_t partE   = constituents[ic].E();
	  Float_t partPt  = constituents[ic].perp();

	  if(constituents[ic].user_index()<0)
	    {
	      fhNeutralPtInJet[type][r]->Fill(jetPt, partPt);
	      if(partPt>leadNePt)
		leadNePt = partPt;
	      if(type<2)
		{
		  AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(constituents[ic].user_index()*(-1)-1);
		  Double_t clsE = cluster->E();
		  if(cluster->Chi2()>0.5) fhTrigNeuPtInJet[type][r]->Fill(jetPt, partPt);
		  if(fStudySubEInHC || fStudyMcOverSubE)
		    {
		      cluster->GetMomentum(gamma, vertex);
		      Double_t sinTheta = TMath::Sqrt(1-TMath::Power(gamma.CosTheta(),2));
		      
		      Double_t subEtmp = cluster->GetDispersion();
		      Double_t mipETmp = cluster->GetEmcCpvDistance();
		      mcSubE += cluster->GetDistanceToBadChannel()*sinTheta;
		      subClsE[0]   += ((mipETmp>clsE)?clsE:mipETmp)*sinTheta;
		      for(Int_t j=1; j<5; j++)
			{
			  subClsE[j]   += ((fraction[j]*subEtmp>clsE)?clsE:fraction[j]*subEtmp)*sinTheta;
			}
		    }
		}
	    }
	  else
	    {
	      if(partPt>leadChPt) {leadChPt = partPt;leadChIndex=ic;}
	      fhChargedPtInJet[type][r]->Fill(jetPt, partPt);
	      chPt += constituents[ic].perp();
	      if(type<2)
		{
		  AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[ic].user_index()-1);
		  Int_t trkType = (Int_t)track->GetTRDQuality();
		  totTrkPt[trkType] += partPt;
		  Int_t clusterIndex = track->GetEMCALcluster();
		  Int_t clusterPos = -1;
		  Bool_t isExist = kFALSE;
		  for(Int_t j=0; j<fClusterArray->GetEntriesFast(); j++)
		    {
		      AliESDCaloCluster *cluster = (AliESDCaloCluster*)fClusterArray->At(j);
		      if( clusterIndex == cluster->GetID() )
			{
			  isExist = kTRUE;
			  clusterPos = j;
			  break;
			}
		    }
		  if(isExist)
		    {
		      AliESDCaloCluster *cluster = (AliESDCaloCluster*)fClusterArray->At(clusterPos);
		      Double_t subEtmp = cluster->GetDispersion();
		      Double_t mipETmp = cluster->GetEmcCpvDistance();
		      for(Int_t k=0; k<5; k++)
			{
			  if(k==0)
			    {
			      if(mipETmp>cluster->E()) subTrkPtClean[k] += partPt;
			      else             subTrkPtAmbig[k] += partPt;
			    }
			  else
			    {
			      if(fraction[k]*subEtmp>cluster->E()) subTrkPtClean[k] += partPt;
			      else                         subTrkPtAmbig[k] += partPt;
			    }
			}
		    }
		}
	    }

	  Float_t dR = TMath::Sqrt( (partEta-jetEta)*(partEta-jetEta) + (partPhi-jetPhi)*(partPhi-jetPhi) );
	  for(Int_t i=0; i<nBins; i++)
	    {
	      if( dR < radiusCut[i] )
		{
		  eFraction[i] += partE;
		  nPart[i]++;
		}
	    }
	}

      fhLeadNePtInJet[type][r]->Fill(jetPt, leadNePt);
      fhLeadChPtInJet[type][r]->Fill(jetPt, leadChPt);
      if(type<2 && leadChIndex>-1)
	{
	  Double_t cz = GetZ(constituents[leadChIndex].px(),constituents[leadChIndex].py(),constituents[leadChIndex].pz(),jetsIncl[ij].px(),jetsIncl[ij].py(),jetsIncl[ij].pz());
	  fhChLeadZVsJetPt[type][r]->Fill(jetPt, cz);
	}

      if((fStudySubEInHC||fStudyMcOverSubE) && type<2)
	{
	  for(Int_t i=0; i<5; i++)
	    {
	      if(fStudySubEInHC) 
		{
		  fhSubClsEVsJetPt[type][r][i]->Fill(jetPt-subClsE[4],subClsE[i]/(jetPt-subClsE[4]));
		  fhHCTrkPtClean[type][r][i]->Fill(jetPt-subClsE[4],subTrkPtClean[i]/(jetPt-subClsE[4]));
		  fhHCTrkPtAmbig[type][r][i]->Fill(jetPt-subClsE[4],subTrkPtAmbig[i]/(jetPt-subClsE[4]));
		}
	      if(type==0 && fIsMC && fStudyMcOverSubE)
		{
		  fHCOverSubE[r][i]->Fill(jetPt-subClsE[4],subClsE[i]-mcSubE);
		  fHCOverSubEFrac[r][i]->Fill(jetPt-subClsE[4],(subClsE[i]-mcSubE)/(jetPt-subClsE[4]));
		}
	    }
	}
      if(type<2 && chPt>0)
	{
	  for(Int_t i=0; i<3; i++)
	    {
	      fRelTrkCon[type][r]->Fill(jetPt,totTrkPt[i]/chPt,i);
	    }
	}


      for(Int_t ibin=0; ibin<nBins; ibin++)
	{
	  Double_t fill1[3] = {jetPt,radiusCut[ibin]-0.005,eFraction[ibin]/jetE};
	  fJetEnergyFraction[type][r]->Fill(fill1);
	  Double_t fill2[3] = {jetPt,radiusCut[ibin]-0.005,nPart[ibin]};
	  fJetNPartFraction[type][r]->Fill(fill2);
	}

      // Get the jet pt containing tracks or clusters above some threshold
      if(type<2)
	{
	  Double_t thres[3] = {15,25,40};
	  Int_t okCh[3] = {0,0,0};
	  Int_t okNe[3] = {0,0,0};
	  Double_t lowPt[2] = {0.,0.};
	  for(UInt_t ic=0; ic<constituents.size(); ic++)
	    {
	      Float_t partPt  = constituents[ic].perp();

	      if(partPt<0.3) lowPt[0] += partPt;
	      else if(partPt<0.5) lowPt[1] += partPt;

	      if(constituents[ic].user_index()>0)
		{
		  for(Int_t it=0; it<3; it++)
		    {
		      if(partPt>thres[it]) okCh[it] = 1;
		    }
		}
	      else
		{
		  for(Int_t icl=0; icl<3; icl++)
		    {
		      if(partPt>thres[icl]) okNe[icl] = 1;
		    }
		}
	    }
	  for(Int_t i=0; i<3; i++)
	    {
	      if(okCh[i]==1)
		fhJetPtWithTrkThres[type][r]->Fill(i,jetPt);
	      if(okNe[i]==1)
		fhJetPtWithClsThres[type][r]->Fill(i,jetPt);
	    }
	  for(Int_t k=0; k<2; k++) fhJetPtVsLowPtCons[type][r][k]->Fill(jetPt,lowPt[k]/jetPt);
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::AnalyzeSecondaryContribution(AliFJWrapper *jetFinder, const Int_t r, const Int_t etaCut)
{
  //
  // Analyze secondaries
  //

  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();
  for(UInt_t ij=0; ij<jetsIncl.size(); ij++)
    {
      Float_t jetEta = jetsIncl[ij].eta();
      Float_t jetPt  = jetsIncl[ij].perp();
      if(TMath::Abs(jetEta)>0.5*etaCut+0.5 || jetPt<1) continue;
      std::vector<fastjet::PseudoJet> constituents = jetFinder->GetJetConstituents(ij);
      Double_t dNKPt = 0, dWPPt = 0, weakPt=0, nkPt=0;
      Int_t type = -1;
      for(UInt_t ic=0; ic<constituents.size(); ic++)
	{
	  Int_t ipart = TMath::Abs(constituents[ic].user_index())-1;
	  AliMCParticle* mcParticle = (AliMCParticle*)fMC->GetTrack(ipart);
	  if(!mcParticle) continue;
	  Int_t pdg = mcParticle->PdgCode();
	  if(pdg==310 || (pdg<=3400 && pdg>=3100)) weakPt += mcParticle->Pt();
	  if(pdg==130 || pdg==2112)  nkPt += mcParticle->Pt();

	  // Weak decaying particles
	  if(pdg==310) type=0;
	  else if (pdg==3122) type=1;
	  else if (pdg<=3400 && pdg>=3100) type=2;
	  else type=-1;
	  if(type>-1)
	    {
	      Int_t binx = fhSecondaryResponse[type]->GetXaxis()->FindFixBin(mcParticle->Pt());
	      TH1F *htmp = (TH1F*)fhSecondaryResponse[type]->ProjectionY(Form("hpro_%d_%d_%d",(Int_t)Entry(),ij,ic),binx,binx);
	      Double_t pro = htmp->GetRandom(); 
	      dWPPt += pro*mcParticle->Pt() - mcParticle->Pt();
	      delete htmp;
	    }

	  // Missing neutron and K0L
	  if(pdg==130 || pdg==2112) dNKPt -= mcParticle->Pt();
	}

      fhNKFracVsJetPt[etaCut][r]->Fill(jetPt,nkPt/jetPt);
      fhWeakFracVsJetPt[etaCut][r]->Fill(jetPt,weakPt/jetPt);

      Double_t jetNewWPPt = jetPt + dWPPt;
      fhJetResponseWP[etaCut][r]->Fill(jetPt,jetNewWPPt);
      fhJetResolutionWP[etaCut][r]->Fill(jetPt,(jetPt-jetNewWPPt)/jetPt);

      Double_t jetNewNKPt = jetPt + dNKPt;
      fhJetResponseNK[etaCut][r]->Fill(jetPt,jetNewNKPt);
      fhJetResolutionNK[etaCut][r]->Fill(jetPt,(jetPt-jetNewNKPt)/jetPt);
    }
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::AnalyzeSecondaryContributionViaMatching(AliFJWrapper *jetFinder, const Int_t r, const Int_t type, const Int_t etaCut)
{
  //
  // Estimate the contribution of missing energies via matching
  // type = 0 -> NK
  // type = 1 -> WP
  //

  if(type!=0 && type!=1) return;

  AliFJWrapper fj("fj","fj");
  fj.CopySettingsFrom(*jetFinder);

  Int_t npart = fMC->GetNumberOfTracks();
  Int_t pType = -1;
  for(Int_t ipart=0; ipart<npart; ipart++)
    {
      AliVParticle* vParticle = fMC->GetTrack(ipart);
      if(!IsGoodMcPartilce(vParticle,ipart)) continue;
      Int_t pdgCode = vParticle->PdgCode();
      Double_t pt = vParticle->Pt();
      if( type==0 && (pdgCode==130 || pdgCode==2112) ) continue; // Reject NK
      if(type==1)
	{
	  if(pdgCode==310) pType=0;
	  else if (pdgCode==3122) pType=1;
	  else if (pdgCode<=3400 && pdgCode>=3100) pType=2;
	  else pType=-1;
	  if(pType>-1)
	    {
	      Int_t binx = fhSecondaryResponse[pType]->GetXaxis()->FindFixBin(vParticle->Pt());
	      TH1F *htmp = (TH1F*)fhSecondaryResponse[pType]->ProjectionY(Form("hpro_%d_%d",(Int_t)Entry(),ipart),binx,binx);
	      Double_t pro = htmp->GetRandom(); 
	      pt = pro * vParticle->Pt();
	      delete htmp;
	    }
	}
      Int_t index = 1;
      if(vParticle->Charge()==0) { index=-1; }     
      fj.AddInputVector(vParticle->Px()*pt/vParticle->Pt(), vParticle->Py()*pt/vParticle->Pt(), vParticle->Pz(), vParticle->P(), (ipart+1)*index);
    }

  fj.Run(); 
  std::vector<fastjet::PseudoJet> jets_incl_mc = fj.GetInclusiveJets(); 
  std::vector<fastjet::PseudoJet> jets_incl_mc_std = jetFinder->GetInclusiveJets(); 

  for(UInt_t ij=0; ij<jets_incl_mc.size(); ij++)
    {
      Float_t jetEta = jets_incl_mc[ij].eta();
      Float_t jetPt  = jets_incl_mc[ij].perp();
      if(TMath::Abs(jetEta)>0.5*etaCut+0.5 || jetPt<5) continue;
      Double_t dEta, dPhi;
      Int_t index = FindSpatialMatchedJet(jets_incl_mc[ij], jetFinder, dEta, dPhi, 1);
      if(index>-1)
	{
	  Double_t jetPtStd = jets_incl_mc_std[index].perp();
	  Double_t dR = TMath::Sqrt(dEta*dEta+dPhi*dPhi);
	  if(type==0) fhJetResolutionNKSM[etaCut][r]->Fill(jetPtStd,(jetPtStd-jetPt)/jetPtStd,dR);
	  if(type==1) fhJetResolutionWPSM[etaCut][r]->Fill(jetPtStd,(jetPtStd-jetPt)/jetPtStd,dR);

	  if(dR<kdRCut[r])
	    {
	      if(type==0) fhJetResponseNKSM[etaCut][r]->Fill(jetPtStd,jetPt);
	      if(type==1) fhJetResponseWPSM[etaCut][r]->Fill(jetPtStd,jetPt);
	    }
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::RunAnalyzeUE(AliFJWrapper *jetFinder, const Int_t type, const Bool_t isMCTruth)
{
  //
  // Run analysis to estimate the underlying event
  // Adapted from Oliver

  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();
  Double_t leadpt=0, leadPhi = -999, leadEta = -999;
  Int_t leadIndex = -1;
  for(UInt_t ij=0; ij<jetsIncl.size(); ij++)
    {
      Double_t jetEta = jetsIncl[ij].eta();
      Double_t jetPt = jetsIncl[ij].perp();
      Double_t jetPhi = jetsIncl[ij].phi();
      if(leadpt<jetPt)
	{
	  leadpt = jetPt;
	  leadPhi = jetPhi;
	  leadEta = jetEta;
	  leadIndex = ij;
	}
    }
  if(leadpt<1 || TMath::Abs(leadEta)>0.5 ) return;

  Int_t backIndex = -1;
  Bool_t isBackToBack = kFALSE;
  for(UInt_t ij=0; ij<jetsIncl.size(); ij++)
    {
      Double_t dPhi = TMath::Abs(jetsIncl[ij].phi()-leadPhi);
      if(dPhi > kPI) dPhi = 2*kPI - dPhi;
      if(dPhi>150*TMath::DegToRad() && jetsIncl[ij].perp()/leadpt > 0.8)
	{
	  backIndex = ij;
	  isBackToBack = kTRUE;
	}
    }

  for(Int_t ij=0; ij<(Int_t)jetsIncl.size(); ij++)
    {
      if(ij==leadIndex || ij==backIndex) continue;
      if(TMath::Abs(jetsIncl[ij].eta())>0.5) continue;
      if(jetsIncl[ij].perp()>15) isBackToBack = kFALSE;
    }


  Double_t axis[2] = {leadPhi+0.5 * kPI, leadPhi-0.5 * kPI};
  if(axis[0]>2*kPI) axis[0] -= 2*kPI;
  if(axis[1]<0) axis[1] += 2*kPI;

  fhUEJetPtNorm[type][0][0]->Fill(leadpt);
  if(isBackToBack)  fhUEJetPtNorm[type][0][1]->Fill(leadpt);

  if(!isMCTruth)
    {
      Int_t ntracks = fTrackArray->GetEntriesFast();
      Double_t vertex[3] = {0, 0, 0};
      fESD->GetVertex()->GetXYZ(vertex) ;
      TLorentzVector gamma;
      Int_t nclusters = fClusterArray->GetEntriesFast();
      
      for(Int_t j=0; j<2; j++)
	{
	  Double_t ueCh = 0, ueChNe = 0.;
	  for(Int_t it=0; it<ntracks; it++)
	    {
	      AliESDtrack *t = (AliESDtrack*) fTrackArray->At(it);
	      if(!t) continue;
	      Double_t dPhi = TMath::Abs(axis[j]-t->Phi());
	      Double_t dEta = TMath::Abs(leadEta-t->Eta());
	      if(dPhi > kPI) dPhi = 2*kPI - dPhi;
	      if(TMath::Sqrt(dPhi*dPhi+dEta*dEta)<0.4) 
		{
		  fhUEJetPtVsConsPt[type][0][0]->Fill(leadpt,t->Pt());
		  if(isBackToBack) fhUEJetPtVsConsPt[type][0][1]->Fill(leadpt,t->Pt());
		  ueCh += t->Pt();
		  if(TMath::Abs(leadEta-0.4)<0.7 && axis[j]>80*TMath::DegToRad()+0.4 && axis[j]<180*TMath::DegToRad()-0.4)
		    {
		      fhUEJetPtVsConsPt[type][1][0]->Fill(leadpt,t->Pt());
		      if(isBackToBack) fhUEJetPtVsConsPt[type][1][1]->Fill(leadpt,t->Pt());
		      ueChNe  += t->Pt();
		    }
		}
	    }
	  fhUEJetPtVsSumPt[type][0][0]->Fill(leadpt,ueCh);
	  if(isBackToBack) fhUEJetPtVsSumPt[type][0][1]->Fill(leadpt,ueCh);
	  
	  if(!(TMath::Abs(leadEta-0.4)<0.7 && axis[j]>80*TMath::DegToRad()+0.4 && axis[j]<180*TMath::DegToRad()-0.4)) continue;
	  fhUEJetPtNorm[type][1][0]->Fill(leadpt);
	  if(isBackToBack) fhUEJetPtNorm[type][1][1]->Fill(leadpt);
	  for(Int_t ic=0; ic<nclusters; ic++)
	    {
	      AliESDCaloCluster * cl = (AliESDCaloCluster *)fClusterArray->At(ic);
	      if(!cl) continue;
	      cl->GetMomentum(gamma, vertex);
	      Double_t clsPhi = gamma.Phi();
	      if(clsPhi<0) clsPhi += 2*kPI;
	      Double_t dPhi = TMath::Abs(axis[j]-clsPhi);
	      Double_t dEta = TMath::Abs(leadEta-gamma.Eta());
	      if(dPhi > kPI) dPhi = 2*kPI - dPhi;
	      if(TMath::Sqrt(dPhi*dPhi+dEta*dEta)<0.4) 
		{
		  ueChNe  += gamma.Pt();
		  fhUEJetPtVsConsPt[type][1][0]->Fill(leadpt,gamma.Pt());
		  if(isBackToBack) fhUEJetPtVsConsPt[type][1][1]->Fill(leadpt,gamma.Pt());
		}
	    }
	  fhUEJetPtVsSumPt[type][1][0]->Fill(leadpt,ueChNe);
	  if(isBackToBack)  fhUEJetPtVsSumPt[type][1][1]->Fill(leadpt,ueChNe); 
	}
    }
  else
    {
      fhUEJetPtNorm[type][1][0]->Fill(leadpt);
      if(isBackToBack) fhUEJetPtNorm[type][1][1]->Fill(leadpt);
      Int_t npart = fMcPartArray->GetSize();
      for(Int_t j=0; j<2; j++)
	{
	  Double_t ueCh = 0, ueChNe = 0., ueChNe2 = 0.;
	  for(Int_t ipos=0; ipos<npart; ipos++)
	    {
	      AliVParticle* vParticle = fMC->GetTrack(fMcPartArray->At(ipos));
	      if(!vParticle) continue;
	      Double_t dPhi = TMath::Abs(axis[j]-vParticle->Phi());
	      Double_t dEta = TMath::Abs(leadEta-vParticle->Eta());
	      if(dPhi > kPI) dPhi = 2*kPI - dPhi;
	      if(TMath::Sqrt(dPhi*dPhi+dEta*dEta)<0.4) 
		{
		  Double_t pt = vParticle->Pt();
		  ueChNe += pt;
		  fhUEJetPtVsConsPt[type][1][0]->Fill(leadpt,pt);
		  if(isBackToBack)  fhUEJetPtVsConsPt[type][1][1]->Fill(leadpt,pt);
		  if( vParticle->Charge()!=0 ) 
		    {
		      fhUEJetPtVsConsPt[type][0][0]->Fill(leadpt,pt);
		      if(isBackToBack) fhUEJetPtVsConsPt[type][0][1]->Fill(leadpt,pt);
		      ueCh += pt;
		      ueChNe2 += pt;
		    }
		  else
		    {
		      if(pt>0.5) ueChNe2 += pt;
		    }
		}
	    }
	  fhUEJetPtVsSumPt[type][0][0]->Fill(leadpt,ueCh);
	  fhUEJetPtVsSumPt[type][1][0]->Fill(leadpt,ueChNe);
	  if(isBackToBack) fhUEJetPtVsSumPt[type][0][1]->Fill(leadpt,ueCh);
	  if(isBackToBack) fhUEJetPtVsSumPt[type][1][1]->Fill(leadpt,ueChNe);
	}
    }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsGoodJet(AliAODJet *jet, Double_t rad)
{
  // 
  // Check if it is a good jet 
  //

  if(jet->Pt()<1) return kFALSE;
  if(TMath::Abs(jet->Eta())>(0.7-rad)) return kFALSE;
  if(jet->Phi() < (80*TMath::DegToRad()+rad) || jet->Phi() > (180*TMath::DegToRad()-rad) ) return kFALSE;
  return kTRUE;
}


//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetLeadingZ(const Int_t jetIndex, AliFJWrapper *jetFinder)
{
  //
  // Get the leading z of the input jet
  //

  Double_t z = 0;
  std::vector<fastjet::PseudoJet> constituents = jetFinder->GetJetConstituents(jetIndex);
  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();
  Int_t index = -1;
  Double_t maxPt = 0;
  for(UInt_t ic=0; ic<constituents.size(); ic++)
    {
      if(constituents[ic].perp()>maxPt)
	{
	  maxPt = constituents[ic].perp();
	  index = ic;
	}
    }
  if(index>-1)
    z = GetZ(constituents[index].px(),constituents[index].py(),constituents[index].pz(),jetsIncl[jetIndex].px(),jetsIncl[jetIndex].py(),jetsIncl[jetIndex].pz());
  return z;
}

//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const
{
  // 
  // Get the z of a constituent inside of a jet
  //

  return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/(jetPx*jetPx+jetPy*jetPy+jetPz*jetPz);
}

//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetMeasuredJetPtResolution(const Int_t jetIndex, AliFJWrapper *jetFinder)
{
  //
  // Get jet energy resoultion due to intrinsic detector effects
  //

  Double_t jetSigma2 = 0;
  std::vector<fastjet::PseudoJet> constituents = jetFinder->GetJetConstituents(jetIndex);
  for(UInt_t ic=0; ic<constituents.size(); ic++)
    {
      if(constituents[ic].user_index()>0)
	{
	  AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents[ic].user_index()-1);
	  jetSigma2 += TMath::Power(track->Pt()*track->Pt()*TMath::Sqrt(track->GetSigma1Pt2()),2);
	}
      else
	{
	  AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(constituents[ic].user_index()*(-1)-1);
	  jetSigma2 += TMath::Power(cluster->GetTOF(),2);
	}
    }
  return TMath::Sqrt(jetSigma2);
}



//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetJetMissPtDueToTrackingEfficiency(const Int_t jetIndex, AliFJWrapper *jetFinder, const Int_t radiusIndex)
{
  //
  // Correct the tracking inefficiency explicitly
  //

  Double_t misspt = 0;
  std::vector<fastjet::PseudoJet> jetsIncl = jetFinder->GetInclusiveJets();
  Double_t jetPt = jetsIncl[jetIndex].perp();
  if(!fhCorrTrkEffPtBin[fTriggerType][radiusIndex])
    {
      printf("Warning: can't get the mean # of tracks per jet with pt=%f in: trigger=%d, radiusIndex=%d\n",jetPt,fTriggerType,radiusIndex);
      return 0;
    }
  Int_t ibin = fhCorrTrkEffPtBin[fTriggerType][radiusIndex]->FindFixBin(jetPt);
  if(!fhCorrTrkEffSample[fTriggerType][radiusIndex][ibin-1] || fhCorrTrkEffSample[fTriggerType][radiusIndex][ibin-1]->Integral()<0.001)
    {
      printf("Warning: no sampling distrubtion for jet with pt=%f\n",jetPt);
      return 0;
    }
  Double_t ntrack = fhCorrTrkEffPtBin[fTriggerType][radiusIndex]->GetBinContent(ibin);
  Int_t nTrk = (Int_t) ntrack;
  Double_t res = ntrack-nTrk;
  Double_t pro1 = fRandomGen->Uniform();
  if(pro1<res) nTrk++;
  for(Int_t itry=0; itry<nTrk; itry++)
    {
      Double_t trkPt = fhCorrTrkEffSample[fTriggerType][radiusIndex][ibin-1]->GetRandom();
      if(trkPt/jetPt>fTrkEffCorrCutZ) continue;
      Double_t eff = GetTrkEff(trkPt);
      Double_t pro = fRandomGen->Uniform();
      if(pro>eff) misspt += trkPt;
    }
  return misspt;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsGoodMcPartilce(const AliVParticle* vParticle, const Int_t ipart)
{
  //
  // Select good primary particles to feed into jet finder
  //

  if(!vParticle) return kFALSE;
  if(!fMC->IsPhysicalPrimary(ipart)) return kFALSE;
  if (TMath::Abs(vParticle->Eta())>2) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskFullppJet::FindSpatialMatchedJet(fastjet::PseudoJet jet, AliFJWrapper *jetFinder, Double_t &dEta, Double_t &dPhi, Double_t maxR)
{
  //
  // Find spatially matched detector and particle jets
  //

  Int_t index=-1;
  dEta=-999, dPhi=-999;
  std::vector<fastjet::PseudoJet> jets_incl = jetFinder->GetInclusiveJets();  
  for(UInt_t ij=0; ij<jets_incl.size(); ij++)
    {
      if(jets_incl[ij].perp()<5) continue;
      if(TMath::Abs(jets_incl[ij].eta())>1) continue;
      Double_t tmpR = TMath::Sqrt( TMath::Power(jet.eta()-jets_incl[ij].eta(), 2) + TMath::Power(jet.phi()-jets_incl[ij].phi(), 2) );
      if(tmpR<maxR)
	{
	  maxR=tmpR;
	  index=ij;
	  dEta = jet.eta()-jets_incl[ij].eta();
	  dPhi = jet.phi()-jets_incl[ij].phi();
	}
    }
  return index;
}

//________________________________________________________________________
Int_t AliAnalysisTaskFullppJet::FindEnergyMatchedJet(AliFJWrapper *jetFinder1, const Int_t index1, AliFJWrapper *jetFinder2, const Double_t fraction)
{
  //
  // Find matched detector and particle jets based on shared constituents
  //

  Int_t matchedIndex = -1;
  std::vector<fastjet::PseudoJet> jetsIncl1 = jetFinder1->GetInclusiveJets();  
  std::vector<fastjet::PseudoJet> jetsIncl2 = jetFinder2->GetInclusiveJets();
  std::vector<fastjet::PseudoJet> constituents1 = jetFinder1->GetJetConstituents(index1);
  Double_t jetPt1 = jetsIncl1[index1].perp();
  if(jetPt1<0) return matchedIndex;

  for(UInt_t ij2=0; ij2<jetsIncl2.size(); ij2++)
    {
      Double_t jetPt2  = jetsIncl2[ij2].perp();
      if(jetPt2<0) return matchedIndex;
      std::vector<fastjet::PseudoJet> constituents2 = jetFinder2->GetJetConstituents(ij2);
      Double_t sharedPt1 = 0., sharedPt2 = 0.;
      for(UInt_t ic2=0; ic2<constituents2.size(); ic2++)
	{
	  Int_t mcLabel = constituents2[ic2].user_index()-1;
	  Double_t consPt2 = constituents2[ic2].perp();
	  for(UInt_t ic1=0; ic1<constituents1.size(); ic1++)
	    {
	      Double_t consPt1 = constituents1[ic1].perp();
	      if(constituents1[ic1].user_index()>0)
		{
		  AliESDtrack *track = (AliESDtrack*) fTrackArray->At(constituents1[ic1].user_index()-1);
		  if(track->GetLabel()==mcLabel)  {sharedPt2 += consPt2; sharedPt1 += consPt1; cout<<"Found a matched track"<<endl;break;}
		}
	      else
		{
		  AliESDCaloCluster *cluster = (AliESDCaloCluster *)fClusterArray->At(constituents1[ic1].user_index()*(-1)-1);
		  if(cluster->GetLabel()==mcLabel) {sharedPt2 += consPt2; sharedPt1 += consPt1; cout<<"Found a matched cluster"<<endl;break;}
		}
	    }
	}
      cout<<sharedPt1/jetPt1<<"  "<<sharedPt2/jetPt2<<endl;
      if(sharedPt1/jetPt1 > fraction && sharedPt2/jetPt2 > fraction)
	{
	  matchedIndex = ij2;
	  break;
	} 
    }
  return matchedIndex;
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::GetESDTrax()
{
  //
  // Get esd tracks.
  //

  Int_t nTrax = fESD->GetNumberOfTracks();
  if(fVerbosity>5) printf("[i] # of tracks in event: %d\n",nTrax);

  for (Int_t i = 0; i < nTrax; ++i) 
    {
      AliESDtrack* esdtrack = fESD->GetTrack(i);
      if (!esdtrack) 
	{
	  AliError(Form("Couldn't get ESD track %d\n", i));
	  continue;
	}

      AliESDtrack *newtrack = GetAcceptTrack(esdtrack);
      if(!newtrack) continue;

      Double_t pt  = newtrack->Pt();
      Double_t eta = newtrack->Eta();
      if (pt < fTrkPtMin[fTriggerType] || pt > fTrkPtMax[fTriggerType] || TMath::Abs(eta) > fTrkEtaMax)
	{
	  delete newtrack;
	  continue;
	}

      if(fSysTrkEff)
	{
	  Double_t rand = fRandomGen->Uniform();
	  if(rand<fVaryTrkEff)
	    {
	      delete newtrack;
	      continue;
	    }
	}
      newtrack->SetIntegratedLength(esdtrack->Pt());
      fTrackArray->Add(newtrack);
    }
  if(fVerbosity>5) printf("[i] # of tracks in event: %d\n", fTrackArray->GetEntries());
  fNTracksPerChunk += fTrackArray->GetEntries();
}

//
//________________________________________________________________________
//
AliESDtrack *AliAnalysisTaskFullppJet::GetAcceptTrack(AliESDtrack *esdtrack)
{
  //
  // Get the hybrid tracks
  //

  AliESDtrack *newTrack = 0x0;

  if(fTrackCutsType==0 || fTrackCutsType==3)
    {
      if(fEsdTrackCuts->AcceptTrack(esdtrack))
	{
	  newTrack = new AliESDtrack(*esdtrack);
	  newTrack->SetTRDQuality(0);
	}
      else if(fHybridTrackCuts1->AcceptTrack(esdtrack))
	{
	  if(esdtrack->GetConstrainedParam())
	    {
	      newTrack = new AliESDtrack(*esdtrack);
	      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
	      newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
	      newTrack->SetTRDQuality(1);		
	    }
	  else 
	    return 0x0;
	}
      else if(fHybridTrackCuts2->AcceptTrack(esdtrack))
	{
	  if(esdtrack->GetConstrainedParam())
	    {
	      newTrack = new AliESDtrack(*esdtrack);
	      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
	      newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
	      newTrack->SetTRDQuality(2);		
	    }
	  else 
	    return 0x0;
	}
      else
	{
	  return 0x0;
	}
    }
  else if(fTrackCutsType==1)
    {
      if(fEsdTrackCuts->AcceptTrack(esdtrack))
	{
	  newTrack = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(fESD),esdtrack->GetID());// use TPC only tracks with non default SPD vertex
	  if(!newTrack) return 0x0;  
	  AliExternalTrackParam exParam;
	  Bool_t relate = newTrack->RelateToVertexTPC(fESD->GetPrimaryVertexSPD(),fESD->GetMagneticField(),kVeryBig,&exParam); //constrain to SPD vertex
	  if( !relate )
	    {
	      delete newTrack;
	      return 0x0;
	    }
	  newTrack->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
	  newTrack->SetTRDQuality(1);
	}
      else
	{
	  return 0x0;
	}
    }
  else if (fTrackCutsType==2)
    {
      if(fEsdTrackCuts->AcceptTrack(esdtrack))
	{
	  newTrack = new AliESDtrack(*esdtrack);
	  newTrack->SetTRDQuality(0);
	}
      else
	return 0x0;
    }
  else
    {
      printf("Unknown track cuts type: %d\n",fTrackCutsType);
      return 0x0;
    }

  return newTrack;
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::GetESDEMCalClusters()
{
  //
  // Get emcal clusters - selected
  //

  if(fSysTrkClsMth)
    {
      if (!TGeoGlobalMagField::Instance()->GetField()) fESD->InitMagneticField();
      fRecoUtil->FindMatches(fESD,0x0,fGeom);
      fRecoUtil->SetClusterMatchedToTrack(fESD);
      fRecoUtil->SetTracksMatchedToCluster(fESD);
    }

  const Int_t nCaloClusters = fESD->GetNumberOfCaloClusters();
  for(Int_t i = 0 ; i < nCaloClusters; i++) 
    {
      AliESDCaloCluster * cl = (AliESDCaloCluster *) fESD->GetCaloCluster(i);
      if(!IsGoodCluster(cl)) continue;
      AliESDCaloCluster *newCluster = new AliESDCaloCluster(*cl);
      if (newCluster->E() < fClsEtMin[fTriggerType] || newCluster->E() > fClsEtMax[fTriggerType]) 
	{delete newCluster; continue;}

      // Absolute scale
      if(fPeriod.Contains("lhc11a",TString::kIgnoreCase)) 
	{
	  Double_t newE = newCluster->E() * 1.02;
	  newCluster->SetE(newE);
	}

      // Trigger efficiency systematic uncertainty
      if(fSysJetTrigEff)
	{
	  Double_t newE = newCluster->E() * (1+fVaryJetTrigEff);
	  newCluster->SetE(newE);
	  if(fTriggerType==0) fhClsE[fTriggerType]->Fill(newCluster->E());
	  if(fTriggerType==1)
	    {
	      if(fPeriod.Contains("lhc12a15a",TString::kIgnoreCase)) fhClsE[0]->Fill(newCluster->E());
	      if(newCluster->Chi2()>0.5) fhClsE[fTriggerType]->Fill(newCluster->E());
	    }
	}

      // Cluster energy scale systematic uncertainty
      if(fSysClusterEScale)
	{
	  Double_t newE = newCluster->E() * (1+fVaryClusterEScale);
	  newCluster->SetE(newE);
	}

      // Cluster energy resolution systematic uncertainty
      if(fSysClusterERes)
        {
          Double_t oldE = newCluster->E();
          Double_t resolution = fClusterEResolution->Eval(oldE);
          Double_t smear = resolution * TMath::Sqrt((1+fVaryClusterERes)*(1+fVaryClusterERes)-1);
          Double_t newE = oldE + fRandomGen->Gaus(0, smear);
          newCluster->SetE(newE);
        }

      // non-linearity systematic uncertainty
      if(fSysNonLinearity)
	{
	  Double_t oldE = newCluster->E();
	  Double_t newE = oldE/fNonLinear->Eval(oldE) * 1.012;
	  newCluster->SetE(newE);
	}

      // clusterizer systematic uncertainty
      if(fSysClusterizer)
	{
	  fhSysClusterE[fTriggerType][0]->Fill(newCluster->E());
	  fhSysNCellVsClsE[fTriggerType][0]->Fill(newCluster->E(),newCluster->GetNCells());
	}

      // hadronic correction
      Double_t clsE = newCluster->E();
      Double_t subE = 0., eRes = 0., mcSubE = 0, MIPE = 0.;
      if(fElectronRejection || fHadronicCorrection)
	subE= SubtractClusterEnergy(newCluster,  eRes, MIPE, mcSubE);

      if(!fStudySubEInHC && !fStudyMcOverSubE)  clsE -= subE;
      if(clsE<0) {delete newCluster; continue;}
      newCluster->SetE(clsE);
      newCluster->SetTOF(eRes);
      newCluster->SetDispersion(subE);
      newCluster->SetEmcCpvDistance(MIPE);
      newCluster->SetDistanceToBadChannel(mcSubE);
      fClusterArray->Add(newCluster);

      // clusterizer systematic uncertainty
      if(fSysClusterizer)
	{
	  fhSysClusterE[fTriggerType][1]->Fill(newCluster->E());
	  fhSysNCellVsClsE[fTriggerType][1]->Fill(newCluster->E(),newCluster->GetNCells());
	}
    }

  if(fVerbosity>5) printf("[i] # of EMCal clusters in event: %d\n", fClusterArray->GetEntries());

}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsGoodCluster(AliESDCaloCluster *cluster)
{
  //
  // Select good clusters
  //

  if(!cluster) return kFALSE;
  if (!cluster->IsEMCAL()) return kFALSE;
  if(fRejectExoticCluster && fRecoUtil->IsExoticCluster(cluster, (AliVCaloCells*)fESD->GetEMCALCells())) return kFALSE;
  if(fRemoveBadChannel && fRecoUtil->ClusterContainsBadChannel(fGeom, cluster->GetCellsAbsId(),cluster->GetNCells())) return kFALSE;
  return kTRUE;
}


//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::SubtractClusterEnergy(AliESDCaloCluster *cluster, Double_t &eRes, Double_t &MIPE, Double_t &mcSubE)
{
  //
  // Hadronic correction
  //

  mcSubE = 0;
  eRes = 0; 
  MIPE = 0;

  eRes += TMath::Power(fClusterEResolution->Eval(cluster->E()),2);
  Double_t subE = 0., sumTrkPt = 0., sumTrkP = 0.;
  Int_t nTrack = 0;
  TArrayI *matched = cluster->GetTracksMatched();
  if(!matched) return 0;
  for(Int_t im=0; im<matched->GetSize(); im++)
    {
      Int_t trkIndex = matched->At(im);
      if(trkIndex<0 || trkIndex>=fESD->GetNumberOfTracks()) continue;
      Bool_t isSelected = kFALSE;
      Int_t index = -1;
      for(Int_t j=0; j<fTrackArray->GetEntriesFast(); j++)
	{
	  AliESDtrack *tr = (AliESDtrack*)fTrackArray->At(j);
	  if( trkIndex == tr->GetID() )
	    {
	      isSelected = kTRUE;
	      index = j;
	      break;
	    }
	}
      if(!isSelected) continue;
      nTrack++;
      AliESDtrack *track = (AliESDtrack*)fTrackArray->At(index);
      Double_t trkP = track->P();
      sumTrkPt += track->Pt(); sumTrkP += track->P();
      if(fSysTrkPtRes) trkP = TMath::Sqrt(TMath::Power(track->Pt()+GetSmearedTrackPt(track),2)+TMath::Power(track->Pz(),2));
      if(IsElectron(track,cluster->E())) //electrons
	{
	  if(fElectronRejection)
	    {
	      subE+= trkP;
	      MIPE+= trkP;
	      eRes += TMath::Power(track->Pt()*track->Pt()*TMath::Sqrt(track->GetSigma1Pt2()),2);
	    }
	}
      else //hadrons
	{
	  if(fHadronicCorrection)
	    {
	      MIPE += (trkP>0.27)?0.27:trkP; //MIP correction
	      if(fFractionHC>2)
		{
		  if(trkP>0.27)
		    {
		      subE+=0.27;
		    }
		  else
		    {
		      subE+=trkP;
		      eRes += TMath::Power(track->Pt()*track->Pt()*TMath::Sqrt(track->GetSigma1Pt2()),2);
		    }
		}
	      else
		{
		  if(trkP>fHCLowerPtCutMIP) subE += 0.27;
		  else subE+=fFractionHC*trkP;
		  eRes += TMath::Power(fFractionHC*track->Pt()*track->Pt()*TMath::Sqrt(track->GetSigma1Pt2()),2);
		}
	    }
	}
    }

  if(fSaveQAHistos) fhNMatchedTrack[fTriggerType]->Fill(nTrack);
  if(nTrack>0)
    {
      Double_t fraction[4] = {0.3,0.5,0.7,1.0};
      for(Int_t j=0; j<4; j++)
	{
	  Double_t subETmp = sumTrkP*fraction[j];
	  if(subETmp>cluster->E()) subETmp = cluster->E();
	  if(fSaveQAHistos) fhSubEVsTrkPt[fTriggerType][j]->Fill(sumTrkPt,subETmp/sumTrkP);
	}
    }

  eRes = TMath::Sqrt(eRes);

  if(fIsMC && nTrack>0)
    {
      Double_t neutralE = 0;
      TArrayI* labels = cluster->GetLabelsArray();
      if(labels)
	{
	  for(Int_t il=0; il<labels->GetSize(); il++)
	    {
	      Int_t ipart = labels->At(il);
	      if(ipart>-1 && ipart<fMC->GetNumberOfTracks())
		{
		  AliVParticle* vParticle = fMC->GetTrack(ipart);
		  if(vParticle->Charge()==0)
		    {
		      neutralE += vParticle->E();
		    }
		}
	    }
	}
      mcSubE = cluster->E() - neutralE;
      if(mcSubE<0)
	mcSubE=0;
    }


  return subE;
}


//
//________________________________________________________________________
//
Double_t AliAnalysisTaskFullppJet::GetExoticEnergyFraction(AliESDCaloCluster *cluster)
{
  // 
  // Exotic fraction: f_cross
  // Adpated from AliEMCalRecoUtils

  if(!cluster) return -1;
  AliVCaloCells *cells = (AliVCaloCells*)fESD->GetEMCALCells();
  if(!cells)  return -1;
  
  // Get highest energy tower
  Int_t iSupMod = -1, absID = -1, ieta = -1, iphi = -1,iTower = -1, iIphi = -1, iIeta = -1; 
  Bool_t share = kFALSE;
  fRecoUtil->GetMaxEnergyCell(fGeom, cells, cluster, absID, iSupMod, ieta, iphi, share);
  fGeom->GetCellIndex(absID,iSupMod,iTower,iIphi,iIeta); 
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);

  Int_t absID1 = fGeom-> GetAbsCellIdFromCellIndexes(iSupMod, iphi+1, ieta);
  Int_t absID2 = fGeom-> GetAbsCellIdFromCellIndexes(iSupMod, iphi-1, ieta);
  Int_t absID3 = fGeom-> GetAbsCellIdFromCellIndexes(iSupMod, iphi, ieta+1);
  Int_t absID4 = fGeom-> GetAbsCellIdFromCellIndexes(iSupMod, iphi, ieta-1);
  
  Float_t  ecell  = 0, ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
  Double_t tcell  = 0, tcell1  = 0, tcell2  = 0, tcell3  = 0, tcell4  = 0;
  Bool_t   accept = 0, accept1 = 0, accept2 = 0, accept3 = 0, accept4 = 0;
  const Int_t bc  = 0;
  
  accept  = fRecoUtil->AcceptCalibrateCell(absID, bc, ecell ,tcell ,cells); 
    
  if(!accept) return -1;
  
  if(ecell < 0.5) return -1;
  
  accept1 = fRecoUtil->AcceptCalibrateCell(absID1,bc, ecell1,tcell1,cells); 
  accept2 = fRecoUtil->AcceptCalibrateCell(absID2,bc, ecell2,tcell2,cells); 
  accept3 = fRecoUtil->AcceptCalibrateCell(absID3,bc, ecell3,tcell3,cells); 
  accept4 = fRecoUtil->AcceptCalibrateCell(absID4,bc, ecell4,tcell4,cells); 
  
  const Double_t exoticCellDiffTime = 1e6;
  if(TMath::Abs(tcell-tcell1)*1.e9 > exoticCellDiffTime) ecell1 = 0 ;
  if(TMath::Abs(tcell-tcell2)*1.e9 > exoticCellDiffTime) ecell2 = 0 ;
  if(TMath::Abs(tcell-tcell3)*1.e9 > exoticCellDiffTime) ecell3 = 0 ;
  if(TMath::Abs(tcell-tcell4)*1.e9 > exoticCellDiffTime) ecell4 = 0 ;

  Float_t eCross = ecell1+ecell2+ecell3+ecell4;

  return 1-eCross/ecell;
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::GetMCInfo()
{
  //
  // Get # of trials per ESD event
  //

  AliStack *stack = fMC->Stack();
  if(stack)
    {
      AliGenEventHeader *head = dynamic_cast<AliGenEventHeader*>(fMC->GenEventHeader());
      if (head == 0x0)
       {
         AliError("Could not get the event header");
         return;
       }

     AliGenPythiaEventHeader *headPy = dynamic_cast<AliGenPythiaEventHeader*>(head);
     if (headPy != 0x0)
       {
         if (headPy->Trials() > 0)
           {
             fhNTrials[0]->Fill(0.5,headPy->Trials());
             
	   }
       }
    }
}



//________________________________________________________________________
void AliAnalysisTaskFullppJet::Terminate(Option_t *) 
{
  //
  // Called once at the end of the query
  //
  Info("Terminate","Terminate");
  AliAnalysisTaskSE::Terminate();

}

//
//________________________________________________________________________
//
Int_t AliAnalysisTaskFullppJet::RunOfflineTrigger() 
{
  //
  // Run trigger offline
  //

  fIsEventTriggerBit = 0;
  Int_t isTrigger = 0;
  Int_t ncl = fESD->GetNumberOfCaloClusters();
  for(Int_t icl=0; icl<ncl; icl++)
    {
      // Check every cluster
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(icl);
      if(!IsGoodCluster(cluster)) continue;
      Double_t pro = GetOfflineTriggerProbability(cluster);
      Double_t rand = fRandomGen->Uniform();
      if(rand<pro)
	{
	  isTrigger = 1;
	  fIsEventTriggerBit = 1;
	  cluster->SetChi2(1);
	}
      else 
	cluster->SetChi2(0);
    }
  return isTrigger;
}


//
//________________________________________________________________________
//
Double_t AliAnalysisTaskFullppJet::GetOfflineTriggerProbability(AliESDCaloCluster *cluster)
{  
  //
  // Get the probablity of the given cluster to trigger the event
  //

  Double_t pro = 0;
  // Check the trigger mask
  AliVCaloCells *cells = (AliVCaloCells*)fESD->GetEMCALCells();
  Int_t iSupMod = -1, absID = -1, ieta = -1, iphi = -1,iTower = -1, iIphi = -1, iIeta = -1; 
  Bool_t share = kFALSE;
  fRecoUtil->GetMaxEnergyCell(fGeom, cells, cluster, absID, iSupMod, ieta, iphi, share); // Get the position of the most energetic cell in the cluster
  fGeom->GetCellIndex(absID,iSupMod,iTower,iIphi,iIeta); 
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
  // convert co global phi eta
  Int_t gphi = iphi + 24*(iSupMod/2);
  Int_t geta = ieta + 48*(iSupMod%2);
  // get corresponding FALTRO
  Int_t fphi = gphi / 2;
  Int_t feta = geta / 2;
  if(fCheckTriggerMask && fTriggerMask->GetBinContent(feta+1,fphi+1)>0.5 && iSupMod>-1) // check the trigger mask
    {
      Double_t clsE = cluster->E();
      if(fSysClusterEScale) clsE = clsE * (1+fVaryClusterEScale); // Used for systematic uncertainty. Not needed for regular analysis
      if(clsE>10) pro = 1; // Probability is 1 at high E
      else
	{
	  Int_t bin = fTriggerCurve[iSupMod]->FindFixBin(clsE);
	  pro = fTriggerCurve[iSupMod]->GetBinContent(bin)/fTriggerEfficiency[iSupMod]->Eval(10); // Read the probability from trigger turn-on curves
	}
    }
  return pro;
}



//
//________________________________________________________________________
//
Int_t AliAnalysisTaskFullppJet::GetClusterSuperModule(AliESDCaloCluster *cluster)
{
  //
  // Return the given cluster supermodule
  //

  Float_t pos[3];
  cluster->GetPosition(pos);
  TVector3 clsVec(pos[0],pos[1],pos[2]);

  Int_t sMod=-1;
  fGeom->SuperModuleNumberFromEtaPhi(clsVec.Eta(), clsVec.Phi(), sMod);
  return sMod;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsElectron(AliESDtrack *track, Double_t clsE) const
{
  //
  // Check if the given track is an electron candidate based on de/dx
  //

  if(track->GetTPCsignal()<=fdEdxMax && track->GetTPCsignal()>=fdEdxMin && (clsE/track->P())<=fEoverPMax && (clsE/track->P())>=fEoverPMin )
    return kTRUE;
  else
    return kFALSE;
}


//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetTrkEff(Double_t inPt)
{
  // 
  // Get tracking efficiency estimated from simulation
  //

  Double_t eff = 1;
  Double_t ptBound[4] = {0, 0.5, 3.8, 300};

  for(Int_t i=0; i<3; i++)
    {
      if( inPt < ptBound[i+1] && inPt >= ptBound[i])
	{
	  eff = fTrkEffFunc[i]->Eval(inPt);
	  break;
	}
    }
  return eff;
}

//________________________________________________________________________
Double_t AliAnalysisTaskFullppJet::GetSmearedTrackPt(AliESDtrack *track)
{
  //
  // Smear track momentum
  //

  Double_t resolution = track->Pt()*track->Pt()*TMath::Sqrt(track->GetSigma1Pt2());
  Double_t smear = resolution*TMath::Sqrt((1+fVaryTrkPtRes)*(1+fVaryTrkPtRes)-1);
  return fRandomGen->Gaus(0, smear);

}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::CheckEventTriggerBit()
{
  //
  // Check if the triggered events have correct trigger bit
  //

 fIsEventTriggerBit = 0;
  // constants
  const Int_t nColsModule = 2;
  const Int_t nRowsModule = 5;
  const Int_t nRowsFeeModule = 24;
  const Int_t nColsFeeModule = 48;
  const Int_t nColsFaltroModule = 24;
  const Int_t nRowsFaltroModule = 12;

  // part 1, trigger extraction -------------------------------------
  Int_t trigtimes[30], globCol, globRow, ntimes;
  Int_t trigger[nColsFaltroModule*nColsModule][nRowsFaltroModule*nRowsModule];

  // erase trigger maps
  for( Int_t i = 0; i < nColsFaltroModule*nColsModule; i++ )
    {
      for( Int_t j = 0; j < nRowsFaltroModule*nRowsModule; j++ )
	{
	  trigger[i][j] = 0;
	}
    }

  Int_t fTrigCutLow = 7, fTrigCutHigh = 10, trigInCut = 0;
  AliESDCaloTrigger * fCaloTrigger = fESD->GetCaloTrigger( "EMCAL" );
  // go through triggers
  if( fCaloTrigger->GetEntries() > 0 )
    {
      // needs reset
      fCaloTrigger->Reset();
      while( fCaloTrigger->Next() )
	{
	  fCaloTrigger->GetPosition( globCol, globRow );  // get position in global 2x2 tower coordinates
	  fCaloTrigger->GetNL0Times( ntimes );   // get dimension of time arrays
	  if( ntimes < 1 ) continue;    // no L0s in this channel. Presence of the channel in the iterator still does not guarantee that L0 was produced!!
	  fCaloTrigger->GetL0Times( trigtimes );  // get timing array
	  if(fCheckTriggerMask && fTriggerMask->GetBinContent(globCol+1,globRow+1)<0.5) continue;
	  trigInCut = 0;
	  for( Int_t i = 0; i < ntimes; i++ )
	    {
	      if( trigtimes[i] > fTrigCutLow && trigtimes[i] < fTrigCutHigh )  // check if in cut
		{
		  trigInCut = 1;
		}
	    }
	  if(trigInCut==1) trigger[globCol][globRow] = 1;
	} // calo trigger entries
    } // has calo trigger entries


  
  // part 2 go through the clusters here -----------------------------------
  Int_t nCell, iCell;
  UShort_t *cellAddrs;
  Int_t absID = -1, nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1, gphi=-1, geta=-1, feta=-1, fphi=-1;
  Bool_t share = kFALSE;
  AliVCaloCells *cells = (AliVCaloCells*)fESD->GetEMCALCells();
  for(Int_t icl=0; icl<fESD->GetNumberOfCaloClusters(); icl++)
    {
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(icl);
      if(!cluster || !cluster->IsEMCAL()) continue;
      cluster->SetChi2(0);
      if(!IsGoodCluster(cluster)) continue;

      //Clusters with most energetic cell in dead region can't be triggered
      fRecoUtil->GetMaxEnergyCell(fGeom, cells, cluster, absID, nSupMod, ieta, iphi, share);
      fGeom->GetCellIndex(absID,nSupMod,nModule, nIphi, nIeta); 
      fGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule, nIphi, nIeta, iphi, ieta);
      gphi = iphi + nRowsFeeModule*(nSupMod/2);
      geta = ieta + nColsFeeModule*(nSupMod%2);
      fphi = gphi / 2;
      feta = geta / 2;

      if(fCheckTriggerMask && fTriggerMask->GetBinContent(feta+1,fphi+1)>0.5)
	{
	  nCell = cluster->GetNCells();  // get cluster cells
	  cellAddrs = cluster->GetCellsAbsId();  // get the cell addresses
	  for( iCell = 0; iCell < nCell; iCell++ )
	    {
	      // get cell position
	      fGeom->GetCellIndex( cellAddrs[iCell], nSupMod, nModule, nIphi, nIeta );
	      fGeom->GetCellPhiEtaIndexInSModule( nSupMod,nModule, nIphi, nIeta, iphi, ieta);
	      
	      // convert co global phi eta
	      gphi = iphi + nRowsFeeModule*(nSupMod/2);
	      geta = ieta + nColsFeeModule*(nSupMod%2);
	      
	      // get corresponding FALTRO
	      fphi = gphi / 2;
	      feta = geta / 2;
	      // try to match with a triggered
	      if( trigger[feta][fphi] == 1)
		{
		  cluster->SetChi2(1);
		  fIsEventTriggerBit = 1;
		  //break;
		}
	    } // cells
	}
    } // clusters
}


//________________________________________________________________________
void AliAnalysisTaskFullppJet::PrintConfig()
{
  //
  // Print configuration
  //

  const char *trackCutName[3] = {"Hybrid tracks","TPCOnly tracks","Golden tracks"};
  const char *triggerType[2]  = {"MB","EMC"};
  const char *decision[2] = {"no","yes"};
  const char *recombination[] = {"E_scheme","pt_scheme","pt2_scheme","Et_scheme","Et2_scheme","BIpt_scheme","BIpt2_scheme"};
  const char *type[2] = {"local","grid"};
  if(fStudySubEInHC || fStudyMcOverSubE)
    {
      printf("\n\n=================================\n");
      printf("======WARNING: HC is ingored!======\n");
      printf("======    NOT for PHYSICS!   ======\n\n");
    }
  printf("Run period: %s\n",fPeriod.Data());
  printf("Reject SPD pileup: %s\n",decision[fRejectPileup]);
  printf("Reject exotic triggered events: %s\n",decision[fRejectExoticTrigger]);
  printf("Is this MC data: %s\n",decision[fIsMC]);
  printf("Only find charged jets in MC: %s\n",decision[fChargedMC]);
  printf("Analyze on local or grid? %s\n",type[fAnaType]);
  printf("Run offline trigger on MC: %s\n",decision[fOfflineTrigger]);
  if(fIsMC)
    printf("Is K0 and n included: %s\n",decision[1-fRejectNK]);
  printf("Constrain tracks in EMCal acceptance: %s\n",decision[fConstrainChInEMCal]);
  printf("Track selection:    %s, |eta| < %2.1f\n",trackCutName[fTrackCutsType], fTrkEtaMax);
  for(Int_t i=0; i<2; i++)
    {
      printf("Track pt cut:       %s -> %2.2f < pT < %2.1f\n",triggerType[i], fTrkPtMin[i], fTrkPtMax[i]);
    }
  for(Int_t i=0; i<2; i++)
    {
      printf("Cluster selection:  %s -> %2.2f < Et < %2.1f\n",triggerType[i],fClsEtMin[i], fClsEtMax[i]);
    }
  printf("Electron selectoin: %2.0f < dE/dx < %2.0f, %1.1f < E/P < %1.1f\n",fdEdxMin, fdEdxMax, fEoverPMin, fEoverPMax);
  printf("Reject exotic cluster: %s\n",decision[fRejectExoticCluster]);
  printf("Remove problematic region in SM4: %s\n",decision[fRemoveBadChannel]);
  printf("Use only good SM (1,2,6,7,8,9) for trigger: %s\n",decision[fUseGoodSM]);
  printf("Reject electron: %s\n", decision[fElectronRejection]);
  printf("Correct hadron: %s\n",decision[fHadronicCorrection]);
  printf("HC fraction: %2.1f up to %2.0f GeV/c\n",fFractionHC,fHCLowerPtCutMIP);
  printf("Find charged jets: %s\n", decision[fFindChargedOnlyJet]);
  printf("Find netural jets: %s\n", decision[fFindNeutralOnlyJet]);
  printf("Find good jets: %s\n",decision[fSpotGoodJet]);
  printf("Jet radius: %s\n",fRadius.Data());
  printf("Jet recombination scheme: %s\n",recombination[fRecombinationScheme]);
  printf("Correct tracking efficiency: %s\n",decision[fCheckTrkEffCorr]);
  printf("Save jet QA histos: %s\n",decision[fSaveQAHistos]);
  printf("Systematics: jet efficiency: %s with variation %1.0f%%\n",decision[fSysJetTrigEff],fVaryJetTrigEff*100);
  printf("Systematics: tracking efficiency: %s with variation %1.0f%%\n",decision[fSysTrkEff],fVaryTrkEff*100);
  printf("Systematics: track pt resolution: %s with variation %1.0f%%\n",decision[fSysTrkPtRes],fVaryTrkPtRes*100);
  printf("Systematics: track-cluster matching: %s with |dEta|<%2.3f, |dPhi|<%2.3f\n",decision[fSysTrkClsMth],fCutdEta,fCutdPhi);
  printf("Systematics: EMCal non-linearity: %s\n",decision[fSysNonLinearity]);
  printf("Systematics: EMCal energy scale: %s with uncertainty of %1.0f%%\n",decision[fSysClusterEScale],fVaryClusterEScale*100);
  printf("Systematics: EMCal energy resolution: %s with uncertainty of %1.0f%%\n",decision[fSysClusterERes],fVaryClusterERes*100);
  printf("Smear lhc12a15a: %s\n",decision[fSmearMC]);
  printf("Run UE analysis: %s\n",decision[fRunUE]);
  printf("Run secondaries: %s\n",decision[fRunSecondaries]);
  printf("=======================================\n\n");
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::CheckExoticEvent()
{
  //
  // Check if the event containts exotic clusters
  //

  Double_t leadingE = 0;
  for(Int_t icl=0; icl<fESD->GetNumberOfCaloClusters(); icl++)
    {
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(icl);
      if(!cluster) continue;
      if (!cluster->IsEMCAL()) continue;
      if(fRecoUtil->IsExoticCluster(cluster, (AliVCaloCells*)fESD->GetEMCALCells()) && cluster->E()>leadingE)
	{
	  leadingE = cluster->E();
	}
    }
  if(leadingE>3) fIsExoticEvent3GeV = kTRUE;
  if(leadingE>5) fIsExoticEvent5GeV = kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::HasPrimaryVertex() const
{
  //
  // Check if the primary vertex exists
  //
  const AliESDVertex* vtx = fESD->GetPrimaryVertex();
  if (!vtx || vtx->GetNContributors()<1) return kFALSE;
  else return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsPrimaryVertexOk() const
{
  //
  // Check if the event vertex is good
  //
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsTPCOnlyVtx() const
{
  //
  // Check if the event only has valid TPC vertex
  //

  const Bool_t isPrimaryVtx = HasPrimaryVertex();

  Bool_t isGoodVtx = kTRUE;
  AliESDVertex* goodvtx = const_cast<AliESDVertex*>(fESD->GetPrimaryVertexTracks());
  if(!goodvtx || goodvtx->GetNContributors()<1)
    goodvtx = const_cast<AliESDVertex*>(fESD->GetPrimaryVertexSPD()); // SPD vertex
  if(!goodvtx || !goodvtx->GetStatus()) isGoodVtx = kFALSE;  // Event rejected

  if( isPrimaryVtx && !isGoodVtx )
    return kTRUE;
  else
    return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskFullppJet::CheckTPCOnlyVtx(const UInt_t trigger)
{
  //
  // Check the fraction of accepted events that have only TPC vertex
  //
  Int_t lTriggerType = -1;
  if (trigger & AliVEvent::kMB)       lTriggerType = 1;
  if (trigger & AliVEvent::kFastOnly) lTriggerType = 0;
  if (trigger & AliVEvent::kEMC1)     lTriggerType = 2;
  if(lTriggerType==-1) return;
  fhEventStatTPCVtx->Fill(0.5+lTriggerType*3);
  if(HasPrimaryVertex()) fhEventStatTPCVtx->Fill(1.5+lTriggerType*3);
  if(IsTPCOnlyVtx()) fhEventStatTPCVtx->Fill(2.5+lTriggerType*3);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFullppJet::IsLEDEvent() const
{
  // 
  // Check if the event is contaminated by LED signal
  //

  AliESDCaloCells *cells = fESD->GetEMCALCells();
  Short_t nCells = cells->GetNumberOfCells();
  Int_t nCellCount[2] = {0,0};
  for(Int_t iCell=0; iCell<nCells; iCell++)
    {
      Int_t cellId = cells->GetCellNumber(iCell);
      Double_t cellE = cells->GetCellAmplitude(cellId);
      Int_t sMod = fGeom->GetSuperModuleNumber(cellId);
      
      if(sMod==3 || sMod==4)
	{
	  if(cellE>0.1)
	    nCellCount[sMod-3]++;
	}
    }
  Bool_t isLED=kFALSE;

  if(fPeriod.CompareTo("lhc11a")==0)
    {
      if (nCellCount[1] > 100)
	isLED = kTRUE;
      Int_t runN = fESD->GetRunNumber();
      if (runN>=146858 && runN<=146860)
	{
	  if(fTriggerType==0 && nCellCount[0]>=21) isLED=kTRUE;
	  if(fTriggerType==1 && nCellCount[0]>=35) isLED=kTRUE;
	}
    }

  return isLED;
}
