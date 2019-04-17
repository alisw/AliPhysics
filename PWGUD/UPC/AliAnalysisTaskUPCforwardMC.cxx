/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODVertex.h"


// my headers
#include "AliAnalysisTaskUPCforwardMC.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskUPCforwardMC;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskUPCforwardMC) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::AliAnalysisTaskUPCforwardMC()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionIncoherentH(0),
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      fVectorGoodRunNumbers(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH(0),
      fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      fVectorCosThetaGenerated(0),
      fVectorCosThetaReconstructed(0),
      fCosThetaGeneratedHelicityFrame(0),
      fCosThetaReconHelicityFrame(0),
      fCounterUPCevent(0),
      fBinMigrationHelicityH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fMCCosThetaHelicityFrameJPsiH(0),
      fMCCosThetaCollinsSoperFrameJPsiH(0),
      fMCPhiHelicityFrameJPsiH(0),
      fMCPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::AliAnalysisTaskUPCforwardMC( const char* name )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionIncoherentH(0),
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      fVectorGoodRunNumbers(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH(0),
      fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      fVectorCosThetaGenerated(0),
      fVectorCosThetaReconstructed(0),
      fCosThetaGeneratedHelicityFrame(0),
      fCosThetaReconHelicityFrame(0),
      fCounterUPCevent(0),
      fBinMigrationHelicityH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fMCCosThetaHelicityFrameJPsiH(0),
      fMCCosThetaCollinsSoperFrameJPsiH(0),
      fMCPhiHelicityFrameJPsiH(0),
      fMCPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0)
{
    FillGoodRunVector(fVectorGoodRunNumbers);

    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::~AliAnalysisTaskUPCforwardMC()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers)
{
  fVectorGoodRunNumbers.clear();
  Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, 296376, 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623 };
  Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, 296977,
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595, 297623, 297624 };
  Int_t sizeOfLHC18q = 0;
  Int_t sizeOfLHC18r = 0;
  for ( Int_t GoodRunNumberLHC18q : listOfGoodRunNumbersLHC18q ) {
        fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18q);
        sizeOfLHC18q++;
  }
  for ( Int_t GoodRunNumberLHC18r : listOfGoodRunNumbersLHC18r ) {
        fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18r);
        sizeOfLHC18r++;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //

  //muon track cuts
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(    AliMuonTrackCuts::kMuEta     |
                                    AliMuonTrackCuts::kMuThetaAbs|
                                    AliMuonTrackCuts::kMuPdca    |
                                    AliMuonTrackCuts::kMuMatchLpt   );
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");


  fOutputList = new TList();          // this is a list which will contain all
                                      // of your histograms at the end of the
                                      // analysis, the contents of this list
                                      // are written to the output file

  fOutputList->SetOwner(kTRUE);       // memory management: the list is owner
                                      // of all objects it contains and will
                                      // delete them if requested

  //_______________________________
  // - Adding histograms
  fNumberMuonsH = new TH1F("fNumberMuonsH", "fNumberMuonsH", 12, -0.5, 11.5);
  fOutputList->Add(fNumberMuonsH);    // don't forget to add it to the list!

  fCounterH = new TH1F("fCounterH", "fCounterH", 24, -0.5, 23.5);
  fOutputList->Add(fCounterH);

  fEtaMuonH = new TH1F("fEtaMuonH", "fEtaMuonH", 90, -2, -5);
  fOutputList->Add(fEtaMuonH);

  fRAbsMuonH = new TH1F("fRAbsMuonH", "fRAbsMuonH", 100, 0, 100);
  fOutputList->Add(fRAbsMuonH);

  fInvariantMassDistributionH = new TH1F("fInvariantMassDistributionH", "fInvariantMassDistributionH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionH);

  fEntriesAgainstRunNumberH = new TH1F("fEntriesAgainstRunNumberH", "fEntriesAgainstRunNumberH", 10000, 290000, 300000);
  fOutputList->Add(fEntriesAgainstRunNumberH);

  /* - Trying to reproduce the histogram for the RunNumbers as they always
     - show it, properly labelled. Inspiration has come from the website:
     - https://root.cern.ch/doc/master/hlabels1_8C.html
     - Let us see if it works properly.
     -
   */
  fEntriesAgainstRunNumberProperlyH = new TH1F("fEntriesAgainstRunNumberProperlyH", "fEntriesAgainstRunNumberProperlyH", 10000, 290000, 300000);
  fEntriesAgainstRunNumberProperlyH->SetStats(0);
  fEntriesAgainstRunNumberProperlyH->SetFillColor(38);
  fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
  fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

  fInvariantMassDistributionCoherentH = new TH1F("fInvariantMassDistributionCoherentH", "fInvariantMassDistributionCoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentH);

  fInvariantMassDistributionIncoherentH = new TH1F("fInvariantMassDistributionIncoherentH", "fInvariantMassDistributionIncoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentH);

  fDimuonPtDistributionH = new TH1F("fDimuonPtDistributionH", "fDimuonPtDistributionH", 2000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionH);

  fTemplatePtDistributionH = new TH1F("fTemplatePtDistributionH", "fTemplatePtDistributionH", 2000, 0, 20);
  fOutputList->Add(fTemplatePtDistributionH);

  fDcaAgainstInvariantMassH = new TH2F("fDcaAgainstInvariantMassH", "fDcaAgainstInvariantMassH", 4000, 0, 40, 2000, -100, 100);
  fOutputList->Add(fDcaAgainstInvariantMassH);

  /* - These histograms have an EXTENDED range (0,20)->(0,40)
     -
   */
  fInvariantMassDistributionExtendedH = new TH1F("fInvariantMassDistributionExtendedH", "fInvariantMassDistributionExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionExtendedH);

  fInvariantMassDistributionCoherentExtendedH = new TH1F("fInvariantMassDistributionCoherentExtendedH", "fInvariantMassDistributionCoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionCoherentExtendedH);

  fInvariantMassDistributionIncoherentExtendedH = new TH1F("fInvariantMassDistributionIncoherentExtendedH", "fInvariantMassDistributionIncoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionIncoherentExtendedH);

  /* - Here starts the list of histograms needed for the analysis of the J/Psi's
     - polarization.
     -
   */
  fAngularDistribOfPositiveMuonRestFrameJPsiH = new TH1F("fAngularDistribOfPositiveMuonRestFrameJPsiH", "fAngularDistribOfPositiveMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfPositiveMuonRestFrameJPsiH);

  fAngularDistribOfNegativeMuonRestFrameJPsiH = new TH1F("fAngularDistribOfNegativeMuonRestFrameJPsiH", "fAngularDistribOfNegativeMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfNegativeMuonRestFrameJPsiH);

  fBinMigrationHelicityH = new TH2F("fBinMigrationHelicityH", "fBinMigrationHelicityH", 1000, -1., 1., 1000, -1., 1.);
  fOutputList->Add(fBinMigrationHelicityH);


  //_______________________________
  // - MC-only plots
  fMCpdgCodesH = new TH1F("fMCpdgCodesH", "fMCpdgCodesH", 3, 0, 3);
  fMCpdgCodesH ->SetStats(0);
  fMCpdgCodesH ->SetFillColor(38);
  fMCpdgCodesH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesH);

  fMCpdgCodesOnlyPrimaryH = new TH1F("fMCpdgCodesOnlyPrimaryH", "fMCpdgCodesOnlyPrimaryH", 3, 0, 3);
  fMCpdgCodesOnlyPrimaryH ->SetStats(0);
  fMCpdgCodesOnlyPrimaryH ->SetFillColor(38);
  fMCpdgCodesOnlyPrimaryH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesOnlyPrimaryH);

  fMCinvariantMassDistrJPsiGeneratedTruthH = new TH1F("fMCinvariantMassDistrJPsiGeneratedTruthH", "fMCinvariantMassDistrJPsiGeneratedTruthH", 2000000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiGeneratedTruthH);

  fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH = new TH1F("fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", "fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", 20000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH);

  fMCphiGeneratedTruthH = new TH1F("fMCphiGeneratedTruthH", "fMCphiGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiGeneratedTruthH);

  fMCetaGeneratedTruthH = new TH1F("fMCetaGeneratedTruthH", "fMCetaGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaGeneratedTruthH);

  fMCpseudorapidityGeneratedTruthH = new TH1F("fMCpseudorapidityGeneratedTruthH", "fMCpseudorapidityGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityGeneratedTruthH);

  fMCptGeneratedTruthH = new TH1F("fMCptGeneratedTruthH", "fMCptGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCptGeneratedTruthH);

  fMCphiDimuonGeneratedTruthH = new TH1F("fMCphiDimuonGeneratedTruthH", "fMCphiDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiDimuonGeneratedTruthH);

  fMCetaDimuonGeneratedTruthH = new TH1F("fMCetaDimuonGeneratedTruthH", "fMCetaDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaDimuonGeneratedTruthH);

  fMCpseudorapidityDimuonGeneratedTruthH = new TH1F("fMCpseudorapidityDimuonGeneratedTruthH", "fMCpseudorapidityDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityDimuonGeneratedTruthH);

  fMCptDimuonGeneratedTruthSingleMuonsH = new TH1F("fMCptDimuonGeneratedTruthSingleMuonsH", "fMCptDimuonGeneratedTruthSingleMuonsH", 2000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthSingleMuonsH);

  fMCptDimuonGeneratedTruthH = new TH1F("fMCptDimuonGeneratedTruthH", "fMCptDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthH);


  /* - Here starts the list of histograms needed for the analysis of the J/Psi's
     - polarization. GENERATED MC TRUTH.
     -
   */
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH = new TH1F("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH", "fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH", 1000, -1., 1.);
  fOutputList->Add(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH);

  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH = new TH1F("fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH", "fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH", 1000, -1., 1.);
  fOutputList->Add(fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH);

  fCheckHelicityRestFrameJPsiH = new TH1F("fCheckHelicityRestFrameJPsiH", "fCheckHelicityRestFrameJPsiH", 100000, -50., 50.);
  fOutputList->Add(fCheckHelicityRestFrameJPsiH);

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin] = new TH1F(
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin]);

    fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin] = new TH1F(
                Form("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH_%d", iRapidityBin),
                Form("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]);

  }

  fCosThetaHelicityFrameJPsiH = new TH1F("fCosThetaHelicityFrameJPsiH", "fCosThetaHelicityFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaHelicityFrameJPsiH);

  fCosThetaCollinsSoperFrameJPsiH = new TH1F("fCosThetaCollinsSoperFrameJPsiH", "fCosThetaCollinsSoperFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaCollinsSoperFrameJPsiH);

  fPhiHelicityFrameJPsiH = new TH1F("fPhiHelicityFrameJPsiH", "fPhiHelicityFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiHelicityFrameJPsiH);

  fPhiCollinsSoperFrameJPsiH = new TH1F("fPhiCollinsSoperFrameJPsiH", "fPhiCollinsSoperFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiCollinsSoperFrameJPsiH);

  fMCCosThetaHelicityFrameJPsiH = new TH1F("fMCCosThetaHelicityFrameJPsiH", "fMCCosThetaHelicityFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fMCCosThetaHelicityFrameJPsiH);

  fMCCosThetaCollinsSoperFrameJPsiH = new TH1F("fMCCosThetaCollinsSoperFrameJPsiH", "fMCCosThetaCollinsSoperFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiH);

  fMCPhiHelicityFrameJPsiH = new TH1F("fMCPhiHelicityFrameJPsiH", "fMCPhiHelicityFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fMCPhiHelicityFrameJPsiH);

  fMCPhiCollinsSoperFrameJPsiH = new TH1F("fMCPhiCollinsSoperFrameJPsiH", "fMCPhiCollinsSoperFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fMCPhiCollinsSoperFrameJPsiH);

  /* - Rapidity-dependent analysis for the helicity of the J/Psi.
     - It is divided into RECONSTRUCTED and GENERATED for both
     - helicity frame and Collins-Soper.
     -
   */
  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  /* - Rapidity-dependent analysis for the helicity of the J/Psi.
     - It is divided into RECONSTRUCTED and GENERATED for both
     - helicity frame and Collins-Soper.
     -
     - The following are the same as above divided in 10 rapidity bins.
     - The reasons for this addition is explained in the file:
     - PWGUD/UPC/AliAnalysisTaskUPCforward.cxx
   */
  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH =
        new TH2F( "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  80, -1, 1,
                  80, -4, 4
                  );
  fOutputList->Add(fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH);

  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH =
        new TH2F( "fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  "fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  80, -1, 1,
                  80, -4, 4
                  );
  fOutputList->Add(fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH);

  //_______________________________
  // - End of the function
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::UserExec(Option_t *)
{
  /* - This iSelectionCounter is used as a token. So at every passing step it is
     - increased by one. All the events are supposed to pass the first step
     - obviously, but only a few get to the end. This is effect is clearly
     - noticeable in fCounterH event with the small trial local sample.
     - Almost 160k possible events at the 0-th step, while only 2k at the 4th step.
   */
  Int_t iSelectionCounter = 0; // no selection applied yet
  fCounterH->Fill(iSelectionCounter); // entering UserExec
  iSelectionCounter++;

  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
      PostData(1, fOutputList);
      return;
  }
  fMCEvent = MCEvent();
  if(!fMCEvent) {
      PostData(1, fOutputList);
      return;
  }
  if(fMCEvent) {
    ProcessMCParticles(fMCEvent);
    fCounterUPCevent += 1;
  }
  /* - We are now checking if there were any tracks. If there were at least one,
     - then the histogram gets filled again. If not we are returning. There
     - would be no point in going further.
   */
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter); // At least one track
  iSelectionCounter++;


  //_______________________________
  // EVENT DATA EXTRACTION
  /* - Eugeny Krishen's event data extraction. I am trying to implement it.
     - The only thing I am a bit worried about is whether it should go before or
     - after the "nTracks<1" check... I will try and switch it if it sounds
     - better. These data are used for the event selection and maybe later on
     - for track selection, but I did not get to that part yet. If after all of
     - this I remember to do so, I will come back to this point and correct this
     - statement. If you find this part, please, keep in mind to check the
     - following.
   */

  /* - Event information:
     - Run Number, maybe to select the GOOD Runs and discard the others;
     - Number of Tracklets, these are in this case the SPD tracklets, so the
     - almost unit vector roughly 2 cm between two pixels of the SPD in different
     - layers.
   */
  fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  // fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  // fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  // fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  // fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

  /* - ZDC: we try to find the ZDC object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without ZDC's information. Or at least, this
     - is my impression of it (filling fCounterH). ZDC information:
     - fZem1Energy:
     - fZem2Energy:
     - fZNAEnergy:
     - fZNCEnergy:
     - fZPAEnergy:
     - fZPCEnergy:
     - fZNATime:
     - fZNCTime:
     - fZNATDC[i]:
     - fZNCTDC[i]:
     - fZPATDC[i]:
     - fZPCTDC[i]:
   */
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();
  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  /* - These lines are the calibration for the ZDC as provided by Evgeny Kryshen.
     -
   */
  // Bool_t calibrated = 0;
  // if ( fRunNum <  295726 ) calibrated = 1;
  // if ( fRunNum == 296509 ) calibrated = 1;
  // if ( fRunNum >  296689 ) calibrated = 1;
  // if ( fRunNum >  296695 ) calibrated = 0;
  // if ( fRunNum == 297219 ) calibrated = 1;
  // if ( fRunNum == 297221 ) calibrated = 1;
  // if ( fRunNum == 297415 ) calibrated = 1;
  //
  // if ( !calibrated ) {
  //   fZNAEnergy *= (2500./190.);
  //   fZNCEnergy *= (2500./190.);
  // }

  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  // Reset event info
  fBBCFlags = 0;
  fBGCFlags = 0;
  fBBAFlags = 0;
  fBGAFlags = 0;
  for (Int_t i=0; i<64; i++) {
    // get array of fired cells
    fBBFlag[i] = dataVZERO->GetBBFlag(i);
    fBGFlag[i] = dataVZERO->GetBGFlag(i);
  }

  for (Int_t i=0; i<32; i++){ // loop over cells
    fBBCFlags+=fBBFlag[i];
    fBGCFlags+=fBGFlag[i];
    fBBAFlags+=fBBFlag[i+32];
    fBGAFlags+=fBGFlag[i+32];
  }


  //_____________________________________
  // RUN SELECTION
  /* - This part is the run selection. We call the std::find() method of the
     - <algorithm> library STL. This returns a kTRUE if you find the value you
     - are looking for inside the vector. So, what happens is that I look for
     - the Run Numbers inside the vector containing them. If I cannot find
     - them I move on to the next event.
     -
   */
   auto findRunNumber = std::find(  std::begin(fVectorGoodRunNumbers),
                                    std::end(fVectorGoodRunNumbers),
                                    fRunNum
                                    );
   if (findRunNumber != std::end(fVectorGoodRunNumbers)) {
        // std::cout << "fVectorGoodRunNumbers DOES     contain: " << fRunNum << std::endl;
   } else {
        PostData(1, fOutputList);
        return;
   }
  // END RUN SELECTION
  //_____________________________________



  /* - We have to get the number of fired V0C cells. So firstly, we get the
     - boolean information about the hit cells for all V0. This is done through
     - the GetBBFlag(i) method, where 0<i<32 stands for the V0C cells and
     - 32<i<64 for the V0A cells. Then I thought the easiest way to check
     - whether the number of fired V0C cells is above 2 is just to add up the
     - boolean numbers for 0<i<32. Let's see.
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
  }

  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();

        // Reset event info
        fBBCFlagsAD = 0;
        fBGCFlagsAD = 0;
        fBBAFlagsAD = 0;
        fBGAFlagsAD = 0;
        for(Int_t i=0; i<16; i++) {
          // get array of fired pads
          fBBFlagAD[i] = dataAD->GetBBFlag(i);
          fBGFlagAD[i] = dataAD->GetBGFlag(i);
        }

        for(Int_t i=0; i<4; i++) { // loop over pairs of pads
          if ( fBBFlagAD[i]   && fBBFlagAD[i+4]  ) fBBCFlagsAD++;
          if ( fBGFlagAD[i]   && fBGFlagAD[i+4]  ) fBGCFlagsAD++;
          if ( fBBFlagAD[i+8] && fBBFlagAD[i+12] ) fBBAFlagsAD++;
          if ( fBGFlagAD[i+8] && fBGFlagAD[i+12] ) fBGAFlagsAD++;
        }
  }
  // END EVENT DATA EXTRACTION
  //_______________________________
  // APPLY TRIGGER MC!
  if(!IsTriggered()) {
    cout << "Ehm" ;
    PostData(1, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;
  //_______________________________
  // EVENT SELECTION
  /* - This is Eugeny Krishen's event selection from the talk in 14/1/2019 for
     - the PWG-UD (UPC oriented) meeting. The event selection requires:
     - CMUP11-B triggers;
     - Maximum 2 V0C cells fired;
     - Empty V0A decision;
     - Empty ADA decision;
     - Empty ADC decision;
     - 0 tracklets in SPD;
     - Exactly 2 unlike-sign muons;
   */
  /* - Empty V0A decision
     - Empty ADA decision
     - Empty ADC decision
   */
  if(fV0ADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADCDecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - 0 tracklets in SPD
   */
  if(fTracklets != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - Maximum 2 V0C cells fired.
   */
  if( fV0TotalNCells > 2 ) {
       PostData(1, fOutputList);
       return;
  }

  /* - We are finally at the starting point. We loop over the tracks and select
     - the good muons. Later on everything should happen in this loop. Let us
     - see what the future has in hold.
     -
     - Saturday: I moved the creation of the AliAODTrack* track outside of the
     - loop as it would have been otherwise created for each single iteration.
     - This could have caused massive memory issues especially to grid. I have
     - added a second AliAODTrack* track[2] to hold the second supposed muon.
     - Now this is ready to send the information to two TLorentzVectors to
     - obtain the invariant mass of the J/Psi through the Mag() method of the
     - class. Hope for the best.
   */
  Int_t nGoodMuons = 0;
  AliAODTrack* track[2];
  track[0]         = 0x0;
  track[1]         = 0x0;
  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    /* - This should be another form of event selection.
       - I am basically requesting the presence of TWO good muons only.
       - Later I will be checking whether of they are likesign or unlikesign.
     */
    if(nGoodMuons > 2) {
         PostData(1, fOutputList);
         return;
    }
    track[nGoodMuons] = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track[nGoodMuons]) return;

    // is it a good muon track?
    if(!track[nGoodMuons]->IsMuonTrack()) {
        // track[nGoodMuons] = 0x0;
        continue;
    }
    if(!fMuonTrackCuts->IsSelected(track[nGoodMuons])) {
        // track[nGoodMuons] = 0x0;
        continue;
    }

    // MUON SELECTION
    /* - This is Eugeny Krishen's MUON selection from the talk in 14/1/2019 for
       - the PWG-UD (UPC oriented) meeting. The event selection requires:
       - Muon trigger matching >=2 (1 GeV/c threshold);
       - (-4) < eta < (-2.5);
       - (17.5 cm) < R_{abs} < (89.5 cm);
       - p \times DCA cut;
    */

    // increase counter
    nGoodMuons++;
  }
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuons != 2 ) {
        PostData(1, fOutputList);
        return;
  }
  /* - Implementing the track cut on the unlike muons
   * -
   */
  if( (track[0]->Charge()) == (track[1]->Charge()) ) {
        PostData(1, fOutputList);
        return;
  }
  for(Int_t iFilling = 0; iFilling < nGoodMuons; iFilling++) {
        fEtaMuonH ->Fill(track[iFilling]->Eta());
        fRAbsMuonH->Fill(track[iFilling]->GetRAtAbsorberEnd());
  }
  // store muons
  fNumberMuonsH->Fill(nGoodMuons);
  fEntriesAgainstRunNumberH->Fill(fRunNum);
  /* - This is the last part of my try to obtain a proper RunNumbers histogram...
     -
   */
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  if (nGoodMuons>0) fCounterH->Fill(iSelectionCounter); // At least one good muon
  iSelectionCounter++;

  /* - Finally the core!!!
   * - What will be happening is that we will instantiate TLorentzVectors to
   * - obtain the invariant mass of the dimuon system. If everything goes fine
   * - after this we should be able to obtain the peak of the J/Psi. But
   * - things never go as expected, so who knows!
   */
  TLorentzVector muons[2];
  TLorentzVector possibleJPsi;
  Double_t       chargeOfMuons[2];
  for(int indexMuon = 0; indexMuon < 2; indexMuon++) {
        muons[indexMuon].SetPtEtaPhiM(   track[indexMuon]->Pt(),
                                         track[indexMuon]->Eta(),
                                         track[indexMuon]->Phi(),
                                         TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                       );
        possibleJPsi += muons[indexMuon];
        chargeOfMuons[indexMuon] = track[indexMuon]->Charge();
  }
  fInvariantMassDistributionH->Fill(possibleJPsi.Mag());
  fInvariantMassDistributionExtendedH->Fill(possibleJPsi.Mag());

  /* - This is a TH2F histogram filled with DCA against the invariant mass of
     - the dimuon pair. This should be plotted as LEGO or as a heat map...
     - Is this the correct way of doing it?? I have to ask to my supervisor!!
     -
   */
  for(int indexMuonForDCA = 0; indexMuonForDCA < 2; indexMuonForDCA++) {
        fDcaAgainstInvariantMassH->Fill(possibleJPsi.Mag(), track[indexMuonForDCA]->DCA());
  }

  /* - Now we are evaluating the pt of the dimuon pair. Generally speaking,
     - if such a pt is less than 0.25 GeV/c then it fills the coherent
     - component, otherwise the incoherent component. At this point we may fill
     - even the dimuon pt distribution histogram and see if it looks like Kay's.
     -
   */
  Double_t ptOfTheDimuonPair = possibleJPsi.Pt();
  if( ptOfTheDimuonPair < 0.25) {
        fInvariantMassDistributionCoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionCoherentExtendedH->Fill(possibleJPsi.Mag());
  } else {
        fInvariantMassDistributionIncoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionIncoherentExtendedH->Fill(possibleJPsi.Mag());
  }
  fDimuonPtDistributionH->Fill(ptOfTheDimuonPair);
  if ( (possibleJPsi.Mag() > 2.8) && (possibleJPsi.Mag() < 3.3) ) fTemplatePtDistributionH->Fill(ptOfTheDimuonPair);


  /* - Filling the J/Psi's polarization plots.
     -
     - Now we are ordering the muons. The first muon will always be positive.
     - This is useful for the histograms...
   */
  TLorentzVector muonsCopy[2];
  TLorentzVector muonsCopy2[2];
  TLorentzVector possibleJPsiCopy;
  if( chargeOfMuons[0] > 0 ){
    muonsCopy[0]     = muons[0];
    muonsCopy[1]     = muons[1];
  } else if( chargeOfMuons[0] < 0 ){
    muonsCopy[0]     = muons[1];
    muonsCopy[1]     = muons[0];
  }
  muonsCopy2[0]      = muonsCopy[0];
  muonsCopy2[1]      = muonsCopy[1];
  possibleJPsiCopy = possibleJPsi;
  for( Int_t iBoosting = 0; iBoosting < 2; iBoosting++ ) {
    // TLorentzVector boostBack = -(possibleJPsiCopy).BoostVector();
    /* - This snippet has beem taken from the website:
       - http://personalpages.to.infn.it/~puccio/htmldoc/src/AliAODDimuon.cxx.html
       -
     */
    TVector3 beta=(-1./possibleJPsiCopy.E())*possibleJPsiCopy.Vect();
    muonsCopy[iBoosting].Boost( beta );
  }
  Double_t cosThetaMuonsRestFrame[2];
  for( Int_t iAngle = 0; iAngle < 2; iAngle++ ) {
    TVector3 muonsCopyVector        = muonsCopy[iAngle].Vect();
    TVector3 possibleJPsiCopyVector = possibleJPsiCopy.Vect();
    Double_t dotProductMuonJPsi     = muonsCopyVector.Dot(possibleJPsiCopyVector);
    cosThetaMuonsRestFrame[iAngle]  = dotProductMuonJPsi/( muonsCopyVector.Mag() * possibleJPsiCopyVector.Mag() );
  }
  /* - If we are in the J/Psi peak, hence 2.8 < M < 3.3 GeV/c, AND if we are
     - in the coherent regime, so if the Pt < 0.25 GeV/c, we fill the plots.
     -
     - In the following note that the rapidity is well computed, so we are
     - dealing with negative values... -4.0 < Y < -2.5 !!!
     -
   */
  if ( (possibleJPsiCopy.Mag() > 2.8) && (possibleJPsiCopy.Mag() < 3.3) && (possibleJPsiCopy.Pt() < 0.25) ) {
    fAngularDistribOfPositiveMuonRestFrameJPsiH->Fill(cosThetaMuonsRestFrame[0]);
    fAngularDistribOfNegativeMuonRestFrameJPsiH->Fill(cosThetaMuonsRestFrame[1]);
    fCheckHelicityRestFrameJPsiH->Fill( muonsCopy[0].Dot(muonsCopy[1]) );

    /* - New part: filling all possible histograms!
       -
     */
    fCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                              muonsCopy2[1],
                                                              possibleJPsiCopy
                                                              )
                                                            );
    fCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                 muonsCopy2[1],
                                                                 possibleJPsiCopy
                                                                 )
                                                               );
    fPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                       muonsCopy2[1],
                                                       possibleJPsiCopy
                                                       )
                                                     );
    fPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                          muonsCopy2[1],
                                                          possibleJPsiCopy
                                                          )
                                                        );
    /* - Now we are filling in terms of rapidity...
       - The easiest way to do so I have envisioned is to simply
       - check everytime if we are below the following threshold
       - in a consecutive sense. This means that if we have not passed
       - the previous check we are at least above it.
       - This readily defines the rapidity bin.
       -
       */
    for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
        if( (possibleJPsiCopy.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
          fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrame[0]);
          /* - New part: filling all possible histograms!
             -
           */
          fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                              muonsCopy2[1],
                                                                                              possibleJPsiCopy
                                                                                              )
                                                                                           );
          fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                                                 muonsCopy2[1],
                                                                                                 possibleJPsiCopy
                                                                                                 )
                                                                                               );
          fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                                                       muonsCopy2[1],
                                                                                       possibleJPsiCopy
                                                                                       )
                                                                                     );
          fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                                                          muonsCopy2[1],
                                                                                          possibleJPsiCopy
                                                                                          )
                                                                                        );
          break;
        }
    }


    /* - Now we are filling in terms of rapidity...
       - The easiest way to do so I have envisioned is to simply
       - check everytime if we are below the following threshold
       - in a consecutive sense. This means that if we have not passed
       - the previous check we are at least above it.
       - This readily defines the rapidity bin.
       -
       -
       - NEW: the following code is for 10 rapidity bins...
       */
    for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
        if( (possibleJPsiCopy.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
          fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                                 muonsCopy2[1],
                                                                                                 possibleJPsiCopy
                                                                                                 )
                                                                                              );
          fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                                                    muonsCopy2[1],
                                                                                                    possibleJPsiCopy
                                                                                                    )
                                                                                               );
          fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                                                          muonsCopy2[1],
                                                                                          possibleJPsiCopy
                                                                                          )
                                                                                     );
          fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                                                             muonsCopy2[1],
                                                                                             possibleJPsiCopy
                                                                                             )
                                                                                        );
          break;
        }
    }


    /* - What we do here is very similar.
       - This time we divide firstly in bins of CosTheta.
       - As many as needed.
       - And then we divide again in terms of Phi.
       - Then we fill.
       - This way we should be able to obtain some kind of map...
       -
     */
    fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                               muonsCopy2[1],
                                                                                               possibleJPsiCopy
                                                                                               ),
                                                                        CosPhiHelicityFrame( muonsCopy2[0],
                                                                                             muonsCopy2[1],
                                                                                             possibleJPsiCopy
                                                                                             )
                                                                        );

  }

  // fVectorCosThetaReconstructed.push_back(cosThetaMuonsRestFrame[0]);
  fCosThetaReconHelicityFrame = cosThetaMuonsRestFrame[0];
  /* - Mind that it could generate segmentation fault without
     - fCounterUPCevent-1, because we are incrementing the counter right after
     - it processes the MC events at Generated level...
     -
     - Comparing old version with vector (unusable on MC due to the huge
     - vector size) to the new easier on resources...
     -
   */
  // if ( fVectorCosThetaGenerated.at(fCounterUPCevent-1) && cosThetaMuonsRestFrame[0] ) {
  //       fBinMigrationHelicityH->Fill( fVectorCosThetaGenerated.at(fCounterUPCevent-1),
  //                                     cosThetaMuonsRestFrame[0]
  //                                   );
  // }
  if ( fCosThetaReconHelicityFrame ) {
        fBinMigrationHelicityH->Fill( fCosThetaGeneratedHelicityFrame,
                                      fCosThetaReconHelicityFrame
                                    );
  }




  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUPCforwardMC::IsTriggered()
{
  /* - This function IS roughly speaking the trigger for the MC.
     - This code has been greatly inspired by David Horak's work.
     - It returns kTRUE if the trigger has been fired.
     -
   */
  Bool_t is0VBAfired = fBBAFlags > 0;
  Bool_t is0VBCfired = fBBCFlags > 0;
  Bool_t is0UBAfired = fBBAFlagsAD > 0;
  Bool_t is0UBCfired = fBBCFlagsAD > 0;
  // cout << "is0VBAfired = ( fBBAFlags   = " << fBBAFlags   << " ) > 0 => " << is0VBAfired << endl;
  // cout << "is0VBCfired = ( fBBCFlags   = " << fBBCFlags   << " ) > 0 => " << is0VBCfired << endl;
  // cout << "is0UBAfired = ( fBBAFlagsAD = " << fBBAFlagsAD << " ) > 0 => " << is0UBAfired << endl;
  // cout << "is0UBCfired = ( fBBCFlagsAD = " << fBBCFlagsAD << " ) > 0 => " << is0UBCfired << endl;
  if (!is0VBAfired && !is0UBAfired && !is0UBCfired ) return kTRUE;
  else return kFALSE;

}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::ProcessMCParticles(AliMCEvent* fMCEventArg)
{
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  mcParticle->Charge()    == 0) continue;
    Double_t pT     = mcParticle->Pt();
    Double_t eta    = mcParticle->Eta();
    Double_t pseudo = mcParticle->Y();
    Double_t phi    = mcParticle->Phi();
    Short_t  charge = mcParticle->Charge();
    Int_t    pdg    = mcParticle->PdgCode();
    fMCpdgCodesH->Fill( Form("%d", pdg) , 1 );
    if ( mcParticle->IsPrimary() ) fMCpdgCodesOnlyPrimaryH->Fill( Form("%d", pdg) , 1 );
    fMCphiGeneratedTruthH           ->Fill(phi);
    fMCetaGeneratedTruthH           ->Fill(eta);
    fMCpseudorapidityGeneratedTruthH->Fill(pseudo);
    fMCptGeneratedTruthH            ->Fill(pT);
  }

  /* - This is loop where we record the MC TRUTH stuff BEFORE both
     - the event selection and the track selection.
     - It is nonetheless strange that we can't directly use the
     - J/Psi and that we have to o through such a roundabout
     - road like reconstructing the invariant mass...
     -
   */
  Int_t nGoodMuonsMC = 0;
  TLorentzVector muonsMC[2];
  TLorentzVector possibleJPsiMC;
  Double_t pT[2];
  Double_t eta[2];
  Double_t pseudo[2];
  Double_t phi[2];
  Short_t  charge[2];
  Int_t    pdg[2];
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    if(nGoodMuonsMC > 2) continue;
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  !mcParticle->IsPrimary()                 ) continue;
    if (   mcParticle->Charge()              == 0  ) continue;
    if (   TMath::Abs(mcParticle->PdgCode()) == 13 ) {
      if ( nGoodMuonsMC < 2 ) {
        // cout << "Ok" << nGoodMuonsMC << endl;
        pT[nGoodMuonsMC]     = mcParticle->Pt();
        eta[nGoodMuonsMC]    = mcParticle->Eta();
        pseudo[nGoodMuonsMC] = mcParticle->Y();
        phi[nGoodMuonsMC]    = mcParticle->Phi();
        charge[nGoodMuonsMC] = mcParticle->Charge();
        pdg[nGoodMuonsMC]    = mcParticle->PdgCode();
        muonsMC[nGoodMuonsMC].SetPtEtaPhiM( pT[nGoodMuonsMC],
                                            eta[nGoodMuonsMC],
                                            phi[nGoodMuonsMC],
                                            TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                           );
        nGoodMuonsMC++;
      }
    }
  }
  // cout << "Tracks in this event:  " << fMCEventArg->GetNumberOfTracks() << endl;
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuonsMC == 2 ) {
    /* - Unlike-sign muons
       -
     */
    if( charge[0] != charge[1] ) {
      TLorentzVector muonsMCcopy[2];
      muonsMCcopy[0] = muonsMC[0];
      muonsMCcopy[1] = muonsMC[1];
      for( Int_t iMuonsMC = 0; iMuonsMC < 2; iMuonsMC++ ) {
        possibleJPsiMC+=muonsMC[iMuonsMC];
        fMCphiDimuonGeneratedTruthH           ->Fill(phi[iMuonsMC]);
        fMCetaDimuonGeneratedTruthH           ->Fill(eta[iMuonsMC]);
        fMCpseudorapidityDimuonGeneratedTruthH->Fill(pseudo[iMuonsMC]);
        fMCptDimuonGeneratedTruthSingleMuonsH ->Fill(pT[iMuonsMC]);
      }
      fMCinvariantMassDistrJPsiGeneratedTruthH->Fill(possibleJPsiMC.Mag());
      fMCptDimuonGeneratedTruthH->Fill(possibleJPsiMC.Pt());
      for( Int_t iBoosting = 0; iBoosting < 2; iBoosting++ ) {
        // TLorentzVector boostBack = -(possibleJPsiCopy).BoostVector();
        /* - This snippet has beem taken from the website:
           - http://personalpages.to.infn.it/~puccio/htmldoc/src/AliAODDimuon.cxx.html
           -
         */
        TVector3 betaMC = (-1./possibleJPsiMC.E())*possibleJPsiMC.Vect();
        muonsMC[iBoosting].Boost( betaMC );
      }
      Double_t cosThetaMuonsRestFrameMC[2];
      for( Int_t iAngle = 0; iAngle < 2; iAngle++ ) {
        TVector3 muonsMCVector            = muonsMC[iAngle].Vect();
        TVector3 possibleJPsiMCVector     = possibleJPsiMC.Vect();
        Double_t dotProductMuonJPsiMC     = muonsMCVector.Dot(possibleJPsiMCVector);
        cosThetaMuonsRestFrameMC[iAngle]  = dotProductMuonJPsiMC/( muonsMCVector.Mag() * possibleJPsiMCVector.Mag() );
      }
      if ( (possibleJPsiMC.Mag() > 2.8) && (possibleJPsiMC.Mag() < 3.3) && (possibleJPsiMC.Pt() < 0.25) ) {
          if( charge[0] > 0 ) {
                  /* - This means that [0] is the positive muon while [1]
                     - is the negative muon!
                     -
                   */
                  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[0]);
                  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[1]);
                  // fVectorCosThetaGenerated.push_back(cosThetaMuonsRestFrameMC[0]);
                  fCosThetaGeneratedHelicityFrame = cosThetaMuonsRestFrameMC[0];
                  /* - New part: filling all possible histograms!
                     -
                   */
                  fMCCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                              muonsMCcopy[1],
                                                                              possibleJPsiMC
                                                                              )
                                                                            );
                  fMCCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                 muonsMCcopy[1],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                               );
                  fMCPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                       muonsMCcopy[1],
                                                                       possibleJPsiMC
                                                                       )
                                                                     );
                  fMCPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                          muonsMCcopy[1],
                                                                          possibleJPsiMC
                                                                          )
                                                                        );
                  /* - Now we are filling in terms of rapidity...
                     - The easiest way to do so I have envisioned is to simply
                     - check everytime if we are below the following threshold
                     - in a consecutive sense. This means that if we have not passed
                     - the previous check we are at least above it.
                     - This readily defines the rapidity bin.
                     -
                     */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
                        fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrameMC[0]);
                        /* - New part: filling all possible histograms!
                           -
                         */
                        fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                              muonsMCcopy[1],
                                                                                                              possibleJPsiMC
                                                                                                              )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                                                 muonsMCcopy[1],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                       muonsMCcopy[1],
                                                                                                       possibleJPsiMC
                                                                                                       )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                                                          muonsMCcopy[1],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                        );
                        break;
                      }
                  }

                  /* - Filling for 10 rapidity bins.
                     -
                   */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
                        fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                                 muonsMCcopy[1],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                                                    muonsMCcopy[1],
                                                                                                                    possibleJPsiMC
                                                                                                                    )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                          muonsMCcopy[1],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                                                             muonsMCcopy[1],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                                        );
                        break;
                      }
                  }
                  /* - What we do here is very similar.
                     - This time we divide firstly in bins of CosTheta.
                     - As many as needed.
                     - And then we divide again in terms of Phi.
                     - Then we fill.
                     - This way we should be able to obtain some kind of map...
                     -
                   */
                  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                               muonsMCcopy[1],
                                                                                                               possibleJPsiMC
                                                                                                               ),
                                                                                        CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                             muonsMCcopy[1],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                        );
          } else  {
                  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[0]);
                  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[1]);
                  // fVectorCosThetaGenerated.push_back(cosThetaMuonsRestFrameMC[1]);
                  fCosThetaGeneratedHelicityFrame = cosThetaMuonsRestFrameMC[1];
                  /* - New part: filling all possible histograms!
                     -
                   */
                  fMCCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                              muonsMCcopy[0],
                                                                              possibleJPsiMC
                                                                              )
                                                                            );
                  fMCCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                 muonsMCcopy[0],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                               );
                  fMCPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                       muonsMCcopy[0],
                                                                       possibleJPsiMC
                                                                       )
                                                                     );
                  fMCPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                          muonsMCcopy[0],
                                                                          possibleJPsiMC
                                                                          )
                                                                        );
                  /* - Now we are filling in terms of rapidity...
                     - The easiest way to do so I have envisioned is to simply
                     - check everytime if we are below the following threshold
                     - in a consecutive sense. This means that if we have not passed
                     - the previous check we are at least above it.
                     - This readily defines the rapidity bin.
                     -
                     */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
                        fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrameMC[1]);
                        /* - New part: filling all possible histograms!
                           -
                         */
                        fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                              muonsMCcopy[0],
                                                                                                              possibleJPsiMC
                                                                                                              )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                                                 muonsMCcopy[0],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                       muonsMCcopy[0],
                                                                                                       possibleJPsiMC
                                                                                                       )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                                                          muonsMCcopy[0],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                        );
                        break;
                      }
                  }

                  /* - Filling for 10 rapidity bins.
                     -
                   */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
                        fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                                 muonsMCcopy[0],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                                                    muonsMCcopy[0],
                                                                                                                    possibleJPsiMC
                                                                                                                    )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                          muonsMCcopy[0],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                                                             muonsMCcopy[0],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                                        );
                        break;
                      }
                  }
                  /* - What we do here is very similar.
                     - This time we divide firstly in bins of CosTheta.
                     - As many as needed.
                     - And then we divide again in terms of Phi.
                     - Then we fill.
                     - This way we should be able to obtain some kind of map...
                     -
                   */
                  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                               muonsMCcopy[0],
                                                                                                               possibleJPsiMC
                                                                                                               ),
                                                                                        CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                             muonsMCcopy[0],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                        );

          }
      }
    }
  }

}
//_____________________________________________________________________________
/* - The following are code snippets adapted from the AliAODDimuon class.
   - The problem is that that class was adapted specifically for the
   - inclusive people's analysis, hence it is not fit for the UPC...
   -
 */
Double_t AliAnalysisTaskUPCforwardMC::CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaCS = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaCS;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                             TLorentzVector muonNegative,
                                                             TLorentzVector possibleJPsi )
{
  /* - This function computes the Helicity cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  /* - Determine the He angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaHE = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaHE;

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();

  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  return   phi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosPhiHelicityFrame(  TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the helicity phi for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
  */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  return   phi;
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
