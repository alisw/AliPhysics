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

// --- ROOT system ---
#include <TH3.h>
#include <TH2F.h>
//#include "Riostream.h"
#include <TCanvas.h>
#include <TPad.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TDatabasePDG.h>

//---- AliRoot system ----
#include "AliAnaPi0.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliMCEvent.h"
#include "AliFiducialCut.h"
#include "AliVEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliNeutralMesonSelection.h"
#include "AliMixedEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"

// --- Detectors ---
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaPi0) ;
/// \endcond

//______________________________________________________
/// Default Constructor. Initialized parameters with default values.
//______________________________________________________
AliAnaPi0::AliAnaPi0() : AliAnaCaloTrackCorrBaseClass(),
fEventsList(0x0),
fUseAngleCut(kFALSE),        fUseAngleEDepCut(kFALSE),     fAngleCut(0),                 fAngleMaxCut(0.),   fUseOneCellSeparation(kFALSE),
fMultiCutAna(kFALSE),        fMultiCutAnaSim(kFALSE),      fMultiCutAnaAcc(kFALSE),
fNPtCuts(0),                 fNAsymCuts(0),                fNCellNCuts(0),               fNPIDBits(0), fNAngleCutBins(0),
fMakeInvPtPlots(kFALSE),     fSameSM(kFALSE),
fFillSMCombinations(kFALSE), fCheckConversion(kFALSE),
fFillBadDistHisto(kFALSE),   fFillSSCombinations(kFALSE),
fFillAngleHisto(kFALSE),     fFillAsymmetryHisto(kFALSE),  
fFillOriginHisto(0),         fFillOriginHistoForMesonsOnly(1), 
fFillArmenterosThetaStar(0), fFillOnlyMCAcceptanceHisto(0),
fFillSecondaryCellTiming(0), fFillOpAngleCutHisto(0),      fCheckAccInSector(0),
fPairWithOtherDetector(0),   fOtherDetectorInputName(""),
fPhotonMom1(),               fPhotonMom1Boost(),           fPhotonMom2(),                fMCPrimMesonMom(),
fMCProdVertex(),

// Histograms
fhReMod(0x0),                fhReSameSideEMCALMod(0x0),    fhReSameSectorEMCALMod(0x0),  fhReDiffPHOSMod(0x0),
fhReSameSectorDCALPHOSMod(0),fhReDiffSectorDCALPHOSMod(0),
fhMiMod(0x0),                fhMiSameSideEMCALMod(0x0),    fhMiSameSectorEMCALMod(0x0),  fhMiDiffPHOSMod(0x0),
fhMiSameSectorDCALPHOSMod(0),fhMiDiffSectorDCALPHOSMod(0),
fhReConv(0x0),               fhMiConv(0x0),                fhReConv2(0x0),  fhMiConv2(0x0),
fhRe1(0x0),                  fhMi1(0x0),                   fhRe2(0x0),                   fhMi2(0x0),
fhRe3(0x0),                  fhMi3(0x0),                   fhReInvPt1(0x0),              fhMiInvPt1(0x0),
fhReInvPt2(0x0),             fhMiInvPt2(0x0),              fhReInvPt3(0x0),              fhMiInvPt3(0x0),
fhRePtNCellAsymCuts(0x0),    fhMiPtNCellAsymCuts(0x0),     fhRePtNCellAsymCutsSM(),
fhRePtNCellAsymCutsOpAngle(0x0),    fhMiPtNCellAsymCutsOpAngle(0x0),                     
fhRePtAsym(0x0),             fhRePtAsymPi0(0x0),           fhRePtAsymEta(0x0),
fhMiPtAsym(0x0),             fhMiPtAsymPi0(0x0),           fhMiPtAsymEta(0x0),
fhEventBin(0),               fhEventMixBin(0),
fhCentrality(0x0),           fhCentralityNoPair(0x0),
fhEventPlaneResolution(0x0),
fhRealOpeningAngle(0x0),     fhRealCosOpeningAngle(0x0),   fhMixedOpeningAngle(0x0),     fhMixedCosOpeningAngle(0x0),
// MC histograms
fhPrimPi0E(0x0),             fhPrimPi0Pt(0x0),             fhPrimPi0PtInCalo(0x0),
fhPrimPi0AccE(0x0),          fhPrimPi0AccPt(0x0),          fhPrimPi0AccPtPhotonCuts(0x0),
fhPrimPi0Y(0x0),             fhPrimPi0AccY(0x0),
fhPrimPi0Yeta(0x0),          fhPrimPi0YetaYcut(0x0),       fhPrimPi0AccYeta(0x0),
fhPrimPi0Phi(0x0),           fhPrimPi0AccPhi(0x0),
fhPrimPi0OpeningAngle(0x0),  fhPrimPi0OpeningAnglePhotonCuts(0x0),
fhPrimPi0OpeningAngleAsym(0x0),fhPrimPi0CosOpeningAngle(0x0),
fhPrimPi0PtCentrality(0),    fhPrimPi0PtEventPlane(0),
fhPrimPi0AccPtCentrality(0), fhPrimPi0AccPtEventPlane(0),
fhPrimEtaE(0x0),             fhPrimEtaPt(0x0),             fhPrimEtaPtInCalo(0x0),        
fhPrimEtaAccE(0x0),          fhPrimEtaAccPt(0x0),          fhPrimEtaAccPtPhotonCuts(0x0),
fhPrimEtaY(0x0),             fhPrimEtaAccY(0x0),
fhPrimEtaYeta(0x0),          fhPrimEtaYetaYcut(0x0),       fhPrimEtaAccYeta(0x0),
fhPrimEtaPhi(0x0),           fhPrimEtaAccPhi(0x0),
fhPrimEtaOpeningAngle(0x0),  fhPrimEtaOpeningAnglePhotonCuts(0x0),
fhPrimEtaOpeningAngleAsym(0x0),fhPrimEtaCosOpeningAngle(0x0),
fhPrimEtaPtCentrality(0),    fhPrimEtaPtEventPlane(0),
fhPrimEtaAccPtCentrality(0), fhPrimEtaAccPtEventPlane(0),
fhPrimChHadronPt(0),
fhPrimPi0PtOrigin(0x0),      fhPrimEtaPtOrigin(0x0),
fhPrimNotResonancePi0PtOrigin(0x0),      fhPrimPi0PtStatus(0x0),
fhMCPi0MassPtRec(0x0),       fhMCPi0MassPtTrue(0x0),
fhMCPi0PtTruePtRec(0x0),     fhMCPi0PtTruePtRecMassCut(0x0),
fhMCEtaMassPtRec(0x0),       fhMCEtaMassPtTrue(0x0),
fhMCEtaPtTruePtRec(0x0),     fhMCEtaPtTruePtRecMassCut(0x0),
fhMCPi0MassPtRecCen(0x0),    fhMCPi0MassPtTrueCen(0x0),
fhMCPi0PtTruePtRecCen(0x0),  fhMCPi0PtTruePtRecMassCutCen(0x0),
fhMCEtaMassPtRecCen(0x0),    fhMCEtaMassPtTrueCen(0x0),
fhMCEtaPtTruePtRecCen(0x0),  fhMCEtaPtTruePtRecMassCutCen(0x0),
fhMCPi0PtTruePtRecDifOverPtTrue(0), fhMCPi0PtRecOpenAngle(0),
fhMCEtaPtTruePtRecDifOverPtTrue(0), fhMCEtaPtRecOpenAngle(0),
fhMCPi0PtTruePtRecDifOverPtTrueMassCut(0), fhMCPi0PtRecOpenAngleMassCut(0),
fhMCEtaPtTruePtRecDifOverPtTrueMassCut(0), fhMCEtaPtRecOpenAngleMassCut(0),
fhMCPi0PtTruePtRecDifOverPtTrueCen(0),     fhMCEtaPtTruePtRecDifOverPtTrueCen(0), 
fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut(0), fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut(0),
fhMCPi0PtOrigin(0),          fhMCEtaPtOrigin(0),
fhMCNotResonancePi0PtOrigin(0),fhMCPi0PtStatus(0x0),
fhMCPi0ProdVertex(0),        fhMCEtaProdVertex(0),
fhPrimPi0ProdVertex(0),      fhPrimEtaProdVertex(0),   fhMCPi0Radius(0), fhMCEtaRadius(0),
fhReMCFromConversion(0),     fhReMCFromNotConversion(0),   fhReMCFromMixConversion(0),
fhCosThStarPrimPi0(0),       fhCosThStarPrimEta(0),
fhEPairDiffTime(0),  
fhReSecondaryCellInTimeWindow(0),  fhMiSecondaryCellInTimeWindow(0),           
fhReSecondaryCellOutTimeWindow(0), fhMiSecondaryCellOutTimeWindow(0)             

{
  InitParameters();
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhArmPrimEta[i] = 0;
    fhArmPrimPi0[i] = 0;
  }
  
  for(Int_t ism = 0; ism < 20; ism++)
  {
    fhRealOpeningAnglePerSM [ism] = 0; 
    fhMixedOpeningAnglePerSM[ism] = 0;    
  }
  
  for(Int_t icut = 0; icut < 10; icut++)
  {
    fhReOpAngleBinMinClusterEtaPhi       [icut] = 0; 
    fhReOpAngleBinMaxClusterEtaPhi       [icut] = 0; 
    fhReOpAngleBinMinClusterColRow       [icut] = 0; 
    fhReOpAngleBinMaxClusterColRow       [icut] = 0; 
    fhReOpAngleBinMinClusterEPerSM       [icut] = 0; 
    fhReOpAngleBinMaxClusterEPerSM       [icut] = 0; 
    fhReOpAngleBinMinClusterTimePerSM    [icut] = 0; 
    fhReOpAngleBinMaxClusterTimePerSM    [icut] = 0;     
    fhReOpAngleBinMinClusterNCellPerSM   [icut] = 0; 
    fhReOpAngleBinMaxClusterNCellPerSM   [icut] = 0; 
    fhReOpAngleBinPairClusterRatioPerSM  [icut] = 0;   
    fhReOpAngleBinPairClusterMass        [icut] = 0;   
    fhReOpAngleBinPairClusterMassPerSM   [icut] = 0;   
//  fhReOpAngleBinPairClusterAbsIdMaxCell[icut] = 0;
    
    fhReOpAngleBinPairClusterMassMCTruePi0[icut] = 0;   
    fhReOpAngleBinPairClusterMassMCTrueEta[icut] = 0;   
    fhPrimEtaAccPtOpAngCuts               [icut] = 0;
    fhPrimPi0AccPtOpAngCuts               [icut] = 0;
    
    fhMiOpAngleBinMinClusterEtaPhi       [icut] = 0; 
    fhMiOpAngleBinMaxClusterEtaPhi       [icut] = 0; 
//  fhMiColRowClusterMinOpAngleBin       [icut] = 0; 
//  fhMiOpAngleBinMaxClusterColRow       [icut] = 0; 
    fhMiOpAngleBinMinClusterEPerSM       [icut] = 0; 
    fhMiOpAngleBinMaxClusterEPerSM       [icut] = 0; 
    fhMiOpAngleBinMinClusterTimePerSM    [icut] = 0; 
    fhMiOpAngleBinMaxClusterTimePerSM    [icut] = 0; 
    fhMiOpAngleBinMinClusterNCellPerSM   [icut] = 0; 
    fhMiOpAngleBinMaxClusterNCellPerSM   [icut] = 0; 
    fhMiOpAngleBinPairClusterRatioPerSM  [icut] = 0;
    fhMiOpAngleBinPairClusterMass        [icut] = 0;   
    fhMiOpAngleBinPairClusterMassPerSM   [icut] = 0;   
//  fhMiOpAngleBinPairClusterAbsIdMaxCell[icut] = 0;
    
    fhPtBinClusterEtaPhi                 [icut] = 0; 
    fhPtBinClusterColRow                 [icut] = 0; 
  }
  
  fhReSS[0] = 0;  fhReSS[1] = 0; fhReSS[2] = 0; 
  
  for(Int_t igen = 0; igen < 10; igen++)
  {
    for(Int_t itag = 0; itag < 10; itag++)
    {    
      fhPairGeneratorsBkgMass               [igen][itag] = 0;
      fhPairGeneratorsBkgMassMCPi0          [igen][itag] = 0;
      fhPairGeneratorsBkgCentMCPi0          [igen][itag] = 0;
      fhPairGeneratorsBkgCentMCPi0MassCut   [igen][itag] = 0;
      fhPairGeneratorsBkgEPrimRecoRatioMCPi0[igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoDiffMCPi0 [igen][itag] = 0; 
      fhPairGeneratorsBkgMassMCEta          [igen][itag] = 0; 
      fhPairGeneratorsBkgCentMCEta          [igen][itag] = 0; 
      fhPairGeneratorsBkgCentMCEtaMassCut   [igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoRatioMCEta[igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoDiffMCEta [igen][itag] = 0;
      fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut [igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[igen][itag] = 0; 
      fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut [igen][itag] = 0; 
      
    }
    
    fhPrimPi0PtPerGenerator   [igen] = 0;  
    fhPrimPi0PtInCaloPerGenerator[igen] = 0;  
    fhPrimPi0AccPtPerGenerator[igen] = 0; 
    fhPrimPi0AccPtPhotonCutsPerGenerator[igen] = 0;
    fhPrimPi0PhiPerGenerator  [igen] = 0;
    fhPrimPi0YPerGenerator    [igen] = 0;
    
    fhPrimEtaPtPerGenerator   [igen] = 0;
    fhPrimEtaPtInCaloPerGenerator[igen] = 0;
    fhPrimEtaAccPtPerGenerator[igen] = 0;
    fhPrimEtaAccPtPhotonCutsPerGenerator[igen] = 0;
    fhPrimEtaPhiPerGenerator  [igen] = 0;
    fhPrimEtaYPerGenerator    [igen] = 0;
  }

  for(Int_t i = 0; i < 17; i++)
  {
    fhMCOrgMass    [i] =0 ;
    fhMCOrgAsym    [i] =0 ;
    fhMCOrgDeltaEta[i] =0 ;
    fhMCOrgDeltaPhi[i] =0 ;
  }
  
  for(Int_t i = 0; i < 3; i++)
  {
    fhMCOrgPi0MassPtConversion[i] =0 ;
    fhMCOrgEtaMassPtConversion[i] =0 ;
  }
}

//_____________________
/// Destructor.
//_____________________
AliAnaPi0::~AliAnaPi0()
{
  // Remove event containers
  
  if(DoOwnMix() && fEventsList)
  {
    for(Int_t ic=0; ic<GetNCentrBin(); ic++)
    {
      for(Int_t iz=0; iz<GetNZvertBin(); iz++)
      {
        for(Int_t irp=0; irp<GetNRPBin(); irp++)
        {
          Int_t bin = GetEventMixBin(ic,iz,irp);
          fEventsList[bin]->Delete() ;
          delete fEventsList[bin] ;
        }
      }
    }
    delete[] fEventsList;
  }
}

//______________________________
/// Init parameters when first called the analysis.
/// Set default parameters.
//______________________________
void AliAnaPi0::InitParameters()
{
  SetInputAODName("CaloTrackParticle");
  
  AddToHistogramsName("AnaPi0_");
  
  fUseAngleEDepCut = kFALSE;
  
  fUseAngleCut = kTRUE;
  fAngleCut    = 0.;
  fAngleMaxCut = DegToRad(80.);  // 80 degrees cut, avoid EMCal/DCal combinations
  fUseOneCellSeparation = kFALSE;
  
  fMultiCutAna    = kFALSE;
  fMultiCutAnaAcc = kFALSE;
  fMultiCutAnaSim = kFALSE;
  
  fNPtCuts = 3;
  fPtCuts[0] = 0.; fPtCuts[1] = 0.3;   fPtCuts[2] = 0.5;
  for(Int_t i = fNPtCuts; i < 10; i++) fPtCuts[i] = 0.;
  for(Int_t i = 0       ; i < 10; i++) fPtCutsMax[i] = 20.;
 
  fNAsymCuts = 2;
  fAsymCuts[0] = 1.;  fAsymCuts[1] = 0.7; //fAsymCuts[2] = 0.6; //  fAsymCuts[3] = 0.1;
  for(Int_t i = fNAsymCuts; i < 10; i++)fAsymCuts[i] = 0.;
  
  fNCellNCuts = 3;
  fCellNCuts[0] = 0; fCellNCuts[1] = 1;   fCellNCuts[2] = 2;
  for(Int_t i = fNCellNCuts; i < 10; i++)fCellNCuts[i]  = 0;
  
  fNPIDBits = 2;
  fPIDBits[0] = 0;   fPIDBits[1] = 2; //  fPIDBits[2] = 4; fPIDBits[3] = 6;// check, no cut,  dispersion, neutral, dispersion&&neutral
  for(Int_t i = fNPIDBits; i < 10; i++)fPIDBits[i] = 0;
  
//  fNAngleCutBins = 7;
//  fAngleCutBinsArray[0] = 0.014; fAngleCutBinsArray[1] = 0.035;   fAngleCutBinsArray[2] = 0.07; fAngleCutBinsArray[3] = 0.5;
//  fAngleCutBinsArray[4] = 0.95 ; fAngleCutBinsArray[5] = 1.03 ;   fAngleCutBinsArray[6] = 1.1 ; fAngleCutBinsArray[7] = 1.4 ;
//  for(Int_t i = fNAngleCutBins+1; i < 11; i++)fAngleCutBinsArray[i] = 1000;
  fNAngleCutBins = 10;
  Float_t cellsize = 0.0143;
  fAngleCutBinsArray[0] =  0;           fAngleCutBinsArray[1] = 1*cellsize;   fAngleCutBinsArray[2] = 0.02; // 1.5*cellsize
  fAngleCutBinsArray[3] =  2*cellsize;  fAngleCutBinsArray[4] = 3*cellsize;   fAngleCutBinsArray[5] = 6*cellsize;  
  fAngleCutBinsArray[6] = 12*cellsize;  fAngleCutBinsArray[7] = 24*cellsize;  fAngleCutBinsArray[8] = 48*cellsize; 
  fAngleCutBinsArray[9] = 96*cellsize;  fAngleCutBinsArray[10]= 128*cellsize; 
  
  fPi0MassWindow[0] = 0.10; fPi0MassWindow[1] = 0.25;
  fEtaMassWindow[0] = 0.45; fEtaMassWindow[1] = 0.65;
  
}

//_______________________________________
/// Save parameters used for analysis.
//_______________________________________
TObjString * AliAnaPi0::GetAnalysisCuts()
{
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliAnaPi0 ---:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of bins in Centrality:  %d;",GetNCentrBin()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of bins in Z vert. pos: %d;",GetNZvertBin()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of bins in Reac. Plain: %d;",GetNRPBin()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Depth of event buffer: %d;",GetNMaxEvMix()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Select pairs with their angle: %d, edep %d, min angle %2.3f, max angle %2.3f, 1cell separation %d;",fUseAngleCut, fUseAngleEDepCut,fAngleCut,fAngleMaxCut,fUseOneCellSeparation) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," Asymmetry cuts: n = %d, asymmetry < ",fNAsymCuts) ;
  for(Int_t i = 0; i < fNAsymCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fAsymCuts[i]);
  parList+=onePar ;
  snprintf(onePar,buffersize," PID selection bits: n = %d, PID bit =",fNPIDBits) ;
  for(Int_t i = 0; i < fNPIDBits; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fPIDBits[i]);
  parList+=onePar ;
  snprintf(onePar,buffersize,"Cuts:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Z vertex position: -%f < z < %f;",GetZvertexCut(),GetZvertexCut()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  if(fMultiCutAna || fMultiCutAnaAcc)
  {
    snprintf(onePar, buffersize," pT cuts: n = %d, pt > ",fNPtCuts) ;
    for(Int_t i = 0; i < fNPtCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fPtCuts[i]);
    parList+=onePar ;
    snprintf(onePar,buffersize, " N cell in cluster cuts: n = %d, nCell > ",fNCellNCuts) ;
    for(Int_t i = 0; i < fNCellNCuts; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fCellNCuts[i]);
    parList+=onePar ;
  }
  
  return new TObjString(parList) ;
}

//_________________________________________
/// Create histograms to be saved in output file and
/// store them in fOutputContainer.
//_________________________________________
TList * AliAnaPi0::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName(GetName());
  
  // Set sizes and ranges
  const Int_t buffersize = 255;
  char key[buffersize] ;
  char title[buffersize] ;
  
  Int_t nptbins   = GetHistogramRanges()->GetHistoPtBins();
  Int_t nphibins  = GetHistogramRanges()->GetHistoPhiBins();
  Int_t netabins  = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax   = GetHistogramRanges()->GetHistoPtMax();
  Float_t phimax  = GetHistogramRanges()->GetHistoPhiMax();
  Float_t etamax  = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin   = GetHistogramRanges()->GetHistoPtMin();
  Float_t phimin  = GetHistogramRanges()->GetHistoPhiMin();
  Float_t etamin  = GetHistogramRanges()->GetHistoEtaMin();
  
  Int_t nmassbins = GetHistogramRanges()->GetHistoMassBins();
  Int_t nasymbins = GetHistogramRanges()->GetHistoAsymmetryBins();
  Float_t massmax = GetHistogramRanges()->GetHistoMassMax();
  Float_t asymmax = GetHistogramRanges()->GetHistoAsymmetryMax();
  Float_t massmin = GetHistogramRanges()->GetHistoMassMin();
  Float_t asymmin = GetHistogramRanges()->GetHistoAsymmetryMin();
//  Int_t ntrmbins  = GetHistogramRanges()->GetHistoTrackMultiplicityBins();
//  Int_t ntrmmax   = GetHistogramRanges()->GetHistoTrackMultiplicityMax();
//  Int_t ntrmmin   = GetHistogramRanges()->GetHistoTrackMultiplicityMin();
  Int_t tdbins    = GetHistogramRanges()->GetHistoDiffTimeBins();
  Float_t tdmax   = GetHistogramRanges()->GetHistoDiffTimeMax();
  Float_t tdmin   = GetHistogramRanges()->GetHistoDiffTimeMin();
  Int_t ntimebins = GetHistogramRanges()->GetHistoTimeBins();         
  Float_t timemax = GetHistogramRanges()->GetHistoTimeMax();         
  Float_t timemin = GetHistogramRanges()->GetHistoTimeMin();

  Int_t   nopanbins = GetHistogramRanges()->GetHistoNOpeningAngleBins();
  Float_t opanmin   = GetHistogramRanges()->GetHistoOpeningAngleMin()  ;
  Float_t opanmax   = GetHistogramRanges()->GetHistoOpeningAngleMax()  ;

  Int_t   nratbins = GetHistogramRanges()->GetHistoRatioBins();
  Float_t ratmin   = GetHistogramRanges()->GetHistoRatioMin() ;
  Float_t ratmax   = GetHistogramRanges()->GetHistoRatioMax() ;

  Int_t   ndifbins = GetHistogramRanges()->GetHistoEDiffBins();
  Float_t difmin   = GetHistogramRanges()->GetHistoEDiffMin() ;
  Float_t difmax   = GetHistogramRanges()->GetHistoEDiffMax() ;
  
  Int_t netabinsopen =  TMath::Nint(netabins*4/(etamax-etamin));
  Int_t nphibinsopen = TMath::Nint(nphibins*TMath::TwoPi()/(phimax-phimin));
  
  Int_t   ncenbin  = GetHistogramRanges()->GetHistoCentralityBins()  ;
  Float_t cenmin   = GetHistogramRanges()->GetHistoCentralityMin()   ;
  Float_t cenmax   = GetHistogramRanges()->GetHistoCentralityMax()   ;
  
  InitHistoRangeArrays();

  TArrayD ptBinsArray   = GetHistogramRanges()->GetHistoPtArr();
  TArrayD massBinsArray = GetHistogramRanges()->GetHistoMassArr();
  TArrayD difBinsArray  = GetHistogramRanges()->GetHistoEDiffArr();
  TArrayD cenBinsArray  = GetHistogramRanges()->GetHistoCentralityArr();
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  //Int_t totalSM = fLastModule-fFirstModule+1;
  //printf("N SM %d, first SM %d, last SM %d, total %d\n",fNModules,fFirstModule,fLastModule, totalSM);
  
  // Cell column-row histograms, see base class for data members setting
  //fNMaxColsFull+2,-1.5,fNMaxColsFull+0.5, fNMaxRowsFull+2,-1.5,fNMaxRowsFull+0.5
  Int_t   ncolcell   = fNMaxColsFull+2;
  Float_t colcellmin = -1.5;
  Float_t colcellmax = fNMaxColsFull+0.5;
  
  Int_t   nrowcell   = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  Float_t rowcellmin = fNMaxRowsFullMin-1.5;
  Float_t rowcellmax = fNMaxRowsFullMax+0.5;

  // Start with pure MC kinematics histograms
  // In case other tasks just need this info like AliAnaPi0EbE
  if ( IsDataMC() && IsGeneratedParticlesAnalysisOn() )
  {
    // Pi0
    
    fhPrimPi0E     = new TH1F("hPrimPi0E","Primary #pi^{0} E, |#it{Y}|<1",
                              nptbins,ptmin,ptmax) ;
    fhPrimPi0E   ->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhPrimPi0E) ;

    fhPrimPi0Pt     = new TH1F("hPrimPi0Pt","Primary #pi^{0} #it{p}_{T} , |#it{Y}|<1",
                               nptbins,ptmin,ptmax) ;
    fhPrimPi0Pt   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimPi0Pt) ;
    
    fhPrimPi0Y      = new TH2F("hPrimPi0Rapidity","Rapidity of primary #pi^{0}",
                               nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimPi0Y   ->SetYTitle("#it{Rapidity}");
    fhPrimPi0Y   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimPi0Y) ;
    
    fhPrimPi0Yeta      = new TH2F("hPrimPi0PseudoRapidity","PseudoRapidity of primary #pi^{0}",
                                  nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimPi0Yeta   ->SetYTitle("#eta");
    fhPrimPi0Yeta   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimPi0Yeta) ;
    
    fhPrimPi0YetaYcut      = new TH2F("hPrimPi0PseudoRapidityYcut","PseudoRapidity of primary #pi^{0}, |#it{Y}|<1",
                                      nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimPi0YetaYcut   ->SetYTitle("#eta");
    fhPrimPi0YetaYcut   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimPi0YetaYcut) ;
    
    fhPrimPi0Phi    = new TH2F("hPrimPi0Phi","#varphi of primary #pi^{0}, |#it{Y}|<1",
                               nptbins,ptmin,ptmax,nphibinsopen,0,360) ;
    fhPrimPi0Phi->SetYTitle("#varphi (deg)");
    fhPrimPi0Phi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimPi0Phi) ;
    
    if ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() )
    {      
      fhPrimPi0PtInCalo     = new TH1F("hPrimPi0PtInCalo","Primary #pi^{0} #it{p}_{T} , in calorimeter acceptance",
                                       nptbins,ptmin,ptmax) ;
      fhPrimPi0PtInCalo   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0PtInCalo) ;

      fhPrimPi0AccE  = new TH1F("hPrimPi0AccE","Primary #pi^{0} #it{E} with both photons in acceptance",
                                nptbins,ptmin,ptmax) ;
      fhPrimPi0AccE->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhPrimPi0AccE) ;
      
      fhPrimPi0AccPt  = new TH1F("hPrimPi0AccPt","Primary #pi^{0} #it{p}_{T} with both photons in acceptance",
                                 nptbins,ptmin,ptmax) ;
      fhPrimPi0AccPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0AccPt) ;
      
      fhPrimPi0AccPtPhotonCuts  = new TH1F("hPrimPi0AccPtPhotonCuts","Primary #pi^{0} #it{p}_{T} with both photons in acceptance",
                                           nptbins,ptmin,ptmax) ;
      fhPrimPi0AccPtPhotonCuts->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0AccPtPhotonCuts) ;
      
      fhPrimPi0AccY   = new TH2F("hPrimPi0AccRapidity","Rapidity of primary #pi^{0} with accepted daughters",
                                 nptbins,ptmin,ptmax,netabins,etamin,etamax) ;
      fhPrimPi0AccY->SetYTitle("Rapidity");
      fhPrimPi0AccY->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0AccY) ;
      
      fhPrimPi0AccYeta      = new TH2F("hPrimPi0AccPseudoRapidity","PseudoRapidity of primary #pi^{0} with accepted daughters",
                                       nptbins,ptmin,ptmax,netabins,etamin,etamax) ;
      fhPrimPi0AccYeta   ->SetYTitle("#eta");
      fhPrimPi0AccYeta   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0AccYeta) ;
      
      fhPrimPi0AccPhi = new TH2F("hPrimPi0AccPhi","#varphi of primary #pi^{0} with accepted daughters",
                                 nptbins,ptmin,ptmax,
                                 nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ;
      fhPrimPi0AccPhi->SetYTitle("#varphi (deg)");
      fhPrimPi0AccPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimPi0AccPhi) ;
    }
    
    // Eta
    
    fhPrimEtaE     = new TH1F("hPrimEtaE","Primary eta E",
                              nptbins,ptmin,ptmax) ;
    fhPrimEtaE   ->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhPrimEtaE) ;

    fhPrimEtaPt     = new TH1F("hPrimEtaPt","Primary #eta #it{p}_{T}",
                               nptbins,ptmin,ptmax) ;
    fhPrimEtaPt   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimEtaPt) ;
    
    fhPrimEtaY      = new TH2F("hPrimEtaRapidity","Rapidity of primary #eta",
                               nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimEtaY->SetYTitle("#it{Rapidity}");
    fhPrimEtaY->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimEtaY) ;
    
    fhPrimEtaYeta      = new TH2F("hPrimEtaPseudoRapidityEta","PseudoRapidity of primary #eta",
                                  nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimEtaYeta->SetYTitle("#it{Rapidity}");
    fhPrimEtaYeta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimEtaYeta) ;
    
    fhPrimEtaYetaYcut      = new TH2F("hPrimEtaPseudoRapidityEtaYcut","PseudoRapidity of primary #eta, |#it{Y}|<1",
                                      nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
    fhPrimEtaYetaYcut->SetYTitle("#it{Pseudorapidity}");
    fhPrimEtaYetaYcut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimEtaYetaYcut) ;
        
    fhPrimEtaPhi    = new TH2F("hPrimEtaPhi","Azimuthal of primary #eta",
                               nptbins,ptmin,ptmax, nphibinsopen,0,360) ;
    fhPrimEtaPhi->SetYTitle("#varphi (deg)");
    fhPrimEtaPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPrimEtaPhi) ;
    
    if ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() )
    {
      fhPrimEtaPtInCalo     = new TH1F("hPrimEtaPtInCalo","Primary #eta #it{p}_{T}, in calorimeter acceptance",
                                 nptbins,ptmin,ptmax) ;
      fhPrimEtaPtInCalo   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaPtInCalo) ;
      
      fhPrimEtaAccE  = new TH1F("hPrimEtaAccE","Primary #eta #it{E} with both photons in acceptance",
                                nptbins,ptmin,ptmax) ;
      fhPrimEtaAccE->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhPrimEtaAccE) ;
      
      fhPrimEtaAccPt  = new TH1F("hPrimEtaAccPt","Primary eta #it{p}_{T} with both photons in acceptance",
                                 nptbins,ptmin,ptmax) ;
      fhPrimEtaAccPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaAccPt) ;
      
      fhPrimEtaAccPtPhotonCuts  = new TH1F("hPrimEtaAccPtPhotonCuts","Primary eta #it{p}_{T} with both photons in acceptance",
                                           nptbins,ptmin,ptmax) ;
      fhPrimEtaAccPtPhotonCuts->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaAccPtPhotonCuts) ;
      
      fhPrimEtaAccPhi = new TH2F("hPrimEtaAccPhi","Azimuthal of primary #eta with accepted daughters",
                                 nptbins,ptmin,ptmax, nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ;
      fhPrimEtaAccPhi->SetYTitle("#varphi (deg)");
      fhPrimEtaAccPhi->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaAccPhi) ;
      
      fhPrimEtaAccY   = new TH2F("hPrimEtaAccRapidity","Rapidity of primary #eta",
                                 nptbins,ptmin,ptmax, netabins,etamin,etamax) ;
      fhPrimEtaAccY->SetYTitle("#it{Rapidity}");
      fhPrimEtaAccY->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaAccY) ;
      
      fhPrimEtaAccYeta  = new TH2F("hPrimEtaAccPseudoRapidity","PseudoRapidity of primary #eta",
                                   nptbins,ptmin,ptmax, netabins,etamin,etamax) ;
      fhPrimEtaAccYeta->SetYTitle("#it{Pseudorapidity}");
      fhPrimEtaAccYeta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPrimEtaAccYeta) ;
    }
      
    // Create histograms only for PbPb or high multiplicity analysis analysis
    if ( IsHighMultiplicityAnalysisOn() )
    {      
      fhPrimPi0PtCentrality     = new TH2F
      ("hPrimPi0PtCentrality",
       "Primary #pi^{0} #it{p}_{T} vs reco centrality, |#it{Y}|<1",
       nptbins,ptmin,ptmax, ncenbin, cenmin, cenmax) ;
      fhPrimPi0PtCentrality   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimPi0PtCentrality   ->SetYTitle("Centrality");
      outputContainer->Add(fhPrimPi0PtCentrality) ;
      
      fhPrimEtaPtCentrality     = new TH2F
      ("hPrimEtaPtCentrality",
       "Primary #eta #it{p}_{T} vs reco centrality, |#it{Y}|<1",
       nptbins,ptmin,ptmax, ncenbin, cenmin, cenmax) ;
      fhPrimEtaPtCentrality   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimEtaPtCentrality   ->SetYTitle("Centrality");
      outputContainer->Add(fhPrimEtaPtCentrality) ;
      
      
      fhPrimPi0PtEventPlane     = new TH2F
      ("hPrimPi0PtEventPlane",
       "Primary #pi^{0} #it{p}_{T} vs reco event plane angle, |#it{Y}|<1",
       nptbins,ptmin,ptmax, 100, 0, TMath::Pi()) ;
      fhPrimPi0PtEventPlane   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimPi0PtEventPlane   ->SetYTitle("Event Plane Angle (rad)");
      outputContainer->Add(fhPrimPi0PtEventPlane) ;
      
      
      fhPrimEtaPtEventPlane     = new TH2F
      ("hPrimEtaPtEventPlane",
       "Primary #eta #it{p}_{T} vs reco event plane angle, |#it{Y}|<1",
       nptbins,ptmin,ptmax, 100, 0, TMath::Pi()) ;
      fhPrimEtaPtEventPlane   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimEtaPtEventPlane   ->SetYTitle("Event Plane Angle (rad)");
      outputContainer->Add(fhPrimEtaPtEventPlane) ;

      if ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() )
      {
        fhPrimPi0AccPtCentrality  = new TH2F
        ("hPrimPi0AccPtCentrality",
         "Primary #pi^{0} with both photons in acceptance #it{p}_{T} vs reco centrality",
         nptbins,ptmin,ptmax, ncenbin, cenmin, cenmax) ;
        fhPrimPi0AccPtCentrality->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPrimPi0AccPtCentrality->SetYTitle("Centrality");
        outputContainer->Add(fhPrimPi0AccPtCentrality) ;
        
        fhPrimEtaAccPtCentrality  = new TH2F
        ("hPrimEtaAccPtCentrality",
         "Primary #eta with both photons in acceptance #it{p}_{T} vs reco centrality",
         nptbins,ptmin,ptmax,  ncenbin, cenmin, cenmax) ;
        fhPrimEtaAccPtCentrality->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPrimEtaAccPtCentrality->SetYTitle("Centrality");
        outputContainer->Add(fhPrimEtaAccPtCentrality) ;

        fhPrimPi0AccPtEventPlane  = new TH2F
        ("hPrimPi0AccPtEventPlane",
         "Primary #pi^{0} with both photons in acceptance #it{p}_{T} vs reco event plane angle",
         nptbins,ptmin,ptmax, 100, 0, TMath::Pi()) ;
        fhPrimPi0AccPtEventPlane->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPrimPi0AccPtEventPlane->SetYTitle("Event Plane Angle (rad)");
        outputContainer->Add(fhPrimPi0AccPtEventPlane) ;

        fhPrimEtaAccPtEventPlane  = new TH2F
        ("hPrimEtaAccPtEventPlane",
         "Primary #eta with both #gamma_{decay} in acceptance #it{p}_{T} vs reco event plane angle",
         nptbins,ptmin,ptmax, 100, 0, TMath::Pi()) ;
        fhPrimEtaAccPtEventPlane->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPrimEtaAccPtEventPlane->SetYTitle("Event Plane Angle (rad)");
        outputContainer->Add(fhPrimEtaAccPtEventPlane) ;
      }
    }
    
    if(fFillAngleHisto && ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() ) )
    {
      fhPrimPi0OpeningAngle  = new TH2F
      ("hPrimPi0OpeningAngle","Angle between all primary #gamma pair vs E_{#pi^{0}}, in acceptance",
       nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
      fhPrimPi0OpeningAngle->SetYTitle("#theta(rad)");
      fhPrimPi0OpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhPrimPi0OpeningAngle) ;

      fhPrimPi0OpeningAnglePhotonCuts  = new TH2F
      ("hPrimPi0OpeningAnglePhotonCuts","Angle between all primary #gamma pair vs E_{#pi^{0}} in acceptance",
       nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
      fhPrimPi0OpeningAnglePhotonCuts->SetYTitle("#theta(rad)");
      fhPrimPi0OpeningAnglePhotonCuts->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhPrimPi0OpeningAnglePhotonCuts) ;
      
      fhPrimPi0OpeningAngleAsym  = new TH2F
      ("hPrimPi0OpeningAngleAsym","Angle between all primary #gamma pair vs #it{Asymmetry}, in acceptance, #it{p}_{T}>5 GeV/#it{c}",
       100,0,1,nopanbins,opanmin,opanmax);
      fhPrimPi0OpeningAngleAsym->SetXTitle("|A|=| (E_{1}-E_{2}) / (E_{1}+E_{2}) |");
      fhPrimPi0OpeningAngleAsym->SetYTitle("#theta(rad)");
      outputContainer->Add(fhPrimPi0OpeningAngleAsym) ;
      
      fhPrimPi0CosOpeningAngle  = new TH2F
      ("hPrimPi0CosOpeningAngle","Cosinus of angle between all primary #gamma pair vs E_{#pi^{0}}, in acceptance",
       nptbins,ptmin,ptmax,100,-1,1);
      fhPrimPi0CosOpeningAngle->SetYTitle("cos (#theta) ");
      fhPrimPi0CosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhPrimPi0CosOpeningAngle) ;
      
      fhPrimEtaOpeningAngle  = new TH2F
      ("hPrimEtaOpeningAngle","Angle between all primary #gamma pair vs E_{#eta}, in acceptance",
       nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
      fhPrimEtaOpeningAngle->SetYTitle("#theta(rad)");
      fhPrimEtaOpeningAngle->SetXTitle("E_{#eta} (GeV)");
      outputContainer->Add(fhPrimEtaOpeningAngle) ;

      fhPrimEtaOpeningAnglePhotonCuts  = new TH2F
      ("hPrimEtaOpeningAnglePhotonCuts","Angle between all primary #gamma pair vs E_{#eta}, in acceptance",
       nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
      fhPrimEtaOpeningAnglePhotonCuts->SetYTitle("#theta(rad)");
      fhPrimEtaOpeningAnglePhotonCuts->SetXTitle("E_{#eta} (GeV)");
      outputContainer->Add(fhPrimEtaOpeningAnglePhotonCuts) ;
      
      fhPrimEtaOpeningAngleAsym  = new TH2F
      ("hPrimEtaOpeningAngleAsym","Angle between all primary #gamma pair vs #it{Asymmetry}, #it{p}_{T}>5 GeV/#it{c}, in acceptance",
       100,0,1,nopanbins,opanmin,opanmax);
      fhPrimEtaOpeningAngleAsym->SetXTitle("|#it{A}|=| (#it{E}_{1}-#it{E}_{2}) / (#it{E}_{1}+#it{E}_{2}) |");
      fhPrimEtaOpeningAngleAsym->SetYTitle("#theta(rad)");
      outputContainer->Add(fhPrimEtaOpeningAngleAsym) ;
      
      fhPrimEtaCosOpeningAngle  = new TH2F
      ("hPrimEtaCosOpeningAngle","Cosinus of angle between all primary #gamma pair vs E_{#eta}, in acceptance",
       nptbins,ptmin,ptmax,100,-1,1);
      fhPrimEtaCosOpeningAngle->SetYTitle("cos (#theta) ");
      fhPrimEtaCosOpeningAngle->SetXTitle("#it{E}_{ #eta} (GeV)");
      outputContainer->Add(fhPrimEtaCosOpeningAngle) ;
    }
  
    // Primary origin
    if ( fFillOriginHisto )
    {
      //K+- & pi+-
      fhPrimChHadronPt=new TH2F("hPrimChHadronPt","Primary K+- #pi+-",nptbins,ptmin,ptmax,2,0,2) ;
      fhPrimChHadronPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimChHadronPt->SetYTitle("Origin");
      fhPrimChHadronPt->GetYaxis()->SetBinLabel(1 ,"K+-");
      fhPrimChHadronPt->GetYaxis()->SetBinLabel(2 ,"pi+-");
      outputContainer->Add(fhPrimChHadronPt) ;
      
      // Pi0
      fhPrimPi0PtOrigin     = new TH2F("hPrimPi0PtOrigin","Primary #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,20,0,20) ;
      fhPrimPi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimPi0PtOrigin->SetYTitle("Origin");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances ");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K0S");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(11 ,"K0L");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(12 ,"K+-");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(13 ,"K*");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(14 ,"#Lambda");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(15 ,"hadron int.");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(16 ,"radius");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(17 ,"#pi+-");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(18 ,"e+-");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(19 ,"#mu+-");
      fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(20 ,"p, n");
      outputContainer->Add(fhPrimPi0PtOrigin) ;
      
      fhPrimNotResonancePi0PtOrigin     = new TH2F("hPrimNotResonancePi0PtOrigin","Primary #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,11,0,11) ;
      fhPrimNotResonancePi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimNotResonancePi0PtOrigin->SetYTitle("Origin");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
      fhPrimNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
      outputContainer->Add(fhPrimNotResonancePi0PtOrigin) ;
      
      fhPrimPi0PtStatus     = new TH2F("hPrimPi0PtStatus","Primary #pi^{0} #it{p}_{T} vs status",nptbins,ptmin,ptmax,101,-50,50) ;
      fhPrimPi0PtStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimPi0PtStatus->SetYTitle("Status");
      outputContainer->Add(fhPrimPi0PtStatus) ;

      // Eta
      fhPrimEtaPtOrigin     = new TH2F("hPrimEtaPtOrigin","Primary #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,7,0,7) ;
      fhPrimEtaPtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimEtaPtOrigin->SetYTitle("Origin");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(5 ,"Other");
      fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(6 ,"#eta prime ");
      outputContainer->Add(fhPrimEtaPtOrigin) ;
      
      // Production vertex
      fhPrimPi0ProdVertex = new TH2F("hPrimPi0ProdVertex","generated #pi^{0} #it{p}_{T} vs production vertex",
                                     200,0.,20.,5000,0,500) ;
      fhPrimPi0ProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimPi0ProdVertex->SetYTitle("#it{R} (cm)");
      outputContainer->Add(fhPrimPi0ProdVertex) ;
      
      fhPrimEtaProdVertex = new TH2F("hPrimEtaProdVertex","generated #eta #it{p}_{T} vs production vertex",
                                     200,0.,20.,5000,0,500) ;
      fhPrimEtaProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPrimEtaProdVertex->SetYTitle("#it{R} (cm)");
      outputContainer->Add(fhPrimEtaProdVertex) ;
    }
    
    if ( fFillArmenterosThetaStar && ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() ) )
    {
      TString ebin[] = {"8 < E < 12 GeV","12 < E < 16 GeV", "16 < E < 20 GeV", "E > 20 GeV" };
      Int_t narmbins = 400;
      Float_t armmin = 0;
      Float_t armmax = 0.4;
      
      for(Int_t i = 0; i < 4; i++)
      {
        fhArmPrimPi0[i] =  new TH2F(Form("hArmenterosPrimPi0EBin%d",i),
                                    Form("Armenteros of primary #pi^{0}, %s",ebin[i].Data()),
                                    200, -1, 1, narmbins,armmin,armmax);
        fhArmPrimPi0[i]->SetYTitle("#it{p}_{T}^{Arm}");
        fhArmPrimPi0[i]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmPrimPi0[i]) ;
        
        fhArmPrimEta[i] =  new TH2F(Form("hArmenterosPrimEtaEBin%d",i),
                                    Form("Armenteros of primary #eta, %s",ebin[i].Data()),
                                    200, -1, 1, narmbins,armmin,armmax);
        fhArmPrimEta[i]->SetYTitle("#it{p}_{T}^{Arm}");
        fhArmPrimEta[i]->SetXTitle("#alpha^{Arm}");
        outputContainer->Add(fhArmPrimEta[i]) ;
      }
      
      // Same as asymmetry ...
      fhCosThStarPrimPi0  = new TH2F
      ("hCosThStarPrimPi0","cos(#theta *) for primary #pi^{0}",nptbins,ptmin,ptmax,200,-1,1);
      fhCosThStarPrimPi0->SetYTitle("cos(#theta *)");
      fhCosThStarPrimPi0->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhCosThStarPrimPi0) ;
      
      fhCosThStarPrimEta  = new TH2F
      ("hCosThStarPrimEta","cos(#theta *) for primary #eta",nptbins,ptmin,ptmax,200,-1,1);
      fhCosThStarPrimEta->SetYTitle("cos(#theta *)");
      fhCosThStarPrimEta->SetXTitle("E_{ #eta} (GeV)");
      outputContainer->Add(fhCosThStarPrimEta) ;
    }
    
    if(fFillOnlyMCAcceptanceHisto)  return outputContainer;
  }

  //
  // Create mixed event containers
  //
  fEventsList = new TList*[GetNCentrBin()*GetNZvertBin()*GetNRPBin()] ;
  
  for(Int_t ic=0; ic<GetNCentrBin(); ic++)
  {
    for(Int_t iz=0; iz<GetNZvertBin(); iz++)
    {
      for(Int_t irp=0; irp<GetNRPBin(); irp++)
      {
        Int_t bin = GetEventMixBin(ic,iz,irp);
        fEventsList[bin] = new TList() ;
        fEventsList[bin]->SetOwner(kFALSE);
      }
    }
  }
      
  fhRe1 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
  fhMi1 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
  if ( fFillBadDistHisto )
  {
    fhRe2 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    fhRe3 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    fhMi2 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    fhMi3 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
  }

  if ( fMakeInvPtPlots )
  {
    fhReInvPt1 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    fhMiInvPt1 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    
    if(fFillBadDistHisto)
    {
      fhReInvPt2 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
      fhReInvPt3 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
      fhMiInvPt2 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
      fhMiInvPt3 = new TH2F*[GetNCentrBin()*fNPIDBits*fNAsymCuts] ;
    }
  }
  
  for(Int_t ic=0; ic<GetNCentrBin(); ic++)
  {
    for(Int_t ipid=0; ipid<fNPIDBits; ipid++)
    {
      for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
      {
        Int_t index = ((ic*fNPIDBits)+ipid)*fNAsymCuts + iasym;
        //printf("cen %d, pid %d, asy %d, Index %d\n",ic,ipid,iasym,index);
        //Distance to bad module 1
        snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
        snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                 ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
        fhRe1[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhRe1[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhRe1[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        //printf("name: %s\n ",fhRe1[index]->GetName());
        outputContainer->Add(fhRe1[index]) ;
        
        if(fFillBadDistHisto)
        {
          //Distance to bad module 2
          snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
          snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                   ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
          fhRe2[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRe2[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhRe2[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhRe2[index]) ;
          
          //Distance to bad module 3
          snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
          snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                   ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
          fhRe3[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRe3[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhRe3[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhRe3[index]) ;
        }
        
        //Inverse pT
        if(fMakeInvPtPlots)
        {
          //Distance to bad module 1
          snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
          snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                   ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
          fhReInvPt1[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhReInvPt1[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhReInvPt1[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhReInvPt1[index]) ;
          
          if(fFillBadDistHisto){
            //Distance to bad module 2
            snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhReInvPt2[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhReInvPt2[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhReInvPt2[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhReInvPt2[index]) ;
            
            //Distance to bad module 3
            snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhReInvPt3[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhReInvPt3[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhReInvPt3[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhReInvPt3[index]) ;
          }
        }
        
        if(DoOwnMix())
        {
          //Distance to bad module 1
          snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
          snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                   ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
          fhMi1[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMi1[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMi1[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhMi1[index]) ;
          if(fFillBadDistHisto){
            //Distance to bad module 2
            snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhMi2[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMi2[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMi2[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMi2[index]) ;
            
            //Distance to bad module 3
            snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhMi3[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMi3[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMi3[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMi3[index]) ;
          }
          
          // Inverse pT
          if(fMakeInvPtPlots)
          {
            // Distance to bad module 1
            snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhMiInvPt1[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMiInvPt1[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMiInvPt1[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMiInvPt1[index]) ;
            if(fFillBadDistHisto){
              // Distance to bad module 2
              snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhMiInvPt2[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMiInvPt2[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              fhMiInvPt2[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMiInvPt2[index]) ;
              
              // Distance to bad module 3
              snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f,dist bad 3",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhMiInvPt3[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMiInvPt3[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              fhMiInvPt3[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMiInvPt3[index]) ;
            }
          }
        }
      }
    }
  }

  if ( !IsDataMC() )
  {
    fhEPairDiffTime = new TH2F("hEPairDiffTime","cluster pair time difference vs #it{p}_{T}",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEPairDiffTime->SetXTitle("#it{p}_{T,pair} (GeV/#it{c})");
    fhEPairDiffTime->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEPairDiffTime);
  }
  
  if ( fFillSecondaryCellTiming )
  {
    fhReSecondaryCellInTimeWindow = new TH2F("hReSecondaryCellInTimeWindow","Real Pair, it{t}_{cell} < 50 ns, w_{cell} > 0",
                                             nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhReSecondaryCellInTimeWindow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReSecondaryCellInTimeWindow->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReSecondaryCellInTimeWindow) ;
    
    fhReSecondaryCellOutTimeWindow = new TH2F("hReSecondaryCellOutTimeWindow","Real Pair, it{t}_{cell} < 50 ns, w_{cell} > 0",
                                              nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhReSecondaryCellOutTimeWindow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReSecondaryCellOutTimeWindow->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReSecondaryCellOutTimeWindow) ;
    
    if ( DoOwnMix() )
    {
      fhMiSecondaryCellInTimeWindow = new TH2F("hMiSecondaryCellInTimeWindow","Real Pair, it{t}_{cell} < 50 ns, w_{cell} > 0",
                                               nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMiSecondaryCellInTimeWindow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiSecondaryCellInTimeWindow->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhMiSecondaryCellInTimeWindow) ;
      
      fhMiSecondaryCellOutTimeWindow = new TH2F("hMiSecondaryCellOutTimeWindow","Real Pair, it{t}_{cell} < 50 ns, w_{cell} > 0",
                                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMiSecondaryCellOutTimeWindow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiSecondaryCellOutTimeWindow->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhMiSecondaryCellOutTimeWindow) ;
      
    }
  }
  
  if ( fCheckConversion )
  {
    fhReConv = new TH2F("hReConv","Real Pair with one recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhReConv->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReConv->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReConv) ;
    
    fhReConv2 = new TH2F("hReConv2","Real Pair with 2 recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhReConv2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReConv2->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReConv2) ;
    
    if ( DoOwnMix() )
    {
      fhMiConv = new TH2F("hMiConv","Mixed Pair with one recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMiConv->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiConv->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhMiConv) ;
      
      fhMiConv2 = new TH2F("hMiConv2","Mixed Pair with 2 recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMiConv2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiConv2->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhMiConv2) ;
    }
  }

  if ( fFillAsymmetryHisto )
  {
    fhRePtAsym = new TH2F("hRePtAsym","#it{Asymmetry} vs #it{p}_{T}, for pairs",
                          nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
    fhRePtAsym->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhRePtAsym->SetYTitle("#it{Asymmetry}");
    outputContainer->Add(fhRePtAsym);
    
    fhRePtAsymPi0 = new TH2F("hRePtAsymPi0",Form("#it{Asymmetry} vs #it{p}_{T}, for pairs %2.2f<M<%2.2f MeV/#it{c}^{2}",
                                                 fPi0MassWindow[0],fPi0MassWindow[1]),
                             nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
    fhRePtAsymPi0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhRePtAsymPi0->SetYTitle("Asymmetry");
    outputContainer->Add(fhRePtAsymPi0);
    
    fhRePtAsymEta = new TH2F("hRePtAsymEta",Form("#it{Asymmetry} vs #it{p}_{T}, for pairs %2.2f<M<%2.2f MeV/#it{c}^{2}",
                                                 fEtaMassWindow[0],fEtaMassWindow[1]),
                             nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
    fhRePtAsymEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhRePtAsymEta->SetYTitle("#it{Asymmetry}");
    outputContainer->Add(fhRePtAsymEta);
    
    if(DoOwnMix())
    {
      fhMiPtAsym = new TH2F("hMiPtAsym","#it{Asymmetry} vs #it{p}_{T}, for mixed pairs",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
      fhMiPtAsym->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiPtAsym->SetYTitle("#it{Asymmetry}");
      outputContainer->Add(fhMiPtAsym);
      
      fhMiPtAsymPi0 = new TH2F("hMiPtAsymPi0",Form("#it{Asymmetry} vs #it{p}_{T}, for mixed pairs %2.2f<M<%2.2f MeV/#it{c}^{2}",
                                                   fPi0MassWindow[0],fPi0MassWindow[1]),
                               nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
      fhMiPtAsymPi0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiPtAsymPi0->SetYTitle("Asymmetry");
      outputContainer->Add(fhMiPtAsymPi0);
      
      fhMiPtAsymEta = new TH2F("hMiPtAsymEta",
                               Form("#it{Asymmetry} vs #it{p}_{T}, for mixed pairs %2.2f<M<%2.2f MeV/#it{c}^{2}",
                                    fEtaMassWindow[0],fEtaMassWindow[1]),
                               nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
      fhMiPtAsymEta->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMiPtAsymEta->SetYTitle("#it{Asymmetry}");
      outputContainer->Add(fhMiPtAsymEta);
    }
  }
  
  if  ( fMultiCutAnaAcc )
  {
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
    {
      fhPtBinClusterEtaPhi[ipt] = new TH2F
      (Form("hPtBin%d_Cluster_EtaPhi",ipt),
       Form("#eta vs #varphi, %2.2f<#it{p}_{T}<%2.2f GeV/#it{c}",fPtCuts[ipt],fPtCuts[ipt+1]),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhPtBinClusterEtaPhi[ipt]->SetYTitle("#varphi (rad)");
      fhPtBinClusterEtaPhi[ipt]->SetXTitle("#eta");
      outputContainer->Add(fhPtBinClusterEtaPhi[ipt]) ;
      
      fhPtBinClusterColRow[ipt] = new TH2F
      (Form("hPtBin%d_Cluster_ColRow",ipt),
       Form("column vs row, %2.2f<#it{p}_{T}<%2.2f GeV/#it{c}",fPtCuts[ipt],fPtCuts[ipt+1]),
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhPtBinClusterColRow[ipt]->SetYTitle("row");
      fhPtBinClusterColRow[ipt]->SetXTitle("column");
      outputContainer->Add(fhPtBinClusterColRow[ipt]) ;
    }
  }

  if ( fMultiCutAna )
  {    
    fhRePtNCellAsymCuts    = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhMiPtNCellAsymCuts    = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];

    if(fFillAngleHisto)
    {
      fhRePtNCellAsymCutsOpAngle = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMiPtNCellAsymCutsOpAngle = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    }
    
    if(fFillSMCombinations)
    {
      for(Int_t iSM = 0; iSM < fNModules; iSM++) 
      {
        if ( iSM < fFirstModule || iSM > fLastModule ) 
        {
          fhRePtNCellAsymCutsSM[iSM] = 0x0;
          if ( fFillAngleHisto ) 
            fhRePtNCellAsymCutsSMOpAngle[iSM] = 0x0;
          continue;
        }
        
        fhRePtNCellAsymCutsSM[iSM] = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        if(fFillAngleHisto) fhRePtNCellAsymCutsSMOpAngle[iSM] = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      }
    }
        
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
    {
      for(Int_t icell=0; icell<fNCellNCuts; icell++)
      {
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
        {
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f ",
                   fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          
          Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
          //printf("ipt %d, icell %d, iassym %d, index %d\n",ipt, icell, iasym, index);
          
          fhRePtNCellAsymCuts[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCuts[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhRePtNCellAsymCuts[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhRePtNCellAsymCuts[index]) ;
          
          if(DoOwnMix())
          {
            snprintf(key,   buffersize,"hMi_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
            snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f",
                     fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
            fhMiPtNCellAsymCuts[index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMiPtNCellAsymCuts[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMiPtNCellAsymCuts[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMiPtNCellAsymCuts[index]) ;
          }
          
          if ( fFillSMCombinations )
          {
            for(Int_t iSM = 0; iSM < fNModules; iSM++)
            {
              //printf("\t sm %d\n",iSM);
              if ( iSM < fFirstModule || iSM > fLastModule ) 
              {
                fhRePtNCellAsymCutsSM[iSM][index] = 0x0;
                continue;
              }
              
              snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d_SM%d",ipt,icell,iasym,iSM) ;
              snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f, SM %d ",
                       fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym],iSM) ;
              fhRePtNCellAsymCutsSM[iSM][index] = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhRePtNCellAsymCutsSM[iSM][index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              fhRePtNCellAsymCutsSM[iSM][index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhRePtNCellAsymCutsSM[iSM][index]) ;
            }
            
            if ( fFillAngleHisto )
            {
              snprintf(key,   buffersize,"hReOpAngle_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
              snprintf(title, buffersize,"Real #theta_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f ",
                       fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
              
              Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
              //printf("Angle ipt %d, icell %d, iassym %d, index %d, %s, %s\n",ipt, icell, iasym, index, key, title);
              
              fhRePtNCellAsymCutsOpAngle[index] = new TH2F(key,title,nptbins,ptmin,ptmax,280,0,1.4) ;
              fhRePtNCellAsymCutsOpAngle[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              fhRePtNCellAsymCutsOpAngle[index]->SetYTitle("#theta_{#gamma,#gamma} (rad)");
              outputContainer->Add(fhRePtNCellAsymCutsOpAngle[index]) ;
              
              if ( DoOwnMix() )
              {
                snprintf(key,   buffersize,"hMiOpAngle_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
                snprintf(title, buffersize,"Mixed #theta_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f",
                         fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
                fhMiPtNCellAsymCutsOpAngle[index] = new TH2F(key,title,nptbins,ptmin,ptmax,280,0,1.4) ;
                fhMiPtNCellAsymCutsOpAngle[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
                fhMiPtNCellAsymCutsOpAngle[index]->SetYTitle("#theta_{#gamma,#gamma} (rad)");
                outputContainer->Add(fhMiPtNCellAsymCutsOpAngle[index]) ;
              }
            
              if ( fFillSMCombinations )
              {
                for(Int_t iSM = 0; iSM < fNModules; iSM++)
                {
                  if ( iSM < fFirstModule || iSM > fLastModule )
                  {
                    fhRePtNCellAsymCutsSMOpAngle[iSM][index] = 0x0;
                    continue;
                  }
                  
                  snprintf(key,   buffersize,"hReOpAngle_pt%d_cell%d_asym%d_SM%d",ipt,icell,iasym,iSM) ;
                  snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for %1.1f< #it{p}_{T} < %1.1f, ncell>%d and asym<%1.2f, SM %d ",
                           fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym],iSM) ;
                  fhRePtNCellAsymCutsSMOpAngle[iSM][index] = new TH2F(key,title,nptbins,ptmin,ptmax,280,0,1.4) ;
                  fhRePtNCellAsymCutsSMOpAngle[iSM][index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
                  fhRePtNCellAsymCutsSMOpAngle[iSM][index]->SetYTitle("#theta_{#gamma,#gamma} (rad)");
                  outputContainer->Add(fhRePtNCellAsymCutsSMOpAngle[iSM][index]) ;
                }
              }
            }
          }
        }
      }
    }
    
  } // multi cuts analysis
  
  if ( fFillSSCombinations )
  {
    fhReSS[0] = new TH2F("hRe_SS_Tight"," 0.01 < #lambda_{0}^{2} < 0.4",
                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
    fhReSS[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReSS[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReSS[0]) ;
    
    
    fhReSS[1] = new TH2F("hRe_SS_Loose"," #lambda_{0}^{2} > 0.4",
                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
    fhReSS[1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReSS[1]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReSS[1]) ;
    
    
    fhReSS[2] = new TH2F("hRe_SS_Both"," cluster_{1} #lambda_{0}^{2} > 0.4; cluster_{2} 0.01 < #lambda_{0}^{2} < 0.4",
                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
    fhReSS[2]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhReSS[2]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
    outputContainer->Add(fhReSS[2]) ;
  }
  
  if ( DoOwnMix() )
  {
    fhEventBin=new TH1I("hEventBin","Number of real pairs per bin(cen,vz,rp)",
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventBin->SetXTitle("bin");
    outputContainer->Add(fhEventBin) ;
    
    
    fhEventMixBin=new TH1I("hEventMixBin","Number of mixed pairs per bin(cen,vz,rp)",
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventMixBin->SetXTitle("bin");
    outputContainer->Add(fhEventMixBin) ;
  }
  
  if ( IsHighMultiplicityAnalysisOn() )
  {
    fhCentrality=new TH1F("hCentralityBin","Number of events in centrality bin",GetNCentrBin(),0.,1.*GetNCentrBin()) ;
    fhCentrality->SetXTitle("Centrality bin");
    outputContainer->Add(fhCentrality) ;
    
    fhCentralityNoPair=new TH1F("hCentralityBinNoPair","Number of events in centrality bin, with no cluster pairs",GetNCentrBin(),0.,1.*GetNCentrBin()) ;
    fhCentralityNoPair->SetXTitle("Centrality bin");
    outputContainer->Add(fhCentralityNoPair) ;
    
    fhEventPlaneResolution=new TH2F("hEventPlaneResolution","Event plane resolution",GetNCentrBin(),0,GetNCentrBin(),100,0.,TMath::TwoPi()) ;
    fhEventPlaneResolution->SetYTitle("Resolution");
    fhEventPlaneResolution->SetXTitle("Centrality Bin");
    outputContainer->Add(fhEventPlaneResolution) ;
  }
  
  if ( fFillAngleHisto )
  {
    fhRealOpeningAngle  = new TH2F
    ("hRealOpeningAngle","Angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
    fhRealOpeningAngle->SetYTitle("#theta(rad)");
    fhRealOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhRealOpeningAngle) ;
    
    fhRealCosOpeningAngle  = new TH2F
    ("hRealCosOpeningAngle","Cosinus of angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,0,1);
    fhRealCosOpeningAngle->SetYTitle("cos (#theta) ");
    fhRealCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhRealCosOpeningAngle) ;
    
    if(DoOwnMix())
    {
      fhMixedOpeningAngle  = new TH2F
      ("hMixedOpeningAngle","Angle between all #gamma pair vs E_{#pi^{0}}, Mixed pairs",nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
      fhMixedOpeningAngle->SetYTitle("#theta(rad)");
      fhMixedOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhMixedOpeningAngle) ;
      
      fhMixedCosOpeningAngle  = new TH2F
      ("hMixedCosOpeningAngle","Cosinus of angle between all #gamma pair vs E_{#pi^{0}}, Mixed pairs",nptbins,ptmin,ptmax,100,0,1);
      fhMixedCosOpeningAngle->SetYTitle("cos (#theta) ");
      fhMixedCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
      outputContainer->Add(fhMixedCosOpeningAngle) ;
    }
    
    if(fFillSMCombinations)
    {
      for(Int_t ism = 0; ism < 20; ism++)
      {
        fhRealOpeningAnglePerSM[ism]  = new TH2F
        (Form("hRealOpeningAngleMod_%d",ism),
         Form("Angle between all #gamma pair vs E_{#pi^{0}}, SM %d",ism),
         nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
        fhRealOpeningAnglePerSM[ism]->SetYTitle("#theta(rad)");
        fhRealOpeningAnglePerSM[ism]->SetXTitle("E_{ #pi^{0}} (GeV)");
        outputContainer->Add(fhRealOpeningAnglePerSM[ism]) ;
        
        if(DoOwnMix())
        {
          fhMixedOpeningAnglePerSM[ism]  = new TH2F
          (Form("hMixedOpeningAngleMod_%d",ism),
           Form("Angle between all #gamma pair vs E_{#pi^{0}}, Mixed pairs, SM %d",ism),
           nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax);
          fhMixedOpeningAnglePerSM[ism]->SetYTitle("#theta(rad)");
          fhMixedOpeningAnglePerSM[ism]->SetXTitle("E_{ #pi^{0}} (GeV)");
          outputContainer->Add(fhMixedOpeningAnglePerSM[ism]) ;
        }
      }
    }
  }
  
  // Histograms filled only if MC data is requested
  if ( IsDataMC() )
  {
    if ( fCheckConversion )
    {
      fhReMCFromConversion = new TH2F("hReMCFromConversion","Invariant mass of 2 clusters originated in conversions",
                                      nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhReMCFromConversion->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhReMCFromConversion->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhReMCFromConversion) ;
      
      fhReMCFromNotConversion = new TH2F("hReMCNotFromConversion","Invariant mass of 2 clusters not originated in conversions",
                                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhReMCFromNotConversion->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhReMCFromNotConversion->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhReMCFromNotConversion) ;
      
      fhReMCFromMixConversion = new TH2F("hReMCFromMixConversion","Invariant mass of 2 clusters one from conversion and the other not",
                                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhReMCFromMixConversion->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhReMCFromMixConversion->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      outputContainer->Add(fhReMCFromMixConversion) ;
    }
    
    if ( fFillOriginHisto && !IsHighMultiplicityAnalysisOn() )
    {
      fhMCPi0PtOrigin     = new TH2F("hMCPi0PtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,20,0,20) ;
      fhMCPi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCPi0PtOrigin->SetYTitle("Origin");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K0S");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(11 ,"K0L");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(12 ,"K+-");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(13 ,"K*");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(14 ,"#Lambda");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(15 ,"hadron int.");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(16 ,"radius");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(17 ,"p, n, #pi+-");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(18 ,"e+-");
      fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(19 ,"#mu+-");
      outputContainer->Add(fhMCPi0PtOrigin) ;
      
      fhMCNotResonancePi0PtOrigin     = new TH2F("hMCNotResonancePi0PtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,11,0,11) ;
      fhMCNotResonancePi0PtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCNotResonancePi0PtOrigin->SetYTitle("Origin");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
      fhMCNotResonancePi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
      outputContainer->Add(fhMCNotResonancePi0PtOrigin) ;
      
      fhMCPi0PtStatus     = new TH2F("hMCPi0PtStatus","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs status",nptbins,ptmin,ptmax,101,-50,50) ;
      fhMCPi0PtStatus->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCPi0PtStatus->SetYTitle("Status");
      outputContainer->Add(fhMCPi0PtStatus) ;
      
      // Eta
      
      fhMCEtaPtOrigin     = new TH2F("hMCEtaPtOrigin","Reconstructed pair from generated #pi^{0} #it{p}_{T} vs origin",nptbins,ptmin,ptmax,7,0,7) ;
      fhMCEtaPtOrigin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCEtaPtOrigin->SetYTitle("Origin");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(5 ,"Other");
      fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(6 ,"#eta prime");
      outputContainer->Add(fhMCEtaPtOrigin) ;
      
      fhMCPi0ProdVertex = new TH2F("hMCPi0ProdVertex","Selected reco pair from generated #pi^{0} #it{p}_{T} vs production vertex",
                                   200,0.,20.,5000,0,500) ;
      fhMCPi0ProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCPi0ProdVertex->SetYTitle("#it{R} (cm)");
      outputContainer->Add(fhMCPi0ProdVertex) ;
      
      fhMCEtaProdVertex = new TH2F("hMCEtaProdVertex","Selected reco pair from generated #eta #it{p}_{T} vs production vertex",
                                   200,0.,20.,5000,0,500) ;
      fhMCEtaProdVertex->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhMCEtaProdVertex->SetYTitle("#it{R} (cm)");
      outputContainer->Add(fhMCEtaProdVertex) ;
      
      fhMCPi0Radius = new TH2F("hPrimPi0Radius",
                               "Production radius of reconstructed pair from generated #pi^{0} corrected by vertex",
                               200,0,20,5000,0,500) ;
      fhMCPi0Radius->SetYTitle("#it{R} (cm)");
      fhMCPi0Radius->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPi0Radius) ;
      
      fhMCEtaRadius = new TH2F("hPrimEtaRadius",
                               "Production radius of reconstructed pair from generated #eta corrected by vertex",
                               200,0,20,5000,0,500) ;
      fhMCEtaRadius->SetYTitle("#it{R} (cm)");
      fhMCEtaRadius->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCEtaRadius) ;
      
      if ( !fFillOriginHistoForMesonsOnly ) 
      {
        // Name of found ancestors of the cluster pairs. Check definitions in FillMCVsRecDataHistograms
        TString ancestorTitle[] = {"Photon, conversion","Electron, conversion",
          "Pi0","Eta","AntiProton","AntiNeutron","Muon & converted stable particles",
          "Resonances","Strings","Initial state interaction","Final state radiations","Colliding protons",
          "Pi0Not2SingleGamma","EtaNot2Gamma","Photon, not both conversion","Electron, not both conversion","not found"};
        
        for(Int_t i = 0; i<17; i++)
        {        
          fhMCOrgMass[i] = new TH2F(Form("hMCOrgMass_%d",i),Form("#it{M} vs #it{p}_{T}, ancestor %s",ancestorTitle[i].Data()),
                                    nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMCOrgMass[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCOrgMass[i]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhMCOrgMass[i]) ;
          
          fhMCOrgAsym[i]= new TH2F(Form("hMCOrgAsym_%d",i),Form("#it{Asymmetry} vs #it{p}_{T}, ancestor %s",ancestorTitle[i].Data()),
                                   nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
          fhMCOrgAsym[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCOrgAsym[i]->SetYTitle("A");
          outputContainer->Add(fhMCOrgAsym[i]) ;
          
          fhMCOrgDeltaEta[i] = new TH2F(Form("hMCOrgDeltaEta_%d",i),Form("#Delta #eta of pair vs #it{p}_{T}, ancestor %s",ancestorTitle[i].Data()),
                                        nptbins,ptmin,ptmax,netabins,-1.4,1.4) ;
          fhMCOrgDeltaEta[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCOrgDeltaEta[i]->SetYTitle("#Delta #eta");
          outputContainer->Add(fhMCOrgDeltaEta[i]) ;
          
          fhMCOrgDeltaPhi[i]= new TH2F(Form("hMCOrgDeltaPhi_%d",i),Form("#Delta #varphi of pair vs #it{p}_{T}, ancestor %s",ancestorTitle[i].Data()),
                                       nptbins,ptmin,ptmax,nphibins,-0.7,0.7) ;
          fhMCOrgDeltaPhi[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMCOrgDeltaPhi[i]->SetYTitle("#Delta #varphi (rad)");
          outputContainer->Add(fhMCOrgDeltaPhi[i]) ;
        }
      }
      
      if ( fCheckConversion )
      {
        //reconstructed gamma clusters coming pi0 (validated in MC true) both not from conversion, with one or two conversions
        fhMCOrgPi0MassPtConversion[0] = new TH2F("hMCOrgPi0MassPtConversion0","Invariant mass of 2 clusters (ancestor #pi^{0}) not originated in conversions",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgPi0MassPtConversion[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgPi0MassPtConversion[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgPi0MassPtConversion[0]) ;
        
        fhMCOrgPi0MassPtConversion[1] = new TH2F("hMCOrgPi0MassPtConversion1","Invariant mass of 2 clusters (ancestor #pi^{0}) one from conversion and the other not",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgPi0MassPtConversion[1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgPi0MassPtConversion[1]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgPi0MassPtConversion[1]) ;
        
        fhMCOrgPi0MassPtConversion[2] = new TH2F("hMCOrgPi0MassPtConversion2","Invariant mass of 2 clusters (ancestor #pi^{0}) originated in conversions",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgPi0MassPtConversion[2]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgPi0MassPtConversion[2]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgPi0MassPtConversion[2]) ;
        
        //reconstructed gamma clusters coming eta (validated in MC true) both not from conversion, with one or two conversions
        fhMCOrgEtaMassPtConversion[0] = new TH2F("hMCOrgEtaMassPtConversion0","Invariant mass of 2 clusters (ancestor #eta) not originated in conversions",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgEtaMassPtConversion[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgEtaMassPtConversion[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgEtaMassPtConversion[0]) ;
        
        fhMCOrgEtaMassPtConversion[1] = new TH2F("hMCOrgEtaMassPtConversion1","Invariant mass of 2 clusters (ancestor #eta) one from conversion and the other not",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgEtaMassPtConversion[1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgEtaMassPtConversion[1]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgEtaMassPtConversion[1]) ;
        
        fhMCOrgEtaMassPtConversion[2] = new TH2F("hMCOrgEtaMassPtConversion2","Invariant mass of 2 clusters (ancestor #eta) originated in conversions",
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCOrgEtaMassPtConversion[2]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhMCOrgEtaMassPtConversion[2]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCOrgEtaMassPtConversion[2]) ;
      }
      
      
      if ( fMultiCutAnaSim )
      {
        fhMCPi0MassPtTrue  = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCPi0MassPtRec   = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCPi0PtTruePtRec = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCEtaMassPtRec   = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCEtaMassPtTrue  = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCEtaPtTruePtRec = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        
        fhMCPi0PtTruePtRecMassCut = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        fhMCEtaPtTruePtRecMassCut = new TH2F*[fNPtCuts*fNAsymCuts*fNCellNCuts];
        
        for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
        {
          for(Int_t icell=0; icell<fNCellNCuts; icell++)
          {
            for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
            {
              Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
              
              fhMCPi0MassPtRec[index] = new TH2F(Form("hMCPi0MassPtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                 Form("Reconstructed #it{M} vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|<%1.2f",
                                                      fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMCPi0MassPtRec[index]->SetXTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
              fhMCPi0MassPtRec[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMCPi0MassPtRec[index]) ;
              
              fhMCPi0MassPtTrue[index] = new TH2F(Form("hMCPi0MassPtTrue_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                  Form("Reconstructed #it{M} vs generated #it{p}_{T} of true #pi^{0} cluster pairs for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|<%1.2f",
                                                       fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                  nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMCPi0MassPtTrue[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCPi0MassPtTrue[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMCPi0MassPtTrue[index]) ;
              
              fhMCPi0PtTruePtRec[index] = new TH2F(Form("hMCPi0PtTruePtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                   Form("Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|>%1.2f",
                                                        fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
              fhMCPi0PtTruePtRec[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCPi0PtTruePtRec[index]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
              outputContainer->Add(fhMCPi0PtTruePtRec[index]) ;
              
              fhMCPi0PtTruePtRecMassCut[index] = new TH2F(Form("hMCPi0PtTruePtRecMassCut_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                          Form("Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2} for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|>%1.2f",
                                                               fPi0MassWindow[0],fPi0MassWindow[1],fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                          nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
              fhMCPi0PtTruePtRecMassCut[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCPi0PtTruePtRecMassCut[index]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
              outputContainer->Add(fhMCPi0PtTruePtRecMassCut[index]) ;
              
              fhMCEtaMassPtRec[index] = new TH2F(Form("hMCEtaMassPtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                 Form("Reconstructed #it{M} vs reconstructed #it{p}_{T} of true #eta cluster pairs for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|<%1.2f",
                                                      fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMCEtaMassPtRec[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCEtaMassPtRec[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMCEtaMassPtRec[index]) ;
              
              fhMCEtaMassPtTrue[index] = new TH2F(Form("hMCEtaMassPtTrue_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                  Form("Reconstructed #it{M} vs generated #it{p}_{T} of true #eta cluster pairs for %1.1f<#it{p}_{T}<%1.1f, #it{N}^{cluster}_{cell}>%d and |#it{A}|<%1.2f",
                                                       fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                  nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMCEtaMassPtTrue[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCEtaMassPtTrue[index]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
              outputContainer->Add(fhMCEtaMassPtTrue[index]) ;
              
              fhMCEtaPtTruePtRec[index] = new TH2F(Form("hMCEtaPtTruePtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                   Form("Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs for %1.1f<#it{p}_{T}<%1.1f, ncell>%d and asym<%1.2f",
                                                        fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
              fhMCEtaPtTruePtRec[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCEtaPtTruePtRec[index]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
              outputContainer->Add(fhMCEtaPtTruePtRec[index]) ;
              
              
              fhMCEtaPtTruePtRecMassCut[index] = new TH2F(Form("hMCEtaPtTruePtRecMassCut_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                          Form("Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2} for %1.1f<#it{p}_{T}<%1.1f, ncell>%d and asym<%1.2f",
                                                               fEtaMassWindow[0],fEtaMassWindow[1],fPtCuts[ipt],fPtCutsMax[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                          nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
              fhMCEtaPtTruePtRecMassCut[index]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
              fhMCEtaPtTruePtRecMassCut[index]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
              outputContainer->Add(fhMCEtaPtTruePtRecMassCut[index]) ;
            }
          }
        }
      }//multi cut ana
      else
      {        
        fhMCPi0MassPtTrue  = new TH2F*[1];
        fhMCPi0MassPtRec   = new TH2F*[1];
        fhMCPi0PtTruePtRec = new TH2F*[1];
        fhMCEtaMassPtTrue  = new TH2F*[1];
        fhMCEtaMassPtRec   = new TH2F*[1];
        fhMCEtaPtTruePtRec = new TH2F*[1];
        
        fhMCPi0PtTruePtRecMassCut = new TH2F*[1];
        fhMCEtaPtTruePtRecMassCut = new TH2F*[1];
        
        fhMCPi0MassPtTrue[0] = new TH2F
        ("hMCPi0MassPtTrue",
         "Reconstructed Mass vs generated #it{p}_{T} of true #pi^{0} cluster pairs",
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhMCPi0MassPtTrue[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCPi0MassPtTrue[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCPi0MassPtTrue[0]) ;
        
        fhMCPi0MassPtRec[0] = new TH2F
        ("hMCPi0MassPtRec",
         "Reconstructed Mass vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs",
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhMCPi0MassPtRec[0]->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCPi0MassPtRec[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCPi0MassPtRec[0]) ;
        
        fhMCPi0PtTruePtRec[0]= new TH2F
        ("hMCPi0PtTruePtRec",
         "Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs",
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
        fhMCPi0PtTruePtRec[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCPi0PtTruePtRec[0]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0PtTruePtRec[0]) ;
        
        fhMCPi0PtTruePtRecMassCut[0]= new TH2F
        ("hMCPi0PtTruePtRecMassCut",
         
         Form("Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}",fPi0MassWindow[0],fPi0MassWindow[1]),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
        fhMCPi0PtTruePtRecMassCut[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCPi0PtTruePtRecMassCut[0]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
        outputContainer->Add(fhMCPi0PtTruePtRecMassCut[0]) ;
        
        fhMCEtaMassPtTrue[0] = new TH2F
        ("hMCEtaMassPtTrue",
         "Reconstructed Mass vs generated #it{p}_{T} of true #eta cluster pairs",
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhMCEtaMassPtTrue[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCEtaMassPtTrue[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCEtaMassPtTrue[0]) ;
        
        fhMCEtaMassPtRec[0] = new TH2F
        ("hMCEtaMassPtRec",
         "Reconstructed Mass vs reconstructed #it{p}_{T} of true #eta cluster pairs",
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhMCEtaMassPtRec[0]->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCEtaMassPtRec[0]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhMCEtaMassPtRec[0]) ;
        
        fhMCEtaPtTruePtRec[0]= new TH2F
        ("hMCEtaPtTruePtRec",
         "Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs",
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
        fhMCEtaPtTruePtRec[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCEtaPtTruePtRec[0]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaPtTruePtRec[0]) ;       
        
        fhMCEtaPtTruePtRecMassCut[0]= new TH2F
        ("hMCEtaPtTruePtRecMassCut",
         Form("Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}",
              fEtaMassWindow[0],fEtaMassWindow[1]), 
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
        fhMCEtaPtTruePtRecMassCut[0]->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
        fhMCEtaPtTruePtRecMassCut[0]->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
        outputContainer->Add(fhMCEtaPtTruePtRecMassCut[0]) ;
        
        fhMCPi0PtTruePtRecDifOverPtTrue = new TH2F
        ("hMCPi0PtTruePtRecDifOverPtTrue",
         "(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} of true #pi^{0} cluster pairs",
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax) ;
        fhMCPi0PtTruePtRecDifOverPtTrue->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCPi0PtTruePtRecDifOverPtTrue->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
        outputContainer->Add(fhMCPi0PtTruePtRecDifOverPtTrue) ;
        
        fhMCPi0PtTruePtRecDifOverPtTrueMassCut = new TH2F
        ("hMCPi0PtTruePtRecDifOverPtTrueMassCut",
         Form("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} of true #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
              fPi0MassWindow[0],fPi0MassWindow[1]),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax) ;
        fhMCPi0PtTruePtRecDifOverPtTrueMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCPi0PtTruePtRecDifOverPtTrueMassCut->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} ");
        outputContainer->Add(fhMCPi0PtTruePtRecDifOverPtTrueMassCut) ;
        
        fhMCEtaPtTruePtRecDifOverPtTrue = new TH2F
        ("hMCEtaPtTruePtRecDifOverPtTrue",
         "(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}  of true #eta cluster pairs",
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax) ;
        fhMCEtaPtTruePtRecDifOverPtTrue->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCEtaPtTruePtRecDifOverPtTrue->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
        outputContainer->Add(fhMCEtaPtTruePtRecDifOverPtTrue) ;
        
        fhMCEtaPtTruePtRecDifOverPtTrueMassCut = new TH2F
        ("hMCEtaPtTruePtRecDifOverPtTrueCutMassCut",
         Form("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} of true #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
              fEtaMassWindow[0],fEtaMassWindow[1]),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax) ;
        fhMCEtaPtTruePtRecDifOverPtTrueMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCEtaPtTruePtRecDifOverPtTrueMassCut->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
        outputContainer->Add(fhMCEtaPtTruePtRecDifOverPtTrueMassCut) ;
      }
      
      if ( fFillAngleHisto )
      {
        fhMCPi0PtRecOpenAngle = new TH2F
        ("hMCPi0PtRecOpenAngle",
         "Opening angle of true #pi^{0} cluster pairs",
         nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax) ;
        fhMCPi0PtRecOpenAngle->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCPi0PtRecOpenAngle->SetYTitle("#theta(rad)");
        outputContainer->Add(fhMCPi0PtRecOpenAngle) ;
        
        fhMCPi0PtRecOpenAngleMassCut = new TH2F
        ("hMCPi0PtRecOpenAngleCutMassCut",
         Form("Opening angle of true #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
              fPi0MassWindow[0],fPi0MassWindow[1]),
         nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax) ;
        fhMCPi0PtRecOpenAngleMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCPi0PtRecOpenAngleMassCut->SetYTitle("#theta(rad)");
        outputContainer->Add(fhMCPi0PtRecOpenAngleMassCut) ;
        
        fhMCEtaPtRecOpenAngle = new TH2F
        ("hMCEtaPtRecOpenAngle",
         "Opening angle of true #eta cluster pairs",
         nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax) ;
        fhMCEtaPtRecOpenAngle->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCEtaPtRecOpenAngle->SetYTitle("#theta(rad)");
        outputContainer->Add(fhMCEtaPtRecOpenAngle) ;
        
        fhMCEtaPtRecOpenAngleMassCut = new TH2F
        ("hMCEtaPtRecOpenAngleCutMassCut",
         Form("Opening angle of true #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
              fEtaMassWindow[0],fEtaMassWindow[1]),
         nptbins,ptmin,ptmax,nopanbins,opanmin,opanmax) ;
        fhMCEtaPtRecOpenAngleMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
        fhMCEtaPtRecOpenAngleMassCut->SetYTitle("#theta(rad)");
        outputContainer->Add(fhMCEtaPtRecOpenAngleMassCut) ;
      }
      
    }
    else if ( fFillOriginHisto && IsHighMultiplicityAnalysisOn() )
    {
      fhMCPi0MassPtTrueCen = new TH3F
      ("hMCPi0MassPtTrueCen",
       "Reconstructed Mass vs generated #it{p}_{T} of true #pi^{0} cluster pairs",
       ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
       massBinsArray.GetSize() - 1, massBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
      fhMCPi0MassPtTrueCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCPi0MassPtTrueCen->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      fhMCPi0MassPtTrueCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0MassPtTrueCen) ;
      
      fhMCPi0MassPtRecCen = new TH3F
      ("hMCPi0MassPtRecCen",
       "Reconstructed Mass vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs",
       ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
       massBinsArray.GetSize() - 1, massBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
      fhMCPi0MassPtRecCen->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCPi0MassPtRecCen->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      fhMCPi0MassPtRecCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0MassPtRecCen) ;
      
      fhMCPi0PtTruePtRecCen= new TH3F
      ("hMCPi0PtTruePtRecCen",
       "Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs",
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCPi0PtTruePtRecCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCPi0PtTruePtRecCen->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
      fhMCPi0PtTruePtRecCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0PtTruePtRecCen) ;
      
      fhMCPi0PtTruePtRecMassCutCen= new TH3F
      ("hMCPi0PtTruePtRecCenMassCut",
       Form("Generated vs reconstructed #it{p}_{T} of true #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}",fPi0MassWindow[0],fPi0MassWindow[1]),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCPi0PtTruePtRecMassCutCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCPi0PtTruePtRecMassCutCen->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
      fhMCPi0PtTruePtRecMassCutCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0PtTruePtRecMassCutCen) ;
      
      fhMCEtaMassPtTrueCen = new TH3F
      ("hMCEtaMassPtTrueCen",
       "Reconstructed Mass vs generated #it{p}_{T} of true #eta cluster pairs",
       ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
       massBinsArray.GetSize() - 1, massBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
      fhMCEtaMassPtTrueCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCEtaMassPtTrueCen->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      fhMCEtaMassPtTrueCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaMassPtTrueCen) ;
      
      fhMCEtaMassPtRecCen = new TH3F
      ("hMCEtaMassPtRecCen",
       "Reconstructed Mass vs reconstructed #it{p}_{T} of true #eta cluster pairs",
       ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
       massBinsArray.GetSize() - 1, massBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
      fhMCEtaMassPtRecCen->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCEtaMassPtRecCen->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
      fhMCEtaMassPtRecCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaMassPtRecCen) ;

      fhMCEtaPtTruePtRecCen= new TH3F
      ("hMCEtaPtTruePtRecCen",
       "Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs",
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCEtaPtTruePtRecCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCEtaPtTruePtRecCen->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
      fhMCEtaPtTruePtRecCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaPtTruePtRecCen) ;       
            
      fhMCEtaPtTruePtRecMassCutCen= new TH3F
      ("hMCEtaPtTruePtRecCenMassCut",
       Form("Generated vs reconstructed #it{p}_{T} of true #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}",
            fEtaMassWindow[0],fEtaMassWindow[1]), 
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCEtaPtTruePtRecMassCutCen->SetXTitle("#it{p}_{T, generated} (GeV/#it{c})");
      fhMCEtaPtTruePtRecMassCutCen->SetYTitle("#it{p}_{T, reconstructed} (GeV/#it{c})");
      fhMCEtaPtTruePtRecMassCutCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaPtTruePtRecMassCutCen) ;
      
      fhMCPi0PtTruePtRecDifOverPtTrueCen = new TH3F
      ("hMCPi0PtTruePtRecDifOverPtTrueCen",
       "(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} of true #pi^{0} cluster pairs",
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       difBinsArray.GetSize() - 1, difBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCPi0PtTruePtRecDifOverPtTrueCen->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCPi0PtTruePtRecDifOverPtTrueCen->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
      fhMCPi0PtTruePtRecDifOverPtTrueCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0PtTruePtRecDifOverPtTrueCen) ;
      
      fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut = new TH3F
      ("hMCPi0PtTruePtRecDifOverPtTrueCenMassCut",
       Form("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} #pi^{0} cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
            fPi0MassWindow[0],fPi0MassWindow[1]),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       difBinsArray.GetSize() - 1, difBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
      fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut) ;
      
      fhMCEtaPtTruePtRecDifOverPtTrueCen = new TH3F
      ("hMCEtaPtTruePtRecDifOverPtTrueCen",
       "(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} of true #eta cluster pairs",
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       difBinsArray.GetSize() - 1, difBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCEtaPtTruePtRecDifOverPtTrueCen->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCEtaPtTruePtRecDifOverPtTrueCen->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} ");
      fhMCEtaPtTruePtRecDifOverPtTrueCen->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaPtTruePtRecDifOverPtTrueCen) ;
      
      fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut = new TH3F
      ("hMCEtaPtTruePtRecDifCutOverPtTrueCenMassCut",
       Form("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener} #eta cluster pairs, %2.2f < rec. mass < %2.2f MeV/#it{c}^{2}", 
            fEtaMassWindow[0],fEtaMassWindow[1]),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       difBinsArray.GetSize() - 1, difBinsArray.GetArray(),
       cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
      fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
      fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut->SetYTitle("(#it{p}_{T, reco} - #it{p}_{T, gener}) / #it{p}_{T, gener}");
      fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut->SetZTitle("Centrality (%)");
      outputContainer->Add(fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut) ;
    } // High Mult, Origin

  } // MC data
  
  if ( fFillSMCombinations )
  {  
    AliDebug(1,Form("*** NMod = %d first %d last %d\n",fNModules, fFirstModule, fLastModule));
    if(fLastModule >= fNModules)
      AliError(Form("Last module number <%d> is larger than total SM number <%d>, please check configuration \n",fLastModule,fNModules));
    
    if(!fPairWithOtherDetector)
    {
      // Init the number of modules, set in the class AliCalorimeterUtils  
      fhReMod  = new TH2F*[fNModules] ;
      
      if(GetCalorimeter() == kPHOS)
      {
        fhReDiffPHOSMod        = new TH2F*[fNModules]   ;
      }
      else
      {
        fhReSameSectorEMCALMod = new TH2F*[fNModules/2] ;
        fhReSameSideEMCALMod   = new TH2F*[fNModules-2] ;
      }
      
      if(DoOwnMix())
      {
        fhMiMod  = new TH2F*[fNModules] ;
        if(GetCalorimeter() == kPHOS)
        {
          fhMiDiffPHOSMod        = new TH2F*[fNModules]   ;
        }
        else
        {
          fhMiSameSectorEMCALMod = new TH2F*[fNModules/2] ;
          fhMiSameSideEMCALMod   = new TH2F*[fNModules-2] ;
        }
      }
      
      // Single super modules
      //      
      for(Int_t imod=0; imod<fNModules; imod++)
      {
        if ( imod < fFirstModule || imod > fLastModule )
        {
          fhReMod[imod] = 0x0;
          if ( DoOwnMix() )
            fhMiMod[imod] = 0x0;
          continue;
        }
        // Module dependent invariant mass
        snprintf(key, buffersize,"hReMod_%d",imod) ;
        snprintf(title, buffersize,"Real #it{M}_{#gamma#gamma} distr. for Module %d",imod) ;
        fhReMod[imod]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhReMod[imod]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhReMod[imod]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhReMod[imod]) ;
        if(DoOwnMix())
        {
          snprintf(key, buffersize,"hMiMod_%d",imod) ;
          snprintf(title, buffersize,"Mixed #it{M}_{#gamma#gamma} distr. for Module %d",imod) ;
          fhMiMod[imod]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiMod[imod]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMiMod[imod]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhMiMod[imod]) ;
        }
      }
      
      // Super modules combinations
      //
      if(GetCalorimeter()==kPHOS)
      { 
        TString pairnamePHOS[] = {"(0-1)","(0-2)","(1-2)","(0-3)","(1-3)","(2-3)","(0-4)","(1-4)","(2-4)","(3-4)"};
        
        for(Int_t imod=0; imod<10; imod++)
        {
          if ( (fNModules == 3 && imod > 2) || (fNModules == 4 && imod > 5) )
          {
            fhReDiffPHOSMod[imod]  = 0x0;
            if ( DoOwnMix() )
              fhMiDiffPHOSMod[imod] = 0x0;
            continue;
          }
          
          snprintf(key, buffersize,"hReDiffPHOSMod_%d",imod) ;
          snprintf(title, buffersize,"Real pairs PHOS, clusters in different Modules: %s",(pairnamePHOS[imod]).Data()) ;
          fhReDiffPHOSMod[imod]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhReDiffPHOSMod[imod]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhReDiffPHOSMod[imod]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhReDiffPHOSMod[imod]) ;
          
          if(DoOwnMix())
          {          
            snprintf(key, buffersize,"hMiDiffPHOSMod_%d",imod) ;
            snprintf(title, buffersize,"Mixed pairs PHOS, clusters in different Modules: %s",(pairnamePHOS[imod]).Data()) ;
            fhMiDiffPHOSMod[imod]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMiDiffPHOSMod[imod]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMiDiffPHOSMod[imod]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMiDiffPHOSMod[imod]) ;
          }
        }
      }
      else
      {
        // EMCAL
        
        // Sectors
        //
        Int_t maxSector = (Int_t) fLastModule /2 ;
        Int_t minSector = (Int_t) fFirstModule/2 ;
        for(Int_t isector=minSector; isector<=maxSector; isector++)
        {
          snprintf(key, buffersize,"hReSameSectorEMCALMod_%d",isector) ;
          snprintf(title, buffersize,"Real pairs EMCAL, clusters in same sector, SM(%d,%d)",isector*2,isector*2+1) ;
          fhReSameSectorEMCALMod[isector]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhReSameSectorEMCALMod[isector]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhReSameSectorEMCALMod[isector]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhReSameSectorEMCALMod[isector]) ;
          
          if(DoOwnMix())
          {         
            snprintf(key, buffersize,"hMiSameSectorEMCALMod_%d",isector) ;
            snprintf(title, buffersize,"Mixed pairs EMCAL, clusters in same sector, SM(%d,%d)",isector*2,isector*2+1) ;
            fhMiSameSectorEMCALMod[isector]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMiSameSectorEMCALMod[isector]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMiSameSectorEMCALMod[isector]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMiSameSectorEMCALMod[isector]) ;
          } 
        }// sectors
        
        // Sides 
        //
        Int_t minSide = 0;
        Int_t maxSide = fNModules-2;
        
        if(fLastModule > 11)
          maxSide = fNModules-4;
        else
          maxSide = fLastModule-1;
        
        if(fFirstModule > 11)
          minSide = 10;
        //printf("** Last %d, First %d, min %d, max %d\n",fLastModule,fFirstModule,minSide,maxSide);
        for(Int_t iside=minSide; iside<maxSide; iside++)
        {
          Int_t ism1 = iside;
          Int_t ism2 = iside+2;
          if(iside > 9) // skip EMCal-DCal combination
          {
            ism1 = iside+2;
            ism2 = iside+4;
          }
          
          //printf("iside %d, sm1 %d, sm2 %d\n",iside,ism1,ism2);
          
          snprintf(key, buffersize,"hReSameSideEMCALMod_%d",iside) ;
          snprintf(title, buffersize,"Real pairs EMCAL, clusters in same side SM(%d,%d)",ism1, ism2) ;
          fhReSameSideEMCALMod[iside]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhReSameSideEMCALMod[iside]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhReSameSideEMCALMod[iside]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhReSameSideEMCALMod[iside]) ;
          
          if(DoOwnMix())
          {         
            snprintf(key, buffersize,"hMiSameSideEMCALMod_%d",iside) ;
            snprintf(title, buffersize,"Mixed pairs EMCAL, clusters in same side SM(%d,%d)",ism1, ism2) ;
            fhMiSameSideEMCALMod[iside]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMiSameSideEMCALMod[iside]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhMiSameSideEMCALMod[iside]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
            outputContainer->Add(fhMiSameSideEMCALMod[iside]) ;
          } // mix
        } // sides
        
      }//EMCAL
    } // Not pair of detectors
    else
    {
      fhReSameSectorDCALPHOSMod = new TH2F*[6] ;
      fhReDiffSectorDCALPHOSMod = new TH2F*[8] ;
      fhMiSameSectorDCALPHOSMod = new TH2F*[6] ;
      fhMiDiffSectorDCALPHOSMod = new TH2F*[8] ;
      
      Int_t dcSameSM[6] = {12,13,14,15,16,17}; // Check eta order
      Int_t phSameSM[6] = {3,  3, 2, 2, 1, 1};
      
      Int_t dcDiffSM[8] = {12,13,14,15,16,17,0,0};
      Int_t phDiffSM[8] = {2,  2, 1, 1, 3, 3,0,0};
      
      for(Int_t icombi = 0; icombi < 8; icombi++)
      {
        snprintf(key, buffersize,"hReDiffSectorDCALPHOS_%d",icombi) ;
        snprintf(title, buffersize,"Real pairs DCAL-PHOS, clusters in different sector, SM(%d,%d)",dcDiffSM[icombi],phDiffSM[icombi]) ;
        fhReDiffSectorDCALPHOSMod[icombi]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhReDiffSectorDCALPHOSMod[icombi]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhReDiffSectorDCALPHOSMod[icombi]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhReDiffSectorDCALPHOSMod[icombi]) ;
        if(DoOwnMix())
        {
          snprintf(key, buffersize,"hMiDiffSectorDCALPHOS_%d",icombi) ;
          snprintf(title, buffersize,"Mixed pairs DCAL-PHOS, clusters in different sector, SM(%d,%d)",dcDiffSM[icombi],phDiffSM[icombi]) ;
          fhMiDiffSectorDCALPHOSMod[icombi]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiDiffSectorDCALPHOSMod[icombi]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMiDiffSectorDCALPHOSMod[icombi]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhMiDiffSectorDCALPHOSMod[icombi]) ;
        }
        
        if ( icombi > 5 ) continue ;
        
        snprintf(key, buffersize,"hReSameSectorDCALPHOS_%d",icombi) ;
        snprintf(title, buffersize,"Real pairs DCAL-PHOS, clusters in same sector, SM(%d,%d)",dcSameSM[icombi],phSameSM[icombi]) ;
        fhReSameSectorDCALPHOSMod[icombi]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhReSameSectorDCALPHOSMod[icombi]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhReSameSectorDCALPHOSMod[icombi]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
        outputContainer->Add(fhReSameSectorDCALPHOSMod[icombi]) ;
        if(DoOwnMix())
        {
          snprintf(key, buffersize,"hMiSameSectorDCALPHOS_%d",icombi) ;
          snprintf(title, buffersize,"Mixed pairs DCAL-PHOS, clusters in same sector, SM(%d,%d)",dcSameSM[icombi],phSameSM[icombi]) ;
          fhMiSameSectorDCALPHOSMod[icombi]  = new TH2F(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiSameSectorDCALPHOSMod[icombi]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhMiSameSectorDCALPHOSMod[icombi]->SetYTitle("#it{M}_{#gamma,#gamma} (GeV/#it{c}^{2})");
          outputContainer->Add(fhMiSameSectorDCALPHOSMod[icombi]) ;
        }
      }
      
    }
  } // SM combinations
  //
  if ( fFillOpAngleCutHisto )
  {
    for(Int_t icut = 0; icut < fNAngleCutBins; icut++)
    {
      fhReOpAngleBinMinClusterEtaPhi[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMin_EtaPhi",icut),
       Form("#eta vs #varphi, cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhReOpAngleBinMinClusterEtaPhi[icut]->SetYTitle("#varphi (rad)");
      fhReOpAngleBinMinClusterEtaPhi[icut]->SetXTitle("#eta");
      outputContainer->Add(fhReOpAngleBinMinClusterEtaPhi[icut]) ;
      
      fhReOpAngleBinMaxClusterEtaPhi[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMax_EtaPhi",icut),
       Form("#eta vs #varphi, cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhReOpAngleBinMaxClusterEtaPhi[icut]->SetYTitle("#varphi (rad)");
      fhReOpAngleBinMaxClusterEtaPhi[icut]->SetXTitle("#eta");
      outputContainer->Add(fhReOpAngleBinMaxClusterEtaPhi[icut]) ;
      
      fhReOpAngleBinMinClusterColRow[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMin_ColRow",icut),
       Form("highest #it{E} cell, column vs row, cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhReOpAngleBinMinClusterColRow[icut]->SetYTitle("row");
      fhReOpAngleBinMinClusterColRow[icut]->SetXTitle("column");
      outputContainer->Add(fhReOpAngleBinMinClusterColRow[icut]) ;
      
      fhReOpAngleBinMaxClusterColRow[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMax_ColRow",icut),
       Form("highest #it{E} cell, column vs row, cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhReOpAngleBinMaxClusterColRow[icut]->SetYTitle("row");
      fhReOpAngleBinMaxClusterColRow[icut]->SetXTitle("column");
      outputContainer->Add(fhReOpAngleBinMaxClusterColRow[icut]) ;
      
      fhReOpAngleBinMinClusterEPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMin_EPerSM",icut),
       Form("cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       nptbins,ptmin,ptmax,20,0,20);
      fhReOpAngleBinMinClusterEPerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMinClusterEPerSM[icut]->SetXTitle("#it{E} (GeV/#it{c})");
      outputContainer->Add(fhReOpAngleBinMinClusterEPerSM[icut]) ;
      
      fhReOpAngleBinMaxClusterEPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMax_EPerSM",icut),
       Form("cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       nptbins,ptmin,ptmax,20,0,20);
      fhReOpAngleBinMaxClusterEPerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMaxClusterEPerSM[icut]->SetXTitle("#it{E} (GeV/#it{c})");
      outputContainer->Add(fhReOpAngleBinMaxClusterEPerSM[icut]) ;
      
      fhReOpAngleBinMinClusterTimePerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMin_TimePerSM",icut),
       Form("cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       ntimebins,timemin,timemax,20,0,20);
      fhReOpAngleBinMinClusterTimePerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMinClusterTimePerSM[icut]->SetXTitle("#it{t} (ns)");
      outputContainer->Add(fhReOpAngleBinMinClusterTimePerSM[icut]) ;
      
      fhReOpAngleBinMaxClusterTimePerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMax_TimePerSM",icut),
       Form("cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       ntimebins,timemin,timemax,20,0,20);
      fhReOpAngleBinMaxClusterTimePerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMaxClusterTimePerSM[icut]->SetXTitle("#it{t} (ns)");
      outputContainer->Add(fhReOpAngleBinMaxClusterTimePerSM[icut]) ;
      
      fhReOpAngleBinMinClusterNCellPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMin_NCellPerSM",icut),
       Form("cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       30,0,30,20,0,20);
      fhReOpAngleBinMinClusterNCellPerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMinClusterNCellPerSM[icut]->SetXTitle("# cells");
      outputContainer->Add(fhReOpAngleBinMinClusterNCellPerSM[icut]) ;
      
      fhReOpAngleBinMaxClusterNCellPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_ClusterMax_NCellPerSM",icut),
       Form("cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       30,0,30,20,0,20);
      fhReOpAngleBinMaxClusterNCellPerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinMaxClusterNCellPerSM[icut]->SetXTitle("# cells");
      outputContainer->Add(fhReOpAngleBinMaxClusterNCellPerSM[icut]) ;
      
      fhReOpAngleBinPairClusterRatioPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_PairCluster_RatioPerSM",icut),
       Form("cluster pair #it{E}_{high}/ #it{E}_{low}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       100,0,1,20,0,20);
      fhReOpAngleBinPairClusterRatioPerSM[icut]->SetYTitle("SM");
      fhReOpAngleBinPairClusterRatioPerSM[icut]->SetXTitle("#it{E}_{low}/ #it{E}_{high}");
      outputContainer->Add(fhReOpAngleBinPairClusterRatioPerSM[icut]) ;
      
      fhReOpAngleBinPairClusterMassPerSM[icut] = new TH2F
      (Form("hReOpAngleBin%d_PairCluster_MassPerSM",icut),
       Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       nmassbins,massmin,massmax,20,0,20);
      fhReOpAngleBinPairClusterMassPerSM[icut]->SetXTitle("#it{M} (GeV/#it{c}^2)");
      fhReOpAngleBinPairClusterMassPerSM[icut]->SetYTitle("SM");
      outputContainer->Add(fhReOpAngleBinPairClusterMassPerSM[icut]) ;
      
      fhReOpAngleBinPairClusterMass[icut] = new TH2F
      (Form("hReOpAngleBin%d_PairCluster_Mass",icut),
       Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
       nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
      fhReOpAngleBinPairClusterMass[icut]->SetYTitle("#it{M} (GeV/#it{c}^2)");
      fhReOpAngleBinPairClusterMass[icut]->SetXTitle("#it{p}_{T} GeV/#it{c}");
      outputContainer->Add(fhReOpAngleBinPairClusterMass[icut]) ;
      
      if ( IsDataMC() )
      {
        if ( fFillOriginHisto )
        {
          fhReOpAngleBinPairClusterMassMCTruePi0[icut] = new TH2F
          (Form("hReOpAngleBin%d_PairCluster_MassMCTruePi0",icut),
           Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f, true mc #pi^{0}",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
           nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
          fhReOpAngleBinPairClusterMassMCTruePi0[icut]->SetYTitle("#it{M} (GeV/#it{c}^2)");
          fhReOpAngleBinPairClusterMassMCTruePi0[icut]->SetXTitle("#it{p}_{T} GeV/#it{c}");
          outputContainer->Add(fhReOpAngleBinPairClusterMassMCTruePi0[icut]) ;
          
          fhReOpAngleBinPairClusterMassMCTrueEta[icut] = new TH2F
          (Form("hReOpAngleBin%d_PairCluster_MassMCTrueEta",icut),
           Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f, true mc #eta",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
           nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
          fhReOpAngleBinPairClusterMassMCTrueEta[icut]->SetYTitle("#it{M} (GeV/#it{c}^2)");
          fhReOpAngleBinPairClusterMassMCTrueEta[icut]->SetXTitle("#it{p}_{T} GeV/#it{c}");
          outputContainer->Add(fhReOpAngleBinPairClusterMassMCTrueEta[icut]) ;
        }
        
        fhPrimPi0AccPtOpAngCuts[icut]  = new TH1F
        (Form("hPrimPi0AccPt_OpAngleBin%d",icut),
         Form("Primary #pi^{0} #it{p}_{T} with both photons in acceptance, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nptbins,ptmin,ptmax) ;
        fhPrimPi0AccPtOpAngCuts[icut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimPi0AccPtOpAngCuts[icut]) ;
        
        fhPrimEtaAccPtOpAngCuts[icut]  = new TH1F
        (Form("hPrimEtaAccPt_OpAngleBin%d",icut),
         Form("Primary #eta #it{p}_{T} with both photons in acceptance, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nptbins,ptmin,ptmax) ;
        fhPrimEtaAccPtOpAngCuts[icut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimEtaAccPtOpAngCuts[icut]) ;
      }
      
      //       fhReOpAngleBinPairClusterAbsIdMaxCell[icut] = new TH2F
      //      (Form("hReOpAngleBin%d_PairCluster_AbsIdCell",icut),
      //       Form("cluster pair Abs Cell ID low #it{E} vs high #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
      //       //17664,0,17664,17664,0,17664);
      //       1689,0,16896,1689,0,16896);
      //      fhReOpAngleBinPairClusterAbsIdMaxCell[icut]->SetYTitle("AbsId-higher");
      //      fhReOpAngleBinPairClusterAbsIdMaxCell[icut]->SetXTitle("AbsId-lower");
      //      outputContainer->Add(fhReOpAngleBinPairClusterAbsIdMaxCell[icut]) ;
      
      if( DoOwnMix() )
      {
        fhMiOpAngleBinMinClusterEtaPhi[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMin_EtaPhi",icut),
         Form("#eta vs #varphi, mixed cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhMiOpAngleBinMinClusterEtaPhi[icut]->SetYTitle("#varphi (rad)");
        fhMiOpAngleBinMinClusterEtaPhi[icut]->SetXTitle("#eta");
        outputContainer->Add(fhMiOpAngleBinMinClusterEtaPhi[icut]) ;
        
        fhMiOpAngleBinMaxClusterEtaPhi[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMax_EtaPhi",icut),
         Form("#eta vs #varphi, mixed cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhMiOpAngleBinMaxClusterEtaPhi[icut]->SetYTitle("#varphi (rad)");
        fhMiOpAngleBinMaxClusterEtaPhi[icut]->SetXTitle("#eta");
        outputContainer->Add(fhMiOpAngleBinMaxClusterEtaPhi[icut]) ;
        
        //        fhMiOpAngleBinMinClusterColRow[icut] = new TH2F
        //        (Form("hMiOpAngleBin%d_ClusterMin_ColRow",icut),
        //         Form("highest #it{E} cell, column vs row, mixed cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
        //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        //        fhMiOpAngleBinMinClusterColRow[icut]->SetYTitle("row");
        //        fhMiOpAngleBinMinClusterColRow[icut]->SetXTitle("column");
        //        outputContainer->Add(fhMiOpAngleBinMinClusterColRow[icut]) ;
        //        
        //        fhMiOpAngleBinMaxClusterColRow[icut] = new TH2F
        //        (Form("hMiOpAngleBin%d_ClusterMax_ColRow",icut),
        //         Form("highest #it{E} cell, column vs row, mixed cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
        //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        //        fhMiOpAngleBinMaxClusterColRow[icut]->SetYTitle("row");
        //        fhMiOpAngleBinMaxClusterColRow[icut]->SetXTitle("column");
        //        outputContainer->Add(fhMiOpAngleBinMaxClusterColRow[icut]) ;
        
        fhMiOpAngleBinMinClusterEPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMin_EPerSM",icut),
         Form("mixed cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nptbins,ptmin,ptmax,20,0,20);
        fhMiOpAngleBinMinClusterEPerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMinClusterEPerSM[icut]->SetXTitle("#it{E} (GeV/#it{c})");
        outputContainer->Add(fhMiOpAngleBinMinClusterEPerSM[icut]) ;
        
        fhMiOpAngleBinMaxClusterEPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMax_EPerSM",icut),
         Form("mixed cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nptbins,ptmin,ptmax,20,0,20);
        fhMiOpAngleBinMaxClusterEPerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMaxClusterEPerSM[icut]->SetXTitle("#it{E} (GeV/#it{c})");
        outputContainer->Add(fhMiOpAngleBinMaxClusterEPerSM[icut]) ;
        
        fhMiOpAngleBinMinClusterTimePerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMin_TimePerSM",icut),
         Form("mixed cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         ntimebins,timemin,timemax,20,0,20);
        fhMiOpAngleBinMinClusterTimePerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMinClusterTimePerSM[icut]->SetXTitle("#it{t} (ns)");
        outputContainer->Add(fhMiOpAngleBinMinClusterTimePerSM[icut]) ;
        
        fhMiOpAngleBinMaxClusterTimePerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMax_TimePerSM",icut),
         Form("mixed cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         ntimebins,timemin,timemax,20,0,20);
        fhMiOpAngleBinMaxClusterTimePerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMaxClusterTimePerSM[icut]->SetXTitle("#it{t} (ns)");
        outputContainer->Add(fhMiOpAngleBinMaxClusterTimePerSM[icut]) ;
        
        fhMiOpAngleBinMinClusterNCellPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMin_NCellPerSM",icut),
         Form("mixed cluster pair lower #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         30,0,30,20,0,20);
        fhMiOpAngleBinMinClusterNCellPerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMinClusterNCellPerSM[icut]->SetXTitle("# cells");
        outputContainer->Add(fhMiOpAngleBinMinClusterNCellPerSM[icut]) ;
        
        fhMiOpAngleBinMaxClusterNCellPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_ClusterMax_NCellPerSM",icut),
         Form("mixed cluster pair higher #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         30,0,30,20,0,20);
        fhMiOpAngleBinMaxClusterNCellPerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinMaxClusterNCellPerSM[icut]->SetXTitle("# cells");
        outputContainer->Add(fhMiOpAngleBinMaxClusterNCellPerSM[icut]) ;
        
        fhMiOpAngleBinPairClusterRatioPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_PairCluster_RatioPerSM",icut),
         Form("mixed cluster pair #it{E}_{high}/ #it{E}_{low}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         100,0,1,20,0,20);
        fhMiOpAngleBinPairClusterRatioPerSM[icut]->SetYTitle("SM");
        fhMiOpAngleBinPairClusterRatioPerSM[icut]->SetXTitle("#it{E}_{low}/ #it{E}_{high}");
        outputContainer->Add(fhMiOpAngleBinPairClusterRatioPerSM[icut]) ;
        
        fhMiOpAngleBinPairClusterMassPerSM[icut] = new TH2F
        (Form("hMiOpAngleBin%d_PairCluster_MassPerSM",icut),
         Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nmassbins,massmin,massmax,20,0,20);
        fhMiOpAngleBinPairClusterMassPerSM[icut]->SetXTitle("#it{M} (GeV/#it{c}^2)");
        fhMiOpAngleBinPairClusterMassPerSM[icut]->SetYTitle("SM");
        outputContainer->Add(fhMiOpAngleBinPairClusterMassPerSM[icut]) ;
        
        fhMiOpAngleBinPairClusterMass[icut] = new TH2F
        (Form("hMiOpAngleBin%d_PairCluster_Mass",icut),
         Form("cluster pair #it{M}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMiOpAngleBinPairClusterMass[icut]->SetYTitle("#it{M} (GeV/#it{c}^2)");
        fhMiOpAngleBinPairClusterMass[icut]->SetXTitle("#it{p}_{T} GeV/#it{c}");
        outputContainer->Add(fhMiOpAngleBinPairClusterMass[icut]) ;
        
        //        fhMiOpAngleBinPairClusterAbsIdMaxCell[icut] = new TH2F
        //        (Form("hMiOpAngleBin%d_PairCluster_AbsIdCell",icut),
        //         Form("mixed cluster pair Abs Cell ID low #it{E} vs high #it{E}, pair %1.4f<#theta<%1.4f",fAngleCutBinsArray[icut],fAngleCutBinsArray[icut+1]),
        //         ntimebins,timemin,timemax,20,0,20);
        //        fhMiOpAngleBinPairClusterRatioPerSM[icut]->SetYTitle("AbsId-higher");
        //        fhMiOpAngleBinPairClusterRatioPerSM[icut]->SetXTitle("AbsId-lower");
        //        outputContainer->Add(fhMiOpAngleBinPairClusterRatioPerSM[icut]) ;
      }
    } // angle bin
  }  
  
  if ( IsDataMC() && IsStudyClusterOverlapsPerGeneratorOn() )
  {
    TString tagBkgNames[] = { 
      "Clean", "HijingBkg", "NotHijingBkg", "HijingAndOtherBkg",
      "Clean_HijingBkg", "Clean_NotHijingBkg", "Clean_HijingAndOtherBkg",
      "HijingBkg_NotHijingBkg", "HijingBkg_HijingAndOtherBkg", "NotHijingBkg_HijingAndOtherBkg" } ;
    TString tagBkgTitle[] = {
      "no overlap", "pair Hijing Bkg", "pair not Hijing bkg", "pair Hijing and other bkg",
      "no overlap and hijing overlap", "no overlap and generator overlap", "no overlap and multiple overlap",
      "hijing overlap and gener overlap", "hijing overlap and multiple overlap", "gener overlap and multiple overlap" } ;
    
    for(Int_t igen = 0; igen < GetNCocktailGenNamesToCheck(); igen++)
    {
      TString add = "_MainGener_";
      if(igen==0)
        add = "";
      else
      {
        fhPrimPi0PtPerGenerator[igen-1]     = new TH1F
        (Form("hPrimPi0Pt%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Primary #pi^{0} #it{p}_{T}, |#it{Y}| < 1, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax) ;
        fhPrimPi0PtPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimPi0PtPerGenerator[igen-1]) ;
        
        fhPrimPi0YPerGenerator[igen-1]      = new TH2F
        (Form("hPrimPi0Rapidity%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Rapidity of primary #pi^{0}, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
        fhPrimPi0YPerGenerator[igen-1]   ->SetYTitle("#it{Rapidity}");
        fhPrimPi0YPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimPi0YPerGenerator[igen-1]) ;
        
        fhPrimPi0PhiPerGenerator[igen-1]    = new TH2F
        (Form("hPrimPi0Phi%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#varphi of primary #pi^{0}, |#it{Y}|<1, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax,nphibinsopen,0,360) ;
        fhPrimPi0PhiPerGenerator[igen-1]->SetYTitle("#varphi (deg)");
        fhPrimPi0PhiPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimPi0PhiPerGenerator[igen-1]) ;
        
        if ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() )
        {
          fhPrimPi0PtInCaloPerGenerator[igen-1]     = new TH1F
          (Form("hPrimPi0PtInCalo%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #pi^{0} #it{p}_{T}, in calorimeter acc. %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimPi0PtInCaloPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimPi0PtInCaloPerGenerator[igen-1]) ;
          
          fhPrimPi0AccPtPerGenerator[igen-1]  = new TH1F
          (Form("hPrimPi0AccPt%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #pi^{0} #it{p}_{T} with both photons in acceptance, generator %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimPi0AccPtPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimPi0AccPtPerGenerator[igen-1]) ;
          
          fhPrimPi0AccPtPhotonCutsPerGenerator[igen-1]  = new TH1F
          (Form("hPrimPi0AccPtPhotonCuts%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #pi^{0} #it{p}_{T} with both photons in acceptance, generator %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimPi0AccPtPhotonCutsPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimPi0AccPtPhotonCutsPerGenerator[igen-1]) ;
        }
        
        //
        
        fhPrimEtaPtPerGenerator[igen-1]     = new TH1F
        (Form("hPrimEtaPt%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Primary #eta #it{p}_{T}, |#it{Y}| < 1, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax) ;
        fhPrimEtaPtPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimEtaPtPerGenerator[igen-1]) ;
        
        Int_t netabinsopen =  TMath::Nint(netabins*4/(etamax-etamin));
        fhPrimEtaYPerGenerator[igen-1]      = new TH2F
        (Form("hPrimEtaRapidity%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Rapidity of primary #eta, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax,netabinsopen,-2, 2) ;
        fhPrimEtaYPerGenerator[igen-1]   ->SetYTitle("#it{Rapidity}");
        fhPrimEtaYPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimEtaYPerGenerator[igen-1]) ;
        
        fhPrimEtaPhiPerGenerator[igen-1]    = new TH2F
        (Form("hPrimEtaPhi%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#varphi of primary #eta, |#it{Y}|<1, generator %s",GetCocktailGenNameToCheck(igen).Data()),
         nptbins,ptmin,ptmax,nphibinsopen,0,360) ;
        fhPrimEtaPhiPerGenerator[igen-1]->SetYTitle("#varphi (deg)");
        fhPrimEtaPhiPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPrimEtaPhiPerGenerator[igen-1]) ;
        
        if ( IsRealCaloAcceptanceOn() || IsFiducialCutOn() )
        {
          fhPrimEtaPtInCaloPerGenerator[igen-1]     = new TH1F
          (Form("hPrimEtaPtInCalo%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #eta #it{p}_{T}, in calorimeter acceptance %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimEtaPtInCaloPerGenerator[igen-1]   ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimEtaPtInCaloPerGenerator[igen-1]) ;
          
          fhPrimEtaAccPtPerGenerator[igen-1]  = new TH1F
          (Form("hPrimEtaAccPt%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #eta #it{p}_{T} with both photons in acceptance, generator %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimEtaAccPtPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimEtaAccPtPerGenerator[igen-1]) ;
          
          fhPrimEtaAccPtPhotonCutsPerGenerator[igen-1]  = new TH1F
          (Form("hPrimEtaAccPtPhotonCuts%s%s",add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Primary #eta #it{p}_{T} with both photons in acceptance, generator %s",GetCocktailGenNameToCheck(igen).Data()),
           nptbins,ptmin,ptmax) ;
          fhPrimEtaAccPtPhotonCutsPerGenerator[igen-1]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPrimEtaAccPtPhotonCutsPerGenerator[igen-1]) ;
        }
        
      }
      
      for(Int_t itag = 0; itag < 10; itag++)
      {
        fhPairGeneratorsBkgMass[igen][itag] = new TH2F
        (Form("h%sGeneratorPairMass%s%s",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Pair Mass with generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhPairGeneratorsBkgMass[igen][itag]->SetYTitle("#it{M} (MeV/#it{c}^{2})");
        fhPairGeneratorsBkgMass[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgMass[igen][itag]) ;
        
        fhPairGeneratorsBkgMassMCPi0[igen][itag] = new TH2F
        (Form("h%sGeneratorPairMass%s%s_MCPi0",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Pair Mass with contribution of true #pi^{0} generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhPairGeneratorsBkgMassMCPi0[igen][itag]->SetYTitle("#it{M} (MeV/#it{c}^{2})");
        fhPairGeneratorsBkgMassMCPi0[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgMassMCPi0[igen][itag]) ;
        
        fhPairGeneratorsBkgMassMCEta[igen][itag] = new TH2F
        (Form("h%sGeneratorPairMass%s%s_MCEta",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("Pair Mass with contribution of true #eta generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhPairGeneratorsBkgMassMCEta[igen][itag]->SetYTitle("#it{M} (MeV/#it{c}^{2})");
        fhPairGeneratorsBkgMassMCEta[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgMassMCEta[igen][itag]) ;
        
        if ( IsHighMultiplicityAnalysisOn() )
        {
          fhPairGeneratorsBkgCentMCPi0[igen][itag] = new TH2F
          (Form("h%sGeneratorPairCent%s%s_MCPi0",
                tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Pair Mass with contribution of true #pi^{0} generator%s, %s ",
                GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
           nptbins,ptmin,ptmax,100,0,100);
          fhPairGeneratorsBkgCentMCPi0[igen][itag]->SetYTitle("Centrality");
          fhPairGeneratorsBkgCentMCPi0[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
          outputContainer->Add(fhPairGeneratorsBkgCentMCPi0[igen][itag]) ;
          
          fhPairGeneratorsBkgCentMCPi0MassCut[igen][itag] = new TH2F
          (Form("h%sGeneratorPairCent%s%s_MCPi0_MassCut",
                tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Pair Mass with contribution of true #pi^{0} generator%s, %s, %2.2f<M<%2.2f",
                GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data(),fPi0MassWindow[0],fPi0MassWindow[1]),
           nptbins,ptmin,ptmax,100,0,100);
          fhPairGeneratorsBkgCentMCPi0MassCut[igen][itag]->SetYTitle("Centrality");
          fhPairGeneratorsBkgCentMCPi0MassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
          outputContainer->Add(fhPairGeneratorsBkgCentMCPi0MassCut[igen][itag]) ;
          
          fhPairGeneratorsBkgCentMCEta[igen][itag] = new TH2F
          (Form("h%sGeneratorPairCent%s%s_MCEta",
                tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Pair Mass with contribution of true #eta generator%s, %s ",
                GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
           nptbins,ptmin,ptmax,100,0,100);
          fhPairGeneratorsBkgCentMCEta[igen][itag]->SetYTitle("Centrality");
          fhPairGeneratorsBkgCentMCEta[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
          outputContainer->Add(fhPairGeneratorsBkgCentMCEta[igen][itag]) ;
          
          fhPairGeneratorsBkgCentMCEtaMassCut[igen][itag] = new TH2F
          (Form("h%sGeneratorPairCent%s%s_MCEta_MassCut",
                tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
           Form("Pair Mass with contribution of true #eta generator%s, %s, %2.2f<M<%2.2f",
                GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data(),fEtaMassWindow[0],fEtaMassWindow[1]),
           nptbins,ptmin,ptmax,100,0,100);
          fhPairGeneratorsBkgCentMCEtaMassCut[igen][itag]->SetYTitle("Centrality");
          fhPairGeneratorsBkgCentMCEtaMassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
          outputContainer->Add(fhPairGeneratorsBkgCentMCEtaMassCut[igen][itag]) ;
        }
        
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoRatio%s%s_MCPi0",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}/#it{E}_{gen} pair with contribution of true #pi^{0} generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nratbins,ratmin,ratmax);
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0[igen][itag]->SetYTitle("#it{E}_{reco}/#it{E}_{gen}");
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoRatioMCPi0[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoRatioMCEta[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoRatio%s%s_MCEta",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}/#it{E}_{gen} pair with contribution of true #eta generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nratbins,ratmin,ratmax);
        fhPairGeneratorsBkgEPrimRecoRatioMCEta[igen][itag]->SetYTitle("#it{E}_{reco}/#it{E}_{gen}");
        fhPairGeneratorsBkgEPrimRecoRatioMCEta[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoRatioMCEta[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoDiff%s%s_MCPi0",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}-#it{E}_{gen} pair with contribution of true #pi^{0} generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax);
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0[igen][itag]->SetYTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoDiffMCPi0[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoDiffMCEta[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoDiff%s%s_MCEta",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}-#it{E}_{gen} pair with contribution of true #eta generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax);
        fhPairGeneratorsBkgEPrimRecoDiffMCEta[igen][itag]->SetYTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
        fhPairGeneratorsBkgEPrimRecoDiffMCEta[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoDiffMCEta[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoRatio%s%s_MCPi0MassCut",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}/#it{E}_{gen} pair with contribution of true #pi^{0} generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nratbins,ratmin,ratmax);
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[igen][itag]->SetYTitle("#it{E}_{reco}/#it{E}_{gen}");
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoRatio%s%s_MCEtaMassCut",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}/#it{E}_{gen} pair with contribution of true #eta generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,nratbins,ratmin,ratmax);
        fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[igen][itag]->SetYTitle("#it{E}_{reco}/#it{E}_{gen}");
        fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoDiff%s%s_MCPi0MassCut",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}-#it{E}_{gen} pair with contribution of true #pi^{0} generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax);
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut[igen][itag]->SetYTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut[igen][itag]) ;
        
        fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut[igen][itag] = new TH2F
        (Form("h%sGeneratorPairEPrimRecoDiff%s%s_MCEtaMassCut",
              tagBkgNames[itag].Data(),add.Data(),GetCocktailGenNameToCheck(igen).Data()),
         Form("#it{E}_{reco}-#it{E}_{gen} pair with contribution of true #eta generator%s, %s ",
              GetCocktailGenNameToCheck(igen).Data(),tagBkgTitle[itag].Data()),
         nptbins,ptmin,ptmax,ndifbins,difmin,difmax);
        fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut[igen][itag]->SetYTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
        fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut[igen][itag]->SetXTitle("#it{E}_{reco} (GeV)");
        outputContainer->Add(fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut[igen][itag]) ;
      }
    }
  }
  
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //  {
  //    printf("Histogram %d, name: %s\n",i, outputContainer->At(i)->GetName());
  //  }
  return outputContainer;
}


//___________________________________________________
/// Print some relevant parameters set for the analysis.
//___________________________________________________
void AliAnaPi0::Print(const Option_t * /*opt*/) const
{
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Number of bins in Centrality:  %d \n",GetNCentrBin()) ;
  printf("Number of bins in Z vert. pos: %d \n",GetNZvertBin()) ;
  printf("Number of bins in Reac. Plain: %d \n",GetNRPBin()) ;
  printf("Depth of event buffer: %d \n",GetNMaxEvMix()) ;
  printf("Pair in same Module: %d \n",fSameSM) ;
  printf("Cuts: \n") ;
  // printf("Z vertex position: -%2.3f < z < %2.3f \n",GetZvertexCut(),GetZvertexCut()) ; //It crashes here, why?
  printf("Number of modules:             %d \n",fNModules) ;
  printf("Select pairs with their angle: %d, edep %d, min angle %2.3f, max angle %2.3f, 1cell %d \n",fUseAngleCut, fUseAngleEDepCut, fAngleCut, fAngleMaxCut, fUseOneCellSeparation) ;
  printf("Asymmetry cuts: n = %d, \n",fNAsymCuts) ;
  printf("\tasymmetry < ");
  for(Int_t i = 0; i < fNAsymCuts; i++) printf("%2.2f ",fAsymCuts[i]);
  printf("\n");
  
  printf("PID selection bits: n = %d, \n",fNPIDBits) ;
  printf("\tPID bit = ");
  for(Int_t i = 0; i < fNPIDBits; i++) printf("%d ",fPIDBits[i]);
  printf("\n");
  
  if(fMultiCutAna || fMultiCutAnaAcc)
  {
    printf("pT cuts: n = %d, \n",fNPtCuts) ;
    printf("\tpT > ");
    for(Int_t i = 0; i < fNPtCuts; i++) printf("%2.2f ",fPtCuts[i]);
    printf("GeV/c\n");
    printf("\tpT < ");
    for(Int_t i = 0; i < fNPtCuts; i++) printf("%2.2f ",fPtCutsMax[i]);
    printf("GeV/c\n");
    
    printf("N cell in cluster cuts: n = %d, \n",fNCellNCuts) ;
    printf("\tnCell > ");
    for(Int_t i = 0; i < fNCellNCuts; i++) printf("%d ",fCellNCuts[i]);
    printf("\n");
    
  }
  printf("------------------------------------------------------\n") ;
}

//________________________________________
/// Fill acceptance histograms if MC data is available.
//________________________________________
void AliAnaPi0::FillAcceptanceHistograms()
{
  if ( !GetMC() ) return;

  Double_t mesonY   = -100 ;
  Double_t mesonE   = -1 ;
  Double_t mesonPt  = -1 ;
  Double_t mesonPhi =  100 ;
  Double_t mesonYeta= -1 ;
  
  Int_t    pdg     = 0 ;
  Int_t    nprim   = GetMC()->GetNumberOfTracks();
  Int_t    nDaught = 0 ;
  Int_t    iphot1  = 0 ;
  Int_t    iphot2  = 0 ;
  
  Float_t cen = GetEventCentrality();
  Float_t ep  = GetEventPlaneAngle();
  
  AliVParticle * primary = 0;
  
  TString genName = "";
  
  for(Int_t i=0 ; i < nprim; i++)
  {
    if ( !GetReader()->AcceptParticleMCLabel( i ) ) continue ;
        
    primary = GetMC()->GetTrack(i) ;
    if(!primary)
    {
      AliWarning("Primaries pointer not available!!");
      continue;
    }
    
    // If too small  skip
    if( primary->E() < 0.4 ) continue;
    
    pdg       = primary->PdgCode();
    // Select only pi0 or eta
    if( pdg != 111 && pdg != 221 && pdg != 321 && pdg != 211) continue ;
    
    nDaught   = primary->GetNDaughters();
    iphot1    = primary->GetDaughterLabel(0) ;
    iphot2    = primary->GetDaughterLabel(1) ;
    
    // Protection against floating point exception
    if ( primary->E() == TMath::Abs(primary->Pz()) || 
        (primary->E() - primary->Pz()) < 1e-3      ||
        (primary->E() + primary->Pz()) < 0           )  continue ; 
    
    //printf("i %d, %s %d  %s %d \n",i, GetMC()->Particle(i)->GetName(), GetMC()->Particle(i)->GetPdgCode(),
    //       prim->GetName(), prim->GetPdgCode());
    
    //Photon kinematics
    primary->Momentum(fMCPrimMesonMom);
    
    mesonY = 0.5*TMath::Log((primary->E()+primary->Pz())/(primary->E()-primary->Pz())) ;
    
    Bool_t inacceptance = kTRUE;
    if ( IsFiducialCutOn() || IsRealCaloAcceptanceOn() ) 
    {
      // Check if pi0 enters the calo
      if(IsRealCaloAcceptanceOn() && !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance( GetCalorimeter(), primary )) inacceptance = kFALSE;
      if(IsFiducialCutOn() && inacceptance && !GetFiducialCut()->IsInFiducialCut(fMCPrimMesonMom.Eta(), fMCPrimMesonMom.Phi(), GetCalorimeter())) inacceptance = kFALSE;
    }
    else inacceptance = kFALSE;
    
    mesonPt  = fMCPrimMesonMom.Pt () ;
    mesonE   = fMCPrimMesonMom.E  () ;
    mesonYeta= fMCPrimMesonMom.Eta() ;
    mesonPhi = fMCPrimMesonMom.Phi() ;
    if( mesonPhi < 0 ) mesonPhi+=TMath::TwoPi();
    mesonPhi *= TMath::RadToDeg();
    
    ////
    Int_t genType = GetNCocktailGenNamesToCheck()-2; // bin 0 is not null 
    Int_t index   = GetReader()->GetCocktailGeneratorAndIndex(i, genName);
    //(GetMC())->GetCocktailGenerator(i,genName);
    
    Float_t weightPt = GetParticlePtWeight(mesonPt, pdg, genName, index) ; 
    
    if(IsStudyClusterOverlapsPerGeneratorOn())
    {
      for(Int_t igen = 1; igen < GetNCocktailGenNamesToCheck(); igen++)
      {       
        if ( GetCocktailGenNameToCheck(igen).Contains(genName) && 
            ( GetCocktailGenIndexToCheck(igen) < 0 || index == GetCocktailGenIndexToCheck(igen)) )
        {
          genType = igen-1;
          break;
        }
      }
    }
    
    //
    if(pdg == 321){//k+-
      if(TMath::Abs(mesonY) < 1.0 && fFillOriginHisto){
    	fhPrimChHadronPt->Fill(mesonPt, 0.5, GetEventWeight()*weightPt);
      }
    }
    else if(pdg == 211){//pi+-
      if(TMath::Abs(mesonY) < 1.0 && fFillOriginHisto){
    	fhPrimChHadronPt->Fill(mesonPt, 1.5, GetEventWeight()*weightPt);
      }
    }
    else if(pdg == 111)
    {
      if(TMath::Abs(mesonY) < 1.0)
      {
        fhPrimPi0E  ->Fill(mesonE ,           GetEventWeight()*weightPt) ;
        fhPrimPi0Pt ->Fill(mesonPt,           GetEventWeight()*weightPt) ;
        fhPrimPi0Phi->Fill(mesonPt, mesonPhi, GetEventWeight()*weightPt) ;
        
        fhPrimPi0YetaYcut->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
        
        if(IsStudyClusterOverlapsPerGeneratorOn())
        {
          fhPrimPi0PtPerGenerator [genType]->Fill(mesonPt,           GetEventWeight()*weightPt) ;
          fhPrimPi0PhiPerGenerator[genType]->Fill(mesonPt, mesonPhi, GetEventWeight()*weightPt) ;
        }
        
        if( IsHighMultiplicityAnalysisOn() )
        {
          fhPrimPi0PtCentrality->Fill(mesonPt, cen, GetEventWeight()*weightPt) ;
          fhPrimPi0PtEventPlane->Fill(mesonPt, ep , GetEventWeight()*weightPt) ;
        }
      }
      
      fhPrimPi0Y   ->Fill(mesonPt, mesonY   , GetEventWeight()*weightPt) ;
      fhPrimPi0Yeta->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
      
      if(inacceptance) fhPrimPi0PtInCalo->Fill(mesonPt,GetEventWeight()*weightPt);
      
      if(IsStudyClusterOverlapsPerGeneratorOn())
      {
        fhPrimPi0YPerGenerator[genType]->Fill(mesonPt, mesonY, GetEventWeight()*weightPt) ;  
        if(inacceptance) fhPrimPi0PtInCaloPerGenerator[genType]->Fill(mesonPt,GetEventWeight()*weightPt);
      }
    }
    else if(pdg == 221)
    {
      if(TMath::Abs(mesonY) < 1.0)
      {
        fhPrimEtaE  ->Fill(mesonE ,           GetEventWeight()*weightPt) ;
        fhPrimEtaPt ->Fill(mesonPt,           GetEventWeight()*weightPt) ;
        fhPrimEtaPhi->Fill(mesonPt, mesonPhi, GetEventWeight()*weightPt) ;
        
        if(IsStudyClusterOverlapsPerGeneratorOn())
        {
          fhPrimEtaPtPerGenerator [genType]->Fill(mesonPt,           GetEventWeight()*weightPt) ;
          fhPrimEtaPhiPerGenerator[genType]->Fill(mesonPt, mesonPhi, GetEventWeight()*weightPt) ;
        }
        
        fhPrimEtaYetaYcut->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
        
        if( IsHighMultiplicityAnalysisOn() )
        {
          fhPrimEtaPtCentrality->Fill(mesonPt, cen, GetEventWeight()*weightPt) ;
          fhPrimEtaPtEventPlane->Fill(mesonPt, ep , GetEventWeight()*weightPt) ;
        }
      }
      
      fhPrimEtaY   ->Fill(mesonPt, mesonY   , GetEventWeight()*weightPt) ;
      fhPrimEtaYeta->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
      
      if(inacceptance) fhPrimEtaPtInCalo->Fill(mesonPt,GetEventWeight()*weightPt);
      
      if(IsStudyClusterOverlapsPerGeneratorOn())
      {
        fhPrimEtaYPerGenerator[genType]->Fill(mesonPt, mesonY, GetEventWeight()*weightPt) ;
        if(inacceptance) fhPrimEtaPtInCaloPerGenerator[genType]->Fill(mesonPt,GetEventWeight()*weightPt);
      }
     }
    
    // Origin of meson
    if(fFillOriginHisto && TMath::Abs(mesonY) < 0.7)
    {
      Int_t momindex = primary->GetMother();
      Int_t mompdg    = -1;
      Int_t momstatus = -1;
      Int_t uniqueID  = primary->GetUniqueID();
      Int_t status    = -1;
      Float_t momR    =  0;
      Float_t momRcorr=  0;
      Double_t vertex[3]={0.,0.,0.};
      GetVertex(vertex);
      
      if(momindex >= 0 && momindex < nprim)
      {
        status = primary->MCStatusCode();
        AliVParticle* mother = GetMC()->GetTrack(momindex);
        mompdg    = TMath::Abs(mother->PdgCode());
        momstatus = mother->MCStatusCode();
        //momR      = mother->R();
        momR      = TMath::Sqrt(mother->Xv()*mother->Xv()+mother->Yv()*mother->Yv());
	momRcorr  = TMath::Sqrt((mother->Xv()-vertex[0])*(mother->Xv()-vertex[0]) +
				(mother->Yv()-vertex[1])*(mother->Yv()-vertex[1]));
	   
        if(pdg == 111)
        {
          fhPrimPi0ProdVertex->Fill(mesonPt, momR  , GetEventWeight()*weightPt);
          fhPrimPi0PtStatus  ->Fill(mesonPt, status, GetEventWeight()*weightPt);

	  if((uniqueID==13 && mompdg!=130 && mompdg!=310 && mompdg!=3122 && mompdg!=321 && mompdg!=211)
	     || mompdg==3212 || mompdg==3222 || mompdg==3112 || mompdg==3322 || mompdg==3212) {
	    //310 - K0S
	    //130 - K0L
	    //3133 - Lambda
	    //321 - K+-
	    //3212 - sigma0
	    //3222,3112 - sigma+-
	    //3322,3312 - cascades
	    fhPrimPi0PtOrigin->Fill(mesonPt, 14.5, GetEventWeight()*weightPt);//hadronic interaction
	  }
          else if(mompdg==310) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 6.5, GetEventWeight()*weightPt);//k0S
	  }
	  else if(mompdg    == 130) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 10.5, GetEventWeight()*weightPt);//k0L
	  }
	  else if(mompdg==321) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 11.5, GetEventWeight()*weightPt);//k+-
	  }
	  else if(mompdg==3122) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 13.5, GetEventWeight()*weightPt);//lambda 
	  }
	  else if(momRcorr>0.1) {//in cm
	    fhPrimPi0PtOrigin->Fill(mesonPt, 15.5, GetEventWeight()*weightPt);//radius too large
	  }
	  else if(mompdg==211) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 16.5, GetEventWeight()*weightPt);//pi+-
	  }
	  else if(mompdg==11) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 17.5, GetEventWeight()*weightPt);//e+-
	  }
	  else if(mompdg==13) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 18.5, GetEventWeight()*weightPt);//mu+-
	  }
	  else if(mompdg==2212 || mompdg==2112) {
	    fhPrimPi0PtOrigin->Fill(mesonPt, 19.5, GetEventWeight()*weightPt);//p,n
	  }
	  
          else if(momstatus  == 21) fhPrimPi0PtOrigin->Fill(mesonPt, 0.5, GetEventWeight()*weightPt);//parton
          else if(mompdg     < 22  && mompdg!=11 && mompdg!=13) fhPrimPi0PtOrigin->Fill(mesonPt, 1.5, GetEventWeight()*weightPt);//quark
          else if(mompdg     > 2100  && mompdg < 2210)
                                    fhPrimPi0PtOrigin->Fill(mesonPt, 2.5, GetEventWeight()*weightPt);// resonances
          else if(mompdg    == 221) fhPrimPi0PtOrigin->Fill(mesonPt, 8.5, GetEventWeight()*weightPt);//eta
          else if(mompdg    == 331) fhPrimPi0PtOrigin->Fill(mesonPt, 9.5, GetEventWeight()*weightPt);//eta prime
          else if(mompdg    == 213) fhPrimPi0PtOrigin->Fill(mesonPt, 4.5, GetEventWeight()*weightPt);//rho
          else if(mompdg    == 223) fhPrimPi0PtOrigin->Fill(mesonPt, 5.5, GetEventWeight()*weightPt);//omega
          else if(mompdg    >  310   && mompdg <= 323)
	                            fhPrimPi0PtOrigin->Fill(mesonPt, 12.5, GetEventWeight()*weightPt);//k*
          else if(momstatus == 11 || momstatus == 12 )
                                    fhPrimPi0PtOrigin->Fill(mesonPt, 3.5, GetEventWeight()*weightPt);//resonances
          else                      fhPrimPi0PtOrigin->Fill(mesonPt, 7.5, GetEventWeight()*weightPt);//other?
          
          //printf("Prim Pi0: index %d, pt %2.2f Prod vertex %3.3f, origin pdg %d, origin status %d, origin UI %d\n",
          //                   momindex, mesonPt,mother->R(),mompdg,momstatus,mother->GetUniqueID());
          
          if(status!=11)
          {
            if     (momstatus  == 21) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 0.5, GetEventWeight()*weightPt);//parton
            else if(mompdg     < 22 ) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 1.5, GetEventWeight()*weightPt);//quark
            else if(mompdg     > 2100  && mompdg < 2210)
                                      fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 2.5, GetEventWeight()*weightPt);// resonances
            else if(mompdg    == 221) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 8.5, GetEventWeight()*weightPt);//eta
            else if(mompdg    == 331) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 9.5, GetEventWeight()*weightPt);//eta prime
            else if(mompdg    == 213) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 4.5, GetEventWeight()*weightPt);//rho
            else if(mompdg    == 223) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 5.5, GetEventWeight()*weightPt);//omega
            else if(mompdg    >= 310   && mompdg <= 323)
                                      fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 6.5, GetEventWeight()*weightPt);//k0S, k+-,k*
            else if(mompdg    == 130) fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 6.5, GetEventWeight()*weightPt);//k0L
            else if(momstatus == 11 || momstatus == 12 )
                                      fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 3.5, GetEventWeight()*weightPt);//resonances
            else                      fhPrimNotResonancePi0PtOrigin->Fill(mesonPt, 7.5, GetEventWeight()*weightPt);//other?
          }
          
        }//pi0
        else
        {
          if     (momstatus == 21 ) fhPrimEtaPtOrigin->Fill(mesonPt, 0.5, GetEventWeight()*weightPt);//parton
          else if(mompdg    < 22  ) fhPrimEtaPtOrigin->Fill(mesonPt, 1.5, GetEventWeight()*weightPt);//quark
          else if(mompdg    > 2100  && mompdg   < 2210)
                                    fhPrimEtaPtOrigin->Fill(mesonPt, 2.5, GetEventWeight()*weightPt);//qq resonances
          else if(mompdg    == 331) fhPrimEtaPtOrigin->Fill(mesonPt, 5.5, GetEventWeight()*weightPt);//eta prime
          else if(momstatus == 11 || momstatus  == 12 )
                                    fhPrimEtaPtOrigin->Fill(mesonPt, 3.5, GetEventWeight()*weightPt);//resonances
          else                      fhPrimEtaPtOrigin->Fill(mesonPt, 4.5, GetEventWeight()*weightPt);//stable, conversions?
          //printf("Other Meson pdg %d, Mother %s, pdg %d, status %d\n",pdg, TDatabasePDG::Instance()->GetParticle(mompdg)->GetName(),mompdg, momstatus );
          
          fhPrimEtaProdVertex->Fill(mesonPt, momR, GetEventWeight()*weightPt);
        }
      } // pi0 has mother
    }
    
    // Check if both photons hit calorimeter or a given fiducial region 
    // only if those settings are specified.
    if ( !IsFiducialCutOn() && !IsRealCaloAcceptanceOn() ) continue ;
    
    if ( nDaught != 2 ) continue; //Only interested in 2 gamma decay
    
    if ( iphot1 < 0 || iphot1 >= nprim || iphot2 < 0 || iphot2 >= nprim ) continue ;
    
    Int_t pdg1 = 0;
    Int_t pdg2 = 0;
    Bool_t inacceptance1 = kTRUE;
    Bool_t inacceptance2 = kTRUE;
        
    AliVParticle * phot1 = GetMC()->GetTrack(iphot1) ;
    AliVParticle * phot2 = GetMC()->GetTrack(iphot2) ;
    
    if(!phot1 || !phot2) continue ;
    
    pdg1 = phot1->PdgCode();
    pdg2 = phot2->PdgCode();
    
    fPhotonMom1.SetPxPyPzE(phot1->Px(),phot1->Py(),phot1->Pz(),phot1->E());
    fPhotonMom2.SetPxPyPzE(phot2->Px(),phot2->Py(),phot2->Pz(),phot2->E());
    
    // Check if photons hit the Calorimeter acceptance
    if(IsRealCaloAcceptanceOn())
    {
      if( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance( GetCalorimeter(), phot1 )) inacceptance1 = kFALSE ;
      if( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance( GetCalorimeter(), phot2 )) inacceptance2 = kFALSE ;
    }
    
    if( pdg1 != 22 || pdg2 !=22) continue ;
    
    // Check if photons hit desired acceptance in the fidutial borders fixed in the analysis
    if(IsFiducialCutOn())
    {
      if( inacceptance1 && !GetFiducialCut()->IsInFiducialCut(fPhotonMom1.Eta(), fPhotonMom1.Phi(), GetCalorimeter()) ) inacceptance1 = kFALSE ;
      if( inacceptance2 && !GetFiducialCut()->IsInFiducialCut(fPhotonMom2.Eta(), fPhotonMom2.Phi(), GetCalorimeter()) ) inacceptance2 = kFALSE ;
    }
    
    if ( fFillArmenterosThetaStar ) FillArmenterosThetaStar(pdg);
    
    Int_t absID1=-1;
    Int_t absID2=-1;

    if(GetCalorimeter()==kEMCAL && inacceptance1 && inacceptance2 && fCheckAccInSector)
    {
      Float_t photonPhi1 = fPhotonMom1.Phi();
      Float_t photonPhi2 = fPhotonMom2.Phi();
      
      if(photonPhi1 < 0) photonPhi1+=TMath::TwoPi();
      if(photonPhi2 < 0) photonPhi2+=TMath::TwoPi();
      
      GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(fPhotonMom1.Eta(),photonPhi1,absID1);
      GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(fPhotonMom2.Eta(),photonPhi2,absID2);
      
      Int_t sm1 = GetEMCALGeometry()->GetSuperModuleNumber(absID1);
      Int_t sm2 = GetEMCALGeometry()->GetSuperModuleNumber(absID2);
      
      Int_t j=0;
      Bool_t sameSector = kFALSE;
      for(Int_t isector = 0; isector < fNModules/2; isector++)
      {
        j=2*isector;
        if((sm1==j && sm2==j+1) || (sm1==j+1 && sm2==j)) sameSector = kTRUE;
      }
      
      if(sm1!=sm2 && !sameSector)
      {
        inacceptance1 = kFALSE;
        inacceptance2 = kFALSE;
      }
      //if(sm1!=sm2)printf("sm1 %d, sm2 %d, same sector %d, in acceptance %d\n",sm1,sm2,sameSector,inacceptance);
      //                  if(GetEMCALGeometry()->Impact(phot1) && GetEMCALGeometry()->Impact(phot2))
      //                    inacceptance = kTRUE;
    }
    
    AliDebug(2,Form("Accepted in %s?: m (%2.2f,%2.2f,%2.2f), p1 (%2.2f,%2.2f,%2.2f), p2 (%2.2f,%2.2f,%2.2f) : in1 %d, in2 %d",
                    GetCalorimeterString().Data(),
                    mesonPt,mesonYeta,mesonPhi,
                    fPhotonMom1.Pt(),fPhotonMom1.Eta(),fPhotonMom1.Phi()*TMath::RadToDeg(),
                    fPhotonMom2.Pt(),fPhotonMom2.Eta(),fPhotonMom2.Phi()*TMath::RadToDeg(),
                    inacceptance1, inacceptance2));
    
    if(inacceptance1 && inacceptance2)
    {
      Float_t  asym  = TMath::Abs((fPhotonMom1.E()-fPhotonMom2.E()) / (fPhotonMom1.E()+fPhotonMom2.E()));
      Double_t angle = fPhotonMom1.Angle(fPhotonMom2.Vect());
      
      Bool_t cutAngle = kFALSE;
      if(fUseAngleCut && (angle < fAngleCut || angle > fAngleMaxCut)) cutAngle = kTRUE;

      Bool_t cutSeparation = kFALSE;
      if(fUseOneCellSeparation)
      {
	if(CheckSeparation(absID1,absID2)) cutSeparation = kTRUE;
      }
      
      AliDebug(2,Form("\t ACCEPTED pdg %d: pt %2.2f, phi %2.2f, eta %2.2f",pdg,mesonPt,mesonPhi,mesonYeta));
      
      if(pdg==111)
      {
        fhPrimPi0AccE   ->Fill(mesonE ,            GetEventWeight()*weightPt) ;
        fhPrimPi0AccPt  ->Fill(mesonPt,            GetEventWeight()*weightPt) ;
        fhPrimPi0AccPhi ->Fill(mesonPt, mesonPhi , GetEventWeight()*weightPt) ;
        fhPrimPi0AccY   ->Fill(mesonPt, mesonY   , GetEventWeight()*weightPt) ;
        fhPrimPi0AccYeta->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
        
        if(IsStudyClusterOverlapsPerGeneratorOn())
          fhPrimPi0AccPtPerGenerator[genType]  ->Fill(mesonPt,GetEventWeight()*weightPt) ;

        if( IsHighMultiplicityAnalysisOn() )
        {
          fhPrimPi0AccPtCentrality->Fill(mesonPt, cen, GetEventWeight()*weightPt) ;
          fhPrimPi0AccPtEventPlane->Fill(mesonPt, ep , GetEventWeight()*weightPt) ;
        }
        
        if(fFillAngleHisto)
        {
          fhPrimPi0OpeningAngle    ->Fill(mesonPt, angle, GetEventWeight()*weightPt);
          if(mesonPt > 5)fhPrimPi0OpeningAngleAsym->Fill(asym, angle, GetEventWeight()*weightPt);
          fhPrimPi0CosOpeningAngle ->Fill(mesonPt, TMath::Cos(angle), GetEventWeight()*weightPt);
        }
        
        if(fPhotonMom1.Pt() > GetMinPt() && fPhotonMom2.Pt() > GetMinPt() && !cutAngle && !cutSeparation)
        {
          fhPrimPi0AccPtPhotonCuts->Fill(mesonPt, GetEventWeight()*weightPt) ;
          
          if(IsStudyClusterOverlapsPerGeneratorOn())
            fhPrimPi0AccPtPhotonCutsPerGenerator[genType]  ->Fill(mesonPt,GetEventWeight()*weightPt) ;

          if(fFillAngleHisto)
            fhPrimPi0OpeningAnglePhotonCuts->Fill(mesonPt, angle, GetEventWeight()*weightPt);
          
          if ( fNAngleCutBins > 0 && fFillOpAngleCutHisto )
          {
            Int_t angleBin = -1;
            for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
            {
              if(angle >= fAngleCutBinsArray[ibin] && 
                 angle <  fAngleCutBinsArray[ibin+1]) angleBin = ibin;
            }
            
            if( angleBin >= 0 && angleBin < fNAngleCutBins)
              fhPrimPi0AccPtOpAngCuts[angleBin]->Fill(mesonPt,GetEventWeight()*weightPt);
          }
        }
      }
      else if(pdg==221)
      {
        fhPrimEtaAccE   ->Fill(mesonE ,            GetEventWeight()*weightPt) ;
        fhPrimEtaAccPt  ->Fill(mesonPt,            GetEventWeight()*weightPt) ;
        fhPrimEtaAccPhi ->Fill(mesonPt, mesonPhi , GetEventWeight()*weightPt) ;
        fhPrimEtaAccY   ->Fill(mesonPt, mesonY   , GetEventWeight()*weightPt) ;
        fhPrimEtaAccYeta->Fill(mesonPt, mesonYeta, GetEventWeight()*weightPt) ;
        
        if(IsStudyClusterOverlapsPerGeneratorOn())
          fhPrimEtaAccPtPerGenerator[genType]  ->Fill(mesonPt,GetEventWeight()*weightPt) ;
        
        if( IsHighMultiplicityAnalysisOn() )
        {
          fhPrimEtaAccPtCentrality->Fill(mesonPt, cen, GetEventWeight()*weightPt) ;
          fhPrimEtaAccPtEventPlane->Fill(mesonPt, ep , GetEventWeight()*weightPt) ;
        }
        
        if(fFillAngleHisto)
        {
          fhPrimEtaOpeningAngle    ->Fill(mesonPt, angle, GetEventWeight()*weightPt);
          if(mesonPt > 5)fhPrimEtaOpeningAngleAsym->Fill(asym, angle, GetEventWeight()*weightPt);
          fhPrimEtaCosOpeningAngle ->Fill(mesonPt, TMath::Cos(angle), GetEventWeight()*weightPt);
        }  
        
        if(fPhotonMom1.Pt() > GetMinPt() && fPhotonMom2.Pt() > GetMinPt() && !cutAngle && !cutSeparation)
        {
          fhPrimEtaAccPtPhotonCuts->Fill(mesonPt, GetEventWeight()*weightPt) ;
          
          if(IsStudyClusterOverlapsPerGeneratorOn())
            fhPrimEtaAccPtPhotonCutsPerGenerator[genType]->Fill(mesonPt,GetEventWeight()*weightPt) ;
          
          if(fFillAngleHisto)
            fhPrimEtaOpeningAnglePhotonCuts->Fill(mesonPt, angle, GetEventWeight()*weightPt);
          
          if ( fNAngleCutBins > 0 && fFillOpAngleCutHisto )
          {
            Int_t angleBin = -1;
            for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
            {
              if(angle >= fAngleCutBinsArray[ibin] && 
                 angle <  fAngleCutBinsArray[ibin+1]) angleBin = ibin;
            }
            
            if( angleBin >= 0 && angleBin < fNAngleCutBins)
              fhPrimEtaAccPtOpAngCuts[angleBin]->Fill(mesonPt,GetEventWeight()*weightPt);
          }
        }
      }
    } // Accepted
  } // loop on primaries
}

//________________________________________________
/// Fill armenteros plots.
//________________________________________________
void AliAnaPi0::FillArmenterosThetaStar(Int_t pdg)
{
  // Get pTArm and AlphaArm
  Float_t momentumSquaredMother = fMCPrimMesonMom.P()*fMCPrimMesonMom.P();
  Float_t momentumDaughter1AlongMother = 0.;
  Float_t momentumDaughter2AlongMother = 0.;
  
  if (momentumSquaredMother > 0.)
  {
    momentumDaughter1AlongMother = (fPhotonMom1.Px()*fMCPrimMesonMom.Px() + fPhotonMom1.Py()*fMCPrimMesonMom.Py()+ fPhotonMom1.Pz()*fMCPrimMesonMom.Pz()) / sqrt(momentumSquaredMother);
    momentumDaughter2AlongMother = (fPhotonMom2.Px()*fMCPrimMesonMom.Px() + fPhotonMom2.Py()*fMCPrimMesonMom.Py()+ fPhotonMom2.Pz()*fMCPrimMesonMom.Pz()) / sqrt(momentumSquaredMother);
  }
  
  Float_t momentumSquaredDaughter1 = fPhotonMom1.P()*fPhotonMom1.P();
  Float_t ptArmSquared = momentumSquaredDaughter1 - momentumDaughter1AlongMother*momentumDaughter1AlongMother;
  
  Float_t pTArm = 0.;
  if (ptArmSquared > 0.)
    pTArm = sqrt(ptArmSquared);
  
  Float_t alphaArm = 0.;
  if(momentumDaughter1AlongMother + momentumDaughter2AlongMother > 0)
    alphaArm = (momentumDaughter1AlongMother -momentumDaughter2AlongMother) / (momentumDaughter1AlongMother + momentumDaughter2AlongMother);
  
  fPhotonMom1Boost = fPhotonMom1;
  fPhotonMom1Boost.Boost(-fMCPrimMesonMom.BoostVector());
  Float_t  cosThStar=TMath::Cos(fPhotonMom1Boost.Vect().Angle(fMCPrimMesonMom.Vect()));
  
  Float_t en   = fMCPrimMesonMom.Energy();
  Int_t   ebin = -1;
  if(en > 8  && en <= 12) ebin = 0;
  if(en > 12 && en <= 16) ebin = 1;
  if(en > 16 && en <= 20) ebin = 2;
  if(en > 20)             ebin = 3;
  if(ebin < 0 || ebin > 3) return ;
  
  if(pdg==111)
  {
    fhCosThStarPrimPi0->Fill(en,       cosThStar, GetEventWeight());
    fhArmPrimPi0[ebin]->Fill(alphaArm, pTArm    , GetEventWeight());
  }
  else
  {
    fhCosThStarPrimEta->Fill(en,      cosThStar, GetEventWeight());
    fhArmPrimEta[ebin]->Fill(alphaArm,pTArm    , GetEventWeight());
  }
  
  //  if(GetDebug() > 2 )
  //  {
  //    Float_t asym = 0;
  //    if(fPhotonMom1.E()+fPhotonMom2.E() > 0) asym = TMath::Abs(fPhotonMom1.E()-fPhotonMom2.E())/(fPhotonMom1.E()+fPhotonMom2.E());
  //
  //    printf("AliAnaPi0::FillArmenterosThetaStar() - E %f, alphaArm %f, pTArm %f, cos(theta*) %f, asymmetry %1.3f\n",
  //         en,alphaArm,pTArm,cosThStar,asym);
  //  }
}

//_______________________________________________________________________________________
/// Do some MC checks on the origin of the pair, is there any common ancestor and if there is one, who?
/// Adjusted for Pythia, need to see what to do for other generators.
/// Array of histograms ordered as follows: 0-Photon, 1-electron, 2-pi0, 3-eta, 4-a-proton, 5-a-neutron, 6-stable particles,
/// 7-other decays, 8-string, 9-final parton, 10-initial parton, intermediate, 11-colliding proton, 12-unrelated
//_______________________________________________________________________________________
void AliAnaPi0::FillMCVersusRecDataHistograms(Int_t ancLabel , Int_t ancPDG, 
                                              Int_t ancStatus, Double_t weightPt,
                                              Int_t iclus1,    Int_t iclus2,
                                              Int_t mctag1,    Int_t mctag2,
                                              Float_t pt1,     Float_t pt2,
                                              Int_t ncell1,    Int_t ncell2,
                                              Double_t mass,   Double_t pt,   Double_t asym,
                                              Double_t deta,   Double_t dphi, Double_t angle)
{  
  Int_t momindex  = -1;
  Int_t mompdg    = -1;
  Int_t momstatus = -1;
  Int_t status    = -1;
  Int_t uniqueID  = -1;
  Float_t prodR = -1;
  Float_t prodRcorr = -1;
  Int_t mcIndex = -1;
  Float_t ptPrim = fMCPrimMesonMom.Pt();
  Float_t ratPrim = 0;
  if ( ptPrim > 0.1 ) ratPrim = (pt-ptPrim)/ptPrim;
  Float_t cent   = GetEventCentrality();
  Double_t vertex[3]={0.,0.,0.};
  GetVertex(vertex);
  if ( ancLabel > -1 )
  {
    AliDebug(1,Form("Common ancestor label %d, pdg %d, name %s, status %d",
                    ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus));
    
    AliVParticle * mom = GetMC()->GetTrack(ancLabel);

    if ( ancPDG == 22 )
    { // gamma, from conversions?
      if ( GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
           GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)    )
        mcIndex = 0;
      else
        mcIndex = 14;
    }
    else if ( TMath::Abs(ancPDG) == 11 )
    {// electron, from conversions?
      if ( GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
           GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)    )
        mcIndex = 1;
      else
        mcIndex = 15;
    }
    else if ( ancPDG==111 &&  mom )
    {//Pi0
      Bool_t ok= kTRUE;
      
      if ( mom->GetNDaughters()!=2 )
      {
        ok = kFALSE;
      }
      else
      {
        AliVParticle * d1 = GetMC()->GetTrack(mom->GetDaughterLabel(0));
        AliVParticle * d2 = GetMC()->GetTrack(mom->GetDaughterLabel(1));

        if      ( !d1 || !d2 )
          ok = kFALSE;
        else if ( d1->PdgCode() != 22 || d2->PdgCode() != 22 )
          ok = kFALSE;
      }
      
      if ( GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCPi0Decay) &&
           GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCPi0Decay) && ok  )
      {
        mcIndex = 2;
        
        if ( IsHighMultiplicityAnalysisOn() )
        {
          fhMCPi0MassPtTrueCen    ->Fill(ptPrim, mass , cent, GetEventWeight()*weightPt);
          fhMCPi0MassPtRecCen     ->Fill(pt    , mass , cent, GetEventWeight()*weightPt);
          fhMCPi0PtTruePtRecCen   ->Fill(ptPrim,   pt , cent, GetEventWeight()*weightPt);
          fhMCPi0PtTruePtRecDifOverPtTrueCen->Fill(pt , ratPrim, cent, GetEventWeight()*weightPt);
          
          if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] ) 
          {
            fhMCPi0PtTruePtRecMassCutCen->Fill(ptPrim, pt, cent, GetEventWeight()*weightPt);
            fhMCPi0PtTruePtRecDifOverPtTrueCenMassCut->Fill(pt, ratPrim, cent, GetEventWeight()*weightPt);
          }
        }
        else if ( fMultiCutAnaSim )
        {
          for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
          {
            for(Int_t icell=0; icell<fNCellNCuts; icell++)
            {
              for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
              {
                Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                if(pt1    >  fPtCuts[ipt]      && pt2    >  fPtCuts[ipt]        &&
                   pt1    <  fPtCutsMax[ipt]   && pt2    <  fPtCutsMax[ipt]     &&
                   asym   <  fAsymCuts[iasym]                                   &&
                   ncell1 >= fCellNCuts[icell] && ncell2 >= fCellNCuts[icell])
                {
                  fhMCPi0MassPtRec [index]->Fill(pt    ,mass, GetEventWeight()*weightPt);
                  fhMCPi0MassPtTrue[index]->Fill(ptPrim,mass, GetEventWeight()*weightPt);
                  
                  if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] )
                    fhMCPi0PtTruePtRecMassCut[index]->Fill(ptPrim, pt, GetEventWeight()*weightPt);
                }//pass the different cuts
              }// pid bit cut loop
            }// icell loop
          }// pt cut loop
        }// Multi cut ana sim
        else
        {
          fhMCPi0PtTruePtRecDifOverPtTrue->Fill(pt, ratPrim, GetEventWeight()*weightPt);
          
          if ( fFillAngleHisto )
            fhMCPi0PtRecOpenAngle->Fill(pt, angle    , GetEventWeight()*weightPt);
          
          if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] )
          {
            fhMCPi0PtTruePtRecDifOverPtTrueMassCut->Fill(pt, ratPrim, GetEventWeight()*weightPt);
            
            if  ( fFillAngleHisto )
              fhMCPi0PtRecOpenAngleMassCut->Fill(pt, angle    , GetEventWeight()*weightPt);
          }
          
          fhMCPi0MassPtTrue [0]->Fill(ptPrim, mass, GetEventWeight()*weightPt);
          fhMCPi0MassPtRec  [0]->Fill(pt    , mass, GetEventWeight()*weightPt);
          fhMCPi0PtTruePtRec[0]->Fill(ptPrim,   pt, GetEventWeight()*weightPt);
          
          if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] )        
          {
            fhMCPi0PtTruePtRecMassCut[0]->Fill(ptPrim, pt, GetEventWeight()*weightPt);
            
            AliVParticle* ancestor = GetMC()->GetTrack(ancLabel);
            status   = ancestor->MCStatusCode();
            momindex = ancestor->GetMother();
            uniqueID = ancestor->GetUniqueID();
            
            Bool_t momOK = kFALSE;
            if ( momindex >= 0 )
            {
              AliVParticle* mother = GetMC()->GetTrack(momindex);
              mompdg    = TMath::Abs(mother->PdgCode());
              momstatus = mother->MCStatusCode();
              prodR = TMath::Sqrt(mother->Xv()*mother->Xv()+mother->Yv()*mother->Yv());
              prodRcorr = TMath::Sqrt((mother->Xv()-vertex[0])*(mother->Xv()-vertex[0]) +
                                      (mother->Yv()-vertex[1])*(mother->Yv()-vertex[1]));
              //prodR = mother->R();
              //uniqueId = mother->GetUniqueID();
              momOK = kTRUE;
            }
            
            if ( momOK )
            {
              //            printf("Reco Pi0: pt %2.2f Prod vertex %3.3f, origin pdg %d, origin status %d, origin UI %d\n",
              //                   pt,prodR,mompdg,momstatus,uniqueId);
              
              fhMCPi0ProdVertex->Fill(pt, prodR , GetEventWeight()*weightPt);
              fhMCPi0PtStatus  ->Fill(pt, status, GetEventWeight()*weightPt);
              fhMCPi0Radius    ->Fill(pt, prodRcorr,GetEventWeight()*weightPt);
              
              if ( (uniqueID==13 && mompdg!=130 && mompdg!=310 && mompdg!=3122 && mompdg!=321)
                  || mompdg==3212 || mompdg==3222 || mompdg==3112 || mompdg==3322 || mompdg==3212) {
                //310 - K0S
                //130 - K0L
                //3133 - Lambda
                //321 - K+-
                //3212 - sigma0
                //3222,3112 - sigma+-
                //3322,3312 - cascades
                fhMCPi0PtOrigin->Fill(pt, 14.5, GetEventWeight()*weightPt);//hadronic interaction
              }
              else if(mompdg==310) {
                fhMCPi0PtOrigin->Fill(pt, 6.5, GetEventWeight()*weightPt);//k0S
              }
              else if(mompdg    == 130) {
                fhMCPi0PtOrigin->Fill(pt, 10.5, GetEventWeight()*weightPt);//k0L
              }
              else if(mompdg==321) {
                fhMCPi0PtOrigin->Fill(pt, 11.5, GetEventWeight()*weightPt);//k+-
              }
              else if(mompdg==3122) {
                fhMCPi0PtOrigin->Fill(pt, 13.5, GetEventWeight()*weightPt);//lambda 
              }
              else if(prodRcorr>0.1) {//in cm
                fhMCPi0PtOrigin->Fill(pt, 15.5, GetEventWeight()*weightPt);//radius too large
              }
              else if(mompdg==2212 || mompdg==2112 || mompdg==211) {
                fhMCPi0PtOrigin->Fill(pt, 16.5, GetEventWeight()*weightPt);//p,n,pi
              }
              else if(mompdg==11) {
                fhMCPi0PtOrigin->Fill(pt, 17.5, GetEventWeight()*weightPt);//e+-
              }
              else if(mompdg==13) {
                fhMCPi0PtOrigin->Fill(pt, 18.5, GetEventWeight()*weightPt);//mu+-
              }
              else if     (momstatus  == 21) {
                fhMCPi0PtOrigin->Fill(pt, 0.5, GetEventWeight()*weightPt);//parton (empty for py8)
              }
              else if(mompdg     < 22 && mompdg!=11 && mompdg!=13) {
                fhMCPi0PtOrigin->Fill(pt, 1.5, GetEventWeight()*weightPt);//quark (no e nor mu)
              }
              else if(mompdg     > 2100  && mompdg   < 2210) {
                fhMCPi0PtOrigin->Fill(pt, 2.5, GetEventWeight()*weightPt);// resonances
              }
              else if(mompdg    == 221) {
                fhMCPi0PtOrigin->Fill(pt, 8.5, GetEventWeight()*weightPt);//eta
              }
              else if(mompdg    == 331) {
                fhMCPi0PtOrigin->Fill(pt, 9.5, GetEventWeight()*weightPt);//eta prime
              }
              else if(mompdg    == 213) {
                fhMCPi0PtOrigin->Fill(pt, 4.5, GetEventWeight()*weightPt);//rho
              }
              else if(mompdg    == 223) {
                fhMCPi0PtOrigin->Fill(pt, 5.5, GetEventWeight()*weightPt);//omega
              }
              else if(mompdg    > 310   && mompdg    <= 323) {
                fhMCPi0PtOrigin->Fill(pt, 12.5, GetEventWeight()*weightPt);//k*
              }
              else if(momstatus == 11 || momstatus  == 12 ) {
                fhMCPi0PtOrigin->Fill(pt, 3.5, GetEventWeight()*weightPt);//resonances
              }
              else {
                fhMCPi0PtOrigin->Fill(pt, 7.5, GetEventWeight()*weightPt);//other? py8 resonances?
              }
              
              if(status!=11) // all the time for py8
              {
                if     (momstatus  == 21) fhMCNotResonancePi0PtOrigin->Fill(pt, 0.5, GetEventWeight()*weightPt);//parton (empty for py8)
                else if(mompdg     < 22 ) fhMCNotResonancePi0PtOrigin->Fill(pt, 1.5, GetEventWeight()*weightPt);//quark
                else if(mompdg     > 2100  && mompdg   < 2210)
                  fhMCNotResonancePi0PtOrigin->Fill(pt, 2.5, GetEventWeight()*weightPt);// resonances
                else if(mompdg    == 221) fhMCNotResonancePi0PtOrigin->Fill(pt, 8.5, GetEventWeight()*weightPt);//eta
                else if(mompdg    == 331) fhMCNotResonancePi0PtOrigin->Fill(pt, 9.5, GetEventWeight()*weightPt);//eta prime
                else if(mompdg    == 213) fhMCNotResonancePi0PtOrigin->Fill(pt, 4.5, GetEventWeight()*weightPt);//rho
                else if(mompdg    == 223) fhMCNotResonancePi0PtOrigin->Fill(pt, 5.5, GetEventWeight()*weightPt);//omega
                else if(mompdg    >= 310   && mompdg    <= 323)
                  fhMCNotResonancePi0PtOrigin->Fill(pt, 6.5, GetEventWeight()*weightPt);//k0S, k+-,k*
                else if(mompdg    == 130) fhMCNotResonancePi0PtOrigin->Fill(pt, 6.5, GetEventWeight()*weightPt);//k0L
                else if(momstatus == 12 )
                  fhMCNotResonancePi0PtOrigin->Fill(pt, 3.5, GetEventWeight()*weightPt);//resonances
                else                      fhMCNotResonancePi0PtOrigin->Fill(pt, 7.5, GetEventWeight()*weightPt);//other? py8 resonances?
              }
            }
            
          }//pi0 mass region
        }
        
        if ( fNAngleCutBins > 0 && fFillOpAngleCutHisto )
        {
          Int_t angleBin = -1;
          for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
          {
            if(angle >= fAngleCutBinsArray[ibin] && 
               angle <  fAngleCutBinsArray[ibin+1]) angleBin = ibin;
          }
          
          if( angleBin >= 0 && angleBin < fNAngleCutBins)
            fhReOpAngleBinPairClusterMassMCTruePi0[angleBin]->Fill(pt, mass, GetEventWeight()*weightPt);
        }
      }
      else
      {
        mcIndex = 12; // other decays than pi0->gamma-gamma or merged
        
//        printf("Et (%2.2f, %2.2f), Mass %f, decay (%d, %d), merged (%d, %d), conv (%d, %d) ok %d, N daug %d\n",
//               pt1, pt2, mass,
//               GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCPi0Decay),
//               GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCPi0Decay),
//               GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCPi0),  
//               GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCPi0),  
//               GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion),  
//               GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion),
//               ok, mom->GetNDaughters());
//        
//        if ( mom->GetNDaughters()==2 )
//        {
//          AliVParticle * d1 = GetMC()->GetTrack(mom->GetDaughterLabel(0));
//          AliVParticle * d2 = GetMC()->GetTrack(mom->GetDaughterLabel(1));
//          printf("\t labels (%d, %d), pdg (%d, %d) \n",
//                 mom->GetDaughterLabel(0), mom->GetDaughterLabel(1),
//                 d1->PdgCode()           , d2->PdgCode()            ); 
//        }
//        
//        Int_t first = 0;
//        AliVCluster * cluster1 = FindCluster(GetEMCALClusters(),iclus1,first);
//        AliVCluster * cluster2 = FindCluster(GetEMCALClusters(),iclus2,iclus1);
//        if(!cluster2 || !cluster1) 
//        { 
//          AliWarning(Form("Cluster1 %p or Cluster 2 %p not found!",cluster1,cluster2));
//          return;
//        }
//
//        printf("xxx Cluster1\n");
//        for(Int_t ilab = 0; ilab < cluster1->GetNLabels(); ilab++ )
//        {
//          printf("label %d\n",ilab);
//          GetMCAnalysisUtils()->PrintAncestry(GetMC(), cluster1->GetLabels()[ilab]);
//        }
//
//        printf("xxx Cluster2\n");
//        for(Int_t ilab = 0; ilab < cluster2->GetNLabels(); ilab++ )
//        {
//          printf("label %d\n",ilab);
//          GetMCAnalysisUtils()->PrintAncestry(GetMC(), cluster2->GetLabels()[ilab]);
//        }
      }

    }
    else if ( ancPDG==221 && mom )
    {//Eta
      Bool_t ok= kTRUE;
      
      if ( mom->GetNDaughters()!=2 ) 
      {
        ok = kFALSE;
      }
      else
      {
        AliVParticle * d1 = GetMC()->GetTrack(mom->GetDaughterLabel(0));
        AliVParticle * d2 = GetMC()->GetTrack(mom->GetDaughterLabel(1));
        if      ( !d1 || !d2 )
          ok = kFALSE;
        else if ( d1->PdgCode() != 22 || d2->PdgCode() != 22 )
          ok = kFALSE;
      }
      
      if ( GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCEtaDecay) &&
           GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCEtaDecay) && ok )
      {        
        mcIndex = 3;
        
        if ( IsHighMultiplicityAnalysisOn() )
        {
          fhMCEtaMassPtTrueCen    ->Fill(ptPrim, mass , cent, GetEventWeight()*weightPt);
          fhMCEtaMassPtRecCen     ->Fill(pt    , mass , cent, GetEventWeight()*weightPt);
          fhMCEtaPtTruePtRecCen   ->Fill(ptPrim,   pt , cent, GetEventWeight()*weightPt);
          fhMCEtaPtTruePtRecDifOverPtTrueCen->Fill(pt , ratPrim, cent, GetEventWeight()*weightPt);
          
          if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] ) 
          {
            fhMCEtaPtTruePtRecMassCutCen->Fill(ptPrim, pt, cent, GetEventWeight()*weightPt);
            fhMCEtaPtTruePtRecDifOverPtTrueCenMassCut->Fill(pt, ratPrim, cent, GetEventWeight()*weightPt);
          }
        }
        else if ( fMultiCutAnaSim )
        {
          for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
          {
            for(Int_t icell=0; icell<fNCellNCuts; icell++)
            {
              for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
              {
                Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                if(pt1    >  fPtCuts[ipt]      && pt2    >  fPtCuts[ipt]        &&
                   pt1    <  fPtCutsMax[ipt]   && pt2    <  fPtCutsMax[ipt]     &&
                   asym   <  fAsymCuts[iasym]                                   &&
                   ncell1 >= fCellNCuts[icell] && ncell2 >= fCellNCuts[icell])
                {
                  fhMCEtaMassPtRec  [index]->Fill(pt    , mass, GetEventWeight()*weightPt);
                  fhMCEtaMassPtTrue [index]->Fill(ptPrim, mass, GetEventWeight()*weightPt);
                  fhMCEtaPtTruePtRec[index]->Fill(ptPrim,   pt, GetEventWeight()*weightPt);
                  
                  if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] )
                    fhMCEtaPtTruePtRecMassCut[index]->Fill(ptPrim, pt, GetEventWeight()*weightPt);
                }//pass the different cuts
              }// pid bit cut loop
            }// icell loop
          }// pt cut loop
        } //Multi cut ana sim
        else
        {      
          fhMCEtaPtTruePtRecDifOverPtTrue->Fill(pt, ratPrim, GetEventWeight()*weightPt);
          
          if ( fFillAngleHisto )
            fhMCEtaPtRecOpenAngle->Fill(pt, angle    , GetEventWeight()*weightPt);
          
          if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] )
          {
            fhMCEtaPtTruePtRecDifOverPtTrueMassCut->Fill(pt, ratPrim, GetEventWeight()*weightPt);
            
            if ( fFillAngleHisto )
              fhMCEtaPtRecOpenAngleMassCut->Fill(pt, angle    , GetEventWeight()*weightPt);
          }
          
          fhMCEtaMassPtTrue [0]->Fill(ptPrim, mass, GetEventWeight()*weightPt);
          fhMCEtaMassPtRec  [0]->Fill(pt    , mass, GetEventWeight()*weightPt);
          fhMCEtaPtTruePtRec[0]->Fill(ptPrim,   pt, GetEventWeight()*weightPt);
          
          if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] ) 
          {
            fhMCEtaPtTruePtRecMassCut[0]->Fill(ptPrim, pt, GetEventWeight()*weightPt);
            
            Bool_t momOK = kFALSE;
            
            AliVParticle* ancestor = GetMC()->GetTrack(ancLabel);
            momindex  = ancestor->GetMother();
            if ( momindex >= 0 ) 
            {
              AliVParticle* mother = GetMC()->GetTrack(momindex);
              mompdg    = TMath::Abs(mother->PdgCode());
              momstatus = mother->MCStatusCode();
              //prodR = mother->R();
              prodR = TMath::Sqrt(mother->Xv()*mother->Xv()+mother->Yv()*mother->Yv());
              prodRcorr = TMath::Sqrt((mother->Xv()-vertex[0])*(mother->Xv()-vertex[0]) +
                                      (mother->Yv()-vertex[1])*(mother->Yv()-vertex[1]));
              momOK = kTRUE;
            }
            
            if ( momOK )
            {
              fhMCEtaProdVertex->Fill(pt, prodR, GetEventWeight()*weightPt);
              fhMCEtaRadius    ->Fill(pt, prodRcorr, GetEventWeight()*weightPt);
              
              
              if     (momstatus == 21 ) { // empty for py8
                fhMCEtaPtOrigin->Fill(pt, 0.5, GetEventWeight()*weightPt);//parton
              }
              else if(mompdg    < 22  ) {
                fhMCEtaPtOrigin->Fill(pt, 1.5, GetEventWeight()*weightPt);//quark, include parton for py8
              }
              else if(mompdg    > 2100  && mompdg  < 2210) {
                fhMCEtaPtOrigin->Fill(pt, 2.5, GetEventWeight()*weightPt);//qq resonances
              }
              else if(mompdg    == 331) {
                fhMCEtaPtOrigin->Fill(pt, 5.5, GetEventWeight()*weightPt);//eta prime
              }
              else if(momstatus == 11 || momstatus == 12 ) { // empty for py8
                fhMCEtaPtOrigin->Fill(pt, 3.5, GetEventWeight()*weightPt);//resonances
              }
              else {
                fhMCEtaPtOrigin->Fill(pt, 4.5, GetEventWeight()*weightPt);//stable, conversions? resonances for py8?
                //printf("Other Meson pdg %d, Mother %s, pdg %d, status %d\n",
                //pdg, TDatabasePDG::Instance()->GetParticle(mompdg)->GetName(),mompdg, momstatus );
              }
              
            } // mom ok
            
          }// eta mass region
        } 
        
        if ( fNAngleCutBins > 0 && fFillOpAngleCutHisto )
        {
          Int_t angleBin = -1;
          for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
          {
            if(angle >= fAngleCutBinsArray[ibin] && 
               angle <  fAngleCutBinsArray[ibin+1]) angleBin = ibin;
          }
          
          if( angleBin >= 0 && angleBin < fNAngleCutBins)
            fhReOpAngleBinPairClusterMassMCTrueEta[angleBin]->Fill(pt, mass, GetEventWeight()*weightPt);
        }
      } // eta gamma gamma decays
      else 
      {
        mcIndex = 13; // other decays than eta->gamma-gamma or merged
      }
    }
    else if(ancPDG==-2212)
    {//AProton
      mcIndex = 4;
    }
    else if(ancPDG==-2112)
    {//ANeutron
      mcIndex = 5;
    }
    else if(TMath::Abs(ancPDG)==13)
    {//muons
      mcIndex = 6;
    }
    else if (TMath::Abs(ancPDG) > 100 && ancLabel > 7)
    {
      if(ancStatus==1) 
      {//Stable particles, converted? not decayed resonances
        mcIndex = 6;
      }
      else
      {//resonances and other decays, more hadron conversions?
        mcIndex = 7;
      }
    }
    else
    {//Partons, colliding protons, strings, intermediate corrections
      if(ancStatus == 11 || ancStatus == 12 || ancStatus == 0 )
      {//String fragmentation
        mcIndex = 8;
      }
      //else if (ancStatus==21)
      //{
      if(ancLabel < 2)
      {//Colliding protons
        mcIndex = 11;
      }//colliding protons
      else if(ancLabel < 6)
      {//partonic initial states interactions, not exactly for py8 it can include few more labels
        mcIndex = 9;
      }
      else if(ancLabel < 8)
      {//Final state partons radiations?, not exactly for py8 it can include few more labels
        mcIndex = 10;
      }
        // else {
        //   printf("AliAnaPi0::FillMCVersusRecDataHistograms() - Check ** Common ancestor label %d, pdg %d, name %s, status %d; \n",
        //          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
        // }
      //}//status 21
      //else {
      //  printf("AliAnaPi0::FillMCVersusRecDataHistograms() - Check *** Common ancestor label %d, pdg %d, name %s, status %d; \n",
      //         ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
      // }
    }////Partons, colliding protons, strings, intermediate corrections
  }//ancLabel > -1
  else
  { //ancLabel <= -1
    //printf("Not related at all label = %d\n",ancLabel);
    AliDebug(1,"Common ancestor not found");

    mcIndex = 16;
  }
  if ( mcIndex < 0 || mcIndex >= 17 ) 
  {
    AliInfo(Form("Wrong ancestor type %d, set it to unknown (12)",mcIndex));
    
    AliInfo(Form("\t Ancestor type not found: label %d, pdg %d, name %s, status %d\n",
           ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus));
    
    mcIndex = 16;
  }
  
  //printf("Pi0TASK: MC index %d\n",mcIndex);
  if ( !fFillOriginHistoForMesonsOnly )
  {
    fhMCOrgMass    [mcIndex]->Fill(pt, mass, GetEventWeight()*weightPt);
    fhMCOrgAsym    [mcIndex]->Fill(pt, asym, GetEventWeight()*weightPt);
    fhMCOrgDeltaEta[mcIndex]->Fill(pt, deta, GetEventWeight()*weightPt);
    fhMCOrgDeltaPhi[mcIndex]->Fill(pt, dphi, GetEventWeight()*weightPt);
  }

  if ( fCheckConversion )
  {
    if ( mcIndex == 2 ) {//pi0 only check conversions
      if(GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
         GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)) {
        //both conversions
        fhMCOrgPi0MassPtConversion[2]->Fill(pt, mass, GetEventWeight()*weightPt);
      } else if(!GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
                !GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)) {
        //no conversion
        fhMCOrgPi0MassPtConversion[0]->Fill(pt, mass, GetEventWeight()*weightPt);
      } else {
        //one conversion and one not
        fhMCOrgPi0MassPtConversion[1]->Fill(pt, mass, GetEventWeight()*weightPt);
      }
    }
    
    if ( mcIndex==3 ) {//eta only, check conversions
      if(GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
         GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)) {
        //both conversions
        fhMCOrgEtaMassPtConversion[2]->Fill(pt, mass, GetEventWeight()*weightPt);
      } else if(!GetMCAnalysisUtils()->CheckTagBit(mctag1,AliMCAnalysisUtils::kMCConversion) &&
                !GetMCAnalysisUtils()->CheckTagBit(mctag2,AliMCAnalysisUtils::kMCConversion)) {
        //no conversion
        fhMCOrgEtaMassPtConversion[0]->Fill(pt, mass, GetEventWeight()*weightPt);
      } else {
        //one conversion and one not
        fhMCOrgEtaMassPtConversion[1]->Fill(pt, mass, GetEventWeight()*weightPt);
      }
    }
  } // conversion checks
  
  if ( IsStudyClusterOverlapsPerGeneratorOn() )
  {      
    // Get the original clusters
    //
    Int_t first = 0;
    AliVCluster * cluster1 = FindCluster(GetEMCALClusters(),iclus1,first);
    AliVCluster * cluster2 = FindCluster(GetEMCALClusters(),iclus2,iclus1);
    if(!cluster2 || !cluster1) 
    { 
      AliWarning(Form("Cluster1 %p or Cluster 2 %p not found!",cluster1,cluster2));
      return;
    }
    
    // Get the generators names of each cluster and background generator tag
    //
    TString genName1 = "", genName2 = "", genNameBkg1 = "", genNameBkg2 = "";
    Int_t indexGen1 = -1, indexGen2 = -1, indexGenBkg1 = -1, indexGenBkg2 = -1;
    Int_t genBkgTag1 = GetCocktailGeneratorBackgroundTag(cluster1, mctag1, genName1, indexGen1, genNameBkg1, indexGenBkg1);
    if     (genBkgTag1 == -1) return;
    else if(genBkgTag1  >  3) printf("Bkg1 generator tag larger than 3; Main %s Bkg %s\n",genName1.Data(),genNameBkg1.Data());

    Int_t genBkgTag2 = GetCocktailGeneratorBackgroundTag(cluster2, mctag2, genName2, indexGen2, genNameBkg2, indexGenBkg2);
    if     (genBkgTag2 == -1) return;
    else if(genBkgTag2  >  3) printf("Bkg2 generator tag larger than 3; Main %s Bkg %s\n",genName2.Data(),genNameBkg2.Data());
    
    // if the clusters do not come from the same generator exclude
    if(genName1!=genName2) return;
    
    Int_t genType = GetNCocktailGenNamesToCheck()-1;
    for(Int_t igen = 1; igen < GetNCocktailGenNamesToCheck(); igen++)
    {
      if ( GetCocktailGenNameToCheck(igen).Contains(genName1) && 
          ( GetCocktailGenIndexToCheck(igen) < 0 || indexGen1 == GetCocktailGenIndexToCheck(igen)))
      {
        genType = igen;
        break;
      }
    }
    
    // Fill histograms
    //
    // assign histogram index number depending on pair combination
    Int_t tag = -1;
    if ( genBkgTag1 == genBkgTag2 )
    {
      tag = genBkgTag1; // 0,1,2,3
    }
    else
    {
      Int_t genBkgMin = -1; 
      Int_t genBkgMax = -1;
      
      if(genBkgTag1 > genBkgTag2)
      {
        genBkgMax = genBkgTag1;
        genBkgMin = genBkgTag2;
      }
      else
      {
        genBkgMax = genBkgTag2;
        genBkgMin = genBkgTag1;
      }
    
      if      ( genBkgMin == 0 )
      {
        if     (genBkgMax == 1 )   tag = 4; // clean+hijing
        else if(genBkgMax == 2 )   tag = 5; // clean+not hijing
        else if(genBkgMax == 3 )   tag = 6; // clean+multiple
      }
      else if ( genBkgMin == 1 ) 
      {
        if      ( genBkgMax == 2 ) tag = 7; // hijing + not hijing
        else if ( genBkgMax == 3 ) tag = 8; // hijing + multiple
      }
      else if ( genBkgMin == 2 )   tag = 9; // not hijing + multiple
    }
    
    if(tag == -1)
    {
      printf("Combination not found, bkg1 tag %d, bkg2 tag %d\n",genBkgTag1,genBkgTag2);
      return;
    }
    
    fhPairGeneratorsBkgMass[genType][tag]->Fill(pt, mass, GetEventWeight()*weightPt);
    fhPairGeneratorsBkgMass      [0][tag]->Fill(pt, mass, GetEventWeight()*weightPt);
    
    //
    if ( ptPrim < 0.1 || pt < 0.5 ) return;
    
    Float_t ratio = pt / ptPrim;
    Float_t diff  = pt - ptPrim;

    if     ( mcIndex==2 ) // Pi0
    {
      fhPairGeneratorsBkgMassMCPi0      [0][tag]->Fill(pt,  mass, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgMassMCPi0[genType][tag]->Fill(pt,  mass, GetEventWeight()*weightPt);

      fhPairGeneratorsBkgEPrimRecoRatioMCPi0[0][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgEPrimRecoDiffMCPi0 [0][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
      
      fhPairGeneratorsBkgEPrimRecoRatioMCPi0[genType][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgEPrimRecoDiffMCPi0 [genType][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
      
      if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] )
      {
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[0][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut [0][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
        
        fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[genType][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut [genType][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
      }
      
      if(IsHighMultiplicityAnalysisOn())
      {
        fhPairGeneratorsBkgCentMCPi0      [0][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgCentMCPi0[genType][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        if ( mass < fPi0MassWindow[1] && mass > fPi0MassWindow[0] )
        {
          fhPairGeneratorsBkgCentMCPi0MassCut      [0][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
          fhPairGeneratorsBkgCentMCPi0MassCut[genType][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        }
      }
    }
    else if( mcIndex==3 ) // Eta
    {
      fhPairGeneratorsBkgMassMCEta      [0][tag]->Fill(pt,  mass, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgMassMCEta[genType][tag]->Fill(pt,  mass, GetEventWeight()*weightPt);

      fhPairGeneratorsBkgEPrimRecoRatioMCEta[0][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgEPrimRecoDiffMCEta [0][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);    
      
      fhPairGeneratorsBkgEPrimRecoRatioMCEta[genType][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
      fhPairGeneratorsBkgEPrimRecoDiffMCEta [genType][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
      
      if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] )
      {
        fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[0][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut [0][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);    
        
        fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[genType][tag]->Fill(pt, ratio, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut [genType][tag]->Fill(pt,  diff, GetEventWeight()*weightPt);
      }
      
      if(IsHighMultiplicityAnalysisOn())
      {
        fhPairGeneratorsBkgCentMCEta      [0][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        fhPairGeneratorsBkgCentMCEta[genType][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        if ( mass < fEtaMassWindow[1] && mass > fEtaMassWindow[0] )
        {
          fhPairGeneratorsBkgCentMCEtaMassCut      [0][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
          fhPairGeneratorsBkgCentMCEtaMassCut[genType][tag]->Fill(pt,  cent, GetEventWeight()*weightPt);
        }
      }
    }
  } // do cluster overlaps from cocktails

  // carefull adding something else here, "returns" can affect
}

//__________________________________________
/// Main method. Process one event and extract photons from AOD branch
/// filled with AliAnaPhoton and fill histos with invariant mass.
/// Do also the mixed event filling and combination (if requested).
//__________________________________________
void AliAnaPi0::MakeAnalysisFillHistograms()
{
  // In case of simulated data, fill acceptance histograms
  if ( IsDataMC() && IsGeneratedParticlesAnalysisOn() )
  {
    FillAcceptanceHistograms();
  
    if ( fFillOnlyMCAcceptanceHisto ) return;
  }
  
//if (GetReader()->GetEventNumber()%10000 == 0)
// printf("--- Event %d ---\n",GetReader()->GetEventNumber());
  
  if ( !GetInputAODBranch() )
  {
    AliFatal(Form("No input aod photons in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
    return;
  }
  
  //
  // Init some variables
  //
  
  // Analysis within the same detector:
  TClonesArray* secondLoopInputData = GetInputAODBranch();
  
  Int_t nPhot      = GetInputAODBranch()->GetEntriesFast() ;
  Int_t nPhot2     = nPhot; // second loop
  Int_t minEntries = 2;
  Int_t last       = 1;     // last entry in first loop to be removed
  
  // Combine 2 detectors:
  if ( fPairWithOtherDetector )
  {
    // Input from CaloTrackCorr tasks
    secondLoopInputData = (TClonesArray *) GetReader()->GetAODBranchList()->FindObject(fOtherDetectorInputName);
    
    // In case input from PCM or other external source
    if(!secondLoopInputData) secondLoopInputData = (TClonesArray *) GetReader()->GetOutputEvent()->FindObject(fOtherDetectorInputName);
    if(!secondLoopInputData) secondLoopInputData = (TClonesArray *) GetReader()->GetInputEvent ()->FindObject(fOtherDetectorInputName);

    if(!secondLoopInputData)
    {
      AliFatal(Form("Input for other detector not found in branch %s!",fOtherDetectorInputName.Data()));
      return ; // coverity
    }
    
    nPhot2     = secondLoopInputData->GetEntriesFast() ; // add one since we remove one for the normal case
    minEntries = 1;
    last       = 0;
  }
  
  AliDebug(1,Form("Photon entries %d", nPhot));
  
  // If less than photon 2 (1) entries in the list, skip this event
  if ( nPhot < minEntries )
  {
    AliDebug(1,Form("nPhotons %d, cent bin %d continue to next event",nPhot, GetEventCentralityBin()));
    
    if ( GetNCentrBin() > 1 && IsHighMultiplicityAnalysisOn() )
        fhCentralityNoPair->Fill(GetEventCentralityBin(), GetEventWeight());
    
    return ;
  }
  else if ( fPairWithOtherDetector && nPhot2 < minEntries )
  {
    AliDebug(1,Form("nPhotons %d, cent bin %d continue to next event",nPhot, GetEventCentralityBin()));
    
    if ( GetNCentrBin() > 1 && IsHighMultiplicityAnalysisOn() )
      fhCentralityNoPair->Fill(GetEventCentralityBin(), GetEventWeight());
    
    return ;
  }
  
  // Init variables

  Int_t ncentr = GetNCentrBin();
  if ( GetNCentrBin() == 0 )
    ncentr = 1;
  Int_t curCentrBin = 0;
  if ( ncentr > 1 )
    curCentrBin = GetEventCentralityBin();

  Int_t module1         = -1;
  Int_t module2         = -1;
  Double_t vert[]       = {0.0, 0.0, 0.0} ; //vertex
  Int_t evtIndex1       = 0 ;
  Int_t currentEvtIndex = -1;
//Int_t curVzBin        = GetEventVzBin();
//Int_t curRPBin        = GetEventRPBin();
  Int_t eventbin        = GetEventMixBin();
    
  if ( eventbin < 0 || eventbin >= GetNCentrBin()*GetNZvertBin()*GetNRPBin() )
  {
    AliWarning(Form("Mix Bin not expected: cen bin %d (<%d), z bin %d (<%d), rp bin %d (<%d), "
                    "total bin %d, Event Centrality %d, z vertex %2.3f, Reaction Plane %2.3f",
                    GetEventCentralityBin(), ncentr,
                    GetEventVzBin(), GetNZvertBin(),
                    GetEventRPBin(), GetNRPBin(),
                    eventbin,GetEventCentrality(),
                    vert[2],GetEventPlaneAngle()));
    return;
  }
  
  // Get shower shape information of clusters
//  TObjArray *clusters = 0;
//  if     (GetCalorimeter()==kEMCAL) clusters = GetEMCALClusters();
//  else if(GetCalorimeter()==kPHOS ) clusters = GetPHOSClusters() ;
  
  //---------------------------------
  // First loop on photons/clusters
  //---------------------------------
  for(Int_t i1 = 0; i1 < nPhot-last; i1++)
  {
    AliCaloTrackParticle * p1 = (AliCaloTrackParticle*) (GetInputAODBranch()->At(i1)) ;

    // Select photons within a pT range
    if ( p1->Pt() < GetMinPt() || p1->Pt()  > GetMaxPt() ) continue ;
    
    //printf("AliAnaPi0::MakeAnalysisFillHistograms() : cluster1 id %d/%d\n",i1,nPhot-1);
    
    // get the event index in the mixed buffer where the photon comes from
    // in case of mixing with analysis frame, not own mixing
    evtIndex1 = GetEventIndex(p1, vert) ;
    if ( evtIndex1 == -1 )
      return ;
    if ( evtIndex1 == -2 )
      continue ;
    
    // Only effective in case of mixed event frame
    if ( TMath::Abs(vert[2]) > GetZvertexCut() ) continue ;   //vertex cut
    
    if ( evtIndex1 != currentEvtIndex )
    {
      // Fill event bin info
      if ( DoOwnMix() ) fhEventBin->Fill(eventbin, GetEventWeight()) ;
      
      if ( IsHighMultiplicityAnalysisOn() )
      {
        fhCentrality->Fill(curCentrBin, GetEventWeight());
        
        if ( GetEventPlane() )
          fhEventPlaneResolution->Fill(curCentrBin, TMath::Cos(2.*GetEventPlane()->GetQsubRes()), GetEventWeight());
      }
      
      currentEvtIndex = evtIndex1 ;
    }
    
    //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 1 Evt %d  Vertex : %f,%f,%f\n",evtIndex1, GetVertex(evtIndex1)[0] ,GetVertex(evtIndex1)[1],GetVertex(evtIndex1)[2]);
    
    //Get the momentum of this cluster
    fPhotonMom1.SetPxPyPzE(p1->Px(),p1->Py(),p1->Pz(),p1->E());
    
    //Get (Super)Module number of this cluster
    module1 =  p1->GetSModNumber();// GetModuleNumber(p1);
    
    //------------------------------------------
    // Recover original cluster
    // Declare variables for absid col-row identification of pair
    // Also fill col-row, eta-phi histograms depending on pt bin
    
    Int_t iclus1 = -1, iclus2 = -1 ;
    Float_t maxCellFraction1 = 0, maxCellFraction2 = 0;
    Int_t absIdMax1 = -1, absIdMax2 = -1;
    Int_t   icol1 = -1, icol2 = -1, icolAbs1 = -1, icolAbs2 = -1;
    Int_t   irow1 = -1, irow2 = -1, irowAbs1 = -1, irowAbs2 = -1;
    Int_t   iRCU1 = -1, iRCU2 = -1;

    if ( fMultiCutAnaAcc || fFillAngleHisto )
    {
      AliVCluster * cluster1 = FindCluster(GetEMCALClusters(),p1->GetCaloLabel(0),iclus1);
      if(!cluster1) AliWarning("Cluster1 not found!");
      
      absIdMax1 = GetCaloUtils()->GetMaxEnergyCell(GetEMCALCells(),cluster1,maxCellFraction1);
      
      GetModuleNumberCellIndexesAbsCaloMap(absIdMax1,GetCalorimeter(), icol1, irow1, iRCU1, icolAbs1, irowAbs1);
      
      if(fMultiCutAnaAcc)
      { 
        for(Int_t ipt = 0; ipt < fNPtCuts; ipt++)
        {
          if( p1->Pt() >   fPtCuts[ipt] && p1->Pt() < fPtCuts[ipt+1] )
          {
            fhPtBinClusterEtaPhi[ipt]->Fill(p1->Eta(),GetPhi(p1->Phi()),GetEventWeight()) ;
            
            fhPtBinClusterColRow[ipt]->Fill(icolAbs1,irowAbs1,GetEventWeight()) ;
          }
        }
      }
    }
    
    //---------------------------------
    // Second loop on photons/clusters
    //---------------------------------
    Int_t first = i1+1;
    if ( fPairWithOtherDetector ) first = 0;
    
    for(Int_t i2 = first; i2 < nPhot2; i2++)
    {
      //AliCaloTrackParticle * p2 = (AliCaloTrackParticle*) (GetInputAODBranch()->At(i2)) ;
      AliCaloTrackParticle * p2 = (AliCaloTrackParticle*) (secondLoopInputData->At(i2)) ;
      
      // Select photons within a pT range
      if ( p2->Pt() < GetMinPt() || p2->Pt()  > GetMaxPt() ) continue ;
      
      //printf("AliAnaPi0::MakeAnalysisFillHistograms() : cluster2 i %d/%d\n",i2,nPhot);
      
      //In case of mixing frame, check we are not in the same event as the first cluster
      Int_t evtIndex2 = GetEventIndex(p2, vert) ;
      if ( evtIndex2 == -1 )
        return ;
      if ( evtIndex2 == -2 )
        continue ;
      if (GetMixedEvent() && (evtIndex1 == evtIndex2))
        continue ;
      
//      //------------------------------------------
//      // Recover original cluster
//      Int_t iclus2 = -1;
//      AliVCluster * cluster2 = FindCluster(clusters,p2->GetCaloLabel(0),iclus2,iclus1+1);
//      // start new loop from iclus1+1 to gain some time
//      if(!cluster2) AliWarning("Cluster2 not found!");
//      
//      // Get the TOF,l0 and ncells from the clusters
//      Float_t tof1  = -1;
//      Float_t l01   = -1;
//      Int_t ncell1  = 0;
//      if(cluster1)
//      {
//        tof1   = cluster1->GetTOF()*1e9;
//        l01    = cluster1->GetM02();
//        ncell1 = cluster1->GetNCells();
//        //printf("cluster1: E %2.2f (%2.2f), l0 %2.2f, tof %2.2f\n",cluster1->E(),p1->E(),l01,tof1);
//      }
//      //else printf("cluster1 not available: calo label %d / %d, cluster ID %d\n",
//      //            p1->GetCaloLabel(0),(GetReader()->GetInputEvent())->GetNumberOfCaloClusters()-1,cluster1->GetID());
//      
//      Float_t tof2  = -1;
//      Float_t l02   = -1;
//      Int_t ncell2  = 0;
//      if(cluster2)
//      {
//        tof2   = cluster2->GetTOF()*1e9;
//        l02    = cluster2->GetM02();
//        ncell2 = cluster2->GetNCells();
//        //printf("cluster2: E %2.2f (%2.2f), l0 %2.2f, tof %2.2f\n",cluster2->E(),p2->E(),l02,tof2);
//      }
//      //else printf("cluster2 not available: calo label %d / %d, cluster ID %d\n",
//      //            p2->GetCaloLabel(0),(GetReader()->GetInputEvent())->GetNumberOfCaloClusters()-1,cluster2->GetID());
//      
//      if(cluster1 && cluster2)
//      {
//        Double_t t12diff = tof1-tof2;
//        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
//      }
      
      Float_t tof1   = p1->GetTime();
      Float_t l01    = p1->GetM02();
      Int_t   ncell1 = p1->GetNCells();
      //printf("cluster1: E %2.2f, l0 %2.2f, tof %2.2f\n",p1->E(),l01,tof1);
      
      Float_t tof2   = p2->GetTime();
      Float_t l02    = p2->GetM02();
      Int_t   ncell2 = p2->GetNCells();
      //printf("cluster2: E %2.2f, l0 %2.2f, tof %2.2f\n",p2->E(),l02,tof2);
      
      if ( !IsDataMC() )
      {
        Double_t t12diff = tof1-tof2;
        fhEPairDiffTime->Fill((fPhotonMom1 + fPhotonMom2).Pt(), t12diff, GetEventWeight());
        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
      }
      
      //------------------------------------------
      
      //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 2 Evt %d  Vertex : %f,%f,%f\n",evtIndex2, GetVertex(evtIndex2)[0] ,GetVertex(evtIndex2)[1],GetVertex(evtIndex2)[2]);
      
      // Get the momentum of this cluster
      fPhotonMom2.SetPxPyPzE(p2->Px(),p2->Py(),p2->Pz(),p2->E());
      
      // Get module number
      module2 = p2->GetSModNumber(); //GetModuleNumber(p2);
      
      //---------------------------------
      // Get pair kinematics
      //---------------------------------
      Double_t m    = (fPhotonMom1 + fPhotonMom2).M() ;
      Double_t pt   = (fPhotonMom1 + fPhotonMom2).Pt();
      Double_t deta = fPhotonMom1.Eta() - fPhotonMom2.Eta();
      Double_t dphi = fPhotonMom1.Phi() - fPhotonMom2.Phi();
      Double_t a    = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
      
      AliDebug(2,Form("E: fPhotonMom1 %f, fPhotonMom2 %f; Pair: pT %f, mass %f, a %f", p1->E(), p2->E(), (fPhotonMom1 + fPhotonMom2).E(),m,a));
      
      //--------------------------------
      // Opening angle selection
      //--------------------------------
      // Check if opening angle is too large or too small compared to what is expected
      Double_t angle   = fPhotonMom1.Angle(fPhotonMom2.Vect());
      if ( fUseAngleEDepCut && 
          !GetNeutralMesonSelection()->IsAngleInWindow((fPhotonMom1+fPhotonMom2).E(),angle+0.05) )
      {
        AliDebug(2,Form("Real pair angle %f (deg) not in E %f window",RadToDeg(angle), (fPhotonMom1+fPhotonMom2).E()));
        continue;
      }
      
      if ( fUseAngleCut && (angle < fAngleCut || angle > fAngleMaxCut) )
      {
        AliDebug(2,Form("Real pair cut %f < angle %f < cut %f (deg)",RadToDeg(fAngleCut), RadToDeg(angle), RadToDeg(fAngleMaxCut)));
        continue;
      }

      if ( fUseOneCellSeparation )
      {
        Bool_t separation = CheckSeparation(p1->GetCellAbsIdMax() ,p2->GetCellAbsIdMax());
        if ( !separation )
        {
          AliDebug(2,Form("Real pair one cell separation required and Yes/No %d", separation));
          continue;
        }
      }
      
      //-----------------------------------
      // In case of MC, get the ancestry and 
      // the weight depending on particle originating the pair if requested
      //-----------------------------------
      Int_t   ancPDG    = 0;
      Int_t   ancStatus = 0;
      Int_t   ancLabel  =-1;
      Float_t weightPt  = 1;
      if ( IsDataMC() )
      {
        ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(p1->GetLabel(), p2->GetLabel(),
                                                             GetMC(), ancPDG, ancStatus,fMCPrimMesonMom, fMCProdVertex);
        if ( ancLabel >= 0 )
        {
          TString genName;
          Int_t index   = GetReader()->GetCocktailGeneratorAndIndex(ancLabel, genName);
          //(GetMC())->GetCocktailGenerator(i,genName);
          
          weightPt = GetParticlePtWeight(fMCPrimMesonMom.Pt(), ancPDG, genName, index) ; 
        }
      }
      
      //-------------------------------------------------------------------------------------------------
      // Fill module dependent histograms, put a cut on assymmetry on the first available cut in the array
      //-------------------------------------------------------------------------------------------------
      if ( a < fAsymCuts[0] && fFillSMCombinations &&
           module1 >=0 && module1<fNModules        && 
           module2 >=0 && module2<fNModules           )
      {
        if ( !fPairWithOtherDetector )
        {
          if ( module1==module2 )
          {
            fhReMod[module1]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if(fFillAngleHisto) fhRealOpeningAnglePerSM[module1]->Fill(pt, angle, GetEventWeight()*weightPt);
          }
          else if ( GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL )
          {
            // Same sector
            Int_t isector1 = module1/2;
            Int_t isector2 = module2/2;
            if ( isector1==isector2 ) 
            {
              fhReSameSectorEMCALMod[isector1]->Fill(pt, m, GetEventWeight()) ;
            }
            // Same side
            else if ( TMath::Abs(isector2-isector1) == 1 )
            {
              Int_t iside1 = module1;
              Int_t iside2 = module2;
              // skip EMCal/DCal combination
              if(module1 > 11) iside1-=2; 
              if(module2 > 11) iside2-=2;
              
              if     ( module1 < module2 && module2-module1==2 ) 
                fhReSameSideEMCALMod[iside1]->Fill(pt, m, GetEventWeight());
              else if( module2 < module1 && module1-module2==2 ) 
                fhReSameSideEMCALMod[iside2]->Fill(pt, m, GetEventWeight());
            } // side
          } // EMCAL
          else if ( GetCalorimeter() == kPHOS )
          { // PHOS
            if((module1==0 && module2==1) || (module1==1 && module2==0)) fhReDiffPHOSMod[0]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==0 && module2==2) || (module1==2 && module2==0)) fhReDiffPHOSMod[1]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==1 && module2==2) || (module1==2 && module2==1)) fhReDiffPHOSMod[2]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==0 && module2==3) || (module1==3 && module2==0)) fhReDiffPHOSMod[3]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==1 && module2==3) || (module1==3 && module2==1)) fhReDiffPHOSMod[4]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==2 && module2==3) || (module1==3 && module2==2)) fhReDiffPHOSMod[5]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==0 && module2==4) || (module1==4 && module2==0)) fhReDiffPHOSMod[6]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==1 && module2==4) || (module1==4 && module2==1)) fhReDiffPHOSMod[7]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==2 && module2==4) || (module1==4 && module2==2)) fhReDiffPHOSMod[8]->Fill(pt, m, GetEventWeight()*weightPt) ;
            if((module1==3 && module2==4) || (module1==4 && module2==3)) fhReDiffPHOSMod[9]->Fill(pt, m, GetEventWeight()*weightPt) ;
          } // PHOS
        }
        else
        {
          Float_t phi1 = GetPhi(fPhotonMom1.Phi());
          Float_t phi2 = GetPhi(fPhotonMom2.Phi());
          Bool_t etaside = 0;
          if(   (p1->GetDetectorTag()==kEMCAL && fPhotonMom1.Eta() < 0) 
             || (p2->GetDetectorTag()==kEMCAL && fPhotonMom2.Eta() < 0)) etaside = 1;
          
          if      (    phi1 > DegToRad(260) && phi2 > DegToRad(260) && phi1 < DegToRad(280) && phi2 < DegToRad(280))  fhReSameSectorDCALPHOSMod[0+etaside]->Fill(pt, m, GetEventWeight()*weightPt);
          else if (    phi1 > DegToRad(280) && phi2 > DegToRad(280) && phi1 < DegToRad(300) && phi2 < DegToRad(300))  fhReSameSectorDCALPHOSMod[2+etaside]->Fill(pt, m, GetEventWeight()*weightPt);
          else if (    phi1 > DegToRad(300) && phi2 > DegToRad(300) && phi1 < DegToRad(320) && phi2 < DegToRad(320))  fhReSameSectorDCALPHOSMod[4+etaside]->Fill(pt, m, GetEventWeight()*weightPt);
          else if (   (phi1 > DegToRad(260) && phi2 > DegToRad(280) && phi1 < DegToRad(280) && phi2 < DegToRad(300)) 
                   || (phi1 > DegToRad(280) && phi2 > DegToRad(260) && phi1 < DegToRad(300) && phi2 < DegToRad(280))) fhReDiffSectorDCALPHOSMod[0+etaside]->Fill(pt, m, GetEventWeight()*weightPt);  
          else if (   (phi1 > DegToRad(280) && phi2 > DegToRad(300) && phi1 < DegToRad(300) && phi2 < DegToRad(320)) 
                   || (phi1 > DegToRad(300) && phi2 > DegToRad(280) && phi1 < DegToRad(320) && phi2 < DegToRad(300))) fhReDiffSectorDCALPHOSMod[2+etaside]->Fill(pt, m, GetEventWeight()*weightPt); 
          else if (   (phi1 > DegToRad(260) && phi2 > DegToRad(300) && phi1 < DegToRad(280) && phi2 < DegToRad(320)) 
                   || (phi1 > DegToRad(300) && phi2 > DegToRad(260) && phi1 < DegToRad(320) && phi2 < DegToRad(280))) fhReDiffSectorDCALPHOSMod[4+etaside]->Fill(pt, m, GetEventWeight()*weightPt); 
          else                                                                                                            fhReDiffSectorDCALPHOSMod[6+etaside]->Fill(pt, m, GetEventWeight()*weightPt);
        }
      } // Fill SM combinations
      
      // In case we want only pairs in same (super) module, check their origin.
      Bool_t ok = kTRUE;
      
      if ( fSameSM )
      {
        if ( !fPairWithOtherDetector )
        {
          if ( module1!=module2 ) ok=kFALSE;
        } 
        else // PHOS and DCal in same sector
        {
          Float_t phi1 = GetPhi(fPhotonMom1.Phi());
          Float_t phi2 = GetPhi(fPhotonMom2.Phi());
          ok=kFALSE;
          if      ( phi1 > DegToRad(260) && phi2 > DegToRad(260) && phi1 < DegToRad(280) && phi2 < DegToRad(280)) ok = kTRUE;
          else if ( phi1 > DegToRad(280) && phi2 > DegToRad(280) && phi1 < DegToRad(300) && phi2 < DegToRad(300)) ok = kTRUE;
          else if ( phi1 > DegToRad(300) && phi2 > DegToRad(300) && phi1 < DegToRad(320) && phi2 < DegToRad(320)) ok = kTRUE;
        }
      } // Pair only in same SM
      
      if ( !ok ) continue;
      
      //
      // Fill histograms with selected cluster pairs
      //
      
      // Check if one of the clusters comes from a conversion
      if ( fCheckConversion )
      {
        if     (p1->IsTagged() && p2->IsTagged()) fhReConv2->Fill(pt, m, GetEventWeight()*weightPt);
        else if(p1->IsTagged() || p2->IsTagged()) fhReConv ->Fill(pt, m, GetEventWeight()*weightPt);
      }
      
      // Fill shower shape cut histograms
      if ( fFillSSCombinations )
      {
        if     ( l01 > 0.01 && l01 < 0.4  &&
                 l02 > 0.01 && l02 < 0.4 )               fhReSS[0]->Fill(pt, m, GetEventWeight()*weightPt); // Tight
        else if( l01 > 0.4  && l02 > 0.4 )               fhReSS[1]->Fill(pt, m, GetEventWeight()*weightPt); // Loose
        else if( l01 > 0.01 && l01 < 0.4  && l02 > 0.4 ) fhReSS[2]->Fill(pt, m, GetEventWeight()*weightPt); // Both
        else if( l02 > 0.01 && l02 < 0.4  && l01 > 0.4 ) fhReSS[2]->Fill(pt, m, GetEventWeight()*weightPt); // Both
      }
      
      //
      // Main invariant mass histograms.
      // Fill histograms for different bad channel distance, centrality, assymmetry cut and pid bit
      //
      for(Int_t ipid=0; ipid<fNPIDBits; ipid++)
      {
        if ( (p1->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton)) && 
             (p2->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton))    )
        {
          for(Int_t iasym=0; iasym < fNAsymCuts; iasym++)
          {
            if ( a < fAsymCuts[iasym] )
            {
              if ( curCentrBin < 0 || curCentrBin >= GetNCentrBin() ) continue;
              
              Int_t index = ((curCentrBin*fNPIDBits)+ipid)*fNAsymCuts + iasym;
              //printf("index %d :(cen %d * nPID %d + ipid %d)*nasym %d + iasym %d - max index %d\n",index,curCentrBin,fNPIDBits,ipid,fNAsymCuts,iasym, curCentrBin*fNPIDBits*fNAsymCuts);
              
              if ( index < 0 || index >= ncentr*fNPIDBits*fNAsymCuts ) continue ;
              
              fhRe1     [index]->Fill(pt, m, GetEventWeight()*weightPt);
              
              if ( fMakeInvPtPlots ) fhReInvPt1[index]->Fill(pt, m, 1./pt * GetEventWeight()*weightPt) ;
              
              if ( fFillBadDistHisto )
              {
                if ( p1->DistToBad()>0 && p2->DistToBad()>0 )
                {
                  fhRe2     [index]->Fill(pt, m, GetEventWeight()*weightPt) ;
                  if ( fMakeInvPtPlots )
                    fhReInvPt2[index]->Fill(pt, m, 1./pt * GetEventWeight()*weightPt) ;
                  
                  if ( p1->DistToBad()>1 && p2->DistToBad()>1 )
                  {
                    fhRe3     [index]->Fill(pt, m, GetEventWeight()*weightPt) ;
                    if ( fMakeInvPtPlots )
                      fhReInvPt3[index]->Fill(pt, m, 1./pt * GetEventWeight()*weightPt) ;
                  }// bad 3
                }// bad2
              }// Fill bad dist histos
            }//assymetry cut
          }// asymmetry cut loop
        }// bad 1
      }// pid bit loop
      
      //
      // Fill histograms with opening angle
      if ( fFillAngleHisto )
      {
        fhRealOpeningAngle   ->Fill(pt, angle, GetEventWeight()*weightPt);
        fhRealCosOpeningAngle->Fill(pt, TMath::Cos(angle), GetEventWeight()*weightPt);
      }
      
      //
      // Fill histograms for different opening angle bins
      if ( fFillOpAngleCutHisto )
      {        
        Int_t angleBin = -1;
        for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
        {
          if ( angle >= fAngleCutBinsArray[ibin] && 
               angle <  fAngleCutBinsArray[ibin+1]  ) angleBin = ibin;
        }
        
        if ( angleBin >= 0 && angleBin < fNAngleCutBins )
        {
          Float_t e1   = fPhotonMom1.E();
          Float_t e2   = fPhotonMom2.E();

          Float_t t1   = tof1;
          Float_t t2   = tof2;
          
          Int_t nc1    = ncell1;
          Int_t nc2    = ncell2;          
          
          Float_t eta1 = fPhotonMom1.Eta(); 
          Float_t eta2 = fPhotonMom2.Eta(); 

          Float_t phi1 = GetPhi(fPhotonMom1.Phi());
          Float_t phi2 = GetPhi(fPhotonMom2.Phi());
          
          Int_t   mod1 = module1;
          Int_t   mod2 = module2;
          
          // Recover original cluster
          AliVCluster * cluster2 = FindCluster(GetEMCALClusters(),p2->GetCaloLabel(0),iclus2);
          if ( !cluster2 ) 
            AliWarning("Cluster2 not found!");

          absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(GetEMCALCells(),cluster2,maxCellFraction2);
          
          if ( e2 > e1 )
          {
            e1   = fPhotonMom2.E();
            e2   = fPhotonMom1.E();

            t1   = tof2;
            t2   = tof1;
            
            nc1  = ncell2;
            nc2  = ncell1;         
            
            eta1 = fPhotonMom2.Eta(); 
            eta2 = fPhotonMom1.Eta(); 
            
            phi1 = GetPhi(fPhotonMom2.Phi());
            phi2 = GetPhi(fPhotonMom1.Phi());
            
            mod1 = module2;
            mod2 = module1;
            
            Int_t tmp = absIdMax2;
            absIdMax2 = absIdMax1;
            absIdMax1 = tmp;
          }

          fhReOpAngleBinMinClusterEPerSM[angleBin]->Fill(e2,mod2,GetEventWeight()*weightPt) ; 
          fhReOpAngleBinMaxClusterEPerSM[angleBin]->Fill(e1,mod1,GetEventWeight()*weightPt) ; 
          
          fhReOpAngleBinMinClusterTimePerSM[angleBin]->Fill(t2,mod2,GetEventWeight()*weightPt) ; 
          fhReOpAngleBinMaxClusterTimePerSM[angleBin]->Fill(t1,mod1,GetEventWeight()*weightPt) ; 
          
          fhReOpAngleBinMinClusterNCellPerSM[angleBin]->Fill(nc2,mod2,GetEventWeight()*weightPt) ; 
          fhReOpAngleBinMaxClusterNCellPerSM[angleBin]->Fill(nc1,mod1,GetEventWeight()*weightPt) ; 

          fhReOpAngleBinPairClusterMass[angleBin]->Fill(pt,m,GetEventWeight()*weightPt) ;
          if ( mod2 == mod1 )  
            fhReOpAngleBinPairClusterMassPerSM[angleBin]->Fill(m,mod1,GetEventWeight()*weightPt) ;
                    
          if ( e1 > 0.01 )
            fhReOpAngleBinPairClusterRatioPerSM[angleBin]->Fill(e2/e1,mod1,GetEventWeight()*weightPt) ;  
          
          fhReOpAngleBinMinClusterEtaPhi[angleBin]->Fill(eta2,phi2,GetEventWeight()*weightPt) ;
          fhReOpAngleBinMaxClusterEtaPhi[angleBin]->Fill(eta1,phi1,GetEventWeight()*weightPt) ;
                    
          GetModuleNumberCellIndexesAbsCaloMap(absIdMax2,GetCalorimeter(), icol2, irow2, iRCU2, icolAbs2, irowAbs2);
          
          //fhReOpAngleBinPairClusterAbsIdMaxCell[angleBin]->Fill(absIdMax1,absIdMax2,GetEventWeight()*weightPt);

          fhReOpAngleBinMinClusterColRow[angleBin]->Fill(icolAbs2,irowAbs2,GetEventWeight()*weightPt) ;
          fhReOpAngleBinMaxClusterColRow[angleBin]->Fill(icolAbs1,irowAbs1,GetEventWeight()*weightPt) ;
        }
      } // fFillOpAngleHisto

      
      // Fill histograms with pair assymmetry
      if ( fFillAsymmetryHisto )
      {
        fhRePtAsym->Fill(pt, a, GetEventWeight()*weightPt);
        if ( m > fPi0MassWindow[0] && m < fPi0MassWindow[1] ) fhRePtAsymPi0->Fill(pt, a, GetEventWeight()*weightPt);
        if ( m > fEtaMassWindow[0] && m < fEtaMassWindow[1] ) fhRePtAsymEta->Fill(pt, a, GetEventWeight()*weightPt);
      }
      
      // Check cell time content in cluster
      if ( fFillSecondaryCellTiming )
      {
        if      ( p1->GetFiducialArea() == 0 && p2->GetFiducialArea() == 0 )
          fhReSecondaryCellInTimeWindow ->Fill(pt, m, GetEventWeight()*weightPt);
        
        else if ( p1->GetFiducialArea() != 0 && p2->GetFiducialArea() != 0 )
          fhReSecondaryCellOutTimeWindow->Fill(pt, m, GetEventWeight()*weightPt);
      }

      //---------
      // MC data
      //---------
      // Do some MC checks on the origin of the pair, is there any common ancestor and if there is one, who?
      if ( IsDataMC() )
      {
        if ( fCheckConversion )
        {
          if(GetMCAnalysisUtils()->CheckTagBit(p1->GetTag(),AliMCAnalysisUtils::kMCConversion) &&
             GetMCAnalysisUtils()->CheckTagBit(p2->GetTag(),AliMCAnalysisUtils::kMCConversion))
          {
            fhReMCFromConversion->Fill(pt, m, GetEventWeight()*weightPt);
          }
          else if(!GetMCAnalysisUtils()->CheckTagBit(p1->GetTag(),AliMCAnalysisUtils::kMCConversion) &&
                  !GetMCAnalysisUtils()->CheckTagBit(p2->GetTag(),AliMCAnalysisUtils::kMCConversion))
          {
            fhReMCFromNotConversion->Fill(pt, m, GetEventWeight()*weightPt);
          }
          else
          {
            fhReMCFromMixConversion->Fill(pt, m, GetEventWeight()*weightPt);
          }
        }
        
        if ( fFillOriginHisto )
        {
          FillMCVersusRecDataHistograms(ancLabel, ancPDG, ancStatus, weightPt,
                                        p1->GetCaloLabel(0), p2->GetCaloLabel(0),
                                        p1->GetTag(),p2->GetTag(),
                                        p1->Pt(), p2->Pt(),
                                        ncell1, ncell2, m, pt, a, deta, dphi, angle);
        }
      }
      
      //-----------------------
      // Multi cuts analysis
      //-----------------------
      if ( fMultiCutAna )
      {        
        // Several pt,ncell and asymmetry cuts
        for(Int_t ipt = 0; ipt < fNPtCuts; ipt++)
        {
          for(Int_t icell = 0; icell < fNCellNCuts; icell++)
          {
            for(Int_t iasym = 0; iasym < fNAsymCuts; iasym++)
            {
              Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
              if(p1->Pt() >   fPtCuts[ipt]      && p2->Pt()  > fPtCuts[ipt]    &&
                 p1->Pt() <   fPtCutsMax[ipt]   && p2->Pt()  < fPtCutsMax[ipt] &&
                 a        <   fAsymCuts[iasym]                                 &&
                 ncell1   >=  fCellNCuts[icell] && ncell2   >= fCellNCuts[icell])
              {
                fhRePtNCellAsymCuts[index]->Fill(pt, m, GetEventWeight()*weightPt) ;
                if(fFillAngleHisto)  fhRePtNCellAsymCutsOpAngle[index]->Fill(pt, angle, GetEventWeight()*weightPt) ;
                
                if(fFillSMCombinations && module1==module2)
                {
                  fhRePtNCellAsymCutsSM[module1][index]->Fill(pt, m, GetEventWeight()*weightPt) ;
                  if(fFillAngleHisto)  fhRePtNCellAsymCutsSMOpAngle[module1][index]->Fill(pt, angle, GetEventWeight()*weightPt) ;
                }
              }
            }// pid bit cut loop
          }// icell loop
        }// pt cut loop
      }// multiple cuts analysis
      
    }// second same event particle
  }// first cluster
  
  //-------------------------------------------------------------
  // Mixing
  //-------------------------------------------------------------
  if ( DoOwnMix() )
  {
    // Recover events in with same characteristics as the current event
    
    // Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
    if ( eventbin < 0 || eventbin >= GetNCentrBin()*GetNZvertBin()*GetNRPBin() ) 
      return ;
    
    TList * evMixList=fEventsList[eventbin] ;
    
    if ( !evMixList )
    {
      AliWarning(Form("Mix event list not available, bin %d",eventbin));
      return;
    }
    
    Int_t nMixed = evMixList->GetSize() ;
    for(Int_t ii=0; ii<nMixed; ii++)
    {
      TClonesArray* ev2= (TClonesArray*) (evMixList->At(ii));
      Int_t nPhot2=ev2->GetEntriesFast() ;
      Double_t m = -999;
      AliDebug(1,Form("Mixed event %d photon entries %d, centrality bin %d",ii, nPhot2, GetEventCentralityBin()));
      
      fhEventMixBin->Fill(eventbin, GetEventWeight()) ;
      
      //---------------------------------
      // First loop on photons/clusters
      //---------------------------------
      for(Int_t i1 = 0; i1 < nPhot; i1++)
      {
        AliCaloTrackParticle * p1 = (AliCaloTrackParticle*) (GetInputAODBranch()->At(i1)) ;
        
        // Select photons within a pT range
        if ( p1->Pt() < GetMinPt() || p1->Pt()  > GetMaxPt() ) continue ;
        
        // Not sure why this line is here
        //if(fSameSM && GetModuleNumber(p1)!=module1) continue;
        
        //Get kinematics of cluster and (super) module of this cluster
        fPhotonMom1.SetPxPyPzE(p1->Px(),p1->Py(),p1->Pz(),p1->E());
        module1 = GetModuleNumber(p1);
        
        //---------------------------------
        // Second loop on other mixed event photons/clusters
        //---------------------------------
        for(Int_t i2 = 0; i2 < nPhot2; i2++)
        {
          AliCaloTrackParticle * p2 = (AliCaloTrackParticle*) (ev2->At(i2)) ;
          
          // Select photons within a pT range
          if ( p2->Pt() < GetMinPt() || p2->Pt()  > GetMaxPt() ) continue ;
          
          // Get kinematics of second cluster and calculate those of the pair
          fPhotonMom2.SetPxPyPzE(p2->Px(),p2->Py(),p2->Pz(),p2->E());
          m           = (fPhotonMom1+fPhotonMom2).M() ;
          Double_t pt = (fPhotonMom1 + fPhotonMom2).Pt();
          Double_t a  = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
          
          // Check if opening angle is too large or too small compared to what is expected
          Double_t angle   = fPhotonMom1.Angle(fPhotonMom2.Vect());
          if ( fUseAngleEDepCut && 
              !GetNeutralMesonSelection()->IsAngleInWindow((fPhotonMom1+fPhotonMom2).E(),angle+0.05) )
          {
            AliDebug(2,Form("Mix pair angle %f (deg) not in E %f window",RadToDeg(angle), (fPhotonMom1+fPhotonMom2).E()));
            continue;
          }
          
          if ( fUseAngleCut && (angle < fAngleCut || angle > fAngleMaxCut) )
          {
            AliDebug(2,Form("Mix pair cut %f < angle %f < cut %f (deg)",RadToDeg(fAngleCut),RadToDeg(angle),RadToDeg(fAngleMaxCut)));
            continue;
          }

          if ( fUseOneCellSeparation )
          {
            Bool_t separation = CheckSeparation(p1->GetCellAbsIdMax() ,p2->GetCellAbsIdMax());
            if ( !separation )
            {
              AliDebug(2,Form("Mix pair one cell separation required and Yes/No %d", separation));
              continue;
            }
          }
          
          AliDebug(2,Form("Mixed Event: pT: fPhotonMom1 %2.2f, fPhotonMom2 %2.2f; Pair: pT %2.2f, mass %2.3f, a %2.3f",p1->Pt(), p2->Pt(), pt,m,a));
          
          // In case we want only pairs in same (super) module, check their origin.
          module2 = GetModuleNumber(p2);
                    
          //-------------------------------------------------------------------------------------------------
          // Fill module dependent histograms, put a cut on assymmetry on the first available cut in the array
          //-------------------------------------------------------------------------------------------------
          if ( a < fAsymCuts[0] && fFillSMCombinations &&
               module1 >=0 && module1<fNModules        && 
               module2 >=0 && module2<fNModules           )
          {
            if ( !fPairWithOtherDetector )
            {
              if ( module1==module2 )
              {
                fhMiMod[module1]->Fill(pt, m, GetEventWeight()) ;
                if(fFillAngleHisto) fhMixedOpeningAnglePerSM[module1]->Fill(pt, angle, GetEventWeight());
              }
              else if ( GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL )
              {
                // Same sector
                Int_t isector1 = module1/2;
                Int_t isector2 = module2/2;
                if ( isector1==isector2 ) 
                {
                  fhMiSameSectorEMCALMod[isector1]->Fill(pt, m, GetEventWeight()) ;
                }
                // Same side
                else if ( TMath::Abs(isector2-isector1) == 1 )
                {
                  Int_t iside1 = module1;
                  Int_t iside2 = module2;
                  // skip EMCal/DCal combination
                  if(module1 > 11) iside1-=2; 
                  if(module2 > 11) iside2-=2;
                  
                  if     ( module1 < module2 && module2-module1==2 ) 
                    fhMiSameSideEMCALMod[iside1]->Fill(pt, m, GetEventWeight());
                  else if( module2 < module1 && module1-module2==2 ) 
                    fhMiSameSideEMCALMod[iside2]->Fill(pt, m, GetEventWeight());
                }
              } // EMCAL
              else if ( GetCalorimeter() == kPHOS )
              { // PHOS
                if((module1==0 && module2==1) || (module1==1 && module2==0)) fhMiDiffPHOSMod[0]->Fill(pt, m, GetEventWeight()) ;
                if((module1==0 && module2==2) || (module1==2 && module2==0)) fhMiDiffPHOSMod[1]->Fill(pt, m, GetEventWeight()) ;
                if((module1==1 && module2==2) || (module1==2 && module2==1)) fhMiDiffPHOSMod[2]->Fill(pt, m, GetEventWeight()) ;
                if((module1==0 && module2==3) || (module1==3 && module2==0)) fhMiDiffPHOSMod[3]->Fill(pt, m, GetEventWeight()) ;
                if((module1==1 && module2==3) || (module1==3 && module2==1)) fhMiDiffPHOSMod[4]->Fill(pt, m, GetEventWeight()) ;
                if((module1==2 && module2==3) || (module1==3 && module2==2)) fhMiDiffPHOSMod[5]->Fill(pt, m, GetEventWeight()) ;
                if((module1==0 && module2==4) || (module1==4 && module2==0)) fhMiDiffPHOSMod[6]->Fill(pt, m, GetEventWeight()) ;
                if((module1==1 && module2==4) || (module1==4 && module2==1)) fhMiDiffPHOSMod[7]->Fill(pt, m, GetEventWeight()) ;
                if((module1==2 && module2==4) || (module1==4 && module2==2)) fhMiDiffPHOSMod[8]->Fill(pt, m, GetEventWeight()) ;
                if((module1==3 && module2==4) || (module1==4 && module2==3)) fhMiDiffPHOSMod[9]->Fill(pt, m, GetEventWeight()) ;
              } // PHOS
            }
            else
            {
              Float_t phi1 = GetPhi(fPhotonMom1.Phi());
              Float_t phi2 = GetPhi(fPhotonMom2.Phi());
              Bool_t etaside = 0;
              if (   (p1->GetDetectorTag()==kEMCAL && fPhotonMom1.Eta() < 0) 
                  || (p2->GetDetectorTag()==kEMCAL && fPhotonMom2.Eta() < 0)) etaside = 1;
              
              if      (    phi1 > DegToRad(260) && phi2 > DegToRad(260) && phi1 < DegToRad(280) && phi2 < DegToRad(280))  fhMiSameSectorDCALPHOSMod[0+etaside]->Fill(pt, m, GetEventWeight());
              else if (    phi1 > DegToRad(280) && phi2 > DegToRad(280) && phi1 < DegToRad(300) && phi2 < DegToRad(300))  fhMiSameSectorDCALPHOSMod[2+etaside]->Fill(pt, m, GetEventWeight());
              else if (    phi1 > DegToRad(300) && phi2 > DegToRad(300) && phi1 < DegToRad(320) && phi2 < DegToRad(320))  fhMiSameSectorDCALPHOSMod[4+etaside]->Fill(pt, m, GetEventWeight());
              else if (   (phi1 > DegToRad(260) && phi2 > DegToRad(280) && phi1 < DegToRad(280) && phi2 < DegToRad(300)) 
                       || (phi1 > DegToRad(280) && phi2 > DegToRad(260) && phi1 < DegToRad(300) && phi2 < DegToRad(280))) fhMiDiffSectorDCALPHOSMod[0+etaside]->Fill(pt, m, GetEventWeight());  
              else if (   (phi1 > DegToRad(280) && phi2 > DegToRad(300) && phi1 < DegToRad(300) && phi2 < DegToRad(320)) 
                       || (phi1 > DegToRad(300) && phi2 > DegToRad(280) && phi1 < DegToRad(320) && phi2 < DegToRad(300))) fhMiDiffSectorDCALPHOSMod[2+etaside]->Fill(pt, m, GetEventWeight()); 
              else if (   (phi1 > DegToRad(260) && phi2 > DegToRad(300) && phi1 < DegToRad(280) && phi2 < DegToRad(320)) 
                       || (phi1 > DegToRad(300) && phi2 > DegToRad(260) && phi1 < DegToRad(320) && phi2 < DegToRad(280))) fhMiDiffSectorDCALPHOSMod[4+etaside]->Fill(pt, m, GetEventWeight()); 
              else                                                                                                            fhMiDiffSectorDCALPHOSMod[6+etaside]->Fill(pt, m, GetEventWeight());
            }            
          } //  different SM combinations
          
          Bool_t ok = kTRUE;          
          if ( fSameSM )
          {
            if ( !fPairWithOtherDetector )
            {
              if ( module1!=module2 ) ok=kFALSE;
            } 
            else // PHOS and DCal in same sector
            {
              Float_t phi1 = GetPhi(fPhotonMom1.Phi());
              Float_t phi2 = GetPhi(fPhotonMom2.Phi());
              ok=kFALSE;
              if      ( phi1 > DegToRad(260) && phi2 > DegToRad(260) && phi1 < DegToRad(280) && phi2 < DegToRad(280)) ok = kTRUE;
              else if ( phi1 > DegToRad(280) && phi2 > DegToRad(280) && phi1 < DegToRad(300) && phi2 < DegToRad(300)) ok = kTRUE;
              else if ( phi1 > DegToRad(300) && phi2 > DegToRad(300) && phi1 < DegToRad(320) && phi2 < DegToRad(320)) ok = kTRUE;
            }
          } // Pair only in same SM
          
          if ( !ok ) continue ;
          
          //
          // Do the final histograms with the selected clusters
          //
          
          // Check if one of the clusters comes from a conversion
          if ( fCheckConversion )
          {
            if     (p1->IsTagged() && p2->IsTagged()) fhMiConv2->Fill(pt, m, GetEventWeight());
            else if(p1->IsTagged() || p2->IsTagged()) fhMiConv ->Fill(pt, m, GetEventWeight());
          }
          
          //
          // Main invariant mass histograms
          // Fill histograms for different bad channel distance, centrality, assymmetry cut and pid bit
          //
          for(Int_t ipid=0; ipid<fNPIDBits; ipid++)
          {
            if ( (p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton)) )
            {
              for(Int_t iasym=0; iasym < fNAsymCuts; iasym++)
              {
                if ( a < fAsymCuts[iasym] )
                {
                  if ( curCentrBin < 0 || curCentrBin >= GetNCentrBin() ) continue;

                  Int_t index = ((curCentrBin*fNPIDBits)+ipid)*fNAsymCuts + iasym;
                  
                  if ( index < 0 || index >= ncentr*fNPIDBits*fNAsymCuts ) continue ;
                  
                  fhMi1[index]->Fill(pt, m, GetEventWeight()) ;
                  if(fMakeInvPtPlots)fhMiInvPt1[index]->Fill(pt, m, 1./pt * GetEventWeight()) ;
                  
                  if ( fFillBadDistHisto )
                  {
                    if ( p1->DistToBad()>0 && p2->DistToBad()>0 )
                    {
                      fhMi2[index]->Fill(pt, m, GetEventWeight()) ;
                      if ( fMakeInvPtPlots )
                        fhMiInvPt2[index]->Fill(pt, m, 1./pt * GetEventWeight()) ;
                      
                      if ( p1->DistToBad()>1 && p2->DistToBad()>1 )
                      {
                        fhMi3[index]->Fill(pt, m, GetEventWeight()) ;
                        if ( fMakeInvPtPlots )
                          fhMiInvPt3[index]->Fill(pt, m, 1./pt * GetEventWeight()) ;
                      }
                    }
                  }// Fill bad dist histo
                  
                }//Asymmetry cut
              }// Asymmetry loop
            }//PID cut
          }// PID loop 

          //-----------------------
          // Multi cuts analysis
          //-----------------------
          Int_t  ncell1 = p1->GetNCells();
          Int_t  ncell2 = p1->GetNCells();
          
          if ( fMultiCutAna )
          {
            // Several pt,ncell and asymmetry cuts
            for(Int_t ipt=0; ipt<fNPtCuts; ipt++)
            {
              for(Int_t icell=0; icell<fNCellNCuts; icell++)
              {
                for(Int_t iasym=0; iasym<fNAsymCuts; iasym++)
                {
                  Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                  
                  if(p1->Pt() >   fPtCuts[ipt]      && p2->Pt() > fPtCuts[ipt]      &&
                     p1->Pt() <   fPtCutsMax[ipt]   && p2->Pt() < fPtCutsMax[ipt]   &&
                     a        <   fAsymCuts[iasym]                                  &&
                     ncell1   >=  fCellNCuts[icell] && ncell2   >= fCellNCuts[icell] 
                     )
                  {
                    //printf("MI ipt %d, iasym%d, icell %d, index %d \n",ipt, iasym, icell, index);
                    //printf("\t %p, %p\n",fhMiPtNCellAsymCuts[index],fhMiPtNCellAsymCutsOpAngle[index]);
                    
                    fhMiPtNCellAsymCuts[index]->Fill(pt, m, GetEventWeight()) ;
                    if ( fFillAngleHisto )
                      fhMiPtNCellAsymCutsOpAngle[index]->Fill(pt, angle, GetEventWeight()) ;
                    
                    //printf("ipt %d, icell%d, iasym %d, name %s\n",ipt, icell, iasym,  fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym]->GetName());
                  }
                }// pid bit cut loop
              }// icell loop
            }// pt cut loop
          } // Multi cut ana
          
          //
          // Fill histograms with opening angle
          if ( fFillAngleHisto )
          {
            fhMixedOpeningAngle   ->Fill(pt, angle, GetEventWeight());
            fhMixedCosOpeningAngle->Fill(pt, TMath::Cos(angle), GetEventWeight());
          }          
          
          //
          // Fill histograms for different opening angle bins
          if ( fFillOpAngleCutHisto )
          {
            Int_t angleBin = -1;
            for(Int_t ibin = 0; ibin < fNAngleCutBins; ibin++)
            {
              if ( angle >= fAngleCutBinsArray[ibin] && 
                   angle <  fAngleCutBinsArray[ibin+1]) angleBin = ibin;
            }
            
            if ( angleBin >= 0 && angleBin < fNAngleCutBins )
            {
              Float_t e1   = fPhotonMom1.E();
              Float_t e2   = fPhotonMom2.E();
              
              Float_t t1   = p1->GetTime();
              Float_t t2   = p2->GetTime();
              
              Int_t nc1    = ncell1;
              Int_t nc2    = ncell2;
              
              Float_t eta1 = fPhotonMom1.Eta(); 
              Float_t eta2 = fPhotonMom2.Eta(); 
              
              Float_t phi1 = GetPhi(fPhotonMom1.Phi());
              Float_t phi2 = GetPhi(fPhotonMom2.Phi());
              
              Int_t   mod1 = module1;
              Int_t   mod2 = module2;
              
              //              // Recover original cluster
              //              Int_t iclus1 = -1, iclus2 = -1 ;
              //              AliVCluster * cluster1 = FindCluster(GetEMCALClusters(),p1->GetCaloLabel(0),iclus1);
              //              AliVCluster * cluster2 = FindCluster(GetEMCALClusters(),p2->GetCaloLabel(0),iclus2);
              //              
              //              Float_t maxCellFraction1 = 0, maxCellFraction2 = 0;
              //              Int_t absIdMax1 = GetCaloUtils()->GetMaxEnergyCell(GetEMCALCells(),cluster1,maxCellFraction1);
              //              Int_t absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(GetEMCALCells(),cluster2,maxCellFraction2);
              
              if ( e2 > e1 )
              {
                e1   = fPhotonMom2.E();
                e2   = fPhotonMom1.E();
                
                t1   = p2->GetTime();
                t2   = p1->GetTime();
                
                nc1  = ncell2;
                nc2  = ncell1;
                
                eta1 = fPhotonMom2.Eta(); 
                eta2 = fPhotonMom1.Eta(); 
                
                phi1 = GetPhi(fPhotonMom2.Phi());
                phi2 = GetPhi(fPhotonMom1.Phi());
                
                mod1 = module2;
                mod2 = module1;
                
                //                Int_t tmp = absIdMax2;
                //                absIdMax2 = absIdMax1;
                //                absIdMax1 = tmp;
              }
              
              fhMiOpAngleBinMinClusterEPerSM[angleBin]->Fill(e2,mod2,GetEventWeight()) ; 
              fhMiOpAngleBinMaxClusterEPerSM[angleBin]->Fill(e1,mod1,GetEventWeight()) ; 
              
              fhMiOpAngleBinMinClusterTimePerSM[angleBin]->Fill(t2,mod2,GetEventWeight()) ; 
              fhMiOpAngleBinMaxClusterTimePerSM[angleBin]->Fill(t1,mod1,GetEventWeight()) ; 
              
              fhMiOpAngleBinMinClusterNCellPerSM[angleBin]->Fill(nc2,mod2,GetEventWeight()) ; 
              fhMiOpAngleBinMaxClusterNCellPerSM[angleBin]->Fill(nc1,mod1,GetEventWeight()) ; 
              
              fhMiOpAngleBinPairClusterMass[angleBin]->Fill(pt,m,GetEventWeight()) ;
              if ( mod2 == mod1 )
                fhMiOpAngleBinPairClusterMassPerSM[angleBin]->Fill(m,mod1,GetEventWeight()) ;
              
              if ( e1 > 0.01 )
                fhMiOpAngleBinPairClusterRatioPerSM[angleBin]->Fill(e2/e1,mod1,GetEventWeight()) ;  
              
              fhMiOpAngleBinMinClusterEtaPhi[angleBin]->Fill(eta2,phi2,GetEventWeight()) ;
              fhMiOpAngleBinMaxClusterEtaPhi[angleBin]->Fill(eta1,phi1,GetEventWeight()) ;
              
              //              Int_t   icol1 = -1, icol2 = -1, icolAbs1 = -1, icolAbs2 = -1;
              //              Int_t   irow1 = -1, irow2 = -1, irowAbs1 = -1, irowAbs2 = -1;
              //              Int_t   iRCU1 = -1, iRCU2 = -1;
              //              GetModuleNumberCellIndexesAbsCaloMap(absIdMax1,GetCalorimeter(), icol1, irow1, iRCU1, icolAbs1, irowAbs1);
              //              GetModuleNumberCellIndexesAbsCaloMap(absIdMax2,GetCalorimeter(), icol2, irow2, iRCU1, icolAbs2, irowAbs2);
              //              
              //              fhMiOpAngleBinPairClusterAbsIdMaxCell[angleBin]->Fill(absIdMax1,absIdMax2,GetEventWeight());
              //
              //              fhMiColRowClusterMinOpAngleBin[angleBin]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
              //              fhMiOpAngleBinMaxClusterColRow[angleBin]->Fill(icolAbs1,irowAbs1,GetEventWeight()) ;
            }
          }
          
          // Fill histograms with pair assymmetry
          if ( fFillAsymmetryHisto )
          {
            fhMiPtAsym->Fill(pt, a, GetEventWeight());
            if ( m > fPi0MassWindow[0] && m < fPi0MassWindow[1] ) fhMiPtAsymPi0->Fill(pt, a, GetEventWeight());
            if ( m > fEtaMassWindow[0] && m < fEtaMassWindow[1] ) fhMiPtAsymEta->Fill(pt, a, GetEventWeight());
          }
          
          // Check cell time content in cluster
          if ( fFillSecondaryCellTiming )
          {
            if      ( p1->GetFiducialArea() == 0 && p2->GetFiducialArea() == 0 )
              fhMiSecondaryCellInTimeWindow ->Fill(pt, m, GetEventWeight());
            
            else if ( p1->GetFiducialArea() != 0 && p2->GetFiducialArea() != 0 )
              fhMiSecondaryCellOutTimeWindow->Fill(pt, m, GetEventWeight());
          }
                  
        }// second cluster loop
      }//first cluster loop
    }//loop on mixed events
    
    //--------------------------------------------------------
    // Add the current event to the list of events for mixing
    //--------------------------------------------------------
    
    //TClonesArray *currentEvent = new TClonesArray(*GetInputAODBranch());
    TClonesArray *currentEvent = new TClonesArray(*secondLoopInputData);
    
    // Add current event to buffer and Remove redundant events
    if ( currentEvent->GetEntriesFast() > 0 )
    {
      evMixList->AddFirst(currentEvent) ;
      currentEvent = 0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
      if ( evMixList->GetSize() >= GetNMaxEvMix() )
      {
        TClonesArray * tmp = (TClonesArray*) (evMixList->Last()) ;
        evMixList->RemoveLast() ;
        delete tmp ;
      }
    }
    else
    { // empty event
      delete currentEvent ;
      currentEvent=0 ;
    }
  }// DoOwnMix
  
  AliDebug(1,"End fill histograms");
}

//________________________________________________________________________
/// It retieves the event index and checks the vertex
///  * in the mixed buffer returns -2 if vertex NOK
///  * for normal events   returns 0 if vertex OK and -1 if vertex NOK
//________________________________________________________________________
Int_t AliAnaPi0::GetEventIndex(AliCaloTrackParticle * part, Double_t * vert)
{
  Int_t evtIndex = -1 ;
    
  if(GetReader()->GetDataType()!=AliCaloTrackReader::kMC)
  {
    if (GetMixedEvent())
    {
      evtIndex = GetMixedEvent()->EventIndexForCaloCluster(part->GetCaloLabel(0)) ;
      GetVertex(vert,evtIndex);
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -2 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
    }
    else
    {
      // Single event
      GetVertex(vert);
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -1 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
      else
        evtIndex = 0 ;
    }
  } // No MC reader
  else
  {
    evtIndex = 0;
    vert[0] = 0. ; 
    vert[1] = 0. ; 
    vert[2] = 0. ; 
  }
  
  return evtIndex ; 
}

//________________________________________________________________________
/// Checks whether two photon cluster maxima have at least one cell separation
/// TRUE - clusters are separated
/// FALSE - are NOT separated
//________________________________________________________________________
Bool_t AliAnaPi0::CheckSeparation(Int_t absID1,Int_t absID2)
{
  if(absID1<0 || absID2<0) return kFALSE;
  Bool_t neighbours = GetCaloUtils()->AreNeighbours(AliFiducialCut::kEMCAL, absID1, absID2);
  //neighbours = kTRUE -> no separation
  
  return (!neighbours);
}
