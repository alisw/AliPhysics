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

#include "TMath.h"

#include "AliHistogramRanges.h"

/// \cond CLASSIMP
ClassImp(AliHistogramRanges) ;
/// \endcond

//_______________________________________________
/// Default constructor. Initialize parameters.
//_______________________________________________
AliHistogramRanges::AliHistogramRanges() : 
TObject(), 
fHistoPtBins(0),              fHistoPtMax(0.),             fHistoPtMin(0.),         fHistoPtArr(0),
fHistoPhiBins(0),             fHistoPhiMax(0.),            fHistoPhiMin(0.),        fHistoPhiArr(0),
fHistoEtaBins(0),             fHistoEtaMax(0.),            fHistoEtaMin(0.),        fHistoEtaArr(0),
fHistoDeltaPhiBins(0),        fHistoDeltaPhiMax(0.),       fHistoDeltaPhiMin(0.),   fHistoDeltaPhiArr(0),
fHistoDeltaEtaBins(0),        fHistoDeltaEtaMax(0.),       fHistoDeltaEtaMin(0.),   fHistoDeltaEtaArr(0),
fHistoMassBins(0),            fHistoMassMax(0.),           fHistoMassMin(0.),       fHistoMassArr(0),
fHistoAsymBins(0),            fHistoAsymMax(0.),           fHistoAsymMin(0.),       fHistoAsymArr(0),
fHistoV0SBins(0),             fHistoV0SMax(0),             fHistoV0SMin(0),         fHistoV0SArr(0),
fHistoV0MBins(0),             fHistoV0MMax(0),             fHistoV0MMin(0),         fHistoV0MArr(0),
fHistoTrMBins(0),             fHistoTrMMax(0),             fHistoTrMMin(0),         fHistoTrMArr(0),
fHistoFinePtBins(1000),       fHistoFinePtMax(5.),         fHistoFinePtMin(0.),     fHistoFinePtArr(0),
fHistoWidePtBins(20),         fHistoWidePtMax(200.),       fHistoWidePtMin(0.),     fHistoWidePtArr(0),
fHistoEOverPBins(100),        fHistoEOverPMax(100.),       fHistoEOverPMin(0.),     fHistoEOverPArr(0),
fHistoNSigmaBins(100),        fHistoNSigmaMax(100.),       fHistoNSigmaMin(-100),   fHistoNSigmaArr(0),
fHistodEdxBins(100),          fHistodEdxMax(100.),         fHistodEdxMin(0.),       fHistodEdxArr(0),
fHistodRBins(100),            fHistodRMax(100.),           fHistodRMin(0.),         fHistodRArr(0),
fHistoTimeBins(100),          fHistoTimeMax(100.),         fHistoTimeMin(0.),       fHistoTimeArr(0),
fHistoNClusCellBins(100),     fHistoNClusCellMax(100),     fHistoNClusCellMin(0),   fHistoNClusCellArr(0),
fHistoNCellsBins(100),        fHistoNCellsMax(100),        fHistoNCellsMin(0),      fHistoNCellsArr(0),
fHistoNClustersBins(100),     fHistoNClustersMax(100),     fHistoNClustersMin(0),   fHistoNClustersArr(0),
fHistoRatioBins(100),         fHistoRatioMax (2.),         fHistoRatioMin (0.),     fHistoRatioArr (0),
fHistoRatio1Bins(100),        fHistoRatio1Max(1.),         fHistoRatio1Min(0.),     fHistoRatio1Arr(0),
fHistoEDiffBins(100),         fHistoEDiffMax(5.),          fHistoEDiffMin(-5.),     fHistoEDiffArr(0),
fHistoHBPBins(100),           fHistoHBPMax(100.),          fHistoHBPMin(0.),        fHistoHBPArr(0),
fHistoVertexDistBins(100),    fHistoVertexDistMax(100.),   fHistoVertexDistMin(0.), fHistoVertexDistArr(0),
fHistoRBins(100),             fHistoRMax(1000),            fHistoRMin(-1000),       fHistoRArr(0),
fHistoXBins(100),             fHistoXMax(1000),            fHistoXMin(-1000),       fHistoXArr(0),
fHistoYBins(100),             fHistoYMax(1000),            fHistoYMin(-1000),       fHistoYArr(0),
fHistoZBins(100),             fHistoZMax(1000),            fHistoZMin(-1000),       fHistoZArr(0),
fHistoSSBins(0),              fHistoSSMax(0),              fHistoSSMin(0),          fHistoSSArr(0),
fHistoDiffTimeBins(0),        fHistoDiffTimeMax(0),        fHistoDiffTimeMin(0),    fHistoDiffTimeArr(0),
fHistoTrackResidualEtaBins(0),fHistoTrackResidualEtaMax(0),fHistoTrackResidualEtaMin(0), fHistoTrackResidualEtaArr(0),
fHistoTrackResidualPhiBins(0),fHistoTrackResidualPhiMax(0),fHistoTrackResidualPhiMin(0), fHistoTrackResidualPhiArr(0),
fHistoNPtSumBins(0),          fHistoPtSumMax(0.),          fHistoPtSumMin(0.),       fHistoPtSumArr(0),
fHistoNPtSumSubBins(0),       fHistoPtSumSubMax(0.),       fHistoPtSumSubMin(0.),    fHistoPtSumSubArr(0),
fHistoNPtInConeBins(0),       fHistoPtInConeMax(0.),       fHistoPtInConeMin(0.),    fHistoPtInConeArr(0),
fHistoOpAngleBins(200),       fHistoOpAngleMax(0.7),       fHistoOpAngleMin(0.),     fHistoOpAngleArr(0),
fHistoCenBins(100),           fHistoCenMax(100),           fHistoCenMin(0),          fHistoCenArr(0),
fHistoNLMBins(20),            fHistoNLMMax(20),            fHistoNLMMin(0),          fHistoNLMArr(0),
fHistoNoverlapBins(20),       fHistoNoverlapMax(20),       fHistoNoverlapMin(0),     fHistoNoverlapArr(0),
fHistoExoticityBins(200),     fHistoExoticityMax(1),       fHistoExoticityMin(-1),   fHistoExoticityArr(0),
fHistoSpherocityBins(100),    fHistoSpherocityMax(1),      fHistoSpherocityMin(0),   fHistoSpherocityArr(0)
{
  InitParameters();
}

//_______________________________________
/// Initialize histogram parameters
//_______________________________________
void AliHistogramRanges::InitParameters()
{   
  fHistoPtBins           = 240 ;  fHistoPtMax           = 120   ; fHistoPtMin           = 0.  ;
  fHistoPhiBins          = 120 ;  fHistoPhiMax          = TMath::TwoPi(); fHistoPhiMin  = 0.  ;
  fHistoEtaBins          = 100 ;  fHistoEtaMax          = 1     ; fHistoEtaMin          = -1  ;
  fHistoDeltaPhiBins     = 163 ;  fHistoDeltaPhiMax     = 4.8   ; fHistoDeltaPhiMin     = -1.7;
  fHistoDeltaEtaBins     =  40 ;  fHistoDeltaEtaMax     = 2     ; fHistoDeltaEtaMin     = -2  ;
  fHistoMassBins         = 200 ;  fHistoMassMax         = 1.    ; fHistoMassMin         = 0.  ;
  fHistoAsymBins         = 10  ;  fHistoAsymMax         = 1.    ; fHistoAsymMin         = 0.  ;
  fHistoV0SBins          = 100 ;  fHistoV0SMax          = 10000 ; fHistoV0SMin          = 0   ;
  fHistoV0MBins          = 100;   fHistoV0MMax          = 10000 ; fHistoV0MMin          = 0   ;
  fHistoTrMBins          = 200 ;  fHistoTrMMax          = 200   ; fHistoTrMMin          = 0   ;
  fHistoEOverPBins       = 100 ;  fHistoEOverPMax       = 10.   ; fHistoEOverPMin       = 0.  ;
  fHistoNSigmaBins       = 300 ;  fHistoNSigmaMax       = 15.   ; fHistoNSigmaMin       =-15.  ;
  fHistodEdxBins         = 200 ;  fHistodEdxMax         = 400.  ; fHistodEdxMin         = 0.  ;  
  fHistodRBins           = 300 ;  fHistodRMax           = 3.15  ; fHistodRMin           = 0.  ;
  fHistoTimeBins         = 200 ;  fHistoTimeMax         = 200   ; fHistoTimeMin         =-200.;//ns
  fHistoNClusCellBins    = 50  ;  fHistoNClusCellMax    = 50    ; fHistoNClusCellMin    = 0   ;
  fHistoNCellsBins       = 300 ;  fHistoNCellsMax       = 300   ; fHistoNCellsMin       = 0   ;
  fHistoNClustersBins    = 50  ;  fHistoNClustersMax    = 50    ; fHistoNClustersMin    = 0   ;
  fHistoRatioBins        = 200 ;  fHistoRatioMax        = 2.    ; fHistoRatioMin        = 0.  ;
  fHistoRatio1Bins       = 100 ;  fHistoRatio1Max       = 1.    ; fHistoRatio1Min       = 0.  ;
  fHistoEDiffBins        = 200 ;  fHistoEDiffMax        = 10    ; fHistoEDiffMin        =-10. ;
  fHistoHBPBins          = 200 ;  fHistoHBPMax          = 10    ; fHistoHBPMin          = 0.  ;
  fHistoVertexDistBins   = 100 ;  fHistoVertexDistMax   = 500.  ; fHistoVertexDistMin   = 0.  ;
  fHistoRBins            = 100 ;  fHistoRMax            = 500   ; fHistoRMin            = -500;//cm
  fHistoXBins            = 100 ;  fHistoXMax            = 500   ; fHistoXMin            = -500;//cm
  fHistoYBins            = 100 ;  fHistoYMax            = 500   ; fHistoYMin            = -500;//cm
  fHistoZBins            = 100 ;  fHistoZMax            = 600   ; fHistoZMin            = -500;//cm
  fHistoSSBins           = 500 ;  fHistoSSMax           = 5     ; fHistoSSMin           = 0   ;  
  fHistoDiffTimeBins     = 400 ;  fHistoDiffTimeMax     = 400   ; fHistoDiffTimeMin     = -400;// ns
  fHistoNPtSumBins       = 100 ;  fHistoPtSumMax        = 100   ; fHistoPtSumMin        = 0.  ;
  fHistoNPtSumSubBins    = 200 ;  fHistoPtSumMax        = 100   ; fHistoPtSumMin        =-100. ;
  fHistoNPtInConeBins    = 100 ;  fHistoPtInConeMax     = 50    ; fHistoPtInConeMin     = 0.  ;
  fHistoOpAngleBins      = 200 ;  fHistoOpAngleMax      = 0.7   ; fHistoOpAngleMin      = 0.  ;
  fHistoCenBins          = 200 ;  fHistoCenMax          = 100   ; fHistoCenMin          = 0.  ;
  fHistoNLMBins          = 20  ;  fHistoNLMMax          = 20    ; fHistoNLMMin          = 0.  ;
  fHistoNoverlapBins     = 10  ;  fHistoNoverlapMax     = 10    ; fHistoNoverlapMin     = 0.  ;
  fHistoExoticityBins    = 200 ;  fHistoExoticityMax    = 1.    ; fHistoExoticityMin    = -1. ;
  fHistoSpherocityBins   = 100 ;  fHistoSpherocityMax   = 1.    ; fHistoSpherocityMin   = 0.  ;

  fHistoTrackResidualEtaBins = 100 ; fHistoTrackResidualEtaMax = 0.15 ; fHistoTrackResidualEtaMin = -0.15;
  fHistoTrackResidualPhiBins = 100 ; fHistoTrackResidualPhiMax = 0.15 ; fHistoTrackResidualPhiMin = -0.15;
}

//________________________________________________________
// Print histogram parameters
//________________________________________________________
void AliHistogramRanges::Print(const Option_t * /*opt*/) const
{  	
  printf("Histograms: %3.1f < pT  < %3.1f,  Nbin = %d\n"             , fHistoPtMin,          fHistoPtMax,          fHistoPtBins);
  printf("Histograms: %3.1f < pT  < %3.1f, fine Nbin = %d\n"         , fHistoFinePtMin,      fHistoFinePtMax,      fHistoFinePtBins);
  printf("Histograms: %3.1f < pT  < %3.1f, wide Nbin = %d\n"         , fHistoWidePtMin,      fHistoWidePtMax,      fHistoWidePtBins);
  printf("Histograms: %3.1f < phi < %3.1f, Nbin = %d\n"              , fHistoPhiMin,         fHistoPhiMax,         fHistoPhiBins);
  printf("Histograms: %3.1f < eta < %3.1f, Nbin = %d\n"              , fHistoEtaMin,         fHistoEtaMax,         fHistoEtaBins);
  printf("Histograms: %3.1f < delta phi < %3.1f, Nbin = %d\n"        , fHistoDeltaPhiMin,    fHistoDeltaPhiMax,    fHistoDeltaPhiBins);
  printf("Histograms: %3.1f < delta eta < %3.1f, Nbin = %d\n"        , fHistoDeltaEtaMin,    fHistoDeltaEtaMax,    fHistoDeltaEtaBins);
  printf("Histograms: %3.1f < mass < %3.1f, Nbin = %d\n"             , fHistoMassMin,        fHistoMassMax,        fHistoMassBins);
  printf("Histograms: %3.1f < asymmetry < %3.1f, Nbin = %d\n"        , fHistoAsymMin,        fHistoAsymMax,        fHistoAsymBins);
  printf("Histograms: %d < V0 Signal < %d, Nbin = %d\n"              , fHistoV0SMin,         fHistoV0SMax,         fHistoV0SBins);
  printf("Histograms: %d < V0 Mult < %d, Nbin = %d\n"                , fHistoV0MMin,         fHistoV0MMax,         fHistoV0MBins);
  printf("Histograms: %d < Track Mult < %d, Nbin = %d\n"             , fHistoTrMMin,         fHistoTrMMax,         fHistoTrMBins);
  printf("Histograms: %3.1f < E/p  < %3.1f, Nbin = %d\n"             , fHistoEOverPMin,      fHistoEOverPMax,      fHistoEOverPBins);
  printf("Histograms: %3.1f < nSigma  < %3.1f, Nbin = %d\n"          , fHistoNSigmaMin,      fHistoNSigmaMax,      fHistoNSigmaBins);
  printf("Histograms: %3.1f < dEdx < %3.1f, Nbin = %d\n"             , fHistodEdxMin,        fHistodEdxMax,        fHistodEdxBins);
  printf("Histograms: %3.1f < dR (track cluster)< %3.1f, Nbin = %d\n", fHistodRMin,          fHistodRMax,          fHistodRBins);
  printf("Histograms: %3.1f < R=sqrt{x^2+y^2}   < %3.1f, Nbin = %d\n", fHistoRMin,           fHistoRMax,           fHistoRBins);
  printf("Histograms: %3.1f < X    < %3.1f, Nbin = %d\n"             , fHistoXMin,           fHistoXMax,           fHistoXBins);
  printf("Histograms: %3.1f < Y    < %3.1f, Nbin = %d\n"             , fHistoYMin,           fHistoYMax,           fHistoYBins);
  printf("Histograms: %3.1f < Z    < %3.1f, Nbin = %d\n"             , fHistoZMin,           fHistoZMax,           fHistoZBins);
  printf("Histograms: %g < Time < %g, Nbin = %d\n"                   , fHistoTimeMin,        fHistoTimeMax,        fHistoTimeBins);
  printf("Histograms: %d < N cells per cluster    < %d, Nbin = %d\n" , fHistoNClusCellMin,   fHistoNClusCellMax,   fHistoNClusCellBins);
  printf("Histograms: %d < N cells   < %d, Nbin = %d\n"              , fHistoNCellsMin,      fHistoNCellsMax,      fHistoNCellsBins);
  printf("Histograms: %d < N clusters   < %d, Nbin = %d\n"           , fHistoNClustersMin,   fHistoNClustersMax,   fHistoNClustersBins);
  printf("Histograms: %3.1f < Ratio < %3.1f, Nbin = %d\n"            , fHistoRatioMin,       fHistoRatioMax,       fHistoRatioBins);
  printf("Histograms: %3.1f < Ratio1 < %3.1f, Nbin = %d\n"           , fHistoRatio1Min,      fHistoRatio1Max,      fHistoRatio1Bins);
  printf("Histograms: %3.1f < En diff< %3.1f, Nbin = %d\n"           , fHistoEDiffMin,       fHistoEDiffMax,       fHistoEDiffBins);
  printf("Histograms: %3.1f < HBP< %3.1f, Nbin = %d\n"               , fHistoHBPMin,         fHistoHBPMax,         fHistoHBPBins);
  printf("Histograms: %3.1f < Vertex Distance < %3.1f,   Nbin = %d\n", fHistoVertexDistMin,  fHistoVertexDistMax,  fHistoVertexDistBins);
  printf("Histograms: %3.1f < pT sum < %3.1f,  Nbin = %d\n"          , fHistoPtSumMin,       fHistoPtSumMax,       fHistoNPtSumBins   );
  printf("Histograms: %3.1f < pT sum-UE < %3.1f,  Nbin = %d\n"       , fHistoPtSumSubMin,    fHistoPtSumSubMax,    fHistoNPtSumSubBins);
  printf("Histograms: %3.1f < pT in cone < %3.1f, Nbin = %d\n"       , fHistoPtInConeMin,    fHistoPtInConeMax,    fHistoNPtInConeBins);
  printf("Histograms: %3.1f < Op. angle < %3.1f, Nbin = %d\n"        , fHistoOpAngleMin,     fHistoOpAngleMax,     fHistoOpAngleBins);
  printf("Histograms: %2.2f < Residual Eta(Z) < %2.2f,   Nbin = %d\n", fHistoTrackResidualEtaMin, fHistoTrackResidualEtaMax,fHistoTrackResidualEtaBins);
  printf("Histograms: %2.2f < Residual Phi(R,X) < %2.2f, Nbin = %d\n", fHistoTrackResidualPhiMin, fHistoTrackResidualPhiMax,fHistoTrackResidualPhiBins);
  printf("Histograms: %2.2f < Centrality < %2.2f, Nbin = %d\n"       , fHistoCenMin,         fHistoCenMax,         fHistoCenBins);
  printf("Histograms: %d < NLM < %d, Nbin = %d\n"                    , fHistoNLMMin,         fHistoNLMMax,         fHistoNLMBins);
  printf("Histograms: %d < Noverlaps < %d, Nbin = %d\n"              , fHistoNoverlapMin,    fHistoNoverlapMax,    fHistoNoverlapBins);
  printf("Histograms: %2.2f < Exoticity < %2.2f, Nbin = %d\n"        , fHistoExoticityMin,   fHistoExoticityMax,   fHistoExoticityBins);
  printf("Histograms: %2.2f < Spherocity < %2.2f, Nbin = %d\n"       , fHistoSpherocityMin,  fHistoSpherocityMax,  fHistoSpherocityBins);

  printf("    \n") ;
} 



