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

//_________________________________________________________________________
// Class containing histogram settings:
//    - Number of bins
//    - Min and max of range
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)
//

#include "TMath.h"

#include "AliHistogramRanges.h"

ClassImp(AliHistogramRanges)


//_______________________________________________
AliHistogramRanges::AliHistogramRanges() : 
TObject(), 
fHistoPtBins(0),              fHistoPtMax(0.),             fHistoPtMin(0.),
fHistoPhiBins(0),             fHistoPhiMax(0.),            fHistoPhiMin(0.),
fHistoEtaBins(0),             fHistoEtaMax(0.),            fHistoEtaMin(0.),
fHistoDeltaPhiBins(0),        fHistoDeltaPhiMax(0.),       fHistoDeltaPhiMin(0.),
fHistoDeltaEtaBins(0),        fHistoDeltaEtaMax(0.),       fHistoDeltaEtaMin(0.),
fHistoMassBins(0),            fHistoMassMax(0.),           fHistoMassMin(0.),
fHistoAsymBins(0),            fHistoAsymMax(0.),           fHistoAsymMin(0.),
fHistoV0SBins(0),             fHistoV0SMax(0),             fHistoV0SMin(0),
fHistoV0MBins(0),             fHistoV0MMax(0),             fHistoV0MMin(0),
fHistoTrMBins(0),             fHistoTrMMax(0),             fHistoTrMMin(0),
fHistoFinePtBins(1000),       fHistoFinePtMax(5.),         fHistoFinePtMin(0.),
fHistoPOverEBins(100),        fHistoPOverEMax(100.),       fHistoPOverEMin(0.),
fHistodEdxBins(100),          fHistodEdxMax(100.),         fHistodEdxMin(0.),
fHistodRBins(100),            fHistodRMax(100.),           fHistodRMin(0.),
fHistoTimeBins(100),          fHistoTimeMax(100.),         fHistoTimeMin(0.),
fHistoNClusCellBins(100),     fHistoNClusCellMax(100),     fHistoNClusCellMin(0),
fHistoNCellsBins(100),        fHistoNCellsMax(100),        fHistoNCellsMin(0),
fHistoNClustersBins(100),     fHistoNClustersMax(100),     fHistoNClustersMin(0),
fHistoRatioBins(100),         fHistoRatioMax(100.),        fHistoRatioMin(0.),
fHistoVertexDistBins(100),    fHistoVertexDistMax(100.),   fHistoVertexDistMin(0.),
fHistoRBins(100),             fHistoRMax(1000),            fHistoRMin(-1000),
fHistoXBins(100),             fHistoXMax(1000),            fHistoXMin(-1000),
fHistoYBins(100),             fHistoYMax(1000),            fHistoYMin(-1000),
fHistoZBins(100),             fHistoZMax(1000),            fHistoZMin(-1000),
fHistoSSBins(0),              fHistoSSMax(0),              fHistoSSMin(0),
fHistoDiffTimeBins(0),        fHistoDiffTimeMax(0),        fHistoDiffTimeMin(0),
fHistoTrackResidualEtaBins(0),fHistoTrackResidualEtaMax(0),fHistoTrackResidualEtaMin(0),
fHistoTrackResidualPhiBins(0),fHistoTrackResidualPhiMax(0),fHistoTrackResidualPhiMin(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
  
}

//_______________________________________
void AliHistogramRanges::InitParameters()
{ 
  //Initialize histogram parameters
  
  fHistoPtBins           = 240 ;  fHistoPtMax           = 120   ; fHistoPtMin           = 0.  ;
  fHistoPhiBins          = 120 ;  fHistoPhiMax          = TMath::TwoPi(); fHistoPhiMin  = 0.  ;
  fHistoEtaBins          = 100 ;  fHistoEtaMax          =  1    ; fHistoEtaMin          = -1  ;
  fHistoDeltaPhiBins     = 140 ;  fHistoDeltaPhiMax     =  3* TMath::TwoPi()/2 ; fHistoDeltaPhiMin  = -  TMath::TwoPi()/2 ;
  fHistoDeltaEtaBins     = 200 ;  fHistoDeltaEtaMax     = TMath::TwoPi()/2 ;     fHistoDeltaEtaMin  = - TMath::TwoPi()/2 ;
  fHistoMassBins         = 200 ;  fHistoMassMax         = 1.    ; fHistoMassMin         = 0.  ;
  fHistoAsymBins         = 10  ;  fHistoAsymMax         = 1.    ; fHistoAsymMin         = 0.  ;
  fHistoV0SBins          = 100 ;  fHistoV0SMax          = 10000 ; fHistoV0SMin          = 0   ;
  fHistoV0MBins          = 100;   fHistoV0MMax          = 10000 ; fHistoV0MMin          = 0   ;
  fHistoTrMBins          = 100 ;  fHistoTrMMax          = 2000  ; fHistoTrMMin          = 0   ;
  fHistoPOverEBins       = 100 ;  fHistoPOverEMax       = 10.   ; fHistoPOverEMin       = 0.  ;
  fHistodEdxBins         = 200 ;  fHistodEdxMax         = 400.  ; fHistodEdxMin         = 0.  ;  
  fHistodRBins           = 300 ;  fHistodRMax           = 3.15  ; fHistodRMin           = 0.  ;
  fHistoTimeBins         = 1000;  fHistoTimeMax         = 1.e3  ; fHistoTimeMin         = 0.  ;//ns
  fHistoNClusCellBins    = 200 ;  fHistoNClusCellMax    = 200   ; fHistoNClusCellMin    = 0   ;
  fHistoNCellsBins       = 300 ;  fHistoNCellsMax       = 300   ; fHistoNCellsMin       = 0   ;
  fHistoNClustersBins    = 50  ;  fHistoNClustersMax    = 50    ; fHistoNClustersMin    = 0   ;
  fHistoRatioBins        = 200 ;  fHistoRatioMax        = 2     ; fHistoRatioMin        = 0.  ;
  fHistoVertexDistBins   = 100 ;  fHistoVertexDistMax   = 500.  ; fHistoVertexDistMin   = 0.  ;
  fHistoRBins            = 100 ;  fHistoRMax            = 500   ; fHistoRMin            = -500;//cm
  fHistoXBins            = 100 ;  fHistoXMax            = 500   ; fHistoXMin            = -500;//cm
  fHistoYBins            = 100 ;  fHistoYMax            = 500   ; fHistoYMin            = -500;//cm
  fHistoZBins            = 100 ;  fHistoZMax            = 600   ; fHistoZMin            = -500;//cm
  fHistoSSBins           = 500 ;  fHistoSSMax           = 5     ; fHistoSSMin           = 0   ;  
  fHistoDiffTimeBins     = 800 ;  fHistoDiffTimeMax     = 400   ; fHistoDiffTimeMin     = -400;// ns
  
  fHistoTrackResidualEtaBins = 100 ; fHistoTrackResidualEtaMax = 0.15 ; fHistoTrackResidualEtaMin = -0.15;
  fHistoTrackResidualPhiBins = 100 ; fHistoTrackResidualPhiMax = 0.15 ; fHistoTrackResidualPhiMin = -0.15;

}

//________________________________________________________
void AliHistogramRanges::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  	
  printf("Histograms: %3.1f < pT  < %3.1f,  Nbin = %d\n"             , fHistoPtMin,          fHistoPtMax,          fHistoPtBins);
  printf("Histograms: %3.1f < phi < %3.1f, Nbin = %d\n"              , fHistoPhiMin,         fHistoPhiMax,         fHistoPhiBins);
  printf("Histograms: %3.1f < eta < %3.1f, Nbin = %d\n"              , fHistoEtaMin,         fHistoEtaMax,         fHistoEtaBins);
  printf("Histograms: %3.1f < delta phi < %3.1f, Nbin = %d\n"        , fHistoDeltaPhiMin,    fHistoDeltaPhiMax,    fHistoDeltaPhiBins);
  printf("Histograms: %3.1f < delta eta < %3.1f, Nbin = %d\n"        , fHistoDeltaEtaMin,    fHistoDeltaEtaMax,    fHistoDeltaEtaBins);
  printf("Histograms: %3.1f < mass < %3.1f, Nbin = %d\n"             , fHistoMassMin,        fHistoMassMax,        fHistoMassBins);
  printf("Histograms: %3.1f < asymmetry < %3.1f, Nbin = %d\n"        , fHistoAsymMin,        fHistoAsymMax,        fHistoAsymBins);
  printf("Histograms: %d < V0 Signal < %d, Nbin = %d\n"              , fHistoV0SMin,         fHistoV0SMax,         fHistoV0SBins);
  printf("Histograms: %d < V0 Mult < %d, Nbin = %d\n"                , fHistoV0MMin,         fHistoV0MMax,         fHistoV0MBins);
  printf("Histograms: %d < Track Mult < %d, Nbin = %d\n"             , fHistoTrMMin,         fHistoTrMMax,         fHistoTrMBins);
  printf("Histograms: %3.1f < p/E  < %3.1f, Nbin = %d\n"             , fHistoPOverEMin,      fHistoPOverEMax,      fHistoPOverEBins);
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
  printf("Histograms: %3.1f < Ratio< %3.1f, Nbin = %d\n"             , fHistoRatioMin,       fHistoRatioMax,       fHistoRatioBins);
  printf("Histograms: %3.1f < Vertex Distance < %3.1f,   Nbin = %d\n", fHistoVertexDistMin,  fHistoVertexDistMax,  fHistoVertexDistBins);
  printf("Histograms: %2.2f < Residual Eta(Z) < %2.2f,   Nbin = %d\n", fHistoTrackResidualEtaMin, fHistoTrackResidualEtaMax,fHistoTrackResidualEtaBins);
  printf("Histograms: %2.2f < Residual Phi(R,X) < %2.2f, Nbin = %d\n", fHistoTrackResidualPhiMin, fHistoTrackResidualPhiMax,fHistoTrackResidualPhiBins);

  printf("    \n") ;
  
} 



