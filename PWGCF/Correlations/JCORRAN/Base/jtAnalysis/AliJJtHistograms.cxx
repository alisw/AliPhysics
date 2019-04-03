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

// Container class for histograms needed in the analysis.

#include "AliJJtHistograms.h"
#include "../AliJCard.h"
#include "../AliJBaseTrack.h"
#include "../AliJPhoton.h"
#include "../AliJTrack.h"
#include <TGrid.h>
#include <TPRegexp.h>

//______________________________________________________________________________
AliJJtHistograms::AliJJtHistograms(AliJCard* cardP) :
  AliJHistogramInterface(cardP),
  fhDphiDetaKlong(),
  fhDphiDetaXlong(),
  fhDphiDetaPta(),
  fhDphiDetaTrackMerge(),
  fhDetaNearMixAcceptance(),
  fhDphiDetaBgKlongEta(),
  fhDphiDetaBgKlongR(),
  fhDphiDetaBgKlongPhi(),
  fhDphiDetaBgXlongEta(),
  fhDphiDetaBgXlongR(),
  fhDphiDetaBgXlongPhi(),
  fhDphiDetaBgPtaEta(),
  fhDphiDetaBgPtaR(),
  fhDphiDetaBgPtaPhi(),
  fhBgAssocKlongEta(),
  fhBgAssocKlongR(),
  fhBgAssocKlongPhi(),
  fhBgAssocXlongEta(),
  fhBgAssocXlongR(),
  fhBgAssocXlongPhi(),
  fhBgAssocPtaEta(),
  fhBgAssocPtaR(),
  fhBgAssocPtaPhi(),
  fhInvariantMassXe(),
  fhInvariantMassKlong(),
  fhInvariantMassPta(),
  fhInvariantMassXeLikeSign(),
  fhInvariantMassKlongLikeSign(),
  fhInvariantMassPtaLikeSign(),
  fhInvariantMassXeUnlikeSign(),
  fhInvariantMassKlongUnlikeSign(),
  fhInvariantMassPtaUnlikeSign(),
  fhIphiTrigg(),
  fhIetaTrigg(),
  fhIphiAssoc(),
  fhIetaAssoc(),
  fhTriggPtBin(),
  fhAssocPtBin(),
  fhJT(),
  fhJTBg(),
  fhJTBgR(),
  fhJTBgPhi(),
  fhJTLikeSign(),
  fhJTBgLikeSign(),
  fhJTBgRLikeSign(),
  fhJTBgPhiLikeSign(),
  fhJTUnlikeSign(),
  fhJTBgUnlikeSign(),
  fhJTBgRUnlikeSign(),
  fhJTBgPhiUnlikeSign(),
  fhJTKlong(),
  fhJTKlongBg(),
  fhJTKlongBgR(),
  fhJTKlongBgPhi(),
  fhJTKlongLikeSign(),
  fhJTKlongBgLikeSign(),
  fhJTKlongBgRLikeSign(),
  fhJTKlongBgPhiLikeSign(),
  fhJTKlongUnlikeSign(),
  fhJTKlongBgUnlikeSign(),
  fhJTKlongBgRUnlikeSign(),
  fhJTKlongBgPhiUnlikeSign(),
  fhJTPta(),
  fhJTPtaBg(),
  fhJTPtaBgR(),
  fhJTPtaBgPhi(),
  fhJTPtaLikeSign(),
  fhJTPtaBgLikeSign(),
  fhJTPtaBgRLikeSign(),
  fhJTPtaBgPhiLikeSign(),
  fhJTPtaUnlikeSign(),
  fhJTPtaBgUnlikeSign(),
  fhJTPtaBgRUnlikeSign(),
  fhJTPtaBgPhiUnlikeSign(),
  fHmgInclusive(NULL),
  fhIetaTriggFromFile(),
  fhIetaAssocFromFile(),
  fhIphiTriggFromFile(),
  fhIphiAssocFromFile(),
  fhLPpt(),
  fhChargedPt(),
  fhChargedPtNoCorr(),
  fhTrackingEfficiency(),
  fhChargedEta(),
  fhLPeta(),
  fhChargedMult(),
  fhCentr(),
  fhiCentr(),
  fhZVert(),
  fhAcceptanceTraditional(),
  fhAcceptanceTraditional2D(),
  fhAcceptance3DNearSide(),
  fhDphiDetaTrackMergeCorrection(),
  fmaxEtaRange(0),
  fenable2DHistos(false),
  fEnableAcceptanceQAHistos(false)
{   // constructor

    fmaxEtaRange = fCard->Get("EtaRange");
  
}

//______________________________________________________________________________
AliJJtHistograms::AliJJtHistograms(const AliJJtHistograms& obj) :
  AliJHistogramInterface(obj),
  fhDphiDetaKlong(obj.fhDphiDetaKlong),
  fhDphiDetaXlong(obj.fhDphiDetaXlong),
  fhDphiDetaPta(obj.fhDphiDetaPta),
  fhDphiDetaTrackMerge(obj.fhDphiDetaTrackMerge),
  fhDetaNearMixAcceptance(obj.fhDetaNearMixAcceptance),
  fhDphiDetaBgKlongEta(obj.fhDphiDetaBgKlongEta),
  fhDphiDetaBgKlongR(obj.fhDphiDetaBgKlongR),
  fhDphiDetaBgKlongPhi(obj.fhDphiDetaBgKlongPhi),
  fhDphiDetaBgXlongEta(obj.fhDphiDetaBgXlongEta),
  fhDphiDetaBgXlongR(obj.fhDphiDetaBgXlongR),
  fhDphiDetaBgXlongPhi(obj.fhDphiDetaBgXlongPhi),
  fhDphiDetaBgPtaEta(obj.fhDphiDetaBgPtaEta),
  fhDphiDetaBgPtaR(obj.fhDphiDetaBgPtaR),
  fhDphiDetaBgPtaPhi(obj.fhDphiDetaBgPtaPhi),
  fhBgAssocKlongEta(obj.fhBgAssocKlongEta),
  fhBgAssocKlongR(obj.fhBgAssocKlongR),
  fhBgAssocKlongPhi(obj.fhBgAssocKlongPhi),
  fhBgAssocXlongEta(obj.fhBgAssocXlongEta),
  fhBgAssocXlongR(obj.fhBgAssocXlongR),
  fhBgAssocXlongPhi(obj.fhBgAssocXlongPhi),
  fhBgAssocPtaEta(obj.fhBgAssocPtaEta),
  fhBgAssocPtaR(obj.fhBgAssocPtaR),
  fhBgAssocPtaPhi(obj.fhBgAssocPtaPhi),
  fhInvariantMassXe(obj.fhInvariantMassXe),
  fhInvariantMassKlong(obj.fhInvariantMassKlong),
  fhInvariantMassPta(obj.fhInvariantMassPta),
  fhInvariantMassXeLikeSign(obj.fhInvariantMassXeLikeSign),
  fhInvariantMassKlongLikeSign(obj.fhInvariantMassKlongLikeSign),
  fhInvariantMassPtaLikeSign(obj.fhInvariantMassPtaLikeSign),
  fhInvariantMassXeUnlikeSign(obj.fhInvariantMassXeUnlikeSign),
  fhInvariantMassKlongUnlikeSign(obj.fhInvariantMassKlongUnlikeSign),
  fhInvariantMassPtaUnlikeSign(obj.fhInvariantMassPtaUnlikeSign),
  fhIphiTrigg(obj.fhIphiTrigg),
  fhIetaTrigg(obj.fhIetaTrigg),
  fhIphiAssoc(obj.fhIphiAssoc),
  fhIetaAssoc(obj.fhIetaAssoc),
  fhTriggPtBin(obj.fhTriggPtBin),
  fhAssocPtBin(obj.fhAssocPtBin),
  fhJT(obj.fhJT),
  fhJTBg(obj.fhJTBg),
  fhJTBgR(obj.fhJTBgR),
  fhJTBgPhi(obj.fhJTBgPhi),
  fhJTLikeSign(obj.fhJTLikeSign),
  fhJTBgLikeSign(obj.fhJTBgLikeSign),
  fhJTBgRLikeSign(obj.fhJTBgRLikeSign),
  fhJTBgPhiLikeSign(obj.fhJTBgPhiLikeSign),
  fhJTUnlikeSign(obj.fhJTUnlikeSign),
  fhJTBgUnlikeSign(obj.fhJTBgUnlikeSign),
  fhJTBgRUnlikeSign(obj.fhJTBgRUnlikeSign),
  fhJTBgPhiUnlikeSign(obj.fhJTBgPhiUnlikeSign),
  fhJTKlong(obj.fhJTKlong),
  fhJTKlongBg(obj.fhJTKlongBg),
  fhJTKlongBgR(obj.fhJTKlongBgR),
  fhJTKlongBgPhi(obj.fhJTKlongBgPhi),
  fhJTKlongLikeSign(obj.fhJTKlongLikeSign),
  fhJTKlongBgLikeSign(obj.fhJTKlongBgLikeSign),
  fhJTKlongBgRLikeSign(obj.fhJTKlongBgRLikeSign),
  fhJTKlongBgPhiLikeSign(obj.fhJTKlongBgPhiLikeSign),
  fhJTKlongUnlikeSign(obj.fhJTKlongUnlikeSign),
  fhJTKlongBgUnlikeSign(obj.fhJTKlongBgUnlikeSign),
  fhJTKlongBgRUnlikeSign(obj.fhJTKlongBgRUnlikeSign),
  fhJTKlongBgPhiUnlikeSign(obj.fhJTKlongBgPhiUnlikeSign),
  fhJTPta(obj.fhJTPta),
  fhJTPtaBg(obj.fhJTPtaBg),
  fhJTPtaBgR(obj.fhJTPtaBgR),
  fhJTPtaBgPhi(obj.fhJTPtaBgPhi),
  fhJTPtaLikeSign(obj.fhJTPtaLikeSign),
  fhJTPtaBgLikeSign(obj.fhJTPtaBgLikeSign),
  fhJTPtaBgRLikeSign(obj.fhJTPtaBgRLikeSign),
  fhJTPtaBgPhiLikeSign(obj.fhJTPtaBgPhiLikeSign),
  fhJTPtaUnlikeSign(obj.fhJTPtaUnlikeSign),
  fhJTPtaBgUnlikeSign(obj.fhJTPtaBgUnlikeSign),
  fhJTPtaBgRUnlikeSign(obj.fhJTPtaBgRUnlikeSign),
  fhJTPtaBgPhiUnlikeSign(obj.fhJTPtaBgPhiUnlikeSign),
  fHmgInclusive(obj.fHmgInclusive),
  fhIetaTriggFromFile(obj.fhIetaTriggFromFile),
  fhIetaAssocFromFile(obj.fhIetaAssocFromFile),
  fhIphiTriggFromFile(obj.fhIphiTriggFromFile),
  fhIphiAssocFromFile(obj.fhIphiAssocFromFile),
  fhLPpt(obj.fhLPpt),
  fhChargedPt(obj.fhChargedPt),
  fhChargedPtNoCorr(obj.fhChargedPtNoCorr),
  fhTrackingEfficiency(obj.fhTrackingEfficiency),
  fhChargedEta(obj.fhChargedEta),
  fhLPeta(obj.fhLPeta),
  fhChargedMult(obj.fhChargedMult),
  fhCentr(obj.fhCentr),
  fhiCentr(obj.fhiCentr),
  fhZVert(obj.fhZVert),
  fhAcceptanceTraditional(obj.fhAcceptanceTraditional),
  fhAcceptanceTraditional2D(obj.fhAcceptanceTraditional2D),
  fhAcceptance3DNearSide(obj.fhAcceptance3DNearSide),
  fhDphiDetaTrackMergeCorrection(obj.fhDphiDetaTrackMergeCorrection),
  fmaxEtaRange(obj.fmaxEtaRange),
  fenable2DHistos(obj.fenable2DHistos),
  fEnableAcceptanceQAHistos(obj.fEnableAcceptanceQAHistos)
{
    // copy constructor
    JUNUSED(obj);
}

//______________________________________________________________________________
AliJJtHistograms& AliJJtHistograms::operator=(const AliJJtHistograms& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}

//______________________________________________________________________________
AliJJtHistograms::~AliJJtHistograms() {
  // destructor
	delete fHMG;
	delete fHmgInclusive;
}

//______________________________________________________________________________
void AliJJtHistograms::CreateCorrelationHistograms()
{
  // Create all the histograms needed in correlation analysis
  fHMG->cd();
  
  int    bins = 240; // 240 is divisible by 2,3,4,612*24=280    -1/3 and  0.5 and 5/3  are bin edges

  double ptbw=10/100.0;  //see hPt histo below, let's make 10 bins per 1GeV/c

  const int nUEBins=20;
  double *uEBinBorders = new double[nUEBins+1];
  
  double uEa = fCard->GetBinBorder(kAssocType, 0), uEb = fCard->GetBinBorder(kAssocType, fCard->GetNoOfBins(kAssocType));
  double logUEbw = (log(uEb)-log(uEa))/nUEBins;
  for(int ij=0;ij<=nUEBins;ij++) uEBinBorders[ij]=uEa*exp(ij*logUEbw);


  if(fCard->GetNoOfBins(kCentrType) > kMaxNoCentrBin ){
    cout<<"ERROR: No of Centrality bins exceed max dim in AliJJtHistograms.cxx "<<endl;
    exit(0);
  }

  //==================================
  //  trigger pt fhistos
  //==================================

  double pTt1 = fPTtBin.GetMin();
  double pTt2 = fPTtBin.GetMax();
  double pTa1 = fPTaBin.GetMin();
  double pTa2 = fPTaBin.GetMax();

  fhIphiTrigg
      << TH1D( "fhIphiTrigg", "",  bins, -kJPi-0.1, kJPi+0.1)
      <<  fCentBin << fPTtBin  << "END";
  fhIetaTrigg
      << TH1D( "hIetaTrigg", "",  80, -fmaxEtaRange, fmaxEtaRange)
      <<  fCentBin << fPTtBin  << "END";// inclusive eta
  fhTriggPtBin
      << TH1D( "hTriggPtBin", "", (int)TMath::Ceil((pTt2-pTt1)/ptbw),pTt1, pTt2)
      <<  fCentBin << fVtxBin << fPTtBin  << "END";

  //=======================================
  //  associated fpt fhistos without etaGaps
  //=======================================
  fhAssocPtBin
      << TH1D( "hAssocPtBin", "", (int)TMath::Ceil((pTa2-pTa1)/ptbw), pTa1, pTa2)
      <<  fCentBin << fPTtBin << fPTaBin  << "END";
  fhIphiAssoc
      << TH1D( "fhIphiAssoc", "",  bins, -kJPi-0.1, kJPi+0.1)
      <<  fCentBin << fPTaBin  << "END";
  fhIetaAssoc
      << TH1D( "hIetaAssoc", "",  80, -fmaxEtaRange, fmaxEtaRange)
      <<  fCentBin << fPTaBin  << "END";
  
  //======================
  // invariant mass histograms
  //======================

  fhInvariantMassXe
      << TH1D("hInvariantMassXe","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhInvariantMassKlong
      << TH1D("hInvariantMassKlong","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";

  fhInvariantMassPta
      << TH1D("hInvariantMassPta","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

  // Like sign pairs for invariant mass histograms
    
  fhInvariantMassXeLikeSign
      << TH1D("hInvariantMassXeLikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhInvariantMassKlongLikeSign
      << TH1D("hInvariantMassKlongLikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";

  fhInvariantMassPtaLikeSign
      << TH1D("hInvariantMassPtaLikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

  // Unlike sign pairs for invariant mass histograms
    
  fhInvariantMassXeUnlikeSign
      << TH1D("hInvariantMassXeUnlikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhInvariantMassKlongUnlikeSign
      << TH1D("hInvariantMassKlongUnlikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";

  fhInvariantMassPtaUnlikeSign
      << TH1D("hInvariantMassPtaUnlikeSign","",1500,0,3)
      <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

  //======================================
  // Histograms for acceptance correction
  //======================================
  
  fhDetaNearMixAcceptance
      << TH1D( "hDEtaNearMixAcceptance", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
      <<  fCentBin << fPTtBin << fPTaBin  << "END";
  
  //=======================
  //jT fhistos
  //=======================
  
  int nJT = 100;
  double jtLow = 0.05, jtHigh = 20;
  
  double logBinsJt[101];
  double logJt = (log(jtHigh)-log(jtLow))/nJT;
  for(int ij=0;ij<=nJT;ij++) logBinsJt[ij]=jtLow*exp(ij*logJt);

  // Histograms in xlong bins
    
  fhJT
      << TH1D( "hJT", "",  nJT, logBinsJt)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhJTBg
      << TH1D( "hJTBg", "",  nJT, logBinsJt)
      <<  fCentBin << fEtaGapBin << fPTtBin << fXEBin  << "END";

  fhJTBgR
      << TH1D( "hJTBgR", "",  nJT, logBinsJt)
      <<  fCentBin << fRGapBin << fPTtBin << fXEBin  << "END";
  
  fhJTBgPhi
      << TH1D( "hJTBgPhi", "",  nJT, logBinsJt)
      <<  fCentBin << fPhiGapBin << fPTtBin << fXEBin  << "END";

  fhJTLikeSign
      << TH1D( "hJTLikeSign", "",  nJT, logBinsJt)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhJTBgLikeSign
      << TH1D( "hJTBgLikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fEtaGapBin << fPTtBin << fXEBin  << "END";

  fhJTBgRLikeSign
      << TH1D( "hJTBgRLikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fRGapBin << fPTtBin << fXEBin  << "END";
  
  fhJTBgPhiLikeSign
      << TH1D( "hJTBgPhiLikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fPhiGapBin << fPTtBin << fXEBin  << "END";

  fhJTUnlikeSign
      << TH1D( "hJTUnlikeSign", "",  nJT, logBinsJt)
      <<  fTypBin << fCentBin << fPTtBin << fXEBin  << "END";

  fhJTBgUnlikeSign
      << TH1D( "hJTBgUnlikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fEtaGapBin << fPTtBin << fXEBin  << "END";

  fhJTBgRUnlikeSign
      << TH1D( "hJTBgRUnlikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fRGapBin << fPTtBin << fXEBin  << "END";
  
  fhJTBgPhiUnlikeSign
      << TH1D( "hJTBgPhiUnlikeSign", "",  nJT, logBinsJt)
      <<  fCentBin << fPhiGapBin << fPTtBin << fXEBin  << "END";

  fhBgAssocXlongEta
      << TH1D( "hBgAssocXlongEta", "",  nUEBins, uEBinBorders)
      <<  fCentBin << fEtaGapBin << fPTtBin << fXEBin  << "END";

  fhBgAssocXlongR
      << TH1D( "hBgAssocXlongR", "",  nUEBins, uEBinBorders)
      <<  fCentBin << fRGapBin << fPTtBin << fXEBin  << "END";
  
  fhBgAssocXlongPhi
      << TH1D( "hBgAssocXlongPhi", "",  nUEBins, uEBinBorders)
      <<  fCentBin << fPhiGapBin << fPTtBin << fXEBin  << "END";

  // Only create deltaEta deltaPhi histograms if they are enabled in the JCard
  if(fCard->Get("EnableDeltaEtaDeltaPhiHistograms")==1){

    fhDphiDetaXlong
        << TH2D( "hDphiDetaXlong", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fTypBin <<  fCentBin << fVtxBin << fPTtBin << fXEBin  << "END";
    
    if(fCard->Get("TrackMergeSystematics")==1){
      
      fhDphiDetaTrackMerge
        << TH2D( "hDphiDetaTrackMerge", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fVtxBin << fPTtBin << fXEBin  << "END";
      
      if(fCard->Get("QualityControlLevel") > 0){
        
        fhDphiDetaTrackMergeCorrection
          << TH2D( "hDphiDetaTrackMergeCorrection", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
          <<  fCentBin << fVtxBin << fPTtBin << fXEBin  << "END";
        
      }
      
    }
    
  }
  
  if(fenable2DHistos){
      
    fhDphiDetaBgXlongEta
        << TH2D( "hDphiDetaBgXlongEta", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fPTtBin << fXEBin  << "END";

    fhDphiDetaBgXlongR
        << TH2D( "hDphiDetaBgXlongR", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fPTtBin << fXEBin  << "END";
  
    fhDphiDetaBgXlongPhi
        << TH2D( "hDphiDetaBgXlongPhi", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fPTtBin << fXEBin  << "END";

  }
  
  // Histograms in klong bins

  // Fill the klong histograms only if they are enabled from JCard
  if(fCard->Get("EnableKlongBins")==1){
  
    fhJTKlong
        << TH1D( "hJTKlong", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";
    
    fhJTKlongBg
        << TH1D( "hJTKlongBg", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fKLongBin  << "END";

    fhJTKlongBgR
        << TH1D( "hJTKlongBgR", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fKLongBin  << "END";
  
    fhJTKlongBgPhi
        << TH1D( "hJTKlongBgPhi", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fKLongBin  << "END";

    fhJTKlongLikeSign
        << TH1D( "hJTKlongLikeSign", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";
    
    fhJTKlongBgLikeSign
        << TH1D( "hJTKlongBgLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fKLongBin  << "END";

    fhJTKlongBgRLikeSign
        << TH1D( "hJTKlongBgRLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fKLongBin  << "END";
  
    fhJTKlongBgPhiLikeSign
        << TH1D( "hJTKlongBgPhiLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fKLongBin  << "END";

    fhJTKlongUnlikeSign
        << TH1D( "hJTKlongUnlikeSign", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fKLongBin  << "END";
    
    fhJTKlongBgUnlikeSign
        << TH1D( "hJTKlongBgUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fKLongBin  << "END";

    fhJTKlongBgRUnlikeSign
        << TH1D( "hJTKlongBgRUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fKLongBin  << "END";
  
    fhJTKlongBgPhiUnlikeSign
        << TH1D( "hJTKlongBgPhiUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fKLongBin  << "END";

    fhBgAssocKlongEta
        << TH1D( "hBgAssocKlongEta", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fEtaGapBin << fPTtBin << fKLongBin  << "END";

    fhBgAssocKlongR
        << TH1D( "hBgAssocKlongR", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fRGapBin << fPTtBin << fKLongBin  << "END";
  
    fhBgAssocKlongPhi
        << TH1D( "hBgAssocKlongPhi", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fPhiGapBin << fPTtBin << fKLongBin  << "END";

    // Only create deltaEta deltaPhi histograms if they are enabled in the JCard
    if(fCard->Get("EnableDeltaEtaDeltaPhiHistograms")==1){

      fhDphiDetaKlong
          << TH2D( "hDphiDetaKlong", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
          <<  fTypBin <<  fCentBin << fVtxBin << fPTtBin << fKLongBin  << "END";
  
    }
  
    if(fenable2DHistos){
      
      fhDphiDetaBgKlongEta
          << TH2D( "hDphiDetaBgKlongEta", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
          <<  fCentBin << fPTtBin << fKLongBin  << "END";

      fhDphiDetaBgKlongR
          << TH2D( "hDphiDetaBgKlongR", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
          <<  fCentBin << fPTtBin << fKLongBin  << "END";
  
      fhDphiDetaBgKlongPhi
          << TH2D( "hDphiDetaBgKlongPhi", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
          <<  fCentBin << fPTtBin << fKLongBin  << "END";

    }
  } // Klong binned histograms
  
  // Histograms in pta bins

  // Fill pTa bins only if they are enebled from JCard
  if(fCard->Get("EnablePtaBins")==1){
  
    fhJTPta
        << TH1D( "hJTPta", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBg
        << TH1D( "hJTPtaBg", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBgR
        << TH1D( "hJTPtaBgR", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fPTaBin  << "END";
  
    fhJTPtaBgPhi
        << TH1D( "hJTPtaBgPhi", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaLikeSign
        << TH1D( "hJTPtaLikeSign", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBgLikeSign
        << TH1D( "hJTPtaBgLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBgRLikeSign
        << TH1D( "hJTPtaBgRLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fPTaBin  << "END";
  
    fhJTPtaBgPhiLikeSign
        << TH1D( "hJTPtaBgPhiLikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaUnlikeSign
        << TH1D( "hJTPtaUnlikeSign", "",  nJT, logBinsJt)
        <<  fTypBin << fCentBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBgUnlikeSign
        << TH1D( "hJTPtaBgUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

    fhJTPtaBgRUnlikeSign
        << TH1D( "hJTPtaBgRUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fRGapBin << fPTtBin << fPTaBin  << "END";
  
    fhJTPtaBgPhiUnlikeSign
        << TH1D( "hJTPtaBgPhiUnlikeSign", "",  nJT, logBinsJt)
        <<  fCentBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

    fhBgAssocPtaEta
        << TH1D( "hBgAssocPtaEta", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

    fhBgAssocPtaR
        << TH1D( "hBgAssocPtaR", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fRGapBin << fPTtBin << fPTaBin  << "END";
  
    fhBgAssocPtaPhi
        << TH1D( "hBgAssocPtaPhi", "",  nUEBins, uEBinBorders)
        <<  fCentBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

    // Only create deltaEta deltaPhi histograms if they are enabled in the JCard
    if(fCard->Get("EnableDeltaEtaDeltaPhiHistograms")==1){
  
      fhDphiDetaPta
          << TH2D( "hDphiDetaPta", "", 400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
          <<  fTypBin <<  fCentBin << fVtxBin << fPTtBin << fPTaBin  << "END";
    
    }
    
    if(fenable2DHistos){
      
      fhDphiDetaBgPtaEta
          << TH2D( "hDphiDetaBgPtaEta", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
          <<  fCentBin << fPTtBin << fPTaBin  << "END";

      fhDphiDetaBgPtaR
          << TH2D( "hDphiDetaBgPtaR", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
          <<  fCentBin << fPTtBin << fPTaBin  << "END";
  
      fhDphiDetaBgPtaPhi
          << TH2D( "hDphiDetaBgPtaPhi", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
          <<  fCentBin << fPTtBin << fPTaBin  << "END";
      
    }
    
  }
  
  // Histograms for checking the acceptance correction is done correctly
  // Not needed in the analysis and should not be activated if not doing QA
  if(fEnableAcceptanceQAHistos){
    fhAcceptanceTraditional
        << TH1D( "hAcceptanceTraditional", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange)
        <<  fCentBin << fPTtBin << fPTaBin << "END";

    fhAcceptanceTraditional2D
        << TH2D( "hAcceptanceTraditional2D", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
        <<  fCentBin << fPTtBin << fPTaBin << "END";
  
    fhAcceptance3DNearSide
        << TH2D( "hAcceptance3DNearSide", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fPTtBin << fXEBin << "END";
    
    fhAcceptanceTraditional2DZ
        << TH2D( "hAcceptanceTraditional2DZ", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -kJPi/2, kJPi/2)
        <<  fCentBin << fVtxBin << fPTtBin << fPTaBin << "END";
    
    fhAcceptance3DNearSideZ
        << TH2D( "hAcceptance3DNearSideZ", "",  400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 640, -kJPi, kJPi)
        <<  fCentBin << fVtxBin << fPTtBin << fXEBin << "END";
  }
  
  delete [] uEBinBorders;
}


//______________________________________________________________________________
void AliJJtHistograms::CreateEventTrackHistos(){
  // Create basic event histograms
  
  CreateDataManagerHistograms();
  
  fHMG->cd();
  int nBINS=150;
  double logBinsX[nBINS+1], limL=0.1, limH=100;
  double logBW = (log(limH)-log(limL))/nBINS;
  for(int ij=0;ij<=nBINS;ij++) logBinsX[ij]=limL*exp(ij*logBW);

  fhLPpt       <<  TH1D("hLPpt","LP pt", nBINS, logBinsX ) << "END";
  fhChargedEta << TH1D("hChargedEta","All eta",100,-1.0,1.0)<< "END";
  fhLPeta      << TH1D("hLPeta","LP eta",100,-1.0,1.0)<< "END";

  fhChargedMult
      << TH1D("hChargedMult","", 300, 0., 3500.)
      << fCentBin << "END";
  fhCentr      << TH1D("hCentr","centrality", 101, -0.5, 100.5) << "END";
  fhiCentr         << TH1D("hiCentr","centrality",10, -0.5, 9.5) << "END";
  fhZVert
      << TH1D("hZVert", "", 100, -30., 30.)
      << fCentBin << "END";
  fhChargedPt
      << TH1D("hChargedPt","", nBINS, logBinsX )
      << fCentBin << "END";
  fhChargedPtNoCorr
      << TH1D("hChargedPtNoCorr","", nBINS, logBinsX )
      << fCentBin << "END";

  fhTrackingEfficiency << TProfile("hTrackingEff","",nBINS, logBinsX)
      << fCentBin << "END";
  
}


//______________________________________________________________________________
void AliJJtHistograms::ReadInclusiveHistos(const char *inclusFileName){
  // Read inclusive histograms
  fHMG->cd();
  
  TPMERegexp sep("::");
  int ncol = sep.Split( inclusFileName );
  TString filename = sep[0];
  
  if (TString(inclusFileName).BeginsWith("alien:"))  TGrid::Connect("alien:");
  TFile *inclusFile = TFile::Open(filename);
  TDirectory * dir =  (TDirectory*) inclusFile;
  if( ncol > 1 ) dir = (TDirectory*)( inclusFile->Get(sep[1]));
  if( !dir ) {
    cout << " ReadInclusiveHistos wrong file name or dirname !!!!" << endl;
    return;
  }
  
  cout<<inclusFileName<<"\t"<<filename<<"\t";
  if( ncol > 1 ) cout<<sep[1];
  cout<<endl;
  dir->Print();
  
  fHmgInclusive = new AliJHistManager("hst",sep[1]);
  fHmgInclusive->LoadConfig();
  
  fhIetaTriggFromFile = fHmgInclusive->GetTH1D("hIetaTrigg");
  fhIetaTriggFromFile.Print();
  fhIetaTriggFromFile[0][0]->Print();
  
  fhIphiTriggFromFile = fHmgInclusive->GetTH1D("fhIphiTrigg");
  fhIphiTriggFromFile.Print();
  fhIetaAssocFromFile = fHmgInclusive->GetTH1D("hIetaAssoc");
  fhIetaAssocFromFile.Print();
  fhIphiAssocFromFile = fHmgInclusive->GetTH1D("fhIphiAssoc");
  fhIphiAssocFromFile.Print();

}
