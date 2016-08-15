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

#include "AliJIaaHistograms.h"
#include "../AliJCard.h"
#include "../AliJBaseTrack.h"
#include "../AliJPhoton.h"
#include "../AliJTrack.h"
#include <TGrid.h>
#include <TPRegexp.h>

//______________________________________________________________________________
AliJIaaHistograms::AliJIaaHistograms(AliJCard* cardP) :
	AliJHistogramInterface(cardP),
	fhDphiAssoc(),
	fhDEtaNear(),
	fhDEtaNearM(),
	fhDetaNearMixAcceptance(),
	fhDphiDetaPta(),
	//fhDphiDetaPtaRgap(),
	fhIphiTrigg(),
	fhIetaTrigg(),
	fhIphiAssoc(),
	fhIetaAssoc(),
	fhTriggPtBin(),
	fhAssocPtBin(),
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
	//fhAcceptanceTraditional(),
	//fhAcceptanceTraditional2D(),
	//fhAcceptance3DNearSide(),
	//fhAcceptanceTraditional2DZ(),
	//fhAcceptance3DNearSideZ(),
	fmaxEtaRange(0)
{   // constructor

    fmaxEtaRange = fCard->Get("EtaRange");
  
}

//______________________________________________________________________________
AliJIaaHistograms::AliJIaaHistograms(const AliJIaaHistograms& obj) :
  AliJHistogramInterface(obj),
  fhDphiAssoc(obj.fhDphiAssoc),
  fhDEtaNear(obj.fhDEtaNear),
  fhDEtaNearM(obj.fhDEtaNearM),
  fhDetaNearMixAcceptance(obj.fhDetaNearMixAcceptance),
  fhDphiDetaPta(obj.fhDphiDetaPta),
  //fhDphiDetaPtaRgap(obj.fhDphiDetaPtaRgap),
  fhIphiTrigg(obj.fhIphiTrigg),
  fhIetaTrigg(obj.fhIetaTrigg),
  fhIphiAssoc(obj.fhIphiAssoc),
  fhIetaAssoc(obj.fhIetaAssoc),
  fhTriggPtBin(obj.fhTriggPtBin),
  fhAssocPtBin(obj.fhAssocPtBin),
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
  //fhAcceptanceTraditional(obj.fhAcceptanceTraditional),
  //fhAcceptanceTraditional2D(obj.fhAcceptanceTraditional2D),
  //fhAcceptance3DNearSide(obj.fhAcceptance3DNearSide),
  //fhAcceptanceTraditional2DZ(obj.fhAcceptanceTraditional2DZ),
  //fhAcceptance3DNearSideZ(obj.fhAcceptance3DNearSideZ),
  fmaxEtaRange(obj.fmaxEtaRange)
{
    // copy constructor
    JUNUSED(obj);
}

//______________________________________________________________________________
AliJIaaHistograms& AliJIaaHistograms::operator=(const AliJIaaHistograms& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}

//______________________________________________________________________________
AliJIaaHistograms::~AliJIaaHistograms() {
  // destructor
	delete fHMG;
	delete fHmgInclusive;
}

//______________________________________________________________________________
void AliJIaaHistograms::CreateCorrelationHistograms()
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
	cout<<"ERROR: No of Centrality bins exceed max dim in AliJIaaHistograms.cxx "<<endl;
    exit(0);
  }

  //==================================
  //  trigger pt histos
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
  //  associated pt histos
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
  

  //======================================
  // correlation histos (and mixed event)
  //======================================
  
  fhDetaNearMixAcceptance
      << TH1D( "hDEtaNearMixAcceptance", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
      <<  fCentBin << fPTtBin << fPTaBin  << "END";
  
  fhDphiAssoc
	  << TH1D( "hDphiAssoc", "",  320, -9./20, 31./20.)
	  <<  fTypBin << fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

//  fhDEtaNear
//	  << TH1D( "hDEtaNear", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
//	  <<  fCentBin << fVtxBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";
//
//  fhDEtaNearM
//	  << TH1D( "hDEtaNearM", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
//	  <<  fCentBin << fVtxBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

  fhDEtaNear
	  << TH1D( "hDEtaNear", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
	  <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

  fhDEtaNearM
	  << TH1D( "hDEtaNearM", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
	  <<  fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

//<< TH2D( "hDphiDetaPta", "", 100*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, -9./20, 31./20.)
  fhDphiDetaPta
	  << TH2D( "hDphiDetaPta", "", 20*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 72, -9./20, 31./20.)
	  <<  fTypBin <<  fCentBin << fPTtBin << fPTaBin  << "END";

//  fhDphiDetaPtaRgap
//      << TH2D( "fhDphiDetaPtaRgap", "", 100*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 100, -0.5, 1.5)
//      << fTypBin << fRGapBin << fCentBin << fPTtBin << fPTaBin  << "END";
  
  delete [] uEBinBorders;
}


//______________________________________________________________________________
void AliJIaaHistograms::CreateEventTrackHistos(){
  // Create basic event histograms
  
  CreateDataManagerHistograms();
  
  fHMG->cd();
  int nBINS=150;
  double logBinsX[nBINS+1], limL=0.1, limH=100;
  double logBW = (log(limH)-log(limL))/nBINS;
  for(int ij=0;ij<=nBINS;ij++) logBinsX[ij]=limL*exp(ij*logBW);

  fhLPpt       << TH1D("hLPpt","LP pt", nBINS, logBinsX ) << "END";
  fhChargedEta << TH1D("hChargedEta","All eta",100,-1.0,1.0)<< "END";
  fhLPeta      << TH1D("hLPeta","LP eta",100,-1.0,1.0)<< "END";
  fhCentr	   << TH1D("hCentr","centrality", 101, -0.5, 100.5) << "END";
  fhiCentr	   << TH1D("hiCentr","centrality",10, -0.5, 9.5) << "END";
  fhZVert      << TH1D("hZVert", "", 100, -30., 30.) << fCentBin << "END";

  fhChargedMult
      << TH1D("hChargedMult","", 300, 0., 3500.)
      << fCentBin << "END";
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
void AliJIaaHistograms::ReadInclusiveHistos(const char *inclusFileName){
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
