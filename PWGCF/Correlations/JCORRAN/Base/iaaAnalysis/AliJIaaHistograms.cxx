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
	//fhDEtaNear(),
	//fhDEtaNearM(),
	fhDphiDetaPta(),
	//fhResonanceCut(),
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
	fhChargedPtPublished(),
	fhTrackingEfficiency(),
	fhChargedEta(),
	fhLPeta(),
	fhChargedMult(),
	fhCentr(),
	fhiCentr(),
	fhZVert(),
	fhEP(),
	fhPhiS(),
	fmaxEtaRange(0)
{
	fmaxEtaRange = fCard->Get("EtaRange");
}

//______________________________________________________________________________
AliJIaaHistograms::AliJIaaHistograms(const AliJIaaHistograms& obj) :
	AliJHistogramInterface(obj),
	//fhDEtaNear(obj.fhDEtaNear),
	//fhDEtaNearM(obj.fhDEtaNearM),
	fhDphiDetaPta(obj.fhDphiDetaPta),
	//fhResonanceCut(obj.fhResonanceCut),
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
	fhChargedPtPublished(obj.fhChargedPtPublished),
	fhTrackingEfficiency(obj.fhTrackingEfficiency),
	fhChargedEta(obj.fhChargedEta),
	fhLPeta(obj.fhLPeta),
	fhChargedMult(obj.fhChargedMult),
	fhCentr(obj.fhCentr),
	fhiCentr(obj.fhiCentr),
	fhZVert(obj.fhZVert),
	fhEP(obj.fhEP),
	fhPhiS(obj.fhPhiS),
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
		<< TH1D( "fhIphiTrigg", "",  bins, -kJPi, kJPi)
		<<  fCentBin << fPTtBin  << "END";
	fhIetaTrigg
		<< TH1D( "hIetaTrigg", "",  80, -fmaxEtaRange, fmaxEtaRange)
		<<  fCentBin << fPTtBin  << "END";// inclusive eta
	fhTriggPtBin // all triggers
		<< TH1D( "hTriggPtBin", "", (int)TMath::Ceil((pTt2-pTt1)/ptbw),pTt1, pTt2)
		<<  fCentBin << fVtxBin << fPTtBin  << "END";

	//=======================================
	//  associated pt histos
	//=======================================
	fhAssocPtBin
		<< TH1D( "hAssocPtBin", "", (int)TMath::Ceil((pTa2-pTa1)/ptbw), pTa1, pTa2)
		<<  fCentBin << fPTtBin << fPTaBin  << "END";
	fhIphiAssoc
		<< TH1D( "fhIphiAssoc", "",  bins, -kJPi, kJPi)
		<<  fCentBin << fPTaBin  << "END";
	fhIetaAssoc
		<< TH1D( "hIetaAssoc", "",  80, -fmaxEtaRange, fmaxEtaRange)
		<<  fCentBin << fPTaBin  << "END";


	//======================================
	// correlation histos (and mixed event)
	//======================================

	//fhDEtaNear
	//		<< TH1D( "hDEtaNear", "",  160, -2*fmaxEtaRange, 2*fmaxEtaRange)
	//        <<  fCentBin << fVtxBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

	//fhDEtaNearM
	//		<< TH1D( "hDEtaNearM", "",  160, -2*fmaxEtaRange, 2*fmaxEtaRange)
	//        <<  fCentBin << fVtxBin << fPhiGapBin << fPTtBin << fPTaBin  << "END";

	fhDphiDetaPta
		<< TH2D( "hDphiDetaPta", "", 160, -2*fmaxEtaRange, 2*fmaxEtaRange, 80, -0.5*kJPi, 1.5*kJPi)
		<<  fTypBin <<  fCentBin << fVtxBin << fPTtBin << fPTaBin  << "END";

	//fhResonanceCut
	//	<< TH2D( "hResonanceCut", "", 160, -2*fmaxEtaRange, 2*fmaxEtaRange, 80, -0.5*kJPi, 1.5*kJPi)
	//	<<  fTypBin <<  fCentBin << fVtxBin << fPTtBin << fPTaBin  << "END";

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

	double xbins_pub[] = {0.15,   0.2,    0.25,   0.3,    0.35,   0.4,    0.45,   0.5,    0.55,   0.6,    0.65,   0.7,    0.75,   0.8,    0.85,   0.9,    0.95,   1,  1.1,    1.2,    1.3,    1.4,    1.5,    1.6,    1.7,    1.8,    1.9,    2,  2.2,    2.4,    2.6,    2.8,    3,  3.2,    3.4,    3.6,    3.8,    4,  4.5,    5,  5.5,    6,  6.5,    7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24,26,  28, 30, 32, 34, 36, 40, 45};
	const int nbins_pub = 64;

	fhLPpt       << TH1D("hLPpt","LP pt", nBINS, logBinsX ) << "END";
	fhChargedEta << TH1D("hChargedEta","All eta",100,-1.0,1.0)<< "END";
	fhLPeta      << TH1D("hLPeta","LP eta",100,-1.0,1.0)<< "END";
	fhCentr	     << TH1D("hCentr","centrality", 101, -0.5, 100.5) << "END";
	fhiCentr	 << TH1D("hiCentr","centrality",10, -0.5, 9.5) << "END";
	fhZVert      << TH1D("hZVert", "", 100, -30., 30.) << fCentBin << "END";

	fhChargedMult
		<< TH1D("hChargedMult","", 300, 0., 3500.)
		<< fCentBin << "END";
	fhZVert
		<< TH1D("hZVert", "", 100, -30., 30.)
		<< fCentBin << "END";
	fhEP
		<< TH1D( "fhEP", "",  240, -kJPi, kJPi)
		<<  fCentBin << "END";
	fhPhiS
		<< TH1D( "fhPhiS", "",  1200, -1., 360.1)
		<<  fCentBin << fPTtBin << "END"; // in degree 0-360
	fhChargedPt
		<< TH1D("hChargedPt","", nBINS, logBinsX )
		<< fCentBin << "END";
	fhChargedPtPublished
		<< TH1D("hChargedPtPublished","",nbins_pub,xbins_pub)
		<< fCentBin << "END";
	fhChargedPtNoCorr
		<< TH1D("hChargedPtNoCorr","", nBINS, logBinsX )
		<< fCentBin << "END";
	fhTrackingEfficiency
		<< TProfile("hTrackingEff","",nBINS, logBinsX)
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
