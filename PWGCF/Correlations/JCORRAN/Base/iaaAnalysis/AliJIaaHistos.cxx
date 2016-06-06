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

#include "AliJIaaHistos.h"
#include  "../AliJCard.h"
#include  "../AliJBaseTrack.h"
#include  "../AliJPhoton.h"
#include  "../AliJTrack.h"
#include <TGrid.h>
#include <TPRegexp.h>

//______________________________________________________________________________
AliJIaaHistos::AliJIaaHistos(AliJCard* cardP) :
    AliJHistogramInterface(cardP),
    fhMixStat(),
    fhPtNear(),
    fhPtFar(),
    fhPhi(),
    fhDphiAssoc(),
    fhDphiDetaPta(),
    fhDetaNearMixAcceptance(),
    fhDeta3DNearMixAcceptance(),
    fhDEtaNear(),
    fhDEtaNearM(),
    fhDEtaFar(),
    fhIphiTrigg(),
    fhIetaTrigg(),
    fhIphiAssoc(),
    fhIetaAssoc(),
    fhTriggPtBin(),
    fhTriggMult(),
    fhLPpt(),
    fhLPpairPt(),
    fhChargedPt(),
    fhChargedPtNoCorr(),
    fhChargedPtJacek(),
    fhChargedPtJacekEta(),
    fhChargedPtFiete(),
    fhVdelta2(),
    fhVdelta3(),
    fhVN(),
    fhTrackingEfficiency(),
    fpV2(),
    fpV3(),
    fpVdeltaNorm(),
    fhChargedEta(),
    fhLPeta(),
    fhAssocMult(),
    fhChargedMult(),
    fhChargedMultCut(),
    fhChargedMultCent(),
    fhV0AMult(),
    fhZVertRawErr(),
    fhZVert(),
    fhCentr(),
    fhiCentr(),
    fhEventPerRun(),
    fmaxEtaRange(0),
    fmaxTriggEtaRange(0),
    fLowRange(0),
    fHighRange(0),
    fenable2DHistos(false),
    fEnableAcceptanceQAHistos(false)
{   // constructor

    fmaxEtaRange = fCard->Get("EtaRange");

    const int nJacek =  73 ;
    double pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
      1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
      10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};

    //const int nJacek =  59 ;
    //double pttJacek[nJacek] = { 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.,
    //    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 9., 10,
    //    11., 12., 13., 14., 15., 16., 18., 20, 22, 24, 26., 28., 30.};

    fNJacek = nJacek;  // Number of bins
    fPttJacek = new double[fNJacek+1];
    for(int i = 0; i <= fNJacek; i++){
      fPttJacek[i] = pttJacek[i];
    }

    const int nEta = 3;
    double eta[nEta+1] = {-0.8,-0.2,0.3,0.8};// eta bins

    fNEta = nEta;  // Number of bins
    fEta = new double[fNEta+1];
    for(int i = 0; i <= fNEta; i++){
      fEta[i] = eta[i];
    }
}

//______________________________________________________________________________
AliJIaaHistos::AliJIaaHistos(const AliJIaaHistos& obj) :
    AliJHistogramInterface(obj),
    fhMixStat(obj.fhMixStat),
    fhPtNear(obj.fhPtNear),
    fhPtFar(obj.fhPtFar),
    fhPhi(obj.fhPhi),
    fhDphiAssoc(obj.fhDphiAssoc),
    fhDphiDetaPta(obj.fhDphiDetaPta),
    fhDetaNearMixAcceptance(obj.fhDetaNearMixAcceptance),
    fhDeta3DNearMixAcceptance(obj.fhDeta3DNearMixAcceptance),
    fhDEtaNear(obj.fhDEtaNear),
    fhDEtaNearM(obj.fhDEtaNearM),
    fhDEtaFar(obj.fhDEtaFar),
    fhIphiTrigg(obj.fhIphiTrigg),
    fhIetaTrigg(obj.fhIetaTrigg),
    fhIphiAssoc(obj.fhIphiAssoc),
    fhIetaAssoc(obj.fhIetaAssoc),
    fhTriggPtBin(obj.fhTriggPtBin),
    fhTriggMult(obj.fhTriggMult),
    fhLPpt(obj.fhLPpt),
    fhLPpairPt(obj.fhLPpairPt),
    fhChargedPt(obj.fhChargedPt),
    fhChargedPtNoCorr(obj.fhChargedPtNoCorr),
    fhChargedPtJacek(obj.fhChargedPtJacek),
    fhChargedPtJacekEta(obj.fhChargedPtJacekEta),
    fhChargedPtFiete(obj.fhChargedPtFiete),
    fhVdelta2(obj.fhVdelta2),
    fhVdelta3(obj.fhVdelta3),
    fhVN(obj.fhVN),
    fhTrackingEfficiency(obj.fhTrackingEfficiency),
    fpV2(obj.fpV2),
    fpV3(obj.fpV3),
    fpVdeltaNorm(obj.fpVdeltaNorm),
    fhChargedEta(obj.fhChargedEta),
    fhLPeta(obj.fhLPeta),
    fhAssocMult(obj.fhAssocMult),
    fhChargedMult(obj.fhChargedMult),
    fhChargedMultCut(obj.fhChargedMultCut),
    fhChargedMultCent(obj.fhChargedMultCent),
    fhV0AMult(obj.fhV0AMult),
    fhZVertRawErr(obj.fhZVertRawErr),
    fhZVert(obj.fhZVert),
    fhCentr(obj.fhCentr),
    fhiCentr(obj.fhiCentr),
    fhEventPerRun(obj.fhEventPerRun),
    fmaxEtaRange(obj.fmaxEtaRange),
    fmaxTriggEtaRange(obj.fmaxTriggEtaRange),
    fLowRange(obj.fLowRange),
    fHighRange(obj.fHighRange),
    fenable2DHistos(obj.fenable2DHistos),
    fEnableAcceptanceQAHistos(obj.fEnableAcceptanceQAHistos)
{
    // copy constructor
    JUNUSED(obj);
}

//______________________________________________________________________________
AliJIaaHistos& AliJIaaHistos::operator=(const AliJIaaHistos& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}

AliJIaaHistos::~AliJIaaHistos() {
    delete fHMG;
    //delete fHmgInclusive;
    //delete []fEta;
}

//______________________________________________________________________________
void AliJIaaHistos::CreateCorrelationHistos()
{
    // Create all the histograms needed in correlation analysis
    fHMG->cd();
  
    int    bins = 240; // 240 is divisible by 2,3,4,612*24=280    -1/3 and  0.5 and 5/3  are bin edges 
    fLowRange = -9./20.;   //lower range for dphi histos
    fHighRange= fLowRange+2;       //upper range for dphi histos;
    double ptbw=10/100.0;  //see hPt histo below, let's make 10 bins per 1GeV/c

    if(fCard->GetNoOfBins(kCentrType) > kMaxNoCentrBin ){
        cout<<"ERROR: No of Centrality bins exceed max dim in AliJHistos.cxx "<<endl;
        exit(0);
    }

    //==================================
    //  trigger pt fhistos 
    //==================================
    double pTt1 = fPTtBin.GetMin();
    double pTt2 = fPTtBin.GetMax();
    double pTa1 = fPTaBin.GetMin();
    double pTa2 = fPTaBin.GetMax();

    fhTriggMult
        << TH1D( "hTriggMult", "",  100, -0.5, 99.5) 
        <<  fCentBin << fPTtBin  << "END";
    fhIphiTrigg
        << TH1D( "fhIphiTrigg", "",  bins, -kJPi-0.1, kJPi+0.1) 
        <<  fCentBin << fPTtBin  << "END";
    fhIetaTrigg
        << TH1D( "hIetaTrigg", "",  80, -fmaxEtaRange, fmaxEtaRange) 
        <<  fCentBin << fPTtBin  << "END";// inclusive eta
    fhTriggPtBin
        << TH1D( "hTriggPtBin", "", (int)TMath::Ceil((pTt2-pTt1)/ptbw),pTt1, pTt2) 
        <<  fCentBin << fVtxBin << fPTtBin  << "END";

    //=====================================
    //  associated pt fhistos with etaGaps
    //=====================================

    fhDphiAssoc
        << TH1D( "hDphiAssoc", "",  bins, fLowRange, fHighRange) 
        <<  fTypBin << fCentBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";

    fhDEtaNear
        << TH1D( "hDEtaNear", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange) 
        <<  fCentBin << fVtxBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";
    fhDEtaNearM
        << TH1D( "hDEtaNearM", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange) 
        <<  fCentBin << fVtxBin << fEtaGapBin << fPTtBin << fPTaBin  << "END";


    //=======================================
    //  associated fpt fhistos without etaGaps
    //=======================================
    fhIphiAssoc
        << TH1D( "fhIphiAssoc", "",  bins, -kJPi-0.1, kJPi+0.1) 
        <<  fCentBin << fPTaBin  << "END";
    fhIetaAssoc
        << TH1D( "hIetaAssoc", "",  80, -fmaxEtaRange, fmaxEtaRange) 
        <<  fCentBin << fPTaBin  << "END";
    fhDEtaFar
        << TH1D( "hDEtaFar", "",  120, -2*fmaxEtaRange, 2*fmaxEtaRange) 
        <<  fTypBin << fCentBin << fPTtBin  << "END";

    fhDphiDetaPta
        << TH2D( "hDphiDetaPta", "", 400*fmaxEtaRange, -2*fmaxEtaRange, 2*fmaxEtaRange, 320, fLowRange, fHighRange)
        <<  fTypBin <<  fCentBin << fPTtBin << fPTaBin  << "END";


    //======================================
    // Histograms for acceptance correction
    //======================================
  
    fhDetaNearMixAcceptance
        << TH1D( "hDEtaNearMixAcceptance", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
        <<  fCentBin << fPTtBin << fPTaBin  << "END";
  
    fhDeta3DNearMixAcceptance
        << TH1D( "hDEta3DNearMixAcceptance", "",  320, -2*fmaxEtaRange, 2*fmaxEtaRange)
        <<  fCentBin << fPTtBin << fXEBin  << "END";
}

//______________________________________________________________________________
void AliJIaaHistos::CreateRunByRunHistos(int runID, int runcounter) const {
    // Todo
    fHMG->cd();
    JUNUSED(runID);
    JUNUSED(runcounter);
} //TODO


//______________________________________________________________________________
void AliJIaaHistos::CreateEventTrackHistos(){
    // comment needed
    fHMG->cd();
    int nBINS=150;
    double logBinsX[nBINS+1], limL=0.1, limH=100;
    double logBW = (log(limH)-log(limL))/nBINS;
    for(int ij=0;ij<=nBINS;ij++) logBinsX[ij]=limL*exp(ij*logBW);

    fhLPpt       <<  TH1D("hLPpt","LP pt", nBINS, logBinsX ) << "END";
    fhLPpairPt   << TH1D("hLPpairPt","LP pair pt", nBINS, logBinsX )<< "END";
    fhChargedEta << TH1D("hChargedEta","All eta",100,-1.0,1.0)<< "END";
    fhLPeta      << TH1D("hLPeta","LP eta",100,-1.0,1.0)<< "END";

    fhAssocMult << TH1D("hAssocMult","Assoc  mlt",100,-0.5,99.5)<< "END";

    fhChargedMult 
        << TH1D("hChargedMult","", 300, 0., 3500.)
        << fCentBin << "END"; 
    fhChargedMultCut 
        << TH1D("hChargedMultCut","",  300, 0., 3500.)
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

    fhVdelta2 << TH1D("hVdelta2","", 100, -0.2, 5)
        << fCentBin << "END";
    fhVdelta3 << TH1D("hVdelta3","", 100, -0.2, 5)
        << fCentBin << "END";
    fhVN << TH1D("hVN","", 100, -kJPi, kJPi)
        << fCentBin << "END";
    fhTrackingEfficiency << TProfile("hTrackingEff","",fNJacek, fPttJacek)
        << fCentBin << "END";
    fhV0AMult << TH1D("hV0Mult","", 2000,0,2000 )
        << fCentBin << "END";
    fpV2 << TProfile("pV2", "v2 with cent", 50, 0, 90) << "END";
    fpV3 << TProfile("pV3", "v3 with cent", 50, 0, 90) << "END";
    fpVdeltaNorm << TProfile("pVdeltaNorm", "mult with cent", 50, 0, 90) << "END";
    fhChargedMultCent << TH2D("hChargedMultCent ", "log(fmult) vs cent", 100, 0, 90, 100, 1, 10) << "END";
    fhZVertRaw << TH1D("hZVertRaw","vertex 0", 120, -30., 30.) << "END";
    fhZVertRawErr << TH1D("hZVertRawErr","vertex 0 Err", 100, 0, 0.1) << "END";
    fhCentr << TH1D("hCentr","centrality", 101, -0.5, 100.5) << "END";
    fhiCentr << TH1D("hiCentr","centrality",10, -0.5, 9.5) << "END";
    fhEventPerRun << TH1D("hEventPerRun","log(eve)/run",200, 0, 30.0) << "END";


    //------------------ for Abs Norm FK --------------------------------
    double   binsVertexMult[] = {0,1,2,3,4,5,10000};
    int   nbinsVertexMult  = sizeof(binsVertexMult)/sizeof(double)-1;
    //double binsVertexZ[]    = {-10,-5,-2,0,2,5,10};
    double binsVertexZ[]    = {-10,-6,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,6,10};
    int   nbinsVertexZ   = sizeof(binsVertexZ)/sizeof(double)-1;
    fhVertexZTriggVtx << TH2D("hVertexZTriggVtx","Vertex counts", nbinsVertexMult, binsVertexMult, nbinsVertexZ, binsVertexZ) << "END";

    // event counter
    fhEvents        
        << TH1D("hEvents","events passing cuts", 100, -0.5, 100-0.5 ) << "END";
    fhEventTrigger  
        << TH1D("hEventTrigger","Trigger count", 50, -0.5, 50.-.5 )
        << "END";
    fhTrackSelection
        << TH1D("hTrackSelection","checking bit convention", 100, -0.5, 100-0.5) << "END";
    // TODO fhEvents->SetXTitle( "0 - all, 1 - SDD selected, 2 - has vertex, 3 - good vertex, 4 - trigger + MB, 5 - trigger, 6 - BX, 7 - centrality" );
}


//______________________________________________________________________________
/*void AliJIaaHistos::ReadInclusiveHistos(const char *inclusFileName){
  // read inclusive histos
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
  
  //dir->cd();    // Not needed in new implementation of AliJHistManager
  
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
  
  // Check if the new acceptance histogram exists in the inclusive file.
  // This check is done for backwards compatibility. If the histogram does not exist,
  // histogram manager will return an empty histogram for the acceptance correction.
  // In AliJCorrelations.cxx, it is required that the histogram has at least 1000 entries,
  // otherwise a simple triangle is used for the acceptance correction. Thus the inclusive
  // code will still work with old data files missing the acceptance histograms.
  if(fHmgInclusive->HistogramExists("hDEtaNearMixAcceptance")){
    fhDEtaNearMixFromFile = fHmgInclusive->GetTH1D("hDEtaNearMixAcceptance");
    fhDEtaNearMixFromFile.Print();
    NormalizeAcceptanceHistos(fhDEtaNearMixFromFile, kAssocType);
  }
  
  if(fHmgInclusive->HistogramExists("hDEta3DNearMixAcceptance")){
    fhDEta3DNearMixFromFile = fHmgInclusive->GetTH1D("hDEta3DNearMixAcceptance");
    fhDEta3DNearMixFromFile.Print();
    NormalizeAcceptanceHistos(fhDEta3DNearMixFromFile, kXeType);
  }

}

void AliJIaaHistos::NormalizeAcceptanceHistos(AliJTH1D &acceptanceHisto, corrType assocType){
  // Method for normalizing and rebinning the inclusive acceptance histograms
  
  // Find the correct binning
  int numCent  = fCard->GetNoOfBins(kCentrType);
  int numPtt   = fCard->GetNoOfBins(kTriggType);
  int numAssoc = fCard->GetNoOfBins(assocType);
  
  // Loop over the input histograms and find the correct normalization
  for (int iCent = 0; iCent < numCent; iCent++) {
    for (int iPtt = 0; iPtt < numPtt; iPtt++){
      for (int iAssoc = 0; iAssoc < numAssoc; iAssoc++){
        
        // Rebin and normalize to the interval [0,1]
        double counts  = acceptanceHisto[iCent][iPtt][iAssoc]->Integral();
        int rebin = 4;
        if(counts<5000) rebin=8;
        if(counts<3000) rebin=10;
        if(counts<1000) rebin=16;
        acceptanceHisto[iCent][iPtt][iAssoc]->Rebin(rebin);
        double maxValue = acceptanceHisto[iCent][iPtt][iAssoc]->GetMaximum();
        if(maxValue > 0) acceptanceHisto[iCent][iPtt][iAssoc]->Scale(1.0/maxValue);
        
      }
    }
  }
  
}
*/
