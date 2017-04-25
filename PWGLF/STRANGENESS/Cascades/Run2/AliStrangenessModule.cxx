//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// General class to do weak decay analysis
//
// WARNING: THIS IS EXPERIMENTAL!
//
// Please send any questions, etc to:
//    david.dobrigkeit.chinellato@cern.ch
//
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TStopwatch.h"
#include "AliVWeakResult.h"
#include "AliV0Result.h"
#include "AliCascadeResult.h"
#include "AliStrangenessModule.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliStrangenessModule);

//________________________________________________________________
AliStrangenessModule::AliStrangenessModule() :
TNamed(),
lDataInput(0x0),
lDataCountersInput(0x0),
lLSDataInput(0x0),
lLSDataCountersInput(0x0),
lMCInput(0x0),
lMCCountersInput(0x0),
lLoMult(0.0),
lHiMult(90.0),
lUseIntegratedMultForFirstFit(kTRUE),
lLoMultIntegrated(0.0),
lHiMultIntegrated(90.0),
lOffsetFromMeanLeft(-0.025),
lOffsetFromMeanRight(+0.025),
lNPtBins(0),
lLoLeftBg(-6),
lHiLeftBg(-4),
lLoPeak(-3),
lHiPeak(3),
lLoRightBg(4),
lHiRightBg(7),
lLSLoMass(1.34),
lLSHiMass(1.36),
lObjectToUseForLS(""),
fF1Mean(0x0),
fF1Sigma(0x0),
lVerbose(kFALSE),
lUseGeant3FlukaCorrection(kFALSE),
lDoOnlyData(kFALSE)
{
    // Dummy Constructor - not to be used!
    for(Long_t ibin = 0; ibin<100; ibin++) lPtBins[ibin] = 0;
    for(Long_t ibin = 0; ibin<100; ibin++) lSigExtTech[ibin] = "linear";
    for(Long_t ibin = 0; ibin<100; ibin++) lSigExtSubLS[ibin] = kFALSE;
}

//________________________________________________________________
AliStrangenessModule::AliStrangenessModule(const char * name, const char * title) :
TNamed(name,title),
lDataInput(0x0),
lDataCountersInput(0x0),
lLSDataInput(0x0),
lLSDataCountersInput(0x0),
lMCInput(0x0),
lMCCountersInput(0x0),
lLoMult(0.0),
lHiMult(90.0),
lUseIntegratedMultForFirstFit(kTRUE),
lLoMultIntegrated(0.0),
lHiMultIntegrated(90.0),
lOffsetFromMeanLeft(-0.025),
lOffsetFromMeanRight(+0.025),
lNPtBins(0),
lLoLeftBg(-6),
lHiLeftBg(-4),
lLoPeak(-3),
lHiPeak(3),
lLoRightBg(4),
lHiRightBg(7),
lLSLoMass(1.34),
lLSHiMass(1.36),
lObjectToUseForLS(""),
fF1Mean(0x0),
fF1Sigma(0x0),
lVerbose(kFALSE),
lUseGeant3FlukaCorrection(kFALSE),
lDoOnlyData(kFALSE)
{
    // Main constructor
    for(Long_t ibin = 0; ibin<100; ibin++) lPtBins[ibin] = 0;
    for(Long_t ibin = 0; ibin<100; ibin++) lSigExtTech[ibin] = "linear";
    for(Long_t ibin = 0; ibin<100; ibin++) lSigExtSubLS[ibin] = kFALSE;

}

//________________________________________________________________
AliStrangenessModule::~AliStrangenessModule(){
    // Nothing to delete
}

//________________________________________________________________
void AliStrangenessModule::SetPtBinning(  Long_t lRecNPtBins, Double_t *lRecPtBins  ){
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    lNPtBins = lRecNPtBins;
    for(int ix = 0;ix<lNPtBins+1;ix++){
        lPtBins[ix] = lRecPtBins[ix];
    }
}

//________________________________________________________________
void AliStrangenessModule::SetSigExtRanges( Double_t lRLoLeftBg, Double_t lRHiLeftBg,  Double_t lRLoPeak,
                                           Double_t lRHiPeak,   Double_t lRLoRightBg, Double_t lRHiRightBg ){
    //Function to set signal extraction ranges (in sigmas!)
    lLoLeftBg  = lRLoLeftBg;
    lHiLeftBg  = lRHiLeftBg;
    lLoPeak    = lRLoPeak;
    lHiPeak    = lRHiPeak;
    lLoRightBg = lRLoRightBg;
    lHiRightBg = lRHiRightBg;
}

//________________________________________________________________
void AliStrangenessModule::SetSigExtTech ( TString lRecSigExtTech ) {
    if( !lRecSigExtTech.Contains("bincounting") &&
       !lRecSigExtTech.Contains("linear") &&
       !lRecSigExtTech.Contains("quadratic") &&
       !lRecSigExtTech.Contains("cubic") &&
       !lRecSigExtTech.Contains("MC") ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(Form(" Sig. ext. mode \"%s\" not recognized!",lRecSigExtTech.Data() ) );
        AliWarning("   Accepted modes: \"bincounting\", \"linear\", \"quadratic\" or \"cubic\"");
        AliWarning("                     WARNING: Will set to linear! ");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        lRecSigExtTech = "linear";
    }
    for(Int_t ibin = 0; ibin<100; ibin++){
        lSigExtTech[ibin] = lRecSigExtTech.Data();
    }
}

//________________________________________________________________
void AliStrangenessModule::SetVariableSigExtTech ( Long_t lRecNPtBins, TString *lRecSigExtTech ) {
    if (lVerbose ) cout<<"Variable signal extraction technique received!"<<endl;
    for(Int_t ibin = 0; ibin<lRecNPtBins; ibin++){
        if( !lRecSigExtTech[ibin].Contains("bincounting") &&
           !lRecSigExtTech[ibin].Contains("linear") &&
           !lRecSigExtTech[ibin].Contains("quadratic") &&
           !lRecSigExtTech[ibin].Contains("cubic") &&
           !lRecSigExtTech[ibin].Contains("MC") ){
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            AliWarning(Form(" Sig. ext. mode \"%s\" not recognized!",lRecSigExtTech[ibin].Data() ) );
            AliWarning("   Accepted modes: \"bincounting\", \"linear\", \"quadratic\" or \"cubic\"");
            AliWarning("                     WARNING: Will set to linear! ");
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            lRecSigExtTech[ibin] = "linear";
        }
        lSigExtTech[ibin] = lRecSigExtTech[ibin].Data();
        if( lVerbose ) cout<<"Bin #"<<ibin<<Form("\t%.1f-%1.f",lPtBins[ibin],lPtBins[ibin+1])<<" Technique: "<<lRecSigExtTech[ibin]<<endl;
    }
}

//________________________________________________________________
void AliStrangenessModule::SetStandardMassFunctionsPbPb2015 () {
    if (lVerbose ) cout<<"Setting up standard mass / mean functions for PbPb 2015 data"<<endl;
    
    TF1 *fMeanFit = new TF1("fMeanFit", "[0]+[1]*TMath::Exp(-[2]*x)+[3]*TMath::Exp(-[4]*x)",0,10);
    fMeanFit->SetParameter(0, 1.32510e+00);
    fMeanFit->SetParameter(1, -3.41173e-03);
    fMeanFit->SetParameter(2, 8.92321e-02);
    fMeanFit->SetParameter(3, -3.30786e-03);
    fMeanFit->SetParameter(4, 2.62183e+00);
    
    TF1 *fSigmaFit = new TF1("fSigmaFit", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)",0,10);
    fSigmaFit->SetParameter(0, 1.40024e-03);
    fSigmaFit->SetParameter(1, 3.52738e-04);
    fSigmaFit->SetParameter(2, 2.20490e-03);
    fSigmaFit->SetParameter(3, 1.83156e+00);
    
    SetMassMeanFunction(fMeanFit);
    SetMassSigmaFunction(fSigmaFit);
}

//________________________________________________________________
TH1D* AliStrangenessModule::DoAnalysis( TString lConfiguration, TString lOutputFile ){
    //Main analysis code
    
    //Report time at the end
    TStopwatch* timer = new TStopwatch();
    timer->Start ( kTRUE );
    
    cout<<"AliStrangenessModule -> Requested to process config: "<<lConfiguration.Data()<<endl;
    
    //"Die hard" fitting, please
    TVirtualFitter::SetMaxIterations( 50000 );
    
    //Check if any LS subtraction requested
    Bool_t lRequestLSSub = kFALSE;
    for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
        if ( lSigExtSubLS[ibin] == kTRUE ) lRequestLSSub = kTRUE;
    }
    
    //Check if input has been provied
    if( !lDataInput ) {
        AliWarning("Data input not defined!"); return 0x0;
    }
    if( !lDataCountersInput ) {
        AliWarning("Data counters list not defined!"); return 0x0;
    }
    if( !lLSDataInput && lRequestLSSub ) {
        AliWarning("Like-sign data input not defined!"); return 0x0;
    }
    if( !lLSDataCountersInput && lRequestLSSub ) {
        AliWarning("Like-sign data counters list not defined!"); return 0x0;
    }
    if( !lMCInput && !lDoOnlyData ) {
        AliWarning("MC input not defined!"); return 0x0;
    }
    if( !lMCCountersInput && !lDoOnlyData  ) {
        AliWarning("MC counters list not defined!"); return 0x0;
    }
    
    //Open File for output: default name if not defined in DoAnalysis call
    if ( lOutputFile.EqualTo("") ) lOutputFile = Form("Results_%s.root",lConfiguration.Data() );
    
    TFile *fFileOut = new TFile(lOutputFile.Data(), "RECREATE");
    TList *fListData = new TList();
    fListData->SetName("cListData");
    fListData->SetOwner(kTRUE);

    TList *fListMC = new TList();
    fListMC->SetName("cListMC");
    fListMC->SetOwner(kTRUE);
    
    TList *fListOutput = new TList();
    fListOutput->SetName("cAnalysisOutput");
    fListOutput->SetOwner(kTRUE);
    
    //In preparation: get number of events from real data
    TH1D *fHistCentrality = (TH1D*) lDataCountersInput->FindObject( "fHistCentrality" )->Clone("fHistCentralityClone");
    fHistCentrality->SetDirectory(0);
    Double_t lNEvents = fHistCentrality->Integral( fHistCentrality->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                                  fHistCentrality->GetXaxis()->FindBin( lHiMult-1e-5 ) );
    cout<<"AliStrangenessModule -> Number of events in multiplicity class "<<lLoMult<<" to "<<lHiMult<<": "<<lNEvents<<endl;
    
    //In preparation bis: get number of events from MC data
    Double_t lNEventsMC = 0;
    if ( !lDoOnlyData ) {
        TH1D *fHistCentralityMC = (TH1D*) lMCCountersInput->FindObject( "fHistCentrality" )->Clone("fHistCentralityCloneMC");
        fHistCentralityMC->SetDirectory(0);
        lNEventsMC = fHistCentralityMC->Integral( fHistCentralityMC->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                                 fHistCentralityMC->GetXaxis()->FindBin( lHiMult-1e-5 ) );
    }
    
    //Save number of events to output
    TH1D *fHistEventCounter = new TH1D("fHistEventCounter", "", 2,0,2);
    fHistEventCounter->SetDirectory(0);
    
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "# Events (Data)");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "# Events (MC)");
    
    fHistEventCounter->SetBinContent(1,lNEvents);
    fHistEventCounter->SetBinContent(2,lNEventsMC);
    
    fListOutput->Add(fHistEventCounter);
    
    //Step 1: Open real data object
    TString lDataName = lConfiguration.Data();
    lDataName.Append("_Data");
    AliVWeakResult *lDataResult = (AliVWeakResult*) lDataInput->FindObject(lConfiguration.Data())->Clone( lDataName.Data() );
    
    //============================================================================
    //Do bookkeeping of base input 
    fListData->Add(lDataResult);
    //============================================================================
    
    if(lVerbose) lDataResult->Print();
    
    //_________________________________________________
    // Process TH3Fs and expand into histograms of interest
    TH3F *f3dHistData = (TH3F*) lDataResult->GetHistogram()->Clone("f3dHistData");
    
    //Check if multiplicity interval requested is possible
    Bool_t lCheckMult = CheckCompatibleMultiplicity ( f3dHistData );
    if ( !lCheckMult ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(" LS: Requested mult interval is inconsistent with the data provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Check if pt interval requested is possible
    Bool_t lCheckPt = CheckCompatiblePt ( f3dHistData );
    if ( !lCheckPt ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(" LS: Requested pt binning is inconsistent with the data provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Project into relevant invariant mass histograms
    TH1D *lHistoData[100];

    //If requested to use integrated for first fit
    TH1D *lHistoDataIntegrated[100];
    
    for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
        lHistoData[ibin] = f3dHistData->ProjectionZ( Form("lHistoData_%.1f_%.1f", lPtBins[ibin], lPtBins[ibin+1]),
                                              f3dHistData->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                              f3dHistData->GetXaxis()->FindBin( lHiMult-1e-5 ),
                                              f3dHistData->GetYaxis()->FindBin( lPtBins[ibin  ]+1e-5 ),
                                              f3dHistData->GetYaxis()->FindBin( lPtBins[ibin+1]-1e-5 )
                                              );
        lHistoData[ibin]->SetDirectory(0); 
        fListData->Add(lHistoData[ibin]);
    }
    
    if ( lUseIntegratedMultForFirstFit ){
        for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
            lHistoDataIntegrated[ibin] = f3dHistData->ProjectionZ( Form("lHistoDataIntegrated_%.1f_%.1f", lPtBins[ibin], lPtBins[ibin+1]),
                                                                  f3dHistData->GetXaxis()->FindBin( lLoMultIntegrated+1e-5 ),
                                                                  f3dHistData->GetXaxis()->FindBin( lHiMultIntegrated-1e-5 ),
                                                                  f3dHistData->GetYaxis()->FindBin( lPtBins[ibin  ]+1e-5 ),
                                                                  f3dHistData->GetYaxis()->FindBin( lPtBins[ibin+1]-1e-5 )
                                                                  );
            lHistoDataIntegrated[ibin]->SetDirectory(0);
            fListData->Add(lHistoDataIntegrated[ibin]);
        }
    }
    
    //If requested, do like-sign subtraction here
    if( lRequestLSSub ){
        //_________________________________________________
        //Project into relevant invariant mass histograms
        TH1D *lHistoLSData[100];
        //_________________________________________________
        // Basic I/O
        if(lVerbose) cout<<"AliStrangenessModule -> Like-sign BG subtraction requested for some bins, doing now..."<<endl; 
        AliVWeakResult *lLSDataResult = (AliVWeakResult*) lLSDataInput->FindObject(lObjectToUseForLS.Data())->Clone( "LikeSignObject" );
        if(lVerbose) lLSDataResult->Print();
        //Compatibility check zero: check if generated with the same cuts
        if( !lLSDataResult->HasSameCuts(lDataResult) ){
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            AliWarning(" You are attempting to correct data with LS results obtained");
            AliWarning(" with different selection criteria! Will stop here. Please ");
            AliWarning(" check your input files!");
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            return 0x0;
        }
        //_________________________________________________
        // Process TH3Fs and expand into histograms of interest
        TH3F *f3dHistLSData = (TH3F*) lLSDataResult->GetHistogram()->Clone("f3dHistLSData");
        //Check if multiplicity interval requested is possible
        Bool_t lCheckMult = CheckCompatibleMultiplicity ( f3dHistLSData );
        if ( !lCheckMult ){
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            AliWarning(" LS: Requested mult interval is inconsistent with the data provided!");
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            return 0x0;
        }
        
        //Check if pt interval requested is possible
        Bool_t lCheckPt = CheckCompatiblePt ( f3dHistLSData );
        if ( !lCheckPt ){
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            AliWarning(" LS: Requested pt binning is inconsistent with the data provided!");
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            return 0x0;
        }
        for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
            lHistoLSData[ibin] = f3dHistLSData->ProjectionZ( Form("lHistoLSData_%.1f_%.1f", lPtBins[ibin], lPtBins[ibin+1]),
                                                        f3dHistLSData->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                                        f3dHistLSData->GetXaxis()->FindBin( lHiMult-1e-5 ),
                                                        f3dHistLSData->GetYaxis()->FindBin( lPtBins[ibin  ]+1e-5 ),
                                                        f3dHistLSData->GetYaxis()->FindBin( lPtBins[ibin+1]-1e-5 )
                                                        );
            lHistoLSData[ibin]->SetDirectory(0);
            fListData->Add(lHistoLSData[ibin]);
        }

        
        
        //_________________________________________________
        // Subtract stuff, please
        for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
            //_____________________________________________
            // Renormalize such that we get similar yields in a specific (settable) region
            Int_t lLoBinForNorm = lHistoData[ibin] -> FindBin( lLSLoMass+1e-5 );
            Int_t lHiBinForNorm = lHistoData[ibin] -> FindBin( lLSHiMass-1e-5 );
            
            Int_t lLoLSBinForNorm = lHistoLSData[ibin] -> FindBin( lLSLoMass+1e-5 );
            Int_t lHiLSBinForNorm = lHistoLSData[ibin] -> FindBin( lLSHiMass-1e-5 );
            
            Double_t lLSScalingRatio = -100;
            Double_t lIntegralData   = lHistoData[ibin]  ->Integral(lLoBinForNorm,  lHiBinForNorm  );
            Double_t lIntegralLSData = lHistoLSData[ibin]->Integral(lLoLSBinForNorm,lHiLSBinForNorm);
            
            if ( TMath::Abs(lIntegralLSData) > 1e-5 ) lLSScalingRatio = lIntegralData/lIntegralLSData;
            
            //Beware error propagation!
            lHistoLSData[ibin] -> Scale( lLSScalingRatio );
            
            //_____________________________________________
            // Proceed with subtraction
            if(lSigExtSubLS[ibin] ) {
                lHistoData[ibin]->Add( lHistoLSData[ibin], -1 );
            }
        }
    }
    
    
    //Perform early fit with gaussian+linear for acquisition of fit range
    Double_t lMeanVsPt[100], lMeanErrVsPt[100];
    Double_t lSigmaVsPt[100], lSigmaErrVsPt[100];
    TH1D* fHistMeanVsPt  = new TH1D("fHistMeanVsPt", "",lNPtBins,lPtBins);
    TH1D* fHistSigmaVsPt = new TH1D("fHistSigmaVsPt","",lNPtBins,lPtBins);
    fHistMeanVsPt->SetDirectory(0);
    fHistSigmaVsPt->SetDirectory(0);
    if(!lVerbose) cout<<"AliStrangenessModule -> Initial fit (data):       ["<<flush;
    Bool_t lExtStatus = kTRUE;
    TH1D *lPointerToRelevantHisto = 0x0;
    
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        lMeanVsPt[ibin] = lDataResult->GetMass(); //Initial guess
        
        //This is the procedure in case mean, sigma functions haven't been set
        if ( fF1Mean==0x0 || fF1Sigma==0x0 ){
        lPointerToRelevantHisto = lHistoData[ibin];
        if ( lUseIntegratedMultForFirstFit ) lPointerToRelevantHisto = lHistoDataIntegrated[ibin];
        
        if(lVerbose) AliWarning(Form("Initial fit on: %s",lHistoData[ibin]->GetName()));
        lExtStatus = PerformInitialFit( lPointerToRelevantHisto, lMeanVsPt[ibin], lMeanErrVsPt[ibin], lSigmaVsPt[ibin], lSigmaErrVsPt[ibin], fListData );
        if( !lVerbose ){
            if( lExtStatus ){ cout<<"="<<flush; } else { cout<<"!"<<flush; }
        }
        }else{
            lMeanVsPt[ibin]  = fF1Mean ->Eval(0.5*(lPtBins[ibin+1]+lPtBins[ibin]));
            lSigmaVsPt[ibin] = fF1Sigma->Eval(0.5*(lPtBins[ibin+1]+lPtBins[ibin]));
            lMeanErrVsPt[ibin]  = 0.0;
            lSigmaErrVsPt[ibin] = 0.0;
            if(lVerbose) AliWarning(Form("Received input mass/mean for %s: mass = %.3f, mean = %.3f ",lHistoData[ibin]->GetName(),
                                         lMeanVsPt[ibin],lSigmaVsPt[ibin]));
        }
        fHistMeanVsPt->SetBinContent(ibin+1, lMeanVsPt[ibin] );
        fHistMeanVsPt->SetBinError(ibin+1, lMeanErrVsPt[ibin] );
        fHistSigmaVsPt->SetBinContent(ibin+1, lSigmaVsPt[ibin] );
        fHistSigmaVsPt->SetBinError(ibin+1, lSigmaErrVsPt[ibin] );
    }
    if(!lVerbose) cout<<"] Done!"<<endl;
    
    if( lVerbose){
        cout<<"---] Peak position, width summary [-------------------------"<<endl;
        for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
            if( lVerbose ) cout<<"Bin #"<<ibin<<Form("\t%.1f-%.1f",lPtBins[ibin],lPtBins[ibin+1])<<" Peak Pos: "<<lMeanVsPt[ibin]<<", width: "<<lSigmaVsPt[ibin] <<endl;
        }
        cout<<"---] End Peak position, width summary [---------------------"<<endl;
    }
    
    //FIXME: Invoke TH1D::Smooth here if asked to do so 
    
    fListData->Add( fHistMeanVsPt );
    fListData->Add( fHistSigmaVsPt );
    
    //Perform actual signal extraction
    Double_t lSignalVsPt[100], lSignalErrVsPt[100];
    TH1D* fHistRawVsPt  = new TH1D("fHistRawVsPt", "",lNPtBins,lPtBins);
    fHistRawVsPt->SetDirectory(0);
    TH1D* fHistBgVsPt  = new TH1D("fHistBgVsPt", "",lNPtBins,lPtBins);
    fHistBgVsPt->SetDirectory(0);
    
    if(!lVerbose) cout<<"AliStrangenessModule -> Extracting signal (data): ["<<flush;
    lExtStatus = kTRUE;
    Double_t lBackground = 0;
    Double_t lBackgroundError = 0;
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        if( lVerbose) cout<<"Extracting yield for bin #"<<ibin<<Form("\t%.1f-%.1f",lPtBins[ibin],lPtBins[ibin+1])<<endl;
        lExtStatus = PerformSignalExtraction( lHistoData[ibin], lSignalVsPt[ibin], lSignalErrVsPt[ibin], lBackground, lBackgroundError, lMeanVsPt[ibin], lSigmaVsPt[ibin], fListData, lSigExtTech[ibin].Data() );
        if( !lVerbose ){
            //Report errors with "!"
            if( lExtStatus ){ cout<<"="<<flush; } else { cout<<"!"<<flush; }
        }
        fHistRawVsPt->SetBinContent(ibin+1, lSignalVsPt[ibin] );
        fHistRawVsPt->SetBinError(ibin+1, lSignalErrVsPt[ibin] );
        fHistBgVsPt->SetBinContent(ibin+1, lBackground );
        fHistBgVsPt->SetBinError(ibin+1, lBackgroundError );
    }
    if( lVerbose){
        cout<<"---] Signal Extraction summary [----------------------------"<<endl;
        for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
            if( lVerbose ) cout<<"Bin #"<<ibin<<Form("\t%.1f-%.1f",lPtBins[ibin],lPtBins[ibin+1])<<" Signal: "<<lSignalVsPt[ibin]<<" +/- "<<lSignalErrVsPt[ibin] <<" via technique: "<<lSigExtTech[ibin]<<endl;
        }
        cout<<"---] End Signal Extraction summary [------------------------"<<endl;
    }
    
    if(!lVerbose) cout<<"] Done!"<<endl;
    
    //Add raw spectra as main analysis output
    fListOutput->Add(fHistRawVsPt);
    fListOutput->Add(fHistBgVsPt);
    
    if(lDoOnlyData){
        fFileOut->cd();
        //Save all objects owned by the TLists
        fListData  ->Write("cListData",       TObject::kSingleKey);
        fListOutput->Write("cAnalysisOutput", TObject::kSingleKey);
        fFileOut->Write();
        fFileOut->Close();
        
        timer->Stop();
        cout<<"AliStrangenessModule -> Done! Total time: "<<timer->RealTime()<<"s (real), "<< timer->CpuTime() <<"s (cpu) (processed data only!)"<<endl;
        delete timer;
        return 0x0;
    }
    
    //Step N: Open MC data object
    TString lMCName = lConfiguration.Data();
    lMCName.Append("_MC");
    AliVWeakResult *lMCResult = (AliVWeakResult*) lMCInput->FindObject(lConfiguration.Data())->Clone( lMCName.Data() ) ;
    if(lVerbose) lMCResult->Print();
    
    //Compatibility check 1: check if generated with the same cuts
    if( !lDataResult->HasSameCuts(lMCResult) ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(" You are attempting to correct data with MC results obtained");
        AliWarning(" with different selection criteria! Will stop here. Please ");
        AliWarning(" check your input files!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //_________________________________________________
    // Process TH3Fs and expand into histograms of interest
    TH3F *f3dHistMC = (TH3F*) lMCResult->GetHistogram()->Clone("f3dHistMC");
    
    //Check if multiplicity interval requested is possible
    Bool_t lCheckMultMC = CheckCompatibleMultiplicity ( f3dHistMC );
    if ( !lCheckMultMC ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning("Requested mult interval is inconsistent with the MC provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Check if pt interval requested is possible
    Bool_t lCheckPtMC = CheckCompatiblePt ( f3dHistMC );
    if ( !lCheckPtMC ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(" Requested pt binning is inconsistent with the MC provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Project into relevant invariant mass histograms
    TH1D *lHistoMC[100];
    
    for( Long_t ibin = 0; ibin<lNPtBins; ibin++){
        lHistoMC[ibin] = f3dHistMC->ProjectionZ( Form("lHistoMC_%.1f_%.1f", lPtBins[ibin], lPtBins[ibin+1]),
                                                    f3dHistMC->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                                    f3dHistMC->GetXaxis()->FindBin( lHiMult-1e-5 ),
                                                    f3dHistMC->GetYaxis()->FindBin( lPtBins[ibin  ]+1e-5 ),
                                                    f3dHistMC->GetYaxis()->FindBin( lPtBins[ibin+1]-1e-5 )
                                                    );
        lHistoMC[ibin]->SetDirectory(0);
        fListMC->Add(lHistoMC[ibin]);
    }
    
    //Perform actual signal extraction
    Double_t lSignalVsPtMC[100], lSignalErrVsPtMC[100];
    TH1D* fHistRawVsPtMC  = new TH1D("fHistRawVsPtMC", "",lNPtBins,lPtBins);
    fHistRawVsPtMC->SetDirectory(0);
    
    if(!lVerbose) cout<<"AliStrangenessModule -> Extracting signal (mc):   ["<<flush;
    lExtStatus = kTRUE;
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        if(lVerbose) AliWarning(Form("MC Signal extraction on: %s",lHistoData[ibin]->GetName()));
        Double_t lBackgroundMC = 0;
        Double_t lBackgroundErrorMC = 0;
        lExtStatus = PerformSignalExtraction( lHistoMC[ibin], lSignalVsPtMC[ibin], lSignalErrVsPtMC[ibin], lBackgroundMC, lBackgroundErrorMC, lMeanVsPt[ibin], lSigmaVsPt[ibin], fListMC, "MC" );
        if( !lVerbose ){
            //Report errors with "!"
            if( lExtStatus ){ cout<<"="<<flush; } else { cout<<"!"<<flush; }
        }
        fHistRawVsPtMC->SetBinContent(ibin+1, lSignalVsPtMC[ibin] );
        fHistRawVsPtMC->SetBinError(ibin+1, lSignalErrVsPtMC[ibin] );
    }
    if(!lVerbose) cout<<"] Done!"<<endl;
    
    //Add raw spectra as main analysis output
    fListOutput->Add(fHistRawVsPtMC);
    
    //Acquire generated histogram
    TString lGenObjName = "fHistGeneratedPtVsYVsCentrality";
    lGenObjName.Append( lDataResult->GetParticleName() ) ;
    TH3F *f3dHistGenMC = (TH3F*) lMCCountersInput->FindObject( lGenObjName.Data() )->Clone("f3dHistGenMC");
    f3dHistGenMC->Sumw2();
    //Project this into a 1D histogram, please
    TH1D* fHistGeneratedOriginal = f3dHistGenMC -> ProjectionX( "fHistGeneratedOriginal",
                                                       f3dHistGenMC->GetYaxis()->FindBin(lDataResult->GetCutMinRapidity()+1e-5),
                                                       f3dHistGenMC->GetYaxis()->FindBin(lDataResult->GetCutMaxRapidity()-1e-5),
                                                       f3dHistGenMC->GetZaxis()->FindBin( lLoMult+1e-5 ),
                                                       f3dHistGenMC->GetZaxis()->FindBin( lHiMult-1e-5 )
                                                       );
    fHistGeneratedOriginal->SetDirectory(0);
    
    //Save this as a very relevant histogram
    fListOutput->Add(fHistGeneratedOriginal);
    
    //Rebin
    TH1D* fHistGenerated = (TH1D*) fHistGeneratedOriginal->Rebin( lNPtBins, "fHistGenerated", lPtBins );
    fHistGenerated->SetDirectory(0);
    
    //Save this as a very relevant histogram
    fListOutput->Add(fHistGeneratedOriginal);
    fListOutput->Add(fHistGenerated);
    
    //Generate Efficiency Histogram
    TH1D* fHistEfficiency = (TH1D*) fHistRawVsPtMC -> Clone ("fHistEfficiency") ;
    fHistEfficiency->SetDirectory(0);
    fHistEfficiency->Divide(fHistGenerated);
    
    //Check if g3/f correction enabled 
    if( lUseGeant3FlukaCorrection ){
        //=====================================================================================
        //Step 1: check which correction function is needed
        TF1 *lFuncG3FCorr = 0x0;
        TString lPartName = lDataResult->GetParticleName();
        if( lPartName.Contains("Lambda")||lPartName.Contains("XiMinus")||lPartName.Contains("OmegaMinus") ){
            //Initialize correction function for protons...
            lFuncG3FCorr = new TF1("lFuncG3FCorr", "1 - [0]*TMath::Exp([1]*x) + [2]",     0.25, 10.0);
            lFuncG3FCorr->SetParameter(0, 4.11235e+03);
            lFuncG3FCorr->SetParameter(1,-3.28947e+01);
            lFuncG3FCorr->SetParameter(2,-9.38341e-03);
        }
        if( lPartName.Contains("AntiLambda")||lPartName.Contains("XiPlus")||lPartName.Contains("OmegaPlus") ){
            //Initialize correction function for antiprotons...
            lFuncG3FCorr = new TF1("lFuncG3FCorr", "1 - [0]*TMath::Exp([1]*x) + [2] + [3]*1/TMath::Power(x, 0.2)*TMath::Log(x)", 0.25, 10.0);
            lFuncG3FCorr->SetParameter(0, 1.77339e+02);
            lFuncG3FCorr->SetParameter(1,-2.20242e+01);
            lFuncG3FCorr->SetParameter(2,-6.53769e-02);
            lFuncG3FCorr->SetParameter(3, 4.43007e-02);
        }
        if( !lFuncG3FCorr ) {
            AliWarning("Something went wrong with the determination of the G3/F correction!"); return 0x0;
        }
        fListOutput->Add(lFuncG3FCorr); //add function to output for completeness
        //=====================================================================================
        //the geant3/fluka logic:
        // --- in geant3: excessive antiparticle yield due to underestimated efficiencies
        // --- correction therefore has to increase efficiencies
        //     -> divide eff by lFuncG3FCorr at appropriate pT
        //=====================================================================================
        
        TH1D *fHistG3FCorrection = new TH1D("fHistG3FCorrection", "", lNPtBins, lPtBins);
        TProfile * fProfProtonComplete = (TProfile*) lMCResult->GetProtonProfileToCopy()->Clone("fProfProtonComplete");
        
        //rebin to match
        TProfile *fProfProton = (TProfile*) fProfProtonComplete->Rebin( lNPtBins, "fProfProton", lPtBins );
        
        fListOutput -> Add(fProfProtonComplete);
        fListOutput -> Add(fProfProton);
        
        for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
            Double_t lProtonMomentum = 1.0;
            //Temporary: fraction of bin center (to be fixed once proton profile is OK)
            lProtonMomentum = fHistG3FCorrection->GetBinCenter(ibin+1) * 0.75;
            
            //Final: use proton profile information from MC (commented for now)
            //lProtonMomentum = fProfProton->GetBinContent(ibin+1);
            fHistG3FCorrection -> SetBinContent(ibin+1, lFuncG3FCorr->Eval(lProtonMomentum) );
            
            //Apply correction on a bin-by-bin basis
            fHistEfficiency->SetBinContent(ibin+1, fHistEfficiency->GetBinContent(ibin+1) / lFuncG3FCorr->Eval(lProtonMomentum) );
            fHistEfficiency->SetBinError  (ibin+1, fHistEfficiency->GetBinError(ibin+1)   / lFuncG3FCorr->Eval(lProtonMomentum) );
        }
        
        //Store histogram with the correction
        fListOutput->Add(fHistG3FCorrection);
    }
    
    fListOutput->Add(fHistEfficiency);
    
    //Generate Corrected Spectrum
    TH1D *fHistSpectra = (TH1D*) fHistRawVsPt->Clone("fHistSpectra");
    fHistSpectra->SetDirectory(0);
    
    //Efficiency correction
    fHistSpectra->Divide(fHistEfficiency);
    
    //Normalize by number of events
    fHistSpectra->Scale(1.0/lNEvents, "width");
    
    //Scale with rapidity window size (default: no scale)
    fHistSpectra->Scale(1.0/(lDataResult->GetCutMaxRapidity()-lDataResult->GetCutMinRapidity()));
    
    //Add to output 
    fListOutput->Add(fHistSpectra);
    
    TH1D *fHistSpectraToReturn = (TH1D*) fHistSpectra->Clone("fHistSpectraToReturn");
    fHistSpectraToReturn->SetDirectory(0);
    
    fFileOut->cd();
    //Save all objects owned by the TLists
    fListData  ->Write("cListData",       TObject::kSingleKey);
    fListMC    ->Write("cListMC",         TObject::kSingleKey);
    fListOutput->Write("cAnalysisOutput", TObject::kSingleKey);
    fFileOut->Write();
    fFileOut->Close();
    
    timer->Stop();
    cout<<"AliStrangenessModule -> Done! Total time: "<<timer->RealTime()<<"s (real), "<< timer->CpuTime() <<"s (cpu)"<<endl;
    delete timer;
}

//________________________________________________________________
Bool_t AliStrangenessModule::CheckCompatibleMultiplicity( TH3F *lHisto ){
    Bool_t lLoMultMatch = kFALSE;
    Bool_t lHiMultMatch = kFALSE;
    for( Long_t ibin = 1; ibin<lHisto->GetNbinsX()+2; ibin++){
        Double_t lLowEdge = lHisto->GetXaxis()->GetBinLowEdge(ibin);
        if( TMath::Abs( lLowEdge - lLoMult) < 1e-5 ) lLoMultMatch = kTRUE;
        if( TMath::Abs( lLowEdge - lHiMult) < 1e-5 ) lHiMultMatch = kTRUE;
    }
    Bool_t lReturnValue = kFALSE;
    if ( lLoMultMatch && lHiMultMatch ) lReturnValue = kTRUE;
    return lReturnValue;
}

//________________________________________________________________
Bool_t AliStrangenessModule::CheckCompatiblePt( TH3F *lHisto ){
    Bool_t lMatchFound = kFALSE;
    Bool_t lReturnValue = kTRUE;
    for( Long_t iptbin = 0; iptbin<lNPtBins+1; iptbin++){
        lMatchFound = kFALSE;
        for( Long_t ibin = 1; ibin<lHisto->GetNbinsY()+2; ibin++){
            Double_t lLowEdge = lHisto->GetYaxis()->GetBinLowEdge(ibin);
            if( TMath::Abs( lLowEdge - lPtBins[iptbin]) < 1e-5 ) lMatchFound = kTRUE;
        }
        if ( ! lMatchFound ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}


//________________________________________________________________
Bool_t AliStrangenessModule::PerformInitialFit( TH1D *lHisto, Double_t &lMean, Double_t &lMeanErr, Double_t &lSigma, Double_t &lSigmaErr, TList *lControlList ){
    //Helper function to perform initial gaussian + linear fit
    
    Bool_t lReturnValue = kTRUE; //everything ok = kTRUE
    
    //Will expect to have mass received as mean value here
    Double_t lMass = lMean;
    TString lName = lHisto->GetName();
    lName.Append("_InitialFit");
    TF1 *fit = new TF1(lName.Data(),"[0]+[1]*x+[2]*TMath::Gaus(x, [3], [4])", lMean+lOffsetFromMeanLeft, lMean+lOffsetFromMeanRight);
    
    //Guess linear parameters
    Double_t lAverageBg = 0.5*(lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean+lOffsetFromMeanLeft) )
                               +lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean+lOffsetFromMeanRight) ));
    Double_t lGuessedSlope = (lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean-0.010) ) -
                              lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean+0.010) ) ) /
                            ( 0.020 ) ;
    
    //constant at zero: lAverageBg - lMean*lGuessedSlope
    fit->SetParameter(0, lAverageBg - lMean*lGuessedSlope);
    fit->SetParameter(1, lGuessedSlope);
    
    //Guess Gaussian Parameters
    fit->SetParameter(2, lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean) ) - lAverageBg );
    fit->SetParLimits(2, 1, 10*(lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean) ) - lAverageBg));
    fit->SetParameter(3, lMean);
    fit->SetParLimits(3, lMean-0.006,lMean+0.006); //ALICE is never off by that much
    fit->SetParameter(4,0.003);
    fit->SetParLimits(4,0.001, 0.015); //ALICE tracking is (usually) never worse than this
    
    //Fit options
    TString lFitOptions = "IREM0S";
    if (!lVerbose) lFitOptions.Append("Q") ;
    
    //Printout options used for the fit
    if ( lVerbose ) cout<<"Fit to histogram "<<lHisto->GetName()<<" will be carried out with options: "<<lFitOptions.Data()<<endl;
    
    //Perform fit
    TFitResultPtr fitResultPtr = lHisto->Fit( lName.Data(), lFitOptions.Data());
    if ( !fitResultPtr->IsValid() ) lReturnValue = kFALSE; //went bad
    
    //Store to control output
    lControlList->Add(fit);
    
    //Provide output variables
    lMean     = fit->GetParameter(3);
    lMeanErr  = fit->GetParError(3);
    lSigma    = fit->GetParameter(4);
    lSigmaErr = fit->GetParError(4);
    
    return lReturnValue;
}

//________________________________________________________________
Bool_t AliStrangenessModule::PerformSignalExtraction( TH1D *lHisto, Double_t &lSignal, Double_t &lSignalErr, Double_t &lBackground, Double_t &lBackgroundError, Double_t lMean, Double_t lSigma, TList *lControlList, TString lOption ){
    //Helper function to perform actual signal extraction
    
    Bool_t lReturnValue = kTRUE; //everything went alright -> kTRUE
    
    TString lFitOptions = "R0S";
    if (!lVerbose) lFitOptions.Append("Q");
    
    //Find bins in which signal extraction is to be performed
    Long_t lBinPeakLo = lHisto->GetXaxis()->FindBin ( lMean + lLoPeak*lSigma );
    Long_t lBinPeakHi = lHisto->GetXaxis()->FindBin ( lMean + lHiPeak*lSigma );
    Long_t lBinLeftBgLo = lHisto->GetXaxis()->FindBin ( lMean + lLoLeftBg*lSigma );
    Long_t lBinLeftBgHi = lHisto->GetXaxis()->FindBin ( lMean + lHiLeftBg*lSigma );
    Long_t lBinRightBgLo = lHisto->GetXaxis()->FindBin ( lMean + lLoRightBg*lSigma );
    Long_t lBinRightBgHi = lHisto->GetXaxis()->FindBin ( lMean + lHiRightBg*lSigma );
    
    //Inclusive on lower and upper limits
    //Get values and use these values for fit ranges: meant to harmonize bin counting wrt fitting
    Double_t lValPeakLo = lHisto->GetBinLowEdge( lBinPeakLo   );
    Double_t lValPeakHi = lHisto->GetBinLowEdge( lBinPeakHi+1 );
    Double_t lValLeftBgLo = lHisto->GetBinLowEdge( lBinLeftBgLo   );
    Double_t lValLeftBgHi = lHisto->GetBinLowEdge( lBinLeftBgHi+1 );
    Double_t lValRightBgLo = lHisto->GetBinLowEdge( lBinRightBgLo   );
    Double_t lValRightBgHi = lHisto->GetBinLowEdge( lBinRightBgHi+1 );
    
    //Get very first guess for linear background
    Double_t lAverageBg = lHisto->Integral(lBinLeftBgLo,lBinLeftBgHi )+lHisto->Integral(lBinRightBgLo,lBinRightBgHi) ;
    lAverageBg = lAverageBg / ( lBinLeftBgHi - lBinLeftBgLo + lBinRightBgHi - lBinRightBgLo + 2);
    
    Double_t lLeftY = lHisto->Integral(lBinLeftBgLo, lBinLeftBgHi  ) / (lBinLeftBgHi - lBinLeftBgLo + 1 );
    Double_t lLeftX = 0.5*(lValLeftBgHi+lValLeftBgLo );
    
    Double_t lRightY = lHisto->Integral(lBinRightBgLo,lBinRightBgHi ) / (lBinRightBgHi - lBinRightBgLo + 1 );
    Double_t lRightX = 0.5*(lValRightBgHi+lValRightBgLo );
    
    Double_t lGuessedSlope = (lRightY-lLeftY)/(lRightX-lLeftX);
    
    Double_t lBgConst = lAverageBg - lMean*lGuessedSlope;
    Double_t lBgSlope = lGuessedSlope;
    
    //Check if this is a low-statistics bin and if so use either "L" or "LL" fit options
    lFitOptions.Append( GetGoodFitOption ( lHisto, lBinLeftBgLo, lBinRightBgHi ) ); 
    
    if ( lVerbose ) cout<<"Fit to histogram "<<lHisto->GetName()<<" will be carried out with options: "<<lFitOptions.Data()<<endl; 
    
    //Clone histogram and create a control histogram with the peak ranges highlighted
    TString lNameHistoPeak = lHisto->GetName();
    lNameHistoPeak.Append("_Peak");
    TH1D *lHistoPeak = (TH1D*) lHisto->Clone(lNameHistoPeak.Data());
    lHistoPeak->SetDirectory(0);
    lHistoPeak->Reset();
    for(Long_t ibin=lBinPeakLo; ibin<lBinPeakHi+1; ibin++){
        lHistoPeak->SetBinContent(ibin, lHisto->GetBinContent(ibin));
        lHistoPeak->SetBinError(ibin, lHisto->GetBinError(ibin));
    }
    lControlList->Add(lHistoPeak);
    
    //Pointer to fit (if it exists)
    TF1 *fit           = 0x0;
    TF1 *fitToSubtract = 0x0;
    
    //N.B. if option "MC" is chosen, background will be assumed zero! 
    Double_t lBgEstimate = 0;
    Double_t lBgEstimateError = 0;
    
    if( lOption.Contains("linear") ){
        //Step 1: Perform Linear Fit to start with
        TString lName = lHisto->GetName();
        lName.Append("_LinearFit");
        fit = new TF1(lName.Data(), this, &AliStrangenessModule::BgPol1,
                      lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma, 4 , "AliStrangenessModule", "BgPol1");
        
        //Start with parameters from initial fit: probably a good initial guess
        fit->FixParameter(0, lMean + lHiLeftBg*lSigma );
        fit->FixParameter(1, lMean + lLoRightBg*lSigma);
        fit->SetParameter(2, lBgConst);
        fit->SetParameter(3, lBgSlope);
        
        //Perform fit
        
        TFitResultPtr fitResultPtr = lHisto->Fit( lName.Data(), lFitOptions.Data());
        if ( !fitResultPtr->IsValid() ) lReturnValue = kFALSE; //went bad
        
        TString lNameToSubtract = lHisto->GetName();
        lNameToSubtract.Append("_FitToSubtract");
        fitToSubtract = new TF1(lNameToSubtract.Data(), "[0]+[1]*x",
                      lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma);
        fitToSubtract->SetParameter( 0, fit->GetParameter(2) );
        fitToSubtract->SetParameter( 1, fit->GetParameter(3) );
        
        lBgEstimate      = fitToSubtract->Integral     ( lValPeakLo, lValPeakHi );
        lBgEstimate     /= lHisto->GetBinWidth(lBinPeakLo); //Transform into counts!
        lBgEstimateError = TMath::Sqrt(lBgEstimate); //fit->IntegralError( lValPeakLo, lValPeakHi );
    }
    
    if ( lOption.Contains("quadratic") ){
        //Step 2: Perform Quadratic Fit to improve results
        TString lNameQuad = lHisto->GetName();
        lNameQuad.Append("_QuadraticFit");
        fit = new TF1(lNameQuad, this, &AliStrangenessModule::BgPol2,
                               lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma, 5 , "AliStrangenessModule", "BgPol2");
        
        //Start with parameters from initial fit: probably a good initial guess
        fit->FixParameter(0, lMean + lHiLeftBg*lSigma );
        fit->FixParameter(1, lMean + lLoRightBg*lSigma);
        fit->SetParameter(2, lBgConst );
        fit->SetParameter(3, lBgSlope );
        fit->SetParameter(4, 0.000 );
        
        //Perform fit - otherwise stick to initial (linear) guess
        TFitResultPtr fitResultPtr = lHisto->Fit( lNameQuad.Data(), lFitOptions.Data());
        if ( !fitResultPtr->IsValid() ) lReturnValue = kFALSE; //went bad
        
        TString lNameToSubtract = lHisto->GetName();
        lNameToSubtract.Append("_FitToSubtract");
        fitToSubtract = new TF1(lNameToSubtract.Data(), "[0]+[1]*x+[2]*x*x",
                                lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma);
        fitToSubtract->SetParameter( 0, fit->GetParameter(2) );
        fitToSubtract->SetParameter( 1, fit->GetParameter(3) );
        fitToSubtract->SetParameter( 2, fit->GetParameter(4) );
        
        lBgEstimate      = fitToSubtract->Integral     ( lValPeakLo, lValPeakHi );
        lBgEstimate     /= lHisto->GetBinWidth(lBinPeakLo); //Transform into counts!
        lBgEstimateError = TMath::Sqrt(lBgEstimate); //fit->IntegralError( lValPeakLo, lValPeakHi );
    }
    
    if ( lOption.Contains("cubic") ){
        //Step 2: Perform Quadratic Fit to improve results
        TString lNameCubic = lHisto->GetName();
        lNameCubic.Append("_CubicFit");
        fit = new TF1(lNameCubic, this, &AliStrangenessModule::BgPol3,
                      lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma, 6 , "AliStrangenessModule", "BgPol3");
        
        //Start with parameters from initial fit: probably a good initial guess
        fit->FixParameter(0, lMean + lHiLeftBg*lSigma );
        fit->FixParameter(1, lMean + lLoRightBg*lSigma);
        fit->SetParameter(2, lBgConst );
        fit->SetParameter(3, lBgSlope );
        fit->SetParameter(4, 0.000 );
        fit->SetParameter(5, 0.000 );
        
        //Perform fit - otherwise stick to initial (linear) guess
        TFitResultPtr fitResultPtr = lHisto->Fit( lNameCubic.Data(), lFitOptions.Data());
        if ( !fitResultPtr->IsValid() ) lReturnValue = kFALSE; //went bad
        
        TString lNameToSubtract = lHisto->GetName();
        lNameToSubtract.Append("_FitToSubtract");
        fitToSubtract = new TF1(lNameToSubtract.Data(), "[0]+[1]*x+[2]*x*x+[3]*x*x*x",
                                lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma);
        fitToSubtract->SetParameter( 0, fit->GetParameter(2) );
        fitToSubtract->SetParameter( 1, fit->GetParameter(3) );
        fitToSubtract->SetParameter( 2, fit->GetParameter(4) );
        fitToSubtract->SetParameter( 3, fit->GetParameter(5) );
        
        lBgEstimate      = fitToSubtract->Integral     ( lValPeakLo, lValPeakHi );
        lBgEstimate     /= lHisto->GetBinWidth(lBinPeakLo); //Transform into counts!
        lBgEstimateError = TMath::Sqrt(lBgEstimate); //fit->IntegralError( lValPeakLo, lValPeakHi );
    }
    
    if ( lOption.Contains("bincounting") || lOption.Contains("MC") ){
        //Check if possible
        if((lLoLeftBg + lHiRightBg > 1e-6 ||
            lHiLeftBg + lLoRightBg > 1e-6 ||
            lLoPeak   + lHiPeak    > 1e-6 ) && !lOption.Contains("MC") ){
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            AliWarning(" Cannot perform bin counting with asymmetric peak and bg regions!");
            AliWarning("   In this case, please prefer the \'linear\' sig. ext. option");
            AliWarning("                  Will not produce spectra");
            AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
            lSignal = 0.0;
            lSignalErr = 0.0;
            return kFALSE;
        }

        //Sum up yields in corresponding bins
        for( Long_t ibin=lBinLeftBgLo; ibin<lBinLeftBgHi+1; ibin++)
            lBgEstimate += lHisto->GetBinContent(ibin);
        for( Long_t ibin=lBinRightBgLo; ibin<lBinRightBgHi+1; ibin++)
            lBgEstimate += lHisto->GetBinContent(ibin);
        
        //Error: sqrt(counts)
        //lBgEstimateError = TMath::Sqrt(lBgEstimate);
        
        //Error: Sum of errors
        lBgEstimateError = 0;
        for( Long_t ibin=lBinLeftBgLo; ibin<lBinLeftBgHi+1; ibin++)
            lBgEstimateError += lHisto->GetBinError(ibin)*lHisto->GetBinError(ibin);
        for( Long_t ibin=lBinRightBgLo; ibin<lBinRightBgHi+1; ibin++)
            lBgEstimateError += lHisto->GetBinError(ibin)*lHisto->GetBinError(ibin);
        lBgEstimateError = TMath::Sqrt(lBgEstimateError);
        
        //Scale according to number of bins
        Double_t lNBinsSummed = lBinLeftBgHi - lBinLeftBgLo + lBinRightBgHi - lBinRightBgLo + 2;
        Double_t lNBinsPeak   = lBinPeakHi - lBinPeakLo + 1;
        Double_t lScalingFactor = lNBinsPeak / lNBinsSummed;
        
        //Scale everything to match the current area of the peak
        lBgEstimate      = lBgEstimate      * lScalingFactor;
        lBgEstimateError = lBgEstimateError * lScalingFactor;
    }
    
    //Save a fit if there is one
    if ( fit           ) lControlList->Add(fit          );
    if ( fitToSubtract ) lControlList->Add(fitToSubtract);
    
    //Provide background information to the outside scope
    lBackground      = lBgEstimate;
    lBackgroundError = lBgEstimateError;
    
    Double_t lPeakPlusBg      = 0;
    Double_t lPeakPlusBgError = 0;
    lPeakPlusBg = lHisto->IntegralAndError(lBinPeakLo,lBinPeakHi,lPeakPlusBgError);
    //lPeakPlusBgError = TMath::Sqrt(lPeakPlusBg);
    
    lSignal = lPeakPlusBg - lBgEstimate;
    lSignalErr = TMath::Sqrt( lPeakPlusBgError*lPeakPlusBgError + lBgEstimateError*lBgEstimateError );
    
    return lReturnValue;
}

//________________________________________________________________
TString AliStrangenessModule::GetGoodFitOption( TH1D *lHisto , Int_t ilow, Int_t ihigh )
{
    //Function for guessing good fit options for the histo lHisto
    //(between bins numbered ilow, ihigh
    TString lResult = "";
    for(Int_t ibin=ilow; ibin<ihigh+1; ibin++){
        if( lHisto->GetBinContent(ibin) < 10 ) lResult = "L";
        //if( lHisto->GetBinContent(ibin) <  20 ) lResult = "LL";
    }
    return lResult;
}

//________________________________________________________________
Double_t AliStrangenessModule::BgPol1(const Double_t *x, const Double_t *par)
{
    //Function for background fitting, rejects peak region
    //Parameter [0] -> Hi LeftBg Boundary
    //Parameter [1] -> Lo RightBg Boundary
    if ( x[0] > par[0] && x[0] < par[1]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[2] + par[3]*x[0];
}

//________________________________________________________________
Double_t AliStrangenessModule::BgPol2(const Double_t *x, const Double_t *par)
{
    //Function for background fitting, rejects peak region
    //Parameter [0] -> Hi LeftBg Boundary
    //Parameter [1] -> Lo RightBg Boundary
    if ( x[0] > par[0] && x[0] < par[1]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[2] + par[3]*x[0] + par[4]*x[0]*x[0];
}

//________________________________________________________________
Double_t AliStrangenessModule::BgPol3(const Double_t *x, const Double_t *par)
{
    //Function for background fitting, rejects peak region
    //Parameter [0] -> Hi LeftBg Boundary
    //Parameter [1] -> Lo RightBg Boundary
    if ( x[0] > par[0] && x[0] < par[1]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[2] + par[3]*x[0] + par[4]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0];
}





