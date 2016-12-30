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
lMCInput(0x0),
lMCCountersInput(0x0),
lLoMult(0.0),
lHiMult(90.0),
lNPtBins(0),
lLoLeftBg(-6),
lHiLeftBg(-4),
lLoPeak(-3),
lHiPeak(3),
lLoRightBg(4),
lHiRightBg(7),
lSigExtTech("linear")
{
    // Dummy Constructor - not to be used!
    for(Long_t ibin = 0; ibin<100; ibin++) lPtBins[ibin] = 0;
}

//________________________________________________________________
AliStrangenessModule::AliStrangenessModule(const char * name, const char * title) :
TNamed(name,title),
lDataInput(0x0),
lDataCountersInput(0x0),
lMCInput(0x0),
lMCCountersInput(0x0),
lLoMult(0.0),
lHiMult(90.0),
lNPtBins(0),
lLoLeftBg(-6),
lHiLeftBg(-4),
lLoPeak(-3),
lHiPeak(3),
lLoRightBg(4),
lHiRightBg(7),
lSigExtTech("linear")
{
    // Main constructor
    for(Long_t ibin = 0; ibin<100; ibin++) lPtBins[ibin] = 0; 
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
TH1D* AliStrangenessModule::DoAnalysis( TString lConfiguration ){
    //Main analysis code
    
    //"Die hard" fitting, please
    TVirtualFitter::SetMaxIterations( 50000 ); 
    
    //Check if input has been provied
    if( !lDataInput ) {
        AliWarning("Data input not defined!"); return 0x0;
    }
    if( !lDataCountersInput ) {
        AliWarning("Data counters list not defined!"); return 0x0;
    }
    if( !lMCInput ) {
        AliWarning("MC input not defined!"); return 0x0;
    }
    if( !lMCCountersInput ) {
        AliWarning("MC counters list not defined!"); return 0x0;
    }
    
    //Open File for output
    TFile *fFileOut = new TFile(Form("Results_%s.root",lConfiguration.Data() ), "RECREATE");
    TList *fListData = new TList();
    fListData->SetName("cListData");
    fListData->SetOwner(kTRUE);

    TList *fListMC = new TList();
    fListMC->SetName("cListMC");
    fListMC->SetOwner(kTRUE);
    
    TList *fListOutput = new TList();
    fListOutput->SetName("cAnalysisOutput");
    fListOutput->SetOwner(kTRUE);
    
    //Step 1: Open real data object
    AliVWeakResult *lDataResult = (AliVWeakResult*) lDataInput->FindObject(lConfiguration.Data());
    TString lDataName = lDataResult->GetName();
    lDataName.Append("_Data");
    lDataResult->SetName( lDataName.Data() );
    lDataResult->Print();
    
    //_________________________________________________
    // Process TH3Fs and expand into histograms of interest
    TH3F *f3dHistData = lDataResult->GetHistogram();
    f3dHistData->SetName("f3dHistData");
    
    //Check if multiplicity interval requested is possible
    Bool_t lCheckMult = CheckCompatibleMultiplicity ( f3dHistData );
    if ( !lCheckMult ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning("Requested mult interval is inconsistent with the data provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Check if pt interval requested is possible
    Bool_t lCheckPt = CheckCompatiblePt ( f3dHistData );
    if ( !lCheckPt ){
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        AliWarning(" Requested pt binning is inconsistent with the data provided!");
        AliWarning("!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!");
        return 0x0;
    }
    
    //Project into relevant invariant mass histograms
    TH1D *lHistoData[100];
    
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
    
    //Perform early fit with gaussian+linear for acquisition of fit range
    Double_t lMeanVsPt[100], lMeanErrVsPt[100];
    Double_t lSigmaVsPt[100], lSigmaErrVsPt[100];
    Double_t lBgConstVsPt[100], lBgSlopeVsPt[100];
    TH1D* fHistMeanVsPt  = new TH1D("fHistMeanVsPt", "",lNPtBins,lPtBins);
    TH1D* fHistSigmaVsPt = new TH1D("fHistSigmaVsPt","",lNPtBins,lPtBins);
    fHistMeanVsPt->SetDirectory(0);
    fHistSigmaVsPt->SetDirectory(0);
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        lMeanVsPt[ibin] = lDataResult->GetMass(); //Initial guess
        AliWarning(Form("Initial fit on: %s",lHistoData[ibin]->GetName()));
        PerformInitialFit( lHistoData[ibin], lMeanVsPt[ibin], lMeanErrVsPt[ibin], lSigmaVsPt[ibin], lSigmaErrVsPt[ibin], lBgConstVsPt[ibin], lBgSlopeVsPt[ibin], fListData );
        fHistMeanVsPt->SetBinContent(ibin+1, lMeanVsPt[ibin] );
        fHistMeanVsPt->SetBinError(ibin+1, lMeanErrVsPt[ibin] );
        fHistSigmaVsPt->SetBinContent(ibin+1, lSigmaVsPt[ibin] );
        fHistSigmaVsPt->SetBinError(ibin+1, lSigmaErrVsPt[ibin] );
    }
    
    //FIXME: Invoke TH1D::Smooth here if asked to do so 
    
    fListData->Add( fHistMeanVsPt );
    fListData->Add( fHistSigmaVsPt );
    
    //Perform actual signal extraction
    Double_t lSignalVsPt[100], lSignalErrVsPt[100];
    TH1D* fHistRawVsPt  = new TH1D("fHistRawVsPt", "",lNPtBins,lPtBins);
    fHistRawVsPt->SetDirectory(0);
    
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        AliWarning(Form("Signal extraction on: %s",lHistoData[ibin]->GetName()));
        PerformSignalExtraction( lHistoData[ibin], lSignalVsPt[ibin], lSignalErrVsPt[ibin], lMeanVsPt[ibin], lSigmaVsPt[ibin], lBgConstVsPt[ibin], lBgSlopeVsPt[ibin], fListData, lSigExtTech.Data() );
        fHistRawVsPt->SetBinContent(ibin+1, lSignalVsPt[ibin] );
        fHistRawVsPt->SetBinError(ibin+1, lSignalErrVsPt[ibin] );
    }
    
    //Add raw spectra as main analysis output
    fListOutput->Add(fHistRawVsPt);
    
    //Step N: Open MC data object
    AliVWeakResult *lMCResult = (AliVWeakResult*) lMCInput->FindObject(lConfiguration.Data());
    TString lMCName = lMCResult->GetName();
    lMCName.Append("_MC");
    lMCResult->SetName( lMCName.Data() );
    lMCResult->Print();
    
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
    TH3F *f3dHistMC = lMCResult->GetHistogram();
    f3dHistMC->SetName("f3dHistMC");
    
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
    
    for(Long_t ibin = 0; ibin<lNPtBins; ibin++){
        AliWarning(Form("MC Signal extraction on: %s",lHistoData[ibin]->GetName()));
        PerformSignalExtraction( lHistoMC[ibin], lSignalVsPtMC[ibin], lSignalErrVsPtMC[ibin], lMeanVsPt[ibin], lSigmaVsPt[ibin], 0, 0, fListMC, "MC" );
        fHistRawVsPtMC->SetBinContent(ibin+1, lSignalVsPtMC[ibin] );
        fHistRawVsPtMC->SetBinError(ibin+1, lSignalErrVsPtMC[ibin] );
    }
    
    //Add raw spectra as main analysis output
    fListOutput->Add(fHistRawVsPtMC);
    
    //Acquire generated histogram
    TString lGenObjName = "fHistGeneratedPtVsYVsCentrality";
    lGenObjName.Append( lDataResult->GetParticleName() ) ;
    TH3F *f3dHistGenMC = (TH3F*) lMCCountersInput->FindObject( lGenObjName.Data() );
    f3dHistGenMC->Sumw2();
    //Project this into a 1D histogram, please
    TH1D* fHistGeneratedOriginal = f3dHistGenMC -> ProjectionX( "fHistGeneratedOriginal",
                                                       f3dHistGenMC->GetYaxis()->FindBin(-0.5+1e-5),
                                                       f3dHistGenMC->GetYaxis()->FindBin(+0.5-1e-5),
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
    fListOutput->Add(fHistEfficiency);
    
    //Get number of events from data
    TH1D *fHistCentrality = (TH1D*) lDataCountersInput->FindObject( "fHistCentrality" );
    fHistCentrality->SetDirectory(0);
    Double_t lNEvents = fHistCentrality->Integral( fHistCentrality->GetXaxis()->FindBin( lLoMult+1e-5 ),
                                                  fHistCentrality->GetXaxis()->FindBin( lHiMult-1e-5 ) );
    AliWarning( Form("Number of events in this multiplicity class: %.0f", lNEvents) );
    
    //Generate Corrected Spectrum
    TH1D *fHistSpectra = (TH1D*) fHistRawVsPt->Clone("fHistSpectra");
    fHistSpectra->SetDirectory(0);
    fHistSpectra->Divide(fHistEfficiency);
    fHistSpectra->Scale(1.0/lNEvents, "width");
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
void AliStrangenessModule::PerformInitialFit( TH1D *lHisto, Double_t &lMean, Double_t &lMeanErr, Double_t &lSigma, Double_t &lSigmaErr, Double_t &lBgConst, Double_t &lBgSlope, TList *lControlList ){
    //Helper function to perform initial gaussian + linear fit
    //Will expect to have mass received as mean value here
    Double_t lMass = lMean;
    TString lName = lHisto->GetName();
    lName.Append("_InitialFit");
    TF1 *fit = new TF1(lName.Data(),"[0]+[1]*x+[2]*TMath::Gaus(x, [3], [4])", lMean-0.025, lMean+0.025);
    //Guess linear parameters
    Double_t lAverageBg = 0.5*(lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean-0.015) )
                             +lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean+0.015) ));
    Double_t lGuessedSlope = (lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean+0.020) ) -
                              lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean-0.020) ) ) / 0.040;
    
    //constant at zero: lAverageBg - lMean*lGuessedSlope
    fit->SetParameter(0, lAverageBg - lMean*lGuessedSlope);
    fit->SetParameter(1, lGuessedSlope);
    
    //Guess Gaussian Parameters
    fit->SetParameter(2, lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean) ) - lAverageBg );
    fit->SetParLimits(2, 10, 10*(lHisto->GetBinContent(lHisto->GetXaxis()->FindBin(lMean) ) - lAverageBg));
    fit->SetParameter(3, lMean);
    fit->SetParLimits(3, lMean-0.006,lMean+0.006); //ALICE is never off by that much
    fit->SetParameter(4,0.003);
    fit->SetParLimits(4,0.001, 0.035); //ALICE tracking is (usually) never worse than this
    
    //Perform fit
    lHisto->Fit( lName.Data(), "IREM0");
    
    //Store to control output
    lControlList->Add(fit);
    
    //Provide output variables

    lBgConst = fit->GetParameter(0);
    lBgSlope = fit->GetParameter(1);
    lMean = fit->GetParameter(3);
    lMeanErr = fit->GetParError(3);
    lSigma = fit->GetParameter(4);
    lSigmaErr = fit->GetParError(4);
}

//________________________________________________________________
void AliStrangenessModule::PerformSignalExtraction( TH1D *lHisto, Double_t &lSignal, Double_t &lSignalErr, Double_t lMean, Double_t lSigma, Double_t lBgConst, Double_t lBgSlope, TList *lControlList, TString lOption ){
    //Helper function to perform actual signal extraction
    
    //Find bins in which signal extraction is to be performed
    Long_t lBinPeakLo = lHisto->GetXaxis()->FindBin ( lMean + lLoPeak*lSigma );
    Long_t lBinPeakHi = lHisto->GetXaxis()->FindBin ( lMean + lHiPeak*lSigma );
    
    //Inclusive on lower and upper limits
    Double_t lValPeakLo = lHisto->GetBinLowEdge( lBinPeakLo   );
    Double_t lValPeakHi = lHisto->GetBinLowEdge( lBinPeakHi+1 );
    
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
        lHisto->Fit( lName.Data(), "R0");
        
        TString lNameToSubtract = lHisto->GetName();
        lNameToSubtract.Append("_FitToSubtract");
        fitToSubtract = new TF1(lNameToSubtract.Data(), "[0]+[1]*x",
                      lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma);
        fitToSubtract->SetParameter( 0, fit->GetParameter(2) );
        fitToSubtract->SetParameter( 1, fit->GetParameter(3) );
        
        lBgEstimate      = fitToSubtract->Integral     ( lValPeakLo, lValPeakHi );
        lBgEstimate      /= lHisto->GetBinWidth(lBinPeakLo); //Transform into counts!
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
        lHisto->Fit( lNameQuad.Data(), "R0");
        
        TString lNameToSubtract = lHisto->GetName();
        lNameToSubtract.Append("_FitToSubtract");
        fitToSubtract = new TF1(lNameToSubtract.Data(), "[0]+[1]*x+[2]*x*x",
                                lMean + lLoLeftBg*lSigma, lMean + lHiRightBg*lSigma);
        fitToSubtract->SetParameter( 0, fit->GetParameter(2) );
        fitToSubtract->SetParameter( 1, fit->GetParameter(3) );
        fitToSubtract->SetParameter( 2, fit->GetParameter(4) );
        
        lBgEstimate      = fitToSubtract->Integral     ( lValPeakLo, lValPeakHi );
        lBgEstimate      /= lHisto->GetBinWidth(lBinPeakLo); //Transform into counts!
        lBgEstimateError = TMath::Sqrt(lBgEstimate); //fit->IntegralError( lValPeakLo, lValPeakHi );
    }
    
    if ( lOption.Contains("bincounting") ){
        //Find bins for LeftBg and RightBg
        Long_t lBinLeftBgLo = lHisto->GetXaxis()->FindBin ( lMean + lLoLeftBg*lSigma );
        Long_t lBinLeftBgHi = lHisto->GetXaxis()->FindBin ( lMean + lHiLeftBg*lSigma );
        Long_t lBinRightBgLo = lHisto->GetXaxis()->FindBin ( lMean + lLoRightBg*lSigma );
        Long_t lBinRightBgHi = lHisto->GetXaxis()->FindBin ( lMean + lHiRightBg*lSigma );
        
        //Sum up yields in corresponding bins
        for( Long_t ibin=lBinLeftBgLo; ibin<lBinLeftBgHi+1; ibin++)
            lBgEstimate += lHisto->GetBinContent(ibin);
        for( Long_t ibin=lBinRightBgLo; ibin<lBinRightBgHi+1; ibin++)
            lBgEstimate += lHisto->GetBinContent(ibin);
        
        //Error: sqrt(counts
        lBgEstimateError = TMath::Sqrt(lBgEstimate);
        
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
    
    Double_t lPeakPlusBg      = 0;
    Double_t lPeakPlusBgError = 0;
    lPeakPlusBg = lHisto->IntegralAndError(lBinPeakLo,lBinPeakHi,lPeakPlusBgError);
    lPeakPlusBgError = TMath::Sqrt(lPeakPlusBg);
    
    lSignal = lPeakPlusBg - lBgEstimate;
    lSignalErr = TMath::Sqrt( lPeakPlusBgError*lPeakPlusBgError + lBgEstimateError*lBgEstimateError );
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



