//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
// This is a container class for the correction factors for the hadronic 
// component of transverse energy
// It is filled by the output of AliAnalysisTaskHadEt from spinning over Monte 
// Carlo data (using AliAnalysisHadEtMonteCarlo)
//It is used by AliAnalysisTaskHadEt while spinning over reconstructed data 
// (using AliAnalysisHadEtReconstructed)
//Please see https://twiki.cern.ch/twiki/bin/view/ALICE/ETCaloAnalysis
#include "AliAnalysisHadEtCorrections.h"
#include "TMath.h"
#include <iostream>
#include "Rtypes.h"
#include "TNamed.h"

using namespace std;

ClassImp(AliAnalysisHadEtCorrections);


AliAnalysisHadEtCorrections::AliAnalysisHadEtCorrections() : TNamed(),
  fEtaCut(0)
  ,fAcceptanceCorrectionFull(0)
  ,fAcceptanceCorrectionEMCAL(0)
  ,fAcceptanceCorrectionPHOS(0)
  ,fNeutralCorrection(0)
  ,fNotHadronicCorrection(0)
  ,fpTcutCorrectionTPC(0)
  ,fpTcutCorrectionITS(0)
  ,fNeutralCorrectionLow(0)
  ,fNotHadronicCorrectionLow(0)
  ,ffpTcutCorrectionTPCLow(0)
  ,ffpTcutCorrectionITSLow(0)
  ,fNeutralCorrectionHigh(0)
  ,fNotHadronicCorrectionHigh(0)
  ,ffpTcutCorrectionTPCHigh(0)
  ,ffpTcutCorrectionITSHigh(0)
  ,fnotIDTPC(0)
  ,fnotIDITS(0)
  ,fnotIDNoID(0)
  ,fEfficiencyPionTPC(0)
  ,fEfficiencyKaonTPC(0)
  ,fEfficiencyProtonTPC(0)
  ,fEfficiencyHadronTPC(0)
  ,fEfficiencyPionITS(0)
  ,fEfficiencyKaonITS(0)
  ,fEfficiencyProtonITS(0)
  ,fEfficiencyHadronITS(0)
  ,fBackgroundTPC(0)
  ,fBackgroundITS(0)
{//default constructor
  //This seems to solve a compiler error
  cout<<"Creating new AliAnalysisHadEtCorrections"<<endl;

}
AliAnalysisHadEtCorrections::~AliAnalysisHadEtCorrections()
{//destructor
  //Clear();
    delete fnotIDTPC;
    delete fnotIDITS;
    delete fnotIDNoID;
    delete fEfficiencyPionTPC;
    delete fEfficiencyKaonTPC;
    delete fEfficiencyProtonTPC;
    delete fEfficiencyHadronTPC;
    delete fEfficiencyPionITS;
    delete fEfficiencyKaonITS;
    delete fEfficiencyProtonITS;
    delete fEfficiencyHadronITS;
    delete fBackgroundTPC;
    delete fBackgroundITS;
//     fnotIDTPC->Clear();
//     fnotIDITS->Clear();
//     fnotIDNoID->Clear();
//     fEfficiencyPionTPC->Clear();
//     fEfficiencyKaonTPC->Clear();
//     fEfficiencyProtonTPC->Clear();
//     fEfficiencyHadronTPC->Clear();
//     fEfficiencyPionITS->Clear();
//     fEfficiencyKaonITS->Clear();
//     fEfficiencyProtonITS->Clear();
//     fEfficiencyHadronITS->Clear();
//     fBackgroundTPC->Clear();
//     fBackgroundITS->Clear();
}
AliAnalysisHadEtCorrections::AliAnalysisHadEtCorrections(const AliAnalysisHadEtCorrections *g): TNamed(),
  fEtaCut(g->fEtaCut)
  ,fAcceptanceCorrectionFull(g->fAcceptanceCorrectionFull)
  ,fAcceptanceCorrectionEMCAL(g->fAcceptanceCorrectionEMCAL)
  ,fAcceptanceCorrectionPHOS(g->fAcceptanceCorrectionPHOS)
  ,fNeutralCorrection(g->fNeutralCorrection)
  ,fNotHadronicCorrection(g->fNotHadronicCorrection)
  ,fpTcutCorrectionTPC(g->fpTcutCorrectionTPC)
  ,fpTcutCorrectionITS(g->fpTcutCorrectionITS)
  ,fNeutralCorrectionLow(g->fNeutralCorrectionLow)
  ,fNotHadronicCorrectionLow(g->fNotHadronicCorrectionLow)
  ,ffpTcutCorrectionTPCLow(g->ffpTcutCorrectionTPCLow)
  ,ffpTcutCorrectionITSLow(g->ffpTcutCorrectionITSLow)
  ,fNeutralCorrectionHigh(g->fNeutralCorrectionHigh)
  ,fNotHadronicCorrectionHigh(g->fNotHadronicCorrectionHigh)
  ,ffpTcutCorrectionTPCHigh(g->ffpTcutCorrectionTPCHigh)
  ,ffpTcutCorrectionITSHigh(g->ffpTcutCorrectionITSHigh)
  ,fnotIDTPC(0)
  ,fnotIDITS(0)
  ,fnotIDNoID(0)
  ,fEfficiencyPionTPC(0)
  ,fEfficiencyKaonTPC(0)
  ,fEfficiencyProtonTPC(0)
  ,fEfficiencyHadronTPC(0)
  ,fEfficiencyPionITS(0)
  ,fEfficiencyKaonITS(0)
  ,fEfficiencyProtonITS(0)
  ,fEfficiencyHadronITS(0)
  ,fBackgroundTPC(0)
  ,fBackgroundITS(0)
{//copy constructor
  //SetName(g->GetName());
  fnotIDTPC = new TH1D(*(g->fnotIDTPC));
  fnotIDITS = new TH1D(*(g->fnotIDITS));
  fnotIDNoID = new TH1D(*(g->fnotIDNoID));
  fEfficiencyPionTPC = new TH1D(*(g->fEfficiencyPionTPC));
  fEfficiencyKaonTPC = new TH1D(*(g->fEfficiencyKaonTPC));
  fEfficiencyProtonTPC = new TH1D(*(g->fEfficiencyProtonTPC));
  fEfficiencyHadronTPC = new TH1D(*(g->fEfficiencyHadronTPC));
  fEfficiencyPionITS = new TH1D(*(g->fEfficiencyPionITS));
  fEfficiencyKaonITS = new TH1D(*(g->fEfficiencyKaonITS));
  fEfficiencyProtonITS = new TH1D(*(g->fEfficiencyProtonITS));
  fEfficiencyHadronITS = new TH1D(*(g->fEfficiencyHadronITS));
  fBackgroundTPC = new TH1D(*(g->fBackgroundTPC));
  fBackgroundITS = new TH1D(*(g->fBackgroundITS));
}


Float_t AliAnalysisHadEtCorrections::GetConstantCorrections(Bool_t totEt, Float_t ptcut, TString type){
  Float_t acceptance = 0.0;
  Float_t neutral = 0.0;
  Float_t ptcorr = 0.0;
  float correction = 0.0;

  //TString *type = new TString(mytype);

  if(type.Contains("Full")) acceptance = fAcceptanceCorrectionFull;
  if(type.Contains("EMCAL")) acceptance = fAcceptanceCorrectionEMCAL;
  if(type.Contains("PHOS")) acceptance = fAcceptanceCorrectionPHOS;

  if(type.Contains("High")){//high bound
    if(totEt) neutral = fNotHadronicCorrectionHigh;
    else{neutral = fNeutralCorrectionHigh;}
    if(ptcut>0.12){ptcorr = ffpTcutCorrectionTPCHigh;}
    else{ptcorr = ffpTcutCorrectionITSHigh;}
    cout<<"Setting correction factor to "<<correction<<endl;
    return correction;
  }
  if(type.Contains("Low")){//high bound
    if(totEt) neutral = fNotHadronicCorrectionLow;
    else{neutral = fNeutralCorrectionLow;}
    if(ptcut>0.12){ptcorr = ffpTcutCorrectionTPCLow;}
    else{ptcorr = ffpTcutCorrectionITSLow;}
    cout<<"Setting correction factor to "<<correction<<endl;
    return correction;
  }

  if(totEt) neutral = fNotHadronicCorrection;
  else{neutral = fNeutralCorrection;}
  if(ptcut>0.12){ptcorr = fpTcutCorrectionTPC;}
  else{ptcorr = fpTcutCorrectionITS;}

  correction = acceptance*neutral*ptcorr;
  cout<<"Setting correction factor for ";
  if(totEt) cout<<"total et";
  else{cout<<"hadronic et";}
  cout<<" with the pt cut off "<<ptcut<<" for "<<type<<" acceptance to "<<correction<<endl;
  //cout<<"Acceptance "<<acceptance<<" neutral "<<neutral<<" ptcorr "<<ptcorr<<endl;
  return correction;

}
// AliAnalysisHadEtCorrections & operator = (const AliAnalysisHadEtCorrections & g) {

//   fEtaCut=g->fEtaCut;
//   fAcceptanceCorrectionFull=g->fAcceptanceCorrectionFull;
//   fAcceptanceCorrectionEMCAL=g->fAcceptanceCorrectionEMCAL;
//   fAcceptanceCorrectionPHOS=g->fAcceptanceCorrectionPHOS;
//   fNeutralCorrection=g->fNeutralCorrection;
//   fNotHadronicCorrection=g->fNotHadronicCorrection;
//   fpTcutCorrectionTPC=g->fpTcutCorrectionTPC;
//   fpTcutCorrectionITS=g->fpTcutCorrectionITS;
//   fNeutralCorrectionLow=g->fNeutralCorrectionLow;
//   fNotHadronicCorrectionLow=g->fNotHadronicCorrectionLow;
//   ffpTcutCorrectionTPCLow=g->ffpTcutCorrectionTPCLow;
//   ffpTcutCorrectionITSLow=g->ffpTcutCorrectionITSLow;
//   fNeutralCorrectionHigh=g->fNeutralCorrectionHigh;
//   fNotHadronicCorrectionHigh=g->fNotHadronicCorrectionHigh;
//   ffpTcutCorrectionTPCHigh=g->ffpTcutCorrectionTPCHigh;
//   ffpTcutCorrectionITSHigh=g->ffpTcutCorrectionITSHigh;

//   fnotIDTPC = g->fnotIDTPC;
//   fnotIDITS = g->fnotIDITS;
//   fnotIDNoID = g->fnotIDNoID;
//   fEfficiencyPionTPC = g->fEfficiencyPionTPC;
//   fEfficiencyKaonTPC = g->fEfficiencyKaonTPC;
//   fEfficiencyProtonTPC = g->fEfficiencyProtonTPC;
//   fEfficiencyHadronTPC = g->fEfficiencyHadronTPC;
//   fEfficiencyPionITS = g->fEfficiencyPionITS;
//   fEfficiencyKaonITS = g->fEfficiencyKaonITS;
//   fEfficiencyProtonITS = g->fEfficiencyProtonITS;
//   fEfficiencyHadronITS = g->fEfficiencyHadronITS;
//   fBackgroundTPC = g->fBackgroundTPC;
//   fBackgroundITS = g->fBackgroundITS;
// }
