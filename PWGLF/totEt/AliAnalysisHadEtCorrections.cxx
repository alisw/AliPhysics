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
#include "TObjArray.h"
#include "AliLog.h"
#include "TH1D.h"

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
							   ,fNotIDConstTPC(0)
							   ,fNotIDConstITS(0)
							   ,fNotIDConstTPCNoID(0)
							   ,fNotIDConstITSNoID(0)
							   ,fNeutralCorrectionLow(0)
							   ,fNotHadronicCorrectionLow(0)
							   ,ffpTcutCorrectionTPCLow(0)
							   ,ffpTcutCorrectionITSLow(0)
							   ,fNeutralCorrectionHigh(0)
							   ,fNotHadronicCorrectionHigh(0)
							   ,ffpTcutCorrectionTPCHigh(0)
							   ,ffpTcutCorrectionITSHigh(0)
							   ,fNotIDConstTPCLow(0)
							   ,fNotIDConstITSLow(0)
							   ,fNotIDConstTPCNoIDLow(0)
							   ,fNotIDConstITSNoIDLow(0)
							   ,fNotIDConstTPCHigh(0)
							   ,fNotIDConstITSHigh(0)
							   ,fNotIDConstTPCNoIDHigh(0)
							   ,fNotIDConstITSNoIDHigh(0)
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
							   ,fEfficiencyTPC(0)
							   ,fEfficiencyITS(0)
							   ,fEfficiencyErrorLow(0)
							   ,fEfficiencyErrorHigh(0)
							   ,fBackgroundErrorLow(0)
							   ,fBackgroundErrorHigh(0)
							   ,fBackgroundTPC(0)
							   ,fBackgroundITS(0)
							   ,fIsEMCal(kTRUE)
							   ,fIsData(kFALSE)
							   ,fDataSet(2009)
							   ,fProduction("ProductionName")
							   ,fProductionDescription("Long production description")
{//default constructor
  Init();
}
void AliAnalysisHadEtCorrections::Init() 
{  //This seems to solve a compiler error
   cout<<"Creating new AliAnalysisHadEtCorrections"<<endl;
   fEfficiencyTPC = new TObjArray();
   fEfficiencyITS = new TObjArray();
}

AliAnalysisHadEtCorrections::~AliAnalysisHadEtCorrections()
{//destructor
  //Clear();
    fnotIDTPC->Clear();
    fnotIDITS->Clear();
    fnotIDNoID->Clear();
    fEfficiencyPionTPC->Clear();
    fEfficiencyKaonTPC->Clear();
    fEfficiencyProtonTPC->Clear();
    fEfficiencyHadronTPC->Clear();
    fEfficiencyPionITS->Clear();
    fEfficiencyKaonITS->Clear();
    fEfficiencyProtonITS->Clear();
    fEfficiencyHadronITS->Clear();
    fBackgroundTPC->Clear();
    fBackgroundITS->Clear();
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
    delete fEfficiencyTPC;
    delete fEfficiencyITS;
    delete fBackgroundTPC;
    delete fBackgroundITS;
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
											      ,fNotIDConstTPC(g->fNotIDConstTPC)
											      ,fNotIDConstITS(g->fNotIDConstITS)
											      ,fNotIDConstTPCNoID(g->fNotIDConstTPCNoID)
											      ,fNotIDConstITSNoID(g->fNotIDConstITSNoID)
											      ,fNeutralCorrectionLow(g->fNeutralCorrectionLow)
											      ,fNotHadronicCorrectionLow(g->fNotHadronicCorrectionLow)
											      ,ffpTcutCorrectionTPCLow(g->ffpTcutCorrectionTPCLow)
											      ,ffpTcutCorrectionITSLow(g->ffpTcutCorrectionITSLow)
											      ,fNeutralCorrectionHigh(g->fNeutralCorrectionHigh)
											      ,fNotHadronicCorrectionHigh(g->fNotHadronicCorrectionHigh)
											      ,ffpTcutCorrectionTPCHigh(g->ffpTcutCorrectionTPCHigh)
											      ,ffpTcutCorrectionITSHigh(g->ffpTcutCorrectionITSHigh)
											      ,fNotIDConstTPCLow(g->fNotIDConstTPCLow)
											      ,fNotIDConstITSLow(g->fNotIDConstITSLow)
											      ,fNotIDConstTPCNoIDLow(g->fNotIDConstTPCNoIDLow)
											      ,fNotIDConstITSNoIDLow(g->fNotIDConstITSNoIDLow)
											      ,fNotIDConstTPCHigh(g->fNotIDConstTPCHigh)
											      ,fNotIDConstITSHigh(g->fNotIDConstITSHigh)
											      ,fNotIDConstTPCNoIDHigh(g->fNotIDConstTPCNoIDHigh)
											      ,fNotIDConstITSNoIDHigh(g->fNotIDConstITSNoIDHigh)
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
											      ,fEfficiencyTPC(0)
											      ,fEfficiencyITS(0)
											      ,fEfficiencyErrorLow(g->fEfficiencyErrorLow)
											      ,fEfficiencyErrorHigh(g->fEfficiencyErrorHigh)
											      ,fBackgroundErrorLow(g->fBackgroundErrorLow)
											      ,fBackgroundErrorHigh(g->fBackgroundErrorHigh)
											      ,fBackgroundTPC(0)
											      ,fBackgroundITS(0)
											      ,fIsEMCal(g->fIsEMCal)
											      ,fIsData(g->fIsData)
											      ,fDataSet(g->fDataSet)
											      ,fProduction(g->fProduction)
											      ,fProductionDescription(g->fProductionDescription)
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
  fEfficiencyTPC = new TObjArray(*(g->fEfficiencyTPC));
  fEfficiencyITS = new TObjArray(*(g->fEfficiencyITS));
  fBackgroundTPC = new TH1D(*(g->fBackgroundTPC));
  fBackgroundITS = new TH1D(*(g->fBackgroundITS));
}


Float_t AliAnalysisHadEtCorrections::GetConstantCorrections(Bool_t totEt, Float_t ptcut, TString type) const {//Get the correction values that are not pt dependent
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
  //cout<<" Acceptance "<<acceptance<<" neutral "<<neutral<<" ptcorr "<<ptcorr<<endl;
  return correction;

}
Float_t AliAnalysisHadEtCorrections::GetSystematicErrorBound(Float_t et,Bool_t isLowBound, Bool_t isHadronic, Bool_t isTPC) const{
  //we calculate factors for each value and then multiply them to get the overall bounds
  //neutral corrections, pt cut, pid, efficiency, background
  float neutral = 1.0;
  float ptcut = 1.0;
  float pid = 1.0;
  float efficiency = 1.0;
  float background = 1.0;
  if(isLowBound){//is lower bound
    if(isHadronic) neutral= fNeutralCorrectionLow/fNeutralCorrection;
    else{neutral = fNotHadronicCorrectionLow/fNotHadronicCorrection;}
    if(isTPC) ptcut = ffpTcutCorrectionTPCLow/fpTcutCorrectionTPC;
    else{ptcut = ffpTcutCorrectionITSLow/fpTcutCorrectionITS;}
    pid = fNotIDConstTPCLow/fNotIDConstTPC;
    efficiency = fEfficiencyErrorLow;
    background = fBackgroundErrorLow;
  }
  else{//is higher bound
    if(isHadronic) neutral= fNeutralCorrectionHigh/fNeutralCorrection;
    else{neutral= fNotHadronicCorrectionHigh/fNotHadronicCorrection;}
    if(isTPC) ptcut = ffpTcutCorrectionTPCHigh/fpTcutCorrectionTPC;
    else{ptcut = ffpTcutCorrectionITSHigh/fpTcutCorrectionITS;}
    pid = fNotIDConstTPCHigh/fNotIDConstTPC;
    efficiency = fEfficiencyErrorHigh;
    background = fBackgroundErrorHigh;
  }
  //cout<<"neutral "<<neutral<<" ptcut "<<ptcut<<" pid "<<pid<<" efficiency "<<efficiency<<" background "<<background<<" overall "<<neutral*ptcut*pid*efficiency*background<<endl;
  return neutral*ptcut*pid*efficiency*background*et;
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
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyPionTPC(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyPionTPC;}
  else{return (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyPionTPC%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyKaonTPC(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyKaonTPC;}
  else{return (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyKaonTPC%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyProtonTPC(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyProtonTPC;}
  else{return (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyProtonTPC%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyHadronTPC(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyHadronTPC;}
  else{return (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyHadronTPC%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyPionITS(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyPionITS;}
  else{return (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyPionITS%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyKaonITS(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyKaonITS;}
  else{return (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyKaonITS%i",cb));}
}
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyProtonITS(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyProtonITS;}
  else{return (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyProtonITS%i",cb));}
}//Proton
TH1D *AliAnalysisHadEtCorrections::GetEfficiencyHadronITS(const int cb){//Get centrality dependent efficiency
  if(cb==-1){return fEfficiencyHadronITS;}
  else{return (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyHadronITS%i",cb));}
}
Float_t AliAnalysisHadEtCorrections::GetTPCEfficiencyCorrectionPion(const float pT, const int cb){//Get the efficiency for reconstructing a pion in the TPC
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyPionTPC){cerr<<"No histogram fEfficiencyPionTPC!"<<endl; return -1.0;}
    eff = fEfficiencyPionTPC->GetBinContent(fEfficiencyPionTPC->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyPionTPC%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyPionTPC%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetTPCEfficiencyCorrectionKaon(const float pT, const int cb){//Get the efficiency for reconstructing a kaon in the TPC
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyKaonTPC){cerr<<"No histogram fEfficiencyKaonTPC!"<<endl; return -1.0;}
    eff = fEfficiencyKaonTPC->GetBinContent(fEfficiencyKaonTPC->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyKaonTPC%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyKaonTPC%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetTPCEfficiencyCorrectionProton(const float pT, const int cb){//Get the efficiency for reconstructing a proton in the TPC
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyProtonTPC){cerr<<"No histogram fEfficiencyProtonTPC!"<<endl; return -1.0;}
    eff = fEfficiencyProtonTPC->GetBinContent(fEfficiencyProtonTPC->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyProtonTPC%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyProtonTPC%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetTPCEfficiencyCorrectionHadron(const float pT, const int cb){//Get the efficiency for reconstructing a hadron in the TPC
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyHadronTPC){cerr<<"No histogram fEfficiencyHadronTPC!"<<endl; return -1.0;}
    eff = fEfficiencyHadronTPC->GetBinContent(fEfficiencyHadronTPC->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyHadronTPC%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyHadronTPC%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetITSEfficiencyCorrectionPion(const float pT, const int cb){//Get the efficiency for reconstructing a pion in the ITS
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyPionITS){cerr<<"No histogram fEfficiencyPionITS!"<<endl; return -1.0;}
    eff = fEfficiencyPionITS->GetBinContent(fEfficiencyPionITS->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyPionITS%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyPionITS%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetITSEfficiencyCorrectionKaon(const float pT, const int cb){//Get the efficiency for reconstructing a kaon in the ITS
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyKaonITS){cerr<<"No histogram fEfficiencyKaonITS!"<<endl; return -1.0;}
    eff = fEfficiencyKaonITS->GetBinContent(fEfficiencyKaonITS->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyKaonITS%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyKaonITS%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetITSEfficiencyCorrectionProton(const float pT, const int cb){//Get the efficiency for reconstructing a proton in the ITS
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyProtonITS){cerr<<"No histogram fEfficiencyProtonITS!"<<endl; return -1.0;}
    eff = fEfficiencyProtonITS->GetBinContent(fEfficiencyProtonITS->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyProtonITS%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyProtonITS%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
Float_t AliAnalysisHadEtCorrections::GetITSEfficiencyCorrectionHadron(const float pT, const int cb){//Get the efficiency for reconstructing a hadron in the ITS
  float eff = -1.0;
  if(cb ==-1){//pp
    if(!fEfficiencyHadronITS){cerr<<"No histogram fEfficiencyHadronITS!"<<endl; return -1.0;}
    eff = fEfficiencyHadronITS->GetBinContent(fEfficiencyHadronITS->FindBin(pT));
    if(eff<=0.0){AliInfo("Efficiency is zero!");  return 0.0;}
  }
  else{
    TH1D *fEfficiency = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyHadronITS%i",cb));
    if(!fEfficiency){cerr<<"No histogram "<<Form("fEfficiencyHadronITS%i",cb)<<endl; return -1.0;}
    eff = fEfficiency->GetBinContent(fEfficiency->FindBin(pT));
    if(eff<=0.0){AliInfo(Form("Pion efficiency is zero for centrality bin %i!",cb));  return 0.0;}
  }
  return 1.0/eff;
}
void AliAnalysisHadEtCorrections::SetEfficiencyPionTPC(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyPionTPC%i",cb));
    if(old){
      fEfficiencyTPC->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyPionTPC%i",cb));
    fEfficiencyTPC->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}
void AliAnalysisHadEtCorrections::SetEfficiencyKaonTPC(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyKaonTPC%i",cb));
    if(old){
      fEfficiencyTPC->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyKaonTPC%i",cb));
    fEfficiencyTPC->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}//Kaon
void AliAnalysisHadEtCorrections::SetEfficiencyProtonTPC(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyProtonTPC%i",cb));
    if(old){
      fEfficiencyTPC->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyProtonTPC%i",cb));
    fEfficiencyTPC->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}//Proton
void AliAnalysisHadEtCorrections::SetEfficiencyHadronTPC(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyHadronTPC%i",cb));
    if(old){
      fEfficiencyTPC->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyHadronTPC%i",cb));
    fEfficiencyTPC->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}
void AliAnalysisHadEtCorrections::SetEfficiencyPionITS(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyPionITS%i",cb));
    if(old){
      fEfficiencyITS->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyPionITS%i",cb));
    fEfficiencyITS->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}
void AliAnalysisHadEtCorrections::SetEfficiencyKaonITS(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyKaonITS%i",cb));
    if(old){
      fEfficiencyITS->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyKaonITS%i",cb));
    fEfficiencyITS->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}//Kaon
void AliAnalysisHadEtCorrections::SetEfficiencyProtonITS(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyProtonITS%i",cb));
    if(old){
      fEfficiencyITS->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyProtonITS%i",cb));
    fEfficiencyITS->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}//Proton
void AliAnalysisHadEtCorrections::SetEfficiencyHadronITS(TH1D *histo, const int cb){//Set centrality dependent efficiency for centrality bin cb
  if(histo){
    //first check to see if the histogram exists already
    TH1D *old = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyHadronITS%i",cb));
    if(old){
      fEfficiencyITS->Remove(old);
      delete old;
    }
    //then if the new histogram exists, add it to the array
    histo->SetName(Form("fEfficiencyHadronITS%i",cb));
    fEfficiencyITS->Add(histo);
  }
  else{cerr<<"Histogram does not exist!"<<endl;}
}


Float_t AliAnalysisHadEtCorrections::GetNotIDCorrectionTPC(const float pT){//get correction for unidentified particles in the TPC
  Float_t val = fnotIDTPC->GetBinContent(fnotIDTPC->FindBin(pT));
  if(val>0.0) return 1.0/(val);
  else{return 0.0;}
}
Float_t AliAnalysisHadEtCorrections::GetNotIDCorrectionITS(const float pT){//Get correction for unidentified particles in the ITS
  Float_t val = fnotIDITS->GetBinContent(fnotIDITS->FindBin(pT));
  if(val>0.0) return 1.0/(val);
  else{return 0.0;}
}
Float_t AliAnalysisHadEtCorrections::GetNotIDCorrectionNoPID(const float pT){//Get correction for particles in the case that there is no particle identification
  Float_t val = fnotIDNoID->GetBinContent(fnotIDNoID->FindBin(pT));
  if(val>0.0) return 1.0/(val);
  else{return 0.0;}
}
Float_t AliAnalysisHadEtCorrections::GetBackgroundCorrectionTPC(const float pT){//Get background correction for TPC tracks
  return (1.0-fBackgroundTPC->GetBinContent(fBackgroundTPC->FindBin(pT)));
}
Float_t AliAnalysisHadEtCorrections::GetBackgroundCorrectionITS(const float pT){//Get background correction for ITS tracks
  return (1.0-fBackgroundITS->GetBinContent(fBackgroundITS->FindBin(pT)));
}
void AliAnalysisHadEtCorrections::Report(){//Gives a report on the status of all corrections
  //This is primarily for cross checking that the results we get from the macro that fills this class, GetCorrections.C, are sane
  cout<<"======================================================================="<<endl;
  cout<<"                   Report from "<<GetName()<<endl;
  cout<<"======================================================================="<<endl;
  cout<<fProductionDescription<<" created from "<<fProduction<<endl;
  cout<<"This for determination of EThad from ";
  if(fIsData) cout<<"data of ";
  else{cout<<"simulation of ";}
  switch(fDataSet){
  case 2009:
    cout<<"p+p collisions at 900 GeV"<<endl;
    break;
  case 2010:
    cout<<"p+p collisions at 7 TeV"<<endl;
    break;
  case 20111:
    cout<<"p+p collisions at 2.76 TeV"<<endl;
    break;
  case 20100:
    cout<<"Pb+Pb collisions at 2.76 TeV"<<endl;
    break;
  default:
    cout<<"an undetermined collision system and energy"<<endl;
  }
  cout<<"This is initialized for the ";
  if(fIsEMCal) cout<<"EMCal";
  else{cout<<"PHOS";}
  cout<<" acceptance"<<endl<<endl;

  cout<<"The acceptance correction for the full  "<<fAcceptanceCorrectionFull<<endl;
  cout<<"                                  EMCal "<<fAcceptanceCorrectionEMCAL<<endl;
  cout<<"                                  PHOS  "<<fAcceptanceCorrectionPHOS<<endl<<endl;

  cout<<Form("The neutral energy correction is %2.4f [%2.4f,%2.4f]",fNeutralCorrection,fNeutralCorrectionLow,fNeutralCorrectionHigh)<<endl;
  cout<<Form("    total                        %2.4f [%2.4f,%2.4f]",fNotHadronicCorrection,fNotHadronicCorrectionLow,fNotHadronicCorrectionHigh)<<endl<<endl;

  cout<<Form("The pT cut correction for 100 MeV is %2.4f [%2.4f,%2.4f]",fpTcutCorrectionITS,ffpTcutCorrectionITSLow,ffpTcutCorrectionITSHigh)<<endl;
  cout<<Form("                          150 MeV    %2.4f [%2.4f,%2.4f]",fpTcutCorrectionTPC,ffpTcutCorrectionTPCLow,ffpTcutCorrectionTPCHigh)<<endl<<endl;

  cout<<Form("The correction for unidentified ITS tracks is %2.4f [%2.4f,%2.4f]",fNotIDConstITS,fNotIDConstITSLow,fNotIDConstITSHigh)<<endl;
  cout<<Form("                                TPC tracks    %2.4f [%2.4f,%2.4f]",fNotIDConstTPC,fNotIDConstTPCLow,fNotIDConstTPCHigh)<<endl<<endl;

  cout<<"Background correction histogram for ITS tracks is";
  if(!fBackgroundITS)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fBackgroundITS->GetEntries()<<" entries"<<endl;}
  cout<<"                                    TPC          ";
  if(!fBackgroundTPC)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fBackgroundTPC->GetEntries()<<" entries"<<endl;}
  cout<<endl;


  cout<<"Efficiency histogram for ITS tracks for hadrons is";
  if(!fEfficiencyHadronITS)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyHadronITS->GetEntries()<<" entries"<<endl;}
  cout<<"                         TPC            hadrons   ";
  if(!fEfficiencyHadronTPC)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyHadronTPC->GetEntries()<<" entries"<<endl;}
  cout<<"                         ITS            pions     ";
  if(!fEfficiencyPionITS)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyPionITS->GetEntries()<<" entries"<<endl;}
  cout<<"                         TPC            pions     ";
  if(!fEfficiencyPionTPC)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyPionTPC->GetEntries()<<" entries"<<endl;}
  cout<<"                         ITS            kaons     ";
  if(!fEfficiencyKaonITS)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyKaonITS->GetEntries()<<" entries"<<endl;}
  cout<<"                         TPC            kaons     ";
  if(!fEfficiencyKaonTPC)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyKaonTPC->GetEntries()<<" entries"<<endl;}
  cout<<"                         ITS            protons   ";
  if(!fEfficiencyProtonITS)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyProtonITS->GetEntries()<<" entries"<<endl;}
  cout<<"                         TPC            protons   ";
  if(!fEfficiencyProtonTPC)cout<<" not set"<<endl;
  else{cout<<" set and has "<<fEfficiencyProtonTPC->GetEntries()<<" entries"<<endl;}
  cout<<endl;

  if(fDataSet==20100){//if Pb+Pb
    cout<<"Efficiency histogram for TPC tracks for hadrons is set for centrality bins ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyHadronTPC%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        pions                              ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyPionTPC%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        kaons                              ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyKaonTPC%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        protons                            ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyTPC->FindObject(Form("fEfficiencyProtonTPC%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                         ITS            hadrons                            ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyHadronITS%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        pions                              ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyPionITS%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        kaons                              ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyKaonITS%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;
    cout<<"                                        protons                            ";
    for(int i = 0;i<=21;i++){
      TH1D *histo = (TH1D*) fEfficiencyITS->FindObject(Form("fEfficiencyProtonITS%i",i));
      if(histo) cout<<i<<" ";
    }
    cout<<endl;


    int nEntries = fEfficiencyTPC->GetEntries();
    int nbadhistograms = 0;
    for(int i=0;i<nEntries;i++){
      TH1D *histo = (TH1D*) fEfficiencyTPC->At(i);
      if(!histo){
	cout<<"Warning:  Histogram in fEfficiencyTPC at "<<i<<" is NULL!"<<endl;
	nbadhistograms++;
      }
      else{
	if(histo->GetEntries()<=1e-2){
	  cout<<"Warning: Histogram "<<histo->GetName()<<" in fEfficiencyTPC is empty!"<<endl;
	  nbadhistograms++;
	}
      }
    }
    cout<<nbadhistograms<<" bad histograms in fEfficiencyTPC"<<endl;

    nEntries = fEfficiencyITS->GetEntries();
    nbadhistograms = 0;
    for(int i=0;i<nEntries;i++){
      TH1D *histo = (TH1D*) fEfficiencyITS->At(i);
      if(!histo){
	cout<<"Warning:  Histogram in fEfficiencyITS at "<<i<<" is NULL!"<<endl;
	nbadhistograms++;
      }
      else{
	if(histo->GetEntries()<=1e-2){
	  cout<<"Warning: Histogram "<<histo->GetName()<<" in fEfficiencyITS is empty!"<<endl;
	  nbadhistograms++;
	}
      }
    }
    cout<<nbadhistograms<<" bad histograms in fEfficiencyITS"<<endl;


  }


  cout<<endl;
  cout<<"======================================================================="<<endl;
  cout<<"======================================================================="<<endl;
}

