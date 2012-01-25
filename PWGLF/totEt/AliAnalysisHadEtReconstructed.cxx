//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for ESD analysis
//  - reconstruction output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliVParticle.h"
#include <iostream>
#include "AliAnalysisHadEtCorrections.h"
#include "TFile.h"
#include "TString.h"
#include "AliAnalysisEtCommon.h"
#include "AliAnalysisHadEt.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliPWG0Helper.h"

using namespace std;

ClassImp(AliAnalysisHadEtReconstructed);


AliAnalysisHadEtReconstructed::AliAnalysisHadEtReconstructed() :
  AliAnalysisHadEt()
  ,fCorrections(0)
  ,fConfigFile("ConfigHadEtAnalysis.C")
  ,fCorrTotEtFullAcceptanceTPC(0)
  ,fCorrTotEtFullAcceptanceITS(0)
  ,fCorrHadEtFullAcceptanceTPC(0)
  ,fCorrHadEtFullAcceptanceITS(0)
  ,fCorrTotEtEMCALAcceptanceTPC(0)
  ,fCorrTotEtEMCALAcceptanceITS(0)
  ,fCorrHadEtEMCALAcceptanceTPC(0)
  ,fCorrHadEtEMCALAcceptanceITS(0)
  ,fCorrTotEtPHOSAcceptanceTPC(0)
  ,fCorrTotEtPHOSAcceptanceITS(0)
  ,fCorrHadEtPHOSAcceptanceTPC(0)
  ,fCorrHadEtPHOSAcceptanceITS(0)
  ,fCorrectedHadEtFullAcceptanceTPCNoPID(0)
  ,fCorrectedHadEtFullAcceptanceITSNoPID(0)
  ,fCorrectedHadEtEMCALAcceptanceTPCNoPID(0)
  ,fCorrectedHadEtEMCALAcceptanceITSNoPID(0)
  ,fCorrectedHadEtPHOSAcceptanceTPCNoPID(0)
  ,fCorrectedHadEtPHOSAcceptanceITSNoPID(0)
  ,fCorrectedHadEtFullAcceptanceTPC(0)
  ,fCorrectedHadEtFullAcceptanceITS(0)
  ,fCorrectedHadEtFullAcceptanceTPCAssumingPion(0)
  ,fCorrectedHadEtFullAcceptanceITSAssumingPion(0)
  ,fCorrectedHadEtFullAcceptanceTPCAssumingProton(0)
  ,fCorrectedHadEtFullAcceptanceITSAssumingProton(0)
  ,fCorrectedHadEtFullAcceptanceTPCAssumingKaon(0)
  ,fCorrectedHadEtFullAcceptanceITSAssumingKaon(0)
  ,fCorrectedHadEtEMCALAcceptanceTPC(0)
  ,fCorrectedHadEtEMCALAcceptanceITS(0)
  ,fCorrectedHadEtPHOSAcceptanceTPC(0)
  ,fCorrectedHadEtPHOSAcceptanceITS(0)
  ,fRawEtFullAcceptanceTPC(0)
  ,fRawEtFullAcceptanceITS(0)
  ,fRawEtEMCALAcceptanceTPC(0)
  ,fRawEtEMCALAcceptanceITS(0)
  ,fRawEtPHOSAcceptanceTPC(0)
  ,fRawEtPHOSAcceptanceITS(0)
  ,fRawEtFullAcceptanceTPCNoPID(0)
  ,fRawEtFullAcceptanceITSNoPID(0)
  ,fRawEtEMCALAcceptanceTPCNoPID(0)
  ,fRawEtEMCALAcceptanceITSNoPID(0)
  ,fRawEtPHOSAcceptanceTPCNoPID(0)
  ,fRawEtPHOSAcceptanceITSNoPID(0)
{
}

AliAnalysisHadEtReconstructed::~AliAnalysisHadEtReconstructed() 
{
  delete fCorrections;
}

Int_t AliAnalysisHadEtReconstructed::AnalyseEvent(AliVEvent* ev, Int_t eventtype)
{ // analyse ESD event
  ResetEventValues();
  if(!ev){
    AliFatal("ERROR: Event does not exist");   
    return 0;
  }

  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev);
  if(!realEvent){  
    AliFatal("ERROR: ESD Event does not exist");
    return 0;
  }
  fCentBin= -1;
  fGoodEvent = kTRUE;//for p+p collisions if we made it this far we have a good event
  if(fDataSet==20100){//If this is Pb+Pb
    AliCentrality *centrality = realEvent->GetCentrality();
    if(fNCentBins<21) fCentBin= centrality->GetCentralityClass10(fCentralityMethod);
    else{ fCentBin= centrality->GetCentralityClass5(fCentralityMethod);}
    if(fCentBin ==-1) fGoodEvent = kFALSE;//but for Pb+Pb events we don't want to count events where we did not find a centrality
  }
  //for PID
  AliESDpid *pID = new AliESDpid();
  pID->MakePID(realEvent);
  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  for(Int_t cutset=0;cutset<2;cutset++){
    bool isTPC = false;
    TString *cutName = NULL;
    TObjArray* list = NULL;
    switch(cutset){
    case 0:
      cutName = strTPCITS;
      list = fEsdtrackCutsITSTPC->GetAcceptedTracks(realEvent);
      isTPC = true;
      break;
    case 1:
      cutName = strITS;
      list = fEsdtrackCutsITS->GetAcceptedTracks(realEvent);
      break;
    case 2:
      cutName = strTPC;
      list = fEsdtrackCutsTPC->GetAcceptedTracks(realEvent);
      break;
    default:
      cerr<<"Error:  cannot fill histograms!"<<endl;
      return -1;
    }
    Int_t nGoodTracks = list->GetEntries();
    for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++)
      {


	AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
	if (!track)
	  {
	    Printf("ERROR: Could not get track %d", iTrack);
	    continue;
	  }
	else{
	  if(TMath::Abs(track->Eta())>fCorrections->GetEtaCut()) continue;
	  Float_t nSigmaPion,nSigmaProton,nSigmaKaon,nSigmaElectron;
	  pID->MakeTPCPID(track);
	  pID->MakeITSPID(track);
	  if(cutset!=1){
	    nSigmaPion = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kPion));
	    nSigmaProton = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kProton));
	    nSigmaKaon = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kKaon));
	    nSigmaElectron = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kElectron));
	  }
	  else{
	    nSigmaPion = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kPion));
	    nSigmaProton = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kProton));
	    nSigmaKaon = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kKaon));
	    nSigmaElectron = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kElectron));
	  }
	  // 	    bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  // 	    bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  // 	    bool isKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
	  // 	    bool isProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);
	  bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  bool isKaon = (nSigmaPion>3.0 && nSigmaProton>3.0 && nSigmaKaon<3.0 && track->Pt()<0.45);
	  bool isProton = (nSigmaPion>3.0 && nSigmaProton<3.0 && nSigmaKaon>3.0 && track->Pt()<0.9);

	  bool unidentified = (!isProton && !isKaon && !isElectron && !isPion);
	  if(cutset==1){//ITS dE/dx identification requires tighter cuts on the tracks and we don't gain much from that so we won't do it
	    unidentified = true;
	    isPion=false;
	    isElectron=false;
	    isKaon=false;
	    isProton=false;
	  }
	  Float_t dEdx = track->GetTPCsignal();
	  if(cutset==1) dEdx = track->GetITSsignal();
	  FillHisto2D(Form("dEdxDataAll%s",cutName->Data()),track->P(),dEdx,1.0);

	  bool inPHOS = IsInPHOS(track);
	  bool inEMCAL = IsInEMCAL(track);

	  Float_t corrBkgd=0.0;
	  Float_t corrNotID=0.0;
	  Float_t corrNoID = fCorrections->GetNotIDCorrectionNoPID(track->Pt());
	  Float_t corrEff = 0.0;
	  Float_t corrEffNoID = 0.0;
	  if(cutset!=1){//TPC
	    corrBkgd = fCorrections->GetBackgroundCorrectionTPC(track->Pt());
	    corrEffNoID = fCorrections->GetTPCEfficiencyCorrectionHadron(track->Pt(),fCentBin);
	    corrNotID = fCorrections->GetNotIDConstCorrectionTPC();
	    corrNoID = fCorrections->GetNotIDConstCorrectionTPCNoID();
	  }
	  if(cutset==1){//ITS
	    corrBkgd = fCorrections->GetBackgroundCorrectionITS(track->Pt());
	    corrEffNoID = fCorrections->GetITSEfficiencyCorrectionHadron(track->Pt(),fCentBin);
	    corrNotID = fCorrections->GetNotIDConstCorrectionITS();
	    corrNoID = fCorrections->GetNotIDConstCorrectionITSNoID();
	  }
	  if(fDataSet==20100){
	    FillHisto2D("fbkgdVsCentralityBin",fCentBin,corrBkgd,1.0);
	    FillHisto2D("fnotIDVsCentralityBin",fCentBin,corrNotID,1.0);
	    FillHisto2D("fpTcutVsCentralityBin",fCentBin,fCorrections->GetpTCutCorrectionTPC(),1.0);
	    if(fCorrHadEtFullAcceptanceTPC>0.0) FillHisto2D("fneutralVsCentralityBin",fCentBin,1.0/fCorrHadEtFullAcceptanceTPC,1.0);
	    if(fCorrections->GetNeutralCorrection()>0.0) FillHisto2D("ConstantCorrectionsVsCentralityBin",fCentBin,1.0/fCorrections->GetNeutralCorrection(),1.0);
	  }
	  Float_t et = 0.0;
	  Float_t etNoID = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	  Float_t etpartialcorrected = 0.0;
	  Float_t etpartialcorrectedPion = 0.0;
	  Float_t etpartialcorrectedKaon = 0.0;
	  Float_t etpartialcorrectedProton = 0.0;
	  Float_t etpartialcorrectedNoID = corrNoID*corrBkgd*corrEffNoID*etNoID;
	  FillHisto2D(Form("EtDataRaw%sNoID",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrectedNoID);

	  if(isPion){
	    FillHisto2D(Form("dEdxDataPion%s",cutName->Data()),track->P(),dEdx,1.0);
	    et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
   	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionPion(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&fDataSet==20100)FillHisto2D("feffPionVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	    if(track->Charge()>0.0){
	      FillHisto2D(Form("EtDataRaw%sPiPlus",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sPiPlus",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	    else{
	      FillHisto2D(Form("EtDataRaw%sPiMinus",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sPiMinus",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	  }
	  if(isKaon){
	    FillHisto2D(Form("dEdxDataKaon%s",cutName->Data()),track->P(),dEdx,1.0);
	    et = Et(track->P(),track->Theta(),fgKPlusCode,track->Charge());
	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionKaon(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&fDataSet==20100)FillHisto2D("feffKaonVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	      
	    if(track->Charge()>0.0){
	      FillHisto2D(Form("EtDataRaw%sKPlus",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sKPlus",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	    else{
	      FillHisto2D(Form("EtDataRaw%sKMinus",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sKMinus",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	  }
	  if(isProton){
	    FillHisto2D(Form("dEdxDataProton%s",cutName->Data()),track->P(),dEdx,1.0);
	    et = Et(track->P(),track->Theta(),fgProtonCode,track->Charge());
	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionProton(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&fDataSet==20100)FillHisto2D("feffProtonVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	      
	    if(track->Charge()>0.0){
	      FillHisto2D(Form("EtDataRaw%sProton",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sProton",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	    else{
	      FillHisto2D(Form("EtDataRaw%sAntiProton",cutName->Data()),track->Pt(),track->Eta(),et);
	      FillHisto2D(Form("EtDataCorrected%sAntiProton",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	  }
	  if(isElectron){
	    FillHisto2D(Form("dEdxDataElectron%s",cutName->Data()),track->P(),dEdx,1.0);
	  }
	  if(unidentified){
	    if(isPion) cerr<<"I should not be here!!  AliAnalysisHadEtReconstructed 273"<<endl; 
	    FillHisto2D(Form("dEdxDataUnidentified%s",cutName->Data()),track->P(),dEdx,1.0);
	    et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	    Float_t etProton = Et(track->P(),track->Theta(),fgProtonCode,track->Charge());
	    Float_t etKaon = Et(track->P(),track->Theta(),fgKPlusCode,track->Charge());
	    if(corrEff>0.0&&fDataSet==20100)FillHisto2D("feffHadronVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	    etpartialcorrected = et*corrBkgd*corrEffNoID*corrNotID;
	    etpartialcorrectedPion = et*corrBkgd*corrEffNoID;
	    etpartialcorrectedProton = etProton*corrBkgd*corrEffNoID;
	    etpartialcorrectedKaon = etKaon*corrBkgd*corrEffNoID;
	    FillHisto2D(Form("EtDataCorrected%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	  }
	  else{
	    etpartialcorrectedPion = etpartialcorrected;
	    etpartialcorrectedKaon = etpartialcorrected;
	    etpartialcorrectedProton = etpartialcorrected;
	  }
	  if(!isTPC){
	    etpartialcorrected = etpartialcorrectedNoID;//Not using PID for ITS
	  }
	  AddEt(et,etNoID,etpartialcorrected,etpartialcorrectedPion,etpartialcorrectedProton,etpartialcorrectedKaon,etpartialcorrectedNoID,track->Pt(),isTPC,inPHOS,inEMCAL);
	}
      }
    delete list;
  }
  Int_t nondiff = (Int_t) AliPWG0Helper::kND;
  Int_t doublediff = (Int_t) AliPWG0Helper::kDD;
  Int_t singlediff = (Int_t) AliPWG0Helper::kSD;
  if(eventtype == nondiff && fGoodEvent){
    FillHisto1D("RecoHadEtFullAcceptanceTPCND",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCND",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSND",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSND",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDND",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDND",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDND",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDND",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
  }
  if(eventtype == doublediff){
    FillHisto1D("RecoHadEtFullAcceptanceTPCDD",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCDD",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSDD",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSDD",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDDD",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDDD",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDDD",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDDD",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
  }
  if(eventtype == singlediff){
    FillHisto1D("RecoHadEtFullAcceptanceTPCSD",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCSD",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSSD",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSSD",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDSD",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDSD",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDSD",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDSD",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
  }
  if(GetCorrectedHadEtFullAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceTPC",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
  if(GetCorrectedTotEtFullAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtFullAcceptanceTPC",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
  if(GetCorrectedHadEtFullAcceptanceTPCAssumingPion()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingPion",GetCorrectedHadEtFullAcceptanceTPCAssumingPion(),1.0);
  if(GetCorrectedHadEtFullAcceptanceTPCAssumingProton()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingProton",GetCorrectedHadEtFullAcceptanceTPCAssumingProton(),1.0);
  if(GetCorrectedHadEtFullAcceptanceTPCAssumingKaon()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingKaon",GetCorrectedHadEtFullAcceptanceTPCAssumingKaon(),1.0);
  if(GetCorrectedHadEtEMCALAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtEMCALAcceptanceTPC",GetCorrectedHadEtEMCALAcceptanceTPC(),1.0);
  if(GetCorrectedTotEtEMCALAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtEMCALAcceptanceTPC",GetCorrectedTotEtEMCALAcceptanceTPC(),1.0);
  if(GetCorrectedHadEtPHOSAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtPHOSAcceptanceTPC",GetCorrectedHadEtPHOSAcceptanceTPC(),1.0);
  if(GetCorrectedTotEtPHOSAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtPHOSAcceptanceTPC",GetCorrectedTotEtPHOSAcceptanceTPC(),1.0);
  if(GetCorrectedHadEtFullAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceTPCNoPID",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedTotEtFullAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtFullAcceptanceTPCNoPID",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedHadEtEMCALAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtEMCALAcceptanceTPCNoPID",GetCorrectedHadEtEMCALAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedTotEtEMCALAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtEMCALAcceptanceTPCNoPID",GetCorrectedTotEtEMCALAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedHadEtPHOSAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtPHOSAcceptanceTPCNoPID",GetCorrectedHadEtPHOSAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedTotEtPHOSAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtPHOSAcceptanceTPCNoPID",GetCorrectedTotEtPHOSAcceptanceTPCNoPID(),1.0);
  if(GetCorrectedHadEtFullAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceITS",GetCorrectedHadEtFullAcceptanceITS(),1.0);
  if(GetCorrectedHadEtFullAcceptanceITSAssumingPion()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceITSAssumingPion",GetCorrectedHadEtFullAcceptanceITSAssumingPion(),1.0);
  if(GetCorrectedHadEtFullAcceptanceITSAssumingProton()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceITSAssumingProton",GetCorrectedHadEtFullAcceptanceITSAssumingProton(),1.0);
  if(GetCorrectedHadEtFullAcceptanceITSAssumingKaon()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceITSAssumingKaon",GetCorrectedHadEtFullAcceptanceITSAssumingKaon(),1.0);
  if(GetCorrectedTotEtFullAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtFullAcceptanceITS",GetCorrectedTotEtFullAcceptanceITS(),1.0);
  if(GetCorrectedHadEtEMCALAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtEMCALAcceptanceITS",GetCorrectedHadEtEMCALAcceptanceITS(),1.0);
  if(GetCorrectedTotEtEMCALAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtEMCALAcceptanceITS",GetCorrectedTotEtEMCALAcceptanceITS(),1.0);
  if(GetCorrectedHadEtPHOSAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtPHOSAcceptanceITS",GetCorrectedHadEtPHOSAcceptanceITS(),1.0);
  if(GetCorrectedTotEtPHOSAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtPHOSAcceptanceITS",GetCorrectedTotEtPHOSAcceptanceITS(),1.0);
  if(GetCorrectedHadEtFullAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtFullAcceptanceITSNoPID",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
  if(GetCorrectedTotEtFullAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtFullAcceptanceITSNoPID",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
  if(GetCorrectedHadEtEMCALAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtEMCALAcceptanceITSNoPID",GetCorrectedHadEtEMCALAcceptanceITSNoPID(),1.0);
  if(GetCorrectedTotEtEMCALAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtEMCALAcceptanceITSNoPID",GetCorrectedTotEtEMCALAcceptanceITSNoPID(),1.0);
  if(GetCorrectedHadEtPHOSAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoHadEtPHOSAcceptanceITSNoPID",GetCorrectedHadEtPHOSAcceptanceITSNoPID(),1.0);
  if(GetCorrectedTotEtPHOSAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoTotEtPHOSAcceptanceITSNoPID",GetCorrectedTotEtPHOSAcceptanceITSNoPID(),1.0);

  if(GetRawEtFullAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtFullAcceptanceTPC",GetRawEtFullAcceptanceTPC(),1.0);
  if(GetRawEtEMCALAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtEMCALAcceptanceTPC",GetRawEtEMCALAcceptanceTPC(),1.0);
  if(GetRawEtPHOSAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtPHOSAcceptanceTPC",GetRawEtPHOSAcceptanceTPC(),1.0);
  if(GetRawEtFullAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtFullAcceptanceTPCNoPID",GetRawEtFullAcceptanceTPCNoPID(),1.0);
  if(GetRawEtEMCALAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtEMCALAcceptanceTPCNoPID",GetRawEtEMCALAcceptanceTPCNoPID(),1.0);
  if(GetRawEtPHOSAcceptanceTPCNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtPHOSAcceptanceTPCNoPID",GetRawEtPHOSAcceptanceTPCNoPID(),1.0);
  if(GetRawEtFullAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtFullAcceptanceITS",GetRawEtFullAcceptanceITS(),1.0);
  if(GetRawEtEMCALAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtEMCALAcceptanceITS",GetRawEtEMCALAcceptanceITS(),1.0);
  if(GetRawEtPHOSAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtPHOSAcceptanceITS",GetRawEtPHOSAcceptanceITS(),1.0);
  if(GetRawEtFullAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtFullAcceptanceITSNoPID",GetRawEtFullAcceptanceITSNoPID(),1.0);
  if(GetRawEtEMCALAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtEMCALAcceptanceITSNoPID",GetRawEtEMCALAcceptanceITSNoPID(),1.0);
  if(GetRawEtPHOSAcceptanceITSNoPID()>0.0 && fGoodEvent)FillHisto1D("RecoRawEtPHOSAcceptanceITSNoPID",GetRawEtPHOSAcceptanceITSNoPID(),1.0);
  if(fCentBin>-1){//if we have Pb+Pb and found a centrality bin
    if(GetCorrectedHadEtFullAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D(Form("RecoHadEtFullAcceptanceTPCCB%i",fCentBin),GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    if(GetCorrectedTotEtFullAcceptanceTPC()>0.0 && fGoodEvent)FillHisto1D(Form("RecoTotEtFullAcceptanceTPCCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    if(GetCorrectedHadEtFullAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D(Form("RecoHadEtFullAcceptanceITSCB%i",fCentBin),GetCorrectedHadEtFullAcceptanceITS(),1.0);
    if(GetCorrectedTotEtFullAcceptanceITS()>0.0 && fGoodEvent)FillHisto1D(Form("RecoTotEtFullAcceptanceITSCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceITS(),1.0);
    if(GetRawEtFullAcceptanceTPC()>0.0 && fGoodEvent)         FillHisto1D(Form("RecoRawEtFullAcceptanceTPCCB%i",fCentBin),GetRawEtFullAcceptanceTPC(),1.0);
    if(GetRawEtFullAcceptanceITS()>0.0 && fGoodEvent)         FillHisto1D(Form("RecoRawEtFullAcceptanceITSCB%i",fCentBin),GetRawEtFullAcceptanceITS(),1.0);
  }
  delete pID;
  delete strTPC;
  delete strITS;
  delete strTPCITS;
  // cout<<"Reconstructed pi/k/p et "<<GetCorrectedPiKPEtFullAcceptanceTPC()<<endl;
  return 1;
}
void AliAnalysisHadEtReconstructed::AddEt(Float_t rawEt, Float_t rawEtNoPID, Float_t corrEt, Float_t corrEtPion, Float_t corrEtProton, Float_t corrEtKaon, Float_t corrEtNoPID, Float_t pt, Bool_t IsTPC, Bool_t InPHOS, Bool_t InEMCAL) {//Adding Et to each of the variables that tracks et event by event
  if(pt>=AliAnalysisHadEt::fgPtTPCCutOff && IsTPC){//TPC tracks
    //adding to the raw Et
    fRawEtFullAcceptanceTPC += rawEt;
    if(InPHOS)fRawEtPHOSAcceptanceTPC += rawEt;
    if(InEMCAL)fRawEtEMCALAcceptanceTPC += rawEt;
    fRawEtFullAcceptanceTPCNoPID += rawEtNoPID;
    if(InPHOS)fRawEtPHOSAcceptanceTPCNoPID += rawEtNoPID;
    if(InEMCAL)fRawEtEMCALAcceptanceTPCNoPID += rawEtNoPID;
    //adding to the corrected Et
    fCorrectedHadEtFullAcceptanceTPC += corrEt;//the pi/k/p et
    fCorrectedHadEtFullAcceptanceTPCAssumingPion += corrEtPion;
    fCorrectedHadEtFullAcceptanceTPCAssumingProton += corrEtProton;
    fCorrectedHadEtFullAcceptanceTPCAssumingKaon += corrEtKaon;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceTPC += corrEt;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceTPC += corrEt;
    fCorrectedHadEtFullAcceptanceTPCNoPID += corrEtNoPID;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceTPCNoPID += corrEtNoPID;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceTPCNoPID += corrEtNoPID;
  }
  if(pt<AliAnalysisHadEt::fgPtTPCCutOff &&pt>=AliAnalysisHadEt::fgPtITSCutOff && !IsTPC){//ITS tracks
    //adding to the raw Et
    fRawEtFullAcceptanceITS += rawEt;
    if(InPHOS)fRawEtPHOSAcceptanceITS += rawEt;
    if(InEMCAL)fRawEtEMCALAcceptanceITS += rawEt;
    fRawEtFullAcceptanceITSNoPID += rawEtNoPID;
    if(InPHOS)fRawEtPHOSAcceptanceITSNoPID += rawEtNoPID;
    if(InEMCAL)fRawEtEMCALAcceptanceITSNoPID += rawEtNoPID;
    //adding to the corrected Et
    fCorrectedHadEtFullAcceptanceITS += corrEt;
    fCorrectedHadEtFullAcceptanceITSAssumingPion += corrEtPion;
    fCorrectedHadEtFullAcceptanceITSAssumingProton += corrEtProton;
    fCorrectedHadEtFullAcceptanceITSAssumingKaon += corrEtKaon;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceITS += corrEt;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceITS += corrEt;
    fCorrectedHadEtFullAcceptanceITSNoPID += corrEtNoPID;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceITSNoPID += corrEtNoPID;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceITSNoPID += corrEtNoPID;
  }
}

Bool_t AliAnalysisHadEtReconstructed::IsInPHOS(AliESDtrack *track){//This function will need to be elaborated on later to include PHOS dead channels
  if(!track){
    cout<<"Error: Track does not exist!!"<<endl;
    return kFALSE;
  }
  return   TMath::Abs(track->Eta()) < fCuts->GetGeometryPhosEtaAccCut()//in eta acceptance
    && track->Phi()*180.0/TMath::Pi() > fCuts->GetGeometryPhosPhiAccMinCut()//greater than the minimum phi
    && track->Phi()*180.0/TMath::Pi() < fCuts->GetGeometryPhosPhiAccMaxCut();//less than the maximum phi
}
Bool_t AliAnalysisHadEtReconstructed::IsInEMCAL(AliESDtrack *track){//This function will need to be elaborated on later to include EMCAL dead channels
  if(!track){
    cout<<"Error: Track does not exist!!"<<endl;
    return kFALSE;
  }
  return   TMath::Abs(track->Eta()) < fCuts->GetGeometryEmcalEtaAccCut()//in eta acceptance
    && track->Phi()*180.0/TMath::Pi() > fCuts->GetGeometryEmcalPhiAccMinCut()//greater than the minimum phi
    && track->Phi()*180.0/TMath::Pi() < fCuts->GetGeometryEmcalPhiAccMaxCut();//less than the maximum phi
}
Bool_t AliAnalysisHadEtReconstructed::CheckGoodVertex(AliVParticle* track)
{ // check vertex

  Float_t bxy = 999.;
  Float_t bz = 999.;
  if(!track){
    AliError("ERROR: no track");
    return kFALSE;
  }
  AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
  if(!esdTrack){
    AliError("ERROR: no track");
    return kFALSE;
  }
  esdTrack->GetImpactParametersTPC(bxy,bz);

  bool status = (TMath::Abs(track->Xv()) < fCuts->GetReconstructedVertexXCut()) && 
    (TMath::Abs(track->Yv()) < fCuts->GetReconstructedVertexYCut()) && 
    (TMath::Abs(track->Zv()) < fCuts->GetReconstructedVertexZCut()) && 
    (TMath::Abs(bxy) < fCuts->GetReconstructedIPxyCut()) && 
    (TMath::Abs(bz) < fCuts->GetReconstructedIPzCut()); 

  return status;
}

void AliAnalysisHadEtReconstructed::Init()
{ // Init
  AliAnalysisHadEt::Init();
  if(fCorrections){
    fCorrTotEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"Full");
    fCorrTotEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"Full");
    fCorrHadEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"Full");
    fCorrHadEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"Full");
    fCorrTotEtEMCALAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"EMCAL");
    fCorrTotEtEMCALAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"EMCAL");
    fCorrHadEtEMCALAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"EMCAL");
    fCorrHadEtEMCALAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"EMCAL");
    fCorrTotEtPHOSAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"PHOS");
    fCorrTotEtPHOSAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"PHOS");
    fCorrHadEtPHOSAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"PHOS");
    fCorrHadEtPHOSAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"PHOS");
  }
  else{
    cout<<"Warning!  You have not set corrections.  Your code will crash.  You have to set the corrections."<<endl;
  }
}

void AliAnalysisHadEtReconstructed::ResetEventValues(){//resetting event by event et's
  AliAnalysisHadEt::ResetEventValues();
  fCorrectedHadEtFullAcceptanceTPCNoPID=0.0;
  fCorrectedHadEtFullAcceptanceITSNoPID=0.0;
  fCorrectedHadEtEMCALAcceptanceTPCNoPID=0.0;
  fCorrectedHadEtEMCALAcceptanceITSNoPID=0.0;
  fCorrectedHadEtPHOSAcceptanceTPCNoPID=0.0;
  fCorrectedHadEtPHOSAcceptanceITSNoPID=0.0;
  fCorrectedHadEtFullAcceptanceTPC=0.0;
  fCorrectedHadEtFullAcceptanceITS=0.0;
  fCorrectedHadEtFullAcceptanceTPCAssumingPion=0.0;
  fCorrectedHadEtFullAcceptanceITSAssumingPion=0.0;
  fCorrectedHadEtFullAcceptanceTPCAssumingProton=0.0;
  fCorrectedHadEtFullAcceptanceITSAssumingProton=0.0;
  fCorrectedHadEtFullAcceptanceTPCAssumingKaon=0.0;
  fCorrectedHadEtFullAcceptanceITSAssumingKaon=0.0;
  fCorrectedHadEtEMCALAcceptanceTPC=0.0;
  fCorrectedHadEtEMCALAcceptanceITS=0.0;
  fCorrectedHadEtPHOSAcceptanceTPC=0.0;
  fCorrectedHadEtPHOSAcceptanceITS=0.0;
  fRawEtFullAcceptanceTPC=0.0;
  fRawEtFullAcceptanceITS=0.0;
  fRawEtEMCALAcceptanceTPC=0.0;
  fRawEtEMCALAcceptanceITS=0.0;
  fRawEtPHOSAcceptanceTPC=0.0;
  fRawEtPHOSAcceptanceITS=0.0;
  fRawEtFullAcceptanceTPCNoPID=0.0;
  fRawEtFullAcceptanceITSNoPID=0.0;
  fRawEtEMCALAcceptanceTPCNoPID=0.0;
  fRawEtEMCALAcceptanceITSNoPID=0.0;
  fRawEtPHOSAcceptanceTPCNoPID=0.0;
  fRawEtPHOSAcceptanceITSNoPID=0.0;

  if(TMath::Abs(fCorrTotEtFullAcceptanceTPC)<1e-3){
    if (fConfigFile.Length()) {
      cout<<"Warning: Rereading fCorrections file..."<<endl;
      gROOT->LoadMacro(fConfigFile);
      fCorrections = (AliAnalysisHadEtCorrections *) gInterpreter->ProcessLine("ConfigHadEtAnalysis()");
      fCorrTotEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"Full");
      fCorrTotEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"Full");
      fCorrHadEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"Full");
      fCorrHadEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"Full");
      fCorrTotEtEMCALAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"EMCAL");
      fCorrTotEtEMCALAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"EMCAL");
      fCorrHadEtEMCALAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"EMCAL");
      fCorrHadEtEMCALAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"EMCAL");
      fCorrTotEtPHOSAcceptanceTPC = fCorrections->GetConstantCorrections(kTRUE,fgPtTPCCutOff,"PHOS");
      fCorrTotEtPHOSAcceptanceITS = fCorrections->GetConstantCorrections(kTRUE,fgPtITSCutOff,"PHOS");
      fCorrHadEtPHOSAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"PHOS");
      fCorrHadEtPHOSAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"PHOS");
    }
    else{cerr<<"Uh-oh!  Unable to open configuration file!"<<endl;}
  }
}
void AliAnalysisHadEtReconstructed::CreateHistograms(){//Creating histograms and adding them to the output TList

  //TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  for(Int_t i=0;i<2;i++){
    TString *cutName = NULL;
    Float_t maxPtdEdx = 10;
    Float_t mindEdx = 35;
    Float_t maxdEdx = 150.0;
    switch(i){
    case 0:
      cutName = strTPCITS;
      break;
    case 1:
      cutName = strITS;
      maxPtdEdx = 5;
      maxdEdx = 500.0;
      break;
      //not deleting this completely since we might want to add it back in later
      //     case 2:
      //       cutName = strTPC;
      //       maxPtdEdx = 5;
      //       maxdEdx = 500.0;
      //       break;
    default:
      cerr<<"Error:  cannot make histograms!"<<endl;
      return;
    }

    CreateEtaPtHisto2D(Form("EtDataRaw%sPiPlus",cutName->Data()),"Raw reconstructed E_{T} from identified #pi^{+}");
    CreateEtaPtHisto2D(Form("EtDataRaw%sPiMinus",cutName->Data()),"Raw reconstructed E_{T} from identified #pi^{-}");
    CreateEtaPtHisto2D(Form("EtDataRaw%sKPlus",cutName->Data()),"Raw reconstructed E_{T} from identified K^{+}");
    //     CreateEtaPtHisto2D(Form("EtDataRaw%sEMinus",cutName->Data()),"Raw reconstructed E_{T} from identified e^{-}");
    //     CreateEtaPtHisto2D(Form("EtDataRaw%sEPlus",cutName->Data()),"Raw reconstructed E_{T} from identified e^{+}");
    CreateEtaPtHisto2D(Form("EtDataRaw%sKMinus",cutName->Data()),"Raw reconstructed E_{T} from identified K^{-}");
    CreateEtaPtHisto2D(Form("EtDataRaw%sProton",cutName->Data()),"Raw reconstructed E_{T} from identified p");
    CreateEtaPtHisto2D(Form("EtDataRaw%sAntiProton",cutName->Data()),"Raw reconstructed E_{T} from identified #bar{p}");
    CreateEtaPtHisto2D(Form("EtDataRaw%sUnidentified",cutName->Data()),"Raw reconstructed E_{T} from unidentified particles using real mass");
    CreateEtaPtHisto2D(Form("EtDataRaw%sNoID",cutName->Data()),"Raw reconstructed E_{T} from unidentified particles using real mass");

    CreateEtaPtHisto2D(Form("EtDataCorrected%sPiPlus",cutName->Data()),"Corrected reconstructed E_{T} from identified #pi^{+}");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sPiMinus",cutName->Data()),"Corrected reconstructed E_{T} from identified #pi^{-}");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sKPlus",cutName->Data()),"Corrected reconstructed E_{T} from identified K^{+}");
    //     CreateEtaPtHisto2D(Form("EtDataCorrected%sEMinus",cutName->Data()),"Corrected reconstructed E_{T} from identified e^{-}");
    //     CreateEtaPtHisto2D(Form("EtDataCorrected%sEPlus",cutName->Data()),"Corrected reconstructed E_{T} from identified e^{+}");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sKMinus",cutName->Data()),"Corrected reconstructed E_{T} from identified K^{-}");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sProton",cutName->Data()),"Corrected reconstructed E_{T} from identified p");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sAntiProton",cutName->Data()),"Corrected reconstructed E_{T} from identified #bar{p}");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sUnidentified",cutName->Data()),"Corrected reconstructed E_{T} from unidentified particles using real mass");
    CreateEtaPtHisto2D(Form("EtDataCorrected%sNoID",cutName->Data()),"Corrected reconstructed E_{T} from unidentified particles using real mass");


    CreateEtaPtHisto2D(Form("EtNData%sPiPlus",cutName->Data()),"Number of reconstructed #pi^{+}");
    CreateEtaPtHisto2D(Form("EtNData%sPiMinus",cutName->Data()),"Number of reconstructed #pi^{-}");
    CreateEtaPtHisto2D(Form("EtNData%sKPlus",cutName->Data()),"Number of reconstructed K^{+}");
    CreateEtaPtHisto2D(Form("EtNData%sKMinus",cutName->Data()),"Number of reconstructed K^{-}");
    CreateEtaPtHisto2D(Form("EtNData%sProton",cutName->Data()),"Number of reconstructed p");
    CreateEtaPtHisto2D(Form("EtNData%sAntiProton",cutName->Data()),"Number of reconstructed #bar{p}");
    CreateEtaPtHisto2D(Form("EtNData%sUnidentified",cutName->Data()),"Number of Reconstructed unidentified particles");

    CreateHisto2D(Form("dEdxDataAll%s",cutName->Data()),"dE/dx for all particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxDataPion%s",cutName->Data()),"dE/dx for #pi^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxDataKaon%s",cutName->Data()),"dE/dx for K^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxDataProton%s",cutName->Data()),"dE/dx for p(#bar{p})","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxDataElectron%s",cutName->Data()),"dE/dx for e^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxDataUnidentified%s",cutName->Data()),"dE/dx for unidentified particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
  }

  Float_t minEt = 0.0;
  Float_t maxEt = 100.0;
  if(fDataSet==20100) maxEt=4000.0;
  Int_t nbinsEt = 200;
  char histoname[200];
  char histotitle[200];
  char xtitle[50];
  TString *ytitle = new TString("Number of events");
  TString *sTPC = new TString("TPC");
  TString *sITS = new TString("ITS");
  TString *sTPCpt = new TString("0.15");
  TString *sITSpt = new TString("0.10");
  TString *sPID = new TString("");
  TString *sNoPID = new TString("NoPID");
  TString *sNoPIDString = new TString(", No PID");
  TString *sHadEt = new TString("HadEt");
  TString *sRawEt = new TString("RawEt");
  TString *sTotEt = new TString("TotEt");
  TString *sTotEtString = new TString("total E_{T}");
  TString *sHadEtString = new TString("hadronic E_{T}");
  TString *sRawEtString = new TString("raw E_{T}");
  TString *sFull = new TString("Full");
  TString *sEMCAL = new TString("EMCAL");
  TString *sPHOS = new TString("PHOS");
  
  for(int tpc = 0;tpc<2;tpc++){
    for(int hadet = 0;hadet<3;hadet++){
      for(int type = 0;type<3;type++){
	for(int pid = 0;pid<2;pid++){
	  TString *detector = NULL;
	  TString *partid = NULL;
	  TString *et = sHadEt;
	  TString *acceptance = NULL;
	  TString *ptstring = NULL;
	  TString *partidstring = NULL;
	  TString *etstring = sHadEtString;
	  if(tpc==1) {detector = sTPC; ptstring = sTPCpt;}
	  else{detector = sITS; ptstring = sITSpt;}
	  if(pid==1){partid = sPID; partidstring = sPID;}
	  else{partid = sNoPID; partidstring = sNoPIDString;}
	  if(hadet==1) {et = sHadEt; etstring = sHadEtString;}
	  if(hadet==0){et = sTotEt; etstring = sTotEtString;}
	  if(hadet==2){et = sRawEt; etstring = sRawEtString;}
	  switch(type){
	  case 0:
	    acceptance = sFull;
	    break;
	  case 1:
	    acceptance = sEMCAL;
	    break;
	  case 2:
	    acceptance = sPHOS;
	    break;
	  default:
	    acceptance = sFull;
	  }
	  snprintf(histoname,200,"Reco%s%sAcceptance%s%s",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	  snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	  snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	  CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	  if(type==0){//full acceptance only
	    snprintf(xtitle,50,"Reconstructed %s",etstring->Data());

	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sND",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed non-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sDD",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed doubly-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sSD",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed singly-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	  }
	  if(fDataSet==20100 && type ==0 &&pid==1){//If this is Pb+Pb and full acceptance with pid
	    Int_t width = 5;
	    if(fNCentBins<21) width = 10;
	    for(Int_t i=0;i<fNCentBins;i++){
	      snprintf(histoname,200,"Reco%s%sAcceptance%s%sCB%i",et->Data(),acceptance->Data(),detector->Data(),partid->Data(),i);
	      snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s for centrality %i-%i",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data(),i*width,(i+1)*width);
	      snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	      CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	    }
	  }
	}
      }
    }
  }
  CreateHisto1D("RecoHadEtFullAcceptanceTPCAssumingPion","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.15 GeV/c assuming pions","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);
  CreateHisto1D("RecoHadEtFullAcceptanceTPCAssumingProton","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.15 GeV/c assuming protons","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);
  CreateHisto1D("RecoHadEtFullAcceptanceTPCAssumingKaon","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.15 GeV/c assuming kaons","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);
  CreateHisto1D("RecoHadEtFullAcceptanceITSAssumingPion","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.10 GeV/c assuming pions","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);
  CreateHisto1D("RecoHadEtFullAcceptanceITSAssumingProton","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.10 GeV/c assuming protons","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);
  CreateHisto1D("RecoHadEtFullAcceptanceITSAssumingKaon","Reconstructing E_{T}^{had} with full acceptance for p_{T}>0.10 GeV/c assuming kaons","Reconstructed E_{T}^{had}","dN_{eve}/dE_{T}^{had}",nbinsEt*2,minEt,maxEt);

  //Cross checks that corrections are applied correctly
  if(fDataSet==20100){
    CreateHisto2D("fbkgdVsCentralityBin","f_{bkgd} vs centrality bin","centrality bin","f_{bkgd}",21,-1.5,19.5,200,0.7,1.05);//
    CreateHisto2D("feffPionVsCentralityBin","Pion efficiency vs centrality bin","centrality bin","pion efficiency",21,-1.5,19.5,200,0,1.2);//
    CreateHisto2D("feffHadronVsCentralityBin","Hadron efficiency vs centrality bin","centrality bin","hadron efficiency",21,-1.5,19.5,200,0,1.2);//
    CreateHisto2D("feffKaonVsCentralityBin","Kaon efficiency vs centrality bin","centrality bin","kaon efficiency",21,-1.5,19.5,200,0,1.2);//
    CreateHisto2D("feffProtonVsCentralityBin","Proton efficiency vs centrality bin","centrality bin","proton efficiency",21,-1.5,19.5,200,0,1.2);//
    CreateHisto2D("fnotIDVsCentralityBin","f_{notID} vs centrality bin","centrality bin","f_{notID}",21,-1.5,19.5,50,0.95,1.05);//
    CreateHisto2D("fpTcutVsCentralityBin","f_{pTcut} vs centrality bin","centrality bin","f_{pTcut}",21,-1.5,19.5,50,0.95,1.05);
    CreateHisto2D("fneutralVsCentralityBin","f_{neutral} vs centrality bin","centrality bin","f_{neutral}",21,-1.5,19.5,50,0.5,1.00);
    CreateHisto2D("ConstantCorrectionsVsCentralityBin","constant corrections vs centrality bin","centrality bin","constant corrections",21,-1.5,19.5,50,0.5,1.00);
  }

  delete sTPC;
  delete sITS;
  delete sTPCpt;
  delete sITSpt;
  delete sPID;
  delete sNoPID;
  delete sNoPIDString;
  delete sHadEt;
  delete sTotEt;
  delete sTotEtString;
  delete sHadEtString;
  delete sFull;
  delete sEMCAL;
  delete sPHOS;

}
