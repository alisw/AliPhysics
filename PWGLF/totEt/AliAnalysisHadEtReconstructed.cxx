//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for ESD analysis
//  - reconstruction output
// implementation file
//
//Created by Christine Nattrass
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
#include "AliMultSelection.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h" 
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliPWG0Helper.h"
#include "TMath.h"

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
  ,fCorrPiKPEtFullAcceptanceTPC(0)
  ,fCorrPiKPEtFullAcceptanceITS(0)
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
  ,fCorrectedHadEtFullAcceptanceTPCNegEta(0)
  ,fCorrectedHadEtFullAcceptanceTPCNoPIDNegEta(0)
  ,fCorrectedHadEtFullAcceptanceTPCPosEta(0)
  ,fCorrectedHadEtFullAcceptanceTPCNoPIDPosEta(0)
  ,fCorrectedHadEtFullAcceptanceTPCLimitedPhi(0)
  ,fCorrectedHadEtFullAcceptanceTPCNoPIDLimitedPhi(0)
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
  ,kIsOfflineV0AND(0)
  ,kDoTriggerChecks(0)
  ,kDoTriggerChecksOnly(0)
  ,useOldCentrality(0)
{
}

AliAnalysisHadEtReconstructed::~AliAnalysisHadEtReconstructed() 
{
  delete fCorrections;
}

Int_t AliAnalysisHadEtReconstructed::AnalyseEvent(AliVEvent* ev, Int_t eventtype)
{ // analyse ESD event
  if(kDoTriggerChecksOnly){return 1;}//In this case we are just after trigger efficiencies and don't care about the ET reconstructed.
  if(kDoTriggerChecks && !kIsOfflineV0AND){return 1;}//In this case we are just after trigger efficiencies and don't care about the ET reconstructed.
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
//   if(fDataSet==20100||fDataSet==2011){//If this is Pb+Pb or pPb
// //     AliCentrality *centrality = realEvent->GetCentrality();
// //     if(fNCentBins<21) fCentBin= centrality->GetCentralityClass10(fCentralityMethod);
// //     else{ fCentBin= centrality->GetCentralityClass5(fCentralityMethod);}
//     AliCentrality *centrality =  realEvent->GetCentrality();
//     fCentBin = GetCentralityBin(fNCentBins, centrality);
//     if(fCentBin ==-1){
//       if(fDataSet==2013){
// 	fCentBin = 19;//For pPb we don't want to throw these events out but there is no CB 19
//       }
//       else{
// 	fGoodEvent = kFALSE;//but for Pb+Pb events we don't want to count events where we did not find a centrality
//       }
//     }
//   }
     //if( fDataSet==2015){
 if(fDataSet==20100||fDataSet==2011 ||  fDataSet==2015){//If this is Pb+Pb or pPb
   AliCentrality *centrality = realEvent->GetCentrality();//if the centrality task exists, use it!
   if(centrality && useOldCentrality){
     //if(centrality){
    fCentBin = GetCentralityBin(fNCentBins, centrality);
   }
   //else{
   AliMultSelection *MultSelection = (AliMultSelection * ) realEvent->FindListObject("MultSelection");
   if(MultSelection && !useOldCentrality){//if the centrality class returns nothing it still exists!
     fCentBin = GetCentralityBin(fNCentBins, MultSelection);
     if(fCentBin ==-1){
       fGoodEvent = kFALSE;//but for Pb+Pb events we don't want to count events where we did not find a centrality
       if(fDataSet==2013){
	 fCentBin = 19;//For pPb we don't want to throw these events out but there is no CB 19
       }
       else{
	 fGoodEvent = kFALSE;//but for Pb+Pb events we don't want to count events where we did not find a centrality
       }
     }
     }
   }
 //}
  //for PID
//   AliESDpid *pID = new AliESDpid();
//   pID->MakePID(realEvent);
  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  //cerr<<"Do track cuts exist reco?  "<<fEsdtrackCutsITSTPC<<", "<<fEsdtrackCutsITS<<", "<<fEsdtrackCutsTPC<<endl;
  if(!fEsdtrackCutsITSTPC || !fEsdtrackCutsITSTPC ||!fEsdtrackCutsITSTPC){
    //gROOT->LoadMacro("ConfigEtReconstructed.C");
    fEsdtrackCutsITSTPC = (AliESDtrackCuts *)fhistoList->FindObject("fEsdTrackCuts");// gInterpreter->ProcessLine("SetTrackCutsITSTPC()");
    fEsdtrackCutsTPC = (AliESDtrackCuts *)  fhistoList->FindObject("fEsdTrackCutsTPCOnly");//  gInterpreter->ProcessLine("SetTrackCutsTPC()");
    fEsdtrackCutsITS = (AliESDtrackCuts *)  fhistoList->FindObject("fEsdTrackCutsITS");//  gInterpreter->ProcessLine("SetTrackCutsITS()");
    //cerr<<"Remade track cuts.  Do track cuts exist reco?  "<<fEsdtrackCutsITSTPC<<", "<<fEsdtrackCutsITS<<", "<<fEsdtrackCutsTPC<<endl;
  }
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
	  if(fCentBin>=0&&cutset==0){
	    FillHisto1D(Form("SpectraCB%i",fCentBin),track->Pt(),1.0); 
	  }
	  Float_t nSigmaPion,nSigmaProton,nSigmaKaon,nSigmaElectron;
	  Float_t nSigmaPionUnsigned,nSigmaProtonUnsigned,nSigmaKaonUnsigned,nSigmaElectronUnsigned;
// 	  pID->MakeTPCPID(track);
// 	  pID->MakeITSPID(track);
	  //if(!fPIDResponse) cout<<"Uh-oh!  No PID Response!"<<endl;
	  if(cutset!=1){
	    nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 
	    nSigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)); 
	    nSigmaKaon =TMath::Abs( fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
	    nSigmaElectron =TMath::Abs( fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));

	    nSigmaPionUnsigned = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion); 
	    nSigmaProtonUnsigned = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton); 
	    nSigmaKaonUnsigned = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon); 
	    nSigmaElectronUnsigned = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 

// 	    nSigmaPion = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kPion));
// 	    nSigmaProton = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kProton));
// 	    nSigmaKaon = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kKaon));
// 	    nSigmaElectron = TMath::Abs(pID->NumberOfSigmasTPC(track,AliPID::kElectron));
	  }
	  else{
	    nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kPion)); 
	    nSigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton)); 
	    nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon)); 
	    nSigmaElectron = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kElectron)); 

	    nSigmaPionUnsigned = fPIDResponse->NumberOfSigmasITS(track, AliPID::kPion); 
	    nSigmaProtonUnsigned = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton); 
	    nSigmaKaonUnsigned = fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon); 
	    nSigmaElectronUnsigned = fPIDResponse->NumberOfSigmasITS(track, AliPID::kElectron); 
// 	    nSigmaPion = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kPion));
// 	    nSigmaProton = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kProton));
// 	    nSigmaKaon = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kKaon));
// 	    nSigmaElectron = TMath::Abs(pID->NumberOfSigmasITS(track,AliPID::kElectron));
	  }
	  //cout<<"Nsigma pion "<<nSigmaPion<<" proton "<<nSigmaProton<<" kaon "<<nSigmaKaon<<" electron "<<nSigmaElectron<<endl;
	  // 	    bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  // 	    bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  // 	    bool isKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
	  // 	    bool isProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);
	  bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  bool isKaon = (nSigmaPion>3.0 && nSigmaProton>3.0 && nSigmaKaon<3.0 && track->Pt()<0.45);
	  bool isProton = (nSigmaPion>3.0 && nSigmaProton<3.0 && nSigmaKaon>3.0 && track->Pt()<0.9);

	  FillHisto2D(Form("dEdxDataNSigmaPionAll%s",cutName->Data()),track->P(),nSigmaPionUnsigned,1.0);
	  FillHisto2D(Form("dEdxDataNSigmaKaonAll%s",cutName->Data()),track->P(),nSigmaKaonUnsigned,1.0);
	  FillHisto2D(Form("dEdxDataNSigmaProtonAll%s",cutName->Data()),track->P(),nSigmaProtonUnsigned,1.0);
	  FillHisto2D(Form("dEdxDataNSigmaElectronAll%s",cutName->Data()),track->P(),nSigmaElectronUnsigned,1.0);
	  
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

	    //if(fCentBin>0){FillHisto1D(Form("SpectraCB%i",fCentBin),track->Pt(),1.0); }

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
	  if(fDataSet==20100||fDataSet==2015||fDataSet==2011){
	    FillHisto2D("rawspectra",fCentBin,track->Pt(),1.0);
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
	    FillHisto2D(Form("dEdxDataNSigmaPion%s",cutName->Data()),track->P(),nSigmaPionUnsigned,1.0);
	    et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
   	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionPion(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&(fDataSet==20100||fDataSet==2015||fDataSet==2011))FillHisto2D("feffPionVsCentralityBin",fCentBin,1.0/corrEff,1.0);
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
	    FillHisto2D(Form("dEdxDataNSigmaKaon%s",cutName->Data()),track->P(),nSigmaKaonUnsigned,1.0);
	    et = Et(track->P(),track->Theta(),fgKPlusCode,track->Charge());
	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionKaon(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&(fDataSet==20100||fDataSet==2015||fDataSet==2011))FillHisto2D("feffKaonVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	      
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
	    FillHisto2D(Form("dEdxDataNSigmaProton%s",cutName->Data()),track->P(),nSigmaProtonUnsigned,1.0);
	    et = Et(track->P(),track->Theta(),fgProtonCode,track->Charge());
	    if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionProton(track->Pt(),fCentBin);}
	    etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	    if(corrEff>0.0&&(fDataSet==20100||fDataSet==2015||fDataSet==2011))FillHisto2D("feffProtonVsCentralityBin",fCentBin,1.0/corrEff,1.0);
	      
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
	    FillHisto2D(Form("dEdxDataNSigmaElectron%s",cutName->Data()),track->P(),nSigmaElectronUnsigned,1.0);
	  }
	  if(unidentified){
	    if(isPion) cerr<<"I should not be here!!  AliAnalysisHadEtReconstructed 273"<<endl; 
	    FillHisto2D(Form("dEdxDataUnidentified%s",cutName->Data()),track->P(),dEdx,1.0);
	    et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	    Float_t etProton = Et(track->P(),track->Theta(),fgProtonCode,track->Charge());
	    Float_t etKaon = Et(track->P(),track->Theta(),fgKPlusCode,track->Charge());
	    if(corrEff>0.0&&(fDataSet==20100||fDataSet==2015||fDataSet==2011))FillHisto2D("feffHadronVsCentralityBin",fCentBin,1.0/corrEff,1.0);
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
	  AddEt(et,etNoID,etpartialcorrected,etpartialcorrectedPion,etpartialcorrectedProton,etpartialcorrectedKaon,etpartialcorrectedNoID,track->Pt(),isTPC,inPHOS,inEMCAL,track->Phi(),track->Eta());
	  if(fCentBin>0) FillHisto2D(Form("ETvsPhiAndEtaCB%i",fCentBin),track->Phi(),track->Eta(),1.0);

	}
      }
    delete list;
  }

  Int_t nondiff = 0;//(Int_t) AliPWG0Helper::kND;
  Int_t doublediff = 0;//(Int_t) AliPWG0Helper::kDD;
  Int_t singlediff = 0;//(Int_t) AliPWG0Helper::kSD;
  if(fDataSet!=20100 && fDataSet!=2015 && fDataSet!=2011){
    nondiff = (Int_t) AliPWG0Helper::kND;
    doublediff = (Int_t) AliPWG0Helper::kDD;
    singlediff = (Int_t) AliPWG0Helper::kSD;
  }
//  cout<<"event type "<<eventtype<<" nondiff event type "<<nondiff<<" data set "<<fDataSet<<" good event "<<fGoodEvent<<endl;
  if((eventtype == nondiff|| fDataSet==20100|| fDataSet==2015 || fDataSet==2011)  && fGoodEvent){
    //cout<<"Filling "<<endl;
    FillHisto1D("RecoHadEtFullAcceptanceTPCND",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCND",GetCorrectedPiKPEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCND",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSND",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSND",GetCorrectedPiKPEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSND",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDND",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCNoPIDND",GetCorrectedPiKPEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDND",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDND",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSNoPIDND",GetCorrectedPiKPEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDND",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);

    FillHisto1D("RecoRawEtFullAcceptanceTPCND",GetRawEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCNoPIDND",GetRawEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSND",GetRawEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSNoPIDND",GetRawEtFullAcceptanceITSNoPID(),1.0);

  }
  if(eventtype == doublediff){
    FillHisto1D("RecoHadEtFullAcceptanceTPCDD",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCDD",GetCorrectedPiKPEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCDD",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSDD",GetCorrectedPiKPEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSDD",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSDD",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDDD",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCNoPIDDD",GetCorrectedPiKPEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDDD",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDDD",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSNoPIDDD",GetCorrectedPiKPEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDDD",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCDD",GetRawEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCNoPIDDD",GetRawEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSDD",GetRawEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSNoPIDDD",GetRawEtFullAcceptanceITSNoPID(),1.0);//ND
  }
  if(eventtype == singlediff){
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCSD",GetCorrectedPiKPEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCSD",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCSD",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSSD",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSSD",GetCorrectedPiKPEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSSD",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPIDSD",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCNoPIDSD",GetCorrectedPiKPEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDSD",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPIDSD",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSNoPIDSD",GetCorrectedPiKPEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPIDSD",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCSD",GetRawEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCNoPIDSD",GetRawEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSSD",GetRawEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSNoPIDSD",GetRawEtFullAcceptanceITSNoPID(),1.0);
  }
  if(fGoodEvent){
    FillHisto1D("RecoPiKPEtFullAcceptanceTPC",GetCorrectedPiKPEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITS",GetCorrectedPiKPEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceTPCNoPID",GetCorrectedPiKPEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoPiKPEtFullAcceptanceITSNoPID",GetCorrectedPiKPEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPC",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPC",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingPion",GetCorrectedHadEtFullAcceptanceTPCAssumingPion(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingProton",GetCorrectedHadEtFullAcceptanceTPCAssumingProton(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCAssumingKaon",GetCorrectedHadEtFullAcceptanceTPCAssumingKaon(),1.0);
    FillHisto1D("RecoHadEtEMCALAcceptanceTPC",GetCorrectedHadEtEMCALAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtEMCALAcceptanceTPC",GetCorrectedTotEtEMCALAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtPHOSAcceptanceTPC",GetCorrectedHadEtPHOSAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtPHOSAcceptanceTPC",GetCorrectedTotEtPHOSAcceptanceTPC(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceTPCNoPID",GetCorrectedHadEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPID",GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtEMCALAcceptanceTPCNoPID",GetCorrectedHadEtEMCALAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtEMCALAcceptanceTPCNoPID",GetCorrectedTotEtEMCALAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtPHOSAcceptanceTPCNoPID",GetCorrectedHadEtPHOSAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoTotEtPHOSAcceptanceTPCNoPID",GetCorrectedTotEtPHOSAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITS",GetCorrectedHadEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSAssumingPion",GetCorrectedHadEtFullAcceptanceITSAssumingPion(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSAssumingProton",GetCorrectedHadEtFullAcceptanceITSAssumingProton(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSAssumingKaon",GetCorrectedHadEtFullAcceptanceITSAssumingKaon(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITS",GetCorrectedTotEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtEMCALAcceptanceITS",GetCorrectedHadEtEMCALAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtEMCALAcceptanceITS",GetCorrectedTotEtEMCALAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtPHOSAcceptanceITS",GetCorrectedHadEtPHOSAcceptanceITS(),1.0);
    FillHisto1D("RecoTotEtPHOSAcceptanceITS",GetCorrectedTotEtPHOSAcceptanceITS(),1.0);
    FillHisto1D("RecoHadEtFullAcceptanceITSNoPID",GetCorrectedHadEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceITSNoPID",GetCorrectedTotEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoHadEtEMCALAcceptanceITSNoPID",GetCorrectedHadEtEMCALAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtEMCALAcceptanceITSNoPID",GetCorrectedTotEtEMCALAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoHadEtPHOSAcceptanceITSNoPID",GetCorrectedHadEtPHOSAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtPHOSAcceptanceITSNoPID",GetCorrectedTotEtPHOSAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNegEta",GetCorrectedTotEtFullAcceptanceTPCNegEta(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCPosEta",GetCorrectedTotEtFullAcceptanceTPCPosEta(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCLimitedPhi",GetCorrectedTotEtFullAcceptanceTPCLimitedPhi(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDNegEta",GetCorrectedTotEtFullAcceptanceTPCNoPIDNegEta(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDPosEta",GetCorrectedTotEtFullAcceptanceTPCNoPIDPosEta(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPCNoPIDLimitedPhi",GetCorrectedTotEtFullAcceptanceTPCNoPIDLimitedPhi(),1.0);

    FillHisto1D("RecoRawEtFullAcceptanceTPC",GetRawEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtEMCALAcceptanceTPC",GetRawEtEMCALAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtPHOSAcceptanceTPC",GetRawEtPHOSAcceptanceTPC(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceTPCNoPID",GetRawEtFullAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtEMCALAcceptanceTPCNoPID",GetRawEtEMCALAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtPHOSAcceptanceTPCNoPID",GetRawEtPHOSAcceptanceTPCNoPID(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITS",GetRawEtFullAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtEMCALAcceptanceITS",GetRawEtEMCALAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtPHOSAcceptanceITS",GetRawEtPHOSAcceptanceITS(),1.0);
    FillHisto1D("RecoRawEtFullAcceptanceITSNoPID",GetRawEtFullAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoRawEtEMCALAcceptanceITSNoPID",GetRawEtEMCALAcceptanceITSNoPID(),1.0);
    FillHisto1D("RecoRawEtPHOSAcceptanceITSNoPID",GetRawEtPHOSAcceptanceITSNoPID(),1.0);
    if(fCentBin>-1){//if we have Pb+Pb and found a centrality bin
      FillHisto1D(Form("RecoHadEtFullAcceptanceTPCCB%i",fCentBin),GetCorrectedHadEtFullAcceptanceTPC(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPC(),1.0);
      FillHisto1D(Form("RecoHadEtFullAcceptanceITSCB%i",fCentBin),GetCorrectedHadEtFullAcceptanceITS(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceITSCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceITS(),1.0);
      FillHisto1D(Form("RecoRawEtFullAcceptanceTPCCB%i",fCentBin),GetRawEtFullAcceptanceTPC(),1.0);
      FillHisto1D(Form("RecoRawEtFullAcceptanceITSCB%i",fCentBin),GetRawEtFullAcceptanceITS(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCNegEtaCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCNegEta(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCPosEtaCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCPosEta(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCLimitedPhiCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCLimitedPhi(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCNoPIDCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCNoPID(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCNoPIDNegEtaCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCNoPIDNegEta(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCNoPIDPosEtaCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCNoPIDPosEta(),1.0);
      FillHisto1D(Form("RecoTotEtFullAcceptanceTPCNoPIDLimitedPhiCB%i",fCentBin),GetCorrectedTotEtFullAcceptanceTPCNoPIDLimitedPhi(),1.0);
    }
  }
//   delete pID;
  delete strTPC;
  delete strITS;
  delete strTPCITS;
  // cout<<"Reconstructed pi/k/p et "<<GetCorrectedPiKPEtFullAcceptanceTPC()<<endl;
  return 1;
}
void AliAnalysisHadEtReconstructed::AddEt(Float_t rawEt, Float_t rawEtNoPID, Float_t corrEt, Float_t corrEtPion, Float_t corrEtProton, Float_t corrEtKaon, Float_t corrEtNoPID, Float_t pt, Bool_t IsTPC, Bool_t InPHOS, Bool_t InEMCAL,Float_t phi, Float_t eta) {//Adding Et to each of the variables that tracks et event by event
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
    if(eta<0){
      fCorrectedHadEtFullAcceptanceTPCNegEta += corrEt;//the pi/k/p et
      fCorrectedHadEtFullAcceptanceTPCNoPIDNegEta += corrEtNoPID;
    }
    else{
      fCorrectedHadEtFullAcceptanceTPCPosEta += corrEt;//the pi/k/p et
      fCorrectedHadEtFullAcceptanceTPCNoPIDPosEta += corrEtNoPID;
    }
    //exclude these regions
    //(0.3<track->Phi()<1.15 && 2<track->Phi()<2.5 && 3.7<track->Phi()<4.6)
    //one could write this in one line but it might be wiser to write code which is less buggy and easier to follow...
    //cout<<"phi "<<phi<<endl;
    if(phi<0.3 || phi>1.15 ){
      if(phi<2 || phi >2.5){
	if(phi<3.7 || phi > 4.6){
	  //cout<<"phi "<<phi<<endl;
	  fCorrectedHadEtFullAcceptanceTPCLimitedPhi += corrEt;//the pi/k/p et
	  fCorrectedHadEtFullAcceptanceTPCNoPIDLimitedPhi += corrEtNoPID;
	}
      }
    }
  }
  //if(pt<AliAnalysisHadEt::fgPtTPCCutOff &&pt>=AliAnalysisHadEt::fgPtITSCutOff && !IsTPC){//ITS tracks
  //If we use standalone tracks - not pure standalone tracks - the only tracks we get are ones that were missed by the TPC+ITS tracking.  Therefore we don't need to add a momentum cut-off
  if(pt<AliAnalysisHadEt::fgPtTPCCutOff && !IsTPC){//ITS tracks
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
    fCorrPiKPEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"PiKP");
    fCorrPiKPEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"PiKP");
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
  fCorrectedHadEtFullAcceptanceTPCNegEta=0;
  fCorrectedHadEtFullAcceptanceTPCNoPIDNegEta=0;
  fCorrectedHadEtFullAcceptanceTPCPosEta=0;
  fCorrectedHadEtFullAcceptanceTPCNoPIDPosEta=0;
  fCorrectedHadEtFullAcceptanceTPCLimitedPhi=0;
  fCorrectedHadEtFullAcceptanceTPCNoPIDLimitedPhi=0;
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
      fCorrPiKPEtFullAcceptanceTPC = fCorrections->GetConstantCorrections(kFALSE,fgPtTPCCutOff,"PiKP");
      fCorrPiKPEtFullAcceptanceITS = fCorrections->GetConstantCorrections(kFALSE,fgPtITSCutOff,"PiKP");
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
  Float_t maxCentbinRange = fNCentBins+0.5;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man) {
    AliFatal("Analysis manager needed");
    return;
  }
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliError("PIDResponse object was not created");


  if(kDoTriggerChecksOnly){return;}//In this case we are just after trigger efficiencies and don't care about the ET reconstructed.
  //TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  Float_t minNSigma = -4;
  Float_t maxNSigma = 4;
  Float_t maxPtdEdxNSigma = 1.0;

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


    CreateHisto2D(Form("dEdxDataNSigmaPion%s",cutName->Data()),"N_{#sigma} dE/dx for PID'd #pi^{#pm}","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaKaon%s",cutName->Data()),"N_{#sigma} dE/dx for PID'd K^{#pm}","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaProton%s",cutName->Data()),"N_{#sigma} dE/dx for PID'd p(#bar{p})","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaElectron%s",cutName->Data()),"N_{#sigma} dE/dx for e^{#pm}","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);

    CreateHisto2D(Form("dEdxDataNSigmaPionAll%s",cutName->Data()),"N_{#sigma} dE/dx for #pi^{#pm} for all","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaKaonAll%s",cutName->Data()),"N_{#sigma} dE/dx for K^{#pm} for all","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaProtonAll%s",cutName->Data()),"N_{#sigma} dE/dx for p(#bar{p}) for all","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);
    CreateHisto2D(Form("dEdxDataNSigmaElectronAll%s",cutName->Data()),"N_{#sigma} dE/dx for e^{#pm} for all","momentum (GeV/c)","N_{#sigma} dE/dx",50,0.0,maxPtdEdxNSigma,100,minNSigma,maxNSigma);


  }

  Float_t minEt = 0.0;
  Float_t maxEt = 100.0;
  Float_t minEtPiKP = 0.0;
  Float_t maxEtPiKP = 100.0;
  if(fDataSet==20100||fDataSet==2011){
    maxEt=4000.0;
    maxEtPiKP = 2500;
  }
  if(fDataSet==2013){
    maxEt=100.0;
    maxEtPiKP = 100.0;
  }
  if(fDataSet==2015){
    maxEt=6000.0;
    maxEtPiKP = 3500;
  }
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
  TString *sPiKPEt = new TString("PiKPEt");
  TString *sRawEt = new TString("RawEt");
  TString *sTotEt = new TString("TotEt");
  TString *sTotEtString = new TString("total E_{T}");
  TString *sHadEtString = new TString("hadronic E_{T}");
  TString *sPiKPEtString = new TString("E_{T}^{#pi,K,p}");
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

	    snprintf(histoname,200,"Reco%s%sAcceptance%s%s",sPiKPEt->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s",sPiKPEtString->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    snprintf(xtitle,50,"Reconstructed %s",sPiKPEtString->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEtPiKP,maxEtPiKP);

	    snprintf(xtitle,50,"Reconstructed %s",sPiKPEt->Data());
	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sND",sPiKPEt->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed non-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",sPiKPEtString->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEtPiKP,maxEtPiKP);
	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sDD",sPiKPEt->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed doubly-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",sPiKPEtString->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEtPiKP,maxEtPiKP);
	    snprintf(histoname,200,"Reco%s%sAcceptance%s%sSD",sPiKPEt->Data(),acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"Reconstructed singly-diffractive events %s with %s acceptance for p_{T}>%s GeV/c%s",sPiKPEtString->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	    CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEtPiKP,maxEtPiKP);//et
	    if(tpc==1 && hadet==0){//full acceptance and TPC and tot ET
	      snprintf(histoname,200,"Reco%s%sAcceptance%s%sNegEta",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	      snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using #eta<0",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	      snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	      CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);

	      snprintf(histoname,200,"Reco%s%sAcceptance%s%sPosEta",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	      snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using #eta>0",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	      snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	      CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);

	      snprintf(histoname,200,"Reco%s%sAcceptance%s%sLimitedPhi",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	      snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using restricted #phi",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	      snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	      CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	    }

	  }
	  //if((fDataSet==20100||fDataSet==2015||fDataSet==2011) && type ==0 &&pid==1){//If this is Pb+Pb and full acceptance with pid
	  if((fDataSet==20100||fDataSet==2015||fDataSet==2011) && type ==0){//If this is Pb+Pb and full acceptance
	    Int_t width = 5;
	    if(fNCentBins<21) width = 10;
	    for(Int_t i=0;i<fNCentBins;i++){
	      snprintf(histoname,200,"Reco%s%sAcceptance%s%sCB%i",et->Data(),acceptance->Data(),detector->Data(),partid->Data(),i);
	      snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s for centrality %i-%i",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data(),i*width,(i+1)*width);
	      snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
	      CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	      if(tpc==1 && hadet==0){//full acceptance and TPC and tot ET

		if(pid==1){//only want to make one for each centrality bin
		  //void CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh);
		  CreateHisto2D(Form("ETvsPhiAndEtaCB%i",i),"E_{T} vs #phi and #eta","#phi","#eta",200,0.0,TMath::Pi()*2,200,-0.7,0.7);

		  snprintf(histoname,200,"SpectraCB%i",i);
		  snprintf(histotitle,200,"Spectra in CB %i",i);
		  snprintf(xtitle,50,"Spectra %s",etstring->Data());
		  CreatePtSpectraHisto1D(histoname,histotitle,"p_{T}","Number of particles");

		}

		snprintf(histoname,200,"Reco%s%sAcceptance%s%sNegEtaCB%i",et->Data(),acceptance->Data(),detector->Data(),partid->Data(),i);
		snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using #eta<0 for centrality %i-%i",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data(),i*width,(i+1)*width);
		snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
		CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
		snprintf(histoname,200,"Reco%s%sAcceptance%s%sPosEtaCB%i",et->Data(),acceptance->Data(),detector->Data(),partid->Data(),i);
		snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using #eta>0 for centrality %i-%i",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data(),i*width,(i+1)*width);
		snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
		CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
		snprintf(histoname,200,"Reco%s%sAcceptance%s%sLimitedPhiCB%i",et->Data(),acceptance->Data(),detector->Data(),partid->Data(),i);
		snprintf(histotitle,200,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s using limited #phi for centrality %i-%i",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data(),i*width,(i+1)*width);
		snprintf(xtitle,50,"Reconstructed %s",etstring->Data());
		CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);

	      }
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
  if(fDataSet==20100 || fDataSet==2015 || fDataSet==2011){
    CreateHisto2D("fbkgdVsCentralityBin","f_{bkgd} vs centrality bin","centrality bin","f_{bkgd}",fNCentBins,-1.5,maxCentbinRange,200,0.7,1.05);//
    CreateHisto2D("feffPionVsCentralityBin","Pion efficiency vs centrality bin","centrality bin","pion efficiency",fNCentBins,-1.5,maxCentbinRange,200,0,1.2);//
    CreateHisto2D("feffHadronVsCentralityBin","Hadron efficiency vs centrality bin","centrality bin","hadron efficiency",fNCentBins,-1.5,maxCentbinRange,200,0,1.2);//
    CreateHisto2D("feffKaonVsCentralityBin","Kaon efficiency vs centrality bin","centrality bin","kaon efficiency",fNCentBins,-1.5,maxCentbinRange,200,0,1.2);//
    CreateHisto2D("feffProtonVsCentralityBin","Proton efficiency vs centrality bin","centrality bin","proton efficiency",fNCentBins,-1.5,maxCentbinRange,200,0,1.2);//
    CreateHisto2D("fnotIDVsCentralityBin","f_{notID} vs centrality bin","centrality bin","f_{notID}",fNCentBins,-1.5,maxCentbinRange,50,0.95,1.05);//
    CreateHisto2D("fpTcutVsCentralityBin","f_{pTcut} vs centrality bin","centrality bin","f_{pTcut}",fNCentBins,-1.5,maxCentbinRange,50,0.95,1.05);
    CreateHisto2D("fneutralVsCentralityBin","f_{neutral} vs centrality bin","centrality bin","f_{neutral}",fNCentBins,-1.5,maxCentbinRange,50,0.5,1.00);
    CreateHisto2D("ConstantCorrectionsVsCentralityBin","constant corrections vs centrality bin","centrality bin","constant corrections",fNCentBins,-1.5,maxCentbinRange,50,0.5,1.00);
    CreateHisto2D("rawspectra","raw spectra vs centrality bin","centrality bin","p_{T}",fNCentBins,-1.5,maxCentbinRange,200,0,20);//
  }

  delete sTPC;
  delete sITS;
  delete sTPCpt;
  delete sITSpt;
  delete sPID;
  delete sNoPID;
  delete sNoPIDString;
  delete sPiKPEt;
  delete sHadEt;
  delete sTotEt;
  delete sTotEtString;
  delete sHadEtString;
  delete sPiKPEtString;
  delete sFull;
  delete sEMCAL;
  delete sPHOS;

}
