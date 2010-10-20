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

Int_t AliAnalysisHadEtReconstructed::AnalyseEvent(AliVEvent* ev)
{ // analyse ESD event
    ResetEventValues();

    AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev);
    //for PID
    AliESDpid *pID = new AliESDpid();
    pID->MakePID(realEvent);

    TString *strTPC = new TString("TPC");
    TString *strITS = new TString("ITS");
    TString *strTPCITS = new TString("TPCITS");
    for(Int_t cutset=0;cutset<2;cutset++){
      bool isTPC = false;
      TString *cutName;
      TObjArray* list;
      switch(cutset){
      case 0:
	cutName = strTPC;
	list = fEsdtrackCutsTPC->GetAcceptedTracks(realEvent);
	isTPC = true;
	break;
      case 1:
	cutName = strITS;
	list = fEsdtrackCutsITS->GetAcceptedTracks(realEvent);
	break;
      case 2:
	cutName = strTPCITS;
	list = fEsdtrackCutsITSTPC->GetAcceptedTracks(realEvent);
	break;
      default:
	cerr<<"Error:  cannot fill histograms!"<<endl;
	return -1;
      }
      Int_t nGoodTracks = list->GetEntries();
      //cout<<nGoodTracks<<" "<<cutName->Data()<<" tracks"<<endl;
      for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++)
	{


	  AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
	  if (!track)
	    {
	      Printf("ERROR: Could not get track %d", iTrack);
	      continue;
	    }
	  else{
	    //if(!isTPC) cout<<"I should add an ITS track..."<<endl;
	    if(TMath::Abs(track->Eta())>fCorrections->GetEtaCut()) continue;
	    Float_t nSigmaPion,nSigmaProton,nSigmaKaon,nSigmaElectron;
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
	    bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	    bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	    bool isKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
	    bool isProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);

	    //bool IsElectron = false;
	    bool unidentified = (!isProton && !isKaon && !isElectron && !isPion);
	    Float_t dEdx = track->GetTPCsignal();
	    if(cutset==1) dEdx = track->GetITSsignal();
	    FillHisto2D(Form("dEdxDataAll%s",cutName->Data()),track->P(),dEdx,1.0);

	    bool inPHOS = IsInPHOS(track);
	    bool inEMCAL = IsInEMCAL(track);
	    //if(!(fCorrections->GetEfficiencyPionTPC())) cerr<<"Uh-oh!  No histogram!"<<endl;

	    Float_t corrBkgd=0.0;
	    Float_t corrNotID=0.0;
	    Float_t corrNoID=0.0;// = fCorrections->GetNotIDCorrectionNoPID(track->Pt());
	    Float_t corrEff = 0.0;
	    Float_t corrEffNoID = 0.0;
	    if(cutset==0){//TPC
	      corrBkgd = fCorrections->GetBackgroundCorrectionTPC(track->Pt());
	      corrEffNoID = fCorrections->GetTPCEfficiencyCorrectionHadron(track->Pt());
	      corrNotID = fCorrections->GetNotIDConstCorrectionTPC();
	      corrNoID = fCorrections->GetNotIDConstCorrectionTPCNoID();
	    }
	    if(cutset==1){//ITS
	      corrBkgd = fCorrections->GetBackgroundCorrectionITS(track->Pt());
	      corrEffNoID = fCorrections->GetITSEfficiencyCorrectionHadron(track->Pt());
	      corrNotID = fCorrections->GetNotIDConstCorrectionITS();
	      corrNoID = fCorrections->GetNotIDConstCorrectionITSNoID();
	    }
	    Float_t et = 0.0;
	    Float_t etNoID = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	    Float_t etpartialcorrected = 0.0;
	    Float_t etpartialcorrectedNoID = corrNoID*corrBkgd*corrEffNoID*etNoID;
	    FillHisto2D(Form("EtDataRaw%sNoID",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrectedNoID);

	    if(isPion){
	      FillHisto2D(Form("dEdxDataPion%s",cutName->Data()),track->P(),dEdx,1.0);
	      et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	      if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionPion(track->Pt());}
	      //else{corrEff = fCorrections->GetITSEfficiencyCorrectionPion(track->Pt());}
	      etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	      
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
	      if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionKaon(track->Pt());}
	      //else{corrEff = fCorrections->GetITSEfficiencyCorrectionKaon(track->Pt());}
	      etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	      
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
	      if(cutset==0){corrEff = fCorrections->GetTPCEfficiencyCorrectionProton(track->Pt());}
	      //else{corrEff = fCorrections->GetITSEfficiencyCorrectionProton(track->Pt());}
	      etpartialcorrected = et*corrBkgd*corrEff*corrNotID;
	      
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
	      FillHisto2D(Form("dEdxDataProton%s",cutName->Data()),track->P(),dEdx,1.0);
	      //et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	    }
	    if(unidentified){
	      //if(!isPion) 
	      FillHisto2D(Form("dEdxDataUnidentified%s",cutName->Data()),track->P(),dEdx,1.0);
	      et = Et(track->P(),track->Theta(),fgPiPlusCode,track->Charge());
	      etpartialcorrected = et*corrBkgd*corrEffNoID*corrNotID;
	      //if(!isPion) 
	      FillHisto2D(Form("EtDataCorrected%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),etpartialcorrected);
	    }
	    if(!isTPC){
	      etpartialcorrected = etpartialcorrectedNoID;//Not using PID for ITS
	      //cout<<"I should add an ITS track..."<<endl;
	    }
	    AddEt(et,etNoID,etpartialcorrected,etpartialcorrectedNoID,track->Pt(),isTPC,inPHOS,inEMCAL);
	    // cout<<"Pt "<<track->Pt()<<" et corrected no id "<<etpartialcorrectedNoID<<"="<<corrNoID<<"*"<<corrBkgd<<"*"<<corrEffNoID<<"*"<<etNoID<<endl;
// 	    if(isTPC){
// 	      //if(!unidentified){
// 	      //cout<<"Pt "<<track->Pt()<<" et corrected   pid "<<etpartialcorrected<<"="<<corrBkgd<<"*"<<corrEff<<"*"<<et;
// 	      cout<<"Pt "<<track->Pt()<<" et corrected  pid "<<etpartialcorrected<<"="<<corrNotID<<"*"<<corrBkgd<<"*"<<corrEff<<"*"<<et;
// 	      if(isPion)cout<<" pion";
// 	      if(isProton)cout<<" proton";
// 	      if(isKaon)cout<<" kaon";
// 	      if(isElectron)cout<<" electron";
// 	      if(unidentified)cout<<" unidentified";
// 	      cout<<endl;
// 		//}
// 	      //else{cout<<"Pt "<<track->Pt()<<" et corrected  !pid "<<etpartialcorrected<<"="<<corrBkgd<<"*"<<corrEff<<"*"<<corrNotID<<"*"<<et<<" Not identified"<<endl;}
// 	    }
	    //if(inEMCAL) cout<<"I should add a track"<<endl;
	  }
	}
    }
//   cout<<"Finishing with Raw/Corrected Et in full, PHOS, EMCAL acceptance of "
//       << GetRawEtFullAcceptanceITS() <<"/"
//       << GetCorrectedHadEtFullAcceptanceITS() <<", "
//       << GetRawEtPHOSAcceptanceITS() <<"/"
//       << GetCorrectedHadEtPHOSAcceptanceITS() <<", "
//       << GetRawEtEMCALAcceptanceITS() <<"/"
//       << GetCorrectedHadEtEMCALAcceptanceITS() <<endl;
//   cout<<"Finishing with Raw/Corrected Et w/o PID in full, PHOS, EMCAL acceptance of "
//       << GetRawEtFullAcceptanceITSNoPID() <<"/"
//       << GetCorrectedHadEtFullAcceptanceITSNoPID() <<", "
//       << GetRawEtPHOSAcceptanceITSNoPID() <<"/"
//       << GetCorrectedHadEtPHOSAcceptanceITSNoPID() <<", "
//       << GetRawEtEMCALAcceptanceITSNoPID() <<"/"
//       << GetCorrectedHadEtEMCALAcceptanceITSNoPID() <<endl;
//   cout<<"Finishing with Raw/Corrected Et in full, PHOS, EMCAL acceptance of "
//       << GetRawEtFullAcceptanceTPC() <<"/"
//       << GetCorrectedHadEtFullAcceptanceTPC() <<", "
//       << GetRawEtPHOSAcceptanceTPC() <<"/"
//       << GetCorrectedHadEtPHOSAcceptanceTPC() <<", "
//       << GetRawEtEMCALAcceptanceTPC() <<"/"
//       << GetCorrectedHadEtEMCALAcceptanceTPC() <<endl;
//   cout<<"Finishing with Raw/Corrected Et w/o PID in full, PHOS, EMCAL acceptance of "
//       << GetRawEtFullAcceptanceTPCNoPID() <<"/"
//       << GetCorrectedHadEtFullAcceptanceTPCNoPID() <<", "
//       << GetRawEtPHOSAcceptanceTPCNoPID() <<"/"
//       << GetCorrectedHadEtPHOSAcceptanceTPCNoPID() <<", "
//       << GetRawEtEMCALAcceptanceTPCNoPID() <<"/"
//       << GetCorrectedHadEtEMCALAcceptanceTPCNoPID() <<endl;
//   cout<<"Correction factors "
//       <<fCorrTotEtFullAcceptanceTPC<<", "<<fCorrTotEtFullAcceptanceITS<<", "<<fCorrHadEtFullAcceptanceTPC<<", "<<fCorrHadEtFullAcceptanceITS<<","
//       <<fCorrTotEtEMCALAcceptanceTPC<<", "<<fCorrTotEtEMCALAcceptanceITS<<", "<<fCorrHadEtEMCALAcceptanceTPC<<", "<<fCorrHadEtEMCALAcceptanceITS<<","
//       <<fCorrTotEtPHOSAcceptanceTPC<<", "<<fCorrTotEtPHOSAcceptanceITS<<", "<<fCorrHadEtPHOSAcceptanceTPC<<", "<<fCorrHadEtPHOSAcceptanceITS<<endl;
    FillHisto1D("RecoHadEtFullAcceptanceTPC",GetCorrectedHadEtFullAcceptanceTPC(),1.0);
    FillHisto1D("RecoTotEtFullAcceptanceTPC",GetCorrectedTotEtFullAcceptanceTPC(),1.0);
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
    delete pID;
    delete strTPC;
    delete strITS;
    delete strTPCITS;
    return 1;
}
void AliAnalysisHadEtReconstructed::AddEt(Float_t rawEt, Float_t rawEtNoPID, Float_t corrEt, Float_t corrEtNoPID, Float_t pt, Bool_t IsTPC, Bool_t InPHOS, Bool_t InEMCAL) {//Adding Et to each of the variables that tracks et event by event
  if(pt>=AliAnalysisHadEt::fgPtTPCCutOff && IsTPC){//TPC tracks
    //adding to the raw Et
    //if(InEMCAL)
    //cout<<"Adding "<<rawEt<<" to the raw TPC Et"<<endl;
    fRawEtFullAcceptanceTPC += rawEt;
    if(InPHOS)fRawEtPHOSAcceptanceTPC += rawEt;
    if(InEMCAL)fRawEtEMCALAcceptanceTPC += rawEt;
    fRawEtFullAcceptanceTPCNoPID += rawEtNoPID;
    if(InPHOS)fRawEtPHOSAcceptanceTPCNoPID += rawEtNoPID;
    if(InEMCAL)fRawEtEMCALAcceptanceTPCNoPID += rawEtNoPID;
    //adding to the corrected Et
    //if(InPHOS) cout<<"Adding "<<corrEt<<" to the corrected Et"<<endl;
    fCorrectedHadEtFullAcceptanceTPC += corrEt;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceTPC += corrEt;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceTPC += corrEt;
    fCorrectedHadEtFullAcceptanceTPCNoPID += corrEtNoPID;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceTPCNoPID += corrEtNoPID;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceTPCNoPID += corrEtNoPID;
    //cout<<"Adding "<<corrEt<<" to the corrected Et, TPC pT "<<pt<<endl;
  }
  if(pt<AliAnalysisHadEt::fgPtTPCCutOff &&pt>=AliAnalysisHadEt::fgPtITSCutOff && !IsTPC){//ITS tracks
    //adding to the raw Et
    fRawEtFullAcceptanceITS += rawEt;
    if(InPHOS)fRawEtPHOSAcceptanceITS += rawEt;
    if(InEMCAL)fRawEtEMCALAcceptanceITS += rawEt;
    fRawEtFullAcceptanceITSNoPID += rawEtNoPID;
    if(InPHOS)fRawEtPHOSAcceptanceITSNoPID += rawEtNoPID;
    if(InEMCAL)fRawEtEMCALAcceptanceITSNoPID += rawEtNoPID;
    //cout<<"Adding "<<rawEt<<" to the raw ITS Et"<<endl;
    //adding to the corrected Et
    fCorrectedHadEtFullAcceptanceITS += corrEt;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceITS += corrEt;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceITS += corrEt;
    fCorrectedHadEtFullAcceptanceITSNoPID += corrEtNoPID;
    if(InPHOS)fCorrectedHadEtPHOSAcceptanceITSNoPID += corrEtNoPID;
    if(InEMCAL)fCorrectedHadEtEMCALAcceptanceITSNoPID += corrEtNoPID;
    //cout<<"Adding "<<corrEt<<" to the corrected Et, ITS pT "<<pt<<endl;
  }
}

Bool_t AliAnalysisHadEtReconstructed::IsInPHOS(AliESDtrack *track){//This function will need to be elaborated on later to include PHOS dead channels
  return   TMath::Abs(track->Eta()) < fCuts->GetGeometryPhosEtaAccCut()//in eta acceptance
    && track->Phi()*180.0/TMath::Pi() > fCuts->GetGeometryPhosPhiAccMinCut()//greater than the minimum phi
    && track->Phi()*180.0/TMath::Pi() < fCuts->GetGeometryPhosPhiAccMaxCut();//less than the maximum phi
}
Bool_t AliAnalysisHadEtReconstructed::IsInEMCAL(AliESDtrack *track){//This function will need to be elaborated on later to include EMCAL dead channels
  //cout<<"Eta: |"<<track->Eta()<<"|<"<< fCuts->GetGeometryEmcalEtaAccCut() <<"; phi: "<<fCuts->GetGeometryEmcalPhiAccMinCut()<<"<"<<track->Phi()*180.0/TMath::Pi()<<"<"<<fCuts->GetGeometryEmcalPhiAccMaxCut()<<endl;
  return   TMath::Abs(track->Eta()) < fCuts->GetGeometryEmcalEtaAccCut()//in eta acceptance
    && track->Phi()*180.0/TMath::Pi() > fCuts->GetGeometryEmcalPhiAccMinCut()//greater than the minimum phi
    && track->Phi()*180.0/TMath::Pi() < fCuts->GetGeometryEmcalPhiAccMaxCut();//less than the maximum phi
}
Bool_t AliAnalysisHadEtReconstructed::CheckGoodVertex(AliVParticle* track)
{ // check vertex

    Float_t bxy = 999.;
    Float_t bz = 999.;
    dynamic_cast<AliESDtrack*>(track)->GetImpactParametersTPC(bxy,bz);

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
  if (fConfigFile.Length() && !fCorrections) {
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
	 cout<<"Rereading fCorrections file..."<<endl;
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

  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  for(Int_t i=0;i<2;i++){
    TString *cutName;
    Float_t maxPtdEdx = 10;
    Float_t mindEdx = 35;
    Float_t maxdEdx = 150.0;
    switch(i){
    case 0:
      cutName = strTPC;
      break;
    case 1:
      cutName = strITS;
      maxPtdEdx = 5;
      maxdEdx = 500.0;
      break;
    case 2:
      cutName = strTPCITS;
      break;
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
	  TString *detector;
	  TString *partid;
	  TString *et = sHadEt;
	  TString *acceptance;
	  TString *ptstring;
	  TString *partidstring;
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
	  sprintf(histoname,"Reco%s%sAcceptance%s%s",et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	  sprintf(histotitle,"Reconstructed %s with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	  sprintf(xtitle,"Reconstructed %s",etstring->Data());
	  CreateHisto1D(histoname,histotitle,xtitle,ytitle->Data(),nbinsEt*2,minEt,maxEt);
	}
      }
    }
  }

  //CreateHisto2D("Efficiency","Efficiency","pT","efficiency",

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
