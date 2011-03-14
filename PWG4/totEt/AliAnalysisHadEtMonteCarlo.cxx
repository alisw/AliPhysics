//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for MC analysis
//  - MC output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#include "AliAnalysisHadEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliVParticle.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include <iostream>
#include "TRandom.h"
#include "AliAnalysisEtCommon.h"
#include "AliCentrality.h"
#include "AliLog.h"

using namespace std;

ClassImp(AliAnalysisHadEtMonteCarlo);


Int_t AliAnalysisHadEtMonteCarlo::fgNumSmearWidths = 4;
Float_t AliAnalysisHadEtMonteCarlo::fgSmearWidths[4] = {0.005,0.006,0.007,0.008};

AliAnalysisHadEtMonteCarlo::AliAnalysisHadEtMonteCarlo():AliAnalysisHadEt()
							,fSimPiKPEt(0)
							,fSimHadEt(0)
							,fSimTotEt(0) 
							,fInvestigateSmearing(0)
							,fInvestigateFull(0)
							,fInvestigateEMCal(0)
							,fInvestigatePHOS(0)
							,fInvestigatePiKP(0)
							,fRequireITSHits(0)
							,fBaryonEnhancement(0)
							,fCentBin(-1)
							,fPtSmearer(0)
{
}
AliAnalysisHadEtMonteCarlo::~AliAnalysisHadEtMonteCarlo(){//destructor
  if(fPtSmearer) delete fPtSmearer;
}

void AliAnalysisHadEtMonteCarlo::ResetEventValues(){//resetting event variables
  AliAnalysisHadEt::ResetEventValues();
    fSimHadEt=0.0;
    fSimTotEt=0.0;
    fSimPiKPEt=0.0;
}
Int_t AliAnalysisHadEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{ // analyse MC and real event info
  FillHisto1D("NEvents",0.5,1);
  if(!ev || !ev2){
    AliFatal("ERROR: Event does not exist");   
    return 0;
  }
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
  if(!mcEvent || !realEvent){  
    AliFatal("ERROR: mcEvent or realEvent does not exist");
    return 0;
  }
  AliStack *stack = mcEvent->Stack();
  fCentBin= -1;
  if(fDataSet==20100){//If this is Pb+Pb
    AliCentrality *centrality = realEvent->GetCentrality();
    if(fNCentBins<20) fCentBin= centrality->GetCentralityClass10(fCentralityMethod);
    else{ fCentBin= centrality->GetCentralityClass5(fCentralityMethod);}
  }
  AnalyseEvent(ev);

  //for PID
  AliESDpid *pID = new AliESDpid();//This is identified as a memory leak in valgrind but I delete this object so I think it may be a problem with AliESDpid.

  //=============================================

  //Roughly following $ALICE_ROOT/PWG0/dNdEta/AlidNdEtaCorrectionTask

  //=============================================TPC&&ITS=============================================
  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  Int_t lastcutset = 1;
  if(fRequireITSHits) lastcutset = 2;
  for(Int_t cutset=0;cutset<=lastcutset;cutset++){
    TString *cutName = NULL;
    TObjArray* list = NULL;
    switch(cutset){
    case 0:
      cutName = strTPC;
      list = fEsdtrackCutsTPC->GetAcceptedTracks(realEvent);
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
    for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++)
      {
	AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
	if (!track)
	  {
	    Printf("ERROR: Could not get track %d", iTrack);
	    continue;
	  }
	else{
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
// 	  bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
// 	  bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
// 	  bool isKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
// 	  bool isProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);
	  bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  bool isKaon = (nSigmaPion>3.0 && nSigmaProton>3.0 && nSigmaKaon<3.0 && track->Pt()<0.45);
	  bool isProton = (nSigmaPion>3.0 && nSigmaProton<3.0 && nSigmaKaon>3.0 && track->Pt()<0.9);

	  bool unidentified = (!isProton && !isKaon && !isElectron);
	  Float_t dEdx = track->GetTPCsignal();
	  if(cutset==1) dEdx = track->GetITSsignal();

	  FillHisto2D(Form("dEdxAll%s",cutName->Data()),track->P(),dEdx,1.0);

	  UInt_t label = (UInt_t)TMath::Abs(track->GetLabel());
	  TParticle  *simPart  = stack->Particle(label);
	  if(!simPart) {
	    Printf("no MC particle\n"); 	 	
	    continue; 	 	
	  }
	  else{//analysis
	    if(stack->IsPhysicalPrimary(label)){
	      if (TMath::Abs(simPart->Eta()) < fCuts->GetCommonEtaCut()){
		Int_t pdgCode =  simPart->GetPDG(0)->PdgCode();
		Int_t mypid = 0;
		if(pdgCode==AliAnalysisHadEt::fgPiPlusCode) mypid = 1;
		if(pdgCode==fgProtonCode) mypid = 2;
		if(pdgCode==fgKPlusCode) mypid = 3;
		if(pdgCode==fgEPlusCode) mypid = 4;
		if(pdgCode==fgPiMinusCode) mypid = 1;
		if(pdgCode==fgAntiProtonCode) mypid = 2;
		if(pdgCode==fgKMinusCode) mypid = 3;
		if(pdgCode==fgEMinusCode) mypid = 4;
		bool filled = false;      
		//============Charged hadrons===================================
		//identified...
		if(isPion){
		  if(pdgCode!=fgPiPlusCode && pdgCode!=fgPiMinusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),1,mypid,1);
		    //if(mypid==0)cerr<<"I was misidentified! I'm not a pion! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxPion%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isProton){
		  if(pdgCode!=fgProtonCode && pdgCode!=fgAntiProtonCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),2,mypid,1);
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedProton",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  if(fBaryonEnhancement){
		    myEt = myEt*ProtonBaryonEnhancement(track->Pt());
		    if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedProtonEnhanced",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		    else{ FillHisto2D(Form("EtReconstructed%sIdentifiedAntiProtonEnhanced",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  }
		  FillHisto2D(Form("dEdxProton%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isKaon){
		  if(pdgCode!=fgKMinusCode && pdgCode!=fgKPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),3,mypid,1);
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedKPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedKMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxKaon%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isElectron){
		  if(pdgCode!=fgEMinusCode && pdgCode!=fgEPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),4,mypid,1);
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedEPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedEMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxElectron%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(unidentified){
		  if(pdgCode!=fgEMinusCode && pdgCode!=fgEPlusCode){
		    float myEtPi = Et(simPart,fgPionMass);
		    float myEtP = Et(simPart,fgProtonMass);
		    float myEtK = Et(simPart,fgKaonMass);
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		    FillHisto2D(Form("EtReconstructed%sUnidentifiedAssumingProton",cutName->Data()),track->Pt(),track->Eta(),myEtP);
		    FillHisto2D(Form("EtReconstructed%sUnidentifiedAssumingKaon",cutName->Data()),track->Pt(),track->Eta(),myEtK);
		    FillHisto2D(Form("EtReconstructed%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),1.0);
		    if(pdgCode == fgPiPlusCode||pdgCode == fgPiMinusCode){
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingPion",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingProton",cutName->Data()),track->Pt(),track->Eta(),myEtP);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingKaon",cutName->Data()),track->Pt(),track->Eta(),myEtK);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedPion",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      FillHisto2D(Form("EtNReconstructed%sUnidentifiedPion",cutName->Data()),track->Pt(),track->Eta(),1.0);
		    }
		    if(pdgCode == fgKPlusCode||pdgCode == fgKMinusCode){
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingPion",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingProton",cutName->Data()),track->Pt(),track->Eta(),myEtP);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingKaon",cutName->Data()),track->Pt(),track->Eta(),myEtK);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedKaon",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      FillHisto2D(Form("EtNReconstructed%sUnidentifiedKaon",cutName->Data()),track->Pt(),track->Eta(),1.0);
		    }
		    if(pdgCode == fgProtonCode||pdgCode == fgAntiProtonCode){
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingPion",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingProton",cutName->Data()),track->Pt(),track->Eta(),myEtP);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingKaon",cutName->Data()),track->Pt(),track->Eta(),myEtK);
		      FillHisto2D(Form("EtReconstructed%sUnidentifiedProton",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      FillHisto2D(Form("EtNReconstructed%sUnidentifiedProton",cutName->Data()),track->Pt(),track->Eta(),1.0);
		      if(fBaryonEnhancement){
			myEt = myEt*ProtonBaryonEnhancement(track->Pt());
			FillHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingPionEnhanced",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
			FillHisto2D(Form("EtReconstructed%sUnidentifiedProtonEnhanced",cutName->Data()),track->Pt(),track->Eta(),myEt);
			FillHisto2D(Form("EtNReconstructed%sUnidentifiedProtonEnhanced",cutName->Data()),track->Pt(),track->Eta(),1.0);
		      }
		    }
		  }
		  FillHisto2D(Form("dEdxUnidentified%s",cutName->Data()),track->P(),dEdx,1.0);
		  FillHisto1D(Form("UnidentifiedPIDs%s",cutName->Data()),mypid,1);
		}
		//...simulated
		if(pdgCode == fgPiPlusCode){		
		  float myEt = Et(simPart);
		  float myEtP = Et(simPart,fgProtonMass);
		  float myEtK = Et(simPart,fgKaonMass);
		  FillHisto2D(Form("EtReconstructed%sPiPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sPiPlusCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sPiPlusAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sPiPlusAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  filled = true;
		}
		if(pdgCode == fgPiMinusCode){
		  float myEt = Et(simPart);
		  float myEtP = Et(simPart,fgProtonMass);
		  float myEtK = Et(simPart,fgKaonMass);
		  FillHisto2D(Form("EtReconstructed%sPiMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sPiMinusCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sPiMinusAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sPiMinusAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  filled = true;
		}
		if(pdgCode == fgKPlusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fgPionMass);
		  float myEtP = Et(simPart,fgProtonMass);
		  FillHisto2D(Form("EtReconstructed%sKPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sKPlusCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKPlusAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sKPlusAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  FillHisto2D(Form("EtReconstructed%sKPlusAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  filled = true;
		}
		if(pdgCode == fgKMinusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fgPionMass);
		  float myEtP = Et(simPart,fgProtonMass);
		  FillHisto2D(Form("EtReconstructed%sKMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sKMinusCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKMinusAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sKMinusAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  FillHisto2D(Form("EtReconstructed%sKMinusAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtP);
		  filled = true;
		}
		if(pdgCode == fgProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fgPionMass);
		  float myEtK = Et(simPart,fgKaonMass);
		  FillHisto2D(Form("EtReconstructed%sProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sProtonCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sProtonAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sProtonAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sProtonAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;

		      if(fBaryonEnhancement){
			float enhancement = ProtonBaryonEnhancement(track->Pt());
			FillHisto2D(Form("EtReconstructed%sProtonEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt*enhancement);
			FillHisto2D(Form("EtNReconstructed%sProtonEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt*enhancement);
			FillHisto2D(Form("EtReconstructed%sProtonAssumingPionEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi*enhancement);
		      }

		}
		if(pdgCode == fgAntiProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fgPionMass);
		  float myEtK = Et(simPart,fgKaonMass);
		  FillHisto2D(Form("EtReconstructed%sAntiProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sAntiProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(fCentBin>=0){//if a centrality bin was defined
		    FillHisto2D(Form("EtNReconstructed%sAntiProtonCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),fCentBin),simPart->Pt(),simPart->Eta(),myEt);
		  }
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingKaon",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtK);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;
		  if(fBaryonEnhancement){
			float enhancement = ProtonBaryonEnhancement(track->Pt());
			FillHisto2D(Form("EtReconstructed%sAntiProtonEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt*enhancement);
		    FillHisto2D(Form("EtNReconstructed%sAntiProtonEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt*enhancement);
		    FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingPionEnhanced",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi*enhancement);
		  }
		}
		if(pdgCode == fgEPlusCode){
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sEPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(!isElectron || unidentified){
		    float myEtPi = Et(simPart,fgPionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  filled = true;
		}
		if(pdgCode == fgEMinusCode){
		  if(!isElectron || unidentified){
		    float myEtPi = Et(simPart,fgPionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sEMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		    filled = true;
		}
	      }
	      
	    }
	    else{//not a primary - we're after V0 daughters!
	      if (TMath::Abs(simPart->Eta()) < fCuts->GetCommonEtaCut()){
		TParticle *mom = stack->Particle(simPart->GetFirstMother());
		if(mom){
		  TParticlePDG *pc = mom->GetPDG(0);
		  if(pc){
		    Int_t pdgCode =  mom->GetPDG(0)->PdgCode();
		    if(pdgCode == fgLambdaCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sLambdaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      Float_t weight = LambdaWeight(mom->Pt());
		      if(fBaryonEnhancement){
			float enhancement = ProtonBaryonEnhancement(track->Pt());
			weight = weight*enhancement;
		      }
		      FillHisto2D(Form("EtReconstructed%sLambdaDaughtersReweighted",cutName->Data()),track->Pt(),track->Eta(),myEt*weight);
		    }
		    if(pdgCode == fgAntiLambdaCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      Float_t weight = AntiLambdaWeight(mom->Pt());
		      if(fBaryonEnhancement){
			float enhancement = ProtonBaryonEnhancement(track->Pt());
			weight = weight*enhancement;
		      }
		      FillHisto2D(Form("EtReconstructed%sAntiLambdaDaughtersReweighted",cutName->Data()),track->Pt(),track->Eta(),myEt*weight);
		    }
		    if(pdgCode == fgK0SCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sK0SDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		      Float_t weight = K0Weight(mom->Pt());
		      FillHisto2D(Form("EtReconstructed%sK0SDaughtersReweighted",cutName->Data()),track->Pt(),track->Eta(),myEt*weight);
		    }
		    if(pdgCode == fgXiCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    }
		    if(pdgCode == fgAntiXiCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    }
		    if(pdgCode == fgOmegaCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    }
		    if(pdgCode == fgXiCode){
		      float myEt = Et(simPart);
		      FillHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    }

		    if(mom->GetFirstMother()>0){
		      TParticle *grandma = stack->Particle(mom->GetFirstMother());
		      if(grandma){
			Int_t pdgCodeMom =  mom->GetPDG(0)->PdgCode();
			if(pdgCodeMom==fgPiPlusCode || pdgCodeMom==fgPiMinusCode || pdgCodeMom==fgProtonCode ||pdgCodeMom==fgAntiProtonCode || pdgCodeMom==fgKPlusCode || pdgCode==fgKMinusCode){
			  Int_t pdgCodeGrandma =  grandma->GetPDG(0)->PdgCode();
		      
			  if(pdgCodeGrandma == fgXiCode){
			    float myEt = Et(simPart);
			    FillHisto2D(Form("EtReconstructed%sXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			  }
			  if(pdgCodeGrandma == fgAntiXiCode){
			    float myEt = Et(simPart);
			    FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			  }
			  if(pdgCodeGrandma == fgOmegaCode){
			    float myEt = Et(simPart);
			    FillHisto2D(Form("EtReconstructed%sOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			  }
			  if(pdgCodeGrandma == fgXiCode){
			    float myEt = Et(simPart);
			    FillHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			  }

			}
		      }
		    }
		  }
		}
	      }
	    }
	  }

	}
      }
    delete list;
  }
  delete pID;
  delete strTPC;
  delete strITS;
  delete strTPCITS;
  //delete AliESDpid;
  return 1;
}
Int_t AliAnalysisHadEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
     ResetEventValues();
     if(!ev){
            AliFatal("ERROR: Event does not exist");   
	    return 0;
     }
     
    // Get us an mc event
    AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
    if(!mcEvent){  
      AliFatal("ERROR: MC Event does not exist");
      return 0;
    }
    // Let's play with the stack!
    AliStack *stack = mcEvent->Stack();

    Int_t nPrim = stack->GetNtrack();

    Float_t fSimPiKPEtPtSmeared = 0;
    Float_t fSimPiKPEtEfficiencySmeared = 0;
    Float_t fSimPiKPEtPtCutSmearedTPC = 0;
    Float_t fSimPiKPEtPtCutSmearedITS = 0;
    Float_t fSimPiKPEtPIDSmeared = 0;
    Float_t fSimPiKPEtPIDSmearedNoID = 0;
    //=================Tracks which may or may not have been reconstructed=================

    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

      TParticle *part = stack->Particle(iPart);//This line is identified as a loss of memory by valgrind, however, the pointer still belongs to the stack, so it's the stack's problem

        if (!part)
	  {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
	  }

        //TParticlePDG *pc = part->GetPDG(0);

        // Check if it is a primary particle
	if (stack->IsPhysicalPrimary(iPart)){//primaries

	  if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())	    {

	    Int_t pdgCode =  part->GetPDG(0)->PdgCode();
	    bool filled = false;
	    //Investigating smearing...
	    //Numbers are realistic correction factors from previous studies
	    if(fInvestigateSmearing){
	      if(pdgCode==fgPiPlusCode ||pdgCode==fgPiMinusCode ||pdgCode==fgKPlusCode ||pdgCode==fgKMinusCode ||pdgCode==fgProtonCode ||pdgCode==fgAntiProtonCode){
		//To investigate Smearing...
		Float_t myet = Et(part);
		fSimPiKPEt += myet;
		Float_t theta = part->Theta();
		Short_t charge = 1;
		Float_t momentum = part->P();
		//pt smearing
		Float_t pSmeared = momentum *  fPtSmearer->Gaus(1,0.005);//Gaussian centered around 1
		fSimPiKPEtPtSmeared += Et(pSmeared,theta,pdgCode,charge);
		//Efficiency smearing
		float efficiency = 2.26545*TMath::Exp(-TMath::Power(9.99977e-01/part->Pt(),7.85488e-02));//simple rough efficiency from fitting curve
		if(fPtSmearer->Binomial(1,efficiency) ==1){
		  fSimPiKPEtEfficiencySmeared += (1.0/efficiency)*myet;
		}
		//pT cut smeared
		if(part->Pt()>0.10){fSimPiKPEtPtCutSmearedITS +=1.00645645*myet;}
		if(part->Pt()>0.15){fSimPiKPEtPtCutSmearedTPC +=1.02000723*myet;}
		//PID smearing
		fSimPiKPEtPIDSmearedNoID += 1.02679314*Et(momentum,theta,fgPiPlusCode,charge);
		if(part->P()<1.0){//then the particle would have been ID'd
		  fSimPiKPEtPIDSmeared += 1.0085942*myet;
		}
		else{//Then it would have been assumed to be a pion
		  fSimPiKPEtPIDSmeared += 1.0085942*Et(momentum,theta,fgPiPlusCode,charge);
		}
	      }
	    }

	    //============Charged hadrons===================================
	    if(pdgCode == fgPiPlusCode){
	      float myEt = Et(part);
	      float myEtP = Et(part,fgProtonMass);
	      float myEtK = Et(part,fgKaonMass);

	      fSimHadEt += myEt;
	      fSimTotEt += myEt;

	      FillHisto2D("EtSimulatedPiPlus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedPiPlus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedPiPlusCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedChargedHadronAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedPiPlusAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedChargedHadronAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedPiPlusAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = 1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	    }
	    if(pdgCode == fgPiMinusCode){
	      float myEt = Et(part);
	      float myEtP = Et(part,fgProtonMass);
	      float myEtK = Et(part,fgKaonMass);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedPiMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedPiMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedPiMinusCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedChargedHadronAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedPiMinusAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedChargedHadronAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedPiMinusAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = -1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	    }
	    if(pdgCode == fgKPlusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fgPionMass);
	      float myEtP = Et(part,fgProtonMass);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedKPlus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKPlus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedKPlusCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKPlusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedChargedHadronAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedKPlusAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = 1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	    }
	    if(pdgCode == fgKMinusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fgPionMass);
	      float myEtP = Et(part,fgProtonMass);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedKMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedKMinusCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKMinusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedChargedHadronAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedKMinusAssumingProton",part->Pt(),part->Eta(),myEtP);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = -1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	    }
	    if(pdgCode == fgProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fgPionMass);
	      float myEtK = Et(part,fgKaonMass);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedProton",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedProton",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedProtonCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedProtonAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedChargedHadronAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedProtonAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = 1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	      if(fBaryonEnhancement){
		float enhancement = ProtonBaryonEnhancement(part->Pt());
		FillHisto2D("EtSimulatedProtonEnhanced",part->Pt(),part->Eta(),myEt*enhancement);
		FillHisto2D("EtNSimulatedProtonEnhanced",part->Pt(),part->Eta(),1.0*enhancement);
		FillHisto2D("EtSimulatedProtonAssumingPionEnhanced",part->Pt(),part->Eta(),myEtPi*enhancement);
	      }
	    }
	    if(pdgCode == fgAntiProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fgPionMass);
	      float myEtK = Et(part,fgKaonMass);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiProton",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedAntiProton",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      if(fCentBin>=0){//if a centrality bin was defined
		FillHisto2D(Form("EtNSimulatedAntiProtonCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
		FillHisto2D(Form("EtNSimulatedChargedHadronCB%i",fCentBin),part->Pt(),part->Eta(),1.0);
	      }
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAntiProtonAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedChargedHadronAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedAntiProtonAssumingKaon",part->Pt(),part->Eta(),myEtK);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Short_t charge = -1;
	      Float_t myEtLow = Et(0.0,part->Theta(),pdgCode,charge);
	      Float_t myEtITS = Et(0.10/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      Float_t myEtTPC = Et(0.15/TMath::Sin(part->Theta()),part->Theta(),pdgCode,charge);
	      FillHisto2D("EtSimulatedChargedHadronAssumingNoPt",part->Pt(),part->Eta(),myEtLow);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut",part->Pt(),part->Eta(),myEtTPC);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPtITSCut",part->Pt(),part->Eta(),myEtITS);
	      filled = true;
	      if(fBaryonEnhancement){
		float enhancement = ProtonBaryonEnhancement(part->Pt());
		FillHisto2D("EtSimulatedAntiProtonEnhanced",part->Pt(),part->Eta(),myEt*enhancement);
		FillHisto2D("EtNSimulatedAntiProtonEnhanced",part->Pt(),part->Eta(),1.0*enhancement);
		FillHisto2D("EtSimulatedAntiProtonAssumingPionEnhanced",part->Pt(),part->Eta(),myEtPi*enhancement);
	      }
	    }
	    //============Other hadrons===================================

	    if(pdgCode == fgNeutronCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgAntiNeutronCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgLambdaCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      //cout<<"I am a simulated lambda! pt "<<part->Pt()<<" eta "<<part->Eta()<<endl;
	      FillHisto2D("EtSimulatedLambda",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Float_t weight = LambdaWeight(part->Pt());
	      if(fBaryonEnhancement){
		float enhancement = ProtonBaryonEnhancement(part->Pt());
		weight = weight*enhancement;
	      }
	      FillHisto2D("EtSimulatedLambdaReweighted",part->Pt(),part->Eta(),myEt*weight);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiMinusCode || daughtercode==fgProtonCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedLambdaDaughters",daughter->Pt(),daughter->Eta(),myEt);
		      FillHisto2D("EtSimulatedLambdaDaughtersReweighted",daughter->Pt(),daughter->Eta(),myEt*weight);
		      //cout<<"Lambda daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"Lambda daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgAntiLambdaCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiLambda",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Float_t weight = AntiLambdaWeight(part->Pt());
	      if(fBaryonEnhancement){
		float enhancement = ProtonBaryonEnhancement(part->Pt());
		weight = weight*enhancement;
	      }
	      FillHisto2D("EtSimulatedAntiLambdaReweighted",part->Pt(),part->Eta(),myEt*weight);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiPlusCode || daughtercode==fgAntiProtonCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedAntiLambdaDaughters",daughter->Pt(),daughter->Eta(),myEt);
		      FillHisto2D("EtSimulatedAntiLambdaDaughtersReweighted",daughter->Pt(),daughter->Eta(),myEt*weight);
		      //cout<<"AntiLambda daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"AntiLambda daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgK0SCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedK0S",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Float_t weight = K0Weight(part->Pt());
	      FillHisto2D("EtSimulatedK0SReweighted",part->Pt(),part->Eta(),myEt*weight);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){

		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiMinusCode || daughtercode==fgPiPlusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedK0SDaughters",daughter->Pt(),daughter->Eta(),myEt);
		      FillHisto2D("EtSimulatedK0SDaughtersReweighted",daughter->Pt(),daughter->Eta(),myEt*weight);
		      //cout<<"K0S daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"K0S daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgK0LCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedK0L",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Float_t weight = K0Weight(part->Pt());
	      FillHisto2D("EtSimulatedK0LReweighted",part->Pt(),part->Eta(),myEt*weight);
	      filled = true;
	    }
	    if(pdgCode == fgOmegaCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedOmega",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->Particle(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){

		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiPlusCode || daughtercode==fgProtonCode || daughtercode==fgKMinusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedOmegaDaughters",daughter->Pt(),daughter->Eta(),myEt);
		    //cout<<"Omega daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"Omega daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgAntiOmegaCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedOmega",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiMinusCode || daughtercode==fgAntiProtonCode || daughtercode==fgKPlusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedAntiOmegaDaughters",daughter->Pt(),daughter->Eta(),myEt);
		      //cout<<"AntiOmega daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"AntiOmega daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    //There are two codes for Sigmas
	    if(pdgCode == fgSigmaCode || pdgCode == -3222){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgAntiSigmaCode || pdgCode == 3222){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgXiCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedXi",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5 || daughterindex>1e5) continue;
		//cerr<<"Daughter index "<<daughterindex<<" npart "<<nPrim<<endl;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){

		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiPlusCode || daughtercode==fgProtonCode || daughtercode==fgPiMinusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedXiDaughters",daughter->Pt(),daughter->Eta(),myEt);
		    //cout<<"Xi daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"Xi daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgAntiXiCode){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiXi",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fgPiPlusCode || daughtercode==fgAntiProtonCode || daughtercode==fgPiMinusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedAntiXiDaughters",daughter->Pt(),daughter->Eta(),myEt);
		      //cout<<"AntiXi daughter is a "<<daughter->GetName()<<endl;
		    }
		  }
		  else{
		    //cout<<"AntiXi daughter is a "<<daughter->GetName()<<endl;
		  }
		}
	      }
	      filled = true;
	    }
	    if(pdgCode == fgXi0Code){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgAntiXi0Code){
	      float myEt = Et(part);
	      fSimHadEt += myEt;
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedAntiXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    //============electrons===================================

	    if(pdgCode == fgEPlusCode){
	      float myEt = Et(part);
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedEPlus",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgEMinusCode){
	      float myEt = Et(part);
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedEMinus",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    //============neutrals===================================
	    if(pdgCode == fgGammaCode){
	      TParticle *mom = stack->Particle(part->GetFirstMother());
	      Int_t pdgCodeMom =  mom->GetPDG(0)->PdgCode();
	      //cout<<"I am a gamma and my mom is "<<mom->GetName()<<endl;
	      //We want to separate the gammas by pi0, eta, omega0 but we don't want to double count energy so we get the et from the gamma daughter
	      if(pdgCodeMom == fgEtaCode){
		float myEt = Et(part);
		fSimTotEt += myEt;
		FillHisto2D("EtSimulatedEta",mom->Pt(),mom->Eta(),myEt);
		filled = true;
	      }
	      if(pdgCodeMom == fgPi0Code){
		float myEt = Et(part);
		fSimTotEt += myEt;
		FillHisto2D("EtSimulatedPi0",mom->Pt(),mom->Eta(),myEt);
		filled = true;
	      }
	      if(pdgCodeMom == fgOmega0Code){
		float myEt = Et(part);
		fSimTotEt += myEt;
		FillHisto2D("EtSimulatedOmega0",mom->Pt(),mom->Eta(),myEt);
		filled = true;
	      }
	      if(!filled){
		float myEt = Et(part);
		fSimTotEt += myEt;
		FillHisto2D("EtSimulatedGamma",part->Pt(),part->Eta(),myEt);
		filled = true;
	      }
	    }
	    if(pdgCode == fgEtaCode){
	      float myEt = Et(part);
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedEta",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgPi0Code){
	      float myEt = Et(part);
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedPi0",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fgOmega0Code){
	      float myEt = Et(part);
	      fSimTotEt += myEt;
	      FillHisto2D("EtSimulatedOmega0",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(!filled){
	      //if( strcmp(pc->ParticleClass(),"Baryon")==0 || strcmp(pc->ParticleClass(),"Meson")==0 ){
	      //cout<<"Did not find a place for "<<part->GetName()<<" "<<pdgCode<<" which is a "<<pc->ParticleClass()<<endl;
		//}
	    }
	  }
	}
    }

    if(fSimTotEt>0.0)FillHisto1D("SimTotEt",fSimTotEt,1.0);
    if(fSimHadEt>0.0)FillHisto1D("SimHadEt",fSimHadEt,1.0);
    if(fSimPiKPEt>0.0)FillHisto1D("SimPiKPEt",fSimPiKPEt,1.0);

    if(fInvestigateSmearing){
      //Smearing histograms
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimPtSmeared",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtPtSmeared)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtPtSmeared",fSimPiKPEtPtSmeared,1.0);
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimEfficiencySmeared",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtEfficiencySmeared)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtEfficiencySmeared",fSimPiKPEtEfficiencySmeared,1.0);
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimPtCutSmearedTPC",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtPtCutSmearedTPC)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtPtCutSmearedTPC",fSimPiKPEtPtCutSmearedTPC,1.0);
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimPtCutSmearedITS",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtPtCutSmearedITS)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtPtCutSmearedITS",fSimPiKPEtPtCutSmearedTPC,1.0);
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimPIDSmeared",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtPIDSmeared)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtPIDSmeared",fSimPiKPEtPIDSmeared,1.0);
      if(fSimPiKPEt>0.0) FillHisto2D("SimPiKPEtMinusSimPIDSmearedNoID",fSimPiKPEt,(fSimPiKPEt-fSimPiKPEtPIDSmearedNoID)/fSimPiKPEt,1.0);
      FillHisto1D("SimPiKPEtPIDSmearedNoID",fSimPiKPEtPIDSmearedNoID,1.0);
    }
    return 1;
    
}

void AliAnalysisHadEtMonteCarlo::Init()
{ // Init
    AliAnalysisHadEt::Init();
    if(!fPtSmearer) fPtSmearer = new TRandom();
}
void AliAnalysisHadEtMonteCarlo::CreateHistograms(){
  //for simulated Et only (no reconstruction)
  CreateEtaPtHisto2D(TString("EtSimulatedPiPlus"),TString("Simulated E_{T} from #pi^{+}"));
  CreateEtaPtHisto2D("EtSimulatedPiMinus","Simulated E_{T} from #pi^{-}");
  CreateEtaPtHisto2D("EtSimulatedKPlus","Simulated E_{T} from K^{+}");
  CreateEtaPtHisto2D("EtSimulatedKMinus","Simulated E_{T} from K^{-}");
  CreateEtaPtHisto2D("EtSimulatedProton","Simulated E_{T} from p");
  CreateEtaPtHisto2D("EtSimulatedAntiProton","Simulated E_{T} from #bar{p}");//Both baryon enhancement and strangeness rescaling
  if(fBaryonEnhancement){
    CreateEtaPtHisto2D("EtSimulatedProtonEnhanced","Simulated E_{T} from p");
    CreateEtaPtHisto2D("EtSimulatedAntiProtonEnhanced","Simulated E_{T} from #bar{p}");
  }
  CreateEtaPtHisto2D("EtSimulatedChargedHadron","Simulated E_{T} from charged hadrons");
  CreateEtaPtHisto2D("EtNSimulatedPiPlus","Number of Simulated #pi^{+}");
  CreateEtaPtHisto2D("EtNSimulatedPiMinus","Number of simulated #pi^{-}");
  CreateEtaPtHisto2D("EtNSimulatedKPlus","Number of simulated K^{+}");
  CreateEtaPtHisto2D("EtNSimulatedKMinus","Number of simulated K^{-}");
  CreateEtaPtHisto2D("EtNSimulatedProton","Number of simulated p");
  CreateEtaPtHisto2D("EtNSimulatedAntiProton","Number of simulated #bar{p}");
  if(fBaryonEnhancement){
    CreateEtaPtHisto2D("EtNSimulatedProtonEnhanced","Number of simulated p");
    CreateEtaPtHisto2D("EtNSimulatedAntiProtonEnhanced","Number of simulated #bar{p}");
  }
  CreateEtaPtHisto2D("EtNSimulatedChargedHadron","Number of simulated charged hadrons");
  if(fDataSet==20100){//If this is Pb+Pb
    Int_t width = 5;
    if(fNCentBins<20) width = 10;
    for(Int_t i=0;i<fNCentBins;i++){
      CreateEtaPtHisto2D(Form("EtNSimulatedPiPlusCB%i",i),Form("Number of Simulated #pi^{+} for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedPiMinusCB%i",i),Form("Number of simulated #pi^{-} for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedKPlusCB%i",i),Form("Number of simulated K^{+} for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedKMinusCB%i",i),Form("Number of simulated K^{-} for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedProtonCB%i",i),Form("Number of simulated p for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedAntiProtonCB%i",i),Form("Number of simulated #bar{p} for %i-%i central",i*width,(i+1)*width));
      CreateEtaPtHisto2D(Form("EtNSimulatedChargedHadronCB%i",i),Form("Number of simulated charged hadrons for %i-%i central",i*width,(i+1)*width));
    }
  }
  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingNoPt","Simulated E_{T} from charged hadrons assuming p_{T}=0");
  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingPtTPCCut","Simulated E_{T} from charged hadrons assuming p_{T}=0.15");
  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingPtITSCut","Simulated E_{T} from charged hadrons assuming p_{T}=0.10");

  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingPion","Simulated E_{T} from charged hadrons assuming they are all pions");
  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingProton","Simulated E_{T} from charged hadrons assuming they are all pions");
  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingKaon","Simulated E_{T} from charged hadrons assuming they are all pions");
  CreateEtaPtHisto2D("EtSimulatedKPlusAssumingPion","Simulated E_{T} from K^{+} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedKMinusAssumingPion","Simulated E_{T} from K^{-} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedProtonAssumingPion","Simulated E_{T} from p assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedAntiProtonAssumingPion","Simulated E_{T} from #bar{p} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedKPlusAssumingProton","Simulated E_{T} from K^{+} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedKMinusAssumingProton","Simulated E_{T} from K^{-} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedPiPlusAssumingProton","Simulated E_{T} from p assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedPiMinusAssumingProton","Simulated E_{T} from #bar{p} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedPiPlusAssumingKaon","Simulated E_{T} from K^{+} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedPiMinusAssumingKaon","Simulated E_{T} from K^{-} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedProtonAssumingKaon","Simulated E_{T} from p assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedAntiProtonAssumingKaon","Simulated E_{T} from #bar{p} assuming #pi mass");
  if(fBaryonEnhancement){
    CreateEtaPtHisto2D("EtSimulatedProtonAssumingPionEnhanced","Simulated E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D("EtSimulatedAntiProtonAssumingPionEnhanced","Simulated E_{T} from #bar{p} assuming #pi mass");
  }

  CreateEtaPtHisto2D("EtSimulatedLambda","Simulated E_{T} from #Lambda");
  CreateEtaPtHisto2D("EtSimulatedAntiLambda","Simulated E_{T} from #bar{#Lambda}");
  CreateEtaPtHisto2D("EtSimulatedK0S","Simulated E_{T} from K^{0}_{S}");
  CreateEtaPtHisto2D("EtSimulatedK0L","Simulated E_{T} from K^{0}_{L}");
  CreateEtaPtHisto2D("EtSimulatedLambdaReweighted","Simulated E_{T} from #Lambda");//These will also be used for baryon enhancement
  CreateEtaPtHisto2D("EtSimulatedAntiLambdaReweighted","Simulated E_{T} from #bar{#Lambda}");
  CreateEtaPtHisto2D("EtSimulatedK0SReweighted","Simulated E_{T} from K^{0}_{S}");
  CreateEtaPtHisto2D("EtSimulatedK0LReweighted","Simulated E_{T} from K^{0}_{L}");
  CreateEtaPtHisto2D("EtSimulatedNeutron","Simulated E_{T} from neutrons");
  CreateEtaPtHisto2D("EtSimulatedAntiNeutron","Simulated E_{T} from #bar{n}");
  CreateEtaPtHisto2D("EtSimulatedEPlus","Simulated E_{T} from e^{+}");
  CreateEtaPtHisto2D("EtSimulatedEMinus","Simulated E_{T} from e^{-}");
  CreateEtaPtHisto2D("EtSimulatedOmega","Simulated E_{T} from #Omega^{-}");
  CreateEtaPtHisto2D("EtSimulatedAntiOmega","Simulated E_{T} from #Omega^{+}");
  CreateEtaPtHisto2D("EtSimulatedXi","Simulated E_{T} from #Xi^{-}");
  CreateEtaPtHisto2D("EtSimulatedAntiXi","Simulated E_{T} from #Xi^{+}");
  CreateEtaPtHisto2D("EtSimulatedSigma","Simulated E_{T} from #Xi^{-}");
  CreateEtaPtHisto2D("EtSimulatedAntiSigma","Simulated E_{T} from #Xi^{+}");
  CreateEtaPtHisto2D("EtSimulatedXi0","Simulated E_{T} from #Xi^{0}");
  CreateEtaPtHisto2D("EtSimulatedAntiXi0","Simulated E_{T} from #Xi^{0}");
  CreateEtaPtHisto2D("EtSimulatedAllHadron","Simulated E_{T} from all hadrons");


  CreateEtaPtHisto2D("EtSimulatedLambdaDaughters","Simulated E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiLambdaDaughters","Simulated E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D("EtSimulatedK0SDaughters","Simulated E_{T} from K^{0}_{S} Daughters");
  CreateEtaPtHisto2D("EtSimulatedLambdaDaughtersReweighted","Simulated E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiLambdaDaughtersReweighted","Simulated E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D("EtSimulatedK0SDaughtersReweighted","Simulated E_{T} from K^{0}_{S} Daughters");
  CreateEtaPtHisto2D("EtSimulatedOmegaDaughters","Simulated E_{T} from #Omega^{-} Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiOmegaDaughters","Simulated E_{T} from #Omega^{+} Daughters");
  CreateEtaPtHisto2D("EtSimulatedXiDaughters","Simulated E_{T} from #Xi^{-} Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiXiDaughters","Simulated E_{T} from #Xi^{+} Daughters");


  CreateEtaPtHisto2D("EtSimulatedGamma","Simulated E_{T} from #gamma");
  CreateEtaPtHisto2D("EtSimulatedEta","Simulated E_{T} from #eta");
  CreateEtaPtHisto2D("EtSimulatedPi0","Simulated E_{T} from #pi^{0}");
  CreateEtaPtHisto2D("EtSimulatedOmega0","Simulated E_{T} from #omega");

  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  Int_t lastcutset = 1;
  if(fRequireITSHits) lastcutset = 2;
  for(Int_t i=0;i<=lastcutset;i++){
    TString *cutName = NULL;
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

    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",cutName->Data()),"Reconstructed E_{T} from identified #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",cutName->Data()),"Reconstructed E_{T} from identified #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKPlus",cutName->Data()),"Reconstructed E_{T} from identified K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEMinus",cutName->Data()),"Reconstructed E_{T} from identified e^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEPlus",cutName->Data()),"Reconstructed E_{T} from identified e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKMinus",cutName->Data()),"Reconstructed E_{T} from identified K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedProton",cutName->Data()),"Reconstructed E_{T} from identified p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",cutName->Data()),"Reconstructed E_{T} from identified #bar{p}");
    if(fBaryonEnhancement){
      CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from identified p");
      CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedAntiProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from identified #bar{p}");
    }
    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentified",cutName->Data()),"Number of Reconstructed unidentified particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",cutName->Data()),"Reconstructed E_{T} from unidentified particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedAssumingProton",cutName->Data()),"Reconstructed E_{T} from unidentified particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified particles assuming pion mass");

    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentifiedKaon",cutName->Data()),"Number of Reconstructed unidentified kaons particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingPion",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingProton",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedKaonAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedKaon",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming kaon mass");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentifiedProton",cutName->Data()),"Number of Reconstructed unidentified proton particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingPion",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingProton",cutName->Data()),"Reconstructed E_{T} from unidentified kaons particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingProton",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProton",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming proton mass");
    if(fBaryonEnhancement){
      CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentifiedProtonEnhanced",cutName->Data()),"Number of Reconstructed unidentified proton particles");
      CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonAssumingPionEnhanced",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming pion mass");
      CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from unidentified proton particles assuming proton mass");
    }
    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentifiedPion",cutName->Data()),"Number of Reconstructed unidentified pions particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingPion",cutName->Data()),"Reconstructed E_{T} from unidentified pions particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingKaon",cutName->Data()),"Reconstructed E_{T} from unidentified pions particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPionAssumingProton",cutName->Data()),"Reconstructed E_{T} from unidentified pions particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedPion",cutName->Data()),"Reconstructed E_{T} from unidentified pions particles assuming pion mass");

    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentified",cutName->Data()),"Reconstructed E_{T} from unidentified particles using real mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),"Reconstructed E_{T} from misidentified electrons");


    CreateEtaPtHisto2D(Form("EtReconstructed%sPiPlus",cutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiMinus",cutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlus",cutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinus",cutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProton",cutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProton",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    if(fBaryonEnhancement){
      CreateEtaPtHisto2D(Form("EtReconstructed%sProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from p");
      CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    }
    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),"Reconstructed E_{T} from charged hadrons");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiPlus",cutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiMinus",cutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKPlus",cutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKMinus",cutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sProton",cutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sAntiProton",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    if(fBaryonEnhancement){
      CreateEtaPtHisto2D(Form("EtNReconstructed%sProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from p");
      CreateEtaPtHisto2D(Form("EtNReconstructed%sAntiProtonEnhanced",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    }
    CreateEtaPtHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),"Reconstructed E_{T} from charged hadrons");
    if(fDataSet==20100){//If this is Pb+Pb
      Int_t width = 5;
      if(fNCentBins<20) width = 10;
      for(Int_t j=0;j<fNCentBins;j++){
	CreateEtaPtHisto2D(Form("EtNReconstructed%sPiPlusCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from #pi^{+} for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sPiMinusCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from #pi^{-} for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sKPlusCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from K^{+} for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sKMinusCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from K^{-} for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sProtonCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from p for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sAntiProtonCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from #bar{p} for %i-%i central",j*width,(j+1)*width));
	CreateEtaPtHisto2D(Form("EtNReconstructed%sChargedHadronCB%i",cutName->Data(),j),Form("Reconstructed E_{T} from charged hadrons for %i-%i central",j*width,(j+1)*width));
      }
    }

    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),"Reconstructed E_{T} from charged hadrons assuming they are all pions");
    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadronAssumingProton",cutName->Data()),"Reconstructed E_{T} from charged hadrons assuming they are all pions");
    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadronAssumingKaon",cutName->Data()),"Reconstructed E_{T} from charged hadrons assuming they are all pions");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlusAssumingPion",cutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinusAssumingPion",cutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingPion",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiPlusAssumingKaon",cutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiMinusAssumingKaon",cutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlusAssumingKaon",cutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinusAssumingKaon",cutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingKaon",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingKaon",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");

    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlusAssumingProton",cutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinusAssumingProton",cutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiMinusAssumingProton",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiPlusAssumingProton",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingProton",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingProton",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");

    if(fBaryonEnhancement){
      CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingPionEnhanced",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
      CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingPionEnhanced",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");
    }

    CreateEtaPtHisto2D(Form("EtReconstructed%sEPlus",cutName->Data()),"Reconstructed E_{T} from e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sEMinus",cutName->Data()),"Reconstructed E_{T} from e^{-}");



    CreateEtaPtHisto2D(Form("EtReconstructed%sLambdaDaughters",cutName->Data()),"Reconstructed E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",cutName->Data()),"Reconstructed E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sK0SDaughters",cutName->Data()),"Reconstructed E_{T} from K^{0}_{S} Daughters");
    CreateEtaPtHisto2D(Form("EtReconstructed%sLambdaDaughtersReweighted",cutName->Data()),"Reconstructed E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiLambdaDaughtersReweighted",cutName->Data()),"Reconstructed E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sK0SDaughtersReweighted",cutName->Data()),"Reconstructed E_{T} from K^{0}_{S} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sOmegaDaughters",cutName->Data()),"Reconstructed E_{T} from #Omega^{-} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",cutName->Data()),"Reconstructed E_{T} from #Omega^{+} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sXiDaughters",cutName->Data()),"Reconstructed E_{T} from #Xi^{-} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiXiDaughters",cutName->Data()),"Reconstructed E_{T} from #Xi^{+} Daughters");

    CreateIntHisto1D(Form("UnidentifiedPIDs%s",cutName->Data()),"PIDs of unidentified particles", "PID", "Number of particles",9, -4,4);
    CreateHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),"PIDs of misidentified particles", "PID real","PID identified",5, -.5,4.5,5, -.5,4.5);
    CreateHisto2D(Form("dEdxAll%s",cutName->Data()),"dE/dx for all particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxPion%s",cutName->Data()),"dE/dx for #pi^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxKaon%s",cutName->Data()),"dE/dx for K^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxProton%s",cutName->Data()),"dE/dx for p(#bar{p})","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxElectron%s",cutName->Data()),"dE/dx for e^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxUnidentified%s",cutName->Data()),"dE/dx for unidentified particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
  }
  delete strTPC;
  delete strITS;
  delete strTPCITS;

  Float_t minEt = 0.0;
  Float_t maxEt = 100.0;
  if(fDataSet==20100) maxEt=5000.0;
  Int_t nbinsEt = 100;
  char histoname[200];
  char histotitle[200];
  char xtitle[50];
  char ytitle[50];
  TString *sTPC = new TString("TPC");
  TString *sITS = new TString("ITS");
  TString *sTPCpt = new TString("0.15");
  TString *sITSpt = new TString("0.10");
  TString *sPID = new TString("");
  TString *sNoPID = new TString("NoPID");
  TString *sNoPIDString = new TString(", No PID");
  TString *sHadEt = new TString("HadEt");
  TString *sTotEt = new TString("TotEt");
  TString *sTotEtString = new TString("total E_{T}");
  TString *sHadEtString = new TString("hadronic E_{T}");
  TString *sFull = new TString("Full");
  TString *sEMCAL = new TString("EMCAL");
  TString *sPHOS = new TString("PHOS");
  float etDiff = 1.5;

  for(int tpc = 0;tpc<lastcutset;tpc++){
    TString *detector = NULL;
    TString *ptstring = NULL;
    if(tpc==1) {detector = sTPC; ptstring = sTPCpt;}
    else{detector = sITS; ptstring = sITSpt;}
    for(int hadet = 0;hadet<2;hadet++){
      TString *et = NULL;
      TString *etstring = NULL;
      if(hadet==1) {et = sHadEt; etstring = sHadEtString;}
      else{et = sTotEt; etstring = sTotEtString;}
      for(int type = 0;type<3;type++){
	if(type==0 && !fInvestigateFull) continue;
	if(type==1 && !fInvestigateEMCal) continue;
	if(type==2 && !fInvestigatePHOS) continue;
	TString *acceptance = NULL;
	switch(type){
	case 0:
	  acceptance = sFull;
	  etDiff = 1.5;
	  break;
	case 1:
	  acceptance = sEMCAL;
	  etDiff = 5;
	  break;
	case 2:
	  acceptance = sPHOS;
	  etDiff = 5;
	  break;
	default:
	  acceptance = sFull;
	}
	for(int pid = 0;pid<2;pid++){
	  TString *partid = NULL;
	  TString *partidstring = NULL;
	  if(pid==1){partid = sPID; partidstring = sPID;}
	  else{partid = sNoPID; partidstring = sNoPIDString;}
	  snprintf(histoname,200,"Sim%sMinusReco%s%sAcceptance%s%s",et->Data(),et->Data(),acceptance->Data(),detector->Data(),partid->Data());
	  snprintf(histotitle,200,"(Simulated %s - reconstructed %s)/(Simulated %s) with %s acceptance for p_{T}>%s GeV/c%s",etstring->Data(),etstring->Data(),etstring->Data(),acceptance->Data(),ptstring->Data(),partidstring->Data());
	  snprintf(ytitle,50,"(Simulated %s - reconstructed %s)/(Simulated %s)",etstring->Data(),etstring->Data(),etstring->Data());
	  snprintf(xtitle,50,"Simulated %s",etstring->Data());
	  CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
	  if(hadet==0 && type==0 && fInvestigatePiKP){//we only want to do this once...  not the most elegant way of coding but hey...
	    snprintf(histoname,200,"SimPiKPMinusRecoPiKP%sAcceptance%s%s",acceptance->Data(),detector->Data(),partid->Data());
	    snprintf(histotitle,200,"(Sim PiKP - reco PiKP)/(Sim PiKP) with %s acceptance for p_{T}>%s GeV/c%s",acceptance->Data(),ptstring->Data(),partidstring->Data());
	    snprintf(ytitle,50,"(Sim PiKP - reco PiKP)/(Sim PiKP)");
	    snprintf(xtitle,50,"Simulated E_{T}^{#pi,K,p}");
	    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
	  }
	}
      }
    }
  }
  CreateHisto1D("SimPiKPEt","Simulated #pi,K,p E_{T}","Simulated #pi,K,p E_{T}","Number of events",nbinsEt,minEt,maxEt);
  CreateHisto1D("SimTotEt","Simulated Total E_{T}","Simulated Total E_{T}","Number of events",nbinsEt*4,minEt,maxEt);
  CreateHisto1D("SimHadEt","Simulated Hadronic E_{T}","Simulated Hadronic E_{T}","Number of events",nbinsEt*4,minEt,maxEt);

  etDiff = 0.15;

  if(fInvestigateSmearing){

    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimPtSmeared");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for 0.5 percent momentum smearing");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff/10.0,etDiff/10.0);
    snprintf(histoname,200,"SimPiKPEtPtSmeared");
    snprintf(histotitle,200,"Simulated E_{T} for 0.5 percent momentum smearing");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);

    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimEfficiencySmeared");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for efficiency smearing");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff*5,etDiff*5);
    snprintf(histoname,200,"SimPiKPEtEfficiencySmeared");
    snprintf(histotitle,200,"Simulated E_{T} for efficiency smearing");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);

    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimPtCutSmearedTPC");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for p_{T}>0.15 GeV/c smearing");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
    snprintf(histoname,200,"SimPiKPEtPtCutSmearedTPC");
    snprintf(histotitle,200,"Simulated E_{T} for p_{T}>0.15 GeV/c smearing");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);


    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimPtCutSmearedITS");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for p_{T}>0.10 GeV/c smearing");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
    snprintf(histoname,200,"SimPiKPEtPtCutSmearedITS");
    snprintf(histotitle,200,"Simulated E_{T} for p_{T}>0.10 GeV/c smearing");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);

    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimPIDSmeared");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for PID smearing");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
    snprintf(histoname,200,"SimPiKPEtPIDSmeared");
    snprintf(histotitle,200,"Simulated E_{T} for PID smearing");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);

    //======================================================================

    snprintf(histoname,200,"SimPiKPEtMinusSimPIDSmearedNoID");
    snprintf(histotitle,200,"Simulated (true-smeared)/true for PID smearing No ID");
    snprintf(ytitle,50,"(true-smeared)/true");
    snprintf(xtitle,50,"true p, K, p E_{T}");
    CreateHisto2D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt,nbinsEt,-etDiff,etDiff);
    snprintf(histoname,200,"SimPiKPEtPIDSmearedNoID");
    snprintf(histotitle,200,"Simulated E_{T} for PID smearing No ID");
    snprintf(ytitle,50,"Number of events");
    snprintf(xtitle,50,"p, K, p E_{T}");
    CreateHisto1D(histoname,histotitle,xtitle,ytitle,nbinsEt,minEt,maxEt);
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
  CreateIntHisto1D("NEvents","Number of events","number of events","Number of events",1,0,1);

}

