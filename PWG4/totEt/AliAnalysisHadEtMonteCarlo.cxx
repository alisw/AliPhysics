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
#include <iostream>
using namespace std;


Int_t AliAnalysisHadEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{ // analyse MC and real event info
  FillHisto1D("NEvents",0.5,1);

  AnalyseEvent(ev);
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
  AliStack *stack = mcEvent->Stack();

  //for PID
  AliESDpid *pID = new AliESDpid();
  pID->MakePID(realEvent);

  //This code taken from https://twiki.cern.ch/twiki/bin/view/ALICE/SelectionOfPrimaryTracksForPp2009DataAnalysis
  //Gets good tracks
  //=============================================
  // Primary vertex
  const AliESDVertex *vertex = realEvent->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = realEvent->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) {
      // NO GOOD VERTEX, SKIP EVENT 
    }
  }
  // apply a cut |zVertex| < CUT, if needed

  

  //fEsdtrackCutsITSTPC->SetEtaRange(-0.8,0.8); // normally, |eta|<0.8
  //=============================================

  //Roughly following $ALICE_ROOT/PWG0/dNdEta/AlidNdEtaCorrectionTask

  //=============================================TPC&&ITS=============================================
  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  for(Int_t cutset=0;cutset<3;cutset++){
    TString *cutName;
    TObjArray* list;
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
	  bool isGood = true;
	  if(cutset==1){//if these are ITS stand alone tracks, apply some specific cuts
	    ULong_t trStatus=track->GetStatus();
	    if(trStatus&AliESDtrack::kTPCin) isGood=false; // reject tracks found in TPC
	    if(trStatus&AliESDtrack::kITSpureSA) isGood=false; // reject "pure standalone" ITS tracks
	    if(!(trStatus&AliESDtrack::kITSrefit)) isGood = false; // require proper refit in ITS 
	  }
	  if(!isGood) continue;
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
	  bool unidentified = (!isProton && !isKaon && !isElectron);
	  Float_t dEdx = track->GetTPCsignal();
	  if(cutset==1) dEdx = track->GetITSsignal();

	  FillHisto2D(Form("dEdxAll%s",cutName->Data()),track->P(),dEdx,1.0);
	  //if(cutset==1) cout<<"filling "<<track->P()<<" "<<dEdx<<endl;

	  UInt_t label = (UInt_t)TMath::Abs(track->GetLabel());
	  TParticle  *simPart  = stack->Particle(label);
	  if(!simPart) {
	    Printf("no MC particle\n"); 	 	
	    continue; 	 	
	  }
	  else{//analysis
	    if(stack->IsPhysicalPrimary(label)){
	      if (TMath::Abs(simPart->Eta()) < fCuts->GetCommonEtaCut())	    {

		Int_t pdgCode =  simPart->GetPDG(0)->PdgCode();
		Int_t mypid = 0;
		if(pdgCode==fPiPlusCode) mypid = 1;
		if(pdgCode==fProtonCode) mypid = 2;
		if(pdgCode==fKPlusCode) mypid = 3;
		if(pdgCode==fEPlusCode) mypid = 4;
		if(pdgCode==fPiMinusCode) mypid = 1;
		if(pdgCode==fAntiProtonCode) mypid = 2;
		if(pdgCode==fKMinusCode) mypid = 3;
		if(pdgCode==fEMinusCode) mypid = 4;
		//cout<<pdgCode->PdgCode()<<" ";
		//fPdgDB->GetSimParticle("pi+")->PdgCode();
		bool filled = false;      
		//============Charged hadrons===================================
		//identified...
		if(isPion){
		  if(pdgCode!=fPiPlusCode && pdgCode!=fPiMinusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),1,mypid,1);
		    //if(mypid==0)cerr<<"I was misidentified! I'm not a pion! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxPion%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isProton){
		  if(pdgCode!=fProtonCode && pdgCode!=fAntiProtonCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),2,mypid,1);
		    // if(mypid==0)cerr<<"I was misidentified!  I'm not a proton! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedProton",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxProton%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isKaon){
		  if(pdgCode!=fKMinusCode && pdgCode!=fKPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),3,mypid,1);
		    //if(mypid==0)cerr<<"I was misidentified!  I'm not a kaon! I am a "<<simPart->GetName()<<" p "<<track->P()<<" nSigmaProton "<<nSigmaProton<<" nSigmaPion "<<nSigmaPion<<" nSigmaKaon "<<nSigmaKaon<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedKPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedKMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxKaon%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(isElectron){
		  if(pdgCode!=fEMinusCode && pdgCode!=fEPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",cutName->Data()),4,mypid,1);
		    //cerr<<"I was misidentified!  I'm not an electron! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedEPlus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedEMinus",cutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxElectron%s",cutName->Data()),track->P(),dEdx,1.0);
		}
		if(unidentified){
		  if(pdgCode!=fEMinusCode && pdgCode!=fEPlusCode){
		    float myEtPi = Et(simPart,fPionMass);
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		    FillHisto2D(Form("EtReconstructed%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sUnidentified",cutName->Data()),track->Pt(),track->Eta(),1.0);
		  }
		  FillHisto2D(Form("dEdxUnidentified%s",cutName->Data()),track->P(),dEdx,1.0);
		  //cout<<"I was not identified.  I am a "<<simPart->GetName()<<" PID "<<pdgCode<<endl;
		  //track what was not identified successfully
		  FillHisto1D(Form("UnidentifiedPIDs%s",cutName->Data()),mypid,1);
		}
		//...simulated
		if(pdgCode == fPiPlusCode){
		  //cout<<"I'm a real primary "<<simPart->GetName()<<"! "<<"my label is "<<simPart->GetFirstMother()<<" track no "<<iTrack<<"/"<<realEvent->GetNumberOfTracks()<<endl;//<<" "<<label<<" "<<pdgCode<<endl;
		
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sPiPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;
		}
		if(pdgCode == fPiMinusCode){
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sPiMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;
		}
		if(pdgCode == fKPlusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fPionMass);
		  FillHisto2D(Form("EtReconstructed%sKPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKPlusAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == fKMinusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fPionMass);
		  FillHisto2D(Form("EtReconstructed%sKMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKMinusAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == fProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fPionMass);
		  FillHisto2D(Form("EtReconstructed%sProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sProtonAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == fAntiProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,fPionMass);
		  FillHisto2D(Form("EtReconstructed%sAntiProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sAntiProton",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",cutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == fEPlusCode){
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sEPlus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(!isElectron || unidentified){
		    float myEtPi = Et(simPart,fPionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  filled = true;
		}
		if(pdgCode == fEMinusCode){
		  if(!isElectron || unidentified){
		    float myEtPi = Et(simPart,fPionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sEMinus",cutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		    filled = true;
		}
		//if(!filled){
		  //TParticlePDG *pc = simPart->GetPDG(0);
		  //if( strcmp(pc->ParticleClass(),"Baryon")==0 || strcmp(pc->ParticleClass(),"Meson")==0 ){
		  //cout<<"Did not find a place for "<<simPart->GetName()<<" "<<pdgCode<<" which is a "<<pc->ParticleClass()<<endl;
		  //}
		  //}
	      }
	      
	    }
	    else{//not a primary - we're after V0 daughters!
	      //cout<<"I'm a secondary "<<simPart->GetName()<<"!";//<<endl;
	      TParticle *mom = stack->Particle(simPart->GetFirstMother());
	      if(mom){
		TParticlePDG *pc = mom->GetPDG(0);
		if(pc){
		  Int_t pdgCode =  mom->GetPDG(0)->PdgCode();
		  if(pdgCode == fLambdaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sLambdaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fAntiLambdaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fK0SCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sK0SDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fXiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fAntiXiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fOmegaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == fXiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
		  }

		  if(mom->GetFirstMother()>0){
		    TParticle *grandma = stack->Particle(mom->GetFirstMother());
		    if(grandma){
		      Int_t pdgCodeMom =  mom->GetPDG(0)->PdgCode();
		      if(pdgCodeMom==fPiPlusCode || pdgCodeMom==fPiMinusCode || pdgCodeMom==fProtonCode ||pdgCodeMom==fAntiProtonCode || pdgCodeMom==fKPlusCode || pdgCode==fKMinusCode){
			//cout<<" my grandmother is "<<grandma->GetName()<<" "<<endl;
			Int_t pdgCodeGrandma =  grandma->GetPDG(0)->PdgCode();
		      
			if(pdgCodeGrandma == fXiCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == fAntiXiCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == fOmegaCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sOmegaDaughters",cutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == fXiCode){
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
    delete list;
  }
  //delete AliESDpid;
  return 1;
}
Int_t AliAnalysisHadEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
     ResetEventValues();
     
    // Get us an mc event
    AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);

    // Let's play with the stack!
    AliStack *stack = mcEvent->Stack();

    Int_t nPrim = stack->GetNtrack();

    //=================Tracks which may or may not have been reconstructed=================

    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

        TParticle *part = stack->Particle(iPart);

        if (!part)
	  {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
	  }

        TParticlePDG *pc = part->GetPDG(0);

        // Check if it is a primary particle
        if (!stack->IsPhysicalPrimary(iPart)){//secondaries...
	}
	else{//primaries
	  // Check for reasonable (for now neutral and singly charged) charge on the particle
	  //note that the charge is stored in units of e/3
	  //if (TMath::Abs(pc->Charge()) != EtMonteCarloCuts::kSingleChargedParticle && pc->Charge() != EtMonteCarloCuts::kNeutralParticle) continue;


	  if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())	    {

	    Int_t pdgCode =  part->GetPDG(0)->PdgCode();
	    //cout<<pdgCode->PdgCode()<<" ";
	    //fPdgDB->GetParticle("pi+")->PdgCode();
	    bool filled = false;
	    //============Charged hadrons===================================
	    if(pdgCode == fPiPlusCode){
	      //cout<<"I'm a simulated primary "<<part->GetName()<<"! "<<"my label is "<<part->GetFirstMother()<<" pt "<<part->Pt()<<endl;
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedPiPlus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedPiPlus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fPiMinusCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedPiMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedPiMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fKPlusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fPionMass);
	      FillHisto2D("EtSimulatedKPlus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKPlus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKPlusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fKMinusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fPionMass);
	      FillHisto2D("EtSimulatedKMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKMinusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fPionMass);
	      FillHisto2D("EtSimulatedProton",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedProton",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedProtonAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fAntiProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,fPionMass);
	      FillHisto2D("EtSimulatedAntiProton",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedAntiProton",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAntiProtonAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    //============Other hadrons===================================

	    if(pdgCode == fNeutronCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fAntiNeutronCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fLambdaCode){
	      float myEt = Et(part);
	      //cout<<"I am a simulated lambda! pt "<<part->Pt()<<" eta "<<part->Eta()<<endl;
	      FillHisto2D("EtSimulatedLambda",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fPiMinusCode || daughtercode==fProtonCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedLambdaDaughters",daughter->Pt(),daughter->Eta(),myEt);
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
	    if(pdgCode == fAntiLambdaCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiLambda",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){
		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fPiPlusCode || daughtercode==fAntiProtonCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedAntiLambdaDaughters",daughter->Pt(),daughter->Eta(),myEt);
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
	    if(pdgCode == fK0SCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedK0S",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      Int_t ndaughters = part->GetNDaughters();
	      for(Int_t idaughter = 0;idaughter<ndaughters;idaughter++){
		Int_t daughterindex = part->GetDaughter(idaughter);
		if(daughterindex<0 || daughterindex>1e5) continue;
		TParticle *daughter = stack->ParticleFromTreeK(daughterindex);
		if(daughter){
		  if(daughter->GetPDG(0)){

		    Int_t daughtercode = daughter->GetPDG(0)->PdgCode();
		    if(daughtercode==fPiMinusCode || daughtercode==fPiPlusCode){
		      myEt = Et(daughter);
		      FillHisto2D("EtSimulatedK0SDaughters",daughter->Pt(),daughter->Eta(),myEt);
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
	    if(pdgCode == fK0LCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedK0L",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fOmegaCode){
	      float myEt = Et(part);
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
		    if(daughtercode==fPiPlusCode || daughtercode==fProtonCode || daughtercode==fKMinusCode){
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
	    if(pdgCode == fAntiOmegaCode){
	      float myEt = Et(part);
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
		    if(daughtercode==fPiMinusCode || daughtercode==fAntiProtonCode || daughtercode==fKPlusCode){
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
	    if(pdgCode == fSigmaCode || pdgCode == -3222){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fAntiSigmaCode || pdgCode == 3222){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fXiCode){
	      float myEt = Et(part);
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
		    if(daughtercode==fPiPlusCode || daughtercode==fProtonCode || daughtercode==fPiMinusCode){
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
	    if(pdgCode == fAntiXiCode){
	      float myEt = Et(part);
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
		    if(daughtercode==fPiPlusCode || daughtercode==fAntiProtonCode || daughtercode==fPiMinusCode){
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
	    if(pdgCode == fXi0Code){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fAntiXi0Code){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    //============electrons===================================

	    if(pdgCode == fEPlusCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedEPlus",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == fEMinusCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedEMinus",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(!filled){
	      if( strcmp(pc->ParticleClass(),"Baryon")==0 || strcmp(pc->ParticleClass(),"Meson")==0 ){
		//cout<<"Did not find a place for "<<part->GetName()<<" "<<pdgCode<<" which is a "<<pc->ParticleClass()<<endl;
	      }
	    }
	  }
	}
    }




//     fTotNeutralEtAcc = fTotNeutralEt;
//     fTotEt = fTotChargedEt + fTotNeutralEt;
//     fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;
    
//     FillHistograms();

    return 1;
    
}

void AliAnalysisHadEtMonteCarlo::Init()
{ // Init
    AliAnalysisHadEt::Init();
}

void AliAnalysisHadEtMonteCarlo::CreateHistograms(){
  //for simulated Et only (no reconstruction)
  CreateEtaPtHisto2D(TString("EtSimulatedPiPlus"),TString("Simulated E_{T} from #pi^{+}"));
  CreateEtaPtHisto2D("EtSimulatedPiMinus","Simulated E_{T} from #pi^{-}");
  CreateEtaPtHisto2D("EtSimulatedKPlus","Simulated E_{T} from K^{+}");
  CreateEtaPtHisto2D("EtSimulatedKMinus","Simulated E_{T} from K^{-}");
  CreateEtaPtHisto2D("EtSimulatedProton","Simulated E_{T} from p");
  CreateEtaPtHisto2D("EtSimulatedAntiProton","Simulated E_{T} from #bar{p}");
  CreateEtaPtHisto2D("EtSimulatedChargedHadron","Simulated E_{T} from charged hadrons");
  CreateEtaPtHisto2D(TString("EtNSimulatedPiPlus"),TString("Number of Simulated #pi^{+}"));
  CreateEtaPtHisto2D("EtNSimulatedPiMinus","Number of simulated #pi^{-}");
  CreateEtaPtHisto2D("EtNSimulatedKPlus","Number of simulated K^{+}");
  CreateEtaPtHisto2D("EtNSimulatedKMinus","Number of simulated K^{-}");
  CreateEtaPtHisto2D("EtNSimulatedProton","Number of simulated p");
  CreateEtaPtHisto2D("EtNSimulatedAntiProton","Number of simulated #bar{p}");
  CreateEtaPtHisto2D("EtNSimulatedChargedHadron","Number of simulated charged hadrons");

  CreateEtaPtHisto2D("EtSimulatedChargedHadronAssumingPion","Simulated E_{T} from charged hadrons assuming they are all pions");
  CreateEtaPtHisto2D("EtSimulatedKPlusAssumingPion","Simulated E_{T} from K^{+} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedKMinusAssumingPion","Simulated E_{T} from K^{-} assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedProtonAssumingPion","Simulated E_{T} from p assuming #pi mass");
  CreateEtaPtHisto2D("EtSimulatedAntiProtonAssumingPion","Simulated E_{T} from #bar{p} assuming #pi mass");

  CreateEtaPtHisto2D("EtSimulatedLambda","Simulated E_{T} from #Lambda");
  CreateEtaPtHisto2D("EtSimulatedAntiLambda","Simulated E_{T} from #bar{#Lambda}");
  CreateEtaPtHisto2D("EtSimulatedK0S","Simulated E_{T} from K^{0}_{S}");
  CreateEtaPtHisto2D("EtSimulatedK0L","Simulated E_{T} from K^{0}_{L}");
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
  CreateEtaPtHisto2D("EtSimulatedOmegaDaughters","Simulated E_{T} from #Omega^{-} Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiOmegaDaughters","Simulated E_{T} from #Omega^{+} Daughters");
  CreateEtaPtHisto2D("EtSimulatedXiDaughters","Simulated E_{T} from #Xi^{-} Daughters");
  CreateEtaPtHisto2D("EtSimulatedAntiXiDaughters","Simulated E_{T} from #Xi^{+} Daughters");

  TString *strTPC = new TString("TPC");
  TString *strITS = new TString("ITS");
  TString *strTPCITS = new TString("TPCITS");
  for(Int_t i=0;i<3;i++){
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

    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",cutName->Data()),"Reconstructed E_{T} from identified #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",cutName->Data()),"Reconstructed E_{T} from identified #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKPlus",cutName->Data()),"Reconstructed E_{T} from identified K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEMinus",cutName->Data()),"Reconstructed E_{T} from identified e^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEPlus",cutName->Data()),"Reconstructed E_{T} from identified e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKMinus",cutName->Data()),"Reconstructed E_{T} from identified K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedProton",cutName->Data()),"Reconstructed E_{T} from identified p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",cutName->Data()),"Reconstructed E_{T} from identified #bar{p}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentified",cutName->Data()),"Number of Reconstructed unidentified particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",cutName->Data()),"Reconstructed E_{T} from unidentified particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentified",cutName->Data()),"Reconstructed E_{T} from unidentified particles using real mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",cutName->Data()),"Reconstructed E_{T} from misidentified electrons");


    CreateEtaPtHisto2D(Form("EtReconstructed%sPiPlus",cutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiMinus",cutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlus",cutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinus",cutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProton",cutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProton",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadron",cutName->Data()),"Reconstructed E_{T} from charged hadrons");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiPlus",cutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiMinus",cutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKPlus",cutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKMinus",cutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sProton",cutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sAntiProton",cutName->Data()),"Reconstructed E_{T} from #bar{p}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sChargedHadron",cutName->Data()),"Reconstructed E_{T} from charged hadrons");

    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",cutName->Data()),"Reconstructed E_{T} from charged hadrons assuming they are all pions");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlusAssumingPion",cutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinusAssumingPion",cutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingPion",cutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",cutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");

    CreateEtaPtHisto2D(Form("EtReconstructed%sEPlus",cutName->Data()),"Reconstructed E_{T} from e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sEMinus",cutName->Data()),"Reconstructed E_{T} from e^{-}");



    CreateEtaPtHisto2D(Form("EtReconstructed%sLambdaDaughters",cutName->Data()),"Reconstructed E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",cutName->Data()),"Reconstructed E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sK0SDaughters",cutName->Data()),"Reconstructed E_{T} from K^{0}_{S} Daughters");
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

  CreateIntHisto1D("NEvents","Number of events","number of events","Number of events",1,0,1);

  //CreateHisto1D("MisidentifiedPIDs","PIDs for particles misidentified that are not a #pi,K,p","PID","number of entries",3000,0.5,3000.5);



//     list->Add(fHistEt);
//     TString histname = "fHistEt" + fHistogramNameSuffix;

//     fHistEt = new TH1F(histname.Data(), "Total E_{T} Distribution", 1000, 0.00, 99);
//     fHistEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
//     fHistEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");
}

