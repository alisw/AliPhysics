//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#include "AliAnalysisHadEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"

#include <iostream>
using namespace std;


Int_t AliAnalysisHadEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2){
  FillHisto1D("NEvents",0.5,1);

  AnalyseEvent(ev);
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
  AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
  AliStack *stack = mcEvent->Stack();

  //for PID
  AliESDpid *PID = new AliESDpid();
  PID->MakePID(realEvent);

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

  

  //esdtrackCutsITSTPC->SetEtaRange(-0.8,0.8); // normally, |eta|<0.8
  //=============================================

  //Roughly following $ALICE_ROOT/PWG0/dNdEta/AlidNdEtaCorrectionTask

  //=============================================TPC&&ITS=============================================
  TString *TPC = new TString("TPC");
  TString *ITS = new TString("ITS");
  TString *TPCITS = new TString("TPCITS");
  for(Int_t cutset=0;cutset<3;cutset++){
    TString *CutName;
    TObjArray* list;
    switch(cutset){
    case 0:
      CutName = TPC;
      list = esdtrackCutsTPC->GetAcceptedTracks(realEvent);
      break;
    case 1:
      CutName = ITS;
      list = esdtrackCutsITS->GetAcceptedTracks(realEvent);
      break;
    case 2:
      CutName = TPCITS;
      list = esdtrackCutsITSTPC->GetAcceptedTracks(realEvent);
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
	  bool IsGood = true;
	  if(cutset==1){//if these are ITS stand alone tracks, apply some specific cuts
	    ULong_t trStatus=track->GetStatus();
	    if(trStatus&AliESDtrack::kTPCin) IsGood=false; // reject tracks found in TPC
	    if(trStatus&AliESDtrack::kITSpureSA) IsGood=false; // reject "pure standalone" ITS tracks
	    if(!(trStatus&AliESDtrack::kITSrefit)) IsGood = false; // require proper refit in ITS 
	  }
	  if(!IsGood) continue;
	  Float_t nSigmaPion,nSigmaProton,nSigmaKaon,nSigmaElectron;
	  if(cutset!=1){
	    nSigmaPion = TMath::Abs(PID->NumberOfSigmasTPC(track,AliPID::kPion));
	    nSigmaProton = TMath::Abs(PID->NumberOfSigmasTPC(track,AliPID::kProton));
	    nSigmaKaon = TMath::Abs(PID->NumberOfSigmasTPC(track,AliPID::kKaon));
	    nSigmaElectron = TMath::Abs(PID->NumberOfSigmasTPC(track,AliPID::kElectron));
	  }
	  else{
	    nSigmaPion = TMath::Abs(PID->NumberOfSigmasITS(track,AliPID::kPion));
	    nSigmaProton = TMath::Abs(PID->NumberOfSigmasITS(track,AliPID::kProton));
	    nSigmaKaon = TMath::Abs(PID->NumberOfSigmasITS(track,AliPID::kKaon));
	    nSigmaElectron = TMath::Abs(PID->NumberOfSigmasITS(track,AliPID::kElectron));
	  }
	  bool IsPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
	  bool IsElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
	  bool IsKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
	  bool IsProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);

	  //bool IsElectron = false;
	  bool Unidentified = (!IsProton && !IsKaon && !IsElectron);
	  Float_t dEdx = track->GetTPCsignal();
	  if(cutset==1) dEdx = track->GetITSsignal();

	  FillHisto2D(Form("dEdxAll%s",CutName->Data()),track->P(),dEdx,1.0);
	  //if(cutset==1) cout<<"filling "<<track->P()<<" "<<dEdx<<endl;

	  UInt_t label = (UInt_t)TMath::Abs(track->GetLabel());
	  TParticle  *simPart  = stack->Particle(label);
	  if(!simPart) {
	    Printf("no MC particle\n"); 	 	
	    continue; 	 	
	  }
	  else{//analysis
	    if(stack->IsPhysicalPrimary(label)){
	      if (TMath::Abs(simPart->Eta()) < fEtaCut)	    {

		Int_t pdgCode =  simPart->GetPDG(0)->PdgCode();
		Int_t mypid = 0;
		if(pdgCode==PiPlusCode) mypid = 1;
		if(pdgCode==ProtonCode) mypid = 2;
		if(pdgCode==KPlusCode) mypid = 3;
		if(pdgCode==EPlusCode) mypid = 4;
		if(pdgCode==PiMinusCode) mypid = 1;
		if(pdgCode==AntiProtonCode) mypid = 2;
		if(pdgCode==KMinusCode) mypid = 3;
		if(pdgCode==EMinusCode) mypid = 4;
		//cout<<pdgCode->PdgCode()<<" ";
		//fPdgDB->GetSimParticle("pi+")->PdgCode();
		bool filled = false;      
		//============Charged hadrons===================================
		//identified...
		if(IsPion){
		  if(pdgCode!=PiPlusCode && pdgCode!=PiMinusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",CutName->Data()),1,mypid,1);
		    //if(mypid==0)cerr<<"I was misidentified! I'm not a pion! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxPion%s",CutName->Data()),track->P(),dEdx,1.0);
		}
		if(IsProton){
		  if(pdgCode!=ProtonCode && pdgCode!=AntiProtonCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",CutName->Data()),2,mypid,1);
		    // if(mypid==0)cerr<<"I was misidentified!  I'm not a proton! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedProton",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxProton%s",CutName->Data()),track->P(),dEdx,1.0);
		}
		if(IsKaon){
		  if(pdgCode!=KMinusCode && pdgCode!=KPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",CutName->Data()),3,mypid,1);
		    //if(mypid==0)cerr<<"I was misidentified!  I'm not a kaon! I am a "<<simPart->GetName()<<" p "<<track->P()<<" nSigmaProton "<<nSigmaProton<<" nSigmaPion "<<nSigmaPion<<" nSigmaKaon "<<nSigmaKaon<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedKPlus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedKMinus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxKaon%s",CutName->Data()),track->P(),dEdx,1.0);
		}
		if(IsElectron){
		  if(pdgCode!=EMinusCode && pdgCode!=EPlusCode){
		    FillHisto2D(Form("MisidentifiedPIDs%s",CutName->Data()),4,mypid,1);
		    //cerr<<"I was misidentified!  I'm not an electron! I am a "<<simPart->GetName()<<endl;
		  }
		  float myEt = Et(simPart);
		  if(track->Charge()>0){ FillHisto2D(Form("EtReconstructed%sIdentifiedEPlus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  else{ FillHisto2D(Form("EtReconstructed%sIdentifiedEMinus",CutName->Data()),track->Pt(),track->Eta(),myEt);}
		  FillHisto2D(Form("dEdxElectron%s",CutName->Data()),track->P(),dEdx,1.0);
		}
		if(Unidentified){
		  if(pdgCode!=EMinusCode && pdgCode!=EPlusCode){
		    float myEtPi = Et(simPart,PionMass);
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",CutName->Data()),track->Pt(),track->Eta(),myEtPi);
		    FillHisto2D(Form("EtReconstructed%sUnidentified",CutName->Data()),track->Pt(),track->Eta(),myEt);
		    FillHisto2D(Form("EtNReconstructed%sUnidentified",CutName->Data()),track->Pt(),track->Eta(),1.0);
		  }
		  FillHisto2D(Form("dEdxUnidentified%s",CutName->Data()),track->P(),dEdx,1.0);
		  //cout<<"I was not identified.  I am a "<<simPart->GetName()<<" PID "<<pdgCode<<endl;
		  //track what was not identified successfully
		  FillHisto1D(Form("UnidentifiedPIDs%s",CutName->Data()),mypid,1);
		}
		//...simulated
		if(pdgCode == PiPlusCode){
		  //cout<<"I'm a real primary "<<simPart->GetName()<<"! "<<"my label is "<<simPart->GetFirstMother()<<" track no "<<iTrack<<"/"<<realEvent->GetNumberOfTracks()<<endl;//<<" "<<label<<" "<<pdgCode<<endl;
		
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sPiPlus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiPlus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;
		}
		if(pdgCode == PiMinusCode){
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sPiMinus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sPiMinus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  filled = true;
		}
		if(pdgCode == KPlusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,PionMass);
		  FillHisto2D(Form("EtReconstructed%sKPlus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKPlus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKPlusAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == KMinusCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,PionMass);
		  FillHisto2D(Form("EtReconstructed%sKMinus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sKMinus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sKMinusAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == ProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,PionMass);
		  FillHisto2D(Form("EtReconstructed%sProton",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sProton",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sProtonAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == AntiProtonCode){
		  float myEt = Et(simPart);
		  float myEtPi = Et(simPart,PionMass);
		  FillHisto2D(Form("EtReconstructed%sAntiProton",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sAntiProton",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  FillHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  FillHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",CutName->Data()),simPart->Pt(),simPart->Eta(),myEtPi);
		  filled = true;
		}
		if(pdgCode == EPlusCode){
		  float myEt = Et(simPart);
		  FillHisto2D(Form("EtReconstructed%sEPlus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
		  if(!IsElectron || Unidentified){
		    float myEtPi = Et(simPart,PionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",CutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  filled = true;
		}
		if(pdgCode == EMinusCode){
		  if(!IsElectron || Unidentified){
		    float myEtPi = Et(simPart,PionMass);
		    FillHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",CutName->Data()),track->Pt(),track->Eta(),myEtPi);
		  }
		  float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sEMinus",CutName->Data()),simPart->Pt(),simPart->Eta(),myEt);
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
		  if(pdgCode == LambdaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sLambdaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == AntiLambdaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == K0SCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sK0SDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == XiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sXiDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == AntiXiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == OmegaCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sOmegaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }
		  if(pdgCode == XiCode){
		    float myEt = Et(simPart);
		    FillHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
		  }

		  if(mom->GetFirstMother()>0){
		    TParticle *grandma = stack->Particle(mom->GetFirstMother());
		    if(grandma){
		      Int_t pdgCodeMom =  mom->GetPDG(0)->PdgCode();
		      if(pdgCodeMom==PiPlusCode || pdgCodeMom==PiMinusCode || pdgCodeMom==ProtonCode ||pdgCodeMom==AntiProtonCode || pdgCodeMom==KPlusCode || pdgCode==KMinusCode){
			//cout<<" my grandmother is "<<grandma->GetName()<<" "<<endl;
			Int_t pdgCodeGrandma =  grandma->GetPDG(0)->PdgCode();
		      
			if(pdgCodeGrandma == XiCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sXiDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == AntiXiCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sAntiXiDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == OmegaCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sOmegaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
			}
			if(pdgCodeGrandma == XiCode){
			  float myEt = Et(simPart);
			  FillHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",CutName->Data()),track->Pt(),track->Eta(),myEt);
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
{
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


	  if (TMath::Abs(part->Eta()) < fEtaCut)	    {

	    Int_t pdgCode =  part->GetPDG(0)->PdgCode();
	    //cout<<pdgCode->PdgCode()<<" ";
	    //fPdgDB->GetParticle("pi+")->PdgCode();
	    bool filled = false;
	    //============Charged hadrons===================================
	    if(pdgCode == PiPlusCode){
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
	    if(pdgCode == PiMinusCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedPiMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedPiMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == KPlusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,PionMass);
	      FillHisto2D("EtSimulatedKPlus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKPlus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKPlusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == KMinusCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,PionMass);
	      FillHisto2D("EtSimulatedKMinus",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedKMinus",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedKMinusAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == ProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,PionMass);
	      FillHisto2D("EtSimulatedProton",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedProton",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtNSimulatedChargedHadron",part->Pt(),part->Eta(),1.0);
	      FillHisto2D("EtSimulatedChargedHadronAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedProtonAssumingPion",part->Pt(),part->Eta(),myEtPi);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == AntiProtonCode){
	      float myEt = Et(part);
	      float myEtPi = Et(part,PionMass);
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

	    if(pdgCode == NeutronCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == AntiNeutronCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiNeutron",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == LambdaCode){
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
		    if(daughtercode==PiMinusCode || daughtercode==ProtonCode){
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
	    if(pdgCode == AntiLambdaCode){
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
		    if(daughtercode==PiPlusCode || daughtercode==AntiProtonCode){
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
	    if(pdgCode == K0SCode){
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
		    if(daughtercode==PiMinusCode || daughtercode==PiPlusCode){
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
	    if(pdgCode == K0LCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedK0L",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == OmegaCode){
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
		    if(daughtercode==PiPlusCode || daughtercode==ProtonCode || daughtercode==KMinusCode){
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
	    if(pdgCode == AntiOmegaCode){
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
		    if(daughtercode==PiMinusCode || daughtercode==AntiProtonCode || daughtercode==KPlusCode){
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
	    if(pdgCode == SigmaCode || pdgCode == -3222){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == AntiSigmaCode || pdgCode == 3222){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiSigma",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == XiCode){
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
		    if(daughtercode==PiPlusCode || daughtercode==ProtonCode || daughtercode==PiMinusCode){
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
	    if(pdgCode == AntiXiCode){
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
		    if(daughtercode==PiPlusCode || daughtercode==AntiProtonCode || daughtercode==PiMinusCode){
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
	    if(pdgCode == Xi0Code){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == AntiXi0Code){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedAntiXi0",part->Pt(),part->Eta(),myEt);
	      FillHisto2D("EtSimulatedAllHadron",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    //============electrons===================================

	    if(pdgCode == EPlusCode){
	      float myEt = Et(part);
	      FillHisto2D("EtSimulatedEPlus",part->Pt(),part->Eta(),myEt);
	      filled = true;
	    }
	    if(pdgCode == EMinusCode){
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
{

    AliAnalysisHadEt::Init();

    fVertexXCut = EtReconstructedCuts::kVertexXCut;
    fVertexYCut = EtReconstructedCuts::kVertexYCut;
    fVertexZCut = EtReconstructedCuts::kVertexZCut;
    fIPxyCut = EtReconstructedCuts::kIPxyCut;
    fIPzCut = EtReconstructedCuts::kIPzCut;
    // Track cuts
    //Bool_t selectPrimaries=kTRUE;
    //esdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selectPrimaries);
    //esdtrackCutsITSTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

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

  TString *TPC = new TString("TPC");
  TString *ITS = new TString("ITS");
  TString *TPCITS = new TString("TPCITS");
  for(Int_t i=0;i<3;i++){
    TString *CutName;
    Float_t maxPtdEdx = 10;
    Float_t mindEdx = 35;
    Float_t maxdEdx = 150.0;
    switch(i){
    case 0:
      CutName = TPC;
      break;
    case 1:
      CutName = ITS;
      maxPtdEdx = 5;
      maxdEdx = 500.0;
      break;
    case 2:
      CutName = TPCITS;
      break;
    default:
      cerr<<"Error:  cannot make histograms!"<<endl;
      return;
    }

    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiPlus",CutName->Data()),"Reconstructed E_{T} from identified #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedPiMinus",CutName->Data()),"Reconstructed E_{T} from identified #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKPlus",CutName->Data()),"Reconstructed E_{T} from identified K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEMinus",CutName->Data()),"Reconstructed E_{T} from identified e^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedEPlus",CutName->Data()),"Reconstructed E_{T} from identified e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedKMinus",CutName->Data()),"Reconstructed E_{T} from identified K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedProton",CutName->Data()),"Reconstructed E_{T} from identified p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sIdentifiedAntiProton",CutName->Data()),"Reconstructed E_{T} from identified #bar{p}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sUnidentified",CutName->Data()),"Number of Reconstructed unidentified particles");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentifiedAssumingPion",CutName->Data()),"Reconstructed E_{T} from unidentified particles assuming pion mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sUnidentified",CutName->Data()),"Reconstructed E_{T} from unidentified particles using real mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sMisidentifiedElectrons",CutName->Data()),"Reconstructed E_{T} from misidentified electrons");


    CreateEtaPtHisto2D(Form("EtReconstructed%sPiPlus",CutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sPiMinus",CutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlus",CutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinus",CutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProton",CutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProton",CutName->Data()),"Reconstructed E_{T} from #bar{p}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadron",CutName->Data()),"Reconstructed E_{T} from charged hadrons");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiPlus",CutName->Data()),"Reconstructed E_{T} from #pi^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sPiMinus",CutName->Data()),"Reconstructed E_{T} from #pi^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKPlus",CutName->Data()),"Reconstructed E_{T} from K^{+}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sKMinus",CutName->Data()),"Reconstructed E_{T} from K^{-}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sProton",CutName->Data()),"Reconstructed E_{T} from p");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sAntiProton",CutName->Data()),"Reconstructed E_{T} from #bar{p}");
    CreateEtaPtHisto2D(Form("EtNReconstructed%sChargedHadron",CutName->Data()),"Reconstructed E_{T} from charged hadrons");

    CreateEtaPtHisto2D(Form("EtReconstructed%sChargedHadronAssumingPion",CutName->Data()),"Reconstructed E_{T} from charged hadrons assuming they are all pions");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKPlusAssumingPion",CutName->Data()),"Reconstructed E_{T} from K^{+} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sKMinusAssumingPion",CutName->Data()),"Reconstructed E_{T} from K^{-} assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sProtonAssumingPion",CutName->Data()),"Reconstructed E_{T} from p assuming #pi mass");
    CreateEtaPtHisto2D(Form("EtReconstructed%sAntiProtonAssumingPion",CutName->Data()),"Reconstructed E_{T} from #bar{p} assuming #pi mass");

    CreateEtaPtHisto2D(Form("EtReconstructed%sEPlus",CutName->Data()),"Reconstructed E_{T} from e^{+}");
    CreateEtaPtHisto2D(Form("EtReconstructed%sEMinus",CutName->Data()),"Reconstructed E_{T} from e^{-}");



    CreateEtaPtHisto2D(Form("EtReconstructed%sLambdaDaughters",CutName->Data()),"Reconstructed E_{T} from #Lambda Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiLambdaDaughters",CutName->Data()),"Reconstructed E_{T} from #bar{#Lambda} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sK0SDaughters",CutName->Data()),"Reconstructed E_{T} from K^{0}_{S} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sOmegaDaughters",CutName->Data()),"Reconstructed E_{T} from #Omega^{-} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiOmegaDaughters",CutName->Data()),"Reconstructed E_{T} from #Omega^{+} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sXiDaughters",CutName->Data()),"Reconstructed E_{T} from #Xi^{-} Daughters");
  CreateEtaPtHisto2D(Form("EtReconstructed%sAntiXiDaughters",CutName->Data()),"Reconstructed E_{T} from #Xi^{+} Daughters");

    CreateIntHisto1D(Form("UnidentifiedPIDs%s",CutName->Data()),"PIDs of unidentified particles", "PID", "Number of particles",9, -4,4);
    CreateHisto2D(Form("MisidentifiedPIDs%s",CutName->Data()),"PIDs of misidentified particles", "PID real","PID identified",5, -.5,4.5,5, -.5,4.5);
    CreateHisto2D(Form("dEdxAll%s",CutName->Data()),"dE/dx for all particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxPion%s",CutName->Data()),"dE/dx for #pi^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxKaon%s",CutName->Data()),"dE/dx for K^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxProton%s",CutName->Data()),"dE/dx for p(#bar{p})","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxElectron%s",CutName->Data()),"dE/dx for e^{#pm}","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
    CreateHisto2D(Form("dEdxUnidentified%s",CutName->Data()),"dE/dx for unidentified particles","momentum (GeV/c)","dE/dx",400,0.0,maxPtdEdx,200,mindEdx,maxdEdx);
  }

  CreateIntHisto1D("NEvents","Number of events","number of events","Number of events",1,0,1);

  //CreateHisto1D("MisidentifiedPIDs","PIDs for particles misidentified that are not a #pi,K,p","PID","number of entries",3000,0.5,3000.5);



//     list->Add(fHistEt);
//     TString histname = "fHistEt" + fHistogramNameSuffix;

//     fHistEt = new TH1F(histname.Data(), "Total E_{T} Distribution", 1000, 0.00, 99);
//     fHistEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
//     fHistEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");
}

