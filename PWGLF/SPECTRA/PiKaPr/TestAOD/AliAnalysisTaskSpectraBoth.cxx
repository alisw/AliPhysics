/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//         AliAnalysisTaskSpectraBoth class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskSpectraBoth.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraBothHistoManager.h"
#include "AliSpectraBothTrackCuts.h"
#include "AliSpectraBothEventCuts.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "AliSpectraBothPID.h"
#include "AliGenEventHeader.h"	
#include <TMCProcess.h>

#include <iostream>




using namespace AliSpectraNameSpaceBoth;
using namespace std;

ClassImp(AliAnalysisTaskSpectraBoth)

//________________________________________________________________________
AliAnalysisTaskSpectraBoth::AliAnalysisTaskSpectraBoth(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0),  fPID(0), fIsMC(0), fNRebin(0),fUseMinSigma(0),fCuts(0),fdotheMCLoopAfterEventCuts(0),fmakePIDQAhisto(1),fMotherWDPDGcode(-1),fUseEtaCut(kFALSE),fIncludecorrectlyidentifiedinMCtemplates(kFALSE)

{
  // Default constructor
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, AliSpectraBothHistoManager::Class());
  DefineOutput(2, AliSpectraBothEventCuts::Class());
  DefineOutput(3, AliSpectraBothTrackCuts::Class());
  DefineOutput(4, AliSpectraBothPID::Class());
  fNRebin=0;
  
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskSpectraBoth::UserCreateOutputObjects()
{
  // create output objects
  fHistMan = new AliSpectraBothHistoManager("SpectraHistos",fNRebin,fmakePIDQAhisto);
  fHistMan->SetIncludecorrectlyidentifiedinMCtemplates(fIncludecorrectlyidentifiedinMCtemplates);	

  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  if (!fPID)       AliFatal("PID object should be set in the steering macro");
  fTrackCuts->SetAliESDtrackCuts(fCuts);
  fEventCuts->InitHisto();
  fTrackCuts->InitHisto();
  if(fTrackCuts->GetYMax()<fTrackCuts->GetYMin()) 
  	fUseEtaCut=kTRUE;	
  Printf(" eta cut %d  Will be used lack of Y cut",fUseEtaCut);


  PostData(1, fHistMan  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fPID      );
}
//________________________________________________________________________
void AliAnalysisTaskSpectraBoth::UserExec(Option_t *)
{
  // main event loop
	Int_t ifAODEvent=AliSpectraBothTrackCuts::kotherobject;
	fAOD = dynamic_cast<AliVEvent*>(InputEvent());
  //	AliESDEvent* esdevent=0x0;
 // 	AliAODEvent* aodevent=0x0;
 
   	TString nameoftrack(fAOD->ClassName()); 
    	if(!nameoftrack.CompareTo("AliESDEvent"))
    	{
		ifAODEvent=AliSpectraBothTrackCuts::kESDobject;
		//esdevent=dynamic_cast<AliESDEvent*>(fAOD);
	}
	else if(!nameoftrack.CompareTo("AliAODEvent"))
	{
		ifAODEvent=AliSpectraBothTrackCuts::kAODobject;
		//aodevent=dynamic_cast<AliAODEvent*>(fAOD);
	}
	else
		AliFatal("Not processing AODs or ESDS") ;
	if(fIsMC)
  	{		
  		if(!fEventCuts->CheckMCProcessType(MCEvent()))
			return ;		
  	}

	if(fdotheMCLoopAfterEventCuts)
  		if(!fEventCuts->IsSelected(fAOD,fTrackCuts,fIsMC,-100,fHistMan->GetEventStatHist()))
			return;//event selection
  	TClonesArray *arrayMC = 0;
  	Int_t npar=0;
  	AliStack* stack=0x0;
	 Double_t mcZ=-100;
 	if (fIsMC)
  	{
		TArrayF mcVertex(3);
  		mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
		AliMCEvent* mcEvent=(AliMCEvent*)MCEvent();
		if (!mcEvent) 
		{
			AliFatal("Error: MC particles branch not found!\n");
		}
		AliHeader* header = mcEvent->Header();
    		if (!header) 
		{
      			AliDebug(AliLog::kError, "Header not available");
      			return;
    		}
	
		AliGenEventHeader* genHeader = header->GenEventHeader();
    		if(genHeader)
		{
			genHeader->PrimaryVertex(mcVertex);
  			mcZ=mcVertex[2];
		}
		  if(ifAODEvent==AliSpectraBothTrackCuts::kAODobject)
		  {
			  arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
			  if (!arrayMC) 
			  {
					AliFatal("Error: MC particles branch not found!\n");
			  }
			  Int_t nMC = arrayMC->GetEntries();
			  for (Int_t iMC = 0; iMC < nMC; iMC++)
			  {
				  AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
				  if(!partMC->Charge()) continue;//Skip neutrals
				  //if(partMC->Eta() > fTrackCuts->GetEtaMin() && partMC->Eta() < fTrackCuts->GetEtaMax()){//charged hadron are filled inside the eta acceptance
				  //Printf("%f     %f-%f",partMC->Eta(),fTrackCuts->GetEtaMin(),fTrackCuts->GetEtaMax());
				  if(partMC->Eta() > fTrackCuts->GetEtaMin() && partMC->Eta() < fTrackCuts->GetEtaMax())
				  {
						fHistMan->GetPtHistogram(kHistPtGen)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary());
				  }
				  else 
				  {
					if(fUseEtaCut)
						continue;
				  }						 				 
				  //rapidity cut
				  if(!fUseEtaCut)
				  {		
				 	 if(partMC->Y() > fTrackCuts->GetYMax()|| partMC->Y() < fTrackCuts->GetYMin() ) 
						continue;
				  }		
				  if(partMC->IsPhysicalPrimary())
				 	 npar++;    
				  // check for true PID + and fill P_t histos 
				  Int_t chargetmp = partMC->Charge() > 0 ? kChPos : kChNeg ;
				  Int_t id = fPID->GetParticleSpecie(partMC);
				  if(id != kSpUndefined) 
				  {
					fHistMan->GetHistogram2D(kHistPtGenTruePrimary,id,chargetmp)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary());
				  }
			  }
		  }
		  if(ifAODEvent==AliSpectraBothTrackCuts::kESDobject)
		  {
			stack = mcEvent->Stack();
			Int_t nMC = stack->GetNtrack();
			for (Int_t iMC = 0; iMC < nMC; iMC++)
			{
				 
				TParticle *partMC = stack->Particle(iMC);
				
				if(!partMC)	
					continue;	
	
				if(!partMC->GetPDG(0))
					continue;
				if(TMath::Abs(partMC->GetPDG(0)->Charge()/3.0)<0.01) 
					continue;//Skip neutrals
			 	if(partMC->Eta() > fTrackCuts->GetEtaMin() && partMC->Eta() < fTrackCuts->GetEtaMax())
				{
					fHistMan->GetPtHistogram(kHistPtGen)->Fill(partMC->Pt(),stack->IsPhysicalPrimary(iMC));
				}
				else 
				{
					if(fUseEtaCut)
						continue;
				}
				if(!fUseEtaCut)
				{ 		
					if(partMC->Y()>fTrackCuts->GetYMax() ||partMC->Y()< fTrackCuts->GetYMin()  ) 
						continue;
				}
				if(stack->IsPhysicalPrimary(iMC))
					 npar++;    
				  // check for true PID + and fill P_t histos 
				Int_t chargetmp = partMC->GetPDG(0)->Charge()/3.0 > 0 ? kChPos : kChNeg ;
				Int_t id = fPID->GetParticleSpecie(partMC);
				if(id != kSpUndefined) 
				{
					fHistMan->GetHistogram2D(kHistPtGenTruePrimary,id,chargetmp)->Fill(partMC->Pt(),stack->IsPhysicalPrimary(iMC));
				}
			  }
		  }
	}
	if(!fdotheMCLoopAfterEventCuts)
  		if(!fEventCuts->IsSelected(fAOD,fTrackCuts,fIsMC,mcZ,fHistMan->GetEventStatHist()))
			return;//event selection
  	//main loop on tracks
	Int_t ntracks=0;
  	//cout<<fAOD->GetNumberOfTracks()<<endl;
  	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
	{
		AliVTrack* track = dynamic_cast<AliVTrack*>(fAOD->GetTrack(iTracks));
      		AliAODTrack* aodtrack=0;
  		AliESDtrack* esdtrack=0;
  		Float_t dca=-999.;
		Float_t dcaz=-999.;
		Short_t ncls=-1; 
		Float_t chi2perndf=-1.0;
  		if(ifAODEvent==AliSpectraBothTrackCuts::kESDobject)
  		{
			esdtrack=dynamic_cast<AliESDtrack*>(track);
			if(!esdtrack)
				continue;
			esdtrack->GetImpactParameters(dca,dcaz);
			ncls=esdtrack->GetTPCNcls();
			if ( ncls > 5) 
       				chi2perndf=(esdtrack->GetTPCchi2()/Float_t(ncls - 5));
  
			else 
       				chi2perndf=-1;
    			
  		}
  		else if (ifAODEvent==AliSpectraBothTrackCuts::kAODobject)
  		{
			aodtrack=dynamic_cast<AliAODTrack*>(track);
			if(!aodtrack)
				continue;		
			dca=aodtrack->DCA();
			dcaz=aodtrack->ZAtDCA();
			ncls=aodtrack->GetTPCNcls();
			chi2perndf=aodtrack->Chi2perNDF();
  		}
  		else
			continue;
    		if (!fTrackCuts->IsSelected(track,kTRUE)) 
			continue;	
			
    		ntracks++;
	
    		
		//calculate DCA for AOD track
    		if(dca==-999.)
		{// track->DCA() does not work in old AOD production
     	 		Double_t d[2], covd[3];
      			AliVTrack* track_clone=(AliVTrack*)track->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
      			Bool_t isDCA = track_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),9999.,d,covd);
      			delete track_clone;
      			if(!isDCA)
				d[0]=-999.;
      			dca=d[0];
			dcaz=d[1];	
    		}
	
     		fHistMan->GetPtHistogram(kHistPtRec)->Fill(track->Pt(),dca);  // PT histo
    		// get identity and charge
    		Bool_t rec[3]={false,false,false};
		Bool_t sel[3]={false,false,false};
    		Int_t idRec  = fPID->GetParticleSpecie(fHistMan,track, fTrackCuts,rec);
		Int_t idGen     =kSpUndefined;
		Bool_t isPrimary           = kFALSE;
		Bool_t isSecondaryMaterial = kFALSE; 
		Bool_t isSecondaryWeak     = kFALSE; 

		for(int irec=kSpPion;irec<kNSpecies;irec++)
    		{
   
			if(fUseMinSigma)
			{
				if(irec>kSpPion)
					break;
			}
			else
			{	
				if(!rec[irec]) 
					idRec = kSpUndefined;
				else	
					idRec=irec;
			}		
   
			Int_t charge = track->Charge() > 0 ? kChPos : kChNeg;	
			// Fill histograms, only if inside y and nsigma acceptance
			if(idRec != kSpUndefined && fTrackCuts->CheckYCut ((BothParticleSpecies_t)idRec))
			{
				fHistMan->GetHistogram2D(kHistPtRecSigma,idRec,charge)->Fill(track->Pt(),dca);
				if(fTrackCuts->GetMakeQAhisto())
				{ 
					fTrackCuts->GetHistoDCAzQA()->Fill(idRec,track->Pt(),dcaz);
					fTrackCuts->GetHistoNclustersQA()->Fill(idRec,track->Pt(),ncls);
					fTrackCuts->GetHistochi2perNDFQA()->Fill(idRec,track->Pt(),chi2perndf);
				}
				sel[idRec]=true;
			}
			//can't put a continue because we still have to fill allcharged primaries, done later
		
			/* MC Part */
			if (arrayMC||stack) 
			{
				//Bool_t isPrimary           = kFALSE;
			 //	Bool_t isSecondaryMaterial = kFALSE; 
			 //	Bool_t isSecondaryWeak     = kFALSE; 
			 //	Int_t idGen     =kSpUndefined;
			 	Int_t pdgcode=0;
				Int_t motherpdg=-1;
				Int_t chargeMC=-2;
				if (ifAODEvent==AliSpectraBothTrackCuts::kAODobject)
				{
					AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
				  	if (!partMC) 
					{ 
						AliError("Cannot get MC particle");
						continue; 
				  	}
				  	// Check if it is primary, secondary from material or secondary from weak decay
				  	isPrimary           = partMC->IsPhysicalPrimary();
					isSecondaryWeak     = partMC->IsSecondaryFromWeakDecay();
					isSecondaryMaterial      = partMC->IsSecondaryFromMaterial();
					//cout<<"AOD tagging "<<isPrimary<<" "<<isSecondaryWeak<<isSecondaryMaterial<<" "<<partMC->GetMCProcessCode()<<endl;

				  	if(!isPrimary&&!isSecondaryWeak&&!isSecondaryMaterial)//old tagging for old AODs 
				  	{
						AliError("old tagging");
						Int_t mfl=-999,codemoth=-999;
						Int_t indexMoth=partMC->GetMother(); // FIXME ignore fakes? TO BE CHECKED, on ESD is GetFirstMother()
						if(indexMoth>=0)
						{//is not fake
					  		AliAODMCParticle* moth = (AliAODMCParticle*) arrayMC->At(indexMoth);
					  		codemoth = TMath::Abs(moth->GetPdgCode());
					  		mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
						}
						//Int_t uniqueID = partMC->GetUniqueID();
						//cout<<"uniqueID: "<<partMC->GetUniqueID()<<"       "<<kPDecay<<endl;
						//cout<<"status: "<<partMC->GetStatus()<<"       "<<kPDecay<<endl;
						// if(uniqueID == kPDecay)Printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
						if(mfl==3) 
							isSecondaryWeak     = kTRUE; // add if(partMC->GetStatus() & kPDecay)? FIXME
						else       
							isSecondaryMaterial = kTRUE;
				  	}
					//cout<<"AOD 2 tagging "<<isPrimary<<" "<<isSecondaryWeak<<isSecondaryMaterial<<" "<<partMC->GetMCProcessCode()<<endl;
					if(isSecondaryWeak)
					{	
						Int_t indexMoth=partMC->GetMother(); // FIXME ignore fakes? TO BE CHECKED, on ESD is GetFirstMother()
						if(indexMoth>=0)
						{
					  		AliAODMCParticle* moth = (AliAODMCParticle*) arrayMC->At(indexMoth);
							if(moth)
					  			motherpdg=TMath::Abs(moth->GetPdgCode());
						}
					}

				  	idGen     = fPID->GetParticleSpecie(partMC);
				  	pdgcode=partMC->GetPdgCode(); 
					chargeMC=partMC->Charge() > 0 ? kChPos : kChNeg ;
				}
				else if (ifAODEvent==AliSpectraBothTrackCuts::kESDobject)
				{
					TParticle *partMC =stack->Particle(TMath::Abs(track->GetLabel()));
					if (!partMC) 
					{ 
						AliError("Cannot get MC particle");
						continue; 
				  	}
				  	isPrimary           = stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()));
					isSecondaryWeak     = stack->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()));
					isSecondaryMaterial      = stack->IsSecondaryFromMaterial(TMath::Abs(track->GetLabel()));
					//cout<<"ESD tagging "<<isPrimary<<" "<<isSecondaryWeak<<isSecondaryMaterial<<endl;
					
					if(isSecondaryWeak)	
					{

						TParticle* moth=stack->Particle(TMath::Abs(partMC->GetFirstMother()));
						if(moth)
							 motherpdg = TMath::Abs(moth->GetPdgCode());

					}
					

				   	idGen     = fPID->GetParticleSpecie(partMC);
				   	pdgcode=partMC->GetPdgCode(); 
					chargeMC = partMC->GetPDG(0)->Charge()/3.0 > 0 ? kChPos : kChNeg ;


				}
				else
					return;
			
			//	  cout<<isPrimary<<" "<<isSecondaryWeak<<" "<<isSecondaryMaterial<<endl;
			//	  cout<<" functions "<<partMC->IsPhysicalPrimary()<<" "<<partMC->IsSecondaryFromWeakDecay()<<" "<<partMC->IsSecondaryFromMaterial()<<endl;
		  
		  		if (isPrimary&&irec==kSpPion)
				{
					fHistMan->GetPtHistogram(kHistPtRecPrimaryAll)->Fill(track->Pt(),dca);  // PT histo of reconstrutsed primaries in defined eta
					if(chargeMC!=charge)
						fHistMan->GetPtHistogram("hHistDoubleCounts")->Fill(track->Pt(),4);

				}	

				if(track->Pt()>fTrackCuts->GetPtTOFMatching(irec)&&(!fTrackCuts->CheckTOFMatchingParticleType(irec)))
					continue;

				//in case of pt depended TOF cut we have to remove particles with pt above their TOF cut but below max TOF cut

	
				if(fUseMinSigma)
				{
					if(idRec == kSpUndefined)
						continue;
					if(!fTrackCuts->CheckYCut ((BothParticleSpecies_t)idRec)) 
						continue;
				}			  
				else
				{
					if(!fTrackCuts->CheckYCut ((BothParticleSpecies_t)irec)) 
						continue;

				}
				
		  	// rapidity cut (reconstructed pt and identity)
		 		// if(!fTrackCuts->CheckYCut ((BothParticleSpecies_t)idRec)) continue;
		  
		  // Get true ID
		  		 	
		  
		 		 if ((idRec == idGen)&&(idGen != kSpUndefined)) 
					fHistMan->GetHistogram2D(kHistPtRecTrue,  idGen, charge)->Fill(track->Pt(),dca); 
		  
		  		if (isPrimary) 
				{
					if(idRec!= kSpUndefined) // any primary 
						fHistMan->GetHistogram2D(kHistPtRecSigmaPrimary, idRec, charge)->Fill(track->Pt(),dca); 
					//if(idGen != kSpUndefined) 
					//{
					 if((idGen != kSpUndefined) &&(irec == idGen)) // genereated primary but does not have to be selected useless for min sigma 
		  				fHistMan->GetHistogram2D(kHistPtRecPrimary,      idGen, charge)->Fill(track->Pt(),dca);
		  			if ((idGen != kSpUndefined) &&(idRec == idGen))  // genereated and selected correctly 
						fHistMan->GetHistogram2D(kHistPtRecTruePrimary,  idGen, charge)->Fill(track->Pt(),dca); 
					//}
		  		}
		 		 //25th Apr - Muons are added to Pions -- FIXME
		 		if ( pdgcode == 13 && idRec == kSpPion) 
				{ 
					fHistMan->GetPtHistogram(kHistPtRecTrueMuonPlus)->Fill(track->Pt(),dca); 
					if(isPrimary)
		  				fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonPlus)->Fill(track->Pt(),dca); 
		  		}
		  		if ( pdgcode == -13 && idRec == kSpPion) 
				{ 
					fHistMan->GetPtHistogram(kHistPtRecTrueMuonMinus)->Fill(track->Pt(),dca); 
					if (isPrimary) 
					{
		  				fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonMinus)->Fill(track->Pt(),dca); 
					}
		  		}
		  		
				if(idRec == kSpUndefined)
					continue;

				
		 		 //here we can use idGen in case of fit approach
		  
		  		// Fill secondaries

				
				if(fMotherWDPDGcode>0) // if the Mother pdg is set we undo the Secondary flag in case it dose not match
				{
					if(motherpdg!=fMotherWDPDGcode)
						isSecondaryWeak=kFALSE;
						
				}
				if(fHistMan->GetIncludecorrectlyidentifiedinMCtemplates())// we have to check if genereted is the same as reconstructed  
				{
					if(idRec!=idGen)
						continue;
				}
		  		if(isSecondaryWeak)
					fHistMan->GetHistogram2D(kHistPtRecSigmaSecondaryWeakDecay, idRec, charge)->Fill(track->Pt(),dca);
		  		if(isSecondaryMaterial)  
					fHistMan->GetHistogram2D(kHistPtRecSigmaSecondaryMaterial , idRec, charge)->Fill(track->Pt(),dca);
		  
			}//end if(arrayMC)
		}
		if(sel[0]&&sel[1]&&sel[2])//pi+k+p
			fHistMan->GetPtHistogram("hHistDoubleCounts")->Fill(track->Pt(),0);
		else if(sel[0]&&sel[1]) //pi+k
			fHistMan->GetPtHistogram("hHistDoubleCounts")->Fill(track->Pt(),1);
		else if(sel[0]&&sel[2]) //pi+k
			fHistMan->GetPtHistogram("hHistDoubleCounts")->Fill(track->Pt(),2);
		else if(sel[1]&&sel[2]) //p+k
			fHistMan->GetPtHistogram("hHistDoubleCounts")->Fill(track->Pt(),3);
		if(fmakePIDQAhisto)
    			fPID->FillQAHistos(fHistMan, track, fTrackCuts,idGen);
	
	
  	} // end loop on tracks

 // cout<< ntracks<<endl;
  fHistMan->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->Fill(npar,ntracks);
    fPID->SetoldT0();
  PostData(1, fHistMan  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fPID      );
}

//_________________________________________________________________
void   AliAnalysisTaskSpectraBoth::Terminate(Option_t *)
{
  // Terminate
}
