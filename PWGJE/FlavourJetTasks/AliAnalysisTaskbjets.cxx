/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

//*****************************//
#include "AliAnalysisTaskbjets.h"
//*****************************//

#include "AliAODTrack.h"
#include <math.h>
#include <AliPIDResponse.h>
#include "AliMultSelection.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
//#include "AliAnalysisTaskParticleInJet.h"
#include "AliAODMCParticle.h"
#include "AliEmcalTrackSelection.h"
#include "AliVTrack.h"


#include <vector> 
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "TChain.h"
#include <map>


class AliAnalysisTaskbjets;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskbjets) // classimp: necessary for root

AliAnalysisTaskbjets::AliAnalysisTaskbjets() : AliAnalysisTaskEmcalJet(), 
    fAOD(0), fHistPt(0), jetconrec(0), jetcongen(0), fNoJetConstituents(0) //fHistjet(),
    //fHistManager()
{
   /* SetNeedEmcalGeom(kFALSE);
    //SetOffTrigger(AliVEvent::kINT7);
    SetOffTrigger(AliVEvent::kAnyINT);
     SetVzRange(-10,10);
     //SetMakeGeneralHistograms(kTRUE);
     
     //fJetCollArray.SetOwner(kTRUE);
     
    DefineOutput(1,  AliEmcalList::Class());
 */
    // ******* default constructor, don't allocate memory here! *******
    // ******* this is used by root for IO purposes, it needs to remain empty *******
}
//_____________________________________________________________________________
AliAnalysisTaskbjets::AliAnalysisTaskbjets(const char* name) : AliAnalysisTaskEmcalJet(name, kTRUE),
    fAOD(0), fHistPt(0), fPIDResponse(0), jetconrec(0), jetcongen(0), fJetRecPt(-99), fJetArea(-99), fNoJetConstituents(0), fJetMass(-99) //, fHistManager(name) 
    
    //, fHistjet(), //, fJetCollArray()
    
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
     
     SetNeedEmcalGeom(kFALSE);
     SetOffTrigger(AliVEvent::kINT7);
     //SetOffTrigger(AliVEvent::kAnyINT);
     SetVzRange(-10,10);
     
     SetMakeGeneralHistograms(kTRUE);
     
     //fJetCollArray.SetOwner(kTRUE);
                                         
                                        // it does its work automatically
    //DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                       // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!

    DefineOutput(1,  AliEmcalList::Class());


}
//_____________________________________________________________________________
AliAnalysisTaskbjets::~AliAnalysisTaskbjets()
{
    // destructor
  /*  if(fOutput) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }*/
}

//____________________________

Bool_t AliAnalysisTaskbjets::IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax) 
   	//decides whether a particle is inside a jet cone 
   	 {
  	if(!part || !jet) AliError(Form("Particle or Jet missing: part=%p, jet=%p\n", part,jet));

  	TVector3 vecJet(jet->Px(), jet->Py(), jet->Pz());
  	TVector3 vecPart(part->Px(), part->Py(), part->Pz());
  	Double_t dR = vecJet.DeltaR(vecPart); //= sqrt(dEta*dEta +dPhi*dPhi)
  	if(dR <= dRMax){
    	//printf("IsParticleInCone:: Accepting as dR=%f, maxDR=%f\n",dR, dRMax);
   	 return kTRUE;
  	}
  	//printf("IsParticleInCone:: Rejecting as dR=%f, maxDR=%f\n",dR, dRMax);
  	return kFALSE;
	}

//_________________________________________________
void AliAnalysisTaskbjets::UserCreateOutputObjects()
{
	//Printf("Check done %i",__LINE__);
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    
  //  AllocateJetHistograms();
    
/*    TIter next(fHistManager.GetListOfHistograms());		//TIter : is the collection abstract base class
    TObject* obj = 0;
    while ((obj = next())) 
    {
    fOutput->Add(obj);
    }
 */   
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if (man) 
	{
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
	 }
    
    //fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutput->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

    // example of a histogram
    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogram
    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 1, 0, 10);
    fHistVtxPos = new TH1F("fHistVtxPos", "fHistVtxPos", 100, -10, 10);
    fHistpteta = new TH2F("fHistpteta","fHistpteta", 100, 0, 10, 100, -10, 10);
    fHistptphi = new TH2F("fHistptphi","fHistptphi", 100, 0, 10, 100, -10, 10);
    fHisttpc = new TH2F("fHisttpc", "fHisttpc", 100, 0.5, 20, 1000, 0, 1000);
    
    
    fHistjet = new TH1F("fHistjet", "fHistjet", 100, 0, 10); 		//histogram for jet-number
    
    
    fHistjetPt = new TH1F("fHistjetPt", "fHistjetPt", 100, 0, 10);		//histogram for jet-pt
    
    //fHisttest = new TH1F("fHisttest", "fHisttest", 100, 0, 10);
    
    fOutput->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
                                        
    fOutput->Add(fHistNEvents);
    fOutput->Add(fHistVtxPos);
    
    //fOutput->Add(fHistjet);
    
    fOutput->Add(fHistpteta);
    fHistpteta->SetXTitle("pT");
    fHistpteta->SetYTitle("#eta");
    
    fOutput->Add(fHistptphi);
    fHistptphi->SetXTitle("pT");
    fHistptphi->SetYTitle("#phi");
    
    fOutput->Add(fHisttpc);
    
    fHisttpc->SetXTitle("p/z (GeV/c)");
    fHisttpc->SetYTitle("TPC dE/dx (arb. units)");
    //fHisttpc->Draw("colz");
    
    fOutput->Add(fHistjet);
    
    fOutput->Add(fHistjetPt);
    
    
    PostData(1, fOutput);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
/*
void AliAnalysisTaskbjets::AllocateJetHistograms()
{
	TString histname;
	TString histtitle;
	TString groupname;
	AliJetContainer* jetCont = 0;
	TIter next(&fJetCollArray);
	while ((jetCont = static_cast<AliJetContainer*>(next())))
	
	//while ((AliJetContainer* jetCont = static_cast<AliJetContainer*>(next())))
	{
		groupname = jetCont->GetName();    
		//Protect against creating the histograms twice
		if (fHistManager.FindObject(groupname))
		{
			AliWarning(TString::Format("%s: Found groupname % in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
			continue;
		}

		fHistManager.CreateHistoGroup(groupname);
		for(Int_t cent = 0; cent < fNcentBins; cent++)
		{
		histname = TString::Format("%s/histJetPt %d" , groupname.Data() , cent);
		histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
		fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);
		
		histname= TString::Format("%s/histNJets %d", groupname.Data(), cent);
		histtitle = TString::Format("%s;number of jets; events", histname.Data());
		if (fForceBeamType != kpp)
			{
			fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
			}
		else 	{
			fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
			} 
		}
	}
}
*/

/*
Bool_t AliAnalysisTaskbjets::FillHistograms()
{
DoJetLoop();
 return kTRUE;
}

void AliAnalysisTaskbjets::DoJetLoop()
{
	TString histname;
	TString groupname;
	AliJetContainer* jetCont = 0;
	TIter next(&fJetCollArray);
	while ((jetCont = static_cast<AliJetContainer*>(next())))
	
	//while ((AliJetContainer* jetCont = static_cast<AliJetContainer*>(next())))
	
	{
		groupname = jetCont->GetName();
		UInt_t count = 0;
		for(auto jet : jetCont->accepted()) 
		{
			if (!jet) continue;
			count++;
		
			
//		histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
//  		      	fHistManager.FillTH1(histname, jet->Pt());
		}
//		histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
 //      	fHistManager.FillTH1(histname, count);
       	
       }
}
*/
//_____________________________________________________________________________

void AliAnalysisTaskbjets::ExecOnce()
{
	AliAnalysisTaskEmcalJet::ExecOnce();
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskbjets::Run()
{
	
//	Printf("Check done %i",__LINE__);
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) {return kFALSE;  }
    
                                     // if the pointer to the event is empty (getting it failed) skip this event
    else {
    
    
    fHistNEvents->Fill(1);
    
    //IsMCJetPartonFast( jet, radius, is_udg);
    AliAODVertex *vertex  = fAOD->GetPrimaryVertex(); 
    
    
    fHistVtxPos->Fill(vertex->GetZ());
    
    
    jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(0));		//Particle-level
    jetcongen->ResetCurrentID();
    
    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(1));		//Detector-level
    jetconrec->ResetCurrentID();
   
    //printf("njets= %i\n", jetconrec->GetNJets());
    
    printf("njets= %i\n", jetcongen->GetNJets());
    
    fHistjet->Fill(1);
    
    
    
    
    
   // Look at the MC Jet Container
   
   /* 
   	* Looking into the jet cone
   	* Plot its pt-spectra
   	* Find out its angular distance wrt to the jet axis and check if this distance is < jet-radius
	
   */	
   /* GetJetContainer(fJetContainerNameMC.Data())->ResetCurrentID();
    AliEmcalJet *mcjet = GetJetContainer(fJetContainerNameMC.Data())->GetNextAcceptJet();
    do{
 	      for(int itrk = 0; itrk < mcjet->GetNumberOfTracks(); itrk++){
    	         mctrack = mcjet->TrackAt(itrk, GetParticleContainer(0)->GetArray());
  	         if(!AcceptParticle(mctrack)) continue;
  	         fHistMgr->FillTH2("hMCjetTrack", mcjet->Pt(), mctrack->Pt());
  	       }
  	     } while ((mcjet = GetJetContainer(fJetContainerNameMC.Data())->GetNextAcceptJet()));
  	  }
   */
   
   AliEmcalJet * jetrec = 0;
   AliEmcalJet * jetmatched = 0;
   Bool_t is_udg;
   Double_t radius = 0.4;
      
   while ((jetrec = jetcongen->GetNextAcceptJet())) 		//start of jet loop
	{	
	
		
		//if(!jetrec) continue;
		//if(jetrec->GetNumberOfTracks()==0) continue;
		
		if (!jetcongen) return kFALSE;
		IsMCJetPartonFast(jetrec, radius, is_udg);
		
//		fHistManager.FillTH1("fh1dNoJetsPerEvent",jetconrec->GetNJets(),1);
		
//		fJetRecPt = jetrec->Pt();	//jet pt-distribution
		fJetMass  = jetrec->M();	// jet-mass
		fJetArea  = jetrec->Area();	// jet-area
		
		printf("Generated: fJetRecPt=%f, fJetMass=%f\n",fJetRecPt, fJetMass);
		printf("Generated: fJetArea=%f\n", fJetArea);
	
    		fHistjetPt->Fill(jetrec->Pt());
    	}
    	
     	
	    
     Float_t centrality(0);
     
     if( (0< centrality) && (centrality < 10) )
     {
    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");      
    }
    
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    for(Int_t i(0); i < iTracks; i++) 
    {                 					// loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                            // if we failed, skip this track
        if(!track->TestFilterBit(1)) continue;

        fHistPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram
        //fHistPt->SetXTitle("");
        //fHistPt->SetYTitle();
        
        Double_t px = track->Px();
        Double_t py = track->Py();
        Double_t pz = track->Pz();
        
        Double_t eta = track->Eta();
        Double_t phi = track->Phi();
        Double_t pt = track->Pt();
        
        if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )
        {
          fHistpteta->Fill(pt, eta);
          
    	  //fHistpteta->Draw("colz");
     	  
     	  
     	  fHistptphi->Fill(pt, phi);
     	  
     	  //fHistptphi->Draw("colz");
     	  
   	 } 
   	 
   
   	 
   	 
           //Double_t kaonSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
	   Double_t pionSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
           //Double_t protonSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

	 if (std::abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 ) 
	 {
  	 //jippy, i'm a pion
	 };

	fHisttpc->Fill(track->P(), track->GetTPCsignal());
	
	
	
                                        // continue until all the tracks are processed
    PostData(1, fOutput);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
    }                                                    // it to a file//
} 
return kTRUE;
}

//Printf("Check done %i",__LINE__);
//*****Function for calculating the angular distance of the track wrt to the jet axis = delta_R = RADIAL SHAPE DISTRIBUTION *****


/*AliEmcalJet * jetrec = 0;
AliEmcalJet * jetmatched = 0;
Int_t pdg=0;*/

//Printf("Check done %i",__LINE__);
/*
Int_t  AliAnalysisTaskbjets::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
{	Printf("LINE NUMBER %d",__LINE__);

	AliVParticle* vp = 0x0;
	    
	if(!jet) return 0;
            if(!(jet->GetNumberOfTracks()>fNoJetConstituents))
            { 
	      Printf("LINE NUMBER %d",__LINE__);            
            
              printf("Throwing away jets with only too few constituent!\n");
              return 0;
            }
         
         
        //while ((jetrec = jetconrec->GetNextAcceptJet()))
        
        while ((jetrec = jetcongen->GetNextAcceptJet()))
        { //start jetloop
        
		Printf("LINE NUMBER %d",__LINE__);
		        
      		if(!jetrec) continue;
	        if(jetrec->GetNumberOfTracks()==0) continue;
	        
	        jet->Print();
	        
	        for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) 
	{//start trackloop jet
	      
	      	Printf("LINE NUMBER %d",__LINE__);
	      	
             	 vp = reinterpret_cast<AliVParticle*>(jet->TrackAt(i)); 
                if (!vp)
                   { 
                   Printf("LINE NUMBER %d",__LINE__);
                   
                  	AliError("AliVParticle associated to constituent not found\n");
                  	continue;
                   }
	}
	Printf("LINE NUMBER %d",__LINE__);
	             
        /*        AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp);
      		
      		if (!vtrack) 
      		{
        		AliError(Form("Could not receive track%d\n", i));
        		continue;
      		}
      			AliAODTrack *trackV = dynamic_cast<AliAODTrack*>(vtrack);
  		
  		Printf("Check done %i",__LINE__);
  		AliAODMCParticle* part = static_cast<AliAODMCParticle*>(vp);
  		{ Printf("Check done %i",__LINE__);
  		if(!part) 
  		{
  			AliError("Finding no Part!\n");
  			return 0;
  		}
  		
  	
		if(!part->IsPrimary())
		{
		continue;
		}
		
		Double_t etajet = part->Eta();		//jet-axis
		Double_t phijet = part->Phi();		
		
		Double_t etatrack = trackV->Eta();		//track-axis
		Double_t phitrack = trackV->Phi();
		
		Int_t pdg = (abs(part->PdgCode()));
		Double_t deta = etatrack - etajet;
		Double_t dphi = phitrack - phijet;
               
                dphi = TVector2::Phi_mpi_pi(dphi);  		//Returns phi angle in the interval [-PI,PI)
                                
                Double_t d = sqrt(deta * deta + dphi * dphi);
                printf("d = %f");
                }
               }
        
   
   Printf("Check done %i",__LINE__);	
   }         
Printf("LINE NUMBER %d",__LINE__);
return 0;
	
	}
}*/
Int_t  AliAnalysisTaskbjets::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
{	

	AliVParticle* vp = 0x0;
	    
	if(!jet) return 0;
            if(!(jet->GetNumberOfTracks()>fNoJetConstituents))
            {
              Printf("Throwing away jets with only too few constituent!\n");
              return 0;
            }
         
            for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) 
           {//start trackloop
              vp = static_cast<AliVParticle*>(jet->Track(i));
              if (!vp){
                AliError("AliVParticle associated to constituent not found");
                continue;
              }
             /* AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp);
              if (!vtrack) {
                AliError(Form("Could not receive track%d\n", i));
                continue;
              }
              AliAODTrack *trackV = dynamic_cast<AliAODTrack*>(vtrack); */
              
              Printf("Track pt=%f\n", vp->Pt());
              
              Printf("LINE NUMBER %d",__LINE__);
              
              AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(vp);
  		 Printf("LINE NUMBER %d",__LINE__);
  		
  	/*	if(!part) {
  		 Printf("LINE NUMBER %d",__LINE__);
  		// AliError("Finding no Part!\n");
  		//	return 0;
  		continue;
  		}
  		
  		Printf("LINE NUMBER %d",__LINE__);
  	
		if(!part->IsPrimary())
		{ Printf("LINE NUMBER %d",__LINE__);
		continue;
		}
	*/	
		Printf("LINE NUMBER %d",__LINE__);
		Double_t etajet = part->Eta();		//jet-axis
		Double_t phijet = part->Phi();		
		
		//Double_t etatrack = trackV->Eta();		//track-axis
		//Double_t phitrack = trackV->Phi();
		
		Double_t etatrack = vp->Eta();		//track-axis
		Double_t phitrack = vp->Phi();
		
		Int_t pdg = (abs(part->GetPdgCode()));
		Double_t deta = etatrack - etajet;
		Double_t dphi = phitrack - phijet;
		
		Printf("pdg = %f\n", pdg);
		              
                dphi = TVector2::Phi_mpi_pi(dphi);  		//Returns phi angle in the interval [-PI,PI)
                                
                Double_t d = sqrt(deta * deta + dphi * dphi);
                printf("d = %f\n", d);
                Printf("LINE NUMBER %d",__LINE__);
                }
 return 0;         
 	}
	
 
//_____________________________________________________________________________

/*void AliAnalysisTaskbjets::Run()
{
 return kTRUE;
}*/
//_____________________________________________________________________________
void AliAnalysisTaskbjets::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
