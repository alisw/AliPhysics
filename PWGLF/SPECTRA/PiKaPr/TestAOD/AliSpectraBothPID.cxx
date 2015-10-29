#include "AliSpectraBothPID.h"
#include "AliAODEvent.h"      
#include "TH1F.h"             
#include "TH2F.h"             
#include "TList.h"            
#include "AliAODTrack.h"      
#include "AliAODMCParticle.h" 
#include "AliPIDResponse.h"   
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliSpectraBothTrackCuts.h"

ClassImp(AliSpectraBothPID)

AliSpectraBothPID::AliSpectraBothPID() : TNamed("PID", "PID object"), fPIDType(kNSigmacircleTPCTOF), fNSigmaPID(3), fPIDResponse(0),fshiftTPC(0),fshiftTOF(0) 
{

}

AliSpectraBothPID::AliSpectraBothPID(BothPIDType_t pidType) : TNamed("PID", "PID object"), fPIDType(pidType), fNSigmaPID(3), fPIDResponse(0), fshiftTPC(0),fshiftTOF(0)
{



}



void AliSpectraBothPID::FillQAHistos(AliSpectraBothHistoManager * hman, AliVTrack * track, AliSpectraBothTrackCuts * trackCuts) 
{

  // fill a bunch of QA histos


  // Get PID response object, if needed
  	if(!fPIDResponse) 
  	{
    		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   		 AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    		fPIDResponse = inputHandler->GetPIDResponse();
  	}
	Bool_t ycut[3]={false,false,false};
   	ycut[0]=trackCuts->CheckYCut ((BothParticleSpecies_t) kSpPion);
  	ycut[1]=trackCuts->CheckYCut ((BothParticleSpecies_t) kSpKaon);
	ycut[2]=trackCuts->CheckYCut ((BothParticleSpecies_t) kSpProton);
 
 
  //Response
  	AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(track);
  
  	hman->GetPIDHistogram(kHistPIDTPC)->Fill(track->GetTPCmomentum(), track->GetTPCsignal()*track->Charge()); // PID histo
  	hman->GetPIDHistogram(kHistPIDTPCPion)->Fill(track->GetTPCmomentum(),fPIDResponse->GetTPCResponse().GetExpectedSignal(track->GetTPCmomentum(),AliPID::kPion)*track->Charge()); // Expected PIDPion histo
  	hman->GetPIDHistogram(kHistPIDTPCKaon)->Fill(track->GetTPCmomentum(),fPIDResponse->GetTPCResponse().GetExpectedSignal(track->GetTPCmomentum(),AliPID::kKaon)*track->Charge()); // Expected PIDKaon histo
  	hman->GetPIDHistogram(kHistPIDTPCProton)->Fill(track->GetTPCmomentum(),fPIDResponse->GetTPCResponse().GetExpectedSignal(track->GetTPCmomentum(),AliPID::kProton)*track->Charge()); // Expected PIDProton histo
  
  	Double_t nsigmaTPCkProton = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton)+fshiftTPC;
  	Double_t nsigmaTPCkKaon = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)+fshiftTPC; 
  	Double_t nsigmaTPCkPion = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)+fshiftTPC; 
  	Double_t nsigmaTOFkProton=100.0,nsigmaTOFkKaon=100.0,nsigmaTOFkPion=100.0;

  	Double_t nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
  	Double_t nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
  	Double_t nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
	if(fPIDType==kNSigmaTOF)
	{
		nsigmaTPCTOFkProton = 100.0;
  		nsigmaTPCTOFkKaon   = 100.0;
  		nsigmaTPCTOFkPion   =100.0;

	}
	Float_t tmppt=TMath::Min(trackCuts->GetPtTOFMatchingPion(),TMath::Min(trackCuts->GetPtTOFMatchingKaon(),trackCuts->GetPtTOFMatchingProton()));
 	if(track->Pt()>=tmppt&&trackCuts->CheckTOFMatchingParticleType(kSpProton)&&trackCuts->CheckTOFMatchingParticleType(kSpKaon)&&trackCuts->CheckTOFMatchingParticleType(kSpPion))
  	{
    
   		 hman->GetPIDHistogram(kHistPIDTOF)->Fill(track->P(),(track->GetTOFsignal()/100)*track->Charge()); // PID histo
    
   		 nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton)+fshiftTOF;
  		 nsigmaTOFkKaon = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)+fshiftTOF; 
   		 nsigmaTOFkPion = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)+fshiftTOF; 
    
    //TOF	 	 
		if(ycut[0])
    		{
    			hman->GetPtHistogram(kHistNSigPionTOF)->Fill(track->P()*track->Charge(),nsigmaTOFkPion );
    			hman->GetPtHistogram(kHistNSigPionPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkPion );
    		}
    		if(ycut[1])
  	 	{	
    			hman->GetPtHistogram(kHistNSigKaonTOF)->Fill(track->P()*track->Charge(),nsigmaTOFkKaon);
   	 		hman->GetPtHistogram(kHistNSigKaonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkKaon );	 
   		}
		if(ycut[2])
		{					
   			hman->GetPtHistogram(kHistNSigProtonTOF)->Fill(track->P()*track->Charge(),nsigmaTOFkProton ); 
    			hman->GetPtHistogram(kHistNSigProtonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkProton );
		}
	
    
   		 nsigmaTPCTOFkProton = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton)/2.0);
    		nsigmaTPCTOFkKaon = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon)/2.0);
    		nsigmaTPCTOFkPion = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion)/2.0);
  
  	}

   	if(ycut[0])
    	{
    		hman->GetPtHistogram(kHistNSigPionTPC)->Fill(track->P()*track->Charge(),nsigmaTPCkPion );
    		hman->GetPtHistogram(kHistNSigPionPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkPion );
  		hman->GetPtHistogram(kHistNSigPionTPCTOF)->Fill(track->P()*track->Charge(),nsigmaTPCTOFkPion);
  		hman->GetPtHistogram(kHistNSigPionPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkPion);
	}
    	if(ycut[1])
  	 {	
    		hman->GetPtHistogram(kHistNSigKaonTPC)->Fill(track->P()*track->Charge(),nsigmaTPCkKaon );
   	 	hman->GetPtHistogram(kHistNSigKaonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkKaon );	 
  		hman->GetPtHistogram(kHistNSigKaonTPCTOF)->Fill(track->P()*track->Charge(),nsigmaTPCTOFkKaon);
  		hman->GetPtHistogram(kHistNSigKaonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkKaon);
		
	}
	if(ycut[2])
	{					
   		hman->GetPtHistogram(kHistNSigProtonTPC)->Fill(track->P()*track->Charge(),nsigmaTPCkProton ); 
    		hman->GetPtHistogram(kHistNSigProtonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkProton );
  		hman->GetPtHistogram(kHistNSigProtonTPCTOF)->Fill(track->P()*track->Charge(),nsigmaTPCTOFkProton);
  		hman->GetPtHistogram(kHistNSigProtonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkProton);
	}

}

Int_t AliSpectraBothPID::GetParticleSpecie(AliAODMCParticle * part) 
{
  // return PID according to MC truth
  	switch(TMath::Abs(part->PdgCode()))
	{
  		case 2212:
   		return kSpProton;
   		 break;
  		case 321:
    		return kSpKaon;
    		break;
  		case 211:
    		return kSpPion;
    		break;
  		default:
    		return kSpUndefined;
  	} 
}

Int_t AliSpectraBothPID::GetParticleSpecie(TParticle* part) 
{
  // return PID according to MC truth
  switch(TMath::Abs(part->GetPdgCode()))
  {
  	case 2212:
    	return kSpProton;
   	 break;
  	case 321:
    	return kSpKaon;
    	break;
  	case 211:
    	return kSpPion;
    	break;
  	default:
    	return kSpUndefined;
  } 
}

Int_t AliSpectraBothPID::GetParticleSpecie(AliSpectraBothHistoManager * hman,AliVTrack      * trk, AliSpectraBothTrackCuts * trackCuts, Bool_t* rec) 
{
  // return PID according to detectors
  // Get PID response object, if needed

	  // guess the particle based on the smaller nsigma
  	rec[kSpPion]=false;
 	 rec[kSpKaon]=false;
  	rec[kSpProton]=false;

 	if(!fPIDResponse) 
 	 {
    		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
   		 fPIDResponse = inputHandler->GetPIDResponse();
  	}

  	if(!fPIDResponse) 
	{
    		AliFatal("Cannot get pid response");
    	return 0;
  	}


  // Compute nsigma for each hypthesis
  	AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(trk);

  // --- TPC
  	Double_t nsigmaTPCkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton)+fshiftTPC);
  	Double_t nsigmaTPCkKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)+fshiftTPC); 
 	 Double_t nsigmaTPCkPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)+fshiftTPC);
  
 

	Double_t nsigmaTOFkProton=0.0;
	Double_t nsigmaTOFkKaon=0.0;
	Double_t nsigmaTOFkPion=0.0;



  	Double_t nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
  	Double_t nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
  	Double_t nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
	if(fPIDType==kNSigmaTOF||fPIDType==kNSigmaTPCorTOF)
	{
		nsigmaTPCTOFkProton = 10000.0;
  		nsigmaTPCTOFkKaon   = 10000.0;
  		nsigmaTPCTOFkPion   =10000.0;
		nsigmaTOFkProton = 10000.0;
  		nsigmaTOFkKaon   = 10000.0;
  		nsigmaTOFkPion   =10000.0;

	}
	

	if(trackCuts->GetUseTypeDependedTOFCut())
  	{
		if(trackCuts->CheckTOFMatchingParticleType(kSpPion)&&trk->Pt()>=trackCuts->GetPtTOFMatchingPion())
		{	
			nsigmaTOFkPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)+fshiftTOF); 
			nsigmaTPCTOFkPion   = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion)/2.);
		}
		else if (trk->Pt()>=trackCuts->GetPtTOFMatchingPion()) // tracks without proper flags
		{
			nsigmaTOFkPion   = 10000.0; 
			nsigmaTPCTOFkPion   =10000.0;

		}	
			
		if(trackCuts->CheckTOFMatchingParticleType(kSpKaon)&&trk->Pt()>=trackCuts->GetPtTOFMatchingKaon())
		{	
			nsigmaTOFkKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)+fshiftTOF); 
			nsigmaTPCTOFkKaon   = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon)/2.);
		}
		else if (trk->Pt()>=trackCuts->GetPtTOFMatchingKaon()) // tracks without proper flags
		{
			nsigmaTOFkKaon   = 10000.0; 
			nsigmaTPCTOFkKaon   =10000.0;

		}	

		if(trackCuts->CheckTOFMatchingParticleType(kSpProton)&&trk->Pt()>=trackCuts->GetPtTOFMatchingProton())
		{	
			nsigmaTOFkProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton)+fshiftTOF); 
			nsigmaTPCTOFkProton   = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton)/2.);
		}
		else if (trk->Pt()>=trackCuts->GetPtTOFMatchingProton()) // tracks without proper flags
		{
			nsigmaTOFkProton   = 10000.0; 
			nsigmaTPCTOFkProton   =10000.0;

		}	


	}
	else if(trk->Pt()>=trackCuts->GetPtTOFMatching())
  	{
    		nsigmaTOFkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton)+fshiftTOF);
   	 	nsigmaTOFkKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)+fshiftTOF); 
    		nsigmaTOFkPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)+fshiftTOF); 
    		// the TOF info is used in combined
    		nsigmaTPCTOFkProton = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton)/2.);
		nsigmaTPCTOFkKaon   = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon)/2.);
		nsigmaTPCTOFkPion   = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion)/2.);
  	}
	 else
	{ 
		if (fPIDType==kNSigmaTOF)
			return kSpUndefined;

	}
	
  	// select the nsigma to be used for the actual PID
  	Double_t nsigmaPion=10000.0; 
	Double_t nsigmaKaon=10000.0; 
	Double_t nsigmaProton=10000.0;

  	switch (fPIDType) 
	{
  		case kNSigmaTPC:
   		nsigmaProton  =  nsigmaTPCkProton;
    		nsigmaKaon	  =  nsigmaTPCkKaon  ;
    		nsigmaPion    =  nsigmaTPCkPion  ;
    		break;
  		case kNSigmaTOF:
    		nsigmaProton  =  nsigmaTOFkProton;
    		nsigmaKaon	  =  nsigmaTOFkKaon  ;
    		nsigmaPion    =  nsigmaTOFkPion  ;
    		break;
  		case kNSigmacircleTPCTOF:
    		nsigmaProton  =  nsigmaTPCTOFkProton;
    		nsigmaKaon	  =  nsigmaTPCTOFkKaon  ;
    		nsigmaPion    =  nsigmaTPCTOFkPion  ;
    		break;
		case kNSigmasquareTPCTOF:
		nsigmaProton  =  TMath::Max(nsigmaTPCkProton,nsigmaTOFkProton);
    		nsigmaKaon	  =  TMath::Max(nsigmaTPCkKaon,nsigmaTOFkKaon);
    		nsigmaPion    =  TMath::Max(nsigmaTPCkPion,nsigmaTOFkPion)  ;
		break;
		case kNSigmaTPCorTOF:
		if(trackCuts->GetUseTypeDependedTOFCut())
		{
			if(trk->Pt()>=trackCuts->GetPtTOFMatchingPion())
				nsigmaPion    =  nsigmaTOFkPion ;
			else
				nsigmaPion = nsigmaTPCkPion ;
			if(trk->Pt()>=trackCuts->GetPtTOFMatchingKaon())
				nsigmaKaon    =  nsigmaTOFkKaon ;
			else
				nsigmaKaon = nsigmaTPCkKaon ;
			if(trk->Pt()>=trackCuts->GetPtTOFMatchingKaon())
				nsigmaProton    =  nsigmaTOFkProton ;
			else
				nsigmaProton = nsigmaTPCkProton ;

			
		}
		else
		{ 
			if(trk->Pt()>=trackCuts->GetPtTOFMatching())
			{
				nsigmaProton  =  nsigmaTOFkProton;
    				nsigmaKaon	  =  nsigmaTOFkKaon  ;
    				nsigmaPion    =  nsigmaTOFkPion  ;

			}
			else
			{
				nsigmaProton  =  nsigmaTPCkProton;
    				nsigmaKaon	  =  nsigmaTPCkKaon  ;
    				nsigmaPion    =  nsigmaTPCkPion  ;

			}		
		}
    		break;

  	}	
  
  	if(nsigmaPion   < fNSigmaPID)
		rec[kSpPion]=true;
  	if(nsigmaKaon   < fNSigmaPID)
		rec[kSpKaon]=true;
 	 if(nsigmaProton   < fNSigmaPID)
		rec[kSpProton]=true;
	
  	if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) 
		return kSpUndefined;//if is the default value for the three
  	if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < fNSigmaPID))
	{
		if(hman->GetPIDHistogram(kHistPIDTPCKaonRec))
    			hman->GetPIDHistogram(kHistPIDTPCKaonRec)->Fill(trk->GetTPCmomentum(), trk->GetTPCsignal()*trk->Charge()); // Reconstructed PIDKaon histo
    		return kSpKaon;
  	}
  	if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < fNSigmaPID))
	{
		if(hman->GetPIDHistogram(kHistPIDTPCPionRec))
    			hman->GetPIDHistogram(kHistPIDTPCPionRec)->Fill(trk->GetTPCmomentum(), trk->GetTPCsignal()*trk->Charge()); // Reconstructed PIDPion histo
    		return kSpPion;
  	}
  	if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID))
	{
		if(hman->GetPIDHistogram(kHistPIDTPCProtonRec))
    			hman->GetPIDHistogram(kHistPIDTPCProtonRec)->Fill(trk->GetTPCmomentum(), trk->GetTPCsignal()*trk->Charge()); // Reconstructed PIDProton histo
    		return kSpProton;
  	}
  	// else, return undefined
  	return kSpUndefined;

}

Long64_t AliSpectraBothPID::Merge(TCollection* list)
{
  // Merging interface.
  // Returns the number of merged objects (including this).

  Printf("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // Actually, we don't do anything here...
  // collections of all histograms
  //  TList collections;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraBothPID* entry = dynamic_cast<AliSpectraBothPID*> (obj);
    if (entry == 0) 
      continue;

    // TH1I * histo = entry->GetHistoCuts();      
    // collections.Add(histo);
    count++;
  }
  
  //  fHistoCuts->Merge(&collections);
  
  delete iter;
  Printf("OK");
  return count+1;
}

