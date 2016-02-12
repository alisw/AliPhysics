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
#include "AliTOFGeometry.h"
#include "AliESDtrack.h"	

ClassImp(AliSpectraBothPID)

AliSpectraBothPID::AliSpectraBothPID() : TNamed("PID", "PID object"), fPIDType(kNSigmacircleTPCTOF), fNSigmaPID(3), fPIDResponse(0),fshiftTPC(0),fshiftTOF(0),foldT0(0),fNoldT0bins(-1) 

{

}

AliSpectraBothPID::AliSpectraBothPID(BothPIDType_t pidType) : TNamed("PID", "PID object"), fPIDType(pidType), fNSigmaPID(3), fPIDResponse(0), fshiftTPC(0),fshiftTOF(0),fNoldT0bins(-1) 
{



}



void AliSpectraBothPID::FillQAHistos(AliSpectraBothHistoManager * hman, AliVTrack * track, AliSpectraBothTrackCuts * trackCuts,Int_t idGen) 
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
 	if(track->Pt()>=tmppt)
  	{
    
   		 hman->GetPIDHistogram(kHistPIDTOF)->Fill(track->P(),(track->GetTOFsignal()/100)*track->Charge()); // PID histo
    			
    //TOF	 	 
		if(ycut[0]&&trackCuts->CheckTOFMatchingParticleType(kSpPion)&&track->Pt()>=trackCuts->GetPtTOFMatchingPion())
    		{
			nsigmaTOFkPion = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)+fshiftTOF;
			nsigmaTPCTOFkPion = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion)/2.0);
    			hman->GetPtHistogram(kHistNSigPionPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkPion );
  			hman->GetPtHistogram(kHistNSigPionPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkPion);
			if(fNoldT0bins>0.0)
			{
				Float_t missmatch=GetMissMatchNsigma(track, AliPID::kPion);
				hman->GetPtHistogram(kHistNSigPionPtTOFmissmatch)->Fill(track->Pt()*track->Charge(),missmatch);
			
			}

			if(idGen==kSpPion)
			{
				hman->GetPtHistogram(kHistNSigTruePionPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkPion );
				hman->GetPtHistogram(kHistNSigTruePionPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkPion);
			}

    		}
    		if(ycut[1]&&trackCuts->CheckTOFMatchingParticleType(kSpKaon)&&track->Pt()>=trackCuts->GetPtTOFMatchingKaon())
  	 	{	
			nsigmaTOFkKaon = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)+fshiftTOF; 
			nsigmaTPCTOFkKaon = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon)/2.0);

   	 		hman->GetPtHistogram(kHistNSigKaonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkKaon );	 
  			hman->GetPtHistogram(kHistNSigKaonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkKaon);
			if(fNoldT0bins>0.0)
			{
				Float_t missmatch=GetMissMatchNsigma(track, AliPID::kKaon);
				hman->GetPtHistogram(kHistNSigKaonPtTOFmissmatch)->Fill(track->Pt()*track->Charge(),missmatch);
			
			}

			if(idGen==kSpKaon)
			{
				hman->GetPtHistogram(kHistNSigTrueKaonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkKaon );
				hman->GetPtHistogram(kHistNSigTrueKaonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkKaon);
			}

   		}
		if(ycut[2]&&trackCuts->CheckTOFMatchingParticleType(kSpProton)&&track->Pt()>=trackCuts->GetPtTOFMatchingProton())
		{		
			nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton)+fshiftTOF;
			nsigmaTPCTOFkProton = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton)/2.0);
			
    			hman->GetPtHistogram(kHistNSigProtonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkProton );
  			hman->GetPtHistogram(kHistNSigProtonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkProton);
			if(fNoldT0bins>0.0)
			{
				Float_t missmatch=GetMissMatchNsigma(track, AliPID::kProton);
				hman->GetPtHistogram(kHistNSigProtonPtTOFmissmatch)->Fill(track->Pt()*track->Charge(),missmatch);
			
			}

			if(idGen==kSpProton)
			{
				hman->GetPtHistogram(kHistNSigTrueProtonPtTOF)->Fill(track->Pt()*track->Charge(),nsigmaTOFkProton );
				hman->GetPtHistogram(kHistNSigTrueProtonPtTPCTOF)->Fill(track->Pt()*track->Charge(),nsigmaTPCTOFkProton);
			}

		}   
  
  	}

   	if(ycut[0])
    	{
		if(idGen==kSpPion)	
    			hman->GetPtHistogram(kHistNSigTruePionPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkPion );
    		hman->GetPtHistogram(kHistNSigPionPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkPion );
	}
    	if(ycut[1])
  	 {
		if(idGen==kSpKaon)
 	   		hman->GetPtHistogram(kHistNSigTrueKaonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkKaon );
   	 	hman->GetPtHistogram(kHistNSigKaonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkKaon );	 
		
	}
	if(ycut[2])
	{
		if(idGen==kSpProton)	
   			hman->GetPtHistogram(kHistNSigTrueProtonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkProton ); 
    		hman->GetPtHistogram(kHistNSigProtonPtTPC)->Fill(track->Pt()*track->Charge(),nsigmaTPCkProton );
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
			if(trk->Pt()>=trackCuts->GetPtTOFMatchingProton())
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
//________________________________________________________________________________________________
void AliSpectraBothPID::SetoldT0()
{
	if(!fPIDResponse) 
  	{
    		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   		 AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    		fPIDResponse = inputHandler->GetPIDResponse();
  	}
	Int_t bin=fPIDResponse->GetTOFResponse().GetMomBin(1.0);
	Int_t nbins=fPIDResponse->GetTOFResponse().GetNmomBins();
	if(fNoldT0bins>0&&fNoldT0bins!=nbins)
		delete foldT0;
	if(!foldT0)
	{
		fNoldT0bins=nbins;
		foldT0= new Float_t[fNoldT0bins];
	}
	//else
	//	Printf("TOF1 %lf",foldT0[bin]);	
	for (int k=0;k<fNoldT0bins;k++)
	{
			foldT0[k]=fPIDResponse->GetTOFResponse().GetT0bin(k);
	}
	//Printf("TOF2 %lf",foldT0[bin]);
}
//___________________________________________________________________________________________________________
Float_t AliSpectraBothPID::GetMissMatchNsigma(AliVTrack* track,AliPID::EParticleType type)
{
	AliESDtrack* trackESD=dynamic_cast<AliESDtrack*>(track);
	Float_t timeexp=-1e8;
	if(trackESD)
	{
		int index=trackESD->GetTOFCalChannel();

		AliTOFGeometry tofGeo;
			Float_t c = TMath::C() * 1.e2 / 1.e12; /* cm/ps */
		Float_t c_1 = 1. / c;
		Int_t det[5];
	  	Float_t length, pos[3];
  
 			 /* compute length and expected time */
	 	tofGeo.GetVolumeIndices(index, det);
		tofGeo.GetPosPar(det, pos);
		length = 0.;
		for (Int_t i = 0; i < 3; i++) 
			length += pos[i] * pos[i];
		length = TMath::Sqrt(length);
  		timeexp = length * c_1;
	}
	else
		timeexp =track->GetTOFsignal();

	Float_t expTime=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,type);
	Float_t expres=fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type));
	Int_t bin=fPIDResponse->GetTOFResponse().GetMomBin(track->P());

	Float_t startTime=foldT0[bin];
	if(expres>0.0)
		return (timeexp-startTime-expTime)/expres;
	else
		return 0.0;
}
