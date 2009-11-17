#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
//#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
//#include <AliESDVertex.h>
class AliLog;
class AliESDVertex;

#include "AliProtonFeedDownAnalysis.h"
#include "AliProtonAnalysisBase.h"

ClassImp(AliProtonFeedDownAnalysis)

//____________________________________________________________________//
AliProtonFeedDownAnalysis::AliProtonFeedDownAnalysis() : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fProtonContainer(0), fAntiProtonContainer(0), fweightfunction(0),fLambda(0),fLambdaweighted(0),fAntiLambda(0),fAntiLambdaweighted(0)
  {
  //Default constructor
 }
//____________________________________________________________________//
/*AliProtonFeedDownAnalysis::AliProtonFeedDownAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(nbinsY), fMinY(fLowY), fMaxY(fHighY),
  fNBinsPt(nbinsPt), fMinPt(fLowPt), fMaxPt(fHighPt),
  fProtonContainer(0), fAntiProtonContainer(0),fweightfunction(0),fLambda(0),fLambdaweighted(0),fAntiLambda(0),fAntiLambdaweighted(0)
  {
	//Default constructor
	
	//setting up the containers
	Int_t iBin[2];
	iBin[0] = nbinsY;
	iBin[1] = nbinsPt;
	Double_t *binLimY = new Double_t[nbinsY+1];
	Double_t *binLimPt = new Double_t[nbinsPt+1];
	//values for bin lower bounds
	for(Int_t i = 0; i <= nbinsY; i++) 
		binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
	for(Int_t i = 0; i <= nbinsPt; i++) 
		binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;
	
	fProtonContainer = new AliCFContainer("containerProtons","container for protons",4,2,iBin);
	fProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
	fProtonContainer->SetBinLimits(1,binLimPt); //pT
	fAntiProtonContainer = new AliCFContainer("containerAntiProtons","container for antiprotons",4,2,iBin);
	fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
	fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
	
	
}*/
//____________________________________________________________________//
AliProtonFeedDownAnalysis::~AliProtonFeedDownAnalysis() 
{
	//Default destructor
	if(fProtonAnalysisBase) delete fProtonAnalysisBase;
	if(fProtonContainer) delete fProtonContainer;
	if(fAntiProtonContainer) delete fAntiProtonContainer;
	if(fweightfunction) delete fweightfunction; 
	if(fLambda) delete fLambda;
	if(fLambdaweighted) delete fLambdaweighted;
	if(fAntiLambda) delete fAntiLambda;
	if(fAntiLambdaweighted) delete fAntiLambdaweighted;
}
//____________________________________________________________________//
void AliProtonFeedDownAnalysis::InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY, Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) 
{
	//Initializes the histograms
	fNBinsY = nbinsY;
	fMinY = fLowY;
	fMaxY = fHighY;
	fNBinsPt = nbinsPt;
	fMinPt = fLowPt;
	fMaxPt = fHighPt;
	
	
	//setting up the containers
	Int_t iBin[2];
	iBin[0] = nbinsY;
	iBin[1] = nbinsPt;
	Double_t *binLimY = new Double_t[nbinsY+1];
	Double_t *binLimPt = new Double_t[nbinsPt+1];
	//values for bin lower bounds
	for(Int_t i = 0; i <= nbinsY; i++) 
		binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
	for(Int_t i = 0; i <= nbinsPt; i++) 
		binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;
	
	fProtonContainer = new AliCFContainer("containerProtons","container for protons",kNSteps,2,iBin);
	fProtonContainer->SetBinLimits(0,binLimY); //rapidity
	fProtonContainer->SetBinLimits(1,binLimPt); //pT
	fAntiProtonContainer = new AliCFContainer("containerAntiProtons","container for antiprotons",kNSteps,2,iBin);
	fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity
	fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
	/*fLambda=new TH1F("Lambda","Lambda",35,0.5,4.0);
	fLambdaweighted=new TH1F("Lambdaweighted","Lambdaweighted",35,0.5,4.0);
	fAntiLambda=new TH1F("AntiLambda","AntiLambda",35,0.5,4.0);
	fAntiLambdaweighted=new TH1F("AntiLambdaweighted","AntiLambdaweighted",35,0.5,4.0);*/
	fLambda=new TH2F("Lambda","Lambda",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
	fLambdaweighted=new TH2F("Lambdaweighted","Lambdaweighted",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
	fAntiLambda=new TH2F("AntiLambda","AntiLambda",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
	fAntiLambdaweighted=new TH2F("AntiLambdaweighted","AntiLambdaweighted",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
	fLambda->GetYaxis()->SetTitle("P_{T} [GeV/c]");
	fLambdaweighted->GetYaxis()->SetTitle("P_{T} [GeV/c]");
	fAntiLambda->GetYaxis()->SetTitle("P_{T} [GeV/c]");
	fAntiLambdaweighted->GetYaxis()->SetTitle("P_{T} [GeV/c]");
	if(fProtonAnalysisBase->GetEtaMode())
	{
		fLambda->GetXaxis()->SetTitle("#eta");
		fLambdaweighted->GetXaxis()->SetTitle("#eta");
		fAntiLambda->GetXaxis()->SetTitle("#eta");
		fAntiLambdaweighted->GetXaxis()->SetTitle("#eta");
	}
	else
	{
		fLambda->GetXaxis()->SetTitle("y");
		fLambdaweighted->GetXaxis()->SetTitle("y");
		fAntiLambda->GetXaxis()->SetTitle("y");
		fAntiLambdaweighted->GetXaxis()->SetTitle("y");
	}	  
}
//_________________________________________________________________________//
void AliProtonFeedDownAnalysis::Analyze(AliESDEvent *esd, const AliESDVertex *vertex,AliStack *stack)
{	
	Int_t nTracks = 0;
	Int_t nIdentifiedProtons = 0, nIdentifiedAntiProtons = 0;
	Int_t nSurvivedProtons = 0, nSurvivedAntiProtons = 0;
	
	//fHistEvents->Fill(0); //number of analyzed events
	Double_t containerInput[2] ;
	Double_t gPt = 0.0, gP = 0.0;
	nTracks = esd->GetNumberOfTracks();
	Float_t weight;
	Int_t trackflag;
	for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) 
	{
		
		AliESDtrack* track = esd->GetTrack(iTracks);
		AliESDtrack trackTPC;
		Int_t label= track->GetLabel();
	/*	Int_t trackflag=LambdaIsMother(label,stack);//1 mother lambda -1 mother anti lambda 0 mother something else 
		if (trackflag!=0)
			weight=GetWeightforProton(label,stack);	
		else
			weight=1.0;	*/
	
		if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) 
		{
			AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
			if(!tpcTrack)
				 continue;
			gPt = tpcTrack->Pt();
			gP = tpcTrack->P();
			if(fProtonAnalysisBase->IsProton(track)) 
			{
				trackflag=LambdaIsMother(label,stack);//1 mother lambda -1 mother anti lambda 0 mother something else 
				if (trackflag!=0)
					weight=GetWeightforProton(label,stack);	
				else
					weight=1.0;	
				if(tpcTrack->Charge() > 0) 
				{
					nIdentifiedProtons += 1;
					if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
						continue;//track cuts
					if(!fProtonAnalysisBase->IsInPhaseSpace(track)) 
						continue; //track outside the analyzed y-Pt
					nSurvivedProtons += 1;
					if(fProtonAnalysisBase->GetEtaMode()) 
					{
						containerInput[0] = tpcTrack->Eta();
					}
					else 
					{
						//fill the container
						containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz());
					}
					containerInput[1] = gPt;
					fProtonContainer->Fill(containerInput,kSelected,weight);   
					//Feed-down check
					if (trackflag==1)
					{
						fProtonContainer->Fill(containerInput,kSelectedfromLambda,weight);
						TParticle *particle  = stack->Particle(label);
						Int_t lmother=particle->GetFirstMother();
						TParticle *mparticle=stack->Particle(lmother);
						Double_t ptlambda= mparticle->Pt();
						Double_t eta_or_y=0.0;
						if(fProtonAnalysisBase->GetEtaMode())
						 	eta_or_y=mparticle->Eta();
						else
							eta_or_y=0.5*TMath::Log((mparticle->Energy()+mparticle->Pz())/(mparticle->Energy()-mparticle->Pz()));
						fLambda->Fill(eta_or_y, ptlambda);
						fLambdaweighted->Fill(eta_or_y, ptlambda,weight);												 
					}	 
					
				}//protons
				else if(tpcTrack->Charge() < 0) 
				{
					nIdentifiedAntiProtons += 1;
					if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
						continue;//track cuts
					if(!fProtonAnalysisBase->IsInPhaseSpace(track)) 
						continue; //track outside the analyzed y-Pt
					nSurvivedAntiProtons += 1;
					if(fProtonAnalysisBase->GetEtaMode()) 
					{
						containerInput[0] = tpcTrack->Eta();
					}
					else 
					{
						//fill the container
						containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz());
					}
					containerInput[1] = gPt;
					fAntiProtonContainer->Fill(containerInput,kSelected,weight);
					//Feed-down check							
					if(trackflag==-1)
					{
						fAntiProtonContainer->Fill(containerInput,kSelectedfromLambda,weight);
						TParticle *particle  = stack->Particle(label);
						Int_t lmother=particle->GetFirstMother();
						TParticle *mparticle=stack->Particle(lmother);
						Double_t ptlambda= mparticle->Pt();
						Double_t eta_or_y=0.0;
						if(fProtonAnalysisBase->GetEtaMode())
						 	eta_or_y=mparticle->Eta();
						else
							eta_or_y=0.5*TMath::Log((mparticle->Energy()+mparticle->Pz())/(mparticle->Energy()-mparticle->Pz()));
						fAntiLambda->Fill(eta_or_y, ptlambda);
						fAntiLambdaweighted->Fill(eta_or_y, ptlambda,weight);	
					}						
				}//antiprotons   
			}//proton check
		}//TPC only tracks
		else if(fProtonAnalysisBase->GetAnalysisMode() == AliProtonAnalysisBase::kGlobal) 
		{
			gPt = track->Pt();
			gP = track->P();
			if(fProtonAnalysisBase->IsProton(track)) 
			{
				trackflag=LambdaIsMother(label,stack);//1 mother lambda -1 mother anti lambda 0 mother something else 
				if (trackflag!=0)
					weight=GetWeightforProton(label,stack);	
				else
					weight=1.0;	
				if(track->Charge() > 0) 
				{
					nIdentifiedProtons += 1;
					if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
						continue;//track cuts
					if(!fProtonAnalysisBase->IsInPhaseSpace(track)) 
						continue; //track outside the analyzed y-Pt
					nSurvivedProtons += 1;
					if(fProtonAnalysisBase->GetEtaMode()) 
					{
						containerInput[0] = track->Eta();
					}
					else 
					{
						//fill the container
						containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz());
					}
					containerInput[1] = gPt;
					fProtonContainer->Fill(containerInput,kSelected,weight);  
					 //Feed-down check
					if (trackflag==1)
					{
						fProtonContainer->Fill(containerInput,kSelectedfromLambda,weight);
						TParticle *particle  = stack->Particle(label);
						Int_t lmother=particle->GetFirstMother();
						TParticle *mparticle=stack->Particle(lmother);
						Double_t ptlambda= mparticle->Pt();
						Double_t eta_or_y=0.0;
						if(fProtonAnalysisBase->GetEtaMode())
						 	eta_or_y=mparticle->Eta();
						else
							eta_or_y=0.5*TMath::Log((mparticle->Energy()+mparticle->Pz())/(mparticle->Energy()-mparticle->Pz()));
						fLambda->Fill(eta_or_y, ptlambda);
						fLambdaweighted->Fill(eta_or_y, ptlambda,weight);		
					}	 
				}//protons
				else if(track->Charge() < 0) 
				{
					nIdentifiedAntiProtons += 1;
					if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
						continue;//track cuts
					if(!fProtonAnalysisBase->IsInPhaseSpace(track))
						continue; //track outside the analyzed y-Pt
					nSurvivedAntiProtons += 1;
					if(fProtonAnalysisBase->GetEtaMode()) 
					{
						containerInput[0] = track->Eta();
					}
					else 
					{
						//fill the container
						containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz());
					}
					containerInput[1] = gPt;
					fAntiProtonContainer->Fill(containerInput,kSelected,weight);   
					if(trackflag==-1)
					{
						fAntiProtonContainer->Fill(containerInput,kSelectedfromLambda,weight);
						TParticle *particle  = stack->Particle(label);
						Int_t lmother=particle->GetFirstMother();
						TParticle *mparticle=stack->Particle(lmother);
						Double_t ptlambda= mparticle->Pt();
						Double_t eta_or_y=0.0;
						if(fProtonAnalysisBase->GetEtaMode())
						 	eta_or_y=mparticle->Eta();
						else
							eta_or_y=0.5*TMath::Log((mparticle->Energy()+mparticle->Pz())/(mparticle->Energy()-mparticle->Pz()));
						fAntiLambda->Fill(eta_or_y, ptlambda);
						fAntiLambdaweighted->Fill(eta_or_y, ptlambda,weight);
					}	
				}//antiprotons
			}//proton check 
		}//combined tracking
	}//track loop 
	
	if(fProtonAnalysisBase->GetDebugMode())
		Printf("Initial number of tracks: %d | Identified (anti)protons: %d - %d | Survived (anti)protons: %d - %d",nTracks,nIdentifiedProtons,nIdentifiedAntiProtons,nSurvivedProtons,nSurvivedAntiProtons);

}
//_________________________________________________________________________//
void AliProtonFeedDownAnalysis::Analyze(AliAODEvent *fAOD)
{
  // Analysis from AOD: to be implemented...
  // in the near future.
  fAOD->Print();

}
//_________________________________________________________________________//
void AliProtonFeedDownAnalysis::Analyze(AliStack *stack)
{
	 Double_t containerInput[2] ;
	 Float_t weight;
 	for (Int_t ipart=0; ipart<stack->GetNtrack(); ipart++) 
	{ 
		TParticle *particle  = stack->Particle(ipart);
		Int_t code=particle->GetPdgCode();
		/*if (code==3122)
		{
			fLambda->Fill(particle->Pt());
			fLambdaweighted->Fill(particle->Pt(),GetWeightforLambda(particle->Pt()));			
		}
		if (code==-3122)
		{
			fAntiLambda->Fill(particle->Pt());
			fAntiLambdaweighted->Fill(particle->Pt(),GetWeightforLambda(particle->Pt()));			
		}*/
		if (TMath::Abs(code)==2212)
		{
			Int_t trackflag=LambdaIsMother(ipart,stack);//1 mother lambda -1 mother anti lambda 0 mother something else 
			if (trackflag!=0)
				weight=GetWeightforProton(ipart,stack);
			else
				weight=1.0;	
			if(particle->GetUniqueID() == 13) //recject hadronic inter.
				continue; 
			if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) 
				continue;
			if(fProtonAnalysisBase->GetEtaMode()) 
			{
				if((particle->Eta()> fMaxY)||(particle->Eta()< fMinY))
					continue; 
			}	
			else
			{	
				if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz())< fMinY)) 
					continue;
			}		
			if(fProtonAnalysisBase->GetEtaMode()) 
			{
				containerInput[0] =particle->Eta();
			}
			else 
			{
				containerInput[0] = fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz());
			}
			//containerInput[0] = particle->Eta() ;
			containerInput[1] = (Float_t)particle->Pt();
			if (particle->GetPdgCode()==2212)
			{
				fProtonContainer->Fill(containerInput,kAll,weight);
				if(ipart<stack->GetNprimary())
				{
					fProtonContainer->Fill(containerInput,kPrimary,weight);
				}
				else
				{
					if(trackflag==1)
						fProtonContainer->Fill(containerInput,kFromLambda,weight);
				}
			}
			if (particle->GetPdgCode()==-2212)
			{
				fAntiProtonContainer->Fill(containerInput,kAll,weight);
				if(ipart<stack->GetNprimary())
				{
					fAntiProtonContainer->Fill(containerInput,kPrimary,weight);
				}
				else
				{								
					if(trackflag==-1)
						fAntiProtonContainer->Fill(containerInput,kFromLambda,weight);
				}			
			}
		}
	}   
}
//______________________________________________________________________________________________
Int_t AliProtonFeedDownAnalysis::LambdaIsMother(Int_t number, AliStack *stack)
{
	 if(number<0)
	 	return 0;
	TParticle *particle  = stack->Particle(number);
	Int_t motherPDG=0;
	Int_t lmother=-1;
	lmother=particle->GetFirstMother();
	if (lmother<0)		
		return 0;
	TParticle *mparticle=stack->Particle(lmother);
	motherPDG=mparticle->GetPdgCode();										
	if(motherPDG==3122)
		return 1;
	if(motherPDG==-3122)	
		return -1;
	return 0;	
} 
//___________________________________________________________________________________________
Float_t AliProtonFeedDownAnalysis::GetWeightforProton(Int_t number,AliStack *stack)
{
	 if(number<0)
	 	return 1.0;
	TParticle *particle  = stack->Particle(number);
	Int_t lmother=particle->GetFirstMother();
	TParticle *mparticle=stack->Particle(lmother);
	return 	fweightfunction->Eval(mparticle->Pt());
}
//______________________________________________________________________________________________
Float_t AliProtonFeedDownAnalysis::GetWeightforLambda(Float_t pt)
{
	return 	fweightfunction->Eval(pt);
}
