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


/* $Id$ */

//_________________________________________________________________________
//  Class for JetFinder Input preparation from simulated data 
//
//*-- Author: Mark Horner (LBL/UCT)
//



#include "AliEMCALJetFinderInputSimPrep.h"

#include <TParticle.h>
#include <TParticlePDG.h>
#include <TPDGCode.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjectTable.h>
#include <TMath.h>


#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliMagF.h"
#include "AliEMCALFast.h"
#include "AliEMCAL.h"
#include "AliEMCALHit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALLoader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliMC.h"



ClassImp(AliEMCALJetFinderInputSimPrep)
	
AliEMCALJetFinderInputSimPrep::AliEMCALJetFinderInputSimPrep() 
  : fEMCALType(kHits),fSmearType(kSmearEffic),fTrackType(kCharged),fFileType(kPythia),fEfficiency(0.90),
    fTimeCut(0.),fEtaMax(0.),fEtaMin(0.),fPhiMax(0.),fPhiMin(0.)
{
	// Default constructor
if (fDebug > 0) Info("AliEMCALJetFinderInputSimPrep","Beginning Constructor");	
  
  fDebug = 0;
  fSmearType = kSmearEffic;
  fEMCALType = kHits;
  fTrackType = kCharged;
  fEfficiency = 0.90;
  
}
AliEMCALJetFinderInputSimPrep::~AliEMCALJetFinderInputSimPrep()
{
if (fDebug > 0) Info("~AliEMCALJetFinderInputSimPrep","Beginning Destructor");	
	
}

void AliEMCALJetFinderInputSimPrep::Reset(AliEMCALJetFinderResetType_t resettype)
{
	// Method to reset data
  if (fDebug > 1) Info("Reset","Beginning Reset");
	switch (resettype){

	case kResetData:
		fInputObject.Reset(resettype);
		break;
	case kResetTracks:
		fInputObject.Reset(resettype);
		break;
	case kResetDigits:
		fInputObject.Reset(resettype);
		break;
	case kResetParameters:
		fSmearType = kSmearEffic;
	       	fEMCALType = kHits;
		fTrackType = kCharged;
		break;
	case kResetAll:
		fSmearType = kSmearEffic;
	       	fEMCALType = kHits;
		fTrackType = kCharged;
		fInputObject.Reset(resettype);
		break;
	case kResetPartons:
	  Warning("FillFromFile", "kResetPartons not implemented") ; 
	  break;
	case kResetParticles:
	  Warning("FillFromFile", "kResetParticles not implemented") ; 
	  break;
	case kResetJets:
	  Warning("FillFromFile", "kResetJets not implemented") ; 
	  break;
	}// end switch

}

Int_t AliEMCALJetFinderInputSimPrep::FillFromFile(TString *filename, AliEMCALJetFinderFileType_t filetype,Int_t EventNumber,TString data)
{
  if (fDebug > 1) Info("FillFromFile","Beginning FillFromFile");
  fFileType = filetype;	
  AliRunLoader *rl=AliRunLoader::Open(filename->Data());
  if (fDebug > 1) Info("FillFromFile","Opened file %s",filename->Data());
   if (data.Contains("S"))
   {
	   rl->LoadKinematics();
	   rl->LoadSDigits();
	   rl->GetEvent(EventNumber);
   }else
   {
	   rl->LoadKinematics();
	   rl->LoadHits();
	   rl->GetEvent(EventNumber);
   }
   if (fDebug > 1) Info("FillFromFile","Got event %i with option \"XH\"",EventNumber);

   if (	fEMCALType == kHits ||
	fEMCALType == kTimeCut )  
   {
	   if (data.Contains("S"))
	   {
		   FillSDigits();
	   }else
	   {
                   FillHits();
	   }
   }
   if ( fTrackType != kNoTracks  )  FillTracks();	
   if ( fFileType  != kData){
	   FillPartons();
   }
   return 0;	
}

void AliEMCALJetFinderInputSimPrep::FillHits()		
{
	// Fill from the hits to input object from simulation
if (fDebug > 1) Info("FillHits","Beginning FillHits");
	
// Access hit information
//  AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
//  Info("AliEMCALJetFinderInputSimPrep","Title of the geometry is %s",pEMCAL->GetTitle());
//    AliEMCALGeometry* geom =  AliEMCALGeometry::GetInstance();
    AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

    //    TTree *treeH = AliEMCALGetter::Instance()->TreeH();
//    Int_t ntracks = (Int_t) treeH->GetEntries();
//
//   Loop over tracks
//
    Double_t etH = 0.0;

/*    for (Int_t track=0; track<ntracks;track++) {
        gAlice->ResetHits();
        nbytes += treeH->GetEvent(track);*/
//
//  Loop over hits
//
        for(Int_t i =0; i < emcalLoader->Hits()->GetEntries() ;i++)
        {
	    const AliEMCALHit *mHit = emcalLoader->Hit(i);	
	    if (fEMCALType == kTimeCut && 
		    (mHit->GetTime()>fTimeCut) ) continue;
            Float_t x      =    mHit->X();         // x-pos of hit
            Float_t y      =    mHit->Y();         // y-pos
            Float_t z      =    mHit->Z();         // z-pos
            Float_t eloss  =    mHit->GetEnergy(); // deposited energy

            Float_t r      =    TMath::Sqrt(x*x+y*y);
            Float_t theta  =    TMath::ATan2(r,z);
            //Float_t eta    =   -TMath::Log(TMath::Tan(theta/2.));
            //Float_t phi    =    TMath::ATan2(y,x);
            etH = eloss*TMath::Sin(theta);
	    if (fDebug > 10) Info("FillHits","Values of hit energy %i",Int_t(1e7*etH));
	    /*  have to be tune for TRD1; May 31,06 
	    if (geom->TowerIndexFromEtaPhi(eta,180.0/TMath::Pi()*phi) == -1 )
	    {
		    if (fDebug >5)  
		    {
			    Error("FillHits","Hit not inside EMCAL!");
			    mHit->Dump();
		    }
		    continue;
	    }
            fInputObject.AddEnergyToDigit(geom->TowerIndexFromEtaPhi(eta,180.0/TMath::Pi()*phi)-1,Int_t(1e7*etH));
	    */
        } // Hit Loop
//    } // Track Loop


}
void AliEMCALJetFinderInputSimPrep::FillTracks()	
{
	// Fill from particles simulating a TPC to input object from simulation
    if (fDebug > 1) Info("FillTracks","Beginning FillTracks");
	
    AliRunLoader *rl = AliRunLoader::GetRunLoader();
    TParticlePDG* pdgP = 0;
    TParticle *mPart;
    Int_t npart = rl->Stack()->GetNprimary();
    if (fDebug > 1) Info("FillTracks","Header says there are %i primaries",npart);
    Float_t bfield,rEMCAL;		 
    
    if (fDebug > 1) Info("FillTracks","Defining the geometry");
    
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
    if (geom == 0 && pEMCAL) 
      geom = AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
    if (geom == 0) 
      geom = AliEMCALGeometry::GetInstance("EMCAL_55_25","");//","SHISH_77_TRD1_2X2_FINAL_110DEG");
    fEtaMax = geom->GetArm1EtaMax();
    fEtaMin = geom->GetArm1EtaMin();
    fPhiMax = TMath::Pi()/180.0*geom->GetArm1PhiMax();
    fPhiMin = TMath::Pi()/180.0*geom->GetArm1PhiMin();

    if (fDebug > 1) Info("FillTracks","Starting particle loop");
	    
    if (fDebug > 1) Info("FillTracks","Particle loop of %i",npart);
    for (Int_t part = 0; part < npart; part++) {
	mPart = rl->Stack()->Particle(part);
	pdgP = mPart->GetPDG();

	if (fDebug > 5) Info("FillTracks","Checking if track is a primary");
	
	if (fFileType == kPythia) {
	    if (mPart->GetStatusCode() != 1) continue;
	} else if (fFileType == kHijing) {
	    if (mPart->GetFirstDaughter() >= 0 && mPart->GetFirstDaughter() < npart) continue;
	}
	
	if (fDebug > 15) Info("FillTracks","Checking if track (eta - %f, phi - %f) is in acceptance",mPart->Eta(),mPart->Phi());
	if (fDebug > 10) Info("FillTracks","Checking if EMCAL acceptance  ( %f < eta < %f, %f < phi < %f) is in acceptance",fEtaMin,fEtaMax,fPhiMin,fPhiMax);

	if ((!fPythiaComparison)&&(mPart->Eta() > fEtaMax || mPart->Eta() < fEtaMin))   continue;
	if ((!fPythiaComparison)&&(mPart->Phi() > fPhiMax || mPart->Phi() < fPhiMin))   continue;
	
/*
	{kProton, kProtonBar, kElectron, kPositron,
         kNuE, kNuEBar, kGamma, kNeutron, kNeutronBar,
         kMuonPlus, kMuonMinus, kK0Long , kPiPlus, kPiMinus,
         kKPlus, kKMinus, kLambda0, kLambda0Bar, kK0Short,
         kSigmaMinus, kSigmaPlus, kSigma0, kPi0, kK0, kK0Bar,
         0,kNuMu,kNuMuBar
*/

	if (fDebug > 15) Info("FillTracks","Checking if track is rejected");

	if ((fSmearType == kEfficiency 	||
             fSmearType == kSmearEffic)	&& 
 	     pdgP->Charge()!=0) {
	    if (AliEMCALFast::RandomReject(fEfficiency)) {
		continue;
            }
	}

	if (gAlice && gAlice->Field()) 
	  bfield = gAlice->Field()->SolenoidField();
	else
	  bfield = 0.4;
	rEMCAL = AliEMCALGeometry::GetInstance()->GetIPDistance();
	Float_t rB = 3335.6 * mPart->Pt() / bfield;  // [cm]  (case of |charge|=1)
	if (2.*rB < rEMCAL) continue;  // track curls away
	
	//if (part%10) gObjectTable->Print();
        switch(fTrackType)
        {

	   case kAllP:  // All Stable particles to be included
		if (fDebug > 5) Info("FillTracks","Storing track");
		if (fSmearType == kSmear ||
		    fSmearType == kSmearEffic ){
			Smear(mPart);/*
			TParticle *tmp = Smear(MPart);
			fInputObject.AddTrack(Smear(MPart));
			delete tmp;*/
		}else{
			fInputObject.AddTrack(*mPart);
		}
	   break;
	   case kEM:   // All Electromagnetic particles
		if (mPart->GetPdgCode() == kElectron  ||
                    mPart->GetPdgCode() == kMuonPlus  || 
                    mPart->GetPdgCode() == kMuonMinus ||
		    mPart->GetPdgCode() == kPositron ){
		      if (fDebug > 5) Info("FillTracks","Storing electron or positron");
	 	      if (fSmearType == kSmear ||
		            fSmearType == kSmearEffic ){
			      Smear(mPart);/*
 			      TParticle *tmp = Smear(MPart); 
 			      fInputObject.AddTrack(tmp);
 			      delete tmp;*/
		      }else{
		          fInputObject.AddTrack(*mPart);
		      }
		}
		if ( mPart->GetPdgCode() == kGamma ){ 
			fInputObject.AddTrack(*mPart);
			if (fDebug > 5) Info("FillTracks","Storing gamma");
		}
			
	   break;
           case kCharged: // All Particles with non-zero charge
		if (pdgP->Charge() != 0) {
			if (fDebug > 5) Info("FillTracks","Storing charged track");
			if (fSmearType == kSmear ||
	                    fSmearType == kSmearEffic ){
				Smear(mPart);/*
				TParticle *tmp = Smear(MPart);
    				fInputObject.AddTrack(tmp);
				delete tmp;*/
	                }else{
	                 fInputObject.AddTrack(*mPart);
	                }
		}
	   break;
	   case kNeutral: // All particles with no charge
		if (pdgP->Charge() == 0){ 
			fInputObject.AddTrack(*mPart);
			if (fDebug > 5) Info("FillTracks","Storing neutral");
		}
	   break;
	   case kHadron: //All hadrons
		if (mPart->GetPdgCode() != kElectron  &&
                    mPart->GetPdgCode() != kPositron  &&
                    mPart->GetPdgCode() != kMuonPlus  &&
                    mPart->GetPdgCode() != kMuonMinus &&
                    mPart->GetPdgCode() != kGamma ) 
		{
			if (fDebug > 5) Info("FillTracks","Storing hadron");
			if (pdgP->Charge() == 0){
		    	    fInputObject.AddTrack(*mPart);
			}else{
				if (fSmearType == kSmear ||
			            fSmearType == kSmearEffic ){
					Smear(mPart);/*
					TParticle *tmp = Smear(MPart);	
    					fInputObject.AddTrack(tmp);
					delete tmp;*/
			        }else{
			            fInputObject.AddTrack(*mPart);
			        }
			}
		}
	   break;
	   case kChargedHadron:  // only charged hadrons
		if (mPart->GetPdgCode() != kElectron  &&
                    mPart->GetPdgCode() != kPositron  &&
                    mPart->GetPdgCode() != kGamma     &&
                    mPart->GetPdgCode() != kMuonPlus  &&
                    mPart->GetPdgCode() != kMuonMinus &&
		    pdgP->Charge() 	!= 0   	   ){
			if (fDebug > 5) Info("FillTracks","Storing charged hadron");
		       	if (fSmearType == kSmear ||
		            fSmearType == kSmearEffic ){
				Smear(mPart);/*
				TParticle *tmp = Smear(MPart);
 				fInputObject.AddTrack(tmp);
				delete tmp;*/
		        }else{
		               fInputObject.AddTrack(*mPart);
		        }
		}
	   break;
	   case kNoTracks:
	   break;
	   case kEMChargedPi0:
		if (pdgP->Charge() != 0 || mPart->GetPdgCode() == kPi0  ||
			mPart->GetPdgCode() == kGamma     )
		{
			if (fDebug > 5) Info("FillTracks","Storing charged track");
			if (fSmearType == kSmear ||
	                    fSmearType == kSmearEffic ){
				Smear(mPart);/*
				TParticle *tmp = Smear(MPart);
    				fInputObject.AddTrack(tmp);
				delete tmp;*/
	                }else{
	                 fInputObject.AddTrack(*mPart);
	                }
		}
	   break;
	   case kNoNeutronNeutrinoKlong:
  	        if ( mPart->GetPdgCode() != kNeutron 	&&
	             mPart->GetPdgCode() != kNeutronBar &&
		     mPart->GetPdgCode() != kK0Long 	&&
		     mPart->GetPdgCode() != kNuE	&&
                     mPart->GetPdgCode() != kNuEBar     &&
                     mPart->GetPdgCode() != kNuMu       &&
                     mPart->GetPdgCode() != kNuMuBar    &&
                     mPart->GetPdgCode() != kNuTau      &&
                     mPart->GetPdgCode() != kNuTauBar   )
     		{
     			if (fDebug > 5) Info("FillTracks","Storing charged track");
     			if (fSmearType == kSmear ||
     					fSmearType == kSmearEffic ){
     				Smear(mPart);/*
     						TParticle *tmp = Smear(MPart);
     						fInputObject.AddTrack(tmp);
     						delete tmp;*/
     			}else{
     				fInputObject.AddTrack(*mPart);
     			}
     		}
	   break;
	   default:
	   break;
	   //	   delete mPart;
        }	//end of switch
//	Info("FillTracks","After Particle Storage");
	//if (part%10) gObjectTable->Print();

   }	//end of particle loop
   //gObjectTable->Print();	
}

void AliEMCALJetFinderInputSimPrep::FillPartons()	
{
	// Fill partons to
	// input object from simulation
if (fDebug > 1) Info("FillPartons","Beginning FillPartons");

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliGenEventHeader* evHeader = rl->GetHeader()->GenEventHeader();
  Float_t triggerJetValues[4];
  AliEMCALParton tempParton;
  
  if (fFileType == kPythia)
  {
	Int_t ntriggerjets = ((AliGenPythiaEventHeader*)evHeader)->NTriggerJets();
	if (fDebug > 1) Info("FillPartons","Number of trigger jets --> %i",ntriggerjets);
       	for (Int_t jetcounter = 0 ; jetcounter < ntriggerjets; jetcounter++){
      	((AliGenPythiaEventHeader*)evHeader)->TriggerJet(jetcounter,triggerJetValues);
       	TLorentzVector tempLParton;
       	tempLParton.SetPxPyPzE(triggerJetValues[0],triggerJetValues[1],triggerJetValues[2],triggerJetValues[3]);
	tempParton.SetEnergy(tempLParton.Energy());
	tempParton.SetEta(tempLParton.Eta());
	tempParton.SetPhi(tempLParton.Phi());
	FillPartonTracks(&tempParton);
	fInputObject.AddParton(&tempParton);
  	}
 	  
  }
}

void AliEMCALJetFinderInputSimPrep::FillParticles()		
{
	// Fill particles to input object from simulation
if (fDebug > 1) Info("FillParticles","Beginning FillParticles");

    Int_t npart = (gAlice->GetHeader())->GetNprimary();
    TParticlePDG* pdgP = 0;
 
    for (Int_t part = 0; part < npart; part++) {
	TParticle *mPart = gAlice->GetMCApp()->Particle(part);
	pdgP = mPart->GetPDG();
	
	if (fDebug > 10) Info("FillParticles","Checking if particle is a primary");
	
	if (fFileType == kPythia) {
	    if (mPart->GetStatusCode() != 1) continue;
	} else if (fFileType == kHijing) {
	    if (mPart->GetFirstDaughter() >= 0 && mPart->GetFirstDaughter() < npart) continue;
	}

	if (fDebug > 10) Info("FillParticles","Checking if particle is in acceptance");
	
	if ((!fPythiaComparison)&&(mPart->Eta() > fEtaMax || mPart->Eta() < fEtaMin))   continue;
	if ((!fPythiaComparison)&&(mPart->Phi() > fPhiMax || mPart->Phi() < fPhiMin))   continue;
	

/*
	{kProton, kProtonBar, kElectron, kPositron,
         kNuE, kNuEBar, kGamma, kNeutron, kNeutronBar,
         kMuonPlus, kMuonMinus, kK0Long , kPiPlus, kPiMinus,
         kKPlus, kKMinus, kLambda0, kLambda0Bar, kK0Short,
         kSigmaMinus, kSigmaPlus, kSigma0, kPi0, kK0, kK0Bar,
         0,kNuMu,kNuMuBar
*/
        switch(fTrackType)
        {

	   case kAllP:  // All Stable particles to be included
		if (fDebug > 5) Info("FillParticles","Storing particle");
		fInputObject.AddParticle(mPart);
	   break;
	   case kEM:   // All Electromagnetic particles
		if (mPart->GetPdgCode() == kElectron ||
		    mPart->GetPdgCode() == kPositron ||
		    mPart->GetPdgCode() == kGamma){
			if (fDebug > 5) Info("FillParticles","Storing electromagnetic particle");
			fInputObject.AddParticle(mPart);
		}
	   break;
           case kCharged: // All Particles with non-zero charge
		if (pdgP->Charge() != 0) {
			if (fDebug > 5) Info("FillParticles","Storing charged particle");
   			fInputObject.AddParticle(mPart);
		}
	   break;
	   case kNeutral: // All particles with no charge
		if (pdgP->Charge() == 0){
			if (fDebug > 5) Info("FillParticles","Storing neutral particle");
		       	fInputObject.AddParticle(mPart);
		}
	   break;
	   case kHadron: //All hadrons
		if (mPart->GetPdgCode() != kElectron &&
                    mPart->GetPdgCode() != kPositron &&
                    mPart->GetPdgCode() != kGamma ) 
		{

		if (fDebug > 5) Info("FillParticles","Storing hadron");
		    fInputObject.AddParticle(mPart);
		}
	   break;
	   case kChargedHadron:  // only charged hadrons
		if (mPart->GetPdgCode() != kElectron &&
                    mPart->GetPdgCode() != kPositron &&
                    mPart->GetPdgCode() != kGamma    &&
		    pdgP->Charge() 	!= 0   	   ){
		if (fDebug > 5) Info("FillParticles","Storing charged hadron");
	               fInputObject.AddParticle(mPart);
		}
	   break;
	   case kNoTracks:
	   break;
	   case kEMChargedPi0:
		if (pdgP->Charge() != 0 || mPart->GetPdgCode() == kPi0  ||
			mPart->GetPdgCode() == kGamma     )
		{
			if (fDebug > 5) Info("FillTracks","Storing kEMChargedPi0 track");
			if (fSmearType == kSmear ||
	                    fSmearType == kSmearEffic ){
				Smear(mPart);/*
				TParticle *tmp = Smear(MPart);
    				fInputObject.AddTrack(tmp);
				delete tmp;*/
	                }else{
	                 fInputObject.AddTrack(*mPart);
	                }
		}
	   break;
	   case kNoNeutronNeutrinoKlong:
  	        if ( mPart->GetPdgCode() != kNeutron 	&&
	             mPart->GetPdgCode() != kNeutronBar &&
		     mPart->GetPdgCode() != kK0Long 	&&
		     mPart->GetPdgCode() != kNuE	&&
                     mPart->GetPdgCode() != kNuEBar     &&
                     mPart->GetPdgCode() != kNuMu       &&
                     mPart->GetPdgCode() != kNuMuBar    &&
                     mPart->GetPdgCode() != kNuTau      &&
                     mPart->GetPdgCode() != kNuTauBar   )
     		{
     			if (fDebug > 5) Info("FillTracks","Storing kNoNeutronNeutrinoKlong track");
     			if (fSmearType == kSmear ||
     					fSmearType == kSmearEffic ){
     				Smear(mPart);/*
     						TParticle *tmp = Smear(MPart);
     						fInputObject.AddTrack(tmp);
     						delete tmp;*/
     			}else{
     				fInputObject.AddTrack(*mPart);
     			}
     		}
	   break;
	   default:
	   break;
        }	//end of switch
    }// end of loop over particles

}
void AliEMCALJetFinderInputSimPrep::FillDigits()  
{
	// Fill digits to input object

}
void AliEMCALJetFinderInputSimPrep::FillSDigits()  
{
Info("FillSDigits","Beginning FillSDigits");
	AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
	 
	// Fill digits to input object
	for (Int_t towerid=0; towerid < emcalLoader->SDigits()->GetEntries(); towerid++)
	{
		fInputObject.AddEnergyToDigit(emcalLoader->SDigit(towerid)->GetId(), emcalLoader->SDigit(towerid)->GetAmp());
	}

}

void AliEMCALJetFinderInputSimPrep::Smear(TParticle *particle)
{
	// Smear particle momentum
if (fDebug > 5) Info("Smear","Beginning Smear");

 	Float_t tmp = AliEMCALFast::SmearMomentum(1,particle->P());
	Float_t tmpt = (tmp/particle->P()) * particle->Pt();
	if (fDebug > 10) Info("Smear","Smearing changed P from %f to %f",particle->Pt(),tmpt);
	TLorentzVector tmpvector;
	tmpvector.SetPtEtaPhiM(tmpt,particle->Eta(), particle->Phi(), particle->GetMass());
// create a new particle with smearing momentum - and no daughter or parent information and the same
// vertex	

	TParticle tmparticle(particle->GetPdgCode(), 1, 0, 0,0,0, 
			tmpvector.Px()	, tmpvector.Py(), tmpvector.Pz(), tmpvector.Energy(), 
			particle->Vx(), particle->Vy(), particle->Vz(), particle->T());
	fInputObject.AddTrack(tmparticle);	
	return;
}

void AliEMCALJetFinderInputSimPrep::FillPartonTracks(AliEMCALParton *parton)
{
	// Populate parton tracks for input distributions
  if (fDebug>1) Info("FillPartonTracks","Beginning FillPartonTracks");
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
	
	Int_t npart = rl->Stack()->GetNprimary();
	Int_t ntracks =0;
	TParticlePDG *getpdg;
	TParticle *tempPart;
	for (Int_t part = 0; part < npart; part++)
	{
		tempPart = rl->Stack()->Particle(part);
		if (tempPart->GetStatusCode() != 1) continue;
		if ((!fPythiaComparison)&&(tempPart->Eta() > fEtaMax || tempPart->Eta() < fEtaMin || 
		    tempPart->Phi() > fPhiMax || tempPart->Phi() < fPhiMin) ){
		 	if (fDebug>10) Info("FillPartonTracks","Excluding particle not pointing at the EMCAL");    
			continue;
		}
		getpdg = tempPart->GetPDG();
		if (getpdg->Charge() == 0.0  ) { 
			if (fDebug>10) Info("FillPartonTracks","Excluding a neutral particle");
			continue;
		}
		if ( (parton->Eta() - tempPart->Eta())*(parton->Eta() - tempPart->Eta()) +
	     	(parton->Phi() - tempPart->Phi())*(parton->Phi() - tempPart->Phi()) < 1.0 ) ntracks++;
	}
	Float_t *energy = new Float_t[ntracks];
	Float_t *eta    = new Float_t[ntracks];
	Float_t *phi    = new Float_t[ntracks];
	Int_t   *pdg    = new Int_t[ntracks];
	ntracks=0;
	for (Int_t part = 0; part < npart; part++)
	{
		tempPart = rl->Stack()->Particle(part);
		if (tempPart->GetStatusCode() != 1) continue;
		if ((!fPythiaComparison)&&(tempPart->Eta() > fEtaMax || tempPart->Eta() < fEtaMin ||
 		    tempPart->Phi() > fPhiMax || tempPart->Phi() < fPhiMin) ){ 
			if (fDebug>10) Info("FillPartonTracks","Excluding particle not pointing at the EMCAL");       
			continue;
		}
	    	if (tempPart->GetStatusCode() != 1) continue;
  		getpdg = tempPart->GetPDG();
		if (getpdg->Charge() == 0.0  ) {
			if (fDebug>10) Info("FillPartonTracks","Excluding a neutral particle");
			continue;
		}
		if (   (parton->Eta() - tempPart->Eta())*(parton->Eta() - tempPart->Eta()) +
	       	(parton->Phi() - tempPart->Phi())*(parton->Phi() - tempPart->Phi()) < 1.0 )
	       	{
	  		energy[ntracks] = tempPart->Pt();
	      		eta[ntracks] = tempPart->Eta();
		  	phi[ntracks] = tempPart->Phi();
	       		pdg[ntracks] = tempPart->GetPdgCode();
	  		ntracks++;
		}
	}
	parton->SetTrackList(ntracks,energy,eta,phi,pdg);
	delete[] energy;
	delete[] eta;
	delete[] phi;
	delete[] pdg;
}

