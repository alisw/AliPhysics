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

#include "Riostream.h"

#include <TNode.h> 
#include <TParticle.h> 
#include <TRandom.h> 
#include <TTUBE.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>

#include "AliConst.h" 
#include "AliPDG.h" 
#include "AliRICHGeometry.h"
#include "AliRICHResponse.h"
#include "AliRICHResponseV0.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHv1.h"
#include "AliRun.h"

ClassImp(AliRICHv1)    
//______________________________________________________________________________
AliRICHv1::AliRICHv1(const char *name, const char *title)
          :AliRICH(name,title)
{// Full version of RICH with hits and diagnostics
  if(GetDebug())Info("named ctor","Start.");
  if(GetDebug())Info("named ctor","Stop.");  
}//named ctor
//______________________________________________________________________________
void AliRICHv1::Init()
{//nothing to do here, all the work in ctor or CreateGeometry()
  if(GetDebug())Info("Init","Start.");
  if(GetDebug())Info("Init","Stop.");    
}//void AliRICHv1::Init()
//______________________________________________________________________________
void AliRICHv1::StepManager()
{//Full Step Manager

    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[22];
    static Float_t ckovData[19];
    TLorentzVector position;
    TLorentzVector momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        localPos[3];
    Float_t        localMom[4];
    Float_t        localTheta,localPhi;
    Float_t        theta,phi;
    Float_t        destep, step;
    Double_t        ranf[2];
    Int_t          nPads=-1;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->GetCurrentTrackNumber()];

    
 
    id=gMC->CurrentVolID(copy);
    idvol = copy-1;
    Float_t cherenkovLoss=0;
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    ckovData[6] = 0;                      // dummy track length
    
    /********************Store production parameters for Cerenkov photons************************/ 

    if (gMC->TrackPid() == 50000050) { //is it a Cerenkov photon? 

      //if (gMC->VolId("GAP ")==gMC->CurrentVolID(copy))
        //{                    
	  Float_t ckovEnergy = current->Energy();
	  //energy interval for tracking
	  if  (ckovEnergy > 5.6e-09 && ckovEnergy < 7.8e-09 )       
	    //if (ckovEnergy > 0)
	    {
	      if (gMC->IsTrackEntering()){        //is track entering?
		//printf("Track entered (1)\n");
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (gMC->IsNewTrack()){                          //is it the first step?
		      //printf("I'm in!\n");
		      Int_t mother = current->GetFirstMother(); 
		      
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->GetCurrentTrackNumber();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      //printf("Produced in FREO\n");
		      fCkovNumber++;
		      fFreonProd=1;
		    }    //first step question
		  }        //freo question
		
		if (gMC->IsNewTrack()){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		    }    //quarz question
		}        //first step question
		
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in META\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			
			gMC->Gmtod(mom,localMom,2);
			Float_t cophi = TMath::Cos(TMath::ATan2(localMom[0], localMom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			gMC->GetRandom()->RndmArray(1,ranf);
			if (ranf[0] > t) {
			  gMC->StopTrack();
			  ckovData[13] = 5;
			  AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
			}
			/**********************************************************************************/
		      }    //gap
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			//printf("Got in CSI\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);

			gMC->Gmtod(mom,localMom,2);
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
			Double_t localRt = TMath::Sqrt(localTc);
			localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])));
			Double_t cotheta = TMath::Abs(cos(localTheta));
			Float_t t = Fresnel(ckovEnergy*1e9,cotheta,1);
			    gMC->GetRandom()->RndmArray(1,ranf);
			    if (ranf[0] < t) {
			      gMC->StopTrack();
			      ckovData[13] = 6;
			      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
				
			      //printf("Added One (2)!\n");
			      //printf("Lost by Fresnel\n");
			    }
			    /**********************************************************************************/
		      }
		  } //track entering?
		  
		  
		  /********************Evaluation of losses************************/
		  /******************still in the old fashion**********************/
		  
		  TArrayI procs;
		  Int_t i1 = gMC->StepProcesses(procs);            //number of physics mechanisms acting on the particle
		  for (Int_t i = 0; i < i1; ++i) {
		    //        Reflection loss 
		    if (procs[i] == kPLightReflection) {        //was it reflected
		      ckovData[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=2;
		    } //reflection question
		     
		    //        Absorption loss 
		    else if (procs[i] == kPLightAbsorption) {              //was it absorbed?
		      //printf("Got in absorption\n");
		      ckovData[13]=20;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=11;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=12;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("META")) 
			ckovData[13]=13;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("GAP ")) 
			ckovData[13]=13;
		      
		      if (gMC->CurrentVolID(copy) == gMC->VolId("SRIC")) 
			ckovData[13]=15;
		      
		      //        CsI inefficiency 
		      if (gMC->CurrentVolID(copy) == gMC->VolId("CSI ")) {
			ckovData[13]=16;
		      }
		      gMC->StopTrack();
		      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (procs[i] == kPStop) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      gMC->StopTrack();
		      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    }	// energy treshold question	    
		  }  //number of mechanisms cycle
		  /**********************End of evaluation************************/
		} //freon production question
	    } //energy interval question
    } //cerenkov photon question
      
    /**************************************End of Production Parameters Storing*********************/ 
    
    
    /*******************************Treat photons that hit the CsI (Ckovs and Feedbacks)************/ 
    
    if (gMC->TrackPid() == 50000050 || gMC->TrackPid() == 50000051) {
      
      if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	  if (gMC->Edep() > 0.){
		gMC->TrackPosition(position);
		gMC->TrackMomentum(momentum);
		pos[0]=position(0);
		pos[1]=position(1);
		pos[2]=position(2);
		mom[0]=momentum(0);
		mom[1]=momentum(1);
		mom[2]=momentum(2);
		mom[3]=momentum(3);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;
		

		gMC->Gmtod(pos,localPos,1);

                                                                    
		gMC->Gmtod(mom,localMom,2);

		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		
		((AliRICHChamber*)fChambers->At(idvol))
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if(idvol<kNCH) {	
		    ckovData[0] = gMC->TrackPid();        // particle type
		    ckovData[1] = pos[0];                 // X-position for hit
		    ckovData[2] = pos[1];                 // Y-position for hit
		    ckovData[3] = pos[2];                 // Z-position for hit
		    ckovData[4] = theta;                      // theta angle of incidence
		    ckovData[5] = phi;                      // phi angle of incidence 
		    ckovData[8] = (Float_t) fNsdigits;      // first sdigit
		    ckovData[9] = -1;                       // last pad hit
		    ckovData[13] = 4;                       // photon was detected
		    ckovData[14] = mom[0];
		    ckovData[15] = mom[1];
		    ckovData[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(kBig);
		    cherenkovLoss  += destep;
		    ckovData[7]=cherenkovLoss;
		    
		    //nPads = Hits2SDigits(localPos[0],localPos[2],cherenkovLoss,idvol,kCerenkov);//for photons in CsI kir
		    		    
		    if (fNsdigits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNsdigits;
		    }


		    ckovData[17] = nPads;
		    AliRICHhit *mipHit =  (AliRICHhit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t mipPx = mipHit->MomX();
			Float_t mipPy = mipHit->MomY();
			Float_t mipPz = mipHit->MomZ();
			
			Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			Float_t rt = TMath::Sqrt(r);
			Float_t mipR = mipPx*mipPx + mipPy*mipPy + mipPz*mipPz;	
			Float_t mipRt = TMath::Sqrt(mipR);
			if ((rt*mipRt) > 0)
			  {
			    coscerenkov = (mom[0]*mipPx + mom[1]*mipPy + mom[2]*mipPz)/(rt*mipRt);
			  }
			else
			  {
			    coscerenkov = 0;
			  }
			Float_t cherenkov = TMath::ACos(coscerenkov);
			ckovData[18]=cherenkov;
		      }
		    //if (sector != -1)
		    //{
		    AddHit(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    //printf("Added One (5)!\n");
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge()){//is MIP?
	if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
	  {
	    gMC->TrackMomentum(momentum);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);
	    hits [19] = mom[0];
	    hits [20] = mom[1];
	    hits [21] = mom[2];
	    fFreonProd=1;
	  }

	if(gMC->VolId("GAP ")==gMC->CurrentVolID(copy)) {//is in GAP?
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;
	    gMC->TrackPosition(position);
	    gMC->TrackMomentum(momentum);
	    pos[0]=position(0);
	    pos[1]=position(1);
	    pos[2]=position(2);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);

	    gMC->Gmtod(pos,localPos,1);
	    gMC->Gmtod(mom,localMom,2);
	    ipart  = gMC->TrackPid();
	    destep = gMC->Edep();step   = gMC->TrackStep();// momentum loss and steplength in last step
	    if(gMC->IsTrackEntering()){	    // record hits when track enters ...

		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		

		Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
		Double_t localRt = TMath::Sqrt(localTc);
		localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])))*kRaddeg;                       
		localPhi     = Float_t(TMath::ATan2(Double_t(localMom[2]),Double_t(localMom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = localPos[0];                 // X-position for hit
		hits[2] = localPos[1];                 // Y-position for hit
		hits[3] = localPos[2];                 // Z-position for hit
		hits[4] = localTheta;                  // theta angle of incidence
		hits[5] = localPhi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNsdigits;    // first sdigit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreonProd;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];
		hits[18] = 0;               // dummy cerenkov angle

		tlength = 0;
		eloss   = 0;
		fFreonProd = 0;
	
		Chamber(idvol).LocaltoGlobal(localPos,hits+1);
	   
		
		//To make chamber coordinates x-y had to pass localPos[0], localPos[2]
		xhit    = localPos[0];
		yhit    = localPos[2];
		// Only if not trigger chamber
		if(idvol<kNCH) {
		    //  Initialize hit position (cursor) in the segmentation model 
		    ((AliRICHChamber*)fChambers->At(idvol))
			->SigGenInit(localPos[0], localPos[2], localPos[1]);
		}
	    }
	    
	    // Calculate the charge induced on a pad (disintegration) in case 
	    // Mip left chamber ...
	    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
		gMC->SetMaxStep(kBig);
		eloss   += destep;
		tlength += step;
		
				
		// Only if not trigger chamber
		if(idvol<kNCH) {
		  if (eloss > 0) 
		    {
		      if(gMC->TrackPid() == kNeutron)
			printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		      //nPads = Hits2SDigits(xhit,yhit,eloss,idvol,kMip); //for MIP kir
		      hits[17] = nPads;
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNsdigits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNsdigits;
		}
		
		new(lhits[fNhits++]) AliRICHhit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);
		eloss = 0; 
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    } else if 
		(((AliRICHChamber*)fChambers->At(idvol))
		 ->SigGenCond(localPos[0], localPos[2], localPos[1]))
	    {
		((AliRICHChamber*)fChambers->At(idvol))
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    //nPads = Hits2SDigits(xhit,yhit,eloss,idvol,kMip);//for N kir
		    hits[17] = nPads;
		    //printf("Npads:%d",NPads);
		  }
		xhit     = localPos[0];
		yhit     = localPos[2]; 
		eloss    = destep;
		tlength += step ;
		//
		// nothing special  happened, add up energy loss
	    } else {        
		eloss   += destep;
		tlength += step ;
	    }
	}//is in GAP?
      }//is MIP?
}//void AliRICHv1::StepManager()
