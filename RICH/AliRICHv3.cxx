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


#include <Riostream.h>

#include <TBRIK.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TNode.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TVirtualMC.h>
#include <TPDGCode.h> //for kNuetron
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPDG.h"
#include "AliRICHGeometry.h"
#include "AliRICHResponseV0.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHv3.h"
#include "AliRun.h"
#include "AliRICHRawCluster.h"
#include "AliRICHDigit.h"
#include "AliRICHRecHit1D.h"
#include "AliMC.h"


ClassImp(AliRICHv3)

//______________________________________________________________
//    Implementation of the RICH version 3 with azimuthal rotation


AliRICHv3::AliRICHv3(const char *sName, const char *sTitle)
      	  :AliRICH(sName,sTitle)
{
// The named ctor currently creates a single copy of 
// AliRICHGeometry AliRICHSegmentationV1 AliRICHResponseV0
// and initialises the corresponding models of all 7 chambers with these stuctures.
// Note: all chambers share the single copy of models. MUST be changed later (???).
  if(GetDebug())Info("named ctor","Start.");

   fCkovNumber=fFreonProd=0;
   
   AliRICHGeometry       *pRICHGeometry    =new AliRICHGeometry;           // ??? to be moved to AlRICHChamber::named ctor
   AliRICHSegmentationV1 *pRICHSegmentation=new AliRICHSegmentationV1;     // ??? to be moved to AlRICHChamber::named ctor
   AliRICHResponseV0     *pRICHResponse    =new AliRICHResponseV0;         // ??? to be moved to AlRICHChamber::named ctor
     
   for (Int_t i=1; i<=kNCH; i++){    
      SetGeometryModel(i,pRICHGeometry);
      SetSegmentationModel(i,pRICHSegmentation);
      SetResponseModel(i,pRICHResponse);
      C(i)->Init(i); // ??? to be removed     
   }
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICHv3::ctor(const char *pcName, const char *pcTitle)

AliRICHv3::~AliRICHv3()
{
// Dtor deletes RICH models. In future (???) AliRICHChamber will be responsible for that.
   if(GetDebug()) cout<<ClassName()<<"::dtor()>\n";
      
   if(fChambers) {
     AliRICHChamber *ch =C(1); 
     if(ch) {
       delete ch->GetGeometryModel();
       delete ch->GetResponseModel();
       delete ch->GetSegmentationModel();
     }
   }
}//AliRICHv3::dtor()
//______________________________________________________________________________
void AliRICHv3::StepManager()
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
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TParticle *current = (TParticle*)(*gAlice->GetMCApp()->Particles())[gAlice->GetMCApp()->GetCurrentTrackNumber()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
 
    id=gMC->CurrentVolID(copy);
    idvol = copy-1;
    Float_t cherenkovLoss=0;
    //gAlice->KeepTrack(gAlice->GetCurrentTrackNumber());
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    //bzero((char *)ckovData,sizeof(ckovData)*19);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    ckovData[6] = 0;                      // dummy track length
    //ckovData[11] = gAlice->GetCurrentTrackNumber();
    
    //printf("\n+++++++++++\nTrack: %d\n++++++++++++\n",gAlice->GetCurrentTrackNumber());

    //AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    /********************Store production parameters for Cerenkov photons************************/ 
//is it a Cerenkov photon? 
    if (gMC->TrackPid() == 50000050) { 

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
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->GetMCApp()->GetCurrentTrackNumber();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      //printf("Produced in FREO\n");
		      fCkovNumber++;
		      fFreonProd=1;
		      //printf("Index: %d\n",fCkovNumber);
		    }    //first step question
		  }        //freo question
		
		if (gMC->IsNewTrack()){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		      //printf("Produced in QUAR\n");
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreonProd);
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Track entered (2)\n");
		    //printf("Current volume (should be META): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("META"),gMC->CurrentVolID(copy));
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
			  AddCerenkov(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
			  //printf("Added One (1)!\n");
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    //printf("Current volume (should be CSI) (1): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
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
			      AddCerenkov(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
				
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
		      //gMC->StopTrack();
		      //AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
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
		      AddCerenkov(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
		      //printf("Added One (3)!\n");
		      //printf("Added cerenkov %d\n",fCkovNumber);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (procs[i] == kPStop) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      gMC->StopTrack();
		      AddCerenkov(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
		      //printf("Added One (4)!\n");
		    }	// energy treshold question	    
		  }  //number of mechanisms cycle
		  /**********************End of evaluation************************/
		} //freon production question
	    } //energy interval question
	//}//inside the proximity gap question
    } //cerenkov photon question
      
    /**************************************End of Production Parameters Storing*********************/ 
    
    
    /*******************************Treat photons that hit the CsI (Ckovs and Feedbacks)************/ 
    
    if (gMC->TrackPid() == 50000050 || gMC->TrackPid() == 50000051) {
      //printf("Cerenkov\n");
      
      //if (gMC->TrackPid() == 50000051)
	//printf("Tracking a feedback\n");
      
      if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	  //printf("Current volume (should be CSI) (2): %s\n",gMC->CurrentVolName());
	  //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
	  //printf("Got in CSI\n");
	  //printf("Tracking a %d\n",gMC->TrackPid());
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

		//Chamber(idvol).GlobaltoLocal(pos,localPos);
                                                                    
		gMC->Gmtod(mom,localMom,2);

		//Chamber(idvol).GlobaltoLocal(mom,localMom);
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		//Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
		//printf("Sector:%d\n",sector);

		/*if (gMC->TrackPid() == 50000051){
		  fFeedbacks++;
		  printf("Feedbacks:%d\n",fFeedbacks);
		}*/	
		
        //PH		((AliRICHChamber*) (*fChambers)[idvol])
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
		    
		    ckovData[17] = Hits2SDigits(localPos[0],localPos[2],cherenkovLoss,idvol,kPhoton);//for photons in CsI 
		    		    
		    if (fNsdigits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNsdigits;
		    }

		    
		    //TClonesArray *Hits = RICH->Hits();
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
		    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
		    AddCerenkov(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,ckovData);
		    //printf("Added One (5)!\n");
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge()){
//If MIP
	/*if (gMC->IsTrackEntering())
	  {                
	    hits[13]=20;//is track entering?
	  }*/
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

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {//is in GAP?
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;

	    //Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
	    //printf("Sector:%d\n",sector);
	    
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
	    
	    //Chamber(idvol).GlobaltoLocal(pos,localPos);
                                                                    
	    gMC->Gmtod(mom,localMom,2);

	    //Chamber(idvol).GlobaltoLocal(mom,localMom);
	    
	    ipart  = gMC->TrackPid();
	    //
	    // momentum loss and steplength in last step
	    destep = gMC->Edep();
	    step   = gMC->TrackStep();
  
	    //
	    // record hits when track enters ...
	    if( gMC->IsTrackEntering()) {
//		gMC->SetMaxStep(fMaxStepGas);
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
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
          //PH		    ((AliRICHChamber*) (*fChambers)[idvol])
		    ((AliRICHChamber*)fChambers->At(idvol))
			->SigGenInit(localPos[0], localPos[2], localPos[1]);
		}
	    }
	    
	    // 
	    // Calculate the charge induced on a pad (disintegration) in case 
	    //
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
		      hits[17] = Hits2SDigits(xhit,yhit,eloss,idvol,kMip); //for MIP 
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNsdigits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNsdigits;
		}
		
		//if(sector !=-1)
		new(lhits[fNhits++]) AliRICHhit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    }else if(((AliRICHChamber*)fChambers->At(idvol))->SigGenCond(localPos[0], localPos[2], localPos[1])){
		((AliRICHChamber*)fChambers->At(idvol))->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    hits[17] = Hits2SDigits(xhit,yhit,eloss,idvol,kMip);//for n
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
    /*************************************************End of MIP treatment**************************************/
}//void AliRICHv3::StepManager()
//__________________________________________________________________________________________________
Int_t AliRICHv3::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
{//calls the charge disintegration method of the current chamber and adds all generated sdigits to the list of digits
   
   Float_t newclust[4][500];
   Int_t clhits[5];
   Int_t iNdigits;
   clhits[0]=fNhits+1;
   
  ((AliRICHChamber*)fChambers->At(idvol))->DisIntegration(eloss, xhit, yhit, iNdigits,newclust, res);
    
  for (Int_t i=0; i<iNdigits; i++) {
    if (Int_t(newclust[0][i]) > 0) {
            clhits[1] = Int_t(newclust[0][i]);//  Cluster Charge
            clhits[2] = Int_t(newclust[1][i]);//  Pad: ix
            clhits[3] = Int_t(newclust[2][i]);//  Pad: iy
            clhits[4] = Int_t(newclust[3][i]);//  Pad: chamber sector
            AddSpecialOld(clhits);
        }
    }
  return iNdigits;
}//Int_t AliRICHv3::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
//__________________________________________________________________________________________________
void AliRICHv3::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
{
  
  Int_t NpadX = 162;                 // number of pads on X
  Int_t NpadY = 162;                 // number of pads on Y
  
  Int_t Pad[162][162];
  for (Int_t i=0;i<NpadX;i++) {
    for (Int_t j=0;j<NpadY;j++) {
      Pad[i][j]=0;
    }
  }
  
  //  Create some histograms

  TH1F *pionspectra1 = new TH1F("pionspectra1","Pion Spectra",200,-4,2);
  TH1F *pionspectra2 = new TH1F("pionspectra2","Pion Spectra",200,-4,2);
  TH1F *pionspectra3 = new TH1F("pionspectra3","Pion Spectra",200,-4,2);
  TH1F *protonspectra1 = new TH1F("protonspectra1","Proton Spectra",200,-4,2);
  TH1F *protonspectra2 = new TH1F("protonspectra2","Proton Spectra",200,-4,2);
  TH1F *protonspectra3 = new TH1F("protonspectra3","Proton Spectra",200,-4,2);
  TH1F *kaonspectra1 = new TH1F("kaonspectra1","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra2 = new TH1F("kaonspectra2","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra3 = new TH1F("kaonspectra3","Kaon Spectra",100,-4,2);
  TH1F *electronspectra1 = new TH1F("electronspectra1","Electron Spectra",100,-4,2);
  TH1F *electronspectra2 = new TH1F("electronspectra2","Electron Spectra",100,-4,2);
  TH1F *electronspectra3 = new TH1F("electronspectra3","Electron Spectra",100,-4,2);
  TH1F *muonspectra1 = new TH1F("muonspectra1","Muon Spectra",100,-4,2);
  TH1F *muonspectra2 = new TH1F("muonspectra2","Muon Spectra",100,-4,2);
  TH1F *muonspectra3 = new TH1F("muonspectra3","Muon Spectra",100,-4,2);
  TH1F *neutronspectra1 = new TH1F("neutronspectra1","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra2 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra3 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *chargedspectra1 = new TH1F("chargedspectra1","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra2 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra3 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *pionptspectrafinal = new TH1F("pionptspectrafinal","Primary Pions Transverse Momenta at HMPID",20,0,5);
  TH1F *pionptspectravertex = new TH1F("pionptspectravertex","Primary Pions Transverse Momenta at vertex",20,0,5);
  TH1F *kaonptspectrafinal = new TH1F("kaonptspectrafinal","Primary Kaons Transverse Momenta at HMPID",20,0,5);
  TH1F *kaonptspectravertex = new TH1F("kaonptspectravertex","Primary Kaons Transverse Momenta at vertex",20,0,5);
  //TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",100,-180,180);
  TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of Theta angle of incidence, all tracks",100,0,50);
  TH1F *hitsTheta500MeV = new TH1F("hitsTheta500MeV","Distribution of Theta angle of incidence, 0.5-1 GeV primary tracks",100,0,50);
  TH1F *hitsTheta1GeV = new TH1F("hitsTheta1GeV","Distribution of Theta angle of incidence, 1-2 GeV primary tracks",100,0,50);
  TH1F *hitsTheta2GeV = new TH1F("hitsTheta2GeV","Distribution of Theta angle of incidence, 2-3 GeV primary tracks",100,0,50);
  TH1F *hitsTheta3GeV = new TH1F("hitsTheta3GeV","Distribution of Theta angle of incidence, >3 GeV primary tracks",100,0,50);
  TH2F *production = new TH2F("production","Mother production vertices",100,-300,300,100,0,600);
   
   
   

//   Start loop over events 

  Int_t pion=0, kaon=0, proton=0, electron=0, positron=0, neutron=0, highneutrons=0, muon=0;
  Int_t chargedpions=0,primarypions=0,highprimarypions=0,chargedkaons=0,primarykaons=0,highprimarykaons=0;
  Int_t photons=0, primaryphotons=0, highprimaryphotons=0;
  TRandom* random=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       
       AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
     
       TTree *treeH = TreeH();
       Int_t ntracks =(Int_t) treeH->GetEntries();
            
// Start loop on tracks in the hits containers
       
       for (Int_t track=0; track<ntracks;track++) {
	   printf ("Processing Track: %d\n",track);
	   gAlice->ResetHits();
	   treeH->GetEvent(track);
	   	   	   
	   for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	       mHit;
	       mHit=(AliRICHhit*)pRICH->NextHit()) 
	     {
	       //Int_t nch  = mHit->fChamber;              // chamber number
	       //Float_t x  = mHit->X();                    // x-pos of hit
	       //Float_t y  = mHit->Z();                    // y-pos
	       //Float_t z  = mHit->Y();
	       //Float_t phi = mHit->Phi();                 //Phi angle of incidence
	       Float_t theta = mHit->Theta();             //Theta angle of incidence
	       Float_t px = mHit->MomX();
	       Float_t py = mHit->MomY();
	       Int_t index = mHit->Track();
	       Int_t particle = (Int_t)(mHit->Particle());    
	       Float_t R;
	       Float_t PTfinal;
	       Float_t PTvertex;

	      TParticle *current = gAlice->GetMCApp()->Particle(index);
	      
	      //Float_t energy=current->Energy(); 

	      R=TMath::Sqrt(current->Vx()*current->Vx() + current->Vy()*current->Vy());
	      PTfinal=TMath::Sqrt(px*px + py*py);
	      PTvertex=TMath::Sqrt(current->Px()*current->Px() + current->Py()*current->Py());
	      
	      

	      if (TMath::Abs(particle) < 10000000)
		{
		  hitsTheta->Fill(theta,(float) 1);
		  if (R<5)
		    {
		      if (PTvertex>.5 && PTvertex<=1)
			{
			  hitsTheta500MeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>1 && PTvertex<=2)
			{
			  hitsTheta1GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>2 && PTvertex<=3)
			{
			  hitsTheta2GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>3)
			{
			  hitsTheta3GeV->Fill(theta,(float) 1);
			}
		    }
		  
		}

	      //if (nch == 3)
		//{
	      
	      if (TMath::Abs(particle) < 50000051)
		{
		  //if (TMath::Abs(particle) == 50000050 || TMath::Abs(particle) == 2112)
		  if (TMath::Abs(particle) == 2112 || TMath::Abs(particle) == 50000050)
		    {
		      //gMC->Rndm(&random, 1);
		      if (random->Rndm() < .1)
			production->Fill(current->Vz(),R,(float) 1);
		      if (TMath::Abs(particle) == 50000050)
			//if (TMath::Abs(particle) > 50000000)
			{
			  photons +=1;
			  if (R<5)
			    {
			      primaryphotons +=1;
			      if (current->Energy()>0.001)
				highprimaryphotons +=1;
			    }
			}	
		      if (TMath::Abs(particle) == 2112)
			{
			  neutron +=1;
			  if (current->Energy()>0.0001)
			    highneutrons +=1;
			}
		    }
		  if (TMath::Abs(particle) < 50000000)
		    {
		      production->Fill(current->Vz(),R,(float) 1);
		    }
		  //mip->Fill(x,y,(float) 1);
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  if (R<5)
		    {
		      pionptspectravertex->Fill(PTvertex,(float) 1);
		      pionptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  if (R<5)
		    {
		      kaonptspectravertex->Fill(PTvertex,(float) 1);
		      kaonptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    pionspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      pionspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  pion +=1;
		  if (TMath::Abs(particle)==211)
		    {
		      chargedpions +=1;
		      if (R<5)
			{
			  primarypions +=1;
			  if (current->Energy()>1)
			    highprimarypions +=1;
			}
		    }	
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    protonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    protonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  proton +=1;
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    kaonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    kaonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  kaon +=1;
		  if (TMath::Abs(particle)==321)
		    {
		      chargedkaons +=1;
		      if (R<5)
			{
			  primarykaons +=1;
			  if (current->Energy()>1)
			    highprimarykaons +=1;
			}
		    }
		}
	      if (TMath::Abs(particle)==11)
		{
		  electronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    electronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    electronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (particle == 11)
		    electron +=1;
		  if (particle == -11)
		    positron +=1;
		}
	      if (TMath::Abs(particle)==13)
		{
		  muonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    muonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    muonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  muon +=1;
		}
	      if (TMath::Abs(particle)==2112)
		{
		  neutronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    neutronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      neutronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  neutron +=1;
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy()-current->GetCalcMass()>1)
		    {
		      chargedspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
			chargedspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (R>250 && R<450)
			chargedspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		}
	      // Fill the histograms
	      //Nh1+=nhits;
	      //h->Fill(x,y,(float) 1);
	      //}
	      //}
	   }          
	   
       }
       
   }
   //   }

   TStyle *mystyle=new TStyle("Plain","mystyle");
   mystyle->SetPalette(1,0);
   mystyle->cd();
   
   //Create canvases, set the view range, show histograms

    TCanvas *c2 = new TCanvas("c2","Angles of incidence",150,150,100,150);
    c2->Divide(2,2);
    //c2->SetFillColor(42);
    
    c2->cd(1);
    hitsTheta500MeV->SetFillColor(5);
    hitsTheta500MeV->Draw();
    c2->cd(2);
    hitsTheta1GeV->SetFillColor(5);
    hitsTheta1GeV->Draw();
    c2->cd(3);
    hitsTheta2GeV->SetFillColor(5);
    hitsTheta2GeV->Draw();
    c2->cd(4);
    hitsTheta3GeV->SetFillColor(5);
    hitsTheta3GeV->Draw();
    
            
   
    TCanvas *c15 = new TCanvas("c15","Mothers Production Vertices",50,50,600,600);
    c15->cd();
    production->SetFillColor(42);
    production->SetXTitle("z (m)");
    production->SetYTitle("R (m)");
    production->Draw();

    TCanvas *c10 = new TCanvas("c10","Pt Spectra",50,50,600,700);
    c10->Divide(2,2);
    c10->cd(1);
    pionptspectravertex->SetFillColor(5);
    pionptspectravertex->SetXTitle("Pt (GeV)");
    pionptspectravertex->Draw();
    c10->cd(2);
    pionptspectrafinal->SetFillColor(5);
    pionptspectrafinal->SetXTitle("Pt (GeV)");
    pionptspectrafinal->Draw();
    c10->cd(3);
    kaonptspectravertex->SetFillColor(5);
    kaonptspectravertex->SetXTitle("Pt (GeV)");
    kaonptspectravertex->Draw();
    c10->cd(4);
    kaonptspectrafinal->SetFillColor(5);
    kaonptspectrafinal->SetXTitle("Pt (GeV)");
    kaonptspectrafinal->Draw();
   
  
   TCanvas *c16 = new TCanvas("c16","Particles Spectra II",150,150,600,350);
   c16->Divide(2,1);
   
   c16->cd(1);
   //TCanvas *c13 = new TCanvas("c13","Electron Spectra",400,10,600,700);
   electronspectra1->SetFillColor(5);
   electronspectra1->SetXTitle("log(GeV)");
   electronspectra2->SetFillColor(46);
   electronspectra2->SetXTitle("log(GeV)");
   electronspectra3->SetFillColor(10);
   electronspectra3->SetXTitle("log(GeV)");
   //c13->SetLogx();
   electronspectra1->Draw();
   electronspectra2->Draw("same");
   electronspectra3->Draw("same");
   
   c16->cd(2);
   //TCanvas *c14 = new TCanvas("c14","Muon Spectra",400,10,600,700);
   muonspectra1->SetFillColor(5);
   muonspectra1->SetXTitle("log(GeV)");
   muonspectra2->SetFillColor(46);
   muonspectra2->SetXTitle("log(GeV)");
   muonspectra3->SetFillColor(10);
   muonspectra3->SetXTitle("log(GeV)");
   //c14->SetLogx();
   muonspectra1->Draw();
   muonspectra2->Draw("same");
   muonspectra3->Draw("same");
   
   //c16->cd(3);
   //TCanvas *c16 = new TCanvas("c16","Neutron Spectra",400,10,600,700);
   //neutronspectra1->SetFillColor(42);
   //neutronspectra1->SetXTitle("log(GeV)");
   //neutronspectra2->SetFillColor(46);
   //neutronspectra2->SetXTitle("log(GeV)");
   //neutronspectra3->SetFillColor(10);
   //neutronspectra3->SetXTitle("log(GeV)");
   //c16->SetLogx();
   //neutronspectra1->Draw();
   //neutronspectra2->Draw("same");
   //neutronspectra3->Draw("same");

   TCanvas *c9 = new TCanvas("c9","Particles Spectra",150,150,600,700);
   //TCanvas *c9 = new TCanvas("c9","Pion Spectra",400,10,600,700);
   c9->Divide(2,2);
   
   c9->cd(1);
   pionspectra1->SetFillColor(5);
   pionspectra1->SetXTitle("log(GeV)");
   pionspectra2->SetFillColor(46);
   pionspectra2->SetXTitle("log(GeV)");
   pionspectra3->SetFillColor(10);
   pionspectra3->SetXTitle("log(GeV)");
   //c9->SetLogx();
   pionspectra1->Draw();
   pionspectra2->Draw("same");
   pionspectra3->Draw("same");
   
   c9->cd(2);
   //TCanvas *c10 = new TCanvas("c10","Proton Spectra",400,10,600,700);
   protonspectra1->SetFillColor(5);
   protonspectra1->SetXTitle("log(GeV)");
   protonspectra2->SetFillColor(46);
   protonspectra2->SetXTitle("log(GeV)");
   protonspectra3->SetFillColor(10);
   protonspectra3->SetXTitle("log(GeV)");
   //c10->SetLogx();
   protonspectra1->Draw();
   protonspectra2->Draw("same");
   protonspectra3->Draw("same");
   
   c9->cd(3);
   //TCanvas *c11 = new TCanvas("c11","Kaon Spectra",400,10,600,700); 
   kaonspectra1->SetFillColor(5);
   kaonspectra1->SetXTitle("log(GeV)");
   kaonspectra2->SetFillColor(46);
   kaonspectra2->SetXTitle("log(GeV)");
   kaonspectra3->SetFillColor(10);
   kaonspectra3->SetXTitle("log(GeV)");
   //c11->SetLogx();
   kaonspectra1->Draw();
   kaonspectra2->Draw("same");
   kaonspectra3->Draw("same");
   
   c9->cd(4);
   //TCanvas *c12 = new TCanvas("c12","Charged Particles Spectra",400,10,600,700);
   chargedspectra1->SetFillColor(5);
   chargedspectra1->SetXTitle("log(GeV)");
   chargedspectra2->SetFillColor(46);
   chargedspectra2->SetXTitle("log(GeV)");
   chargedspectra3->SetFillColor(10);
   chargedspectra3->SetXTitle("log(GeV)");
   //c12->SetLogx();
   chargedspectra1->Draw();
   chargedspectra2->Draw("same");
   chargedspectra3->Draw("same");
   


   printf("*****************************************\n");
   printf("* Particle                   *  Counts  *\n");
   printf("*****************************************\n");

   printf("* Pions:                     *   %4d   *\n",pion);
   printf("* Charged Pions:             *   %4d   *\n",chargedpions);
   printf("* Primary Pions:             *   %4d   *\n",primarypions);
   printf("* Primary Pions (p>1GeV/c):  *   %4d   *\n",highprimarypions);
   printf("* Kaons:                     *   %4d   *\n",kaon);
   printf("* Charged Kaons:             *   %4d   *\n",chargedkaons);
   printf("* Primary Kaons:             *   %4d   *\n",primarykaons);
   printf("* Primary Kaons (p>1GeV/c):  *   %4d   *\n",highprimarykaons);
   printf("* Muons:                     *   %4d   *\n",muon);
   printf("* Electrons:                 *   %4d   *\n",electron);
   printf("* Positrons:                 *   %4d   *\n",positron);
   printf("* Protons:                   *   %4d   *\n",proton);
   printf("* All Charged:               *   %4d   *\n",(chargedpions+chargedkaons+muon+electron+positron+proton));
   printf("*****************************************\n");
   //printf("* Photons:                   *   %3.1f   *\n",photons); 
   //printf("* Primary Photons:           *   %3.1f   *\n",primaryphotons);
   //printf("* Primary Photons (p>1MeV/c):*   %3.1f   *\n",highprimaryphotons);
   //printf("*****************************************\n");
   //printf("* Neutrons:                  *   %3.1f   *\n",neutron);
   //printf("* Neutrons (p>100keV/c):     *   %3.1f   *\n",highneutrons);
   //printf("*****************************************\n");

   if (gAlice->TreeD())
     {
       gAlice->TreeD()->GetEvent(0);
   
       Float_t occ[7]; 
       Float_t sum=0;
       Float_t mean=0; 
       printf("\n*****************************************\n");
       printf("* Chamber   * Digits      * Occupancy   *\n");
       printf("*****************************************\n");
       
       for (Int_t ich=0;ich<7;ich++)
	 {
	   TClonesArray *Digits = DigitsAddress(ich);    //  Raw clusters branch
	   Int_t ndigits = Digits->GetEntriesFast();
	   occ[ich] = Float_t(ndigits)/(160*144);
	   sum += Float_t(ndigits)/(160*144);
	   printf("*   %d      *    %d      *   %3.1f%%     *\n",ich,ndigits,occ[ich]*100);
	 }
       mean = sum/7;
       printf("*****************************************\n");
       printf("* Mean occupancy          *   %3.1f%%     *\n",mean*100);
       printf("*****************************************\n");
     }
 
  printf("\nEnd of analysis\n");
   
}//void AliRICHv3::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
//__________________________________________________________________________________________________

