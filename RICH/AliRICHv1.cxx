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

/*
  $Log$
  Revision 1.9  2000/05/31 08:19:38  jbarbosa
  Fixed bug in StepManager

  Revision 1.8  2000/05/26 17:30:08  jbarbosa
  Cerenkov angle now stored within cerenkov data structure.

  Revision 1.7  2000/05/18 10:31:36  jbarbosa
  Fixed positioning of spacers inside freon.
  Fixed positioning of proximity gap
  inside methane.
  Fixed cut on neutral particles in the StepManager.

  Revision 1.6  2000/04/28 11:51:58  morsch
   Dimensions of arrays hits and Ckov_data corrected.

  Revision 1.5  2000/04/19 13:28:46  morsch
  Major changes in geometry (parametrised), materials (updated) and
  step manager (diagnostics) (JB, AM)

*/



//////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version //
//////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliRICHv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliPDG.h" 
#include "TGeant3.h"

ClassImp(AliRICHv1)
    
//___________________________________________
AliRICHv1::AliRICHv1() : AliRICHv0()
{
    //fChambers = 0;
}

//___________________________________________
AliRICHv1::AliRICHv1(const char *name, const char *title)
    : AliRICHv0(name,title)
{
    fCkov_number=0;
    fFreon_prod=0;

    fChambers = new TObjArray(7);
    for (Int_t i=0; i<7; i++) {
    
	(*fChambers)[i] = new AliRICHChamber();  
	
    }
}


//___________________________________________
void AliRICHv1::StepManager()
{
    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[18];
    static Float_t Ckov_data[19];
    TLorentzVector Position;
    TLorentzVector Momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        Localpos[3];
    Float_t        Localmom[4];
    Float_t        Localtheta,Localphi;
    Float_t        theta,phi;
    Float_t        destep, step;
    Float_t        ranf[2];
    Int_t          NPads;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t big=1.e10;
       
    TClonesArray &lhits = *fHits;
    TGeant3 *geant3 = (TGeant3*) gMC;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->CurrentTrack()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkov_loss=0;
    //gAlice->KeepTrack(gAlice->CurrentTrack());
    
    gMC->TrackPosition(Position);
    pos[0]=Position(0);
    pos[1]=Position(1);
    pos[2]=Position(2);
    Ckov_data[1] = pos[0];                 // X-position for hit
    Ckov_data[2] = pos[1];                 // Y-position for hit
    Ckov_data[3] = pos[2];                 // Z-position for hit
    //Ckov_data[11] = gAlice->CurrentTrack();

    //AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    /********************Store production parameters for Cerenkov photons************************/ 
//is it a Cerenkov photon? 
    if (gMC->TrackPid() == 50000050) {          

      //if (gMC->VolId("GAP ")==gMC->CurrentVolID(copy))
        //{                    
	  Float_t Ckov_energy = current->Energy();
	  //energy interval for tracking
	  if  (Ckov_energy > 5.6e-09 && Ckov_energy < 7.8e-09 )       
	    //if (Ckov_energy > 0)
	    {
	      if (gMC->IsTrackEntering()){                                     //is track entering?
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (geant3->Gctrak()->nstep<1){                          //is it the first step?
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      Ckov_data[10] = mother;
		      Ckov_data[11] = gAlice->CurrentTrack();
		      Ckov_data[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      fCkov_number++;
		      fFreon_prod=1;
		      //printf("Index: %d\n",fCkov_number);
		    }    //first step question
		  }        //freo question
		
		if (geant3->Gctrak()->nstep<1){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      Ckov_data[12] = 2;
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreon_prod);
	      }   //track entering question
	      
	      if (Ckov_data[12] == 1)                                        //was it produced in Freon?
		//if (fFreon_prod == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Got in");
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in\n");
			gMC->TrackMomentum(Momentum);
			mom[0]=Momentum(0);
			mom[1]=Momentum(1);
			mom[2]=Momentum(2);
			mom[3]=Momentum(3);
			// Z-position for hit
			
			
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			gMC->Rndm(ranf, 1);
			//printf("grid calculation:%f\n",t);
			if (ranf[0] > t) {
			  //geant3->StopTrack();
			  Ckov_data[13] = 5;
			  AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			gMC->TrackMomentum(Momentum);
			mom[0]=Momentum(0);
			mom[1]=Momentum(1);
			mom[2]=Momentum(2);
			mom[3]=Momentum(3);
			
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = Fresnel(Ckov_energy*1e9,cophi,1);
			gMC->Rndm(ranf, 1);
			if (ranf[0] < t) {
			  //geant3->StopTrack();
			  Ckov_data[13] = 6;
			  AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
			  //printf("Lost by Fresnel\n");
			}
			/**********************************************************************************/
		      }
		  } //track entering?
		  
		  
		  /********************Evaluation of losses************************/
		  /******************still in the old fashion**********************/
		  
		  Int_t i1 = geant3->Gctrak()->nmec;            //number of physics mechanisms acting on the particle
		  for (Int_t i = 0; i < i1; ++i) {
		    //        Reflection loss 
		    if (geant3->Gctrak()->lmec[i] == 106) {        //was it reflected
		      Ckov_data[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			Ckov_data[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			Ckov_data[13]=2;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		    } //reflection question
		    
		    
		    //        Absorption loss 
		    else if (geant3->Gctrak()->lmec[i] == 101) {              //was it absorbed?
		      Ckov_data[13]=20;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			Ckov_data[13]=11;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			Ckov_data[13]=12;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("META")) 
			Ckov_data[13]=13;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("GAP ")) 
			Ckov_data[13]=13;
		      
		      if (gMC->CurrentVolID(copy) == gMC->VolId("SRIC")) 
			Ckov_data[13]=15;
		      
		      //        CsI inefficiency 
		      if (gMC->CurrentVolID(copy) == gMC->VolId("CSI ")) {
			Ckov_data[13]=16;
		      }
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		      //printf("Added cerenkov %d\n",fCkov_number);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (geant3->Gctrak()->lmec[i] == 30) {                 //is it below energy treshold?
		      Ckov_data[13]=21;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
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
	if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	    
	  if (gMC->Edep() > 0.){
		gMC->TrackPosition(Position);
		gMC->TrackMomentum(Momentum);
		pos[0]=Position(0);
		pos[1]=Position(1);
		pos[2]=Position(2);
		mom[0]=Momentum(0);
		mom[1]=Momentum(1);
		mom[2]=Momentum(2);
		mom[3]=Momentum(3);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		gMC->Gmtod(pos,Localpos,1);                                                                    
		gMC->Gmtod(mom,Localmom,2);
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		//Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(Localpos[0], Localpos[2]);
		//printf("Sector:%d\n",sector);

		/*if (gMC->TrackPid() == 50000051){
		  fFeedbacks++;
		  printf("Feedbacks:%d\n",fFeedbacks);
		}*/	
		
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		if(idvol<7) {	
		    Ckov_data[0] = gMC->TrackPid();        // particle type
		    Ckov_data[1] = pos[0];                 // X-position for hit
		    Ckov_data[2] = pos[1];                 // Y-position for hit
		    Ckov_data[3] = pos[2];                 // Z-position for hit
		    Ckov_data[4] = theta;                      // theta angle of incidence
		    Ckov_data[5] = phi;                      // phi angle of incidence 
		    Ckov_data[8] = (Float_t) fNPadHits;      // first padhit
		    Ckov_data[9] = -1;                       // last pad hit
		    Ckov_data[13] = 4;                       // photon was detected
		    Ckov_data[14] = mom[0];
		    Ckov_data[15] = mom[1];
		    Ckov_data[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(big);
		    cherenkov_loss  += destep;
		    Ckov_data[7]=cherenkov_loss;
		    
		    NPads = MakePadHits(Localpos[0],Localpos[2],cherenkov_loss,idvol,cerenkov);
		    if (fNPadHits > (Int_t)Ckov_data[8]) {
			Ckov_data[8]= Ckov_data[8]+1;
			Ckov_data[9]= (Float_t) fNPadHits;
		    }

		    Ckov_data[17] = NPads;
		    //printf("Npads:%d",NPads);
		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHHit *mipHit =  (AliRICHHit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t Mip_px = mipHit->fMomX;
			Float_t Mip_py = mipHit->fMomY;
			Float_t Mip_pz = mipHit->fMomZ;
			
			Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			Float_t rt = TMath::Sqrt(r);
			Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
			Float_t Mip_rt = TMath::Sqrt(Mip_r);
			if ((rt*Mip_rt) > 0)
			  {
			    coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt);
			  }
			else
			  {
			    coscerenkov = 0;
			  }
			Float_t cherenkov = TMath::ACos(coscerenkov);
			Ckov_data[18]=cherenkov;
		      }
		    //if (sector != -1)
		    //{
		    AddHit(gAlice->CurrentTrack(),vol,Ckov_data);
		    AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge())
    //else if (1 == 1)
      {
//If MIP
	/*if (gMC->IsTrackEntering())
	  {                
	    hits[13]=20;//is track entering?
	  }*/
	if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
	  {
	    fFreon_prod=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;

	    //Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(Localpos[0], Localpos[2]);
	    //printf("Sector:%d\n",sector);
	    
	    gMC->TrackPosition(Position);
	    gMC->TrackMomentum(Momentum);
	    pos[0]=Position(0);
	    pos[1]=Position(1);
	    pos[2]=Position(2);
	    mom[0]=Momentum(0);
	    mom[1]=Momentum(1);
	    mom[2]=Momentum(2);
	    mom[3]=Momentum(3);
	    gMC->Gmtod(pos,Localpos,1);                                                                    
	    gMC->Gmtod(mom,Localmom,2);
	    
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
		

		Double_t Localtc = Localmom[0]*Localmom[0]+Localmom[2]*Localmom[2];
		Double_t Localrt = TMath::Sqrt(Localtc);
		Localtheta   = Float_t(TMath::ATan2(Localrt,Double_t(Localmom[1])))*kRaddeg;                       
		Localphi     = Float_t(TMath::ATan2(Double_t(Localmom[2]),Double_t(Localmom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = Localpos[0];                 // X-position for hit
		hits[2] = Localpos[1];                 // Y-position for hit
		hits[3] = Localpos[2];                 // Z-position for hit
		hits[4] = Localtheta;                  // theta angle of incidence
		hits[5] = Localphi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNPadHits;    // first padhit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreon_prod;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];

		tlength = 0;
		eloss   = 0;
		fFreon_prod = 0;
	
		Chamber(idvol).LocaltoGlobal(Localpos,hits+1);
	   
		
		//To make chamber coordinates x-y had to pass LocalPos[0], LocalPos[2]
		xhit    = Localpos[0];
		yhit    = Localpos[2];
		// Only if not trigger chamber
		if(idvol<7) {
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
		    ((AliRICHChamber*) (*fChambers)[idvol])
			->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		}
	    }
	    
	    // 
	    // Calculate the charge induced on a pad (disintegration) in case 
	    //
	    // Mip left chamber ...
	    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
		gMC->SetMaxStep(big);
		eloss   += destep;
		tlength += step;
		
				
		// Only if not trigger chamber
		if(idvol<7) {
		  if (eloss > 0) 
		    {
		      if(gMC->TrackPid() == kNeutron)
			printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		      NPads = MakePadHits(xhit,yhit,eloss,idvol,mip);
		      hits[17] = NPads;
		      //printf("Npads:%d",NPads);
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNPadHits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNPadHits;
		}
		
		//if(sector !=-1)
		new(lhits[fNhits++]) AliRICHHit(fIshunt,gAlice->CurrentTrack(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    } else if 
		(((AliRICHChamber*) (*fChambers)[idvol])
		 ->SigGenCond(Localpos[0], Localpos[2], Localpos[1]))
	    {
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    NPads = MakePadHits(xhit,yhit,eloss,idvol,mip);
		    hits[17] = NPads;
		    //printf("Npads:%d",NPads);
		  }
		xhit     = Localpos[0];
		yhit     = Localpos[2]; 
		eloss    = destep;
		tlength += step ;
		//
		// nothing special  happened, add up energy loss
	    } else {        
		eloss   += destep;
		tlength += step ;
	    }
	}
      }
    /*************************************************End of MIP treatment**************************************/
   //}
}
