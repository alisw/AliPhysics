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
  Revision 1.3  2000/06/13 13:06:38  jbarbosa
  Fixed compiling error for HP (multiple declaration)

  Revision 1.2  2000/06/12 15:36:16  jbarbosa
  Cleaned up version.

  Revision 1.1  2000/06/09 15:00:31  jbarbosa
  New full version. All parameters configurable.

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
#include <TParticle.h> 

#include "AliRICHv1.h"
#include "AliRICHHit.h"
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

// Default constructor fo AliRICHvv1 (full version)

    //fChambers = 0;
}

//___________________________________________
AliRICHv1::AliRICHv1(const char *name, const char *title)
    : AliRICHv0(name,title)
{

// Full version of RICH with hits and diagnostics

    fCkovNumber=0;
    fFreonProd=0;
  
    fChambers = new TObjArray(kNCH);
    for (Int_t i=0; i<kNCH; i++) {
    
	(*fChambers)[i] = new AliRICHChamber();  
	
    }
}

void AliRICHv1::Init()
{

  printf("*********************************** RICH_INIT ***********************************\n");
  printf("*                                                                               *\n");
  printf("*                    AliRICHv1 Configurable version started                     *\n");
  printf("*                                                                               *\n");

  
  AliRICHSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHResponse*  response;


    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=1; i<kNCH; i++) {
	//printf ("i:%d",i);
	( (AliRICHChamber*) (*fChambers)[i])->Init();  
    }  
    
    //
    // Set the chamber (sensitive region) GEANT identifier
    
    ((AliRICHChamber*)(*fChambers)[0])->SetGid(1);  
    ((AliRICHChamber*)(*fChambers)[1])->SetGid(2);  
    ((AliRICHChamber*)(*fChambers)[2])->SetGid(3);  
    ((AliRICHChamber*)(*fChambers)[3])->SetGid(4);  
    ((AliRICHChamber*)(*fChambers)[4])->SetGid(5);  
    ((AliRICHChamber*)(*fChambers)[5])->SetGid(6);  
    ((AliRICHChamber*)(*fChambers)[6])->SetGid(7); 

    Float_t pos1[3]={0,471.8999,165.2599};
    Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90));

    Float_t pos2[3]={171,470,0};
    Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],new TRotMatrix("rot994","rot994",90,-20,90,70,0,0));

    Float_t pos3[3]={0,500,0};
    Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],new TRotMatrix("rot995","rot995",90,0,90,90,0,0));
    
    Float_t pos4[3]={-171,470,0};
    Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], new TRotMatrix("rot996","rot996",90,20,90,110,0,0));  

    Float_t pos5[3]={161.3999,443.3999,-165.3};
    Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70));

    Float_t pos6[3]={0., 471.9, -165.3,};
    Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90));

    Float_t pos7[3]={-161.399,443.3999,-165.3};
    Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110));
    
    segmentation=Chamber(0).GetSegmentationModel(0);
    geometry=Chamber(0).GetGeometryModel();
    response=Chamber(0).GetResponseModel();
    
     
    printf("*                            Pads            : %3dx%3d                          *\n",segmentation->Npx(),segmentation->Npy());
    printf("*                            Pad size        : %5.2f x%5.2f mm2                 *\n",segmentation->Dpx(),segmentation->Dpy()); 
    printf("*                            Gap Thickness   : %5.1f mm                         *\n",geometry->GetGapThickness());
    printf("*                            Radiator Width  : %5.1f mm                         *\n",geometry->GetQuartzWidth());
    printf("*                            Radiator Length : %5.1f mm                         *\n",geometry->GetQuartzLength());
    printf("*                            Freon Thickness : %5.1f mm                         *\n",geometry->GetFreonThickness());
    printf("*                            Charge Slope    : %5.1f ADC                        *\n",response->ChargeSlope());
    printf("*                            Feedback Prob.  : %5.2f %%                          *\n",response->AlphaFeedback());
    printf("*                                                                               *\n");
    printf("*                                   Success!                                    *\n");
    printf("*                                                                               *\n");
    printf("*********************************************************************************\n");

}

//___________________________________________
void AliRICHv1::StepManager()
{

// Full Step Manager

    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[18];
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
    Float_t        ranf[2];
    Int_t          nPads;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TGeant3 *geant3 = (TGeant3*) gMC;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->CurrentTrack()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkovLoss=0;
    //gAlice->KeepTrack(gAlice->CurrentTrack());
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    //ckovData[11] = gAlice->CurrentTrack();

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
	      if (gMC->IsTrackEntering()){                                     //is track entering?
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (geant3->Gctrak()->nstep<1){                          //is it the first step?
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->CurrentTrack();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      fCkovNumber++;
		      fFreonProd=1;
		      //printf("Index: %d\n",fCkovNumber);
		    }    //first step question
		  }        //freo question
		
		if (geant3->Gctrak()->nstep<1){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreonProd);
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Got in");
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			// Z-position for hit
			
			
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			gMC->Rndm(ranf, 1);
			//printf("grid calculation:%f\n",t);
			if (ranf[0] > t) {
			  //geant3->StopTrack();
			  ckovData[13] = 5;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = Fresnel(ckovEnergy*1e9,cophi,1);
			gMC->Rndm(ranf, 1);
			if (ranf[0] < t) {
			  //geant3->StopTrack();
			  ckovData[13] = 6;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
		      ckovData[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=2;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		    } //reflection question
		    
		    
		    //        Absorption loss 
		    else if (geant3->Gctrak()->lmec[i] == 101) {              //was it absorbed?
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
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		      //printf("Added cerenkov %d\n",fCkovNumber);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (geant3->Gctrak()->lmec[i] == 30) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
		gMC->Gmtod(pos,localPos,1);                                                                    
		gMC->Gmtod(mom,localMom,2);
		
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
		
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if(idvol<kNCH) {	
		    ckovData[0] = gMC->TrackPid();        // particle type
		    ckovData[1] = pos[0];                 // X-position for hit
		    ckovData[2] = pos[1];                 // Y-position for hit
		    ckovData[3] = pos[2];                 // Z-position for hit
		    ckovData[4] = theta;                      // theta angle of incidence
		    ckovData[5] = phi;                      // phi angle of incidence 
		    ckovData[8] = (Float_t) fNPadHits;      // first padhit
		    ckovData[9] = -1;                       // last pad hit
		    ckovData[13] = 4;                       // photon was detected
		    ckovData[14] = mom[0];
		    ckovData[15] = mom[1];
		    ckovData[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(kBig);
		    cherenkovLoss  += destep;
		    ckovData[7]=cherenkovLoss;
		    
		    nPads = MakePadHits(localPos[0],localPos[2],cherenkovLoss,idvol,kCerenkov);
		    if (fNPadHits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNPadHits;
		    }

		    ckovData[17] = nPads;
		    //printf("nPads:%d",nPads);
		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHHit *mipHit =  (AliRICHHit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t mipPx = mipHit->fMomX;
			Float_t mipPy = mipHit->fMomY;
			Float_t mipPz = mipHit->fMomZ;
			
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
		    AddHit(gAlice->CurrentTrack(),vol,ckovData);
		    AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
	    fFreonProd=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
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
	    gMC->Gmtod(mom,localMom,2);
	    
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
		hits[8] = (Float_t) fNPadHits;    // first padhit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreonProd;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];

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
		    ((AliRICHChamber*) (*fChambers)[idvol])
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
		      nPads = MakePadHits(xhit,yhit,eloss,idvol,kMip);
		      hits[17] = nPads;
		      //printf("nPads:%d",nPads);
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
		 ->SigGenCond(localPos[0], localPos[2], localPos[1]))
	    {
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    nPads = MakePadHits(xhit,yhit,eloss,idvol,kMip);
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
	}
      }
    /*************************************************End of MIP treatment**************************************/
   //}
}
