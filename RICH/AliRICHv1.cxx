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


#include <TParticle.h> 
#include <TRandom.h> 
#include <TVirtualMC.h>
#include <TPDGCode.h>

#include "AliRICHv1.h"
#include "AliRICHParam.h"
#include <AliConst.h>
#include <AliPDG.h>
#include <AliRun.h>

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
  static Int_t   iCurrentChamber;
  static Int_t   vol[2];
  Int_t          ipart;
  static Float_t hits[22];
  static Float_t ckovData[19];
  TLorentzVector x4,p4;
  Float_t        pos[3],mom[4],localPos[3],localMom[4];
  Float_t        theta,phi,localTheta,localPhi;
  Float_t        destep, step;
  Double_t        ranf[2];
  Int_t          nPads=kBad;
  Float_t        coscerenkov;
  static Float_t eloss,  tlength;
  const  Float_t kBig=1.e10;
       
  TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->GetCurrentTrackNumber()];

    
 
  id=gMC->CurrentVolID(copy); iCurrentChamber=copy;
  Float_t cherenkovLoss=0;
    
  gMC->TrackPosition(x4);    pos[0]=x4(0);    pos[1]=x4(1);    pos[2]=x4(2);
  ckovData[1] = pos[0];   ckovData[2] = pos[1];   ckovData[3] = pos[2]; 
  ckovData[6] = 0;                      // dummy track length
    
    /********************Store production parameters for Cerenkov photons************************/ 

  if(gMC->TrackPid()==kCerenkov){//C
    Float_t ckovEnergy = current->Energy();
    if(ckovEnergy > 5.6e-09 && ckovEnergy < 7.8e-09 ){//C+E
      if(gMC->IsTrackEntering()){//C+E+enter
        if(gMC->VolId("FRE1")==gMC->CurrentVolID(copy)||gMC->VolId("FRE2")==gMC->CurrentVolID(copy)){//C+E+enter+FRE
          if(gMC->IsNewTrack()){//C+E+enter+FRE+new
            Int_t mother=current->GetFirstMother(); 
            ckovData[10]=mother;
            ckovData[11]=gAlice->GetCurrentTrackNumber();
            ckovData[12]=1;             //Media where photon was produced 1->Freon, 2->Quarz
            fCkovNumber++;
            fFreonProd=1;
          }//C+E+enter+FRE+new
        }//C+E+enter+FRE
        if(gMC->IsNewTrack()&&gMC->VolId("QUAR")==gMC->CurrentVolID(copy))  ckovData[12]=2;		
      }//C+E+enter	      
      if(ckovData[12]==1){//C+E+produced in Freon
        if(gMC->IsTrackEntering()){                                     //is track entering?
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy)){                //is it in gap?      
                  gMC->TrackMomentum(p4);   mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);   mom[3]=p4(3);
			
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
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy)){             //is it in csi?      
		      
			gMC->TrackMomentum(p4);	mom[0]=p4(0);	mom[1]=p4(1);	mom[2]=p4(2);	mom[3]=p4(3);

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
      }//C+E+produced in Freon
    }//C+E
  }//C
//*************************************End of Production Parameters Storing*********************
    
  if(gMC->TrackPid()==kCerenkov||gMC->TrackPid()==kFeedback){//CF      
    if(gMC->VolId("CSI ")==gMC->CurrentVolID(copy)){//CF+CSI
      if(gMC->Edep()>0.){//CF+CSI+DE
        gMC->TrackPosition(x4);   pos[0]=x4(0);   pos[1]=x4(1);   pos[2]=x4(2);
        gMC->TrackMomentum(p4);   mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);   mom[3]=p4(3);
        
	Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
	Double_t rt = TMath::Sqrt(tc);
	theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
	phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		
	gMC->CurrentVolOffID(2,copy);	vol[0]=copy;	iCurrentChamber=vol[0];
	
        gMC->Gmtod(pos,localPos,1);     gMC->Gmtod(mom,localMom,2);

	Param()->SigGenInit(localPos[0],localPos[2]);
        ckovData[0]=gMC->TrackPid();   
        ckovData[1]=pos[0];      ckovData[2]=pos[1];    ckovData[3]=pos[2];            
        ckovData[4]=theta;       ckovData[5]=phi;       //theta-phi angles of incidence 
        ckovData[8]=(Float_t) fNsdigits;      // first sdigit
        ckovData[9]=-1;                       // last pad hit
        ckovData[13]=4;                       // photon was detected
        ckovData[14]=mom[0];     ckovData[15]=mom[1];    ckovData[16]=mom[2];
		    
	destep = gMC->Edep();
	gMC->SetMaxStep(kBig);
	cherenkovLoss  += destep;
	ckovData[7]=cherenkovLoss;
		    
	GenerateFeedbacks(iCurrentChamber,cherenkovLoss);//CF+CSI+DE
		    		    
	if (fNsdigits > (Int_t)ckovData[8]) {
	  ckovData[8]= ckovData[8]+1;
          ckovData[9]= (Float_t) fNsdigits;
        }

        ckovData[17] = nPads;
        AliRICHhit *mipHit =  (AliRICHhit*) (fHits->UncheckedAt(0));
        if(mipHit){
          mom[0] = current->Px();   mom[1] = current->Py();   mom[2] = current->Pz();
          Float_t mipPx = mipHit->MomX();   Float_t mipPy = mipHit->MomY();   Float_t mipPz = mipHit->MomZ();
			
          Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
          Float_t rt = TMath::Sqrt(r);
          Float_t mipR = mipPx*mipPx + mipPy*mipPy + mipPz*mipPz;	
          Float_t mipRt = TMath::Sqrt(mipR);
          if((rt*mipRt) > 0)
            coscerenkov = (mom[0]*mipPx + mom[1]*mipPy + mom[2]*mipPz)/(rt*mipRt);
          else
            coscerenkov = 0;
          ckovData[18]=TMath::ACos(coscerenkov);//Cerenkov angle
        }
        AddHit(gAlice->GetCurrentTrackNumber(),vol,ckovData);//PHOTON HIT CF+CSI+DE
        AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
      }//CF+CSI+DE
    }//CF+CSI
  }/*CF*/else if(gMC->TrackCharge()){//MIP
    if(gMC->VolId("FRE1")==gMC->CurrentVolID(copy)||gMC->VolId("FRE2")==gMC->CurrentVolID(copy)){//MIP+FRE
      gMC->TrackMomentum(p4);	    mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);   mom[3]=p4(3);
      hits[19]=mom[0];	    hits [20] = mom[1];	    hits [21] = mom[2];	    fFreonProd=1;
    }//MIP+FRE
    if(gMC->VolId("GAP ")==gMC->CurrentVolID(copy)){//MIP+GAP
      gMC->CurrentVolOffID(3,copy);    vol[0]=copy;    iCurrentChamber=vol[0];
      gMC->TrackPosition(x4);   pos[0]=x4(0);   pos[1]=x4(1);   pos[2]=x4(2);
      gMC->TrackMomentum(p4);   mom[0]=p4(0);   mom[1]=p4(1);   mom[2]=p4(2);    mom[3]=p4(3);
      gMC->Gmtod(pos,localPos,1);	    gMC->Gmtod(mom,localMom,2);
      ipart =gMC->TrackPid();
      destep = gMC->Edep();step   = gMC->TrackStep();// momentum loss and steplength in last step
      if(gMC->IsTrackEntering()){//MIP+GAP+Enter  record hit when mip enters ...
        Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
        Double_t rt = TMath::Sqrt(tc);
        theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
        phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
	Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
	Double_t localRt = TMath::Sqrt(localTc);
	localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])))*kRaddeg;                       
	localPhi     = Float_t(TMath::ATan2(Double_t(localMom[2]),Double_t(localMom[0])))*kRaddeg;    
        hits[0] = Float_t(ipart);         // particle type
        hits[1] = localPos[0];    hits[2] = localPos[1];     hits[3] = localPos[2];                 
        hits[4] = localTheta;     hits[5] = localPhi;              // theta-phi angles of incidence 
        hits[8] = (Float_t) fNsdigits;    // first sdigit
        hits[9] = -1;                     // last pad hit
        hits[13] = fFreonProd;           // did id hit the freon?
        hits[14] = mom[0];        hits[15] = mom[1];        hits[16] = mom[2];
        hits[18] = 0;               // dummy cerenkov angle
        tlength = 0;        eloss   = 0;        fFreonProd = 0;
        C(iCurrentChamber)->LocaltoGlobal(localPos,hits+1);
          
        Param()->SigGenInit(localPos[0], localPos[2]);
      }/*MIP+GAP+Enter*/else if(gMC->IsTrackExiting()||gMC->IsTrackStop()||gMC->IsTrackDisappeared()){//MIP+GAP+Exit
        gMC->SetMaxStep(kBig);
	eloss   += destep;
	tlength += step;
        if (eloss > 0) {
          if(gMC->TrackPid() == kNeutron) printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
          GenerateFeedbacks(iCurrentChamber,eloss);//MIP+GAP+Exit
          hits[17] = nPads;
        }
        hits[6]=tlength;        hits[7]=eloss;
        if(fNsdigits > (Int_t)hits[8]) {
          hits[8]= hits[8]+1;
          hits[9]= (Float_t) fNsdigits;
        }
	AddHit(gAlice->GetCurrentTrackNumber(),vol,hits);//MIP HIT MIP+GAP+Exit
	eloss = 0; 
      }/*MIP+GAP+Exit*/else if(Param()->SigGenCond(localPos[0], localPos[2])){//MIP+GAP+Spec
        Param()->SigGenInit(localPos[0], localPos[2]);
        if(eloss>0){
          if(gMC->TrackPid() == kNeutron) printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
          GenerateFeedbacks(iCurrentChamber,eloss);//MIP+GAP+Spec
          hits[17] = nPads;
        }
        
	eloss    = destep;	tlength += step ;		
      }/*MIP+GAP+Spec*/else{//MIP+GAP+nothing special
        eloss   += destep;
        tlength += step ;
      }//MIP+GAP+nothing special
    }//MIP+GAP
  }//MIP
}//void AliRICHv1::StepManager()
