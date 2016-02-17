// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************


#include "AliHMPIDvG4.h"         // class header
#include "AliHMPIDParam.h"      // StepManager()
#include "AliHMPIDHit.h"        // Hits2SDigs(),StepManager()
#include "AliHMPIDDigit.h"      // Digits2Raw(), Raw2SDigits()
#include "AliHMPIDRawStream.h"  // Digits2Raw(), Raw2SDigits()
#include "AliRawReader.h"       // Raw2SDigits()
#include "AliTrackReference.h"  //
#include <TVirtualMC.h>         // StepManager() for TVirtualMC::GetMC()
#include <TPDGCode.h>           // StepHistory() 
#include <AliStack.h>           // StepManager(),Hits2SDigits()78.6
#include <AliLoader.h>          // Hits2SDigits()
#include <AliRunLoader.h>       // Hits2SDigits()
#include <AliMC.h>              // StepManager()      
#include <AliRun.h>             // CreateMaterials()    
#include <AliMagF.h>            // CreateMaterials()
#include "AliGeomManager.h"     // AddAlignableVolumes()
#include <AliCDBEntry.h>        // CreateMaterials()
#include <AliCDBManager.h>      // CreateMaterials()
#include <TF1.h>                // DefineOpticalProperties()
#include <TF2.h>                // DefineOpticalProperties()
#include <TGeoCompositeShape.h> // CradleBaseVolume()
#include <TGeoGlobalMagField.h> //
#include <TGeoPhysicalNode.h>   // AddAlignableVolumes()
#include <TGeoXtru.h>           // CradleBaseVolume()
#include <TLorentzVector.h>     // IsLostByFresnel() 
#include <TString.h>            // StepManager()
#include <TTree.h>              // 

ClassImp(AliHMPIDvG4)    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDvG4::GenFee(Float_t qtot)
{
// Generate FeedBack photons for the current particle. To be invoked from StepManager().
// eloss=0 means photon so only pulse height distribution is to be analysed.
  TLorentzVector x4;
  TVirtualMC::GetMC()->TrackPosition(x4); 
  Int_t iNphotons=TVirtualMC::GetMC()->GetRandom()->Poisson(0.02*qtot);  //# of feedback photons is proportional to the charge of hit
  AliDebug(1,Form("N photons=%i",iNphotons));
  if (iNphotons > fMaxFeed) return;
  Int_t j;
  Float_t cthf, phif, enfp = 0, sthf, e1[3], e2[3], e3[3], vmod, uswop,dir[3], phi,pol[3], mom[4];
//Generate photons
  for(Int_t i=0;i<iNphotons;i++){//feedbacks loop
    Double_t ranf[2];
    TVirtualMC::GetMC()->GetRandom()->RndmArray(2,ranf);    //Sample direction
    cthf=ranf[0]*2-1.0;
    if(cthf<0) continue;
    sthf = TMath::Sqrt((1. - cthf) * (1. + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    
    if(Double_t randomNumber=TVirtualMC::GetMC()->GetRandom()->Rndm()<=0.57)
      enfp = 7.5e-9;
    else if(randomNumber<=0.7)
      enfp = 6.4e-9;
    else
      enfp = 7.9e-9;
    

    dir[0] = sthf * TMath::Sin(phif);    dir[1] = cthf;    dir[2] = sthf * TMath::Cos(phif);
    TVirtualMC::GetMC()->Gdtom(dir, mom, 2);
    mom[0]*=enfp;    mom[1]*=enfp;    mom[2]*=enfp;
    mom[3] = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    
    // Polarisation
    e1[0]=      0;    e1[1]=-dir[2];    e1[2]= dir[1];
    e2[0]=-dir[1];    e2[1]= dir[0];    e2[2]=      0;
    e3[0]= dir[1];    e3[1]=      0;    e3[2]=-dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;  for(j=0;j<3;j++) vmod+=e1[j]*e1[j];  vmod=TMath::Sqrt(1/vmod);  for(j=0;j<3;j++) e1[j]*=vmod;    
    vmod=0;  for(j=0;j<3;j++) vmod+=e2[j]*e2[j];  vmod=TMath::Sqrt(1/vmod);  for(j=0;j<3;j++) e2[j]*=vmod;
    
    phi = TVirtualMC::GetMC()->GetRandom()->Rndm()* 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    TVirtualMC::GetMC()->Gdtom(pol, pol, 2);
    Int_t outputNtracksStored;    
    gAlice->GetMCApp()->PushTrack(1,                             //transport
                     gAlice->GetMCApp()->GetCurrentTrackNumber(),//parent track 
                     50000051,                                   //PID
		     mom[0],mom[1],mom[2],mom[3],                //track momentum  
                     x4.X(),x4.Y(),x4.Z(),x4.T(),                //track origin 
                     pol[0],pol[1],pol[2],                       //polarization
		     kPFeedBackPhoton,                           //process ID   
                     outputNtracksStored,                        //on return how many new photons stored on stack
                     1.0);                                       //weight
  }//feedbacks loop
  AliDebug(1,"Stop.");
}//GenerateFeedbacks()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AliHMPIDvG4::StepManager()
{
// Full Step Manager.
// Arguments: none
//   Returns: none           
//  StepHistory(); return; //uncomment to print tracks history
 //  StepCount(); return;     //uncomment to count photons
  
   TString volname = TVirtualMC::GetMC()->CurrentVolName();

//Treat photons    
    if((TVirtualMC::GetMC()->TrackPid()==50000050||TVirtualMC::GetMC()->TrackPid()==50000051)&&volname.Contains("Hpad")){ //photon (Ckov or feedback) hits on module PC (Hpad)
    if(TVirtualMC::GetMC()->Edep()>0){                                                                           //photon survided QE test i.e. produces electron
      if(IsLostByFresnel()){ TVirtualMC::GetMC()->StopTrack(); return;}                                          //photon lost due to fersnel reflection on PC       
      Int_t   tid=     TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();                                 //take TID
      Int_t   pid=     TVirtualMC::GetMC()->TrackPid();                                                          //take PID
      Float_t etot=    TVirtualMC::GetMC()->Etot();                                                              //total hpoton energy, [GeV] 
      Double_t x[3];   TVirtualMC::GetMC()->TrackPosition(x[0],x[1],x[2]);                                       //take MARS position at entrance to PC
      Float_t hitTime= (Float_t)TVirtualMC::GetMC()->TrackTime();                                                //hit formation time       
      TString tmpname = volname; tmpname.Remove(0,4); Int_t idch = tmpname.Atoi();               //retrieve the chamber number
      Float_t xl,yl;   AliHMPIDParam::Instance()->Mars2Lors(idch,x,xl,yl);                       //take LORS position 
      new((*fHits)[fNhits++])AliHMPIDHit(idch,etot,pid,tid,xl,yl,hitTime,x);                             //HIT for photon, position at P, etot will be set to Q
      if(fDoFeed) {
	Int_t nfeed = etot * 0.02;
	if (nfeed > fMaxFeed) {
	  printf("Nfeed: %5d eloss: %13.3f pid: %5d tid: %5d \n", nfeed, etot, pid, tid);
	  StepHistory();
	}
	GenFee(etot);                                                                  //generate feedback photons etot is modified in hit ctor to Q of hit
      }
    }//photon hit PC and DE >0 
  }//photon hit PC
   
  
  //Treat charged particles  
  static Float_t eloss;                                                                           //need to store mip parameters between different steps    
  static Double_t in[3];                                                                          

  if(TVirtualMC::GetMC()->IsTrackEntering() && TVirtualMC::GetMC()->TrackCharge() && volname.Contains("Hpad")) //Trackref stored when entering in the pad volume
    AddTrackReference(TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber(), AliTrackReference::kHMPID);       //for acceptance calculations
   

  if(TVirtualMC::GetMC()->TrackCharge() && volname.Contains("Hcel")){                                           //charged particle in amplification gap (Hcel)
    if(TVirtualMC::GetMC()->IsTrackEntering()||TVirtualMC::GetMC()->IsNewTrack()) {                                               //entering or newly created
      eloss=0;                                                                                    //reset Eloss collector                         
      TVirtualMC::GetMC()->TrackPosition(in[0],in[1],in[2]);                                                      //take position at the entrance
    }else if(TVirtualMC::GetMC()->IsTrackExiting()||TVirtualMC::GetMC()->IsTrackStop()||TVirtualMC::GetMC()->IsTrackDisappeared()){               //exiting or disappeared
      eloss              +=TVirtualMC::GetMC()->Edep();                                                           //take into account last step Eloss
      Int_t tid=          TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();                               //take TID
      Int_t pid=          TVirtualMC::GetMC()->TrackPid();                                                        //take PID
      Double_t out[3];    TVirtualMC::GetMC()->TrackPosition(out[0],out[1],out[2]);                               //take MARS position at exit
      Float_t hitTime= (Float_t)TVirtualMC::GetMC()->TrackTime();                                                         //hit formation time       
      out[0]=0.5*(out[0]+in[0]);                                                                  //
      out[1]=0.5*(out[1]+in[1]);                                                                  //take hit position at the anod plane
      out[2]=0.5*(out[2]+in[2]);
      TString tmpname = volname;  tmpname.Remove(0,4);  Int_t idch = tmpname.Atoi();              //retrieve the chamber number
      Float_t xl,yl;AliHMPIDParam::Instance()->Mars2Lors(idch,out,xl,yl);                         //take LORS position
      if(eloss>0) {
        new((*fHits)[fNhits++])AliHMPIDHit(idch,eloss,pid,tid,xl,yl,hitTime,out);                           //HIT for MIP, position near anod plane, eloss will be set to Q 
        if(fDoFeed){
	  Int_t nfeed = 0.02 * eloss;
	  if (nfeed > fMaxFeed) {
	    printf("Nfeed: %5d eloss: %13.3f pid: %5d tid: %5d \n", nfeed, eloss, pid, tid);
	    StepHistory();
	  }
	  GenFee(eloss);                                                                  //generate feedback photons 
	}
        eloss=0;
      }
    }else                                                                                         //just going inside
      eloss          += TVirtualMC::GetMC()->Edep();                                                              //collect this step eloss 
  }//MIP in GAP
 
}//StepManager()
