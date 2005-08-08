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


#include "AliRICHv1.h"     //class header
#include "AliRICHParam.h"
#include "AliRICHChamber.h"
#include <TParticle.h> 
#include <TRandom.h> 
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include <AliStack.h>      //Hits2SDigits()
#include <AliConst.h>
#include <AliPDG.h>
#include <AliRun.h>
#include <AliMC.h>
#include <AliRawDataHeader.h> //Digits2Raw()

ClassImp(AliRICHv1)    
//__________________________________________________________________________________________________
void AliRICHv1::StepManager()
{
// Full Step Manager.
// 3- Ckovs absorbed on Collection electrods
// 5- Ckovs absorbed on Cathode wires
// 6- Ckovs absorbed on Anod wires
         
  Int_t          copy;
  static Int_t   iCurrentChamber;
  static Int_t idRRAD = gMC->VolId("RRAD");
  static Int_t idRRWI = gMC->VolId("RRWI");
  static Int_t idRICH = gMC->VolId("RICH");
  static Int_t idRPC  = gMC->VolId("RPC ");
  static Int_t idRGAP = gMC->VolId("RGAP");
//history of Cerenkovs
  if(gMC->TrackPid()==kCerenkov){
    if( gMC->IsNewTrack()   && gMC->CurrentVolID(copy)==idRRAD) fCounters(0)++;// 0- Ckovs produced in radiator
    if(!gMC->IsTrackAlive() && gMC->CurrentVolID(copy)==idRRAD) fCounters(1)++;// 1- Ckovs absorbed in radiator
    if(!gMC->IsTrackAlive() && gMC->CurrentVolID(copy)==idRRWI) fCounters(2)++;// 2- Ckovs absorbed in radiator window
    if(!gMC->IsTrackAlive() && gMC->CurrentVolID(copy)==idRICH) fCounters(4)++;// 4- Ckovs absorbed in CH4
  }
          
//Treat photons    
  static TLorentzVector cerX4;
  if((gMC->TrackPid()==kCerenkov||gMC->TrackPid()==kFeedback)&&gMC->CurrentVolID(copy)==idRPC){//photon in PC
    if(gMC->Edep()>0){//photon in PC +DE
      if(IsLostByFresnel()){ 
        if(gMC->TrackPid()==kCerenkov) fCounters(7)++;// 7- Ckovs reflected from CsI
        gMC->StopTrack();
        return;
      }        
      gMC->TrackPosition(cerX4); gMC->CurrentVolOffID(2,iCurrentChamber);//RICH-RPPF-RPC
	
      HitAdd(iCurrentChamber,gAlice->GetMCApp()->GetCurrentTrackNumber(),cerX4.Vect(),cerX4.Vect());//HIT for PHOTON in conditions CF+CSI+DE
      fCounters(8)++;//4- Ckovs converted to electron on CsI
      GenerateFeedbacks(iCurrentChamber);
    }//photon in PC and DE >0 
  }//photon in PC
  
//Treat charged particles  
  static Float_t eloss;
  static TLorentzVector mipInX4,mipOutX4;
  if(gMC->TrackCharge() && gMC->CurrentVolID(copy)==idRGAP){//MIP in GAP
    gMC->CurrentVolOffID(1,iCurrentChamber);//RICH-RGAP
    if(gMC->IsTrackEntering()||gMC->IsNewTrack()) {//MIP in GAP entering or newly created
      eloss=0;                                                           
      gMC->TrackPosition(mipInX4);
    }else if(gMC->IsTrackExiting()||gMC->IsTrackStop()||gMC->IsTrackDisappeared()){//MIP in GAP exiting or disappeared
      eloss+=gMC->Edep();//take into account last step dEdX
      gMC->TrackPosition(mipOutX4);  
      HitAdd(iCurrentChamber,gAlice->GetMCApp()->GetCurrentTrackNumber(),mipInX4.Vect(),mipOutX4.Vect(),eloss);//HIT for MIP: MIP in GAP Exiting
      GenerateFeedbacks(iCurrentChamber,eloss);//MIP+GAP+Exit
    }else//MIP in GAP going inside
      eloss   += gMC->Edep();
  }//MIP in GAP
}//StepManager()
//__________________________________________________________________________________________________
Bool_t AliRICHv1::IsLostByFresnel()
{
// Calculate probability for the photon to be lost by Fresnel reflection.
  TLorentzVector p4;
  Double_t mom[3],localMom[3];
  gMC->TrackMomentum(p4);   mom[0]=p4(1);   mom[1]=p4(2);   mom[2]=p4(3);
  localMom[0]=0; localMom[1]=0; localMom[2]=0;
  gMC->Gmtod(mom,localMom,2);
  Double_t localTc    = localMom[0]*localMom[0]+localMom[2]*localMom[2];
  Double_t localTheta = TMath::ATan2(TMath::Sqrt(localTc),localMom[1]);
  Double_t cotheta = TMath::Abs(TMath::Cos(localTheta));
  if(gMC->GetRandom()->Rndm() < Fresnel(p4.E()*1e9,cotheta,1)){
    AliDebug(1,"Photon lost");
    return kTRUE;
  }else
    return kFALSE;
}//IsLostByFresnel()
//__________________________________________________________________________________________________
void AliRICHv1::GenerateFeedbacks(Int_t iChamber,Float_t eloss)
{
// Generate FeedBack photons for the current particle. To be invoked from StepManager().
// eloss=0 means photon so only pulse height distribution is to be analysed. This one is done in AliRICHParam::TotQdc()
  
  TLorentzVector x4;
  gMC->TrackPosition(x4);  
  TVector2 x2=C(iChamber)->Mrs2Pc(x4);//hit position on photocathode plane
  TVector2 xspe=x2;
  Int_t sector=P()->Loc2Sec(xspe);  if(sector==-1) return; //hit in dead zone, nothing to produce
  Int_t iTotQdc=P()->TotQdc(x2,eloss);
  Int_t iNphotons=gMC->GetRandom()->Poisson(P()->C(iChamber)->AlphaFeedback(sector)*iTotQdc);    
  AliDebug(1,Form("N photons=%i",iNphotons));
  Int_t j;
  Float_t cthf, phif, enfp = 0, sthf, e1[3], e2[3], e3[3], vmod, uswop,dir[3], phi,pol[3], mom[4];
//Generate photons
  for(Int_t i=0;i<iNphotons;i++){//feedbacks loop
    Double_t ranf[2];
    gMC->GetRandom()->RndmArray(2,ranf);    //Sample direction
    cthf=ranf[0]*2-1.0;
    if(cthf<0) continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    
    if(Double_t randomNumber=gMC->GetRandom()->Rndm()<=0.57)
      enfp = 7.5e-9;
    else if(randomNumber<=0.7)
      enfp = 6.4e-9;
    else
      enfp = 7.9e-9;
    

    dir[0] = sthf * TMath::Sin(phif);    dir[1] = cthf;    dir[2] = sthf * TMath::Cos(phif);
    gMC->Gdtom(dir, mom, 2);
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
    
    phi = gMC->GetRandom()->Rndm()* 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    Int_t outputNtracksStored;    
    gAlice->GetMCApp()->PushTrack(1,                 //do not transport
                     gAlice->GetMCApp()->GetCurrentTrackNumber(),//parent track 
                     kFeedback,                      //PID
		     mom[0],mom[1],mom[2],mom[3],    //track momentum  
                     x4.X(),x4.Y(),x4.Z(),x4.T(),    //track origin 
                     pol[0],pol[1],pol[2],           //polarization
		     kPFeedBackPhoton,
                     outputNtracksStored,
                     1.0);    
  }//feedbacks loop
  AliDebug(1,"Stop.");
}//GenerateFeedbacks()
//__________________________________________________________________________________________________
void AliRICHv1::Hits2SDigits()
{
// Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
  AliDebug(1,"Start.");
  for(Int_t iEventN=0;iEventN<GetLoader()->GetRunLoader()->GetAliRun()->GetEventsPerRun();iEventN++){//events loop
    GetLoader()->GetRunLoader()->GetEvent(iEventN);//get next event
  
    if(!GetLoader()->TreeH()) GetLoader()->LoadHits();    GetLoader()->GetRunLoader()->LoadHeader(); 
    if(!GetLoader()->GetRunLoader()->TreeK())             GetLoader()->GetRunLoader()->LoadKinematics();//from
    if(!GetLoader()->TreeS()) GetLoader()->MakeTree("S"); MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop 
        AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);//get current hit                
        TVector2 x2 = C(pHit->C())->Mrs2Anod(0.5*(pHit->InX3()+pHit->OutX3()));//hit position in the anod plane
        Int_t iTotQdc=P()->TotQdc(x2,pHit->Eloss());//total charge produced by hit, 0 if hit in dead zone
        if(iTotQdc==0) continue;
        //
        //need to quantize the anod....
        TVector padHit=AliRICHParam::Loc2Pad(x2);
        TVector2 padHitXY=AliRICHParam::Pad2Loc(padHit);
        TVector2 anod;
        if((x2.Y()-padHitXY.Y())>0) anod.Set(x2.X(),padHitXY.Y()+AliRICHParam::PitchAnod()/2);
        else anod.Set(x2.X(),padHitXY.Y()-AliRICHParam::PitchAnod()/2);
        //end to quantize anod
        //
        TVector area=P()->Loc2Area(anod);//determine affected pads, dead zones analysed inside
        AliDebug(1,Form("hitanod(%6.2f,%6.2f)->area(%3.0f,%3.0f)-(%3.0f,%3.0f) QDC=%4i",anod.X(),anod.Y(),area[0],area[1],area[2],area[3],iTotQdc));
        TVector pad(2);
        for(pad[1]=area[1];pad[1]<=area[3];pad[1]++)//affected pads loop
          for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){                    
            Double_t padQdc=iTotQdc*P()->FracQdc(anod,pad);
            AliDebug(1,Form("current pad(%3.0f,%3.0f) with QDC  =%6.2f",pad[0],pad[1],padQdc));
            if(padQdc>0.1) SDigitAdd(pHit->C(),pad,padQdc,GetLoader()->GetRunLoader()->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }//affected pads loop 
      }//hits loop
    }//prims loop
    GetLoader()->TreeS()->Fill();
    GetLoader()->WriteSDigits("OVERWRITE");
    SDigitsReset();
  }//events loop  
  GetLoader()->UnloadHits(); GetLoader()->GetRunLoader()->UnloadHeader(); GetLoader()->GetRunLoader()->UnloadKinematics();
  GetLoader()->UnloadSDigits();  
  AliDebug(1,"Stop.");
}//Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICHv1::Digits2Raw()
{
//Creates raw data files in DDL format. Invoked by AliSimulation
//loop over events is done outside in AliSimulation
//Arguments: none
//  Returns: none    
  AliDebug(1,"Start.");
  GetLoader()->LoadDigits();
  GetLoader()->TreeD()->GetEntry(0);
  
  ofstream file[AliRICHDigit::kNddls];   //output streams
  Int_t    cnt[AliRICHDigit::kNddls];        //data words counters for DDLs
  AliRawDataHeader header; //empty DDL header
  UInt_t w32=0;            //32 bits data word 
  
  for(Int_t i=0;i<AliRICHDigit::kNddls;i++){//open all 14 DDL in parallel
    file[i].open(Form("RICH_%4i.ddl",AliRICHDigit::kDdlOffset+i));  //first is 0x700
    file[i].write((char*)&header,sizeof(header));                //write dummy header as place holder, actual will be written later when total size of DDL is known
    cnt[i]=0; //reset counters
  }
  
  for(Int_t iChN=1;iChN<=kNchambers;iChN++){ //digits are stored on chamber by chamber basis   
    for(Int_t iDigN=0;iDigN<Digits(iChN)->GetEntriesFast();iDigN++){//digits loop for a given chamber
      AliRICHDigit *pDig=(AliRICHDigit*)Digits(iChN)->At(iDigN);
      Int_t ddl=pDig->Dig2Raw(w32);  //ddl is 0..13 
      file[ddl].write((char*)&w32,sizeof(w32));  cnt[ddl]++;//write formated digit to the propriate file (as decided in Dig2Raw) and increment corresponding counter
    }//digits loop for a given chamber
  }//chambers loop    
  for(Int_t i=0;i<AliRICHDigit::kNddls;i++){
    header.fSize=sizeof(header)+cnt[i]*sizeof(w32); //now calculate total number of bytes for each DDL file  
    header.SetAttribute(0); 
    file[i].seekp(0); file[i].write((char*)&header,sizeof(header));//rewrite DDL header with fSize field properly set
    file[i].close();                                               //close DDL file
  }
  GetLoader()->UnloadDigits();
  AliDebug(1,"Stop.");      
}//Digits2Raw()
//__________________________________________________________________________________________________
