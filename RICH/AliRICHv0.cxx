//**************************************************************************
//  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//                                                                         *
//  Author: The ALICE Off-line Project.                                    *
//  Contributors are mentioned in the code where appropriate.              *
//                                                                         *
//  Permission to use, copy, modify and distribute this software and its   *
//  documentation strictly for non-commercial purposes is hereby granted   *
//  without fee, provided that the above copyright notice appears in all   *
//  copies and that both the copyright notice and this permission notice   *
//  appear in the supporting documentation. The authors make no claims     *
//  about the suitability of this software for any purpose. It is          *
//  provided "as is" without express or implied warranty.                  *
//**************************************************************************

#include "AliRICHv0.h"
#include "AliRICHChamber.h" 
#include <AliRun.h>
#include <AliMC.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <TGeoManager.h>

ClassImp(AliRICHv0)

void AliRICHv0::StepManager()
{
//This StepManager is a provision for different test-learn activities on the current MC layer  
  const char *sParticle;
  switch(gMC->TrackPid()){
    case kProton:
      sParticle="proton";break;
    case kNeutron:
      sParticle="neutron";break;
    case kGamma:
      sParticle="gamma";break;
    case kCerenkov:
      sParticle="photon";break;
    case kPi0:
      sParticle="Pi0";break;  
    case kElectron:
      sParticle="electron";break;  
    default:
      sParticle="not known";break;
  }

  Info("","event=%i hunt=%i tid=%i pid=%i(%s) m=%f q=%3.1f dEdX=%9.3f",
                            gMC->CurrentEvent(),
                            fIshunt,
                            gAlice->GetMCApp()->GetCurrentTrackNumber(),
                            gMC->TrackPid(),
                            sParticle,
                            gMC->TrackMass(),
                            gMC->TrackCharge(),
                            gMC->Edep());
  Info("","Flags:alive(%i) disap(%i) enter(%i) exit(%i) inside(%i) out(%i) stop(%i) new(%i)",
                            gMC->IsTrackAlive(),
                            gMC->IsTrackDisappeared(),
                            gMC->IsTrackEntering(),
                            gMC->IsTrackExiting(),
                            gMC->IsTrackInside(),
                            gMC->IsTrackOut(),
                            gMC->IsTrackStop(),
                            gMC->IsNewTrack());
  Int_t copy0,copy1,copy2,copy3;
  Int_t vid0=gMC->CurrentVolID(copy0);
  Int_t vid1=gMC->CurrentVolOffID(1,copy1);
  Int_t vid2=gMC->CurrentVolOffID(2,copy2);
  Int_t vid3=gMC->CurrentVolOffID(3,copy3);
  Info("","vid0=%i(%s)c%i vid1=%i(%s)c%i vid2=%i(%s)c%i vid3=%i(%s)c%i   %s-%s-%s-%s",
                      vid0,gMC->VolName(vid0),copy0, 
                      vid1,gMC->VolName(vid1),copy1, 
                      vid2,gMC->VolName(vid2),copy2, 
                      vid3,gMC->VolName(vid3),copy3, 
                      gMC->CurrentVolName(),
                      gMC->CurrentVolOffName(1),
                      gMC->CurrentVolOffName(2),
                      gMC->CurrentVolOffName(3));
  
  Float_t a,z,den,rad,abs; a=z=den=rad=abs=kBad;
  Int_t mid=gMC->CurrentMaterial(a,z,den,rad,abs);
  Info("","mid=%i a=%7.2f z=%7.2f den=%7.2f rad=%7.2f abs=%7.2f",mid,a,z,den,rad,abs);
  
  TLorentzVector x4;
  gMC->TrackPosition(x4);
  Float_t glo[3],loc[3];
  glo[0]=x4.X();glo[1]=x4.Y();glo[2]=x4.Z();  
  gMC->Gmtod(glo,loc,1);
  Info("","glo(%+8.3f,%+8.3f,%+8.3f) r=%8.3f theta=%8.3f phi=%8.3f",
                      glo[0],glo[1],glo[2],x4.Rho(),x4.Theta()*TMath::RadToDeg(),x4.Phi()*TMath::RadToDeg());  
  Info("","loc(%+8.3f,%+8.3f,%8.3f) by gMC->Gmtod()",         loc[0],loc[1],loc[2]);  
  if(gMC->VolId("CSI ")==gMC->CurrentVolID(copy0)){
    Int_t iChamber;
    gMC->CurrentVolOffID(2,iChamber);
    TVector2 x2=C(iChamber)->Glob2Loc(x4);
    Info("","loc(%+8.3f,%+8.3f) by Glob2Loc",      x2.X(),x2.Y());  
  }
  Info("","end of current step\n");
}//StepManager()

void AliRICHv0::CreateGeometry()
{
  if(GetDebug())Info("CreateGeometry","Start v0.");  
  
  Double_t cm=1,mm=0.1;
  //place radioactive source compartment if needed
  Double_t containerLen=21.8*mm/2                        ,containerR=38*mm/2;  
  Double_t screwLen=15*mm/2                              ,screwR=5*mm/2;  
  Double_t srcLen=10*mm/2                                ,srcR=2*mm/2;  
  Double_t perpexLen=20*mm/2                             ,perpexR=34*mm/2;  
  Double_t perpexWholeLen=10*mm/2                        ,perpexWholeR=4*mm/2;  
  Double_t alBottomLen=containerLen-perpexLen            ,alWholeR=5*mm/2;      
//volumes    
  Double_t par[3];

//  Double_t anodWireD=20*mkm, cathWireD=100*mkm, collWireD=100*mkm;
        
  Double_t pcZ2=0.5*mm;
  par[0]=68.8*cm;            par[1]=70.86*cm;                par[2]=13*cm;           gMC->Gsvolu("RICH","BOX ",(*fIdtmed)[kCH4],par,3); //RICH    
  par[0]=P()->PcSizeX()/2;   par[1]=P()->PcSizeY()/2;        par[2]=pcZ2;            gMC->Gsvolu("CSI ","BOX ",(*fIdtmed)[kCSI],par,3); //CSI
  par[0]=P()->PcSizeX()/2;   par[1]=P()->PcSizeY()/2;        par[2]=P()->GapAmp();   gMC->Gsvolu("GAP ","BOX ",(*fIdtmed)[kCH4],par,3); //GAP
//  par[0]=0;                  par[1]=;                        par[2]=P()->PcSizeX()/2;gMC->Gsvolu("WIAN","TUBE",(*fIdtmed)[kW],par,3);//Anod wire
//  par[0]=0;                  par[1]=srcR;                    par[2]=P()->PcSizeX()/2;gMC->Gsvolu("WICA","TUBE",(*fIdtmed)[kCu],par,3);//Cathod wire
//  par[0]=0;                  par[1]=srcR;                    par[2]=P()->PcSizeX()/2;gMC->Gsvolu("WICO","TUBE",(*fIdtmed)[kCu],par,3);//Collect wire
    
  par[0]=0      ;par[1]=containerR  ;par[2]=containerLen;  gMC->Gsvolu("CONT","TUBE",(*fIdtmed)[kCH4]   ,par,3); //container 
  par[0]=perpexR;par[1]=containerR  ;par[2]=containerLen;  gMC->Gsvolu("ALWA","TUBE",(*fIdtmed)[kAl]    ,par,3); //Al cylindric wall
  par[0]=0      ;par[1]=perpexR     ;par[2]=alBottomLen;   gMC->Gsvolu("ALBO","TUBE",(*fIdtmed)[kAl]    ,par,3); //Al bottom 
  par[0]=0      ;par[1]=alWholeR    ;par[2]=alBottomLen;   gMC->Gsvolu("ALWH","TUBE",(*fIdtmed)[kCH4]   ,par,3); //Whole in Al bottom
  par[0]=0      ;par[1]=perpexR     ;par[2]=perpexLen;     gMC->Gsvolu("PEPL","TUBE",(*fIdtmed)[kPerpex],par,3); //Perpex plug
  par[0]=0      ;par[1]=perpexWholeR;par[2]=perpexWholeLen;gMC->Gsvolu("PEWH","TUBE",(*fIdtmed)[kCH4]   ,par,3); //Whole in Perpex
  par[0]=0      ;par[1]=screwR      ;par[2]=screwLen;      gMC->Gsvolu("SCRE","TUBE",(*fIdtmed)[kSteel] ,par,3); //Screw
  par[0]=0      ;par[1]=srcR        ;par[2]=srcLen;        gMC->Gsvolu("SOUR","TUBE",(*fIdtmed)[kSteel] ,par,3); //Source itself
//nodes    
  gMC->Gspos("RICH",1,"ALIC",0,0,-9                                ,0,"ONLY");     //RICH in ALIC
  gMC->Gspos("CSI ",1,"RICH",0,0,-(P()->GapProx()+pcZ2)           ,0,"ONLY");     //CsI  in ALIC
  gMC->Gspos("GAP ",1,"RICH",0,0,-(P()->GapProx()-P()->GapAmp()/2),0,"ONLY");     //GAP in ALIC
    
  gMC->Gspos("CONT",1,"RICH",0,0,1*cm                             ,0,"ONLY");     //Sr90 container in RICH
  gMC->Gspos("ALWA",1,"CONT",0,0,0                                ,0,"ONLY");     //Al wall
  gMC->Gspos("ALBO",1,"CONT",0,0,-containerLen+alBottomLen        ,0,"ONLY");     //Al bottom
  gMC->Gspos("PEPL",1,"CONT",0,0,containerLen-perpexLen           ,0,"ONLY");     //Perpex plug
   
  gMC->Gspos("ALWH",1,"ALBO",6*mm,0,0                             ,0,"ONLY");     //Whole in Al bottom
    
  gMC->Gspos("PEWH",1,"PEPL",6*mm,0,-perpexLen+perpexWholeLen     ,0,"ONLY");     //Whole in Perpex plug
  gMC->Gspos("SOUR",1,"PEPL",6*mm,0, perpexLen-srcLen             ,0,"ONLY");     //Source in Perpex plug
  gMC->Gspos("SCRE",1,"PEPL",0,0,perpexLen-screwLen               ,0,"ONLY");     //Screw in Perpex plug
  if(GetDebug())Info("CreateGeometry","Stop v0.");  
}//CreateGeometry()
