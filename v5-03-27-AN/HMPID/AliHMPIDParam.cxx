//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************
#include "AliHMPIDParam.h"  //class header
#include "AliHMPIDDigit.h"  //ctor
#include "AliLog.h"         //general
#include <AliRunLoader.h>   //Stack()
#include <AliStack.h>       //Stack()
#include "AliCDBManager.h"  //ctor
#include "AliCDBEntry.h"    //ctor
#include <TLatex.h>         //TestTrans()  
#include <TView.h>          //TestTrans()
#include <TPolyMarker3D.h>  //TestTrans()
#include <TRotation.h>
#include <TParticle.h>      //Stack()    
#include <TGeoPhysicalNode.h> //ctor
#include <TGeoBBox.h>
#include <TF1.h>                 //ctor

ClassImp(AliHMPIDParam)


// Mathieson constant definition
const Double_t AliHMPIDParam::fgkD     = 0.222500;  // ANODE-CATHODE distance 0.445/2
//                                                                                          K3 = 0.66 along the wires (anode-cathode/wire pitch=0.5625)
const Double_t AliHMPIDParam::fgkSqrtK3x = TMath::Sqrt(0.66);
const Double_t AliHMPIDParam::fgkK2x     = TMath::PiOver2()*(1 - 0.5*fgkSqrtK3x);
const Double_t AliHMPIDParam::fgkK1x     = 0.25*fgkK2x*fgkSqrtK3x/TMath::ATan(fgkSqrtK3x);
const Double_t AliHMPIDParam::fgkK4x     = fgkK1x/(fgkK2x*fgkSqrtK3x);
//                                                                                          K3 = 0.87 along the wires (anode-cathode/wire pitch=0.5625)
const Double_t AliHMPIDParam::fgkSqrtK3y = TMath::Sqrt(0.87);
const Double_t AliHMPIDParam::fgkK2y     = TMath::PiOver2()*(1 - 0.5*fgkSqrtK3y);
const Double_t AliHMPIDParam::fgkK1y     = 0.25*fgkK2y*fgkSqrtK3y/TMath::ATan(fgkSqrtK3y);
const Double_t AliHMPIDParam::fgkK4y     = fgkK1y/(fgkK2y*fgkSqrtK3y);
//
  

Float_t AliHMPIDParam::fgkMinPcX[]={0.,0.,0.,0.,0.,0.};
Float_t AliHMPIDParam::fgkMaxPcX[]={0.,0.,0.,0.,0.,0.};
Float_t AliHMPIDParam::fgkMinPcY[]={0.,0.,0.,0.,0.,0.};
Float_t AliHMPIDParam::fgkMaxPcY[]={0.,0.,0.,0.,0.,0.};

Bool_t AliHMPIDParam::fgMapPad[160][144][7];

Float_t AliHMPIDParam::fgCellX=0.;
Float_t AliHMPIDParam::fgCellY=0.;

Float_t AliHMPIDParam::fgPcX=0;
Float_t AliHMPIDParam::fgPcY=0;

Float_t AliHMPIDParam::fgAllX=0;
Float_t AliHMPIDParam::fgAllY=0;

Bool_t AliHMPIDParam::fgInstanceType=kTRUE;  

AliHMPIDParam* AliHMPIDParam::fgInstance=0x0;        //singleton pointer               

Int_t AliHMPIDParam::fgNSigmas  = 4;
Int_t AliHMPIDParam::fgThreshold= 4;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDParam::AliHMPIDParam(Bool_t noGeo):
  TNamed("HmpidParam","default version"),
  fX(0), fY(0), fRefIdx(1.28947),fPhotEMean(6.675),fTemp(25)                          //just set a refractive index for C6F14 at ephot=6.675 eV @ T=25 C
{
// Here all the intitializition is taken place when AliHMPIDParam::Instance() is invoked for the first time.
// In particular, matrices to be used for LORS<->MARS trasnformations are initialized from TGeo structure.    
// Note that TGeoManager should be already initialized from geometry.root file  

  AliCDBManager *pCDB = AliCDBManager::Instance();
  if(!pCDB) {
     AliWarning("No Nmean C6F14 from OCDB. Default is taken from ctor.");
  } else {
    AliCDBEntry *pNmeanEnt =pCDB->Get("HMPID/Calib/Nmean"); //contains TObjArray of 42 TF1 + 1 EPhotMean
    if(!pNmeanEnt) {
      AliWarning("No Nmean C6F14 from OCDB. Default is taken from ctor.");
    } else {
      TObjArray *pNmean = (TObjArray*)pNmeanEnt->GetObject();
      if(pNmean->GetEntries()==43) {                                               //for backward compatibility
        Double_t tmin,tmax;
        ((TF1*)pNmean->At(42))->GetRange(tmin,tmax);
        fPhotEMean = ((TF1*)pNmean->At(42))->Eval(tmin);                          //photon eMean from OCDB
        AliInfo(Form("EPhotMean = %f eV successfully loaded from OCDB",fPhotEMean));
      } else {
        AliWarning("For backward compatibility EPhotMean is taken from ctor.");
      }
    }
  }

  fRefIdx = MeanIdxRad(); //initialization of the running ref. index of freon
  
  Float_t dead=2.6;// cm of the dead zones between PCs-> See 2CRC2099P1


  if(noGeo==kTRUE) fgInstanceType=kFALSE;                                                   //instance from ideal geometry, no actual geom is present
    
  if(noGeo==kFALSE && !gGeoManager)  
  {
    TGeoManager::Import("geometry.root");
    if(!gGeoManager) AliFatal("!!!!!!No geometry loaded!!!!!!!");
  }
  
  fgCellX=0.8;fgCellY=0.84;
  
  if(!noGeo==kTRUE){
    TGeoVolume *pCellVol = gGeoManager->GetVolume("Hcel");
    if(pCellVol) {
      TGeoBBox *bcell = (TGeoBBox *)pCellVol->GetShape();
      fgCellX=2.*bcell->GetDX(); fgCellY = 2.*bcell->GetDY();  // overwrite the values with the read ones
    }
  }    
  fgPcX=80.*fgCellX; fgPcY = 48.*fgCellY;
  fgAllX=2.*fgPcX+dead;
  fgAllY=3.*fgPcY+2.*dead;

  fgkMinPcX[1]=fgPcX+dead; fgkMinPcX[3]=fgkMinPcX[1];  fgkMinPcX[5]=fgkMinPcX[3];
  fgkMaxPcX[0]=fgPcX; fgkMaxPcX[2]=fgkMaxPcX[0];  fgkMaxPcX[4]=fgkMaxPcX[2];
  fgkMaxPcX[1]=fgAllX; fgkMaxPcX[3]=fgkMaxPcX[1];  fgkMaxPcX[5]=fgkMaxPcX[3];

  fgkMinPcY[2]=fgPcY+dead; fgkMinPcY[3]=fgkMinPcY[2];  
  fgkMinPcY[4]=2.*fgPcY+2.*dead; fgkMinPcY[5]=fgkMinPcY[4];
  fgkMaxPcY[0]=fgPcY; fgkMaxPcY[1]=fgkMaxPcY[0];  
  fgkMaxPcY[2]=2.*fgPcY+dead; fgkMaxPcY[3]=fgkMaxPcY[2]; 
  fgkMaxPcY[4]=fgAllY; fgkMaxPcY[5]=fgkMaxPcY[4];   
    
  fX=0.5*SizeAllX();
  fY=0.5*SizeAllY();
  
      
  for(Int_t ich=kMinCh;ich<=kMaxCh;ich++) {
    for(Int_t padx=0;padx<160;padx++) {
       for(Int_t pady=0;pady<144;pady++) {
         fgMapPad[padx][pady][ich] = kTRUE;             //init all the pads are active at the beginning....
       }
     }
   }
     

  for(Int_t i=kMinCh;i<=kMaxCh;i++)
    if(gGeoManager && gGeoManager->IsClosed()) {
      TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(Form("/HMPID/Chamber%i",i));
      if (!pne) {
        AliErrorClass(Form("The symbolic volume %s does not correspond to any physical entry!",Form("HMPID_%i",i)));
        fM[i]=new TGeoHMatrix;
        IdealPosition(i,fM[i]);
      } else {
        TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
        if(pnode) fM[i]=new TGeoHMatrix(*(pnode->GetMatrix()));
        else {
          fM[i]=new TGeoHMatrix;
          IdealPosition(i,fM[i]);
        }
      }
    } else{
      fM[i]=new TGeoHMatrix;
      IdealPosition(i,fM[i]);
    } 
  fgInstance=this; 
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDParam::Print(Option_t* opt) const
{
// print some usefull (hopefully) info on some internal guts of HMPID parametrisation 
  
  for(Int_t i=0;i<7;i++) fM[i]->Print(opt);
}//Print()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDParam::IdealPosition(Int_t iCh, TGeoHMatrix *pMatrix)
{
// Construct ideal position matrix for a given chamber
// Arguments: iCh- chamber ID; pMatrix- pointer to precreated unity matrix where to store the results
//   Returns: none
  const Double_t kAngHor=19.5;        //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;          //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;          //  common HMPID rotation with respect to x axis  30   grad     
  const Double_t kTrans[3]={490,0,0}; //  center of the chamber is on window-gap surface
  pMatrix->RotateY(90);               //  rotate around y since initial position is in XY plane -> now in YZ plane
  pMatrix->SetTranslation(kTrans);    //  now plane in YZ is shifted along x 
  switch(iCh){
    case 0:                pMatrix->RotateY(kAngHor);  pMatrix->RotateZ(-kAngVer);  break; //right and down 
    case 1:                                            pMatrix->RotateZ(-kAngVer);  break; //down              
    case 2:                pMatrix->RotateY(kAngHor);                               break; //right 
    case 3:                                                                         break; //no rotation
    case 4:                pMatrix->RotateY(-kAngHor);                              break; //left   
    case 5:                                            pMatrix->RotateZ(kAngVer);   break; //up
    case 6:                pMatrix->RotateY(-kAngHor); pMatrix->RotateZ(kAngVer);   break; //left and up 
  }
  pMatrix->RotateZ(kAngCom);     //apply common rotation  in XY plane    
   
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDParam::Stack(Int_t evt,Int_t tid)
{
// Prints some useful info from stack
// Arguments: evt - event number. if not -1 print info only for that event
//            tid - track id. if not -1 then print it and all it's mothers if any   
//   Returns: mother tid of the given tid if any
  AliRunLoader *pAL=AliRunLoader::Open(); 
  if(pAL->LoadHeader()) return -1;
  if(pAL->LoadKinematics()) return -1;
  
  Int_t mtid=-1;
  Int_t iNevt=pAL->GetNumberOfEvents();
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(evt!=-1 && evt!=iEvt) continue; //in case one needs to print the requested event, ignore all others
    pAL->GetEvent(iEvt);    
    AliStack *pStack=pAL->Stack();  
    if(tid==-1){                        //print all tids for this event
      for(Int_t i=0;i<pStack->GetNtrack();i++) pStack->Particle(i)->Print();
          Printf("totally %i tracks including %i primaries for event %i out of %i event(s)",
          pStack->GetNtrack(),pStack->GetNprimary(),iEvt,iNevt);
    }else{                              //print only this tid and it;s mothers
      if(tid<0 || tid>pStack->GetNtrack()) {Printf("Wrong tid, valid tid range for event %i is 0-%i",iEvt,pStack->GetNtrack());break;}
      TParticle *pTrack=pStack->Particle(tid); mtid=pTrack->GetFirstMother();
      TString str=pTrack->GetName();
      while((tid=pTrack->GetFirstMother()) >= 0){
        pTrack=pStack->Particle(tid);
        str+=" from ";str+=pTrack->GetName();
      } 
    }//if(tid==-1)      
  }//events loop
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return mtid;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDParam::StackCount(Int_t pid,Int_t evt)
{
// Counts total number of particles of given sort (including secondary) for a given event
  AliRunLoader *pAL=AliRunLoader::Open(); 
  pAL->GetEvent(evt);    
  if(pAL->LoadHeader()) return 0;
  if(pAL->LoadKinematics()) return 0;
  AliStack *pStack=pAL->Stack();
  
  Int_t iCnt=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) if(pStack->Particle(i)->GetPdgCode()==pid) iCnt++;
  
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return iCnt;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDParam::Sigma2(Double_t trkTheta,Double_t trkPhi,Double_t ckovTh, Double_t ckovPh)
{
// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  TVector3 v(-999,-999,-999);
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*GetRefIdx());
  
  if(trkBeta > 1) trkBeta = 1;                 //protection against bad measured thetaCer  
  if(trkBeta < 0) trkBeta = 0.0001;            //

  v.SetX(SigLoc (trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));
  v.SetY(SigGeom(trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));
  v.SetZ(SigCrom(trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));

  return v.Mag2();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDParam::SigLoc(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosf     = TMath::Cos(trkPhi);
  Double_t sinfd    = TMath::Sin(phiDelta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);        // formula (after 8 in the text)
  if (k<0) return 1e10;
  Double_t mu =sint*sinf+tantheta*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
  Double_t e  =sint*cosf+tantheta*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

  Double_t kk = betaM*TMath::Sqrt(k)/(GapThick()*alpha);                            // formula (6) and (7)
  Double_t dtdxc = kk*(k*(cosfd*cosf-cost*sinfd*sinf)-(alpha*mu/(betaM*betaM))*sint*sinfd); // formula (6)           
  Double_t dtdyc = kk*(k*(cosfd*sinf+cost*sinfd*cosf)+(alpha* e/(betaM*betaM))*sint*sinfd); // formula (7)            pag.4

  Double_t errX = 0.2,errY=0.25;                                                            //end of page 7
  return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDParam::SigCrom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t dtdn = cost*GetRefIdx()*betaM*betaM/(alpha*tantheta);                    // formula (12)
            
//  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
  Double_t f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

  return f*dtdn;
}//SigCrom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDParam::SigGeom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of geometric error (due to lack of knowledge of creation point in radiator) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Formulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    

  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t costheta = TMath::Cos(thetaC);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                  // formula (11)
  
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);         // formula (after 8 in the text)
  if (k<0) return 1e10;

  Double_t eTr = 0.5*RadThick()*betaM*TMath::Sqrt(k)/(GapThick()*alpha);     // formula (14)
  Double_t lambda = (1.-sint*sinf)*(1.+sint*sinf);                                                  // formula (15)

  Double_t c1 = 1./(1.+ eTr*k/(alpha*alpha*costheta*costheta));                              // formula (13.a)
  Double_t c2 = betaM*TMath::Power(k,1.5)*tantheta*lambda/(GapThick()*alpha*alpha);  // formula (13.b)
  Double_t c3 = (1.+eTr*k*betaM*betaM)/((1+eTr)*alpha*alpha);                                // formula (13.c)
  Double_t c4 = TMath::Sqrt(k)*tantheta*(1-lambda)/(GapThick()*betaM);               // formula (13.d)
  Double_t dtdT = c1 * (c2+c3*c4);
  Double_t trErr = RadThick()/(TMath::Sqrt(12.)*cost);

  return trErr*dtdT;
}//SigGeom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
