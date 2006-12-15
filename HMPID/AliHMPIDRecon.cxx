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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDRecon                                                         //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom    //
// for single chamber                                                   //
//////////////////////////////////////////////////////////////////////////

#include "AliHMPIDRecon.h"  //class header
#include "AliHMPIDCluster.h" //CkovAngle()
#include <TRotation.h>     //TracePhoton()
#include <TH1D.h>          //HoughResponse()
#include <TClonesArray.h>  //CkovAngle()

#include <TTree.h>         //Display()
#include <TFile.h>         //Display()
#include <AliESD.h>        //Display()
#include <TPolyMarker.h>   //Display()
#include <TLatex.h>        //Display()
#include <TCanvas.h>       //Display()


const Double_t AliHMPIDRecon::fgkRadThick=1.5;
const Double_t AliHMPIDRecon::fgkWinThick=0.5;
const Double_t AliHMPIDRecon::fgkGapThick=8.0;
const Double_t AliHMPIDRecon::fgkRadIdx  =1.292;
const Double_t AliHMPIDRecon::fgkWinIdx  =1.5787;
const Double_t AliHMPIDRecon::fgkGapIdx  =1.0005;


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecon::AliHMPIDRecon():TTask("RichRec","RichPat"),
  fPhotCnt(-1),
  fCkovSigma2(0),
  fIsWEIGHT(kFALSE),
  fDTheta(0.001),
  fWindowWidth(0.045),
  fTrkDir(TVector3(0,0,1)),fTrkPos(TVector2(30,40))  
{
// main ctor
  for (Int_t i=0; i<3000; i++) {
    fPhotFlag[i] =  0;
    fPhotCkov[i] = -1;
    fPhotPhi [i] = -1;
    fPhotWei [i] =  0;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::CkovAngle(TClonesArray *pCluLst,Int_t &iNaccepted)
{
// Pattern recognition method based on Hough transform
// Arguments: pCluLst  - list of clusters for this chamber   
//   Returns:          - track ckov angle, [rad], 
  
  if(pCluLst->GetEntries()>200) fIsWEIGHT = kTRUE; // offset to take into account bkg in reconstruction
  else                          fIsWEIGHT = kFALSE;

  // Photon Flag:  Flag = 0 initial set; Flag = 1 good candidate (charge compatible with photon); Flag = 2 photon used for the ring;

  fPhotCnt=0;                                                      
  for (Int_t iClu=0; iClu<pCluLst->GetEntriesFast();iClu++){//clusters loop
    AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(iClu);                       //get pointer to current cluster    
    if(pClu->Q()>100) continue;                                                             //avoid MIP clusters from bkg
    
    fPhotCkov[fPhotCnt]=FindPhotCkov(pClu->X(),pClu->Y());                                  //find ckov angle for this  photon candidate
    fPhotCnt++;         //increment counter of photon candidates
  }//clusters loop

  iNaccepted=FlagPhot(HoughResponse()); //flag photons according to individual theta ckov with respect to most probable track theta ckov
  if(iNaccepted<1) return -11; 
  else             return FindRingCkov(pCluLst->GetEntries());  //find best Theta ckov for ring i.e. track
}//ThetaCerenkov()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FindPhotCkov(Double_t cluX,Double_t cluY)
{
// Finds Cerenkov angle  for this photon candidate
// Arguments: cluX,cluY - position of cadidate's cluster  
//   Returns: Cerenkov angle 

  TVector2 pos(cluX,cluY); Double_t cluR=(pos-fTrkPos).Mod();  Double_t phi=FindPhotPhi(cluX,cluY);      
  Double_t ckov1=0,ckov2=0.75;
  const Double_t kTol=0.05; 
  Int_t iIterCnt = 0;
  while(1){
    if(iIterCnt>=50) return -1;
    Double_t ckov=0.5*(ckov1+ckov2);
    Double_t dist=cluR-TracePhot(ckov,phi,pos); iIterCnt++;   //get distance between trial point and cluster position
    if     (dist> kTol) ckov1=ckov;                           //cluster @ larger ckov 
    else if(dist<-kTol) ckov2=ckov;                           //cluster @ smaller ckov
    else                return ckov;                          //precision achived         
  }
}//FindPhotTheta()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FindPhotPhi(Double_t cluX,Double_t cluY)
{
// Finds phi angle og photon candidate by considering the cluster's position  of this candudate w.r.t track position
  
  Double_t emiss=0; 
  return fPhotPhi[fPhotCnt]=TMath::ATan2(cluY-fTrkPos.Y()-emiss*TMath::Tan(fTrkDir.Theta())*TMath::Sin(fTrkDir.Phi()),
                                         cluX-fTrkPos.X()-emiss*TMath::Tan(fTrkDir.Theta())*TMath::Cos(fTrkDir.Phi()));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FindRingArea(Double_t ckovAng)const
{
// Find area inside the cerenkov ring which lays inside PCs
// Arguments: ckovThe - cernkov    
//   Returns: area of the ring in cm^2 for given theta ckov
   
  
  TVector2 pos1,pos2;
  
  const Int_t kN=100;
  Double_t area=0;
  for(Int_t i=0;i<kN;i++){
    TracePhot(ckovAng,Double_t(TMath::TwoPi()*i    /kN),pos1);//trace this photon 
    TracePhot(ckovAng,Double_t(TMath::TwoPi()*(i+1)/kN),pos2);//trace this photon 
    area+=(pos1-fTrkPos)*(pos2-fTrkPos);
      
  }
  return area;
}//FindRingArea()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FindRingCkov(Int_t)
{
// Loops on all Ckov candidates and estimates the best Theta Ckov for a ring formed by those candidates. Also estimates an error for that Theat Ckov
// collecting errors for all single Ckov candidates thetas. (Assuming they are independent)  
// Arguments: iNclus- total number of clusters in chamber for background estimation
//    Return: best estimation of track Theta ckov

  Double_t wei = 0.;
  Double_t weightThetaCerenkov = 0.;

  Double_t ckovMin=9999.,ckovMax=0.;
  Double_t sigma2 = 0;   //to collect error squared for this ring
  
  for(Int_t i=0;i<fPhotCnt;i++){//candidates loop
    if(fPhotFlag[i] == 2){
      if(fPhotCkov[i]<=0) continue;//?????????????????Flag photos = 2 may imply CkovEta = 0?????????????? 
      if(fPhotCkov[i]<ckovMin) ckovMin=fPhotCkov[i];  //find max and min Theta ckov from all candidates within probable window
      if(fPhotCkov[i]>ckovMax) ckovMax=fPhotCkov[i]; 
      weightThetaCerenkov += fPhotCkov[i]*fPhotWei[i];   wei += fPhotWei[i];                 //collect weight as sum of all candidate weghts   
      
      sigma2 += 1./Sigma2(fPhotCkov[i],fPhotPhi[i]);
    }
  }//candidates loop
  
  if(sigma2>0) fCkovSigma2=1./sigma2;
  else         fCkovSigma2=1e10;  
  
  if(wei != 0.) weightThetaCerenkov /= wei; else weightThetaCerenkov = 0.;  
  return weightThetaCerenkov;
}//FindCkovRing()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRecon::FlagPhot(Double_t ckov)
{
// Flag photon candidates if their individual ckov angle is inside the window around ckov angle returned by  HoughResponse()
// Arguments: ckov- value of most probable ckov angle for track as returned by HoughResponse()
//   Returns: number of photon candidates happened to be inside the window

  
  Int_t steps = (Int_t)((ckov )/ fDTheta); //how many times we need to have fDTheta to fill the distance between 0  and thetaCkovHough

  Double_t tmin = (Double_t)(steps - 1)*fDTheta;
  Double_t tmax = (Double_t)(steps)*fDTheta;
  Double_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;  tmax = tavg + 0.5*fWindowWidth;

  Int_t iInsideCnt = 0; //count photons which Theta ckov inside the window
  for(Int_t i=0;i<fPhotCnt;i++){//photon candidates loop
    if(fPhotCkov[i] >= tmin && fPhotCkov[i] <= tmax)	{ 
      fPhotFlag[i]=2;	  
      iInsideCnt++;
    }
  }
  return iInsideCnt;
}//FlagPhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::TracePhot(Double_t ckovThe,Double_t ckovPhi,TVector2 &pos)const
{
// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles, [rad]  (warning: not photon theta and phi)     
//   Returns: distance between photon point on PC and track projection  
  TRotation mtheta;   mtheta.RotateY(fTrkDir.Theta());
  TRotation mphi;       mphi.RotateZ(fTrkDir.Phi());  
  TRotation mrot=mphi*mtheta;
  
  TVector3  posCkov(fTrkPos.X(),fTrkPos.Y(),-0.5*fgkRadThick-fgkWinThick-fgkGapThick);   //RAD: photon position is track position @ middle of RAD 
  TVector3  dirCkov;   dirCkov.SetMagThetaPhi(1,ckovThe,ckovPhi);                        //initially photon is directed according to requested ckov angle
                                               dirCkov=mrot*dirCkov;                     //now we know photon direction in LORS
                       dirCkov.SetPhi(ckovPhi);   
  if(dirCkov.Theta() > TMath::ASin(1./fgkRadIdx)) return -999;//total refraction on WIN-GAP boundary
  
  Propagate(dirCkov,posCkov,-fgkWinThick-fgkGapThick); //go to RAD-WIN boundary  remeber that z=0 is PC plane
  Refract  (dirCkov,         fgkRadIdx,fgkWinIdx    ); //RAD-WIN refraction
  Propagate(dirCkov,posCkov,-fgkGapThick           );  //go to WIN-GAP boundary
  Refract  (dirCkov,         fgkWinIdx,fgkGapIdx    ); //WIN-GAP refraction
  Propagate(dirCkov,posCkov,0                     );   //go to PC
  
  pos.Set(posCkov.X(),posCkov.Y());
  return (pos-fTrkPos).Mod();
}//TracePhoton()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::Propagate(const TVector3 &dir,TVector3 &pos,Double_t z)const
{
// Finds an intersection point between a line and XY plane shifted along Z.
// Arguments:  dir,pos   - vector along the line and any point of the line
//             z         - z coordinate of plain 
//   Returns:  none
//   On exit:  pos is the position if this intesection if any
  static TVector3 nrm(0,0,1); 
         TVector3 pnt(0,0,z);
  
  TVector3 diff=pnt-pos;
  Double_t sint=(nrm*diff)/(nrm*dir);
  pos+=sint*dir;
}//Propagate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::Refract(TVector3 &dir,Double_t n1,Double_t n2)const
{
// Refract direction vector according to Snell law
// Arguments: 
//            n1 - ref idx of first substance
//            n2 - ref idx of second substance
//   Returns: none
//   On exit: dir is new direction
  Double_t sinref=(n1/n2)*TMath::Sin(dir.Theta());
  if(sinref>1.)    dir.SetXYZ(-999,-999,-999);
  else             dir.SetTheta(TMath::ASin(sinref));
}//Refract()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::HoughResponse()
{
//
//
//       
  Double_t kThetaMax=0.75;
  Int_t nChannels = (Int_t)(kThetaMax/fDTheta+0.5);
  TH1D *phots   = new TH1D("Rphot"  ,"phots"         ,nChannels,0,kThetaMax);
  TH1D *photsw  = new TH1D("RphotWeighted" ,"photsw" ,nChannels,0,kThetaMax);
  TH1D *resultw = new TH1D("resultw","resultw"       ,nChannels,0,kThetaMax);
  Int_t nBin = (Int_t)(kThetaMax/fDTheta);
  Int_t nCorrBand = (Int_t)(fWindowWidth/(2*fDTheta));
  
  for (Int_t i=0; i< fPhotCnt; i++){//photon cadidates loop
    Double_t angle = fPhotCkov[i];  if(angle<0||angle>kThetaMax) continue;
    phots->Fill(angle);
    Int_t bin = (Int_t)(0.5+angle/(fDTheta));
    Double_t weight=1.;
    if(fIsWEIGHT){
      Double_t lowerlimit = ((Double_t)bin)*fDTheta - 0.5*fDTheta;  Double_t upperlimit = ((Double_t)bin)*fDTheta + 0.5*fDTheta;   
      Double_t diffArea = FindRingArea(upperlimit)-FindRingArea(lowerlimit);
      if(diffArea>0) weight = 1./diffArea;
    }
    photsw->Fill(angle,weight);
    fPhotWei[i]=weight;
  }//photon candidates loop 
   
  for (Int_t i=1; i<=nBin;i++){
    Int_t bin1= i-nCorrBand;
    Int_t bin2= i+nCorrBand;
    if(bin1<1) bin1=1;
    if(bin2>nBin)bin2=nBin;
    Double_t sumPhots=phots->Integral(bin1,bin2);
    if(sumPhots<3) continue;                            // if less then 3 photons don't trust to this ring
    Double_t sumPhotsw=photsw->Integral(bin1,bin2);
    resultw->Fill((Double_t)((i+0.5)*fDTheta),sumPhotsw);
  } 
// evaluate the "BEST" theta ckov as the maximum value of histogramm
  Double_t *pVec = resultw->GetArray();
  Int_t locMax = TMath::LocMax(nBin,pVec);
  phots->Delete();photsw->Delete();resultw->Delete(); // Reset and delete objects
  
  return (Double_t)(locMax*fDTheta+0.5*fDTheta); //final most probable track theta ckov   
}//HoughResponse()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::Sigma2(Double_t ckovTh, Double_t ckovPh)const
{
// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  TVector3 v(-999,-999,-999);
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*fgkRadIdx);

  v.SetX(SigLoc (ckovTh,ckovPh,trkBeta));
  v.SetY(SigGeom(ckovTh,ckovPh,trkBeta));
  v.SetZ(SigCrom(ckovTh,ckovPh,trkBeta));

  return v.Mag2();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::SigLoc(Double_t thetaC, Double_t phiC,Double_t betaM)const
{
// Analithical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  Double_t phiDelta = phiC - fTrkDir.Phi();

  Double_t alpha =TMath::Cos(fTrkDir.Theta())-TMath::Tan(thetaC)*TMath::Cos(phiDelta)*TMath::Sin(fTrkDir.Theta());
  Double_t k = 1.-fgkRadIdx*fgkRadIdx+alpha*alpha/(betaM*betaM);
  if (k<0) return 1e10;

  Double_t mu =TMath::Sin(fTrkDir.Theta())*TMath::Sin(fTrkDir.Phi())+TMath::Tan(thetaC)*(TMath::Cos(fTrkDir.Theta())*TMath::Cos(phiDelta)*TMath::Sin(fTrkDir.Phi())+TMath::Sin(phiDelta)*TMath::Cos(fTrkDir.Phi()));
  Double_t e  =TMath::Sin(fTrkDir.Theta())*TMath::Cos(fTrkDir.Phi())+TMath::Tan(thetaC)*(TMath::Cos(fTrkDir.Theta())*TMath::Cos(phiDelta)*TMath::Cos(fTrkDir.Phi())-TMath::Sin(phiDelta)*TMath::Sin(fTrkDir.Phi()));

  Double_t kk = betaM*TMath::Sqrt(k)/(8*alpha);
  Double_t dtdxc = kk*(k*(TMath::Cos(phiDelta)*TMath::Cos(fTrkDir.Phi())-TMath::Cos(fTrkDir.Theta())*TMath::Sin(phiDelta)*TMath::Sin(fTrkDir.Phi()))-(alpha*mu/(betaM*betaM))*TMath::Sin(fTrkDir.Theta())*TMath::Sin(phiDelta));
  Double_t dtdyc = kk*(k*(TMath::Cos(phiDelta)*TMath::Sin(fTrkDir.Phi())+TMath::Cos(fTrkDir.Theta())*TMath::Sin(phiDelta)*TMath::Cos(fTrkDir.Phi()))+(alpha* e/(betaM*betaM))*TMath::Sin(fTrkDir.Theta())*TMath::Sin(phiDelta));

  return  TMath::Sqrt(0.2*0.2*dtdxc*dtdxc + 0.25*0.25*dtdyc*dtdyc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::SigCrom(Double_t thetaC, Double_t phiC,Double_t betaM)const
{
// Analithical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  Double_t phiDelta = phiC - fTrkDir.Phi();
  Double_t alpha =TMath::Cos(fTrkDir.Theta())-TMath::Tan(thetaC)*TMath::Cos(phiDelta)*TMath::Sin(fTrkDir.Theta());

  Double_t dtdn = TMath::Cos(fTrkDir.Theta())*fgkRadIdx*betaM*betaM/(alpha*TMath::Tan(thetaC));
            
  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);

  return f*dtdn;
}//SigCrom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::SigGeom(Double_t thetaC, Double_t phiC,Double_t betaM)const
{
// Analithical calculation of geometric error (due to lack of knowledge of creation point in radiator) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Formulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    

  Double_t phiDelta = phiC - fTrkDir.Phi();
  Double_t alpha =TMath::Cos(fTrkDir.Theta())-TMath::Tan(thetaC)*TMath::Cos(phiDelta)*TMath::Sin(fTrkDir.Theta());

  Double_t k = 1.-fgkRadIdx*fgkRadIdx+alpha*alpha/(betaM*betaM);
  if (k<0) return 1e10;

  Double_t eTr = 0.5*1.5*betaM*TMath::Sqrt(k)/(8*alpha);
  Double_t lambda = 1.-TMath::Sin(fTrkDir.Theta())*TMath::Sin(fTrkDir.Theta())*TMath::Sin(phiC)*TMath::Sin(phiC);

  Double_t c = 1./(1.+ eTr*k/(alpha*alpha*TMath::Cos(thetaC)*TMath::Cos(thetaC)));
  Double_t i = betaM*TMath::Tan(thetaC)*lambda*TMath::Power(k,1.5);
  Double_t ii = 1.+eTr*betaM*i;

  Double_t err = c * (i/(alpha*alpha*8) +  ii*(1.-lambda) / ( alpha*alpha*8*betaM*(1.+eTr)) );
  Double_t trErr = 1.5/(TMath::Sqrt(12.)*TMath::Cos(fTrkDir.Theta()));

  return trErr*err;
}//SigGeom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::Display() 
{
// Display digits, reconstructed tracks intersections and HMPID rings if available 
// Arguments: none
//   Returns: none    
  TFile *pEsdFl=TFile::Open("AliESDs.root");     if(!pEsdFl || !pEsdFl->IsOpen()) return;//open AliESDs.root                                                                    
  TTree *pEsdTr=(TTree*) pEsdFl->Get("esdTree"); if(!pEsdTr)                      return;//get ESD tree
                                                                 
  AliESD *pEsd=new AliESD;  pEsdTr->SetBranchAddress("ESD", &pEsd);
  
  TPolyMarker  *pDigMap[7]; //digits map
  TPolyMarker  *pTrkMap[7]; Int_t aTrkCnt[7]; //TRKxPC intersection map
  
  for(Int_t i=0;i<7;i++){
                  pDigMap[i]=new TPolyMarker(); pDigMap[i]->SetMarkerStyle(25); pDigMap[i]->SetMarkerSize(0.5); pDigMap[i]->SetMarkerColor(kGreen); 
    aTrkCnt[i]=0; pTrkMap[i]=new TPolyMarker(); pTrkMap[i]->SetMarkerStyle(4);  pTrkMap[i]->SetMarkerSize(0.5); pTrkMap[i]->SetMarkerColor(kRed); 
  }

  AliHMPIDRecon rec;
  
  TLatex t;
  TCanvas *pC = new TCanvas("HMPIDDisplay","HMPID Display",0,0,1226,900);  pC->Divide(3,3);
  
  for(Int_t iEvt=0;iEvt<pEsdTr->GetEntries();iEvt++) {                //events loop
    pC->cd(3);  t.DrawText(0.2,0.4,Form("Event %i",iEvt));        //print current event number
    pEsdTr->GetEntry(iEvt);                                       //get ESD for this event   
    for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//ESD tracks loop
      AliESDtrack *pTrk = pEsd->GetTrack(iTrk);             //
      Float_t th,ph,x,y; pTrk->GetHMPIDtrk(x,y,th,ph); if(x<0) continue;
      Int_t ch=pTrk->GetHMPIDcluIdx()/1000000;
      pTrkMap[ch]->SetPoint(aTrkCnt[ch]++,x,y);
    }//ESD tracks loop
    
//     al->GetEvent(iEvt);   rl->TreeD()->GetEntry(0); //get digits list
    for(Int_t iCh=0;iCh<7;iCh++) {//chambers loop
//       for(Int_t iDig=0;iDig < r->DigLst(iCh)->GetEntries();iDig++) {      //digits loop
//         AliHMPIDDigit *pDig = (AliHMPIDDigit*)r->DigLst(iCh)->At(iDig);     
//         pDigMap[iCh]->SetPoint(iDig,pDig->LorsX(),pDig->LorsY());
//       }                                                             //digits loop
// 
//       
      if(iCh==6) pC->cd(1); if(iCh==5) pC->cd(2);
      if(iCh==4) pC->cd(4); if(iCh==3) pC->cd(5); if(iCh==2) pC->cd(6);
                            if(iCh==1) pC->cd(8); if(iCh==0) pC->cd(9);
                          
      AliHMPIDDigit::DrawPc();  pTrkMap[iCh]->Draw(); pDigMap[iCh]->Draw();
    }//chambers loop
//    pC->Update();
//    pC->Modified();
//    if(iEvt<iEvtTo) {gPad->WaitPrimitive();pC->Clear();}
    
    
    
  }//events loop
  delete pEsd;  pEsdFl->Close();//close AliESDs.root
//  rl->UnloadDigits();
}//Display()
