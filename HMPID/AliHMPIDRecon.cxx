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

#include "AliHMPIDRecon.h"   //class header
#include "AliHMPIDParam.h"   //CkovAngle()
#include "AliHMPIDCluster.h" //CkovAngle()
#include <TMinuit.h>         //FitEllipse()
#include <TRotation.h>       //TracePhot()
#include <TH1D.h>            //HoughResponse()
#include <TClonesArray.h>    //CkovAngle()
#include <AliESDtrack.h>     //CkovAngle()

const Double_t AliHMPIDRecon::fgkRadThick=1.5;
const Double_t AliHMPIDRecon::fgkWinThick=0.5;
const Double_t AliHMPIDRecon::fgkGapThick=8.0;
const Double_t AliHMPIDRecon::fgkWinIdx  =1.5787;
const Double_t AliHMPIDRecon::fgkGapIdx  =1.0005;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecon::AliHMPIDRecon():TTask("RichRec","RichPat"),
  fRadNmean(1.292),  
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
//hidden algorithm
  fMipX=fMipY=fThTrkFit=fPhTrkFit=fCkovFit=fMipQ=fRadX=fRadY=-999;
  fIdxMip=fNClu=0;
  fCkovSig2=0;
  for (Int_t i=0; i<100; i++) {
    fXClu[i] = fYClu[i] = 0;
    fClCk[i] = kTRUE;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::CkovAngle(AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean,Double_t qthre)
{
// Pattern recognition method based on Hough transform
// Arguments:   pTrk     - track for which Ckov angle is to be found
//              pCluLst  - list of clusters for this chamber   
//   Returns:            - track ckov angle, [rad], 
    
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  
  if(pCluLst->GetEntries()>pParam->MultCut()) fIsWEIGHT = kTRUE; // offset to take into account bkg in reconstruction
  else                                        fIsWEIGHT = kFALSE;

  Float_t xRa,yRa,th,ph;
  pTrk->GetHMPIDtrk(xRa,yRa,th,ph);        //initialize this track: th and ph angles at middle of RAD 
  SetTrack(xRa,yRa,th,ph);

  fRadNmean=nmean;

  Float_t dMin=999,mipX=-1,mipY=-1;Int_t chId=-1,mipId=-1,mipQ=-1;                                                                           
  fPhotCnt=0;                                                      
  for (Int_t iClu=0; iClu<pCluLst->GetEntriesFast();iClu++){//clusters loop
    AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(iClu);                       //get pointer to current cluster    
    chId=pClu->Ch();
    if(pClu->Q()>qthre){                                                                      //charge compartible with MIP clusters      
      Float_t dX=fPc.X()-pClu->X(),dY=fPc.Y()-pClu->Y(),d =TMath::Sqrt(dX*dX+dY*dY);          //distance between current cluster and intersection point
      if( d < dMin) {mipId=iClu; dMin=d;mipX=pClu->X();mipY=pClu->Y();mipQ=(Int_t)pClu->Q();} //current cluster is closer, overwrite data for min cluster
    }else{                                                                                    //charge compatible with photon cluster
      Double_t thetaCer,phiCer;
      if(FindPhotCkov(pClu->X(),pClu->Y(),thetaCer,phiCer)){                                  //find ckov angle for this  photon candidate
        fPhotCkov[fPhotCnt]=thetaCer;                                                         //actual theta Cerenkov (in TRS)
        fPhotPhi [fPhotCnt]=phiCer;                                                           //actual phi   Cerenkov (in TRS): -pi to come back to "unusual" ref system (X,Y,-Z)
	//PH        Printf("photon n. %i reconstructed theta = %f",fPhotCnt,fPhotCkov[fPhotCnt]);
        fPhotCnt++;                                                                           //increment counter of photon candidates
      }
    }
  }//clusters loop
  if(fPhotCnt<=3) pTrk->SetHMPIDsignal(kNoPhotAccept);                                        //no reconstruction with <=3 photon candidates
  Int_t iNacc=FlagPhot(HoughResponse());                                                      //flag photons according to individual theta ckov with respect to most probable
  pTrk->SetHMPIDmip(mipX,mipY,mipQ,iNacc);                                                    //store mip info 

  if(mipId==-1)              {pTrk->SetHMPIDsignal(kMipQdcCut);  return;}                     //no clusters with QDC more the threshold at all
  if(dMin>pParam->DistCut()) {pTrk->SetHMPIDsignal(kMipDistCut); return;}                     //closest cluster with enough charge is still too far from intersection
  pTrk->SetHMPIDcluIdx(chId,mipId);                                                           //set index of cluster
  if(iNacc<1){
    pTrk->SetHMPIDsignal(kNoPhotAccept);                                                      //no photon candidates is accepted
  }
  else {
    pTrk->SetHMPIDsignal(FindRingCkov(pCluLst->GetEntries()));                                //find best Theta ckov for ring i.e. track
    pTrk->SetHMPIDchi2(fCkovSigma2);                                                          //errors squared
  }

}//CkovAngle()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::FindPhotCkov(Double_t cluX,Double_t cluY,Double_t &thetaCer,Double_t &phiCer)
{
// Finds Cerenkov angle  for this photon candidate
// Arguments: cluX,cluY - position of cadidate's cluster  
// Returns: Cerenkov angle 

  TVector3 dirCkov;
  
  Double_t zRad= -0.5*fgkRadThick-0.5*fgkWinThick;                   //z position of middle of RAD
  TVector3 rad(fTrkPos.X(),fTrkPos.Y(),zRad);                        //impact point at middle of RAD
  TVector3  pc(cluX,cluY,0.5*fgkWinThick+fgkGapIdx);                 //mip at PC
  Double_t cluR = TMath::Sqrt((cluX-fTrkPos.X())*(cluX-fTrkPos.X())+
                              (cluY-fTrkPos.Y())*(cluY-fTrkPos.Y()));//ref. distance impact RAD-CLUSTER   
  Double_t phi=(pc-rad).Phi();                                       //phi of photon
    
  Double_t ckov1=0;
  Double_t ckov2=0.75+fTrkDir.Theta();                        //start to find theta cerenkov in DRS
  const Double_t kTol=0.01;
  Int_t iIterCnt = 0;
  while(1){
    if(iIterCnt>=50) return kFALSE;
    Double_t ckov=0.5*(ckov1+ckov2);
    dirCkov.SetMagThetaPhi(1,ckov,phi);
    TVector2 posC=TraceForward(dirCkov);                      //trace photon with actual angles
    Double_t dist=cluR-(posC-fTrkPos).Mod();                  //get distance between trial point and cluster position
    if(posC.X()==-999) dist = - 999;                          //total reflection problem
    iIterCnt++;                                               //counter step
    if     (dist> kTol) ckov1=ckov;                           //cluster @ larger ckov
    else if(dist<-kTol) ckov2=ckov;                           //cluster @ smaller ckov
    else{                                                     //precision achived: ckov in DRS found
      dirCkov.SetMagThetaPhi(1,ckov,phi);                     //
      RecPhot(dirCkov,thetaCer,phiCer);                       //find ckov (in TRS:the effective Cherenkov angle!)
      return kTRUE;
    }
  }
}//FindPhotTheta()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliHMPIDRecon::TraceForward(TVector3 dirCkov)const
{
  //Trace forward a photon from (x,y) up to PC
  // Arguments: dirCkov photon vector in LORS
  //   Returns: pos of traced photon at PC
  TVector2 pos(-999,-999);
  Double_t thetaCer = dirCkov.Theta();
  if(thetaCer > TMath::ASin(1./fRadNmean))  return pos;         //total refraction on WIN-GAP boundary
  Double_t zRad= -0.5*fgkRadThick-0.5*fgkWinThick;              //z position of middle of RAD
  TVector3  posCkov(fTrkPos.X(),fTrkPos.Y(),zRad);              //RAD: photon position is track position @ middle of RAD 
  Propagate(dirCkov,posCkov,           -0.5*fgkWinThick);       //go to RAD-WIN boundary  
  Refract  (dirCkov,         fRadNmean,fgkWinIdx);              //RAD-WIN refraction
  Propagate(dirCkov,posCkov,            0.5*fgkWinThick);       //go to WIN-GAP boundary
  Refract  (dirCkov,         fgkWinIdx,fgkGapIdx);              //WIN-GAP refraction
  Propagate(dirCkov,posCkov,0.5*fgkWinThick+fgkGapThick);       //go to PC
  pos.Set(posCkov.X(),posCkov.Y());
  return pos;
}//TraceForward()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::RecPhot(TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer)
{
  //Theta Cerenkov reconstruction 
  // Arguments: (x,y) of initial point in LORS, dirCkov photon vector in LORS
  //   Returns: thetaCer theta cerenkov reconstructed
//  TVector3 dirTrk;
//  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi());
//  Double_t thetaCer = TMath::ACos(dirCkov*dirTrk);
  TRotation mtheta;   mtheta.RotateY(- fTrkDir.Theta());
  TRotation mphi;       mphi.RotateZ(- fTrkDir.Phi());
  TRotation mrot=mtheta*mphi;
  TVector3 dirCkovTRS;
  dirCkovTRS=mrot*dirCkov;
  phiCer  = dirCkovTRS.Phi();                                          //actual value of the phi of the photon
  thetaCer= dirCkovTRS.Theta();                                        //actual value of thetaCerenkov of the photon
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FindRingArea(Double_t ckovAngMin,Double_t ckovAngMax)const
{
// Find area between 2 cerenkov angles in the PC acceptance
// Arguments: ckovAngMin - cerenkov angle Min    
// Arguments: ckovAngMax - cerenkov angle Max
//   Returns: area of the ring in cm^2 for given theta ckov
   
  const Int_t kN=100;
  Int_t np=0;
  Double_t xP[2*kN],yP[2*kN];
  
  Double_t area=0;
//---  find points from first ring
  for(Int_t i=0;i<kN;i++){
    TVector2 pos=TracePhot(ckovAngMin,Double_t(TMath::TwoPi()*(i+1)/kN));//trace the next photon
    if(pos.X()==-999) continue;                                       //no area: open ring
    if(AliHMPIDParam::IsInside(pos.X(),pos.Y(),0)) continue;
    xP[np] = pos.X();
    yP[np] = pos.Y();
    np++;                                                              
  }
//--- find points from last ring
  for(Int_t i=kN-1;i>=0;i--){
    TVector2 pos=TracePhot(ckovAngMax,Double_t(TMath::TwoPi()*(i+1)/kN));//trace the next photon
    if(pos.X()==-999) continue;                                       
    if(AliHMPIDParam::IsInside(pos.X(),pos.Y(),0)) continue;
    xP[np] = pos.X();
    yP[np] = pos.Y();
    np++;                                                              
  }
//--calculate delta area from array of points...
//  for(Int_t i=0;i<np;i++) {
  area = 1;    
  area*=0.5;
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
      if(fPhotCkov[i]<ckovMin) ckovMin=fPhotCkov[i];                         //find max and min Theta ckov from all candidates within probable window
      if(fPhotCkov[i]>ckovMax) ckovMax=fPhotCkov[i]; 
      weightThetaCerenkov += fPhotCkov[i]*fPhotWei[i];
      wei += fPhotWei[i];                                                    //collect weight as sum of all candidate weghts   
      
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

// Photon Flag:  Flag = 0 initial set; 
//               Flag = 1 good candidate (charge compatible with photon); 
//               Flag = 2 photon used for the ring;
  
  Int_t steps = (Int_t)((ckov )/ fDTheta); //how many times we need to have fDTheta to fill the distance between 0  and thetaCkovHough

  Double_t tmin = (Double_t)(steps - 1)*fDTheta;
  Double_t tmax = (Double_t)(steps)*fDTheta;
  Double_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;  tmax = tavg + 0.5*fWindowWidth;

  Int_t iInsideCnt = 0; //count photons which Theta ckov inside the window
  for(Int_t i=0;i<fPhotCnt;i++){//photon candidates loop
    fPhotFlag[i] = 0;
    if(fPhotCkov[i] >= tmin && fPhotCkov[i] <= tmax)	{ 
      fPhotFlag[i]=2;	  
      iInsideCnt++;
    }
  }
  return iInsideCnt;
}//FlagPhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliHMPIDRecon::TracePhot(Double_t ckovThe,Double_t ckovPhi)const
{
// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in DRS, [rad]    
//   Returns: distance between photon point on PC and track projection  
  TRotation mtheta;   mtheta.RotateY(fTrkDir.Theta());
  TRotation mphi;       mphi.RotateZ(fTrkDir.Phi());  
  TRotation mrot=mphi*mtheta;
  TVector3  dirCkov,dirCkovTors;   

  dirCkovTors.SetMagThetaPhi(1,ckovThe,ckovPhi);                    //initially photon is directed according to requested ckov angle
  dirCkov=mrot*dirCkovTors;                                         //now we know photon direction in LORS
  return TraceForward(dirCkov);
}//TracePhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::Propagate(const TVector3 dir,TVector3 &pos,Double_t z)const
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
  if(TMath::Abs(sinref)>1.) dir.SetXYZ(-999,-999,-999);
  else             dir.SetTheta(TMath::ASin(sinref));
}//Refract()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::HoughResponse()
{
//
//    fIdxMip = mipId;

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
      Double_t diffArea = FindRingArea(lowerlimit,upperlimit);
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
  delete phots;delete photsw;delete resultw; // Reset and delete objects
  
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
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*fRadNmean);
  
  if(trkBeta > 1) trkBeta = 1;                 //protection against bad measured thetaCer  
  if(trkBeta < 0) trkBeta = 0.0001;            //

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

  Double_t sint     = TMath::Sin(fTrkDir.Theta());
  Double_t cost     = TMath::Cos(fTrkDir.Theta());
  Double_t sinf     = TMath::Sin(fTrkDir.Phi());
  Double_t cosf     = TMath::Cos(fTrkDir.Phi());
  Double_t sinfd    = TMath::Sin(phiDelta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t k = 1.-fRadNmean*fRadNmean+alpha*alpha/(betaM*betaM);                            // formula (after 8 in the text)
  if (k<0) return 1e10;
  Double_t mu =sint*sinf+tantheta*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
  Double_t e  =sint*cosf+tantheta*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

  Double_t kk = betaM*TMath::Sqrt(k)/(fgkGapThick*alpha);                                   // formula (6) and (7)
  Double_t dtdxc = kk*(k*(cosfd*cosf-cost*sinfd*sinf)-(alpha*mu/(betaM*betaM))*sint*sinfd); // formula (6)           
  Double_t dtdyc = kk*(k*(cosfd*sinf+cost*sinfd*cosf)+(alpha* e/(betaM*betaM))*sint*sinfd); // formula (7)            pag.4

  Double_t errX = 0.2,errY=0.25;                                                            //end of page 7
  return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
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

  Double_t sint     = TMath::Sin(fTrkDir.Theta());
  Double_t cost     = TMath::Cos(fTrkDir.Theta());
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t dtdn = cost*fRadNmean*betaM*betaM/(alpha*tantheta);                              // formula (12)
            
//  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
  Double_t f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

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

  Double_t sint     = TMath::Sin(fTrkDir.Theta());
  Double_t cost     = TMath::Cos(fTrkDir.Theta());
  Double_t sinf     = TMath::Sin(fTrkDir.Phi());
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t costheta = TMath::Cos(thetaC);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                  // formula (11)
  
  Double_t k = 1.-fRadNmean*fRadNmean+alpha*alpha/(betaM*betaM);                             // formula (after 8 in the text)
  if (k<0) return 1e10;

  Double_t eTr = 0.5*fgkRadThick*betaM*TMath::Sqrt(k)/(fgkGapThick*alpha);                   // formula (14)
  Double_t lambda = 1.-sint*sint*sinf*sinf;                                                  // formula (15)

  Double_t c1 = 1./(1.+ eTr*k/(alpha*alpha*costheta*costheta));                              // formula (13.a)
  Double_t c2 = betaM*TMath::Power(k,1.5)*tantheta*lambda/(fgkGapThick*alpha*alpha);         // formula (13.b)
  Double_t c3 = (1.+eTr*k*betaM*betaM)/((1+eTr)*alpha*alpha);                                // formula (13.c)
  Double_t c4 = TMath::Sqrt(k)*tantheta*(1-lambda)/(fgkGapThick*betaM);                      // formula (13.d)
  Double_t dtdT = c1 * (c2+c3*c4);
  Double_t trErr = fgkRadThick/(TMath::Sqrt(12.)*cost);

  return trErr*dtdT;
}//SigGeom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// From here HTA....
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::CkovHiddenTrk(AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean, Double_t qthre)
{
// Pattern recognition method without any infos from tracking:HTA (Hidden Track Algorithm)...
// The method finds in the chmber the cluster with the highest charge
// compatibile with a MIP, then the strategy is applied
// Arguments:  pTrk     - pointer to ESD track
//             pCluLs   - list of clusters for a given chamber 
//             nmean    - mean freon ref. index
//   Returns:           - 0=ok,1=not fitted 
  
  fRadNmean=nmean;

  if(pCluLst->GetEntriesFast()>100) return kFALSE;                                            //boundary check for CluX,CluY...
  Float_t mipX=-1,mipY=-1;Int_t mipId=-1,mipQ=-1;                                                                           
  Double_t qRef = 0;
  Int_t nCh=0;
  for (Int_t iClu=0;iClu<pCluLst->GetEntriesFast();iClu++){                                   //clusters loop
    AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(iClu);                       //get pointer to current cluster    
    nCh = pClu->Ch();
    fXClu[iClu] = pClu->X();fYClu[iClu] = pClu->Y();                                          //store x,y for fitting procedure
    fClCk[iClu] = kTRUE;                                                                      //all cluster are accepted at this stage to be reconstructed
    if(pClu->Q()>qRef){                                                                       //searching the highest charge to select a MIP      
      qRef = pClu->Q();
      mipId=iClu; mipX=pClu->X();mipY=pClu->Y();mipQ=(Int_t)pClu->Q();
    }                                                                                    
  }//clusters loop

  fNClu = pCluLst->GetEntriesFast();
  if(qRef>qthre){                                                                     //charge compartible with MIP clusters
    fIdxMip = mipId;
    fClCk[mipId] = kFALSE;
    fMipX = mipX; fMipY=mipY; fMipQ = qRef;
    if(!DoRecHiddenTrk(pCluLst)) {
      pTrk->SetHMPIDsignal(kNoPhotAccept);
      return kFALSE;
    }                                                                           //Do track and ring reconstruction,if problems returns 1
    pTrk->SetHMPIDtrk(fRadX,fRadY,fThTrkFit,fPhTrkFit);                                        //store track intersection info
    pTrk->SetHMPIDmip(fMipX,fMipY,(Int_t)fMipQ,fNClu);                                         //store mip info 
    pTrk->SetHMPIDcluIdx(nCh,fIdxMip);                                                         //set cham number and index of cluster
    pTrk->SetHMPIDsignal(fCkovFit);                                                            //find best Theta ckov for ring i.e. track
    pTrk->SetHMPIDchi2(fCkovSig2);                                                             //errors squared
//    Printf(" n clusters tot %i accepted %i",pCluLst->GetEntriesFast(),fNClu);
    return kTRUE;
  }
  
  return kFALSE;
}//CkovHiddenTrk()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::DoRecHiddenTrk(TClonesArray *pCluLst)
{
// Pattern recognition method without any infos from tracking...
// First a preclustering filter to avoid part of the noise
// Then only ellipsed-rings are fitted (no possibility, 
// for the moment, to reconstruct very inclined tracks)
// Finally a fitting with (th,ph) free, starting by very close values
// previously evaluated.
// Arguments:   none
//   Returns:   none
  Double_t phiRec;
  if(!CluPreFilter(pCluLst)) {return kFALSE;}
  if(!FitEllipse(phiRec)) {return kFALSE;}
  Int_t nClTmp1 = pCluLst->GetEntriesFast()-1;  //minus MIP...
  Int_t nClTmp2 = 0;
  while(nClTmp1 != nClTmp2){
    SetNClu(pCluLst->GetEntriesFast());
    if(!FitFree(phiRec)) {return kFALSE;}
    nClTmp2 = NClu();
    if(nClTmp2!=nClTmp1) {nClTmp1=nClTmp2;nClTmp2=0;}
  }
  fNClu = nClTmp2;
  return kTRUE;
}//DoRecHiddenTrk()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::CluPreFilter(TClonesArray *pCluLst)
{
// Filter of bkg clusters
// based on elliptical-shapes...
//
  if(pCluLst->GetEntriesFast()>50||pCluLst->GetEntriesFast()<4) return kFALSE; 
  else return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::FitEllipse(Double_t &phiRec)
{
//Fit a set of clusters with an analitical conical section function:
  //
  // Ax^2 + B*y^2 + 2Hxy + 2Gx + 2Fy + 1 = 0   ---> conical section
  //
  //  H*H - A*B > 0 hyperbola
  //            < 0 ellipse
  //            = 0 parabola
  //
  // tan 2alfa = 2H/(A-B)  alfa=angle of rotation
  //
  // coordinate of the centre of the conical section:
  //   x = x' + a
  //   y = y' + b
  //
  //       HF - BG
  //  a = ---------
  //       AB - H^2
  //
  //       HG - AF
  //  b = --------
  //       AB - H^2
  Double_t cA,cB,cF,cG,cH;
  Double_t aArg=-1;      Int_t iErrFlg;                                                //tmp vars for TMinuit

  if(!gMinuit) gMinuit = new TMinuit(5);                                               //init MINUIT with this number of parameters (5 params)
  gMinuit->mncler();                                                                   // reset Minuit list of paramters
  gMinuit->SetObjectFit((TObject*)this);  gMinuit->SetFCN(AliHMPIDRecon::FunMinEl);    //set fit function
  gMinuit->mnexcm("SET PRI",&aArg,1,iErrFlg);                                          //suspend all printout from TMinuit 
  gMinuit->mnexcm("SET NOW",&aArg,0,iErrFlg);                                          //suspend all warning printout from TMinuit
  
  Double_t d1,d2,d3;
  TString sName;

  gMinuit->mnparm(0," A ",1,0.01,0,0,iErrFlg);
  gMinuit->mnparm(1," B ",1,0.01,0,0,iErrFlg);
  gMinuit->mnparm(2," H ",1,0.01,0,0,iErrFlg);
  gMinuit->mnparm(3," G ",1,0.01,0,0,iErrFlg);
  gMinuit->mnparm(4," F ",1,0.01,0,0,iErrFlg);

  gMinuit->mnexcm("SIMPLEX",&aArg,0,iErrFlg);
  gMinuit->mnexcm("MIGRAD" ,&aArg,0,iErrFlg);   
  gMinuit->mnpout(0,sName,cA,d1,d2,d3,iErrFlg);
  gMinuit->mnpout(1,sName,cB,d1,d2,d3,iErrFlg);
  gMinuit->mnpout(2,sName,cH,d1,d2,d3,iErrFlg);
  gMinuit->mnpout(3,sName,cG,d1,d2,d3,iErrFlg);
  gMinuit->mnpout(4,sName,cF,d1,d2,d3,iErrFlg);
  delete gMinuit;

  Double_t i2 = cA*cB-cH*cH;                                       //quartic invariant : i2 > 0  ellipse, i2 < 0 hyperbola
  if(i2<=0) return kFALSE;
  Double_t aX = (cH*cF-cB*cG)/i2;                                  //x centre of the canonical section 
  Double_t bY = (cH*cG-cA*cF)/i2;                                  //y centre of the canonical section 
  Double_t alfa1 = TMath::ATan(2*cH/(cA-cB));                      //alpha = angle of rotation of the conical section
  if(alfa1<0) alfa1+=TMath::Pi(); 
  alfa1*=0.5;
//  Double_t alfa2 = alfa1+TMath::Pi();
  Double_t phiref = TMath::ATan2(bY-fMipY,aX-fMipX);               //evaluate in a unique way the angle of rotation comparing it
  if(phiref<0) phiref+=TMath::TwoPi();                             //with the vector that points to the centre from the mip 
  if(i2<0) phiref+=TMath::Pi();
  if(phiref>TMath::TwoPi()) phiref-=TMath::TwoPi();

//  Printf(" alfa1 %f",alfa1*TMath::RadToDeg());
//  Printf(" alfa2 %f",alfa2*TMath::RadToDeg());
//  Printf(" firef %f",phiref*TMath::RadToDeg());
//  if(TMath::Abs(alfa1-phiref)<TMath::Abs(alfa2-phiref)) phiRec = alfa1; else phiRec = alfa2;  
  
//  Printf("FitEllipse: phi reconstructed %f",phiRec*TMath::RadToDeg());
  phiRec=phiref;
  return kTRUE;
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::FitFree(Double_t phiRec)
{
// Fit performed by minimizing RMS/sqrt(n) of the
// photons reconstructed. First phi is fixed and theta
// is fouond, then (th,ph) of the track
// as free parameters
// Arguments:    PhiRec phi of the track
//   Returns:    none
  Double_t aArg=-1;  Int_t iErrFlg;                                                    //tmp vars for TMinuit
  if(!gMinuit) gMinuit = new TMinuit(2);                                               //init MINUIT with this number of parameters (5 params)
  gMinuit->mncler();                                                                   // reset Minuit list of paramters
  gMinuit->SetObjectFit((TObject*)this);  gMinuit->SetFCN(AliHMPIDRecon::FunMinPhot);  //set fit function
  gMinuit->mnexcm("SET PRI",&aArg,1,iErrFlg);                                          //suspend all printout from TMinuit 
  gMinuit->mnexcm("SET NOW",&aArg,0,iErrFlg);                                          //suspend all warning printout from TMinuit
  
  Double_t d1,d2,d3;
  TString sName;
  Double_t th,ph;
  
  gMinuit->mnparm(0," theta ",  0.01,0.01,0,TMath::PiOver2(),iErrFlg);
  gMinuit->mnparm(1," phi   ",phiRec,0.01,0,TMath::TwoPi()  ,iErrFlg);
  
  gMinuit->FixParameter(1);
  gMinuit->mnexcm("SIMPLEX" ,&aArg,0,iErrFlg);   
  gMinuit->mnexcm("MIGRAD"  ,&aArg,0,iErrFlg);
  gMinuit->Release(1);  
  gMinuit->mnexcm("MIGRAD"  ,&aArg,0,iErrFlg);
  
  gMinuit->mnpout(0,sName,th,d1,d2,d3,iErrFlg);
  gMinuit->mnpout(1,sName,ph,d1,d2,d3,iErrFlg);   

  Double_t outPar[2] = {th,ph}; Double_t g; Double_t f;Int_t flag = 3;
  gMinuit->Eval(2, &g, f, outPar,flag);  

  SetTrkFit(th,ph);
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDRecon::FunConSect(Double_t *c,Double_t x,Double_t y)
{
  return c[0]*x*x+c[1]*y*y+2*c[2]*x*y+2*c[3]*x+2*c[4]*y+1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::FunMinEl(Int_t &/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t /* */)
{
  AliHMPIDRecon *pRec=(AliHMPIDRecon*)gMinuit->GetObjectFit();
  Double_t minFun = 0;
  Int_t np = pRec->NClu();
  for(Int_t i=0;i<np;i++) {
    if(i==pRec->IdxMip()) continue;
    Double_t el = pRec->FunConSect(par,pRec->XClu(i),pRec->YClu(i));
    minFun +=el*el;
  }
  f = minFun;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::FunMinPhot(Int_t &/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t iflag)
{
  AliHMPIDRecon *pRec=(AliHMPIDRecon*)gMinuit->GetObjectFit();
  Double_t sizeCh = 0.5*fgkRadThick+fgkWinThick+fgkGapThick;
  Double_t thTrk = par[0]; 
  Double_t phTrk = par[1]; 
  Double_t xrad = pRec->MipX() - sizeCh*TMath::Tan(thTrk)*TMath::Cos(phTrk);
  Double_t yrad = pRec->MipY() - sizeCh*TMath::Tan(thTrk)*TMath::Sin(phTrk);
  pRec->SetRadXY(xrad,yrad);
  pRec->SetTrack(xrad,yrad,thTrk,phTrk);

  Double_t meanCkov =0;
  Double_t meanCkov2=0;
  Double_t thetaCer,phiCer;
  Int_t nClAcc = 0;
  Int_t nClTot=pRec->NClu();
    
  for(Int_t i=0;i<nClTot;i++) {
    if(!(pRec->ClCk(i))) continue;
    pRec->FindPhotCkov(pRec->XClu(i),pRec->YClu(i),thetaCer,phiCer);  
    meanCkov  += thetaCer;
    meanCkov2 += thetaCer*thetaCer;
    nClAcc++;
  }
  if(nClAcc==0) {f=999;return;}
  meanCkov/=nClAcc;
  Double_t rms = (meanCkov2 - meanCkov*meanCkov*nClAcc)/nClAcc;
  if(rms<0) Printf(" rms2 = %f, strange!!!",rms);
  rms = TMath::Sqrt(rms);
  f = rms/TMath::Sqrt(nClAcc);
  
  
  if(iflag==3) {
    Printf("FunMinPhot before: photons candidates %i used %i",nClTot,nClAcc);
    nClAcc = 0;
    Double_t meanCkov1=0;
    Double_t meanCkov2=0;
    for(Int_t i=0;i<nClTot;i++) {
      if(!(pRec->ClCk(i))) continue;
      pRec->FindPhotCkov(pRec->XClu(i),pRec->YClu(i),thetaCer,phiCer);  
      if(TMath::Abs(thetaCer-meanCkov)<2*rms) {
        meanCkov1 += thetaCer;
        meanCkov2 += thetaCer*thetaCer;
        nClAcc++;
      } else pRec->SetClCk(i,kFALSE);
    }
    meanCkov1/=nClAcc;
    Double_t rms2 = (meanCkov2 - meanCkov*meanCkov*nClAcc)/nClAcc;
    Printf("FunMinPhot after: photons candidates %i used %i thetaCer %f",nClTot,nClAcc,meanCkov1);
    pRec->SetCkovFit(meanCkov1);
    pRec->SetCkovSig2(rms2);
    pRec->SetNClu(nClAcc);
  }
}//FunMinPhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ended Hidden track algorithm....
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
