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
#include "AliHMPIDCluster.h" //CkovAngle()
#include <TRotation.h>       //TracePhot()
#include <TH1D.h>            //HoughResponse()
#include <TClonesArray.h>    //CkovAngle()
#include <AliESDtrack.h>     //CkovAngle()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecon::AliHMPIDRecon():TTask("RichRec","RichPat")
{
//..
//init of data members
//..
  
  fPhotCnt  = -1;
  fPhotFlag = 0x0;
  fPhotCkov = 0x0;
  fPhotPhi  = 0x0;
  fPhotWei  = 0x0;
  fCkovSigma2 = 0;
  fIsWEIGHT = kFALSE;
  fDTheta   = 0.001;
  fWindowWidth = 0.045;
  fTrkDir = TVector3(0,0,1); // init just for test
  fTrkPos = TVector2(30,40); // init just for test
  
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  fParam = pParam;
  
  fParam->SetRefIdx(fParam->MeanIdxRad()); // initialization of ref index to a default one
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::InitVars(Int_t n)
{
//..
//Init some variables
//..
  if(n<0) return;
  fPhotFlag = new Int_t[n];
  fPhotCkov = new Double_t[n];
  fPhotPhi  = new Double_t[n];
  fPhotWei  = new Double_t[n];
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::DeleteVars()
{
//..
//Delete variables
//..
  delete fPhotFlag;
  delete fPhotCkov;
  delete fPhotPhi;
  delete fPhotWei;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::CkovAngle(AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean,Double_t qthre)
{
// Pattern recognition method based on Hough transform
// Arguments:   pTrk     - track for which Ckov angle is to be found
//              pCluLst  - list of clusters for this chamber   
//   Returns:            - track ckov angle, [rad], 
    
  Int_t nClusTot = pCluLst->GetEntries();
  if(nClusTot>fParam->MultCut()) fIsWEIGHT = kTRUE; // offset to take into account bkg in reconstruction
  else                           fIsWEIGHT = kFALSE;

  InitVars(nClusTot);
  
  Float_t xRa,yRa,th,ph;
  pTrk->GetHMPIDtrk(xRa,yRa,th,ph);        //initialize this track: th and ph angles at middle of RAD 
  SetTrack(xRa,yRa,th,ph);

  fParam->SetRefIdx(nmean);

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
  fMipPos.Set(mipX,mipY);
  if(fPhotCnt<=3) pTrk->SetHMPIDsignal(kNoPhotAccept);                                        //no reconstruction with <=3 photon candidates
  Int_t iNacc=FlagPhot(HoughResponse());                                                      //flag photons according to individual theta ckov with respect to most probable
  pTrk->SetHMPIDmip(mipX,mipY,mipQ,iNacc);                                                    //store mip info 

  if(mipId==-1)              {pTrk->SetHMPIDsignal(kMipQdcCut);  return;}                     //no clusters with QDC more the threshold at all
  if(dMin>fParam->DistCut()) {pTrk->SetHMPIDsignal(kMipDistCut); return;}                     //closest cluster with enough charge is still too far from intersection
  pTrk->SetHMPIDcluIdx(chId,mipId);                                                           //set index of cluster
  if(iNacc<1){
    pTrk->SetHMPIDsignal(kNoPhotAccept);                                                      //no photon candidates is accepted
  }
  else {
    pTrk->SetHMPIDsignal(FindRingCkov(pCluLst->GetEntries()));                                //find best Theta ckov for ring i.e. track
    pTrk->SetHMPIDchi2(fCkovSigma2);                                                          //errors squared
  }

  DeleteVars();
}//CkovAngle()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRecon::FindPhotCkov(Double_t cluX,Double_t cluY,Double_t &thetaCer,Double_t &phiCer)
{
// Finds Cerenkov angle  for this photon candidate
// Arguments: cluX,cluY - position of cadidate's cluster  
// Returns: Cerenkov angle 

  TVector3 dirCkov;
  
  Double_t zRad= -0.5*fParam->RadThick()-0.5*fParam->WinThick();     //z position of middle of RAD
  TVector3 rad(fTrkPos.X(),fTrkPos.Y(),zRad);                        //impact point at middle of RAD
  TVector3  pc(cluX,cluY,0.5*fParam->WinThick()+fParam->GapIdx());   //mip at PC
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
  if(thetaCer > TMath::ASin(1./fParam->GetRefIdx())) return pos;          //total refraction on WIN-GAP boundary
  Double_t zRad= -0.5*fParam->RadThick()-0.5*fParam->WinThick();          //z position of middle of RAD
  TVector3  posCkov(fTrkPos.X(),fTrkPos.Y(),zRad);                        //RAD: photon position is track position @ middle of RAD 
  Propagate(dirCkov,posCkov,           -0.5*fParam->WinThick());          //go to RAD-WIN boundary  
  Refract  (dirCkov,         fParam->GetRefIdx(),fParam->WinIdx());       //RAD-WIN refraction
  Propagate(dirCkov,posCkov,            0.5*fParam->WinThick());          //go to WIN-GAP boundary
  Refract  (dirCkov,         fParam->WinIdx(),fParam->GapIdx());          //WIN-GAP refraction
  Propagate(dirCkov,posCkov,0.5*fParam->WinThick()+fParam->GapThick());   //go to PC
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
Double_t AliHMPIDRecon::FindRingArea(Double_t ckovAng)const
{
// Find area covered in the PC acceptance
// Arguments: ckovAng - cerenkov angle     
//   Returns: area of the ring in cm^2 for given theta ckov
   
  const Int_t kN=50;
  TVector2 pos1;
  Double_t area=0;
  Bool_t first=kFALSE;
  for(Int_t i=0;i<kN;i++){
   if(!first) {
     pos1=TracePhot(ckovAng,Double_t(TMath::TwoPi()*(i+1)/kN));                                     //find a good trace for the first photon
     if(pos1.X()==-999) continue;                                                                   //no area: open ring 	          
     if(!fParam->IsInside(pos1.X(),pos1.Y(),0)) pos1 = IntWithEdge(fMipPos,pos1);                   // ffind the very first intersection...
     first=kTRUE;
     continue;
   }
   TVector2 pos2=TracePhot(ckovAng,Double_t(TMath::TwoPi()*(i+1)/kN));                              //trace the next photon
   if(pos2.X()==-999) continue;                                                                     //no area: open ring 	     
   if(!fParam->IsInside(pos2.X(),pos2.Y(),0)) {
     pos2 = IntWithEdge(fMipPos,pos2);
   }
   area+=TMath::Abs((pos1-fMipPos).X()*(pos2-fMipPos).Y()-(pos1-fMipPos).Y()*(pos2-fMipPos).X());   //add area of the triangle... 	     
   pos1 = pos2;
  }
//---  find points from ring
  area*=0.5;
  return area;
}//FindRingArea()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliHMPIDRecon::IntWithEdge(TVector2 p1,TVector2 p2)const
{
// It finds the intersection of the line for 2 points traced as photons
// and the edge of a given PC
// Arguments: 2 points obtained tracing the photons
//   Returns: intersection point with detector (PC) edges

  Double_t xmin = (p1.X()<p2.X())? p1.X():p2.X(); 
  Double_t xmax = (p1.X()<p2.X())? p2.X():p1.X(); 
  Double_t ymin = (p1.Y()<p2.Y())? p1.Y():p2.Y(); 
  Double_t ymax = (p1.Y()<p2.Y())? p2.Y():p1.Y(); 
  
  Double_t m = TMath::Tan((p2-p1).Phi());
  TVector2 pint;
  //intersection with low  X
  pint.Set((Double_t)(p1.X() + (0-p1.Y())/m),0.);
  if(pint.X()>=0 && pint.X()<=fParam->SizeAllX() &&
     pint.X()>=xmin && pint.X()<=xmax            &&
     pint.Y()>=ymin && pint.Y()<=ymax) return pint;
  //intersection with high X  
  pint.Set((Double_t)(p1.X() + (fParam->SizeAllY()-p1.Y())/m),(Double_t)(fParam->SizeAllY()));
  if(pint.X()>=0 && pint.X()<=fParam->SizeAllX() &&
     pint.X()>=xmin && pint.X()<=xmax            &&
     pint.Y()>=ymin && pint.Y()<=ymax) return pint;
  //intersection with left Y  
  pint.Set(0.,(Double_t)(p1.Y() + m*(0-p1.X())));
  if(pint.Y()>=0 && pint.Y()<=fParam->SizeAllY() &&
     pint.Y()>=ymin && pint.Y()<=ymax            &&
     pint.X()>=xmin && pint.X()<=xmax) return pint;
  //intersection with righ Y  
  pint.Set((Double_t)(fParam->SizeAllX()),(Double_t)(p1.Y() + m*(fParam->SizeAllX()-p1.X())));
  if(pint.Y()>=0 && pint.Y()<=fParam->SizeAllY() &&
     pint.Y()>=ymin && pint.Y()<=ymax            &&
     pint.X()>=xmin && pint.X()<=xmax) return pint;
  return p1;
}//IntWithEdge()
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
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
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
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*fParam->GetRefIdx());
  
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
  Double_t k = 1.-fParam->GetRefIdx()*fParam->GetRefIdx()+alpha*alpha/(betaM*betaM);        // formula (after 8 in the text)
  if (k<0) return 1e10;
  Double_t mu =sint*sinf+tantheta*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
  Double_t e  =sint*cosf+tantheta*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

  Double_t kk = betaM*TMath::Sqrt(k)/(fParam->GapThick()*alpha);                            // formula (6) and (7)
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
  Double_t dtdn = cost*fParam->GetRefIdx()*betaM*betaM/(alpha*tantheta);                    // formula (12)
            
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
  
  Double_t k = 1.-fParam->GetRefIdx()*fParam->GetRefIdx()+alpha*alpha/(betaM*betaM);         // formula (after 8 in the text)
  if (k<0) return 1e10;

  Double_t eTr = 0.5*fParam->RadThick()*betaM*TMath::Sqrt(k)/(fParam->GapThick()*alpha);     // formula (14)
  Double_t lambda = 1.-sint*sint*sinf*sinf;                                                  // formula (15)

  Double_t c1 = 1./(1.+ eTr*k/(alpha*alpha*costheta*costheta));                              // formula (13.a)
  Double_t c2 = betaM*TMath::Power(k,1.5)*tantheta*lambda/(fParam->GapThick()*alpha*alpha);  // formula (13.b)
  Double_t c3 = (1.+eTr*k*betaM*betaM)/((1+eTr)*alpha*alpha);                                // formula (13.c)
  Double_t c4 = TMath::Sqrt(k)*tantheta*(1-lambda)/(fParam->GapThick()*betaM);               // formula (13.d)
  Double_t dtdT = c1 * (c2+c3*c4);
  Double_t trErr = fParam->RadThick()/(TMath::Sqrt(12.)*cost);

  return trErr*dtdT;
}//SigGeom()
