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
#include <AliESDfriendTrack.h>     //CkovAngle()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecon::AliHMPIDRecon():
  TTask("RichRec","RichPat"),
  fPhotCnt(-1),
  fPhotFlag(0x0),
  fPhotClusIndex(0x0),    
  fPhotCkov(0x0),
  fPhotPhi(0x0),
  fPhotWei(0x0),
  fCkovSigma2(0),
  fIsWEIGHT(kFALSE),
  fDTheta(0.001),
  fWindowWidth(0.045),
  fRingArea(0),
  fRingAcc(0),
  fTrkDir(0,0,1),  // Just for test
  fTrkPos(30,40),  // Just for test
  fMipPos(0,0),
  fPc(0,0),
  fParam(AliHMPIDParam::Instance())
{
//..
//init of data members
//..
  
  fParam->SetRefIdx(fParam->MeanIdxRad()); // initialization of ref index to a default one
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::InitVars(Int_t n)
{
//..
//Init some variables
//..
  if(n<=0) return;
  fPhotFlag = new Int_t[n];
  fPhotClusIndex  = new Int_t[n];
  fPhotCkov = new Double_t[n];
  fPhotPhi  = new Double_t[n];
  fPhotWei  = new Double_t[n];
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::DeleteVars()const
{
//..
//Delete variables
//..
  delete [] fPhotFlag;
  delete [] fPhotClusIndex;
  delete [] fPhotCkov;
  delete [] fPhotPhi;
  delete [] fPhotWei;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::CkovAngle(AliESDtrack *pTrk,TClonesArray *pCluLst,Int_t index,Double_t nmean,Float_t xRa,Float_t yRa)
{
// Pattern recognition method based on Hough transform
// Arguments:   pTrk     - track for which Ckov angle is to be found
//              pCluLst  - list of clusters for this chamber   
//   Returns:            - track ckov angle, [rad], 
       
  const Int_t nMinPhotAcc = 3;                      // Minimum number of photons required to perform the pattern recognition
  
  Int_t nClusTot = pCluLst->GetEntries();
  if(nClusTot>fParam->MultCut()) fIsWEIGHT = kTRUE; // offset to take into account bkg in reconstruction
  else                           fIsWEIGHT = kFALSE;

  InitVars(nClusTot);
  
  Float_t xPc,yPc,th,ph;
  pTrk->GetHMPIDtrk(xPc,yPc,th,ph);        //initialize this track: th and ph angles at middle of RAD 
  SetTrack(xRa,yRa,th,ph);

  fParam->SetRefIdx(nmean);

  Float_t mipX=-1,mipY=-1;
  Int_t chId=-1,mipQ=-1,sizeClu = -1;
  
  fPhotCnt=0;
  
  for (Int_t iClu=0; iClu<pCluLst->GetEntriesFast();iClu++){//clusters loop
    AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(iClu);                       //get pointer to current cluster    
    if(iClu == index) {                                                                       // this is the MIP! not a photon candidate: just store mip info
      mipX = pClu->X();
      mipY = pClu->Y();
      mipQ=(Int_t)pClu->Q();
      sizeClu=pClu->Size();
      continue;                                                             
    }
    chId=pClu->Ch();
    if(pClu->Q()>2*fParam->QCut()) continue;
    Double_t thetaCer,phiCer;
    if(FindPhotCkov(pClu->X(),pClu->Y(),thetaCer,phiCer)){                                    //find ckov angle for this  photon candidate
      fPhotCkov[fPhotCnt]=thetaCer;                                                           //actual theta Cerenkov (in TRS)
      fPhotPhi [fPhotCnt]=phiCer;
      fPhotClusIndex[fPhotCnt]=iClu;                                                             //actual phi   Cerenkov (in TRS): -pi to come back to "unusual" ref system (X,Y,-Z)
      fPhotCnt++;                                                                             //increment counter of photon candidates
    }
  }//clusters loop

  pTrk->SetHMPIDmip(mipX,mipY,mipQ,fPhotCnt);                                                 //store mip info in any case 
  pTrk->SetHMPIDcluIdx(chId,index+1000*sizeClu);                                              //set index of cluster
  
  if(fPhotCnt<=nMinPhotAcc) {                                                                 //no reconstruction with <=3 photon candidates
    pTrk->SetHMPIDsignal(kNoPhotAccept);                                                      //set the appropriate flag
    return;
  }
    
  fMipPos.Set(mipX,mipY);
    
//PATTERN RECOGNITION STARTED: 
  
  Int_t iNrec=FlagPhot(HoughResponse(),pCluLst,pTrk);                                                      //flag photons according to individual theta ckov with respect to most probable
  
  pTrk->SetHMPIDmip(mipX,mipY,mipQ,iNrec);                                                    //store mip info 

  if(iNrec<1){
    pTrk->SetHMPIDsignal(kNoPhotAccept);                                                      //no photon candidates are accepted
    return;
  }
  
  Double_t thetaC = FindRingCkov(pCluLst->GetEntries());                                    //find the best reconstructed theta Cherenkov
//    FindRingGeom(thetaC,2);
  pTrk->SetHMPIDsignal(thetaC);                                                             //store theta Cherenkov
  pTrk->SetHMPIDchi2(fCkovSigma2);                                                          //store errors squared
  
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
      Lors2Trs(dirCkov,thetaCer,phiCer);                       //find ckov (in TRS:the effective Cherenkov angle!)
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
void AliHMPIDRecon::Lors2Trs(TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer)const
{
  //Theta Cerenkov reconstruction 
  // Arguments: dirCkov photon vector in LORS
  //   Returns: thetaCer of photon in TRS
  //              phiCer of photon in TRS
//  TVector3 dirTrk;
//  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi());
//  Double_t thetaCer = TMath::ACos(dirCkov*dirTrk);
  TRotation mtheta;   mtheta.RotateY(-fTrkDir.Theta());
  TRotation mphi;       mphi.RotateZ(-fTrkDir.Phi());
  TRotation mrot=mtheta*mphi;
  TVector3 dirCkovTRS;
  dirCkovTRS=mrot*dirCkov;
  phiCer  = dirCkovTRS.Phi();                                          //actual value of the phi of the photon
  thetaCer= dirCkovTRS.Theta();                                        //actual value of thetaCerenkov of the photon
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::Trs2Lors(TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer)const
{
  //Theta Cerenkov reconstruction 
  // Arguments: dirCkov photon vector in TRS
  //   Returns: thetaCer of photon in LORS
  //              phiCer of photon in LORS
  TRotation mtheta;   mtheta.RotateY(fTrkDir.Theta());
  TRotation mphi;       mphi.RotateZ(fTrkDir.Phi());
  TRotation mrot=mphi*mtheta;
  TVector3 dirCkovLORS;
  dirCkovLORS=mrot*dirCkov;
  phiCer  = dirCkovLORS.Phi();                                          //actual value of the phi of the photon
  thetaCer= dirCkovLORS.Theta();                                        //actual value of thetaCerenkov of the photon
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecon::FindRingGeom(Double_t ckovAng,Int_t level)
{
// Find area covered in the PC acceptance
// Arguments: ckovAng - cerenkov angle
//            level   - precision in finding area and portion of ring accepted (multiple of 50)
//   Returns: area of the ring in cm^2 for given theta ckov
   
  Int_t kN=50*level;
  Int_t nPoints = 0;
  Double_t area=0;
  
  Bool_t first=kFALSE;
  TVector2 pos1;
  
  for(Int_t i=0;i<kN;i++){
   if(!first) {
      pos1=TracePhot(ckovAng,Double_t(TMath::TwoPi()*(i+1)/kN));                                    //find a good trace for the first photon
     if(pos1.X()==-999) continue;                                                                   //no area: open ring 	          
     if(!fParam->IsInside(pos1.X(),pos1.Y(),0)) {
       pos1 = IntWithEdge(fMipPos,pos1);                                                            // find the very first intersection...
     } else {
       if(!AliHMPIDParam::IsInDead(pos1.X(),pos1.Y())) nPoints++;                                   //photon is accepted if not in dead zone
     }
     first=kTRUE;
     continue;
   }
   TVector2 pos2=TracePhot(ckovAng,Double_t(TMath::TwoPi()*(i+1)/kN));                              //trace the next photon
   if(pos2.X()==-999) continue;                                                                     //no area: open ring 	     
   if(!fParam->IsInside(pos2.X(),pos2.Y(),0)) {
     pos2 = IntWithEdge(fMipPos,pos2);
   } else {
     if(!AliHMPIDParam::IsInDead(pos2.X(),pos2.Y())) nPoints++;                                     //photon is accepted if not in dead zone
   }
   area+=TMath::Abs((pos1-fMipPos).X()*(pos2-fMipPos).Y()-(pos1-fMipPos).Y()*(pos2-fMipPos).X());   //add area of the triangle... 	     
   pos1 = pos2;
  }
//---  find area and length of the ring;
  fRingAcc = (Double_t)nPoints/(Double_t)kN;
  area*=0.5;
  fRingArea = area;
}//FindRingGeom()
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
      
      sigma2 += 1./fParam->Sigma2(fTrkDir.Theta(),fTrkDir.Phi(),fPhotCkov[i],fPhotPhi[i]);
    }
  }//candidates loop
  
  if(sigma2>0) fCkovSigma2=1./sigma2;
  else         fCkovSigma2=1e10;  
  
  if(wei != 0.) weightThetaCerenkov /= wei; else weightThetaCerenkov = 0.;
  return weightThetaCerenkov;
}//FindCkovRing()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRecon::FlagPhot(Double_t ckov,TClonesArray *pCluLst, AliESDtrack *pTrk)
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
    if(fPhotCkov[i] >= tmin && fPhotCkov[i] <= tmax) { 
      fPhotFlag[i]=2;
      AddObjectToFriends(pCluLst,i,pTrk);
      iInsideCnt++;
    }
  }
      
  return iInsideCnt;
  
}//FlagPhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  AliHMPIDRecon::AddObjectToFriends(TClonesArray *pCluLst, Int_t photonIndex, AliESDtrack *pTrk)
{
// Add AliHMPIDcluster object to ESD friends
    
  AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(fPhotClusIndex[photonIndex]);	  
  AliHMPIDCluster *pClus = new AliHMPIDCluster(*pClu);
  pClus->SetChi2(fPhotCkov[photonIndex]);  
  pTrk->AddCalibObject(pClus);   
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliHMPIDRecon::TracePhot(Double_t ckovThe,Double_t ckovPhi)const
{
// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
//   Returns: distance between photon point on PC and track projection  
  
  Double_t theta,phi;
  TVector3  dirTRS,dirLORS;   
  dirTRS.SetMagThetaPhi(1,ckovThe,ckovPhi);                     //photon in TRS
  Trs2Lors(dirTRS,theta,phi);
  dirLORS.SetMagThetaPhi(1,theta,phi);                          //photon in LORS
  return TraceForward(dirLORS);                                 //now foward tracing
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
      FindRingGeom(lowerlimit);
      Double_t areaLow  = GetRingArea();
      FindRingGeom(upperlimit);
      Double_t areaHigh = GetRingArea();
      Double_t diffArea = areaHigh - areaLow;
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
    if((Double_t)((i+0.5)*fDTheta)>0.7) continue;
    resultw->Fill((Double_t)((i+0.5)*fDTheta),sumPhotsw);
  } 
// evaluate the "BEST" theta ckov as the maximum value of histogramm
  Double_t *pVec = resultw->GetArray();
  Int_t locMax = TMath::LocMax(nBin,pVec);
  delete phots;delete photsw;delete resultw; // Reset and delete objects
  
  return (Double_t)(locMax*fDTheta+0.5*fDTheta); //final most probable track theta ckov   
}//HoughResponse()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t AliHMPIDRecon::FindRingExt(Double_t ckov,Int_t ch,Double_t xPc,Double_t yPc,Double_t thRa,Double_t phRa)
{
// To find the acceptance of the ring even from external inputs. 
//    
//       
  Double_t xRa = xPc - (fParam->RadThick()+fParam->WinThick()+fParam->GapThick())*TMath::Cos(phRa)*TMath::Tan(thRa); //just linear extrapolation back to RAD
  Double_t yRa = yPc - (fParam->RadThick()+fParam->WinThick()+fParam->GapThick())*TMath::Sin(phRa)*TMath::Tan(thRa);
  
  Int_t nStep = 500;
  Int_t nPhi = 0;  

  Int_t ipc,ipadx,ipady;
    
  if(ckov>0){
    SetTrack(xRa,yRa,thRa,phRa);
    for(Int_t j=0;j<nStep;j++){
      TVector2 pos; pos=TracePhot(ckov,j*TMath::TwoPi()/(Double_t)(nStep-1));
      if(fParam->IsInDead(pos.X(),pos.Y())) continue;
      fParam->Lors2Pad(pos.X(),pos.Y(),ipc,ipadx,ipady);
      ipadx+=(ipc%2)*fParam->kPadPcX;
      ipady+=(ipc/2)*fParam->kPadPcY;
      if(fParam->IsDeadPad(ipadx,ipady,ch)) continue;
      nPhi++;
    }//point loop
  return ((Double_t)nPhi/(Double_t)nStep); 
  }//if
  return -1;
}  
