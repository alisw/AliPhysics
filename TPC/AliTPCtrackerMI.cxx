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


//-------------------------------------------------------
//          Implementation of the TPC tracker
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//  AliTPC parallel tracker - 
//  How to use?  - 
//  run AliTPCFindClusters.C macro - clusters neccessary for tracker are founded
//  run AliTPCFindTracksMI.C macro - to find tracks
//  tracks are written to AliTPCtracks.root file
//  for comparison also seeds are written to the same file - to special branch
//-------------------------------------------------------


/* $Id$ */



#include "Riostream.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>

#include "AliComplexCluster.h"
#include "AliESD.h"
#include "AliHelix.h"
#include "AliRunLoader.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "AliTPCclusterMI.h"
#include "AliTPCpolyTrack.h"
#include "AliTPCreco.h" 
#include "AliTPCtrackerMI.h"
#include "TStopwatch.h"
//

ClassImp(AliTPCseed)
ClassImp(AliTPCtrackerMI)


class AliTPCFastMath {
public:
  AliTPCFastMath();  
  static Double_t FastAsin(Double_t x);   
 private: 
  static Double_t fgFastAsin[20000];  //lookup table for fast asin computation
};

Double_t AliTPCFastMath::fgFastAsin[20000];
AliTPCFastMath gAliTPCFastMath; // needed to fill the LUT

AliTPCFastMath::AliTPCFastMath(){
  //
  // initialized lookup table;
  for (Int_t i=0;i<10000;i++){
    fgFastAsin[2*i] = TMath::ASin(i/10000.);
    fgFastAsin[2*i+1] = (TMath::ASin((i+1)/10000.)-fgFastAsin[2*i]);
  }
}

Double_t AliTPCFastMath::FastAsin(Double_t x){
  //
  // return asin using lookup table
  if (x>0){
    Int_t index = int(x*10000);
    return fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1];
  }
  x*=-1;
  Int_t index = int(x*10000);
  return -(fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1]);
}




Int_t AliTPCtrackerMI::UpdateTrack(AliTPCseed * track, Int_t accept){
  //
  //update track information using current cluster - track->fCurrentCluster


  AliTPCclusterMI* c =track->fCurrentCluster;
  if (accept>0) track->fCurrentClusterIndex1 |=0x8000;  //sign not accepted clusters

  UInt_t i = track->fCurrentClusterIndex1;

  Int_t sec=(i&0xff000000)>>24; 
  //Int_t row = (i&0x00ff0000)>>16; 
  track->fRow=(i&0x00ff0000)>>16;
  track->fSector = sec;
  //  Int_t index = i&0xFFFF;
  if (sec>=fParam->GetNInnerSector()) track->fRow += fParam->GetNRowLow(); 
  track->SetClusterIndex2(track->fRow, i);  
  //track->fFirstPoint = row;
  //if ( track->fLastPoint<row) track->fLastPoint =row;
  //  if (track->fRow<0 || track->fRow>160) {
  //  printf("problem\n");
  //}
  if (track->fFirstPoint>track->fRow) 
    track->fFirstPoint = track->fRow;
  if (track->fLastPoint<track->fRow) 
    track->fLastPoint  = track->fRow;
  

  track->fClusterPointer[track->fRow] = c;  
  //

  Float_t angle2 = track->GetSnp()*track->GetSnp();
  angle2 = TMath::Sqrt(angle2/(1-angle2)); 
  //
  //SET NEW Track Point
  //
  //  if (debug)
  {
    AliTPCTrackerPoint   &point =*(track->GetTrackPoint(track->fRow));
    //
    point.SetSigmaY(c->GetSigmaY2()/track->fCurrentSigmaY2);
    point.SetSigmaZ(c->GetSigmaZ2()/track->fCurrentSigmaZ2);
    point.SetErrY(sqrt(track->fErrorY2));
    point.SetErrZ(sqrt(track->fErrorZ2));
    //
    point.SetX(track->GetX());
    point.SetY(track->GetY());
    point.SetZ(track->GetZ());
    point.SetAngleY(angle2);
    point.SetAngleZ(track->GetTgl());
    if (point.fIsShared){
      track->fErrorY2 *= 4;
      track->fErrorZ2 *= 4;
    }
  }  

  Double_t chi2 = track->GetPredictedChi2(track->fCurrentCluster);
  //
  track->fErrorY2 *= 1.3;
  track->fErrorY2 += 0.01;    
  track->fErrorZ2 *= 1.3;   
  track->fErrorZ2 += 0.005;      
    //}
  if (accept>0) return 0;
  if (track->GetNumberOfClusters()%20==0){
    //    if (track->fHelixIn){
    //  TClonesArray & larr = *(track->fHelixIn);    
    //  Int_t ihelix = larr.GetEntriesFast();
    //  new(larr[ihelix]) AliHelix(*track) ;    
    //}
  }
  track->fNoCluster =0;
  return track->Update(c,chi2,i);
}



Int_t AliTPCtrackerMI::AcceptCluster(AliTPCseed * seed, AliTPCclusterMI * cluster, Float_t factor, 
                                      Float_t cory, Float_t corz)
{
  //
  // decide according desired precision to accept given 
  // cluster for tracking
  Double_t sy2=ErrY2(seed,cluster)*cory;
  Double_t sz2=ErrZ2(seed,cluster)*corz;
  //sy2=ErrY2(seed,cluster)*cory;
  //sz2=ErrZ2(seed,cluster)*cory;
  
  Double_t sdistancey2 = sy2+seed->GetSigmaY2();
  Double_t sdistancez2 = sz2+seed->GetSigmaZ2();
  
  Double_t rdistancey2 = (seed->fCurrentCluster->GetY()-seed->GetY())*
    (seed->fCurrentCluster->GetY()-seed->GetY())/sdistancey2;
  Double_t rdistancez2 = (seed->fCurrentCluster->GetZ()-seed->GetZ())*
    (seed->fCurrentCluster->GetZ()-seed->GetZ())/sdistancez2;
  
  Double_t rdistance2  = rdistancey2+rdistancez2;
  //Int_t  accept =0;
  
  if (rdistance2>16) return 3;
  
  
  if ((rdistancey2>9.*factor || rdistancez2>9.*factor) && cluster->GetType()==0)  
    return 2;  //suspisiouce - will be changed
  
  if ((rdistancey2>6.25*factor || rdistancez2>6.25*factor) && cluster->GetType()>0)  
    // strict cut on overlaped cluster
    return  2;  //suspisiouce - will be changed
  
  if ( (rdistancey2>1.*factor || rdistancez2>6.25*factor ) 
       && cluster->GetType()<0){
    seed->fNFoundable--;
    return 2;    
  }
  return 0;
}




//_____________________________________________________________________________
AliTPCtrackerMI::AliTPCtrackerMI(const AliTPCParam *par): 
AliTracker(), fkNIS(par->GetNInnerSector()/2), fkNOS(par->GetNOuterSector()/2)
{
  //---------------------------------------------------------------------
  // The main TPC tracker constructor
  //---------------------------------------------------------------------
  fInnerSec=new AliTPCSector[fkNIS];         
  fOuterSec=new AliTPCSector[fkNOS];
 
  Int_t i;
  for (i=0; i<fkNIS; i++) fInnerSec[i].Setup(par,0);
  for (i=0; i<fkNOS; i++) fOuterSec[i].Setup(par,1);

  fN=0;  fSectors=0;

  fSeeds=0;
  fNtracks = 0;
  fParam = par;  
  Int_t nrowlow = par->GetNRowLow();
  Int_t nrowup = par->GetNRowUp();

  
  for (Int_t i=0;i<nrowlow;i++){
    fXRow[i]     = par->GetPadRowRadiiLow(i);
    fPadLength[i]= par->GetPadPitchLength(0,i);
    fYMax[i]     = fXRow[i]*TMath::Tan(0.5*par->GetInnerAngle());
  }

  
  for (Int_t i=0;i<nrowup;i++){
    fXRow[i+nrowlow]      = par->GetPadRowRadiiUp(i);
    fPadLength[i+nrowlow] = par->GetPadPitchLength(60,i);
    fYMax[i+nrowlow]      = fXRow[i+nrowlow]*TMath::Tan(0.5*par->GetOuterAngle());
  }
  fSeeds=0;
  //
  fInput    = 0;
  fOutput   = 0;
  fSeedTree = 0;
  fTreeDebug =0;
  fNewIO     =0;
  fDebug     =0;
  fEvent     =0;
}

//_____________________________________________________________________________
AliTPCtrackerMI::~AliTPCtrackerMI() {
  //------------------------------------------------------------------
  // TPC tracker destructor
  //------------------------------------------------------------------
  delete[] fInnerSec;
  delete[] fOuterSec;
  if (fSeeds) {
    fSeeds->Delete(); 
    delete fSeeds;
  }
}

void AliTPCtrackerMI::SetIO()
{
  //
  fNewIO   =  kTRUE;
  fInput   =  AliRunLoader::GetTreeR("TPC", kFALSE,AliConfig::GetDefaultEventFolderName());
  
  fOutput  =  AliRunLoader::GetTreeT("TPC", kTRUE,AliConfig::GetDefaultEventFolderName());
  if (fOutput){
    AliTPCtrack *iotrack= new AliTPCtrack;
    fOutput->Branch("tracks","AliTPCtrack",&iotrack,32000,100);
    delete iotrack;
  }
}


void AliTPCtrackerMI::SetIO(TTree * input, TTree * output, AliESD * event)
{

  // set input
  fNewIO = kFALSE;
  fInput    = 0;
  fOutput   = 0;
  fSeedTree = 0;
  fTreeDebug =0;
  fInput = input;
  if (input==0){
    return;
  }  
  //set output
  fOutput = output;
  if (output){
    AliTPCtrack *iotrack= new AliTPCtrack;
    //    iotrack->fHelixIn   = new TClonesArray("AliHelix");
    //iotrack->fHelixOut  = new TClonesArray("AliHelix");    
    fOutput->Branch("tracks","AliTPCtrack",&iotrack,32000,100);
    delete iotrack;
  }
  if (output && (fDebug&2)){
    //write the full seed information if specified in debug mode
    //
    fSeedTree =  new TTree("Seeds","Seeds");
    AliTPCseed * vseed = new AliTPCseed;
    //
    TClonesArray * arrtr = new TClonesArray("AliTPCTrackPoint",160);
    arrtr->ExpandCreateFast(160);
    TClonesArray * arre = new TClonesArray("AliTPCExactPoint",160);
    //
    vseed->fPoints = arrtr;
    vseed->fEPoints = arre;
    //    vseed->fClusterPoints = arrcl;
    fSeedTree->Branch("seeds","AliTPCseed",&vseed,32000,99);
    delete arrtr;
    delete arre;    
    fTreeDebug = new TTree("trackDebug","trackDebug");
    TClonesArray * arrd = new TClonesArray("AliTPCTrackPoint2",0);
    fTreeDebug->Branch("debug",&arrd,32000,99);
  }


  //set ESD event  
  fEvent  = event;  
}

void AliTPCtrackerMI::FillESD(TObjArray* arr)
{
  //
  //
  //fill esds using updated tracks
  if (fEvent){
    // write tracks to the event
    // store index of the track
    Int_t nseed=arr->GetEntriesFast();
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
      if (!pt) continue; 
      pt->PropagateTo(fParam->GetInnerRadiusLow());
      if ( (pt->GetNumberOfClusters()>70)&& (Float_t(pt->GetNumberOfClusters())/Float_t(pt->fNFoundable))>0.55) {
	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);	
	//iotrack.SetTPCindex(i);
	fEvent->AddTrack(&iotrack);
      }        
    }
  }
}

void AliTPCtrackerMI::WriteTracks(TTree * tree)
{
  //
  // write tracks from seed array to selected tree
  //
  fOutput  = tree;
  if (fOutput){
    AliTPCtrack *iotrack= new AliTPCtrack;
    fOutput->Branch("tracks","AliTPCtrack",&iotrack,32000,100);
  }
  WriteTracks();
}

void AliTPCtrackerMI::WriteTracks()
{
  //
  // write tracks to the given output tree -
  // output specified with SetIO routine
  if (!fSeeds)  return;
  if (!fOutput){
    SetIO();
  }

  if (fOutput){
    AliTPCtrack *iotrack= 0;
    Int_t nseed=fSeeds->GetEntriesFast();
    //for (Int_t i=0; i<nseed; i++) {
    //  iotrack= (AliTPCtrack*)fSeeds->UncheckedAt(i);
    //  if (iotrack) break;      
    //}    
    //TBranch * br = fOutput->Branch("tracks","AliTPCtrack",&iotrack,32000,100);
    TBranch * br = fOutput->GetBranch("tracks");
    br->SetAddress(&iotrack);
    //
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i);    
      if (!pt) continue;    
      AliTPCtrack * track = new AliTPCtrack(*pt);
      iotrack = track;
      pt->fLab2 =i; 
      //      br->SetAddress(&iotrack);
      fOutput->Fill();
      delete track;
      iotrack =0;
    }
    //fOutput->GetDirectory()->cd();
    //fOutput->Write();
  }
  // delete iotrack;
  //
  if (fSeedTree){
    //write the full seed information if specified in debug mode
      
    AliTPCseed * vseed = new AliTPCseed;
    //
    TClonesArray * arrtr = new TClonesArray("AliTPCTrackPoint",160);
    arrtr->ExpandCreateFast(160);
    //TClonesArray * arrcl = new TClonesArray("AliTPCclusterMI",160);
    //arrcl->ExpandCreateFast(160);
    TClonesArray * arre = new TClonesArray("AliTPCExactPoint",160);
    //
    vseed->fPoints = arrtr;
    vseed->fEPoints = arre;
    //    vseed->fClusterPoints = arrcl;
    //TBranch * brseed = seedtree->Branch("seeds","AliTPCseed",&vseed,32000,99);
    TBranch * brseed = fSeedTree->GetBranch("seeds");
    
    Int_t nseed=fSeeds->GetEntriesFast();
    
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i);    
      if (!pt) continue;     
      pt->fPoints = arrtr;
      //      pt->fClusterPoints = arrcl;
      pt->fEPoints       = arre;
      pt->RebuildSeed();
      vseed = pt;
      brseed->SetAddress(&vseed);
      fSeedTree->Fill();
      pt->fPoints  = 0;
      pt->fEPoints = 0;
      //      pt->fClusterPoints = 0;
    }
    fSeedTree->Write();
    if (fTreeDebug) fTreeDebug->Write();
  }

}
  



Double_t AliTPCtrackerMI::ErrY2(AliTPCseed* seed, AliTPCclusterMI * cl){
  //
  //
  //seed->SetErrorY2(0.1);
  //return 0.1;
  //calculate look-up table at the beginning
  static Bool_t  ginit = kFALSE;
  static Float_t gnoise1,gnoise2,gnoise3;
  static Float_t ggg1[10000];
  static Float_t ggg2[10000];
  static Float_t ggg3[10000];
  static Float_t glandau1[10000];
  static Float_t glandau2[10000];
  static Float_t glandau3[10000];
  //
  static Float_t gcor01[500];
  static Float_t gcor02[500];
  static Float_t gcorp[500];
  //

  //
  if (ginit==kFALSE){
    for (Int_t i=1;i<500;i++){
      Float_t rsigma = float(i)/100.;
      gcor02[i] = TMath::Max(0.78 +TMath::Exp(7.4*(rsigma-1.2)),0.6);
      gcor01[i] = TMath::Max(0.72 +TMath::Exp(3.36*(rsigma-1.2)),0.6);
      gcorp[i]  = TMath::Max(TMath::Power((rsigma+0.5),1.5),1.2);
    }

    //
    for (Int_t i=3;i<10000;i++){
      //
      //
      // inner sector
      Float_t amp = float(i);
      Float_t padlength =0.75;
      gnoise1 = 0.0004/padlength;
      Float_t nel     = 0.268*amp;
      Float_t nprim   = 0.155*amp;
      ggg1[i]          = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.001*nel/(padlength*padlength))/nel;
      glandau1[i]      = (2.+0.12*nprim)*0.5* (2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau1[i]>1) glandau1[i]=1;
      glandau1[i]*=padlength*padlength/12.;      
      //
      // outer short
      padlength =1.;
      gnoise2   = 0.0004/padlength;
      nel       = 0.3*amp;
      nprim     = 0.133*amp;
      ggg2[i]      = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
      glandau2[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau2[i]>1) glandau2[i]=1;
      glandau2[i]*=padlength*padlength/12.;
      //
      //
      // outer long
      padlength =1.5;
      gnoise3   = 0.0004/padlength;
      nel       = 0.3*amp;
      nprim     = 0.133*amp;
      ggg3[i]      = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
      glandau3[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau3[i]>1) glandau3[i]=1;
      glandau3[i]*=padlength*padlength/12.;
      //
    }
    ginit = kTRUE;
  }
  //
  //
  //
  Int_t amp = int(TMath::Abs(cl->GetQ()));  
  if (amp>9999) {
    seed->SetErrorY2(1.);
    return 1.;
  }
  Float_t snoise2;
  Float_t z = TMath::Abs(fParam->GetZLength()-TMath::Abs(seed->GetZ()));
  Int_t ctype = cl->GetType();  
  Float_t padlength= GetPadPitchLength(seed->fRow);
  Float_t angle2 = seed->GetSnp()*seed->GetSnp();
  angle2 = angle2/(1-angle2); 
  //
  //cluster "quality"
  Int_t rsigmay = int(100.*cl->GetSigmaY2()/(seed->fCurrentSigmaY2));
  Float_t res;
  //
  if (fSectors==fInnerSec){
    snoise2 = gnoise1;
    res     = ggg1[amp]*z+glandau1[amp]*angle2;     
    if (ctype==0) res *= gcor01[rsigmay];
    if ((ctype>0)){
      res+=0.002;
      res*= gcorp[rsigmay];
    }
  }
  else {
    if (padlength<1.1){
      snoise2 = gnoise2;
      res     = ggg2[amp]*z+glandau2[amp]*angle2; 
      if (ctype==0) res *= gcor02[rsigmay];      
      if ((ctype>0)){
	res+=0.002;
	res*= gcorp[rsigmay];
      }
    }
    else{
      snoise2 = gnoise3;      
      res     = ggg3[amp]*z+glandau3[amp]*angle2; 
      if (ctype==0) res *= gcor02[rsigmay];
      if ((ctype>0)){
	res+=0.002;
	res*= gcorp[rsigmay];
      }
    }
  }  

  if (ctype<0){
    res+=0.005;
    res*=2.4;  // overestimate error 2 times
  }
  res+= snoise2;
 
  if (res<2*snoise2)
    res = 2*snoise2;
  
  seed->SetErrorY2(res);
  return res;


}



Double_t AliTPCtrackerMI::ErrZ2(AliTPCseed* seed, AliTPCclusterMI * cl){
  //
  //
  //seed->SetErrorY2(0.1);
  //return 0.1;
  //calculate look-up table at the beginning
  static Bool_t  ginit = kFALSE;
  static Float_t gnoise1,gnoise2,gnoise3;
  static Float_t ggg1[10000];
  static Float_t ggg2[10000];
  static Float_t ggg3[10000];
  static Float_t glandau1[10000];
  static Float_t glandau2[10000];
  static Float_t glandau3[10000];
  //
  static Float_t gcor01[1000];
  static Float_t gcor02[1000];
  static Float_t gcorp[1000];
  //

  //
  if (ginit==kFALSE){
    for (Int_t i=1;i<1000;i++){
      Float_t rsigma = float(i)/100.;
      gcor02[i] = TMath::Max(0.81 +TMath::Exp(6.8*(rsigma-1.2)),0.6);
      gcor01[i] = TMath::Max(0.72 +TMath::Exp(2.04*(rsigma-1.2)),0.6);
      gcorp[i]  = TMath::Max(TMath::Power((rsigma+0.5),1.5),1.2);
    }

    //
    for (Int_t i=3;i<10000;i++){
      //
      //
      // inner sector
      Float_t amp = float(i);
      Float_t padlength =0.75;
      gnoise1 = 0.0004/padlength;
      Float_t nel     = 0.268*amp;
      Float_t nprim   = 0.155*amp;
      ggg1[i]          = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.001*nel/(padlength*padlength))/nel;
      glandau1[i]      = (2.+0.12*nprim)*0.5* (2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau1[i]>1) glandau1[i]=1;
      glandau1[i]*=padlength*padlength/12.;      
      //
      // outer short
      padlength =1.;
      gnoise2   = 0.0004/padlength;
      nel       = 0.3*amp;
      nprim     = 0.133*amp;
      ggg2[i]      = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
      glandau2[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau2[i]>1) glandau2[i]=1;
      glandau2[i]*=padlength*padlength/12.;
      //
      //
      // outer long
      padlength =1.5;
      gnoise3   = 0.0004/padlength;
      nel       = 0.3*amp;
      nprim     = 0.133*amp;
      ggg3[i]      = fParam->GetDiffT()*fParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
      glandau3[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
      if (glandau3[i]>1) glandau3[i]=1;
      glandau3[i]*=padlength*padlength/12.;
      //
    }
    ginit = kTRUE;
  }
  //
  //
  //
  Int_t amp = int(TMath::Abs(cl->GetQ()));  
  if (amp>9999) {
    seed->SetErrorY2(1.);
    return 1.;
  }
  Float_t snoise2;
  Float_t z = TMath::Abs(fParam->GetZLength()-TMath::Abs(seed->GetZ()));
  Int_t ctype = cl->GetType();  
  Float_t padlength= GetPadPitchLength(seed->fRow);
  //
  Float_t angle2 = seed->GetSnp()*seed->GetSnp();
  //  if (angle2<0.6) angle2 = 0.6;
  angle2 = seed->GetTgl()*seed->GetTgl()*(1+angle2/(1-angle2)); 
  //
  //cluster "quality"
  Int_t rsigmaz = int(100.*cl->GetSigmaZ2()/(seed->fCurrentSigmaZ2));
  Float_t res;
  //
  if (fSectors==fInnerSec){
    snoise2 = gnoise1;
    res     = ggg1[amp]*z+glandau1[amp]*angle2;     
    if (ctype==0) res *= gcor01[rsigmaz];
    if ((ctype>0)){
      res+=0.002;
      res*= gcorp[rsigmaz];
    }
  }
  else {
    if (padlength<1.1){
      snoise2 = gnoise2;
      res     = ggg2[amp]*z+glandau2[amp]*angle2; 
      if (ctype==0) res *= gcor02[rsigmaz];      
      if ((ctype>0)){
	res+=0.002;
	res*= gcorp[rsigmaz];
      }
    }
    else{
      snoise2 = gnoise3;      
      res     = ggg3[amp]*z+glandau3[amp]*angle2; 
      if (ctype==0) res *= gcor02[rsigmaz];
      if ((ctype>0)){
	res+=0.002;
	res*= gcorp[rsigmaz];
      }
    }
  }  

  if (ctype<0){
    res+=0.002;
    res*=1.3;
  }
  if ((ctype<0) &&amp<70){
    res+=0.002;
    res*=1.3;  
  }
  res += snoise2;
  if (res<2*snoise2)
     res = 2*snoise2;
  if (res>3) res =3;
  seed->SetErrorZ2(res);
  return res;
}



/*
Double_t AliTPCtrackerMI::ErrZ2(AliTPCseed* seed, AliTPCclusterMI * cl){
  //
  //
  //seed->SetErrorZ2(0.1);
  //return 0.1;

  Float_t snoise2;
  Float_t z = TMath::Abs(fParam->GetZLength()-TMath::Abs(seed->GetZ()));
  //
  Float_t rsigmaz = cl->GetSigmaZ2()/(seed->fCurrentSigmaZ2);
  Int_t ctype = cl->GetType();
  Float_t amp = TMath::Abs(cl->GetQ());
  
  Float_t nel;
  Float_t nprim;
  //
  Float_t landau=2 ;    //landau fluctuation part
  Float_t gg=2;         // gg fluctuation part
  Float_t padlength= GetPadPitchLength(seed->GetX());
 
  if (fSectors==fInnerSec){
    snoise2 = 0.0004/padlength;
    nel     = 0.268*amp;
    nprim   = 0.155*amp;
    gg      = (2+0.001*nel/(padlength*padlength))/nel;
    landau  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
    if (landau>1) landau=1;
  }
  else {
    snoise2 = 0.0004/padlength;
    nel     = 0.3*amp;
    nprim   = 0.133*amp;
    gg      = (2+0.0008*nel/(padlength*padlength))/nel;
    landau  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
    if (landau>1) landau=1;
  }
  Float_t sdiff = gg*fParam->GetDiffT()*fParam->GetDiffT()*z;

  //
  Float_t angle2 = seed->GetSnp()*seed->GetSnp();
  angle2 = TMath::Sqrt((1-angle2));
  if (angle2<0.6) angle2 = 0.6;
  //angle2 = 1;

  Float_t angle = seed->GetTgl()/angle2;
  Float_t angular = landau*angle*angle*padlength*padlength/12.;
  Float_t res = sdiff + angular;

  
  if ((ctype==0) && (fSectors ==fOuterSec))
    res *= 0.81 +TMath::Exp(6.8*(rsigmaz-1.2));

  if ((ctype==0) && (fSectors ==fInnerSec))
    res *= 0.72 +TMath::Exp(2.04*(rsigmaz-1.2));
  
  if ((ctype>0)){
    res+=0.005;
    res*= TMath::Power(rsigmaz+0.5,1.5);  //0.31+0.147*ctype;
  }
  if (ctype<0){
    res+=0.002;
    res*=1.3;
  }
  if ((ctype<0) &&amp<70){
    res+=0.002;
    res*=1.3;  
  }
  res += snoise2;
  if (res<2*snoise2)
     res = 2*snoise2;

  seed->SetErrorZ2(res);
  return res;
}
*/



void AliTPCseed::Reset(Bool_t all)
{
  //
  //
  SetNumberOfClusters(0);
  fNFoundable = 0;
  SetChi2(0);
  ResetCovariance();
  /*
  if (fTrackPoints){
    for (Int_t i=0;i<8;i++){
      delete [] fTrackPoints[i];
    }
    delete fTrackPoints;
    fTrackPoints =0;
  }
  */

  if (all){   
    for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
    for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;
  }

}


void AliTPCseed::Modify(Double_t factor)
{

  //------------------------------------------------------------------
  //This function makes a track forget its history :)  
  //------------------------------------------------------------------
  if (factor<=0) {
    ResetCovariance();
    return;
  }
  fC00*=factor;
  fC10*=0;  fC11*=factor;
  fC20*=0;  fC21*=0;  fC22*=factor;
  fC30*=0;  fC31*=0;  fC32*=0;  fC33*=factor;
  fC40*=0;  fC41*=0;  fC42*=0;  fC43*=0;  fC44*=factor;
  SetNumberOfClusters(0);
  fNFoundable =0;
  SetChi2(0);
  fRemoval = 0;
  fCurrentSigmaY2 = 0.000005;
  fCurrentSigmaZ2 = 0.000005;
  fNoCluster     = 0;
  //fFirstPoint = 160;
  //fLastPoint  = 0;
}




Int_t  AliTPCseed::GetProlongation(Double_t xk, Double_t &y, Double_t & z) const
{
  //-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=xk.
  // doesn't change internal state of the track
  //-----------------------------------------------------------------
  
  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1;

  if (TMath::Abs(fP4*xk - fP2) >= 0.999) {   
    return 0;
  }

  //  Double_t y1=fP0, z1=fP1;
  Double_t c1=fP4*x1 - fP2, r1=sqrt(1.- c1*c1);
  Double_t c2=fP4*x2 - fP2, r2=sqrt(1.- c2*c2);
  
  y = fP0;
  z = fP1;
  //y += dx*(c1+c2)/(r1+r2);
  //z += dx*(c1+c2)/(c1*r2 + c2*r1)*fP3;
  
  Double_t dy = dx*(c1+c2)/(r1+r2);
  Double_t dz = 0;
  //
  Double_t delta = fP4*dx*(c1+c2)/(c1*r2 + c2*r1);
  /*
  if (TMath::Abs(delta)>0.0001){
    dz = fP3*TMath::ASin(delta)/fP4;
  }else{
    dz = dx*fP3*(c1+c2)/(c1*r2 + c2*r1);
  }
  */
  dz =  fP3*AliTPCFastMath::FastAsin(delta)/fP4;
  //
  y+=dy;
  z+=dz;
  

  return 1;  
}


//_____________________________________________________________________________
Double_t AliTPCseed::GetPredictedChi2(const AliTPCclusterMI *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  //Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters();
    if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}


//_________________________________________________________________________________________


Int_t AliTPCseed::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the sector - for given sector according z
  //-----------------------------------------------------------------
  AliTPCseed *t=(AliTPCseed*)o;

  if (fSort == 0){
    if (t->fRelativeSector>fRelativeSector) return -1;
    if (t->fRelativeSector<fRelativeSector) return 1;
    Double_t z2 = t->GetZ();
    Double_t z1 = GetZ();
    if (z2>z1) return 1;
    if (z2<z1) return -1;
    return 0;
  }
  else {
    Float_t f2 =1;
    f2 = 1-20*TMath::Sqrt(t->fC44)/(TMath::Abs(t->GetC())+0.0066);
    if (t->fBConstrain) f2=1.2;

    Float_t f1 =1;
    f1 = 1-20*TMath::Sqrt(fC44)/(TMath::Abs(GetC())+0.0066);

    if (fBConstrain)   f1=1.2;
 
    if (t->GetNumberOfClusters()*f2 <GetNumberOfClusters()*f1) return -1;
    else return +1;
  }
}

void AliTPCtrackerMI::RotateToLocal(AliTPCseed *seed)
{
  //rotate to track "local coordinata
  Float_t x = seed->GetX();
  Float_t y = seed->GetY();
  Float_t ymax = x*TMath::Tan(0.5*fSectors->GetAlpha());
  
  if (y > ymax) {
    seed->fRelativeSector= (seed->fRelativeSector+1) % fN;
    if (!seed->Rotate(fSectors->GetAlpha())) 
      return;
  } else if (y <-ymax) {
    seed->fRelativeSector= (seed->fRelativeSector-1+fN) % fN;
    if (!seed->Rotate(-fSectors->GetAlpha())) 
      return;
  }   

}




//_____________________________________________________________________________
Int_t AliTPCseed::Update(const AliTPCclusterMI *c, Double_t chisq, UInt_t /*index*/) {
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;

  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  Double_t cur=fP4 + k40*dy + k41*dz, eta=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.9) {
    return 0;
  }

  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = eta;
  fP3 += k30*dy + k31*dz;
  fP4  = cur;

  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k40*c03+k41*c13; 

  fC44-=k40*c04+k41*c14; 

  Int_t n=GetNumberOfClusters();
  //  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);

  return 1;
}



//_____________________________________________________________________________
Double_t AliTPCtrackerMI::F1old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  if ( xr*xr+yr*yr<=0.00000000000001) return 100;
  return -xr*yr/sqrt(xr*xr+yr*yr); 
}



//_____________________________________________________________________________
Double_t AliTPCtrackerMI::F1(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (det==0) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t AliTPCtrackerMI::F2(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (det==0) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  x0+=x1;
  x0*=c2;  
  return x0;
}



//_____________________________________________________________________________
Double_t AliTPCtrackerMI::F2old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
Double_t AliTPCtrackerMI::F3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


Double_t AliTPCtrackerMI::F3n(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2, Double_t c) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------

  //  Double_t angle1;
  
  //angle1    =  (z1-z2)*c/(TMath::ASin(c*x1-ni)-TMath::ASin(c*x2-ni));
  //
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  //  Double_t   angle2    =  TMath::ASin(d*c*0.5);
  //  Double_t   angle2    =  AliTPCFastMath::FastAsin(d*c*0.5);
  Double_t   angle2    = (d*c*0.5>0.1)? TMath::ASin(d*c*0.5): AliTPCFastMath::FastAsin(d*c*0.5);

  angle2  = (z1-z2)*c/(angle2*2.);
  return angle2;
}

Bool_t   AliTPCtrackerMI::GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z)
{//-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=x2.
  //-----------------------------------------------------------------
  
  Double_t dx=x2-x1;

  if (TMath::Abs(x[4]*x1 - x[2]) >= 0.999) {   
    return kFALSE;
  }

  Double_t c1=x[4]*x1 - x[2], r1=sqrt(1.- c1*c1);
  Double_t c2=x[4]*x2 - x[2], r2=sqrt(1.- c2*c2);  
  y = x[0];
  z = x[1];
  
  Double_t dy = dx*(c1+c2)/(r1+r2);
  Double_t dz = 0;
  //
  Double_t delta = x[4]*dx*(c1+c2)/(c1*r2 + c2*r1);
  
  if (TMath::Abs(delta)>0.01){
    dz = x[3]*TMath::ASin(delta)/x[4];
  }else{
    dz = x[3]*AliTPCFastMath::FastAsin(delta)/x[4];
  }
  
  //dz = x[3]*AliTPCFastMath::FastAsin(delta)/x[4];

  y+=dy;
  z+=dz;
  
  return kTRUE;  
}

Int_t  AliTPCtrackerMI::LoadClusters (TTree *tree)
{
  //
  //
  fInput = tree;
  return LoadClusters();
}

Int_t  AliTPCtrackerMI::LoadClusters()
{
  //
  // load clusters to the memory
  AliTPCClustersRow *clrow= new AliTPCClustersRow;
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(0);
  clrow->GetArray()->ExpandCreateFast(10000);
  //
  //  TTree * tree = fClustersArray.GetTree();

  TTree * tree = fInput;
  TBranch * br = tree->GetBranch("Segment");
  br->SetAddress(&clrow);
  //
  Int_t j=Int_t(tree->GetEntries());
  for (Int_t i=0; i<j; i++) {
    br->GetEntry(i);
    //  
    Int_t sec,row;
    fParam->AdjustSectorRow(clrow->GetID(),sec,row);
    //
    AliTPCRow * tpcrow=0;
    Int_t left=0;
    if (sec<fkNIS*2){
      tpcrow = &(fInnerSec[sec%fkNIS][row]);    
      left = sec/fkNIS;
    }
    else{
      tpcrow = &(fOuterSec[(sec-fkNIS*2)%fkNOS][row]);
      left = (sec-fkNIS*2)/fkNOS;
    }
    if (left ==0){
      tpcrow->fN1 = clrow->GetArray()->GetEntriesFast();
      tpcrow->fClusters1 = new AliTPCclusterMI[tpcrow->fN1];
      for (Int_t i=0;i<tpcrow->fN1;i++) 
	tpcrow->fClusters1[i] = *(AliTPCclusterMI*)(clrow->GetArray()->At(i));
    }
    if (left ==1){
      tpcrow->fN2 = clrow->GetArray()->GetEntriesFast();
      tpcrow->fClusters2 = new AliTPCclusterMI[tpcrow->fN2];
      for (Int_t i=0;i<tpcrow->fN2;i++) 
	tpcrow->fClusters2[i] = *(AliTPCclusterMI*)(clrow->GetArray()->At(i));
    }
  }
  //
  delete clrow;
  LoadOuterSectors();
  LoadInnerSectors();
  return 0;
}


void AliTPCtrackerMI::UnloadClusters()
{
  //
  // unload clusters from the memory
  //
  Int_t nrows = fOuterSec->GetNRows();
  for (Int_t sec = 0;sec<fkNOS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCRow*  tpcrow = &(fOuterSec[sec%fkNOS][row]);
      //      if (tpcrow){
      //	if (tpcrow->fClusters1) delete []tpcrow->fClusters1; 
      //	if (tpcrow->fClusters2) delete []tpcrow->fClusters2; 
      //}
      tpcrow->ResetClusters();
    }
  //
  nrows = fInnerSec->GetNRows();
  for (Int_t sec = 0;sec<fkNIS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCRow*  tpcrow = &(fInnerSec[sec%fkNIS][row]);
      //if (tpcrow){
      //	if (tpcrow->fClusters1) delete []tpcrow->fClusters1; 
      //if (tpcrow->fClusters2) delete []tpcrow->fClusters2; 
      //}
      tpcrow->ResetClusters();
    }

  return ;
}


//_____________________________________________________________________________
Int_t AliTPCtrackerMI::LoadOuterSectors() {
  //-----------------------------------------------------------------
  // This function fills outer TPC sectors with clusters.
  //-----------------------------------------------------------------
  Int_t nrows = fOuterSec->GetNRows();
  UInt_t index=0;
  for (Int_t sec = 0;sec<fkNOS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCRow*  tpcrow = &(fOuterSec[sec%fkNOS][row]);  
      Int_t sec2 = sec+2*fkNIS;
      //left
      Int_t ncl = tpcrow->fN1;
      while (ncl--) {
	AliTPCclusterMI *c= &(tpcrow->fClusters1[ncl]);
	index=(((sec2<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //right
      ncl = tpcrow->fN2;
      while (ncl--) {
	AliTPCclusterMI *c= &(tpcrow->fClusters2[ncl]);
	index=((((sec2+fkNOS)<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //
      // write indexes for fast acces
      //
      for (Int_t i=0;i<510;i++)
	tpcrow->fFastCluster[i]=-1;
      for (Int_t i=0;i<tpcrow->GetN();i++){
        Int_t zi = Int_t((*tpcrow)[i]->GetZ()+255.);
	tpcrow->fFastCluster[zi]=i;  // write index
      }
      Int_t last = 0;
      for (Int_t i=0;i<510;i++){
	if (tpcrow->fFastCluster[i]<0)
	  tpcrow->fFastCluster[i] = last;
	else
	  last = tpcrow->fFastCluster[i];
      }
    }  
  fN=fkNOS;
  fSectors=fOuterSec;
  return 0;
}


//_____________________________________________________________________________
Int_t  AliTPCtrackerMI::LoadInnerSectors() {
  //-----------------------------------------------------------------
  // This function fills inner TPC sectors with clusters.
  //-----------------------------------------------------------------
  Int_t nrows = fInnerSec->GetNRows();
  UInt_t index=0;
  for (Int_t sec = 0;sec<fkNIS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCRow*  tpcrow = &(fInnerSec[sec%fkNIS][row]);
      //
      //left
      Int_t ncl = tpcrow->fN1;
      while (ncl--) {
	AliTPCclusterMI *c= &(tpcrow->fClusters1[ncl]);
	index=(((sec<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //right
      ncl = tpcrow->fN2;
      while (ncl--) {
	AliTPCclusterMI *c= &(tpcrow->fClusters2[ncl]);
	index=((((sec+fkNIS)<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //
      // write indexes for fast acces
      //
      for (Int_t i=0;i<510;i++)
	tpcrow->fFastCluster[i]=-1;
      for (Int_t i=0;i<tpcrow->GetN();i++){
        Int_t zi = Int_t((*tpcrow)[i]->GetZ()+255.);
	tpcrow->fFastCluster[zi]=i;  // write index
      }
      Int_t last = 0;
      for (Int_t i=0;i<510;i++){
	if (tpcrow->fFastCluster[i]<0)
	  tpcrow->fFastCluster[i] = last;
	else
	  last = tpcrow->fFastCluster[i];
      }

    }  
   
  fN=fkNIS;
  fSectors=fInnerSec;
  return 0;
}



//_________________________________________________________________________
AliTPCclusterMI *AliTPCtrackerMI::GetClusterMI(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t sec=(index&0xff000000)>>24; 
  Int_t row=(index&0x00ff0000)>>16; 
  Int_t ncl=(index&0x00007fff)>>00;

  const AliTPCRow * tpcrow=0;
  AliTPCclusterMI * clrow =0;
  if (sec<fkNIS*2){
    tpcrow = &(fInnerSec[sec%fkNIS][row]);
    if (sec<fkNIS) 
      clrow = tpcrow->fClusters1;
    else
      clrow = tpcrow->fClusters2;
  }
  else{
    tpcrow = &(fOuterSec[(sec-fkNIS*2)%fkNOS][row]);
    if (sec-2*fkNIS<fkNOS)
      clrow = tpcrow->fClusters1;
    else
      clrow = tpcrow->fClusters2;
  }
  if (tpcrow==0) return 0;
  if (tpcrow->GetN()<=ncl) return 0;
  //  return (AliTPCclusterMI*)(*tpcrow)[ncl];      
  return &(clrow[ncl]);      
  
}



Int_t AliTPCtrackerMI::FollowToNext(AliTPCseed& t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  //
  Double_t  x= GetXrow(nr), ymax=GetMaxY(nr);
  AliTPCclusterMI *cl=0;
  Int_t tpcindex= t.GetClusterIndex2(nr);
  //
  // update current shape info every 5 pad-row
  //  if ( (nr%5==0) || t.GetNumberOfClusters()<2 || (t.fCurrentSigmaY2<0.0001) ){
    GetShape(&t,nr);    
    //}
  //  
  if (fIteration>0 && tpcindex>=-1){  //if we have already clusters 
    //        
    if (tpcindex==-1) return 0; //track in dead zone
    if (tpcindex>0){     //
      cl = t.fClusterPointer[nr];
      if ( (cl==0) ) cl = GetClusterMI(tpcindex);
      t.fCurrentClusterIndex1 = tpcindex; 
    }
    if (cl){      
      Int_t relativesector = ((tpcindex&0xff000000)>>24)%18;  // if previously accepted cluster in different sector
      Float_t angle = relativesector*fSectors->GetAlpha()+fSectors->GetAlphaShift();
      //
      if (angle<-TMath::Pi()) angle += 2*TMath::Pi();
      if (angle>=TMath::Pi()) angle -= 2*TMath::Pi();
      
      if (TMath::Abs(angle-t.GetAlpha())>0.001){
	Double_t rotation = angle-t.GetAlpha();
	t.fRelativeSector= relativesector;
	t.Rotate(rotation); 	
      }
      t.PropagateTo(x);
      //
      t.fCurrentCluster = cl; 
      t.fRow = nr;
      Int_t accept = AcceptCluster(&t,t.fCurrentCluster,1.);
      if ((tpcindex&0x8000)==0) accept =0;
      if (accept<3) { 
	//if founded cluster is acceptible
	if (cl->IsUsed(11)) {  // id cluster is shared inrease uncertainty
	  t.fErrorY2 += 0.03;
	  t.fErrorZ2 += 0.03; 
	  t.fErrorY2 *= 3;
	  t.fErrorZ2 *= 3; 
	}
	t.fNFoundable++;
	UpdateTrack(&t,accept);
	return 1;
      }    
    }
  }
  if (fIteration>1) return 0;  // not look for new cluster during refitting
  //
  UInt_t index=0;
  if (TMath::Abs(t.GetSnp())>0.95 || TMath::Abs(x*t.GetC()-t.GetEta())>0.95) return 0;
  Double_t  y=t.GetYat(x);
  if (TMath::Abs(y)>ymax){
    if (y > ymax) {
      t.fRelativeSector= (t.fRelativeSector+1) % fN;
      if (!t.Rotate(fSectors->GetAlpha())) 
	return 0;
    } else if (y <-ymax) {
      t.fRelativeSector= (t.fRelativeSector-1+fN) % fN;
      if (!t.Rotate(-fSectors->GetAlpha())) 
	return 0;
    }
    //return 1;
  }
  //
  if (!t.PropagateTo(x)) {
    if (fIteration==0) t.fRemoval = 10;
    return 0;
  }
  y=t.GetY(); 
  Double_t z=t.GetZ();
  //
  const AliTPCRow &krow=GetRow(t.fRelativeSector,nr);
  if ( (t.GetSigmaY2()<0) || t.GetSigmaZ2()<0) return 0;
  Double_t  roady  =1.;
  Double_t  roadz = 1.;
  //
  if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
    t.fInDead = kTRUE;
    t.SetClusterIndex2(nr,-1); 
    return 0;
  } 
  else
    {
      if (TMath::Abs(z)<(1.05*x+10)) t.fNFoundable++;
      else
	return 0;
    }   
  //calculate 
  if (krow) {
    //    cl = krow.FindNearest2(y+10.,z,roady,roadz,index);    
    cl = krow.FindNearest2(y,z,roady,roadz,index);    
    if (cl) t.fCurrentClusterIndex1 = krow.GetIndex(index);       
  }  
  if (cl) {
    t.fCurrentCluster = cl; 
    t.fRow = nr;
    if (fIteration==2&&cl->IsUsed(10)) return 0; 
    Int_t accept = AcceptCluster(&t,t.fCurrentCluster,1.);
    if (fIteration==2&&cl->IsUsed(11)) {
      t.fErrorY2 += 0.03;
      t.fErrorZ2 += 0.03; 
      t.fErrorY2 *= 3;
      t.fErrorZ2 *= 3; 
    }
    /*    
    if (t.fCurrentCluster->IsUsed(10)){
      //
      //     

      t.fNShared++;
      if (t.fNShared>0.7*t.GetNumberOfClusters()) {
	t.fRemoval =10;
	return 0;
      }
    }
    */
    if (accept<3) UpdateTrack(&t,accept);

  } else {  
    if ( fIteration==0 && t.fNFoundable*0.5 > t.GetNumberOfClusters()) t.fRemoval=10;
    
  }
  return 1;
}

Int_t AliTPCtrackerMI::FollowToNextFast(AliTPCseed& t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  //
  Double_t  x= GetXrow(nr), ymax=GetMaxY(nr);
  Double_t y,z; 
  if (!t.GetProlongation(x,y,z)) {
    t.fRemoval = 10;
    return 0;
  }
  //
  //
  if (TMath::Abs(y)>ymax){
    
    if (y > ymax) {
      t.fRelativeSector= (t.fRelativeSector+1) % fN;
      if (!t.Rotate(fSectors->GetAlpha())) 
	return 0;
    } else if (y <-ymax) {
      t.fRelativeSector= (t.fRelativeSector-1+fN) % fN;
      if (!t.Rotate(-fSectors->GetAlpha())) 
	return 0;
    }
    if (!t.PropagateTo(x)) {
      return 0;
    } 
    t.GetProlongation(x,y,z);
  }
  //
  // update current shape info every 3 pad-row
  if ( (nr%6==0) || t.GetNumberOfClusters()<2 || (t.fCurrentSigmaY2<0.0001) ){
    //    t.fCurrentSigmaY = GetSigmaY(&t);
    //t.fCurrentSigmaZ = GetSigmaZ(&t);
    GetShape(&t,nr);
  }
  //  
  AliTPCclusterMI *cl=0;
  UInt_t index=0;
  
  
  //Int_t nr2 = nr;
  const AliTPCRow &krow=GetRow(t.fRelativeSector,nr);
  if ( (t.GetSigmaY2()<0) || t.GetSigmaZ2()<0) return 0;
  Double_t  roady  =1.;
  Double_t  roadz = 1.;
  //
  Int_t row = nr;
  if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
    t.fInDead = kTRUE;
    t.SetClusterIndex2(row,-1); 
    return 0;
  } 
  else
    {
      if (TMath::Abs(z)>(1.05*x+10)) t.SetClusterIndex2(row,-1);
    }   
  //calculate 
  
  if ((cl==0)&&(krow)) {
    //    cl = krow.FindNearest2(y+10,z,roady,roadz,index);    
    cl = krow.FindNearest2(y,z,roady,roadz,index);    

    if (cl) t.fCurrentClusterIndex1 = krow.GetIndex(index);       
  }  

  if (cl) {
    t.fCurrentCluster = cl; 
    //    Int_t accept = AcceptCluster(&t,t.fCurrentCluster,1.);        
    //if (accept<3){
      t.SetClusterIndex2(row,index);
      t.fClusterPointer[row] = cl;
      //}
  }
  return 1;
}



Int_t AliTPCtrackerMI::UpdateClusters(AliTPCseed& t,  Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  t.fCurrentCluster  = 0;
  t.fCurrentClusterIndex1 = 0;   
   
  Double_t xt=t.GetX();
  Int_t     row = GetRowNumber(xt)-1; 
  Double_t  ymax= GetMaxY(nr);

  if (row < nr) return 1; // don't prolongate if not information until now -
  if (TMath::Abs(t.GetSnp())>0.9 && t.GetNumberOfClusters()>40. && fIteration!=2) {
    t.fRemoval =10;
    return 0;  // not prolongate strongly inclined tracks
  } 
  if (TMath::Abs(t.GetSnp())>0.95) {
    t.fRemoval =10;
    return 0;  // not prolongate strongly inclined tracks
  }

  Double_t x= GetXrow(nr);
  Double_t y,z;
  //t.PropagateTo(x+0.02);
  //t.PropagateTo(x+0.01);
  if (!t.PropagateTo(x)){
    return 0;
  }
  //
  y=t.GetY();
  z=t.GetZ();

  if (TMath::Abs(y)>ymax){
    if (y > ymax) {
      t.fRelativeSector= (t.fRelativeSector+1) % fN;
      if (!t.Rotate(fSectors->GetAlpha())) 
	return 0;
    } else if (y <-ymax) {
      t.fRelativeSector= (t.fRelativeSector-1+fN) % fN;
      if (!t.Rotate(-fSectors->GetAlpha())) 
	return 0;
    }
    //    if (!t.PropagateTo(x)){
    //  return 0;
    //}
    return 1;
    //y = t.GetY();    
  }
  //

  AliTPCRow &krow=GetRow(t.fRelativeSector,nr);

  if (TMath::Abs(TMath::Abs(y)-ymax)<krow.fDeadZone){
    t.fInDead = kTRUE;
    t.SetClusterIndex2(nr,-1); 
    return 0;
  } 
  else
    {
      if (TMath::Abs(t.GetZ())<(1.05*t.GetX()+10)) t.fNFoundable++;
      else
	return 0;      
    }

  // update current
  if ( (nr%6==0) || t.GetNumberOfClusters()<2){
    //    t.fCurrentSigmaY = GetSigmaY(&t);
    //t.fCurrentSigmaZ = GetSigmaZ(&t);
    GetShape(&t,nr);
  }
    
  AliTPCclusterMI *cl=0;
  UInt_t index=0;
  //
  Double_t roady = 1.;
  Double_t roadz = 1.;
  //

  if (!cl){
    index = t.GetClusterIndex2(nr);    
    if ( (index>0) && (index&0x8000)==0){
      cl = t.fClusterPointer[nr];
      if ( (cl==0) && (index>0)) cl = GetClusterMI(index);
      t.fCurrentClusterIndex1 = index;
      if (cl) {
	t.fCurrentCluster  = cl;
	return 1;
      }
    }
  }

  if (krow) {    
    //cl = krow.FindNearest2(y+10,z,roady,roadz,index);      
    cl = krow.FindNearest2(y,z,roady,roadz,index);      
  }

  if (cl) t.fCurrentClusterIndex1 = krow.GetIndex(index);   
  t.fCurrentCluster  = cl;

  return 1;
}


Int_t AliTPCtrackerMI::FollowToNextCluster(AliTPCseed & t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------

  //update error according neighborhoud

  if (t.fCurrentCluster) {
    t.fRow = nr; 
    Int_t accept = AcceptCluster(&t,t.fCurrentCluster,1.);
    
    if (t.fCurrentCluster->IsUsed(10)){
      //
      //
      //  t.fErrorZ2*=2;
      //  t.fErrorY2*=2;
      t.fNShared++;
      if (t.fNShared>0.7*t.GetNumberOfClusters()) {
	t.fRemoval =10;
	return 0;
      }
    }   
    if (fIteration>0) accept = 0;
    if (accept<3)  UpdateTrack(&t,accept);  
 
  } else {
    if (fIteration==0){
      if ( ( (t.GetSigmaY2()+t.GetSigmaZ2())>0.16)&& t.GetNumberOfClusters()>18) t.fRemoval=10;      
      if (  t.GetChi2()/t.GetNumberOfClusters()>6 &&t.GetNumberOfClusters()>18) t.fRemoval=10;      

      if (( (t.fNFoundable*0.5 > t.GetNumberOfClusters()) || t.fNoCluster>15)) t.fRemoval=10;
    }
  }
  return 1;
}



//_____________________________________________________________________________
Int_t AliTPCtrackerMI::FollowProlongation(AliTPCseed& t, Int_t rf, Int_t step) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  //
  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  //
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN;
    
  Int_t first = GetRowNumber(xt)-1;
  for (Int_t nr= first; nr>=rf; nr-=step) {    
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;
    if (FollowToNext(t,nr)==0) 
      if (!t.IsActive()) 
	return 0;
    
  }   
  return 1;
}


//_____________________________________________________________________________
Int_t AliTPCtrackerMI::FollowProlongationFast(AliTPCseed& t, Int_t rf, Int_t step) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  //
  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN;
    
  for (Int_t nr=GetRowNumber(xt)-1; nr>=rf; nr-=step) {
    
    if (FollowToNextFast(t,nr)==0) 
      if (!t.IsActive()) return 0;
    
  }   
  return 1;
}





Int_t AliTPCtrackerMI::FollowBackProlongation(AliTPCseed& t, Int_t rf) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  //  Double_t xt=t.GetX();  
  //
  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN;
    
  Int_t first = t.fFirstPoint;
  //
  if (first<0) first=0;
  for (Int_t nr=first; nr<=rf; nr++) {
    //if ( (t.GetSnp()<0.9))
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;
    FollowToNext(t,nr);                                                             
  }   
  return 1;
}




   
Float_t AliTPCtrackerMI::OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t & sum2)
{
  //
  //
  sum1=0;
  sum2=0;
  Int_t sum=0;
  //
  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;  

  Float_t dy2 =TMath::Abs((s1->GetY() - s2->GetY()));
  dy2*=dy2;
  Float_t distance = TMath::Sqrt(dz2+dy2);
  if (distance>4.) return 0; // if there are far away  - not overlap - to reduce combinatorics
 
  //  Int_t offset =0;
  Int_t firstpoint = TMath::Min(s1->fFirstPoint,s2->fFirstPoint);
  Int_t lastpoint = TMath::Max(s1->fLastPoint,s2->fLastPoint);
  if (lastpoint>160) 
    lastpoint =160;
  if (firstpoint<0) 
    firstpoint = 0;
  if (firstpoint>lastpoint) {
    firstpoint =lastpoint;
    //    lastpoint  =160;
  }
    
  
  for (Int_t i=firstpoint-1;i<lastpoint+1;i++){
    if (s1->GetClusterIndex2(i)>0) sum1++;
    if (s2->GetClusterIndex2(i)>0) sum2++;
    if (s1->GetClusterIndex2(i)==s2->GetClusterIndex2(i) && s1->GetClusterIndex2(i)>0) {
      sum++;
    }
  }
  if (sum<5) return 0;

  Float_t summin = TMath::Min(sum1+1,sum2+1);
  Float_t ratio = (sum+1)/Float_t(summin);
  return ratio;
}

void  AliTPCtrackerMI::SignShared(AliTPCseed * s1, AliTPCseed * s2)
{
  //
  //
  if (TMath::Abs(s1->GetC()-s2->GetC())>0.004) return;
  if (TMath::Abs(s1->GetTgl()-s2->GetTgl())>0.6) return;

  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;
  Float_t dy2 =(s1->GetY() - s2->GetY());
  dy2*=dy2;
  Float_t distance = dz2+dy2;
  if (distance>325.) return ; // if there are far away  - not overlap - to reduce combinatorics
  
  //
  Int_t sumshared=0;
  //
  Int_t firstpoint = TMath::Max(s1->fFirstPoint,s2->fFirstPoint);
  Int_t lastpoint = TMath::Min(s1->fLastPoint,s2->fLastPoint);
  //
  if (firstpoint>=lastpoint-5) return;;

  for (Int_t i=firstpoint;i<lastpoint;i++){
    //    if ( (s1->GetClusterIndex2(i)&0xFFFF8FFF)==(s2->GetClusterIndex2(i)&0xFFFF8FFF) && s1->GetClusterIndex2(i)>0) {
    if ( (s1->GetClusterIndex2(i))==(s2->GetClusterIndex2(i)) && s1->GetClusterIndex2(i)>0) {
      sumshared++;
    }
  }
  if (sumshared>4){
    // sign clusters
    //
    for (Int_t i=firstpoint;i<lastpoint;i++){
      //      if ( (s1->GetClusterIndex2(i)&0xFFFF8FFF)==(s2->GetClusterIndex2(i)&0xFFFF8FFF) && s1->GetClusterIndex2(i)>0) {
      if ( (s1->GetClusterIndex2(i))==(s2->GetClusterIndex2(i)) && s1->GetClusterIndex2(i)>0) {
	AliTPCTrackerPoint *p1  = s1->GetTrackPoint(i);
	AliTPCTrackerPoint *p2  = s2->GetTrackPoint(i);; 
	if (s1->IsActive()&&s2->IsActive()){
	  p1->fIsShared = kTRUE;
	  p2->fIsShared = kTRUE;
	}	
      }
    }
  }
  //  
  if (sumshared>10){
    for (Int_t i=0;i<4;i++){
      if (s1->fOverlapLabels[3*i]==0){
	s1->fOverlapLabels[3*i] = s2->GetLabel();
	s1->fOverlapLabels[3*i+1] = sumshared;
	s1->fOverlapLabels[3*i+2] = s2->GetUniqueID();
	break;
      }	
    }
    for (Int_t i=0;i<4;i++){
      if (s2->fOverlapLabels[3*i]==0){
	s2->fOverlapLabels[3*i] = s1->GetLabel();
	s2->fOverlapLabels[3*i+1] = sumshared;
	s2->fOverlapLabels[3*i+2] = s1->GetUniqueID();
	break;
      }	
    }    
  }
  
}

void  AliTPCtrackerMI::SignShared(TObjArray * arr)
{
  //
  //sort trackss according sectors
  //  
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    //if (pt) RotateToLocal(pt);
    pt->fSort = 0;
  }
  arr->UnSort();
  arr->Sort();  // sorting according z
  arr->Expand(arr->GetEntries());
  //
  //
  Int_t nseed=arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    for (Int_t j=0;j<=12;j++){
      pt->fOverlapLabels[j] =0;
    }
  }
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->fRemoval>10) continue;
    for (Int_t j=i+1; j<nseed; j++){
      AliTPCseed *pt2=(AliTPCseed*)arr->UncheckedAt(j);
      //      if (pt2){
      if (pt2->fRemoval<=10) {
	if ( TMath::Abs(pt->fRelativeSector-pt2->fRelativeSector)>0) break;
	SignShared(pt,pt2);
      }
    }  
  }
}

void  AliTPCtrackerMI::RemoveDouble(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removalindex)
{
  //
  //sort trackss according sectors
  //
  if (fDebug&1) {
    Info("RemoveDouble","Number of tracks before double removal- %d\n",arr->GetEntries());
  }
  //
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    pt->fSort = 0;
  }
  arr->UnSort();
  arr->Sort();  // sorting according z
  arr->Expand(arr->GetEntries());
  //
  //reset overlap labels
  //
  Int_t nseed=arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    pt->SetUniqueID(i);
    for (Int_t j=0;j<=12;j++){
      pt->fOverlapLabels[j] =0;
    }
  }
  //
  //sign shared tracks
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->fRemoval>10) continue;
    Float_t deltac = pt->GetC()*0.1;
    for (Int_t j=i+1; j<nseed; j++){
      AliTPCseed *pt2=(AliTPCseed*)arr->UncheckedAt(j);
      //      if (pt2){
      if (pt2->fRemoval<=10) {
	if ( TMath::Abs(pt->fRelativeSector-pt2->fRelativeSector)>0) break;
	if (TMath::Abs(pt->GetC()  -pt2->GetC())>deltac) continue;
	if (TMath::Abs(pt->GetTgl()-pt2->GetTgl())>0.05) continue;
	//
	SignShared(pt,pt2);
      }
    }
  }
  //
  // remove highly shared tracks
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->fRemoval>10) continue;
    //
    Int_t sumshared =0;
    for (Int_t j=0;j<4;j++){
      sumshared = pt->fOverlapLabels[j*3+1];      
    }
    Float_t factor = factor1;
    if (pt->fRemoval>0) factor = factor2;
    if (sumshared/pt->GetNumberOfClusters()>factor){
      for (Int_t j=0;j<4;j++){
	if (pt->fOverlapLabels[3*j]==0) continue;
	if (pt->fOverlapLabels[3*j+1]<5) continue; 
	if (pt->fRemoval==removalindex) continue;      
	AliTPCseed * pt2 = (AliTPCseed*)arr->UncheckedAt(pt->fOverlapLabels[3*j+2]);
	if (!pt2) continue;
	if (pt2->GetSigma2C()<pt->GetSigma2C()){
	  //	  pt->fRemoval = removalindex;
	  delete arr->RemoveAt(i);	  
	  break;
	}
      }      
    }
  }
  arr->Compress();
  if (fDebug&1) {
    Info("RemoveDouble","Number of tracks after double removal- %d\n",arr->GetEntries());
  }
}






void AliTPCtrackerMI::SortTracks(TObjArray * arr, Int_t mode) const
{
  //
  //sort tracks in array according mode criteria
  Int_t nseed = arr->GetEntriesFast();    
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    pt->fSort = mode;
  }
  arr->UnSort();
  arr->Sort();
}

void AliTPCtrackerMI::RemoveUsed(TObjArray * arr, Float_t factor1,  Float_t factor2, Int_t removalindex)
{

  //Loop over all tracks and remove "overlaps"
  //
  //
  Int_t nseed = arr->GetEntriesFast();  
  Int_t good =0;

  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      delete arr->RemoveAt(i);
    }
    else{
      pt->fSort =1;
      pt->fBSigned = kFALSE;
    }
  }
  arr->Compress();
  nseed = arr->GetEntriesFast();
  arr->UnSort();
  arr->Sort();
  //
  //unsign used
  UnsignClusters();
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    Int_t found,foundable,shared;
    if (pt->IsActive()) 
      pt->GetClusterStatistic(0,160,found, foundable,shared,kFALSE);
    else
      pt->GetClusterStatistic(0,160,found, foundable,shared,kTRUE); 
    //
    Double_t factor = factor2;
    if (pt->fBConstrain) factor = factor1;

    if ((Float_t(shared)/Float_t(found))>factor){
      pt->Desactivate(removalindex);
      continue;
    }

    good++;
    for (Int_t i=0; i<160; i++) {
      Int_t index=pt->GetClusterIndex2(i);
      if (index<0 || index&0x8000 ) continue;
      AliTPCclusterMI *c= pt->fClusterPointer[i];        
      if (!c) continue;
      //      if (!c->IsUsed(10)) c->Use(10);
      //if (pt->IsActive()) 
      c->Use(10);  
      //else
      //	c->Use(5);
    }
    
  }
  fNtracks = good;
  if (fDebug>0){
    Info("RemoveUsed","\n*****\nNumber of good tracks after shared removal\t%d\n",fNtracks);
  }
}

void AliTPCtrackerMI::UnsignClusters() 
{
  //
  // loop over all clusters and unsign them
  //
  
  for (Int_t sec=0;sec<fkNIS;sec++){
    for (Int_t row=0;row<fInnerSec->GetNRows();row++){
      AliTPCclusterMI *cl = fInnerSec[sec][row].fClusters1;
      for (Int_t icl =0;icl< fInnerSec[sec][row].fN1;icl++)
	//	if (cl[icl].IsUsed(10)) 	
	cl[icl].Use(-1);
      cl = fInnerSec[sec][row].fClusters2;
      for (Int_t icl =0;icl< fInnerSec[sec][row].fN2;icl++)
	//if (cl[icl].IsUsed(10)) 	
	  cl[icl].Use(-1);      
    }
  }
  
  for (Int_t sec=0;sec<fkNOS;sec++){
    for (Int_t row=0;row<fOuterSec->GetNRows();row++){
      AliTPCclusterMI *cl = fOuterSec[sec][row].fClusters1;
      for (Int_t icl =0;icl< fOuterSec[sec][row].fN1;icl++)
	//if (cl[icl].IsUsed(10)) 	
	  cl[icl].Use(-1);
      cl = fOuterSec[sec][row].fClusters2;
      for (Int_t icl =0;icl< fOuterSec[sec][row].fN2;icl++)
	//if (cl[icl].IsUsed(10)) 	
	cl[icl].Use(-1);      
    }
  }
  
}



void AliTPCtrackerMI::SignClusters(TObjArray * arr, Float_t fnumber, Float_t fdensity)
{
  //
  //sign clusters to be "used"
  //
  // snumber and sdensity sign number of sigmas - bellow mean value to be accepted
  // loop over "primaries"
  
  Float_t sumdens=0;
  Float_t sumdens2=0;
  Float_t sumn   =0;
  Float_t sumn2  =0;
  Float_t sumchi =0;
  Float_t sumchi2 =0;

  Float_t sum    =0;

  TStopwatch timer;
  timer.Start();

  Int_t nseed = arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    if (!(pt->IsActive())) continue;
    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->fNFoundable);
    if ( (dens>0.7) && (pt->GetNumberOfClusters()>70)){
      sumdens += dens;
      sumdens2+= dens*dens;
      sumn    += pt->GetNumberOfClusters();
      sumn2   += pt->GetNumberOfClusters()*pt->GetNumberOfClusters();
      Float_t chi2 = pt->GetChi2()/pt->GetNumberOfClusters();
      if (chi2>5) chi2=5;
      sumchi  +=chi2;
      sumchi2 +=chi2*chi2;
      sum++;
    }
  }

  Float_t mdensity = 0.9;
  Float_t meann    = 130;
  Float_t meanchi  = 1;
  Float_t sdensity = 0.1;
  Float_t smeann    = 10;
  Float_t smeanchi  =0.4;
  

  if (sum>20){
    mdensity = sumdens/sum;
    meann    = sumn/sum;
    meanchi  = sumchi/sum;
    //
    sdensity = sumdens2/sum-mdensity*mdensity;
    sdensity = TMath::Sqrt(sdensity);
    //
    smeann   = sumn2/sum-meann*meann;
    smeann   = TMath::Sqrt(smeann);
    //
    smeanchi = sumchi2/sum - meanchi*meanchi;
    smeanchi = TMath::Sqrt(smeanchi);
  }


  //REMOVE  SHORT DELTAS or tracks going out of sensitive volume of TPC
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (pt->fBSigned) continue;
    if (pt->fBConstrain) continue;    
    //if (!(pt->IsActive())) continue;
    /*
    Int_t found,foundable,shared;    
    pt->GetClusterStatistic(0,160,found, foundable,shared);
    if (shared/float(found)>0.3) {
      if (shared/float(found)>0.9 ){
	//delete arr->RemoveAt(i);
      }
      continue;
    }
    */
    Bool_t isok =kFALSE;
    if ( (pt->fNShared/pt->GetNumberOfClusters()<0.5) &&pt->GetNumberOfClusters()>60)
      isok = kTRUE;
    if ((TMath::Abs(1/pt->GetC())<100.) && (pt->fNShared/pt->GetNumberOfClusters()<0.7))
      isok =kTRUE;
    if  (TMath::Abs(pt->GetZ()/pt->GetX())>1.1)
      isok =kTRUE;
    if ( (TMath::Abs(pt->GetSnp()>0.7) && pt->GetD(0,0)>60.))
      isok =kTRUE;
    
    if (isok)     
      for (Int_t i=0; i<160; i++) {	
	Int_t index=pt->GetClusterIndex2(i);
	if (index<0) continue;
	AliTPCclusterMI *c= pt->fClusterPointer[i];
	if (!c) continue;
	//if (!(c->IsUsed(10))) c->Use();  
	c->Use(10);  
      }
  }
  
  
  //
  Double_t maxchi  = meanchi+2.*smeanchi;

  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    //if (!(pt->IsActive())) continue;
    if (pt->fBSigned) continue;
    Double_t chi     = pt->GetChi2()/pt->GetNumberOfClusters();
    if (chi>maxchi) continue;

    Float_t bfactor=1;
    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->fNFoundable);
   
    //sign only tracks with enoug big density at the beginning
    
    if ((pt->GetDensityFirst(40)<0.75) && pt->GetNumberOfClusters()<meann) continue; 
    
    
    Double_t mindens = TMath::Max(double(mdensity-sdensity*fdensity*bfactor),0.65);
    Double_t minn    = TMath::Max(Int_t(meann-fnumber*smeann*bfactor),50);
   
    //    if (pt->fBConstrain) mindens = TMath::Max(mdensity-sdensity*fdensity*bfactor,0.65);
    if ( (pt->fRemoval==10) && (pt->GetSnp()>0.8)&&(dens>mindens))
      minn=0;

    if ((dens>mindens && pt->GetNumberOfClusters()>minn) && chi<maxchi ){
      //Int_t noc=pt->GetNumberOfClusters();
      pt->fBSigned = kTRUE;
      for (Int_t i=0; i<160; i++) {

	Int_t index=pt->GetClusterIndex2(i);
	if (index<0) continue;
	AliTPCclusterMI *c= pt->fClusterPointer[i];
	if (!c) continue;
	//	if (!(c->IsUsed(10))) c->Use();  
	c->Use(10);  
      }
    }
  }
  //  gLastCheck = nseed;
  //  arr->Compress();
  if (fDebug>0){
    timer.Print();
  }
}


void  AliTPCtrackerMI::StopNotActive(TObjArray * arr, Int_t row0, Float_t th0, Float_t th1, Float_t th2) const
{
  // stop not active tracks
  // take th1 as threshold for number of founded to number of foundable on last 10 active rows
  // take th2 as threshold for number of founded to number of foundable on last 20 active rows 
  Int_t nseed = arr->GetEntriesFast();  
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (!(pt->IsActive())) continue;
    StopNotActive(pt,row0,th0, th1,th2);
  }
}



void  AliTPCtrackerMI::StopNotActive(AliTPCseed * seed, Int_t row0, Float_t th0, Float_t th1,
 Float_t th2) const
{
  // stop not active tracks
  // take th1 as threshold for number of founded to number of foundable on last 10 active rows
  // take th2 as threshold for number of founded to number of foundable on last 20 active rows 
  Int_t sumgood1  = 0;
  Int_t sumgood2  = 0;
  Int_t foundable = 0;
  Int_t maxindex = seed->fLastPoint;  //last foundable row
  if (seed->fNFoundable*th0 > seed->GetNumberOfClusters()) {
    seed->Desactivate(10) ;
    return;
  }

  for (Int_t i=row0; i<maxindex; i++){
    Int_t index = seed->GetClusterIndex2(i);
    if (index!=-1) foundable++;
    //if (!c) continue;
    if (foundable<=30) sumgood1++;
    if (foundable<=50) {
      sumgood2++;
    }
    else{ 
      break;
    }        
  }
  if (foundable>=30.){ 
     if (sumgood1<(th1*30.)) seed->Desactivate(10);
  }
  if (foundable>=50)
    if (sumgood2<(th2*50.)) seed->Desactivate(10);
}


Int_t AliTPCtrackerMI::RefitInward(AliESD *event)
{
  //
  // back propagation of ESD tracks
  //
  //return 0;
  fEvent = event;
  ReadSeeds(event,2);
  fIteration=2;
  //PrepareForProlongation(fSeeds,1);
  PropagateForward2(fSeeds);
  Int_t ntracks=0;
  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed * seed = (AliTPCseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    seed->PropagateTo(fParam->GetInnerRadiusLow());
    AliESDtrack *esd=event->GetTrack(i);
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    if (seed->GetNumberOfClusters()>60){
      esd->UpdateTrackParams(seed,AliESDtrack::kTPCrefit); 
      ntracks++;
    }
    else{
      //printf("problem\n");
    }
  }
  Info("RefitInward","Number of refitted tracks %d",ntracks);
  fEvent =0;
  //WriteTracks();
  return 0;
}


Int_t AliTPCtrackerMI::PropagateBack(AliESD *event)
{
  //
  // back propagation of ESD tracks
  //

  fEvent = event;
  fIteration = 1;
  ReadSeeds(event,0);
  PropagateBack(fSeeds);
  Int_t nseed = fSeeds->GetEntriesFast();
  Int_t ntracks=0;
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed * seed = (AliTPCseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    AliESDtrack *esd=event->GetTrack(i);
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    if (seed->GetNumberOfClusters()>60){
      esd->UpdateTrackParams(seed,AliESDtrack::kTPCout);
      ntracks++;
    }
  }
  Info("PropagateBack","Number of back propagated tracks %d",ntracks);
  fEvent =0;
  //WriteTracks();
  return 0;
}


void AliTPCtrackerMI::DeleteSeeds()
{
  //
  //delete Seeds
  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed * seed = (AliTPCseed*)fSeeds->At(i);
    if (seed) delete fSeeds->RemoveAt(i);
  }
  delete fSeeds;
  fSeeds =0;
}

void AliTPCtrackerMI::ReadSeeds(AliESD *event, Int_t direction)
{
  //
  //read seeds from the event
  
  Int_t nentr=event->GetNumberOfTracks();
  if (fDebug>0){
    Info("ReadSeeds", "Number of ESD tracks: %d\n", nentr);
  }
  if (fSeeds) 
    DeleteSeeds();
  if (!fSeeds){   
    fSeeds = new TObjArray(nentr);
  }
  UnsignClusters();
  //  Int_t ntrk=0;  
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);
    ULong_t status=esd->GetStatus();    
    AliTPCtrack t(*esd);
    AliTPCseed *seed = new AliTPCseed(t,t.GetAlpha());
    if ((status==AliESDtrack::kTPCin)&&(direction==1)) seed->ResetCovariance(); 
    if ( direction ==2 &&(status & AliESDtrack::kTRDrefit) == 0 ) seed->ResetCovariance();
    if ( direction ==2 && ((status & AliESDtrack::kTPCout) == 0) ) {
      fSeeds->AddAt(0,i);
      delete seed;
      continue;    
    }
    if ( direction ==2 &&(status & AliESDtrack::kTRDrefit) > 0 )  {
      Double_t par0[5],par1[5],x;
      esd->GetInnerExternalParameters(x,par0);
      esd->GetExternalParameters(x,par1);
      Double_t delta1 = TMath::Abs(par0[4]-par1[4])/(0.000000001+TMath::Abs(par0[4]+par1[4]));
      Double_t delta2 = TMath::Abs(par0[3]-par1[3]);
      //reset covariance if suspicious 
      if ( (delta1>0.1) || (delta2>0.006))
	seed->ResetCovariance();
    }

    //
    //
    // rotate to the local coordinate system
   
    fSectors=fInnerSec; fN=fkNIS;
    
    Double_t alpha=seed->GetAlpha() - fSectors->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/fSectors->GetAlpha())%fN;
    alpha =ns*fSectors->GetAlpha() + fSectors->GetAlphaShift();
    if (alpha<-TMath::Pi()) alpha += 2*TMath::Pi();
    if (alpha>=TMath::Pi()) alpha -= 2*TMath::Pi();
    alpha-=seed->GetAlpha();  
    if (!seed->Rotate(alpha)) {
      delete seed;
      continue;
    }
    seed->fEsd = esd;
    //
    //seed->PropagateTo(fSectors->GetX(0));
    //
    //    Int_t index = esd->GetTPCindex();
    //AliTPCseed * seed2= (AliTPCseed*)fSeeds->At(index);
    //if (direction==2){
    //  AliTPCseed * seed2  = ReSeed(seed,0.,0.5,1.);
    //  if (seed2) {
    //	delete seed;
    //	seed = seed2;
    //  }
    //}
    //
    // sign clusters
    for (Int_t irow=0;irow<160;irow++){
      Int_t index = seed->GetClusterIndex2(irow);    
      if (index>0){ 
	//
	AliTPCclusterMI * cl = GetClusterMI(index);
	seed->fClusterPointer[irow] = cl;
	if (cl){
	  if ((index & 0x8000)==0){
	    cl->Use(10);  // accepted cluster	  
	  }else{
	    cl->Use(6);   // close cluster not accepted
	  }	
      	}else{
	   Info("ReadSeeds","Not found cluster");
	}
      }
    }
    fSeeds->AddAt(seed,i);
  }
}



//_____________________________________________________________________________
void AliTPCtrackerMI::MakeSeeds3(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				 Float_t deltay, Int_t ddsec) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  // SEEDING WITH VERTEX CONSTRAIN 
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut
  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  AliTPCseed * seed = new AliTPCseed;
  Double_t alpha=fSectors->GetAlpha(), shift=fSectors->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  //  Double_t x1 =fOuterSec->GetX(i1);
  //Double_t xx2=fOuterSec->GetX(i2);
  
  Double_t x1 =GetXrow(i1);
  Double_t xx2=GetXrow(i2);

  Double_t x3=GetX(), y3=GetY(), z3=GetZ();

  Int_t imiddle = (i2+i1)/2;    //middle pad row index
  Double_t xm = GetXrow(imiddle); // radius of middle pad-row
  const AliTPCRow& krm=GetRow(sec,imiddle); //middle pad -row
  //
  Int_t ns =sec;   

  const AliTPCRow& kr1=GetRow(ns,i1);
  Double_t ymax  = GetMaxY(i1)-kr1.fDeadZone-1.5;  
  Double_t ymaxm = GetMaxY(imiddle)-kr1.fDeadZone-1.5;  

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1-x3)*(x1-x3)+(ymax+5-y3)*(ymax+5-y3));
  if (dvertexmax*0.5*cuts[0]>0.85){
    cuts[0] = 0.85/(dvertexmax*0.5+1.);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut

  //  Int_t ddsec = 1;
  if (deltay>0) ddsec = 0; 
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;
    Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();    
    //if (TMath::Abs(y1)>ymax) continue;

    if (deltay>0 && TMath::Abs(ymax-TMath::Abs(y1))> deltay ) continue;  // seed only at the edge

    // find possible directions    
    Float_t anglez = (z1-z3)/(x1-x3); 
    Float_t extraz = z1 - anglez*(x1-xx2);  // extrapolated z      
    //
    //
    //find   rotation angles relative to line given by vertex and point 1
    Double_t dvertex2 = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
    Double_t dvertex  = TMath::Sqrt(dvertex2);
    Double_t angle13  = TMath::ATan((y1-y3)/(x1-x3));
    Double_t cs13     = cos(-angle13), sn13 = sin(-angle13);            
    
    //
    // loop over 2 sectors
    Int_t dsec1=-ddsec;
    Int_t dsec2= ddsec;
    if (y1<0)  dsec2= 0;
    if (y1>0)  dsec1= 0;
    
    Double_t dddz1=0;  // direction of delta inclination in z axis
    Double_t dddz2=0;
    if ( (z1-z3)>0)
      dddz1 =1;    
    else
      dddz2 =1;
    //
    for (Int_t dsec = dsec1; dsec<=dsec2;dsec++){
      Int_t sec2 = sec + dsec;
      // 
      //      AliTPCRow&  kr2  = fOuterSec[(sec2+fkNOS)%fkNOS][i2];
      //AliTPCRow&  kr2m = fOuterSec[(sec2+fkNOS)%fkNOS][imiddle];
      AliTPCRow&  kr2  = GetRow((sec2+fkNOS)%fkNOS,i2);
      AliTPCRow&  kr2m = GetRow((sec2+fkNOS)%fkNOS,imiddle);
      Int_t  index1 = TMath::Max(kr2.Find(extraz-0.6-dddz1*TMath::Abs(z1)*0.05)-1,0);
      Int_t  index2 = TMath::Min(kr2.Find(extraz+0.6+dddz2*TMath::Abs(z1)*0.05)+1,kr2);

      // rotation angles to p1-p3
      Double_t cs13r     = cos(-angle13+dsec*alpha)/dvertex, sn13r = sin(-angle13+dsec*alpha)/dvertex;            
      Double_t x2,   y2,   z2; 
      //
      //      Double_t dymax = maxangle*TMath::Abs(x1-xx2);

      //
      Double_t dxx0 =  (xx2-x3)*cs13r;
      Double_t dyy0 =  (xx2-x3)*sn13r;
      for (Int_t js=index1; js < index2; js++) {
	const AliTPCclusterMI *kcl = kr2[js];
	if (kcl->IsUsed(10)) continue;	
	//
	//calcutate parameters
	//	
	Double_t yy0 =  dyy0 +(kcl->GetY()-y3)*cs13r;
	// stright track
	if (TMath::Abs(yy0)<0.000001) continue;
	Double_t xx0 =  dxx0 -(kcl->GetY()-y3)*sn13r;
	Double_t y0  =  0.5*(xx0*xx0+yy0*yy0-xx0)/yy0;
	Double_t r02 = (0.25+y0*y0)*dvertex2;	
	//curvature (radius) cut
	if (r02<r2min) continue;		
       
	nin0++;	
	//
	Double_t c0  = 1/TMath::Sqrt(r02);
	if (yy0>0) c0*=-1.;	
	       
       
	//Double_t dfi0   = 2.*TMath::ASin(dvertex*c0*0.5);
	//Double_t dfi1   = 2.*TMath::ASin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);
	Double_t dfi0   = 2.*AliTPCFastMath::FastAsin(dvertex*c0*0.5);
	Double_t dfi1   = 2.*AliTPCFastMath::FastAsin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);  
	//
	//
	Double_t z0  =  kcl->GetZ();  
	Double_t zzzz2    = z1-(z1-z3)*dfi1/dfi0;
	if (TMath::Abs(zzzz2-z0)>0.5) continue;       
	nin1++;              
	//	
	Double_t dip    = (z1-z0)*c0/dfi1;        
	Double_t x0 = (0.5*cs13+y0*sn13)*dvertex*c0;
	//
	y2 = kcl->GetY(); 
	if (dsec==0){
	  x2 = xx2; 
	  z2 = kcl->GetZ();	  
	}
	else
	  {
	    // rotation	
	    z2 = kcl->GetZ();  
	    x2= xx2*cs-y2*sn*dsec;
	    y2=+xx2*sn*dsec+y2*cs;
	  }
	
	x[0] = y1;
	x[1] = z1;
	x[2] = x0;
	x[3] = dip;
	x[4] = c0;
	//
	//
	// do we have cluster at the middle ?
	Double_t ym,zm;
	GetProlongation(x1,xm,x,ym,zm);
	UInt_t dummy; 
	AliTPCclusterMI * cm=0;
	if (TMath::Abs(ym)-ymaxm<0){	  
	  cm = krm.FindNearest2(ym,zm,1.0,0.6,dummy);
	  if ((!cm) || (cm->IsUsed(10))) {	  
	    continue;
	  }
	}
	else{	  
	  // rotate y1 to system 0
	  // get state vector in rotated system 
	  Double_t yr1  = (-0.5*sn13+y0*cs13)*dvertex*c0;
	  Double_t xr2  =  x0*cs+yr1*sn*dsec;
	  Double_t xr[5]={kcl->GetY(),kcl->GetZ(), xr2, dip, c0};
	  //
	  GetProlongation(xx2,xm,xr,ym,zm);
	  if (TMath::Abs(ym)-ymaxm<0){
	    cm = kr2m.FindNearest2(ym,zm,1.0,0.6,dummy);
	    if ((!cm) || (cm->IsUsed(10))) {	  
	      continue;
	    }
	  }
	}
       

	Double_t dym = 0;
	Double_t dzm = 0;
	if (cm){
	  dym = ym - cm->GetY();
	  dzm = zm - cm->GetZ();
	}
	nin2++;


	//
	//
        Double_t sy1=kr1[is]->GetSigmaY2()*2., sz1=kr1[is]->GetSigmaZ2()*2.;
        Double_t sy2=kcl->GetSigmaY2()*2.,     sz2=kcl->GetSigmaZ2()*2.;
	//Double_t sy3=400*3./12., sy=0.1, sz=0.1;
	Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	//Double_t sy3=25000*x[4]*x[4]*60+0.5, sy=0.1, sz=0.1;

	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
	
	Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
	
        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
	
        UInt_t index=kr1.GetIndex(is);
	AliTPCseed *track=new(seed) AliTPCseed(index, x, c, x1, ns*alpha+shift);
	
	track->fIsSeeding = kTRUE;
	track->fSeed1 = i1;
	track->fSeed2 = i2;
	track->fSeedType=3;

       
	//if (dsec==0) {
	  FollowProlongation(*track, (i1+i2)/2,1);
	  Int_t foundable,found,shared;
	  track->GetClusterStatistic((i1+i2)/2,i1, found, foundable, shared, kTRUE);
	  if ((found<0.55*foundable)  || shared>0.5*found || (track->GetSigmaY2()+track->GetSigmaZ2())>0.5){
	    seed->Reset();
	    seed->~AliTPCseed();
	    continue;
	  }
	  //}
	
	nin++;
	FollowProlongation(*track, i2,1);
	
	
	//Int_t rc = 1;
	track->fBConstrain =1;
	//	track->fLastPoint = i1+fInnerSec->GetNRows();  // first cluster in track position
	track->fLastPoint = i1;  // first cluster in track position
	track->fFirstPoint = track->fLastPoint;
	
	if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->fNFoundable*0.6 || 
	    track->fNShared>0.4*track->GetNumberOfClusters() ) {
	  seed->Reset();
	  seed->~AliTPCseed();
	  continue;
	}
	nout1++;
        // Z VERTEX CONDITION
	Double_t zv;
        zv = track->GetZ()+track->GetTgl()/track->GetC()*
	  ( asin(-track->GetEta()) - asin(track->GetX()*track->GetC()-track->GetEta()));
	if (TMath::Abs(zv-z3)>cuts[2]) {
	  FollowProlongation(*track, TMath::Max(i2-20,0));
	  zv = track->GetZ()+track->GetTgl()/track->GetC()*
	    ( asin(-track->GetEta()) - asin(track->GetX()*track->GetC()-track->GetEta()));
	  if (TMath::Abs(zv-z3)>cuts[2]){
	    FollowProlongation(*track, TMath::Max(i2-40,0));
	    zv = track->GetZ()+track->GetTgl()/track->GetC()*
	      ( asin(-track->GetEta()) - asin(track->GetX()*track->GetC()-track->GetEta()));
	    if (TMath::Abs(zv-z3)>cuts[2] &&(track->GetNumberOfClusters() > track->fNFoundable*0.7)){
	      // make seed without constrain
	      AliTPCseed * track2 = MakeSeed(track,0.2,0.5,1.);
	      FollowProlongation(*track2, i2,1);
	      track2->fBConstrain = kFALSE;
	      track2->fSeedType = 1;
	      arr->AddLast(track2); 
	      seed->Reset();
	      seed->~AliTPCseed();
	      continue;		
	    }
	    else{
	      seed->Reset();
	      seed->~AliTPCseed();
	      continue;
	    
	    }
	  }
	}
	
	track->fSeedType =0;
	arr->AddLast(track); 
	seed = new AliTPCseed; 	
	nout2++;
	// don't consider other combinations
	if (track->GetNumberOfClusters() > track->fNFoundable*0.8)
	  break;
      }
    }
  }
  if (fDebug>3){
    Info("MakeSeeds3","\nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2);
  }
  delete seed;
}


void AliTPCtrackerMI::MakeSeeds5(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				 Float_t deltay) {
  


  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut


  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;
  Int_t nout3 =0;
  Double_t x[5], c[15];
  //
  // make temporary seed
  AliTPCseed * seed = new AliTPCseed;
  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  //

  // first 3 padrows
  Double_t x1 = GetXrow(i1-1);
  const    AliTPCRow& kr1=GetRow(sec,i1-1);
  Double_t y1max  = GetMaxY(i1-1)-kr1.fDeadZone-1.5;  
  //
  Double_t x1p = GetXrow(i1);
  const    AliTPCRow& kr1p=GetRow(sec,i1);
  //
  Double_t x1m = GetXrow(i1-2);
  const    AliTPCRow& kr1m=GetRow(sec,i1-2);

  //
  //last 3 padrow for seeding
  AliTPCRow&  kr3  = GetRow((sec+fkNOS)%fkNOS,i1-7);
  Double_t    x3   =  GetXrow(i1-7);
  //  Double_t    y3max= GetMaxY(i1-7)-kr3.fDeadZone-1.5;  
  //
  AliTPCRow&  kr3p  = GetRow((sec+fkNOS)%fkNOS,i1-6);
  Double_t    x3p   = GetXrow(i1-6);
  //
  AliTPCRow&  kr3m  = GetRow((sec+fkNOS)%fkNOS,i1-8);
  Double_t    x3m   = GetXrow(i1-8);

  //
  //
  // middle padrow
  Int_t im = i1-4;                           //middle pad row index
  Double_t xm         = GetXrow(im);         // radius of middle pad-row
  const AliTPCRow& krm=GetRow(sec,im);   //middle pad -row
  //  Double_t ymmax = GetMaxY(im)-kr1.fDeadZone-1.5;  
  //
  //
  Double_t deltax  = x1-x3;
  Double_t dymax   = deltax*cuts[1];
  Double_t dzmax   = deltax*cuts[3];
  //
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;
    Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();    
    //
    if (deltay>0 && TMath::Abs(y1max-TMath::Abs(y1))> deltay ) continue;  // seed only at the edge    
    // 
    Int_t  index1 = TMath::Max(kr3.Find(z1-dzmax)-1,0);
    Int_t  index2 = TMath::Min(kr3.Find(z1+dzmax)+1,kr3);
    //    
    Double_t y3,   z3;
    //
    //
    UInt_t index;
    for (Int_t js=index1; js < index2; js++) {
      const AliTPCclusterMI *kcl = kr3[js];
      if (kcl->IsUsed(10)) continue;
      y3 = kcl->GetY(); 
      // apply angular cuts
      if (TMath::Abs(y1-y3)>dymax) continue;
      x3 = x3; 
      z3 = kcl->GetZ();	
      if (TMath::Abs(z1-z3)>dzmax) continue;
      //
      Double_t angley = (y1-y3)/(x1-x3);
      Double_t anglez = (z1-z3)/(x1-x3);
      //
      Double_t erry = TMath::Abs(angley)*(x1-x1m)*0.5+0.5;
      Double_t errz = TMath::Abs(anglez)*(x1-x1m)*0.5+0.5;
      //
      Double_t yyym = angley*(xm-x1)+y1;
      Double_t zzzm = anglez*(xm-x1)+z1;

      const AliTPCclusterMI *kcm = krm.FindNearest2(yyym,zzzm,erry,errz,index);
      if (!kcm) continue;
      if (kcm->IsUsed(10)) continue;
      
      erry = TMath::Abs(angley)*(x1-x1m)*0.4+0.5;
      errz = TMath::Abs(anglez)*(x1-x1m)*0.4+0.5;
      //
      //
      //
      Int_t used  =0;
      Int_t found =0;
      //
      // look around first
      const AliTPCclusterMI *kc1m = kr1m.FindNearest2(angley*(x1m-x1)+y1,
						      anglez*(x1m-x1)+z1,
						      erry,errz,index);
      //
      if (kc1m){
	found++;
	if (kc1m->IsUsed(10)) used++;
      }
      const AliTPCclusterMI *kc1p = kr1p.FindNearest2(angley*(x1p-x1)+y1,
						      anglez*(x1p-x1)+z1,
						      erry,errz,index);
      //
      if (kc1p){
	found++;
	if (kc1p->IsUsed(10)) used++;
      }
      if (used>1)  continue;
      if (found<1) continue; 

      //
      // look around last
      const AliTPCclusterMI *kc3m = kr3m.FindNearest2(angley*(x3m-x3)+y3,
						      anglez*(x3m-x3)+z3,
						      erry,errz,index);
      //
      if (kc3m){
	found++;
	if (kc3m->IsUsed(10)) used++;
      }
      else 
	continue;
      const AliTPCclusterMI *kc3p = kr3p.FindNearest2(angley*(x3p-x3)+y3,
						      anglez*(x3p-x3)+z3,
						      erry,errz,index);
      //
      if (kc3p){
	found++;
	if (kc3p->IsUsed(10)) used++;
      }
      else 
	continue;
      if (used>1)  continue;
      if (found<3) continue;       
      //
      Double_t x2,y2,z2;
      x2 = xm;
      y2 = kcm->GetY();
      z2 = kcm->GetZ();
      //
                  	
      x[0]=y1;
      x[1]=z1;
      x[4]=F1(x1,y1,x2,y2,x3,y3);
      //if (TMath::Abs(x[4]) >= cuts[0]) continue;
      nin0++;
      //
      x[2]=F2(x1,y1,x2,y2,x3,y3);
      nin1++;
      //
      x[3]=F3n(x1,y1,x2,y2,z1,z2,x[4]);
      //if (TMath::Abs(x[3]) > cuts[3]) continue;
      nin2++;
      //
      //
      Double_t sy1=0.1,  sz1=0.1;
      Double_t sy2=0.1,  sz2=0.1;
      Double_t sy3=0.1,  sy=0.1, sz=0.1;
      
      Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
      Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
      Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
      Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
      Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
      Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
      
      Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
      Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
      Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
      Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
      
      c[0]=sy1;
      c[1]=0.;       c[2]=sz1; 
      c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
      c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
      c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
      c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
      c[13]=f30*sy1*f40+f32*sy2*f42;
      c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
      
      //	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
      
      UInt_t index=kr1.GetIndex(is);
      AliTPCseed *track=new(seed) AliTPCseed(index, x, c, x1, sec*alpha+shift);
      
      track->fIsSeeding = kTRUE;

      nin++;      
      FollowProlongation(*track, i1-7,1);
      if (track->GetNumberOfClusters() < track->fNFoundable*0.75 || 
	  track->fNShared>0.6*track->GetNumberOfClusters() || ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.6){
	seed->Reset();
	seed->~AliTPCseed();
	continue;
      }
      nout1++;
      nout2++;	
      //Int_t rc = 1;
      FollowProlongation(*track, i2,1);
      track->fBConstrain =0;
      track->fLastPoint = i1+fInnerSec->GetNRows();  // first cluster in track position
      track->fFirstPoint = track->fLastPoint;
      
      if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	  track->GetNumberOfClusters()<track->fNFoundable*0.7 || 
	  track->fNShared>2. || track->GetChi2()/track->GetNumberOfClusters()>6 || ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.5 ) {
	seed->Reset();
	seed->~AliTPCseed();
	continue;
      }
   
      {
	FollowProlongation(*track, TMath::Max(i2-10,0),1);
	AliTPCseed * track2 = MakeSeed(track,0.2,0.5,0.9);
	FollowProlongation(*track2, i2,1);
	track2->fBConstrain = kFALSE;
	track2->fSeedType = 4;
	arr->AddLast(track2); 
	seed->Reset();
	seed->~AliTPCseed();
      }
      
   
      //arr->AddLast(track); 
      //seed = new AliTPCseed; 	
      nout3++;
    }
  }
  
  if (fDebug>3){
    Info("MakeSeeds5","\nSeeding statiistic:\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2,nout3);
  }
  delete seed;
}


//_____________________________________________________________________________
void AliTPCtrackerMI::MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t */*cuts[4]*/,
				 Float_t deltay, Bool_t /*bconstrain*/) {
  //-----------------------------------------------------------------
  // This function creates track seeds - without vertex constraint
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut        - not applied
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut    - not applied 
  // cuts[3]   - fP3 cut
  Int_t nin0=0;
  Int_t nin1=0;
  Int_t nin2=0;
  Int_t nin3=0;
  //  Int_t nin4=0;
  //Int_t nin5=0;

  

  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  Int_t row0 = (i1+i2)/2;
  Int_t drow = (i1-i2)/2;
  const AliTPCRow& kr0=fSectors[sec][row0];
  AliTPCRow * kr=0;

  AliTPCpolyTrack polytrack;
  Int_t nclusters=fSectors[sec][row0];
  AliTPCseed * seed = new AliTPCseed;

  Int_t sumused=0;
  Int_t cused=0;
  Int_t cnused=0;
  for (Int_t is=0; is < nclusters; is++) {  //LOOP over clusters
    Int_t nfound =0;
    Int_t nfoundable =0;
    for (Int_t iter =1; iter<2; iter++){   //iterations
      const AliTPCRow& krm=fSectors[sec][row0-iter];
      const AliTPCRow& krp=fSectors[sec][row0+iter];      
      const AliTPCclusterMI * cl= kr0[is];
      
      if (cl->IsUsed(10)) {
	cused++;
      }
      else{
	cnused++;
      }
      Double_t x = kr0.GetX();
      // Initialization of the polytrack
      nfound =0;
      nfoundable =0;
      polytrack.Reset();
      //
      Double_t y0= cl->GetY();
      Double_t z0= cl->GetZ();
      Float_t erry = 0;
      Float_t errz = 0;
      
      Double_t ymax = fSectors->GetMaxY(row0)-kr0.fDeadZone-1.5;
      if (deltay>0 && TMath::Abs(ymax-TMath::Abs(y0))> deltay ) continue;  // seed only at the edge
      
      erry = (0.5)*cl->GetSigmaY2()/TMath::Sqrt(cl->GetQ())*6;	    
      errz = (0.5)*cl->GetSigmaZ2()/TMath::Sqrt(cl->GetQ())*6;      
      polytrack.AddPoint(x,y0,z0,erry, errz);

      sumused=0;
      if (cl->IsUsed(10)) sumused++;


      Float_t roady = (5*TMath::Sqrt(cl->GetSigmaY2()+0.2)+1.)*iter;
      Float_t roadz = (5*TMath::Sqrt(cl->GetSigmaZ2()+0.2)+1.)*iter;
      //
      x = krm.GetX();
      AliTPCclusterMI * cl1 = krm.FindNearest(y0,z0,roady,roadz);
      if (cl1 && TMath::Abs(ymax-TMath::Abs(y0))) {
	erry = (0.5)*cl1->GetSigmaY2()/TMath::Sqrt(cl1->GetQ())*3;	    
	errz = (0.5)*cl1->GetSigmaZ2()/TMath::Sqrt(cl1->GetQ())*3;
	if (cl1->IsUsed(10))  sumused++;
	polytrack.AddPoint(x,cl1->GetY(),cl1->GetZ(),erry,errz);
      }
      //
      x = krp.GetX();
      AliTPCclusterMI * cl2 = krp.FindNearest(y0,z0,roady,roadz);
      if (cl2) {
	erry = (0.5)*cl2->GetSigmaY2()/TMath::Sqrt(cl2->GetQ())*3;	    
	errz = (0.5)*cl2->GetSigmaZ2()/TMath::Sqrt(cl2->GetQ())*3;
	if (cl2->IsUsed(10)) sumused++;	 
	polytrack.AddPoint(x,cl2->GetY(),cl2->GetZ(),erry,errz);
      }
      //
      if (sumused>0) continue;
      nin0++;
      polytrack.UpdateParameters();
      // follow polytrack
      roadz = 1.2;
      roady = 1.2;
      //
      Double_t yn,zn;
      nfoundable = polytrack.GetN();
      nfound     = nfoundable; 
      //
      for (Int_t ddrow = iter+1; ddrow<drow;ddrow++){
	Float_t maxdist = 0.8*(1.+3./(ddrow));
	for (Int_t delta = -1;delta<=1;delta+=2){
	  Int_t row = row0+ddrow*delta;
	  kr = &(fSectors[sec][row]);
	  Double_t xn = kr->GetX();
	  Double_t ymax = fSectors->GetMaxY(row)-kr->fDeadZone-1.5;
	  polytrack.GetFitPoint(xn,yn,zn);
	  if (TMath::Abs(yn)>ymax) continue;
	  nfoundable++;
	  AliTPCclusterMI * cln = kr->FindNearest(yn,zn,roady,roadz);
	  if (cln) {
	    Float_t dist =  TMath::Sqrt(  (yn-cln->GetY())*(yn-cln->GetY())+(zn-cln->GetZ())*(zn-cln->GetZ()));
	    if (dist<maxdist){
	      /*
	      erry = (dist+0.3)*cln->GetSigmaY2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));	    
	      errz = (dist+0.3)*cln->GetSigmaZ2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));
	      if (cln->IsUsed(10)) {
		//	printf("used\n");
		sumused++;
		erry*=2;
		errz*=2;
	      }
	      */
	      erry=0.1;
	      errz=0.1;
	      polytrack.AddPoint(xn,cln->GetY(),cln->GetZ(),erry, errz);
	      nfound++;
	    }
	  }
	}
	if ( (sumused>3) || (sumused>0.5*nfound) || (nfound<0.6*nfoundable))  break;     
	polytrack.UpdateParameters();
      }           
    }
    if ( (sumused>3) || (sumused>0.5*nfound))  {
      //printf("sumused   %d\n",sumused);
      continue;
    }
    nin1++;
    Double_t dy,dz;
    polytrack.GetFitDerivation(kr0.GetX(),dy,dz);
    AliTPCpolyTrack track2;
    
    polytrack.Refit(track2,0.5+TMath::Abs(dy)*0.3,0.4+TMath::Abs(dz)*0.3);
    if (track2.GetN()<0.5*nfoundable) continue;
    nin2++;

    if ((nfound>0.6*nfoundable) &&( nfoundable>0.4*(i1-i2))) {
      //
      // test seed with and without constrain
      for (Int_t constrain=0; constrain<=0;constrain++){
	// add polytrack candidate

	Double_t x[5], c[15];
	Double_t x1,x2,x3,y1,y2,y3,z1,z2,z3;
	track2.GetBoundaries(x3,x1);	
	x2 = (x1+x3)/2.;
	track2.GetFitPoint(x1,y1,z1);
	track2.GetFitPoint(x2,y2,z2);
	track2.GetFitPoint(x3,y3,z3);
	//
	//is track pointing to the vertex ?
	Double_t x0,y0,z0;
	x0=0;
	polytrack.GetFitPoint(x0,y0,z0);

	if (constrain) {
	  x2 = x3;
	  y2 = y3;
	  z2 = z3;
	  
	  x3 = 0;
	  y3 = 0;
	  z3 = 0;
	}
	x[0]=y1;
	x[1]=z1;
	x[4]=F1(x1,y1,x2,y2,x3,y3);
		
	//	if (TMath::Abs(x[4]) >= cuts[0]) continue;  //
	x[2]=F2(x1,y1,x2,y2,x3,y3);
	
	//if (TMath::Abs(x[4]*x1-x[2]) >= cuts[1]) continue;
	//x[3]=F3(x1,y1,x2,y2,z1,z2);
	x[3]=F3n(x1,y1,x3,y3,z1,z3,x[4]);
	//if (TMath::Abs(x[3]) > cuts[3]) continue;

	
	Double_t sy =0.1, sz =0.1;
	Double_t sy1=0.02, sz1=0.02;
	Double_t sy2=0.02, sz2=0.02;
	Double_t sy3=0.02;

	if (constrain){
	  sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	}
	
	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;

	Double_t f30=(F3(x1,y1+sy,x3,y3,z1,z3)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x3,y3,z1+sz,z3)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x3,y3+sy,z1,z3)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x3,y3,z1,z3+sz)-x[3])/sz;

	
	c[0]=sy1;
	c[1]=0.;       c[2]=sz1;
	c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
	c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
	c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
	c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
	c[13]=f30*sy1*f40+f32*sy2*f42;
	c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//Int_t row1 = fSectors->GetRowNumber(x1);
	Int_t row1 = GetRowNumber(x1);

	UInt_t index=0;
	//kr0.GetIndex(is);
	AliTPCseed *track=new (seed) AliTPCseed(index, x, c, x1, sec*alpha+shift);
	track->fIsSeeding = kTRUE;
	Int_t rc=FollowProlongation(*track, i2);	
	if (constrain) track->fBConstrain =1;
	else
	  track->fBConstrain =0;
	track->fLastPoint = row1+fInnerSec->GetNRows();  // first cluster in track position
	track->fFirstPoint = track->fLastPoint;

	if (rc==0 || track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->fNFoundable*0.6 || 
	    track->fNShared>0.4*track->GetNumberOfClusters()) {
	  //delete track;
	  seed->Reset();
	  seed->~AliTPCseed();
	}
	else {
	  arr->AddLast(track);
	  seed = new AliTPCseed;
	}
	nin3++;
      }
    }  // if accepted seed
  }
  if (fDebug>3){
    Info("MakeSeeds2","\nSeeding statiistic:\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin3);
  }
  delete seed;
}


AliTPCseed *AliTPCtrackerMI::MakeSeed(AliTPCseed *track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using track points
  Int_t p0 = int(r0*track->GetNumberOfClusters());     // point 0 
  Int_t p1 = int(r1*track->GetNumberOfClusters());
  Int_t p2 = int(r2*track->GetNumberOfClusters());   // last point
  Int_t pp2=0;
  Double_t  x0[3],x1[3],x2[3];
  x0[0]=-1;
  x0[0]=-1;
  x0[0]=-1;

  // find track position at given ratio of the length
  Int_t  sec0, sec1, sec2;
  sec0=0;
  sec1=0;
  sec2=0;
  Int_t index=-1;
  Int_t clindex;
  for (Int_t i=0;i<160;i++){
    if (track->fClusterPointer[i]){
      index++;
      AliTPCTrackerPoint   *trpoint =track->GetTrackPoint(i);
      if ( (index<p0) || x0[0]<0 ){
	if (trpoint->GetX()>1){
	  clindex = track->GetClusterIndex2(i);
	  if (clindex>0){	
	    x0[0] = trpoint->GetX();
	    x0[1] = trpoint->GetY();
	    x0[2] = trpoint->GetZ();
	    sec0  = ((clindex&0xff000000)>>24)%18;
	  }
	}
      }

      if ( (index<p1) &&(trpoint->GetX()>1)){
	clindex = track->GetClusterIndex2(i);
	if (clindex>0){
	  x1[0] = trpoint->GetX();
	  x1[1] = trpoint->GetY();
	  x1[2] = trpoint->GetZ();
	  sec1  = ((clindex&0xff000000)>>24)%18;
	}
      }
      if ( (index<p2) &&(trpoint->GetX()>1)){
	clindex = track->GetClusterIndex2(i);
	if (clindex>0){
	  x2[0] = trpoint->GetX();
	  x2[1] = trpoint->GetY();
	  x2[2] = trpoint->GetZ(); 
	  sec2  = ((clindex&0xff000000)>>24)%18;
	  pp2 = i;
	}
      }
    }
  }
  
  Double_t alpha, cs,sn, xx2,yy2;
  //
  alpha = (sec1-sec2)*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= x1[0]*cs-x1[1]*sn;
  yy2= x1[0]*sn+x1[1]*cs;
  x1[0] = xx2;
  x1[1] = yy2;
  //
  alpha = (sec0-sec2)*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= x0[0]*cs-x0[1]*sn;
  yy2= x0[0]*sn+x0[1]*cs;
  x0[0] = xx2;
  x0[1] = yy2;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=x2[1];
  x[1]=x2[2];
  x[4]=F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  //  if (x[4]>1) return 0;
  x[2]=F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  x[3]=F3n(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2],x[4]);
  //if (TMath::Abs(x[3]) > 2.2)  return 0;
  //if (TMath::Abs(x[2]) > 1.99) return 0;
  //  
  Double_t sy =0.1,  sz =0.1;
  //
  Double_t sy1=0.02+track->GetSigmaY2(), sz1=0.02+track->GetSigmaZ2();
  Double_t sy2=0.01+track->GetSigmaY2(), sz2=0.01+track->GetSigmaZ2();
  Double_t sy3=0.01+track->GetSigmaY2();
  //
  Double_t f40=(F1(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[4])/sy;
  Double_t f42=(F1(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[4])/sy;
  Double_t f43=(F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[4])/sy;
  Double_t f20=(F2(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[2])/sy;
  Double_t f22=(F2(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[2])/sy;
  Double_t f23=(F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(x2[0],x2[1]+sy,x0[0],x0[1],x2[2],x0[2])-x[3])/sy;
  Double_t f31=(F3(x2[0],x2[1],x0[0],x0[1],x2[2]+sz,x0[2])-x[3])/sz;
  Double_t f32=(F3(x2[0],x2[1],x0[0],x0[1]+sy,x2[2],x0[2])-x[3])/sy;
  Double_t f34=(F3(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  //  Int_t row1 = fSectors->GetRowNumber(x2[0]);
  AliTPCseed *seed=new  AliTPCseed(0, x, c, x2[0], sec2*fSectors->GetAlpha()+fSectors->GetAlphaShift());
  //  Double_t y0,z0,y1,z1, y2,z2;
  //seed->GetProlongation(x0[0],y0,z0);
  // seed->GetProlongation(x1[0],y1,z1);
  //seed->GetProlongation(x2[0],y2,z2);
  //  seed =0;
  seed->fLastPoint  = pp2;
  seed->fFirstPoint = pp2;
  

  return seed;
}


AliTPCseed *AliTPCtrackerMI::ReSeed(AliTPCseed *track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using founded clusters 
  //
  // Find the number of clusters
  Int_t nclusters = 0;
  for (Int_t irow=0;irow<160;irow++){
    if (track->GetClusterIndex(irow)>0) nclusters++;
  }
  //
  Int_t ipos[3];
  ipos[0] = TMath::Max(int(r0*nclusters),0);             // point 0 cluster
  ipos[1] = TMath::Min(int(r1*nclusters),nclusters-1);   // 
  ipos[2] = TMath::Min(int(r2*nclusters),nclusters-1);   // last point
  //
  //
  Double_t  xyz[3][3];
  Int_t     row[3],sec[3]={0,0,0};
  //
  // find track row position at given ratio of the length
  Int_t index=-1;
  for (Int_t irow=0;irow<160;irow++){    
    if (track->GetClusterIndex2(irow)<0) continue;
    index++;
    for (Int_t ipoint=0;ipoint<3;ipoint++){
      if (index<=ipos[ipoint]) row[ipoint] = irow;
    }        
  }
  //
  //Get cluster and sector position
  for (Int_t ipoint=0;ipoint<3;ipoint++){
    Int_t clindex = track->GetClusterIndex2(row[ipoint]);
    AliTPCclusterMI * cl = GetClusterMI(clindex);
    if (cl==0) {
      //Error("Bug\n");
      //      AliTPCclusterMI * cl = GetClusterMI(clindex);
      return 0;
    }
    sec[ipoint]     = ((clindex&0xff000000)>>24)%18;
    xyz[ipoint][0]  = GetXrow(row[ipoint]);
    xyz[ipoint][1]  = cl->GetY();
    xyz[ipoint][2]  = cl->GetZ();
  }
  //
  //
  // Calculate seed state vector and covariance matrix

  Double_t alpha, cs,sn, xx2,yy2;
  //
  alpha = (sec[1]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[1][0]*cs-xyz[1][1]*sn;
  yy2= xyz[1][0]*sn+xyz[1][1]*cs;
  xyz[1][0] = xx2;
  xyz[1][1] = yy2;
  //
  alpha = (sec[0]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[0][0]*cs-xyz[0][1]*sn;
  yy2= xyz[0][0]*sn+xyz[0][1]*cs;
  xyz[0][0] = xx2;
  xyz[0][1] = yy2;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=xyz[2][1];
  x[1]=xyz[2][2];
  x[4]=F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[2]=F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[3]=F3n(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1;
  //
  Double_t sy1=0.2, sz1=0.2;
  Double_t sy2=0.2, sz2=0.2;
  Double_t sy3=0.2;
  //
  Double_t f40=(F1(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f42=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f43=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[4])/sy;
  Double_t f20=(F2(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f22=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f23=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(xyz[2][0],xyz[2][1]+sy,xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f31=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2]+sz,xyz[0][2])-x[3])/sz;
  Double_t f32=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1]+sy,xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f34=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  //  Int_t row1 = fSectors->GetRowNumber(xyz[2][0]);
  AliTPCseed *seed=new  AliTPCseed(0, x, c, xyz[2][0], sec[2]*fSectors->GetAlpha()+fSectors->GetAlphaShift());
  seed->fLastPoint  = row[2];
  seed->fFirstPoint = row[2];  
  return seed;
}

Int_t  AliTPCtrackerMI::CheckKinkPoint(AliTPCseed*seed, Float_t th)
{
  //
  //
  // 
  for (Int_t i=0;i<12;i++) seed->fKinkPoint[i]=0;
  //
  if (TMath::Abs(seed->GetC())>0.01) return 0;
  //

  Float_t x[160], y[160], erry[160], z[160], errz[160];
  Int_t sec[160];
  Float_t xt[160], yt[160], zt[160];
  Int_t i1 = 200;
  Int_t i2 = 0;
  Int_t secm   = -1;
  Int_t padm   = -1;
  Int_t middle = seed->GetNumberOfClusters()/2;
  //
  //
  // find central sector, get local cooordinates
  Int_t count = 0;
  for (Int_t i=seed->fFirstPoint;i<=seed->fLastPoint;i++) {
    sec[i]= seed->GetClusterSector(i)%18;
    x[i]  = GetXrow(i);  
    if (sec[i]>=0) {
      AliTPCclusterMI * cl = seed->fClusterPointer[i];
      //      if (cl==0)        cl = GetClusterMI(seed->GetClusterIndex2(i));
      if (cl==0) {
	sec[i] = -1;
	continue;
      }
      //
      //
      if (i>i2)  i2 = i;  //last  point with cluster
      if (i2<i1) i1 = i;  //first point with cluster
      y[i] = cl->GetY();
      z[i] = cl->GetZ();
      AliTPCTrackerPoint * point = seed->GetTrackPoint(i);
      xt[i] = x[i];
      yt[i] = point->GetY();
      zt[i] = point->GetZ();
  
      if (point->GetX()>0){
	erry[i] = point->GetErrY();
	errz[i] = point->GetErrZ();	
      }

      count++;
      if (count<middle) {
	secm = sec[i];  //central sector
	padm = i;       //middle point with cluster
      }
    }
  }
  //
  // rotate position to global coordinate system connected to  sector at last the point
  //
  for (Int_t i=i1;i<=i2;i++){
    //    
    if (sec[i]<0) continue;
    Double_t alpha = (sec[i2]-sec[i])*fSectors->GetAlpha();
    Double_t cs = TMath::Cos(alpha);
    Double_t sn = TMath::Sin(alpha);    
    Float_t xx2= x[i]*cs+y[i]*sn;
    Float_t yy2= -x[i]*sn+y[i]*cs;
    x[i] = xx2;
    y[i] = yy2;    
    //
    xx2= xt[i]*cs+yt[i]*sn;
    yy2= -xt[i]*sn+yt[i]*cs;
    xt[i] = xx2;
    yt[i] = yy2;    

  }
  //get "state" vector
  Double_t xh[5],xm = x[padm];  
  xh[0]=yt[i2];
  xh[1]=zt[i2];
  xh[4]=F1(xt[i2],yt[i2],xt[padm],yt[padm],xt[i1],yt[i1]);  
  xh[2]=F2(xt[i2],yt[i2],xt[padm],yt[padm],xt[i1],yt[i1]);
  xh[3]=F3n(xt[i2],yt[i2],xt[i1],yt[i1],zt[i2],zt[i1],xh[4]);
  //
  //
  for (Int_t i=i1;i<=i2;i++){
    Double_t yy,zz;
    if (sec[i]<0) continue;    
    GetProlongation(x[i2], x[i],xh,yy,zz);
    if (TMath::Abs(y[i]-yy)>4||TMath::Abs(z[i]-zz)>4){
      //Double_t xxh[5];
      //xxh[4]=F1old(x[i2],y[i2],x[padm],y[padm],x[i1],y[i1]);  
      //xxh[2]=F2old(x[i2],y[i2],x[padm],y[padm],x[i1],y[i1]);
      Error("AliTPCtrackerMI::CheckKinkPoint","problem\n");
    }
    y[i] = y[i] - yy;
    z[i] = z[i] - zz;
  }
  Float_t dyup[160],dydown[160], dzup[160], dzdown[160];
  Float_t yup[160], ydown[160],  zup[160],  zdown[160];
 
  AliTPCpolyTrack ptrack1,ptrack2;
  //
  // derivation up
  for (Int_t i=i1;i<=i2;i++){
    AliTPCclusterMI * cl = seed->fClusterPointer[i];
    if (!cl) continue;
    if (cl->GetType()<0) continue;
    if (cl->GetType()>10) continue;

    if (sec[i]>=0){
      ptrack1.AddPoint(x[i]-xm,y[i],z[i],0.1,0.1);
    }
    if (ptrack1.GetN()>4.){
      ptrack1.UpdateParameters();
      Double_t ddy,ddz;
      ptrack1.GetFitDerivation(x[i]-xm,ddy,ddz);
      Double_t yy,zz;
      ptrack1.GetFitPoint(x[i]-xm,yy,zz);

      dyup[i] = ddy;
      dzup[i] = ddz;
      yup[i]  = yy;
      zup[i]  = zz;

    }
    else{
      dyup[i]=0.;  //not enough points
    }
  }
  //
  // derivation down
  for (Int_t i=i2;i>=i1;i--){
    AliTPCclusterMI * cl = seed->fClusterPointer[i];
    if (!cl) continue;
    if (cl->GetType()<0) continue;
    if (cl->GetType()>10) continue;
    if (sec[i]>=0){
      ptrack2.AddPoint(x[i]-xm,y[i],z[i],0.1,0.1);
    }
    if (ptrack2.GetN()>4){
      ptrack2.UpdateParameters();
      Double_t ddy,ddz;
      ptrack2.GetFitDerivation(x[i]-xm,ddy,ddz);
      Double_t yy,zz;
      ptrack2.GetFitPoint(x[i]-xm,yy,zz);

      dydown[i] = ddy;
      dzdown[i] = ddz;
      ydown[i]  = yy;
      zdown[i]  = zz;
    }
    else{
      dydown[i]=0.;  //not enough points
    }
  }
  //
  //
  // find maximal difference of the derivation
  for (Int_t i=0;i<12;i++) seed->fKinkPoint[i]=0;


  for (Int_t i=i1+10;i<i2-10;i++){
    if ( (TMath::Abs(dydown[i])<0.00000001)  ||  (TMath::Abs(dyup[i])<0.00000001) ||i<30)continue;
    //    printf("%f\t%f\t%f\t%f\t%f\n",x[i],dydown[i],dyup[i],dzdown[i],dzup[i]);
    //
    Float_t ddy = TMath::Abs(dydown[i]-dyup[i]);
    Float_t ddz = TMath::Abs(dzdown[i]-dzup[i]);    
    if ( (ddy+ddz)> th){
      seed->fKinkPoint[0] = i;
      seed->fKinkPoint[1] = ddy;
      seed->fKinkPoint[2] = ddz;
      th = ddy+ddz;      
    }
  }

  if (fTreeDebug){
    //
    //write information to the debug tree
    TBranch * br = fTreeDebug->GetBranch("debug");
    TClonesArray * arr = new TClonesArray("AliTPCTrackPoint2");
    arr->ExpandCreateFast(i2-i1);
    br->SetAddress(&arr);
    //
    AliTPCclusterMI cldummy;
    cldummy.SetQ(0);
    AliTPCTrackPoint2 pdummy;
    pdummy.GetTPoint().fIsShared = 10;
    //
    Double_t alpha = sec[i2]*fSectors->GetAlpha();
    Double_t cs    = TMath::Cos(alpha);
    Double_t sn    = TMath::Sin(alpha);    

    for (Int_t i=i1;i<i2;i++){
      AliTPCTrackPoint2 *trpoint = (AliTPCTrackPoint2*)arr->UncheckedAt(i-i1);
      //cluster info
      AliTPCclusterMI * cl0 = seed->fClusterPointer[i];
      //      
      AliTPCTrackerPoint * point = seed->GetTrackPoint(i);
      
      if (cl0){
	Double_t x = GetXrow(i);
	trpoint->GetTPoint() = *point;
	trpoint->GetCPoint() = *cl0;
	trpoint->GetCPoint().SetQ(TMath::Abs(cl0->GetQ()));
	trpoint->fID    = seed->GetUniqueID();
	trpoint->fLab   = seed->GetLabel();
	//
	trpoint->fGX =  cs *x + sn*point->GetY();
	trpoint->fGY = -sn *x + cs*point->GetY() ;
	trpoint->fGZ = point->GetZ();
	//
	trpoint->fDY = y[i];
	trpoint->fDZ = z[i];
	//
	trpoint->fDYU = dyup[i];
	trpoint->fDZU = dzup[i];
	//
	trpoint->fDYD = dydown[i];
	trpoint->fDZD = dzdown[i];
	//
	if (TMath::Abs(dyup[i])>0.00000000001 &&TMath::Abs(dydown[i])>0.00000000001){
	  trpoint->fDDY = dydown[i]-dyup[i];
	  trpoint->fDDZ = dzdown[i]-dzup[i];
	}else{
	  trpoint->fDDY = 0.;
	  trpoint->fDDZ = 0.;
	}	
      }
      else{
	*trpoint = pdummy;
	trpoint->GetCPoint()= cldummy;
	trpoint->fID = -1;
      }
      //     
    }
    fTreeDebug->Fill();
  }
  
  
  return 0;
  
}





AliTPCseed*  AliTPCtrackerMI::ReSeed(AliTPCseed *t)
{
  //
  // reseed - refit -  track
  //
  Int_t first = 0;
  //  Int_t last  = fSectors->GetNRows()-1;
  //
  if (fSectors == fOuterSec){
    first = TMath::Max(first, t->fFirstPoint-fInnerSec->GetNRows());
    //last  = 
  }
  else
    first = t->fFirstPoint;
  //
  AliTPCseed * seed = MakeSeed(t,0.1,0.5,0.9);
  FollowBackProlongation(*t,fSectors->GetNRows()-1);
  t->Reset(kFALSE);
  FollowProlongation(*t,first);
  return seed;
}







//_____________________________________________________________________________
Int_t AliTPCtrackerMI::ReadSeeds(const TFile *inp) {
  //-----------------------------------------------------------------
  // This function reades track seeds.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  TFile *in=(TFile*)inp;
  if (!in->IsOpen()) {
     cerr<<"AliTPCtrackerMI::ReadSeeds(): input file is not open !\n";
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     cerr<<"AliTPCtrackerMI::ReadSeeds(): ";
     cerr<<"can't get a tree with track seeds !\n";
     return 2;
  }
  AliTPCtrack *seed=new AliTPCtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     fSeeds->AddLast(new AliTPCseed(*seed,seed->GetAlpha()));
  }
  
  delete seed;
  delete seedTree; 
  savedir->cd();
  return 0;
}

Int_t AliTPCtrackerMI::Clusters2Tracks (AliESD *esd)
{
  //
  if (fSeeds) DeleteSeeds();
  fEvent = esd;
  Clusters2Tracks();
  if (!fSeeds) return 1;
  FillESD(fSeeds);
  return 0;
  //
}


//_____________________________________________________________________________
Int_t AliTPCtrackerMI::Clusters2Tracks() {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 
  TStopwatch timer;

  fIteration = 0;
  fSeeds = Tracking();

  if (fDebug>0){
    Info("Clusters2Tracks","Time for tracking: \t");timer.Print();timer.Start();
  }
  //activate again some tracks
  for (Int_t i=0; i<fSeeds->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<20) {
      delete fSeeds->RemoveAt(i);
      continue;
    }
    if (pt->fRemoval==10) {
      if (pt->GetDensityFirst(20)>0.8 || pt->GetDensityFirst(30)>0.8 || pt->GetDensityFirst(40)>0.7)
	pt->Desactivate(10);  // make track again active
      else{
	pt->Desactivate(20); 	
	delete fSeeds->RemoveAt(i);
      }
    } 
  }
  RemoveDouble(fSeeds,0.2,0.6,11);
  RemoveUsed(fSeeds,0.5,0.5,6);

  //
  Int_t nseed=fSeeds->GetEntriesFast();
  Int_t found = 0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<15) {
      delete fSeeds->RemoveAt(i);
      continue;
    }
    CookLabel(pt,0.1); //For comparison only
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->fRemoval==10) )){
      found++;      
      if (fDebug>0) cerr<<found<<'\r';      
      pt->fLab2 = i;
    }
    else
      delete fSeeds->RemoveAt(i);
  }

  
  //RemoveOverlap(fSeeds,0.99,7,kTRUE);  
  SignShared(fSeeds);  
  //RemoveUsed(fSeeds,0.9,0.9,6);
  // 
  nseed=fSeeds->GetEntriesFast();
  found = 0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<15) {
      delete fSeeds->RemoveAt(i);
      continue;
    }
    t.SetUniqueID(i);
    t.CookdEdx(0.02,0.6);
    //    CheckKinkPoint(&t,0.05);
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->fRemoval==10) )){
      found++;
      if (fDebug>0){
	cerr<<found<<'\r';      
      }
      pt->fLab2 = i;
    }
    else
      delete fSeeds->RemoveAt(i);
    //AliTPCseed * seed1 = ReSeed(pt,0.05,0.5,1);
    //if (seed1){
    //  FollowProlongation(*seed1,0);
    //  Int_t n = seed1->GetNumberOfClusters();
    //  printf("fP4\t%f\t%f\n",seed1->GetC(),pt->GetC());
    //  printf("fN\t%d\t%d\n", seed1->GetNumberOfClusters(),pt->GetNumberOfClusters());
    //
    //}
    //AliTPCseed * seed2 = ReSeed(pt,0.95,0.5,0.05);
    
  }

  SortTracks(fSeeds, 1);
  
  /*    
  fIteration = 1;
  PrepareForBackProlongation(fSeeds,5.);
  PropagateBack(fSeeds);
  printf("Time for back propagation: \t");timer.Print();timer.Start();
  
  fIteration = 2;
  
  PrepareForProlongation(fSeeds,5.);
  PropagateForward2(fSeeds);
   
  printf("Time for FORWARD propagation: \t");timer.Print();timer.Start();
  // RemoveUsed(fSeeds,0.7,0.7,6);
  //RemoveOverlap(fSeeds,0.9,7,kTRUE);
   
  nseed=fSeeds->GetEntriesFast();
  found = 0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<15) {
      delete fSeeds->RemoveAt(i);
      continue;
    }
    t.CookdEdx(0.02,0.6);
    //    CookLabel(pt,0.1); //For comparison only
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->fRemoval==10) )){
      cerr<<found++<<'\r';      
    }
    else
      delete fSeeds->RemoveAt(i);
    pt->fLab2 = i;
  }
  */
 
  //  fNTracks = found;
  if (fDebug>0){
    Info("Clusters2Tracks","Time for overlap removal, track writing and dedx cooking: \t"); timer.Print();timer.Start();
  }
  //
  //  cerr<<"Number of found tracks : "<<"\t"<<found<<endl;  
  Info("Clusters2Tracks","Number of found tracks %d",found);  
  savedir->cd();
  //  UnloadClusters();
  //  
  return 0;
}

void AliTPCtrackerMI::Tracking(TObjArray * arr)
{
  //
  // tracking of the seeds
  //

  fSectors = fOuterSec;
  ParallelTracking(arr,150,63);
  fSectors = fOuterSec;
  ParallelTracking(arr,63,0);
}

TObjArray * AliTPCtrackerMI::Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4], Float_t dy, Int_t dsec)
{
  //
  //
  //tracking routine
  TObjArray * arr = new TObjArray;
  // 
  fSectors = fOuterSec;
  TStopwatch timer;
  timer.Start();
  for (Int_t sec=0;sec<fkNOS;sec++){
    if (seedtype==3) MakeSeeds3(arr,sec,i1,i2,cuts,dy, dsec);
    if (seedtype==4) MakeSeeds5(arr,sec,i1,i2,cuts,dy);    
    if (seedtype==2) MakeSeeds2(arr,sec,i1,i2,cuts,dy);
  }
  if (fDebug>0){
    Info("Tracking","\nSeeding - %d\t%d\t%d\t%d\n",seedtype,i1,i2,arr->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  Tracking(arr);  
  if (fDebug>0){
    timer.Print();
  }

  return arr;
}

TObjArray * AliTPCtrackerMI::Tracking()
{
  //
  //
  TStopwatch timer;
  timer.Start();
  Int_t nup=fOuterSec->GetNRows()+fInnerSec->GetNRows();

  TObjArray * seeds = new TObjArray;
  TObjArray * arr=0;
  
  Int_t gap =20;
  Float_t cuts[4];
  cuts[0] = 0.002;
  cuts[1] = 1.5;
  cuts[2] = 3.;
  cuts[3] = 3.;
  Float_t fnumber  = 3.0;
  Float_t fdensity = 3.0;
  
  //  
  //find primaries  
  cuts[0]=0.0066;
  for (Int_t delta = 0; delta<18; delta+=6){
    //
    cuts[0]=0.0070;
    cuts[1] = 1.5;
    arr = Tracking(3,nup-1-delta,nup-1-delta-gap,cuts,-1,1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity); 
    //
    for (Int_t i=2;i<6;i+=2){
      // seed high pt tracks
      cuts[0]=0.0022;
      cuts[1]=0.3;
      arr = Tracking(3,nup-i-delta,nup-i-delta-gap,cuts,-1,0);
      SumTracks(seeds,arr);   
      SignClusters(seeds,fnumber,fdensity);        
    }
  }
  fnumber  = 4;
  fdensity = 4.;
  //  RemoveUsed(seeds,0.9,0.9,1);
  //  UnsignClusters();
  //  SignClusters(seeds,fnumber,fdensity);    

  //find primaries  
  cuts[0]=0.0077;
  for (Int_t delta = 20; delta<120; delta+=10){
    //
    // seed high pt tracks
    cuts[0]=0.0060;
    cuts[1]=0.3;
    cuts[2]=6.;
    arr = Tracking(3,nup-delta,nup-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            

    cuts[0]=0.003;
    cuts[1]=0.3;
    cuts[2]=6.;
    arr = Tracking(3,nup-delta-5,nup-delta-5-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            
  }

  cuts[0] = 0.01;
  cuts[1] = 2.0;
  cuts[2] = 3.;
  cuts[3] = 2.0;
  fnumber  = 2.;
  fdensity = 2.;
  
  if (fDebug>0){
    Info("Tracking()","\n\nPrimary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  //  RemoveUsed(seeds,0.75,0.75,1);
  //UnsignClusters();
  //SignClusters(seeds,fnumber,fdensity);
  
  // find secondaries

  cuts[0] = 0.3;
  cuts[1] = 1.5;
  cuts[2] = 3.;
  cuts[3] = 1.5;

  arr = Tracking(4,nup-1,nup-1-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-2,nup-2-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-3,nup-3-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //


  for (Int_t delta = 3; delta<30; delta+=5){
    //
    cuts[0] = 0.3;
    cuts[1] = 1.5;
    cuts[2] = 3.;
    cuts[3] = 1.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
    //
    arr = Tracking(4,nup-3-delta,nup-5-delta-gap,cuts,4);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity); 
    //
  } 
  fnumber  = 1;
  fdensity = 1;
  //
  // change cuts
  fnumber  = 2.;
  fdensity = 2.;
  cuts[0]=0.0080;

  // find secondaries
  for (Int_t delta = 30; delta<70; delta+=10){
    //
    cuts[0] = 0.3;
    cuts[1] = 1.5;
    cuts[2] = 3.;
    cuts[3] = 1.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
    //
    arr = Tracking(4,nup-5-delta,nup-5-delta-gap,cuts,5 );
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  }
 
  if (fDebug>0){
    Info("Tracking()","\n\nSecondary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }

  return seeds;
  //
      
}


void AliTPCtrackerMI::SumTracks(TObjArray *arr1,TObjArray *arr2) const
{
  //
  //sum tracks to common container
  //remove suspicious tracks
  Int_t nseed = arr2->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt=(AliTPCseed*)arr2->UncheckedAt(i);    
    if (pt){
      
      // NORMAL ACTIVE TRACK
      if (pt->IsActive()){
	arr1->AddLast(arr2->RemoveAt(i));
	continue;
      }
      //remove not usable tracks
      if (pt->fRemoval!=10){
	delete arr2->RemoveAt(i);
	continue;
      }
      // REMOVE VERY SHORT  TRACKS
      if (pt->GetNumberOfClusters()<20){ 
	delete arr2->RemoveAt(i);
	continue;
      }
      // ENABLE ONLY ENOUGH GOOD STOPPED TRACKS
      if (pt->GetDensityFirst(20)>0.8 || pt->GetDensityFirst(30)>0.8 || pt->GetDensityFirst(40)>0.7)
	arr1->AddLast(arr2->RemoveAt(i));
      else{      
	delete arr2->RemoveAt(i);
      }
    }
  }
  delete arr2;  
}



void  AliTPCtrackerMI::ParallelTracking(TObjArray * arr, Int_t rfirst, Int_t rlast)
{
  //
  // try to track in parralel

  Int_t nseed=arr->GetEntriesFast();
  //prepare seeds for tracking
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i), &t=*pt; 
    if (!pt) continue;
    if (!t.IsActive()) continue;
    // follow prolongation to the first layer
    if ( (fSectors ==fInnerSec) || (t.fFirstPoint-fParam->GetNRowLow()>rfirst+1) )  
      FollowProlongation(t, rfirst+1);
  }


  //
  for (Int_t nr=rfirst; nr>=rlast; nr--){ 
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;
    // make indexes with the cluster tracks for given       

    // find nearest cluster
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i), &t=*pt;       
      if (!pt) continue;
      if (!pt->IsActive()) continue;
      //      if ( (fSectors ==fOuterSec) && (pt->fFirstPoint-fParam->GetNRowLow())<nr) continue;
      if (pt->fRelativeSector>17) {
	continue;
      }
      UpdateClusters(t,nr);
    }
    // prolonagate to the nearest cluster - if founded
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i); 
      if (!pt) continue;
      if (!pt->IsActive()) continue; 
      // if ((fSectors ==fOuterSec) && (pt->fFirstPoint-fParam->GetNRowLow())<nr) continue;
      if (pt->fRelativeSector>17) {
	continue;
      }
      FollowToNextCluster(*pt,nr);
    }
  }    
}

void AliTPCtrackerMI::PrepareForBackProlongation(TObjArray * arr,Float_t fac) const
{
  //
  //
  // if we use TPC track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
      //
      //rotate to current local system at first accepted  point    
      Int_t index  = pt->GetClusterIndex2(pt->fFirstPoint); 
      Int_t sec    = (index&0xff000000)>>24;
      sec = sec%18;
      Float_t angle1 = fInnerSec->GetAlpha()*sec+fInnerSec->GetAlphaShift();
      if (angle1>TMath::Pi()) 
	angle1-=2.*TMath::Pi();
      Float_t angle2 = pt->GetAlpha();
      
      if (TMath::Abs(angle1-angle2)>0.001){
	pt->Rotate(angle1-angle2);
	//angle2 = pt->GetAlpha();
	//pt->fRelativeSector = pt->GetAlpha()/fInnerSec->GetAlpha();
	//if (pt->GetAlpha()<0) 
	//  pt->fRelativeSector+=18;
	//sec = pt->fRelativeSector;
      }
	
    }
    
  }


}
void AliTPCtrackerMI::PrepareForProlongation(TObjArray * arr, Float_t fac) const
{
  //
  //
  // if we use TPC track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
      pt->fFirstPoint = pt->fLastPoint; 
    }
    
  }


}

Int_t AliTPCtrackerMI::PropagateBack(TObjArray * arr)
{
  //
  // make back propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) { 
      //AliTPCseed *pt2 = new AliTPCseed(*pt);
      fSectors = fInnerSec;
      //FollowBackProlongation(*pt,fInnerSec->GetNRows()-1);
      //fSectors = fOuterSec;
      FollowBackProlongation(*pt,fInnerSec->GetNRows()+fOuterSec->GetNRows()-1);     
      //if (pt->GetNumberOfClusters()<(pt->fEsd->GetTPCclusters(0)) ){
      //	Error("PropagateBack","Not prolonged track %d",pt->GetLabel());
      //	FollowBackProlongation(*pt2,fInnerSec->GetNRows()+fOuterSec->GetNRows()-1);
      //}
    }      
  }
  return 0;
}


Int_t AliTPCtrackerMI::PropagateForward2(TObjArray * arr)
{
  //
  // make forward propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  //
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) { 
      FollowProlongation(*pt,0);
    }
  }
  return 0;
}


Int_t AliTPCtrackerMI::PropagateForward()
{
  //
  // propagate track forward
  //UnsignClusters();
  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)fSeeds->UncheckedAt(i);
    if (pt){
      AliTPCseed &t = *pt;
      Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
      if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
      if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
      t.fRelativeSector = Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN;
    }
  }
  
  fSectors = fOuterSec;
  ParallelTracking(fSeeds,fOuterSec->GetNRows()+fInnerSec->GetNRows()-1,fInnerSec->GetNRows());
  fSectors = fInnerSec;
  ParallelTracking(fSeeds,fInnerSec->GetNRows()-1,0);
  //WriteTracks();
  return 1;
}






Int_t AliTPCtrackerMI::PropagateBack(AliTPCseed * pt, Int_t row0, Int_t row1)
{
  //
  // make back propagation, in between row0 and row1
  //
  
  if (pt) { 
    fSectors = fInnerSec;
    Int_t  r1;
    //
    if (row1<fSectors->GetNRows()) 
      r1 = row1;
    else 
      r1 = fSectors->GetNRows()-1;

    if (row0<fSectors->GetNRows()&& r1>0 )
      FollowBackProlongation(*pt,r1);
    if (row1<=fSectors->GetNRows())
      return 0;
    //
    r1 = row1 - fSectors->GetNRows();
    if (r1<=0) return 0;
    if (r1>=fOuterSec->GetNRows()) return 0;
    fSectors = fOuterSec;
    return FollowBackProlongation(*pt,r1);
  }        
  return 0;
}




void  AliTPCtrackerMI::GetShape(AliTPCseed * seed, Int_t row)
{
  //
  //
  Float_t sd2 = TMath::Abs((fParam->GetZLength()-TMath::Abs(seed->GetZ())))*fParam->GetDiffL()*fParam->GetDiffL();
  //  Float_t padlength =  fParam->GetPadPitchLength(seed->fSector);
  Float_t padlength =  GetPadPitchLength(row);
  //
  Float_t sresy = (seed->fSector < fParam->GetNSector()/2) ? 0.2 :0.3;
  Float_t angulary  = seed->GetSnp();
  angulary = angulary*angulary/(1-angulary*angulary);
  seed->fCurrentSigmaY2 = sd2+padlength*padlength*angulary/12.+sresy*sresy;  
  //
  Float_t sresz = fParam->GetZSigma();
  Float_t angularz  = seed->GetTgl();
  seed->fCurrentSigmaZ2 = sd2+padlength*padlength*angularz*angularz*(1+angulary)/12.+sresz*sresz;
  /*
  Float_t wy = GetSigmaY(seed);
  Float_t wz = GetSigmaZ(seed);
  wy*=wy;
  wz*=wz;
  if (TMath::Abs(wy/seed->fCurrentSigmaY2-1)>0.0001 || TMath::Abs(wz/seed->fCurrentSigmaZ2-1)>0.0001 ){
    printf("problem\n");
  }
  */
}


Float_t  AliTPCtrackerMI::GetSigmaY(AliTPCseed * seed)
{
  //
  //  
  Float_t sd2 = TMath::Abs((fParam->GetZLength()-TMath::Abs(seed->GetZ())))*fParam->GetDiffL()*fParam->GetDiffL();
  Float_t padlength =  fParam->GetPadPitchLength(seed->fSector);
  Float_t sres = (seed->fSector < fParam->GetNSector()/2) ? 0.2 :0.3;
  Float_t angular  = seed->GetSnp();
  angular = angular*angular/(1-angular*angular);
  //  angular*=angular;
  //angular  = TMath::Sqrt(angular/(1-angular));
  Float_t res = TMath::Sqrt(sd2+padlength*padlength*angular/12.+sres*sres);
  return res;
}
Float_t  AliTPCtrackerMI::GetSigmaZ(AliTPCseed * seed)
{
  //
  //
  Float_t sd2 = TMath::Abs((fParam->GetZLength()-TMath::Abs(seed->GetZ())))*fParam->GetDiffL()*fParam->GetDiffL();
  Float_t padlength =  fParam->GetPadPitchLength(seed->fSector);
  Float_t sres = fParam->GetZSigma();
  Float_t angular  = seed->GetTgl();
  Float_t res = TMath::Sqrt(sd2+padlength*padlength*angular*angular/12.+sres*sres);
  return res;
}




//__________________________________________________________________________
void AliTPCtrackerMI::CookLabel(AliTPCseed *t, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  if (noc<10){
    //printf("\nnot founded prolongation\n\n\n");
    //t->Dump();
    return ;
  }
  Int_t lb[160];
  Int_t mx[160];
  AliTPCclusterMI *clusters[160];
  //
  for (Int_t i=0;i<160;i++) {
    clusters[i]=0;
    lb[i]=mx[i]=0;
  }

  Int_t i;
  Int_t current=0;
  for (i=0; i<160 && current<noc; i++) {
     
     Int_t index=t->GetClusterIndex2(i);
     if (index<=0) continue; 
     if (index&0x8000) continue;
     //     
     //clusters[current]=GetClusterMI(index);
     if (t->fClusterPointer[i]){
       clusters[current]=t->fClusterPointer[i];     
       current++;
     }
  }
  noc = current;

  Int_t lab=123456789;
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i];
    if (!c) continue;
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<noc; j++) if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }

  Int_t max=0;
  for (i=0; i<noc; i++) if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i]; 
    if (!c) continue;
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }

  if ((1.- Float_t(max)/noc) > wrong) lab=-lab;

  else {
     Int_t tail=Int_t(0.10*noc);
     max=0;
     Int_t ind=0;
     for (i=1; i<=160&&ind<tail; i++) {
       //       AliTPCclusterMI *c=clusters[noc-i];
       AliTPCclusterMI *c=clusters[i];
       if (!c) continue;
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
       ind++;
     }
     if (max < Int_t(0.5*tail)) lab=-lab;
  }

  t->SetLabel(lab);

  //  delete[] lb;
  //delete[] mx;
  //delete[] clusters;
}


Int_t  AliTPCtrackerMI::AliTPCSector::GetRowNumber(Double_t x) const 
{
  //return pad row number for this x
  Double_t r;
  if (fN < 64){
    r=fRow[fN-1].GetX();
    if (x > r) return fN;
    r=fRow[0].GetX();
    if (x < r) return -1;
    return Int_t((x-r)/fPadPitchLength + 0.5);}
  else{    
    r=fRow[fN-1].GetX();
    if (x > r) return fN;
    r=fRow[0].GetX();
    if (x < r) return -1;
    Double_t r1=fRow[64].GetX();
    if(x<r1){       
      return Int_t((x-r)/f1PadPitchLength + 0.5);}
    else{
      return (Int_t((x-r1)/f2PadPitchLength + 0.5)+64);} 
  }
}

//_________________________________________________________________________
void AliTPCtrackerMI::AliTPCSector::Setup(const AliTPCParam *par, Int_t f) {
  //-----------------------------------------------------------------------
  // Setup inner sector
  //-----------------------------------------------------------------------
  if (f==0) {
     fAlpha=par->GetInnerAngle();
     fAlphaShift=par->GetInnerAngleShift();
     fPadPitchWidth=par->GetInnerPadPitchWidth();
     fPadPitchLength=par->GetInnerPadPitchLength();
     fN=par->GetNRowLow();
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiLow(i));
       fRow[i].fDeadZone =1.5;  //1.5 cm of dead zone
     }
  } else {
     fAlpha=par->GetOuterAngle();
     fAlphaShift=par->GetOuterAngleShift();
     fPadPitchWidth  = par->GetOuterPadPitchWidth();
     fPadPitchLength = par->GetOuter1PadPitchLength();
     f1PadPitchLength = par->GetOuter1PadPitchLength();
     f2PadPitchLength = par->GetOuter2PadPitchLength();

     fN=par->GetNRowUp();
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiUp(i)); 
       fRow[i].fDeadZone =1.5;  // 1.5 cm of dead zone
     }
  } 
}

AliTPCtrackerMI::AliTPCRow::AliTPCRow() {
  //
  // default constructor
  fN=0;
  fN1=0;
  fN2=0;
  fClusters1=0;
  fClusters2=0;
}

AliTPCtrackerMI::AliTPCRow::~AliTPCRow(){
  //

}



//_________________________________________________________________________
void 
AliTPCtrackerMI::AliTPCRow::InsertCluster(const AliTPCclusterMI* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMaxClusterPerRow) {
    cerr<<"AliTPCRow::InsertCluster(): Too many clusters !\n"; return;
  }
  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTPCclusterMI*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  fIndex[i]=index; fClusters[i]=c; fN++;
}

void AliTPCtrackerMI::AliTPCRow::ResetClusters() {
   //
   // reset clusters
   fN  = 0; 
   fN1 = 0;
   fN2 = 0;
   //delete[] fClusterArray; 
   if (fClusters1) delete []fClusters1; 
   if (fClusters2) delete []fClusters2; 
   //fClusterArray=0;
   fClusters1 = 0;
   fClusters2 = 0;
}


//___________________________________________________________________
Int_t AliTPCtrackerMI::AliTPCRow::Find(Double_t z) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster 
  //-----------------------------------------------------------------------
  if (fN==0) return 0;
  if (z <= fClusters[0]->GetZ()) return 0;
  if (z > fClusters[fN-1]->GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}



//___________________________________________________________________
AliTPCclusterMI * AliTPCtrackerMI::AliTPCRow::FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;

  AliTPCclusterMI *cl =0;
  for (Int_t i=Find(z-roadz); i<fN; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( (c->GetY()-y) >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
      }
  }
  return cl;      
}

AliTPCclusterMI * AliTPCtrackerMI::AliTPCRow::FindNearest2(Double_t y, Double_t z, Double_t roady, Double_t roadz,UInt_t & index) const 
{
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;
  Int_t iz1 = TMath::Max(fFastCluster[Int_t(z-roadz+254.5)]-1,0);
  Int_t iz2 = TMath::Min(fFastCluster[Int_t(z+roadz+255.5)]+1,fN);

  AliTPCclusterMI *cl =0;
  //FindNearest3(y,z,roady,roadz,index);
  //  for (Int_t i=Find(z-roadz); i<fN; i++) {
  for (Int_t i=iz1; i<iz2; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( c->GetY()-y >  roady ) continue;
      if ( y-c->GetY() >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
	index =i;
	//roady = TMath::Sqrt(maxdistance);
      }
  }
  return cl;      
}



AliTPCclusterMI * AliTPCtrackerMI::AliTPCRow::FindNearest3(Double_t y, Double_t z, Double_t roady, Double_t roadz,UInt_t & index) const 
{
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;
  //  Int_t iz = Int_t(z+255.);
  AliTPCclusterMI *cl =0;
  for (Int_t i=Find(z-roadz); i<fN; i++) {
    //for (Int_t i=fFastCluster[iz-2]; i<fFastCluster[iz+2]; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( c->GetY()-y >  roady ) continue;
      if ( y-c->GetY() >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
	index =i;
	//roady = TMath::Sqrt(maxdistance);
      }
  }
  return cl;      
}




AliTPCseed::AliTPCseed():AliTPCtrack(){
  //
  fRow=0; 
  fRemoval =0; 
  for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
  for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;

  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  =0;
  fRemoval = 0;
  fSort =0;
  fFirstPoint =0;
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
}

AliTPCseed::AliTPCseed(const AliTPCtrack &t):AliTPCtrack(t){
  //
  //copy constructor
  fPoints = 0;
  fEPoints = 0;
  fNShared  =0; 
  //  fTrackPoints =0;
  fRemoval =0;
  fSort =0;
  for (Int_t i=0;i<160;i++) {
    fClusterPointer[i] = 0;
    Int_t index = t.GetClusterIndex(i);
    if (index>=-1){ 
      SetClusterIndex2(i,index);
    }
    else{
      SetClusterIndex2(i,-3); 
    }    
  }
  fFirstPoint =0;
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
}

AliTPCseed::AliTPCseed(const AliKalmanTrack &t, Double_t a):AliTPCtrack(t,a){
  //
  //copy constructor
  fRow=0;
  for (Int_t i=0;i<160;i++) {
    fClusterPointer[i] = 0;
    Int_t index = t.GetClusterIndex(i);
    SetClusterIndex2(i,index);
  }
  
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0; 
  fNShared  =0; 
  //  fTrackPoints =0;
  fRemoval =0;
  fSort = 0;
  fFirstPoint =0;
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;

}

AliTPCseed::AliTPCseed(UInt_t index, const Double_t xx[5], const Double_t cc[15], 
					Double_t xr, Double_t alpha):      
  AliTPCtrack(index, xx, cc, xr, alpha) {
  //
  //
  //constructor
  fRow =0;
  for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
  for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  = 0;
  //  fTrackPoints =0;
  fRemoval =0;
  fSort =0;
  fFirstPoint =0;
  //  fHelixIn = new TClonesArray("AliHelix",0);
  //fHelixOut = new TClonesArray("AliHelix",0);
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
}

AliTPCseed::~AliTPCseed(){
  //
  // destructor
  if (fPoints) delete fPoints;
  fPoints =0;
  if (fEPoints) delete fEPoints;
  fEPoints = 0;
  fNoCluster =0;
}

AliTPCTrackerPoint * AliTPCseed::GetTrackPoint(Int_t i)
{
  //
  // 
  return &fTrackPoints[i];
}

void AliTPCseed::RebuildSeed()
{
  //
  // rebuild seed to be ready for storing
  AliTPCclusterMI cldummy;
  cldummy.SetQ(0);
  AliTPCTrackPoint pdummy;
  pdummy.GetTPoint().fIsShared = 10;
  for (Int_t i=0;i<160;i++){
    AliTPCclusterMI * cl0 = fClusterPointer[i];
    AliTPCTrackPoint *trpoint = (AliTPCTrackPoint*)fPoints->UncheckedAt(i);     
    if (cl0){
      trpoint->GetTPoint() = *(GetTrackPoint(i));
      trpoint->GetCPoint() = *cl0;
      trpoint->GetCPoint().SetQ(TMath::Abs(cl0->GetQ()));
    }
    else{
      *trpoint = pdummy;
      trpoint->GetCPoint()= cldummy;
    }
    
  }

}


Double_t AliTPCseed::GetDensityFirst(Int_t n)
{
  //
  //
  // return cluster for n rows bellow first point
  Int_t nfoundable = 1;
  Int_t nfound      = 1;
  for (Int_t i=fLastPoint-1;i>0&&nfoundable<n; i--){
    Int_t index = GetClusterIndex2(i);
    if (index!=-1) nfoundable++;
    if (index>0) nfound++;
  }
  if (nfoundable<n) return 0;
  return Double_t(nfound)/Double_t(nfoundable);

}


void AliTPCseed::GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2)
{
  // get cluster stat.  on given region
  //
  found       = 0;
  foundable   = 0;
  shared      =0;
  for (Int_t i=first;i<last; i++){
    Int_t index = GetClusterIndex2(i);
    if (index!=-1) foundable++;
    if (fClusterPointer[i]) {
      found++;
    }
    else 
      continue;

    if (fClusterPointer[i]->IsUsed(10)) {
      shared++;
      continue;
    }
    if (!plus2) continue; //take also neighborhoud
    //
    if ( (i>0) && fClusterPointer[i-1]){
      if (fClusterPointer[i-1]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    if ( fClusterPointer[i+1]){
      if (fClusterPointer[i+1]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    
  }
  //if (shared>found){
    //Error("AliTPCseed::GetClusterStatistic","problem\n");
  //}
}

//_____________________________________________________________________________
void AliTPCseed::CookdEdx(Double_t low, Double_t up,Int_t i1, Int_t i2, Bool_t onlyused) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[200];
  Float_t angular[200];
  Float_t weight[200];
  Int_t index[200];
  //Int_t nc = 0;
  //  TClonesArray & arr = *fPoints; 
  Float_t meanlog = 100.;
  
  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  for (Int_t of =0; of<4; of++){    
    for (Int_t i=of+i1;i<i2;i+=4)
      {
	Int_t index = fIndex[i];
	if (index<0||index&0x8000) continue;

	//AliTPCTrackPoint * point = (AliTPCTrackPoint *) arr.At(i);
	AliTPCTrackerPoint * point = GetTrackPoint(i);
	//AliTPCTrackerPoint * pointm = GetTrackPoint(i-1);
	//AliTPCTrackerPoint * pointp = 0;
	//if (i<159) pointp = GetTrackPoint(i+1);

	if (point==0) continue;
	AliTPCclusterMI * cl = fClusterPointer[i];
	if (cl==0) continue;	
	if (onlyused && (!cl->IsUsed(10))) continue;
	if (cl->IsUsed(11)) {
	  fNShared++;
	  continue;
	}
	Int_t   type   = cl->GetType();
	//if (point->fIsShared){
	//  fNShared++;
	//  continue;
	//}
	//if (pointm) 
	//  if (pointm->fIsShared) continue;
	//if (pointp) 
	//  if (pointp->fIsShared) continue;

	if (type<0) continue;
	//if (type>10) continue;       
	//if (point->GetErrY()==0) continue;
	//if (point->GetErrZ()==0) continue;

	//Float_t ddy = (point->GetY()-cl->GetY())/point->GetErrY();
	//Float_t ddz = (point->GetZ()-cl->GetZ())/point->GetErrZ();
	//if ((ddy*ddy+ddz*ddz)>10) continue; 


	//	if (point->GetCPoint().GetMax()<5) continue;
	if (cl->GetMax()<5) continue;
	Float_t angley = point->GetAngleY();
	Float_t anglez = point->GetAngleZ();

	Float_t rsigmay2 =  point->GetSigmaY();
	Float_t rsigmaz2 =  point->GetSigmaZ();
	/*
	Float_t ns = 1.;
	if (pointm){
	  rsigmay +=  pointm->GetTPoint().GetSigmaY();
	  rsigmaz +=  pointm->GetTPoint().GetSigmaZ();
	  ns+=1.;
	}
	if (pointp){
	  rsigmay +=  pointp->GetTPoint().GetSigmaY();
	  rsigmaz +=  pointp->GetTPoint().GetSigmaZ();
	  ns+=1.;
	}
	rsigmay/=ns;
	rsigmaz/=ns;
	*/

	Float_t rsigma = TMath::Sqrt(rsigmay2*rsigmaz2);

	Float_t ampc   = 0;     // normalization to the number of electrons
	if (i>64){
	  //	  ampc = 1.*point->GetCPoint().GetMax();
	  ampc = 1.*cl->GetMax();
	  //ampc = 1.*point->GetCPoint().GetQ();	  
	  //	  AliTPCClusterPoint & p = point->GetCPoint();
	  //	  Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.6)) - TMath::Abs(p.GetY()/0.6)+0.5);
	  // Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);
	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	  Float_t zfactor = (1.05-0.0004*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 
	}
	else{ 
	  //ampc = 1.0*point->GetCPoint().GetMax(); 
	  ampc = 1.0*cl->GetMax(); 
	  //ampc = 1.0*point->GetCPoint().GetQ(); 
	  //AliTPCClusterPoint & p = point->GetCPoint();
	  // Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.4)) - TMath::Abs(p.GetY()/0.4)+0.5);
	  //Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);

	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	Float_t zfactor = (1.02-0.000*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 

	}
	ampc *= 2.0;     // put mean value to channel 50
	//ampc *= 0.58;     // put mean value to channel 50
	Float_t w      =  1.;
	//	if (type>0)  w =  1./(type/2.-0.5); 
	//	Float_t z = TMath::Abs(cl->GetZ());
	if (i<64) {
	  ampc /= 0.6;
	  //ampc /= (1+0.0008*z);
	} else
	  if (i>128){
	    ampc /=1.5;
	    //ampc /= (1+0.0008*z);
	  }else{
	    //ampc /= (1+0.0008*z);
	  }
	
	if (type<0) {  //amp at the border - lower weight
	  // w*= 2.;
	  
	  continue;
	}
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	amp[nc[of]]        = ampc;
	angular[nc[of]]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	weight[nc[of]]     = w;
	nc[of]++;
      }
    
    TMath::Sort(nc[of],amp,index,kFALSE);
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;
    //meanlog = amp[index[Int_t(nc[of]*0.33)]];
    meanlog = 50;
    for (Int_t i=int(nc[of]*low+0.5);i<int(nc[of]*up+0.5);i++){
      Float_t ampl      = amp[index[i]]/angular[index[i]];
      ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
      //
      sumw    += weight[index[i]]; 
      sumamp  += weight[index[i]]*ampl;
      sumamp2 += weight[index[i]]*ampl*ampl;
      norm[of]    += angular[index[i]]*weight[index[i]];
    }
    if (sumw<1){ 
      SetdEdx(0);  
    }
    else {
      norm[of] /= sumw;
      mean[of]  = sumamp/sumw;
      sigma[of] = sumamp2/sumw-mean[of]*mean[of];
      if (sigma[of]>0.1) 
	sigma[of] = TMath::Sqrt(sigma[of]);
      else
	sigma[of] = 1000;
      
    mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
    //mean  *=(1-0.02*(sigma/(mean*0.17)-1.));
    //mean *=(1-0.1*(norm-1.));
    }
  }

  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //  mean[0]*= (1-0.05*(sigma[0]/(0.01+mean[1]*0.18)-1));
  //  mean[1]*= (1-0.05*(sigma[1]/(0.01+mean[0]*0.18)-1));

  
  //  dedx = (mean[0]* TMath::Sqrt((1.+nc[0]))+ mean[1]* TMath::Sqrt((1.+nc[1])) )/ 
  //  (  TMath::Sqrt((1.+nc[0]))+TMath::Sqrt((1.+nc[1])));

  Int_t norm2 = 0;
  Int_t norm3 = 0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i];
      fSdEdx    += sigma[i]*(nc[i]-2);
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i];
      norm3     += nc[i]-2;
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return;
  }
  //  Float_t dedx1 =dedx;
  /*
  dedx =0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      mean[i]   = mean[i]*(1-0.12*(sigma[i]/(fSdEdx)-1.));
      dedx      += mean[i] *nc[i];
    }
    fDEDX[i]  = mean[i];                
  }
  dedx /= norm2;
  */

  
  SetdEdx(dedx);
    
  //mi deDX



  //Very rough PID
  Double_t p=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt()));

  if (p<0.6) {
    if (dedx < 39.+ 12./(p+0.25)/(p+0.25)) { SetMass(0.13957); return;}
    if (dedx < 39.+ 12./p/p) { SetMass(0.49368); return;}
    SetMass(0.93827); return;
  }

  if (p<1.2) {
    if (dedx < 39.+ 12./(p+0.25)/(p+0.25)) { SetMass(0.13957); return;}
    SetMass(0.93827); return;
  }

  SetMass(0.13957); return;

}



/*



void AliTPCseed::CookdEdx2(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[200];
  Float_t angular[200];
  Float_t weight[200];
  Int_t index[200];
  Bool_t inlimit[200];
  for (Int_t i=0;i<200;i++) inlimit[i]=kFALSE;
  for (Int_t i=0;i<200;i++) amp[i]=10000;
  for (Int_t i=0;i<200;i++) angular[i]= 1;;
  

  //
  Float_t meanlog = 100.;
  Int_t indexde[4]={0,64,128,160};

  Float_t amean     =0;
  Float_t asigma    =0;
  Float_t anc       =0;
  Float_t anorm     =0;

  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  //  for (Int_t of =0; of<3; of++){    
  //  for (Int_t i=indexde[of];i<indexde[of+1];i++)
  for (Int_t i =0; i<160;i++)
    {
	AliTPCTrackPoint * point = GetTrackPoint(i);
	if (point==0) continue;
	if (point->fIsShared){
	  fNShared++;	  
	  continue;
	}
	Int_t   type   = point->GetCPoint().GetType();
	if (type<0) continue;
	if (point->GetCPoint().GetMax()<5) continue;
	Float_t angley = point->GetTPoint().GetAngleY();
	Float_t anglez = point->GetTPoint().GetAngleZ();
	Float_t rsigmay =  point->GetCPoint().GetSigmaY();
	Float_t rsigmaz =  point->GetCPoint().GetSigmaZ();
	Float_t rsigma = TMath::Sqrt(rsigmay*rsigmaz);

	Float_t ampc   = 0;     // normalization to the number of electrons
	if (i>64){
	  ampc =  point->GetCPoint().GetMax();
	}
	else{ 
	  ampc = point->GetCPoint().GetMax(); 
	}
	ampc *= 2.0;     // put mean value to channel 50
	//	ampc *= 0.565;     // put mean value to channel 50

	Float_t w      =  1.;
	Float_t z = TMath::Abs(point->GetCPoint().GetZ());
	if (i<64) {
	  ampc /= 0.63;
	} else
	  if (i>128){
	    ampc /=1.51;
	  }	  	
	if (type<0) {  //amp at the border - lower weight	  	  
	  continue;
	}
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	angular[i]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	amp[i]        = ampc/angular[i];
	weight[i]     = w;
	anc++;
    }

  TMath::Sort(159,amp,index,kFALSE);
  for (Int_t i=int(anc*low+0.5);i<int(anc*up+0.5);i++){      
    inlimit[index[i]] = kTRUE;  // take all clusters
  }
  
  //  meanlog = amp[index[Int_t(anc*0.3)]];
  meanlog =10000.;
  for (Int_t of =0; of<3; of++){    
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;    
   for (Int_t i=indexde[of];i<indexde[of+1];i++)
      {
	if (inlimit[i]==kFALSE) continue;
	Float_t ampl      = amp[i];
	///angular[i];
	ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
	//
	sumw    += weight[i]; 
	sumamp  += weight[i]*ampl;
	sumamp2 += weight[i]*ampl*ampl;
	norm[of]    += angular[i]*weight[i];
	nc[of]++;
      }
   if (sumw<1){ 
     SetdEdx(0);  
   }
   else {
     norm[of] /= sumw;
     mean[of]  = sumamp/sumw;
     sigma[of] = sumamp2/sumw-mean[of]*mean[of];
     if (sigma[of]>0.1) 
       sigma[of] = TMath::Sqrt(sigma[of]);
     else
       sigma[of] = 1000;      
     mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
   }
  }
    
  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //
  Int_t norm2 = 0;
  Int_t norm3 = 0;
  Float_t www[3] = {12.,14.,17.};
  //Float_t www[3] = {1.,1.,1.};

  for (Int_t i =0;i<3;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i]*www[i]/sigma[i];
      fSdEdx    += sigma[i]*(nc[i]-2)*www[i]/sigma[i];
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i]*www[i]/sigma[i];
      norm3     += (nc[i]-2)*www[i]/sigma[i];
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return;
  }
  //  Float_t dedx1 =dedx;
  
  dedx =0;
  Float_t norm4 = 0;
  for (Int_t i =0;i<3;i++){
    if (nc[i]>2&&nc[i]<1000&&sigma[i]>3){
      //mean[i]   = mean[i]*(1+0.08*(sigma[i]/(fSdEdx)-1.));
      dedx      += mean[i] *(nc[i])/(sigma[i]);
      norm4     += (nc[i])/(sigma[i]);
    }
    fDEDX[i]  = mean[i];                
  }
  if (norm4>0) dedx /= norm4;
  

  
  SetdEdx(dedx);
    
  //mi deDX

}

*/
