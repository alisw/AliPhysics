#include <iostream.h>
#include <TMath.h>
 
#include "AliITSRawCluster.h"

ClassImp(AliITSRawCluster)
 
ClassImp(AliITSRawClusterSDD)
//--------------------------------------
AliITSRawClusterSDD::AliITSRawClusterSDD(Int_t wing, Float_t Anode,Float_t Time,Float_t Charge,Float_t PeakAmplitude,Float_t Asigma, Float_t Tsigma,Float_t DriftPath,Float_t AnodeOffset,Int_t Samples) {
  // constructor
  fWing = wing;
  fAnode = Anode;
  fTime = Time;
  fQ = Charge;
  fPeakAmplitude = PeakAmplitude;
  fNanodes = 1;
  fNsamples = Samples;
  Int_t sign = 1;
  Int_t i;
  for(i=0;i<fWing; i++) sign*=(-1);
  fX = DriftPath*sign/10000.;
  fZ = AnodeOffset/10000.;
}

//----------------------------------------
void AliITSRawClusterSDD::Add(AliITSRawClusterSDD* clJ) {
  // add
  fAnode = (fAnode*fQ + clJ->A()*clJ->Q())/(fQ+clJ->Q());
  fTime = (fTime*fQ + clJ->T()*clJ->Q())/(fQ+clJ->Q());
  fX = (fX*fQ + clJ->X()*clJ->Q())/(fQ+clJ->Q());
  fZ = (fZ*fQ + clJ->Z()*clJ->Q())/(fQ+clJ->Q());
  fQ += clJ->Q();
  fNsamples += (Int_t) (clJ->Samples());
  (fNanodes)++;
  if(clJ->PeakAmpl() > fPeakAmplitude) fPeakAmplitude = clJ->PeakAmpl();
  
  return;
} 
//--------------------------------------
Bool_t AliITSRawClusterSDD::Brother(AliITSRawClusterSDD* cluster,Float_t danode,Float_t dtime) {
  // brother
  Bool_t brother = kTRUE;
  if(fWing != cluster->W()) return brother = kFALSE;
  if(TMath::Abs(fTime-cluster->T()) > dtime) return brother = kFALSE;
  if(TMath::Abs(fAnode-cluster->A()) > danode) return brother = kFALSE;
  return brother;
}

//--------------------------------------
void AliITSRawClusterSDD::Print() {
  // print
  cout << ", Anode " << fAnode << ", Time: " << fTime << ", Charge: " << fQ;
  cout << ", Samples: " << fNsamples;
  cout << ", X: " << fX << ", Z: " << fZ << endl;
}
//--------------------------------------


ClassImp(AliITSRawClusterSPD)
  //--------------------------------------

  AliITSRawClusterSPD::AliITSRawClusterSPD(Float_t clz,Float_t clx,Float_t Charge,Int_t ClusterSizeZ,Int_t ClusterSizeX,Int_t xstart,Int_t xstop,Int_t xstartf,Int_t xstopf,Float_t zstart,Float_t zstop,Int_t zend) {
  // constructor
  
  fZ = clz;
  fX = clx;
  fQ = Charge;
  fNClZ = ClusterSizeZ;
  fNClX = ClusterSizeX;
  fXStart = xstart;
  fXStop = xstop;
  fXStartf = xstartf;
  fXStopf = xstopf;
  fZStart = zstart;
  fZStop = zstop;
  fZend = zend;
}

//--------------------------------------
void AliITSRawClusterSPD::Add(AliITSRawClusterSPD* clJ) {
  // Recolculate the new center of gravity coordinate and cluster sizes
  // in both directions after grouping of clusters
  
  if(this->fZStop < clJ->ZStop()) this->fZStop = clJ->ZStop();  
  
  this->fZ = (this->fZ + clJ->Z())/2.;
  this->fX = (this->fX + clJ->X())/2.;
  this->fQ = this->fQ + clJ->Q();
  
  this->fXStart = clJ->XStart(); // for a comparison with the next
  this->fXStop = clJ->XStop();   // z column
  
  if(this->fXStartf > clJ->XStartf()) this->fXStartf = clJ->XStartf();
  if(this->fXStopf < clJ->XStopf()) this->fXStopf = clJ->XStopf();
  if(this->fZend < clJ->Zend()) this->fZend = clJ->Zend();
  this->fNClX = this->fXStopf - this->fXStartf + 1; 
  (this->fNClZ)++;
  
  return;
} 

//--------------------------------------
Bool_t AliITSRawClusterSPD::Brother(AliITSRawClusterSPD* cluster,Float_t dz,Float_t dx) {
  // fXStart, fXstop and fZend information is used now instead of dz and dx
  // to check an absent (or a present) of the gap between two pixels in 
  // both x and z directions. The increasing order of fZend is used.
  
  Bool_t brother = kFALSE;  
  Bool_t test2 = kFALSE;  
  Bool_t test3 = kFALSE;  
  
  // Diagonal clusters are included:
  if(fXStop >= (cluster->XStart() -1) && fXStart <= (cluster->XStop()+1)) test2 = kTRUE;
  
  // Diagonal clusters are excluded:   
  //       if(fXStop >= cluster->XStart() && fXStart <= cluster->XStop()) test2 = kTRUE;
  if(cluster->Zend() == (fZend + 1)) test3 = kTRUE; 
  if(test2 && test3) {
    // cout<<"test 2,3 0k, brother = true "<<endl;
    return brother = kTRUE;
  }
  return brother;
}

//--------------------------------------
void AliITSRawClusterSPD::Print() 
{
  // print
  cout << ", Z: " << fZ << ", X: " << fX << ", Charge: " << fQ<<endl;
  cout << " Z cluster size: " << fNClZ <<", X cluster size "<< fNClX <<endl;
  cout <<" XStart, XStop, XStartf,XStopf,Zend ="<<fXStart<<","<<fXStop<<","<<fXStartf<<","<<fXStopf<<","<<fZend<<endl;
  
}


ClassImp(AliITSRawClusterSSD)
  //--------------------------------------
  AliITSRawClusterSSD::AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn) {  
  // constructor
  //fProbability   = Prob;
  fMultiplicity  = Sp;
  fMultiplicityN = Sn;
  
}
