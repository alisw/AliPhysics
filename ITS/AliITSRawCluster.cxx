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
 
/* $Id$ */

#include <Riostream.h>
#include <TMath.h>
 
#include "AliITSRawCluster.h"

ClassImp(AliITSRawCluster)
ClassImp(AliITSRawClusterSDD)

//______________________________________________________________________
AliITSRawClusterSDD::AliITSRawClusterSDD(Int_t wing,
					 Float_t Anode,Float_t Time,
					 Float_t Charge,Float_t PeakAmplitude,
					 Int_t PeakPosition,
					 Float_t Asigma,Float_t Tsigma,
					 Float_t DriftPath,
					 Float_t AnodeOffset,
					 Int_t Samples,Int_t Tstart,
					 Int_t Tstop,Int_t Tstartf,
					 Int_t Tstopf,Int_t Anodes, 
					 Int_t Astart, Int_t Astop){
    // constructor

    fWing          = wing;
    fAnode         = Anode;
    fTime          = Time;
    fQ             = Charge;
    fPeakAmplitude = PeakAmplitude;
    fPeakPosition  = PeakPosition;
    fAsigma        = Asigma;
    fTsigma        = Tsigma;
    fNanodes       = Anodes;
    fTstart        = Tstart;
    fTstop         = Tstop;
    fTstartf       = Tstartf;
    fTstopf        = Tstopf;
    fAstart        = Astart;
    fAstop         = Astop;
    fMultiplicity  = Samples;
    fSumAmplitude  = 0;

    Int_t sign = 1;
    for(Int_t i=0;i<fWing; i++) sign *= (-1);
    fX = DriftPath*sign/10000.;
    fZ = AnodeOffset/10000.;
}
//______________________________________________________________________
AliITSRawClusterSDD::AliITSRawClusterSDD(const AliITSRawClusterSDD & source):
    AliITSRawCluster(source){
    // copy constructor

    fWing          = source.fWing;
    fAnode         = source.fAnode;
    fTime          = source.fTime;
    fQ             = source.fQ;
    fPeakAmplitude = source.fPeakAmplitude;
    fPeakPosition  = source.fPeakPosition;
    fAsigma        = source.fAsigma;
    fTsigma        = source.fTsigma;
    fNanodes       = source.fNanodes;
    fTstart        = source.fTstart;
    fTstop         = source.fTstop;
    fTstartf       = source.fTstartf;
    fTstopf        = source.fTstopf;
    fAstart        = source.fAstart;
    fAstop         = source.fAstop;

    fMultiplicity  = source.fMultiplicity;
    fSumAmplitude  = source.fSumAmplitude;
    fX             = source.fX;
    fZ             = source.fZ;
}
//______________________________________________________________________
void AliITSRawClusterSDD::Add(AliITSRawClusterSDD* clJ) {
    // add

    fAnode = (fAnode*fQ + clJ->A()*clJ->Q())/(fQ+clJ->Q());
    fTime  = ( fTime*fQ + clJ->T()*clJ->Q())/(fQ+clJ->Q());
    fX     = (    fX*fQ + clJ->X()*clJ->Q())/(fQ+clJ->Q());
    fZ     = (    fZ*fQ + clJ->Z()*clJ->Q())/(fQ+clJ->Q());
    fQ += clJ->Q();
    if(fSumAmplitude == 0) fSumAmplitude += fPeakAmplitude;
    /*
      fAnode = (fAnode*fSumAmplitude+clJ->A()*clJ->PeakAmpl())/
               (fSumAmplitude+clJ->PeakAmpl());
      fTime = (fTime*fSumAmplitude +clJ->T()*clJ->PeakAmpl())/
              (fSumAmplitude+clJ->PeakAmpl());
      fX = (fX*fSumAmplitude +clJ->X()*clJ->PeakAmpl())/
           (fSumAmplitude+clJ->PeakAmpl());
      fZ = (fZ*fSumAmplitude +clJ->Z()*clJ->PeakAmpl())/
           (fSumAmplitude+clJ->PeakAmpl());
    */
    fSumAmplitude += clJ->PeakAmpl();

    fTstart = clJ->Tstart();
    fTstop  = clJ->Tstop();
    if(fTstartf > clJ->Tstartf()) fTstartf = clJ->Tstartf();
    if( fTstopf < clJ->Tstopf() ) fTstopf  = clJ->Tstopf();
    if(  fAstop < clJ->Astop()  ) fAstop   = clJ->Astop();

    fMultiplicity += (Int_t) (clJ->Samples());
    (fNanodes)++;
    if(clJ->PeakAmpl() > fPeakAmplitude) {
	fPeakAmplitude = clJ->PeakAmpl();
	fPeakPosition = clJ->PeakPos();
    } // end if

    return;
}
//______________________________________________________________________
Bool_t AliITSRawClusterSDD::Brother(AliITSRawClusterSDD* cluster,
				    Float_t danode,Float_t dtime) {

    Bool_t brother = kFALSE;
    Bool_t test2 = kFALSE;
    Bool_t test3 = kFALSE;
    Bool_t test4 = kFALSE;
    Bool_t test5 = kFALSE;
  
    if(fWing != cluster->W()) return brother;

    if(fTstopf >= cluster->Tstart() &&
       fTstartf <= cluster->Tstop()) test2 = kTRUE;
    if(cluster->Astop() == (fAstop+1)) test3 = kTRUE;

    if(TMath::Abs(fTime-cluster->T()) < dtime) test4 = kTRUE;
    if(TMath::Abs(fAnode-cluster->A()) < danode) test5 = kTRUE;

    if((test2 && test3) || (test4 && test5) ) {
	return brother = kTRUE;
    } // end if
  
    return brother;
}
//______________________________________________________________________
void AliITSRawClusterSDD::PrintInfo() {
    // print

    cout << ", Anode " << fAnode << ", Time: " << fTime << ", Charge: " << fQ;
    cout << ", Samples: " << fMultiplicity;
    cout << ", X: " << fX << ", Z: " << fZ << "tstart " << fTstart 
	 << "tstop "<< fTstop <<endl;
}
//======================================================================
ClassImp(AliITSRawClusterSPD)
//______________________________________________________________________
AliITSRawClusterSPD::AliITSRawClusterSPD(Float_t clz,Float_t clx,
					 Float_t Charge,Int_t ClusterSizeZ,
					 Int_t ClusterSizeX,Int_t xstart,
					 Int_t xstop,
					 Float_t zstart,Float_t zstop,
					 Int_t zend,Int_t module) {
    // constructor

    fZ       = clz;
    fX       = clx;
    fQ       = Charge;
    fNClZ    = ClusterSizeZ;
    fNClX    = ClusterSizeX;
    fXStart  = xstart;
    fXStop   = xstop;
    fZStart  = zstart;
    fZStop   = zstop;
    fZend    = zend;
    fModule  = module;
}
//______________________________________________________________________
void AliITSRawClusterSPD::Add(AliITSRawClusterSPD* clJ) {
    // Recolculate the new center of gravity coordinate and cluster sizes
    // in both directions after grouping of clusters

    if(this->fZStop < clJ->ZStop()) this->fZStop = clJ->ZStop();
    this->fZ      = this->fZ + clJ->Z();
    this->fX      = (this->fX + clJ->X())/2.;
    this->fQ      = this->fQ + clJ->Q();
    this->fXStart = clJ->XStart(); // for a comparison with the next
    this->fXStop  = clJ->XStop();  // z column
    if(this->fZend < clJ->Zend())       this->fZend    = clJ->Zend();
    this->fNClX   = this->fXStop - this->fXStart + 1; 
    (this->fNClZ)++;

    return;
}
//______________________________________________________________________
Bool_t AliITSRawClusterSPD::Brother(AliITSRawClusterSPD* cluster,
				    Float_t dz,Float_t dx) {
    // fXStart, fXstop and fZend information is used now instead of dz and dx
    // to check an absent (or a present) of the gap between two pixels in 
    // both x and z directions. The increasing order of fZend is used.
    Bool_t brother = kFALSE;  
    Bool_t test2 = kFALSE;  
    Bool_t test3 = kFALSE;

    dx = dz = 0; // to remove unused variable warning.
    // Diagonal clusters are included:
    if(fXStop >= (cluster->XStart() -1) && 
       fXStart <= (cluster->XStop()+1)) test2 = kTRUE;

    // Diagonal clusters are excluded:   
    // if(fXStop >= cluster->XStart() &&
    //    fXStart <= cluster->XStop()) test2 = kTRUE;
    if(cluster->Zend() == (fZend + 1)) test3 = kTRUE; 
    if(test2 && test3) {
	// cout<<"test 2,3 0k, brother = true "<<endl;
	return brother = kTRUE;
    } // end if
    return brother;
}
//______________________________________________________________________
void AliITSRawClusterSPD::PrintInfo(){
    // print

    cout << ", Z: " << fZ << ", X: " << fX << ", Charge: " << fQ<<endl;
    cout << " Z cluster size: " << fNClZ <<", X cluster size "<< fNClX <<endl;
    cout <<" XStart, XStop,Zend, Module ="<<fXStart<<","
	 <<fXStop<<","<<fZend << "," << fModule<<endl;
}
//======================================================================
ClassImp(AliITSRawClusterSSD)
//______________________________________________________________________
AliITSRawClusterSSD::AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn) {  
    // constructor

    Prob = 0.0; // added to remove unused variable warning.
    //fProbability   = Prob;
    fMultiplicity  = Sp;
    fMultiplicityN = Sn;
}
