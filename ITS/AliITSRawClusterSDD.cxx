/**************************************************************************
 * Copyright(c) 2000-2004, ALICE Experiment at CERN, All rights reserved. *
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
#include <Riostream.h>
#include <TMath.h>
 
#include "AliITSRawClusterSDD.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SDD                          //
//                                                //
////////////////////////////////////////////////////

ClassImp(AliITSRawClusterSDD)
//______________________________________________________________________
AliITSRawClusterSDD::AliITSRawClusterSDD():
fX(0),
fZ(0),
fQ(0),
fWing(0),
fAnode(0),
fTime(0),
fAsigma(0),
fTsigma(0),
fPeakAmplitude(0),
fSumAmplitude(0),
fPeakPosition(-1),
fNanodes(1),
fTstart(0),
fTstop(0),
fTstartf(0),
fTstopf(0),
fAstart(0),
fAstop(0)
{
	// default constructor
  fMultiplicity = 0;
}

//______________________________________________________________________
AliITSRawClusterSDD::AliITSRawClusterSDD(Int_t wing,
					 Float_t anode,Float_t time,
					 Float_t charge,Float_t peakAmplitude,
					 Int_t peakPosition,
					 Float_t asigma,Float_t tsigma,
					 Float_t driftPath,
					 Float_t anodeOffset,
					 Int_t samples,Int_t tstart,
					 Int_t tstop,Int_t tstartf,
					 Int_t tstopf,Int_t anodes, 
					 Int_t astart, Int_t astop):
fX(0),
fZ(0),
fQ(charge),
fWing(wing),
fAnode(anode),
fTime(time),
fAsigma(asigma),
fTsigma(tsigma),
fPeakAmplitude(peakAmplitude),
fSumAmplitude(0),
fPeakPosition(peakPosition),
fNanodes(anodes),
fTstart(tstart),
fTstop(tstop),
fTstartf(tstartf),
fTstopf(tstopf),
fAstart(astart),
fAstop(astop){
    // constructor

    Int_t sign = 1;
    for(Int_t i=0;i<fWing; i++) sign *= (-1);
    fX = driftPath*sign/10000.;
    fZ = anodeOffset/10000.;
    fMultiplicity = samples;
}
//______________________________________________________________________
AliITSRawClusterSDD::AliITSRawClusterSDD(const AliITSRawClusterSDD & source):
    AliITSRawCluster(source),
fX(source.fX),
fZ(source.fZ),
fQ(source.fQ),
fWing(source.fWing),
fAnode(source.fAnode),
fTime(source.fTime),
fAsigma(source.fAsigma),
fTsigma(source.fTsigma),
fPeakAmplitude(source.fPeakAmplitude),
fSumAmplitude(source.fSumAmplitude),
fPeakPosition(source.fPeakPosition),
fNanodes(source.fNanodes),
fTstart(source.fTstart),
fTstop(source.fTstop),
fTstartf(source.fTstartf),
fTstopf(source.fTstopf),
fAstart(source.fAstart),
fAstop(source.fAstop){
    // copy constructor

    fMultiplicity  = source.fMultiplicity;
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
  // compare this with "cluster"
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
void AliITSRawClusterSDD::PrintInfo() const {
    // print

    cout << ", Anode " << fAnode << ", Time: " << fTime << ", Charge: " << fQ;
    cout << ", Samples: " << fMultiplicity;
    cout << ", X: " << fX << ", Z: " << fZ << "tstart " << fTstart 
	 << "tstop "<< fTstop <<endl;
}
