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
 
#include "AliITSRawClusterSPD.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SPD                          //
//                                                //
////////////////////////////////////////////////////

ClassImp(AliITSRawClusterSPD)

//______________________________________________________________________
AliITSRawClusterSPD::AliITSRawClusterSPD():
fX(0.0),
fZ(0.0),
fQ(0.0),
fNClZ(0),
fNClX(0),
fXStart(0),
fXStop(0),
fZStart(0.0),
fZStop(0.0),
fZend(0),
fNTracks(0),
fModule(0)
{
  // Default constructor
  fTracks[0]=fTracks[1]=fTracks[2]=-3;
}

//______________________________________________________________________
AliITSRawClusterSPD::AliITSRawClusterSPD(Float_t clz,Float_t clx,
					 Float_t Charge,Int_t ClusterSizeZ,
					 Int_t ClusterSizeX,Int_t xstart,
					 Int_t xstop,
					 Float_t zstart,Float_t zstop,
					 Int_t zend,Int_t module):
fX(clx),
fZ(clz),
fQ(Charge),
fNClZ(ClusterSizeZ),
fNClX(ClusterSizeX),
fXStart(xstart),
fXStop(xstop),
fZStart(zstart),
fZStop(zstop),
fZend(zend),
fNTracks(0),
fModule(module) {
    // constructor

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
				    Float_t dz,Float_t dx) const {
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
void AliITSRawClusterSPD::PrintInfo() const{
    // print

    cout << ", Z: " << fZ << ", X: " << fX << ", Charge: " << fQ<<endl;
    cout << " Z cluster size: " << fNClZ <<", X cluster size "<< fNClX <<endl;
    cout <<" XStart, XStop,Zend, Module ="<<fXStart<<","
	 <<fXStop<<","<<fZend << "," << fModule<<endl;
}

//______________________________________________________________________
void AliITSRawClusterSPD::SetTracks(Int_t track0, Int_t track1, Int_t track2) {
  // set tracks in cluster (no more than three ones)
  fTracks[0]=track0;
  fTracks[1]=track1;
  fTracks[2]=track2;
}
