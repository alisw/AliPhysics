#ifndef ALIITSRAWCLUSTERSPD_H
#define ALIITSRAWCLUSTERSPD_H

#include "AliITSRawCluster.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SPD                          //
//                                                //
////////////////////////////////////////////////////

class AliITSRawClusterSPD : public AliITSRawCluster {
 public:
    AliITSRawClusterSPD();
    AliITSRawClusterSPD(Float_t clz,Float_t clx,Float_t Charge,
			Int_t ClusterSizeZ,Int_t ClusterSizeX,
			Int_t xstart,Int_t xstop,Float_t zstart,
			Float_t zstop,Int_t zend,Int_t module);
    virtual ~AliITSRawClusterSPD() {// destructor
    }
    void Add(AliITSRawClusterSPD* clJ); 
    Bool_t Brother(AliITSRawClusterSPD* cluster,Float_t dz,Float_t dx) const;
    void PrintInfo() const;
    // Getters
    Float_t Q() const {// Q
	return fQ ;}
    Float_t Z() const {// Z
	return fZ ;}
    Float_t X() const {// X
	return fX ;}
    Int_t NclZ() const {// NclZ
	return fNClZ ;}
    Int_t NclX() const {// NclX
	return fNClX ;}
    Int_t   XStart() const {//XStart
	return fXStart;}
    Int_t   XStop() const {//XStop
	return fXStop;}
    Int_t   XStartf() const {//XStartf
	return fXStart;}
    Int_t   XStopf() const {//XStopf
	return fXStop;}
    Float_t ZStart() const {//ZStart
	return fZStart;}
    Float_t ZStop() const {//ZStop
	return fZStop;}
    Int_t   Zend() const {//Zend
	return fZend;}
    Int_t   NTracks() const {//NTracks
	return fNTracks;}
    Int_t Module() const {//Returns module where this cluster came from
	return fModule;}
    void GetTracks(Int_t &track0,Int_t &track1,Int_t &track2) const {track0=fTracks[0]; track1=fTracks[1]; track2=fTracks[2];}
    void   SetTracks(Int_t track0, Int_t track1, Int_t track2);
    void   SetNTracks(Int_t ntracks) {
	// set ntracks
	fNTracks=ntracks;
    }
 protected:
    Float_t   fX;           // X of cluster
    Float_t   fZ;           // Z of cluster
    Float_t   fQ;           // Q of cluster
    Int_t     fNClZ;        // Cluster size in Z direction
    Int_t     fNClX;        // Cluster size in X direction
    Int_t     fXStart;      // number of first pixel in cluster
    Int_t     fXStop;       // number of last pixel in cluster
    Float_t   fZStart;      // number of first pixel in cluster
    Float_t   fZStop;       // number of last pixel in cluster
    Int_t     fZend;        // Zend
    Int_t     fNTracks;     // number of tracks created a cluster
    Int_t     fTracks[3];   // tracks created a cluster
    Int_t     fModule;      // Module number for this culuster
  
  ClassDef(AliITSRawClusterSPD,2)  // AliITSRawCluster class for SPD
};

#endif
