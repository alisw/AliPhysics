#ifndef ALIITSRECPOINTU_H
#define ALIITSRECPOINTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////////////////////
//  Class to set the local coordinates in ITS Upgrade recpoint //
/////////////////////////////////////////////////////////////////

#include <AliCluster.h>


class AliITSRecPointU : public AliCluster {

public :
  AliITSRecPointU();
  virtual ~AliITSRecPointU() {}; // distructor
  AliITSRecPointU(const AliITSRecPointU& pt);
  AliITSRecPointU& operator=(const AliITSRecPointU &source);

  void SetLocalCoord(Float_t x, Float_t z) {fXloc=x; fZloc=z;}
  void SetModule(Int_t i){fModule=i;} 
  void SetNTracksIdMC(Int_t nLabels) {fNTracksIdMC = nLabels;}
  void AddTrackID(Int_t tid); 

  Int_t GetModule() const {return fModule;}
  Int_t GetNTracksIdMC() const {return fNTracksIdMC;}
  Int_t GetTrackID(Int_t ipart) const {if(ipart<0 || ipart >=kMaxLab) return -1; else return fTrackIdMC[ipart];}

  enum {kMaxLab=24}; // maximum number of MC labels associated to the cluster
  void CleanLabels() {SetNTracksIdMC(0); for(Int_t i=0; i<kMaxLab ; i++) fTrackIdMC[i]=-3; }

  virtual void Print(Option_t* option = "") const;

  Float_t GetDetLocalX() const {return fXloc;} // gets fX
  Float_t GetDetLocalZ() const {return fZloc;} // gets fZ

 protected:

  Float_t   fXloc ;        //X of cluster (local coordinates)
  Float_t   fZloc ;        //Z of cluster (local coordinates)

  Int_t fModule;         // segmentation element within the same layer
  Int_t fNTracksIdMC;     // total number of associated MC labels (could be more than 3!)
  Int_t fTrackIdMC[kMaxLab];  // MC track labels 

  ClassDef(AliITSRecPointU,2)  // AliITSRecPointU class

};
#endif
