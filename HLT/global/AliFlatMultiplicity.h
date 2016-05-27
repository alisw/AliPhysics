#ifndef ALIFLATMultiplicity_H
#define ALIFLATMultiplicity_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Author : Mikolaj Krzewicki <mkrzewic@cern.ch>     */

/**
 * Flat structure representing a ESD Multiplicity 
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliMultiplicity.h"
#include <algorithm>

class AliFlatMultiplicity
{
 public:
  // -- Constructor / Destructors
  
  AliFlatMultiplicity();
 ~AliFlatMultiplicity() {}

  // constructor and method for reinitialisation of virtual table
  AliFlatMultiplicity( AliVConstructorReinitialisationFlag );
  void Reinitialize() const {} // no virtual table - do nothing

  void SetFromMultiplicity(const AliMultiplicity &v );
  void GetMultiplicity( AliMultiplicity &v ) const;
  
  void SetNumberOfTracklets(Int_t n) {fNtracks = n;}
  void SetITSClusters(const UInt_t* c) {std::copy(c,c+6,fITSClusters);}
  Int_t GetNumberOfTracklets() const {return fNtracks;}
  UInt_t GetITSClusters(Int_t layer) const {return (layer<6)?fITSClusters[layer]:0;}

  static size_t GetSize() { return sizeof(AliFlatMultiplicity); }
  
private: 
  Int_t fNtracks;            // Number of tracklets
  UInt_t fITSClusters[6];    // Number of ITS cluster per layer
};

inline AliFlatMultiplicity::AliFlatMultiplicity() :
  fNtracks(0)
{   
  // Default constructor 
  for( int i=0; i<6; i++ ) fITSClusters[i] = 0;
}

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatMultiplicity::AliFlatMultiplicity( AliVConstructorReinitialisationFlag ){}  // do nothing
#pragma GCC diagnostic warning "-Weffc++" 

inline void AliFlatMultiplicity::SetFromMultiplicity(const AliMultiplicity &v )
{
  fNtracks = v.GetNumberOfTracklets();
  fITSClusters[0] = v.GetNumberOfITSClusters(0);
  fITSClusters[1] = v.GetNumberOfITSClusters(1);
  fITSClusters[2] = v.GetNumberOfITSClusters(2);
  fITSClusters[3] = v.GetNumberOfITSClusters(3);
  fITSClusters[4] = v.GetNumberOfITSClusters(4);
  fITSClusters[5] = v.GetNumberOfITSClusters(5);
}

inline void AliFlatMultiplicity::GetMultiplicity( AliMultiplicity &v ) const
{
  v.SetNumberOfTracklets(fNtracks);
  v.SetITSClusters(0,fITSClusters[0]);
  v.SetITSClusters(1,fITSClusters[1]);
  v.SetITSClusters(2,fITSClusters[2]);
  v.SetITSClusters(3,fITSClusters[3]);
  v.SetITSClusters(4,fITSClusters[4]);
  v.SetITSClusters(5,fITSClusters[5]);
}

#endif
