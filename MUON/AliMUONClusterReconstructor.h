#ifndef ALIMUONCLUSTERRECONSTRUCTOR_H
#define ALIMUONCLUSTERRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONClusterReconstructor
/// \brief Steering class for clusterizing MUON tracking chambers

#include <TObject.h>

class AliMUONVClusterFinder;
class AliMUONGeometryTransformer;
class AliMUONVClusterStore;
class AliMUONVDigitStore;

class AliMUONClusterReconstructor : public TObject 
{
 public:
  AliMUONClusterReconstructor(AliMUONVClusterFinder* finder = 0x0,
                              const AliMUONGeometryTransformer* transformer = 0x0
                              ); 

  virtual ~AliMUONClusterReconstructor(void); // Destructor

 
  // Cluster Finding 
  virtual void Digits2Clusters(const AliMUONVDigitStore& digitStore, AliMUONVClusterStore& clusterStore);
  
 protected:
  /// Not implemented
  AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs); // copy constructor
  /// Not implemented
  AliMUONClusterReconstructor& operator=(const AliMUONClusterReconstructor& rhs); // assignment operator

  void ClusterizeOneDE(Int_t detElemId, const AliMUONVDigitStore& digitStore);

private:
  AliMUONVClusterFinder* fClusterFinder; //!< the object doing the real clustering job (not owner)  

  const AliMUONGeometryTransformer* fTransformer; //!< to go from local to global (not owner)
  
  AliMUONVClusterStore* fClusterStore; //!< not owner
  
  Int_t fNCluster; //!< number of clusters in the cluster store (used to define the cluster ID)
  
  ClassDef(AliMUONClusterReconstructor,0) // Clustering steering
};
	
#endif
