#ifndef ALIMUONCLUSTERRECONSTRUCTOR_H
#define ALIMUONCLUSTERRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONClusterReconstructor
/// \brief MUON cluster reconstructor in ALICE

#include <TObject.h>

class AliMUONClusterFinderVS;
class AliMUONData;
class TClonesArray;
class AliMUONVClusterFinder;
class AliMUONGeometryTransformer;

class AliMUONClusterReconstructor : public TObject 
{
 public:
  AliMUONClusterReconstructor(AliMUONData* data = 0x0,
                              AliMUONVClusterFinder* finder = 0x0,
                              const AliMUONGeometryTransformer* transformer = 0x0
                              ); 
  virtual ~AliMUONClusterReconstructor(void); // Destructor

 
  // Cluster Finding & Trigger
  virtual void   Digits2Clusters(Int_t chBeg = 0);
  virtual void   Trigger2Trigger() ;

  // Reco Model
  AliMUONClusterFinderVS* GetRecoModel() {return fRecModel;}

  void SetRecoModel(AliMUONClusterFinderVS* rec);

 protected:

  AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs); // copy constructor
  AliMUONClusterReconstructor& operator=(const AliMUONClusterReconstructor& rhs); // assignment operator

  void ClusterizeOneDE(Int_t detElemId);
  void ClusterizeOneDEV2(Int_t detElemId);
  
private:
  AliMUONVClusterFinder* fClusterFinder; //!< the object doing the real job (not owner)  
  AliMUONData*            fMUONData;           //!< Data container for MUON subsystem 
  AliMUONClusterFinderVS* fRecModel;           //!< cluster recontruction model

  TClonesArray* fDigitsCath0; //!< digits for cathode 0 of the current DE
  TClonesArray* fDigitsCath1; //!< digits for cathode 1 of the current DE
  
  const AliMUONGeometryTransformer* fTransformer; //!< to go from local to global (not owner)
    
  ClassDef(AliMUONClusterReconstructor,0) // MUON cluster reconstructor in ALICE
};
	
#endif
