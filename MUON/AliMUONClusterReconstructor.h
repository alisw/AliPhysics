#ifndef ALIMUONCLUSTERRECONSTRUCTOR_H
#define ALIMUONCLUSTERRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONClusterReconstructor
/// \brief MUON cluster reconstructor in ALICE
///
/////////////////////////////////////
/// MUON event reconstructor in ALICE
/////////////////////////////////////

#include <TObject.h>
#include "AliMUONClusterFinderVS.h" //AZ

class AliLoader;
class AliMUON;
class AliMUONRawCluster;
//AZ class AliMUONClusterFinderVS;
class AliMUONData;
class AliRawReader;

class AliMUONClusterReconstructor : public TObject 
{
 public:
  AliMUONClusterReconstructor(AliLoader* loader, AliMUONData* data = 0x0); // Constructor
  virtual ~AliMUONClusterReconstructor(void); // Destructor

 
  // Cluster Finding & Trigger
  virtual void   Digits2Clusters();
  virtual void   Trigger2Trigger() ;

  // pointer to data container
  AliMUONData*   GetMUONData() {return fMUONData;}
  // Reco Model
  AliMUONClusterFinderVS* GetRecoModel() {return fRecModel;}
  //  AliMUONClusterFinderAZ* GetRecoModel() {return fRecModel;}
  //AZ void   SetRecoModel(AliMUONClusterFinderVS* rec) {fRecModel = rec;}
  void   SetRecoModel(AliMUONClusterFinderVS* rec) {if (fRecModel) delete fRecModel; fRecModel = rec;} //AZ
  //  void   SetRecoModel(AliMUONClusterFinderAZ* rec) {fRecModel = rec;}


 protected:
  AliMUONClusterReconstructor();                  // Default constructor
  AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs); // copy constructor
  AliMUONClusterReconstructor& operator=(const AliMUONClusterReconstructor& rhs); // assignment operator

 private:

  AliMUONData*            fMUONData;           //! Data container for MUON subsystem 
  AliMUONClusterFinderVS* fRecModel;           //! cluster recontruction model
  //AliMUONClusterFinderAZ* fRecModel;           //! cluster recontruction model

  // alice loader
  AliLoader* fLoader;


  ClassDef(AliMUONClusterReconstructor,0) // MUON cluster reconstructor in ALICE
};
	
#endif
