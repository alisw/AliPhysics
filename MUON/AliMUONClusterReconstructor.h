#ifndef ALIMUONCLUSTERRECONSTRUCTOR_H
#define ALIMUONCLUSTERRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

////////////////////////////////////
// MUON event reconstructor in ALICE
////////////////////////////////////

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

  // Interface with AliMUONData
  virtual void       SetTreeAddress(){};
    
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

  // print level
  Int_t GetPrintLevel(void) const {return fPrintLevel;}
  void SetPrintLevel(Int_t printLevel) {fPrintLevel = printLevel;}

 protected:
  AliMUONClusterReconstructor();                  // Default constructor
  AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs); // copy constructor
  AliMUONClusterReconstructor& operator=(const AliMUONClusterReconstructor& rhs); // assignment operator

 private:
  static const Int_t fgkDefaultPrintLevel;     // Default print level

  AliMUONData*            fMUONData;           //! Data container for MUON subsystem 
  AliMUONClusterFinderVS* fRecModel;           //! cluster recontruction model
  //AliMUONClusterFinderAZ* fRecModel;           //! cluster recontruction model

 // print level
  Int_t fPrintLevel;

  // debug
  Int_t fDebug;
  
  // alice loader
  AliLoader* fLoader;


  ClassDef(AliMUONClusterReconstructor,0) // MUON cluster reconstructor in ALICE
};
	
#endif
