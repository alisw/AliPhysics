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

class AliLoader;
class AliMUON;
class AliMUONChamber;
class AliMUONRawCluster;
class AliMUONClusterFinderVS;
class AliMUONData;


class AliMUONClusterReconstructor : public TObject 
{
 public:
  AliMUONClusterReconstructor(AliLoader* loader); // Constructor
  virtual ~AliMUONClusterReconstructor(void); // Destructor

  // Interface with AliMUONData
  virtual void       SetTreeAddress(){};
    
  // Cluster Finding & Trigger
  virtual void   Digits2Clusters();


  // void EventDump(void);  // dump reconstructed event
  
  // Set Reconstruction Model
  virtual void   SetReconstructionModel(Int_t id, AliMUONClusterFinderVS* reconst);
 
  AliMUONData*   GetMUONData() {return fMUONData;}

  Int_t GetPrintLevel(void) const {return fPrintLevel;}
  void SetPrintLevel(Int_t printLevel) {fPrintLevel = printLevel;}

 protected:
  AliMUONClusterReconstructor();                  // Default constructor
  AliMUONClusterReconstructor (const AliMUONClusterReconstructor& rhs); // copy constructor
  AliMUONClusterReconstructor& operator=(const AliMUONClusterReconstructor& rhs); // assignment operator

 private:
  static const Int_t fgkDefaultPrintLevel;     // Default print level

  Int_t                   fNCh;                // Number of chambers   
  Int_t                   fNTrackingCh;        // Number of tracking chambers*
  AliMUONData*            fMUONData;           //! Data container for MUON subsystem 
  AliMUON*                fMUON;               //! pointer to MUON  
  TObjArray*              fChambers;           //! List of Tracking Chambers

 // print level
  Int_t fPrintLevel;

  // debug
  Int_t fDebug;
  
  // alice loader
  AliLoader* fLoader;


  ClassDef(AliMUONClusterReconstructor,0) // MUON cluster reconstructor in ALICE
};
	
#endif
