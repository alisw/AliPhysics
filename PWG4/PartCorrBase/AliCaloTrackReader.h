#ifndef ALICALOTRACKREADER_H
#define ALICALOTRACKREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors.
// Not all MC particles/tracks/clusters are kept, some kinematical restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader: Fills Kinematics data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TClonesArrays (PHOS, EMCAL, CTS) 
//                          
// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include "TObject.h" 
class TClonesArray ; 
class TLorentzVector ;
class TString ;
#include "TArrayF.h"  

//--- ANALYSIS system ---
class AliStack ; 
class AliHeader ; 
class AliGenEventHeader ; 
#include "AliESDEvent.h" 
#include "AliAODEvent.h" 
#include "AliMCEvent.h" 
class AliLog ;
#include "AliFidutialCut.h"

class AliCaloTrackReader : public TObject {

public: 

  AliCaloTrackReader() ; // ctor
  AliCaloTrackReader(const AliCaloTrackReader & g) ; // cpy ctor
  AliCaloTrackReader & operator = (const AliCaloTrackReader & g) ;//cpy assignment
  virtual ~AliCaloTrackReader() ;//virtual dtor

  enum inputDataType {kESD, kAOD, kMC};
  
  virtual void InitParameters();
  virtual void Print(const Option_t * opt) const;

  virtual Int_t GetDebug()         const { return fDebug ; }
  virtual void  SetDebug(Int_t d)        { fDebug = d ; }
  virtual Int_t GetDataType()      const { return fDataType ; }
  virtual void  SetDataType(Int_t data ) { fDataType = data ; }

  virtual Int_t GetEventNumber() const {return fEventNumber ; }

  //Minimum pt setters and getters 
  virtual Float_t  GetEMCALPtMin() const { return fEMCALPtMin  ; }
  virtual Float_t  GetPHOSPtMin()  const { return fPHOSPtMin  ; }
  virtual Float_t  GetCTSPtMin()   const { return fCTSPtMin  ; }

  virtual void SetEMCALPtMin(Float_t  pt) { fEMCALPtMin = pt ; }
  virtual void SetPHOSPtMin(Float_t  pt)  { fPHOSPtMin = pt ; }
  virtual void SetCTSPtMin(Float_t  pt)   { fCTSPtMin = pt ; }
  
  //Input setters and getters

  Bool_t IsCTSSwitchedOn()  const { return fFillCTS ; }
  void SwitchOnCTS()    {fFillCTS = kTRUE ; }
  void SwitchOffCTS()   {fFillCTS = kFALSE ; }

  Bool_t IsEMCALSwitchedOn() const { return fFillEMCAL ; }
  void SwitchOnEMCAL()  {fFillEMCAL = kTRUE ; }
  void SwitchOffEMCAL() {fFillEMCAL = kFALSE ; }

  Bool_t IsPHOSSwitchedOn()  const { return fFillPHOS ; }
  void SwitchOnPHOS()   {fFillPHOS = kTRUE ; }
  void SwitchOffPHOS()  {fFillPHOS = kFALSE ; }

  Bool_t IsEMCALCellsSwitchedOn() const { return fFillEMCALCells ; }
  void SwitchOnEMCALCells()  {fFillEMCALCells = kTRUE ; }
  void SwitchOffEMCALCells() {fFillEMCALCells = kFALSE ; }

  Bool_t IsPHOSCellsSwitchedOn()  const { return fFillPHOSCells ; }
  void SwitchOnPHOSCells()   {fFillPHOSCells = kTRUE ; }
  void SwitchOffPHOSCells()  {fFillPHOSCells = kFALSE ; }

  virtual void FillInputEvent(Int_t iEntry)  ;
  virtual void FillInputCTS()   {;}
  virtual void FillInputEMCAL() {;}
  virtual void FillInputPHOS()  {;}
  virtual void FillInputEMCALCells() {;}
  virtual void FillInputPHOSCells()  {;}

  virtual TClonesArray* GetAODCTS()   const {return fAODCTS ;}
  virtual TClonesArray* GetAODEMCAL() const {return fAODEMCAL ;}
  virtual TClonesArray* GetAODPHOS()  const {return fAODPHOS ;}
  virtual TNamed* GetEMCALCells()     const {return fEMCALCells ;}
  virtual TNamed* GetPHOSCells()      const {return fPHOSCells ;}

  virtual AliStack*          GetStack()  const ;
  virtual AliHeader*         GetHeader() const ;
  virtual AliGenEventHeader* GetGenEventHeader() const ;
  virtual AliESDEvent* GetESD() const {return fESD;}
  virtual AliAODEvent* GetAOD() const {return fAOD;}
  virtual AliMCEvent*  GetMC()  const {return fMC;}
  virtual AliVEvent*   GetInputEvent()        const {return (new AliESDEvent());}
  virtual void         GetVertex(Double_t * ) const {;}

  virtual void SetESD( AliESDEvent* esd) {fESD = esd;}
  virtual void SetAOD(AliAODEvent* aod)  {fAOD = aod;}
  virtual void SetMC(AliMCEvent* mc)     {fMC  = mc;}

  virtual void ResetLists();

  virtual AliFidutialCut * GetFidutialCut() const {return  fFidutialCut ;}
  virtual void SetFidutialCut(AliFidutialCut * fc) { fFidutialCut = fc ;}

  virtual void SetInputEvent(TObject* /*esd*/, TObject* /*aod*/, TObject* /*mc*/) {;}

 protected:
  Int_t			   fEventNumber; // Event number
  Int_t            fDataType ;   // Select MC:Kinematics, Data:ESD/AOD, MCData:Both
  Int_t            fDebug;       // Debugging level
  AliFidutialCut * fFidutialCut; // Acceptance cuts
		
  Float_t        fCTSPtMin;      // pT  Threshold on charged particles 
  Float_t        fEMCALPtMin;    // pT Threshold on emcal clusters
  Float_t        fPHOSPtMin;     // pT  Threshold on phos clusters

  TClonesArray * fAODCTS ;        //! temporal array with tracks
  TClonesArray * fAODEMCAL ;      //! temporal array with EMCAL CaloClusters
  TClonesArray * fAODPHOS ;       //! temporal array with PHOS CaloClusters
  TNamed *       fEMCALCells ;    //! temporal array with EMCAL CaloCells, ESD or AOD
  TNamed *       fPHOSCells ;     //! temporal array with PHOS CaloCells, ESD or AOD

  AliESDEvent *  fESD;            //! pointer to esd
  AliAODEvent *  fAOD;            //! pointer to aod
  AliMCEvent  *  fMC;             //! Monte Carlo Event Handler  

  Bool_t         fFillCTS;        // use data from CTS
  Bool_t         fFillEMCAL;      // use data from EMCAL
  Bool_t         fFillPHOS;       // use data from PHOS
  Bool_t         fFillEMCALCells; // use data from EMCAL
  Bool_t         fFillPHOSCells;  // use data from PHOS

  ClassDef(AliCaloTrackReader,1)
} ;


#endif //ALICALOTRACKREADER_H



