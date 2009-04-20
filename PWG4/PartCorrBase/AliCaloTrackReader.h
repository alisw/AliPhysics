#ifndef ALICALOTRACKREADER_H
#define ALICALOTRACKREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors.
// Not all MC particles/tracks/clusters are kept, some kinematical restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TRefArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader: Fills Kinematics data in 3 TRefArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TRefArrays (PHOS, EMCAL, CTS) 
//                          
// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include "TObject.h" 
class TRefArray ; 
class TLorentzVector ;
class TString ;
class TRefArray;
class TArrayF;  

//--- ANALYSIS system ---
class AliStack ; 
class AliHeader ; 
class AliGenEventHeader ; 
class AliVEvent;
class AliAODEvent;  
class AliMCEvent;
class AliFidutialCut;

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

  virtual TRefArray* GetAODCTS()   const {return fAODCTS ;}
  virtual TRefArray* GetAODEMCAL() const {return fAODEMCAL ;}
  virtual TRefArray* GetAODPHOS()  const {return fAODPHOS ;}
  virtual TNamed* GetEMCALCells()  const {return fEMCALCells ;}
  virtual TNamed* GetPHOSCells()   const {return fPHOSCells ;}

  virtual AliStack*    GetStack()      const ;
  virtual AliHeader*   GetHeader()     const ;
  virtual AliGenEventHeader* GetGenEventHeader() const ;
  virtual AliVEvent*   GetInputEvent()  const {return fInputEvent;}
  virtual AliAODEvent* GetOutputEvent() const {return fOutputEvent;}
  virtual AliMCEvent*  GetMC()          const {return fMC;}
  virtual void         GetVertex(Double_t * ) const {;}

  virtual void SetInputEvent(AliVEvent* input)  {fInputEvent  = input;}
  virtual void SetOutputEvent(AliAODEvent* aod) {fOutputEvent = aod;}
  virtual void SetMC(AliMCEvent* mc)            {fMC  = mc;}

  virtual void ResetLists();

  virtual AliFidutialCut * GetFidutialCut() const {return  fFidutialCut ;}
  virtual void SetFidutialCut(AliFidutialCut * fc) { fFidutialCut = fc ;}

  virtual void SetInputOutputMCEvent(AliVEvent* /*esd*/, AliAODEvent* /*aod*/, AliMCEvent* /*mc*/) {;}

 protected:
  Int_t	           fEventNumber; // Event number
  Int_t            fDataType ;   // Select MC:Kinematics, Data:ESD/AOD, MCData:Both
  Int_t            fDebug;       // Debugging level
  AliFidutialCut * fFidutialCut; // Acceptance cuts
		
  Float_t        fCTSPtMin;      // pT Threshold on charged particles 
  Float_t        fEMCALPtMin;    // pT Threshold on emcal clusters
  Float_t        fPHOSPtMin;     // pT Threshold on phos clusters

  TRefArray *    fAODCTS ;        //! temporal referenced array with tracks
  TRefArray *    fAODEMCAL ;      //! temporal referenced array with EMCAL CaloClusters
  TRefArray *    fAODPHOS ;       //! temporal referenced array with PHOS CaloClusters
  TNamed *       fEMCALCells ;    //! temporal array with EMCAL CaloCells, ESD or AOD
  TNamed *       fPHOSCells ;     //! temporal array with PHOS CaloCells, ESD or AOD

  AliVEvent   *  fInputEvent;     //! pointer to esd or aod input
  AliAODEvent *  fOutputEvent;    //! pointer to aod output
  AliMCEvent  *  fMC;             //! Monte Carlo Event Handler  

  Bool_t         fFillCTS;        // use data from CTS
  Bool_t         fFillEMCAL;      // use data from EMCAL
  Bool_t         fFillPHOS;       // use data from PHOS
  Bool_t         fFillEMCALCells; // use data from EMCAL
  Bool_t         fFillPHOSCells;  // use data from PHOS

  ClassDef(AliCaloTrackReader,3)
} ;


#endif //ALICALOTRACKREADER_H



