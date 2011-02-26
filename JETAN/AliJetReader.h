#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          Magali Estienne <magali.estienne@subatech.in2p3.fr>  

#include <TObject.h>
#include <TChain.h>
#include <TArrayI.h>

class TTree;
class TChain;
class TTask;
class TClonesArray;
class TRefArray;
class AliEMCALGeoUtils;
class AliJetReaderHeader;
class AliESDEvent;
class AliHeader;
class AliJetUnitArray;
class AliJetHadronCorrection;
class AliJet;
class AliJetFillUnitArray;
class AliOADBContainer;


class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  virtual TClonesArray*       GetMomentumArray()     const {return fMomentumArray;}
  virtual TRefArray*          GetReferences()        const {return 0;}   
  virtual TClonesArray        *GetUnitArray() const {return fUnitArray;}  
  virtual AliJetReaderHeader* GetReaderHeader()      const {return fReaderHeader;}
  virtual AliHeader           *GetAliHeader() const  {return fAliHeader;}
  virtual Int_t               GetSignalFlag(Int_t i) const {return fSignalFlag[i];}
  virtual Int_t               GetCutFlag(Int_t i)    const {return fCutFlag[i];}
  virtual Int_t               GetArrayInitialised()  const {return fArrayInitialised;}
  virtual Int_t               GetNumCandidate()      const {return fNumCandidate;} 
  virtual Int_t               GetNumCandidateCut()   const {return fNumCandidateCut;}
  
  // Setters
  virtual Bool_t FillMomentumArray() {return kTRUE;}
  virtual Bool_t ReadEventLoader(Int_t) {return kTRUE;}
  virtual void   InitUnitArray() {}
  virtual void   InitParameters() {}
  virtual void   CreateTasks(TChain* /*tree*/) {}
  virtual Bool_t ExecTasks(Bool_t /*procid*/, TRefArray* /*refArray*/) {return kFALSE;}
  // Correction of hadronic energy 
  virtual void   SetHadronCorrector(AliJetHadronCorrection*) {;} 
  virtual void   SetApplyMIPCorrection(Bool_t /*val*/){;}
  virtual void   SetApplyFractionHadronicCorrection(Bool_t /*val*/){;}
  virtual void   SetFractionHadronicCorrection(Double_t /*val*/){;}
  virtual void   SetApplyElectronCorrection(Int_t /*flag*/) {;}
  virtual void   SetReaderHeader(AliJetReaderHeader* header) 
  {fReaderHeader = header;}
  virtual void   SetESD(AliESDEvent* esd) { fESD = esd;}

  // Others
  virtual void   OpenInputFiles() {}
  virtual void   SetInputEvent(const TObject* /*esd*/, const TObject* /*aod*/, const TObject* /*mc*/) {;}
  virtual void   ConnectTree(TTree* /*tree*/, TObject* /*data*/) {}
  virtual Bool_t GetGenJets(AliJet* /*genJets*/) {return kFALSE;}
  
  void ClearArray();
  
  virtual const TString GetJetanOADBPath()  {return fJetanOADBpath.Data();}
  void SetJetanOADBPath(TString name) {fJetanOADBpath = name;}
 
  virtual void SetDebug(Int_t debug = 0) {fDebug = debug;}
  
 protected:
  AliJetReader(const AliJetReader& rJetReader);
  AliJetReader& operator = (const AliJetReader& rhsr);
  Bool_t SetEMCALGeometry();
  

  TString                         fJetanOADBpath;          //! path to official OADB, to be set by the task
  static AliEMCALGeoUtils         *fGeom;                  //! EMCAL Geometry 
  TChain                          *fChain;                 // chain for reconstructed tracks
  TChain                          *fTree;                  // tree for reconstructed tracks
  TClonesArray                    *fMomentumArray;         // array of particle momenta
  TClonesArray                    *fArrayMC;               //! array of mc particles
  TTask                           *fFillUnitArray;         //! task list for filling the UnitArray
  AliESDEvent                     *fESD;                   // pointer to esd
  AliJetReaderHeader              *fReaderHeader;          // pointer to header
  AliHeader                       *fAliHeader;             // AliHeader
  TArrayI                          fSignalFlag;            // to flag if a particle comes from pythia or
                                                           // from the underlying event
  TArrayI                          fCutFlag;               // to flag if a particle passed the pt cut or not
  TClonesArray                    *fUnitArray;             // array of digit position and energy 
  Bool_t                           fArrayInitialised;      // To check that array of units is initialised  
  AliJetFillUnitArray             *fFillUAFromTracks;      // For charged particle task
  AliJetFillUnitArray             *fFillUAFromEMCalDigits; // For neutral particle task
  Int_t                            fNumCandidate;          // Number of entries different from zero in unitarray
  Int_t                            fNumCandidateCut;       // Number of entries different from zero in unitarray
                                                           // which pass pt cut
  AliJetHadronCorrection          *fHadronCorrector;       //! Pointer to hadronic correction 
  Int_t                            fHCorrection;           //  Hadron correction flag 
  Int_t                            fECorrection;           //  Electron correction flag 
  Bool_t                           fEFlag;                 //  Electron correction flag 
  Int_t                            fDebug;                //! Debug option

  ClassDef(AliJetReader,1)
};
 
#endif
