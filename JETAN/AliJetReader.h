#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          Magali Estienne <magali.estienne@IReS.in2p3.fr>  

#include <TObject.h>
#include <TChain.h>
#include <TArrayI.h>

class TTree;
class TTask;
class TChain;
class TClonesArray;
class TRefArray;
class AliJetReaderHeader;
class AliESDEvent;
class AliHeader;
class AliJetUnitArray;
class AliJetHadronCorrectionv1;
class AliJet;
class AliJetFillUnitArrayTracks;
class AliJetFillUnitArrayEMCalDigits;

class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  virtual TClonesArray*        GetMomentumArray()    const {return fMomentumArray;}
  virtual TRefArray*           GetReferences()       const {return 0;}   
  virtual TClonesArray        *GetUnitArray()        const {return fUnitArray;}  
  virtual TRefArray           *GetRefArray()         const {return fRefArray;}
  virtual TClonesArray        *GetUnitArrayNoCuts()  const {return fUnitArrayNoCuts;} 

  virtual AliJetReaderHeader* GetReaderHeader()      const {return fReaderHeader;}
  virtual Int_t               GetSignalFlag(Int_t i) const {return fSignalFlag[i];}
  virtual Int_t               GetCutFlag(Int_t i)    const {return fCutFlag[i];}
  virtual Int_t               GetArrayInitialised()  const {return fArrayInitialised;}
  virtual Int_t               GetNumCandidate()      const {return fNumCandidate;} 
  virtual Int_t               GetNumCandidateCut()   const {return fNumCandidateCut;}
  
  // Setters
  virtual Bool_t FillMomentumArray(Int_t) {return kTRUE;}
  virtual void   FillUnitArrayFromTPCTracks(Int_t) {}     // temporarily not used
  virtual void   FillUnitArrayFromEMCALHits() {}          // temporarily not used
  virtual void   FillUnitArrayFromEMCALDigits(Int_t) {}   // temporarily not used
  virtual void   FillUnitArrayFromEMCALClusters(Int_t) {} // temporarily not used
  virtual void   InitUnitArray() {}
  virtual void   InitParameters() {}
  virtual void   CreateTasks() {}
  //  virtual void   ExecTasks(Int_t) {}
  virtual Bool_t   ExecTasks(Int_t) {return kTRUE;}
  /*   // Correction of hadronic energy 
       virtual void SetHadronCorrector(AliEMCALHadronCorrectionv1* corr) {fHadronCorrector = corr;} 
       virtual void SetHadronCorrection(Int_t flag = 1) {fHCorrection = flag;} */
  virtual void   SetReaderHeader(AliJetReaderHeader* header) 
      {fReaderHeader = header;}
  virtual void   SetESD(AliESDEvent* esd) { fESD = esd;}
  //  virtual Int_t  SetNumCandidate(Int_t cand) {fNumCandidate = cand;} 
  //  virtual Int_t  SetNumCandidateCut(Int_t candcut) {fNumCandidateCut = candcut;}
  

  // Others
  virtual void   OpenInputFiles() {}
  virtual void   SetInputEvent(TObject* /*esd*/, TObject* /*aod*/, TObject* /*mc*/) {;}
  virtual void   ConnectTree(TTree* /*tree*/, TObject* /*data*/) {}
  virtual Bool_t GetGenJets(AliJet* /*genJets*/) {return kFALSE;}
  
  void ClearArray();
 
 protected:
  AliJetReader(const AliJetReader& rJetReader);
  AliJetReader& operator = (const AliJetReader& rhsr);
  TChain                  *fChain;            // chain for reconstructed tracks
  TClonesArray            *fMomentumArray;    // array of particle momenta
  TClonesArray            *fArrayMC;          //! array of mc particles
  TTask                   *fFillUnitArray;    //! task list for filling the UnitArray
  AliESDEvent             *fESD;              // pointer to esd
  AliJetReaderHeader      *fReaderHeader;     // pointer to header
  TArrayI                  fSignalFlag;       // to flag if a particle comes from pythia or
                                              // from the underlying event
  TArrayI                  fCutFlag;          // to flag if a particle passed the pt cut or not
  TClonesArray            *fUnitArray;        // array of digit position and energy 
  TRefArray               *fRefArray;         // array of digit position and energy 
  TClonesArray            *fUnitArrayNoCuts;  // array of digit position and energy 
  Bool_t                   fArrayInitialised; // To check that array of units is initialised  
  AliJetFillUnitArrayTracks        *fFillUAFromTracks;
  AliJetFillUnitArrayEMCalDigits   *fFillUAFromEMCalDigits;
  Int_t                    fNumCandidate;
  Int_t                    fNumCandidateCut;
  ClassDef(AliJetReader,1)
};
 
#endif
