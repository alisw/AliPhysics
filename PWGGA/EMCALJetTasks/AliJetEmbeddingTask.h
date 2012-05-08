#ifndef ALIJETEMBEDDINGTASK_H
#define ALIJETEMBEDDINGTASK_H

// $Id$

class TClonesArray;
class AliEMCALGeometry;

#include "AliAnalysisTaskSE.h"

class AliJetEmbeddingTask : public AliAnalysisTaskSE {
 public:
  AliJetEmbeddingTask();
  AliJetEmbeddingTask(const char *name); 
  virtual ~AliJetEmbeddingTask();

  void         Init();
  void         UserExec(Option_t* /*option*/);
  void         Terminate(const Option_t* /*option*/) {;}

  void         SetClusName(const char *n)            { fCaloName     = n;    }
  void         SetTracksName(const char *n)          { fTracksName   = n;    }
  void         SetEtaRange(Float_t min, Float_t max) { fEtaMin       = min;  fEtaMax = max; }
  void         SetPhiRange(Float_t min, Float_t max) { fPhiMin       = min;  fPhiMax = max; }
  void         SetPtRange(Float_t min, Float_t max)  { fPtMin        = min;  fPtMax  = max;  }
  void         SetCopyArray(Bool_t copy)             { fCopyArray    = copy; }
  void         SetNClusters(Int_t n)                 { fNEmbClusters = n;    }
  void         SetNTracks(Int_t n)                   { fNEmbTracks   = n;    }
  void         SetGeometryName(const char *n)        { fGeomName     = n;    }

 protected:

  virtual void           Embed();                 // do embedding

  TString                fGeomName;               // EMCal geometry name
  TString                fTracksName;             // name of track collection
  TString                fOutTracksName;          // name of output track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fOutCaloName;            // name of output cluster collection
  Float_t                fEtaMin;                 // eta minimum value
  Float_t                fEtaMax;                 // eta maximum value
  Float_t                fPhiMin;                 // phi minimum value
  Float_t                fPhiMax;                 // phi maximum value
  Float_t                fPtMin;                  // pt minimum value
  Float_t                fPtMax;                  // pt maximum value
  Bool_t                 fCopyArray;              // whether or not the array will be copied to a new one before embedding
  Int_t                  fNEmbClusters;           // how many clusters are being embedded
  Int_t                  fNEmbTracks;             // how many tracks are being embedded
  AliEMCALGeometry      *fGeom;                   //!pointer to EMCal geometry
  TClonesArray          *fClusters;               //!cluster collection
  TClonesArray          *fOutClusters;            //!output cluster collection
  TClonesArray          *fTracks;                 //!track collection
  TClonesArray          *fOutTracks;              //!output track collection

 private:
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 1) // Jet embedding task
};
#endif
