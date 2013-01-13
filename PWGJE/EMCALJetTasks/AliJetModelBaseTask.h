#ifndef ALIJETMODELBASETASK_H
#define ALIJETMODELBASETASK_H

// $Id$

class TClonesArray;
class AliEMCALGeometry;
class AliVCluster;
class AliPicoTrack;

#include <TH1F.h>
#include <TF1.h>

#include "AliAnalysisTaskSE.h"

class AliJetModelBaseTask : public AliAnalysisTaskSE {
 public:
  AliJetModelBaseTask();
  AliJetModelBaseTask(const char *name); 
  virtual ~AliJetModelBaseTask();

  void                   UserExec(Option_t* /*option*/);

  void                   SetClusName(const char *n)            { fCaloName     = n;    }
  void                   SetCopyArray(Bool_t copy)             { fCopyArray    = copy; }
  void                   SetEtaRange(Float_t min, Float_t max) { fEtaMin       = min;  fEtaMax = max; }
  void                   SetGeometryName(const char *n)        { fGeomName     = n;    }
  void                   SetMarkMC(Bool_t m)                   { fMarkMC       = m;    }
  void                   SetNClusters(Int_t n)                 { fNClusters    = n;    }
  void                   SetNTracks(Int_t n)                   { fNTracks      = n;    }
  void                   SetPhiRange(Float_t min, Float_t max) { fPhiMin       = min;  fPhiMax = max; }
  void                   SetPtRange(Float_t min, Float_t max)  { fPtMin        = min;  fPtMax  = max;  }
  void                   SetPtSpectrum(TH1 *f)                 { fPtSpectrum   = f;    }
  void                   SetPtSpectrum(TF1 *f)                 { fPtSpectrum   = new TH1F("ptSpectrum","ptSpectrum",250,f->GetXmin(),f->GetXmax()); 
                                                                 fPtSpectrum->Add(f); }
  void                   SetSuffix(const char *s)              { fSuffix       = s;    }
  void                   SetTracksName(const char *n)          { fTracksName   = n;    }

 protected:
  AliVCluster           *AddCluster(Double_t e = -1, Double_t eta = -999, Double_t phi = -1);   // add a cluster; if values are -1 generate random parameters
  AliVCluster           *AddCluster(Double_t e, Int_t absId);                                   // add a cluster with given energy and position
  AliPicoTrack          *AddTrack(Double_t pt = -1, Double_t eta = -999, Double_t phi = -1);    // add a track; if values are -1 generate random parameters
  void                   CopyClusters();
  void                   CopyTracks();
  virtual Bool_t         ExecOnce();
  void                   GetRandomCell(Double_t &eta, Double_t &phi, Int_t &absId);             // generate a random cell in the calorimeter
  Double_t               GetRandomEta();                                                        // generate a random eta value in the given range
  Double_t               GetRandomPhi();                                                        // generate a random phi value in the given range
  Double_t               GetRandomPt();                                                         // generate a random pt value in the given range
  virtual void           Run();                                                                 // do jet model action

  TString                fGeomName;               // EMCal geometry name
  TString                fTracksName;             // name of track collection
  TString                fOutTracksName;          // name of output track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fOutCaloName;            // name of output cluster collection
  TString                fSuffix;                 // suffix to add in the name of new collections
  Float_t                fEtaMin;                 // eta minimum value
  Float_t                fEtaMax;                 // eta maximum value
  Float_t                fPhiMin;                 // phi minimum value
  Float_t                fPhiMax;                 // phi maximum value
  Float_t                fPtMin;                  // pt minimum value
  Float_t                fPtMax;                  // pt maximum value
  Bool_t                 fCopyArray;              // whether or not the array will be copied to a new one before modelling
  Int_t                  fNClusters;              // how many clusters are being processed
  Int_t                  fNTracks;                // how many tracks are being processed
  Bool_t                 fMarkMC;                 // whether or not mark new tracks/cluster as MC
  TH1                   *fPtSpectrum;             // pt spectrum parametrization to extract random pt values
  Bool_t                 fIsInit;                 //=true if initialized
  AliEMCALGeometry      *fGeom;                   //!pointer to EMCal geometry
  TClonesArray          *fClusters;               //!cluster collection
  TClonesArray          *fOutClusters;            //!output cluster collection
  TClonesArray          *fTracks;                 //!track collection
  TClonesArray          *fOutTracks;              //!output track collection

 private:
  AliJetModelBaseTask(const AliJetModelBaseTask&);            // not implemented
  AliJetModelBaseTask &operator=(const AliJetModelBaseTask&); // not implemented

  ClassDef(AliJetModelBaseTask, 4) // Jet modelling task
};
#endif
