#ifndef ALIEMCALPARTICLEMAKER_H
#define ALIEMCALPARTICLEMAKER_H

// $Id$

class TClonesArray;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalParticleMaker : public AliAnalysisTaskEmcal {
 public:
  AliEmcalParticleMaker();
  AliEmcalParticleMaker(const char *name);
  virtual ~AliEmcalParticleMaker();

  const char *       GetClusOutName()   const           { return fCaloOutName;        }
  const char *       GetTracksOutName() const           { return fTracksOutName;      }
  void               SetTracksOutName(const char *name) { fTracksOutName      = name; }
  void               SetClusOutName(const char *name)   { fCaloOutName        = name; }

 protected:  
  void               ExecOnce();
  Bool_t             Run();

  TString            fTracksOutName;    // name of output track array
  TString            fCaloOutName;      // name of output clusters array
  TClonesArray      *fTracksOut;        //!track array out
  TClonesArray      *fCaloClustersOut;  //!track clusters out

 private:
  AliEmcalParticleMaker(const AliEmcalParticleMaker&);            // not implemented
  AliEmcalParticleMaker &operator=(const AliEmcalParticleMaker&); // not implemented

  ClassDef(AliEmcalParticleMaker, 2); // Task to make array of EMCAL particle
};
#endif
