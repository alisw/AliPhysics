#ifndef ALIJETPARTICLESREADERHLT_H
#define ALIJETPARTICLESREADERHLT_H

//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// File reader for HLT tracks ESD                                          //
//                                                                         //
// loizides@ikf.uni-frankfurt.de                                           //
/////////////////////////////////////////////////////////////////////////////

//#include <TString.h>
//#include <AliESDtrack.h>
#include "AliJetParticlesReaderESD.h"

//class TFile;
//class AliESD;

class AliJetParticlesReaderHLT: public AliJetParticlesReaderESD
{
  public:
  AliJetParticlesReaderHLT(Bool_t bMapper, const Char_t* esdfilename = "AliESDs.root") ;
  AliJetParticlesReaderHLT(Bool_t bMapper, TObjArray* dirs, const Char_t* esdfilename = "AliESDs.root");

  virtual ~AliJetParticlesReaderHLT();

  protected:
  Int_t    ReadESD(AliESD* esd); //read esd file/objects

  Bool_t fTrackerType; //track type: 0==Conformal Tracker, 1==Hough Tracker

  ClassDef(AliJetParticlesReaderHLT,1) //
};

#endif
