#ifndef ALIJETPARTICLESREADERHLT_H
#define ALIJETPARTICLESREADERHLT_H

//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// File reader for HLT tracks ESD                                          //
//                                                                         //
// loizides@ikf.uni-frankfurt.de                                           //
/////////////////////////////////////////////////////////////////////////////

#include "AliJetParticlesReaderESD.h"

class AliJetParticlesReaderHLT: public AliJetParticlesReaderESD
{
  public:
  AliJetParticlesReaderHLT(Bool_t bMapper, const Char_t* esdfilename = "AliESDs.root") ;
  AliJetParticlesReaderHLT(Bool_t bMapper, TObjArray* dirs, const Char_t* esdfilename = "AliESDs.root");

  virtual ~AliJetParticlesReaderHLT();

  void SetMinHits(Int_t i){fMinHits=i;}
  void SetMinWeight(Int_t i){fMinWeight=i;}

  protected:
  Int_t    ReadESD(AliESD* esd); //read esd file/objects

  Bool_t fTrackerType; //track type: 0==Conformal Tracker, 1==Hough Tracker
  Int_t fMinHits;      //minimum hits required
  Int_t fMinWeight;    //minimum weight required 

  ClassDef(AliJetParticlesReaderHLT,1) //
};

#endif
