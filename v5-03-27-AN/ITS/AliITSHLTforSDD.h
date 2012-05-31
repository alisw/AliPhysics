#ifndef ALIITSHLTFORSDD_H
#define ALIITSHLTFORSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Class to store the HLT status used to define the SDD raw data format //
// (when HLT is in mode C SDD data are compressed,                      // 
// see AliITSCompressRawDataSDD)                                        //
// Origin: F.Prino, Torino, prino@to.infn.it                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include<TObject.h>

class TString;

class AliITSHLTforSDD : public TObject {

 public:
  AliITSHLTforSDD();
  AliITSHLTforSDD(TString hltMode);
  virtual ~AliITSHLTforSDD(){};

  void SetHLTmodeC(Bool_t isHLTmodC){fHLTmodeC=isHLTmodC;}
  Bool_t IsHLTmodeC() const {return fHLTmodeC;}


 protected:
  Bool_t fHLTmodeC;  // flag for the HLT status

  ClassDef(AliITSHLTforSDD,1);
};
#endif
