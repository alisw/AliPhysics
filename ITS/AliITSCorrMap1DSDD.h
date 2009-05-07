#ifndef ALIITSCORRMAP1DSDD_H
#define ALIITSCORRMAP1DSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for SDD maps in 1D used to correct for                  //
// voltage divider shape and doping fluctuations                 //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSsegmentationSDD.h"
#include<TNamed.h>
#include "AliLog.h"
class TH1F;

class AliITSCorrMap1DSDD : public AliITSCorrMapSDD {

 public:
  AliITSCorrMap1DSDD();
  AliITSCorrMap1DSDD(Char_t *mapname);
  AliITSCorrMap1DSDD(Char_t *mapname, Int_t nbinsdr);
  virtual ~AliITSCorrMap1DSDD(){};

  virtual void ResetMap();
  virtual void Set1DMap(TH1F* hmap);
  virtual void SetCellContent(Int_t /*iAn*/, Int_t iTb, Float_t devMicron){
    if(CheckDriftBounds(iTb)) fCorrMap[iTb]=(Short_t)(devMicron*10.+0.5);
  }

  virtual Float_t GetCellContent(Int_t /*iAn*/, Int_t iTb) const {
    if(CheckDriftBounds(iTb)) return (Float_t)fCorrMap[iTb]/10.;
    else return 0.;
  }

 protected:
  Short_t fCorrMap[kMaxNDriftPts];           // map of deviations
                                       // stored as Short_t: integer 
                                       // values from -32000 to 32000
                                       // in the range -3.2 - 3.2 mm

  ClassDef(AliITSCorrMap1DSDD,1);
};
#endif
