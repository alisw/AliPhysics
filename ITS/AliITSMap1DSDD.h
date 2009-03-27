#ifndef ALIITSMAP1DSDD_H
#define ALIITSMAP1DSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

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

class AliITSMap1DSDD : public AliITSMapSDD {

 public:
  AliITSMap1DSDD();
  AliITSMap1DSDD(Char_t *mapname);
  AliITSMap1DSDD(Char_t *mapname, Int_t nbinsdr);
  virtual ~AliITSMap1DSDD(){};

  virtual void ResetMap();
  virtual void Set1DMap(TH1F* hmap);
  virtual void SetCellContent(Int_t /*iAn*/, Int_t iTb, Float_t devMicron){
    if(CheckDriftBounds(iTb)) fMap[iTb]=(Short_t)(devMicron*10.+0.5);
  }

  virtual Float_t GetCellContent(Int_t /*iAn*/, Int_t iTb) const {
    if(CheckDriftBounds(iTb)) return (Float_t)fMap[iTb]/10.;
    else return 0.;
  }

 protected:
  Short_t fMap[kMaxNDriftPts];           // map of deviations
                                       // stored as Short_t: integer 
                                       // values from -32000 to 32000
                                       // in the range -3.2 - 3.2 mm

  ClassDef(AliITSMap1DSDD,1);
};
#endif
