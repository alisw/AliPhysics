#ifndef ALIITSCORRMAP2DSDD_H
#define ALIITSCORRMAP2DSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for SDD maps in 2D used to correct for                  //
// voltage divider shape and doping fluctuations                 //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSsegmentationSDD.h"
#include<TNamed.h>
#include "AliLog.h"
class TH1F;
class TH2F;

class AliITSCorrMap2DSDD : public AliITSCorrMapSDD {

 public:
  AliITSCorrMap2DSDD();
  AliITSCorrMap2DSDD(Char_t *mapname);
  AliITSCorrMap2DSDD(Char_t *mapname, Int_t nbinsan, Int_t nbinsdr);
  virtual ~AliITSCorrMap2DSDD(){};

  virtual void ResetMap();
  virtual void Set2DMap(TH2F* hmap);
  virtual void SetCellContent(Int_t iAn, Int_t iTb, Float_t devMicron){
    if(CheckAnodeBounds(iAn) && CheckDriftBounds(iTb)) fCorrMap[iAn][iTb]=(Short_t)(devMicron*10.+0.5);
  }

  virtual Float_t GetCellContent(Int_t iAn, Int_t iTb) const {
   if(CheckAnodeBounds(iAn) && CheckDriftBounds(iTb)) return (Float_t)fCorrMap[iAn][iTb]/10.;
    else return 0.;
  }

 protected:
  Short_t fCorrMap[kMaxNAnodePts][kMaxNDriftPts];   // map of deviations
                                                // stored as Short_t: integer 
                                                // values from -32000 to 32000
                                                // in the range -3.2 - 3.2 mm

  ClassDef(AliITSCorrMap2DSDD,1);
};
#endif
