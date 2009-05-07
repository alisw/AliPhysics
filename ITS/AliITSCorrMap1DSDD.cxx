/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the base class for SDD map 2D corrections   //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSCorrMap1DSDD.h"

ClassImp(AliITSCorrMap1DSDD)
//______________________________________________________________________
AliITSCorrMap1DSDD::AliITSCorrMap1DSDD():
AliITSCorrMapSDD()
{
  // default constructor
  ResetMap();
  SetNBinsAnode(1);
}
//______________________________________________________________________
AliITSCorrMap1DSDD::AliITSCorrMap1DSDD(Char_t *mapname):
AliITSCorrMapSDD(mapname)
{
  // standard constructor
  ResetMap();
  SetNBinsAnode(1);
}
//______________________________________________________________________
AliITSCorrMap1DSDD::AliITSCorrMap1DSDD(Char_t *mapname, Int_t nbinsdr):
AliITSCorrMapSDD(mapname)
{
  ResetMap();
  SetNBinsAnode(1);
  SetNBinsDrift(nbinsdr);
}
//______________________________________________________________________
void AliITSCorrMap1DSDD::ResetMap(){
  // Sets contents to zero
  for(Int_t iDr=0;iDr<kMaxNDriftPts; iDr++){
    fCorrMap[iDr]=0;
  }
}
//______________________________________________________________________
void AliITSCorrMap1DSDD::Set1DMap(TH1F* hmap){
  // Fill map staring from 1D histo of rediduals vs. x
  if(hmap->GetNbinsX()!=fNDriftPts){ 
    AliError(Form("N. of histo bins (%d) not matching N. of map cells (%d)\n",hmap->GetNbinsX(),fNDriftPts));
    return;
  }
  for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
    SetCellContent(0,iDr,hmap->GetBinContent(iDr+1));
  }
}

