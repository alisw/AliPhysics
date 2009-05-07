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
#include "TH2F.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSCorrMap2DSDD.h"

ClassImp(AliITSCorrMap2DSDD)
//______________________________________________________________________
AliITSCorrMap2DSDD::AliITSCorrMap2DSDD():
AliITSCorrMapSDD()
{
  // default constructor
  ResetMap();
}
//______________________________________________________________________
AliITSCorrMap2DSDD::AliITSCorrMap2DSDD(Char_t *mapname):
AliITSCorrMapSDD(mapname)
{
  // standard constructor
  ResetMap();
}
//______________________________________________________________________
AliITSCorrMap2DSDD::AliITSCorrMap2DSDD(Char_t *mapname, Int_t nbinsan, Int_t nbinsdr):
AliITSCorrMapSDD(mapname)
{
  // standard constructor
  ResetMap();
  SetNBinsAnode(nbinsan);
  SetNBinsDrift(nbinsdr);
}
//______________________________________________________________________
void AliITSCorrMap2DSDD::ResetMap(){
  // Sets contents to zero
  for(Int_t iAn=0;iAn<kMaxNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<kMaxNDriftPts; iDr++){
      fCorrMap[iAn][iDr]=0;
    }
  }
}

//______________________________________________________________________
void AliITSCorrMap2DSDD::Set2DMap(TH2F* hmap){
  // Fill map staring from 2D histo 
  // with anodes on x axis and drift dist. on y axis
  if(hmap->GetNbinsX()!=fNAnodePts || hmap->GetNbinsY()!=fNDriftPts){ 
    AliError(Form("N. of histo bins (%dX%d) not matching N. of map cells (%dX%d)\n",hmap->GetNbinsX(),hmap->GetNbinsY(),fNAnodePts,fNDriftPts));
    return;
  }
  for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
      SetCellContent(iAn,iDr,hmap->GetBinContent(iAn+1,iDr+1));
    }
  }
}
