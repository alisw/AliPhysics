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
// Implementation of the base class for SDD map corrections      //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "TH2F.h"
#include "AliITSCorrMapSDD.h"

const Int_t AliITSCorrMapSDD::fgkNAnodePtsDefault = 1;
const Int_t AliITSCorrMapSDD::fgkNDriftPtsDefault = 72;

ClassImp(AliITSCorrMapSDD)
//______________________________________________________________________
AliITSCorrMapSDD::AliITSCorrMapSDD():
TNamed("defaultmap",""),
fNAnodePts(fgkNAnodePtsDefault),
fNDriftPts(fgkNDriftPtsDefault)
{
  // default constructor  
}
//______________________________________________________________________
AliITSCorrMapSDD::AliITSCorrMapSDD(Char_t *mapname):
TNamed(mapname,""),
fNAnodePts(fgkNAnodePtsDefault),
fNDriftPts(fgkNDriftPtsDefault)
{
  // standard constructor
}
//______________________________________________________________________
Float_t AliITSCorrMapSDD::GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg){
  // returns correction in cm starting from local coordinates on the module
  const Double_t kMicronTocm = 1.0e-4; 
  Int_t nAnodes=seg->Npz();
  Int_t nAnodesHybrid=seg->NpzHalf();
  Int_t bina =(Int_t) seg->GetAnodeFromLocal(x,z);
  if(bina>nAnodes)  AliError("Wrong anode anumber!");
  if(bina>=nAnodesHybrid) bina-=nAnodesHybrid;
  Float_t stept = seg->Dx()*kMicronTocm/(Float_t)fNDriftPts;
  Int_t bint = TMath::Abs((Int_t)(x/stept));
  if(bint==fNDriftPts) bint-=1;
  if(bint>=fNDriftPts) AliError("Wrong bin number along drift direction!");
  return kMicronTocm*GetCellContent(bina,bint);
}
//______________________________________________________________________
TH2F* AliITSCorrMapSDD::GetMapHisto() const{
  // Returns a TH2F histogram with map of residuals
  Char_t hname[50];
  sprintf(hname,"h%s",GetName());
  TH2F* hmap=new TH2F(hname,"",fNAnodePts,-0.5,255.5,fNDriftPts,0.,35.);
  for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
      hmap->SetBinContent(iAn+1,iDr+1,GetCellContent(iAn,iDr));
    }
  }
  return hmap;
}
//______________________________________________________________________
TH1F* AliITSCorrMapSDD::GetResidualDistr(Float_t dmin, Float_t dmax) const{
  // Returns a TH1F histogram with distribution of residual
  Char_t hname[50];
  sprintf(hname,"hd%s",GetName());
  TH1F* hd=new TH1F(hname,"",100,dmin,dmax);
  for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
      hd->Fill(GetCellContent(iAn,iDr));
    }
  }
  return hd;
}
