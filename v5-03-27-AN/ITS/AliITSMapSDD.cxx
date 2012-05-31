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
#include "AliITSMapSDD.h"
#include "AliITSCorrMap1DSDD.h"

ClassImp(AliITSMapSDD)
//______________________________________________________________________
AliITSMapSDD::AliITSMapSDD():TNamed("defaultmap","")
{
  // default constructor
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      fMap[iAn][iDr]=0;
    }
  }
}
//______________________________________________________________________
AliITSMapSDD::AliITSMapSDD(Char_t *mapname):TNamed(mapname,"")
{
  // standard constructor
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      fMap[iAn][iDr]=0;
    }
  }
}

//______________________________________________________________________
void AliITSMapSDD::SetMap(TH2F* hmap){
  // Fill map staring from 2D histo 
  // with anodes on x axis and drift dist. on y axis
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      SetCellContent(iAn,iDr,hmap->GetBinContent(iAn+1,iDr+1));
    }
  }
}
//______________________________________________________________________
Float_t AliITSMapSDD::GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg){
  // returns correction in cm starting from local coordinates on the module
  const Double_t kMicronTocm = 1.0e-4; 
  Int_t nAnodes=seg->Npz();
  Int_t nAnodesHybrid=seg->NpzHalf();
  Int_t bina =(Int_t) seg->GetAnodeFromLocal(x,z);
  if(bina>nAnodes)  AliError("Wrong anode anumber!");
  if(bina>=nAnodesHybrid) bina-=nAnodesHybrid;
  Float_t stept = seg->Dx()*kMicronTocm/(Float_t)fgkNDrifPts;
  Int_t bint = TMath::Abs((Int_t)(x/stept));
  if(bint==fgkNDrifPts) bint-=1;
  if(bint>=fgkNDrifPts) AliError("Wrong bin number along drift direction!");
  return kMicronTocm*GetCellContent(bina,bint);
}
//______________________________________________________________________
TH2F* AliITSMapSDD::GetMapHisto() const{
  // Returns a TH2F histogram with map of residuals
  TString hname;
  hname.Form("h%s",GetName());
  TH2F* hmap=new TH2F(hname.Data(),"",fgkNAnodPts,-0.5,255.5,fgkNDrifPts,0.,35.);
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      hmap->SetBinContent(iAn+1,iDr+1,GetCellContent(iAn,iDr));
    }
  }
  return hmap;
}
//______________________________________________________________________
TH1F* AliITSMapSDD::GetResidualDistr(Float_t dmin, Float_t dmax) const{
  // Returns a TH1F histogram with distribution of residual
  TString hname;
  hname.Form("hd%s",GetName());
  TH1F* hd=new TH1F(hname.Data(),"",100,dmin,dmax);
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      hd->Fill(GetCellContent(iAn,iDr));
    }
  }
  return hd;
}
//______________________________________________________________________
AliITSCorrMapSDD* AliITSMapSDD::ConvertToNewFormat() const{
  // convert correction map to new format  
  Char_t* name=(Char_t*)GetName();
  AliITSCorrMapSDD* newmap=new AliITSCorrMap1DSDD(name,fgkNDrifPts);
  for(Int_t i=0; i<fgkNDrifPts; i++){
    newmap->SetCellContent(0,i,GetCellContent(0,i));
  }
  return newmap;
}

