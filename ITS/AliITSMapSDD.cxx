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
      fMap[iAn][iDr]=hmap->GetBinContent(iAn+1,iDr+1);
    }
  }
}
//______________________________________________________________________
TH2F* AliITSMapSDD::GetMapHisto() const{
  // Returns a TH2F histogram with map of residuals
  Char_t hname[50];
  sprintf(hname,"h%s",GetName());
  TH2F* hmap=new TH2F(hname,"",fgkNAnodPts,-0.5,255.5,fgkNDrifPts,0.,35.);
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      hmap->SetBinContent(iAn+1,iDr+1,fMap[iAn][iDr]);
    }
  }
  return hmap;
}
//______________________________________________________________________
TH1F* AliITSMapSDD::GetResidualDistr(Float_t dmin, Float_t dmax) const{
  // Returns a TH1F histogram with distribution of residual
  Char_t hname[50];
  sprintf(hname,"hd%s",GetName());
  TH1F* hd=new TH1F(hname,"",100,dmin,dmax);
  for(Int_t iAn=0;iAn<fgkNAnodPts; iAn++){
    for(Int_t iDr=0;iDr<fgkNDrifPts; iDr++){
      hd->Fill(fMap[iAn][iDr]);
    }
  }
  return hd;
}
