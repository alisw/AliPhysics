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
AliITSCorrMapSDD::AliITSCorrMapSDD():TNamed("defaultmap",""),
  fNAnodePts(fgkNAnodePtsDefault),
  fNDriftPts(fgkNDriftPtsDefault),
  fXt1(0.),
  fXt2(0.),
  fXm1(0.),
  fXm2(0.),
  fDrLen(0.)
{
  // default constructor  
}
//______________________________________________________________________
AliITSCorrMapSDD::AliITSCorrMapSDD(Char_t *mapname):
  TNamed(mapname,""),
  fNAnodePts(fgkNAnodePtsDefault),
  fNDriftPts(fgkNDriftPtsDefault),
  fXt1(0.),
  fXt2(0.),
  fXm1(0.),
  fXm2(0.),
  fDrLen(0.)
{
  // standard constructor
}
//______________________________________________________________________
void AliITSCorrMapSDD::ComputeGridPoints(Float_t z, Float_t x, AliITSsegmentationSDD *seg, Bool_t isReco){
  // extracts points from the discrete grid with the correction map

  const Double_t kMicronTocm = 1.0e-4; 
  Int_t nAnodes=seg->Npz();
  Int_t nAnodesHybrid=seg->NpzHalf();
  Int_t bina =(Int_t) seg->GetAnodeFromLocal(x,z);
  if(bina>nAnodes)  AliError("Wrong anode anumber!");
  if(bina>=nAnodesHybrid) bina-=nAnodesHybrid;
  Float_t stept = seg->Dx()*kMicronTocm/(Float_t)fNDriftPts;
  fDrLen= seg->Dx()*kMicronTocm-TMath::Abs(x);
  Int_t bint = TMath::Abs((Int_t)(fDrLen/stept));
  if(bint==fNDriftPts) bint-=1;
  if(bint>=fNDriftPts){
    AliError("Wrong bin number along drift direction!");
    bint=fNDriftPts-1;
  }
  fXt1=stept*bint;
  fXm1=fXt1-GetCellContent(bina,bint)*kMicronTocm;
  if((bint+1)<fNDriftPts){
    fXt2=stept*(bint+1);
    fXm2=fXt2-GetCellContent(bina,bint+1)*kMicronTocm;
  }else{
    fXt2=stept*(bint-1);
    fXm2=fXt2-GetCellContent(bina,bint-1)*kMicronTocm;
  }
  if(isReco){
    if(fXm1<fDrLen && fXm2>fDrLen) return;
    if(bint==0 || bint==(fNDriftPts-1)) return;
    if(fXm1>fDrLen){
      for(Int_t itry=1; itry<=10; itry++){
	Float_t xmtest=(bint-itry)*stept-GetCellContent(bina,bint-itry)*kMicronTocm;
	if(xmtest<fDrLen){
	  fXt1=stept*(bint-itry);
	  fXt2=fXt1+stept;
	  fXm1=fXt1-GetCellContent(bina,bint-itry)*kMicronTocm;
	  fXm2=fXt2-GetCellContent(bina,bint+1-itry)*kMicronTocm;
	  return;
	}
      }
    }
    if(fXm2<fDrLen){
      for(Int_t itry=1; itry<=10; itry++){
	Float_t xmtest=(bint+1+itry)*stept-GetCellContent(bina,bint+1+itry)*kMicronTocm;
	if(xmtest>fDrLen){
	  fXt1=stept*(bint+itry);
	  fXt2=fXt1+stept;
	  fXm1=fXt1-GetCellContent(bina,bint+itry)*kMicronTocm;
	  fXm2=fXt2-GetCellContent(bina,bint+1+itry)*kMicronTocm;
	  return;
	}
      }
    }
  }
}
//______________________________________________________________________
Float_t AliITSCorrMapSDD::GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg){
  // returns correction in cm starting from local coordinates on the module
  ComputeGridPoints(z,x,seg,kTRUE);
  Float_t m=(fXt2-fXt1)/(fXm2-fXm1);
  Float_t q=fXt1-m*fXm1;
  Float_t xcorr=m*fDrLen+q;
  // fDrLen is the measured drift distance, xcorr is the corresponding true
  return (xcorr-fDrLen); 
}
//______________________________________________________________________
Float_t AliITSCorrMapSDD::GetShiftForSimulation(Float_t z, Float_t x, AliITSsegmentationSDD *seg){
  // returns shift to be appiled in digitizarion (in cm) starting from local coordinates on the module
  ComputeGridPoints(z,x,seg,kFALSE);
  Float_t m=(fXm2-fXm1)/(fXt2-fXt1);
  Float_t q=fXm1-m*fXt1;
  Float_t xshifted=m*fDrLen+q;
  // fDrLen is the true drift distance, xshifted is the one with map shift
  return (fDrLen-xshifted);
}
//______________________________________________________________________
TH2F* AliITSCorrMapSDD::GetMapHisto() const{
  // Returns a TH2F histogram with map of residuals
  TString hname;
  hname.Form("h%s",GetName());
  TH2F* hmap=new TH2F(hname.Data(),"",fNAnodePts,-0.5,255.5,fNDriftPts,0.,35.);
  for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
      hmap->SetBinContent(iAn+1,iDr+1,GetCellContent(iAn,iDr));
    }
  }
  return hmap;
}
//______________________________________________________________________
TH1F* AliITSCorrMapSDD::GetMapProfile() const{
  // Returns a TH1F with the projection of the map along drift coordinate
  TString hname;
  hname.Form("p%s",GetName());
  TH1F* hprof=new TH1F(hname.Data(),"",fNDriftPts,0.,35.);
  for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
    Float_t meanval=0.;
    for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
      meanval+=GetCellContent(iAn,iDr);
    }
    hprof->SetBinContent(iDr+1,meanval/fNAnodePts);
  }
  return hprof;
  
}
//______________________________________________________________________
TH1F* AliITSCorrMapSDD::GetResidualDistr(Float_t dmin, Float_t dmax) const{
  // Returns a TH1F histogram with distribution of residual
  TString hname;
  hname.Form("hd%s",GetName());
  TH1F* hd=new TH1F(hname.Data(),"",100,dmin,dmax);
  for(Int_t iAn=0;iAn<fNAnodePts; iAn++){
    for(Int_t iDr=0;iDr<fNDriftPts; iDr++){
      hd->Fill(GetCellContent(iAn,iDr));
    }
  }
  return hd;
}
