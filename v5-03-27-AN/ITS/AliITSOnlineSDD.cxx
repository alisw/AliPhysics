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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the base class for SDD detector algorithms  //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSOnlineSDD.h"
#include "TH2F.h"

ClassImp(AliITSOnlineSDD)
//______________________________________________________________________
  AliITSOnlineSDD::AliITSOnlineSDD():TObject(),fDDL(0),fCarlos(0),fSide(0),fFirstGoodTB(0),fLastGoodTB(0)
{
  // default constructor
  SetFirstGoodTB();
  SetLastGoodTB();
}
//______________________________________________________________________
  AliITSOnlineSDD::AliITSOnlineSDD(Int_t nddl, Int_t ncarlos, Int_t sid):TObject(),fDDL(0),fCarlos(0),fSide(0),fFirstGoodTB(0),fLastGoodTB(0)
{
  // standard constructor
  SetDDL(nddl);
  SetCarlos(ncarlos);
  SetDetectorSide(sid);
  SetFirstGoodTB();
  SetLastGoodTB();
}
//______________________________________________________________________
TH2F* AliITSOnlineSDD::ApplyZeroSuppression(TH2F* hRaw, Float_t basl, Int_t tL, Int_t tH){
  // apply zero suppression 

  Int_t nx=hRaw->GetNbinsX();
  Float_t xmin=hRaw->GetXaxis()->GetXmin();
  Float_t xmax=hRaw->GetXaxis()->GetXmax();
  Int_t ny=hRaw->GetNbinsY();
  Float_t ymin=hRaw->GetYaxis()->GetXmin();
  Float_t ymax=hRaw->GetYaxis()->GetXmax();
  TH2F* hZeroSupp=new TH2F(Form("%s_ZeroSupp",hRaw->GetName()),"",nx,xmin,xmax,ny,ymin,ymax);
  for(Int_t ix=1; ix<=nx;ix++){
    for(Int_t iy=1; iy<=ny;iy++){
      //                     N
      // Get "quintuple":   WCE
      //                     S
      Float_t cC=hRaw->GetBinContent(ix,iy)-basl;
      Int_t nLow=0, nHigh=0;      
      if(cC<tL){
	hZeroSupp->SetBinContent(ix,iy,0.);
	continue;
      }
      nLow++; // cC is greater than tL
      if(cC>tH) nHigh++;
      Float_t wW=0.;
      if(iy>0) wW=hRaw->GetBinContent(ix,iy-1)-basl;
      if(wW>tL) nLow++;
      if(wW>tH) nHigh++;
      Float_t eE=0.;
      if(iy<ny) eE=hRaw->GetBinContent(ix,iy+1)-basl;
      if(eE>tL) nLow++;
      if(eE>tH) nHigh++;
      Float_t nN=0.;
      if(ix<nx) nN=hRaw->GetBinContent(ix+1,iy)-basl;
      if(nN>tL) nLow++;
      if(nN>tH) nHigh++;
      Float_t sS=0.;
      if(ix>0) sS=hRaw->GetBinContent(ix-1,iy)-basl;
      if(sS>tL) nLow++;
      if(sS>tH) nHigh++;
      if(nLow>=2 && nHigh>=1){
	hZeroSupp->SetBinContent(ix,iy,cC);
      }else{
	hZeroSupp->SetBinContent(ix,iy,0.);
      }
    }
  }
  return hZeroSupp;
}
