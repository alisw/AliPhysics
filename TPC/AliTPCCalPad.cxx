/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for parameters which saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include <TObjArray.h>
#include <TAxis.h>
#include <TGraph.h>
ClassImp(AliTPCCalPad)

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad():TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = 0;
  }

}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTPCCalPad constructor
  //
  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = new AliTPCCalROC(isec);
  }
}


//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const AliTPCCalPad &c):TNamed(c)
{
  //
  // AliTPCCalPad copy constructor
  //

  ((AliTPCCalPad &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(TObjArray * array):TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = (AliTPCCalROC *)array->At(isec);
  }

}


///_____________________________________________________________________________
AliTPCCalPad::~AliTPCCalPad()
{
  //
  // AliTPCCalPad destructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      delete fROC[isec];
      fROC[isec] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTPCCalPad &AliTPCCalPad::operator=(const AliTPCCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      fROC[isec]->Copy(*((AliTPCCalPad &) c).fROC[isec]);
    }
  }
  TObject::Copy(c);
}

TGraph  *  AliTPCCalPad::MakeGraph(Int_t type, Float_t ratio){
  //
  //   type=1 - mean
  //        2 - median
  //        3 - LTM
  Int_t npoints = 0;
  for (Int_t i=0;i<72;i++) if (fROC[i]) npoints++;
  TGraph * graph = new TGraph(npoints);
  npoints=0;   
  for (Int_t isec=0;isec<72;isec++){
    if (!fROC[isec]) continue;
    if (type==0)  graph->SetPoint(npoints,isec,fROC[isec]->GetMean());      
    if (type==1)  graph->SetPoint(npoints,isec,fROC[isec]->GetMedian());
    if (type==2)  graph->SetPoint(npoints,isec,fROC[isec]->GetLTM(0,ratio));    
    npoints++;
  }

  graph->GetXaxis()->SetTitle("Sector"); 
  if (type==0) {
    graph->GetYaxis()->SetTitle("Mean");   
    graph->SetMarkerStyle(22);    
  }
  if (type==1) {
    graph->GetYaxis()->SetTitle("Median");   
    graph->SetMarkerStyle(22);    
  }
  if (type==2) {
      graph->GetYaxis()->SetTitle(Form("Mean%f",ratio));      
      graph->SetMarkerStyle(24);
  }

  return graph;
}
