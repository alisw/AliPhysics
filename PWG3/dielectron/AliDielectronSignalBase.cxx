/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron SignalBase                                  //
//                                                                       //
//                                                                       //
/*
Base class for signal extraction from a histogram or an array of histograms
The histogram is assumed to be an inv. mass spectrum,
the array of histograms is assumed to be an array with inv. mass histograms
resulting from single and mixed events, as defined in AliDielectron.cxx

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "TPaveText.h"
#include "AliDielectronSignalBase.h"

ClassImp(AliDielectronSignalBase)

AliDielectronSignalBase::AliDielectronSignalBase() :
  TNamed(),
  fValues(4),
  fErrors(4),
  fIntMin(2.99),
  fIntMax(3.19)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalBase::AliDielectronSignalBase(const char* name, const char* title) :
  TNamed(name, title),
  fValues(4),
  fErrors(4),
  fIntMin(2.99),
  fIntMax(3.19)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronSignalBase::~AliDielectronSignalBase()
{
  //
  // Default Destructor
  //
  
}

//______________________________________________
TPaveText* AliDielectronSignalBase::DrawStats(Double_t x1/*=0.*/, Double_t y1/*=0.*/, Double_t x2/*=0.*/, Double_t y2/*=0.*/)
{
  //
  // Draw extracted values in a TPaveText
  // with the corners x1,y2,x2,y2
  //
  if (TMath::Abs(x1)<1e-20&&TMath::Abs(x2)<1e-20){
    x1=.6;
    x2=.9;
    y1=.7;
    y2=.9;
  }
  TPaveText *t=new TPaveText(x1,y1,x2,y2,"brNDC");
  t->SetFillColor(kWhite);
  t->SetBorderSize(1);
  t->SetTextAlign(12);
  t->AddText(Form("Singal : %.2f #pm %.2f",GetSignal(),GetSignalError()));
  t->AddText(Form("Backgnd: %.2f #pm %.2f",GetBackground(),GetBackgroundError()));
  t->AddText(Form("Signif.: %.2f #pm %.2f",GetSignificance(),GetSignificanceError()));
  t->AddText(Form("SoB    : %.2f #pm %.2f",GetSignalOverBackground(),GetSignalOverBackgroundError()));
  t->Draw();

  return t;
}

