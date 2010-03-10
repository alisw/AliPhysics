/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  ISSD class for the HLT EVE display

#include "AliHLTEveISSD.h"
#include "TEvePointSet.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"

ClassImp(AliHLTEveISSD);

AliHLTEveISSD::AliHLTEveISSD() : 
  AliHLTEveITS("ISSD"),
  f2DCanvas(NULL),
  f2DHistoCount(0)
{
  // Constructor.
}

AliHLTEveISSD::~AliHLTEveISSD()
{
  //Destructor
  if(f2DCanvas)
    delete f2DCanvas;
  f2DCanvas = NULL;
}

void AliHLTEveISSD::SetUpPointSet(TEvePointSet * ps) {
  //See header file for documentation
  ps->SetMainColor(kBlue);
  ps->SetMarkerStyle((Style_t)kFullDotMedium);
}


void AliHLTEveISSD::AddHistogramsToCanvas(AliHLTHOMERBlockDesc* block, TCanvas *canvas, Int_t &cdCount ) {
  //See header file for documentation


  if(!f2DCanvas) f2DCanvas = CreateCanvas("ISSD 2D QA", "ISSD 2D QA");


  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if(histo){
      ++cdCount;
      canvas->cd(cdCount);
      histo->Draw();
    }

  } else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista) {
      ++f2DHistoCount;
      f2DCanvas->cd(f2DHistoCount);
      hista->Draw("COLZ");
    }
  
  } else if ( ! block->GetClassName().CompareTo("TObjArray")) {
    TIter next((TObjArray*)(block->GetTObject()));
    TObject *object;
    while (( object = (TObject*) next())) {
      TString string;
      string = "TH1F";
      TString string2;
      string2 = "TH2F";

      if ( !(string.CompareTo(object->ClassName())) ) {
	TH1F* histo = reinterpret_cast<TH1F*>(object);
	++cdCount;
	canvas->cd(cdCount);
	histo->Draw();
	
      
      } else if ( !(string2.CompareTo(object->ClassName()) ) ) {
	TH2F* histo = reinterpret_cast<TH2F*>(object);
	++f2DHistoCount;
	f2DCanvas->cd(f2DHistoCount);
	histo->Draw("COLZ");
      }
    }
  }

  canvas->cd();  f2DCanvas->cd();
  
}


void AliHLTEveISSD::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  if(f2DCanvas) f2DCanvas->Update();
  if(fPointSet) fPointSet->ElementChanged();

}

void AliHLTEveISSD::ResetElements() {
  //See header file for documentation
  fHistoCount = 0;
  f2DHistoCount = 0;
  if(fPointSet) fPointSet->Reset();
}
