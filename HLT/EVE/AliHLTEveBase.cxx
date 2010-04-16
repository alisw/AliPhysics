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

/// @file   AliHLTEveBase.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  Base class for the HLT eve display detector elements


#include "AliHLTEveBase.h"
#include "AliHLTHOMERBlockDesc.h"
//#include "TCollection.h"
#include "AliEveHOMERManager.h"
#include "TCanvas.h"
#include "TEveWindow.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliHLTEveBase);

AliHLTEveBase::AliHLTEveBase() : 
  fEventManager(NULL), 
  fCanvas(NULL),
  fHistoCount(0),
  fMaxHistos(0), 
  fDetector("")
{
  // Constructor.
}

AliHLTEveBase::~AliHLTEveBase()
{
  //Destructor

  if(fCanvas) 
    delete fCanvas;
  fCanvas = NULL;

  fEventManager = NULL;
}



TCanvas * AliHLTEveBase::CreateCanvas(TString  tabTitle, TString  canvasTitle ) {
   //See header file for documentation

  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(fEventManager->GetEveManager()->GetBrowser()->GetTabRight());
  slot->StartEmbedding();
  TCanvas * canvas = new TCanvas(canvasTitle.Data(),canvasTitle.Data(), 600, 400);
  slot->StopEmbedding(tabTitle.Data());

  return canvas;
}


void AliHLTEveBase::AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount ) {
  //See header file for documentation

   
  if ( ! block->GetClassName().CompareTo("TObjArray")) {
    TIter next((TObjArray*)(block->GetTObject()));
    TObject *object;
    
    while (( object = (TObject*) next())) {
      if (cdCount == fMaxHistos) {
	cout << "Too many histograms from detector, increase division size!!!!!!!!!!!! Detector: " << GetDetector() << endl;
	break;
      }
      TH2F* histo = dynamic_cast<TH2F*>(object);
      if(histo){
	canvas->cd(++cdCount);
	histo->Draw("COLZ");
      } else {
	TH1F* hist = dynamic_cast<TH1F*>(object);
	if (hist) {
	  canvas->cd(++cdCount);
	  hist->Draw();
	} else {
	  cout <<"AliHLTEveCaloBase::AddHistogramsTocCanvas: Histogram neither TH1F nor TH2F"<<endl;
	}
      }
    }
  }
    
  else if ( ! block->GetClassName().CompareTo("TH1F")) {

    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    ++cdCount;
    if(cdCount > fMaxHistos){
    
      cout << "Too many histograms, divide canvas more or create additional. Or ask svein to fix it! Detector: "<< GetDetector()<<endl;
      cdCount = 1;
    }  canvas->cd(cdCount);
    histo->Draw();
    
  } 
  
  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *histo = reinterpret_cast<TH2F*>(block->GetTObject());
    if (histo) {
      ++cdCount;

      if(cdCount > fMaxHistos){ 
	cdCount = 1;
	cout << "Too many histograms, divide canvas more or create additional. Or ask svein to fix it! Detector :" << GetDetector()<<endl;
      }

      canvas->cd(cdCount);
      histo->Draw("COLZ");
    }
  }

  canvas->cd();
}


