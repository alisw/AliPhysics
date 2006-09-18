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
//  This is the base class for the creation and filling of monitor           //
//  histograms.                                                              //
//  Derived classes have to implement the methods CreateHistos and           //
//  FillHistos. CreateHistos has to create the monitor histograms, put       //
//  them into a new folder and add this folder to the given root folder.     //
//  FillHistos has to fill the data of the current event into the created    //
//  monitor histograms.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>

#include "AliLog.h"
#include "AliMonitorTrend.h"

#include "AliMonitor.h"


ClassImp(AliMonitor) 


//_____________________________________________________________________________
AliMonitor::AliMonitor():
  TObject(),
  fFolder(NULL)
{
}

//_____________________________________________________________________________
void AliMonitor::CreateBranches(TTree*)
{
// add branches to the monitor tree
// by default no branches are added
// this method can be overwritten by derived classes

}


//_____________________________________________________________________________
AliMonitorHisto* AliMonitor::CreateHisto1(const char* name, const char* title,
			      Int_t xBins, Double_t xMin, Double_t xMax,
			      const char* xTitle, const char* yTitle,
			      AliMonitorHisto::ENorm norm)
{
// create a 1 dimensional monitor histogram and add it to fFolder

  TH1F* histo = new TH1F(name, title, xBins, xMin, xMax);
  histo->SetMarkerStyle(kFullCircle);
  histo->GetXaxis()->SetTitle(xTitle);
  histo->GetYaxis()->SetTitle(yTitle);
  AliMonitorHisto* result = new AliMonitorHisto(histo, norm);
  fFolder->Add(result);
  return result;
}

//_____________________________________________________________________________
AliMonitorHisto* AliMonitor::CreateHisto2(const char* name, const char* title,
			      Int_t xBins, Double_t xMin, Double_t xMax,
			      Int_t yBins, Double_t yMin, Double_t yMax,
			      const char* xTitle, const char* yTitle,
			      const char* zTitle,
			      AliMonitorHisto::ENorm norm)
{
// create a 2 dimensional monitor histogram and add it to fFolder

  TH2F* histo = new TH2F(name, title, xBins, xMin, xMax, yBins, yMin, yMax);
  histo->SetOption("BOX");
  histo->GetXaxis()->SetTitle(xTitle);
  histo->GetYaxis()->SetTitle(yTitle);
  histo->GetZaxis()->SetTitle(zTitle);
  AliMonitorHisto* result = new AliMonitorHisto(histo, norm);
  fFolder->Add(result);
  return result;
}

//_____________________________________________________________________________
AliMonitorTrend* AliMonitor::CreateTrend(const char* name, const char* title,
					 const char* label, 
					 Double_t min, Double_t max)
{
// create a trend monitor histogram and add it to fFolder

  AliMonitorTrend* result = new AliMonitorTrend(name, title, label, min, max);
  fFolder->Add(result);
  return result;
}

