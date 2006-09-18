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
//  This is the base class for plots used to monitor the quality of the      //
//  recorded data.                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorPlot.h"


ClassImp(AliMonitorPlot) 


Bool_t  AliMonitorPlot::fgDrawRef = kTRUE;
Float_t AliMonitorPlot::fgThreshold = 5.0;

Color_t AliMonitorPlot::fgColorData = kBlack;
Color_t AliMonitorPlot::fgColorRef = kBlue;
Color_t AliMonitorPlot::fgColorCompare = kRed;


//_____________________________________________________________________________
AliMonitorPlot::AliMonitorPlot() :
  TNamed(),
  fDescription(),
  fNumberOfEvents(0)
{
// default contructor

}

//_____________________________________________________________________________
AliMonitorPlot::AliMonitorPlot(const AliMonitorPlot& plot) :
  TNamed(plot),
  fDescription(plot.fDescription),
  fNumberOfEvents(plot.fNumberOfEvents)
{
// copy constructor

}

//_____________________________________________________________________________
AliMonitorPlot& AliMonitorPlot::operator =(const AliMonitorPlot& plot)
{
// assignment operator

  TNamed::operator =(plot);
  fNumberOfEvents = plot.fNumberOfEvents;
  return *this;
}

//_____________________________________________________________________________
AliMonitorPlot::AliMonitorPlot(const char* name, const char* title) :
  TNamed(name, title),
  fDescription(),
  fNumberOfEvents(0)
{
// constructor setting name and title

}



//_____________________________________________________________________________
void AliMonitorPlot::SetDrawRef(Bool_t drawRef)
{
// set the flag for drawing the reference plot

  fgDrawRef = drawRef;
}

//_____________________________________________________________________________
void AliMonitorPlot::SetThreshold(Float_t threshold)
{
// set the threshold in standard deviations for the comparison 
// to the reference plot.
// no comparison is performed, if the threshold is <= 0

  fgThreshold = threshold;
}


//_____________________________________________________________________________
void AliMonitorPlot::DrawEvent(Int_t number)
{
// draw the normalized monitor plot together with the reference
// plot and the comparison plot (if available)
// for the "number"th last event

  if (!GetEvent(number)) return;
  if (fgDrawRef) ComparePlot();
  DrawPlot();
}

//_____________________________________________________________________________
void AliMonitorPlot::DrawSum(Int_t number)
{
// draw the normalized monitor plot together with the reference
// plot and the comparison plot (if available)
// for the sum of the last "number" events

  if (!GetSum(number)) return;
  if (fgDrawRef) ComparePlot();
  DrawPlot();
}

//_____________________________________________________________________________
void AliMonitorPlot::DrawRun()
{
// draw the normalized monitor plot together with the reference
// plot and the comparison plot (if available)
// for all monitored events of the current run

  if (!GetRun()) return;
  if (fgDrawRef) ComparePlot();
  DrawPlot();
}


//_____________________________________________________________________________
Bool_t AliMonitorPlot::CompareEvent(Int_t number)
{
// compare the normalized monitor plot for the "number"th last event
// to the reference plot

  if (!GetEvent(number)) return kTRUE;
  return ComparePlot();
}

//_____________________________________________________________________________
Bool_t AliMonitorPlot::CompareSum(Int_t number)
{
// compare the normalized monitor plot for the sum of the last 
// "number" events to the reference plot

  if (!GetSum(number)) return kTRUE;
  return ComparePlot();
}

//_____________________________________________________________________________
Bool_t AliMonitorPlot::CompareRun()
{
// compare the normalized monitor plot for all monitored events 
// of the current run to the reference plot

  if (!GetRun()) return kTRUE;
  return ComparePlot();
}

