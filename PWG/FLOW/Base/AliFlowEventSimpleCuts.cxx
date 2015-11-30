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

// AliFlowEventSimpleCuts:
// An event cut class for the flow framework
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include "AliFlowEventSimpleCuts.h"
#include "AliFlowEventSimple.h"

ClassImp(AliFlowEventSimpleCuts)

//-----------------------------------------------------------------------
AliFlowEventSimpleCuts::AliFlowEventSimpleCuts():
  TNamed(),
  fCutCentralityPercentile(kFALSE),
  fUseNewCentralityFramework(kFALSE),
  fCentralityPercentileMax(100.),
  fCentralityPercentileMin(0.)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowEventSimpleCuts::AliFlowEventSimpleCuts(const char* name, const char* title):
  TNamed(name, title),
  fCutCentralityPercentile(kFALSE),
  fUseNewCentralityFramework(kFALSE),
  fCentralityPercentileMax(100.),
  fCentralityPercentileMin(0.)
{
  //constructor 
}

////-----------------------------------------------------------------------
AliFlowEventSimpleCuts::AliFlowEventSimpleCuts(const AliFlowEventSimpleCuts& that):
  TNamed(that),
  fCutCentralityPercentile(that.fCutCentralityPercentile),
  fUseNewCentralityFramework(that.fUseNewCentralityFramework),
  fCentralityPercentileMax(that.fCentralityPercentileMax),
  fCentralityPercentileMin(that.fCentralityPercentileMin)
{
  //copy ctor
}

////-----------------------------------------------------------------------
AliFlowEventSimpleCuts::~AliFlowEventSimpleCuts()
{
  //dtor
}

////-----------------------------------------------------------------------
AliFlowEventSimpleCuts& AliFlowEventSimpleCuts::operator=(const AliFlowEventSimpleCuts& that)
{
  //assignment
  if (this==&that) return *this;

  fCutCentralityPercentile=that.fCutCentralityPercentile;
  fUseNewCentralityFramework=that.fUseNewCentralityFramework;
  fCentralityPercentileMax=that.fCentralityPercentileMax;
  fCentralityPercentileMin=that.fCentralityPercentileMin;
  return *this;
}

//----------------------------------------------------------------------- 
Bool_t AliFlowEventSimpleCuts::IsSelected(TObject* obj, TObject* /*mcobj*/)
{
  //check cuts
  AliFlowEventSimple* event = dynamic_cast<AliFlowEventSimple*>(obj);
  if (!event) return kFALSE;  //when passed wrong type of object
  Double_t centrality = event->GetCentrality();
  if (centrality>=fCentralityPercentileMin && centrality<fCentralityPercentileMax)
    return kTRUE;
  return kFALSE;
}
