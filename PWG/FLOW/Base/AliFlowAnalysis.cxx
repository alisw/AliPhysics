/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

#define ALIFLOWANALYSISBASE_CXX
 
#include "AliFlowAnalysis.h"
#include "AliFlowEventSimpleCuts.h"
#include "AliFlowEventSimple.h"

//////////////////////////////////////////////////////////////////////////////
// AliFlowAnalysis:
// Description: Base class for flow analysis methods
// author: M. Krzewicki mikolaj.krzewicki@cern.ch
//////////////////////////////////////////////////////////////////////////////

ClassImp(AliFlowAnalysis)

//-----------------------------------------------------------------------
AliFlowAnalysis::AliFlowAnalysis(const char* name):
TNamed(name,name),
fEventCuts(NULL),
fPOItype(1)
{
  //ctor
}

//-----------------------------------------------------------------------
AliFlowAnalysis::~AliFlowAnalysis() 
{
  //destructor
  delete fEventCuts;
}

//-----------------------------------------------------------------------
void AliFlowAnalysis::ProcessEvent(AliFlowEventSimple* event)
{
  //process the event
  //only call the analysis if we pass the event cuts (if any)
  Bool_t pass=kTRUE;
  if (fEventCuts) pass=fEventCuts->IsSelected(event,(TObject*)NULL);
  if (pass) Make(event);
}
