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

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSClusterizer.h"
#include "AliPHOSGetter.h" 
#include "AliPHOSQualAssDataMaker.h" 

// --- AliRoot header files ---
#include "AliPHOSQualAssDataMaker.h"

ClassImp(AliPHOSClusterizer)

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer():
  TTask("",""),
  fEventFolderName(""),
  fFirstEvent(0),
  fLastEvent(-1),
  fRawReader(0),
  fQADM(0)
{
 // ctor
}

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer(const TString alirunFileName, 
				       const TString eventFolderName):
  TTask("PHOS"+AliConfig::Instance()->GetReconstructionerTaskName(), 
	alirunFileName), fEventFolderName(eventFolderName),
  fFirstEvent(0),
  fLastEvent(-1),
  fRawReader(0),
  fQADM(0)
{
  // ctor
  fQADM = new AliPHOSQualAssDataMaker() ;  
  //initiaizes the quality assurance data maker
  fQADM ->Init(AliQualAss::kRECPOINTS) ;    
}

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer(const AliPHOSClusterizer & clusterizer) :
  TTask(clusterizer),fEventFolderName(clusterizer.GetEventFolderName()),
  fFirstEvent(clusterizer.GetFirstEvent()),fLastEvent(clusterizer.GetLastEvent()),
  fRawReader(clusterizer.GetRawReader()),
  fQADM(clusterizer.GetQualAssDataMaker())
{
  //Copy constructor
  //initiaizes the quality assurance data maker
  GetQualAssDataMaker()->Init(AliQualAss::kRECPOINTS) ;    
}
//____________________________________________________________________________
AliPHOSClusterizer::~AliPHOSClusterizer()
{
  // dtor
         
 //Remove this from the parental task before destroying
  if(AliPHOSGetter::Instance()->PhosLoader())
    AliPHOSGetter::Instance()->PhosLoader()->CleanReconstructioner();
  delete fQADM ; 
}


