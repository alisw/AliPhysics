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
//                Dielectron ClusterCuts                                  //
//                                                                       //
//                                                                       //
/*
 Detailed description
 
 
 */
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <TMath.h>

#include "AliDielectronClusterCuts.h"
#include "AliVCluster.h"

ClassImp(AliDielectronClusterCuts)

//______________________________________________
AliDielectronClusterCuts::AliDielectronClusterCuts() : AliAnalysisCuts(),
  fCaloType(kAny),
  fMinNCells(-999),
  fRejectExotics(kFALSE),
  fRequireTrackMatch(kFALSE),
  fM02Min(-999.),
  fM02Max(-999.),
  fM20Min(-999.),
  fM20Max(-999.),
  fTrackDxMin(-999.),
  fTrackDxMax(-999.),
  fTrackDzMin(-999.),
  fTrackDzMax(-999.)
{
  //
  // default constructor
  //
}

//______________________________________________
AliDielectronClusterCuts::AliDielectronClusterCuts(const char* name, const char* title) : AliAnalysisCuts(name, title),
  fCaloType(kAny),
  fMinNCells(-999),
  fRejectExotics(kFALSE),
  fRequireTrackMatch(kFALSE),
  fM02Min(-999.),
  fM02Max(-999.),
  fM20Min(-999.),
  fM20Max(-999.),
  fTrackDxMin(-999.),
  fTrackDxMax(-999.),
  fTrackDzMin(-999.),
  fTrackDzMax(-999.)
{
  //
  // named constructor
  //
}

//______________________________________________
AliDielectronClusterCuts::~AliDielectronClusterCuts()
{
  //
  // default Destructor
  //
}

//______________________________________________
Bool_t AliDielectronClusterCuts::IsSelected(TObject* cluster)
{
  //
  // apply configured cuts
  //
  AliVCluster *vcluster = dynamic_cast<AliVCluster*>(cluster);
  if (!vcluster) return kFALSE;
  
  Bool_t accept = kTRUE;

  // calo typ
  if (fCaloType!=kAny) {
    if (fCaloType==kEMCal)  accept *= vcluster->IsEMCAL();
    if (fCaloType==kPHOS)   accept *= vcluster->IsPHOS();
  }
  
  // min number of cells
  if (fMinNCells!=-999) accept *= (vcluster->GetNCells()>=fMinNCells);
  
  // track match
  if (fRequireTrackMatch) accept *= (vcluster->GetNTracksMatched()>0);
  
  // M02
  if (fM02Min!=-999.) accept *= (vcluster->GetM02()>=fM02Min);
  if (fM02Max!=-999.) accept *= (vcluster->GetM02()<=fM02Max);
  
  // M20
  if (fM20Min!=-999.) accept *= (vcluster->GetM20()>=fM20Min);
  if (fM20Max!=-999.) accept *= (vcluster->GetM20()<=fM20Max);

  // track Dx
  if (fTrackDxMin!=-999.) accept *= (vcluster->GetTrackDx()>=fTrackDxMin);
  if (fTrackDxMax!=-999.) accept *= (vcluster->GetTrackDx()<=fTrackDxMax);

  // track Dz
  if (fTrackDzMin!=-999.) accept *= (vcluster->GetTrackDz()>=fTrackDzMin);
  if (fTrackDzMax!=-999.) accept *= (vcluster->GetTrackDz()<=fTrackDzMax);

  // exotics
  if (fRejectExotics) accept *= !(vcluster->GetIsExotic());
  
  return accept;
}
