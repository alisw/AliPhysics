/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TMath.h>

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliEMCalTriggerExtraCuts.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts)

using namespace PWGJE::EMCALJetTasks;

/**
 * Constructor
 */
AliEMCalTriggerExtraCuts::AliEMCalTriggerExtraCuts() :
  fMinCrossedRowsTPC(0),
  fRequestBitmap(16)
{
}

/**
 * Apply track selection
 * \param o The object to check (must be of type AliVTrack)
 * \return True if the track is selected, false otherwise
 */
Bool_t AliEMCalTriggerExtraCuts::IsSelected(TObject *o){
  AliVTrack *rectrack(NULL);
  if(!(rectrack = dynamic_cast<AliVTrack *>(o))){
    AliError("Object not of type AliVTrack");
    return false;
  }
  Bool_t isSelected(true);

  if(fRequestBitmap.TestBitNumber(kTPCCrossedRows)){
    //test
    Float_t ncrossed = GetTPCCrossedRows(rectrack);
    AliDebug(2, Form("Crossed rows %f, min %d\n", ncrossed, fMinCrossedRowsTPC));
    if(ncrossed < fMinCrossedRowsTPC)
      isSelected = false;
    AliDebug(1, Form("Selected : %s", isSelected ? "Yes" : "No"));
  }

  if(fRequestBitmap.TestBitNumber(kTPCTrackLength)){
    Double_t tracklength = CalculateTPCTrackLength(rectrack),
             cut = 0.85*(130-5*TMath::Abs(1./rectrack->Pt()));
    AliDebug(2, Form("track length %f, min %f\n", tracklength, cut));
    if(tracklength < cut)
      isSelected = false;
    AliDebug(1, Form("Selected : %s", isSelected ? "Yes" : "No"));
  }

  return isSelected;
}

/**
 * Temporary, until interface is unified
 * @return
 */
Float_t AliEMCalTriggerExtraCuts::GetTPCCrossedRows(const AliVTrack *const trk) const {
  if(trk->IsA() == AliESDtrack::Class()){
    return (static_cast<const AliESDtrack *>(trk))->GetTPCCrossedRows();
  } else if(trk->IsA() == AliAODTrack::Class()){
    return static_cast<Float_t>((static_cast<const AliAODTrack * >(trk))->GetTPCNCrossedRows());
  }
  return 0.;
}

/**
 * Virtual implementation of the calculation of the track length in the active TPC volume
 * Code provided by Philipp Luettig
 *
 * \param trk Track to check
 * \return Track length in the active TPC volume
 */
Double_t AliEMCalTriggerExtraCuts::CalculateTPCTrackLength(AliVTrack *trk) const{
  Double_t result = 0;
  Short_t sign = trk->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[21];

  memset(cv, 0, sizeof(Double_t) * 21);
  memset(pxpypz, 0, sizeof(Double_t) * 50);
  memset(xyz, 0, sizeof(Double_t) * 50);

  Double_t bMagZ = trk->GetEvent()->GetMagneticField();

  trk->GetXYZ(xyz);
  trk->GetPxPyPz(pxpypz);
  trk->GetCovarianceXYZPxPyPz(cv);

  AliExternalTrackParam par(xyz, pxpypz, cv, sign);
  result = AliESDtrack::GetLengthInActiveZone(&par,3,236, bMagZ ,0,0);
  return result;
}
