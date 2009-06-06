// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************


#include "AliHLTITSDetector.h"

#include "AliITSgeomTGeo.h"
#include "AliITSCalibration.h"
#include "AliITSsegmentation.h"
#include "AliITSDetTypeRec.h"

//------------------------------------------------------------------------
AliHLTITSDetector::AliHLTITSDetector(const AliHLTITSDetector& det):
fR(det.fR),
fRmisal(det.fRmisal),
fPhi(det.fPhi),
fSinPhi(det.fSinPhi),
fCosPhi(det.fCosPhi),
fYmin(det.fYmin),
fYmax(det.fYmax),
fZmin(det.fZmin),
fZmax(det.fZmax),
fIsBad(det.fIsBad),
fNChips(det.fNChips),
fChipIsBad(det.fChipIsBad)
{
  //Copy constructor
}
//------------------------------------------------------------------------
void AliHLTITSDetector::ReadBadDetectorAndChips(Int_t ilayer,Int_t idet,
						const AliITSDetTypeRec *detTypeRec)
{
  //--------------------------------------------------------------------
  // Read bad detectors and chips from calibration objects in AliITSDetTypeRec
  //--------------------------------------------------------------------

  // In AliITSDetTypeRec, detector numbers go from 0 to 2197
  // while in the tracker they start from 0 for each layer
  for(Int_t il=0; il<ilayer; il++) 
    idet += AliITSgeomTGeo::GetNLadders(il+1)*AliITSgeomTGeo::GetNDetectors(il+1);

  Int_t detType;
  if (ilayer==0 || ilayer==1) {        // ----------  SPD
    detType = 0;
  } else if (ilayer==2 || ilayer==3) { // ----------  SDD
    detType = 1;
  } else if (ilayer==4 || ilayer==5) { // ----------  SSD
    detType = 2;
  } else {
    printf("AliITStrackerHLT::AliHLTITSDetector::InitBadFromOCDB: Wrong layer number %d\n",ilayer);
    return;
  }

  // Get calibration from AliITSDetTypeRec
  AliITSCalibration *calib = (AliITSCalibration*)detTypeRec->GetCalibrationModel(idet);
  calib->SetModuleIndex(idet);
  AliITSCalibration *calibSPDdead = 0;
  if(detType==0) calibSPDdead = (AliITSCalibration*)detTypeRec->GetSPDDeadModel(idet); // TEMPORARY
  if (calib->IsBad() ||
      (detType==0 && calibSPDdead->IsBad())) // TEMPORARY
    {
      SetBad();
      //      printf("lay %d bad %d\n",ilayer,idet);
    }

  // Get segmentation from AliITSDetTypeRec
  AliITSsegmentation *segm = (AliITSsegmentation*)detTypeRec->GetSegmentationModel(detType);

  // Read info about bad chips
  fNChips = segm->GetMaximumChipIndex()+1;
  //printf("ilayer %d  detType %d idet %d fNChips %d %d  GetNumberOfChips %d\n",ilayer,detType,idet,fNChips,segm->GetMaximumChipIndex(),segm->GetNumberOfChips());
  if(fChipIsBad) { delete [] fChipIsBad; fChipIsBad=NULL; }
  fChipIsBad = new Bool_t[fNChips];
  for (Int_t iCh=0;iCh<fNChips;iCh++) {
    fChipIsBad[iCh] = calib->IsChipBad(iCh);
    if (detType==0 && calibSPDdead->IsChipBad(iCh)) fChipIsBad[iCh] = kTRUE; // TEMPORARY
    //if(fChipIsBad[iCh]) {printf("lay %d det %d bad chip %d\n",ilayer,idet,iCh);}
  }

  return;
}


