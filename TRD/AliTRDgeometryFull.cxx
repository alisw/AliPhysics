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

/*
$Log$
Revision 1.13  2004/02/27 15:40:18  cblume
Fix definition of rotation matrices

Revision 1.12  2003/09/18 09:06:07  cblume
Geometry update, Removal of compiler warnings

Revision 1.10  2002/11/21 22:38:47  alibrary
Removing AliMC and AliMCProcess

Revision 1.9  2002/10/31 17:45:35  cblume
New chamber geometry

Revision 1.8  2002/02/11 14:21:16  cblume
Update of the geometry. Get rid of MANY

Revision 1.7  2001/05/11 07:56:12  hristov
Consistent declarations needed on Alpha

Revision 1.6  2001/02/14 18:22:26  cblume
Change in the geometry of the padplane

Revision 1.5  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.1.4.6  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.1.4.5  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.4  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.3  2000/09/22 14:43:41  cblume
Allow the pad/timebin-dimensions to be changed after initialization

Revision 1.4  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.2  2000/05/08 14:46:44  cblume
Include options SetPHOShole() and SetRICHhole()

Revision 1.1.4.1  2000/04/27 12:46:04  cblume
Corrected bug in full geometry

Revision 1.1  2000/02/28 19:01:15  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry for the spaceframe without holes                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TVirtualMC.h"

#include "AliTRDgeometryFull.h"
#include "AliTRDparameter.h"

ClassImp(AliTRDgeometryFull)

//_____________________________________________________________________________
AliTRDgeometryFull::AliTRDgeometryFull():AliTRDgeometry()
{
  //
  // AliTRDgeometryFull default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometryFull::~AliTRDgeometryFull()
{
  //
  // AliTRDgeometryFull destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometryFull::Init()
{
  //
  // Initializes the geometry parameter
  //

  Int_t icham;
  Int_t iplan;

  fPHOShole = kFALSE;
  fRICHhole = kFALSE;

  // The outer lengths of the chambers for the sectors with holes for the PHOS
  Float_t lengthPH[kNplan][kNcham] = { { 124.0, 117.0,   0.0, 117.0, 124.0 }
				     , { 131.0, 124.0,   0.0, 124.0, 131.0 }
				     , { 138.0, 131.0,   0.0, 131.0, 138.0 }
				     , { 145.0, 138.0,   0.0, 138.0, 145.0 }
				     , { 147.0, 140.0,   0.0, 140.0, 147.0 }
				     , { 147.0, 140.0,   0.0, 140.0, 147.0 } };

  // The outer lengths of the chambers for the sectors with holes for the RICH
  Float_t lengthRH[kNplan][kNcham] = { {  87.5,   0.0,   0.0,   0.0,  87.5 }
				     , { 101.5,   0.0,   0.0,   0.0, 101.5 }
				     , { 115.5,   0.0,   0.0,   0.0, 115.5 }
				     , { 129.5,   0.0,   0.0,   0.0, 129.5 }
				     , { 133.5,   0.0,   0.0,   0.0, 133.5 }
				     , { 133.5,   0.0,   0.0,   0.0, 133.5 } };

  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      fClengthPH[iplan][icham] = lengthPH[iplan][icham];
      fClengthRH[iplan][icham] = lengthRH[iplan][icham];
    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometryFull::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the TRD geometry without hole
  //
  //
  // Names of the TRD volumina (xx = detector number):
  //
  //      Lower part of the readout chambers (gas volume + radiator)
  //
  //        UAxx    Aluminum frames             (Al)
  //        UBxx    G10 frames                  (C)
  //        UCxx    Inner volumes               (Air)
  //
  //      Upper part of the readout chambers (readout plane + fee)
  //
  //        UDxx    G10 frames                  (C)
  //        UExx    Inner volumes of the G10    (Air)
  //        UFxx    Aluminum frames             (Al)
  //        UGxx    Inner volumes of the Al     (Air)
  //
  //      Inner material layers
  //
  //        UHxx    Radiator                    (Rohacell)
  //        UIxx    Entrance window             (Mylar)
  //        UJxx    Drift volume                (Xe/CO2)
  //        UKxx    Amplification volume        (Xe/CO2)
  //        ULxx    Pad plane                   (Cu)
  //        UMxx    Support structure           (Rohacell)
  //

  const Int_t kNdet    = kNplan * kNcham;

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;

  Float_t xpos, ypos, zpos;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Char_t  cTagV[5];
  Char_t  cTagM[5];

  AliTRDgeometry::CreateGeometry(idtmed);

  // The TRD mother volume for one sector (Air), full length in z-direction
  // Provides material for side plates of super module
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("UTR1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  // The TRD mother volume for one sector (Al), leaving hole for PHOS
  if (fPHOShole) {
    gMC->Gsvolu("UTR2","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }
  // The TRD mother volume for one sector (Al), leaving hole for RICH
  if (fRICHhole) {
    gMC->Gsvolu("UTR3","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }  

  // 
  // The side plates of the super module (Al)
  parTrd[0] = fgkSwidth1/2. - fgkSMgapT;
  parTrd[1] = fgkSwidth2/2. - fgkSMgapT;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("UTS1","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  // The TRD mother volume for one sector (Al), leaving hole for PHOS
  if (fPHOShole) {
    gMC->Gsvolu("UTS2","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  }
  // The TRD mother volume for one sector (Al), leaving hole for RICH
  if (fRICHhole) {
    gMC->Gsvolu("UTS3","TRD1",idtmed[1301-1],parTrd,kNparTrd);
  }  

  // The inner part of the TRD mother volume for one sector (Air), 
  // full length in z-direction
  parTrd[0] = fgkSwidth1/2. - fgkSMgapT - fgkSMpltT;
  parTrd[1] = fgkSwidth2/2. - fgkSMgapT - fgkSMpltT;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("UTI1","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  // The TRD mother volume for one sector (Air), leaving hole for PHOS
  if (fPHOShole) {
    gMC->Gsvolu("UTI2","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }
  // The TRD mother volume for one sector (Air), leaving hole for RICH
  if (fRICHhole) {
    gMC->Gsvolu("UTI3","TRD1",idtmed[1302-1],parTrd,kNparTrd);
  }  

  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {  

      Int_t iDet = GetDetectorSec(iplan,icham);

      // The lower part of the readout chambers (gas volume + radiator) 
      // The aluminum frames 
      sprintf(cTagV,"UA%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCraH/2. + fgkCdrH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The G10 frames 
      sprintf(cTagV,"UB%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. - fgkCalT; 
      parCha[1] = -1.;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part (air)
      sprintf(cTagV,"UC%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. - fgkCalT - fgkCclsT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCclfT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          // The aluminum frames 
          sprintf(cTagV,"UA%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2.;
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCraH/2. + fgkCdrH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
          // The G10 frames 
          sprintf(cTagV,"UB%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. - fgkCalT; 
          parCha[1] = -1.;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
          // The inner part (air)
          sprintf(cTagV,"UC%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. - fgkCalT - fgkCclsT; 
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.- fgkCclfT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          // The aluminum frames 
          sprintf(cTagV,"UA%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2.;
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCraH/2. + fgkCdrH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
          // The G10 frames 
          sprintf(cTagV,"UB%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. - fgkCalT; 
          parCha[1] = -1.;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
          // The inner part (air)
          sprintf(cTagV,"UC%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. - fgkCalT - fgkCclsT; 
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.- fgkCclfT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
	}
      }

      // The upper part of the readout chambers (readout plane)
      // The G10 frames
      sprintf(cTagV,"UD%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCamH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part of the G10 frame (air)
      sprintf(cTagV,"UE%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCcuT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCcuT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
      // The aluminum frames
      sprintf(cTagV,"UF%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCroH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The inner part of the aluminum frames
      sprintf(cTagV,"UG%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCauT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCauT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          sprintf(cTagV,"UD%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW;
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCamH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
          // The inner part of the G10 frame (air)
          sprintf(cTagV,"UE%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCcuT; 
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.- fgkCcuT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
          // The aluminum frames
          sprintf(cTagV,"UF%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW;
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCroH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
          // The inner part of the aluminum frames
          sprintf(cTagV,"UG%02d",iDet+kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCauT; 
          parCha[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.- fgkCauT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          sprintf(cTagV,"UD%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW;
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCamH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
          // The inner part of the G10 frame (air)
          sprintf(cTagV,"UE%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCcuT; 
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.- fgkCcuT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
          // The aluminum frames
          sprintf(cTagV,"UF%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW;
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.;
          parCha[2] = fgkCroH/2.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
          // The inner part of the aluminum frames
          sprintf(cTagV,"UG%02d",iDet+2*kNdet);
          parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCauT; 
          parCha[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.- fgkCauT;
          parCha[2] = -1.;
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
	}
      }

      // The material layers inside the chambers
      parCha[0] = -1.;
      parCha[1] = -1.;
      // Rohacell layer (radiator)
      parCha[2] = fgkRaThick/2;
      sprintf(cTagV,"UH%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
      // Mylar layer (entrance window + HV cathode) 
      parCha[2] = fgkMyThick/2;
      sprintf(cTagV,"UI%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1308-1],parCha,kNparCha);
      // Xe/Isobutane layer (drift volume) 
      parCha[2] = fgkDrThick/2.;
      sprintf(cTagV,"UJ%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);
      // Xe/Isobutane layer (amplification volume)
      parCha[2] = fgkAmThick/2.;
      sprintf(cTagV,"UK%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
      // Cu layer (pad plane)
      parCha[2] = fgkCuThick/2;
      sprintf(cTagV,"UL%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
      // G10 layer (support structure / honeycomb)
      parCha[2] = fgkSuThick/2;
      sprintf(cTagV,"UM%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          // Rohacell layer (radiator)
          parCha[2] = fgkRaThick/2;
          sprintf(cTagV,"UH%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
          // Mylar layer (entrance window + HV cathode) 
          parCha[2] = fgkMyThick/2;
          sprintf(cTagV,"UI%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1308-1],parCha,kNparCha);
          // Xe/Isobutane layer (drift volume) 
          parCha[2] = fgkDrThick/2.;
          sprintf(cTagV,"UJ%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);
          // Xe/Isobutane layer (amplification volume)
          parCha[2] = fgkAmThick/2.;
          sprintf(cTagV,"UK%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
          // Cu layer (pad plane)
          parCha[2] = fgkCuThick/2;
          sprintf(cTagV,"UL%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
          // G10 layer (support structure / honeycomb)
          parCha[2] = fgkSuThick/2;
          sprintf(cTagV,"UM%02d",iDet+kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          // Rohacell layer (radiator)
          parCha[2] = fgkRaThick/2;
          sprintf(cTagV,"UH%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
          // Mylar layer (entrance window + HV cathode) 
          parCha[2] = fgkMyThick/2;
          sprintf(cTagV,"UI%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1308-1],parCha,kNparCha);
          // Xe/Isobutane layer (drift volume) 
          parCha[2] = fgkDrThick/2.;
          sprintf(cTagV,"UJ%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);
          // Xe/Isobutane layer (amplification volume)
          parCha[2] = fgkAmThick/2.;
          sprintf(cTagV,"UK%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
          // Cu layer (pad plane)
          parCha[2] = fgkCuThick/2;
          sprintf(cTagV,"UL%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
          // G10 layer (support structure / honeycomb)
          parCha[2] = fgkSuThick/2;
          sprintf(cTagV,"UM%02d",iDet+2*kNdet);
          gMC->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
	}
      }

      // Position the layers in the chambers
      xpos = 0;
      ypos = 0;
      // Lower part
      // Rohacell layer (radiator)
      zpos = fgkRaZpos;
      sprintf(cTagV,"UH%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Mylar layer (entrance window + HV cathode)   
      zpos = fgkMyZpos;
      sprintf(cTagV,"UI%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Xe/Isobutane layer (drift volume) 
      zpos = fgkDrZpos;
      sprintf(cTagV,"UJ%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Upper part
      // Xe/Isobutane layer (amplification volume)
      zpos = fgkAmZpos;
      sprintf(cTagV,"UK%02d",iDet);
      sprintf(cTagM,"UE%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Readout part
      // Cu layer (pad plane)
      zpos = fgkCuZpos; 
      sprintf(cTagV,"UL%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10 layer (support structure)
      zpos = fgkSuZpos;
      sprintf(cTagV,"UM%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          // Lower part
          // Rohacell layer (radiator)
          zpos = fgkRaZpos;
          sprintf(cTagV,"UH%02d",iDet+kNdet);
          sprintf(cTagM,"UC%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Mylar layer (entrance window + HV cathode)   
          zpos = fgkMyZpos;
          sprintf(cTagV,"UI%02d",iDet+kNdet);
          sprintf(cTagM,"UC%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Xe/Isobutane layer (drift volume) 
          zpos = fgkDrZpos;
          sprintf(cTagV,"UJ%02d",iDet+kNdet);
          sprintf(cTagM,"UC%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Upper part
          // Xe/Isobutane layer (amplification volume)
          zpos = fgkAmZpos;
          sprintf(cTagV,"UK%02d",iDet+kNdet);
          sprintf(cTagM,"UE%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Readout part
          // Cu layer (pad plane)
          zpos = fgkCuZpos; 
          sprintf(cTagV,"UL%02d",iDet+kNdet);
          sprintf(cTagM,"UG%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // G10 layer (support structure)
          zpos = fgkSuZpos;
          sprintf(cTagV,"UM%02d",iDet+kNdet);
          sprintf(cTagM,"UG%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          // Lower part
          // Rohacell layer (radiator)
          zpos = fgkRaZpos;
          sprintf(cTagV,"UH%02d",iDet+2*kNdet);
          sprintf(cTagM,"UC%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Mylar layer (entrance window + HV cathode)   
          zpos = fgkMyZpos;
          sprintf(cTagV,"UI%02d",iDet+2*kNdet);
          sprintf(cTagM,"UC%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Xe/Isobutane layer (drift volume) 
          zpos = fgkDrZpos;
          sprintf(cTagV,"UJ%02d",iDet+2*kNdet);
          sprintf(cTagM,"UC%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Upper part
          // Xe/Isobutane layer (amplification volume)
          zpos = fgkAmZpos;
          sprintf(cTagV,"UK%02d",iDet+2*kNdet);
          sprintf(cTagM,"UE%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // Readout part
          // Cu layer (pad plane)
          zpos = fgkCuZpos; 
          sprintf(cTagV,"UL%02d",iDet+2*kNdet);
          sprintf(cTagM,"UG%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // G10 layer (support structure)
          zpos = fgkSuZpos;
          sprintf(cTagV,"UM%02d",iDet+2*kNdet);
          sprintf(cTagM,"UG%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
	}
      }

      // Position the inner volumes of the chambers in the frames
      xpos      = 0.0;
      ypos      = 0.0;
      zpos      = 0.0;
      // The inside of the lower G10 frame
      sprintf(cTagV,"UC%02d",iDet);
      sprintf(cTagM,"UB%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The lower G10 frame inside the aluminum frame
      sprintf(cTagV,"UB%02d",iDet);
      sprintf(cTagM,"UA%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper G10 frame
      sprintf(cTagV,"UE%02d",iDet);
      sprintf(cTagM,"UD%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper aluminum frame
      sprintf(cTagV,"UG%02d",iDet);
      sprintf(cTagM,"UF%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          // The inside of the lower G10 frame
          sprintf(cTagV,"UC%02d",iDet+kNdet);
          sprintf(cTagM,"UB%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The lower G10 frame inside the aluminum frame
          sprintf(cTagV,"UB%02d",iDet+kNdet);
          sprintf(cTagM,"UA%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The inside of the upper G10 frame
          sprintf(cTagV,"UE%02d",iDet+kNdet);
          sprintf(cTagM,"UD%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The inside of the upper aluminum frame
          sprintf(cTagV,"UG%02d",iDet+kNdet);
          sprintf(cTagM,"UF%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          // The inside of the lower G10 frame
          sprintf(cTagV,"UC%02d",iDet+2*kNdet);
          sprintf(cTagM,"UB%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The lower G10 frame inside the aluminum frame
          sprintf(cTagV,"UB%02d",iDet+2*kNdet);
          sprintf(cTagM,"UA%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The inside of the upper G10 frame
          sprintf(cTagV,"UE%02d",iDet+2*kNdet);
          sprintf(cTagM,"UD%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
          // The inside of the upper aluminum frame
          sprintf(cTagV,"UG%02d",iDet+2*kNdet);
          sprintf(cTagM,"UF%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      
	}
      }

      // Position the frames of the chambers in the TRD mother volume
      xpos  = 0.;
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
      for (Int_t ic = 0; ic < icham; ic++) {
        ypos += fClength[iplan][ic];        
      }
      ypos += fClength[iplan][icham]/2.;
      zpos  = fgkCraH/2. + fgkCdrH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
      // The lower aluminum frame, radiator + drift region
      sprintf(cTagV,"UA%02d",iDet);
      gMC->Gspos(cTagV,1,"UTI1",xpos,ypos,zpos,0,"ONLY");
      // The upper G10 frame, amplification region
      sprintf(cTagV,"UD%02d",iDet);
      zpos += fgkCamH/2. + fgkCraH/2. + fgkCdrH/2.;
      gMC->Gspos(cTagV,1,"UTI1",xpos,ypos,zpos,0,"ONLY");
      // The upper aluminum frame
      sprintf(cTagV,"UF%02d",iDet);
      zpos += fgkCroH/2. + fgkCamH/2.;
      gMC->Gspos(cTagV,1,"UTI1",xpos,ypos,zpos,0,"ONLY");
      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          xpos  = 0.;
          ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
          for (Int_t ic = 0; ic < icham; ic++) {
            ypos += fClength[iplan][ic];        
          }
          if (icham > 2) {
            ypos += fClength[iplan][icham];
            ypos -= fClengthPH[iplan][icham]/2.;
	  }
          else {
            ypos += fClengthPH[iplan][icham]/2.;
	  }
          zpos  = fgkCraH/2. + fgkCdrH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
          // The lower aluminum frame, radiator + drift region
          sprintf(cTagV,"UA%02d",iDet+kNdet);
          gMC->Gspos(cTagV,1,"UTI2",xpos,ypos,zpos,0,"ONLY");
          // The upper G10 frame, amplification region
          sprintf(cTagV,"UD%02d",iDet+kNdet);
          zpos += fgkCamH/2. + fgkCraH/2. + fgkCdrH/2.;
          gMC->Gspos(cTagV,1,"UTI2",xpos,ypos,zpos,0,"ONLY");
          // The upper aluminum frame
          sprintf(cTagV,"UF%02d",iDet+kNdet);
          zpos += fgkCroH/2. + fgkCamH/2.;
          gMC->Gspos(cTagV,1,"UTI2",xpos,ypos,zpos,0,"ONLY");
	}
      }
      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          xpos  = 0.;
          ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
          for (Int_t ic = 0; ic < icham; ic++) {
            ypos += fClength[iplan][ic];        
          }
          if (icham > 2) {
            ypos += fClength[iplan][icham];
            ypos -= fClengthRH[iplan][icham]/2.;
	  }
          else {
            ypos += fClengthRH[iplan][icham]/2.;
	  }
          zpos  = fgkCraH/2. + fgkCdrH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
          // The lower aluminum frame, radiator + drift region
          sprintf(cTagV,"UA%02d",iDet+2*kNdet);
          gMC->Gspos(cTagV,1,"UTI3",xpos,ypos,zpos,0,"ONLY");
          // The upper G10 frame, amplification region
          sprintf(cTagV,"UD%02d",iDet+2*kNdet);
          zpos += fgkCamH/2. + fgkCraH/2. + fgkCdrH/2.;
          gMC->Gspos(cTagV,1,"UTI3",xpos,ypos,zpos,0,"ONLY");
          // The upper aluminum frame
          sprintf(cTagV,"UF%02d",iDet+2*kNdet);
          zpos += fgkCroH/2. + fgkCamH/2.;
          gMC->Gspos(cTagV,1,"UTI3",xpos,ypos,zpos,0,"ONLY");
	}
      }

    }
  }

  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTI1",1,"UTS1",xpos,ypos,zpos,0,"ONLY");
  if (fPHOShole) {
    gMC->Gspos("UTI2",2,"UTS2",xpos,ypos,zpos,0,"ONLY");
  }
  if (fRICHhole) {
    gMC->Gspos("UTI3",3,"UTS3",xpos,ypos,zpos,0,"ONLY");
  }

  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTS1",1,"UTR1",xpos,ypos,zpos,0,"ONLY");
  if (fPHOShole) {
    gMC->Gspos("UTS2",2,"UTR2",xpos,ypos,zpos,0,"ONLY");
  }
  if (fRICHhole) {
    gMC->Gspos("UTS3",3,"UTR3",xpos,ypos,zpos,0,"ONLY");
  }

  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTR1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  if (fPHOShole) {
    gMC->Gspos("UTR2",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  }
  else {
    gMC->Gspos("UTR1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  }
  if (fRICHhole) {
    gMC->Gspos("UTR3",3,"BTR3",xpos,ypos,zpos,0,"ONLY");
  }
  else {
    gMC->Gspos("UTR1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");
  }

  // Create the volumes of the super module frame
  CreateFrame(idtmed);

  // Create the volumes of the services
  CreateServices(idtmed);

}

//_____________________________________________________________________________
void AliTRDgeometryFull::CreateFrame(Int_t *idtmed)
{
  //
  // Create the geometry of the frame of the supermodule
  //
  // Names of the TRD services volumina
  //
  //        USRL    Support rails for the chambers (Al)
  //        USxx    Support cross bars between the chambers (Al)
  //

  Int_t   iplan = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  Char_t  cTagV[5];

  //
  // The chamber support rails
  //

  const Float_t kSRLwid  = 2.0;
  const Float_t kSRLhgt  = 2.3;
  const Float_t kSRLdst  = 0.6;
  const Int_t   kNparSRL = 3;
  Float_t parSRL[kNparSRL];
  parSRL[0] = kSRLwid/2.;
  parSRL[1] = fgkSlenTR1/2.;
  parSRL[2] = kSRLhgt/2.;
  gMC->Gsvolu("USRL","BOX ",idtmed[1301-1],parSRL,kNparSRL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  for (iplan = 0; iplan < kNplan; iplan++) {
    
    xpos  = fCwidth[iplan]/2. + kSRLwid/2. + kSRLdst;
    ypos  = 0.0;
    zpos  = fgkCraH + fgkCdrH - fgkSheight/2. - kSRLhgt/2. 
          + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("USRL",iplan+1         ,"UTI1", xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("USRL",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      gMC->Gspos("USRL",iplan+1+2*kNplan,"UTI2", xpos,ypos,zpos,0,"ONLY");
      gMC->Gspos("USRL",iplan+1+3*kNplan,"UTI2",-xpos,ypos,zpos,0,"ONLY");
    }
    if (fRICHhole) {
      gMC->Gspos("USRL",iplan+1+4*kNplan,"UTI3", xpos,ypos,zpos,0,"ONLY");
      gMC->Gspos("USRL",iplan+1+5*kNplan,"UTI3",-xpos,ypos,zpos,0,"ONLY");
    }

  }

  //
  // The cross bars between the chambers
  //

  const Float_t kSCBwid  = 1.0;
  const Int_t   kNparSCB = 3;
  Float_t parSCB[kNparSCB];
  parSCB[1] = kSCBwid/2.;
  parSCB[2] = fgkCH/2.;

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  for (iplan = 0; iplan < kNplan; iplan++) {

    parSCB[0] = fCwidth[iplan]/2. + kSRLdst/2.;

    sprintf(cTagV,"US0%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  =   fgkSlenTR1/2. - kSCBwid/2.;
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }
    if (fRICHhole) {
      gMC->Gspos(cTagV,3,"UTI3", xpos,ypos,zpos,0,"ONLY");
    }

    sprintf(cTagV,"US1%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  = fClength[iplan][2]/2. + fClength[iplan][1];
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }
    if (fRICHhole) {
      ypos += fClength[iplan][0] - fClengthRH[iplan][0];
      gMC->Gspos(cTagV,3,"UTI3", xpos,ypos,zpos,0,"ONLY");
    }

    sprintf(cTagV,"US2%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  = fClength[iplan][2]/2.;
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      ypos += fClength[iplan][1] - fClengthPH[iplan][1];
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }

    sprintf(cTagV,"US3%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  = - fClength[iplan][2]/2.;
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      ypos -= fClength[iplan][3] - fClengthPH[iplan][3];
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }

    sprintf(cTagV,"US4%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  = - fClength[iplan][2]/2. - fClength[iplan][1];
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }
    if (fRICHhole) {
      ypos -= fClength[iplan][4] - fClengthRH[iplan][4];
      gMC->Gspos(cTagV,3,"UTI3", xpos,ypos,zpos,0,"ONLY");
    }

    sprintf(cTagV,"US5%01d",iplan);
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB);
    xpos  = 0.0;
    ypos  = - fgkSlenTR1/2. + kSCBwid/2.;
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
    gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY");
    if (fPHOShole) {
      gMC->Gspos(cTagV,2,"UTI2", xpos,ypos,zpos,0,"ONLY");
    }
    if (fRICHhole) {
      gMC->Gspos(cTagV,3,"UTI3", xpos,ypos,zpos,0,"ONLY");
    }

  }

}

//_____________________________________________________________________________
void AliTRDgeometryFull::CreateServices(Int_t *idtmed)
{
  //
  // Create the geometry of the services
  //
  // Names of the TRD services volumina
  //
  //        UTCL    Cooling arterias (Al)
  //        UTCW    Cooling arterias (Water)
  //        UUxx    Volumes for the services at the chambers (Air)
  //        UTPW    Power bars       (Cu)
  //        UTCP    Cooling pipes    (Al)
  //        UTCH    Cooling pipes    (Water)
  //        UTPL    Power lines      (Cu)
  //        UMCM    Readout MCMs     (G10/Cu/Si)
  //

  const Int_t kNdet = kNplan * kNcham;

  Int_t   iplan = 0;
  Int_t   icham = 0;

  Float_t xpos  = 0.0;
  Float_t ypos  = 0.0;
  Float_t zpos  = 0.0;

  Char_t  cTagV[5];

  // The rotation matrices
  const Int_t kNmatrix = 3;
  Int_t   matrix[kNmatrix];
  gMC->Matrix(matrix[0],100.0,  0.0, 90.0, 90.0, 10.0,  0.0);
  gMC->Matrix(matrix[1], 80.0,  0.0, 90.0, 90.0, 10.0,180.0);
  gMC->Matrix(matrix[2],  0.0,  0.0, 90.0, 90.0, 90.0,  0.0);

  AliTRDparameter *parameter = new AliTRDparameter("par","TRD parameter");

  //
  // The cooling arterias
  //

  // Width of the cooling arterias
  const Float_t kCOLwid  =  0.5; 
  // Height of the cooling arterias
  const Float_t kCOLhgt  =  5.5;
  // Positioning of the cooling 
  const Float_t kCOLposx =  1.6;
  const Float_t kCOLposz = -0.2;
  // Thickness of the walls of the cooling arterias
  const Float_t kCOLthk  =  0.1;
  const Int_t   kNparCOL = 3;
  Float_t parCOL[kNparCOL];
  parCOL[0]  = kCOLwid/2.;
  parCOL[1]  = fgkSlenTR1/2.;
  parCOL[2]  = kCOLhgt/2.;
  gMC->Gsvolu("UTCL","BOX ",idtmed[1324-1],parCOL,kNparCOL);
  parCOL[0] -= kCOLthk;
  parCOL[1]  = fgkSlenTR1/2.;
  parCOL[2] -= kCOLthk;
  gMC->Gsvolu("UTCW","BOX ",idtmed[1314-1],parCOL,kNparCOL);

  xpos  = 0.0;
  ypos  = 0.0;
  zpos  = 0.0;
  gMC->Gspos("UTCW",1,"UTCL", xpos,ypos,zpos,0,"ONLY");

  for (iplan = 1; iplan < kNplan; iplan++) {
    
    xpos  = fCwidth[iplan]/2. + kCOLwid/2. + kCOLposx;
    ypos  = 0.0;
    zpos  = kCOLhgt/2. - fgkSheight/2. + kCOLposz + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("UTCL",iplan+1         ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
    gMC->Gspos("UTCL",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY");
    if (fPHOShole) {
      gMC->Gspos("UTCL",iplan+1+2*kNplan,"UTI2", xpos,ypos,zpos,matrix[0],"ONLY");
      gMC->Gspos("UTCL",iplan+1+3*kNplan,"UTI2",-xpos,ypos,zpos,matrix[1],"ONLY");
    }
    if (fRICHhole) {
      gMC->Gspos("UTCL",iplan+1+4*kNplan,"UTI3", xpos,ypos,zpos,matrix[0],"ONLY");
      gMC->Gspos("UTCL",iplan+1+5*kNplan,"UTI3",-xpos,ypos,zpos,matrix[1],"ONLY");
    }

  }

  //
  // The power bars
  //

  const Float_t kPWRwid  =  0.6;
  const Float_t kPWRhgt  =  4.5;
  const Float_t kPWRposx =  1.05;
  const Float_t kPWRposz =  0.9;
  const Int_t   kNparPWR = 3;
  Float_t parPWR[kNparPWR];
  parPWR[0] = kPWRwid/2.;
  parPWR[1] = fgkSlenTR1/2.;
  parPWR[2] = kPWRhgt/2.;
  gMC->Gsvolu("UTPW","BOX ",idtmed[1325-1],parPWR,kNparPWR);

  for (iplan = 1; iplan < kNplan; iplan++) {
    
    xpos  = fCwidth[iplan]/2. + kPWRwid/2. + kPWRposx;
    ypos  = 0.0;
    zpos  = kPWRhgt/2. - fgkSheight/2. + kPWRposz + iplan * (fgkCH + fgkVspace);
    gMC->Gspos("UTPW",iplan+1         ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY");
    gMC->Gspos("UTPW",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY");
    if (fPHOShole) {
      gMC->Gspos("UTPW",iplan+1+2*kNplan,"UTI2", xpos,ypos,zpos,matrix[0],"ONLY");
      gMC->Gspos("UTPW",iplan+1+3*kNplan,"UTI2",-xpos,ypos,zpos,matrix[1],"ONLY");
    }
    if (fRICHhole) {
      gMC->Gspos("UTPW",iplan+1+4*kNplan,"UTI3", xpos,ypos,zpos,matrix[0],"ONLY");
      gMC->Gspos("UTPW",iplan+1+5*kNplan,"UTI3",-xpos,ypos,zpos,matrix[1],"ONLY");
    }

  }

  //
  // The volumes for the services at the chambers
  //

  const Int_t kNparServ = 3;
  Float_t parServ[kNparServ];

  for (icham = 0; icham < kNcham; icham++) {
    //for (iplan = 0; iplan < kNplan; iplan++) {
    // Take out upper plane until TRD mothervolume is adjusted
    for (iplan = 0; iplan < kNplan-1; iplan++) {

      Int_t iDet = GetDetectorSec(iplan,icham);

      sprintf(cTagV,"UU%02d",iDet);
      parServ[0] = fCwidth[iplan]/2.;
      parServ[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parServ[2] = fgkVspace/2.;
      gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ);
      xpos  = 0.;
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
      for (Int_t ic = 0; ic < icham; ic++) {
        ypos += fClength[iplan][ic];        
      }
      ypos += fClength[iplan][icham]/2.;
      zpos  = fgkCH + fgkVspace/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
      gMC->Gspos(cTagV,1,"UTI1",xpos,ypos,zpos,0,"ONLY");

      if (fPHOShole) {
        if (fClengthPH[iplan][icham] > 0.0) {
          sprintf(cTagV,"UU%02d",iDet+kNdet);
          parServ[0] = fCwidth[iplan]/2.;
          parServ[1] = fClengthPH[iplan][icham]/2. - fgkHspace/2.;
          parServ[2] = fgkVspace/2.;
          gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ);
          xpos  = 0.;
          ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
          for (Int_t ic = 0; ic < icham; ic++) {
            ypos += fClength[iplan][ic];        
          }
          if (icham > 2) {
            ypos += fClength[iplan][icham];
            ypos -= fClengthPH[iplan][icham]/2.;
	  }
          else {
            ypos += fClengthPH[iplan][icham]/2.;
	  }
          zpos  = fgkCH + fgkVspace/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
          gMC->Gspos(cTagV,1,"UTI2",xpos,ypos,zpos,0,"ONLY");
	}
      }

      if (fRICHhole) {
        if (fClengthRH[iplan][icham] > 0.0) {
          sprintf(cTagV,"UU%02d",iDet+2*kNdet);
          parServ[0] = fCwidth[iplan]/2.;
          parServ[1] = fClengthRH[iplan][icham]/2. - fgkHspace/2.;
          parServ[2] = fgkVspace/2.;
          gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ);
          xpos  = 0.;
          ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
          for (Int_t ic = 0; ic < icham; ic++) {
            ypos += fClength[iplan][ic];        
          }
          if (icham > 2) {
            ypos += fClength[iplan][icham];
            ypos -= fClengthRH[iplan][icham]/2.;
	  }
          else {
            ypos += fClengthRH[iplan][icham]/2.;
	  }
          zpos  = fgkCH + fgkVspace/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
          gMC->Gspos(cTagV,1,"UTI3",xpos,ypos,zpos,0,"ONLY");
	}
      }

    }
  }

  //
  // The cooling pipes inside the service volumes
  //

  const Int_t kNparTube = 3;
  Float_t parTube[kNparTube];
  // The aluminum pipe for the cooling
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  gMC->Gsvolu("UTCP","TUBE",idtmed[1324-1],parTube,0);
  // The cooling water
  parTube[0] =  0.0;
  parTube[1] =  0.2/2.;
  parTube[2] = -1.;
  gMC->Gsvolu("UTCH","TUBE",idtmed[1314-1],parTube,kNparTube);
  // Water inside the cooling pipe
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("UTCH",1,"UTCP",xpos,ypos,zpos,0,"ONLY");

  // Position the cooling pipes in the mother volume
  const Int_t kNpar = 3;
  Float_t par[kNpar];
  for (icham = 0; icham < kNcham;   icham++) {
    //for (iplan = 0; iplan < kNplan; iplan++) {
    // Take out upper plane until TRD mothervolume is adjusted
    for (iplan = 0; iplan < kNplan-1; iplan++) { 
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 100;
      Int_t   nMCMrow = parameter->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos   = 0.0;
        ypos   = (0.5 + iMCMrow) * ySize - 1.9 
               - fClength[iplan][icham]/2. + fgkHspace/2.;
        zpos   = 0.0;                 
        par[0] = 0.0;
        par[1] = 0.3/2.; // Thickness of the cooling pipes
        par[2] = fCwidth[iplan]/2.;
        gMC->Gsposp("UTCP",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",par,kNpar);
      }
      if (fPHOShole) {
        sprintf(cTagV,"UU%02d",iDet+kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          xpos   = 0.0;
          ypos   = (0.5 + iMCMrow) * ySize - 1.9 
                 - fClengthPH[iplan][icham]/2. + fgkHspace/2.;
          zpos   = 0.0;                 
          if (ypos < (fClengthPH[iplan][icham]/2. - fgkHspace/2.)) {
            par[0] = 0.0;
            par[1] = 0.3/2.; // Thickness of the cooling pipes
            par[2] = fCwidth[iplan]/2.;
            gMC->Gsposp("UTCP",iCopy+iMCMrow+nMCMrow,cTagV,xpos,ypos,zpos
                              ,matrix[2],"ONLY",par,kNpar);
	  }
	}
      }
      if (fRICHhole) {
        sprintf(cTagV,"UU%02d",iDet+2*kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          xpos   = 0.0;
          ypos   = (0.5 + iMCMrow) * ySize - 1.9 
                 - fClengthRH[iplan][icham]/2. + fgkHspace/2.;
          zpos   = 0.0;                 
          if (ypos < (fClengthRH[iplan][icham]/2. - fgkHspace/2.)) {
            par[0] = 0.0;
            par[1] = 0.3/2.; // Thickness of the cooling pipes
            par[2] = fCwidth[iplan]/2.;
            gMC->Gsposp("UTCP",iCopy+iMCMrow+2*nMCMrow,cTagV,xpos,ypos,zpos
                              ,matrix[2],"ONLY",par,kNpar);
	  }
	}
      }
    }
  }

  //
  // The power lines
  //

  // The copper power lines
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  gMC->Gsvolu("UTPL","TUBE",idtmed[1305-1],parTube,0);

  // Position the power lines in the mother volume
  for (icham = 0; icham < kNcham;   icham++) {
    //for (iplan = 0; iplan < kNplan; iplan++) {
    // Take out upper plane until TRD mothervolume is adjusted
    for (iplan = 0; iplan < kNplan-1; iplan++) { 
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 100;
      Int_t   nMCMrow = parameter->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        xpos   = 0.0;
        ypos   = (0.5 + iMCMrow) * ySize - 1.0 
               - fClength[iplan][icham]/2. + fgkHspace/2.;
        zpos   = -0.4;
        par[0] = 0.0;
        par[1] = 0.2/2.; // Thickness of the power lines
        par[2] = fCwidth[iplan]/2.;
        gMC->Gsposp("UTPL",iCopy+iMCMrow,cTagV,xpos,ypos,zpos
                          ,matrix[2],"ONLY",par,kNpar);
      }
      if (fPHOShole) {
        sprintf(cTagV,"UU%02d",iDet+kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          xpos   = 0.0;
          ypos   = (0.5 + iMCMrow) * ySize - 1.0 
                 - fClengthPH[iplan][icham]/2. + fgkHspace/2.;
          zpos   = -0.4;                 
          if (ypos < (fClengthPH[iplan][icham]/2. - fgkHspace/2.)) {
            par[0] = 0.0;
            par[1] = 0.2/2.; // Thickness of the power lines
            par[2] = fCwidth[iplan]/2.;
            gMC->Gsposp("UTPL",iCopy+iMCMrow+nMCMrow,cTagV,xpos,ypos,zpos
                              ,matrix[2],"ONLY",par,kNpar);
	  }
	}
      }
      if (fRICHhole) {
        sprintf(cTagV,"UU%02d",iDet+2*kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          xpos   = 0.0;
          ypos   = (0.5 + iMCMrow) * ySize - 1.0 
                 - fClengthRH[iplan][icham]/2. + fgkHspace/2.;
          zpos   = -0.4;                 
          if (ypos < (fClengthRH[iplan][icham]/2. - fgkHspace/2.)) {
            par[0] = 0.0;
            par[1] = 0.2/2.; // Thickness of the power lines
            par[2] = fCwidth[iplan]/2.;
            gMC->Gsposp("UTPL",iCopy+iMCMrow+2*nMCMrow,cTagV,xpos,ypos,zpos
                              ,matrix[2],"ONLY",par,kNpar);
	  }
	}
      }
    }
  }

  //
  // The MCMs
  //

  // The mother volume for the MCMs (air)
  const Int_t kNparMCM = 3;
  Float_t parMCM[kNparMCM];
  parMCM[0] = 3.0/2.;
  parMCM[1] = 3.0/2.;
  parMCM[2] = 0.14/2.;
  gMC->Gsvolu("UMCM","BOX",idtmed[1302-1],parMCM,kNparMCM);

  // The MCM carrier G10 layer
  parMCM[0] = 3.0/2.;
  parMCM[1] = 3.0/2.;
  parMCM[2] = 0.1/2.;
  gMC->Gsvolu("UMC1","BOX",idtmed[1319-1],parMCM,kNparMCM);
  // The MCM carrier Cu layer
  parMCM[0] = 3.0/2.;
  parMCM[1] = 3.0/2.;
  parMCM[2] = 0.0162/2.;
  gMC->Gsvolu("UMC2","BOX",idtmed[1318-1],parMCM,kNparMCM);
  // The silicon of the chips
  parMCM[0] = 3.0/2.;
  parMCM[1] = 3.0/2.;
  parMCM[2] = 0.003/2.;
  gMC->Gsvolu("UMC3","BOX",idtmed[1320-1],parMCM,kNparMCM);

  // Put the MCM material inside the MCM mother volume
  xpos  =  0.0;
  ypos  =  0.0;
  zpos  = -0.07      + 0.1/2.;
  gMC->Gspos("UMC1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.1/2.    + 0.0162/2.;
  gMC->Gspos("UMC2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.00162/2 + 0.003/2.;
  gMC->Gspos("UMC3",1,"UMCM",xpos,ypos,zpos,0,"ONLY");

  // Position the MCMs in the mother volume
  for (icham = 0; icham < kNcham;   icham++) {
    //for (iplan = 0; iplan < kNplan; iplan++) {
    // Take out upper plane until TRD mothervolume is adjusted
    for (iplan = 0; iplan < kNplan-1; iplan++) { 
      Int_t   iDet    = GetDetectorSec(iplan,icham);
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 1000;
      Int_t   nMCMrow = parameter->GetRowMax(iplan,icham,0);
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW) 
                      / ((Float_t) nMCMrow);
      Int_t   nMCMcol = 8;
      Float_t xSize   = (GetChamberWidth(iplan) - 2.* fgkCpadW)
	              / ((Float_t) nMCMcol);
      sprintf(cTagV,"UU%02d",iDet);
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
        for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
          xpos   = (0.5 + iMCMcol) * xSize + 1.0 
                 - fCwidth[iplan]/2.;
          ypos   = (0.5 + iMCMrow) * ySize + 1.0 
                 - fClength[iplan][icham]/2. + fgkHspace/2.;
          zpos   = -0.4;
          par[0] = 0.0;
          par[1] = 0.2/2.; // Thickness of the power lines
          par[2] = fCwidth[iplan]/2.;
          gMC->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol,cTagV
                           ,xpos,ypos,zpos,0,"ONLY");
	}
      }
      if (fPHOShole) {
        sprintf(cTagV,"UU%02d",iDet+kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
            xpos   = (0.5 + iMCMcol) * xSize + 1.0 
                   - fCwidth[iplan]/2.;
            ypos   = (0.5 + iMCMrow) * ySize + 1.0 
                   - fClengthPH[iplan][icham]/2. + fgkHspace/2.;
            zpos   = -0.4;
            if (ypos < (fClengthPH[iplan][icham]/2. - fgkHspace/2.)) {
              par[0] = 0.0;
              par[1] = 0.2/2.; // Thickness of the power lines
              par[2] = fCwidth[iplan]/2.;
              gMC->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol+10*nMCMrow,cTagV
                               ,xpos,ypos,zpos,0,"ONLY");
	    }
	  }
        }
      }
      if (fRICHhole) {
        sprintf(cTagV,"UU%02d",iDet+2*kNdet);
        for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
          for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
            xpos   = (0.5 + iMCMcol) * xSize + 1.0 
                   - fCwidth[iplan]/2.;
            ypos   = (0.5 + iMCMrow) * ySize + 1.0 
                   - fClengthRH[iplan][icham]/2. + fgkHspace/2.;
            zpos   = -0.4;
            if (ypos < (fClengthRH[iplan][icham]/2. - fgkHspace/2.)) {
              par[0] = 0.0;
              par[1] = 0.2/2.; // Thickness of the power lines
              par[2] = fCwidth[iplan]/2.;
              gMC->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol+20*nMCMrow,cTagV
                               ,xpos,ypos,zpos,0,"ONLY");
	    }
	  }
        }
      }

    }
  }

  delete parameter;

}

//_____________________________________________________________________________
void AliTRDgeometryFull::SetOldGeometry()
{
  //
  // Use the old chamber lengths
  //

  Int_t icham;
  Int_t iplan;

  AliTRDgeometry::SetOldGeometry();

  Float_t lengthPH[kNplan][kNcham] = { { 123.5, 116.5,   0.0, 116.5, 123.5 }
				     , { 131.0, 124.0,   0.0, 124.0, 131.0 }
				     , { 134.5, 131.5,   0.0, 131.5, 134.5 }
				     , { 142.0, 139.0,   0.0, 139.0, 142.0 }
				     , { 142.0, 146.0,   0.0, 146.0, 142.0 }
                                     , { 134.5, 153.5,   0.0, 153.5, 134.5 } };

  Float_t lengthRH[kNplan][kNcham] = { {  86.5,   0.0,   0.0,   0.0,  86.5 }
				     , { 101.5,   0.0,   0.0,   0.0, 101.5 }
				     , { 112.5,   0.0,   0.0,   0.0, 112.5 }
				     , { 127.5,   0.0,   0.0,   0.0, 127.5 }
				     , { 134.5,   0.0,   0.0,   0.0, 134.5 }
                                     , { 134.5,   0.0,   0.0,   0.0, 134.5 } };
                                                                               
  for (icham = 0; icham < kNcham; icham++) {
    for (iplan = 0; iplan < kNplan; iplan++) {
      fClengthPH[iplan][icham] = lengthPH[iplan][icham];
      fClengthRH[iplan][icham] = lengthRH[iplan][icham];
    }
  }

}
