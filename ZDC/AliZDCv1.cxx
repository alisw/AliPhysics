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
Revision 1.5  2000/10/02 21:28:20  fca
Removal of useless dependecies via forward declarations

Revision 1.3.2.1  2000/08/24 09:25:47  hristov
Patch by P.Hristov: Bug in ZDC geometry corrected by E.Scomparin

Revision 1.4  2000/08/24 09:23:59  hristov
Bug in ZDC geometry corrected by E.Scomparin

Revision 1.3  2000/07/12 06:59:16  fca
Fixing dimension of hits array

Revision 1.2  2000/07/11 11:12:34  fca
Some syntax corrections for non standard HP aCC

Revision 1.1  2000/07/10 13:58:01  fca
New version of ZDC from E.Scomparin & C.Oppedisano

Revision 1.7  2000/01/19 17:17:40  fca

Revision 1.6  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter                                                  //
//  This class contains the basic functions for the ZDC                      //
//  Functions specific to one particular geometry are                        //
//  contained in the derived classes                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TNode.h>
#include <TMath.h>
#include <TSystem.h>

#include "stdio.h"
#include "AliZDCv1.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliCallf77.h"
#include "AliConst.h"
#include "AliPDG.h"
 
 
ClassImp(AliZDCv1)
 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter version 1                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliZDCv1::AliZDCv1() : AliZDC()
{
  //
  // Default constructor for Zero Degree Calorimeter
  //
  fMedSensF1  = 0;
  fMedSensF2  = 0;
  fMedSensZN  = 0;
  fMedSensZP  = 0;
  fMedSensGR  = 0;
  fMedSensZEM = 0;
  fMedSensPI  = 0;
  fNoShower   = 0;
}
 
//_____________________________________________________________________________
AliZDCv1::AliZDCv1(const char *name, const char *title)
  : AliZDC(name,title)
{
  //
  // Standard constructor for Zero Degree Calorimeter 
  //
  fMedSensF1  = 0;
  fMedSensF2  = 0;
  fMedSensZN  = 0;
  fMedSensZP  = 0;
  fMedSensGR  = 0;
  fMedSensZEM = 0;
  fMedSensPI  = 0;
  fNoShower   = 0;
}
 
//_____________________________________________________________________________
void AliZDCv1::CreateGeometry()
{
  //
  // Create the geometry for the Zero Degree Calorimeter version 1
  //* Initialize COMMON block ZDC_CGEOM
  //*

  CreateBeamLine();
  CreateZDC();
}
  
//_____________________________________________________________________________
void AliZDCv1::CreateBeamLine()
{
  
  Float_t angle;
  Float_t zq, conpar[9], elpar[3], tubpar[3];
  Int_t im1, im2;
  Float_t zd1, zd2;
  
  
  Int_t *idtmed = fIdtmed->GetArray();
  
  // -- Mother of the ZDC 
  
  conpar[0] = 0.;
  conpar[1] = 360.;
  conpar[2] = 2.;
  conpar[3] = 805.;
  conpar[4] = 0.;
  conpar[5] = 55.;
  conpar[6] = 13060.;
  conpar[7] = 0.;
  conpar[8] = 55.;
  gMC->Gsvolu("ZDC ", "PCON", idtmed[10], conpar, 9);
  gMC->Gspos("ZDC ", 1, "ALIC", 0., 0., 0., 0, "ONLY");

  // -- FIRST SECTION OF THE BEAM PIPE (from compensator dipole to 
  //    beginning of D1) 
  
  zd1 = 1921.6;
  
  tubpar[0] = 6.3/2.;
  tubpar[1] = 6.7/2.;
  tubpar[2] = 3916.7/2.;
  gMC->Gsvolu("P001", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P001", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  //-- SECOND SECTION OF THE BEAM PIPE (FROM THE END OF D1 TO THE BEGINNING OF
  //    D2) 
  
  //-- FROM MAGNETIC BEGINNING OG D1 TO MAGNETIC END OF D1 + 23.5 cm
  //-- Elliptic pipe
  
  zd1 = 6310.8-472.5;
  
  elpar[0] = 6.84/2.;
  elpar[1] = 5.86/2.;
  elpar[2] = 945./2.;
//  gMC->Gsvolu("E001", "ELTU", idtmed[5], elpar, 3);
//  gMC->Gspos("E001", 1, "ZDC ", 0., 0., elpar[2] + zd1, 0, "ONLY");
//  
  elpar[0] = 6.44/2.;
  elpar[1] = 5.46/2.;
  elpar[2] = 945./2.;
//  gMC->Gsvolu("E002", "ELTU", idtmed[10], elpar, 3);
//  gMC->Gspos("E002", 1, "E001", 0., 0., 0., 0, "ONLY");

  zd1 += 2.*elpar[2];
  
  elpar[0] = 6.84/2.;
  elpar[1] = 5.86/2.;
  elpar[2] = 13.5/2.;
//  gMC->Gsvolu("E003", "ELTU", idtmed[5], elpar, 3);
//  gMC->Gspos("E002", 1, "ZDC ", 0., 0., elpar[2] + zd1, 0, "ONLY");
  
  elpar[0] = 6.44/2.;
  elpar[1] = 5.46/2.;
  elpar[2] = 13.5/2.;
//  gMC->Gsvolu("E004", "ELTU", idtmed[10], elpar, 3);
//  gMC->Gspos("E004", 1, "E003", 0., 0., 0., 0, "ONLY");

  zd1 += 2.*elpar[2];
  
  conpar[0] = 25./2.;
  conpar[1] = 6.44/2.;
  conpar[2] = 6.84/2.;
  conpar[3] = 10./2.;
  conpar[4] = 10.4/2.;
  gMC->Gsvolu("C001", "CONE", idtmed[5], conpar, 5);
  gMC->Gspos("C001", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");

  zd1 += 2.*conpar[0];
  
  tubpar[0] = 10./2.;
  tubpar[1] = 10.4/2.;
  tubpar[2] = 50./2.;
  gMC->Gsvolu("P002", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P002", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 10./2.;
  tubpar[1] = 10.4/2.;
  tubpar[2] = 10./2.;
  gMC->Gsvolu("P003", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P003", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 10./2.;
  tubpar[1] = 10.4/2.;
  tubpar[2] = 3.16/2.;
  gMC->Gsvolu("P004", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P004", 1, "ZDC ", 0., 0., tubpar[0] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 10.0/2.;
  tubpar[1] = 10.4/2;
  tubpar[2] = 190./2.;
  gMC->Gsvolu("P005", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P005", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 30./2.;
  conpar[1] = 10./2.;
  conpar[2] = 10.4/2.;
  conpar[3] = 20.6/2.;
  conpar[4] = 21./2.;
  gMC->Gsvolu("P006", "CONE", idtmed[5], conpar, 5);
  gMC->Gspos("P006", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 20.6/2.;
  tubpar[1] = 21./2.;
  tubpar[2] = 450./2.;
  gMC->Gsvolu("P007", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P007", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 13.6/2.;
  conpar[1] = 20.6/2.;
  conpar[2] = 21./2.;
  conpar[3] = 25.4/2.;
  conpar[4] = 25.8/2.;
  gMC->Gsvolu("P008", "CONE", idtmed[5], conpar, 5);
  gMC->Gspos("P008", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 25.4/2.;
  tubpar[1] = 25.8/2.;
  tubpar[2] = 205.8/2.;
  gMC->Gsvolu("P009", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P009", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 50./2.;
  tubpar[1] = 50.4/2.;
  tubpar[2] = 505.4/2.;
  gMC->Gsvolu("P010", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P010", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 50./2.;
  tubpar[1] = 50.4/2.;
  tubpar[2] = 700./2.;
  gMC->Gsvolu("P011", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P011", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 50./2.;
  tubpar[1] = 50.4/2.;
  tubpar[2] = 778.5/2.;
  gMC->Gsvolu("P012", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P012", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 14.18/2.;
  conpar[1] = 50./2.;
  conpar[2] = 50.4/2.;
  conpar[3] = 55./2.;
  conpar[4] = 55.4/2.;
  gMC->Gsvolu("P013", "CONE", idtmed[5], conpar, 5);
  gMC->Gspos("P013", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 55./2.;
  tubpar[1] = 55.4/2.;
  tubpar[2] = 730./2.;
  gMC->Gsvolu("P014", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P014", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 36.86/2.;
  conpar[1] = 55./2.;
  conpar[2] = 55.4/2.;
  conpar[3] = 68./2.;
  conpar[4] = 68.4/2.;
  gMC->Gsvolu("P015", "CONE", idtmed[5], conpar, 5);
  gMC->Gspos("P015", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 68./2.;
  tubpar[1] = 68.4/2.;
  tubpar[2] = 927.3/2.;
  gMC->Gsvolu("P016", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("P016", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 0./2.;
  tubpar[1] = 68.4/2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("P017", "TUBE", idtmed[8], tubpar, 3);
  gMC->Gspos("P017", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 0./2.;
  tubpar[1] = 5./2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("Q017", "TUBE", idtmed[10], tubpar, 3);
  
  //-- Position Q017 inside P017
  gMC->Gspos("Q017", 1, "P017", -7.7, 0., 0., 0, "ONLY");
  
  tubpar[0] = 0./2.;
  tubpar[1] = 7./2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("R017", "TUBE", idtmed[10], tubpar, 3);
  
  //-- Position R017 inside P017
  gMC->Gspos("R017", 1, "P017", 7.7, 0., 0., 0, "ONLY");
  
  //-- BEAM PIPE BETWEEN END OF CONICAL PIPE AND BEGINNING OF D2
  
  tubpar[0] = 5./2.;
  tubpar[1] = 5.4/2.;
  tubpar[2] = 678./2.;
  gMC->Gsvolu("P018", "TUBE", idtmed[5], tubpar, 3);

  tubpar[0] = 7./2.;
  tubpar[1] = 7.4/2.;
  tubpar[2] = 678./2.;
  gMC->Gsvolu("P019", "TUBE", idtmed[5], tubpar, 3);
  
  // -- ROTATE PIPES 

  AliMatrix(im1, 90.-0.071, 0., 90., 90., .071, 180.);
  angle = .071*kDegrad;
  gMC->Gspos("P018", 1, "ZDC ", TMath::Sin(angle) * 645. / 2. - 9.7 + 
	       TMath::Sin(angle) * 945. / 2., 0., tubpar[2] + zd1, im1, "ONLY");
  AliMatrix(im2, 90.+0.071, 0., 90., 90., .071, 0.);
  gMC->Gspos("P019", 1, "ZDC ", 9.7 - TMath::Sin(angle) * 645. / 2., 0., 
	       tubpar[2] + zd1, im2, "ONLY");
  
  // --  END OF BEAM PIPE VOLUME DEFINITION. MAGNET DEFINITION FOLLOWS 
  //     (LHC OPTICS 6) 
  
  // -- COMPENSATOR DIPOLE (MBXW) 
  //     GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 4.5;
  tubpar[2] = 340./2.;
  gMC->Gsvolu("MBXW", "TUBE", idtmed[11], tubpar, 3);
  gMC->Gspos("MBXW", 1, "ZDC ", 0., 0., tubpar[2] + 805., 0, "ONLY");
  
  // --  YOKE (IRON WITHOUT MAGNETIC FIELD) 
  
  tubpar[0] = 4.5;
  tubpar[1] = 55.;
  tubpar[2] = 340./2.;
  gMC->Gsvolu("YMBX", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("YMBX", 1, "ZDC ", 0., 0., tubpar[2] + 805., 0, "ONLY");
  
  // -- COMPENSATOR DIPOLE (MCBWA) 
  //     GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 4.5;
  tubpar[2] = 170./2.;
  gMC->Gsvolu("MCBW", "TUBE", idtmed[11], tubpar, 3);
  gMC->Gspos("MCBW", 1, "ZDC ", 0., 0., tubpar[2] + 1921.6, 0, "ONLY");
  
  // --  YOKE (IRON WITHOUT MAGNETIC FIELD) 
  
  tubpar[0] = 4.5;
  tubpar[1] = 55.;
  tubpar[2] = 170./2.;
  gMC->Gsvolu("YMCB", "TUBE", idtmed[5], tubpar, 3);
  gMC->Gspos("YMCB", 1, "ZDC ", 0., 0., tubpar[2] + 1921.6, 0, "ONLY");
  
  // -- INNER TRIPLET 
  
  zq = 2300.;
  
  // -- DEFINE MQXL AND MQX QUADRUPOLE ELEMENT 
  
  //     MQXL 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 3.5;
  tubpar[2] = 630./2.;
  gMC->Gsvolu("MQXL", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 630./2.;
  gMC->Gsvolu("YMQL", "TUBE", idtmed[5], tubpar, 3);
  
  gMC->Gspos("MQXL", 1, "ZDC ", 0., 0., tubpar[2] + zq, 0, "ONLY");
  gMC->Gspos("YMQL", 1, "ZDC ", 0., 0., tubpar[2] + zq, 0, "ONLY");
  
  gMC->Gspos("MQXL", 2, "ZDC ", 0., 0., tubpar[2] + zq + 2430., 0, "ONLY");
  gMC->Gspos("YMQL", 2, "ZDC ", 0., 0., tubpar[2] + zq + 2430., 0, "ONLY");
  
  // --  MQX 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 3.5;
  tubpar[2] = 550./2.;
  gMC->Gsvolu("MQX ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 550./2.;
  gMC->Gsvolu("YMQ ", "TUBE", idtmed[5], tubpar, 3);
  
  gMC->Gspos("MQX ", 1, "ZDC ", 0., 0., tubpar[2] + zq + 880.,  0, "ONLY");
  gMC->Gspos("YMQ ", 1, "ZDC ", 0., 0., tubpar[2] + zq + 880.,  0, "ONLY");
  
  gMC->Gspos("MQX ", 2, "ZDC ", 0., 0., tubpar[2] + zq + 1530., 0, "ONLY");
  gMC->Gspos("YMQ ", 2, "ZDC ", 0., 0., tubpar[2] + zq + 1530., 0, "ONLY");
  
  // -- SEPARATOR DIPOLE D1 
  
  zd1 = 5838.3;
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 7.5/2.;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("D1  ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 0.;
  tubpar[1] = 110./2;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("YD1 ", "TUBE", idtmed[5], tubpar, 3);
  
  gMC->Gspos("YD1 ", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  gMC->Gspos("D1  ", 1, "YD1 ", 0., 0., 0., 0, "ONLY");
  
  // -- DIPOLE D2 
  
  zd2 = 12147.6;
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 7.5/2.;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("D2  ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 0.;
  tubpar[1] = 55.;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("YD2 ", "TUBE", idtmed[5], tubpar, 3);
  
  gMC->Gspos("YD2 ", 1, "ZDC ", 0., 0., tubpar[2] + zd2, 0, "ONLY");
  
  gMC->Gspos("D2  ", 1, "YD2 ", -9.7, 0., 0., 0, "ONLY");
  gMC->Gspos("D2  ", 2, "YD2 ",  9.7, 0., 0., 0, "ONLY");
  
  // -- END OF MAGNET DEFINITION 
}
  
//_____________________________________________________________________________
void AliZDCv1::CreateZDC()
{
  
  Int_t *idtmed = fIdtmed->GetArray();
  Int_t irot1, irot2;
  Float_t DimPb[6], DimVoid[6];

  
  //-- Create calorimeters geometry
  
  //--> Neutron calorimeter (ZN) 
  
  gMC->Gsvolu("ZNEU", "BOX ", idtmed[1], fDimZN, 3); // Passive material  
  gMC->Gsvolu("ZNF1", "TUBE", idtmed[3], fFibZN, 3); // Active material
  gMC->Gsvolu("ZNF2", "TUBE", idtmed[4], fFibZN, 3); 
  gMC->Gsvolu("ZNF3", "TUBE", idtmed[4], fFibZN, 3); 
  gMC->Gsvolu("ZNF4", "TUBE", idtmed[3], fFibZN, 3); 
  gMC->Gsvolu("ZNG1", "BOX ", idtmed[12], fGrvZN, 3); // Empty grooves 
  gMC->Gsvolu("ZNG2", "BOX ", idtmed[12], fGrvZN, 3); 
  gMC->Gsvolu("ZNG3", "BOX ", idtmed[12], fGrvZN, 3); 
  gMC->Gsvolu("ZNG4", "BOX ", idtmed[12], fGrvZN, 3); 
  
  // Divide ZNEU in towers (for hits purposes) 
  
  gMC->Gsdvn("ZNTX", "ZNEU", fTowZN[0], 1); // x-tower 
  gMC->Gsdvn("ZN1 ", "ZNTX", fTowZN[1], 2); // y-tower
  
  //-- Divide ZN1 in minitowers 
  //  fDivZN[0]= NUMBER OF FIBERS PER TOWER ALONG X-AXIS, 
  //  fDivZN[1]= NUMBER OF FIBERS PER TOWER ALONG Y-AXIS
  //  (4 fibres per minitower) 
  
  gMC->Gsdvn("ZNSL", "ZN1 ", fDivZN[1], 2); // Slices 
  gMC->Gsdvn("ZNST", "ZNSL", fDivZN[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks (4 grooves per stick)
  Float_t dx = fDimZN[0] / fDivZN[0] / 4.;
  Float_t dy = fDimZN[1] / fDivZN[1] / 4.;
  
  gMC->Gspos("ZNG1", 1, "ZNST", 0.-dx, 0.+dy, 0., 0, "ONLY");
  gMC->Gspos("ZNG2", 1, "ZNST", 0.+dx, 0.+dy, 0., 0, "ONLY");
  gMC->Gspos("ZNG3", 1, "ZNST", 0.-dx, 0.-dy, 0., 0, "ONLY");
  gMC->Gspos("ZNG4", 1, "ZNST", 0.+dx, 0.-dy, 0., 0, "ONLY");
  
  // --- Position the fibers in the grooves 
  gMC->Gspos("ZNF1", 1, "ZNG1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZNF2", 1, "ZNG2", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZNF3", 1, "ZNG3", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZNF4", 1, "ZNG4", 0., 0., 0., 0, "ONLY");
  
  // --- Position the neutron calorimeter in ZDC 
  gMC->Gspos("ZNEU", 1, "ZDC ", fPosZN[0], fPosZN[1], fPosZN[2] + fDimZN[2], 0, "ONLY");
  

  //--> Proton calorimeter (ZP)  
  
  gMC->Gsvolu("ZPRO", "BOX ", idtmed[2], fDimZP, 3); // Passive material
  gMC->Gsvolu("ZPF1", "TUBE", idtmed[3], fFibZP, 3); // Active material
  gMC->Gsvolu("ZPF2", "TUBE", idtmed[4], fFibZP, 3); 
  gMC->Gsvolu("ZPF3", "TUBE", idtmed[4], fFibZP, 3); 
  gMC->Gsvolu("ZPF4", "TUBE", idtmed[3], fFibZP, 3); 
  gMC->Gsvolu("ZPG1", "BOX ", idtmed[12], fGrvZP, 3); // Empty grooves 
  gMC->Gsvolu("ZPG2", "BOX ", idtmed[12], fGrvZP, 3); 
  gMC->Gsvolu("ZPG3", "BOX ", idtmed[12], fGrvZP, 3); 
  gMC->Gsvolu("ZPG4", "BOX ", idtmed[12], fGrvZP, 3); 
    
  //-- Divide ZPRO in towers(for hits purposes) 
  
  gMC->Gsdvn("ZPTX", "ZPRO", fTowZP[0], 1); // x-tower 
  gMC->Gsdvn("ZP1 ", "ZPTX", fTowZP[1], 2); // y-tower
  
  
  //-- Divide ZP1 in minitowers 
  //  fDivZP[0]= NUMBER OF FIBERS ALONG X-AXIS PER MINITOWER, 
  //  fDivZP[1]= NUMBER OF FIBERS ALONG Y-AXIS PER MINITOWER
  //  (4 fiber per minitower) 
  
  gMC->Gsdvn("ZPSL", "ZP1 ", fDivZP[1], 2); // Slices 
  gMC->Gsdvn("ZPST", "ZPSL", fDivZP[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks (4 grooves per stick)
  dx = fDimZP[0] / fTowZP[0] / fDivZP[0] / 2.;
  dy = fDimZP[1] / fTowZP[1] / fDivZP[1] / 2.;
  
  gMC->Gspos("ZPG1", 1, "ZPST", 0.-dx, 0.+dy, 0., 0, "ONLY");
  gMC->Gspos("ZPG2", 1, "ZPST", 0.+dx, 0.+dy, 0., 0, "ONLY");
  gMC->Gspos("ZPG3", 1, "ZPST", 0.-dx, 0.-dy, 0., 0, "ONLY");
  gMC->Gspos("ZPG4", 1, "ZPST", 0.+dx, 0.-dy, 0., 0, "ONLY");
  
  // --- Position the fibers in the grooves 
  gMC->Gspos("ZPF1", 1, "ZPG1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZPF2", 1, "ZPG2", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZPF3", 1, "ZPG3", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("ZPF4", 1, "ZPG4", 0., 0., 0., 0, "ONLY");
  

  // --- Position the proton calorimeter in ZDC 
  gMC->Gspos("ZPRO", 1, "ZDC ", fPosZP[0], fPosZP[1], fPosZP[2] + fDimZP[2], 0, "ONLY");
    
  
  
  //--> EM calorimeter (ZEM)  
  
  gMC->Gsvolu("ZEM ", "PARA", idtmed[10], fDimZEM, 6);
  
  gMC->Matrix(irot1,0.,0.,90.,90.,90.,180.); 		    // Rotation matrix 1  
  gMC->Matrix(irot2,180.,0.,90.,fDimZEM[3]+90.,90.,fDimZEM[3]); // Rotation matrix 2
//  printf("irot1 = %d, irot2 = %d \n", irot1, irot2);
  
  gMC->Gsvolu("ZEMF", "TUBE", idtmed[3], fFibZEM, 3); // Active material

  gMC->Gsdvn("ZETR", "ZEM ", fDivZEM[2], 1); 	     // Tranches 
  
  DimPb[0] = fDimZEMPb;			// Lead slices 
  DimPb[1] = fDimZEM[2];
  DimPb[2] = fDimZEM[1];
  DimPb[3] = 90.-fDimZEM[3];
  DimPb[4] = 0.;
  DimPb[5] = 0.;
  gMC->Gsvolu("ZEL0", "PARA", idtmed[6], DimPb, 6);
  gMC->Gsvolu("ZEL1", "PARA", idtmed[6], DimPb, 6);
  gMC->Gsvolu("ZEL2", "PARA", idtmed[6], DimPb, 6);
  
  // --- Position the lead slices in the tranche 
  Float_t zTran = fDimZEM[0]/fDivZEM[2]; 
  Float_t zTrPb = -zTran+fDimZEMPb;
  gMC->Gspos("ZEL0", 1, "ZETR", zTrPb, 0., 0., 0, "ONLY");
  gMC->Gspos("ZEL1", 1, "ZETR", fDimZEMPb, 0., 0., 0, "ONLY");
  
  // --- Vacuum zone (to be filled with fibres)
  DimVoid[0] = (zTran-2*fDimZEMPb)/2.;
  DimVoid[1] = fDimZEM[2];
  DimVoid[2] = fDimZEM[1];
  DimVoid[3] = 90.-fDimZEM[3];
  DimVoid[4] = 0.;
  DimVoid[5] = 0.;
  gMC->Gsvolu("ZEV0", "PARA", idtmed[10], DimVoid,6);
  gMC->Gsvolu("ZEV1", "PARA", idtmed[10], DimVoid,6);
  
  // --- Divide the vacuum slice into sticks along x axis
  gMC->Gsdvn("ZES0", "ZEV0", fDivZEM[0], 3); 
  gMC->Gsdvn("ZES1", "ZEV1", fDivZEM[0], 3); 
  
  // --- Positioning the fibers into the sticks
  gMC->Gspos("ZEMF", 1,"ZES0", 0., 0., 0., irot2, "ONLY");
  gMC->Gspos("ZEMF", 1,"ZES1", 0., 0., 0., irot2, "ONLY");
  
  // --- Positioning the vacuum slice into the tranche
  Float_t DisplFib = fDimZEM[1]/fDivZEM[0];
  gMC->Gspos("ZEV0", 1,"ZETR", -DimVoid[0], 0., 0., 0, "ONLY");
  gMC->Gspos("ZEV1", 1,"ZETR", -DimVoid[0]+zTran, 0., DisplFib, 0, "ONLY");

  // --- Positioning the ZEM into the ZDC - rotation for 90 degrees  
  gMC->Gspos("ZEM ", 1,"ZDC ", fPosZEM[0], fPosZEM[1], fPosZEM[2], irot1, "ONLY");
  
  // --- Adding last slice at the end of the EM calorimeter 
  Float_t zLastSlice = fPosZEM[2]+fDimZEMPb+fDimZEM[0];
  gMC->Gspos("ZEL2", 1,"ZDC ", fPosZEM[0], fPosZEM[1], zLastSlice, irot1, "ONLY");
  
}
 
//_____________________________________________________________________________
void AliZDCv1::DrawModule()
{
  //
  // Draw a shaded view of the Zero Degree Calorimeter version 1
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ZDC ","SEEN",0);
  gMC->Gsatt("P001","SEEN",1);
  gMC->Gsatt("E001","SEEN",1);
  gMC->Gsatt("E002","SEEN",1);
  gMC->Gsatt("E003","SEEN",1);
  gMC->Gsatt("E004","SEEN",1);
  gMC->Gsatt("C001","SEEN",1);
  gMC->Gsatt("P002","SEEN",1);
  gMC->Gsatt("P003","SEEN",1);
  gMC->Gsatt("P004","SEEN",1);
  gMC->Gsatt("P005","SEEN",1);
  gMC->Gsatt("P006","SEEN",1);
  gMC->Gsatt("P007","SEEN",1);
  gMC->Gsatt("P008","SEEN",1);
  gMC->Gsatt("P009","SEEN",1);
  gMC->Gsatt("P010","SEEN",1);
  gMC->Gsatt("P011","SEEN",1);
  gMC->Gsatt("P012","SEEN",1);
  gMC->Gsatt("P013","SEEN",1);
  gMC->Gsatt("P014","SEEN",1);
  gMC->Gsatt("P015","SEEN",1);
  gMC->Gsatt("P016","SEEN",1);
  gMC->Gsatt("P017","SEEN",1);
  gMC->Gsatt("Q017","SEEN",1);
  gMC->Gsatt("R017","SEEN",1);
  gMC->Gsatt("P018","SEEN",1);
  gMC->Gsatt("P019","SEEN",1);
  gMC->Gsatt("MBXW","SEEN",1);
  gMC->Gsatt("YMBX","SEEN",1);
  gMC->Gsatt("MCBW","SEEN",1);
  gMC->Gsatt("YMCB","SEEN",1);
  gMC->Gsatt("MQXL","SEEN",1);
  gMC->Gsatt("YMQL","SEEN",1);
  gMC->Gsatt("MQX ","SEEN",1);
  gMC->Gsatt("YMQ ","SEEN",1);
  gMC->Gsatt("D1  ","SEEN",1);
  gMC->Gsatt("YD1 ","SEEN",1);
  gMC->Gsatt("D2  ","SEEN",1);
  gMC->Gsatt("YD2 ","SEEN",1);
  gMC->Gsatt("ZNEU","SEEN",0);
  gMC->Gsatt("ZNF1","SEEN",0);
  gMC->Gsatt("ZNF2","SEEN",0);
  gMC->Gsatt("ZNF3","SEEN",0);
  gMC->Gsatt("ZNF4","SEEN",0);
  gMC->Gsatt("ZNG1","SEEN",0);
  gMC->Gsatt("ZNG2","SEEN",0);
  gMC->Gsatt("ZNG3","SEEN",0);
  gMC->Gsatt("ZNG4","SEEN",0);
  gMC->Gsatt("ZNTX","SEEN",0);
  gMC->Gsatt("ZN1 ","COLO",2); 
  gMC->Gsatt("ZN1 ","SEEN",1);
  gMC->Gsatt("ZNSL","SEEN",0);
  gMC->Gsatt("ZNST","SEEN",0);
  gMC->Gsatt("ZPRO","SEEN",0);
  gMC->Gsatt("ZPF1","SEEN",0);
  gMC->Gsatt("ZPF2","SEEN",0);
  gMC->Gsatt("ZPF3","SEEN",0);
  gMC->Gsatt("ZPF4","SEEN",0);
  gMC->Gsatt("ZPG1","SEEN",0);
  gMC->Gsatt("ZPG2","SEEN",0);
  gMC->Gsatt("ZPG3","SEEN",0);
  gMC->Gsatt("ZPG4","SEEN",0);
  gMC->Gsatt("ZPTX","SEEN",0);
  gMC->Gsatt("ZP1 ","COLO",2); 
  gMC->Gsatt("ZP1 ","SEEN",1);
  gMC->Gsatt("ZPSL","SEEN",0);
  gMC->Gsatt("ZPST","SEEN",0);
  gMC->Gsatt("ZEM ","COLO",2); 
  gMC->Gsatt("ZEM ","SEEN",1);
  gMC->Gsatt("ZEMF","SEEN",0);
  gMC->Gsatt("ZETR","SEEN",0);
  gMC->Gsatt("ZEL0","SEEN",0);
  gMC->Gsatt("ZEL1","SEEN",0);
  gMC->Gsatt("ZEL2","SEEN",0);
  gMC->Gsatt("ZEV0","SEEN",0);
  gMC->Gsatt("ZEV1","SEEN",0);
  gMC->Gsatt("ZES0","SEEN",0);
  gMC->Gsatt("ZES1","SEEN",0);
  
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 100, -100, 100, 12000, 16000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 488, 220, .07, .07);
  gMC->Gdhead(1111, "Zero Degree Calorimeter Version 1");
  gMC->Gdman(18, 4, "MAN");
}

//_____________________________________________________________________________
void AliZDCv1::CreateMaterials()
{
  //
  // Create Materials for the Zero Degree Calorimeter
  //
  // Origin    : E. Scomparin 
  
  Int_t *idtmed = fIdtmed->GetArray();
  
  Float_t dens, ubuf[1], wmat[2], a[2], z[2], epsil=0.001, stmin=0.01;
  Int_t   i, isvolActive, isvol, inofld;
  Float_t fieldm = gAlice->Field()->Max();
  Float_t tmaxfd=gAlice->Field()->Max();
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t deemax=-1;
  Float_t stemax;
  
  // --- Store in UBUF r0 for nuclear radius calculation R=r0*A**1/3 

  // --- Tantalum -> ZN passive material
  ubuf[0] = 1.1;
  AliMaterial(1, "TANT", 180.95, 73., 16.65, .4, 11.9, ubuf, 1);
    
  // --- Tungsten 
//  ubuf[0] = 1.11;
//  AliMaterial(1, "TUNG", 183.85, 74., 19.3, .35, 10.3, ubuf, 1);
  
  // --- Brass (CuZn)  -> ZP passive material
  dens = 8.48;
  a[0] = 63.546;
  a[1] = 65.39;
  z[0] = 29.;
  z[1] = 30.;
  wmat[0] = .63;
  wmat[1] = .37;
  AliMixture(2, "BRASS               ", a, z, dens, 2, wmat);
  
  // --- SiO2 
  dens = 2.64;
  a[0] = 28.086;
  a[1] = 15.9994;
  z[0] = 14.;
  z[1] = 8.;
  wmat[0] = 1.;
  wmat[1] = 2.;
  AliMixture(3, "SIO2                ", a, z, dens, -2, wmat);
  
  
  // --- Lead 
  ubuf[0] = 1.12;
  AliMaterial(5, "LEAD", 207.19, 82., 11.35, .56, 18.5, ubuf, 1);

  // --- Copper 
//  ubuf[0] = 1.1;
//  AliMaterial(7, "COPP", 63.54, 29., 8.96, 1.4, 0., ubuf, 1);
  
  // --- Iron (energy loss taken into account)
  ubuf[0] = 1.1;
  AliMaterial(6, "IRON", 55.85, 26., 7.87, 1.76, 0., ubuf, 1);
  
  // --- Iron (no energy loss)
  ubuf[0] = 1.1;
  AliMaterial(7, "IRON", 55.85, 26., 7.87, 1.76, 0., ubuf, 1);
  
  // --- Vacuum (no magnetic field) 
  AliMaterial(10, "VOID", 1e-16, 1e-16, 1e-16, 1e16, 1e16, ubuf,0);
  
  // --- Vacuum (with magnetic field) 
  AliMaterial(11, "VOIM", 1e-16, 1e-16, 1e-16, 1e16, 1e16, ubuf,0);
  
  // --- Air (no magnetic field)
  AliMaterial(12, "Air    $", 14.61, 7.3, .001205, 30420., 67500., ubuf, 0);
  
  // ---  Definition of tracking media: 
  
  // --- Tantalum = 1 ; 
  // --- Brass = 2 ; 
  // --- Fibers (SiO2) = 3 ; 
  // --- Fibers (SiO2) = 4 ; 
  // --- Lead = 5 ; 
  // --- Iron (with energy loss) = 6 ; 
  // --- Iron (without energy loss) = 7 ; 
  // --- Vacuum (no field) = 10 
  // --- Vacuum (with field) = 11 
  // --- Air (no field) = 12 
  
  
  // --- Tracking media parameters 
  epsil  = .01;
  stemax = 1.;
  isvol  = 0;
  isvolActive = 1;
  inofld = 0;
  fieldm = 0.;
  
  AliMedium(1, "ZTANT", 1, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
//  AliMedium(1, "ZW", 1, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "ZBRASS", 2, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "ZSIO2", 3, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "ZQUAR", 3, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(6, "ZLEAD", 5, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
//  AliMedium(7, "ZCOPP", 7, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(5, "ZIRON", 6, isvolActive, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(8, "ZIRONN", 7, isvol, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(10, "ZVOID", 10, isvol, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(12, "ZAIR", 12, 0, inofld, fieldm, tmaxfd, stemax,deemax, epsil, stmin);
  
  fieldm = 45.;
  AliMedium(11, "ZVOIM", 11, isvol, isxfld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  
  // Thresholds for showering in the ZDCs 
  
  i = 1;
  gMC->Gstpar(idtmed[i], "CUTGAM", .001);
  gMC->Gstpar(idtmed[i], "CUTELE", .001);
  gMC->Gstpar(idtmed[i], "CUTNEU", .01);
  gMC->Gstpar(idtmed[i], "CUTHAD", .01);
  i = 2;
  gMC->Gstpar(idtmed[i], "CUTGAM", .001);
  gMC->Gstpar(idtmed[i], "CUTELE", .001);
  gMC->Gstpar(idtmed[i], "CUTNEU", .01);
  gMC->Gstpar(idtmed[i], "CUTHAD", .01);
  i = 6;
  gMC->Gstpar(idtmed[i], "CUTGAM", .001);
  gMC->Gstpar(idtmed[i], "CUTELE", .001);
  gMC->Gstpar(idtmed[i], "CUTNEU", .01);
  gMC->Gstpar(idtmed[i], "CUTHAD", .01);
  
  // Avoid too detailed showering along the beam line 
  
  i = 5;
  gMC->Gstpar(idtmed[i], "CUTGAM", .1);
  gMC->Gstpar(idtmed[i], "CUTELE", .1);
  gMC->Gstpar(idtmed[i], "CUTNEU", 1.);
  gMC->Gstpar(idtmed[i], "CUTHAD", 1.);
  
  // Avoid interaction in fibers (only energy loss allowed) 
  i = 3;
  gMC->Gstpar(idtmed[i], "DCAY", 0.);
  gMC->Gstpar(idtmed[i], "MULS", 0.);
  gMC->Gstpar(idtmed[i], "PFIS", 0.);
  gMC->Gstpar(idtmed[i], "MUNU", 0.);
  gMC->Gstpar(idtmed[i], "LOSS", 1.);
  gMC->Gstpar(idtmed[i], "PHOT", 0.);
  gMC->Gstpar(idtmed[i], "COMP", 0.);
  gMC->Gstpar(idtmed[i], "PAIR", 0.);
  gMC->Gstpar(idtmed[i], "BREM", 0.);
  gMC->Gstpar(idtmed[i], "DRAY", 0.);
  gMC->Gstpar(idtmed[i], "ANNI", 0.);
  gMC->Gstpar(idtmed[i], "HADR", 0.);
  i = 4;
  gMC->Gstpar(idtmed[i], "DCAY", 0.);
  gMC->Gstpar(idtmed[i], "MULS", 0.);
  gMC->Gstpar(idtmed[i], "PFIS", 0.);
  gMC->Gstpar(idtmed[i], "MUNU", 0.);
  gMC->Gstpar(idtmed[i], "LOSS", 1.);
  gMC->Gstpar(idtmed[i], "PHOT", 0.);
  gMC->Gstpar(idtmed[i], "COMP", 0.);
  gMC->Gstpar(idtmed[i], "PAIR", 0.);
  gMC->Gstpar(idtmed[i], "BREM", 0.);
  gMC->Gstpar(idtmed[i], "DRAY", 0.);
  gMC->Gstpar(idtmed[i], "ANNI", 0.);
  gMC->Gstpar(idtmed[i], "HADR", 0.);
  
  // Avoid interaction in void 
  i = 10;
  gMC->Gstpar(idtmed[i], "DCAY", 0.);
  gMC->Gstpar(idtmed[i], "MULS", 0.);
  gMC->Gstpar(idtmed[i], "PFIS", 0.);
  gMC->Gstpar(idtmed[i], "MUNU", 0.);
  gMC->Gstpar(idtmed[i], "LOSS", 0.);
  gMC->Gstpar(idtmed[i], "PHOT", 0.);
  gMC->Gstpar(idtmed[i], "COMP", 0.);
  gMC->Gstpar(idtmed[i], "PAIR", 0.);
  gMC->Gstpar(idtmed[i], "BREM", 0.);
  gMC->Gstpar(idtmed[i], "DRAY", 0.);
  gMC->Gstpar(idtmed[i], "ANNI", 0.);
  gMC->Gstpar(idtmed[i], "HADR", 0.);

  //
  fMedSensF1  = idtmed[3];  // Sensitive volume: fibres type 1
  fMedSensF2  = idtmed[4];  // Sensitive volume: fibres type 2
  fMedSensZN  = idtmed[1];  // Sensitive volume: ZN passive material
  fMedSensZP  = idtmed[2];  // Sensitive volume: ZP passive material
  fMedSensZEM = idtmed[6];  // Sensitive volume: ZEM passive material
  fMedSensGR  = idtmed[12]; // Sensitive volume: air into the grooves
  fMedSensPI  = idtmed[5];  // Sensitive volume: beam pipes
} 

//_____________________________________________________________________________
void AliZDCv1::Init()
{
 InitTables();
}

//_____________________________________________________________________________
void AliZDCv1::InitTables()
{
  Int_t k, j;
  //Initialize parameters for light tables and read them
  fNalfan = 90;
  fNalfap = 90;
  fNben = 18;
  fNbep = 28;
  
  char *lightfName1,*lightfName2,*lightfName3,*lightfName4,
       *lightfName5,*lightfName6,*lightfName7,*lightfName8;
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;

  lightfName1 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620362207s");
  if((fp1 = fopen(lightfName1,"r")) == NULL){
     printf("Cannot open file fp1 \n");
     return;
  }
  lightfName2 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620362208s");
  if((fp2 = fopen(lightfName2,"r")) == NULL){
     printf("Cannot open file fp2 \n");
     return;
  }  
  lightfName3 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620362209s");
  if((fp3 = fopen(lightfName3,"r")) == NULL){
     printf("Cannot open file fp3 \n");
     return;
  }
  lightfName4 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620362210s");
  if((fp4 = fopen(lightfName4,"r")) == NULL){
     printf("Cannot open file fp4 \n");
     return;
  }
//  printf(" --- Reading light tables for ZN \n");
  for(k=0; k<fNalfan; k++){
     for(j=0; j<fNben; j++){
       fscanf(fp1,"%f",&fTablen[0][k][j]);
       fscanf(fp2,"%f",&fTablen[1][k][j]);
       fscanf(fp3,"%f",&fTablen[2][k][j]);
       fscanf(fp4,"%f",&fTablen[3][k][j]);
     } 
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  
  lightfName5 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620552207s");
  if((fp5 = fopen(lightfName5,"r")) == NULL){
     printf("Cannot open file fp5 \n");
     return;
  }
  lightfName6 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620552208s");
  if((fp6 = fopen(lightfName6,"r")) == NULL){
     printf("Cannot open file fp6 \n");
     return;
  }
  lightfName7 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620552209s");
  if((fp7 = fopen(lightfName7,"r")) == NULL){
     printf("Cannot open file fp7 \n");
     return;
  }
  lightfName8 = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/light22620552210s");
  if((fp8 = fopen(lightfName8,"r")) == NULL){
     printf("Cannot open file fp8 \n");
     return;
  }
//  printf(" --- Reading light tables for ZP and ZEM \n");
  for(k=0; k<fNalfap; k++){
     for(j=0; j<fNbep; j++){
       fscanf(fp5,"%f",&fTablep[0][k][j]);
       fscanf(fp6,"%f",&fTablep[1][k][j]);
       fscanf(fp7,"%f",&fTablep[2][k][j]);
       fscanf(fp8,"%f",&fTablep[3][k][j]);
     } 
  }
  fclose(fp5);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);
}
//_____________________________________________________________________________
void AliZDCv1::StepManager()
{
  //
  // Routine called at every step in the Zero Degree Calorimeters
  //

  Int_t j;

  Int_t vol[2], ibeta=0, ialfa, ibe;
  Float_t x[3], xdet[3], destep, hits[10], m, ekin, um[3], ud[3], be, radius, out;
  TLorentzVector s, p;
  const char *knamed;


  if((gMC->GetMedium() == fMedSensZN) || (gMC->GetMedium() == fMedSensZP) ||
     (gMC->GetMedium() == fMedSensGR) || (gMC->GetMedium() == fMedSensF1) ||
     (gMC->GetMedium() == fMedSensF2) || (gMC->GetMedium() == fMedSensZEM) ||
     (gMC->GetMedium() == fMedSensPI)){
  
  // If particle interacts with beam pipe -> return
    if(gMC->GetMedium() == fMedSensPI){ 
  
  // If option NoShower is set -> StopTrack
      if(fNoShower==1) gMC->StopTrack();
      return;
    }
  
  //Particle coordinates 
    gMC->TrackPosition(s);
    for(j=0; j<=2; j++){
       x[j] = s[j];
    }
    hits[0] = x[0];
    hits[1] = x[1];
    hits[2] = x[2];

  // Determine in which ZDC the particle is
    knamed = gMC->CurrentVolName();
    if(!strncmp(knamed,"ZN",2))vol[0]=1;
    if(!strncmp(knamed,"ZP",2))vol[0]=2;
    if(!strncmp(knamed,"ZE",2))vol[0]=3;
  
  // Determine in which quadrant the particle is
    
    //Quadrant in ZN
    if(vol[0]==1){
      xdet[0] = x[0]-fPosZN[0];
      xdet[1] = x[1]-fPosZN[1];
      if((xdet[0]<=0.) && (xdet[1]>=0.)) vol[1]=1;
      if((xdet[0]>0.) && (xdet[1]>0.)) vol[1]=2;
      if((xdet[0]<0.) && (xdet[1]<0.)) vol[1]=3;
      if((xdet[0]>0.) && (xdet[1]<0.)) vol[1]=4;
    }
    
    //Quadrant in ZP
    if(vol[0]==2){
      xdet[0] = x[0]-fPosZP[0];
      xdet[1] = x[1]-fPosZP[1];
      if(xdet[0]>fDimZP[0])xdet[0]=fDimZP[0]-0.01;
      if(xdet[0]<-fDimZP[0])xdet[0]=-fDimZP[0]+0.01;
      Float_t xqZP = xdet[0]/(fDimZP[0]/2);
      for(int i=1; i<=4; i++){
         if(xqZP>=(i-3) && xqZP<(i-2)){
 	   vol[1] = i;
 	   break;
 	}
     }
    }
    
    //ZEM has only 1 quadrant
    if(vol[0] == 3){
      vol[1] = 1;
      xdet[0] = x[0]-fPosZEM[0];
      xdet[1] = x[1]-fPosZEM[1];
//      printf("x %f %f xdet %f %f\n",x[0],x[1],xdet[0],xdet[1]);
    }

    if(vol[1]>4){
    printf("\n-> Det. %d Quad. %d \n", vol[0], vol[1]);
    printf("x %f %f xdet %f %f\n",x[0],x[1],xdet[0],xdet[1]);}

  // Store impact point and kinetic energy of the ENTERING particle
    
//    Int_t Curtrack = gAlice->CurrentTrack();
//    Int_t Prim = gAlice->GetPrimary(Curtrack);
//    printf ("Primary: %d, Current Track: %d \n", Prim, Curtrack); 
    
//    if(Curtrack==Prim){
      if(gMC->IsTrackEntering()){
        //Particle energy
        gMC->TrackMomentum(p);
//   	 printf("p[0] = %f, p[1] = %f, p[2] = %f, p[3] = %f \n", 
//                 p[0], p[1], p[2], p[3]);
        hits[3] = p[3];

        // Impact point on ZDC  
        hits[4] = xdet[0];
        hits[5] = xdet[1];
	hits[6] = 0;
        hits[7] = 0;
        hits[8] = 0;
        hits[9] = 0;

//	  Int_t PcID = gMC->TrackPid();
//	  printf("Pc ID -> %d\n",PcID);
	AddHit(gAlice->CurrentTrack(), vol, hits);
	
	if(fNoShower==1){
	gMC->StopTrack();
	return;
	}
      }
//    }
             
      // Charged particles -> Energy loss
      if((destep=gMC->Edep())){
         if(gMC->IsTrackStop()){
           gMC->TrackMomentum(p);
	   m = gMC->TrackMass();
	   ekin = p[3]-m;
	   if(ekin<0.) printf("ATTENTION!!!!!!!!!!!!!!! ->	ekin = %f <0 (?)",ekin);
	   hits[9] = ekin;
	   hits[7] = 0.;
	   hits[8] = 0.;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	   }
	 else{
	   hits[9] = destep;
	   hits[7] = 0.;
	   hits[8] = 0.;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	   }
//	 printf("	-> Charged particle -> Dep. E = %f eV \n",hits[8]);
	 }
//	 printf(" \n");
  }


  // *** Light production in fibres 
  if((gMC->GetMedium() == fMedSensF1) || (gMC->GetMedium() == fMedSensF2)){

     //Select charged particles
     if((destep=gMC->Edep())){
//       printf("		-> CHARGED particle!!! \n");

       // Particle velocity
       gMC->TrackMomentum(p);
       Float_t ptot=TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
       Float_t beta =  ptot/p[3];
//       Int_t pcID = gMC->TrackPid();
//       printf("	Pc %d in quadrant %d -> beta = %f \n", pcID, vol[1], beta);
       if(beta<0.67) return;
       if((beta>=0.67) && (beta<=0.75)) ibeta = 0;
       if((beta>0.75)  && (beta<=0.85)) ibeta = 1;
       if((beta>0.85)  && (beta<=0.95)) ibeta = 2;
//       if((beta>0.95)  && (beta<=1.00)) ibeta = 3;
       if(beta>0.95) ibeta = 3;
 
       // Angle between particle trajectory and fibre axis
       // 1 -> Momentum directions
       um[0] = p[0]/ptot;
       um[1] = p[1]/ptot;
       um[2] = p[2]/ptot;
       gMC->Gmtod(um,ud,2);
       // 2 -> Angle < limit angle
       Double_t alfar = TMath::ACos(ud[2]);
       Double_t alfa = alfar*kRaddeg;
       if(alfa>=110.) return;
       ialfa = Int_t(1.+alfa/2.);
 
       // Distance between particle trajectory and fibre axis
       gMC->TrackPosition(s);
       for(j=0; j<=2; j++){
   	  x[j] = s[j];
       }
       gMC->Gmtod(x,xdet,1);
       if(TMath::Abs(ud[0])>0.00001){
         Float_t dcoeff = ud[1]/ud[0];
         be = TMath::Abs((xdet[1]-dcoeff*xdet[0])/TMath::Sqrt(dcoeff*dcoeff+1.));
       }
       else{
         be = TMath::Abs(ud[0]);
       }
 
       if((vol[0]==1)) radius = fFibZN[1];
       if((vol[0]==2)) radius = fFibZP[1];
       ibe = Int_t(be*1000.+1);
 
       //Looking into the light tables 
       Float_t charge = gMC->TrackCharge();
       
       // (1)  ZN
       if((vol[0]==1)) {
         if(ibe>fNben) ibe=fNben;
         out =  charge*charge*fTablen[ibeta][ialfa][ibe];
	 if(gMC->GetMedium() == fMedSensF1){
	   hits[7] = out;  	//fLightPMQ
	   hits[8] = 0;
	   hits[9] = 0;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	 }
	 else{
	   hits[7] = 0;
	   hits[8] = out;	//fLightPMC
	   hits[9] = 0;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	 }
       } 
       
       // (2) ZP
       if((vol[0]==2)) {
         if(ibe>fNbep) ibe=fNbep;
         out =  charge*charge*fTablep[ibeta][ialfa][ibe];
	 if(gMC->GetMedium() == fMedSensF1){
	   hits[7] = out;  	//fLightPMQ
	   hits[8] = 0;
	   hits[9] = 0;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	 }
	 else{
	   hits[7] = 0;
	   hits[8] = out;	//fLightPMC
	   hits[9] = 0;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
	 }
       } 
       // (3) ZEM
       if((vol[0]==3)) {
         if(ibe>fNbep) ibe=fNbep;
         out =  charge*charge*fTablep[ibeta][ialfa][ibe];
	   hits[7] = out;  	//fLightPMQ
	   hits[8] = 0;
	   hits[9] = 0;
	   AddHit(gAlice->CurrentTrack(), vol, hits);
       } 
     }
       
   }
}
