/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ZDC reconstruction parameters                                  //
// Origin: Chiara.Oppedisano@to.infn.it                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#include <TF1.h>
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParamPbPb.h"

ClassImp(AliZDCRecoParamPbPb)

//_____________________________________________________________________________
AliZDCRecoParamPbPb::AliZDCRecoParamPbPb() :
  AliZDCRecoParam(),
  fZNCen(0), 
  fZNPer(0), 
  fZPCen(0), 
  fZPPer(0), 
  fZDCCen(0),
  fZDCPer(0),
  fbCen(0),  
  fbPer(0),  
  fZEMn(0),  
  fZEMp(0),  
  fZEMsp(0), 
  fZEMb(0),  
  fZEMEndValue(0),
  fZEMCutFraction(0),
  fDZEMSup(0),
  fDZEMInf(0),
  fEZN1MaxValue(0),
  fEZP1MaxValue(0),
  fEZDC1MaxValue(0),
  fEZN2MaxValue(0),
  fEZP2MaxValue(0),
  fEZDC2MaxValue(0)
{
  //
  //Default constructor
}
//_____________________________________________________________________________
AliZDCRecoParamPbPb::~AliZDCRecoParamPbPb()
{
// destructor

  delete fZNCen;
  delete fZNPer;
  delete fZPCen;
  delete fZPPer;
  delete fZDCCen;
  delete fZDCPer;
  delete fbCen;
  delete fbPer;
  delete fZEMn;
  delete fZEMp;
  delete fZEMsp;
  delete fZEMb;

}

//_____________________________________________________________________________
AliZDCRecoParamPbPb *AliZDCRecoParamPbPb::GetPbPbRecoParam()
{
  //
  // Makes default reconstruction parameters for low flux environment
  //
  AliZDCRecoParamPbPb *param = new AliZDCRecoParamPbPb();
 
  param->SetfZNCen("(-2.287920+sqrt(2.287920*2.287920-4*(-0.007629)*(11.921710-x)))/(2*(-0.007629))",0.,164.); 
  param->SetfZNPer("(-37.812280-sqrt(37.812280*37.812280-4*(-0.190932)*(-1709.249672-x)))/(2*(-0.190932))",0.,164.); 
  param->SetfZPCen("(-1.321353+sqrt(1.321353*1.321353-4*(-0.007283)*(3.550697-x)))/(2*(-0.007283))",0.,60.); 
  param->SetfZPPer("(-42.643308-sqrt(42.643308*42.643308-4*(-0.310786)*(-1402.945615-x)))/(2*(-0.310786))",0.,60.); 
  param->SetfZDCCen("(-1.934991+sqrt(1.934991*1.934991-4*(-0.004080)*(15.111124-x)))/(2*(-0.004080))",0.,225.); 
  param->SetfZDCPer("(-34.380639-sqrt(34.380639*34.380639-4*(-0.104251)*(-2612.189017-x)))/(2*(-0.104251))",0.,225.);
  param->SetfbCen("-0.056923+0.079703*x-0.0004301*x*x+0.000001366*x*x*x",0.,220.);
  param->SetfbPer("17.943998-0.046846*x+0.000074*x*x",0.,220.);
  param->SetfZEMn("121.7-0.1934*x+0.00007565*x*x",0.,1200.);
  param->SetfZEMp("80.05-0.1315*x+0.00005327*x*x",0.,1200.);
  param->SetfZEMsp("201.7-0.325*x+0.0001292*x*x",0.,1200.); 
  param->SetfZEMb("13.83-0.02851*x+5.101e-5*x*x-7.305e-8*x*x*x+5.101e-11*x*x*x*x-1.25e-14*x*x*x*x*x",0.,1200.);
  
  param->SetZEMEndValue(1200.);
  param->SetZEMCutFraction(0.1);
  param->SetDZEMSup(0.04);
  param->SetDZEMInf(0.05);
  param->SetEZN1MaxValue(161.);
  param->SetEZP1MaxValue(59.);
  param->SetEZDC1MaxValue(220.);
  param->SetEZN2MaxValue(161.);
  param->SetEZP2MaxValue(59.);
  param->SetEZDC2MaxValue(161.);
  
  return param;

}

//_____________________________________________________________________________
void AliZDCRecoParamPbPb::PrintParameters() const 
{
  //
  // print reconstruction parameters
  //
  printf("\n\n\t AliZDCRecoParamPbPb -> parameters set for reconstruction\n");
  printf("\t Functions for reconstruction of centrality varibles (Pb-Pb):\n");
  
  fZNCen->Print(""); 
  fZNPer->Print(""); 
  fZPCen->Print(""); 
  fZPPer->Print(""); 
  fZDCCen->Print("");
  fZDCPer->Print("");
  fbCen->Print("");  
  fbPer->Print("");  
  fZEMn->Print("");  
  fZEMp->Print("");  
  fZEMsp->Print(""); 
  fZEMb->Print("");  
  
  printf("\n ####### Parameters from EZDC vs. ZEM correlation #######  \n");
  printf("\tZEMEndPoint = %1.2f, ZEMCutFraction = %1.2f \n"
    "  DZEMInf = %1.2f, DZEMSup = %1.2f\n",
    fZEMEndValue, fZEMCutFraction, fDZEMInf, fDZEMSup);
  printf("\n ####### Parameters from EZDC vs. Nspec correlation #######        \n");
  printf("\tEZN1MaxValue = %1.2f, EZP1MaxValue = %1.2f, EZDC1MaxValue = %1.2f \n"
    "\tEZN2MaxValue = %1.2f, EZP2MaxValue = %1.2f, EZDC2MaxValue = %1.2f \n\n",
    fEZN1MaxValue, fEZP1MaxValue, fEZDC1MaxValue,
    fEZN2MaxValue, fEZP2MaxValue, fEZDC2MaxValue);
}
