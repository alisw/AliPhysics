//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICHChamber.h"
#include <TRotMatrix.h>

ClassImp(AliRICHChamber)	
//______________________________________________________________________________
AliRICHChamber::AliRICHChamber(Int_t iChamber):TNamed()
{
//main ctor. Defines all geometry parameters for the given module.
  SetToZenith();//put to up position   
  switch(iChamber){
    case 1:
      RotateX(-AliRICHParam::AngleYZ());   
      RotateZ(-AliRICHParam::AngleXY());      
      fName="RICHc1";fTitle="RICH chamber 1";
      break;      
    case 2:
      RotateZ(-AliRICHParam::AngleXY());      
      fName="RICHc2";fTitle="RICH chamber 2";
      break;      
    case 3:
      RotateX(-AliRICHParam::AngleYZ());
      fName="RICHc3";fTitle="RICH chamber 3";
      break;      
    case 4:          
      fName="RICHc4";fTitle="RICH chamber 4";  //no turns
      break;      
    case 5:
      RotateX(AliRICHParam::AngleYZ());
      fName="RICHc5";fTitle="RICH chamber 5";
      break;      
    case 6:
      RotateZ(AliRICHParam::AngleXY());      
      fName="RICHc6";fTitle="RICH chamber 6";
      break;      
    case 7:
      RotateX(AliRICHParam::AngleYZ());            
      RotateZ(AliRICHParam::AngleXY());      
      fName="RICHc7";fTitle="RICH chamber 7";
      break;      
    default:
      Fatal("named ctor","Wrong chamber number %i, check CreateChamber ctor",iChamber);
  }//switch(iModuleN)
  RotateZ(AliRICHParam::AngleRot());//apply common rotation  
  fpRotMatrix=new TRotMatrix("rot"+fName,"rot"+fName, Rot().ThetaX()*TMath::RadToDeg(), Rot().PhiX()*TMath::RadToDeg(),
                                                      Rot().ThetaY()*TMath::RadToDeg(), Rot().PhiY()*TMath::RadToDeg(),
                                                      Rot().ThetaZ()*TMath::RadToDeg(), Rot().PhiZ()*TMath::RadToDeg());
}//main ctor
//__________________________________________________________________________________________________
void AliRICHChamber::Print(Option_t *opt) const
{
// Debug printout
//  printf("%s r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f  %6.2f,%6.2f %6.2f,%6.2f %6.2f,%6.2f\n",fName.Data(),
//                     Rho(), ThetaD(),PhiD(),   X(),    Y(),    Z(),
//                     ThetaXd(),PhiXd(),ThetaYd(),PhiYd(),ThetaZd(),PhiZd());
  fCenterV3.Print(opt);fPcX3.Print(opt);
}//Print()
//__________________________________________________________________________________________________
