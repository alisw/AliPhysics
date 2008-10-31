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

#define AliFlowAnalysisWithLYZEventPlane_cxx

// AliFlowLYZEventPlane:
//
// Class to calculate the event plane and event weight from the LYZ method
// It needs input from the standard LYZ first and second run
// author: N. van der Kolk (kolk@nikhef.nl)
 
#include "Riostream.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"
#include "TList.h"

#include "AliFlowVector.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowLYZEventPlane.h"

class AliFlowTrackSimple;

ClassImp(AliFlowLYZEventPlane)

  //-----------------------------------------------------------------------
  
AliFlowLYZEventPlane::AliFlowLYZEventPlane():
  fSecondRunList(0),
  fWR(0),            
  fPsi(0),
  fSecondReDtheta(0),
  fSecondImDtheta(0),
  fFirstr0theta(0)
{
  // Constructor.
  
}
//-----------------------------------------------------------------------


AliFlowLYZEventPlane::~AliFlowLYZEventPlane() 
{
  //destructor
  delete fSecondRunList;
}
 

//-----------------------------------------------------------------------
void AliFlowLYZEventPlane::Init() 
{
  //Declare histograms & get input files
  cout<<"---Lee Yang Zeros Event Plane Method---"<<endl;

  //input histograms
  if (fSecondRunList) {
    fSecondReDtheta = (TProfile*)fSecondRunList->FindObject("Second_FlowPro_ReDtheta_LYZ");
    fSecondImDtheta = (TProfile*)fSecondRunList->FindObject("Second_FlowPro_ImDtheta_LYZ");
    fFirstr0theta   = (TProfile*)fSecondRunList->FindObject("First_FlowPro_r0theta_LYZ");

    //warnings
    if (!fSecondReDtheta) {cout<<"fSecondReDtheta is NULL!"<<endl; }
    if (!fSecondImDtheta) {cout<<"fSecondImDtheta is NULL!"<<endl; }
    if (!fFirstr0theta)   {cout<<"fFirstr0theta is NULL!"<<endl; }
  }


}

//-----------------------------------------------------------------------
void AliFlowLYZEventPlane::CalculateRPandW(AliFlowVector aQ)
{
  //declare variables
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Double_t dCosTerm = 0;
  Double_t dSinTerm = 0;
  TComplex cDtheta;
  TComplex cRatio;
  
  if (aQ.Mod()==0.) { cout<<"Q vector is NULL"<<endl; }
  else {
    for (Int_t theta=0;theta<iNtheta;theta++)	{
      Double_t dTheta = ((float)theta/iNtheta)*TMath::Pi()/2;  
      //Calculate Qtheta 
      Double_t dQtheta = aQ.X()*cos(2*dTheta)+aQ.Y()*sin(2*dTheta);  //get Qtheta from Q vector
            
      //get R0
      Double_t dR0 = fFirstr0theta->GetBinContent(theta+1); 
      //cerr<<"dR0 = "<<dR0<<endl;

      //get Dtheta
      Double_t dReDtheta = fSecondReDtheta->GetBinContent(theta+1);
      //cerr<<"dReDtheta = "<<dReDtheta<<endl;
      Double_t dImDtheta = fSecondImDtheta->GetBinContent(theta+1);
      //cerr<<"dImDtheta = "<<dImDtheta<<endl;
      cDtheta(dReDtheta,dImDtheta);
    
      TComplex cExpo(0.,dR0*dQtheta);                  //Complex number: 0 +(i r0 Qtheta)
      if (cDtheta.Rho()!=0.) { cRatio =(TComplex::Exp(cExpo))/cDtheta; } //(e^(i r0 Qtheta))/Dtheta
      else { cRatio(0.,0.); }
      
      //sum over theta
      dCosTerm += cRatio.Re() * TMath::Cos(2*dTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } cos(2 theta)  
      dSinTerm += cRatio.Re() * TMath::Sin(2*dTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } sin(2 theta)
      
    }//loop over theta

    //average over theta
    dCosTerm /= iNtheta;  
    dSinTerm /= iNtheta;
    //cerr<<"cosTerm & sinTerm are: "<<dCosTerm<<" & "<<dSinTerm<<endl;
    
    //calculate fWR
    fWR = TMath::Sqrt(dCosTerm*dCosTerm + dSinTerm*dSinTerm);
        
    //calculate fRP
    fPsi = 0.5*TMath::ATan2(dSinTerm,dCosTerm);   //takes care of the signs correctly!
    if (fPsi < 0.) { fPsi += TMath::Pi(); }       //to shift distribution from (-pi/2 to pi/2) to (0 to pi)
  }

}

