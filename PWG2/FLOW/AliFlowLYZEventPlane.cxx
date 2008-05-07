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
 
#include "Riostream.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"

#include "AliFlowLYZConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowLYZEventPlane.h"

class AliFlowTrackSimple;

// AliFlowLYZEventPlane:
//
// Class to calculate the event plane and event weight from the LYZ method
//
// author: N. van der Kolk (kolk@nikhef.nl)

ClassImp(AliFlowLYZEventPlane)

  //-----------------------------------------------------------------------
  
AliFlowLYZEventPlane::AliFlowLYZEventPlane():
  fFirstRunFile(0),
  fSecondRunFile(0),
  fFirstRunFileName("noname.ESD"),
  fSecondRunFileName("noname.ESD"),
  fWR(0),            
  fPsi(0),
  fSecondReDtheta(0),
  fSecondImDtheta(0),
  fFirstr0theta(0)
{
  // Constructor.
  fQ.Set(0.,0.);           // flow vector
}
//-----------------------------------------------------------------------


AliFlowLYZEventPlane::~AliFlowLYZEventPlane() 
{
  //destructor
   
}
 

//-----------------------------------------------------------------------
void AliFlowLYZEventPlane::Init() 
{
  //Declare histograms & get input files
  cout<<"---Lee Yang Zeros Event Plane Method---"<<endl;

  //input histograms
  if (fSecondRunFile->IsZombie()){ //check if file exists
    cout << "Error opening file, run first regular LYZ second run" << endl;
    exit(-1);
  } else if (fSecondRunFile->IsOpen()){
    cout<<"----secondRunFile is open----"<<endl;
    fSecondReDtheta = ( TProfile*)fSecondRunFile->Get("Second_FlowPro_ReDtheta_LYZ");
    fSecondImDtheta = ( TProfile*)fSecondRunFile->Get("Second_FlowPro_ImDtheta_LYZ");
  }
  if (fFirstRunFile->IsZombie()){ //check if file exists
    cout << "Error opening file, run first regular LYZ first run" << endl;
    exit(-1);
  } else if (fFirstRunFile->IsOpen()){
    cout<<"----firstRunFile is open----"<<endl<<endl;
    fFirstr0theta = (TProfile*)fFirstRunFile->Get("First_FlowPro_r0theta_LYZ");
  }



}


//-----------------------------------------------------------------------
void AliFlowLYZEventPlane::CalculateRPandW(TVector2 fQ)
{
  //declare variables
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Double_t cosTerm = 0;
  Double_t sinTerm = 0;
  TComplex fDtheta;
  TComplex ratio;

  for (Int_t theta=0;theta<fNtheta;theta++)	{
    Double_t fTheta = ((float)theta/fNtheta)*TMath::Pi()/2;  
    //Calculate Qtheta 
    Double_t fQtheta = fQ.X()*cos(2*fTheta)+fQ.Y()*sin(2*fTheta);  //get Qtheta from Q vector
            
    //get R0
    Double_t fR0 = fFirstr0theta->GetBinContent(theta+1); 
    //cerr<<"fR0 = "<<fR0<<endl;

    //get Dtheta
    Double_t fReDtheta = fSecondReDtheta->GetBinContent(theta+1);
    //cerr<<"fReDtheta = "<<fReDtheta<<endl;
    Double_t fImDtheta = fSecondImDtheta->GetBinContent(theta+1);
    //cerr<<"fImDtheta = "<<fImDtheta<<endl;
    fDtheta(fReDtheta,fImDtheta);
    
    TComplex fExpo(0.,fR0*fQtheta);                  //Complex number: 0 +(i r0 Qtheta)
    if (fDtheta.Rho()!=0.) { ratio =(TComplex::Exp(fExpo))/fDtheta; } //(e^(i r0 Qtheta))/Dtheta
    else { ratio(0.,0.); }
      
    //sum over theta
    cosTerm += ratio.Re() * TMath::Cos(2*fTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } cos(2 theta)  
    sinTerm += ratio.Re() * TMath::Sin(2*fTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } sin(2 theta)
      
  }//loop over theta

  //average over theta
  cosTerm /= fNtheta;  
  sinTerm /= fNtheta;
  //cerr<<"cosTerm & sinTerm are: "<<cosTerm<<" & "<<sinTerm<<endl;
    
  //calculate fWR
  fWR = TMath::Sqrt(cosTerm*cosTerm + sinTerm*sinTerm);
        
  //calculate fRP
  fPsi = 0.5*TMath::ATan2(sinTerm,cosTerm);   //takes care of the signs correctly!
  if (fPsi < 0.) { fPsi += TMath::Pi(); }     //to shift distribution from (-pi/2 to pi/2) to (0 to pi)
  



}

//-----------------------------------------------------------------------   
TVector2 AliFlowLYZEventPlane::GetQ(AliFlowEventSimple* fEvent) 
{
  //get the Q vector
  TVector2 fQ = fEvent->GetQ();
  
  return fQ;
  
}
