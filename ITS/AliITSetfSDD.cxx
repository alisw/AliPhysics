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
Revision 1.1.2.2  2000/06/12 18:09:36  barbera
fixed posible compilation errors on HP unix

Revision 1.1.2.1  2000/01/12 20:21:30  nilsen
missed AliITSetfSDD files ealier

$Name$
$Author$
$Id$
*/

#include <iostream.h>
#include <TMath.h>
#include "AliITSetfSDD.h"

////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Piergiorgio Cerello
// November 23 1999
//
// AliITSmapSDD is the map of SDDs.
//
//Begin_Html
/*
<img src="picts/ITS/AliITShit_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the relasionships between the ITS hit class and the rest of Aliroot.
</font>
<pre>
*/
//End_Htm
//extern Int_t power(Int_t b, Int_t e);

//_____________________________________________________________________________
ClassImp(AliITSetfSDD)

Int_t ppower(Int_t b, Int_t e) {
  Int_t power = 1,i;
  for(i=0; i<e; i++) power *= b;
  return power;
}

AliITSetfSDD::AliITSetfSDD(Double_t timeclock)
{
  // sampling time in ns
  fSamplingTime = 1000./timeclock;
  fT0 = 0.;
  fDf = ppower(10,9)/(fMaxNofSamples*fSamplingTime);
  fA0 = 9000.;

  Int_t i,j,k;
  for(i=0; i<fMaxNofPoles; i++) {
    fZero_M[i] = 0.;
    fZero_R[i] = 0.;
    fZero_I[i] = 0.;
    fPole_M[i] = 0.;
    fPole_R[i] = 0.;
    fPole_I[i] = 0.;
  }
  fPole_M[0] = 1.;
  fPole_R[0] = -2100000.;
  fPole_I[0] = fPole_R[0];
  fPole_M[1] = 1.;
  fPole_R[1] = -2100000.;
  fPole_I[1] = -fPole_R[1];

   // Compute Transfer Function

  Double_t PI = acos(-1.);
  for(i=0; i<=fMaxNofSamples/2; i++) {
    Double_t frequency = fDf*i;
    Double_t VM = fA0;
    Double_t VA = 0.;
    for(k=0; k<fMaxNofPoles; k++) {
      if(fZero_M[k]) {
        Double_t VZR = -fZero_R[k];
        Double_t VZI = frequency - fZero_I[k];
        Double_t VZM = TMath::Sqrt(VZR*VZR+VZI*VZI);
        Double_t VZA = TMath::ATan2(VZI,VZR);
	//	cout << "VZM: " << VZM << ", VZA: " << VZA << endl;
	//	cout << "VZR: " << VZR << ", VZI: " << VZI << endl;
        for(j=1; j<= (Int_t) fZero_M[k]; j++) {
          VM *= VZM;
          VA += VZA;
          if(VA >= PI) VA -= (2.*PI);
          if(VA <= -PI) VA += (2.*PI);
	  //cout << "VM: " << VM << ", VA: " << VA << endl;
        }
      }

      if(fPole_M[k]) {
        Double_t VPR = -fPole_R[k];
        Double_t VPI = frequency - fPole_I[k];
	Double_t VPM = TMath::Sqrt(VPR*VPR+VPI*VPI);
	Double_t VPA = TMath::ATan2(VPI,VPR);
	//cout << "VPM: " << VPM << ", VPA: " << VPA << endl;
	//cout << "VPR: " << VPR << ", VPI: " << VPI << endl;
        for(j=1; j<= (Int_t) fPole_M[k]; j++) {
          VM /= VPM;
          VA -= VPA;
          if(VA >= PI) VA -= (2.*PI);
          if(VA <= -PI) VA += (2.*PI);
	  //cout << "VM: " << VM << ", VA: " << VA << endl;
        }
      }
      Double_t VR = VM*cos(VA);
      Double_t VI = VM*sin(VA);
      //cout << "VM: " << VM << ", VA: " << VA << endl;
      //cout << "VR: " << VR << ", VI: " << VI << endl;
      fTf_R[i] = VR*ppower(10,9);
      fTf_I[i] = VI*ppower(10,9);
      //cout << "fTf_R[" << i << "] = " << fTf_R[i] << endl;
      //cout << "fTf_I[" << i << "] = " << fTf_I[i] << endl;
      if(i) {
        fTf_R[fMaxNofSamples-i] = fTf_R[i];
        fTf_I[fMaxNofSamples-i] = -fTf_I[i];
      }
    }
  }
  
  // Compute Fourier Weights

  for(i=0; i<=fMaxNofSamples/2; i++) {
    fW_R[i] = cos(-2.*PI*i/fMaxNofSamples);
    fW_I[i] = sin(-2.*PI*i/fMaxNofSamples);
    if(i) {
      fW_R[fMaxNofSamples-i] = fW_R[i];
      fW_I[fMaxNofSamples-i] = -fW_I[i];
    }
  }

}

void AliITSetfSDD::Print()
{
  Int_t i;
  cout << "Sampling Time " << fSamplingTime << endl;
  cout << "Number of Time Samples " << fMaxNofSamples << endl;
  cout << "fT0 " << fT0 << endl;
  cout << "fDf " << fDf << endl;
  cout << "fA0 " << fA0 << endl;

  cout << "Zero's and Pole's" << endl;
  cout << "fZero_M " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fZero_M[i] << endl;
  cout << "fZero_R " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fZero_R[i] << endl;
  cout << "fZero_I " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fZero_I[i] << endl;
  cout << "fPole_M " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fPole_M[i] << endl;
  cout << "fPole_R " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fPole_R[i] << endl;
  cout << "fPole_I " << endl;
  for(i=0; i<fMaxNofPoles; i++) cout << fPole_I[i] << endl;

  cout << "Transfer function" << endl;
  cout << "Real Part" << endl;
  for(i=0; i<fMaxNofSamples; i++) cout << fTf_R[i] << endl;
  cout << "Imaginary Part " << endl;
  for(i=0; i<fMaxNofSamples; i++) cout << fTf_I[i] << endl;

  cout << "Fourier Weights" << endl;
  cout << "Real Part" << endl;
  for(i=0; i<fMaxNofSamples; i++) cout << fW_R[i] << endl;
  cout << "Imaginary Part " << endl;
  for(i=0; i<fMaxNofSamples; i++) cout << fW_I[i] << endl;
}








