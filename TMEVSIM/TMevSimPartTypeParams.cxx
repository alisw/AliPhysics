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
*/

#include "TMevSimPartTypeParams.h"

ClassImp(TMevSimPartTypeParams) 

  
//______________________________________________________________________________
  TMevSimPartTypeParams::TMevSimPartTypeParams() 
{ 
// Default constructor
  
  for (Int_t i = 0; i < 4; i++)
    for (Int_t j = 0; j < NFLOWTERMS; j++) {
      fVnMean[j][i] = 0.0;
      fVnStDev[j][i] = 0.0;
    }
}

//__________________________________________________________________________
TMevSimPartTypeParams::TMevSimPartTypeParams(Int_t agpid, Int_t amultmean, Int_t amultvc, 
					     Float_t atempmean, Float_t atempstdev, Float_t asigmamean,
					     Float_t asigmastdev, Float_t aexpvelmean, Float_t aexpvelstdev) 
{
// Construct the particle type parametrs class. Use the values provide
// by the user to initialize the internal variables. For the meaning of
// the parametrs see the TMevSim class documentation.
   
   fGPid = agpid;
   fMultMean = amultmean;
   fMultVarianceControl = amultvc;
   fTempMean = atempmean;
   fTempStDev = atempstdev;
   fSigmaMean = asigmamean;
   fSigmaStDev = asigmastdev;
   fExpVelMean = aexpvelmean;
   fExpVelStDev = aexpvelstdev;
   for (Int_t i = 0; i < 4; i++)
     for (Int_t j = 0; j < NFLOWTERMS; j++) {
	fVnMean[j][i] = 0.0;
	fVnStDev[j][i] = 0.0;
     }
   
}

//______________________________________________________________________________
TMevSimPartTypeParams::~TMevSimPartTypeParams() {}

//______________________________________________________________________________
TMevSimPartTypeParams::TMevSimPartTypeParams (const TMevSimPartTypeParams& pars) {
// The copy constructor
  
   this->fGPid = pars.GetGPid();
   this->fMultMean = pars.GetMultMean();
   this->fMultVarianceControl = pars.GetMultVarianceControl();
   this->fTempMean = pars.GetTempMean();
   this->fTempStDev = pars.GetTempStDev();
   this->fSigmaMean = pars.GetSigmaMean();
   this->fSigmaStDev = pars.GetSigmaStDev();
   this->fExpVelMean = pars.GetExpVelMean();
   this->fExpVelStDev = pars.GetExpVelStDev();
   for (Int_t i = 0; i < 4; i++)
     for (Int_t j = 0; j < NFLOWTERMS; j++) {
       this->fVnMean[j][i] = pars.GetVnMeanComponent(j, i);
       this->fVnStDev[j][i] = pars.GetVnStDevComponent(j, i);
     }
}

//______________________________________________________________________________
TMevSimPartTypeParams& TMevSimPartTypeParams::operator=(const TMevSimPartTypeParams& pars) {
// The assignment operator
   
   this->fGPid = pars.GetGPid();
   this->fMultMean = pars.GetMultMean();
   this->fMultVarianceControl = pars.GetMultVarianceControl();
   this->fTempMean = pars.GetTempMean();
   this->fTempStDev = pars.GetTempStDev();
   this->fSigmaMean = pars.GetSigmaMean();
   this->fSigmaStDev = pars.GetSigmaStDev();
   this->fExpVelMean = pars.GetExpVelMean();
   this->fExpVelStDev = pars.GetExpVelStDev();
   for (Int_t i = 0; i < 4; i++)
     for (Int_t j = 0; j < NFLOWTERMS; j++) {
       this->fVnMean[j][i] = GetVnMeanComponent(j, i);
       this->fVnStDev[j][i] = GetVnStDevComponent(j, i);
     }
   return (*this);
}

//______________________________________________________________________________
void    TMevSimPartTypeParams::SetGPid(Int_t gpid) {
   fGPid = gpid;
}
//______________________________________________________________________________
Int_t   TMevSimPartTypeParams::GetGPid() const {
   return fGPid;
}
//______________________________________________________________________________
void    TMevSimPartTypeParams::SetMultMean(Int_t multmean) {
   fMultMean = multmean;
}
//______________________________________________________________________________
Int_t TMevSimPartTypeParams::GetMultMean() const {
   return fMultMean;
}
//______________________________________________________________________________
void    TMevSimPartTypeParams::SetMultVarianceControl(Int_t multvc) {
   fMultVarianceControl = multvc;
}
//______________________________________________________________________________
Int_t   TMevSimPartTypeParams::GetMultVarianceControl() const {
   return fMultVarianceControl;
}
//______________________________________________________________________________
void    TMevSimPartTypeParams::SetTempParams(Float_t tempmean, Float_t tempsdev) { 
   fTempMean = tempmean;
   fTempStDev = tempsdev;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetTempMean() const { 
   return fTempMean;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetTempStDev() const { 
   return fTempStDev;
}
//______________________________________________________________________________
void    TMevSimPartTypeParams::SetSigmaPrams(Float_t sigmamean, Float_t sigmastdev) {  
   fSigmaMean = sigmamean;
   fSigmaStDev = sigmastdev;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetSigmaMean() const { 
   return fSigmaMean;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetSigmaStDev() const { 
   return fSigmaStDev;
}
//______________________________________________________________________________
void    TMevSimPartTypeParams::SetExpVelParams(Float_t expvelmean, Float_t expvelstdev) {  
   fExpVelMean = expvelmean;
   fExpVelStDev = expvelstdev;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetExpVelMean() const { 
   return fExpVelMean;
}
//______________________________________________________________________________
Float_t TMevSimPartTypeParams::GetExpVelStDev() const { 
   return fExpVelStDev;
}
//______________________________________________________________________________

void TMevSimPartTypeParams::SetVnMeanComponent
(Int_t nComponent, Float_t mean1, Float_t mean2,  Float_t mean3, Float_t mean4)  {  


  if (nComponent < 0 || nComponent > NFLOWTERMS )
    Error("SetVnMeanComponent", "Wrong Vn component n = %d (must be [%d , %d])", 
	  nComponent, 0, NFLOWTERMS);

  fVnMean[nComponent][0] = mean1;
  fVnMean[nComponent][1] = mean2;
  fVnMean[nComponent][2] = mean3;
  fVnMean[nComponent][3] = mean4;
}
//______________________________________________________________________________

void TMevSimPartTypeParams::SetVnStDevComponent
(Int_t nComponent, Float_t stdev1, Float_t stdev2,Float_t stdev3, Float_t stdev4) {  

  if (nComponent < 0 || nComponent > NFLOWTERMS )
    Error("SetVnStDevComponent", "Wrong Vn component n = %d (must be [%d , %d])",
	  nComponent, 0, NFLOWTERMS);

  fVnStDev[nComponent][0] = stdev1;
  fVnStDev[nComponent][1] = stdev2;
  fVnStDev[nComponent][2] = stdev3;
  fVnStDev[nComponent][3] = stdev4;
}

//______________________________________________________________________________

void TMevSimPartTypeParams::SetVnMeanComponent (Int_t nComponent, Float_t mean[4]) {
  
  if (nComponent < 0 || nComponent > NFLOWTERMS )
    Error("SetVnMeanComponent", "Wrong Vn component (%d) (must be [%d,%d])",
	nComponent, 0, NFLOWTERMS);
  
  Int_t N = 4;
  for (Int_t i=0; i<N; i++) 
    fVnMean[nComponent][i] = mean[i];

}

//______________________________________________________________________________

void TMevSimPartTypeParams::SetVnStDevComponent (Int_t nComponent, Float_t stdev[4]) {

  if (nComponent < 0 || nComponent > NFLOWTERMS )
    Error("SetVnStDevComponent", "Wrong Vn component n = %d (must be [%d , %d])",
	  nComponent, 0, NFLOWTERMS);
  
  Int_t N = 4;
  for (Int_t i=0; i<N; i++) 
    fVnStDev[nComponent][i] = stdev[i];  
}

//______________________________________________________________________________

Float_t     TMevSimPartTypeParams::GetVnMeanComponent(Int_t nComponent, Int_t nMean) const {

  if ((nComponent < 0) || (nComponent>NFLOWTERMS) || (nMean < 0) || (nMean > 3)) return 0.0;
  return fVnMean[nComponent][nMean];
}
//______________________________________________________________________________
Float_t     TMevSimPartTypeParams::GetVnStDevComponent(Int_t nComponent, Int_t nStDev) const {

  if ((nComponent < 0) || (nComponent>NFLOWTERMS) || (nStDev < 0) || (nStDev > 3)) return 0.0;
  return fVnMean[nComponent][nStDev];
}

//______________________________________________________________________________





