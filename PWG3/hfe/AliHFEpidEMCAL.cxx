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
//
// Class for EMCAL PID
// Implements the abstract base class AliHFEpidBase
// IsInitialized() does the PID decision with energy and 
// momentum matching (e/p) 
// 
// Authors:
//  Shingo Sakai 
//   
//
//
#include <TMath.h>
#include "AliESDInputHandler.h"
#include "AliESDpid.h"

#include "AliHFEdetPIDqa.h"
#include "AliHFEpidEMCAL.h"
#include "AliHFEpidQAmanager.h"

#include "AliHFEemcalPIDqa.h"
//#include "AliVCluster.h"
//#include "AliVCaloCells.h"
//#include "AliVEvent.h"
#include "AliLog.h"
#include "AliPID.h"
//#include "AliESDEvent.h"
//#include "AliESDtrack.h"
#include "AliHFEpidEMCAL.h"

ClassImp(AliHFEpidEMCAL)

//___________________________________________________________________
AliHFEpidEMCAL::AliHFEpidEMCAL():
  AliHFEpidBase()
  , fPID(NULL)
  , feopMim(0.8)
  , feopMax(1.4)
{
  //
  // Constructor
  //
} 

//___________________________________________________________________
AliHFEpidEMCAL::AliHFEpidEMCAL(const Char_t *name):
  AliHFEpidBase(name)
  , fPID(NULL)
  , feopMim(0.8)
  , feopMax(1.4)
{
  //
  // Constructor
  //
}

//___________________________________________________________________
AliHFEpidEMCAL::AliHFEpidEMCAL(const AliHFEpidEMCAL &c):
  AliHFEpidBase("")
  , fPID(NULL)
  , feopMim(0.8)
  , feopMax(1.4)
{  
  // 
  // Copy operator
  //

  c.Copy(*this);
}
//___________________________________________________________________
AliHFEpidEMCAL &AliHFEpidEMCAL::operator=(const AliHFEpidEMCAL &ref){
  //
  // Assignment operator
  //

  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}
//___________________________________________________________________
AliHFEpidEMCAL::~AliHFEpidEMCAL(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
}
//___________________________________________________________________
void AliHFEpidEMCAL::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidEMCAL &target = dynamic_cast<AliHFEpidEMCAL &>(ref);

  target.fPID = fPID;          
  target.feopMax = feopMax;
  target.feopMim = feopMim;

  AliHFEpidBase::Copy(ref);
}
//___________________________________________________________________
Bool_t AliHFEpidEMCAL::InitializePID(Int_t /*run*/){
  //
  // InitializePID: EMCAL experts have to implement code here
  //
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidEMCAL::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{ //Function to a return a code indicating whether or not an electron candidate is selected
  //
  //
  if(track==NULL)return 0;

  if(!fkPIDResponse) return 0;
  AliDebug(2, "PID object available");
  // EMCal not fESDpid  (s.s Feb. 11)
	
  const AliVTrack *trk = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  if (trk == NULL) return 0;
  //AliHFEpidObject::AnalysisType_t anaType = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  if(!(trk->GetStatus() & AliESDtrack::kEMCALmatch)) return 0;
  AliDebug(2, "Track Has EMCAL PID");
  
  //if(pidqa) 
  pidqa->ProcessTrack(track, AliHFEpid::kEMCALpid, AliHFEdetPIDqa::kBeforePID);
  // not QA for EMCal, will be added (s.s Feb. 11)


     double feopMimCut = 0.0; 
     double feopMaxCut = 0.0; 

     feopMimCut = feopMim;
     feopMaxCut = feopMax;

     if(feopMax >900)  // not use fix cut
       {
        feopMimCut = CalEopCutMim(track->GetRecTrack(),0);
        feopMaxCut = CalEopCutMax(track->GetRecTrack(),0); 
       }
     if(feopMax < -900)  // not use fix cut for MC
       {
        feopMimCut = CalEopCutMim(track->GetRecTrack(),1);
        feopMaxCut = CalEopCutMax(track->GetRecTrack(),1); 
       }

  //printf("eop cuts org; %g; %g \n",feopMim,feopMax);
  //printf("eop cuts ; %g; %g \n",feopMimCut,feopMaxCut);

  Double_t eop = MomentumEnergyMatchV2(track->GetRecTrack()); // get eop (What is GetRecTrack ?)
  AliDebug(2, Form("Energy - Momentum Matching e/p : %f", eop));
  AliDebug(2, Form("E/p cut ; %g ; %g \n", feopMim,feopMax));
  Int_t pdg = 0;
  //if(eop>feopMim && eop<feopMax){
  //if(eop>feopMim && eop<feopMax){
  if(eop>feopMimCut && eop<feopMaxCut){
    pdg = 11;
    //printf("eop cuts ; %g; %g ; %g\n",feopMimCut,feopMaxCut,eop);
    //if(pidqa) 
    pidqa->ProcessTrack(track, AliHFEpid::kEMCALpid, AliHFEdetPIDqa::kAfterPID);
  }
  else
  {
   pdg = 211; // EMCal doesn't separate pi,k.p by e/p. return pion code as same as TRD
  } 

    AliDebug(1, Form("eID %g ; %d \n",eop,pdg));  
    //printf("eID %g ; %d \n",eop,pdg);  

  return pdg;
}


//__________________________________________________________________________
Double_t AliHFEpidEMCAL::CalEopCutMax(const AliVParticle *const track, Int_t flageop) const
{
  double maxCut = 0.0;
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack==NULL)return maxCut;
  double Pt = esdtrack->Pt();

  if(flageop<0.5)
    {
     double meanP[3] = {0.991,1.0819,0.235}; 
     double sigP[3] = {6.52e-02,2.04e-01,3.34e-01}; 
     double mean = meanP[0]*tanh(meanP[1]+meanP[2]*Pt); 
     double sig = sigP[0]/tanh(sigP[1]+sigP[2]*Pt); 
     maxCut = mean+3.0*sig; 
    }
  else
    {
     double meanP[3] = {0.99,1.299,0.35}; 
     double sigP[3] = {4.11e-02,1.588e-01,2.664e-01}; 
     double mean = meanP[0]*tanh(meanP[1]+meanP[2]*Pt); 
     double sig = sigP[0]/tanh(sigP[1]+sigP[2]*Pt); 
     maxCut = mean+3.0*sig; 
    }
     return maxCut;
}

Double_t AliHFEpidEMCAL::CalEopCutMim(const AliVParticle *const track, Int_t flageop) const
{
  double mimCut = 0.0;
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack==NULL)return mimCut;
  double Pt = esdtrack->Pt();

  if(flageop<0.5) // real
     { 
      //printf("real"); 
     double meanP[3] = {0.991,1.0819,0.235}; 
     double sigP[3] = {6.52e-02,2.04e-01,3.34e-01}; 
     double mean = meanP[0]*tanh(meanP[1]+meanP[2]*Pt); 
     double sig = sigP[0]/tanh(sigP[1]+sigP[2]*Pt); 
     mimCut = mean-2.0*sig;
    }
   else // MC
    { 
      //printf("MC"); 
     double meanP[3] = {0.99,1.299,0.35}; 
     double sigP[3] = {4.11e-02,1.588e-01,2.664e-01}; 
     double mean = meanP[0]*tanh(meanP[1]+meanP[2]*Pt); 
     double sig = sigP[0]/tanh(sigP[1]+sigP[2]*Pt); 
     mimCut = mean-2.0*sig;
    }
     return mimCut;
}

//___________________________________________________________________________
Double_t AliHFEpidEMCAL::MomentumEnergyMatchV2(const AliVParticle *const track) const
{ // Returns e/p & Rmatch

  Double_t matchclsE = 9999.9;
  double feop = -9999.9;

  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack==NULL)return feop;
  const AliESDEvent *evt = esdtrack->GetESDEvent();

   Int_t icl = esdtrack->GetEMCALcluster();

   AliVCluster *cluster = (AliVCluster*) evt->GetCaloCluster(icl);
   if(!cluster->IsEMCAL()) {return feop;}

       matchclsE = cluster->E();
  

   if(matchclsE<9999.0) feop = matchclsE/esdtrack->P();

   return feop;

}


/*
//____________________________________________________________________________________
Double_t AliHFEpidEMCAL::MomentumEnergyMatchV1(const AliVParticle *track) const
{ // Returns e/p if an electron is matched

  Float_t  clsPos[3];
  Double_t trkPos[3];
  Double_t matchclsE = 9999.9;

  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  AliESDEvent *evt = esdtrack->GetESDEvent();
  Double_t  magF = evt->GetMagneticField();
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  //printf("magF ; %g ; %g \n", magF,magSign);
 
  if (!TGeoGlobalMagField::Instance()->GetField()) {
	  printf("Loading field map...\n");
	  //AliMagF* field = new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG);
	  AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG); // for 10d
	  TGeoGlobalMagField::Instance()->SetField(field);
	            }

  AliEMCALTrack *emctrack = new AliEMCALTrack(*esdtrack);
  Double_t fieldB[3]; 
  emctrack->GetBxByBz(fieldB);
  //printf("%g %g %g \n", fieldB[0], fieldB[1], fieldB[2]);

  for(Int_t icl=0; icl<evt->GetNumberOfCaloClusters(); icl++)
   {

   AliVCluster *cluster = (AliVCluster*) evt->GetCaloCluster(icl);
   if(!cluster->IsEMCAL()) continue;

   cluster->GetPosition(clsPos);
   if(!emctrack->PropagateToGlobal(clsPos[0],clsPos[1],clsPos[2],0.,0.) )  continue;
   emctrack->GetXYZ(trkPos);

   TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
   TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);

   Double_t delEmcphi = clsPosVec.Phi()-trkPosVec.Phi();  // track cluster matching
   Double_t delEmceta = clsPosVec.Eta()-trkPosVec.Eta();  // track cluster matching

   double rmatch = sqrt(pow(delEmcphi,2)+pow(delEmceta,2));

   if(rmatch<0.02)
      {
       matchclsE = cluster->E();
      }
  }
   delete emctrack;

   double feop = -9999.9;
   if(matchclsE<9999) feop = matchclsE/esdtrack->P();

   //   if(feop!=-9999.9)printf("%f\n",feop) ; 
   return feop;

}
*/

