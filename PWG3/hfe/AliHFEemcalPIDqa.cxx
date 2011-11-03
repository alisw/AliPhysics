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
// Class AliHFEemcalPIDqa
//
// Monitoring EMCAL PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   TPC Signal distribution 
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//
//   Shingo Sakai <ssakai@lbl.gov>
#include <TClass.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include <TMath.h>
#include "AliESDInputHandler.h"
//#include "AliVCluster.h"
//#include "AliVCaloCells.h"
//#include "AliVEvent.h"
#include "AliMagF.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"
//#include "AliVTrack.h"
//#include "AliESDtrack.h"
#include "AliHFEcollection.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidEMCAL.h"
#include "AliHFEtools.h"
#include "AliHFEemcalPIDqa.h"

ClassImp(AliHFEemcalPIDqa)

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for EMCAL")
  , fHistos(NULL)
{
  //
  // Default constructor
  //
}

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa(const AliHFEemcalPIDqa &o):
    AliHFEdetPIDqa(o)
  , fHistos()
{
  //
  // Copy constructor
  //
  o.Copy(*this);
}

//_________________________________________________________
AliHFEemcalPIDqa &AliHFEemcalPIDqa::operator=(const AliHFEemcalPIDqa &o){
  //
  // Do assignment
  //
  AliHFEdetPIDqa::operator=(o);
  if(&o != this) o.Copy(*this);
  return *this;
}

//_________________________________________________________
AliHFEemcalPIDqa::~AliHFEemcalPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}

//_________________________________________________________
void AliHFEemcalPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEemcalPIDqa &target = dynamic_cast<AliHFEemcalPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);
}

//_________________________________________________________
Long64_t AliHFEemcalPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEemcalPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEemcalPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEemcalPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("emcalqahistos", "Collection of EMCAL QA histograms");

  // Make common binning
  const Int_t kCentralityBins = 11;
  const Double_t kMinP = 0.;
  const Double_t kMaxP = 20.;
  const Double_t kTPCSigMim = 40.;
  const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, Centrality)
  Int_t nBins[6] = {AliPID::kSPECIES + 1, 20, 20, 400, kCentralityBins, 2};
  Double_t min[6] = {-1, kMinP, kMinP, kTPCSigMim,  0, 0.};
  Double_t max[6] = {AliPID::kSPECIES, kMaxP, kMaxP, kTPCSigMax,  11., 2.};
  fHistos->CreateTHnSparse("EMCAL_TPCdedx", "EMCAL signal; species; p [GeV/c]; pT [GeV/c] ; TPC signal [a.u.]; Centrality; PID Step", 6, nBins, min, max);

  //2nd histogram: EMCAL signal - E/p and Rmatch
  Int_t nBins2[6] = {AliPID::kSPECIES + 1, 40, 40, 1000, 250, 2};
  Double_t min2[6] = {-1, kMinP, kMinP, 0,  0, 0.};
  Double_t max2[6] = {AliPID::kSPECIES, kMaxP, kMaxP, 10,  1., 2.};
  fHistos->CreateTHnSparse("EMCAL_Signal", "EMCAL true signal; species; p [GeV/c]; pT [GeV/c] ; E/p; Rmatch; PID Step", 6, nBins2, min2, max2);
    
}




//_________________________________________________________
void AliHFEemcalPIDqa::ProcessTrack(const AliHFEpidObject *track,AliHFEdetPIDqa::EStep_t step){
  //
  // Fill TPC histograms
  //
  //AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t centrality = track->GetCentrality();

  //const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  //const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(vtrack);
  const AliESDtrack *esdtrack = static_cast<const AliESDtrack *>(track->GetRecTrack());

  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;

  Double_t contentSignal[6];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = track->GetRecTrack()->Pt();
  contentSignal[3] = esdtrack->GetTPCsignal(); //?
  contentSignal[4] = centrality;
  contentSignal[5] = step == AliHFEdetPIDqa::kBeforePID ? 0. : 1.;

  TVector3 emcsignal = MomentumEnergyMatchV2(esdtrack);


  Double_t contentSignal2[6];
  contentSignal2[0] = species;
  contentSignal2[1] = track->GetRecTrack()->P();
  contentSignal2[2] = track->GetRecTrack()->Pt();
  contentSignal2[3] = emcsignal(0);
  contentSignal2[4] = emcsignal(1);
  contentSignal2[5] = contentSignal[5];

  //printf("ProcessTrack ; Print Content %g; %g; %g; %g \n",contentSignal[0],contentSignal[1],contentSignal[2],contentSignal[3]); 
  fHistos->Fill("EMCAL_TPCdedx", contentSignal);
  fHistos->Fill("EMCAL_Signal", contentSignal2);
}

//_________________________________________________________
TH1 *AliHFEemcalPIDqa::GetHistogram(const char *name){
  // 
  // Get the histogram
  //
  if(!fHistos) return NULL;
  TObject *histo = fHistos->Get(name);
  if(!histo->InheritsFrom("TH1")) return NULL;
  return dynamic_cast<TH1 *>(histo);
}


//___________________________________________________________________________
TVector3 AliHFEemcalPIDqa::MomentumEnergyMatchV2(const AliESDtrack *esdtrack) const
{ // Returns e/p & Rmatch

  Float_t  clsPos[3];
  Double_t trkPos[3];
  Double_t matchclsE = 9999.9;
  TVector3 refVec(-9999,-9999,-9999);

  const AliESDEvent *evt = esdtrack->GetESDEvent();

   //Int_t icl = esdtrack->GetEMCALcluster();
   Int_t icl = (const_cast<AliESDtrack *>(esdtrack))->GetEMCALcluster();

   AliVCluster *cluster = (AliVCluster*) evt->GetCaloCluster(icl);
   if(!cluster->IsEMCAL()) return refVec;

   cluster->GetPosition(clsPos);
   esdtrack->GetXYZ(trkPos);

   TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
   TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);


   Double_t delEmcphi = clsPosVec.Phi()-trkPosVec.Phi();  // track cluster matching
   Double_t delEmceta = clsPosVec.Eta()-trkPosVec.Eta();  // track cluster matching

   double rmatch = sqrt(pow(delEmcphi,2)+pow(delEmceta,2));

   matchclsE = cluster->E();

   //double feop = -9999.9;
   //if(matchclsE<9999) 
   double feop = matchclsE/esdtrack->P();

   //   if(feop!=-9999.9)printf("%f\n",feop) ; 
   TVector3 emcsignal(feop,rmatch,0);

   return emcsignal;

}




/*
//___________________________________________________________________________
TVector3 AliHFEemcalPIDqa::MomentumEnergyMatchV1(const AliESDtrack *esdtrack) const
{ // Returns e/p & Rmatch

  Float_t  clsPos[3];
  Double_t trkPos[3];
  Double_t Rmatch;
  Double_t matchclsE = 9999.9;
  TVector3 refVec(-9999,-9999,-9999);

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
   if(!cluster->IsEMCAL()) return refVec;

   cluster->GetPosition(clsPos);
   if(!emctrack->PropagateToGlobal(clsPos[0],clsPos[1],clsPos[2],0.,0.) )  return refVec;
   emctrack->GetXYZ(trkPos);

   TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
   TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);

   Double_t delEmcphi = clsPosVec.Phi()-trkPosVec.Phi();  // track cluster matching
   Double_t delEmceta = clsPosVec.Eta()-trkPosVec.Eta();  // track cluster matching

   double rmatch = sqrt(pow(delEmcphi,2)+pow(delEmceta,2));

   if(rmatch<0.02)
      {
       matchclsE = cluster->E();
       Rmatch = rmatch;
      }
  }
   delete emctrack;

   //double feop = -9999.9;
   //if(matchclsE<9999) 
   double feop = matchclsE/esdtrack->P();

   //   if(feop!=-9999.9)printf("%f\n",feop) ; 
   TVector3 emcsignal(feop, Rmatch, 0);
 
   return emcsignal;

}
*/

