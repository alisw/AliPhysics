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
#include "AliPIDResponse.h"

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
#include "AliTrackerBase.h"

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
  const Double_t kMaxP = 50.;
  const Double_t kTPCSigMim = 40.;
  const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, phi, eta,  Sig,  e/p,  ,match, cell, M02, M20, Disp, Centrality, select)
  Int_t nBins[16] = {AliPID::kSPECIES + 1, 500, 500,          400, 300,   100,   600,  300, 100,   40,   300, 300, 300, kCentralityBins, 6, 2};
  Double_t min[16] = {-1,               kMinP, kMinP,  kTPCSigMim,  0.,  -1.0,  -8.0,    0,   0,    0,   0.0, 0.0, 0.0,   0, -3.0, 0.};
  Double_t max[16] = {AliPID::kSPECIES, kMaxP, kMaxP,  kTPCSigMax, 6.0,   1.0,   4.0,  3.0, 0.1,   40,   3.0, 3.0, 3.0,   11., 3.0, 2.};
  fHistos->CreateTHnSparse("EMCAL_TPCdedx", "EMCAL signal; species; p [GeV/c]; pT [GeV/c] ; TPC signal [a.u.]; phi ; eta ; nSig ; E/p ; Rmatch ; Ncell ; M02 ; M20 ; Disp ; Centrality;charge ;PID Step; ", 16, nBins, min, max);

  //2nd histogram: EMCAL signal - E/p 
  //Int_t nBins2[7] = {AliPID::kSPECIES + 1, 500, 500, 500, 100, 600, 2};
  //Double_t min2[7] = {-1, kMinP, kMinP, 0,  0, -8., 0.};
  //Double_t max2[7] = {AliPID::kSPECIES, kMaxP, kMaxP, 5,  0.1, 4., 2.};
  //fHistos->CreateTHnSparse("EMCAL_Signal", "EMCAL true signal; species; p [GeV/c]; pT [GeV/c] ; E/p; Rmatch; TPCnsigma; PID Step", 7, nBins2, min2, max2);
    
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

  // Get TPC nSigma (based on AliHFEtpcPIDqa)
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid)); 
  const AliPIDResponse *pidResponse = tpcpid ? tpcpid->GetPIDResponse() : NULL;
  Float_t nSigmatpc = pidResponse ? pidResponse->NumberOfSigmasTPC(static_cast<const AliVTrack *>(track->GetRecTrack()), AliPID::kElectron) : 0.; 
  //printf("nSigmatpc = %f\n",nSigmatpc);
  AliDebug(2, Form("nSigmatpc = %f\n",nSigmatpc));
  //

  double charge = esdtrack->Charge();
  //printf("charge %f\n",charge); 


  //printf("phi %f, eta %f\n;",phi,eta);

  // Get E/p
  Double_t showerEMC[4];
  TVector3 emcsignal = MomentumEnergyMatchV2(esdtrack,showerEMC);
  AliDebug(2, Form("shower shape ; %f %f %f %f \n",showerEMC[0],showerEMC[1],showerEMC[2],showerEMC[3]));

  /*
  Double_t contentSignal2[7];
  contentSignal2[0] = species;
  contentSignal2[1] = track->GetRecTrack()->P();
  contentSignal2[2] = track->GetRecTrack()->Pt();
  contentSignal2[3] = emcsignal(0);//e over p
  contentSignal2[4] = emcsignal(1);//residual
  contentSignal2[5] = nSigmatpc;
  contentSignal2[6] = step == AliHFEdetPIDqa::kBeforePID ? 0. : 1.;
  */
  // QA array
  Double_t contentSignal[16];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = track->GetRecTrack()->Pt();
  contentSignal[3] = esdtrack->GetTPCsignal(); //?
  double phi  = track->GetRecTrack()->Phi();
  double eta  = track->GetRecTrack()->Eta();
  contentSignal[4] = phi;
  contentSignal[5] = eta;
  contentSignal[6] = nSigmatpc;
  contentSignal[7] = emcsignal(0);
  contentSignal[8] = emcsignal(2);
  contentSignal[9] = showerEMC[0];
  contentSignal[10] = showerEMC[1];
  contentSignal[11] = showerEMC[2];
  contentSignal[12] = showerEMC[3];
  contentSignal[13] = centrality;
  contentSignal[14] = charge;
  contentSignal[15] = step == AliHFEdetPIDqa::kBeforePID ? 0. : 1.;


  //printf("ProcessTrack ; Print Content %g; %g; %g; %g \n",contentSignal[0],contentSignal[1],contentSignal[2],contentSignal[3]); 
  fHistos->Fill("EMCAL_TPCdedx", contentSignal);
  //fHistos->Fill("EMCAL_Signal", contentSignal2);
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
TVector3 AliHFEemcalPIDqa::MomentumEnergyMatchV2(const AliESDtrack *esdtrack, Double_t *shower) const
{ // Returns e/p & Rmatch
  /*
  Float_t  clsPos[3];
  Double_t trkPos[3];
  Double_t fMass=0.139;
  Double_t fStep=10; //This is taken from EMCAL tender, hardcoded!
  */
  TVector3 refVec(-9999,-9999,-9999);
  Double_t matchclsE = 9999.9;
  for(int i=0; i<4; i++)shower[i] = -9999;

  const AliESDEvent *evt = esdtrack->GetESDEvent();

  //Proper trackmatching, added by Tomas
  //Get matched cluster
  Int_t icl = esdtrack->GetEMCALcluster(); //From tender
  AliVCluster *cluster = (AliVCluster*) evt->GetCaloCluster(icl);
  if(cluster==NULL) return refVec;
  if(!cluster->IsEMCAL()) return refVec;
  printf("EMCal cluster pointer ? %p\n",(void*)cluster);
  printf("EMCal cluster ? %d\n",cluster->IsEMCAL());

  // from ESDs
  double delphi = cluster->GetTrackDx(); 
  double deleta = cluster->GetTrackDz(); 
  double rmatch1 = sqrt(pow(delphi,2)+pow(deleta,2));

  matchclsE = cluster->E();
  
  //double feop = -9999.9;
  //if(matchclsE<9999) 
  double feop = matchclsE/esdtrack->P();
  
	shower[0] = cluster->GetNCells(); // number of cells in cluster
	shower[1] = cluster->GetM02(); // long axis
	shower[2] = cluster->GetM20(); // short axis
	shower[3] = cluster->GetDispersion(); // dispersion

  //   if(feop!=-9999.9)printf("%f\n",feop) ; 
  TVector3 emcsignal(feop,0.0,rmatch1);
  
  return emcsignal;
  
}
