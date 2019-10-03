//Analysis Task to find the thermal photon distribution.  Made By Tyler Browning, Aug 26, 2014.

//#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
//#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"

#include <iostream>
#include <vector>

#include "TNtupleD.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "TList.h"
#include "TString.h"
#include "AliLog.h"

#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"

#include "AliAnalysisTaskThermalGAFlow.h"
#include "AliAnalysisTaskThermalGAFlowMC.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

using namespace std;

ClassImp(AliAnalysisTaskThermalGAFlowMC);

//_______________________________Constructor_______________________________//

AliAnalysisTaskThermalGAFlowMC::AliAnalysisTaskThermalGAFlowMC( const char *name) :
 AliAnalysisTaskThermalGAFlow(name),
 fStack(0x0)
{
}

//_______________________________Deconstructor_____________________________//
AliAnalysisTaskThermalGAFlowMC::~AliAnalysisTaskThermalGAFlowMC(){
if (fAnalist) {delete fAnalist;  fAnalist = 0;}
}

//_______________________________Create Output (Histogram TList)___________//

void AliAnalysisTaskThermalGAFlowMC::UserCreateOutputObjects() {
AliAnalysisTaskThermalGAFlow::UserCreateOutputObjects(); //Do setup and add normal Data hists
//now add MC hists

  Int_t Nbins = fNptbins; //("Clus_Pt", "Single Photon Pt Spectrum", 150, 0.5, 20.5); and exp scaling
  const Int_t somebins = Nbins+1;
  Double_t Ptbins[somebins];
  Double_t a = TMath::Power(10, 1);
  for(Int_t ibin = 0; ibin < somebins; ibin++){
  Ptbins[ibin] = 0.5 - (20/(a - 1)) + ((20/(a - 1)) * TMath::Power(10, ibin/(Nbins*1.0)));
  }

  TH1F *Clus_PrimeERatio = new TH1F("Clus_PrimeERatio", "Ratio of Highest Energy Prime to total cluster Energy", 50, 0, 1);
  Clus_PrimeERatio->GetXaxis()->SetTitle("PrimeERatio");
  Clus_PrimeERatio->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_PrimeERatio);

  TH1F *Clus_PtGamma = new TH1F("Clus_PtGamma", "Pt of the MC Gammas", Nbins, Ptbins);
  Clus_PtGamma->GetXaxis()->SetTitle("Pt");
  Clus_PtGamma->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_PtGamma);

  TH2F *Pion_M = new TH2F("Pion_M", "Pion Mass", 400,0,1, Nbins, Ptbins);
  Pion_M->GetXaxis()->SetTitle("M (GeV/c^2)");
  Pion_M->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(Pion_M);

  TH2F *Pion_MCalc = new TH2F("Pion_MCalc", "Pion Mass", 400,0,1, Nbins, Ptbins);
  Pion_MCalc->GetXaxis()->SetTitle("M (GeV/c^2)");
  Pion_MCalc->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(Pion_MCalc);

  TH1F *Pion_Pz = new TH1F("Pion_Pz", "Pz of the MC Gammas", Nbins, Ptbins);
  Pion_Pz->GetXaxis()->SetTitle("Pz");
  Pion_Pz->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Pion_Pz);

  TH1F *Pion_Pt = new TH1F("Pion_Pt", "Pt of the MC Gammas", Nbins, Ptbins);
  Pion_Pt->GetXaxis()->SetTitle("Pt");
  Pion_Pt->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Pion_Pt);

  TH1F *Pion_E = new TH1F("Pion_E", "E of the MC Gammas", Nbins, Ptbins);
  Pion_E->GetXaxis()->SetTitle("E");
  Pion_E->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Pion_E);

  TH2F *PionDaughter_M = new TH2F("PionDaughter_M", "Pion Mass", 400,0,1, Nbins, Ptbins);
  PionDaughter_M->GetXaxis()->SetTitle("M (GeV/c^2)");
  PionDaughter_M->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(PionDaughter_M);

  TH1I *Stat_Purity = new TH1I("Stat_Purity", "Purity of each cut number", Nbins, Ptbins);
  Stat_Purity->GetXaxis()->SetTitle("Counts");
  Stat_Purity->GetYaxis()->SetTitle("Cut Number");
  fAnalist->Add(Stat_Purity);

  TH1F *Stat_Flavor_dir = new TH1F("Stat_Flavor_dir", "Flavor of photon spectrum-dir", Nbins, Ptbins);
  Stat_Flavor_dir->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_dir->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_dir);

  TH1F *Stat_Flavor_pi = new TH1F("Stat_Flavor_pi", "Flavor of photon spectrum-pi", Nbins, Ptbins);
  Stat_Flavor_pi->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_pi->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_pi);

  TH1F *Stat_Flavor_eta = new TH1F("Stat_Flavor_eta", "Flavor of photon spectrum-eta", Nbins, Ptbins);
  Stat_Flavor_eta->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_eta->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_eta);

  TH1F *Stat_Flavor_f0 = new TH1F("Stat_Flavor_f0", "Flavor of photon spectrum-f0", Nbins, Ptbins);
  Stat_Flavor_f0->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_f0->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_f0);

  TH1F *Stat_Flavor_omega = new TH1F("Stat_Flavor_omega", "Flavor of photon spectrum-omega", Nbins, Ptbins);
  Stat_Flavor_omega->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_omega->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_omega);

  TH1F *Stat_Flavor_kaon = new TH1F("Stat_Flavor_kaon", "Flavor of photon spectrum-kaon", Nbins, Ptbins);
  Stat_Flavor_kaon->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_kaon->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_kaon);

  TH1F *Stat_Flavor_total = new TH1F("Stat_Flavor_total", "Flavor of photon spectrum-total", Nbins, Ptbins);
  Stat_Flavor_total->GetXaxis()->SetTitle("Pt");
  Stat_Flavor_total->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Stat_Flavor_total);

PostData(1,fAnalist);
}

//_______________________________Main Loop_________________________________//
void AliAnalysisTaskThermalGAFlowMC::UserExec( Option_t * option){
//printf("Anything.\n");
FillHist("Stat_Efficiency", 15);
PostData(1,fAnalist);
AliAnalysisTaskThermalGAFlow::CaptainsLog(0);

//Step 0: Configure the environment - take out your toys
AliAnalysisTaskThermalGAFlow::ConfigureEvent();
FillHist("Stat_Efficiency", 16);
PostData(1,fAnalist);
if(fRunNumber != fEvent->GetRunNumber()){fRunNumber = fEvent->GetRunNumber(); InitializeGeometry();}
FillHist("Stat_Efficiency", 17);
PostData(1,fAnalist);
AliAnalysisTaskThermalGAFlow::ScanBasicEventParameters();
AliAnalysisTaskThermalGAFlow::CaptainsLog(1);

//Step 1: Internal Triggering and QA - set up your toys
if(!AliAnalysisTaskThermalGAFlow::ApplyEventCuts()){PostData(1,fAnalist);}
else {
FillHist("Stat_Efficiency", 18);
AliAnalysisTaskThermalGAFlow::ScanClusters();
AliAnalysisTaskThermalGAFlow::SaveEvent();
AliAnalysisTaskThermalGAFlow::CaptainsLog(2);

//Step 2: Extract the diphoton mass - play with your toys
//printf("fSextant = %p\n", fSextant);
AliAnalysisTaskThermalGAFlow::MesonExclusion();
AliAnalysisTaskThermalGAFlow::CaptainsLog(3);
//Step 3: Compute thermal photon spectrum - have fun and be happy
AliAnalysisTaskThermalGAFlow::Hypnotic();

//Step 4: Compute thermal photon flow - put away your toys when you're done with them 
//Invalid - In MC, there is no flow.

//delete fSextant;
fSextant = 0x0;

//Step 5: Extract extra information since we're in MC - use your imagination!
fStack = GetMCStack();
ScanClustersMC();  //Then extract the needed MC information.
ExtractRealPions(); //Extract the real pion distribution, no cuts.
TasteFlavor();

PostData(1,fAnalist);
}


//if(!ApplyEventCuts()){PostData(1,fAnalist);}
//else{
//fStack = GetMCStack();
//AliAnalysisTaskThermalGAFlow::UserExec(option); //Do the Same thing as Data
//ScanClustersMC();  //Then extract the needed MC information.
//ExtractRealPions(); //Extract the real pion distribution, no cuts.

//PostData(1,fAnalist);
//}

}

//_________________Final step,  Called once at the end of the Task_________//

void AliAnalysisTaskThermalGAFlowMC::Terminate( Option_t * option){
AliAnalysisTaskThermalGAFlow::Terminate(option);
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                                                          //
//                               Functions                                  //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


//___________________Retrieve the Stack from the Kinematics Tree___________//
AliStack* AliAnalysisTaskThermalGAFlowMC::GetMCStack()
{
fStack = 0x0;
AliVEventHandler* vHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
AliMCEventHandler* mcHandler = 0x0;
if(vHandler) {mcHandler = dynamic_cast<AliMCEventHandler*> (vHandler);}
//printf("%p -- %p -- %p\n", mcHandler, static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent(), mcHandler->MCEvent());
if(mcHandler) {fStack = mcHandler->MCEvent()->Stack();}
if(!fStack){AliError("I couldn't find the Monte Carlo Stack!");}
return fStack;
}

//___________________Extract Useful information from the MC Clusters________//
void AliAnalysisTaskThermalGAFlowMC::ScanClustersMC()
{
TLorentzVector p;

if(fesd){
for(Int_t icluster = 0; icluster < fesd->GetNumberOfCaloClusters(); icluster++) {
  AliESDCaloCluster *cluster = fesd->GetCaloCluster(icluster);
//  GuessPure(cluster); //Note:  Might mess up some other plots.  Try to only guess pure when you need to.  Post analysis macro should give full purity anyways.
  if(ApplyClusterCuts(cluster)){
  if(IsPhotonShower(cluster)){
    cluster->GetMomentum(p, fvertexx);
    FillHist("Clus_PtGamma", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
  }}
}}

}

//___________________Try to Guess the purity_______________________________//
void AliAnalysisTaskThermalGAFlowMC::GuessPure(AliESDCaloCluster *cluster){
//if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 0);} //This can't work - you haven't done basic cuts yet, so labels are meaningless.

//Note:  Might mess up some other plots.  Try to only guess pure when you need to.  Post analysis macro should give full purity anyways.
if(WeakCuts(cluster)){
 if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 2);}

 if(cluster->GetTOF() > -0.0000001 && cluster->GetTOF() < 0.0000001){
  if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 3);}

  if(cluster->GetDistanceToBadChannel() > 2.2){
   if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 4);}

   if(ApplyCPCuts(cluster)){
    if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 5);}

    if(ApplyCoreShapeCut(cluster)){
     if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 6);}

     if(ApplyDispMatrixCut(cluster)){
      if(IsPhotonShower(cluster)){FillHist("Stat_Purity", 7);}   
}}}}}}
}

//___________________Check to see if a cluster was from a MC photon________//
Bool_t AliAnalysisTaskThermalGAFlowMC::IsPhotonShower(AliESDCaloCluster *cluster)
{
//Note:  I'm not really sure, but I think primaries are MC particles that contributed energy to the shower.  Labels, then, is the primaries sorted by deposited energy - Label[0] would be the highest contribution.
TParticle *p = fStack->Particle(cluster->GetLabelAt(0));
Int_t pCode = p->GetPdgCode(); //This is the MC PID code.
Float_t primaryERatio = p->Energy() / cluster->E(); //Prime Energy Ratio.  PHOS_PbPb used 0.9.  I'm going to use something smaller to start.
FillHist("Clus_PrimeERatio", primaryERatio);
if(pCode == 22 && primaryERatio > 0.80){return 1;} 
//22 is the photon code... I think.  11 is electron, -11 is position... again... I think. (This is verified - See "Monte Carlo Particle Numbering Scheme" June 2006, Garren et al.)
//pi0 = 111, eta = 221, f0 = 10221
return 0;
}

//__________________Extract the real Pion Spectrum________________________//
void AliAnalysisTaskThermalGAFlowMC::ExtractRealPions()
{
const TObjArray *pArray = fStack->Particles();
Int_t Npart = pArray->GetEntriesFast();
TParticle *D0, *D1;
TParticle *part;
TLorentzVector p0, p1, p; 
for(Int_t ipart = 0; ipart < Npart; ipart++){
 part = static_cast<TParticle*> (pArray->At(ipart));
// FillHist("Stat_Efficiency", 12);
 if(part->GetPdgCode() == 111){
// FillHist("Stat_Efficiency", 13);
  if(part->GetDaughter(0) > -1 && part->GetDaughter(1) > -1){
// FillHist("Stat_Efficiency", 14);
   D0 = static_cast<TParticle*> (pArray->At(part->GetDaughter(0)));
   D1 = static_cast<TParticle*> (pArray->At(part->GetDaughter(1)));
   if(D0->GetPdgCode() == 22 && D1->GetPdgCode() == 22){
// FillHist("Stat_Efficiency", 15);
    D0->Momentum(p0);
    D1->Momentum(p1);
    p = p0+p1;
    FillHist("PionDaughter_M", sqrt(p.E()*p.E()-p.Px()*p.Px()-p.Py()*p.Py()-p.Pz()*p.Pz()), sqrt(p.Px()*p.Px()+p.Py()*p.Py()));
}
  }
 part->Momentum(p);
 FillHist("Pion_MCalc", sqrt(p.E()*p.E()-p.Px()*p.Px()-p.Py()*p.Py()-p.Pz()*p.Pz()), sqrt(p.Px()*p.Px()+p.Py()*p.Py())); 
 FillHist("Pion_M", part->GetCalcMass(), part->Pt());
 FillHist("Pion_Pz", part->Pz());
 FillHist("Pion_Pt", part->Pt());
 FillHist("Pion_E", part->Energy());
 }
}
//pi0 = 111, eta = 221, f0 = 10221
}

//____________________Break down photons into catagories based on their mothers_______//
void AliAnalysisTaskThermalGAFlowMC::TasteFlavor()
{
const TObjArray *pArray = fStack->Particles();
Int_t Npart = pArray->GetEntriesFast();
TParticle *part, *mother;
for(Int_t ipart = 0; ipart < Npart; ipart++){
 part = static_cast<TParticle*> (pArray->At(ipart));
 if(part->GetPdgCode() == (22)){
  if(part->IsPrimary()){FillHist("Stat_Flavor_dir", part->Pt());} //direct gamma
  else{
   mother = static_cast<TParticle*> (pArray->At(part->GetMother(0)));
   if(mother->GetPdgCode() == 111 && mother->IsPrimary()){FillHist("Stat_Flavor_pi", part->Pt());} //primary pi0
   if(mother->GetPdgCode() == 221){FillHist("Stat_Flavor_eta", part->Pt());} //eta
   if(mother->GetPdgCode() == 10221){FillHist("Stat_Flavor_f0", part->Pt());} //f0(1710)
   if(mother->GetPdgCode() == 30223){FillHist("Stat_Flavor_omega", part->Pt());}  //omega(1650)
   if(mother->GetPdgCode() == 111 && mother->GetMother(0)>-1 && (static_cast<TParticle*>(pArray->At(mother->GetMother(0))))->GetPdgCode() == 311){FillHist("Stat_Flavor_kaon", part->Pt());} //K0 short (or long, just to be safe)
  }
  FillHist("Stat_Flavor_total", part->Pt());
}}
}

//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//                                                                                  //
//                                Helper Functions                                  //
//                                                                                  //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////


//____________________________Keep a record of events so the user knows what's going on_____//
void AliAnalysisTaskThermalGAFlowMC::CaptainsLog(Int_t s){
if((fDebug == 1)){printf("Step %d Completed.\n", s);}
}

//______________________________Method to Fill Histo as Seg Fualt workaround_________//
void AliAnalysisTaskThermalGAFlowMC::FillHist(const char * key, Double_t x, Double_t y) const{
  TObject * temp = fAnalist->FindObject(key);
  if(!temp) {
  AliInfo(Form("I can't find the histogram <%s> ", key));
  return;
  }
  if(temp->IsA() == TClass::GetClass("TH2F")){
  ((TH2F*)temp)->Fill(x,y);
  return;
  }
  if(temp->IsA() == TClass::GetClass("TH2I")){
  ((TH2I*)temp)->Fill(x,y);
  return;
  }
}
void AliAnalysisTaskThermalGAFlowMC::FillHist(const char * key, Double_t x) const{
  TObject * temp = fAnalist->FindObject(key);
  if(!temp) {
    AliInfo(Form("I can't find the histogram <%s>", key));
    return;
  }
  if(temp->IsA() == TClass::GetClass("TH1F")) {
    ((TH1F*)temp)->Fill(x);
    return;
  }
  if(temp->IsA() == TClass::GetClass("TH1I")){
  ((TH1I*)temp)->Fill(x);
  return;
  }
}
void AliAnalysisTaskThermalGAFlowMC::FillHist(const char * key, Double_t x, Double_t y, Double_t z) const{
  TObject * temp = fAnalist->FindObject(key);
  if(!temp) {
    AliInfo(Form("I can't find the histogram <%s>", key));
    return;
  }
  if(temp->IsA() == TClass::GetClass("TH3F")) {
    ((TH3F*)temp)->Fill(x,y,z);
   return;
  }
}
