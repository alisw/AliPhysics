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

  Int_t Nbins = 150; //("Clus_Pt", "Single Photon Pt Spectrum", 150, 0.5, 20.5); and exp scaling
  Double_t Ptbins[151];
  Double_t a = TMath::Power(10, 1);
  for(Int_t ibin = 0; ibin < 151; ibin++){
  Ptbins[ibin] = 0.5 - (20/(a - 1)) + ((20/(a - 1)) * TMath::Power(10, ibin/150.));
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

  TH1I *Stat_Purity = new TH1I("Stat_Purity", "Purity of each cut number", 20, 0, 20);
  Stat_Purity->GetXaxis()->SetTitle("Counts");
  Stat_Purity->GetYaxis()->SetTitle("Cut Number");
  fAnalist->Add(Stat_Purity);

PostData(1,fAnalist);
}

//_______________________________Main Loop_________________________________//
void AliAnalysisTaskThermalGAFlowMC::UserExec( Option_t * option){
fStack = GetMCStack();
AliAnalysisTaskThermalGAFlow::UserExec(option); //Do the Same thing as Data
ScanClustersMC();  //Then extract the needed MC information.
ExtractRealPions(); //Extract the real pion distribution, no cuts.

PostData(1,fAnalist);
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
  GuessPure(cluster);
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

if(WeakCuts(cluster)){
 if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 1);}

 if(cluster->GetTOF() > -0.0000001 && cluster->GetTOF() < 0.0000001){
  if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 2);}

  if(cluster->GetDistanceToBadChannel() > 4.4){
   if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 3);}

   if(ApplyCPCuts(cluster)){
    if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 4);}

    if(ApplyCoreShapeCut(cluster)){
     if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 5);}

     if(ApplyDispMatrixCut(cluster)){
      if(!IsPhotonShower(cluster)){FillHist("Stat_Purity", 6);}   
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
if(pCode == 22 && primaryERatio > 0.1){return 1;} 
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
TLorentzVector p0, p1, p; 
for(Int_t ipart = 0; ipart < Npart; ipart++){
 TParticle *part = static_cast<TParticle*> (pArray->At(ipart));
 if(part->GetPdgCode() == (111/* || 221*/)){
  if(part->GetDaughter(0) >= 1 && part->GetDaughter(1) >= 1){
   D0 = static_cast<TParticle*> (pArray->At(part->GetDaughter(0)));
   D1 = static_cast<TParticle*> (pArray->At(part->GetDaughter(1)));
//   if(D0->GetPdgCode() == 22 && D1->GetPdgCode() == 22){
    D0->Momentum(p0);
    D1->Momentum(p1);
    p = p0+p1;
    FillHist("PionDaughter_M", sqrt(p.E()*p.E()-p.Px()*p.Px()-p.Py()*p.Py()-p.Pz()*p.Pz()), sqrt(p.Px()*p.Px()+p.Py()*p.Py()));
//}
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
