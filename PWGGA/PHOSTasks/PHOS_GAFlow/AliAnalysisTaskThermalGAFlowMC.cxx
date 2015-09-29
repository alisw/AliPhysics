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

  TH1F *Clus_PrimeERatio = new TH1F("Clus_PrimeERatio", "Ratio of Highest Energy Prime to total cluster Energy", 50, 0, 1);
  Clus_PrimeERatio->GetXaxis()->SetTitle("PrimeERatio");
  Clus_PrimeERatio->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_PrimeERatio);

  TH1F *Clus_PtGamma = new TH1F("Clus_PtGamma", "Pt of the MC Gammas", 150, 0.5 , 20.5);
  Clus_PtGamma->GetXaxis()->SetTitle("Pt");
  Clus_PtGamma->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_PtGamma);

PostData(1,fAnalist);
}

//_______________________________Main Loop_________________________________//
void AliAnalysisTaskThermalGAFlowMC::UserExec( Option_t * option){
fStack = GetMCStack();
AliAnalysisTaskThermalGAFlow::UserExec(option); //Do the Same thing as Data
ScanClustersMC();  //Then extract the needed MC information.

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
printf("%p -- %p \n", mcHandler, static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent());
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
  if(ApplyClusterCuts(cluster)){
  if(IsPhotonShower(cluster)){
    cluster->GetMomentum(p, fvertexx);
    FillHist("Clus_PtGamma", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
  }}
}}

}

//___________________Check to see if a cluster was from a MC photon________//
Bool_t AliAnalysisTaskThermalGAFlowMC::IsPhotonShower(AliESDCaloCluster *cluster)
{
//Note:  I'm not really sure, but I think primaries are MC particles that contributed energy to the shower.  Labels, then, is the primaries sorted by deposited energy - Label[0] would be the highest contribution.
TParticle *p = fStack->Particle(cluster->GetLabelAt(0));
Int_t pCode = p->GetPdgCode(); //This is the MC PID code.
Float_t primaryERatio = p->Energy() / cluster->E(); //Prime Energy Ratio.  PHOS_PbPb used 0.9.  I'm going to use something smaller to start.
FillHist("Clus_PrimeERatio", primaryERatio);
if(pCode == 22 && primaryERatio > 0.5){return 1;} //22 is the photon code... I think.  11 is electron, -11 is position... again... I think.
return 0;
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
