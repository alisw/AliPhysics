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

using namespace std;

ClassImp(AliAnalysisTaskThermalGAFlow);
//_______________________________Constructor__________________________//

AliAnalysisTaskThermalGAFlow::AliAnalysisTaskThermalGAFlow( const char *name) :
 AliAnalysisTaskSE(name),
 fRunNumber(-1),
 fAnalist(0),
 fCompass(0),
 fPHOSgeomU(0x0),
 fPHOSgeom(0),
 fEvent(0x0),
 fesd(0x0),
 faod(0x0),
 fvertexx(),
 fcentrality(0),
 fNonLinearCorr(0x0),

//Task Arguments
 fDebug(0),

// CUTS! Note:  Should be overwritten in addtask, but still parameterized here for old tasks.//
 fMinCells(3),
 fMinE(0.3),
 fMinTrackDr(0.3),
 fMaxVertexx(10),
 fMinCentrality(-1),
 fMaxCentrality(100),
 fCoreRadius(3.5),
 fMinCoreEnergyRatio(0.4)
// END CUTS! //

{
  AliInfo("*** CONSTRUCTOR ***");
  DefineOutput(1, TList::Class());
}

//_______________________________Deconstructor____________________________//
AliAnalysisTaskThermalGAFlow::~AliAnalysisTaskThermalGAFlow(){
if (fAnalist) {delete fAnalist;  fAnalist = 0;}
}

//_______________________________Create Output (Histogram TList)___________//

void AliAnalysisTaskThermalGAFlow::UserCreateOutputObjects() {
  AliInfo("*** UserCreateOutputObjects ***");
  AliAnalysisTaskSE::UserCreateOutputObjects();
  OpenFile(1);

  fAnalist = new TList();
  fAnalist->SetOwner();

  TH2I *Cell_N0Clus_MB = new TH2I("Cell_N0Clus_MB", "Clusters per Cell M0 (MB)", 64,0,64,56,0,56);
  Cell_N0Clus_MB->GetXaxis()->SetTitle("Row");
  Cell_N0Clus_MB->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N0Clus_MB);
//  fCompass.push_back();

  TH2I *Cell_N0Clus = new TH2I("Cell_N0Clus", "Clusters per Cell M0", 64,0,64,56,0,56);
  Cell_N0Clus->GetXaxis()->SetTitle("Row");
  Cell_N0Clus->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N0Clus);

  TH2I *Cell_N1Clus_MB = new TH2I("Cell_N1Clus_MB", "Clusters per Cell M1 (MB)", 64,0,64,56,0,56);
  Cell_N1Clus_MB->GetXaxis()->SetTitle("Row");
  Cell_N1Clus_MB->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N1Clus_MB);

  TH2I *Cell_N1Clus = new TH2I("Cell_N1Clus", "Clusters per Cell M1", 64,0,64,56,0,56);
  Cell_N1Clus->GetXaxis()->SetTitle("Row");
  Cell_N1Clus->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N1Clus);

  TH2I *Cell_N2Clus_MB = new TH2I("Cell_N2Clus_MB", "Clusters per Cell M2 (MB)", 64,0,64,56,0,56);
  Cell_N2Clus_MB->GetXaxis()->SetTitle("Row");
  Cell_N2Clus_MB->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N2Clus_MB);

  TH2I *Cell_N2Clus = new TH2I("Cell_N2Clus", "Clusters per Cell M2", 64,0,64,56,0,56);
  Cell_N2Clus->GetXaxis()->SetTitle("Row");
  Cell_N2Clus->GetYaxis()->SetTitle("Col");
  fAnalist->Add(Cell_N2Clus);

  TH1F *Clus_Trackdr = new TH1F("Clus_Trackdr", "2d extrapolation to nearest track", 15000, -1, 1500);
  Clus_Trackdr->GetXaxis()->SetTitle("dr (cm)");
  Clus_Trackdr->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_Trackdr);

  TH1F *Clus_TOF = new TH1F("Clus_TOF", "Time of Flight", 2000, -0.00001, 0.00001);
  Clus_TOF->GetXaxis()->SetTitle("TOF (s)");
  Clus_TOF->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_TOF);

  TH1F *Clus_Pt_MB = new TH1F("Clus_Pt_MB", "Single Photon Pt Spectrum (MB)", 150, 0.5, 6.5);
  Clus_Pt_MB->GetXaxis()->SetTitle("Pt (GeV/c)");
  Clus_Pt_MB->GetYaxis()->SetTitle("dN/dPt");
  fAnalist->Add(Clus_Pt_MB);

  TH1F *Clus_Pt = new TH1F("Clus_Pt", "Single Photon Pt Spectrum", 150, 0.5, 6.5);
  Clus_Pt->GetXaxis()->SetTitle("Pt (GeV/c)");
  Clus_Pt->GetYaxis()->SetTitle("dN/dPt");
  fAnalist->Add(Clus_Pt);

  TH1F *Event_NClus_MB = new TH1F("Event_NClus_MB", "N Clusters per Event (MB)", 200, 0, 200);
  Event_NClus_MB->GetXaxis()->SetTitle("N Clusters");
  Event_NClus_MB->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_NClus_MB);

  TH1F *Event_NClus = new TH1F("Event_NClus", "N Clusters per Event", 200, 0, 200);
  Event_NClus->GetXaxis()->SetTitle("N Clusters");
  Event_NClus->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_NClus);

  TH1F *Event_Vertexx = new TH1F("Event_Vertexx", "Vertex Position Along Beam Axis", 100, -50, 50);
  Event_Vertexx->GetXaxis()->SetTitle("Vertex X");
  Event_Vertexx->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_Vertexx);

  TH1F *Event_Cent = new TH1F("Event_Cent", "Event V0 Centrality percentile", 51, -2, 100);
  Event_Cent->GetXaxis()->SetTitle("Centrality");
  Event_Cent->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_Cent);

  TH2F *Event_N0ClusVsCent = new TH2F("Event_N0ClusVsCent","NClus vs Cent", 51, -2, 100, 200, 0, 200);
  Event_N0ClusVsCent->GetXaxis()->SetTitle("Cent");
  Event_N0ClusVsCent->GetYaxis()->SetTitle("NClus");
  fAnalist->Add(Event_N0ClusVsCent);

  TH2F *Event_N1ClusVsCent = new TH2F("Event_N1ClusVsCent","NClus vs Cent", 51, -2, 100, 200, 0, 200);
  Event_N1ClusVsCent->GetXaxis()->SetTitle("Cent");
  Event_N1ClusVsCent->GetYaxis()->SetTitle("NClus");
  fAnalist->Add(Event_N1ClusVsCent);

  TH2F *Event_N2ClusVsCent = new TH2F("Event_N2ClusVsCent","NClus vs Cent", 51, -2, 100, 200, 0, 200);
  Event_N2ClusVsCent->GetXaxis()->SetTitle("Cent");
  Event_N2ClusVsCent->GetYaxis()->SetTitle("NClus");
  fAnalist->Add(Event_N2ClusVsCent);

  TH2F *Event_NClusVsCent = new TH2F("Event_NClusVsCent","NClus vs Cent", 51, -2, 100, 200, 0, 200);
  Event_NClusVsCent->GetXaxis()->SetTitle("Cent");
  Event_NClusVsCent->GetYaxis()->SetTitle("NClus");
  fAnalist->Add(Event_NClusVsCent);

  TH2F *Cluster_EnergyVsNCells = new TH2F("Cluster_EnergyVsNCells", "ClusE", 60, 0, 60, 100, 0, 20);
  Cluster_EnergyVsNCells->GetXaxis()->SetTitle("NCells");
  Cluster_EnergyVsNCells->GetYaxis()->SetTitle("Cluster E");
  fAnalist->Add(Cluster_EnergyVsNCells);

  TH2F *Cluster_CoreEVsNCells = new TH2F("Cluster_CoreEVsNCells", "CoreE", 60, 0, 60, 100, 0 , 20);
  Cluster_CoreEVsNCells->GetXaxis()->SetTitle("NCells");
  Cluster_CoreEVsNCells->GetYaxis()->SetTitle("Core Energy");
  fAnalist->Add(Cluster_CoreEVsNCells);

  TH2F *Cluster_CoreERatioVsNCells = new TH2F("Cluster_CoreERatioVsNCells", "Ratio", 60, 0, 60, 100, 0, 1);
  Cluster_CoreERatioVsNCells->GetXaxis()->SetTitle("NCells");
  Cluster_CoreERatioVsNCells->GetYaxis()->SetTitle("Ratio");
  fAnalist->Add(Cluster_CoreERatioVsNCells);

  TH2F *Cluster_AvgShape = new TH2F("Cluster_AvgShape", "Shape", 40, 0, 15, 100, 0, 20);
  Cluster_AvgShape->GetXaxis()->SetTitle("Distance To Center");
  Cluster_AvgShape->GetYaxis()->SetTitle("Cell Energy (before nonlinear corr)");
  fAnalist->Add(Cluster_AvgShape);

  PostData(1, fAnalist);
}

//______________________________Main Loop, Called for each Event_____________//

void AliAnalysisTaskThermalGAFlow::UserExec( Option_t *){

//Step 0: Configure the environment - take out your toys
ConfigureEvent();
if(fRunNumber != fEvent->GetRunNumber()){fRunNumber = fEvent->GetRunNumber(); InitializeGeometry();}
ScanBasicEventParameters();
CaptainsLog(0);

//Step 1: Internal Triggering and QA - set up your toys
if(!ApplyEventCuts()){PostData(1,fAnalist);}
else {
ScanClusters();
CaptainsLog(1);

//Step 2: Extract Inclusive Photon Spectrum

CaptainsLog(2);

PostData(1,fAnalist);
}
}

//_____________________________Final step,  Called once at the end of the Task_________//

void AliAnalysisTaskThermalGAFlow::Terminate( Option_t *){

}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                                                          //
//                               Functions                                  //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


//_____________________________Get the event and put it in the proper format_____________//
void AliAnalysisTaskThermalGAFlow::ConfigureEvent(){
fEvent = InputEvent();
fesd = dynamic_cast<AliESDEvent*> (fEvent);
faod = dynamic_cast<AliAODEvent*> (fEvent);

if(!fEvent){
  AliError("Holdup, I couldn't find the event.");
  PostData(1,fAnalist);
}
if(!fesd && !faod){
  AliError("Holdup, The event isn't AOD or ESD.");
  PostData(1,fAnalist);
}

const char *nonlinform = "1.015*(0.0241+1.0504*x+0.000249*x*x)";
fNonLinearCorr = new TF1("nonlin", nonlinform, 0, 40); //Set the Nonlinearity Correction
}

//____________________________Load up the detector Geometry________________//
void AliAnalysisTaskThermalGAFlow::InitializeGeometry(){
fPHOSgeom = AliPHOSGeometry::GetInstance("IHEP");
fPHOSgeomU = new AliPHOSGeoUtils("name", "title");
fPHOSgeomU->SetMisalMatrix(fesd->GetPHOSMatrix(0),0);
fPHOSgeomU->SetMisalMatrix(fesd->GetPHOSMatrix(1),1);
fPHOSgeomU->SetMisalMatrix(fesd->GetPHOSMatrix(2),2);
}

//____________________________Extract the basic perameters for the event______//
void AliAnalysisTaskThermalGAFlow::ScanBasicEventParameters(){ 
  const AliVVertex *primaryVertex = fEvent->GetPrimaryVertex();
  if(primaryVertex){primaryVertex->GetXYZ(fvertexx);}
  else {printf("No primary vertex. Setting it to Zero Vector."); fvertexx[0] = 0; fvertexx[1] = 0; fvertexx[2] = 0;}
  if(fEvent->GetCentrality()){fcentrality = fesd->GetCentrality()->GetCentralityPercentile("V0M");}
  else {printf("No V0C centrality found.  Setting it to -2."); fcentrality = -2;} 
}

//___________________________Extract useful information from the clusters______________//
void AliAnalysisTaskThermalGAFlow::ScanClusters(){

Int_t nPHOSClusters = 0;
Int_t nGoodClusters = 0;
Int_t x[4];
TLorentzVector p;

Int_t nClus0 = 0;
Int_t nClus1 = 0;
Int_t nClus2 = 0;

if(fesd){

for(Int_t icluster = 0; icluster < fesd->GetNumberOfCaloClusters(); icluster++) {

  AliESDCaloCluster *cluster = fesd->GetCaloCluster(icluster);
  if(cluster->IsPHOS()){
    GetPos(cluster, x);
    cluster->GetMomentum(p, fvertexx);
    nPHOSClusters = nPHOSClusters+1;

    if(x[0] == 0){FillHist("Cell_N0Clus_MB", x[2], x[3]); nClus0 = nClus0+1;}
    if(x[0] == 1){FillHist("Cell_N1Clus_MB", x[2], x[3]); nClus1 = nClus1+1;}
    if(x[0] == 2){FillHist("Cell_N2Clus_MB", x[2], x[3]); nClus2 = nClus2+1;}

    FillHist("Clus_Pt_MB", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
  }

 FillHist("Event_N0ClusVsCent", fcentrality, nClus0);
 FillHist("Event_N1ClusVsCent", fcentrality, nClus1);
 FillHist("Event_N2ClusVsCent", fcentrality, nClus2);
 FillHist("Event_NClusVsCent", fcentrality, nClus0+nClus1+nClus2);

  if(ApplyClusterCuts(cluster)){
    if(ApplyCPCuts(cluster)){
      if(ApplyShapeCuts(cluster)){

    nGoodClusters = nGoodClusters+1;
    if(x[0] == 0){FillHist("Cell_N0Clus", x[2], x[3]);}
    if(x[0] == 1){FillHist("Cell_N1Clus", x[2], x[3]);}
    if(x[0] == 2){FillHist("Cell_N2Clus", x[2], x[3]);}
    FillHist("Clus_Pt", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
    FillHist("Clus_TOF", cluster->GetTOF());

    SavePhoton(cluster);
  }}}

}}
if(faod){
for(Int_t icluster = 0; icluster < faod->GetNumberOfCaloClusters(); icluster++) {
  AliAODCaloCluster *cluster = faod->GetCaloCluster(icluster);
  GetPos(cluster, x);
  cluster->GetMomentum(p, fvertexx);
  if(cluster->IsPHOS()){
    nPHOSClusters = nPHOSClusters+1;
    if(x[0] == 0){FillHist("Cell_N0Clus_MB", x[2], x[3]);}
    if(x[0] == 1){FillHist("Cell_N1Clus_MB", x[2], x[3]);}
    if(x[0] == 2){FillHist("Cell_N2Clus_MB", x[2], x[3]);}
    FillHist("Clus_Pt_MB", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
  }
  if(ApplyClusterCuts(cluster)){
    nGoodClusters = nGoodClusters+1;
    if(x[0] == 0){FillHist("Cell_N0Clus", x[2], x[3]);}
    if(x[0] == 1){FillHist("Cell_N1Clus", x[2], x[3]);}
    if(x[0] == 2){FillHist("Cell_N2Clus", x[2], x[3]);}
    FillHist("Clus_Pt", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
    FillHist("Clus_TOF", cluster->GetTOF());
  }
}}

FillHist("Event_NClus_MB", nPHOSClusters);
FillHist("Event_NClus", nGoodClusters);

}

//___________________________Fill Module Specific Data______________________________//
void AliAnalysisTaskThermalGAFlow::ScanMod(Int_t* x, TLorentzVector p){
if(x[0] == 0){
 FillHist("Cell_N0Clus_MB", x[2], x[3]);
// FillHist("Event_N0ClusVsCent", fcentrality, fesd->GetNumberOfCaloClusters());
}
if(x[0] == 1){
 FillHist("Cell_N1Clus_MB", x[2], x[3]);
// FillHist("Event_N1ClusVsCent", fcentrality, fesd->GetNumberOfCaloClusters());
}
if(x[0] == 2){
 FillHist("Cell_N2Clus_MB", x[2], x[3]);
// FillHist("Event_N2ClusVsCent", fcentrality, fesd->GetNumberOfCaloClusters());
}

}

//___________________________Apply Event Cuts________________________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyEventCuts(){
if(fDebug == 1){printf("Vert: %f ;; Cent: %f \n", fvertexx[2], fcentrality);}

if(!(TMath::Abs(fvertexx[2]) < fMaxVertexx)){return 0;}
if(!((fcentrality > fMinCentrality) && (fcentrality < fMaxCentrality))){return 0;}
if(!(fesd->GetNumberOfCaloClusters() > 3)){return 0;}

//Temporary Mod 1 Corruption Cut
AliESDCaloCells* PHOSCells = fesd->GetPHOSCells();
Int_t x[4] = {2, 0, 0, 0}; //mod, (PbW04 = 0 or CPV = -1), row, col
Int_t absID = 0;
Int_t CorruptFlag = 1;
for(Int_t row = 0; row<=64; row++){
for(Int_t col = 0; col<=28; col++){
  x[2] = row;
  x[3] = col;
  if(fPHOSgeomU->RelToAbsNumbering(x, absID)){
    if(PHOSCells->GetCellAmplitude(absID) != 0){CorruptFlag = 0;}
  }
}}
if(CorruptFlag){return 0;}
//End Corruption Cut

FillHist("Event_Vertexx", fvertexx[2]);
FillHist("Event_Cent",fcentrality) ;

return 1;

}

//___________________________Apply Cluster Cuts______________________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyClusterCuts(AliESDCaloCluster *cluster){
if(!(cluster->IsPHOS())){return 0;}
if(!(cluster->GetNCells() >= fMinCells)){return 0;}
if(!(cluster->E() >= fMinE)){return 0;}
if(!(cluster->GetTOF() > -0.0000002 && cluster->GetTOF() < 0.0000002)){return 0;}
return 1;
}

Bool_t AliAnalysisTaskThermalGAFlow::ApplyClusterCuts(AliAODCaloCluster *cluster){
if(!(cluster->IsPHOS())){return 0;}
if(!(cluster->GetNCells() >= fMinCells)){return 0;}
if(!(cluster->E() >= fMinE)){return 0;}
return 1;
}

//___________________________Apply Charged Particle Cuts____________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyCPCuts(AliESDCaloCluster *cluster){

Double_t dx = cluster->GetTrackDx();
Double_t dz = cluster->GetTrackDz();

//  if(dx != 0){printf("Not Zero\n");}
//  printf("%f\n", dx);
//AliESDtrack* matchtrack = fesd->GetTrack(cluster->GetTracksMatched()->At(1));
//Int_t tracknum = -7;
//TArrayI* trackarr = cluster->GetTracksMatched();
//Int_t* arr = trackarr->GetArray();
//if(arr){printf("I got the array, it's %i\n", arr[0]);}
//tracknum = (cluster->GetTracksMatched())->At(0);
//printf("Tracknum: %d\n", tracknum);
//if(matchtrack){printf("Got a track, at least.\n");}


  Double_t trackDr = sqrt(dx*dx + dz*dz);
//  printf("trackDr: %f ; trackDx: %f . \n", trackDr, cluster->GetTrackDx());
  FillHist("Clus_Trackdr", trackDr);
  if(!(trackDr > fMinTrackDr)){return 0;}
  return 1;
}

//___________________________Apply Cluster Shape Cuts______________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyShapeCuts(AliESDCaloCluster *cluster){
Int_t nCells = cluster->GetNCells();
Int_t realID[4];
Float_t x, z, xc, zc;
Double32_t xzematrix[nCells][3]; //Nx3 matrix consisting of the full position and energy information for a given even.
Double32_t coreE = 0;
Double_t cellDistance;
AliVCaloCells *cells = dynamic_cast<AliVCaloCells*> (fEvent->GetPHOSCells());
Double_t cellETot = 0;

GetPos(cluster, realID);
realID[0] = realID[0] + 1; //Have to undo the GetPos correction to give a valid RelId again.
realID[1] = realID[1] + 1;
realID[2] = realID[2] + 1;
realID[3] = realID[3] + 1;
fPHOSgeomU->RelPosInModule(realID, xc, zc);

for(Int_t cellIndex = 0; cellIndex < nCells; cellIndex++){
  if(!fPHOSgeomU->AbsToRelNumbering(cluster->GetCellAbsId(cellIndex), realID)){
    printf("Failed to map Absolute ID to Real ID in Abs ID %i \n", cluster->GetCellAbsId(cellIndex));
    return 0;
  }  

  cellETot = cells->GetCellAmplitude(cluster->GetCellAbsId(cellIndex));
  fPHOSgeomU->RelPosInModule(realID, x, z);

  xzematrix[cellIndex][0] = x - xc;
  xzematrix[cellIndex][1] = z - zc;
  xzematrix[cellIndex][2] = cellETot * cluster->GetCellAmplitudeFraction(cellIndex);

  cellDistance = TMath::Sqrt(xzematrix[cellIndex][0]*xzematrix[cellIndex][0]+xzematrix[cellIndex][1]*xzematrix[cellIndex][1]);

  FillHist("Cluster_AvgShape", cellDistance, xzematrix[cellIndex][2]);
  
  if(cellDistance < fCoreRadius){
  coreE = coreE + xzematrix[cellIndex][2];
  }
} 

coreE = fNonLinearCorr->Eval(coreE); //Don't forget to apply the nonlinearity correction.

FillHist("Cluster_EnergyVsNCells", nCells, cluster->E());
FillHist("Cluster_CoreEVsNCells", nCells, coreE);
FillHist("Cluster_CoreERatioVsNCells", nCells, coreE/cluster->E());

if(coreE / cluster->E() > fMinCoreEnergyRatio){return true;}
else {return false;}
}

//___________________________Save Photon Candidate_________________________________//
void AliAnalysisTaskThermalGAFlow::SavePhoton(AliESDCaloCluster* cluster){
if(cluster){/*printf("I got a photon!\n");*/}
}

//___________________________Get the Position of the center of the cluster_________//
Int_t *AliAnalysisTaskThermalGAFlow::GetPos(AliESDCaloCluster* cluster, Int_t x[4]){
  Float_t cellGlobalPosArr[3];
  cluster->GetPosition(cellGlobalPosArr);
  TVector3 globalPos(cellGlobalPosArr);
  if(!fPHOSgeomU->GlobalPos2RelId(globalPos, x)){
    printf("PHOS is not at the global pos given: %f, %f, %f \n", globalPos[0], globalPos[1], globalPos[2]);}
  x[0] = x[0] - 1; //For some reason, RelId counts from 1 instead of the normal 0. This corrects it. 
  x[1] = x[1] - 1;
  x[2] = x[2] - 1;
  x[3] = x[3] - 1;
  return x;
}

Int_t *AliAnalysisTaskThermalGAFlow::GetPos(AliAODCaloCluster* cluster, Int_t x[4]){
  if(!cluster){Printf("No AOD Cluster.");}
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 0; 
  printf("Unable to handle AOD global position information.  Returning Zero Vector.");
  return x;
}



//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//                                                                                  //
//                                Helper Functions                                  //
//                                                                                  //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////


//____________________________Keep a record of events so the user knows what's going on_____//
void AliAnalysisTaskThermalGAFlow::CaptainsLog(Int_t s){
if(fDebug == 1){printf("Step %d Completed.\n", s);}
}

//______________________________Method to Fill Histo as Seg Fualt workaround_________//
void AliAnalysisTaskThermalGAFlow::FillHist(const char * key, Double_t x, Double_t y) const{
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
void AliAnalysisTaskThermalGAFlow::FillHist(const char * key, Double_t x) const{
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

void AliAnalysisTaskThermalGAFlow::FillHist(const char * key, Double_t x, Double_t y, Double_t z) const{
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
