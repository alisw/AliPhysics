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
//#include <vector>

#include "TNtupleD.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "TList.h"
#include "TString.h"
#include "AliLog.h"

#include "AliCaloPhoton.h"
#include "AliEventplane.h"
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
 fAtlas(0),
 fCompass(0x0),
 fSextant(0x0),
 fphoton(0x0),
 fmixphoton(0x0),
 fPHOSgeomU(0x0),
 fPHOSgeom(0),
 fEvent(0x0),
 fesd(0x0),
 faod(0x0),
 fvertexx(),
 fcentrality(0),
 fevPlane(0),
 fevPlaneA(0),
 fevPlaneC(0),
 fevScaler(0),
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
 fMinCoreEnergyRatio(0.4),
 fMinLambdaDisp(0.3),
 fMinCPVStd(0),

 fMixVertxbins(1),
 fMixCentbins(1),
 fMixEvbins(1),
 fNptbins(150)
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

  Int_t CompassBinning = fMixVertxbins*fMixCentbins*fMixEvbins; //Vertex Binning * Centrality Binning * ev Plane Binning
  fCompass = new TObjArray(CompassBinning);
  fCompass->SetOwner();

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

  TH1F *Clus_Trackdr = new TH1F("Clus_Trackdr", "2d extrapolation to nearest track", 1010, -1, 100);
  Clus_Trackdr->GetXaxis()->SetTitle("dr (cm)");
  Clus_Trackdr->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_Trackdr);

  TH1F *Clus_TrackMatchStd = new TH1F("Clus_TrackMatchStd", "Track Match Std", 100, 0, 20);
  Clus_TrackMatchStd->GetXaxis()->SetTitle("TrackMatchStd");
  Clus_TrackMatchStd->GetYaxis()->SetTitle("dN/do");
  fAnalist->Add(Clus_TrackMatchStd);

  TH1F *Clus_TOF = new TH1F("Clus_TOF", "Time of Flight", 2000, -0.00001, 0.00001);
  Clus_TOF->GetXaxis()->SetTitle("TOF (s)");
  Clus_TOF->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Clus_TOF);

  Int_t Nbins = fNptbins; //("Clus_Pt", "Single Photon Pt Spectrum", 150, 0.5, 20.5); and exp scaling
  const Int_t somebins = Nbins+1;
  Double_t Ptbins[somebins];
  Double_t a = TMath::Power(10, 1);
  for(Int_t ibin = 0; ibin < somebins; ibin++){
  Ptbins[ibin] = 0.5 - (20/(a - 1)) + ((20/(a - 1)) * TMath::Power(10, ibin/(Nbins*1.0)));
  }

  TH1F *Clus_Pt_MB = new TH1F("Clus_Pt_MB", "Single Photon Pt Spectrum (MB)", Nbins, Ptbins);
  Clus_Pt_MB->GetXaxis()->SetTitle("Pt (GeV/c)");
  Clus_Pt_MB->GetYaxis()->SetTitle("dN/dPt");
  fAnalist->Add(Clus_Pt_MB);

  TH1F *Clus_Pt = new TH1F("Clus_Pt", "Single Photon Pt Spectrum", Nbins, Ptbins);
  Clus_Pt->GetXaxis()->SetTitle("Pt (GeV/c)");
  Clus_Pt->GetYaxis()->SetTitle("dN/dPt");
  fAnalist->Add(Clus_Pt);

  TH1F *Clus_BCDistance = new TH1F("Clus_BCDistance", "Distance between maximum and nearest bad channel", 1000, 0, 100);
  Clus_BCDistance->GetXaxis()->SetTitle("Distance");
  Clus_BCDistance->GetYaxis()->SetTitle("dN/dD");
  fAnalist->Add(Clus_BCDistance);

  TH1F *Clus_DispLambda = new TH1F("Clus_DispLambda", "Dispersion Matrix Eigenvalue", 200, 0, 20);
  Clus_DispLambda->GetXaxis()->SetTitle("Dispersion Matrix Eigenvalue");
  Clus_DispLambda->GetYaxis()->SetTitle("dN/dD");
  fAnalist->Add(Clus_DispLambda);

  TH1F *Event_Vertexx = new TH1F("Event_Vertexx", "Vertex Position Along Beam Axis", 100, -50, 50);
  Event_Vertexx->GetXaxis()->SetTitle("Vertex X");
  Event_Vertexx->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_Vertexx);

  TH1F *Event_Cent = new TH1F("Event_Cent", "Event V0 Centrality percentile", 51, -2, 100);
  Event_Cent->GetXaxis()->SetTitle("Centrality");
  Event_Cent->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_Cent);

  TH1F *Event_EP = new TH1F("Event_EP", "Event Plane Angle of rotation", 100, 0, TMath::Pi());
  Event_EP->GetXaxis()->SetTitle("Event Plane Psi");
  Event_EP->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_EP);

  TH1F *Event_EPSP = new TH1F("Event_EPSP", "Event Plane Scaler Product", 100, 0, 1);
  Event_EPSP->GetXaxis()->SetTitle("|cos(Psi_V0A - Psi_V0C)|");
  Event_EPSP->GetYaxis()->SetTitle("Counts");
  fAnalist->Add(Event_EPSP);

  TH2F *Event_EPAC = new TH2F("Event_EPAC", "Event Plane A and C Sides", 100, 0, TMath::Pi(), 100, 0, TMath::Pi()); 
  Event_EPAC->GetXaxis()->SetTitle("V0A EP");
  Event_EPAC->GetYaxis()->SetTitle("V0C_EP");
  fAnalist->Add(Event_EPAC);

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

  TH1I *Stat_Efficiency = new TH1I("Stat_Efficiency", "Efficiency of each cut number", 20, 0, 20);
  Stat_Efficiency->GetXaxis()->SetTitle("Counts");
  Stat_Efficiency->GetYaxis()->SetTitle("Cut Number");
  fAnalist->Add(Stat_Efficiency);

  TH2F *InclGamma_PhiVsPt = new TH2F("InclGamma_PhiVsPt", "Inclusive Gamma Azimuthal yield", 200, -1*TMath::Pi(), TMath::Pi(), Nbins, Ptbins);
  InclGamma_PhiVsPt->GetXaxis()->SetTitle("Phi");
  InclGamma_PhiVsPt->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(InclGamma_PhiVsPt);

//  TH2F *DecGamma_PhiVsPt = new TH2F("DecGamma_PhiVsPt", "Decay Gamma Azimuthal yield", 200, -1*TMath::Pi(), TMath::Pi(), Nbins, Ptbins);
//  DecGamma_PhiVsPt->GetXaxis()->SetTitle("Phi");
//  DecGamma_PhiVsPt->GetYaxis()->SetTitle("Pt");
//  fAnalist->Add(DecGamma_PhiVsPt);

  TH2F *InclGamma_DPhiVsPt_A = new TH2F("InclGamma_DPhiVsPt_A", "Inclusive Gamma Azimuthal yield wrt evplane_A", 100, 0, TMath::Pi(), Nbins, Ptbins);
  InclGamma_DPhiVsPt_A->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsPt_A->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(InclGamma_DPhiVsPt_A);

//  TH2F *DecGamma_DPhiVsPt_A = new TH2F("DecGamma_DPhiVsPt_A", "Decay Gamma Azimuthal yield wrt explane_A", 100, 0, TMath::Pi(), Nbins, Ptbins);
//  DecGamma_DPhiVsPt_A->GetXaxis()->SetTitle("DPhi");
//  DecGamma_DPhiVsPt_A->GetYaxis()->SetTitle("Pt");
//  fAnalist->Add(DecGamma_DPhiVsPt_A);

  TH2F *InclGamma_DPhiVsPt_C = new TH2F("InclGamma_DPhiVsPt_C", "Inclusive Gamma Azimuthal yield wrt evplane_C", 100, 0, TMath::Pi(), Nbins, Ptbins);
  InclGamma_DPhiVsPt_C->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsPt_C->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(InclGamma_DPhiVsPt_C);

//  TH2F *DecGamma_DPhiVsPt_C = new TH2F("DecGamma_DPhiVsPt_C", "Decay Gamma Azimuthal yield wrt explane_C", 100, 0, TMath::Pi(), Nbins, Ptbins);
//  DecGamma_DPhiVsPt_C->GetXaxis()->SetTitle("DPhi");
//  DecGamma_DPhiVsPt_C->GetYaxis()->SetTitle("Pt");
//  fAnalist->Add(DecGamma_DPhiVsPt_C);

  TH2F *InclGamma_DPhiVsPt_AC = new TH2F("InclGamma_DPhiVsPt_AC", "Inclusive Gamma Azimuthal yield wrt evplane_AC (average of A and C)", 100, 0, TMath::Pi(), Nbins, Ptbins);
  InclGamma_DPhiVsPt_AC->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsPt_AC->GetYaxis()->SetTitle("Pt");
  fAnalist->Add(InclGamma_DPhiVsPt_AC);

//  TH2F *DecGamma_DPhiVsPt_AC = new TH2F("DecGamma_DPhiVsPt_AC", "Decay Gamma Azimuthal yield wrt explane_AC", 100, 0, TMath::Pi(), Nbins, Ptbins);
//  DecGamma_DPhiVsPt_AC->GetXaxis()->SetTitle("DPhi");
//  DecGamma_DPhiVsPt_AC->GetYaxis()->SetTitle("Pt");
//  fAnalist->Add(DecGamma_DPhiVsPt_AC);

  TH2F *InclGamma_DPhiVsCent_A = new TH2F("InclGamma_DPhiVsCent_A", "Inclusive Gamma Azimuthal yield wrt evplane_A", 100, 0, 25, 0, 100);
  InclGamma_DPhiVsCent_A->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsCent_A->GetYaxis()->SetTitle("Cent");
  fAnalist->Add(InclGamma_DPhiVsCent_A);

  TH2F *InclGamma_DPhiVsCent_C = new TH2F("InclGamma_DPhiVsCent_C", "Inclusive Gamma Azimuthal yield wrt evplane_C", 100, 0, TMath::Pi(), 25, 0, 100);
  InclGamma_DPhiVsCent_C->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsCent_C->GetYaxis()->SetTitle("Cent");
  fAnalist->Add(InclGamma_DPhiVsCent_C);

  TH2F *InclGamma_DPhiVsCent_AC = new TH2F("InclGamma_DPhiVsCent_AC", "Inclusive Gamma Azimuthal yield wrt evplane_AC", 100, 0, 25, 0 , 100);
  InclGamma_DPhiVsCent_AC->GetXaxis()->SetTitle("DPhi");
  InclGamma_DPhiVsCent_AC->GetYaxis()->SetTitle("Cent");
  fAnalist->Add(InclGamma_DPhiVsCent_AC);

//  TH2F *Meson_RawM = new TH2F("Meson_RawM", "Raw Diphoton Mass", 280,0,.7,500,0,25);
//  Meson_RawM->GetXaxis()->SetTitle("M (GeV/c^2)");
//  Meson_RawM->GetYaxis()->SetTitle("Pt (Px&Py in Lorentz Vector)");
//  fAnalist->Add(Meson_RawM);

//  TH2F *Meson_XM = new TH2F("Meson_XM", "Mixed Event Diphoton Mass", 280,0,.7,500,0,25);
//  Meson_XM->GetXaxis()->SetTitle("M (GeV/c^2)");
//  Meson_XM->GetYaxis()->SetTitle("Pt (Px&Py in Lorentz Vector)");
//  fAnalist->Add(Meson_XM);

  TH2F *Meson_RawM = new TH2F("Meson_RawM", "Raw Diphoton Mass", 400,0,1, Nbins, Ptbins);
  Meson_RawM->GetXaxis()->SetTitle("M (GeV/c^2)");
  Meson_RawM->GetYaxis()->SetTitle("Pt (Px&Py in Lorentz Vector)");
  fAnalist->Add(Meson_RawM);

  TH2F *Meson_XM = new TH2F("Meson_XM", "Mixed Event Diphoton Mass", 400,0,1, Nbins, Ptbins);
  Meson_XM->GetXaxis()->SetTitle("M (GeV/c^2)");
  Meson_XM->GetYaxis()->SetTitle("Pt (Px&Py in Lorentz Vector)");
  fAnalist->Add(Meson_XM);

  TH2F *Cluster_M = new TH2F("Cluster_M", "Di-Cluster Mass", 400,0,1, Nbins, Ptbins);
  Cluster_M->GetXaxis()->SetTitle("M (GeV/c^2)");
  Cluster_M->GetYaxis()->SetTitle("Pt (Px&Py in Lorentz Vector)");
  fAnalist->Add(Cluster_M);

//Begin Pt slices of mass - incomplete
/*  TH1F *fHistM_slice;
  Int_t maslicebin = 500;
  Double_t maslicemin = 0;
  Double_t maslicemax = 0.7;
  maslicenum = 31;
  maslicewidth = 0.2;

  for (Int_t slice = 0; slice < maslicenum; slice++){
    char names[10], titles[19];
    Int_t min = slice;
    Int_t max = (slice+1);
    sprintf(names, "fHistM_[%i]", slice);
    sprintf(titles[slice], "(%i*%i)<Pt(diphoton)<(%i+1)*%i", min, maslicewidth, max, maslicewidth);

    fHistM_slice = new TH1F(names[names], titles[slice], maslicebin, maslicemin, maslicemax);
    fHistM_slice->GetXaxis()->SetTitle("M (GeV/c^2)");
    fHistM_slice->GetYaxis()->SetTitle("dN/dM");
    fHistM_slice->SetMarkerStyle(kFullCircle);
    fAnalist->Add(fHistM_slice);
  }
*/
//End Pt slices of mass

  PostData(1, fAnalist);
}

//______________________________Main Loop, Called for each Event_____________//

void AliAnalysisTaskThermalGAFlow::UserExec( Option_t *){
//printf("Anything.\n");
FillHist("Stat_Efficiency", 15);
PostData(1,fAnalist);
CaptainsLog(0);

//Step 0: Configure the environment - take out your toys
ConfigureEvent();
FillHist("Stat_Efficiency", 16);
PostData(1,fAnalist);
if(fRunNumber != fEvent->GetRunNumber()){fRunNumber = fEvent->GetRunNumber(); InitializeGeometry();}
FillHist("Stat_Efficiency", 17);
PostData(1,fAnalist);
ScanBasicEventParameters();
CaptainsLog(1);

//Step 1: Internal Triggering and QA - set up your toys
if(!ApplyEventCuts()){PostData(1,fAnalist);}
else {
FillHist("Stat_Efficiency", 18);
ScanClusters();
SaveEvent();
CaptainsLog(2);

//Step 2: Extract the diphoton mass - play with your toys
//printf("fSextant = %p\n", fSextant);
MesonExclusion();
CaptainsLog(3);
//Step 3: Compute thermal photon spectrum - have fun and be happy
Hypnotic();

//Step 4: Compute thermal photon flow - put away your toys when you're done with them 
IncFlow();

//delete fSextant;
fSextant = 0x0;
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

//Now, get the Reaction plane.  3 methods to do this:  V0A, V0C, and TPC.  V0A and V0C look for a fourier transform in the zero degree calorimeter on the A and C side, respectively.  TPC looks for the fourier transform in the TPC.  Studies show (see https://indico.cern.ch/event/405255/) that the TPC method is not as reliable as V0A and V0C, which are both about equally accurate.  So we take the average of V0C and V0A for the best results.

Double_t RP;
Double_t RPscaler;

AliEventplane *eventplane = fEvent->GetEventplane();
if(!eventplane){AliError("Event Plane isn't reconstructed.  Setting to zero. \n"); fevPlane = 0; return;}
Double_t V0A_RP = eventplane->GetEventplane("V0A", fEvent);
Double_t V0C_RP = eventplane->GetEventplane("V0C", fEvent);

if(fDebug > 3 && (V0A_RP > TMath::Pi() || V0A_RP < 0)){AliInfo(Form("It seems the V0A (%f) is outside of the valid range [0,pi].  Correcting.", V0A_RP));}
if(fDebug > 3 && (V0C_RP > TMath::Pi() || V0C_RP < 0)){AliInfo(Form("It seems the V0C (%f) is outside of the valid range [0,pi].  Correcting.", V0C_RP));}

while(V0A_RP < 0){V0A_RP = V0A_RP + TMath::Pi();}
while(V0C_RP < 0){V0C_RP = V0C_RP + TMath::Pi();}
while(V0A_RP > TMath::Pi()){V0A_RP = V0A_RP - TMath::Pi();}
while(V0C_RP > TMath::Pi()){V0C_RP = V0C_RP - TMath::Pi();}

RP = (V0A_RP + V0C_RP)/2;;
RPscaler = TMath::Abs(TMath::Cos((V0A_RP - V0C_RP)));

FillHist("Event_EPAC", V0A_RP, V0C_RP);

fevPlane = RP;
fevPlaneA = V0A_RP;
fevPlaneC = V0C_RP;

fevScaler = RPscaler;

//printf("EP: %f \n", fevPlane);
}

//___________________________Apply Event Cuts________________________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyEventCuts(){
if(fDebug == 1){printf("Vert: %f ;; Cent: %f \n", fvertexx[2], fcentrality);}

if(!(TMath::Abs(fvertexx[2]) < fMaxVertexx)){return 0;}
if(!((fcentrality > fMinCentrality) && (fcentrality < fMaxCentrality))){return 0;}
if(!(fesd->GetNumberOfCaloClusters() > 2)){return 0;}
if(!(fevScaler > 0.7)){return 0;}

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

//if(CorruptFlag){return 0;}  //Go ahead and disable corruption cut.  It should be fixed by now?
//End Corruption Cut

FillHist("Event_Vertexx", fvertexx[2]);
FillHist("Event_Cent",fcentrality);
FillHist("Event_EP", fevPlane);
FillHist("Event_EPSP", fevScaler); 

return 1;

}

//___________________________Extract useful information from the clusters______________//
void AliAnalysisTaskThermalGAFlow::ScanClusters(){

//printf("Scanning photons...\n");

Int_t nPHOSClusters = 0;
Int_t nGoodClusters = 0;
Int_t x[4];
TLorentzVector p;

Int_t nClus0 = 0;
Int_t nClus1 = 0;
Int_t nClus2 = 0;

if(fesd){
FillHist("Stat_Efficiency", 19);
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
    FillHist("Clus_BCDistance", cluster->GetDistanceToBadChannel());
  }

  FillHist("Event_N0ClusVsCent", fcentrality, nClus0);
  FillHist("Event_N1ClusVsCent", fcentrality, nClus1);
  FillHist("Event_N2ClusVsCent", fcentrality, nClus2);
  FillHist("Event_NClusVsCent", fcentrality, nClus0+nClus1+nClus2);

  if(ApplyClusterCuts(cluster)){
    nGoodClusters = nGoodClusters+1;
    if(x[0] == 0){FillHist("Cell_N0Clus", x[2], x[3]);}
    if(x[0] == 1){FillHist("Cell_N1Clus", x[2], x[3]);}
    if(x[0] == 2){FillHist("Cell_N2Clus", x[2], x[3]);}
    FillHist("Clus_Pt", sqrt((p.Px()*p.Px()+p.Py()*p.Py())));
    FillHist("Clus_TOF", cluster->GetTOF());

//    printf("Photon Passed.\n");

    SavePhoton(cluster);
  }

}}
}

//___________________________Apply Cluster Cuts______________________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyClusterCuts(AliESDCaloCluster *cluster){
//Basic QA Cuts - mainly targetting bad channcels.
if(!(cluster->IsPHOS())){return 0;}
  FillHist("Stat_Efficiency", 0);
if(!(cluster->GetNCells() >= fMinCells)){return 0;}
  FillHist("Stat_Efficiency", 1);
if(!(cluster->E() >= fMinE)){return 0;}
  FillHist("Stat_Efficiency", 2);
if(!(cluster->GetTOF() > -0.0000001 && cluster->GetTOF() < 0.0000001)){return 0;}
  FillHist("Stat_Efficiency", 3);
if(!(cluster->GetDistanceToBadChannel() > 2.2)){return 0;}
  FillHist("Stat_Efficiency", 4);
//CPV
if(!(ApplyCPCuts(cluster))){return 0;}
  FillHist("Stat_Efficiency", 5);
//Shape
if(!(ApplyCoreShapeCut(cluster))){return 0;}
  FillHist("Stat_Efficiency", 6);
if(!(ApplyDispMatrixCut(cluster))){return 0;}
  FillHist("Stat_Efficiency", 7);
return 1;
}

//___________________________Apply Charged Particle Veto____________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyCPCuts(AliESDCaloCluster *cluster){
//Double_t dx = cluster->GetTrackDx();
//Double_t dz = cluster->GetTrackDz();
//  Double_t trackDr = sqrt(dx*dx + dz*dz);

Double_t trackDr = cluster->GetEmcCpvDistance();

FillHist("Clus_Trackdr", trackDr);
//if(!(trackDr > fMinTrackDr)){return 0;}

//TRY TO ADD PHOS_PbPb METHOD.
  //Parameterization of LHC10h period -- June 10 2015 PHOS meeting Yuri Kharlov said is same as 11h
  //_true if neutral_

TArrayI * itracks = cluster->GetTracksMatched() ;

if(!(itracks->GetSize() > 0)){return 1;} 
Int_t iTr = itracks->At(0);
if(!(iTr >= 0 && iTr < fEvent->GetNumberOfTracks())){return 0;}
//if(!(iTr >= 0 && iTr < fEvent->GetNumberOfTracks())){printf("iTr: %i, ntracks: %i\n", iTr, fEvent->GetNumberOfTracks()); return 0;}

Double_t dx = cluster->GetTrackDx();
Double_t dz = cluster->GetTrackDz();
AliVParticle* track = fEvent->GetTrack(iTr);
Double_t pt = track->Pt();
Int_t charge = track->Charge();

//What comes next is a byzentine series of computation that I don't understand.  The PWGGA seems to think that this extracts the std deviation of the CPV.  ...We'll see.
Double_t meanX=0;
Double_t meanZ=0;
Double_t mf = 0;
Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60);
mf = fesd->GetMagneticField(); //Positive for ++ and negative for --

if(mf<0){ //field --
  meanZ = -0.468318 ;
  if(charge>0){meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01));}
  else{meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01));}
}
else{ //Field ++
  meanZ= -0.468318;
  if(charge>0){meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01));}
  else{meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01));}
}

Double_t rz=(dz-meanZ)/sz ;
Double_t rx=(dx-meanX)/sx ;
Double_t r = TMath::Sqrt(rx*rx+rz*rz); //I understand this to be the number of std of projected track trajectory is from the center of the cluster.

FillHist("Clus_TrackMatchStd", r);
//printf("r = %f", r);
if(r < fMinCPVStd){return 0;}
return 1;
}

//___________________________Apply Core Shape Cut______________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyCoreShapeCut(AliESDCaloCluster *cluster){
//Core Energy Ratio Cut will to check for multiple maxima.
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

if(!(coreE / cluster->E() > fMinCoreEnergyRatio)){return 0;}
return 1;
}

//__________________________Apply Dispersion Matrix Cut___________________________//
Bool_t AliAnalysisTaskThermalGAFlow::ApplyDispMatrixCut(AliESDCaloCluster *cluster){
//Dispersion matrix eigenvalue will check ellipsidity (is that a word?)
Double_t dispLambda = cluster->GetM02();

FillHist("Clus_DispLambda", dispLambda);

if(!(dispLambda > fMinLambdaDisp)){return 0;}
return 1;
}

//___________________________Save Photon Candidate_________________________________//
void AliAnalysisTaskThermalGAFlow::SavePhoton(AliESDCaloCluster* cluster){

//printf("Saving A Photon...\n");

if(fSextant == 0x0){fSextant = new TObjArray(1); fSextant->SetOwner();}

Float_t cellGlobalPosArr[3];
TLorentzVector p;

cluster->GetPosition(cellGlobalPosArr);
cluster->GetMomentum(p, fvertexx);

AliCaloPhoton* photon = new AliCaloPhoton(p.Px(), p.Py(), p.Pz(), p.E());
photon->SetMomV2(&p);
photon->SetEMCx(cellGlobalPosArr[0]);
photon->SetEMCy(cellGlobalPosArr[1]);
photon->SetEMCz(cellGlobalPosArr[2]);

if(photon == 0x0){printf("The photon is null. That shouldn't happen."); return;}
fSextant->AddLast(photon);

//AliCaloPhoton *y = static_cast<AliCaloPhoton*> (fSextant->At(fSextant->GetLast()));
//printf("pPx: %f -- y: %f: -- yPx: %f\n", p.Px(), photon->GetMomV2()->Px() ,y->GetMomV2()->Px()  );

//printf("Photon Saved! \n");

/*
Need to save:
EVENT LEVEL - 
1) Centrality Information
2) Vertexx information
3) Event Plane Information
PHOTON LEVEL - Classic Array
4) Position information
5) Momentum information
6) From decay? (MC only)

They do it this way:  They make a hierarchy of arrays and lists going like:
TObjArray* fCaloPhotonsPHOSLists (multidimensional array of event types - each dimension is Cent, vert, and ep info - each element is a TList* arrayList - This level is hidden behind a method - GetCaloPhotonsPHOSList())

TList* arrayList (List of events - each element is an event)

TObjArray* mixPHOS (Array of photons - each element is a photon)

AliCaloPhoton* ph (Custom made object storing everything you could want to know (and then some) about the photon.  Documentation in $ALICE_PHYSICS/inst/include/AliCaloPhoton.h.  I think this much data is unnessiary - they have already cut on these things (and so have I).  They only use it to make post cut histograms for QA so far as I can tell.  Consider creating my own object or using this and simply not filling all the information or simply using my own standard (cpp) array of doubles.
*/
}

//____________________________Save the event for mixing_____________________________//
void AliAnalysisTaskThermalGAFlow::SaveEvent(){

//printf("Saving Event...\n");

fAtlas = ReferenceAtlas(fvertexx[2], fcentrality, fevPlane);

if(fSextant == 0x0){return;} //If there's no valid photons... don't bother saving.

//<Changed March 18, 2016>
if(fAtlas->At(100) == 0){fAtlas->AddLast(fSextant);} //if the Atlas is less than 100 entries long, work as normal
else{fAtlas->Delete(); fAtlas->AddLast(fSextant);} //if Atlas is more than 100 entries long, freak out and destory the continent, then start over
//</Changed March 18, 2016>

//printf("Event Saved! \n");
}

//____________________________Create Mixed Events___________________________________//
void AliAnalysisTaskThermalGAFlow::MesonExclusion(){
//Note:  Can only be done after Event and constituent photons have been saved... obviously.
//printf("Check\n");
if(fSextant == 0x0){return;}
//printf("Me!!\n");
//Use the mixed event method.

AliCaloPhoton *y1, *y2;
TList* atlas = ReferenceAtlas(fvertexx[2], fcentrality, fevPlane);

for(Int_t i1 = 0; i1 < fSextant->GetEntriesFast(); i1++){
 y1 = static_cast<AliCaloPhoton*> (fSextant->At(i1));
 for(Int_t i2 = i1+1; i2 < fSextant->GetEntriesFast(); i2++){
  y2 = static_cast<AliCaloPhoton*> (fSextant->At(i2));
//  printf("y1: %p, y2: %p, piarr[1]: %f, piarr[2]: %f \n", y1, y2, PionMass(y1, y2), PionPt(y1, y2));
  FillHist("Meson_RawM", PionMass(y1, y2), PionPt(y1, y2));
 }
 for(TObjLink* XLink = atlas->FirstLink(); XLink != atlas->LastLink(); XLink = XLink->Next()){
  //Don't match to the last link.  That's the Sextant.
  TObjArray *XEvent = static_cast<TObjArray*> (XLink->GetObject());
  for(Int_t ix2 = 0; ix2 < XEvent->GetEntriesFast(); ix2++){ 
   y2 = static_cast<AliCaloPhoton*> (XEvent->At(ix2));
   FillHist("Meson_XM", PionMass(y1, y2), PionPt(y1, y2));
  }
 }
}
}

//____________________________Get the Event List you need___________________________//
TList* AliAnalysisTaskThermalGAFlow::ReferenceAtlas(Double_t vertx, Double_t cent, Double_t evplane){
TList* atlas;
Int_t index;
Int_t vertxbins = fMixVertxbins;
Int_t centbins = fMixCentbins;
Int_t evplanebins = fMixEvbins;

Int_t vertx_index = std::floor( ((vertx / (2*fMaxVertexx)) + 0.5) * vertxbins );
Int_t cent_index = std::floor( (((cent - fMinCentrality) / (fMaxCentrality - fMinCentrality)) * centbins) );
Int_t evplane_index = std::floor( (evplane/TMath::Pi()) * evplanebins );

//printf("vx: %f\n", vertx);
//printf("vx_i: %i |  c_i: %i |  ev_i: %i \n", vertx_index, cent_index, evplane_index);

if(vertx_index == vertxbins){vertx_index = vertx_index - 1;}
if(cent_index == centbins){cent_index = cent_index - 1;}
if(evplane_index == evplanebins){evplane_index = evplane_index - 1;}

index = vertx_index + cent_index*vertxbins + evplane_index*vertxbins*centbins;

//printf("index: %i\n", index);

atlas = static_cast<TList*> (fCompass->At(index));

//printf("Compass at: %p | cast to: %p | atlas: %p \n", fCompass, fCompass->At(index), atlas);

if(atlas){return atlas;}
else {
if(fDebug > 1){AliError("The Altas Reference failed. Creating new Atlas....");}
atlas = new TList();
atlas->SetOwner();
fCompass->AddAt(atlas, index);
return atlas;}
}

//____________________________Compute Invariant Mass________________________________//
Double_t AliAnalysisTaskThermalGAFlow::PionMass(AliCaloPhoton* y1, AliCaloPhoton* y2){
const TLorentzVector *p1 = y1->GetMomV2();
const TLorentzVector *p2 = y2->GetMomV2();
TLorentzVector ppi = *p1 + *p2;

Double_t Mpi = sqrt(ppi.E()*ppi.E()-ppi.Px()*ppi.Px()-ppi.Py()*ppi.Py()-ppi.Pz()*ppi.Pz());
return Mpi;
}

//____________________________Compute Pion Pt________________________________//
Double_t AliAnalysisTaskThermalGAFlow::PionPt(AliCaloPhoton* y1, AliCaloPhoton* y2){
const TLorentzVector *p1 = y1->GetMomV2();
const TLorentzVector *p2 = y2->GetMomV2();
TLorentzVector ppi = *p1 + *p2;

Double_t Ptshared = sqrt((ppi.Px()*ppi.Px()+ppi.Py()*ppi.Py()));
return Ptshared;
}

//____________________________Compute Flow___________________________________//
void AliAnalysisTaskThermalGAFlow::IncFlow(){
AliCaloPhoton *y;
Float_t phi, Dphi, DphiA, DphiC, pt;

if(fSextant == 0x0){return;}

for(Int_t i = 0; i < fSextant->GetEntriesFast(); i++){
 y = static_cast<AliCaloPhoton*> (fSextant->At(i));
 phi = y->Phi(); //This *should* be the phi in the global coordinate system computed by taking the arctan of the position vector's y (vertical height) and x (left-right) in the ALICE coordinate system.  However, it might also be y and x in momentum space *or* x and z in PHOS coordinate system (that's what Dmitri seems to use - but surely that's not right?  Whatever - we'll see what happens when we do it.  If it returns 0.0, then it implies x or z is zero.)

 pt = sqrt((y->Px()*y->Px()+y->Py()*y->Py())); 

 Dphi = phi - fevPlane;
 DphiA = phi - fevPlaneA;
 DphiC = phi - fevPlaneC;

while(Dphi < 0){Dphi = Dphi + TMath::Pi();}
while(DphiA < 0){DphiA = DphiA + TMath::Pi();}
while(DphiC < 0){DphiC = DphiC + TMath::Pi();}
while(Dphi > TMath::Pi()){Dphi = Dphi - TMath::Pi();}
while(DphiA > TMath::Pi()){DphiA = DphiA - TMath::Pi();}
while(DphiC > TMath::Pi()){DphiC = DphiC - TMath::Pi();}

// FillHist("InclGamma_PhiVsPt", Dphi, sqrt((y->Px()*y->Px()+y->Py()*y->Py()))); //Remember to Fourier expand in each pt bin in post analysis! 
//printf("Dphi: %f, Pt: %f\n", -1*phi, sqrt((y->Px()*y->Px()+y->Py()*y->Py())));

 FillHist("InclGamma_PhiVsPt", phi, pt); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsPt_AC", Dphi, pt); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsCent_AC", Dphi, fcentrality); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsPt_A", DphiA, pt); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsCent_A", DphiA, fcentrality); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsPt_C", DphiC, pt); //Remember to Fourier expand in each pt bin in post analysis! 
 FillHist("InclGamma_DPhiVsCent_C", DphiC, fcentrality); //Remember to Fourier expand in each pt bin in post analysis! 

//Also, might need to correct for EP resolution - recommended by Dmitri in "Elliptic and triangular flow of inclusive and direct photons...".  Otherwise... that *should* do it for InclGamma phi. ...It's really pretty straitforward.  OR at least it should be.
}

}

//___________________________Compute Decay Flow_______________________________//


//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//                                                                                  //
//                                Helper Functions                                  //
//                                                                                  //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

//___________________________Hypnotic Protocol_____________________________________//
void AliAnalysisTaskThermalGAFlow::Hypnotic(){

TLorentzVector p1, p2, p;
AliESDCaloCluster *cluster1, *cluster2;

for(Int_t i1cluster = 0; i1cluster < fesd->GetNumberOfCaloClusters(); i1cluster++) {
 cluster1 = fesd->GetCaloCluster(i1cluster);
 if(WeakCuts(cluster1)){
  for(Int_t i2cluster = i1cluster+1; i2cluster < fesd->GetNumberOfCaloClusters(); i2cluster++)  {
   cluster2 = fesd->GetCaloCluster(i2cluster);
   if(WeakCuts(cluster2)){
    cluster1->GetMomentum(p1, fvertexx);
    cluster2->GetMomentum(p2, fvertexx);
    p=p1+p2;
    FillHist("Cluster_M", sqrt(p.E()*p.E()-p.Px()*p.Px()-p.Py()*p.Py()-p.Pz()*p.Pz()), sqrt((p.Px()*p.Px()+p.Py()*p.Py()))); 
}}}}

}

//_________________________WeakCuts______________________________________//
Bool_t AliAnalysisTaskThermalGAFlow::WeakCuts(AliESDCaloCluster *cluster){
if(!(cluster->IsPHOS())){return 0;}
if(!(cluster->GetNCells() >= fMinCells)){return 0;}
if(!(cluster->E() >= fMinE)){return 0;}
return 1;
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

//____________________________Keep a record of events so the user knows what's going on_____//
void AliAnalysisTaskThermalGAFlow::CaptainsLog(Int_t s){
if(fDebug > 0){
if(s == 0){AliInfo(Form("Beginning New Event. Standby."));}
if(s == 3){AliInfo(Form("...Finished Event.\n"));}
if(s != 0 && s != 3){AliInfo(Form("Step %i Completed.", s));}
}

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
