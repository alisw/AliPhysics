//--- For C++ ----
#include <TApplication.h>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <THashList.h>
#include <TDirectoryFile.h>
#include <TStyle.h>
#include <TEllipse.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <cstdlib>
#include "TStopwatch.h"
#include "TLegend.h"
#include "TLatex.h"

using namespace std;

#include "AliExternalTrackParam.h"
#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliESDcascade.h"
#include "AliESDVertex.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
//#include "AliLog.h"

Int_t lNotStationary = 0;
Long_t lNegX = 0;
Long_t lPosX = 0;

#include "TMath.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TTreeStream.h"

//stuff for mat corr
#include <TChain.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "TGeoManager.h"
#include <TRegexp.h>

#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliLog.h"
#include "AliTrackerBase.h"

void sandbox01() {
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //
    //  Sandbox macro
    //  =============
    //
    //      This code generates candidates from
    //      two-track combinations using a TTree in which
    //      full AliExternalTrackParam/ESDtrack information is stored
    //      for all V0 daughter tracks
    //
    //      Allows for faster methodology studies than runnning
    //      on the grid, as reconstruction can be replayed from
    //      'scratch' (in the offline weak decay reco scheme)
    //
    //      IMPORTANT: This uses an output from the task
    //
    //          AliAnalysisTaskStrangenessVsMultiplicityMCRun2
    //
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    gStyle->SetOptStat(0);
    
    //Open file, get stuff, please
    TFile *file = new TFile("AnalysisResultsFindable.root", "READ") ;
    TTree *fTreeV0 = (TTree*) file->FindObjectAny("fTreeV0");
    
    TGaxis::SetMaxDigits(2);
    
    Float_t fPVx, fPVy, fPVz, fMagF, fCentrality;
    
    AliExternalTrackParam *fpTrack = new AliExternalTrackParam();
    AliExternalTrackParam *fnTrack = new AliExternalTrackParam();
    
    //Variable Definition=============================================
    //Kinematic
    Float_t lPt, lRap, lPtMC, lNegEta, lPosEta;
    //DCA Variables
    Float_t lDcaV0Daughters;
    Float_t lDcaPosToPrimVertex,  lDcaNegToPrimVertex;
    //Cosine of Pointing Angle variable
    Float_t lV0CosinePointingAngle;
    Float_t lDCAV0ToPrimVertex;
    //Decay Radius and distance over total momentum
    Float_t lV0Radius, lDistOverTotMom;
    //Least Number of TPC Clusters
    Int_t lLeastNbrCrossedRows;
    Float_t lLeastNbrCrossedRowsOverFindable;
    //TPC dE/dx acquired with AliPIDResponse class
    Float_t lNSigmasPosProton,lNSigmasNegProton,lNSigmasPosPion,lNSigmasNegPion;
    Float_t lArmPt,lArmAlpha;
    //Multiplicity Variable
    Double_t lMultiplicity = 0; //for ease of handling later
    Int_t lPID = 3122, lPIDPositive, lPIDNegative;
    Int_t fTreeVariableLeastCrossedRows;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // Reconstructed vars
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // Sandbox vars
    fTreeV0->SetBranchAddress("fTreeVariablePosTrack", &fpTrack);
    fTreeV0->SetBranchAddress("fTreeVariableNegTrack", &fnTrack);
    fTreeV0->SetBranchAddress("fTreeVariablePrimVertexX", &fPVx);
    fTreeV0->SetBranchAddress("fTreeVariablePrimVertexY", &fPVy);
    fTreeV0->SetBranchAddress("fTreeVariablePrimVertexZ", &fPVz);
    fTreeV0->SetBranchAddress("fTreeVariableMagneticField", &fMagF);
    fTreeV0->SetBranchAddress("fTreeVariableCentrality", &fCentrality);
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    fTreeV0->SetBranchAddress("fTreeVariablePtMC",&lPtMC);
    fTreeV0->SetBranchAddress("fTreeVariablePID",&lPID);
    fTreeV0->SetBranchAddress("fTreeVariablePIDPositive",&lPIDPositive);
    fTreeV0->SetBranchAddress("fTreeVariablePIDNegative",&lPIDNegative);
    fTreeV0->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastCrossedRows);
    //--- TPC dEdx Variables ------------------------------------------
    fTreeV0->SetBranchAddress("fTreeVariableNSigmasPosProton",&lNSigmasPosProton);
    fTreeV0->SetBranchAddress("fTreeVariableNSigmasNegProton",&lNSigmasNegProton);
    fTreeV0->SetBranchAddress("fTreeVariableNSigmasPosPion",&lNSigmasPosPion);
    fTreeV0->SetBranchAddress("fTreeVariableNSigmasNegPion",&lNSigmasNegPion);
    //--- Topological selection variables -----------------------------
    fTreeV0->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
    fTreeV0->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    Int_t lRunNumber;
    fTreeV0->SetBranchAddress("fTreeVariableRun",&lRunNumber);
    
    //Base processing
    cout<<"Tree detected with "<<fTreeV0->GetEntries() << " entries."<<endl;
    
    TFile *fileout = new TFile("output.root", "RECREATE");
    
    //Tree for debug purpose: two simultaneous recos, please
    TTree *fTreeDebug = new TTree("fTreeDebug", "");
    
    //The topo basics
    Float_t vHypertritonMass, vCosPA, vDCAV0Dau, vDCAPosToPV, vDCANegToPV, vV0Radius;
    
    //The typical vars
    Float_t vPt, vPtMC;
    Int_t vPID; 
    
    //Attaching to tree
    fTreeDebug->Branch("vHypertritonMass",  &vHypertritonMass,  "vHypertritonMass/F");
    fTreeDebug->Branch("vCosPA",            &vCosPA,            "vCosPA/F"          );
    fTreeDebug->Branch("vDCAV0Dau",         &vDCAV0Dau,         "vDCAV0Dau/F"       );
    fTreeDebug->Branch("vDCAPosToPV",       &vDCAPosToPV,       "vDCAPosToPV/F"     );
    fTreeDebug->Branch("vDCANegToPV",       &vDCANegToPV,       "vDCANegToPV/F"     );
    fTreeDebug->Branch("vV0Radius",         &vV0Radius,         "vV0Radius/F"     );
    
    fTreeDebug->Branch("vPID",              &vPID,              "vPID/I"            );
    fTreeDebug->Branch("vPt",               &vPt,               "vPt/F"             );
    fTreeDebug->Branch("vPtMC",             &vPtMC,             "vPtMC/F"           );

    Bool_t lPerfectSignal = kFALSE, lPerfectDaughters = kFALSE;
    fTreeDebug->Branch("lPerfectSignal",&lPerfectSignal,"lPerfectSignal/O");
    fTreeDebug->Branch("lPerfectDaughters",&lPerfectDaughters,"lPerfectDaughters/O");
    
    //Calls to GetDCAV0Dau will be done directly with improved weak decay finder
    AliAnalysisTaskWeakDecayVertexer *wdvtask = new AliAnalysisTaskWeakDecayVertexer("wdv");
    wdvtask -> SetUseImprovedFinding();

    //Hypertriton mass calculation variables
    TLorentzVector Hypertriton;
    TLorentzVector posHe3, negPi; //Lorentz vector of helium-3 and pion in the LAB
    
    //Performance metrics, counter
    Long_t ipassed = 0;
    TStopwatch* timer = new TStopwatch();
    timer->Start ( kTRUE );
    //_______________________________________________________________________
    for(Long_t i=0;i<fTreeV0->GetEntries(); i++) {
        fTreeV0->GetEntry(i);
        if ( i%10000==0 )
            cout<<"At entry "<<i<<" out of "<<fTreeV0->GetEntries()<<" ("<<100.*((Double_t)i)/fTreeV0->GetEntries()<<"%), "<<ipassed<<" candidates... "<<endl;
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Base selections
        if(
           lDcaNegToPrimVertex <= -1 ||
           lDcaPosToPrimVertex <= -1 ||
           fTreeVariableLeastCrossedRows <= 80
           ) continue;
        ipassed++;
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Determine if true hypertriton signal or not
        lPerfectSignal = kFALSE;
        lPerfectDaughters = kFALSE;
        if ( lPID == 1010010030 && lPIDPositive == 1000020030 && lPIDNegative == -211 ) lPerfectSignal = kTRUE;
        if ( lPIDPositive == 1000020030 && lPIDNegative == -211 ) lPerfectDaughters = kTRUE;
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Store vars that are already there
        vDCAPosToPV = lDcaPosToPrimVertex; vDCANegToPV = lDcaNegToPrimVertex;
        vPID = lPID; vPtMC = lPtMC;
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Replay V0 finding sequence
        AliExternalTrackParam nt(*fnTrack), pt(*fpTrack), *ntp=&nt, *ptp=&pt;
        Double_t xp, xn;
        vDCAV0Dau = wdvtask->GetDCAV0Dau(ptp, ntp, xp, xn, fMagF);
        ntp->PropagateTo(xn,fMagF);
        ptp->PropagateTo(xp,fMagF);
        
        //use these tracks to create a V0
        AliESDv0 v0vertex(nt,0,pt,1);
        v0vertex.SetDcaV0Daughters(vDCAV0Dau);
        
        //Do V0 refit
        v0vertex.Refit();

        //Get V0 vertex variables
        vCosPA  = v0vertex.GetV0CosineOfPointingAngle(fPVx,fPVy,fPVz);
        Double_t lVtx[3]; v0vertex.GetXYZ(lVtx[0],lVtx[1],lVtx[2]);
        vV0Radius = TMath::Sqrt(lVtx[0]*lVtx[0]+lVtx[1]*lVtx[1]);
        vPt = v0vertex.Pt();
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Calculate Hypertriton mass
        Hypertriton.Clear();
        posHe3.Clear();
        negPi.Clear();
        posHe3.SetXYZM(2*ptp->Px(),2*ptp->Py(),2*ptp->Pz(),AliPID::ParticleMass(AliPID::kHe3));
        negPi.SetXYZM(ntp->Px(),ntp->Py(),ntp->Pz(),AliPID::ParticleMass(AliPID::kPion));
        Hypertriton=posHe3+negPi;
        vHypertritonMass = Hypertriton.M();
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

        //Everything needs to be set properly here!
        fTreeDebug->Fill();
    }
    fileout->Write();
}
