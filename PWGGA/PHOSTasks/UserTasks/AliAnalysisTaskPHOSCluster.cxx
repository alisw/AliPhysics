//-------------------------------------------------------------------------
////     AnalyisTask for PHOS clusters
////     Runs on ESDs and AODs
////     Authors: Alexej Kraiker (partly based on analysis task PHOSNeutralMesons by Malte Hecker and Fabian Pliquett)
////     Date: 11/07/2017
////-------------------------------------------------------------------------

#include "AliAnalysisTaskPHOSCluster.h"

#include <vector>

#include "TRefArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TChain.h"
#include "TParticle.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"

//AliRoot include files
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
//#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliOADBContainer.h"
#include "AliPHOSCalibData.h"
#include "AliVCluster.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliPHOSGeometry.h"
//#include "AliEMCALGeometry.h"




using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPHOSCluster);


AliAnalysisTaskPHOSCluster::AliAnalysisTaskPHOSCluster()
  :AliAnalysisTaskSE(),
  fOutput(0),
  fAnyEv(0),
  fEventCounter(0),
  fClusterTree(0),
  fClusterCellsInformation(0),
  fZVertex(0),
  fMinCells(0),
  fPHOSGeo(0),
  fPHOSCalibData(0),
  fUtils(0),
  fcaloClusters(0x0),
  fCluster(0x0),
  fClusterCells(0x0),
  fPHOSHitmap(0x0),
  fClusterEnergyVsNoC(0x0),
  fClusterEnergyVsM02(0x0),
  fClusterEnergyVsM20(0x0),
  fEventCuts(0x0),
  fClusterCuts(0x0),
  fFillHitmapCellByCell(0)
{
// Dummy constructor always needed for I/O
}


AliAnalysisTaskPHOSCluster::AliAnalysisTaskPHOSCluster(const char* name)
  :AliAnalysisTaskSE(name),
  fOutput(0),
  fAnyEv(0),
  fEventCounter(0),
  fClusterTree(0),
  fClusterCellsInformation(0),
  fZVertex(0),
  fMinCells(0),
  fPHOSGeo(0),
  fPHOSCalibData(0),
  fUtils(0),
  fcaloClusters(0x0),
  fCluster(0x0),
  fClusterCells(0x0),
  fPHOSHitmap(0x0),
  fClusterEnergyVsNoC(0x0),
  fClusterEnergyVsM02(0x0),
  fClusterEnergyVsM20(0x0),
  fEventCuts(0x0),
  fClusterCuts(0x0),
  fFillHitmapCellByCell(0)
{
  DefineOutput(1, TList::Class());
}


//_______________________________________ Destructor _______________________________________
// Destructor. Clean-up the output list, but not the histograms that are put inside
// (the list is owner and will clean-up these histograms). Protect in PROOF case.
AliAnalysisTaskPHOSCluster::~AliAnalysisTaskPHOSCluster() {

  if (fClusterTree && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){ //delete tree before deleting the fOutput!
    delete fClusterTree;
    fClusterTree = 0x0;
  }
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput;
    fOutput = 0x0;
  }
}

//_______________________________________ UserCreateOutputObjects _______________________________________
// Create histogramms and output
// Called once (on the worker node)
void AliAnalysisTaskPHOSCluster::UserCreateOutputObjects() {
  cout << "AliAnalysisTaskPHOSCluster - Input settings:" << endl;
  cout << Form("Z Vertex cut: %.3f", fZVertex) << endl;
  cout << Form("MinCells cut: %d", fMinCells) << endl;


  if (!fOutput) fOutput = new TList();
  fOutput->SetOwner();


  fClusterCellsInformation = new AliCaloClusterContent();

  fClusterTree   = new TTree("fClusterTree", "AliPHOSCluster");
  fClusterTree ->Branch("fClusterCellsInformation", &fClusterCellsInformation);
  fClusterTree ->SetAutoSave(100000);


  fcaloClusters = new TRefArray();



// Creation of the histograms
  fClusterCuts = new TH1F("fClusterCuts", "N o Cluster after cuts", 3, 0., 3.);
  fClusterCuts ->GetXaxis()->SetBinLabel(1, "All clusters");
  fClusterCuts ->GetXaxis()->SetBinLabel(2, "PHOS Cluster");
  fClusterCuts ->GetXaxis()->SetBinLabel(3, "Min. Cells");
  fClusterCuts ->GetYaxis()->SetTitle("# Clusters");

  fEventCuts   = new TH1F("fEventCuts", "N o Events after cuts", 4, 0., 4.);
  fEventCuts   ->GetXaxis()->SetBinLabel(1, "All Events");
  fEventCuts   ->GetXaxis()->SetBinLabel(2, "Primary Vertex");
  fEventCuts   ->GetXaxis()->SetBinLabel(3, "PileUp");
  fEventCuts   ->GetXaxis()->SetBinLabel(4, "ZVertex");
  fEventCuts   ->GetYaxis()->SetTitle("# Events");

  fClusterEnergyVsNoC = new TH2F("fClusterEnergyVsNoC", "Cluster energy vs number of cells per cluster", 500, 0., 5., 100, 0., 100.);
  fClusterEnergyVsNoC ->GetXaxis()->SetTitle("#it{E}_{Cluster} (GeV)");
  fClusterEnergyVsNoC ->GetYaxis()->SetTitle("# Cells");

  fClusterEnergyVsM02 = new TH2F("fClusterEnergyVsM02", "Cluster energy vs cluster M02", 500, 0., 5., 400, 0., 60.);
  fClusterEnergyVsM02 ->GetXaxis()->SetTitle("#it{E}_{Cluster} (GeV)");
  fClusterEnergyVsM02 ->GetYaxis()->SetTitle("M02");

  fClusterEnergyVsM20 = new TH2F("fClusterEnergyVsM20", "Cluster energy vs cluster M20", 500, 0., 5., 400, 0., 20);
  fClusterEnergyVsM20 ->GetXaxis()->SetTitle("#it{E}_{Cluster} (GeV)");
  fClusterEnergyVsM20 ->GetYaxis()->SetTitle("M20");

  fPHOSHitmap         = new TH2F("fPHOSHitmap", "Hitmap of PHOS", 4*64, 0., 4*64., 56, 0., 56.);
  fPHOSHitmap         ->GetXaxis()->SetTitle("#it{#phi} (rad)");
  fPHOSHitmap         ->GetYaxis()->SetTitle("#it{#eta}");




  fOutput ->Add(fClusterTree);

  fOutput ->Add(fClusterCuts);
  fOutput ->Add(fEventCuts);
  fOutput ->Add(fClusterEnergyVsNoC);
  fOutput ->Add(fClusterEnergyVsM02);
  fOutput ->Add(fClusterEnergyVsM20);
  fOutput ->Add(fPHOSHitmap);


  PostData(1,fOutput);
}

//_______________________________________ UserExec _______________________________________

void AliAnalysisTaskPHOSCluster::UserExec(Option_t *){

  // ========== START Load Event ========== //
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
  if (!aodH && !esdH)  AliError("Could not get ESD or AODInputHandler");
  if (esdH){
    fAnyEv = dynamic_cast<AliESDEvent*> (esdH->GetEvent());
  }
  else if (aodH){
    fAnyEv = dynamic_cast<AliAODEvent*> (aodH->GetEvent());
  }
  else{
    AliFatal("Neither ESD nor AOD event found");
    return;
  }
  fEventCuts->Fill(0.5);  // count all events
  // ========== END Load Event ========== //


  // ========== START Remove events with no vertex ========== //
  if (esdH) {
    if (!(dynamic_cast<AliESDEvent*>(fAnyEv)->GetPrimaryVertex()->GetStatus())){
      return;
    }
  }
  else if (aodH) {
    //use NContributors instead of GetStatus to check if there is a primary vertex
    if (!((dynamic_cast<AliAODEvent*>(fAnyEv)->GetPrimaryVertex()->GetNContributors())>0)) {
      return;
    }
  }
  fEventCuts ->Fill(1.5); // count all events with vertex
  // ========== END Remove events with no vertex ========== //


  // Remove PileUp Events
  if(fAnyEv->IsPileupFromSPD(3,.8,3.,2.,5.)) return;
  fEventCuts ->Fill(2.5);

//  if (esdH) dynamic_cast<AliESDEvent*>(fAnyEv)->GetVertex()->GetXYZ(vertex);
//  else if (aodH) dynamic_cast<AliAODEvent*>(fAnyEv)->GetVertex(0)->GetXYZ(vertex); // function GetVertex needs argument for AODs


  // ZVertex Cut
  if(fabs(fAnyEv->GetPrimaryVertex()->GetZ()) > fZVertex) return;  // remove events
  fEventCuts ->Fill(3.5); // count events after ZVertex cut


  // ========== START Setting PHOS matrix ========== //
  Int_t runNumber = 0;
  runNumber = fAnyEv->GetRunNumber();
  if (fPHOSGeo==0) {

    fPHOSGeo = AliPHOSGeometry::GetInstance() ;

    if(!fPHOSGeo){ //Geometry not yet constructed with Tender
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP","");

      AliOADBContainer geomContainer("phosGeo");
      geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
      TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
      for(Int_t mod=0; mod<5; mod++) {
        if(!matrixes->At(mod)) continue;
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        printf("....TASK.....Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
        ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
      }
    }
  }

  if(fEventCounter == 0) { // Only done for the first Event
    cout << "Estimating reconstruction pass" <<  endl;
    Int_t recoPass = -1;
    TTree * t = am->GetTree();
    if(t){
      TFile * f = t->GetCurrentFile() ;
      if(f){
        TString fname(f->GetName());
        if(fname.Contains("pass1")) recoPass=1;
        else
          if(fname.Contains("pass2")) recoPass=2;
          else
            if(fname.Contains("pass3")) recoPass=3;
            else
              if(fname.Contains("pass4")) recoPass=4;
		  cout << "Reconstruction pass is: " << recoPass << endl;
        cout << "Run number is: " << runNumber << endl;
      }
    }
    if(recoPass<0){
      AliError("Can not find pass number from file name, is set to -1");
    }
    //Load recalibration data
    AliOADBContainer calibContainer("phosRecalibration");
    calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
    TObjArray *recalib = (TObjArray*)calibContainer.GetObject(runNumber,"PHOSRecalibration");
    if(!recalib){
      AliWarning(Form("Can not read calibrations for run %d.\n",runNumber)) ;
    }
    else{
      cout << "Doing PHOS calibration" << endl;
      fPHOSCalibData = (AliPHOSCalibData*)recalib->At(recoPass-1) ;
      if(!fPHOSCalibData) {
        AliWarning(Form("Can not find calibration for run %d, pass %d \n",runNumber, recoPass)) ;
      }
    }
  } //if(fEventCounter == 0)

  // ========== END Setting PHOS matrix for ESDs ========== //



  // ========== START Loading and saving clusters ========== //

  Int_t nCluster = 0;
  nCluster = fAnyEv ->GetNumberOfCaloClusters();



  for(Int_t iCurrentCluster=0; iCurrentCluster<nCluster; iCurrentCluster++){
    if(esdH)      fCluster    = dynamic_cast<AliESDCaloCluster*>(fAnyEv->GetCaloCluster(iCurrentCluster)); // pointer to cluster
    else if(aodH) fCluster    = dynamic_cast<AliAODCaloCluster*>(fAnyEv->GetCaloCluster(iCurrentCluster));
    if(!fCluster){
      AliError(Form("Cluster %d is neither ESD or AOD Cluster!", iCurrentCluster));
      continue;
    }

    fClusterCuts ->Fill(0.5); // count ALL clusters (EmCal AND PHOS)
    if(fCluster->IsPHOS()){ // analyse only PHOS cluster
      fClusterCuts->Fill(1.5);

      if(fMinCells > fCluster->GetNCells()) continue; // remove clusters with less than 2 cells
      fClusterCuts->Fill(2.5);


		Float_t pos[3] = {0,0,0};  // get cluster coordinates
		fCluster->GetPosition(pos);
		TVector3 vpos(pos);

		Int_t    modNrClusterPos, relId[4], cellX, cellZ;
		fPHOSGeo->GlobalPos2RelId(vpos,relId);  // calculate relative id
		modNrClusterPos  = relId[0];
		cellX = relId[2];
		cellZ = relId[3];


//		if(fRecalibrateModuleWise) { // apply recalibration of cluster energy
//		  Double_t recalib[3] = {fRecalFactorMod1, fRecalFactorMod2, fRecalFactorMod3 };  //recalFactors are set in AddTask!
//		  fCluster->SetE(fCluster->E()*recalib[modNrClusterPos-1]);
//		}

      UShort_t* CellsID = fCluster ->GetCellsAbsId();  // get all cells of cluster
      UShort_t  NCells  = fCluster ->GetNCells();
		fClusterCells = fAnyEv->GetPHOSCells();   // get PHOS cells
      Int_t cellRelID[4];


      if(fFillHitmapCellByCell){
        for(Int_t icell=0; icell<NCells; icell++){  // loop over all cluster cells
          fPHOSGeo ->AbsToRelNumbering(CellsID[icell], cellRelID);
          fPHOSHitmap ->Fill((cellRelID[0]-1)*64+cellRelID[2], cellRelID[3], fClusterCells->GetCellAmplitude(CellsID[icell]));  // fill histogramm cell by cell
        }
	   }
		else fPHOSHitmap ->Fill((modNrClusterPos-1)*64 + cellX, cellZ, fCluster->E());  // fill histogramm cluster by cluster

		// Fill histogramms
		if(fCluster->GetM02()    >= 0) fClusterEnergyVsM02 ->Fill(fCluster->E(), fCluster->GetM02());
		if(fCluster->GetM20()    >= 0) fClusterEnergyVsM20 ->Fill(fCluster->E(), fCluster->GetM20());
		if(fCluster->GetNCells() >= 0) fClusterEnergyVsNoC ->Fill(fCluster->E(), fCluster->GetNCells());

		fClusterCellsInformation->SetClusterAndCells(fCluster, fClusterCells, fPHOSGeo);  // save information about cluster and its cells
		fClusterTree ->Fill();

//    if(fClusterCellsInformation->IsFilled()) fClusterCellsInformation->Reset();  // Reset cluster object
      fClusterCellsInformation->Reset();  // Reset cluster object

    } // if(fCluster->IsPHOS())
  } // cluster loop
  // ========== END Loading and saving clusters ========== //

  fEventCounter++;
  PostData(1,fOutput);
}




//TODO histo mit events wie viele mit und ohne cluster
