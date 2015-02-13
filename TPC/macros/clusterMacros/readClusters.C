/// \file readClusters.C
/// \brief Macro to read TPC clusters from TPC.RecPoints.root and AliESDfriends.root and fill THnSparse with the content
///
/// Used for the HLT-TPC cluster verification
///
/// Usage:
/// ~~~
/// aliroot -b -l -q readClusters.C'("<Folder w/ sim/rec output>", "<Simulation Id>", "<Simulation Version>",
///                                  <minimal folder Id>, <maximal folder Id>, <use ESDfriends 0 or 1>)'
/// ~~~
///
/// Example:
/// ~~~
/// aliroot -b -l -q readClusters.C'("/lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia","Pythia","20a",30000,30001,1)'
/// ~~~
///
/// Will read :
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30000/offline20a/AliESDfriends.root
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30000/HLThw20a/AliESDfriends.root
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30000/HLThwRedux20a/AliESDfriends.root
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30001/offline20a/AliESDfriends.root
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30001/HLThw20a/AliESDfriends.root
/// * /lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia/30001/HLThwRedux20a/AliESDfriends.root
///
/// Will write :
/// * $CWD/results/results_friends_Pythia_20a.root
///
/// \author Jochen Thaeder <jochen@thaeder.de>

/*
   #include "TROOT.h"
   #include "TStyle.h"
   #include "TFile.h"
   #include "TCanvas.h"
   #include "THnSparse.h"
   #include "TColor.h"
   #include "TTree.h"
   #include "TAxis.h"
   #include "TDirectoryFile.h"

   #include "TDirectoryFile.h"
   #include "AliESDfriend.h"
   #include "AliGeomManager.h"
   #include "AliCDBManager.h"
   #include "AliTPCClustersRow.h"
   #include "AliTPCclusterMI.h"
   #include "AliTPCseed.h"
*/

void FillFriends(AliESDfriend* esdFriend, THnSparseF *sparse);
void FillRecPoints(TTree* tree, THnSparseF *sparse);
void FillClusters(THnSparseF *sparse,  AliTPCclusterMI* cl);
void SetupStyle();
THnSparseF* CreateTHnSparse(Char_t* name);

// ==================================================================================
void readClusters( Char_t *folder = "/lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06_Pythia",
		   Char_t *type = "Pythia", Char_t *version = "8",
		   Int_t minId = 30000, Int_t maxId = 30020,
		   Bool_t useFriends=kFALSE) {

  /// --------------------------------------------------------
  /// -- Setup
  /// --------------------------------------------------------

  // -- Setup style
  SetupStyle();

  // -- Setup OCDB and Geometry
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  printf("AliCDBconnect: #### Loading geometry...\n");
  AliGeomManager::LoadGeometry();
  if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD") ) {
    printf("Problem with align objects");
  }

  // --------------------------------------------------------
  // -- Create THnSparse
  // --------------------------------------------------------

  THnSparseF* spo    = CreateTHnSparse("spo");
  THnSparseF* sphhw  = CreateTHnSparse("sphhw");
  THnSparseF* sphhwR = CreateTHnSparse("sphhwR");

  // --------------------------------------------------------
  // -- Open input files / loop over them / fill clusters
  // --------------------------------------------------------

  TFile* fo    = NULL;
  TFile* fhhw  = NULL;
  TFile* fhhwR = NULL;

  const Char_t *fileName[] = {"TPC.RecPoints.root", "AliESDfriends.root"};

  // -- Loop over the folder -> files
  for (Int_t fileIdx = minId; fileIdx <= maxId; ++fileIdx) {

    printf(" -- Open new set of files : %d -- \n", fileIdx);

    fo    = TFile::Open(Form("%s/%05d/offline%s/%s",    folder, fileIdx, version, fileName[useFriends]));
    fhhw  = TFile::Open(Form("%s/%05d/HLThw%s/%s",      folder, fileIdx, version, fileName[useFriends]));
    fhhwR = TFile::Open(Form("%s/%05d/HLThwRedux%s/%s", folder, fileIdx, version, fileName[useFriends]));

    if (!fo)
      continue;

    // -------------------------------------------------------
    // -- Get ptr to friends
    // -------------------------------------------------------

    TTree* toF    = NULL;
    TTree* thhwF  = NULL;
    TTree* thhwRF = NULL;
    AliESDfriend* fFo    = NULL;
    AliESDfriend* fFhhw  = NULL;
    AliESDfriend* fFhhwR = NULL;

    if (useFriends) {
      toF    = dynamic_cast<TTree*> (fo->Get("esdFriendTree"));
      thhwF  = dynamic_cast<TTree*> (fhhw->Get("esdFriendTree"));
      thhwRF = dynamic_cast<TTree*> (fhhwR->Get("esdFriendTree"));

      fFo    = new AliESDfriend();
      fFhhw  = new AliESDfriend();
      fFhhwR = new AliESDfriend();

      toF->GetBranch("ESDfriend.")->SetAddress(&fFo);
      thhwF->GetBranch("ESDfriend.")->SetAddress(&fFhhw);
      thhwRF->GetBranch("ESDfriend.")->SetAddress(&fFhhwR);
    }

    // -------------------------------------------------------
    // -- Lop over the events in file
    // -------------------------------------------------------
    Int_t nEvents = (useFriends) ? toF->GetEntries() : fo->GetNkeys();

    for (Int_t event=0 ; event < nEvents; ++event) {
      if (useFriends) {
	if (fo) {
	  fo->cd();
	  toF->GetEntry(event);
	  FillFriends(fFo, spo);
	}
	if (fhhw) {
	  fhhw->cd();
	  thhwF->GetEntry(event);
	  FillFriends(fFhhw, sphhw);
	}
	if (fhhwR) {
	  fhhwR->cd();
	  thhwRF->GetEntry(event);
	  FillFriends(fFhhwR, sphhwR);
	}
      }
      else {
	if (fo) {
	  TDirectoryFile *do = fo->Get(Form("Event%d",event));
	  FillRecPoints(dynamic_cast<TTree*>(do->Get("TreeR")),spo);
	}
	if (fhhw) {
	  TDirectoryFile *dhhw = fhhw->Get(Form("Event%d",event));
	  FillRecPoints(dynamic_cast<TTree*>(dhhw->Get("TreeR")),sphhw);
	}
	if (fhhwR) {
	  TDirectoryFile *dhhwR = fhhwR->Get(Form("Event%d",event));
	  FillRecPoints(dynamic_cast<TTree*>(dhhwR->Get("TreeR")),sphhwR);
	}
      }
    } // for (Int_t event=0 ; event < nEvents; ++event) {

    // -------------------------------------------------------
    // -- Close input files
    // -------------------------------------------------------
    printf(" -- Close files -- \n");

    if (fo)    { fo->Close();    fo = NULL; }
    if (fhhw)  { fhhw->Close();  fhhw = NULL; }
    if (fhhwR) { fhhwR->Close(); fhhwR = NULL; }

  } // for (Int_t fileIdx = minId; fileIdx <= maxId; ++fileIdx) {

  // --------------------------------------------------------
  // -- Write results
  // --------------------------------------------------------

  gSystem->Exec("if [ ! -d ./results ] ; then mkdir -p results ; fi");

  const Char_t *source[] = {"recPoints", "friends"};
  TFile *resultFile = TFile::Open(Form("results/results_%s_%s_%s.root", source[useFriends], type, version),"RECREATE");

  if (spo)    spo->Write();
  if (sphhw)  sphhw->Write();
  if (sphhwR) sphhwR->Write();

  resultFile->Close();

  return;
}

// ==================================================================================
void FillFriends(AliESDfriend* esdFriend, THnSparseF *sparse) {
  /// -- Fill friends per event (out of TPC.RecPoints)

  for (Int_t nTrk = 0; nTrk < esdFriend->GetNumberOfTracks(); nTrk++) {

    AliESDfriendTrack* track = esdFriend->GetTrack(nTrk);
    if (!track)
	continue;

    TObject *calibObject;
    AliTPCseed *seed = NULL;
    for (Int_t l=0; (calibObject=track->GetCalibObject(l)) ; ++l)
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;

    if (!seed)
      continue;

    for (Int_t iRow = 0; iRow < 160 ; ++iRow)
      FillClusters(sparse, seed->GetClusterFast(iRow));
  }
}

// ==================================================================================
void FillRecPoints(TTree* tree, THnSparseF *sparse) {
  /// -- Fill friends per event (out of AliESDfriend.root)

  AliTPCClustersRow *row = new AliTPCClustersRow();
  tree->GetBranch("Segment")->SetAddress(&row);

  for ( Int_t entry =0 ; entry < tree->GetEntriesFast() ; ++entry) {
    tree->GetEntry(entry);
    TClonesArray* c = row->GetArray();
    if (!c) continue;

    for (Int_t idx = 0; idx < c->GetEntriesFast(); ++idx)
      FillClusters(sparse, static_cast<AliTPCclusterMI*>(c->UncheckedAt(idx)));

  } // for ( Int_t entry =0 ; entry <  tree->GetEntriesFast() ; ++entry) {
}

// ==================================================================================
void FillClusters(THnSparseF *sparse,  AliTPCclusterMI* cl) {
  /// -- Fill cluster into THnSparse

  if (!cl)
    return;

  Float_t xyz[3];
  cl->GetGlobalXYZ(xyz);

  Double_t arr[15];

  arr[0]  = cl->GetX();
  arr[1]  = cl->GetY();
  arr[2]  = cl->GetZ();
  arr[3]  = xyz[0];
  arr[4]  = xyz[1];
  arr[5]  = xyz[2];
  /*
      arr[3]  = 0.;   TODO
      arr[4]  = 0.;   TODO
      arr[5]  = 0.;   TODO
  */
  arr[6]  = cl->GetSigmaY2();
  arr[7]  = cl->GetSigmaZ2();
  arr[8]  = cl->GetY()/cl->GetX();
  arr[9]  = cl->GetTimeBin();
  arr[10] = cl->GetRow();

  if (cl->GetDetector() > 35)
    arr[10] += 63;

  arr[11] = cl->GetPad();
  arr[12] = cl->GetDetector();
  arr[13] = cl->GetQ();
  arr[14] = cl->GetMax();

  sparse->Fill(arr);
}

// ==================================================================================
void SetupStyle() {
  /// -- setup style

  gROOT->SetStyle("Plain");

  gStyle->SetHatchesSpacing(0.8);
  gStyle->SetHatchesLineWidth(1);

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetFillStyle(1001);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(10);

  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);

  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);

  gStyle->SetLegendBorderSize(0);

  Int_t font = 42;

  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);

  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);

  gStyle->SetTickLength(0.02,"xy");
  gStyle->SetEndErrorSize(3);

  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");

  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.3,"xyz");
  gStyle->SetTitleSize(0.04,"xyz");
  gStyle->SetTitleSize(0.04);

  gStyle->SetMarkerSize(1.2);
  gStyle->SetPalette(1,0);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetLineWidth(1);

  return;
}

// ==================================================================================
THnSparseF* CreateTHnSparse(Char_t* name) {
  /// -- Create new THnSparse

  // -- Setup THnSparse binning
  Int_t    bin[15] = { 60, 101,  100,    10,   10,   10, 100, 100, 100,  100,  160,   75,  76,  100, 100 };
  Double_t min[15] = { 80.,-80.,-350.,-350.,-350.,-350.,  0.,  0.,-0.2,   0., -0.5,   0.,  0.,   0.,   0.};
  Double_t max[15] = {280., 80., 350., 350., 350., 350.,  1.,  1., 0.2,1000., 159.5., 150., 75., 900., 900.};


  THnSparseF *sp = new THnSparseF(name, "x:y:z:globalX:globalY:globalZ:sigmaY2:sigmaZ2:xRy:timebin:row:pad:Det:Q:Max",
				 15, bin, min, max);

  sp->GetAxis(0)->SetTitle("Local X");
  sp->GetAxis(1)->SetTitle("Local Y");
  sp->GetAxis(2)->SetTitle("Local Z");
  sp->GetAxis(3)->SetTitle("Global X");
  sp->GetAxis(4)->SetTitle("Global Y");
  sp->GetAxis(5)->SetTitle("Global Z");
  sp->GetAxis(6)->SetTitle("#sigma Y^{2}");
  sp->GetAxis(7)->SetTitle("#sigma Z^{2}");
  sp->GetAxis(8)->SetTitle("Local X/Local Y");
  sp->GetAxis(9)->SetTitle("Timebin");
  sp->GetAxis(10)->SetTitle("Row");
  sp->GetAxis(11)->SetTitle("Pad");
  sp->GetAxis(12)->SetTitle("Sector");
  sp->GetAxis(13)->SetTitle("Q_{tot}");
  sp->GetAxis(14)->SetTitle("Q_{max}");

  return sp;
}
