//
// This macro draws time resolution for true and false particle hypothesis
// on a given reference plane in a given momentum range.
// True hypothesis is painted in BLUE while false in RED
//
// Input parameters:
// nRefPlane - indicate the reference plane.
//   1: in TPC  
//   2: out TPC
//   3: out TRD
// 
// minP, maxP - total momentum cuts. momentum from gAlice are used for cuts
// 
// The macro is using SORTED reference points in a file "trackRef.root".
// Track references can be sorted using CopyReferences.C
// 
// Before using check filenames and names of trees invlolved
// 
// ITS tracks (nRefPlane == 1) are compared at inner TPC ref Plane
// make shure they are properly propagated to this plane.
//
// Sylwester Radomski, e-mail: S.Radomski@gsi.de
// Feb 14, 2002 
// 

CompareTimeHypothesis (Int_t nRefPlane, Double_t minP, Double_t maxP) {

  // Check input data Consistancy
  
  if (nRefPlane < 1 || nRefPlane > 3) {
    cout << "Wrong Reference Plane Id ( " << nRefPalne << ") [1-3] " << endl;
    return;
  }

  if (maxP < minP) {
    cout << "MaxP lesser than MinP" << endl;
    return;
  }

  // Styling
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);

  // magnetic field 0.4 T
  Double_t cc = 0.299792458;
  AliKalmanTrack::SetConvConst(100/cc/0.4);

  const char* gFileName = "galice.root";
  const char* refFileName = "trackRef.root";
  
  const char* trackFileName[] = {"AliTPCtracks.root", "AliTPCBackTracks.root", "AliTRDtracks.root"};
  const char *treeName[] = {"TreeT_ITSb_0", "seedsTPCtoTRD_0", "TRDb_0"};
  const char *refTreeName[] = {"TPC", "TPC", "TRD"};
  
  // Histograms

  const Int_t nTypes = 5; 
  Int_t codes[] = {11, 13, 211, 321, 2212};
  
  //Double_t cutoff[] = {20, 20, 20, 100, 100};
  Double_t cutoff[] = {30, 30, 30, 100, 200};

  const char *names[3*nTypes] = {
    "Electrons", "Muons", "Pions", "Kaons", "Protons",
    "11s", "13s", "211s", "321s", "2212s",
    "11o", "13o", "211o", "321o", "2212o"
  };

  THStack *stack[nTypes];
  TH1D *hDiffSame[nTypes], *hDiffOdd[nTypes];

  for(Int_t i=0; i<nTypes; i++) {

    stack[i] = new THStack(names[i], names[i]);
    hDiffSame[i] = new TH1D(names[nTypes+i], names[nTypes+i], 50, -cutoff[i], cutoff[i]);
    hDiffOdd[i] = new TH1D(names[2*nTypes+i], names[2*nTypes+i], 50, -cutoff[i], cutoff[i]);
 
    hDiffSame[i]->SetFillColor(kBlue);
    hDiffOdd[i]->SetFillColor(kRed);
    
    //hDiffSame[i]->SetFillStyle(3003);
    //hDiffOdd[i]->SetFillStyle(3008);

    stack[i]->Add(hDiffOdd[i]);
    stack[i]->Add(hDiffSame[i]);
  }
  
  TH1D *hLength = new TH1D("length", "Length Resolution", 100, -0.5, 0.5);

  // Particles

  TFile *refFile = new TFile(gFileName, "READ");
  gAlice = (AliRun*)refFile->Get("gAlice");
  gAlice->GetEvent(0);


  // Reference tracks

  refFile = new TFile(refFileName, "READ");
  TTree *treeRef = (TTree*)refFile->Get("TreeTR0_Sorted"); 
  TBranch *trkRef = treeRef->GetBranch(refTreeName[nRefPlane-1]);
  TClonesArray *trkRefs = new TClonesArray("AliTrackReference", 100);
  trkRef->SetAddress(&trkRefs);
  
  Int_t ntracks = trkRef->GetEntries();
  cout << ntracks << endl;


  // Tracks
  TFile *trackFile = new TFile(trackFileName[nRefPlane-1]);
  TTree *tracks = (TTree*)trackFile->Get(treeName[nRefPlane-1]);
  TBranch *trkTracks = tracks->GetBranch("tracks");
  AliTPCtrack *track = 0;
  trkTracks->SetAddress(&track);

  for(Int_t t=0; t<tracks->GetEntries(); t++) {

    cout << t << "\r";
    
    trkTracks->GetEvent(t);
    
    Int_t lab = track->GetLabel();
    if (lab < 0 || lab > 100000) continue;
    
    Int_t pdg = gAlice->Particle(lab)->GetPdgCode();
    treeRef->GetEvent(lab);
    
    Int_t nEntrTRK = trkRefs->GetEntries();
    
    // Primaries within momentum range
    if (gAlice->Particle(lab)->P() < minP ) continue; 
    if (gAlice->Particle(lab)->P() > maxP ) continue;
    if (gAlice->Particle(lab)->GetFirstMother() != -1) continue;    
    if (!(nEntrTRK)) continue;
      
    Double_t timeTrue, timeTrack;
    Double_t trueLength;
    
    Int_t index =  nEntrTRK-1;
    if (nRefPlane == 1) index = 0;

    timeTrue = ((AliTrackReference*)(*trkRefs)[index])->GetTime() * 1e12;

    //timeTrue = gRandom->Gaus(timeTrue, 10);

    trueLength = ((AliTrackReference*)(*trkRefs)[index])->GetLength();
    hLength->Fill(trueLength - track->GetIntegratedLength());

    for(Int_t i=0; i<nTypes; i++) {

      timeTrack = track->GetIntegratedTime(codes[i]);
      if (abs(pdg) == codes[i]) hDiffSame[i]->Fill(timeTrue - timeTrack);
      else hDiffOdd[i]->Fill(timeTrue - timeTrack);      
    }
  }  
  
  TCanvas *c = new TCanvas();
  c->SetWindowSize(640, 800);
  c->Divide(2,3);

  c->cd(1);
  hLength->Draw();
  hLength->SetXTitle("True Length - Measured Length [cm]");

  for(Int_t i=0; i<nTypes; i++) {
    c->cd(2+i);
    gPad->SetGridx();
    gPad->SetGridy();
    stack[i]->Draw();  
    stack[i]->GetXaxis()->SetTitle("time true - measured [ps]");    
  }

  if (treeRef) delete treeRef;
  if (tracks) delete tracks;
}
