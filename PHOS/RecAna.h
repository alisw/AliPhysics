//////////////////////////////////////////////////////////
//  A Demo Macro that shows how to analyze the 
//  reconstructed Tree
//
//  Y. Schutz (SUBATECH)
//////////////////////////////////////////////////////////


#ifndef RecAna_h
#define RecAna_h

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#endif

const Int_t kMaxPHOSTS = 9;
const Int_t kMaxPHOSRP = 9;

class RecAna {
  public:
  AliPHOSv0 * fPHOS ; 
  TTree          *fTree;    //pointer to the analyzed TTree or TChain
  TTree          *fCurrent; //pointer to the current TTree
  //Declaration of leaves types
  TObjArray       *PHOSEmcRP;
  TObjArray       *PHOSPpsdRP;
  Int_t           PHOSTS_;
  Int_t           PHOSTS_fEmcRecPoint[kMaxPHOSTS];
  Int_t           PHOSTS_fPpsdLowRecPoint[kMaxPHOSTS];
  Int_t           PHOSTS_fPpsdUpRecPoint[kMaxPHOSTS];
  UInt_t          PHOSTS_fUniqueID[kMaxPHOSTS];
  UInt_t          PHOSTS_fBits[kMaxPHOSTS];
  Int_t           PHOSRP_;
  Int_t           PHOSRP_fPHOSTrackSegment[kMaxPHOSRP];
  Int_t           PHOSRP_fIndexInList[kMaxPHOSRP];
  Int_t           PHOSRP_fPrimary[kMaxPHOSRP];
  Int_t           PHOSRP_fType[kMaxPHOSRP];
  Int_t           PHOSRP_fPdgCode[kMaxPHOSRP];
  Int_t           PHOSRP_fStatusCode[kMaxPHOSRP];
  Int_t           PHOSRP_fMother[2][kMaxPHOSRP];
  Int_t           PHOSRP_fDaughter[2][kMaxPHOSRP];
  Float_t         PHOSRP_fWeight[kMaxPHOSRP];
  Double_t        PHOSRP_fCalcMass[kMaxPHOSRP];
  Double_t        PHOSRP_fPx[kMaxPHOSRP];
  Double_t        PHOSRP_fPy[kMaxPHOSRP];
  Double_t        PHOSRP_fPz[kMaxPHOSRP];
  Double_t        PHOSRP_fE[kMaxPHOSRP];
  Double_t        PHOSRP_fVx[kMaxPHOSRP];
  Double_t        PHOSRP_fVy[kMaxPHOSRP];
  Double_t        PHOSRP_fVz[kMaxPHOSRP];
  Double_t        PHOSRP_fVt[kMaxPHOSRP];
  Double_t        PHOSRP_fPolarTheta[kMaxPHOSRP];
  Double_t        PHOSRP_fPolarPhi[kMaxPHOSRP];
  UInt_t          PHOSRP_fUniqueID[kMaxPHOSRP];
  UInt_t          PHOSRP_fBits[kMaxPHOSRP];
  Short_t         PHOSRP_fLineColor[kMaxPHOSRP];
  Short_t         PHOSRP_fLineStyle[kMaxPHOSRP];
  Short_t         PHOSRP_fLineWidth[kMaxPHOSRP];
  
  //List of branches
  TBranch        *b_PHOSEmcRP;
  TBranch        *b_PHOSPpsdRP;
  TBranch        *b_PHOSTS_;
  TBranch        *b_PHOSTS_fEmcRecPoint;
  TBranch        *b_PHOSTS_fPpsdLowRecPoint;
  TBranch        *b_PHOSTS_fPpsdUpRecPoint;
  TBranch        *b_PHOSTS_fUniqueID;
  TBranch        *b_PHOSTS_fBits;
  TBranch        *b_PHOSRP_;
  TBranch        *b_PHOSRP_fPHOSTrackSegment;
  TBranch        *b_PHOSRP_fIndexInList;
  TBranch        *b_PHOSRP_fPrimary;
  TBranch        *b_PHOSRP_fType;
  TBranch        *b_PHOSRP_fPdgCode;
  TBranch        *b_PHOSRP_fStatusCode;
  TBranch        *b_PHOSRP_fMother;
  TBranch        *b_PHOSRP_fDaughter;
  TBranch        *b_PHOSRP_fWeight;
  TBranch        *b_PHOSRP_fCalcMass;
  TBranch        *b_PHOSRP_fPx;
  TBranch        *b_PHOSRP_fPy;
  TBranch        *b_PHOSRP_fPz;
  TBranch        *b_PHOSRP_fE;
  TBranch        *b_PHOSRP_fVx;
  TBranch        *b_PHOSRP_fVy;
  TBranch        *b_PHOSRP_fVz;
  TBranch        *b_PHOSRP_fVt;
  TBranch        *b_PHOSRP_fPolarTheta;
  TBranch        *b_PHOSRP_fPolarPhi;
  TBranch        *b_PHOSRP_fUniqueID;
  TBranch        *b_PHOSRP_fBits;
  TBranch        *b_PHOSRP_fLineColor;
  TBranch        *b_PHOSRP_fLineStyle;
  TBranch        *b_PHOSRP_fLineWidth;
  
  RecAna() {};
  RecAna(char * filename);
  RecAna(TTree *tree) {};
  ~RecAna() {;}
  Int_t GetEntry(Int_t entry = 0);
  Int_t GetEvent(Int_t evt);
  Int_t LoadTree(Int_t entry = 0);
  void  Init(TTree *tree);
  void  Loop();
  void  Notify();
  void  Show(Int_t entry = -1);
};

#endif

#ifdef RecAna_cxx
RecAna::RecAna(char * filename)
{
  // connect the file used to generate this class and read the Tree.
  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!f) 
    f = new TFile(filename);
  
  // setup the gAlice evironment
  
  gAlice = (AliRun*) f->Get("gAlice") ;   
  
  // get the PHOS detectector, the geometry instance and the index to object converter instance
  
  fPHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS") ;     
  AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
  AliPHOSIndexToObject::GetInstance(fPHOS) ;
}

Int_t RecAna::GetEntry(Int_t entry)
{
  // Read contents of entry, always = 0.
  
  if (!fTree) 
    return 0;
  
  return fTree->GetEntry(entry);
}

Int_t RecAna::GetEvent(Int_t evt)
{
  // get the selected event
  
  gAlice->GetEvent(evt);
  
  // connect to the Reconstruction tree
  
  fTree = gAlice->TreeR();     
  
  // set the branches 
  
  Init(fTree);
  
  // gets the data
  
  GetEntry(0); 
  cout << "macro EmcRecpoints = " << fPHOS->EmcRecPoints() << endl ; 
}

Int_t RecAna::LoadTree(Int_t entry)
{
  // Set the environment to read one entry, always = 0.
  
  if (!fTree) 
    return -5;
  Int_t centry = fTree->LoadTree(entry);
  if (centry < 0) 
    return centry;
  if (fTree->GetTree() != fCurrent) {
    fCurrent = fTree->GetTree();
    Notify();
  }
  return centry;
}

void RecAna::Init(TTree *tree)
{
  //   Set branch addresses
  if (tree == 0) return;
  fTree    = tree;
  fCurrent = 0;
  
  fTree->SetBranchAddress("PHOSEmcRP",&PHOSEmcRP);
  fTree->SetBranchAddress("PHOSPpsdRP",&PHOSPpsdRP);
  fTree->SetBranchAddress("PHOSTS_",&PHOSTS_);
  fTree->SetBranchAddress("PHOSTS.fEmcRecPoint",PHOSTS_fEmcRecPoint);
  fTree->SetBranchAddress("PHOSTS.fPpsdLowRecPoint",PHOSTS_fPpsdLowRecPoint);
  fTree->SetBranchAddress("PHOSTS.fPpsdUpRecPoint",PHOSTS_fPpsdUpRecPoint);
  fTree->SetBranchAddress("PHOSTS.fUniqueID",PHOSTS_fUniqueID);
  fTree->SetBranchAddress("PHOSTS.fBits",PHOSTS_fBits);
  fTree->SetBranchAddress("PHOSRP_",&PHOSRP_);
  fTree->SetBranchAddress("PHOSRP.fPHOSTrackSegment",PHOSRP_fPHOSTrackSegment);
  fTree->SetBranchAddress("PHOSRP.fIndexInList",PHOSRP_fIndexInList);
  fTree->SetBranchAddress("PHOSRP.fPrimary",PHOSRP_fPrimary);
  fTree->SetBranchAddress("PHOSRP.fType",PHOSRP_fType);
  fTree->SetBranchAddress("PHOSRP.fPdgCode",PHOSRP_fPdgCode);
  fTree->SetBranchAddress("PHOSRP.fStatusCode",PHOSRP_fStatusCode);
  fTree->SetBranchAddress("PHOSRP.fMother[2]",PHOSRP_fMother);
  fTree->SetBranchAddress("PHOSRP.fDaughter[2]",PHOSRP_fDaughter);
  fTree->SetBranchAddress("PHOSRP.fWeight",PHOSRP_fWeight);
  fTree->SetBranchAddress("PHOSRP.fCalcMass",PHOSRP_fCalcMass);
  fTree->SetBranchAddress("PHOSRP.fPx",PHOSRP_fPx);
  fTree->SetBranchAddress("PHOSRP.fPy",PHOSRP_fPy);
  fTree->SetBranchAddress("PHOSRP.fPz",PHOSRP_fPz);
  fTree->SetBranchAddress("PHOSRP.fE",PHOSRP_fE);
  fTree->SetBranchAddress("PHOSRP.fVx",PHOSRP_fVx);
  fTree->SetBranchAddress("PHOSRP.fVy",PHOSRP_fVy);
  fTree->SetBranchAddress("PHOSRP.fVz",PHOSRP_fVz);
  fTree->SetBranchAddress("PHOSRP.fVt",PHOSRP_fVt);
  fTree->SetBranchAddress("PHOSRP.fPolarTheta",PHOSRP_fPolarTheta);
  fTree->SetBranchAddress("PHOSRP.fPolarPhi",PHOSRP_fPolarPhi);
  fTree->SetBranchAddress("PHOSRP.fUniqueID",PHOSRP_fUniqueID);
  fTree->SetBranchAddress("PHOSRP.fBits",PHOSRP_fBits);
  fTree->SetBranchAddress("PHOSRP.fLineColor",PHOSRP_fLineColor);
  fTree->SetBranchAddress("PHOSRP.fLineStyle",PHOSRP_fLineStyle);
  fTree->SetBranchAddress("PHOSRP.fLineWidth",PHOSRP_fLineWidth);
}

void RecAna::Notify()
{
  //   called by LoadTree when loading a new file
  //   get branch pointers
  b_PHOSEmcRP = fTree->GetBranch("PHOSEmcRP");
  b_PHOSPpsdRP = fTree->GetBranch("PHOSPpsdRP");
  b_PHOSTS_ = fTree->GetBranch("PHOSTS_");
  b_PHOSTS_fEmcRecPoint = fTree->GetBranch("PHOSTS.fEmcRecPoint");
  b_PHOSTS_fPpsdLowRecPoint = fTree->GetBranch("PHOSTS.fPpsdLowRecPoint");
  b_PHOSTS_fPpsdUpRecPoint = fTree->GetBranch("PHOSTS.fPpsdUpRecPoint");
  b_PHOSTS_fUniqueID = fTree->GetBranch("PHOSTS.fUniqueID");
  b_PHOSTS_fBits = fTree->GetBranch("PHOSTS.fBits");
  b_PHOSRP_ = fTree->GetBranch("PHOSRP_");
  b_PHOSRP_fPHOSTrackSegment = fTree->GetBranch("PHOSRP.fPHOSTrackSegment");
  b_PHOSRP_fIndexInList = fTree->GetBranch("PHOSRP.fIndexInList");
  b_PHOSRP_fPrimary = fTree->GetBranch("PHOSRP.fPrimary");
  b_PHOSRP_fType = fTree->GetBranch("PHOSRP.fType");
  b_PHOSRP_fPdgCode = fTree->GetBranch("PHOSRP.fPdgCode");
  b_PHOSRP_fStatusCode = fTree->GetBranch("PHOSRP.fStatusCode");
  b_PHOSRP_fMother = fTree->GetBranch("PHOSRP.fMother[2]");
  b_PHOSRP_fDaughter = fTree->GetBranch("PHOSRP.fDaughter[2]");
  b_PHOSRP_fWeight = fTree->GetBranch("PHOSRP.fWeight");
  b_PHOSRP_fCalcMass = fTree->GetBranch("PHOSRP.fCalcMass");
  b_PHOSRP_fPx = fTree->GetBranch("PHOSRP.fPx");
  b_PHOSRP_fPy = fTree->GetBranch("PHOSRP.fPy");
  b_PHOSRP_fPz = fTree->GetBranch("PHOSRP.fPz");
  b_PHOSRP_fE = fTree->GetBranch("PHOSRP.fE");
  b_PHOSRP_fVx = fTree->GetBranch("PHOSRP.fVx");
  b_PHOSRP_fVy = fTree->GetBranch("PHOSRP.fVy");
  b_PHOSRP_fVz = fTree->GetBranch("PHOSRP.fVz");
  b_PHOSRP_fVt = fTree->GetBranch("PHOSRP.fVt");
  b_PHOSRP_fPolarTheta = fTree->GetBranch("PHOSRP.fPolarTheta");
  b_PHOSRP_fPolarPhi = fTree->GetBranch("PHOSRP.fPolarPhi");
  b_PHOSRP_fUniqueID = fTree->GetBranch("PHOSRP.fUniqueID");
  b_PHOSRP_fBits = fTree->GetBranch("PHOSRP.fBits");
  b_PHOSRP_fLineColor = fTree->GetBranch("PHOSRP.fLineColor");
  b_PHOSRP_fLineStyle = fTree->GetBranch("PHOSRP.fLineStyle");
  b_PHOSRP_fLineWidth = fTree->GetBranch("PHOSRP.fLineWidth");
}

void RecAna::Show(Int_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fTree) return;
  fTree->Show(entry);
}
#endif // #ifdef RecAna_cxx

