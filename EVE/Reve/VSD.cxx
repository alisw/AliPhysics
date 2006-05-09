// $Header$

#include "VSD.h"
#include <TFile.h>

using namespace Reve;

//______________________________________________________________________
// VSD
//

ClassImp(VSD)

VSD::VSD()
{
  mFile      = 0;
  mDirectory = 0;
  InitTreeVars();
}

VSD::VSD(const Text_t* , const Text_t* )
{
  mFile      = 0;
  mDirectory = 0;
  InitTreeVars();
}

void VSD::InitTreeVars()
{
  fBuffSize = 128*1024;

  mTreeK  = 0;
  //mTreeTR = 0;
  mTreeH  = 0;
  mTreeC  = 0;
  mTreeR  = 0;
  mTreeKK = 0;
  mTreeV0 = 0;
  mTreeGI = 0;

  mpK  = &mK;
  mpH  = &mH;
  mpC  = &mC;
  mpV0 = &mV0;
  mpKK = &mKK;
  mpR  = &mR;
  mpGI = &mGI;
}

/**************************************************************************/
/**************************************************************************/

void VSD::SetDirectory(TDirectory* dir)
{
  mDirectory = dir;
}

/**************************************************************************/
/**************************************************************************/

void VSD::CreateTrees()
{
  mDirectory->cd();
  // TR missing ...
  mTreeK  = new TTree("Kinematics", "Simulated tracks.");
  mTreeH  = new TTree("Hits",       "Combined detector hits.");
  mTreeC  = new TTree("Clusters",   "Reconstructed clusters.");
  mTreeR  = new TTree("RecTracks",  "Reconstructed tracks.");
  mTreeKK = new TTree("RecKinks",   "Reconstructed kinks.");
  mTreeV0 = new TTree("RecV0s",     "Reconstructed V0s.");
  mTreeGI = new TTree("GenInfo",    "Objects prepared for cross query.");
}

void VSD::DeleteTrees()
{
  delete mTreeK;      mTreeK      = 0;
  // delete mTreeTR;     mTreeTR     = 0;
  delete mTreeH;      mTreeH      = 0;
  delete mTreeC;      mTreeC      = 0;
  delete mTreeR;      mTreeR      = 0;
  delete mTreeV0;     mTreeV0     = 0;
  delete mTreeKK;     mTreeKK     = 0;
  delete mTreeGI;     mTreeGI     = 0;
}

void VSD::CreateBranches()
{
  // TR missing ...
  if(mTreeK)
    mTreeK ->Branch("K",  "Reve::MCTrack",  &mpK,  fBuffSize);
  if(mTreeH)
    mTreeH ->Branch("H",  "Reve::Hit",      &mpH,  fBuffSize);
  if(mTreeC)
    mTreeC ->Branch("C",  "Reve::Cluster",  &mpC,  fBuffSize);
  if(mTreeR)
    mTreeR ->Branch("R",  "Reve::RecTrack", &mpR,  fBuffSize);
  if(mTreeKK)
    mTreeKK->Branch("KK", "Reve::RecKink",  &mpKK, fBuffSize);
  if(mTreeV0)
    mTreeV0->Branch("V0", "Reve::RecV0",    &mpV0, fBuffSize);

  if(mTreeGI) {
    mTreeGI->Branch("GI", "Reve::GenInfo",  &mpGI, fBuffSize);
    mTreeGI->Branch("K.", "Reve::MCTrack",  &mpK);
    mTreeGI->Branch("R.", "Reve::RecTrack", &mpR);
  }
}

void VSD::SetBranchAddresses()
{
  // TR missing ...
  if(mTreeK)
    mTreeK ->SetBranchAddress("K",  &mpK);
  if(mTreeH)
    mTreeH ->SetBranchAddress("H",  &mpH);
  if(mTreeC)
    mTreeC ->SetBranchAddress("C",  &mpC);
  if(mTreeR)
    mTreeR ->SetBranchAddress("R",  &mpR);
  if(mTreeKK)
    mTreeKK->SetBranchAddress("KK", &mpKK);
  if(mTreeV0)
    mTreeV0->SetBranchAddress("V0", &mpV0);

  if(mTreeGI) {
    mTreeGI->SetBranchAddress("GI", &mpGI);
    mTreeGI->SetBranchAddress("K.", &mpK);
    mTreeGI->SetBranchAddress("R.", &mpR);
  }
}

void VSD::WriteTrees()
{
  // Does nothing here ...
}

/**************************************************************************/
/**************************************************************************/

void VSD::LoadTrees()
{
  static const Exc_t eH("VSD::LoadTrees ");
  
  if(mDirectory == 0)
    throw(eH + "directory not set.");

  printf("Reading kinematics.\n");
  mTreeK = (TTree*) mDirectory->Get("Kinematics");
  if(mTreeK == 0) {
    printf("%s Kinematics not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

  printf("Reading hits.\n");  
  mTreeH = (TTree*) mDirectory->Get("Hits");
  if(mTreeH == 0) {
    printf("%s Hits not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

  printf("Reading clusters.\n");
  mTreeC = (TTree*) mDirectory->Get("Clusters");
  if(mTreeC == 0) {
    printf("%s Clusters not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

  printf("Reading reconstructed tracks.\n");
  mTreeR = (TTree*) mDirectory->Get("RecTracks");
  if(mTreeR == 0) {
    printf("%s RecTracks not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

  printf("Reading reconstructed kinks. \n");
  mTreeKK =  (TTree*) mDirectory->Get("RecKinks");
  if(mTreeKK == 0) {
    printf("%s Kinks not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

  printf("Reading Reconstructed V0s.\n");
  mTreeV0 =  (TTree*) mDirectory->Get("RecV0s");
  if(mTreeV0 == 0) {
    printf("%s V0 not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }
 
  printf("Reading GenInfo.\n");
  mTreeGI = (TTree*)mDirectory->Get("GenInfo");
  if(mTreeGI == 0) {
    printf("%s GenInfo not available in mDirectory %s.\n", 
           eH.Data(), mDirectory->GetName());
  }

}

void VSD::LoadVSD(const Text_t* vsd_file_name, const Text_t* dir_name)
{
  static const Exc_t eH("VSD::LoadVSD ");

  mFile = TFile::Open(vsd_file_name);
  if(mFile == 0)
    throw(eH + "can not open VSD file '" + vsd_file_name + "'.");

  mDirectory = (TDirectory*) mFile->Get(dir_name);
  if(mDirectory == 0)
    throw(eH + "directory '" + dir_name + "' not found in VSD file '" + vsd_file_name + "'.");
  printf("%p\n", (void*)mDirectory);
  LoadTrees();
  SetBranchAddresses();
}

/**************************************************************************/
/**************************************************************************/
