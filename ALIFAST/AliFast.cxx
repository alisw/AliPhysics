
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast                                                              //
//                                                                      //
// Main class to control the AliFast program.                           //
//                                                                      //
// This class is a general framework for programs that needs to:        //
//    - Initialise some parameters                                      //
//    - Loop on events                                                  //
//    - Print results and save histograms, etc                          //
//                                                                      //
// The event processor AliFast::Make loops on a list of Makers          //
// where each maker performs some task on the event data and generates  //
// results.                                                             //
// New Makers can be inserted by a user without modifying this class.   //
// Note that the order in which the Makers are called is the order      //
// of insertion in the list of Makers.                                  //
// Each Maker is responsible for creating its branch of the Tree.       //
// The following table shows the list of makers currently implemented   //
// The default option to Save the Maker info in the Tree is mentioned.  //
//                                                                      //
//    Maker name        Save in Tree                                    //
//    ==========        ============                                    //
//    MCMaker             NO                                            //
//    TrackMaker          NO                                            //
//                                                                      //
// Makers must derive from the base class AliFMaker.                    //
// AliFMaker provides a common interface to all Makers.                 //
// Each Maker is responsible for defining its own parameters and        //
// histograms.                                                          //
// Each Maker has its own list of histograms.                           //
// Each Maker has an associated companion class corresponding to the    //
// type of physics object reconstructed by the Maker.                   //
// For example, AliFClusterMaker creates AliFCluster objects.           //
//              AliFTriggerMaker creates one single AliFTrigger object. //
// The pointer supporting the created object(s) is defined in AliFMaker //
//   fFruits may point to a single object (eg. AliFTrigger) or to a    //
//           TClonesArray of objects (eg. AliFCluster).                 //
//                                                                      //
// The function AliFast::Maketree must be called after the creation     //
// of the AliFast object to create a Root Tree.                         //
//                                                                      //
// An example of main program/macro to use AliFast is given below:      //
//========================================================================
//void umain(Int_t nevents=100)
//{
//   gROOT->Reset();
//   gSystem->Load("libalifast.so");  // dynamically link the compiled shared library
//
//   // Open the root output file
//   TFile file("alifast.root","recreate","AliFast root file",2);
//   
//   AliFast alifast("alifast");     // create main object to run alifast
//
//   User user;           // create an object of the User class defined in user.C
//
//   alifast.Init();      // Initialise event (maker histograms,etc)
//   alifast.MakeTree();  // Create the Root tree
//
//   gROOT->LoadMacro("user.C");  // compile/interpret user file
//
//   for (Int_t i=0; i<nevents; i++) {
//      if (i%100 == 0) printf("In loop:%d\n",i);
//      alifast.Make(i);       // Generate and reconstruct event
//      user.FillHistograms(); // User has possibility to decide if store event here!
//      alifast.FillTree();
//      alifast.Clear();       // Clear all event lists
//   }
//   alifast.Finish();
//
//   // save objects in Root file
//   alifast.Write();  //save main alifast object (and run parameters)
//}
//========================================================================
//                                                                      //
// This example illustrates how to:                                     //
//    - Load a shared library                                           //
//    - Open a Root file                                                //
//    - Initialise AliFast                                              //
//    - Load some user code (interpreted)                               //
//      This user code may redefine some Maker parameters               //
//    - Make a loop on events                                           //
//    - Save histograms and the main AliFast object and its Makers      //
//                                                                      //
//========================================================================
//  An example of a User class is given below:                          //
//========================================================================
//
//#ifndef user_H
//#define user_H
//
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// User                                                                 //
//                                                                      //
// Example of a user class to perform user specific tasks when running  //
// the ALIfast program.                                                 //
//                                                                      //
// This class illustrates:                                              //
//   - How to set run parameters                                        //
//   - How to create and fill histograms                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//
//class TH1F;
//class AliFast;
//class AliFClusterMaker;
//class AliFPhotonMaker;
//
//class User {
//
//private:
//   TH1F             *fhist1;       //pointer to histogram
//   TH1F             *fhist2;       //pointer to histogram
//   TH1F             *fhist3;       //pointer to histogram
//public:
//               User();
//   void        FillHistograms();
//   void        SetRunParameters();
//
//#endif
//};
//
//_________________________________________________________________________
//User::User() 
//{
//   SetRunParameters();  //change default parameters
//
//         Create a few histograms
//   fhist1 = new TH1F("hist1","Number of tracks per event",100,0,100);
//   fhist2 = new TH1F("hist2","Number of clusters",100,0,100);
//   fhist3 = new TH1F("hist3","Number of isolated muons",20,0,20);
//}
//
//_________________________________________________________________________
//void User::FillHistograms()
//{
////   fhist1.Fill(event->GetNtracks());
////   fhist2.Fill(event->GetNclusters));
////   fhist3.Fill(event->GetNIsoMuons());
//}
//
//_________________________________________________________________________
//void User::SetRunParameters()
//{
//  // change Alifast default parameters
//
//   gAliFast->SetSmearMuonOpt(0);
//   gAliFast->ClusterMaker()->SetGranBarrelEta(0.12);
//   gAliFast->PhotonMaker()->SetMinPT(6.);
//   gAliFast->TriggerMaker()->SetMuoEtaCoverage(2.8);
//
//}
//======================end of User class=================================
//
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TBrowser.h>
#include <TClonesArray.h>
#include "AliRun.h"
#include "AliFast.h"
//#include "AliFMCMaker.h"
#include "AliFTrackMaker.h"
#include "AliFHistBrowser.h"
#include "AliFBigBang.h"
#include "AliFVirtualDisplay.h"


R__EXTERN AliRun * gAlice;
AliFast *gAliFast;

ClassImp(AliFast)


//_____________________________________________________________________________
  //AliFast::AliFast() : TNamed("alifast","The ALICE fast simulation")
  AliFast::AliFast() : AliRun("alifast","The ALICE fast simulation")
{

   fTree          = 0;
   fMakers        = 0;
   fMode          = 0;
   //   fMCMaker       = 0;
   fTrackMaker    = 0;
   fDisplay       = 0;
   fDet           = new AliFDet("Detector","Make AliFast detector");
   gAliFast     = this;
   gAlice       = (AliRun*)this;
}

//_____________________________________________________________________________
//AliFast::AliFast(const char *name, const char *title): TNamed(name,title)
AliFast::AliFast(const char *name, const char *title): AliRun(name,title)
{

   gAliFast      = this;
   fVersion     = 001;       //AliFAST  version number and release date
   fVersionDate = 150399;
   fTree        = 0;
   fMode        = 0;
   fDisplay     = 0;
   
   SetDefaultParameters();

   gROOT->GetListOfBrowsables()->Add(this,"AliFast");

// create the support list for the various lists of AliFast objects
   fMakers  = new TList();

// create "standard" makers and add them to the list of makers (in AliFMaker constructor
// Note that the order in which makers are added to the list of makers is important
// makers will be processed in this order !!

   //fMCMaker       = new AliFMCMaker("MCMaker","Make MC events");
   fTrackMaker    = new AliFTrackMaker("TrackMaker","Make AliFast tracks");
//create detector
   fDet           = new AliFDet("Detector","Make AliFast detector");
}

//_____________________________________________________________________________
AliFast::~AliFast()
{
//   fMakers->Delete();
//   delete fMakers;
}


//______________________________________________________________________________
void AliFast::Browse(TBrowser *b)
{

  if( b == 0) return;

  if (fTree) b->Add(fTree,fTree->GetName());
  // fca
  b->Add(&fHistBrowser, "Histograms");
  b->Add(&fBigBang, "BigBang");
  // fca

  TIter next(fMakers);
  AliFMaker *maker;
  while ((maker = (AliFMaker*)next())) {
     b->Add(maker,maker->GetName());
   }
}

//_____________________________________________________________________________
void AliFast::Clear(Option_t *option)
{
//    Reset lists of event objects
   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->Clear(option);
   }
   //fca   if (fDisplay) fDisplay->Clear();
}

//_____________________________________________________________________________
void AliFast::Draw(Option_t *option)
{
//    Insert current event in graphics pad list

    // Check if the Event Display object has been created
   if (!fDisplay) {
      Error("Draw","You must create an AliFDisplay object first");
      return;
   }

   //fca   fDisplay->Draw(option);
}

//_____________________________________________________________________________
void  AliFast::GetTreeEvent(Int_t event)
{
//    Read event from Tree
   if (fTree) fTree->GetEvent(event);
   fEvent = event;  
}

//_____________________________________________________________________________
void AliFast::Init()
{
//  Initialise detector
   AliFDet *detector=gAliFast->Detector();  
   detector->InitDetParam(); 

//    Initialise makers
   TIter next(fMakers);
   AliFMaker *maker;
   TObject *objfirst, *objlast;
   while ((maker = (AliFMaker*)next())) {
     // save last created histogram in current Root directory
      objlast = gDirectory->GetList()->Last();

     // Initialise maker
      maker->Init();

     // Add Maker histograms in Maker list of histograms
      if (objlast) objfirst = gDirectory->GetList()->After(objlast);
      else         objfirst = gDirectory->GetList()->First();
      while (objfirst) {
         maker->Histograms()->Add(objfirst);
         objfirst = gDirectory->GetList()->After(objfirst);
      }
   }

}

//_____________________________________________________________________________
void AliFast::Paint(Option_t *option)
{
//    Paint AliFast objects

  //fca   fDisplay->Paint(option);
}

//_____________________________________________________________________________
void AliFast::PrintInfo()
{
//     Gives information about versions etc.
   printf("\n\n");
   printf("**************************************************************\n");
   printf("*             AliFast version:00%1d     last update  %6d    *\n",
                                                       fVersion, fVersionDate);
   printf("**************************************************************\n");
   printf("*                                                            *\n");
   printf("*     Simulates and reconstructs  events on particle level   *\n");
   printf("*     Package by: Yiota Foka and Elzbieta Richter-Was        *\n");
   printf("*           Design based on ATLFast++                        *\n");
   printf("*         by R. Brun and E. Richter-Was                      *\n");
   printf("**************************************************************\n");
   printf("\n\n");

//     Print info for detector geometry
   AliFDet *detector=gAliFast->Detector();
   detector->PrintDetInfo(); 

//     Print info for all defined Makers
   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->PrintInfo();
   }
}

//_____________________________________________________________________________
void AliFast::FillTree()
{
//  Fill the ROOT tree, looping on all active branches

  // Clean generated particles (depending on option Save)
  //MCMaker()->CleanParticles();

  // Now ready to fill the Root Tree
   if(fTree) fTree->Fill();
}

//_____________________________________________________________________________
void AliFast::InitChain(TChain *chain)
{
//  Initialize branch addresses for all makers in a TChain

   if (chain == 0) return;

   fTree = chain;

   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->SetChainAddress(chain);
   }
}

//_____________________________________________________________________________
void AliFast::MakeTree(const char* name, const char*title)
{
//  Create a ROOT tree
//  Loop on all makers to create the Root branch (if any)

   if (fTree) return;

   fTree = new TTree(name,title);

   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->MakeBranch();
   }
}

//_____________________________________________________________________________
void AliFast::SetDefaultParameters()
{

//    Setters for flags and switches
   SetLuminosity();
   SetBfield();
   SetSmearing();
   SetSUSYcodeLSP();
   SetTrackFinding();
}

//_____________________________________________________________________________
void AliFast::Make(Int_t i)
{
   fEvent = i;

//   Loop on all makers
   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->Make();
   }

}

//_____________________________________________________________________________
void AliFast::FillClone()
{
   // Fill Makers fruits clones
   
   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->FillClone();
   }
}

//_____________________________________________________________________________
void AliFast::Finish()
{
//    Terminate a run
//   place to make operations on histograms, normalization,etc.

   TIter next(fMakers);
   AliFMaker *maker;
   while ((maker = (AliFMaker*)next())) {
      maker->Finish();
   }
}

//_____________________________________________________________________________
void AliFast::SortDown(Int_t n1, Float_t *a, Int_t *index, Bool_t down)
{
//  sort the n1 elements of array a.
//  In output the array index contains the indices of the sorted array.
//  if down is false sort in increasing order (default is decreasing order)
//   This is a translation of the CERNLIB routine sortzv (M101)
//   based on the quicksort algorithm

   Int_t i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_t i22 = 0;
   Float_t ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n1;i++) index[i]--;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}
//______________________________________________________________________________
void AliFast::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliFast.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      if (!gAliFast) gAliFast = this;
      gROOT->GetListOfBrowsables()->Add(this,"AliFast");
      
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      
      AliFast::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);

      fTree = (TTree*)gDirectory->Get("T");
   } else {
       AliFast::Class()->WriteBuffer(R__b,this);
   }
}





