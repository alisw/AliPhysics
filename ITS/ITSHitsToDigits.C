#include "iostream.h"

void ITSHitsToDigits (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nsignal  =25, Int_t size=-1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice=0;
   }


// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   printf("file %p\n",file);
   if (file) file->Close(); 
   if (!file) file = new TFile("galice.root","UPDATE");
   file->ls();

   printf ("I'm after Map \n");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   printf ("I'm after gAlice \n");
   
   AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
   if (!ITS) return;


   // Set the simulation models

   AliITSgeom *geom = ITS->GetITSgeom();

   // SDD
   // SDD compression param: 2 fDecrease, 2fTmin, 2fTmax or disable, 2 fTolerance
                
   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) {
         res1=new AliITSresponseSDD();
         ITS->SetResponseModel(1,res1);
   }

   //res1->SetChargeLoss(0.);
   Float_t baseline;
   Float_t noise;
   res1->GetNoiseParam(noise,baseline);
   Float_t noise_after_el = res1->GetNoiseAfterElectronics();
   cout << "noise_after_el: " << noise_after_el << endl; 
   Float_t fCutAmp;
   fCutAmp = baseline;
   fCutAmp += (2.*noise_after_el);  // noise
   cout << "Cut amplitude: " << fCutAmp << endl;
   Int_t cp[8]={0,0,fCutAmp,fCutAmp,0,0,0,0};
   res1->SetCompressParam(cp);
   //   res1->SetElectronics(2);  // 1 = Pascal, 2 = OLA

   res1->Print();

   //cout << "SDD segmentation" << endl;

   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) {
       seg1 = new AliITSsegmentationSDD(geom,res1);
       ITS->SetSegmentationModel(1,seg1);
   }
   seg1->Print();

   //cout << "SDD segmentation" << endl;
   AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
   ITS->SetSimulationModel(1,sim1);
   sim1->Print();
   
   

   // SPD

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
   AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
   ITS->SetSimulationModel(0,sim0);
   // test
   //printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   //printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());
   //printf("SPD pitches %d %d \n",seg0->Dpz(0),seg0->Dpx(0));
   // end test


   // SSD

   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
   res2->SetSigmaSpread(3.,2.);
   AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
   ITS->SetSimulationModel(2,sim2);


//
// Event Loop
//


   // create the TreeD 

   Int_t nparticles=gAlice->GetEvent(0);
   printf("Create TreeD \n");
   if(!gAlice->TreeD()) gAlice->MakeTree("D");
   //make branch
   ITS->MakeBranch("D");

   Int_t nbgr_ev=0;
	
	
	cout<<"Digitizing ITS...\n";
   TStopwatch timer;
	
   for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {
       cout << "nev         " <<nev<<endl;
       if(nev>0) {
	 nparticles = gAlice->GetEvent(nev);
	 gAlice->SetEvent(nev);
	 if(!gAlice->TreeD()) gAlice-> MakeTree("D");
	 ITS->MakeBranch("D");
       }
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       Int_t nbgr_ev=0;
       if(nsignal) nbgr_ev=Int_t(nev/nsignal);
       timer.Start();
       ITS->HitsToDigits(nev,nbgr_ev,size," ","All"," ");
       timer.Stop(); timer.Print();
   } // event loop 

   delete sim0;
   delete sim1;
   delete sim2;


   file->Close();
}














