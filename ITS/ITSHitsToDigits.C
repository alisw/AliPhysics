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
   Float_t baseline = 10.;
   Float_t noise = 1.75;

                
   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) {
         res1=new AliITSresponseSDD();
         ITS->SetResponseModel(1,res1);
   }
   res1->SetMagicValue(900.);

   Float_t maxadc = res1->MaxAdc();    
   Float_t topValue = res1->MagicValue();
   Float_t norm = maxadc/topValue;

   Float_t fCutAmp = baseline + 2.*noise;
   fCutAmp *= norm;
   Int_t cp[8]={0,0,fCutAmp,fCutAmp,0,0,0,0}; //1D

   //res1->SetZeroSupp("2D");
   res1->SetZeroSupp("1D");
   res1->SetNoiseParam(noise,baseline);
   res1->SetDo10to8(kTRUE);
   res1->SetMinVal(4);
   res1->SetCompressParam(cp);
   res1->SetDiffCoeff(3.6,40.);

   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) {
       seg1 = new AliITSsegmentationSDD(geom,res1);
       ITS->SetSegmentationModel(1,seg1);
   }

   AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
   sim1->SetDoFFT(1);
   sim1->SetCheckNoise(kFALSE);

   ITS->SetSimulationModel(1,sim1);
   
   

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
   gAlice->MakeTree("D");
   printf("TreeD %p\n",gAlice->TreeD());
   //make branch
   ITS->MakeBranch("D");

   Int_t nbgr_ev=0;
   for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {
       cout << "nev         " <<nev<<endl;
       if(nev>0) {
	 nparticles = gAlice->GetEvent(nev);
	 gAlice->SetEvent(nev);
	 gAlice-> MakeTree("D");
	 ITS->MakeBranch("D");
       }
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       Int_t nbgr_ev=0;
       if(nsignal) nbgr_ev=Int_t(nev/nsignal);
       ITS->HitsToDigits(nev,nbgr_ev,size," ","All"," ");
   } // event loop 

   delete sim0;
   delete sim1;
   delete sim2;


   file->Close();
}














