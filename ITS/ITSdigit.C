#include "iostream.h"

void ITSdigit (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nsignal  =25, Int_t size=-1) 
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
   }


// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   printf("file %p\n",file);
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
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

   AliITSgeom *geom = ITS->GetITSgeom();

//
// Event Loop
//

   // SDD


   // SDD compression param: 2 fDecrease, 2fTmin, 2fTmax or disable, 2 fTolerance
   //Int_t cp[8]={5,7,17,17,20,20,0,0}; 
   //Int_t cp[8]={0,0,21,21,21,21,0,0};

   Int_t cp[8]={0,0,0,0,0,0,0,0};

   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) {
         res1=new AliITSresponseSDD();
         ITS->SetResponseModel(1,res1);
   }
   res1->SetZeroSupp("2D");
   //res1->SetZeroSupp("1D");
   res1->SetParamOptions("same","same");
   //res1->SetFilenames(" ","$(ALICE_ROOT)/ITS/base.dat","$(ALICE_ROOT)/ITS/2D.dat ");
   //res1->SetNoiseParam(3.,20.);
   res1->SetNoiseParam(0.,0.);
   res1->SetCompressParam(cp);
   res1->SetMinVal(4);

   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) {
       seg1 = new AliITSsegmentationSDD(geom,res1);
       ITS->SetSegmentationModel(1,seg1);
   }

   AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
   ITS->SetSimulationModel(1,sim1);

   // SPD

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
   AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
   ITS->SetSimulationModel(0,sim0);
   // test
   printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());
   printf("SPD pitches %d %d \n",seg0->Dpz(0),seg0->Dpx(0));
   // end test


   // SSD

   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
   res2->SetSigmaSpread(3.,2.);
   AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
   ITS->SetSimulationModel(2,sim2);


   // tests
   printf("sim0 sim1 sim2 %p %p %p\n",sim0,sim1,sim2);
   Float_t n,b;
   res1->GetNoiseParam(n,b);
    printf("SDD: noise baseline %f %f zs option %s data type %s\n",n,b,res1->ZeroSuppOption(),res1->DataType());
   printf("SDD: DriftSpeed %f TopValue %f\n",res1->DriftSpeed(),res1->MagicValue());
   printf("SDD: DiffCoeff %f Qref %f\n",res1->DiffCoeff(),res1->Qref());
   // end tests


   Int_t nbgr_ev=0;

   for (int nev=evNumber1; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " <<nev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       Int_t nbgr_ev=Int_t(nev/nsignal);
       //printf("nbgr_ev %d\n",nbgr_ev);
       ITS->HitsToDigits(nev,nbgr_ev,evNumber2,size," ","All"," ");
       //ITS->HitsToDigits(nev,nbgr_ev,evNumber2,size," ","SSD"," ");
   } // event loop 

   file->Close();
}














