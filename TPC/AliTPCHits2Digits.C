
void AliTPCHits2Digits(const  char * name= "pokusD_")
{
 
  // Dynamically link some shared libs
        if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
      loadlibs();
     }
  
   //names of trees
   const char * inFile = "galice.root";
   //   const * char ident= "TreeD1par_";

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
   if (file) file->Close();
   file = new TFile(inFile,"UPDATE");
// Get AliRun object from file or create it if not on file

   if(gAlice){
     delete gAlice;
     gAlice=0;
   }
      if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
        }
   gAlice->GetEvent(0);
   AliTPC *TPC = (AliTPC*)gAlice->GetModule("TPC");      
   TPC->Dump();
   //adjust parameters
  
  AliTPCD *paramd = TPC->GetDigParam();
  paramd->Dump();
  paramd->SetName("Param1");  
  paramd->MakeTree();
  //set pointers to parameters
  paramd->Dump();
  AliTPCParam &param = paramd->GetParam();
  AliTPCPRF2D &prf = paramd->GetPRF2D();
  AliTPCRF1D  & rf  = paramd->GetRF();
  
  param.SetPadLength(2.05);
  param.SetPadWidth(0.35);
  param.SetPadPitchLength(2.05);
  param.SetPadPitchWidth(0.35);
  param.SetNWires(5);
  param.SetZeroSup(5);
  param.SetDiffT(0.022);
  param.SetDiffL(0.022);
  param.SetNoise(500);
  param.SetGasGain(1.e4);
  param.SetChipGain(24); 
  param.SetSectorAngles(40.,0.,20.,10.);
  param.SetInnerRadiusLow(83.7);
  param.SetInnerRadiusUp(132.9);
  param.SetOuterRadiusLow(146.9);
  param.SetOuterRadiusUp(249.4);     
  param.Update();

    //Set z (time) response function

  rf.SetOffset(3.*param.GetZSigma());
  rf.SetGauss(param.GetZSigma(),param.GetZWidth(),0.4);
  rf.Update();
  //Set two dimensional pad response function
  TFile f("TPC/AliTPCprf2d.root");
  //  prf.Read("prf_205035_Gati_062074_d03");
  prf.Read("prf_205035_Gati_062074_d03");
  f.Close();
  
  printf("**********Digit object dump start********************\n");
  paramd->Dump();
  printf("**********AliTPCParam**************************\n");
  param.Dump();
  printf("**********Time response function***************\n");
  rf.Dump();
  printf("**********Pad response function params*********\n");
  prf.Dump();
  printf("**********Digit object dump end********************\n");

   TPC->Hits2DigitsSector(0);     
         
 
   file->cd();
   TPC->GetDigParam()->Write();
};


