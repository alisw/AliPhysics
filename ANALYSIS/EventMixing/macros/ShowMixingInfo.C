Int_t ShowMixingInfo(TString filename="MixInfo.root") {

   Int_t num=0;

   if (gSystem->Load("libTree.so") < 0) {num++; return num;}
   if (gSystem->Load("libGeom.so") < 0) {num++; return num;}
   if (gSystem->Load("libVMC.so") < 0) {num++; return num;}
   if (gSystem->Load("libMinuit.so") < 0) {num++; return num;}
   if (gSystem->Load("libPhysics.so") < 0) {num++; return num;}
   if (gSystem->Load("libSTEERBase.so") < 0) {num++; return num;}
   if (gSystem->Load("libESD.so") < 0) {num++; return num;}
   if (gSystem->Load("libAOD.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSIS.so") < 0) {num++; return num;}
   if (gSystem->Load("libOADB.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSISalice.so") < 0) {num++; return num;}
   if (gSystem->Load("libEventMixing.so") < 0) {num++; return num;}

   TFile *f = TFile::Open(filename.Data(),"READ");


   TList *fOutputList = (TList*)f->Get("cMixInfoList");
   if (!fOutputList) {
      AliError("fOutputList not available");
      return;
   }
   fOutputList->Print();
   fMixInfo = (AliMixInfo *) fOutputList->FindObject("mixInfo");
   if (fMixInfo) {
      fMixInfo->Draw("HIST");
      AliMixEventPool *evPool = (AliMixEventPool *) fMixInfo->GetEventPool("mixEventPool");
      if (evPool) evPool->Print();
   }


   return 0;
}
