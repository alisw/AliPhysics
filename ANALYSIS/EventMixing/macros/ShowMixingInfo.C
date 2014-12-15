Int_t ShowMixingInfo(TString filename="MixInfo.root") {

   Int_t num=0;

   if (gSystem->Load("libTree") < 0) {num++; return num;}
   if (gSystem->Load("libGeom") < 0) {num++; return num;}
   if (gSystem->Load("libVMC") < 0) {num++; return num;}
   if (gSystem->Load("libMinuit") < 0) {num++; return num;}
   if (gSystem->Load("libPhysics") < 0) {num++; return num;}
   if (gSystem->Load("libSTEERBase") < 0) {num++; return num;}
   if (gSystem->Load("libESD") < 0) {num++; return num;}
   if (gSystem->Load("libAOD") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSIS") < 0) {num++; return num;}
   if (gSystem->Load("libOADB") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSISalice") < 0) {num++; return num;}
   if (gSystem->Load("libEventMixing") < 0) {num++; return num;}

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
