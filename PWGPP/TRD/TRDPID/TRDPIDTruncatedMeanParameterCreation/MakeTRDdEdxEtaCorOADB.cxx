void MakeTRDdEdxEtaCorOADB(Int_t index=-1,const Bool_t kMC=0)
{
  //
  //make OADB object
  //currently only default
  //
 gSystem->Load("libOADB"); 

 //Load files
 TFile *inTPCtgl = TFile::Open("output/LHC18d_TPCtglMap_LHC18d_from4to6Tracklets__NclsCorr.root");
 TFile *inNcls = TFile::Open("output/LHC18d_ClusterMap_LHC18d_from4to6Tracklets_.root");
 TFile *inCent  =  TFile::Open("output/LHC18d_CentMap_LHC18d_from4to6Tracklets__NclsCorr.root"); 

 TObjArray array13b_f(5); // original (1)
 array13b_f.SetName("Corrections"); // old EtaCor

 char name[30];
 TH2D *mapEta = (TH2D*)inTPCtgl->Get("mapTPCtgl");
 mapEta->SetName("TRDEtaMap");
 array13b_f.Add(mapEta);

 TH2D *mapNcls[3];
 for (int i=0; i<3; i++)
   {
     mapNcls[i] = (TH2D*)inNcls->Get(Form("NclsMap_nCh%i", i+4));
     mapNcls[i]->SetName(Form("TRDNclsMap_Nch%i", i+4));
     array13b_f.Add(mapNcls[i]);
   }

 if (inCent!=0) {
   cout << "Adding Centrality map" << endl;
   TH2D *TRDCentralityMap = (TH2D*)inCent->Get("CentMap");
   TRDCentralityMap->SetName("TRDCentralityMap");
   array13b_f.Add(TRDCentralityMap);
 }
 else cout << "No Centrality map" << endl;

 
 TString containerName = "TRDCorrectionMaps";
 AliOADBContainer *cont = new AliOADBContainer(containerName.Data());
  
 // current convention as other detetor is only to provide for data
 // no MC is seen in current (31/07/2014)  AliROOT except for HMPID

 //const TString filePathNamePackage="TRDTest.root";
 const TString filePathNamePackage=Form("/home/ikp/alice/AliPhysics/OADB/COMMON/PID/%s/TRDdEdxCorrectionMaps.root", "data");
  
 Int_t statusCont = cont->InitFromFile(filePathNamePackage.Data(), cont->GetName());
 if (statusCont) {
   printf("No OADBContainer for the current settings found - creating a new one...\n");
 }


 Int_t fRunRange[2]={286350, 285978}; // LHC15n
 // Int_t fRunRange[2]={244340, 244628}; // LHC15n
 // Int_t fRunRange[2]= {244824,246994}; // LHC15o
 //  Int_t fRunRange[2] ={195344, 197388};// LHC13 bc
 
 Int_t CheckForExisitingData=-1;
 CheckForExisitingData=cont->GetIndexForRun(fRunRange[0]); //fRunRange[0]);
 cout << "ExistingData" << endl;
 cout << CheckForExisitingData << endl;

 //cont.CleanDefaultList();               // Remove old default objects at first
 //cont.AddDefaultObject((TObject*) &array13b_f);
 //cont.AddDefaultObject((TObject*) &map);
 
 if (CheckForExisitingData!=-1) {
   cout << "UPDATE" << endl;
   cont->UpdateObject(cont->GetIndexForRun(fRunRange[0]), &array13b_f, fRunRange[0], fRunRange[1]);
 }
 else {
   cout << "APPEND" << endl;
   cont->AppendObject(&array13b_f,fRunRange[0],fRunRange[1]);
  //cont.AppendObject(&map,fRunRange[0],fRunRange[1]);
 }

 TFile* f = new TFile(filePathNamePackage.Data(), "update");
 f->Delete(cont->GetName());
 cont->Write(0, TObject::kOverwrite);
 f->Purge();
 f->Close();
 
 printf("TRDdEdxCorrectionMaps done\n");


 //non-default params

 // cont.WriteToFile("TRDdEdxCorrectionMaps.root");
 
 //printf("TRDdEdxCorrectionMaps done\n");
}
