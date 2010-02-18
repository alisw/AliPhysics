#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TArray.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <TGrid.h>
#include <TAlienCollection.h>
#include <TGridCollection.h>
#include <TNtuple.h>
#include <TGeoManager.h>
#include "AliGeomManager.h"
#include "AliAlignmentTracks.h"
#include "AliTrackFitter.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackResidualsFast.h"
#include "AliTrackResidualsChi2.h"
#include "AliTrackResidualsLinear.h"
#endif

TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName,Int_t nFiles);

//**************************************************************************************************************************************************
void HmpAlignNew(Int_t runNum=105160, const Int_t iterations=1,const Int_t minTrackPoint=30,TString fileintro="AliTrackPoints.root",TString geometryfile="geometry.root"){

  Int_t   nFile2Xml = 1;         //on Grid - how many ESD files to chain in xml
  TString xmlName;

  xmlName=Form("run000%d.xml",runNum);

  TChain *chain = new TChain("esdTree");

  TGrid::Connect("alien://");
  chain = CreateChainFromCollection(xmlName.Data(),"esdTree",nFile2Xml);   

  if(!gGeoManager) AliGeomManager::LoadGeometry(geometryfile.Data());

  AliAlignmentTracks *AliAlTrack=new AliAlignmentTracks();

  AliAlTrack->AddESD(chain);
  
  AliAlTrack->ProcessESD(); 

  TStopwatch *timer=new TStopwatch();
  timer->Start();
 
  AliAlTrack->SetPointsFilename(fileintro.Data());

  AliAlTrack->BuildIndex();

//  AliTrackFitter *fitter=new AliTrackFitterKalman();
  AliTrackFitter *fitter=new AliTrackFitterRieman();
  fitter->SetMinNPoints(minTrackPoint);

  AliAlTrack->SetTrackFitter(fitter);

//  AliTrackResiduals *res = new AliTrackResidualsFast();
  AliTrackResiduals *res = new AliTrackResidualsChi2();
  res->SetMinNPoints(1);
  AliAlTrack->SetMinimizer(res);

  for(Int_t i = 0; i<7; i++){
    
  TArrayI tmpId(1); 
  tmpId.AddAt(AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,i),0);  
  const TArrayI *volId = &tmpId;
   
  AliAlTrack->AlignVolumes(volId, NULL, AliGeomManager::kTPC1,AliGeomManager::kTPC2,iterations);  
 }

  timer->Stop();
  timer->Print();
  return;
}
//*************************************************************************************************************************************************************
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName,Int_t nFiles)
{
// Create a chain from an alien collection.                                                                          
  
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  

  
  TAlienCollection *myCollection  = TAlienCollection::Open(xmlfile);

   
   if (!myCollection) {
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
     return NULL ;
   }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  Int_t iCount = 0;
  while ( myCollection->Next() ){
    if(nFiles!=0)iCount++;
    if(iCount > nFiles)break;
    chain->Add(myCollection->GetTURL("")) ;
    Printf("Adding %s",myCollection->GetTURL(""));
  }
  chain->ls();
  return chain;
}    
