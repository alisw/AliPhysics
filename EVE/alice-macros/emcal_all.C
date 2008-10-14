//************************************************************************
// A macro to read and visualize EMCAL data
// The macro: 
// - can read hits, digits and clusters information from AliRunLoader:
//      emcal_data->LoadHits(ht); 
//      emcal_data->LoadDigits(dt);
//      emcal_data->LoadRecPoints(rt);
// - can read hits, digits and clusters information from AliEMCALLoader:
//      rl->GetEvent(evtNum);
//      emcal_data->LoadHitsFromEMCALLoader(emcl);       // Does not work
//      emcal_data->LoadDigitsFromEMCALLoader(emcl);     
//      emcal_data->LoadRecPointsFromEMCALLoader(emcl); 
// - can read hits, digits and clusters information from ESDs
//      emcal_data->LoadDigitsFromESD();
//      emcal_data->LoadClustersFromESD();
// - will read hits, digits and clusters information from raw
//   => To be implemented
//
//************************************************************************
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//************************************************************************

#include <Riostream.h>
#include <TMath.h>


class AliEveEMCALData;
AliEveEMCALData     *emcal_data       = 0;

void emcal_all(const UInt_t evtNum = 0, Bool_t digFile = 0, 
	       const UInt_t eventsToProcess = 5, TString dirName = "./", 
	       const TString esdTreeName = "esdTree", const char *  pattern = ".")
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  if (rl == 0x0)
    cout<<"Can not instatiate the Run Loader"<<endl;
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  AliEMCALLoader *emcl = dynamic_cast<AliEMCALLoader*> (rl->GetDetectorLoader("EMCAL"));
  Int_t evtID = AliEveEventManager::GetMaster()->GetEventId();
  if(evtID != evtNum) AliEveEventManager::GetMaster()->GotoEvent(evtNum);

  //Load Hits
  rl->LoadHits("EMCAL");
  //Load Digits
  rl->LoadDigits("EMCAL");
  //Load RecPoints
  rl->LoadRecPoints("EMCAL");
  
  TTree* ht = rl->GetTreeH("EMCAL",false);
  TTree* dt = rl->GetTreeD("EMCAL",false);
  TTree *rt = rl->GetTreeR("EMCAL",false);

  gGeoManager = gEve->GetDefaultGeometry();
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
  emcal_data = new AliEveEMCALData(rl,node,m);
  emcal_data->SetESD(esd);
  // RunLoader information
//  emcal_data->LoadHits(ht); // Does not work with my aliroot version ?
//  emcal_data->LoadDigits(dt);
//  emcal_data->LoadRecPoints(rt);

  // To be uncommented if use of emcalLoader
  rl->GetEvent(evtNum);
  emcal_data->LoadHitsFromEMCALLoader(emcl);       
//  emcal_data->LoadDigitsFromEMCALLoader(emcl);     
//  emcal_data->LoadRecPointsFromEMCALLoader(emcl); 

  // To be uncommented to read esds
  emcal_data->LoadDigitsFromESD();
  emcal_data->LoadRecPointsFromESD();

  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  l->SetMainColor(Color_t(2));
  gEve->AddElement(l);

  for (Int_t sm=0; sm<12; sm++)
    {
      cout << "\n Form: " << Form("SM %d", sm) << endl;
      AliEveEMCALSModule* esm = new AliEveEMCALSModule(sm,Form("SM %d Element \n", sm),"test");
      //      esm->SetSModuleID(sm);
      esm->SetDataSource(emcal_data);
      esm->ComputeBBox();
      esm->UpdateQuads();
      l->AddElement(esm);
    }

  gEve->Redraw3D(kTRUE);

  gEve->EnableRedraw();


}
