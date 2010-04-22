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

  Int_t iLoader             = 1;
  Int_t iESD                = 1;
  Int_t iAOD                = 0;
  Int_t iHits               = 1;
  Int_t iDigits             = 1;
  Int_t iClusters           = 1;

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  // runloader check already in AssertRunLoader function 
  AliESDEvent* esd = 0x0;
  if(iESD) esd = AliEveEventManager::AssertESD();
  // esd check already in AssertESD function 
  AliEMCALLoader *emcl = dynamic_cast<AliEMCALLoader*> (rl->GetDetectorLoader("EMCAL"));
  Int_t evtID = AliEveEventManager::GetMaster()->GetEventId();
  if(evtID != evtNum) AliEveEventManager::GetMaster()->GotoEvent(evtNum);

  TTree* ht = 0x0; 
  TTree* dt = 0x0; 
  TTree* rt = 0x0; 
  if(iLoader)
    {
      //Load Hits
      if(iHits) {
	if(!rl->LoadHits("EMCAL"))
	  ht = rl->GetTreeH("EMCAL",false);
	else {printf("Please make sure a have a EMCal.Hits.root file \n"); return;}
      }
      //Load Digits
      if(iDigits) {
	if(!rl->LoadDigits("EMCAL"))
	  dt = rl->GetTreeD("EMCAL",false);
	else {printf("Please make sure a have a EMCal.Digits.root file \n"); return;}
      }
      //Load RecPoints
      if(iClusters) {
	if(!rl->LoadRecPoints("EMCAL"))
	  rt = rl->GetTreeR("EMCAL",false);
	else {printf("Please make sure a have a EMCal.RecPoints.root file \n"); return;}
      }
    }

  //  gGeoManager = gEve->GetDefaultGeometry();
  AliEveEventManager::AssertGeometry();
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
  emcal_data = new AliEveEMCALData(rl,node,m);
  if(iESD) emcal_data->SetESD(esd);

  // Get information from RunLoader
  if(iLoader)
    {
      if(iHits)     emcal_data->LoadHits(ht); // Does not work with my aliroot version 
      if(iDigits)   emcal_data->LoadDigits(dt);
      if(iClusters) emcal_data->LoadRecPoints(rt);
      
      rl->GetEvent(evtNum);

      if(iHits)     emcal_data->LoadHitsFromEMCALLoader(emcl);       
      if(iDigits)   emcal_data->LoadDigitsFromEMCALLoader(emcl);     
      if(iClusters) emcal_data->LoadRecPointsFromEMCALLoader(emcl); 
    }

  // Get information from ESDs
  if(iESD)
    {
      rl->GetEvent(evtNum);
      if(iDigits) emcal_data->LoadDigitsFromESD();
      if(iClusters) emcal_data->LoadRecPointsFromESD();
    }

  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  l->SetMainColor(Color_t(2));
  gEve->AddElement(l);

  for (Int_t sm=0; sm<12; sm++)
    {
      AliEveEMCALSModule* esm = new AliEveEMCALSModule(sm,Form("SM %d Element \n", sm),"test");
      //      esm->SetSModuleID(sm);
      esm->SetDataSource(emcal_data);
      esm->UpdateQuads();
      l->AddElement(esm);
    }

  gEve->Redraw3D(kTRUE);

  gEve->EnableRedraw();


}
