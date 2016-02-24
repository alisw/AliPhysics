/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//************************************************************************
///
/// \file emcal_all.C
/// \brief Visualize EMCAL data
///
/// A macro to read and visualize EMCAL data in different formats with the help of
/// the classes managing the EMCAL data visulization:
/// * AliEveEMCALData,
/// * AliEveEMCALSModule,
/// * AliEveEMCALSModuleData.
///  
/// This macro: 
/// * can read hits, digits and clusters information from AliRunLoader:
///     * emcal_data->LoadHits();
///     * emcal_data->LoadDigits();
///     * emcal_data->LoadRecPoints();
///     * Same methods with the AliEMCALLoader implemented in the class AliEveEMCALData, (emcal_data->Load(DataType)FromEMCALLoader()) but it does not work.
/// * can read hits (note tested in May 2015 implementation), digits and clusters information from ESDs (easily extendable to AODs, but not implemented)
///     * emcal_data->LoadDigitsFromESD();
///     * emcal_data->LoadClustersFromESD();
/// * raw data not implemented
///
/// \param iLoader: Bool, do the analysis of reconstructed real/simulated data (not working online for the moment)
/// \param iESD: Bool, do the analysis of reconstructed data from ESDs
/// \param iHits: Bool, do the analysis of the generated hits from simulation, iLoader must be on
/// \param iDigits: Bool, do the analysis of digits or ESD Cells (depends on iLoader and iESD setting)
/// \param iClusters: Bool, do the analysis of reconstructed clusters, RecPoints or ESD CaloClusters (depends on iLoader and iESD setting)
///
/// Include it in the macro visscan_init.C in this way, for the different cases:
///   * exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "EMCAL DIGITS", "emcal_all.C", "emcal_all","1,0,0,1,0",kTRUE));
///   * exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "EMCAL REC POINTS", "emcal_all.C", "emcal_all","1,0,0,0,1",kTRUE));
///   * exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "EMCAL ESD CELLS/DIGITS", "emcal_all.C", "emcal_all","0,1,0,1,0",kTRUE));
///   * exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "EMCAL ESD CLUSTERS", "emcal_all.C", "emcal_all","0,1,0,0,1",kTRUE));
/// (the last parameter of the visscan_init macro indicates that this line is active or not).
///
/// \author Magali Estienne <magali.estienne@cern.ch>, SUBATECH. EMCal implementation, June 2008
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. DCal implementation + doxygen, May 2015.
//************************************************************************

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TStyle.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliEMCAL.h>
#include <AliEMCALLoader.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliTrackPointArray.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#include <AliEveEMCALData.h>
#include <AliEveEMCALSModule.h>

#else
class AliEveEMCALData;
#endif

AliEveEMCALData * emcal_data = 0;

void emcal_all
(
 Bool_t iLoader    = 1,
 Bool_t iESD       = 0,
 Bool_t iHits      = 0,
 Bool_t iDigits    = 1,
 Bool_t iClusters  = 0
 )
{
//  printf("------------------------------------------------------------------------------------\n");
//  printf("emcal_all.C - Selected options: Loaders %d, ESDs %d; Hits %d, Digits %d, Clusters %d\n",iLoader,iESD,iHits,iDigits,iClusters);
//  printf("------------------------------------------------------------------------------------\n");
  
  //
  // Get the data mangers, AliRunLoader or AliESDEvent and geometry.
  //
  AliRunLoader   * rl   = 0x0;
  if(iLoader)
  {
    rl   = AliEveEventManager::AssertRunLoader();
    // runloader check already in AssertRunLoader function 

    //   Int_t evtID = AliEveEventManager::Instance()->GetEventId();
    //   rl->GetEvent(evtID);
  }
  
  AliESDEvent* esd = 0x0;
  if(iESD) esd = AliEveEventManager::Instance()->AssertESD();
  // esd check already in AssertESD function 
      
  //  gGeoManager = gEve->GetDefaultGeometry();
  AliEveEventManager::Instance()->AssertGeometry();
  
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  //TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
  
  //printf("*** nodes %d\n",node->GetNdaughters());
  
  //
  // Initialize the EMCAL data manager
  //
  emcal_data = new AliEveEMCALData(rl,node);//,m);
  
  //printf("*** AliEveEMCALData %p\n",emcal_data);
  
  if(iESD) emcal_data->SetESD(esd);

  //printf("*** AliEveEMCALData set ESD\n");

  //
  // Get the EMCAL information from RunLoader
  //
  if(iLoader)
  {
    //printf("*** Execute Loader methods \n");

    if ( iHits    ) emcal_data->LoadHits(); 

    if ( iDigits  ) emcal_data->LoadDigits();

    if ( iClusters) emcal_data->LoadRecPoints();
  }
  
  //
  // Get the EMCAL information from ESDs
  //
  if(iESD)
  {
    //if(iLoader) rl ->GetEvent(evtNum);

    //printf("*** Execute ESD methods \n");
    
    if(iDigits)   emcal_data->LoadDigitsFromESD();
    
    if(iClusters) emcal_data->LoadRecPointsFromESD();
  }
  
  //printf("*** Data reading executed\n");

  //
  // EVE stuff
  //
  gStyle->SetPalette(1, 0);
  
  gEve->DisableRedraw();
  
  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  l->SetMainColor(Color_t(2));
  gEve->AddElement(l);
  
  //printf("*** Loop SM data, push data \n");

  //
  // Pass the recovered EMCAL data per super-module to EVE
  //
  for (Int_t sm = 0; sm < node->GetNdaughters(); sm++)
  {
    AliEveEMCALSModule* esm = new AliEveEMCALSModule(sm,Form("SM %d Element \n", sm),"EveEMCAL");
    // When/where is this created object cleaned?
    
    esm->SetDataSource(emcal_data);
    
    esm->UpdateQuads(iHits, iDigits, iClusters);
    
    //l->AddElement(esm); // comment, it crashes, replace by:
    
    if ( iDigits   ) gEve->AddElement(esm->GetDigitQuadSet()  , l);
    
    if ( iClusters ) gEve->AddElement(esm->GetClusterQuadSet(), l);
    
    if ( iHits )     gEve->AddElement(esm->GetHitPointSet()   , l);
    
    esm->DropData(); // Not sure it is needed, it works locally with it, but better clean the arrays.    
  }

  //
  // Draw
  //
  gEve->Redraw3D(kTRUE);
  
  gEve->EnableRedraw();
}
