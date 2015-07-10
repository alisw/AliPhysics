
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//************************************************************************
///
/// \file emcal_esdcells.C
/// \brief Visualize EMCAL ESD cells
///
/// A macro to read and visualize EMCAL digits. 
/// Standalone, it does not used the goodies of the classes AliEveEMCALXXX.
/// It could be used as a simple testing tool for further development in the classes.
///
/// Include it in the macro visscan_init.C in this way:
/// exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "EMCAL ESD CELLS", "emcal_esdcells.C", "emcal_esdcells","",kTRUE));
/// (the last parameter of the visscan_init macro indicates that this line is active or not).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. DCal implementation + doxygen, May 2015.
//************************************************************************


#ifndef __CINT__

#include <TEveManager.h>
#include <TEveBoxSet.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TEveTrans.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>

#include <EveBase/AliEveEventManager.h>

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALDigit.h>
#include <AliLog.h>

#endif

Bool_t fMatrixEMSet = kFALSE;

void emcal_esdcellsV3()
{   
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();

  AliEveEventManager::GetMaster()->AssertGeometry();
  
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  if (!node) return;
  
  Int_t nModules = node->GetNdaughters();

  // Get the EMCAL geometry
  //
  AliEMCALGeometry * geom  = AliEMCALGeometry::GetInstance();  
  if (!geom) 
  {
    printf("xxx Set default geo as Run2 xxx\n");
    geom  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  }
    
  // Set the geometry alignment matrices from the ones stored in the data.
  if(!fMatrixEMSet)
  {
    printf("LOAD EMCAL MATRICES\n");
    
    for(Int_t mod = 0; mod < (geom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
    { 
      printf("Load EMCAL ESD matrix %d, %p\n",mod,esd->GetEMCALMatrix(mod));
      
      if( esd->GetEMCALMatrix(mod) ) 
        geom->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
      else // set default identity matrix
        geom->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over super modules	     
    
    fMatrixEMSet = kTRUE;
  }
  
  // Check that the EMCAL geo and the nodes from EMCAL have the same number of entries
  //
  if(nModules != geom->GetNumberOfSuperModules())
    printf("*** === EMCAL DIGITS - N Daughter Nodes %d - N super mod %d === ***\n", 
           node->GetNdaughters(), geom->GetNumberOfSuperModules());
    
  // Get the EMCAL bounding boxes for the super modules.
  // 4 kind of SM: 10 Full EMCal, 2 1/3 EMCal, 6 DCal (2/3 EMCal) and 2 1/3 EMCal in DCal region.
  //
  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
  TEveFrameBox* frame_big = new TEveFrameBox();
  frame_big->SetFrameColorRGBA(200,200,0,50);
  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());
  
  TEveFrameBox* frame_sml  = 0x0;
  TEveFrameBox* frame_dcl  = 0x0;
  TEveFrameBox* frame_smld = 0x0;

  if (nModules > 10) 
  {
    TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();
    frame_sml = new TEveFrameBox();
    frame_sml->SetFrameColorRGBA(200,200,0,50);
    frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());
  }

  if (nModules > 12) 
  {
    TGeoBBox* dbbox = (TGeoBBox*) node->GetDaughter(12)->GetVolume()->GetShape();
    frame_dcl = new TEveFrameBox();
    frame_dcl->SetFrameColorRGBA(200,200,0,50);
    frame_dcl->SetAABoxCenterHalfSize(0, 0, 0, dbbox->GetDX(), dbbox->GetDY(), dbbox->GetDZ());
    
    TGeoBBox* sdbbox = (TGeoBBox*) node->GetDaughter(18)->GetVolume()->GetShape();
    frame_smld = new TEveFrameBox();
    frame_smld->SetFrameColorRGBA(200,200,0,50);
    frame_smld->SetAABoxCenterHalfSize(0, 0, 0, sdbbox->GetDX(), sdbbox->GetDY(), sdbbox->GetDZ());
  }
  
  // Define EVE stuff
  //
  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  gEve->AddElement(l);

  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 1024);

  // Here we will store the EMCAL data that will be treated by EVE
  // per each super-module.
  // Pass the SM bounding boxes (frames).
  //
  const Int_t nSM = nModules;
  TEveBoxSet* smodules[nSM];
  memset(smodules,0,nModules*sizeof(TEveBoxSet*));
  
  for (Int_t sm = 0; sm < nModules; ++sm)
  {
    TEveBoxSet* q = new TEveBoxSet(Form("SM %d", sm+1));
    q->SetOwnIds(kTRUE);
    
    // Type of object to be displayed, rectangle with cell size
    //q->Reset(TEveBoxSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    //q->SetDefWidth (geom->GetPhiTileSize());
    //q->SetDefHeight(geom->GetEtaTileSize());

    q->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);

    q->RefMainTrans().SetFrom(*node->GetDaughter(sm)->GetMatrix());

//    if     (sm < 10) q->SetFrame(frame_big );
//    else if(sm < 12) q->SetFrame(frame_sml );
//    else if(sm < 18) q->SetFrame(frame_dcl );
//    else if(sm < 20) q->SetFrame(frame_smld);
    
    q->SetPalette(pal);

    gEve->AddElement(q, l);
    smodules[sm] = q;
  }

  // EMCAL data reading
    
  if(!esd)
  {
    printf("emcal_esdcells: ESD event not available\n");
    return;
  }
  
  AliESDCaloCells &cells= *(esd->GetEMCALCells());
  
  if(!esd->GetEMCALCells()) 
  {
    printf("emcal_esdcells: ESDCaloCells not available\n");
    return;
  }
  
  Int_t ncell = cells.GetNumberOfCells() ;  
  
  printf("Number of ESD CaloCells %d\n",ncell);
  
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  Float_t amp   = -1 ;
  Float_t time  = -1 ;
  Int_t id      = -1 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;
  Int_t ieta    =  0 ;

  
  // Extract digit information from the ESDs
  for (Int_t icell=  0; icell <  ncell; icell++) 
  {
    id  = cells.GetCellNumber(icell);
    amp = cells.GetAmplitude (icell); // GeV
        
     //Geometry methods
    geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
    //Gives SuperModule and Tower numbers
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    //Gives label of cell in eta-phi position per each supermodule
    
    
    // It should not happen, but in case the OCDB file is not the
    // correct one.
    if(iSupMod >= nModules) continue;
    
    
    //printf("CaloCell %d, ID %d, energy %2.3f, SM %d/%d, icol %d, irow %d\n",
    //       icell,id,amp, iSupMod, nModules,iphi,ieta);
    
    //
    geom->RelPosCellInSModule(id, x, y, z);
    
    // Push the data to the visualization tools
    TEveBoxSet* q = smodules[iSupMod];
    if (q) 
    {
      //q->AddQuad(y, z);
      //q->QuadValue(amp);
          
      q->AddBox(15, y,  z, amp, 6.0, 6.0);
            
      q->DigitValue(static_cast<Int_t>(amp));
    }

  }
    
  // Send the data to EVE?
  for (Int_t sm = 0; sm < nModules; ++sm)
  {
    smodules[sm]->RefitPlex();
  }
  
  gEve->Redraw3D(kTRUE);
}
