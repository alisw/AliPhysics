// $Id$
// Category: interfaces
//
// Author: D. Adamova, I. Hrivnacova
//==============================================================
//
//----------------TG4GeometryGUI.cxx--------------------------//
//----------------AG4 Geometry Browser----------------------//
//
//=================================================================

  
		 
#include "TG4GeometryGUI.h"
#include "TG4GuiVolume.h"
#include "TG4GUI.h"
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <TGListTree.h>
#include <TObjArray.h>
		 
ClassImp(TG4GeometryGUI)    

TG4GeometryGUI::TG4GeometryGUI()
{
 // Constructor		 
    fPanel   =   new TG4GUI(gClient->GetRoot(), 650, 500);
 
    ReadGeometryTree();
    
}

void TG4GeometryGUI::ReadGeometryTree()
{
// Linking logical volumes to gui volumes

//  Icons for folders
    const TGPicture* kFolder     = gClient->GetPicture("folder_t.xpm");
    const TGPicture* kOpenFolder = gClient->GetPicture("ofolder_t.xpm");
 
    TG4GuiVolume*  volume;
    G4LogicalVolumeStore* volumeStore;
    G4LogicalVolume* top;
         
    TGListTreeItem* itemi;  
    
    volumeStore = G4LogicalVolumeStore::GetInstance();
    top = (*volumeStore )[0];  

    G4String vname = top->GetName();
    volume = new TG4GuiVolume( vname, top);//--->TObject
    itemi = fPanel->AddItem(volume,0,vname, kOpenFolder, kFolder);
 
    RegisterLogicalVolume( top, itemi);
 
}	   

void TG4GeometryGUI::RegisterLogicalVolume(G4LogicalVolume* lv,
                                           TGListTreeItem* itemv) 
{
// Filling  up gui volumes objArray  

typedef G4std::set <G4String, G4std::less<G4String> > TG4StringSet;
TG4StringSet     lVolumeNames;     //set of names of solids  

//  Icons for folders
    const TGPicture* kFolder     = gClient->GetPicture("folder_t.xpm");
    const TGPicture* kOpenFolder = gClient->GetPicture("ofolder_t.xpm");
 
  TG4GuiVolume*  volume;
  G4LogicalVolume* lDaughter;
  G4VPhysicalVolume* pDaughter;
  TGListTreeItem* itemi;
 
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters == 0) return;

//---------->only for nofDaughters > 0  
  // open composition
  G4String lvName = lv->GetName();
  TObjArray* guiVolumes = new TObjArray();
  
 
  G4int ii;
  for (ii=0; ii<nofDaughters; ii++) {
 
    pDaughter = lv->GetDaughter(ii);
    lDaughter = pDaughter->GetLogicalVolume();
    G4String vname = lDaughter->GetName();

//-------> process lv only if it was not yet processed
    if ( (lVolumeNames.find(vname)) == (lVolumeNames.end()) ) {

      volume = new TG4GuiVolume( vname, lDaughter);
      itemi = fPanel->AddItem(volume, itemv, vname, kOpenFolder, kFolder);
      volume->SetItem( itemi);

//-----> store the name of logical volume in the set
      lVolumeNames.insert(lVolumeNames.begin(),vname); 
      guiVolumes->AddLast( volume );
   };
 
 }; 
  

//-----> process daughters
    for (ii=0; ii<guiVolumes->GetEntriesFast(); ii++) {

      TG4GuiVolume* guiVolume = (TG4GuiVolume*)(guiVolumes->At(ii));
      G4LogicalVolume* lvd = guiVolume->GetLogicalVolume();
      TGListTreeItem* itemvd = guiVolume->GetItem();
      RegisterLogicalVolume(lvd, itemvd);
  }
  
  delete guiVolumes;   
}  

