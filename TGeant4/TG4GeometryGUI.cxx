// $Id$
// Category: interfaces
//
// Author: D. Adamova
//==============================================================
//
//----------------TG4GeometryGUI.cxx--------------------------//
//----------------AG4 Geometry Browser----------------------//
//
//=================================================================

  
		 
#include "TG4GeometryGUI.h"
#include "TG4GuiVolume.h"
#include "TG4MainFrame.h"
#include "TG4ListTreeFrame.h"
#include "TG4VolumesFrames.h"
#include "TG4MaterialsFrames.h"
#include "TG4Globals.h"

#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <TGListTree.h>
#include <TObjArray.h>


		 
ClassImp(TG4GeometryGUI)    

TG4GeometryGUI::TG4GeometryGUI()
{
//---> Constructor		 
    fPanel   =   new TG4MainFrame(gClient->GetRoot(), 650, 500);

    G4cout << "\n***********************************************" << G4endl;
    G4cout << "***********************************************" << G4endl;
    G4cout << " Welcome to the Geometry Browser for AliGeant4 " << G4endl;
    G4cout << "\n***********************************************" << G4endl;
    G4cout << "***********************************************\n" << G4endl;
    
    ReadGeometryTree();
    
    ReadMaterials();  

 }
 
TG4GeometryGUI::TG4GeometryGUI(const TG4GeometryGUI& gg) 
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4GeometryGUI copy constructor.");
}

TG4GeometryGUI& TG4GeometryGUI::operator=(const TG4GeometryGUI& gg)
{
  // check assignement to self
  if (this == &gg) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4GeometryGUI singleton.");
    
  return *this;  
}    

TG4GeometryGUI::~TG4GeometryGUI(){
//---> liquidator

  G4cout << "\n Now in TG4GeometryGUI destructor \n" << G4endl;
  delete fPanel;
// ----> guiVolumes not taken care of yet

}

void TG4GeometryGUI::ReadGeometryTree()
{
//--->Linking logical volumes to gui volumes

//  Icons for folders
    const TGPicture* kFolder     = gClient->GetPicture("folder_t.xpm");
    const TGPicture* kOpenFolder = gClient->GetPicture("ofolder_t.xpm");
    const TGPicture* kDocument   = gClient->GetPicture("doc_t.xpm");
    TG4GuiVolume*  volume;

    G4LogicalVolumeStore* volumeStore;
    G4LogicalVolume* top;
         
    TGListTreeItem* itemi;
    
    volumeStore = G4LogicalVolumeStore::GetInstance();
    top = (*volumeStore )[0];  

    G4String vname = top->GetName();
    volume = new TG4GuiVolume( vname, top);//--->TObject
    itemi = fPanel->GetListTreeFrame()
      ->AddItem(volume,0,vname, kOpenFolder, kFolder);
 
    RegisterLogicalVolume( top, itemi);
    
//    delete volume; ---> inactivated to get UserData  for the
//                   ---> ListTree items in the MainFrame
 
}	   

void TG4GeometryGUI::RegisterLogicalVolume(G4LogicalVolume* lv,
                                           TGListTreeItem* itemv) 
{
//--->Filling  up gui volumes objArray  

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
      itemi = fPanel->GetListTreeFrame()
        ->AddItem(volume, itemv, vname, kOpenFolder, kFolder);
      
      itemi->SetUserData(volume);
      
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
  };
  
  delete guiVolumes;
//  delete volume;---> if deleted, wrong logical volumes assigned 
//                ---> to the ListTree items in the MainFrame display
}  

void TG4GeometryGUI::ReadMaterials() const
{
//-----> Puts logical volumes and materials names 
//-----> into ComboBoxes 
    TG4VolumesFrames* vFrame = fPanel->GetVolumesFrames();
    vFrame->SetVolumesComboEntries();

    TG4MaterialsFrames* mFrame = fPanel->GetMaterialsFrames();
    mFrame->SetMaterialsComboEntries();

}

