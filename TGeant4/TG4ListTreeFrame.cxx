// $Id$
// Category: interfaces
//
// Author: D. Adamova
//========================================================
//
//------------TG4ListTreeFrame.cxx--------------------------------//
//--------- Frame for the ListTree container---//
//
//========================================================= 
 
#include "TG4ListTreeFrame.h"
#include "TG4GuiVolume.h"
#include "TG4Globals.h"

#include <TGListTree.h>
#include <TGCanvas.h>

#include <G4LogicalVolume.hh>
#include <G4UImanager.hh>

ClassImp(TG4ListTreeFrame)

TG4ListTreeFrame::TG4ListTreeFrame( TGCompositeFrame* Parent, TGMainFrame* ActionFrame)
{
//------>canvas for the ListTree
   fCanvasWindow = new TGCanvas( Parent, 400, 240);   
   ULong_t back= TGFrame::GetWhitePixel(); 
   fCanvasWindow->ChangeBackground(back);
 
//----->ListTree for the logical volumes
   fVolumesListTree= new TGListTree(fCanvasWindow->GetViewPort(), 10, 10, kHorizontalFrame);
   fVolumesListTree->Associate(ActionFrame);
//------->container for the ListTree
   fCanvasWindow->SetContainer(fVolumesListTree);

   TGLayoutHints *lCanvasLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   Parent->AddFrame(fCanvasWindow, lCanvasLayout);
   
}

TG4ListTreeFrame::TG4ListTreeFrame(const TG4ListTreeFrame& ltf) 
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4ListTreeFrame copy constructor.");
}

TG4ListTreeFrame& TG4ListTreeFrame::operator=(const TG4ListTreeFrame& ltf)
{
  // check assignement to self
  if (this == &ltf) return *this;

  TG4Globals::Exception(
    "Attempt to assign  singleton.");
    
  return *this;  
}    

TG4ListTreeFrame::~TG4ListTreeFrame()
{
//---> liquidator 

   G4cout << "\n Now in  TG4ListTreeFrame destructor \n" << G4endl;
   delete fCanvasWindow;
   delete fVolumesListTree;   

}   


TGListTree* TG4ListTreeFrame::GetVolumesListTree() const
{
//---> get VolumesListTree for use in the MainFrame
  return fVolumesListTree;
}

void  TG4ListTreeFrame::DrawSelectedVolume(TGListTreeItem* item)
{
//---> draw the volume belonging to the item

  TG4GuiVolume* volume=((TG4GuiVolume*) item->GetUserData());
  G4LogicalVolume* lvolume = volume->GetLogicalVolume();	     
//-------->>>>>>inserted part for getting parent and number of daughters
  TGListTreeItem* iParent = item -> GetParent();
  G4int noDght, ii;
  G4int icopy = 0;
  G4String parentName, lName;
    
  if (iParent) {
    
      TG4GuiVolume* vParent   = ((TG4GuiVolume*) iParent->GetUserData());
      G4LogicalVolume* lvParent = vParent->GetLogicalVolume();
      parentName = lvParent->GetName();
      noDght = lvParent->GetNoDaughters();
//------>>>> end of getting parent and number of daughters

   //------->>>>inserted for getting no on copies
      ii = 0;
      while ( ii < noDght){
         if ( lvParent->GetDaughter(ii)->GetLogicalVolume() ==  lvolume )
              icopy++;  
     
             ii++;
          };
   //------>>>> end of getting no of copies

      lName = lvolume->GetName();
  }
//------->>>>> case of iParent = 0        
  else {
      lName = "ALIC";
      icopy = 1;
  };
     
  G4cout << "For logical volume   " << lName << "   "
         << "the number of copies is :!   "  << icopy << endl;
//------> looping over number of copies
   ii=0;	      
   while ( ii < icopy) { 
      G4String lCommand = "/vis/scene/add/volume  " ;
      lCommand += lName;
      lCommand += "  ";
      TG4Globals::AppendNumberToString( lCommand, ii);
      lCommand += "  -1";
      G4cout << "!!!!ADD VOLUME COMMAND IS!!!!   " << lCommand << G4endl;
      G4UImanager::GetUIpointer()->ApplyCommand( lCommand );
      G4UImanager::GetUIpointer()->ApplyCommand(
	          "/control/execute vis_cont.mac");
      ii++;
      };
      
}
	      
 
 
