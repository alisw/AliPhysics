// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//========================================================
//
//------------TG4VolumesFrames.cxx--------------------------------//
//--------- Frames for the the display of volumes properties---//
//
//========================================================= 
 
#include "TG4VolumesFrames.h"
#include "TG4Globals.h" 

#include <TGTextBuffer.h>
#include <TGTextEntry.h>
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGTab.h>
#include <TGFrame.h>

#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VSolid.hh>
 
 
 ClassImp(TG4VolumesFrames)

TG4VolumesFrames::TG4VolumesFrames( TGTab* Tab, TGMainFrame* ActionFrame)
{ 
//---> creates the volumes properties display frame
//---> and plunges it into the main frame
   TGCompositeFrame* parent = Tab->AddTab("Volumes Properties");
   fCapFrame = new TGCompositeFrame(parent, 60, 20, kHorizontalFrame);
   ULong_t back= TGFrame::GetBlackPixel(); 
   fCapFrame->ChangeBackground(back);
   fVolSubframe1 = new TGCompositeFrame(fCapFrame, 60, 20, kVerticalFrame);
   fVolFrameLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

// ComboBox for lvolumes
   fVolumesCombo = new TGComboBox(fVolSubframe1, 100);
   TGLayoutHints* lLayoutHints3 = 
             new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);   
   Text_t* lComboLabelText = " Pick up a volume here ";
   fComboLabel = new TGLabel( fVolSubframe1, lComboLabelText);
   fVolSubframe1->AddFrame(fComboLabel, lLayoutHints3);
   fVolSubframe1->AddFrame(fVolumesCombo, fVolFrameLayout);

   fVolumesCombo->Resize(200, 20);
   fVolumesCombo->Associate(ActionFrame); 

   
// text labels with lvolumes properties 

   Text_t* labelText[3]  = 
   {"Shape's Name", 
    "Material    ",
    "User Limits " }; 

// Entries for lvolumes properties
   TGLayoutHints* lLayoutHints4   = 
       new TGLayoutHints(kLHintsTop | kLHintsExpandY, 5, 5, 5, 5);
   TGLayoutHints* lLayoutHints5 = 
       new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
   fVolSubframe2 = new TGCompositeFrame(fCapFrame, 60, 20, kVerticalFrame);

   { // local scope for i
     for (Int_t i=0; i<3; i++) {
       Int_t idT=i+1;   
       fHframe[i] = new TGHorizontalFrame(fVolSubframe2, 500, 100, kFixedWidth);
       fVolSubframe2->AddFrame(fHframe[i], lLayoutHints4);
       fVolTextBuff[i] = new TGTextBuffer(200);
       fVolTextEntry[i] = new TGTextEntry(fHframe[i], fVolTextBuff[i], 300);
       fLabel[i] = new TGLabel(fHframe[i], labelText[i]);
       fHframe[i]->AddFrame(fLabel[i], lLayoutHints5);
       fHframe[i]->AddFrame(fVolTextEntry[i], lLayoutHints5);
       fVolTextEntry[i]->Associate(ActionFrame); 
     }
   } 

// making up the Volumes frame   
     fCapFrame->AddFrame(fVolSubframe1,fVolFrameLayout);  
     fCapFrame->AddFrame(fVolSubframe2,fVolFrameLayout);
// going to the main frame     
     parent->AddFrame(fCapFrame, fVolFrameLayout);
}

TG4VolumesFrames::TG4VolumesFrames(const TG4VolumesFrames& vf) 
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4VolumesFrames copy constructor.");
}

TG4VolumesFrames& TG4VolumesFrames::operator=(const TG4VolumesFrames& vf)
{
  // check assignement to self
  if (this == &vf) return *this;

  TG4Globals::Exception(
    "Attempt to assign  singleton.");
    
  return *this;  
}    

TG4VolumesFrames::~TG4VolumesFrames()
{
//---> liquidator 

   G4cout << "\n Now in  TG4VolumesFrames destructor \n" << G4endl;
   delete fVolSubframe1;
   delete fVolFrameLayout;
   delete fVolumesCombo;
   delete fComboLabel;
   delete fVolSubframe2;
   delete fCapFrame;
  
   Int_t i;
   for (i=0; i<3; i++) {
     delete fHframe[i];
     delete fVolTextBuff[i];
     delete fVolTextEntry[i];
     delete fLabel[i];
   }

}

void TG4VolumesFrames::SetVolumesComboEntries() 
{
//--->//---> puts names of lvolumes into the combo box entries

    G4LogicalVolumeStore* lComboEntries = G4LogicalVolumeStore::GetInstance();

    G4int ig = lComboEntries->size();
    G4String name;
    
    for (int ii=0; ii < ig; ii++)
        { name = ((*lComboEntries )[ii])->GetName() ;
	  AddLogicalVolumeName( name, ii+1);
	};
	
    name = "  " ;
    AddLogicalVolumeName( name, ig+1);
    	
    
}

void TG4VolumesFrames::AddLogicalVolumeName( const char* name, Int_t index) const
{
//-----> adds an lvolume name to the combo box  

   fVolumesCombo->AddEntry( name, index);
   fVolumesCombo->Select(index);
   fVolumesCombo->Resize(200, 20);
}

void TG4VolumesFrames::DisplayVolumeCharacteristics()
{
//-----> shows informations about a logical volume 

   G4LogicalVolumeStore* lComboEntries = G4LogicalVolumeStore::GetInstance();
   G4int ientr = lComboEntries->size();
   G4int index = fVolumesCombo->GetSelected();
   
   G4cout << "\nThe clicked-on volumes entry has the index:  " << index << G4endl;

   if( index < ientr+1 ) {
   
     G4int ii = index-1;
     G4LogicalVolume* lVolume = (*lComboEntries )[ii];
     G4Material* lvMaterial = ((*lComboEntries )[ii])->GetMaterial();

     G4cout << lVolume->GetName() << "  " 
            << lVolume->GetSolid()->GetEntityType() << "  "
            << lvMaterial->GetName() << "  "
	    << lVolume->GetUserLimits() << "  "
	    << G4endl;
	  
	  
   char buff[100];
   
    sprintf(buff, lVolume->GetSolid()->GetEntityType());
    fVolTextBuff[0]->Clear();
    fVolTextBuff[0]->AddText(0, buff);
    gClient->NeedRedraw(fVolTextEntry[0]);
    
    sprintf(buff, lvMaterial->GetName());
    fVolTextBuff[1]->Clear();
    fVolTextBuff[1]->AddText(0, buff);
    gClient->NeedRedraw(fVolTextEntry[1]);

    
    sprintf(buff, "User limits undefined" );
    if (lVolume->GetUserLimits())
       sprintf(buff, "User limits defined" ); 
    fVolTextBuff[2]->Clear();
    fVolTextBuff[2]->AddText(0, buff);
    gClient->NeedRedraw(fVolTextEntry[2]);
   };
   
   if( index == ientr+1 ) {
   
      for ( G4int ii=0; ii<3; ii++) {
        fVolTextBuff[ii]->Clear();          
        gClient->NeedRedraw(fVolTextEntry[ii]);
	};
     };
}

