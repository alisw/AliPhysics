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
#include "TG4MaterialsFrames.h"
#include "TG4MainFrame.h"
#include "TG4Editor.h"
#include "TG4Globals.h" 
#include "TG4Limits.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"

#include <TGTextBuffer.h>
#include <TGTextEntry.h>
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGTab.h>
//#include <TGFrame.h>
//#include <TGButton.h>

#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VSolid.hh>
#include <G4UserLimits.hh>
#include <G4Track.hh>
 
 
 ClassImp(TG4VolumesFrames)

TG4VolumesFrames::TG4VolumesFrames( TGTab* Tab, TG4MainFrame* ActionFrame)
  : fPanel(ActionFrame)
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

   fVolumesCombo->Resize(200, 30);
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
     for (Int_t i=0; i<2; i++) {
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
     
// --->a group frame for displaying user's limits with a text entry added

  fGrFrame =
  new TGGroupFrame(fVolSubframe2, "----- User Limits Showcase -----", kVerticalFrame);
  fHframe[2] = new TGHorizontalFrame( fGrFrame, 410, 100, kFixedWidth);
  TGLayoutHints* lLayoutHints6 = new TGLayoutHints(kLHintsTop | kLHintsLeft ,
                           4, 0, 20, 2);
  fLabel[2] = new TGLabel(fHframe[2], labelText[2]);
  fHframe[2]->AddFrame(fLabel[2],lLayoutHints6);
  fVolTextBuff[2] = new TGTextBuffer(200);
  fVolTextEntry[2] = new TGTextEntry(fHframe[2], fVolTextBuff[2], 300);
  fVolTextEntry[2]->Resize(1000, fVolTextEntry[2]->GetDefaultHeight());
  fHframe[2]->AddFrame(fVolTextEntry[2], lLayoutHints6); 
  fVolTextEntry[2]->Associate(ActionFrame);

  fGrFrame->AddFrame(fHframe[2], lLayoutHints4);

// ---> adding to the group frame another frame with text buttons 
// for calling up a full display of user's limits properties

  fGrHFrame= new TGHorizontalFrame( fGrFrame, 500, 100);
  fbtsumm = new TGTextButton(fGrHFrame, "&User Limits", 301);
  fbtcuts = new TGTextButton(fGrHFrame, "&Show Cuts", 302);
  fbtcontrols = new TGTextButton(fGrHFrame, "&Show Controls", 303);
  TGLayoutHints* lLayoutHints7   = 
       new TGLayoutHints(kLHintsTop | kLHintsLeft, 50, 0, 5, 5);
  fGrHFrame->SetLayoutManager(new TGMatrixLayout(fGrHFrame, 0, 3, 15));
  fGrHFrame->AddFrame( fbtsumm );
  fbtsumm->Resize(90, fbtsumm->GetDefaultHeight());
  fGrHFrame->AddFrame( fbtcuts );
  fbtcuts->Resize(90, fbtcuts->GetDefaultHeight());
  fGrHFrame->AddFrame( fbtcontrols );
  fbtcontrols->Resize(90, fbtcontrols->GetDefaultHeight());

  fbtsumm->Associate(ActionFrame);
  fbtcuts->Associate(ActionFrame);
  fbtcontrols->Associate(ActionFrame);
  
  fGrFrame->AddFrame(fGrHFrame, lLayoutHints4);

//---> adding the group frame to the subrame 2
  fVolSubframe2->AddFrame( fGrFrame, lLayoutHints7);
  fGrFrame->Resize(fGrFrame->GetDefaultSize());

// ---> text for the user's limits display window when no volume specified yet 
  fDisplBuffLimits = new TGTextBuffer(1000);
  fDisplBuffLimits->Clear(); 
  fDisplBuffLimits->AddText(0, "\n\n***  No volume specified, "
  "no limits displayed *** ");
  
// ---> text for the user's cuts display window when no volume specified yet 
  fDisplBuffCuts = new TGTextBuffer(1000);
  fDisplBuffCuts->Clear(); 
  fDisplBuffCuts->AddText(0, "\n\n***  No volume specified, "
  "no cuts displayed *** ");
    
// ---> text for the user's controls display window when no volume specified yet 
  fDisplBuffControls = new TGTextBuffer(1000);
  fDisplBuffControls->Clear(); 
  fDisplBuffControls->AddText(0, "\n\n***  No volume specified, "
  "no controls displayed *** ");
      

// ---> making up the Volumes frame   
     fCapFrame->AddFrame(fVolSubframe1,fVolFrameLayout);  
     fCapFrame->AddFrame(fVolSubframe2,fVolFrameLayout);

// --->  going to the main frame     
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
   delete fGrFrame;
   delete fGrHFrame;
   delete fbtsumm;
   delete fbtcuts;
   delete fbtcontrols;
  
   Int_t i;
   for (i=0; i<3; i++) {
     delete fHframe[i];
     delete fVolTextBuff[i];
     delete fVolTextEntry[i];
     delete fLabel[i];
   }
     delete fDisplBuffLimits;
     delete fDisplBuffCuts;
     delete fDisplBuffControls;
     
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
   G4int imat = 0;
   
   G4cout << "\nThe clicked-on volumes entry has the index:  " << index << G4endl;
  
   TG4Limits* lLimits;
   
   if( index < ientr+1 ) {
   
     G4int ii = index-1;
     G4LogicalVolume* lVolume = (*lComboEntries )[ii];
     G4Material* lvMaterial = ((*lComboEntries )[ii])->GetMaterial();
     lLimits = (TG4Limits*)lVolume->GetUserLimits(); 
     TString lDisplayLimits = GetLimitsDisplay(lLimits);
     TString lDisplayCuts = GetCutsDisplay(lLimits);
     TString lDisplayControls = GetControlsDisplay(lLimits);
     
//---> fills up the buffer for popup frame display
     fDisplBuffLimits->Clear(); 
     fDisplBuffLimits->AddText(0,lDisplayLimits);
     
     fDisplBuffCuts->Clear(); 
     fDisplBuffCuts->AddText(0,lDisplayCuts);
     
     fDisplBuffControls->Clear(); 
     fDisplBuffControls->AddText(0,lDisplayControls);
     

     G4cout << lVolume->GetName() << "  " 
            << lVolume->GetSolid()->GetEntityType() << "  "
            << lvMaterial->GetName() << "  "
  	    << lVolume->GetUserLimits()->GetType() <<  G4endl;

//---> putting text in the text entries      
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
       sprintf(buff, lLimits->GetName()); 
    fVolTextBuff[2]->Clear();
    fVolTextBuff[2]->AddText(0, buff);
    gClient->NeedRedraw(fVolTextEntry[2]);
    
    imat = lvMaterial->GetIndex();    
    
   };
   
   if( index == ientr+1 ) {
   
      for ( G4int ii=0; ii<3; ii++) {
        fVolTextBuff[ii]->Clear();          
        gClient->NeedRedraw(fVolTextEntry[ii]);
	};
     fDisplBuffLimits->Clear(); 
     fDisplBuffLimits->AddText(0, "\n\n***  No volume specified, "
     "no limits displayed *** ");
     
     fDisplBuffCuts->Clear(); 
     fDisplBuffCuts->AddText(0, "\n\n***  No volume specified, "
     "no cuts displayed *** ");
     
     fDisplBuffControls->Clear(); 
     fDisplBuffControls->AddText(0, "\n\n***  No volume specified, "
     "no controls displayed *** ");
     
      
     };

//---> setting appropriate  display in the MaterialsFrames    
   TG4MaterialsFrames* mFrames = fPanel->GetMaterialsFrames();
   mFrames->DisplayMaterialCharacteristics( imat + 1 );

}

void TG4VolumesFrames::DisplayUserLimits()
{
//-----> displays User Limits associated with the logical volume 

  const char* cdisplay = fDisplBuffLimits->GetString();
  TG4Editor* ed = new TG4Editor( fCapFrame, 450, 300);
                      ed->LoadBuffer(cdisplay);
		      ed->Popup();      

}

void TG4VolumesFrames::DisplayCuts()
{
//-----> displays Cuts associated with the logical volume 

  const char* cdisplay = fDisplBuffCuts->GetString();
  TG4Editor* ed = new TG4Editor( fCapFrame, 450, 300);
                      ed->LoadBuffer(cdisplay);
		      ed->Popup();      

}

void TG4VolumesFrames::DisplayControls()
{
//-----> displays Controls associated with the logical volume 

  const char* cdisplay = fDisplBuffControls->GetString();
  TG4Editor* ed = new TG4Editor( fCapFrame, 450, 300);
                      ed->LoadBuffer(cdisplay);
		      ed->Popup();      

}

//_______________________________________________________________________
TString TG4VolumesFrames::GetLimitsDisplay(G4UserLimits* limits) const
{
// Returns text for the user limits display in a separate frame
// ---
  G4String display;
  G4Track dummy;
  char buff[200];

  display = "\n\n**************************************";
  display += "\n**************************************\n\n";

  const char* lname = ((TG4Limits*)limits)->GetName();
  sprintf( buff, "\"  %s \"  limits: \n\n  Max step length (mm):   %g ",
            lname, limits->GetMaxAllowedStep(dummy)/mm );
  display += buff;
  
  sprintf( buff, "\n\n  Max track length (mm):  %g  \n\n  Max time (s)         :  %g ",
             limits->GetUserMaxTrackLength(dummy)/mm, limits->GetUserMaxTime(dummy)/s);             
  display += buff;
   
  sprintf( buff, "\n\n  Min kin. energy (MeV):  %g \n\n  Min range (mm):         %g ",
            limits->GetUserMinEkine(dummy)/MeV, limits->GetUserMinRange(dummy)/mm );
  display += buff;
  
  display += "\n\n**************************************";
  display += "\n**************************************\n\n";	     

  const char* tmp = display;
  
  return TString(tmp);

}
//===================inserted================================================
//_______________________________________________________________________
TString TG4VolumesFrames::GetCutsDisplay(G4UserLimits* limits) const
{
// Returns text for the cuts display in a separate frame
// ---
  G4String display;

  const TG4G3CutVector* lCutVector = ((TG4Limits*)limits) -> GetCutVector();

  if ( ((TG4Limits*)limits) -> IsCut() )
     display = lCutVector -> Format(); 
  else
     display =  "\n\n****  No special cuts.  **** "; 

  const char* tmp = display;
  
  return TString(tmp);
  
}
//_________________________________________________________________________
TString TG4VolumesFrames::GetControlsDisplay(G4UserLimits* limits) const
{
// Returns text for the controls display in a separate frame
 
 G4String display;

 const TG4G3ControlVector* lControlVector = ((TG4Limits*)limits) -> GetControlVector();
 
 if ( ((TG4Limits*)limits) -> IsControl() )
    display = lControlVector -> Format();
 else 
    display =  "\n\n****  No special controls.  **** "; 

  const char* tmp = display;
  
  return TString(tmp);

}

//================end of inserted==================================================  
