// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//========================================================
//
//------------TG4MaterialsFrames.cxx--------------------------------//
//--------- Frames for the the display of materials properties---//
//
//========================================================= 
 
#include "TG4MaterialsFrames.h" 
#include "TG4Globals.h"

#include <TGTextBuffer.h>
#include <TGTextEntry.h>
#include <TGComboBox.h>
#include <TGLabel.h>
 
#include <G4Material.hh>
#include <G4Element.hh>



 ClassImp(TG4MaterialsFrames)

TG4MaterialsFrames::TG4MaterialsFrames( TGCompositeFrame* Parent, TGMainFrame* ActionFrame )
{ 
//---> creates the materials properties display frame
//---> and plunges it into the main frame
   fCapFrame = new TGCompositeFrame(Parent, 60, 20, kHorizontalFrame);
   fMatSubframe1 = new TGCompositeFrame(fCapFrame, 60, 20, kVerticalFrame);
   fMatFrameLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

   // ComboBox for materials
   fMaterialsCombo = new TGComboBox(fMatSubframe1, 200);
   TGLayoutHints* lLayoutHints3 = 
             new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);   
   Text_t* lComboLabelText = " Pick up a material here ";
   fComboLabel = new TGLabel( fMatSubframe1, lComboLabelText);
   fMatSubframe1->AddFrame(fComboLabel, lLayoutHints3);
   fMatSubframe1->AddFrame(fMaterialsCombo, fMatFrameLayout);


   fMaterialsCombo->Resize(200, 20);
   fMaterialsCombo->Associate(ActionFrame); 


   
// text labels with material properties 
   Text_t* labelText[8]  = 
   {"Index             ", 
    "Number of elements",
    "Elements list     ",
    "Atomic mass       ",
    "Density           ", 
    "State             ", 
    "Radiation Length  ",
    "Abs. Length       " }; 

// Entries for material properties
   TGLayoutHints* lLayoutHints4   = 
       new TGLayoutHints(kLHintsTop | kLHintsExpandY, 5, 5, 5, 5);
   TGLayoutHints* lLayoutHints5 = 
       new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
   fMatSubframe2 = new TGCompositeFrame(fCapFrame, 60, 20, kVerticalFrame);
                
   { // local scope for i
     for (Int_t i=0; i<8; i++) {
       Int_t idT=i+1;   
       fHframe[i] = new TGHorizontalFrame(fMatSubframe2, 500, 100, kFixedWidth);
       fMatTextBuff[i] = new TGTextBuffer(200);
       fMatTextEntry[i] = new TGTextEntry(fHframe[i], fMatTextBuff[i], 300);
       fLabel[i] = new TGLabel(fHframe[i], labelText[i]);
       fHframe[i]->AddFrame(fLabel[i], lLayoutHints5);
       fHframe[i]->AddFrame(fMatTextEntry[i], lLayoutHints5);
       fMatSubframe2->AddFrame(fHframe[i], lLayoutHints4);
       fMatTextEntry[i]->Associate(ActionFrame);
     }
   } 

// making up Materials frame   
     fCapFrame->AddFrame(fMatSubframe1,fMatFrameLayout);  
     fCapFrame->AddFrame(fMatSubframe2,fMatFrameLayout);

// going to the main frame     
     Parent->AddFrame(fCapFrame, fMatFrameLayout);   
   
}

TG4MaterialsFrames::TG4MaterialsFrames(const TG4MaterialsFrames& mf) 
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4MaterialsFrames copy constructor.");
}

TG4MaterialsFrames& TG4MaterialsFrames::operator=(const TG4MaterialsFrames& mf)
{
  // check assignement to self
  if (this == &mf) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4MaterialsFrames singleton.");
    
  return *this;  
}    

TG4MaterialsFrames::~TG4MaterialsFrames()
{
  //---> liquidator
  
   G4cout << "\n Now in  TG4MaterialsFrames destructor \n"<< G4endl;
   delete fMatSubframe1;
   delete fMatFrameLayout;
   delete fMaterialsCombo;
   delete fComboLabel;
   delete fMatSubframe2;
   delete fCapFrame;

   Int_t i;
   for (i=0; i<8; i++) {
     delete fHframe[i];
     delete fMatTextBuff[i];
     delete fMatTextEntry[i];
     delete fLabel[i];
   }

}

void TG4MaterialsFrames::SetMaterialsComboEntries()
{
//---> puts names of materials into the combo box entries

   const G4MaterialTable* lComboEntries = G4Material::GetMaterialTable();

    G4int ig = lComboEntries->entries();
    G4String name;
    
    for (int ii=0; ii < ig; ii++)
        { name = ((*lComboEntries )[ii])->GetName() ;
          AddMaterialName( name, ii+1);
	};

    name = "  " ;
    AddMaterialName( name, ig+1);
    		   
}

void TG4MaterialsFrames::AddMaterialName( const char* name, Int_t index) const
{
  
//-----> adds a material name to the combo box  

   fMaterialsCombo->AddEntry( name, index);
   fMaterialsCombo->Select(index);
   fMaterialsCombo->Resize(200, 20);
}

void TG4MaterialsFrames::DisplayMaterialCharacteristics()
{
  
//-----> shows informations about materials listed in G4MaterialTable 

   const G4MaterialTable* lComboEntries = G4Material::GetMaterialTable();
   G4int ientr = lComboEntries->entries();
   G4int index = fMaterialsCombo->GetSelected();
   
   G4cout << "\nThe clicked-on material has the index:  " << index << G4endl;
   
   if( index < ientr+1 ) {
   
     G4int ii = index-1;
     G4Material* lvMaterial = (*lComboEntries )[ii];
     const G4ElementVector* allElements = lvMaterial->GetElementVector();

     G4cout << lvMaterial->GetName() << "  "
	    << lvMaterial->GetNumberOfElements() << "  "
	    << (*allElements )[0]->GetName() << "...  "
	    << lvMaterial->GetDensity() << "  "
	    << lvMaterial->GetState() << "  "
	    << lvMaterial->GetRadlen() << "  "
	    << G4endl;
	  
   char buff[200];
   
    sprintf(buff, "%10i",index );
    fMatTextBuff[0]->Clear();
    fMatTextBuff[0]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[0]);
    
    G4int noe = lvMaterial->GetNumberOfElements();
    sprintf(buff, "%10i", noe);
    fMatTextBuff[1]->Clear();
    fMatTextBuff[1]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[1]);

    G4String stringOfElements = "  ";
    for (G4int ie=0; ie < noe; ie++) {
    stringOfElements += (*allElements )[ie]->GetName();
    stringOfElements += "  ";
    };
    sprintf(buff, stringOfElements);
    fMatTextBuff[2]->Clear();
    fMatTextBuff[2]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[2]);

    sprintf(buff, " Multi element material" );
    if( noe < 2 )
        sprintf(buff, "%10.2e", lvMaterial->GetA() );
    fMatTextBuff[3]->Clear();
    fMatTextBuff[3]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[3]);
    
    sprintf(buff, "%10.2e", lvMaterial->GetDensity() );
    fMatTextBuff[4]->Clear();
    fMatTextBuff[4]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[4]);
    
    sprintf(buff,"%10i", lvMaterial->GetState());
    fMatTextBuff[5]->Clear();
    fMatTextBuff[5]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[5]);
    
    sprintf(buff,"%10.2e", lvMaterial->GetRadlen());
    fMatTextBuff[6]->Clear();
    fMatTextBuff[6]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[6]);
    
    sprintf(buff, "           " );
    fMatTextBuff[7]->Clear();
    fMatTextBuff[7]->AddText(0, buff);
    gClient->NeedRedraw(fMatTextEntry[7]);
    
   };
   
   if( index == ientr+1 ) {
   
      for ( G4int ii=0; ii<8; ii++) {
        fMatTextBuff[ii]->Clear();          
        gClient->NeedRedraw(fMatTextEntry[ii]);
	};
     };      
   
}

