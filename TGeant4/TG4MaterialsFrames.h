// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//======================================================
//
//------------TG4MaterialsFrames.h--------------------------------//
//---------Frames for the the display of materials properties---//
//
//=======================================================

#ifndef TG4_MATERIALSFRAMES_H
#define TG4_MATERIALSFRAMES_H

#include <TObject.h>
#include <TGFrame.h>

class TGLabel;
class TGTextBuffer;
class TGTextEntry;
class TGComboBox;

class TG4MaterialsFrames : public TObject {

public:   

    TG4MaterialsFrames(  TGCompositeFrame* Parent, TGMainFrame* ActionFrame);
    virtual ~TG4MaterialsFrames();
    
    void DisplayMaterialCharacteristics();
    void SetMaterialsComboEntries() ;

protected:

    TG4MaterialsFrames(const TG4MaterialsFrames& mf) ;
    TG4MaterialsFrames& operator=(const TG4MaterialsFrames& mf) ;


private:
    TGCompositeFrame*   fCapFrame; // the top frame for materials display
    TGCompositeFrame*   fMatSubframe1; // frame for the combo box
    TGCompositeFrame*   fMatSubframe2; // frame for the text entries
    TGLayoutHints*      fMatFrameLayout; // layout hints for SubFrames



    TGHorizontalFrame*  fHframe[8]; // horizontal frames for text entries
    TGLabel*            fLabel[8];  // labels for text entries 
    TGTextBuffer*       fMatTextBuff[8]; //text buffs for mat properties
    TGTextEntry*        fMatTextEntry[8]; //text entries for mat properties
    TGComboBox*         fMaterialsCombo; // materials  combo box
    TGLabel*            fComboLabel; // label for the combo box
 

    void AddMaterialName( const char* name, Int_t index) const;

    ClassDef(TG4MaterialsFrames,0) 
         // class for the composition of the materials display frame   
  };
  
#endif
    
    
     
    
    
 
