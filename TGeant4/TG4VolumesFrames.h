// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//======================================================
//
//------------TG4VolumesFrames.h--------------------------------//
//---------Frames for the the display of volumes properties---//
//
//=======================================================

#ifndef TG4_VOLUMESFRAMES_H
#define TG4_VOLUMESFRAMES_H

#include <TObject.h>
#include <TGFrame.h>

class TGLabel;
class TGTextBuffer;
class TGTextEntry;
class TGComboBox;


class TG4VolumesFrames : public TObject {

public:   

    TG4VolumesFrames( TGCompositeFrame* parent, TGMainFrame* actionFrame);
    virtual ~TG4VolumesFrames();
    
    void SetVolumesComboEntries();
    void DisplayVolumeCharacteristics();

protected:

    TG4VolumesFrames(const TG4VolumesFrames& vf) ;
    TG4VolumesFrames& operator=(const TG4VolumesFrames& vf) ;

private:

    TGCompositeFrame*   fCapFrame; // the top frame for volumes properties display
    TGCompositeFrame*   fVolSubframe1; // frame for the combo box
    TGCompositeFrame*   fVolSubframe2; //  frame for the text entries        
    TGLayoutHints*      fVolFrameLayout; // layout hints for SubFrames
    TGHorizontalFrame*  fHframe[3];     // horizontal frames for text entries
    TGLabel*            fLabel[3];      // labels for text entries
    TGTextBuffer*       fVolTextBuff[3]; //text buffs for vols propertie
    TGTextEntry*        fVolTextEntry[3]; //text entries for vols properties
    TGComboBox*         fVolumesCombo; // volumes  combo box
    TGLabel*            fComboLabel;   // label for combo box


    void AddLogicalVolumeName( const char* name, Int_t index) const;

    ClassDef(TG4VolumesFrames,0)
         // class for the composition of the volumes display frame    
  };
  
#endif
    
    
     
    
    
 
