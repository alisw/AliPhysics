// $Id$
// Category: interfaces
//
// Author: D. Adamova
//========================================================
//
//---------------TG4Editor.cxx---------------------------//
//------- A supplementary service class for--------------//
//-----------AG4 Geometry Browser------------------------//
//
//=========================================================



#ifndef TG4_EDITOR_H
#define TG4_EDITOR_H

#include <TGFrame.h> 

class TGTextEdit;
class TGTextButton;
class TGLayoutHints;

class TG4Editor : public TGTransientFrame {
 
public:
   TG4Editor(const TGWindow* main, UInt_t w, UInt_t h);
   virtual ~TG4Editor();

   void   LoadBuffer(const char* buffer);
   void   SetTitle();
   void   Popup();
   void   CloseWindow();
   Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
   
//--->inline
   TGTextEdit* GetEditor() const { return fEdit; }   

protected:
   TG4Editor(const TG4Editor& ge); 
   TG4Editor& operator=(const TG4Editor& ge) ;

   
private:
   TGTextEdit*       fEdit;   // text edit widget
   TGTextButton*     fOK;     // OK button
   TGLayoutHints*    fL1;     // layout of TGTextEdit
   TGLayoutHints*    fL2;     // layout of OK button
   
   ClassDef(TG4Editor,0)   // service Editor window for GUI
};

#endif
