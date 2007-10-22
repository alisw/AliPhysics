#ifndef ROOT_TGCommandPlugin
#define ROOT_TGCommandPlugin

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

class TGLabel;
class TGComboBox;
class TGTextEntry;
class TGTextBuffer;
class TGTextView;

class TGCommandPlugin : public TGMainFrame {

protected:
   Int_t              fPid;               // current process id
   TGHorizontalFrame *fHf;
   TGLabel           *fLabel;             // "command :" label
   TGComboBox        *fComboCmd;          // commands combobox
   TGTextEntry       *fCommand;           // command text entry widget
   TGTextBuffer      *fCommandBuf;        // command text buffer
   TGTextView        *fStatus;

public:

   TGCommandPlugin(const TGWindow *p, UInt_t w, UInt_t h);
   virtual ~TGCommandPlugin();

   void     CheckRemote(const char * /*str*/);
   void     HandleCommand();

   ClassDef(TGCommandPlugin, 0)
};

#endif
