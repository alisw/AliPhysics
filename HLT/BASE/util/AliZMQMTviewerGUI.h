#include "TGFrame.h"
#include <string>
#include "TCanvas.h"
#include "TRootCanvas.h"

class TThread;
class TH1F;
class TRootEmbeddedCanvas;
class TGTextButton;
class TGTextEntry;
class TCanvas;
struct AliZMQhistViewer;
class TVirtualPad;
class TPad;
class TGStatusBar;
class TGTextButton;

typedef struct {
   Bool_t    fRun;
   TThread  *fThread;
   TCanvas  *fCanvas;
   AliZMQhistViewer* viewer;
} ThreadArgs_t;

class AliZMQMTviewerGUIview {
  public:
  AliZMQMTviewerGUIview(const char* name, const char* title, Int_t a, Int_t b, Int_t c, Int_t d)
    : fCanvas(name,title,a,b,c,d), fDrawnObjects() {
      TRootCanvas *canvas = (TRootCanvas *)fCanvas.GetCanvasImp();
      canvas->Connect("CloseWindow()","AliZMQMTviewerGUIview",this,"CleanUp()");
    }
  TCanvas fCanvas;
  TList fDrawnObjects;
  void CleanUp();
  virtual ~AliZMQMTviewerGUIview() {}
  ClassDef(AliZMQMTviewerGUIview, 0)
};

class AliZMQMTviewerGUI : public TGMainFrame {
private:
   TRootEmbeddedCanvas* fECanvas;
   TGStatusBar*         fStatusBar;   
   TGTextButton*        fPullButton;
   TGTextEntry*         fSelectionEntry;
   TGTextEntry*         fUnSelectionEntry;
   TThread*             fThread;
   TCanvas*             fCanvas;
   ThreadArgs_t         fThreadArgs;
   AliZMQhistViewer*    fViewer;
   TTimer*              fTimer;
   std::string          fWindowTitle;
   void*                fZMQviewerConfig;
   int                  fInitStatus;
   TPRegexp*            fSelection;
   TPRegexp*            fUnSelection;

   AliZMQMTviewerGUI(const AliZMQMTviewerGUI&);
   AliZMQMTviewerGUI& operator=(const AliZMQMTviewerGUI&);
public:
   AliZMQMTviewerGUI(const TGWindow *p, UInt_t w, UInt_t h, int argc, char** argv);
   virtual ~AliZMQMTviewerGUI() { }

   void     CloseWindow();
   Bool_t   HandleTimer(TTimer *t);
   void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);   
   void PadSelected(TVirtualPad* pad, TObject *object, Int_t event);
   void PadPicked(TPad* pad, TObject *object, Int_t event);
   void UpdateCanvas();
   void ReconfigureViewer(std::string string);
   int GetInitStatus() const {return fInitStatus;}
   void DoPullButton();
   void DoSelectionEntry();

   // thread method
   static void* RunViewer(void *ptr = 0);

   static const char* fUSAGE;

   ClassDef(AliZMQMTviewerGUI, 0)
};




