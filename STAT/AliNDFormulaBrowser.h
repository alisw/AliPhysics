#include "TObject.h"
#include <map>

class TGLayoutHints;
class TGTransientFrame;
class TGVerticalFrame;
class TRootEmbeddedCanvas;
class TGGroupFrame;
class TGTextEntry;
class TGLabel;
class TFormula;
class TGWindow;

class AliNDFormulaBrowser : public TObject{

  //  RQ_OBJECT("AliNDFormulaBrowser")

public:
  TGTransientFrame  *fMain;
  TGVerticalFrame   *fVframeFormula;
  TGLayoutHints     *fBly, *fBfly1;
  //
  //                                     // fromula description
  TFormula          *fFormula;           // pointer to the formula
  TGGroupFrame      *fFormulaFrame;      // formula frame
  TGTextEntry       *fDrawFormula;       // draw formula  
  TGLabel           *fDrawFormulaLabel;  // draw label
  TGTextEntry       *fDrawOption;        // draw option  
  TGTextEntry       *fFormulaEval;       // formula value
  TGTextEntry       *fFormulaNPoints;    // formula number of points
  TGTextEntry       *fFormulaNLines;     // formula number of lines
  //
  //                                     // parameters description
  TVectorD          *fFormulaParams;     // parameter vector
  TObjArray         *fParamFrame;        // frame
  TObjArray         *fSliders;           // array of sliders
  TObjArray         *fCurrentEntry;      // array of current entries 
  TObjArray         *fMinEntry;          // array of min values
  TObjArray         *fMaxEntry;          // array of max values 
  TObjArray         *fParamLabels;       // array of param lables (description)  
  TRootEmbeddedCanvas *fCanvas;          // canvas
  static TLatex *   fgkLatex;            // Latex style
  static TLegend *  fgkLegend;           // Legend 
  static void SetDefaultStyle();         // SetDefaultStyle
  static std::map<std::string,TFormula*> fgkFormulaMap;   // map of registered formulas 
public:
  AliNDFormulaBrowser(const TGWindow *p, const TGWindow *main, TFormula *formula,  UInt_t width, UInt_t height);
  virtual ~AliNDFormulaBrowser();
  // slots
  void CloseWindow();
  void DoText(const char *text);
  void DoSlider(Int_t pos = 0);
  void UpdateFormula(); 
  void UpdateCanvas(); 
  ClassDef(AliNDFormulaBrowser,0);
};
