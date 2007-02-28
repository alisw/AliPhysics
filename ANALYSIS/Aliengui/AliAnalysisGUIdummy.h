#ifndef ALIANALYSISGUIDUMMY_H
#define ALIANALYSISGUIDUMMY_H

#include <cstdlib>
#include <TGFrame.h>
#include <AliLog.h>

class TGWindow;

//___________________________________________________________________________
class AliAnalysisGUI : public TGMainFrame {
  
 public:
  AliAnalysisGUI(const TGWindow *, UInt_t , UInt_t ) {
    AliError("No XML support in Root! Exit...");
    exit(1);
  }
  ~AliAnalysisGUI(){}

  ClassDef(AliAnalysisGUI, 0); // AliAnalysisGUI
};

#endif


