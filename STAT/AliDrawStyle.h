
#include "TObject.h"
#include <map>
#include <vector>
#include <string>
#include "TString.h"
class TPRegexp; 
class TStyle;

class AliDrawStyle : public TObject{
public:
  static void ApplyStyle(const char* styleName);
  static void SetDefaults();
  static TString GetLatexAlice(const char * symbol);
  static void AddLatexSymbol(const char * symbolName, const char * symbolTitle);
  static Int_t GetMarkerStyle(const char *style, Int_t index);
  static Int_t GetMarkerColor(const char *style, Int_t index);
  static Int_t GetFillColor(const char *style, Int_t index); 
  static void PrintLatexSymbols(Option_t *option,TPRegexp& regExp);
  static void PrintStyles(Option_t *option, TPRegexp& regExp);
protected:
  static std::map<TString, TString> fLatexAlice;
  static std::map<TString, TStyle*>  fStyleAlice;
  static std::map<TString, std::vector<int>> fMarkerStyles;
  static std::map<TString, std::vector<int>> fMarkerColors;
  static std::map<TString, std::vector<int>> fFillColors;
  //
  static void  RegisterDefaultLatexSymbols();
  static void  RegisterDefaultStyle();
  static void  RegisterDefaultMarkers();

  ClassDef(AliDrawStyle,1);
};
