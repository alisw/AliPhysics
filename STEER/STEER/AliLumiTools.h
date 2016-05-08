#ifndef ALILUMITOOLS_H
#define ALILUMITOOLS_H

#include <TString.h>
#include <TGraph.h>

class AliLumiTools : public TObject
{
 public:
  enum {kLumiCTP,kLumiDIP,kNLumiTypes};
 public:
  AliLumiTools() {}
  virtual ~AliLumiTools() {}
  //
  static TGraph* GetLumiFromCTP(Int_t run=-1, const char * ocdbPathDef="raw://", TString refClassName="", Double_t refSigma=-1);
  static TGraph* GetLumiFromDIP(Int_t run=-1, const char * ocdbPathDef="raw://");
  static TGraph* GetLumiGraph(int tp,Int_t run=-1, const char * ocdbPathDef="raw://");
  //
 protected:
  static Bool_t GetLumiCTPRefClass(int run, TString& refClass, double &refSigma);

 private:
  AliLumiTools(const AliLumiTools& src);
  AliLumiTools& operator=(const AliLumiTools& src) const;
  //
ClassDef(AliLumiTools,0)
};

#endif
