#ifndef ALILUMITOOLS_H
#define ALILUMITOOLS_H

#include <TString.h>
#include <TGraph.h>

class AliLumiTools : public TObject
{
 public:
  enum {kLumiCTP, // use CTP scalers
	kLumiDIP, // extract from Alice delivered lumi record (based on the same T0)
	kLumiDIPInst, // use Alice lumi from T0, communicated to LHC and extracted from the DIP file
	kLumiDIPDel,  // use Alice delivered lumi, communicated to LHC and extracted from the DIP file
	kNLumiTypes};
 public:
  AliLumiTools() {}
  virtual ~AliLumiTools() {}
  //
  static TGraph* GetLumiFromCTP(Int_t run=-1, const char * ocdbPathDef="raw://", TString refClassName="", Double_t refSigma=-1);
  static TGraph* GetLumiFromDIPDel(Int_t run=-1, const char * ocdbPathDef="raw://");
  static TGraph* GetLumiFromDIPInst(Int_t run=-1, const char * ocdbPathDef="raw://");
  static TGraph* GetLumiFromDIP(Int_t run=-1, const char * ocdbPathDef="raw://") { return GetLumiFromDIPDel(run, ocdbPathDef); }
  static TGraph* GetLumiGraph(int tp,Int_t run=-1, const char * ocdbPathDef="raw://");
  //
  static Double_t GetMuEstimate()    {return fgMuEst;}
  static Double_t GetXSectEstimate() {return fgXSecEst;}
  //
  static Float_t  GetScaleDnDeta2pp13TeV(int run=-1,const char * ocdbPathDef="raw://"); 
  //
  static Double_t GetScaleFactor() {return fgScaleFactor;}
  static void     SetScaleFactor(double v=1) {fgScaleFactor = v;}
  //
 protected:
  static Bool_t GetLumiCTPRefClass(int run, TString& refClass, double &refSigma, double &refEff, const char* ocdbPathDef);
  static TObject* GetCDBObjectForRun(int run, const char* path, const char* ocdbPathDef);

 protected:
  static Double_t fgMuEst;     // mu estimate
  static Double_t fgXSecEst;   // X-section to be used for conversion to IR
  static Double_t fgScaleFactor; // user settable scaling factor to tweak the lumi

 private:  // dummy
  AliLumiTools(const AliLumiTools& src);
  AliLumiTools& operator=(const AliLumiTools& src) const;
  //
ClassDef(AliLumiTools,0)
};

#endif
