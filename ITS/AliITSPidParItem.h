#ifndef AliITSPIDPARITEM_H
#define AliITSPIDPARITEM_H
//////////////////////////////////////////////////////////
//Class for PID in the ITS                              //
//Origin: Elena Bruna bruna@to.infn.it                  //
//////////////////////////////////////////////////////////

class AliITSPidParItem : public TObject{
  
 public:
  AliITSPidParItem();  
  AliITSPidParItem(Float_t center,Float_t width,Double_t *buff);  
  virtual ~AliITSPidParItem(){;}
  Float_t GetMomentumCenter() const {return fPCenter;}
  Float_t GetWidthMom() const {return fPWidth;}
  void GetParameters(Double_t *buff) const;
  void PrintParameters() const;
  void GetProtonPar(Double_t *buffp) const;
  void GetKaonPar(Double_t *buffk) const;
  void GetPionPar(Double_t *buffpi) const;
  void GetPar0(Double_t *buff0) const;
  void GetPar1(Double_t *buff1) const;
  void GetPar2(Double_t *buff2) const;
  void GetPar3(Double_t *buff3) const;
  void GetChisquare(Double_t *buffchi) const;
  void GetNDF(Double_t *buffndf) const;
  void GetProParFun(Double_t *pfun) const;
  void GetKaoParFun(Double_t *kfun) const;
  void GetPiParFun(Double_t *pifun) const;
  void GetRangeLim(Double_t *range) const;
  void GetProtonParErr(Double_t *bufferp)const;
  void GetKaonParErr(Double_t *bufferk)const;
  void GetPionParErr(Double_t *bufferpi) const;
  static TF1* CookFunIts(TString namefun,Double_t *par,Double_t rangei,Double_t rangef,TString comment);
  static Double_t Langaufun(Double_t *x, Double_t *par);
  static Double_t Langaufun2(Double_t *x, Double_t *par);
  
 private:
  Float_t fPCenter;     //value for center
  Float_t fPWidth;      //value for width
  Double_t fBuff[39];   //array for PID parameters
  
  ClassDef(AliITSPidParItem,1);
};
#endif
