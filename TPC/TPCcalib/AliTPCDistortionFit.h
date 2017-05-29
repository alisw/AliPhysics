#ifndef ALITPCDISTORTIONFIT_H
#define ALITPCDISTORTIONFIT_H

#include <map>
class AliTPCChebCorr;
class TPRegexp;

class AliTPCDistortionFit : public TNamed{
public:
  // static interface to distortion map - AliTPCChebCorr - used for TTree and TFormula evaluation  
  static Int_t RegisterFitters();
  static Double_t LineFieldLocal(const Double_t *x, const Double_t *param);  
  static Double_t NLinesFieldLocal(const Double_t *x, const Double_t *param);  
  static Int_t LoadDistortionMaps(Int_t run, const char *storage=NULL);
  static TString   LoadDistortionTree(const char *chinput);
  static void  SetMetadata(TTree * tree, TString friendName,Int_t run);
  static void  PrintMap(TPRegexp *filter=NULL);
  static void MakeFitExample1(Int_t run, const char * chinput);
  //
  static void   RegisterMap(TString mapName,  AliTPCChebCorr *map);
  static Int_t  GetNameHash(string name){  return fgkMapNameHash[name];}
  static string GetHashName(Int_t index){  return fgkMapHashName[index];}
  static const AliTPCChebCorr * GetCheb(Int_t index){  return fgkMapCheb[index];}
  static const AliTPCChebCorr * GetCheb(string name){  return fgkMapCheb[fgkMapNameHash[name]];}
  // static eval
  static Double_t Eval(Int_t hashIndex, int sector, int row, float y2x, float z2x, int dimOut);
  static Double_t EvalEfield(Int_t hashIndex, Int_t hashref, int sector, int row, float y2x, float z2x, int dimOut, Int_t dir, Double_t wt, Double_t dz, Double_t Ez);
  static Double_t EvalRho(Int_t hashIndex, Int_t hashref, int sector, int row, float y2x, float z, Int_t dir, Double_t wt, Double_t dz,  Double_t dr, Double_t drphi,  Double_t Ez);
  // static eval with interpolation
  static Double_t EvalSector(Int_t hashIndex, Double_t fsector, int row, float z2x, int dimOut, int interpol, Double_t mndiv);  
  static Double_t EvalSector(Int_t hashIndex, Int_t hashRef,  Double_t fsector, int row, float z2x, int dimOut);  
  static Double_t EvalSectorBckg(Int_t hashIndex, Int_t hashRef, Double_t fsector, int row, float z2x, int dimOut, Double_t dRow, Double_t dSec);  
  static Double_t EvalEfieldSector(Int_t hashIndex, Int_t hashref, Float_t fsector, int row, float z2x, int dimOut, Int_t polarity, Double_t wt, Double_t dz, Double_t Ez, Int_t interpolationType, Double_t mBinSize);
  static Double_t EvalRhoSector(Int_t hashIndex, Int_t hashref, Float_t fsector, int row, float z2x, Int_t polarity, Double_t wt, Double_t dz,  Double_t dr, Double_t drphi,  Double_t Ez, Int_t interpolationType, Double_t mBinSize);
protected:
  static  map<string,int> fgkMapNameHash;              // map  name->index
  static  map<int,string> fgkMapHashName;              // map  index->name
  static  map<int,const AliTPCChebCorr *>fgkMapCheb;   // map  index->AliTPCChebCorr *>fgkMapCheb
public:
  static  TTree * fgkDistortionTree;                   // registered distortion tree
private:  
  ClassDef(AliTPCDistortionFit, 1)   // 
};


#endif

