#ifndef ALIMCTREETOOLS_H
#define ALIMCTREETOOLS_H


class AliMCTreeTools : public TNamed
{
  public:
  Int_t fVerobse;
  static TChain *fKineChain;
  static TChain *fTRChain;
  static void SetWDir(TString dir);
  static void ClearCache();
  static void InitStack(Int_t iEvent);
  static Double_t GetValueAt(Int_t iEvent, Int_t itrack, Double_t x, Double_t y, Double_t z, Int_t returnValue, Int_t interpolationType, Int_t verbose, Double_t rotation=0);
  static Double_t FindNearestReference(Int_t iEvent, Int_t itrack, Double_t x, Double_t y, Double_t z, Int_t returnValue, Int_t verbose=0);
  static TChain *  MakeKineChain(const char *inputList="kinematics.list");
  ClassDef(AliMCTreeTools, 1) 
};

#endif
