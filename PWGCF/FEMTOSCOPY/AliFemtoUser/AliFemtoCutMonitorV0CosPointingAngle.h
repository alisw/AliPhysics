/// \file AliFemtoCutMonitorV0CosPointingAngle
// Purpose is to use 2d fCosPointingAnglewParentInfo with MC data, project out different
// contributions, and do template fit to experimental data to determine ratio of primary Lambda to
// Lambdas from secondary sources

#ifndef AliFemtoCutMonitorV0CosPointingAngle_H
#define AliFemtoCutMonitorV0CosPointingAngle_H

class TH2F;
#include "AliFemtoCutMonitorV0.h"
#include "AliFemtoModelHiddenInfo.h"

class AliFemtoCutMonitorV0CosPointingAngle : public AliFemtoCutMonitorV0 {
public:
  AliFemtoCutMonitorV0CosPointingAngle(int aPrimaryPID, bool aBuildCosPointingAnglewParentInfo=false, int aNbins=100, double aCosMin=0.99, double aCosMax=1.0);
  AliFemtoCutMonitorV0CosPointingAngle(const char *aName, int aPrimaryPID, bool aBuildCosPointingAnglewParentInfo=false, int aNbins=100, double aCosMin=0.99, double aCosMax=1.0);
  AliFemtoCutMonitorV0CosPointingAngle(const AliFemtoCutMonitorV0CosPointingAngle &aCut);
  virtual ~AliFemtoCutMonitorV0CosPointingAngle();

  AliFemtoCutMonitorV0CosPointingAngle& operator=(const AliFemtoCutMonitorV0CosPointingAngle& aCut);

  virtual AliFemtoString Report();
  double GetMotherBin(AliFemtoModelHiddenInfo *aInfo);
  virtual void Fill(const AliFemtoV0* aV0);
  void Write();

  virtual TList *GetOutputList();


  struct ParentPIDInfo
  {
    TString parentName="";
    int parentPID=-1;
  };


  void CreateDefaultParentPIDInfoVec();
  void AddParentToPIDInfoVec(TString aName, int aPID);
  void CreateNewParentPIDInfoVec(vector<TString> &aNames, vector<int> &aPIDs);
  void SetBuildCosPointingAnglewParentInfo(bool aBuild, TString aName="", int aNbins=100, double aCosMin=0.99, double aCosMax=1.0);

protected:
  int fNbinsCPA;
  double fMinCPA, fMaxCPA;
  int fPrimaryPID;
  bool fBuildCosPointingAnglewParentInfo;
  vector<ParentPIDInfo> fParentPIDInfoVec;
  TH2F *fCosPointingAnglewParentInfo;
  TH2F *fDcaV0ToPrimVertexwParentInfo;

#ifdef __ROOT__
  ClassDef(AliFemtoCutMonitorV0CosPointingAngle, 1);
#endif
};

#endif
