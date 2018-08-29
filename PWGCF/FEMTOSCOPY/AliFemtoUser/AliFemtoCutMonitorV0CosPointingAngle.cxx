/// \file AliFemtoCutMonitorV0CosPointingAngle.cxx

#include "AliFemtoCutMonitorV0CosPointingAngle.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include "AliFemtoModelHiddenInfo.h"

AliFemtoCutMonitorV0CosPointingAngle::AliFemtoCutMonitorV0CosPointingAngle(int aPrimaryPID, bool aBuildCosPointingAnglewParentInfo, int aNbins, double aCosMin, double aCosMax):
  AliFemtoCutMonitorV0(),
  fNbinsCPA(aNbins),
  fMinCPA(aCosMin), 
  fMaxCPA(aCosMax), 
  fPrimaryPID(aPrimaryPID),
  fBuildCosPointingAnglewParentInfo(true),
  fParentPIDInfoVec(0),
  fCosPointingAnglewParentInfo(nullptr)
{
  //if fBuildCosPointingAnglewParentInfo=false, the following line will simply rebin fCosPointingAngle
  //which is good, as the default binning isn't very useful
  if(fBuildCosPointingAnglewParentInfo) CreateDefaultParentPIDInfoVec();
  SetBuildCosPointingAnglewParentInfo(fBuildCosPointingAnglewParentInfo, "", fNbinsCPA, fMinCPA, fMaxCPA);
}

AliFemtoCutMonitorV0CosPointingAngle::AliFemtoCutMonitorV0CosPointingAngle(const char *aName, int aPrimaryPID, bool aBuildCosPointingAnglewParentInfo, int aNbins, double aCosMin, double aCosMax):
  AliFemtoCutMonitorV0(aName),
  fNbinsCPA(aNbins),
  fMinCPA(aCosMin), 
  fMaxCPA(aCosMax), 
  fPrimaryPID(aPrimaryPID),
  fBuildCosPointingAnglewParentInfo(true),
  fParentPIDInfoVec(0),
  fCosPointingAnglewParentInfo(nullptr)
{
  //if fBuildCosPointingAnglewParentInfo=false, the following line will simple rebin fCosPointingAngle
  //which is good, as the default binning isn't very useful
  if(fBuildCosPointingAnglewParentInfo) CreateDefaultParentPIDInfoVec();
  SetBuildCosPointingAnglewParentInfo(fBuildCosPointingAnglewParentInfo, TString(aName), fNbinsCPA, fMinCPA, fMaxCPA);
}

AliFemtoCutMonitorV0CosPointingAngle::AliFemtoCutMonitorV0CosPointingAngle(const AliFemtoCutMonitorV0CosPointingAngle &aCut):
  AliFemtoCutMonitorV0(aCut),
  fNbinsCPA(aCut.fNbinsCPA),
  fMinCPA(aCut.fMinCPA), 
  fMaxCPA(aCut.fMaxCPA), 
  fPrimaryPID(aCut.fPrimaryPID),
  fBuildCosPointingAnglewParentInfo(aCut.fBuildCosPointingAnglewParentInfo),
  fParentPIDInfoVec(aCut.fParentPIDInfoVec),
  fCosPointingAnglewParentInfo(nullptr)
{
  /// copy constructor

  if (fCosPointingAnglewParentInfo) delete fCosPointingAnglewParentInfo;
  fCosPointingAnglewParentInfo = new TH2F(*aCut.fCosPointingAnglewParentInfo);
  fCosPointingAnglewParentInfo->Sumw2();

}

AliFemtoCutMonitorV0CosPointingAngle::~AliFemtoCutMonitorV0CosPointingAngle()
{
  /// Destructor
  if(fBuildCosPointingAnglewParentInfo) delete fCosPointingAnglewParentInfo;
}

AliFemtoCutMonitorV0CosPointingAngle& AliFemtoCutMonitorV0CosPointingAngle::operator=(const AliFemtoCutMonitorV0CosPointingAngle& aCut)
{
  /// assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitorV0::operator=(aCut);

  fNbinsCPA = aCut.fNbinsCPA;
  fMinCPA = aCut.fMinCPA;
  fMaxCPA = aCut.fMaxCPA;
  fPrimaryPID = aCut.fPrimaryPID;
  fBuildCosPointingAnglewParentInfo = aCut.fBuildCosPointingAnglewParentInfo;
  fParentPIDInfoVec = aCut.fParentPIDInfoVec;

  if (fCosPointingAnglewParentInfo) delete fCosPointingAnglewParentInfo;
  fCosPointingAnglewParentInfo = new TH2F(*aCut.fCosPointingAnglewParentInfo);
  fCosPointingAnglewParentInfo->Sumw2();

  return *this;
}

AliFemtoString AliFemtoCutMonitorV0CosPointingAngle::Report(){
  /// Prepare report from the execution

  string stemp = "*** AliFemtoCutMonitorV0CosPointingAngle report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

double AliFemtoCutMonitorV0CosPointingAngle::GetMotherBin(int aPID, int aMotherPID)
{
  //Note: subtract 0.5 from return value to place it within bin, instead of at boundary
  if(aPID==fPrimaryPID)
  {
    if(aMotherPID==0) return 1-0.5;
    for(unsigned int i=1; i<fParentPIDInfoVec.size(); i++)
    {
      if(aMotherPID==fParentPIDInfoVec[i].parentPID) return (i+1)-0.5;
    }
    return fParentPIDInfoVec.size()+1-0.5;
  }
  else return fParentPIDInfoVec.size()+2-0.5;
}


void AliFemtoCutMonitorV0CosPointingAngle::Fill(const AliFemtoV0* aV0)
{
  /// Fill momentum resolution histograms for the particle
  AliFemtoCutMonitorV0::Fill(aV0);

  if(fBuildCosPointingAnglewParentInfo)
  {
    AliFemtoModelHiddenInfo *tInfo = (AliFemtoModelHiddenInfo*)aV0->GetHiddenInfo();
    if(tInfo!=NULL) {
      Int_t partID = TMath::Abs(tInfo->GetPDGPid());
      Int_t motherID = TMath::Abs(tInfo->GetMotherPdgCode());

      double tMotherBin = GetMotherBin(partID, motherID);
      fCosPointingAnglewParentInfo->Fill(aV0->CosPointingAngle(), tMotherBin);
    }
  }
}

void AliFemtoCutMonitorV0CosPointingAngle::Write()
{
  /// Write out the relevant histograms
  AliFemtoCutMonitorV0::Write();
  if(fBuildCosPointingAnglewParentInfo) fCosPointingAnglewParentInfo->Write();
}

TList *AliFemtoCutMonitorV0CosPointingAngle::GetOutputList()
{
  /// Get the list of histograms to write
  TList *tOutputList = AliFemtoCutMonitorV0::GetOutputList();
  if(fBuildCosPointingAnglewParentInfo) tOutputList->Add(fCosPointingAnglewParentInfo);

  return tOutputList;
}


void AliFemtoCutMonitorV0CosPointingAngle::CreateDefaultParentPIDInfoVec()
{
  fParentPIDInfoVec.clear();
  ParentPIDInfo tPIDInfo;

  tPIDInfo.parentName = TString("Primary");
  tPIDInfo.parentPID = fPrimaryPID;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Lambda");
  tPIDInfo.parentPID = 3122;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("K^{0}_{S}");
  tPIDInfo.parentPID = 310;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Sigma^{0}");
  tPIDInfo.parentPID = 3212;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Xi^{0}");
  tPIDInfo.parentPID = 3322;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Xi^{-}");
  tPIDInfo.parentPID = 3312;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Sigma^{*+}");
  tPIDInfo.parentPID = 3224;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Sigma^{*-}");
  tPIDInfo.parentPID = 3114;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("#Sigma^{*0}");
  tPIDInfo.parentPID = 3214;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("K^{*0}");
  tPIDInfo.parentPID = 313;
  fParentPIDInfoVec.push_back(tPIDInfo);

  tPIDInfo.parentName = TString("K^{*+}");
  tPIDInfo.parentPID = 323;
  fParentPIDInfoVec.push_back(tPIDInfo);
}

void AliFemtoCutMonitorV0CosPointingAngle::AddParentToPIDInfoVec(TString aName, int aPID)
{
  ParentPIDInfo tPIDInfo;
  tPIDInfo.parentName = aName;
  tPIDInfo.parentPID = aPID;

  fParentPIDInfoVec.push_back(tPIDInfo);
}


void AliFemtoCutMonitorV0CosPointingAngle::CreateNewParentPIDInfoVec(vector<TString> &aNames, vector<int> &aPIDs)
{
  fParentPIDInfoVec.clear();
  ParentPIDInfo tPIDInfo;

  tPIDInfo.parentName = TString("Primary");
  tPIDInfo.parentPID = fPrimaryPID;
  fParentPIDInfoVec.push_back(tPIDInfo);

  for(unsigned int i=0; i<aNames.size(); i++) 
  {
    tPIDInfo.parentName = aNames[i];
    tPIDInfo.parentPID = aPIDs[i];
    fParentPIDInfoVec.push_back(tPIDInfo);
  }
}


void AliFemtoCutMonitorV0CosPointingAngle::SetBuildCosPointingAnglewParentInfo(bool aBuild, TString aName, int aNbins, double aCosMin, double aCosMax)
{
  fBuildCosPointingAnglewParentInfo = aBuild;
  if(fBuildCosPointingAnglewParentInfo)
  {
    TString tName = TString::Format("CosPointingAnglewParentInfo%s", aName.Data());
    int tNbinsPID = fParentPIDInfoVec.size();
    fCosPointingAnglewParentInfo = new TH2F(tName, tName, 
                                            aNbins, aCosMin, aCosMax,
                                            tNbinsPID+2, 0, tNbinsPID+2);
    for(int i=0; i<tNbinsPID; i++)
    {
      fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(i+1, fParentPIDInfoVec[i].parentName);
    }
    fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+1, "Other");
    fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+2, "Fake");

    fCosPointingAnglewParentInfo->Sumw2();
  }

  fCosPointingAngle = new TH1F(fCosPointingAngle->GetName(), fCosPointingAngle->GetTitle(),
                               aNbins, aCosMin, aCosMax);
  fCosPointingAngle->Sumw2();
}


