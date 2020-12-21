/// \file AliFemtoCutMonitorV0CosPointingAngle.cxx

#include "AliFemtoCutMonitorV0CosPointingAngle.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>


AliFemtoCutMonitorV0CosPointingAngle::AliFemtoCutMonitorV0CosPointingAngle(int aPrimaryPID, bool aBuildCosPointingAnglewParentInfo, int aNbins, double aCosMin, double aCosMax):
  AliFemtoCutMonitorV0(),
  fNbinsCPA(aNbins),
  fMinCPA(aCosMin), 
  fMaxCPA(aCosMax), 
  fPrimaryPID(aPrimaryPID),
  fBuildCosPointingAnglewParentInfo(true),
  fParentPIDInfoVec(0),
  fCosPointingAnglewParentInfo(nullptr),
  fDcaV0ToPrimVertexwParentInfo(nullptr)
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
  fCosPointingAnglewParentInfo(nullptr),
  fDcaV0ToPrimVertexwParentInfo(nullptr)
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
  fCosPointingAnglewParentInfo(nullptr),
  fDcaV0ToPrimVertexwParentInfo(nullptr)
{
  /// copy constructor

  if (fCosPointingAnglewParentInfo) delete fCosPointingAnglewParentInfo;
  fCosPointingAnglewParentInfo = new TH2F(*aCut.fCosPointingAnglewParentInfo);
  fCosPointingAnglewParentInfo->Sumw2();

  if (fDcaV0ToPrimVertexwParentInfo) delete fDcaV0ToPrimVertexwParentInfo;
  fDcaV0ToPrimVertexwParentInfo = new TH2F(*aCut.fDcaV0ToPrimVertexwParentInfo);
  fDcaV0ToPrimVertexwParentInfo->Sumw2();

}

AliFemtoCutMonitorV0CosPointingAngle::~AliFemtoCutMonitorV0CosPointingAngle()
{
  /// Destructor
  if(fBuildCosPointingAnglewParentInfo) 
  {
    delete fCosPointingAnglewParentInfo;
    delete fDcaV0ToPrimVertexwParentInfo;
  }
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

  if (fDcaV0ToPrimVertexwParentInfo) delete fDcaV0ToPrimVertexwParentInfo;
  fDcaV0ToPrimVertexwParentInfo = new TH2F(*aCut.fDcaV0ToPrimVertexwParentInfo);
  fDcaV0ToPrimVertexwParentInfo->Sumw2();

  return *this;
}

AliFemtoString AliFemtoCutMonitorV0CosPointingAngle::Report(){
  /// Prepare report from the execution

  string stemp = "*** AliFemtoCutMonitorV0CosPointingAngle report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

double AliFemtoCutMonitorV0CosPointingAngle::GetMotherBin(AliFemtoModelHiddenInfo *aInfo)
{
  int tPID = TMath::Abs(aInfo->GetPDGPid());
  int tMotherPID = TMath::Abs(aInfo->GetMotherPdgCode());

/*
  //Note: subtract 0.5 from return value to place it within bin, instead of at boundary
  if(aInfo->GetOrigin()==2) return fParentPIDInfoVec.size()+3-0.5;
  if(tPID==fPrimaryPID)
  {
    if(tMotherPID==0) return 1-0.5;
    for(unsigned int i=1; i<fParentPIDInfoVec.size(); i++)
    {
      if(tMotherPID==fParentPIDInfoVec[i].parentPID) return (i+1)-0.5;
    }
    return fParentPIDInfoVec.size()+1-0.5;
  }
  else return fParentPIDInfoVec.size()+2-0.5;
*/

  //Note: subtract 0.5 from return value to place it within bin, instead of at boundary
  if(aInfo->GetOrigin()==2) return fParentPIDInfoVec.size()+3-0.5;  //Material
  if(tPID==fPrimaryPID)  //real, not fake
  {
    if(tMotherPID==0)
    {
      if(aInfo->GetOrigin()==0) return 1-0.5;  //this is most certainly primary
      else return -1;                          //Disagreement between tMotherPID and aInfo->GetOrigin()
    }
    for(unsigned int i=1; i<fParentPIDInfoVec.size(); i++)
    {
      if(tMotherPID==fParentPIDInfoVec[i].parentPID) 
      {
        if(tMotherPID==3212) return (i+1)-0.5;  //Sigma0, not sure what aInfo->GetOrigin() should be, as this is EM decay, not weak
        else
        {
          if(aInfo->GetOrigin()==1) return (i+1)-0.5;  //This is secondary
          else return -1;                              //Disagreement between tMotherPID and aInfo->GetOrigin()
        }
      }
    }
    return fParentPIDInfoVec.size()+1-0.5;  //Other
  }
  else return fParentPIDInfoVec.size()+2-0.5;  //Fake
}


void AliFemtoCutMonitorV0CosPointingAngle::Fill(const AliFemtoV0* aV0)
{
  /// Fill momentum resolution histograms for the particle
  AliFemtoCutMonitorV0::Fill(aV0);

  if(fBuildCosPointingAnglewParentInfo)
  {
    AliFemtoModelHiddenInfo *tInfo = (AliFemtoModelHiddenInfo*)aV0->GetHiddenInfo();
    if(tInfo!=NULL) {
      double tMotherBin = GetMotherBin(tInfo);
      fCosPointingAnglewParentInfo->Fill(aV0->CosPointingAngle(), tMotherBin);
      fDcaV0ToPrimVertexwParentInfo->Fill(aV0->DcaV0ToPrimVertex(), tMotherBin);
    }
  }
}

void AliFemtoCutMonitorV0CosPointingAngle::Write()
{
  /// Write out the relevant histograms
  AliFemtoCutMonitorV0::Write();
  if(fBuildCosPointingAnglewParentInfo) 
  {
    fCosPointingAnglewParentInfo->Write();
    fDcaV0ToPrimVertexwParentInfo->Write();
  }
}

TList *AliFemtoCutMonitorV0CosPointingAngle::GetOutputList()
{
  /// Get the list of histograms to write
  TList *tOutputList = AliFemtoCutMonitorV0::GetOutputList();
  if(fBuildCosPointingAnglewParentInfo) 
  {
    tOutputList->Add(fCosPointingAnglewParentInfo);
    tOutputList->Add(fDcaV0ToPrimVertexwParentInfo);
  }

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
    TString tNameCPA = TString::Format("CosPointingAnglewParentInfo%s", aName.Data());
    int tNbinsPID = fParentPIDInfoVec.size();
    fCosPointingAnglewParentInfo = new TH2F(tNameCPA, tNameCPA, 
                                            aNbins, aCosMin, aCosMax,
                                            tNbinsPID+3, 0, tNbinsPID+3);

    TString tNameDCA = TString::Format("DcaV0ToPrimVertexwParentInfo%s", aName.Data());
    fDcaV0ToPrimVertexwParentInfo = new TH2F(tNameDCA, tNameDCA, 
                                             fDcaV0ToPrimVertex->GetNbinsX(), 
                                             fDcaV0ToPrimVertex->GetXaxis()->GetBinLowEdge(1), 
                                             fDcaV0ToPrimVertex->GetXaxis()->GetBinUpEdge(fDcaV0ToPrimVertex->GetNbinsX()),
                                             tNbinsPID+3, 0, tNbinsPID+3);

    for(int i=0; i<tNbinsPID; i++)
    {
      fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(i+1, fParentPIDInfoVec[i].parentName);
      fDcaV0ToPrimVertexwParentInfo->GetYaxis()->SetBinLabel(i+1, fParentPIDInfoVec[i].parentName);
    }
    fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+1, "Other");
    fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+2, "Fake");
    fCosPointingAnglewParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+3, "Material");
    fCosPointingAnglewParentInfo->Sumw2();

    fDcaV0ToPrimVertexwParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+1, "Other");
    fDcaV0ToPrimVertexwParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+2, "Fake");
    fDcaV0ToPrimVertexwParentInfo->GetYaxis()->SetBinLabel(tNbinsPID+3, "Material");
    fDcaV0ToPrimVertexwParentInfo->Sumw2();
  }

  fCosPointingAngle = new TH1F(fCosPointingAngle->GetName(), fCosPointingAngle->GetTitle(),
                               aNbins, aCosMin, aCosMax);
  fCosPointingAngle->Sumw2();
}


