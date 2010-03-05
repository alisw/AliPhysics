#ifndef ALITPCCALIBQACHECKER_H
#define ALITPCCALIBQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//                  QA checking class                                                  //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TString.h>
#include <TH1.h>

class TTree;
class TIterator;
class TGraph;
class TObjArray;
class TVirtualPad;

class AliTPCCalibQAChecker : public TNamed {
public:
  enum { kNQualityFlags=5 };
  enum QualityFlag_t { kNULLFLAG=-1, kINFO, kWARNING, kERROR, kFATAL, kNFLAG };
  enum AlarmType_t   { kMean=0, kBinAny, kBinAll, kNentries};

  AliTPCCalibQAChecker();
  AliTPCCalibQAChecker(const char* name, const char *title);
  
  virtual ~AliTPCCalibQAChecker();
  
  void Process();

  void SetTreeChecker(TTree* &tree)       {fTreePtr=&tree;}
  void SetHistChecker(TH1* &hist)         {fHistPtr=&hist;}
  void SetGraphChecker(TGraph* &graph)    {fGraphPtr=&graph;}
  void SetNumberChecker(Double_t & number) {fNumberPtr=&number;}

  const AliTPCCalibQAChecker* GetSubChecker(const char* name, Bool_t recursive=kTRUE) const;
  AliTPCCalibQAChecker* NextSubChecker();
  Int_t GetNumberOfSubCheckers(Bool_t recursive=kTRUE) const;
  Bool_t HasSubCheckers() const {return GetNumberOfSubCheckers(kFALSE)>0;}
  
  void AddSubChecker(AliTPCCalibQAChecker *alarm);

  //getters
  void GetAlarmThreshold(Double_t &min, Double_t &max, QualityFlag_t quality=kERROR) const {min=fThresMin[quality]; max=fThresMax[quality];}
  //
  const char* GetDrawString() { return fStrDraw.Data(); }
  const char* GetCutsString() { return fStrCuts.Data(); }
  const char* GetDrawOptString() {return fStrDrawOpt.Data(); }
  
  //tree related
  void SetDrawRepresentation(const char *draw, const char* drawOpt="") {fStrDrawRep=draw; fStrDrawRepOpt=drawOpt;}
  void SetDrawAlarm(const char *draw, const char* drawOpt="")          {fStrDraw=draw; fStrDrawOpt=drawOpt;}
  void SetCutString(const char *cutString )                            {fStrCuts=cutString;}

  //general thresholds for the different qualities
  void SetAlarmThreshold(const Double_t min, const Double_t max, const QualityFlag_t quality=kERROR);
  void ResetAlarmThreshold(const QualityFlag_t quality);
  void ResetAlarmThresholds();

  //descriptions
  void SetQualityDescription(const char* text, const QualityFlag_t quality=kERROR);
  
  //alarm type
  void SetAlarmType(AlarmType_t type) {fAlarmType=type;}
  

  QualityFlag_t GetQuality()      const {return fQualityLevel;}
  Color_t       GetQualityColor() const {return AliTPCCalibQAChecker::QualityColor(fQualityLevel);}
  const char*   GetQualityName()  const {return AliTPCCalibQAChecker::QualityName(fQualityLevel);}
  const char*   GetQualityDescription() const { return QualityDescription(fQualityLevel);}
  
  static const char* QualityName(const AliTPCCalibQAChecker::QualityFlag_t quality);
  static Color_t QualityColor(const AliTPCCalibQAChecker::QualityFlag_t quality);
  const char* QualityDescription(const QualityFlag_t quality) const;
  
  virtual void Draw(Option_t *option="");
  virtual void Print(Option_t *option="") const;  

 private:
  //alarm decision variables
  TTree   **fTreePtr;                    //! Pointer to the Tree pointer
  TH1     **fHistPtr;                    //! Pointer to the hist pointer
  TGraph  **fGraphPtr;                   //! Pointer to the graph pointer
  Double_t *fNumberPtr;                  //! Pointer to number
  TH1      *fHist;                       //! Hist pointer for tree processing

  TIterator *fIterSubCheckers;           //! iterator over sub checkers

  TObjArray *fArrSubCheckers;      //array with checkers to process
  TObjArray *fArrAlarmDescriptions; //array with alarm descriptions
  
  TString fStrDrawRep;             //draw string for representation histogram to visualise
  TString fStrDrawRepOpt;          //draw option for representation histogram
  
  TString fStrDraw;                //draw string for alarm histogram
  TString fStrDrawOpt;             //draw option for alarm histogram

  TString fStrCuts;                //cut string
  
  AlarmType_t fAlarmType;          //type of the alarm
  QualityFlag_t fQualityLevel;     //quality level
  
  TObject* fHistRep;                   //visualised histogram

  Double_t fThresMin[kNQualityFlags];//minimum thresholds
  Double_t fThresMax[kNQualityFlags];//maximum thresholds
  
  void CreateRepresentationHist();
  void ResetRepresentationHist() {if (fHistRep) {delete fHistRep; fHistRep=0x0;}}
  //general processing
  void ProcessTree();
  void ProcessHist();
  void ProcessGraph();
  void ProcessNumber();
  void ProcessSub();
  //special processing
  void ProcessEntries();
  void ProcessMean();
  void ProcessBin();
  //
  void CreateAlarmHist();
  void ResetAlarmHist();
  //
  Int_t DrawInPad(TVirtualPad *pad, Int_t sub=1);
  void DrawSubNodes(Option_t *option);
  void DrawRepresentationHist(const Option_t *option);
  void AddQualityLines(TH1 *hist);
  //
  AliTPCCalibQAChecker(const AliTPCCalibQAChecker &cfg);
  AliTPCCalibQAChecker& operator = (const AliTPCCalibQAChecker &cfg);
  
  QualityFlag_t GetQuality(Double_t value) const;
  QualityFlag_t GetQuality(Int_t n, const Double_t *arr) const;
  
  ClassDef(AliTPCCalibQAChecker,1);
};

//
//inline functions
//

//_________________________________________________________________________
inline AliTPCCalibQAChecker::QualityFlag_t AliTPCCalibQAChecker::GetQuality(Double_t value) const
{
  //
  // check quality of a value
  //

  QualityFlag_t quality=kINFO;
  //loop over Quality levels
  for (Int_t i=(Int_t)kINFO; i<kNQualityFlags; ++i){
    if (fThresMin[i]>=fThresMax[i]) continue;
    if (value<fThresMin[i]||value>fThresMax[i]) quality=(QualityFlag_t)i;
  }
  return quality;
}
//_________________________________________________________________________
inline AliTPCCalibQAChecker::QualityFlag_t AliTPCCalibQAChecker::GetQuality(Int_t n, const Double_t *arr) const
{
  //
  // check quality of an array
  //
  
  QualityFlag_t quality=kINFO;
  //loop over Quality levels
  for (Int_t i=(Int_t)kINFO; i<kNQualityFlags; ++i){
    if (fThresMin[i]>=fThresMax[i]) continue;
    for (Int_t ientry=0; ientry<n; ++ientry){
      Double_t value=arr[ientry];
      if (value<fThresMin[i]||value>fThresMax[i]) quality=(QualityFlag_t)i;
    }
  }
  return quality;
}

#endif
