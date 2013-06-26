#ifndef AliToyMCEventGenerator_H
#define AliToyMCEventGenerator_H

#include <TString.h>

class TFile;
class TTree;

class AliTPCParam;
class AliTPCSpaceCharge3D;
class AliTrackPointArray;

class AliToyMCTrack;
class AliToyMCEvent;

class AliToyMCEventGenerator : public TObject {
 public:
   enum EGasType {
     kNeCO2_9010=0
  };

  enum EEpsilon {
    kEps5=0,
    kEps10,
    kEps20
  };

  enum ECollRate {
    k50kHz = 0
  };
  
  AliToyMCEventGenerator();
  AliToyMCEventGenerator(const AliToyMCEventGenerator &gen);
  virtual ~AliToyMCEventGenerator();

  virtual AliToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(AliToyMCTrack &trackIn, Double_t t0=0);
  void CreateSpacePoints(AliToyMCTrack &trackIn,
                        AliTrackPointArray &arrUdist,
                        AliTrackPointArray &arrDist);
  void SetPoint(Float_t xyz[3], AliTrackPoint &point);
  void ConvertTrackPointsToLocalClusters(AliTrackPointArray &arrPoints, AliToyMCTrack &tr, Double_t t0, Int_t type);
  Bool_t SetupCluster(AliTPCclusterMI &tempCl, Float_t xyz[3], Int_t sec, Double_t t0);
  
  void SetOutputFileName(const char* file) { fOutputFileName=file; }
  const char* GetOutputFileName()    const { return fOutputFileName.Data(); }

  void SetSpaceCharge(EEpsilon epsilon, EGasType gasType=kNeCO2_9010, ECollRate collRate=k50kHz);
  void SetSpaceChargeFile(const char* file) { fSpaceChargeFile=file; }
  
  Int_t GetSector(Float_t xyz[3]);

  void InitSpaceCharge();

  void SetStepCorrection(Bool_t step=kTRUE) { fUseStepCorrection=step;   }
  Bool_t GetStepCorrection() const          { return fUseStepCorrection; }

  void SetUseMaterialBudget(Bool_t use) { fUseMaterialBudget=use;    }
  Bool_t GetUseMaterialBudget() const   { return fUseMaterialBudget; }
  
 protected:
  AliTPCParam *fTPCParam;
  AliToyMCEvent *fEvent;
  
  Bool_t ConnectOutputFile();
  Bool_t CloseOutputFile();
  void FillTree();

  
 private:
  AliToyMCEventGenerator& operator= (const AliToyMCEventGenerator& );
   
  AliTPCSpaceCharge3D *fSpaceCharge;

  TString fSpaceChargeFile;
  TString fOutputFileName;
  TFile   *fOutFile;
  TTree   *fOutTree;

  Bool_t fUseStepCorrection;
  Bool_t fUseMaterialBudget;
  
  ClassDef(AliToyMCEventGenerator, 1)
     
};

#endif

