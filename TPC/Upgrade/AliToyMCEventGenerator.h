#ifndef AliToyMCEventGenerator_H
#define AliToyMCEventGenerator_H

#include <TString.h>

class TFile;
class TTree;
class TObjArray;

class AliTPCParam;
class AliTPCCorrection;
class AliTrackPointArray;
class AliTrackPoint;
class AliTPCclusterMI;

class AliToyMCTrack;
class AliToyMCEvent;

class AliToyMCEventGenerator : public TObject {
 public:
   enum EGasType {
     kNeCO2_9010=0,
     kNeCO2N2_90105
   };

  enum EEpsilon {
    kEps5=0,
    kEps10,
    kEps20,
    kEps25,
    kEps30,
    kEps35,
    kEps40
  };

  enum ECollRate {
    k50kHz = 0
  };

  enum ECorrection {
    kLookup=0,
    kSpaceChargeFile
  };
  
  AliToyMCEventGenerator();
  AliToyMCEventGenerator(const AliToyMCEventGenerator &gen);
  virtual ~AliToyMCEventGenerator();

  virtual AliToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(AliToyMCTrack &trackIn, Double_t t0=0);
  void MakeITSClusters(AliToyMCTrack &trackIn/*, Double_t t0*/);
  void MakeTRDClusters(AliToyMCTrack &trackIn/*, Double_t t0*/);
  void MakeTPCClusters(AliToyMCTrack &trackIn, Double_t t0);
  void CreateSpacePoints(AliToyMCTrack &trackIn,
                        AliTrackPointArray &arrUdist,
                        AliTrackPointArray &arrDist);
  void SetPoint(Float_t xyz[3], Float_t sigmaY, Float_t sigmaZ, AliTrackPoint &point);
  void ConvertTrackPointsToLocalClusters(AliTrackPointArray &arrPoints, AliToyMCTrack &tr, Double_t t0, Int_t type);
  Bool_t SetupCluster(AliTPCclusterMI &tempCl, Float_t xyz[3], Int_t sec, Double_t t0);
  
  void SetOutputFileName(const char* file) { fOutputFileName=file; }
  const char* GetOutputFileName()    const { return fOutputFileName.Data(); }

  void SetSpaceCharge(EEpsilon epsilon, EGasType gasType=kNeCO2_9010, ECollRate collRate=k50kHz, ECorrection corrType=kLookup);
  void SetSpaceChargeFile(const char* file) { fCorrectionFile=file; }
  
  Int_t GetSector(Float_t xyz[3]);

  void InitSpaceCharge();

  void SetStepCorrection(Bool_t step=kTRUE) { fUseStepCorrection=step;   }
  Bool_t GetStepCorrection() const          { return fUseStepCorrection; }

  void SetUseMaterialBudget(Bool_t use) { fUseMaterialBudget=use;    }
  Bool_t GetUseMaterialBudget() const   { return fUseMaterialBudget; }

  void SetIsLaser(Bool_t use) { fIsLaser=use;    }
  Bool_t GetIsLaser() const   { return fIsLaser; }

  void   SetSCListFile(const char* file) { fSCListFile=file;              }
  const char* GetSCListFile() const      { return fSCListFile.Data();     }
  void   SetPrereadSCList(Bool_t b)      { fPrereadSCList=b;              }
  Bool_t GetPrereadSCList() const        { return fPrereadSCList;         }
  Bool_t HasSCList() const               { return  !fSCListFile.IsNull(); }

  void SetCalculateScaling(Bool_t  val) { fCalculateScaling = val; }
  Bool_t  GetCalculateScaling() const { return fCalculateScaling; }

  static Float_t GetSCScalingFactor(AliTPCCorrection *corr, AliTPCCorrection *averageCorr, Float_t &chi2);
  static void SetCorrectionFromFile(TString file, AliTPCCorrection* &corr);
  
 protected:
  AliTPCParam *fTPCParam;                //! TPC params
  AliToyMCEvent *fEvent;                 //! Toy event
  
  Bool_t ConnectOutputFile();
  Bool_t CloseOutputFile();
  void FillTree();
  void IterateSC(Int_t ipos=-1);
  void SetSCScalingFactor();
  
  UInt_t fCurrentTrack;                  // unique track id within the current event generation

  
 private:
  AliToyMCEventGenerator& operator= (const AliToyMCEventGenerator& );
   
  AliTPCCorrection *fTPCCorrection;      //! distortion correction
  AliTPCCorrection *fTPCCorrectionAv;    //! average distortion correction
  
  TObjArray   *fSCList;                  //! list with
  TString fSCListFile;                   // file with a list of space charge files
  TString fCorrectionFile;               // name of a sinfle SC file
  TString fOutputFileName;               // name of the output file
  TFile   *fOutFile;                     //! output file
  TTree   *fOutTree;                     //! output tree

  Bool_t fUseStepCorrection;             // use integralDz method?
  Bool_t fUseMaterialBudget;             // use material budget in tracking?
  Bool_t fIsLaser;                       // is a laser event?
  Bool_t fPrereadSCList;                 // preread all SC files from the SC list
  Bool_t fCalculateScaling;              // calculate scaling factor

  void InitSpaceChargeList();
  
  ClassDef(AliToyMCEventGenerator, 2)
  
};

#endif

