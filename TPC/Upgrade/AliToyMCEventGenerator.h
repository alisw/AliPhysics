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
  AliToyMCEventGenerator();
  AliToyMCEventGenerator(const AliToyMCEventGenerator &gen);
  virtual ~AliToyMCEventGenerator();

  virtual AliToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(AliToyMCTrack &trackIn, Double_t t0);
  void CreateSpacePoints(AliToyMCTrack &trackIn,
                        AliTrackPointArray &arrUdist,
                        AliTrackPointArray &arrDist);
  void SetPoint(Float_t xyz[3], AliTrackPoint &point);
  void ConvertTrackPointsToLocalClusters(AliTrackPointArray &arrPoints, AliToyMCTrack &tr, Double_t t0, Int_t type);
  Bool_t SetupCluster(AliTPCclusterMI &tempCl, Float_t xyz[3], Int_t sec, Double_t t0);
  
  void SetOutputFileName(const char* file) { fOutputFileName=file; }
  const char* GetOutputFileName()    const { return fOutputFileName.Data(); }

  Int_t GetSector(Float_t xyz[3]);
  
 protected:
  AliTPCParam *fTPCParam;
  AliToyMCEvent *fEvent;
  
  Bool_t ConnectOutputFile();
  Bool_t CloseOutputFile();
  void FillTree();
  
 private:
  AliToyMCEventGenerator& operator= (const AliToyMCEventGenerator& );
   
  AliTPCSpaceCharge3D *fSpaceCharge;

  TString fOutputFileName;
  TFile   *fOutFile;
  TTree   *fOutTree;
  
  ClassDef(AliToyMCEventGenerator, 1)
     
};

#endif

