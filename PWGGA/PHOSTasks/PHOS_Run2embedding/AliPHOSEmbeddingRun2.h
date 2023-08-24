#ifndef AliPHOSEmbeddingRun2_h
#define AliPHOSEmbeddingRun2_h

// Class to perform embedding on the AOD level
// Author: D.Peressounko

#include "AliAnalysisTaskSE.h"
class TChain;
class TClonesArray;
class TH2F;

class AliPHOSClusterizerv1;
class AliPHOSReconstructor;
class AliAODEvent;
class AliAODtrack;
class AliAODCaloCells;
class AliPHOSGeometry;
class AliPHOSCalibData;

class AliPHOSEmbeddingRun2 : public AliAnalysisTaskSE
{
 public:
  AliPHOSEmbeddingRun2(const char* name = "AliPHOSEmbeddingRun2");
  virtual ~AliPHOSEmbeddingRun2() {}

  // Standard methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*) {}

  // Chain with signal AOD for embedding
  void SetSignalChain(TChain* signal) { fAODChain = signal; }
  void SetPrivateOADBPath(TString path = "$ALICE_PHYSICS/OADB/PHOS/PHOSMCCalibrations.root")
  {
    fPathPrivateOADBMC = path;
  } // simply copy from AliPhysics which you want to use for analysis is enough.
  void SetSignalCalibration(double corr) { fSignalECorrection = corr; }
  void SetSelectedEvProp(float p = 0.2) { fSelEventsPart = p; } // proportion of events which will pass Ev.Selection

 private:
  AliPHOSEmbeddingRun2(const AliPHOSEmbeddingRun2&);            // not implemented
  AliPHOSEmbeddingRun2& operator=(const AliPHOSEmbeddingRun2&); // not implemented

  void RunSignalSimulation(int nevents);
  void Init();
  void InitMF(); // Mag.Field initialization for track matching

  bool GetNextSignalEvent(void); // returns true is successfully read

  void CopyRecalibrateSignal();
  void CopyRecalibrateBackground();
  void MakeEmbedding();

  // void CopyRecalibrateDigits() ;
  void MakeDigits(const AliAODEvent* bgevent, const AliAODEvent* signal);

  // Add new branches
  void ConvertEmbeddedClusters();
  void ConvertEmbeddedCells();
  void ConvertMCParticles();

  double CorrectNonlinearity(double en);
  double CorrectNonlinearityMC(double en);
  bool IsGoodChannel(int mod, int ix, int iz);
  double CoreEnergy(AliVCluster* clu);
  void EvalLambdas(AliVCluster* clu, double& m02, double& m20);
  double TestCoreLambda(double pt, double l1, double l2);
  double TestFullLambda(double pt, double l1, double l2);
  // void DistanceToBadChannel(int mod, TVector3 * locPos, double &minDist) ;
  double EvalTOF(AliVCluster* clu, AliVCaloCells* cells);

  int FindTrackMatching(int mod, TVector3* locpos, double& dx, double& dz, double& pttrack, int& charge);
  float TestCPV(double dx, double dz, double pt, int charge);
  float TestCPVRun2(double dx, double dz, double pt, int charge);

 private:
  bool fAddNoiseMC = false;
  float fNoiseMC = 0.001;
  float fZScut = 0.;
  float fSelEventsPart = 0.2; // only fSelEventsPart will pass selection, simulate smaller number of events

  TChain* fAODChain;    //! Signal
  AliAODEvent* fSignal; //!

  TTree* fDigitsTree;       //! Digits
  TTree* fClustersTree;     //! Clusters
  TTree* fTreeOut;          //! Output AOD
  TClonesArray* fDigitsArr; //!

  TClonesArray* fmcparticles;      //!
  TClonesArray* fSignalClusters;   //!
  TClonesArray* fEmbeddedClusters; //!
  AliAODCaloCells* fEmbeddedCells; //!
  AliAODCaloCells* fCellsPHOS;     //! Old PHOS cells
  TObjArray* fEmcRecPoints;        //!
  TObjArray* fCpvRecPoints;        //!

  AliPHOSGeometry* fPHOSGeo;                //!
  AliPHOSClusterizerv1* fClusterizer;       //!
  AliPHOSReconstructor* fPHOSReconstructor; //!

  TH2I* fPHOSBadMap[6];   //! Bad channels map
  float fRunByRunCorr[5]; // Per module run-by-run correction
  int fL1phase[15];       // L1phases for PHOS DDLs (run2 only)

  TH2F* fOldPHOSCalibration[5];            //! Calibration coeff. used in Bg production
  AliPHOSCalibData* fCalibMC = nullptr;    //! Calibration of signal
  AliPHOSCalibData* fCalibData = nullptr;  //! Calibration of Data
  AliPHOSCalibData* fgCalibData = nullptr; //!
  TString fPathPrivateOADBMC;              // path to private OADB.
  double fSignalECorrection;               // Correction for the Signal clibration
  int fNSignal;                            //! Number of signal evetns processed
  int fNCaloClustersOld;                   //! Number of CaloClusters already in Bg
  int fRunNumber;                          //! Current run
  int fMF;                                 //! Magnetic field
  bool fInitialized;                       //!
  TVector3 fVtx;                           //! vertex in current event

  ClassDef(AliPHOSEmbeddingRun2, 7); // PHOS analysis task
};

#endif
