#ifndef ALIMUONALIGNMENTTASK_H
#define ALIMUONALIGNMENTTASK_H

/// \ingroup ""
/// \class AliMUONAlignmentTask
/// \brief Task to align the muon spectrometer
///
//  Author Javier Castillo, CEA/Saclay - Irfu/SPhN

class TList;
class TGraphErrors;
class AliESDEvent;
class AliMUONAlignment;
class AliMUONGeoemetryTransformer;

#include "AliAnalysisTask.h"

class AliMUONAlignmentTask : public AliAnalysisTask {
 public:
  //  AliMUONAlignmentTask(const char *name = "AliMUONAlignmentTask");
  AliMUONAlignmentTask(const char *name = "AliMUONAlignmentTask", const char *geofilename = "geometry.root");
  AliMUONAlignmentTask(const AliMUONAlignmentTask& obj);
  AliMUONAlignmentTask& operator=(const AliMUONAlignmentTask& other); 
  virtual ~AliMUONAlignmentTask();
  
  virtual void   LocalInit();
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
 private:
  AliESDEvent *fESD;                      //!< ESD object

  AliMUONAlignment *fAlign;               ///< The MUON alignment object
    TString fGeoFilename;                 ///< Geometry file name
  AliMUONGeometryTransformer *fTransform; ///< MUON geometry transformer
    
  Int_t fTrackTot;             ///< Number of track read 
  Int_t fTrackOk;              ///< Number of track read 

  Double_t fParameters[4*156]; ///< Array of alignment parameters
  Double_t fErrors[4*156];     ///< Array of alignment parameters errors
  Double_t fPulls[4*156];      ///< Array of alignment parameters pulls

  TGraphErrors *fMSDEx ;   ///< Graph of translations along x
  TGraphErrors *fMSDEy ;   ///< Graph of translations along y
  TGraphErrors *fMSDEz ;   ///< Graph of translations along z
  TGraphErrors *fMSDEp;    ///< Graph of rotation about z 

  TList   *fList;          ///< list of graphs
   
  ClassDef(AliMUONAlignmentTask, 1); // example of analysis
};

#endif

