#ifndef ALITPCClusterHistograms_H
#define ALITPCClusterHistograms_H

/* $Id$ */

//
// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the reconstructed clusters.
//

#include <TNamed.h>

class TH3;
class TH2F;
class TH1F;
class TCanvas;
class TProfile;
class TProfile2D;
class TObjArray;
class TString;

class AliTPCclusterMI;
class AliTPCseed;

class AliTPCClusterHistograms : public TNamed
{
public:

  AliTPCClusterHistograms();
  AliTPCClusterHistograms(Int_t detector, const Char_t* comment="", Int_t timeStart=-1, Int_t timeStop=-1, Bool_t edgeSuppression = kFALSE);
  
  AliTPCClusterHistograms(const AliTPCClusterHistograms& c);
  virtual ~AliTPCClusterHistograms();
  AliTPCClusterHistograms& operator=(const AliTPCClusterHistograms& corrMatrix);

  virtual Long64_t Merge(TCollection* list);

  virtual void SaveHistograms();

  void FillTrack(const AliTPCseed* seed);
  void FillCluster(AliTPCclusterMI* clusterMI, Int_t time=-1);

  Bool_t IsClusterOnEdge(AliTPCclusterMI* clusterMI);

  TCanvas* DrawHistograms(const Char_t* opt="");

  static const char* FormDetectorName(Int_t detector, Bool_t edgeSuppression = kFALSE, const char* comment = 0);

protected:

  Int_t       fTimeStart;               // begin time of run(s)
  Int_t       fTimeStop;                // end time of runs(s)

  TH2F*       fhQmaxVsRow;              //        QmaxVsRow
  TH2F*       fhQtotVsRow;              //	  QtotVsRow
						  					    
  TProfile*   fhQtotProfileVsRow;       //	  QtotProfileVsRow
  TProfile*   fhQmaxProfileVsRow;       //	  QmaxProfileVsRow

  TH2F*       fhNClustersYVsRow;        // 
  TH2F*       fhNClustersZVsRow;        // 
						  					    
  TH2F*       fhSigmaYVsRow;            //	  SigmaYVsRow
  TH2F*       fhSigmaZVsRow;            //	  SigmaZVsRow
						  					    
  TProfile2D* fhQmaxProfileYVsRow;      //	  QmaxProfileYVsRow
  TProfile2D* fhQtotProfileYVsRow;      // 	  QtotProfileYVsRow
  TProfile2D* fhSigmaYProfileYVsRow;    //	  SigmaYProfileYVsRow
  TProfile2D* fhSigmaZProfileYVsRow;    //	  SigmaZProfileYVsRow
						  					    
  TProfile2D* fhQmaxProfileZVsRow;      //	  QmaxProfileZVsRow
  TProfile2D* fhQtotProfileZVsRow;      //	  QtotProfileZVsRow
  TProfile2D* fhSigmaYProfileZVsRow;    //	  SigmaYProfileZVsRow
  TProfile2D* fhSigmaZProfileZVsRow;    //	  SigmaZProfileZVsRow
						  					    
  TProfile*   fhQtotVsTime;             //	  QtotVsTime
  TProfile*   fhQmaxVsTime;             //	  QmaxVsTime
						  					    
  TH1F*       fhTrackQtotPerCluster;    //	  TrackQtotPerCluster
						  					    
  TH2F*       fhTrackQtotPerClusterVsSnp; //	  TrackQtotPerClusterVsSnp
  TH2F*       fhTrackQtotPerClusterVsTgl; //	  TrackQtotPerClusterVsTgl
						  					    
  TProfile*   fhTrackMeanQtotPerClusterVsSnp; //  TrackMeanQtotPerClusterVsSnp
  TProfile*   fhTrackMeanQtotPerClusterVsTgl; //  TrackMeanQtotPerClusterVsTgl
 

  Int_t       fDetector;                // number of detector
  Bool_t      fIsIROC;                  // true = IROC, false = OROC
  Bool_t      fEdgeSuppression;         // if set edges are not taken into account for histograms

  ClassDef(AliTPCClusterHistograms,1)
};

#endif

