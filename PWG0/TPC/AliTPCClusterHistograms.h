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
  void         SetCommentToHistograms(const Char_t* text) {fCommentToHistograms = TString(text);}

  void FillEvent(Int_t time, Int_t nTracks);
  void FillTrack(const AliTPCseed* seed);
  void FillCluster(AliTPCclusterMI* clusterMI, Int_t time=-1);
  
  void    StartEvent();
  void    FinishEvent(Int_t timeStamp);
  
  Int_t   GetNClusters() {return fNClustersInEvent;}
  Float_t GetQtotInEvent() {return fQtotInEvent;}
  Float_t GetMaxQtotInEvent() {return fMaxQtotInEvent;}

  Bool_t  KeepThisEvent(TString& why);
  const Char_t* WhyKeepEvent() {return fWhyKeepEvent.Data();}

  Bool_t  IsClusterOnEdge(AliTPCclusterMI* clusterMI);
  Float_t DistanceToEdge(AliTPCclusterMI* clusterMI);
 
  TCanvas* DrawHistograms(const Char_t* opt="");

  static const char* FormDetectorName(Int_t detector, Bool_t edgeSuppression = kFALSE, const char* comment = 0);

protected:

  Int_t       fTimeStart;               // begin time of run(s)
  Int_t       fTimeStop;                // end time of runs(s)

  TH2F*       fhQmaxVsRow;              //        QmaxVsRow
  TH2F*       fhQtotVsRow;              //	  QtotVsRow
						  					    
  TProfile*   fhQtotProfileVsRow;       //	  QtotProfileVsRow
  TProfile*   fhQmaxProfileVsRow;       //	  QmaxProfileVsRow

  TH2F*       fhQtotVsDistanceToEdge;        // qtot vs distance to edge
  TProfile*   fhQtotProfileVsDistanceToEdge; // qtot vs distance to edge

  TH2F*       fhNClustersYVsRow;        //        n clusters y vs row
  TH2F*       fhNClustersZVsRow;        //        n clusters z vs row
						  					    
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
						  					    
  TProfile*   fhMeanQtotVsTime;         //	  mean qtot vs time
  TH2F*       fhQtotVsTime;             //	  qtot vs time
  TProfile*   fhMeanNClustersVsTime;    //        mean number of clusters vs time
  TH2F*       fhNClustersVsTime;        //        number of clusters vs time
						  					    
  TH1F*       fhTrackQtotPerCluster;    //	  TrackQtotPerCluster
						  					    
  TH2F*       fhTrackQtotPerClusterVsPhi; //	  TrackQtotPerClusterVsPhi
  TH2F*       fhTrackQtotPerClusterVsTheta; //	  TrackQtotPerClusterVsTheta
						  					    
  TProfile*   fhTrackMeanQtotPerClusterVsPhi; //  TrackMeanQtotPerClusterVsPhi
  TProfile*   fhTrackMeanQtotPerClusterVsTheta; //  TrackMeanQtotPerClusterVsTheta

  TProfile*   fhMeanNTracksVsTime;        // mean number of tracks vs time   
  TH1F*       fhNEventsVsTime;            // number of events vs time

  Int_t       fDetector;                // number of detector
  Bool_t      fIsIROC;                  // true = IROC, false = OROC
  Bool_t      fEdgeSuppression;         // if set edges are not taken into account for histograms

  Int_t       fNClustersInEvent;        // number of clusters in event
  Float_t     fQtotInEvent;             // total qtot in event
  Float_t     fMaxQtotInEvent;          // max qtot in event

  Bool_t      fKeepEvent;               // keep this event
  TString     fWhyKeepEvent;            // why

  TString     fCommentToHistograms;     // comments to histograms

  ClassDef(AliTPCClusterHistograms,1)
};

#endif

