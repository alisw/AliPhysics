#ifndef ALITPCClusterHistograms_H
#define ALITPCClusterHistograms_H

/* $Id$ */

//
// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the reconstructed clusters.
//

#include "TNamed.h"

class TH3;
class TH2F;
class TCanvas;
class TProfile;
class TProfile2D;
class TString;

class AliTPCclusterMI;

class AliTPCClusterHistograms : public TNamed
{
public:

  AliTPCClusterHistograms();
  AliTPCClusterHistograms(Int_t detector, const Char_t* comment="", Int_t timeStart=-1, Int_t timeStop=-1);
  
  AliTPCClusterHistograms(const AliTPCClusterHistograms& c);
  virtual ~AliTPCClusterHistograms();
  AliTPCClusterHistograms& operator=(const AliTPCClusterHistograms& corrMatrix);

  virtual Long64_t Merge(TCollection* list);

  virtual void SaveHistograms();

  void FillCluster(AliTPCclusterMI* clusterMI, Int_t time=-1);
  //  void FillCluster(AliTPCclusterMI* clusterMI, Int_t time=-1, Float_t trackangle);

  TCanvas* DrawHistograms(const Char_t* opt="");

protected:

  Int_t       fTimeStart;               // begin time of run(s)
  Int_t       fTimeStop;                // end time of runs(s)

  TH2F*       fhQmaxVsRow;              //
  TH2F*       fhQtotVsRow;              //
  TH2F*       fhSigmaYVsRow;            //
  TH2F*       fhSigmaZVsRow;            //

  TProfile2D* fhQmaxProfileYVsRow;      //
  TProfile2D* fhQtotProfileYVsRow;      // 
  TProfile2D* fhSigmaYProfileYVsRow;    //
  TProfile2D* fhSigmaZProfileYVsRow;    //

  TProfile2D* fhQmaxProfileZVsRow;      //
  TProfile2D* fhQtotProfileZVsRow;      //
  TProfile2D* fhSigmaYProfileZVsRow;    //
  TProfile2D* fhSigmaZProfileZVsRow;    //

  TProfile*   fhQtotVsTime;             //
  TProfile*   fhQmaxVsTime;             //

  ClassDef(AliTPCClusterHistograms,1)
};

#endif

