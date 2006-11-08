#ifndef ALITPCClusterHistograms_H
#define ALITPCClusterHistograms_H

/* $Id$ */


#include <TNamed.h>
//#include <../TPC/AliTPCclusterMI.h>

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
  //AliTPCClusterHistograms(const Char_t* name, const Char_t* title);
  AliTPCClusterHistograms(Int_t detector, const Char_t* comment="", Int_t timeStart=-1, Int_t timeStop=-1);
  
  AliTPCClusterHistograms(const AliTPCClusterHistograms& c);
  virtual ~AliTPCClusterHistograms();
  AliTPCClusterHistograms& operator=(const AliTPCClusterHistograms& corrMatrix);

  virtual Long64_t Merge(TCollection* list);

  virtual void SaveHistograms();

  void FillCluster(AliTPCclusterMI* clusterMI, Int_t time=-1);

  TCanvas* DrawHistograms(const Char_t* opt="");

protected:

  Int_t       fTimeStart;
  Int_t       fTimeStop;

  TH2F*       fhQmaxVsRow;
  TH2F*       fhQtotVsRow;
  TH2F*       fhSigmaYVsRow;
  TH2F*       fhSigmaZVsRow;

  TProfile2D* fhQmaxProfileYVsRow;
  TProfile2D* fhQtotProfileYVsRow;
  TProfile2D* fhSigmaYProfileYVsRow;
  TProfile2D* fhSigmaZProfileYVsRow;

  TProfile2D* fhQmaxProfileZVsRow;
  TProfile2D* fhQtotProfileZVsRow;
  TProfile2D* fhSigmaYProfileZVsRow;
  TProfile2D* fhSigmaZProfileZVsRow;

  TProfile*   fhQtotVsTime;
  TProfile*   fhQmaxVsTime;

  ClassDef(AliTPCClusterHistograms,1)
};

#endif

