#ifndef ALITPCClusterHistograms_H
#define ALITPCClusterHistograms_H

/* $Id$ */


#include <TNamed.h>
//#include <../TPC/AliTPCclusterMI.h>

class TH3;
class TH2F;
class TProfile2D;
class AliTPCclusterMI;

class AliTPCClusterHistograms : public TNamed
{
public:

  AliTPCClusterHistograms();
  AliTPCClusterHistograms(const Char_t* name, const Char_t* title);
  AliTPCClusterHistograms(const AliTPCClusterHistograms& c);
  virtual ~AliTPCClusterHistograms();
  AliTPCClusterHistograms& operator=(const AliTPCClusterHistograms& corrMatrix);

  virtual Long64_t Merge(TCollection* list);

  virtual void SaveHistograms();

  void FillCluster(AliTPCclusterMI* clusterMI);

protected:

  TH2F*       fhQmaxVsRow;
  TH2F*       fhQtotVsRow;
  TH2F*       fhSigmaYVsRow;
  TH2F*       fhSigmaZVsRow;

  TProfile2D* fhQmaxProfileYVsRow;
  TProfile2D* fhQtotProfileYVsRow;
  TProfile2D* fhSigmaYProfileYVsRow;
  TProfile2D* fhSigmaZProfileYVsRow;


  ClassDef(AliTPCClusterHistograms,1)
};

#endif

