#ifndef ALITRDCALIBRATION_H
#define ALITRDCALIBRATION_H

// macro for extremely simple analysis


#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TList;
class TObject;
class TH1F;
class TProfile2D;
class TH2I;
class TTree;
class TObjArray;
class AliTRDtrackV1;
class AliTRDCalibraFillHisto;
class AliTRDCalibraVLFDebug;
class AliTRDCalibraPRFDebug;
class AliTRDcluster;
class AliTRDtrackInfo;
class AliTRDcalibration : public AliTRDrecoTask 
{
public:
  AliTRDcalibration();
  virtual ~AliTRDcalibration(){;}
  
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetLow(Int_t low)        {flow=low;};
  void SetHigh(Int_t high)      {fhigh=high;};
  void SetFillZero(Bool_t fillZero)     {ffillZero =  fillZero;   };

  
private:
  AliTRDtrackInfo  *fTrackInfo;                  // track info

  AliTRDtrackV1 *ftrdTrack;                      //trdtrack
  AliTRDcluster *fcl;                            //cluster
  
  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;  //calibration analyse object
  TH1F        *fNbTRDTrackUsed;                  //nb ESD tracks used for calibration
  TH1F        *fNbTimeBin;                       //nb Time Bin
  TH1F        *fNbClusters;                      //nb Clusters
  TProfile2D  *fPHSum;                           //sum PH
  TH2I        *fCHSum;                           //sum CH
  
  Int_t       flow;                              //lower limit
  Int_t       fhigh;                             //higher limit
  Int_t       fNbTimeBins;                       //number of timebins 
  Bool_t      ffillZero;                         //fill zero


  AliTRDcalibration(const AliTRDcalibration&); 
  AliTRDcalibration& operator=(const AliTRDcalibration&); 

  ClassDef(AliTRDcalibration, 1) // calibration task
};

#endif
