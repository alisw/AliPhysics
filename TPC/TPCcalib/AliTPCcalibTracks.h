#ifndef AliTPCCALIBTRACKS_H
#define AliTPCCALIBTRACKS_H


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     Class to analyse tracks for calibration                               //
//    to be used as a component in AliTPCSelectorTracks                     //
//     In the constructor you have to specify name and title                 //
//     to get the Object out of a file.                                      //
//     The parameter 'clusterParam', a AliTPCClusterParam object             //
//      (needed for TPC cluster error and shape parameterization)            //
//     Normally you get this object out of the file 'TPCClusterParam.root'   //
//     In the parameter 'cuts' the cuts are specified, that decide           //
//     weather a track will be accepted for calibration or not.              //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include <AliTPCcalibBase.h>
#include "THnSparse.h"

class TF2;
class TH3F;
class TH1F;
class TH1I;
class TH2I;
class TH2D;
class TCollection;
class TTreeSRedirector;
class TLinearFitter;
class AliTPCClusterParam;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
//class AliESDtrack;
class AliTPCclusterMI;
class AliTPCcalibTracksCuts;
class AliTPCCalPad;
class TChain;
class TTree;
class TMutex;
//class AliESDEvent;
class AliVEvent;
class AliVTrack;

using namespace std;

class AliTPCcalibTracks : public AliTPCcalibBase {
public :
   AliTPCcalibTracks();                         // default constructor
  AliTPCcalibTracks(const AliTPCcalibTracks&calibTracks); // copy constructor
  AliTPCcalibTracks(const Text_t *name, const Text_t *title, AliTPCClusterParam *clusterParam, AliTPCcalibTracksCuts* cuts, Int_t logLevel = 0);
  AliTPCcalibTracks & operator=(const AliTPCcalibTracks& calibTracks);
  
  virtual ~AliTPCcalibTracks();                // destructor
  
  virtual void            Process(AliTPCseed *track);  // to be called by the Selector
  //void     Process(AliESDEvent *event) {AliTPCcalibBase::Process(event);}
  void     Process(AliVEvent *event) {AliTPCcalibBase::Process(event);}
  //void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  virtual Long64_t Merge(TCollection *li);
  void    AddHistos(AliTPCcalibTracks* calib);
  void     MakeResPlotsQTree(Int_t minEntries = 100, const char* pathName = "plots");
   void     MakeReport(Int_t stat, const char* pathName = "plots");     // calls all functions that procude pictures, results are written to pathName, stat is the minimal statistic threshold
  //   
  Int_t           AcceptTrack(AliTPCseed * track);
  void            FillResolutionHistoLocal(AliTPCseed * track);  // the MAIN-FUNCTION, called for each track to fill the histograms, called by Process(...)
  
  
  void     SetStyle() const;
  
  TObjArray* GetfArrayQDY() const {return fArrayQDY;}
  TObjArray* GetfArrayQDZ() const {return fArrayQDZ;}
  TObjArray* GetfArrayQRMSY() const {return fArrayQRMSY;}
  TObjArray* GetfArrayQRMSZ() const {return fArrayQRMSZ;}
  TObjArray* GetfResolY() const {return fResolY;}
  TObjArray* GetfResolZ() const {return fResolZ;}
  TObjArray* GetfRMSY() const {return fRMSY;}
  TObjArray* GetfRMSZ() const {return fRMSZ;}
  TH1I*      GetfRejectedTracksHisto() const {return fRejectedTracksHisto;}
  TH2I*      GetfClusterCutHisto() const  {return fClusterCutHisto;}
  AliTPCCalPad*          GetfCalPadClusterPerPad() const {return fCalPadClusterPerPad; }
  AliTPCCalPad*          GetfCalPadClusterPerPadRaw() const {return fCalPadClusterPerPadRaw;}
  AliTPCcalibTracksCuts* GetCuts() {return fCuts;}
  void MakeHistos();  //make THnSparse
  int UpdateClusterParam( AliTPCClusterParam *cParam, Bool_t MirrorZ=1, Bool_t MirrorPad=1, Bool_t MirrorAngle=1, Int_t MinStat=10 );

  static void MakeSummaryTree(THnSparse *hisInput, TTreeSRedirector *pcstream, Int_t ptype);
  static int GetTHnStat( const  THnBase *H, THnBase *&Mean, THnBase *&Sigma, THnBase *&Entr );
  static int CreateWaveCorrection( const  THnBase *DeltaY, THnBase *&MeanY, THnBase *&SigmaY, THnBase *&EntrY,
				   Bool_t MirrorZ=1, Bool_t MirrorPad=1, Bool_t MirrorAngle=1, Int_t MinStat=10 );
 
  static void SetMergeEntriesCut(Double_t entriesCut){fgkMergeEntriesCut = entriesCut;}

protected:         
  
private:

   static Int_t   GetBin(Float_t q, Int_t pad);
   static Int_t   GetBin(Int_t  iq, Int_t pad);
   static Float_t GetQ(Int_t bin);
   static Float_t GetPad(Int_t bin);
   AliTPCClusterParam *fClusterParam; // pointer to cluster parameterization
   AliTPCROC *fROC;          //!
public:
  THnSparse  *fHisDeltaY;    // THnSparse - delta Y 
  THnSparse  *fHisDeltaZ;    // THnSparse - delta Z 
  THnSparse  *fHisRMSY;      // THnSparse - rms Y 
  THnSparse  *fHisRMSZ;      // THnSparse - rms Z 
  THnSparse  *fHisQmax;      // THnSparse - qmax 
  THnSparse  *fHisQtot;      // THnSparse - qtot 

private:
  Double_t fPtDownscaleRatio;       // pt downscaling ratio (use subsample of data)
  Double_t fQDownscaleRatio;        // Q downscaling ratio (use subsample of dta)

   TObjArray *fArrayQDY;    // q binned delta Y histograms
   TObjArray *fArrayQDZ;    // q binned delta Z histograms 
   TObjArray *fArrayQRMSY;  // q binned delta Y histograms
   TObjArray *fArrayQRMSZ;  // q binned delta Z histograms 

   TObjArray *fResolY;      // array of resolution histograms Y
   TObjArray *fResolZ;      // array of resolution histograms Z
   TObjArray *fRMSY;        // array of RMS histograms Y
   TObjArray *fRMSZ;        // array of RMS histograms Z
   AliTPCcalibTracksCuts *fCuts; // object with cuts, that is passed to the constructor
   TH1I      *fRejectedTracksHisto; // histogram of rejecteced tracks, the number coresponds to the failed cut
   TH2I      *fClusterCutHisto;     // histogram showing in which padRow the clusters were cutted by which criterium
   AliTPCCalPad *fCalPadClusterPerPad;    // AliTPCCalPad showing the number of clusters per Pad
   AliTPCCalPad *fCalPadClusterPerPadRaw; // AliTPCCalPad showing the number of clusters per Pad before cuts on clusters are applied
   static Double_t            fgkMergeEntriesCut;//maximal number of entries for merging  -can be modified via setter

  ClassDef(AliTPCcalibTracks,2)
   
};

#endif
