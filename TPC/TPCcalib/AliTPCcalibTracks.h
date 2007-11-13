#ifndef AliTPCCALIBTRACKS_H
#define AliTPCCALIBTRACKS_H


#include <TNamed.h>
#include <TH2D.h>
#include <TF2.h>
#include <TSystem.h>
#include <TCollection.h>
#include <iostream>
#include <TLinearFitter.h>
class AliTPCClusterParam;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
class AliESDtrack;
class AliTPCclusterMI;
class TH3F;
class TH1F;
class TH1I;
class TTreeSRedirector;
class AliTPCcalibTracksCuts;
class AliTPCCalPadRegion;
class AliTPCCalPad;
class TChain;
class TMutex;

using namespace std;

class AliTPCcalibTracks : public TNamed {
public :
   AliTPCcalibTracks();                         // default constructor
   AliTPCcalibTracks(AliTPCcalibTracks* ct);    // copy constructor
   AliTPCcalibTracks(const Text_t *name, const Text_t *title, AliTPCClusterParam *clusterParam, AliTPCcalibTracksCuts* cuts, Int_t logLevel = 0);
   virtual ~AliTPCcalibTracks();                // destructor
   
   static void     AddInfo(TChain *chain, char *fileName);        // add clusterParametrization as user info to the chain
   void            Process(AliTPCseed *track, AliESDtrack *esd);  // to be called by the Selector
   
   Int_t           AcceptTrack(AliTPCseed * track);
   void            FillResolutionHistoLocal(AliTPCseed * track);  // the MAIN-FUNCTION, called for each track to fill the histograms, called by Process(...)
   static  TH2D*   MakeDiff(TH2D * hfit, TF2 * func);

   Long64_t Merge(TCollection *li);
   static AliTPCcalibTracks* TestMerge(AliTPCcalibTracks *ct, AliTPCClusterParam *clusterParam, Int_t nCalTracks = 50);
   
   void     Draw(Option_t* opt);                                  // draw some exemplaric histograms for fast result-check
   void     MakeReport(Int_t stat, char* pathName = "plots");     // calls all functions that procude pictures, results are written to pathName, stat is the minimal statistic threshold
   void     MakeAmpPlots(Int_t stat, char* pathName = "plots");
   void     MakeDeltaPlots(char* pathName = "plots");
   void     MakeChargeVsDriftLengthPlotsOld(char* pathName = "plots");
   void     MakeChargeVsDriftLengthPlots(char* pathName = "plots");
   void     FitResolutionNew(char* pathName = "plots");
   void     FitRMSNew(char* pathName = "plots");
   void     MakeResPlotsQTree(Int_t minEntries = 100, char* pathName = "plots"); 
   void     MakeResPlotsQTreeThread(Int_t minEntries = 100, char* pathName = "plots"); 
   static void* MakeResPlotsQTreeThreadFunction(void* arg); 
   void     SetStyle();
   
//protected:   
   TObjArray* GetfArrayAmpRow() {return fArrayAmpRow;}
   TObjArray* GetfArrayAmp() {return fArrayAmp;}
   TObjArray* GetfArrayQDY() {return fArrayQDY;}
   TObjArray* GetfArrayQDZ() {return fArrayQDZ;}
   TObjArray* GetfArrayQRMSY() {return fArrayQRMSY;}
   TObjArray* GetfArrayQRMSZ() {return fArrayQRMSZ;}
   TObjArray* GetfArrayChargeVsDriftlength() {return fArrayChargeVsDriftlength;}
   TH1F*      GetfDeltaY() {return fDeltaY;}
   TH1F*      GetfDeltaZ() {return fDeltaZ;}
   TObjArray* GetfResolY() {return fResolY;}
   TObjArray* GetfResolZ() {return fResolZ;}
   TObjArray* GetfRMSY() {return fRMSY;}
   TObjArray* GetfRMSZ() {return fRMSZ;}
   TH1I*      GetfHclus() {return fHclus;}
   TH1I*      GetfRejectedTracksHisto() {return fRejectedTracksHisto;}
   TH1I*      GetfHclusterPerPadrow() {return fHclusterPerPadrow;}
   TH1I*      GetfHclusterPerPadrowRaw() {return fHclusterPerPadrowRaw;}
   TH2I*      GetfClusterCutHisto() {return fClusterCutHisto;}
   AliTPCCalPad*          GetfCalPadClusterPerPad() {return fCalPadClusterPerPad; }
   AliTPCCalPad*          GetfCalPadClusterPerPadRaw() {return fCalPadClusterPerPadRaw;}
   AliTPCCalPadRegion*    GetCalPadRegionchargeVsDriftlength() {return fcalPadRegionChargeVsDriftlength;}
   AliTPCcalibTracksCuts* GetCuts() {return fCuts;}
   void       SetLogLevel(Int_t level) {fDebugLevel = level;}
   Int_t      GetLogLevel() {return fDebugLevel;}

   
   
private:
   static Int_t   GetBin(Float_t q, Int_t pad);
   static Int_t   GetBin(Int_t  iq, Int_t pad);
   static Float_t GetQ(Int_t bin);
   static Float_t GetPad(Int_t bin);
   void FillResolutionHistoLocalDebugPart(AliTPCseed *track, AliTPCclusterMI *cluster0, Int_t irow, Float_t  angley, Float_t  anglez, Int_t nclFound, Int_t kDelta);
   AliTPCClusterParam *fClusterParam; //! pointer to cluster parameterization
   TTreeSRedirector   *fDebugStream;  //! debug stream for
   AliTPCROC *fROC;          //!
   TObjArray *fArrayAmpRow; // array with amplitudes versus row for given sector 
   TObjArray *fArrayAmp;    // array with amplitude for sectors
   TObjArray *fArrayQDY;    // q binned delta Y histograms
   TObjArray *fArrayQDZ;    // q binned delta Z histograms 
   TObjArray *fArrayQRMSY;  // q binned delta Y histograms
   TObjArray *fArrayQRMSZ;  // q binned delta Z histograms 
   TObjArray *fArrayChargeVsDriftlength; // array of arrays of TProfiles with charge vs. driftlength for each padsize and sector
   AliTPCCalPadRegion *fcalPadRegionChargeVsDriftlength;
   TH1F      *fDeltaY;      // integrated delta y histo
   TH1F      *fDeltaZ;      // integrated delta z histo
   TObjArray *fResolY;      // array of resolution histograms Y
   TObjArray *fResolZ;      // array of resolution histograms Z
   TObjArray *fRMSY;        // array of RMS histograms Y
   TObjArray *fRMSZ;        // array of RMS histograms Z
   AliTPCcalibTracksCuts *fCuts; // object with cuts, that is passed to the constructor
   TH1I      *fHclus;       // number of clusters per track
   TH1I      *fRejectedTracksHisto; // histogram of rejecteced tracks, the number coresponds to the failed cut
   TH1I      *fHclusterPerPadrow;   // histogram showing the number of clusters per padRow
   TH1I      *fHclusterPerPadrowRaw;// histogram showing the number of clusters per padRow before cuts on clusters are applied
   TH2I      *fClusterCutHisto;     // histogram showing in which padRow the clusters were cutted by which criterium
   AliTPCCalPad *fCalPadClusterPerPad;    // AliTPCCalPad showing the number of clusters per Pad
   AliTPCCalPad *fCalPadClusterPerPadRaw; // AliTPCCalPad showing the number of clusters per Pad before cuts on clusters are applied
   Int_t      fDebugLevel;  // log level - debug output: -1: silence, 0: default, 1: things like constructor called, 5: write fDebugStream, 6: waste your screen
   TLinearFitter *fFitterLinY1;   //!
   TLinearFitter *fFitterLinZ1;   //! 
   TLinearFitter *fFitterLinY2;   //! 
   TLinearFitter *fFitterLinZ2;   //!
   TLinearFitter *fFitterParY;    //! 
   TLinearFitter *fFitterParZ;    //!
   static Int_t fgLoopCounter;   //! only for MakeResPlotsQTreeThread, display status   
   static TMutex *fgWriteMutex;  //!
   static TMutex *fgFitResMutex; //!
   static TMutex *fgFitRmsMutex; //!
   struct TthreadParameterStruct {
      Int_t logLevel;
      Int_t minEntries;
      Int_t dim;
      Int_t pad;
      TH3F *(*resArray)[2][3][11];
      TH3F *(*rmsArray)[2][3][11];
      TString *fileName;
      TTreeSRedirector *fTreeResol;
   };   

protected:         
   ClassDef(AliTPCcalibTracks,1)
};

#endif
