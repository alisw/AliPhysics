#ifndef AliTPCCALIBTRACKS_H
#define AliTPCCALIBTRACKS_H


#include <TNamed.h>
#include <TH2D.h>
#include <TF2.h>
#include <TSystem.h>
#include <TCollection.h>
#include <iostream>
using namespace std;

// #include "AliTPCClusterParam.h"


class AliTPCClusterParam;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
class AliESDtrack;
class TH3F;
class TH1F;
class TH1I;
class TTreeSRedirector;
class AliTPCcalibTracksCuts;
class TChain;


class AliTPCcalibTracks : public TNamed {
public :
   AliTPCcalibTracks();
   AliTPCcalibTracks(AliTPCcalibTracks* ct);
   AliTPCcalibTracks(const Text_t *name, const Text_t *title, AliTPCClusterParam *clusterParam, AliTPCcalibTracksCuts* cuts);
   virtual ~AliTPCcalibTracks();
   
   static void     AddInfo(TChain * chain, char* fileName);
   static void     AddCuts(TChain * chain, char* ctype);
   void            Process(AliTPCseed *track, AliESDtrack *esd);
   
   Float_t         TPCBetheBloch(Float_t bg);
   Bool_t          AcceptTrack(AliTPCseed * track); 
   void            FillHistoCluster(AliTPCseed * track);
   void            FillResolutionHistoLocal(AliTPCseed * track);
   void            AlignUpDown(AliTPCseed * track, AliESDtrack *esd);
   static  Int_t   GetBin(Float_t q, Int_t pad);
   static  Int_t   GetBin(Int_t  iq, Int_t pad);
   static  Float_t GetQ(Int_t bin);
   static  Float_t GetPad(Int_t bin){return bin%3;}

   Long64_t Merge(TCollection *li);
   static AliTPCcalibTracks* TestMerge(AliTPCcalibTracks *ct, AliTPCClusterParam *clusterParam, Int_t nCalTracks = 50);
   void     MakeReport(Int_t stat, char* pathName = "plots");
   void     MakeAmpPlots(Int_t stat, char* pathName = "plots");
   void     MakeDeltaPlots(char* pathName = "plots");
   void     FitResolutionNew(char* pathName = "plots");
   void     FitRMSNew(char* pathName = "plots");
   void     MakeResPlotsQ(Int_t minEntries = 1,  Bool_t bDraw=kFALSE, char* pathName = "plots"); 
   void     MakeResPlotsQTree(Int_t minEntries = 1, char* pathName = "plots"); 
   void     Draw(Option_t* opt);
   void     SetStyle();

   static TH2D      *MakeDiff(TH2D * hfit, TF2 * func);
   static TObjArray *FitProjections(TH3F * hfit, Int_t val = 2, Int_t minEntry=500, Bool_t bDraw=kFALSE);
   
//protected:   
   TObjArray* GetfArrayAmpRow() {return fArrayAmpRow;}
   TObjArray* GetfArrayAmp() {return fArrayAmp;}
   TObjArray* GetfArrayQDY() {return fArrayQDY;}
   TObjArray* GetfArrayQDZ() {return fArrayQDZ;}
   TObjArray* GetfArrayQRMSY() {return fArrayQRMSY;}
   TObjArray* GetfArrayQRMSZ() {return fArrayQRMSZ;}
   TH1F*      GetfDeltaY() {return fDeltaY;}
   TH1F*      GetfDeltaZ() {return fDeltaZ;}
   TObjArray* GetfResolY() {return fResolY;}
   TObjArray* GetfResolZ() {return fResolZ;}
   TObjArray* GetfRMSY() {return fRMSY;}
   TObjArray* GetfRMSZ() {return fRMSZ;}
   TH1I*      GetfHclus() {return fHclus;}
   AliTPCcalibTracksCuts* GetCuts() {return fCuts;}

   
   
private:
   AliTPCClusterParam *fClusterParam; //! pointer to cluster parameterization
   TTreeSRedirector   *fDebugStream;  //! debug stream for
   AliTPCROC *fROC;          //!
   TObjArray *fArrayAmpRow; // array with amplitudes versus row for given sector 
   TObjArray *fArrayAmp;    // array with amplitude for sectors
   TObjArray *fArrayQDY;    // q binned delta Y histograms
   TObjArray *fArrayQDZ;    // q binned delta Z histograms 
   TObjArray *fArrayQRMSY;  // q binned delta Y histograms
   TObjArray *fArrayQRMSZ;  // q binned delta Z histograms 
   TH1F      *fDeltaY;      // integrated delta y histo
   TH1F      *fDeltaZ;      // integrated delta z histo
   TObjArray *fResolY;      // array of resolution histograms Y
   TObjArray *fResolZ;      // array of resolution histograms Z
   TObjArray *fRMSY;        // array of RMS histograms Y
   TObjArray *fRMSZ;        // array of RMS histograms Z
   AliTPCcalibTracksCuts *fCuts; // object with cuts, that is passed to the constructor
   TH1I *fHclus;             // number of clusters per track

protected:         
   ClassDef(AliTPCcalibTracks,1)
};

#endif




#ifndef AliTPCCALIBTRACKSCUTS_H
#define AliTPCCALIBTRACKSCUTS_H

class AliTPCcalibTracksCuts: public TNamed {
   //////////////////////////////////////////////////////
   //                                                  //
   //     Class to specify cuts for track analysis     //
   //     with AliTPCcalibTracks                       //
   //                                                  //
   //////////////////////////////////////////////////////

public:
   AliTPCcalibTracksCuts(Int_t minClusters = 20, Float_t minRatio = 0.4, Float_t max1pt = 0.5,
      Float_t edgeXZCutNoise = 0.13, Float_t edgeThetaCutNoise = 0.018):
         TNamed("calibTracksCuts", "calibTracksCuts") {
      // 
      // Constuctor for AliTPCcalibTracksCuts
      // specify the cuts to be set on the processed tracks
      // 
      fMinClusters = minClusters;
      fMinRatio = minRatio;
      fMax1pt = max1pt;
      fEdgeYXCutNoise = edgeXZCutNoise;
      fEdgeThetaCutNoise = edgeThetaCutNoise;
   }
   virtual ~AliTPCcalibTracksCuts(){cout << "AliTPCcalibTracksCuts destructor called, nothing happend." << endl;}
   void SetMinClusters(Int_t minClusters){fMinClusters = minClusters;}
   void SetMinRatio(Float_t minRatio){fMinRatio = minRatio;}
   void SetMax1pt(Float_t max1pt){fMax1pt = max1pt;}
   void SetEdgeXYCutNoise(Float_t edgeCutNoise){fEdgeYXCutNoise = edgeCutNoise;}
   void SetEdgeThetaCutNoise(Float_t edgeCutNoise){fEdgeThetaCutNoise = edgeCutNoise;}
   const Int_t   GetMinClusters(){return fMinClusters;}
   const Float_t GetMinRatio(){return fMinRatio;}
   const Float_t GetMax1pt(){return fMax1pt;}
   const Float_t GetEdgeYXCutNoise(){return fEdgeYXCutNoise;}
   const Float_t GetEdgeThetaCutNoise(){return fEdgeThetaCutNoise;}
   virtual void Print(Option_t* option = ""){
      option = option;  // to avoid compiler warnings
      cout << "<AliTPCcalibTracksCuts>: The following cuts are specified: " << endl;
      cout << "fMinClusters: " << fMinClusters << endl;
      cout << "fMinRatio: " << fMinRatio << endl;
      cout << "fMax1pt: " << fMax1pt << endl;
      cout << "fEdgeYXCutNoise: " << fEdgeYXCutNoise << endl;
      cout << "fEdgeThetaCutNoise: " << fEdgeThetaCutNoise << endl;
   }  // Prints out the specified cuts
   
private:
   Int_t   fMinClusters;         // number of clusters
   Float_t fMinRatio;            // kMinRratio = 0.4
   Float_t fMax1pt;              // kMax1pt = 0.5
   Float_t fEdgeYXCutNoise;      // kEdgeYXCutNoise = 0.13
   Float_t fEdgeThetaCutNoise;   // kEdgeThetaCutNoise = 0.018

protected:         
   ClassDef(AliTPCcalibTracksCuts,1)
};


#endif


