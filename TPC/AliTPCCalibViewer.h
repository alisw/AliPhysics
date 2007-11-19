#ifndef ALITPCCALIBVIEWER_H
#define ALITPCCALIBVIEWER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalibViewer.h,v */

///////////////////////////////////////////////////
//                                               //
//  TPC calibration viewer/visualization class   //
//  use Tree for visualization                   //
///////////////////////////////////////////////////

#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TFriendElement.h"

#include "AliMathBase.h"
class TLegend;


class AliTPCCalibViewer : public TObject {
public:
   AliTPCCalibViewer();
   AliTPCCalibViewer(const AliTPCCalibViewer &c);
   AliTPCCalibViewer(TTree* tree);
   AliTPCCalibViewer(char* fileName, char* treeName = "calPads");
   AliTPCCalibViewer &operator = (const AliTPCCalibViewer & param);
   virtual ~AliTPCCalibViewer();
   virtual void Delete(Option_t* option = "");
   
   virtual void     Draw(Option_t* opt="") { fTree->Draw(opt); }
   virtual Long64_t Draw(const char* varexp, const TCut& selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0) { return fTree->Draw(varexp, selection, option, nentries, firstentry); };
   virtual Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0) { return fTree->Draw(varexp, selection, option, nentries, firstentry); };

   Int_t EasyDraw(const char* drawCommand, const char* sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw(const char* drawCommand, Int_t sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw1D(const char* drawCommand, const char* sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw1D(const char* drawCommand, Int_t sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   
   Int_t  DrawHisto1D(const char* drawCommand,       Int_t sector, const char* cuts = 0, const char *sigmas = "2;4;6", Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE) const; // draws 1d histograms and superimposes mean, median, ltm and several sigma cuts
   Int_t  DrawHisto1D(const char* drawCommand, const char* sector, const char* cuts = 0, const char *sigmas = "2;4;6", Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE) const; // draws 1d histograms and superimposes mean, median, ltm and several sigma cuts
   Int_t     SigmaCut(const char* drawCommand,       Int_t sector, const char* cuts = 0, Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE, Bool_t pm = kFALSE, const char *sigmas = "", Float_t sigmaStep = -1) const;    // draws fraction of used pads over different sigma cuts
   Int_t     SigmaCut(const char* drawCommand, const char* sector, const char* cuts = 0, Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE, Bool_t pm = kFALSE, const char *sigmas = "", Float_t sigmaStep = -1) const;    // draws fraction of used pads over different sigma cuts
   Int_t    Integrate(const char* drawCommand, const char* sector, const char* cuts = 0, Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE, const char *sigmas = "", Float_t sigmaStep = -1) const;    // draws an integrated histogram
   Int_t    Integrate(const char* drawCommand,       Int_t sector, const char* cuts = 0, Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE, const char *sigmas = "", Float_t sigmaStep = -1) const;    // draws an integrated histogram
   
   AliTPCCalPad* GetCalPad(const char* desiredData, char* cuts = "", char* calPadName = "NoName") const;     // returns an AliTPCCalPad object containing the specified data with cuts applied
   AliTPCCalROC* GetCalROC(const char* desiredData, UInt_t sector, char* cuts = "") const;  // returns an AliTPCCalROC object containing the specified data for sector with cuts applied
   
   TObjArray* GetArrayOfCalPads();
   TObjArray* GetListOfVariables(Bool_t printList = kFALSE);
   TObjArray* GetListOfNormalizationVariables(Bool_t printList = kFALSE) const;
   
   static void MakeTreeWithObjects(const char * fileName, TObjArray * array, const char * mapFileName = 0);
   static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
   
   TFriendElement* AddReferenceTree(const char* filename, const char* treename = "calPads", const char* refname = "R");
   TFriendElement* AddFriend(const char* treename, const char* filename) {return fTree->AddFriend(treename, filename);};
   TFriendElement* AddFriend(TTree* tree, const char* alias, Bool_t warn=kFALSE) {return fTree->AddFriend(tree, alias, warn);};
   TFriendElement* AddFriend(const char* treename, TFile* file) {return fTree->AddFriend(treename, file);};
   TTree * GetTree() { return fTree;}

   TString* Fit(const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, TVectorD &fitParam, TMatrixD &covMatrix);

  // 
  // Array tools
  // 
  static Double_t GetLTM(Int_t n, Double_t *array, Double_t *sigma = 0, Double_t fraction = 0.9);
  static Int_t GetBin(Float_t value, Int_t nbins, Double_t binLow, Double_t binUp);
  static TH1F* SigmaCut(Int_t n, Float_t *array, Float_t mean, Float_t sigma, Int_t nbins, Float_t binLow, Float_t binUp, Float_t sigmaMax, Float_t sigmaStep = -1, Bool_t pm = kFALSE);
  static TH1F* SigmaCut(TH1F *histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep = -1, Bool_t pm = kFALSE);
  static TH1F* Integrate(TH1F *histogram, Float_t mean = 0, Float_t sigma = 0, Float_t sigmaMax = 0, Float_t sigmaStep = -1);
  static TH1F* Integrate(Int_t n, Float_t *array, Int_t nbins, Float_t binLow, Float_t binUp, Float_t mean = 0, Float_t sigma = 0, Float_t sigmaMax = 0, Float_t sigmaStep = -1);
  
  static TH1F* SigmaCut(Int_t n, Double_t *array, Double_t mean, Double_t sigma, Int_t nbins, Double_t *xbins, Double_t sigmaMax);
   
   
         
protected:
   TTree* fTree;     // tree containing visualization data (e.g. written by AliTPCCalPad::MakeTree(...)
   TFile* fFile;     // file that contains a calPads tree (e.g. written by AliTPCCalPad::MakeTree(...)
   TObjArray* fListOfObjectsToBeDeleted;  //Objects, that will be deleted when the destructor ist called
   Bool_t fTreeMustBeDeleted;  // decides weather the tree must be deleted in destructor or not 
   
   void DrawLines(TH1F *cutHistoMean, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const;
  
   
   ClassDef(AliTPCCalibViewer,1)    //  TPC calibration viewer class
};

#endif
