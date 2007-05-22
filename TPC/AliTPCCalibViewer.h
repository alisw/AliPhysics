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


class AliTPCCalibViewer : public TObject {
public:
   AliTPCCalibViewer();
   AliTPCCalibViewer(const AliTPCCalibViewer &c);
   AliTPCCalibViewer(TTree* tree);
   AliTPCCalibViewer(char* fileName, char* treeName = "calPads");
   AliTPCCalibViewer &operator = (const AliTPCCalibViewer & param);
   virtual ~AliTPCCalibViewer();
   
   virtual void     Draw(Option_t* opt="") { fTree->Draw(opt); }
   virtual Long64_t Draw(const char* varexp, const TCut& selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0) { return fTree->Draw(varexp, selection, option, nentries, firstentry); };
   virtual Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0) { return fTree->Draw(varexp, selection, option, nentries, firstentry); };

   Int_t EasyDraw(const char* drawCommand, const char* sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw(const char* drawCommand, Int_t sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw1D(const char* drawCommand, const char* sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t EasyDraw1D(const char* drawCommand, Int_t sector, const char* cuts = 0, const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   // easy drawing of data, use '~' for abbreviation of '.fElements'
   Int_t DrawHisto1D(const char* type, Int_t sector, TVectorF& nsigma, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE) const; // draws 1d histograms and superimposes mean, median, ltm and several sigma cuts
   void SigmaCut(const char* type, Int_t sector, Float_t sigmaMax = 5, Float_t sigmaStep = 0.5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, Bool_t plotLTM = kTRUE) const;    // draws fraction of used pads over different sigma cuts
   
   AliTPCCalPad* GetCalPad(const char* desiredData, char* cuts = "", char* calPadName = "NoName") const;     // returns an AliTPCCalPad object containing the specified data with cuts applied
   AliTPCCalROC* GetCalROC(const char* desiredData, UInt_t sector, char* cuts = "") const;  // returns an AliTPCCalROC object containing the specified data for sector with cuts applied
   
   TObjArray* GetArrayOfCalPads();
   TObjArray* GetListOfVariables(Bool_t printList = kFALSE);
   TObjArray* GetListOfNormalizationVariables(Bool_t printList = kFALSE);
   
   static void MakeTreeWithObjects(const char * fileName, TObjArray * array, const char * mapFileName = 0);
   static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
   
   TFriendElement* AddReferenceTree(const char* filename, const char* treename = "calPads", const char* refname = "R");
   TFriendElement* AddFriend(const char* treename, const char* filename) {return fTree->AddFriend(treename, filename);};
   TFriendElement* AddFriend(TTree* tree, const char* alias, Bool_t warn=kFALSE) {return fTree->AddFriend(tree, alias, warn);};
   TFriendElement* AddFriend(const char* treename, TFile* file) {return fTree->AddFriend(treename, file);};
   
protected:
   TTree* fTree;     // tree containing visualization data (e.g. written by AliTPCCalPad::MakeTree(...)
   TFile* fFile;     // file that contains a calPads tree (e.g. written by AliTPCCalPad::MakeTree(...)
   TObjArray* fListOfObjectsToBeDeleted;  //Objects, that will be deleted when the destructor ist called
   
   ClassDef(AliTPCCalibViewer,1)    //  TPC calibration viewer class
};

#endif
