#ifndef ALIBASECALIBVIEWER_H
#define ALIBASECALIBVIEWER_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class for the AliTPCCalibViewer and AliTRDCalibViewer               //
//  used for the calibration monitor                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include "TFriendElement.h"

#include "AliMathBase.h"

class TLegend;

class AliBaseCalibViewer : public TObject {
 public:
   AliBaseCalibViewer();
   AliBaseCalibViewer(const AliBaseCalibViewer &c);
   AliBaseCalibViewer(TTree* tree);
   AliBaseCalibViewer(const Char_t* fileName, const Char_t* treeName = "tree");
   AliBaseCalibViewer &operator = (const AliBaseCalibViewer & param);
   virtual ~AliBaseCalibViewer();
   virtual void Delete(Option_t* option = "");

   TString& GetAbbreviation()  { return fAbbreviation;  }
   TString& GetAppendString()  { return fAppendString; }
   void SetAbbreviation(const Char_t* abr) { fAbbreviation = abr; }
   void SetAppendString(const Char_t* str) { fAppendString = str; }

   //virtual void GetTimeInfoOCDB(const Char_t* runList, const Char_t* outFile,
                                //Int_t firstRun, Int_t lastRun, UInt_t infoFlags,
                                //const Char_t* ocdbStorage) = 0;

   virtual void     Draw(Option_t* opt="") { fTree->Draw(opt); }
   virtual Long64_t Draw(const Char_t* varexp, const TCut& selection, Option_t* option = "", 
			 Long64_t nentries = 1000000000, Long64_t firstentry = 0) { 
     return fTree->Draw(varexp, selection, option, nentries, firstentry); 
   };
   virtual Long64_t Draw(const Char_t* varexp, const Char_t* selection, Option_t* option = "", 
			 Long64_t nentries = 1000000000, Long64_t firstentry = 0) { 
     return fTree->Draw(varexp, selection, option, nentries, firstentry); 
   };

   virtual const char* AddAbbreviations(Char_t* c, Bool_t printDrawCommand = kFALSE) = 0;
   // easy drawing of data, use '~' for abbreviation of '.fElements'
   virtual Int_t EasyDraw(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts = 0, 
			  const Char_t* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const = 0;   
   // easy drawing of data, use '~' for abbreviation of '.fElements'
   virtual Int_t EasyDraw(const Char_t* drawCommand, Int_t sector, const Char_t* cuts = 0, 
			  const Char_t* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const = 0;   
   // easy drawing of data, use '~' for abbreviation of '.fElements'
   virtual Int_t EasyDraw1D(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts = 0, 
			    const Char_t* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const = 0;   
   // easy drawing of data, use '~' for abbreviation of '.fElements'
   virtual Int_t EasyDraw1D(const Char_t* drawCommand, Int_t sector, const Char_t* cuts = 0, 
			    const Char_t* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const = 0;   
   // formats title and axis labels of histo, removes '.fElements'
   void FormatHistoLabels(TH1 *histo) const;   
   // draws 1d histograms and superimposes mean, median, ltm and several sigma cuts
   Int_t  DrawHisto1D(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts = 0, 
		      const Char_t *sigmas = "2;4;6", Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, 
		      Bool_t plotLTM = kTRUE) const; 
   // draws fraction of used pads over different sigma cuts
   Int_t     SigmaCut(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts = 0, 
		      Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, 
		      Bool_t plotLTM = kTRUE, Bool_t pm = kFALSE, const Char_t *sigmas = "", 
		      Float_t sigmaStep = -1) const;    
   // draws an integrated histogram
   Int_t    Integrate(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts = 0, 
		      Float_t sigmaMax = 5, Bool_t plotMean = kTRUE, Bool_t plotMedian = kTRUE, 
		      Bool_t plotLTM = kTRUE, const Char_t *sigmas = "", Float_t sigmaStep = -1) const;    

   virtual TObjArray* GetListOfVariables(Bool_t printList = kFALSE) = 0;
   virtual TObjArray* GetListOfNormalizationVariables(Bool_t printList = kFALSE) const = 0;

   TFriendElement* AddReferenceTree(const Char_t* filename, const Char_t* treename = "tree", const Char_t* refname = "R");
   TFriendElement* AddFriend(const Char_t* treename, const Char_t* filename) 
   {return fTree->AddFriend(treename, filename);};
   TFriendElement* AddFriend(TTree* tree, const Char_t* alias, Bool_t warn=kFALSE) 
   {return fTree->AddFriend(tree, alias, warn);};
   TFriendElement* AddFriend(const Char_t* treename, TFile* file) 
   {return fTree->AddFriend(treename, file);};
   TTree * GetTree() const { return fTree;}

   TString* Fit(const Char_t* drawCommand, const Char_t* formula, const Char_t* cuts, 
		Double_t & chi2, TVectorD &fitParam, TMatrixD &covMatrix);
   static Double_t GetLTM(Int_t n, Double_t *array, Double_t *sigma = 0, Double_t fraction = 0.9);
   static Int_t GetBin(Float_t value, Int_t nbins, Double_t binLow, Double_t binUp);
   static TH1F* SigmaCut(Int_t n, Float_t *array, Float_t mean, Float_t sigma, Int_t nbins, 
			 Float_t binLow, Float_t binUp, Float_t sigmaMax, Float_t sigmaStep = -1, Bool_t pm = kFALSE);
   static TH1F* SigmaCut(TH1F *histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, 
			 Float_t sigmaStep = -1, Bool_t pm = kFALSE);
   static TH1F* Integrate(TH1F *histogram, Float_t mean = 0, Float_t sigma = 0, 
			  Float_t sigmaMax = 0, Float_t sigmaStep = -1);
   static TH1F* Integrate(Int_t n, Float_t *array, Int_t nbins, Float_t binLow, Float_t binUp, 
			  Float_t mean = 0, Float_t sigma = 0, Float_t sigmaMax = 0, Float_t sigmaStep = -1);
   static TH1F* SigmaCut(Int_t n, Double_t *array, Double_t mean, Double_t sigma, 
			 Int_t nbins, Double_t *xbins, Double_t sigmaMax);

 protected:
   TTree* fTree;     // tree containing visualization data (e.g. written by AliTPCCalPad::MakeTree(...)
   TFile* fFile;     // file that contains a calPads tree (e.g. written by AliTPCCalPad::MakeTree(...)
   TObjArray* fListOfObjectsToBeDeleted;  //Objects, that will be deleted when the destructor ist called
   Bool_t fTreeMustBeDeleted;  // decides weather the tree must be deleted in destructor or not 
   TString fAbbreviation;       // the abreviation for '.fElements'
   TString fAppendString;      // '.fElements', stored in a TStrig
   
   void DrawLines(TH1F *cutHistoMean, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const;
   void DrawLines(TGraph *graph, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const;
   
   
   ClassDef(AliBaseCalibViewer,1)    //  Base calibration viewer class
};

#endif
