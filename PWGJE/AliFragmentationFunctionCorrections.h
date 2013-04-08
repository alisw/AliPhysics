// *****************************************************************************
// * Task for corrections to output from AliAnalysisTaskFragmentationFunctions *
//  ****************************************************************************

#ifndef ALIFRAGMENTATIONFUNCTIONCORRECTIONS_H
#define ALIFRAGMENTATIONFUNCTIONCORRECTIONS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class ThnSparse;

class AliFragmentationFunctionCorrections : public TObject {

 public:
  
 //----------------------------------------
  class AliFragFuncCorrHistos : public TObject
  {
    
    public:
    
    AliFragFuncCorrHistos();
    AliFragFuncCorrHistos(const AliFragFuncCorrHistos& copy);
    AliFragFuncCorrHistos& operator=(const AliFragFuncCorrHistos &o);
    virtual ~AliFragFuncCorrHistos();
    AliFragFuncCorrHistos(const char* label,Int_t arraySize);
    void AddCorrHistos(Int_t slice,TH1F* histPt=0,TH1F* histZ=0,TH1F* histXi=0);
    void ReplaceCorrHistos(Int_t slice,TH1F* histPt=0,TH1F* histZ=0,TH1F* histXi=0);

    TH1F* GetTrackPt(const Int_t slice);
    TH1F* GetZ(const Int_t slice);
    TH1F* GetXi(const Int_t slice);

    TString GetLabel() { return fCorrLabel; }

    private:

    Int_t fArraySize;

    TH1F** fh1CorrFFTrackPt;  //! corrected FF histos 
    TH1F** fh1CorrFFZ;        //! corrected FF histos 
    TH1F** fh1CorrFFXi;       //! corrected FF histos 

    TString fCorrLabel;    //! correction label 

    ClassDef(AliFragFuncCorrHistos, 1);
  };

  AliFragmentationFunctionCorrections(); 
  AliFragmentationFunctionCorrections(const  AliFragmentationFunctionCorrections &copy);
  AliFragmentationFunctionCorrections& operator=(const  AliFragmentationFunctionCorrections &o);
  virtual ~AliFragmentationFunctionCorrections();
  
  virtual void SetDebugLevel(Int_t debug){ fDebug = debug; }
  
  void DeleteHistoArray(TH1F** hist) const;
  void DeleteTHnSparseArray(THnSparse** hist) const;
  TH1F** BookHistoArray();
  THnSparse** BookTHnSparseArray();
  void AddCorrectionLevel(const char* label = "");
  void AddCorrectionLevelBgr(const char* label = "");
  void AddCorrectionLevelSinglePt(const char* label = "");
    
  void SetJetPtSlices(Float_t* bins, const Int_t nJetPtSlices);

  void SetHistoBins(const Int_t jetPtSlice, const Int_t sizeBins, Double_t* bins,Int_t type);
  void SetHistoBins(const Int_t jetPtSlice, const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth,Int_t type);
  TArrayD* GetHistoBins(const Int_t jetPtSlice,  const Int_t type);

  void SetHistoBinsSinglePt(const Int_t sizeBins, Double_t* bins);
  void SetHistoBinsSinglePt(const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth);

  // set histo bins for inclusive pt spectra 

  void NormalizeTH1(TH1* hist, const Float_t nJets);
  void NormalizeFF();
  void NormalizeBgr();
  void ReadRawFF(TString strfile, TString strID, TString strFFID = "RecCuts");
  void ReadRawFF(TString strfile, TString strdir, TString strlist, TString strFFID);
  void ReadRawBgr(TString strfile, TString strID, TString strBgrID = "Perp", TString strFFID = "RecCuts");
  void ReadRawBgr(TString strfile, TString strdir, TString strlist, TString strBgrID, TString strFFID);
  void ReadRawBgrEmbedding(TString strfile, TString strID, TString strFFID);
  void ReadRawBgrEmbedding(TString strfile, TString strdir, TString strlist, TString strFFID);

  void WriteOutput(TString strfile, TString strdir = "", Bool_t updateOutfile = kTRUE);

  THnSparse* TH1toSparse(const TH1F* hist, TString strName, TString strTit, const Bool_t fillConst = kFALSE);

  TH1F* Unfold(THnSparse* hnHist, const THnSparse* hnResponse, const THnSparse* hnEff, const Int_t nIter, 
	       const Bool_t useCorrelatedErrors = kTRUE, const THnSparse* hnPrior = 0x0);

  void UnfoldHistos(const Int_t nIter, const Bool_t useCorrelatedErrors, const Int_t type);

  void UnfoldPt(const Int_t nIter=5, const Bool_t useCorrelatedErrors=kTRUE);
  void UnfoldZ(const Int_t nIter=5, const Bool_t useCorrelatedErrors=kTRUE);
  void UnfoldXi(const Int_t nIter=5, const Bool_t useCorrelatedErrors=kTRUE);

  TH1F* ApplyResponse(const TH1F* hist, THnSparse* hnResponse);
  
  void ReadEfficiency(TString strfile, TString strdir = "", TString strlist = "");
  void ReadBgrEfficiency(TString strfile, TString strdir = "", TString strlist = "");

  void EffCorr(); 
  void EffCorrBgr();

  void XiShift(const Int_t corrLevel); 

  void SubtractBgr(Double_t sysErr = 0);

  void WriteSingleTrackEff(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "", TString strPostfix = "");
  void WriteSingleTrackEff(TString strInfile, TString strdir, TString strlist, TString strOutfile, Bool_t updateOutfile = kTRUE, TString strOutDir = "", 
			   TString strPostfix = "");
 
  void WriteSingleTrackSecCorr(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "");
  void WriteSingleTrackSecCorr(TString strInfile, TString strdir, TString strlist, TString strOutfile, Bool_t updateOutfile = kTRUE, TString strOutDir = "");
  
  void WriteSingleResponse(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "");
  void WriteSingleResponse(TString strInfile, TString strdir, TString strlist, TString strOutfile, Bool_t updateOutfile = kTRUE, TString strOutDir = "");
 
  void WriteJetTrackEff(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE);
  void WriteJetTrackEff(TString strInfile, TString strdir, TString strlist, TString strOutfile, Bool_t updateOutfile = kTRUE);

  void WriteJetSecCorr(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "");
  void WriteBgrJetSecCorr(TString strInfile, TString strBgrID, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "");

  void WriteJetSecCorr(TString strInfile, TString strdir, TString strlist, TString strOutfile, Bool_t updateOutfile = kTRUE, 
		       TString strOutDir = "",Bool_t writeBgr=kFALSE,TString strBgrID="");
 
  void WriteJetResponse(TString strInfile, TString strID, TString strOutfile,Bool_t updateOutfile = kTRUE, TString strOutDir = "");
  void WriteJetResponse(TString strInfile, TString strdir, TString strlist,TString strOutfile, Bool_t updateOutfile, TString strOutDir = "");

  void ReadResponse(TString strfile, TString strdir="", TString strlist="");
  void ReadPriors(TString strfile,const Int_t type); 

  void ProjectSingleResponseMatrix(TString strOutfile, Bool_t updateOutfile, TString strOutDir = "");
  void ProjectJetResponseMatrices(TString strOutfile);

  void RebinHisto(const Int_t jetPtSlice, const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth, const Int_t type);

  void WriteJetSpecResponse(TString strInfile, TString strdir, TString strlist /*, TString strOutfile*/ );

  void ReadSingleTrackEfficiency(TString strfile, TString strdir="", TString strlist="", TString strname="hSingleTrackEffPt");
  void ReadSingleTrackResponse(TString strfile, TString strdir="", TString strlist="", TString strname="fhnResponseSinglePt");
  void ReadSingleTrackSecCorr(TString strfile, TString strdir="", TString strlist="", TString strname="hSingleTrackSecCorrPt");
  void ReadSingleTrackCorrection(TString strfile, TString strdir, TString strlist, TString strname, const Int_t type);
  
  void ReadRawPtSpec(TString strInfile, TString strID);
  void ReadRawPtSpec(TString strfile, TString strdir, TString strlist);
  void ReadRawPtSpecQATask(TString strfile, TString strdir, TString strlist); // spectra from Martas QA task
  void EffCorrSinglePt();
  void UnfoldSinglePt(const Int_t nIter, const Bool_t useCorrelatedErrors);
  void SecCorrSinglePt();
  void dNdz2dNdxi();

  void WriteBinShiftCorr(TString strInfile, TString strIDGen,  TString strIDRec,  
			 TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim = kTRUE,  
			 TString strOutDir = "");

  void WriteBgrBinShiftCorr(TString strInfile, TString strBgrID, TString strIDGen,  TString strIDRec,  
			    TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim = kTRUE,  
			    TString strOutDir = "");

  void WriteBinShiftCorr(TString strInfile, TString strdirGen, TString strlistGen, 
			 TString strdirRec, TString strlistRec, 
			 TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim = kTRUE, 
			 TString strOutDir = "",Bool_t writeBgr = kFALSE, TString strBgrID = "");

  void ReadBgrBinShiftCorr(TString strfile,  TString strBgrID, TString strdir="", TString strlist="");
  void ReadBinShiftCorr(TString strfile, TString strdir="", TString strlist="", Bool_t readBgr = kFALSE, TString strBgrID="");

  void ReadFoldingCorr(TString strfile, TString strdir="", TString strlist="");

  void BbBCorr();
  void BbBCorrBgr();

  void FoldingCorr();

  void ReadBgrJetSecCorr(TString strfile, TString strBgrID, TString strdir="", TString strlist="", Bool_t useScaledStrangeness=kTRUE);
  void ReadJetSecCorr(TString strfile, TString strdir="", TString strlist="", Bool_t useScaledStrangeness = kTRUE, 
		      Bool_t readBgr=kFALSE, TString strBgrID="");

  void JetSecCorr();
  void JetSecCorrBgr();


  void WriteBinShiftCorrSinglePt(TString strInfile, TString strIDGen,  TString strIDRec,  
				 TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, TString strOutDir);

  void WriteBinShiftCorrSinglePt(TString strInfile, TString strdirGen, TString strlistGen, 
				 TString strdirRec, TString strlistRec, 
				 TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, TString strOutDir = "");
  
  void ReadBinShiftCorrSinglePt(TString strfile, TString strdir = "", TString strlist = "", Bool_t useRecPrim = kTRUE);

  void BbBCorrSinglePt();


  enum {kFlagPt=0,kFlagZ,kFlagXi,kFlagSinglePt};
  enum {kFlagEfficiency=0,kFlagResponse,kFlagSecondaries};


 private:

  static const Int_t fgMaxNCorrectionLevels = 10;  //! max number of corrections 
  
  Int_t fDebug;              //! Debug level
  Int_t fNJetPtSlices;       //! n slices in jet pt
  TArrayF* fJetPtSlices;     //! array to hold slices in jet pt 

  TArrayF* fNJets;           //! nJets per jet pt slice - non-int e.g. for xsec/nTrials scaled spectra
  TArrayF* fNJetsBgr;        //! nJets bgr per jet pt slice - non-int  e.g. for xsec/nTrials scaled spectra
 
  Int_t fNHistoBinsSinglePt;  //! nBins inclusive pt spec histos - leave undefinded to use original binning
  TArrayD* fHistoBinsSinglePt; //! inclusive pt spec histo bins

  Int_t* fNHistoBinsPt;      //! nBins FF histos in each jet pt slice - leave undefinded for any slice to use original binning
  Int_t* fNHistoBinsZ;       //! nBins FF histos in each jet pt slice - leave undefinded for any slice to use original binning
  Int_t* fNHistoBinsXi;      //! nBins FF histos in each jet pt slice - leave undefinded for any slice to use original binning

  TArrayD** fHistoBinsPt;    //! FF histo bins 
  TArrayD** fHistoBinsZ;     //! FF histo bins
  TArrayD** fHistoBinsXi;    //! FF histo bins

  Int_t fNCorrectionLevels;        //! corrections applied: efficiency, secondaries, resolution unfolding, bgr subtraction
  AliFragFuncCorrHistos** fCorrFF; //! array of fragmentation functions, dimensions: jet pt bins, correction steps

  Int_t fNCorrectionLevelsBgr;      //! corrections applied: efficiency, secondaries, resolution unfolding, bgr subtraction
  AliFragFuncCorrHistos** fCorrBgr; //! array of bgr fragmentation functions, dimensions: jet pt bins, correction steps

  Int_t fNCorrectionLevelsSinglePt;      //! corrections applied: efficiency, secondaries, resolution unfolding, bgr subtraction
  AliFragFuncCorrHistos** fCorrSinglePt; //! array to keep single track pt spectra, 1D in jet pt bins dimension 



  // xi shift
  TH1F** fh1FFXiShift;          //! FF: track xi, corrected for shift in jet energy

  // eff correction
  TH1F*  fh1EffSinglePt;       //!  efficiency all tracks

  TH1F** fh1EffPt;             //! reconstruction efficiency track pt
  TH1F** fh1EffZ;              //! reconstruction efficiency z
  TH1F** fh1EffXi;             //! reconstruction efficiency xi

  TH1F** fh1EffBgrPt;          //! reconstruction efficiency bgr track pt
  TH1F** fh1EffBgrZ;  	       //! reconstruction efficiency bgr z
  TH1F** fh1EffBgrXi;	       //! reconstruction efficiency bgr xi

  // bin-by-bin correction
  TH1F*  fh1BbBCorrSinglePt;   //!  BbB corr track pt 

  TH1F** fh1BbBPt;             //! bin-by-bin correction track pt
  TH1F** fh1BbBZ;              //! bin-by-bin correction z
  TH1F** fh1BbBXi;             //! bin-by-bin correction xi

  TH1F** fh1BbBBgrPt;          //! bin-by-bin correction UE track pt
  TH1F** fh1BbBBgrZ;           //! bin-by-bin correction UE z
  TH1F** fh1BbBBgrXi;          //! bin-by-bin correction UE xi

  TH1F** fh1FoldingCorrPt;   //! corr factor rec/folded 
  TH1F** fh1FoldingCorrZ;    //! corr factor rec/folded 
  TH1F** fh1FoldingCorrXi;   //! corr factor rec/folded 


  // secondaries correction
  TH1F** fh1SecCorrPt;             //! secondaries correction track pt
  TH1F** fh1SecCorrZ;              //! secondaries correction z
  TH1F** fh1SecCorrXi;             //! secondaries correction xi

  TH1F** fh1SecCorrBgrPt;          //! secondaries correction track pt
  TH1F** fh1SecCorrBgrZ;           //! secondaries correction z
  TH1F** fh1SecCorrBgrXi;          //! reconstruction efficiency xi

  // unfolding

  TH1F** fh1FFTrackPtBackFolded;  //! FF: track pt backfolded (unfolded + smeared with response matrix)
  TH1F** fh1FFZBackFolded;        //! FF: track z, backfolded (unfolded + smeared with response matrix)
  TH1F** fh1FFXiBackFolded;       //! FF: track xi,backfolded (unfolded + smeared with response matrix)

  TH1F** fh1FFRatioTrackPtFolded;  //! ratio FF: track pt unfolded / original 
  TH1F** fh1FFRatioZFolded;        //! ratio FF: track z  unfolded / original
  TH1F** fh1FFRatioXiFolded;       //! ratio FF: track xi unfolded / original

  TH1F** fh1FFRatioTrackPtBackFolded;  //! ratio FF: track pt backfolded / original
  TH1F** fh1FFRatioZBackFolded;        //! ratio FF: track z  backfolded / original
  TH1F** fh1FFRatioXiBackFolded;       //! ratio FF: track xi backfolded / original

  TH1F*  fh1SingleTrackPtBackFolded;      //! inclusive track pt backfolded (unfolded + smeared with response matrix)
  TH1F*  fh1RatioSingleTrackPtFolded;     //! ratio inclusive track pt unfolded / original 
  TH1F*  fh1RatioSingleTrackPtBackFolded; //! ratio inblusive track pt backfolded / original

  THnSparse*  fhnResponseSinglePt;  //!  response matrix pt gen vs rec all tracks
  THnSparse** fhnResponsePt;        //!  response matrix pt gen vs rec 
  THnSparse** fhnResponseZ;         //!  response matrix z  gen vs rec 
  THnSparse** fhnResponseXi;        //!  response matrix xi gen vs rec 

  TH1F** fh1FFTrackPtPrior;  //! FF: track pt prior 
  TH1F** fh1FFZPrior;        //! FF: track z  prior
  TH1F** fh1FFXiPrior;       //! FF: track xi prior


  // secondaries 
  TH1F*  fh1SecCorrSinglePt;       //!  secondaries correction all tracks


  ClassDef(AliFragmentationFunctionCorrections, 1);
};

#endif
