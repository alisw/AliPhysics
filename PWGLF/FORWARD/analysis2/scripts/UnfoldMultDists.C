/**
 * @file   UnfoldMultDists.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Nov 12 09:25:52 2013
 * 
 * @brief  A class to do unfolding 
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TClass.h>
#include <TRegexp.h>
#include <TMath.h>
#include <TParameter.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include <fstream>

/**
 * Class to do unfolding of raw histograms produced by AliForwardMultDists 
 * 
 * @ingroup pwglf_forward_multdist
 */
struct Unfolder
{
  /**
   * Colours used 
   * 
   */
  enum {
    kColorMeasured = kOrange-2, 
    kColorTruth    = kBlue-3, 
    kColorAccepted = kMagenta-3,
    kColorTrgVtx   = kBlack, 
    kColorUnfolded = kOrange+2,
    kColorCorrected= kRed+2, 
    kColorError    = kBlue-10,
    kColorALICE    = kPink+1, 
    kColorCMS      = kGreen+2
  };

  /** 
   * Constructor 
   */
  Unfolder() {}
  /** 
   * Get a top collection from a file
   * 
   * @param fileName Name of file 
   * @param results  Wheter it's the results collection or not 
   * 
   * @return Collection or null
   */
  static TCollection* GetTop(const TString& fileName, Bool_t results=false)
  {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file) { 
      Error("GetTop", "Failed to open %s", fileName.Data());
      return 0;
    }
    TCollection* ret = 0;
    TString cName(Form("ForwardMult%s", results ? "Results" : "Sums"));
    file->GetObject(cName, ret);
    if (!ret) 
      Error("GetTop", "Failed to get collection %s from %s", 
	    cName.Data(), fileName.Data());
    file->Close();
    return ret;
  }
  /** 
   * Get an object from a collection 
   * 
   * @param c    Collection
   * @param name Name of object
   * @param cl   Possible class to check against
   * @param verbose  Be verbose
   * 
   * @return Pointer to object or null
   */
  static TObject* GetObject(TCollection* c, const TString& name, 
			    TClass* cl, Bool_t verbose=true)
  {
    TObject* o = c->FindObject(name);
    if (!o) { 
      if (verbose)
	Warning("GetObject", "%s not found in %s", name.Data(), c->GetName());
      return 0;
    }
    if (cl && !o->IsA()->InheritsFrom(cl)) {
      if (verbose) 
	Warning("GetCollection", "%s is not a %s but a %s", 
		name.Data(), cl->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  /** 
   * Get a collection 
   * 
   * @param c        Parent collection
   * @param name     Name of object to findf
   * @param verbose  Be verbose
   * 
   * @return 
   */
  static TCollection* GetCollection(TCollection*   c, 
				    const TString& name, 
				    Bool_t         verbose=-true)
  {
    return static_cast<TCollection*>(GetObject(c, name, 
					       TCollection::Class(),
					       verbose));
  }
  /** 
   * Get a 1D histogram from a collection
   * 
   * @param c    Collection
   * @param name Nanme of histogram
   * @param verbose  Be verbose
   * 
   * @return Pointer to object or null
   */
  static TH1* GetH1(TCollection* c, const TString& name, Bool_t verbose=true) 
  {
    return static_cast<TH1*>(GetObject(c, name, TH1::Class(), verbose));
  }
  /** 
   * Get a 2D histogram from a collection
   * 
   * @param c    Collection
   * @param name Nanme of histogram
   * @param verbose  Be verbose
   * 
   * @return Pointer to object or null
   */
  static TH2* GetH2(TCollection* c, const TString& name, Bool_t verbose=true) 
  {
    return static_cast<TH2*>(GetObject(c, name, TH2::Class(), verbose));
  }
  /** 
   * Get an unsigned short parameter from the collection 
   * 
   * @param c    Collection
   * @param name Parameter name 
   * @param v    Value
   * 
   * @return Value 
   */
  static void GetParameter(TCollection* c, const TString& name, UShort_t& v)
  {
    TObject* o = GetObject(c, name, TParameter<int>::Class(), true);
    v = (!o ? 0 : o->GetUniqueID());
  }
  /** 
   * Get an unsigned short parameter from the collection 
   * 
   * @param c    Collection
   * @param name Parameter name 
   * @param v    Value
   * 
   * @return Value 
   */
  static void GetParameter(TCollection* c, const TString& name, ULong_t& v)
  {
    TObject* o = GetObject(c, name, TParameter<long>::Class(), true);
    v = (!o ? 0 : o->GetUniqueID());
  }
  /** 
   * Get an unsigned short parameter from the collection 
   * 
   * @param c    Collection
   * @param name Parameter name 
   * @param v    Value
   * 
   * @return Value 
   */
  static void GetParameter(TCollection* c, const TString& name, Double_t& v)
  {
    TObject* o = GetObject(c, name, TParameter<double>::Class(), true);
    v = (!o ? 0 : static_cast<TParameter<double>*>(o)->GetVal());
  }
  /** 
   * Get the method identifier 
   * 
   * @param method Method 
   * 
   * @return Method identifier 
   */    
  static UInt_t MethodId(TString& method) 
  {
    struct Method { 
      UInt_t  id;
      TString name;
    };
    const Method methods[] = { {RooUnfold::kNone,    "None"},
			       {RooUnfold::kBayes,   "Bayes"},
			       {RooUnfold::kSVD,     "SVD"},
			       {RooUnfold::kBinByBin,"BinByBin"},
			       {RooUnfold::kTUnfold, "TUnfold"},
			       {RooUnfold::kInvert,  "Invert"},
			       {RooUnfold::kDagostini,"Dagostini"}, 
			       {0xDeadBeef,           "unknown"} };
    const Method* pMethod = methods;
    while (pMethod->id != 0xDeadBeef) {
      if (method.BeginsWith(pMethod->name, TString::kIgnoreCase)) {
	method = pMethod->name;
	break;
      }
      pMethod++;
    }
    if (pMethod->id == 0xDeadBeef) 
      Error("MethodId", "Unknown unfolding method: %s", method.Data());

    return pMethod->id;
  }
  
  /** 
   * Run the unfolding and correction task. 
   *
   * The @a measuredFile is assumed to have the structure 
   *
   * @verbatim
   *     /- ForwardMultSums         (TCollection)
   *          |- [type]             (TCollection)
   *          |    |- [bin]         (TCollection) 
   *          |    |    `- rawDist  (TH1)
   *          |    |- [bin]
   *          |    ...
   *          |- [type]
   *          ...
   * @endverbatim 
   * 
   * and @a corrFile is assumed to have the structure 
   *
   * @verbatim
   *     /- ForwardMultResults            (TCollection)
   *          |- [type]                   (TCollection)
   *          |    |- [bin]               (TCollection) 
   *          |    |    |- truth          (TH1)
   *          |    |    |- truthAccepted  (TH1)
   *          |    |    |- triggerVertex  (TH1)
   *          |    |    `- response       (TH2)
   *          |    |- [bin]
   *          |    ...
   *          |- [type]
   *          ...
   * @endverbatim 
   *
   * where @c [type] is one of <i>symmetric</i>, <i>positive</i>,
   * <i>negative</i>, or <i>other</i>, and [bin] is the @f$ \eta@f$
   * bin named like
   *
   * @verbatim
   *   [bin]          := [eta_spec] _ [eta_spec]
   *   [eta_spec]     := [sign_char] [integer_part] d [decimal_part]
   *   [sign_part]    := p          positive eta 
   *                  |  m          negative eta 
   *   [integer_part] := [number]
   *   [decimal_part] := [number] [number]
   *   [number]       := 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 
   * @endverbatim
   *
   * That is, the bin @f$ -3\le\eta\ge3@f$ is labeled
   * <b>m3d00_p3d00</b>, @f$ 0\le\eta\ge2.5@f$ is <b>p0d00_p2d50</b> 
   *
   * @a measuredFile and @a corrFile can point to the same file.  If
   * @a corrFile is not specified, it is assumed that @a measuredFile
   * has the expected @a corrFile @e in @e addition to the
   * expected content of that file.
   * 
   * @param measuredFile Name of file containing measured data
   * @param corrFile     Name of file containing correction data
   * @param method       Unfolding method to use 
   * @param regParam     Regularization parameter 
   */  
  void Run(const TString& measuredFile, const TString& corrFile,
	   const TString& method="Bayes", Double_t regParam=4)
  {
    // Get the input collections 
    if (measuredFile.IsNull()) {
      Error("Run", "No measurements given");
      return;
    }
    TCollection* mTop = GetTop(measuredFile, false);
    TCollection* cTop = GetTop(corrFile.IsNull() ? measuredFile : corrFile, 
			       true);
    if (!mTop || !cTop) return;

    // Get some info from the input collection 
    UShort_t sys;
    UShort_t sNN; 
    ULong_t  trig; 
    Double_t minZ; 
    Double_t maxZ;
    GetParameter(mTop, "sys",     sys);	  
    GetParameter(mTop, "sNN",     sNN);	  
    GetParameter(mTop, "trigger", trig);
    GetParameter(mTop, "minIpZ",  minZ); 
    GetParameter(mTop, "maxIpZ",  maxZ); 
    if (sys == 0 || sNN == 0) 
      Warning("Run", "System (%d) and/or collision energy (%d) unknown", 
	      sys, sNN);
    
    // Open the output file 
    TFile* out = TFile::Open("forward_unfolded.root", "RECREATE");
    if (!out) { 
      Error("Run", "Failed to open output file");
      return;
    }

    // Decode method option and store in file 
    TString meth(method);
    UInt_t  mId = MethodId(meth);
    if (mId == 0xDeadBeef) return;

    // Store information 
    SaveInformation(out,meth,mId,regParam,sys,sNN,trig,minZ,maxZ,
		    corrFile.IsNull());

    // Load other data 
    TString savPath(gROOT->GetMacroPath());
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2/scripts",
                             gROOT->GetMacroPath()));
    // Always recompile 
    if (!gROOT->GetClass("OtherPNch"))
      gROOT->LoadMacro("OtherPNchData.C++");
    gROOT->SetMacroPath(savPath);

    // Loop over the input 
    const char*  inputs[] = { "symmetric", "positive", "negative", 0 };
    const char** pinput   = inputs;
    while (*pinput) { 
      TCollection* mInput = GetCollection(mTop, *pinput, false);
      TCollection* cInput = GetCollection(cTop, *pinput, false);
      if (mInput && cInput)
	ProcessType(mInput, cInput, mId, regParam, out,
		    sys, sNN);
      pinput++;
    }      

    out->Write();
    // out->ls();
    out->Close();

    SaveSummarize();
  }
  /** 
   * Append an & to a string and the next term.
   * 
   * @param trg  Output string
   * @param what Term
   */
  static void AppendAnd(TString& trg, const TString& what)
  {
    if (!trg.IsNull()) trg.Append(" & ");
    trg.Append(what);
  }
  /** 
   * Store some information on the output
   * 
   * @param dir      Where to store
   * @param method   Method used
   * @param mId      Method identifier 
   * @param regParam Regularization parameter 
   * @param sys      Collision system
   * @param sNN      Center of mass energy 
   * @param trigger  Trigger mask 
   * @param minIpZ   Least z coordinate of interaction point
   * @param maxIpZ   Largest z coordinate of interaction point
   * @param self     Self-consistency check
   */
  void SaveInformation(TDirectory* dir, 
		       const TString& method,
		       Int_t          mId,
		       Double_t       regParam, 
		       UShort_t       sys, 
		       UShort_t       sNN,
		       UInt_t         trigger,
		       Double_t       minIpZ, 
		       Double_t       maxIpZ, 
		       Bool_t         self) const
  {
    dir->cd();

    TParameter<bool>* pM = new TParameter<bool>("self", self);
    pM->SetBit(BIT(19));
    pM->Write();

    TNamed* outMeth = new TNamed("method", method.Data());
    outMeth->SetUniqueID(mId);
    outMeth->Write();

    TParameter<double>* pR = new TParameter<double>("regParam", regParam);
    pR->SetBit(BIT(19));
    pR->Write();
    
    TString tS = (sys == 1 ? "pp" : sys == 2 ? "PbPb" : sys == 3 ? "pPb" : "?");
    TNamed* pS = new TNamed("sys", tS.Data()); pS->SetUniqueID(sys);
    pS->Write();
    
    TString tE;
    if      (sNN <  1000) tE = Form("%dGeV", sNN);
    else if (sNN <  3000) tE = Form("%4.2fTeV", float(sNN)/1000);
    else                  tE = Form("%dTeV", sNN/1000);
    TNamed* pE = new TNamed("sNN", tE.Data()); pE->SetUniqueID(sNN);
    pE->Write();
    
    TString tT;
    /** 
     * Bits of the trigger pattern
     */
    enum { 
      /** In-elastic collision */
      kInel        = 0x0001, 
      /** In-elastic collision with at least one SPD tracklet */
      kInelGt0     = 0x0002, 
      /** Non-single diffractive collision */
      kNSD         = 0x0004, 
      /** Empty bunch crossing */
      kEmpty       = 0x0008, 
      /** A-side trigger */
      kA           = 0x0010, 
      /** B(arrel) trigger */
      kB           = 0x0020, 
      /** C-side trigger */
      kC           = 0x0080,  
      /** Empty trigger */
      kE           = 0x0100,
      /** pileup from SPD */
      kPileUp      = 0x0200,    
      /** true NSD from MC */
      kMCNSD       = 0x0400,    
      /** Offline MB triggered */
      kOffline     = 0x0800,
      /** At least one SPD cluster */ 
      kNClusterGt0 = 0x1000,
      /** V0-AND trigger */
      kV0AND       = 0x2000, 
      /** Satellite event */
      kSatellite   = 0x4000
    };
    if ((trigger & kInel)        != 0x0) AppendAnd(tT, "INEL");
    if ((trigger & kInelGt0)     != 0x0) AppendAnd(tT, "INEL>0");
    if ((trigger & kNSD)         != 0x0) AppendAnd(tT, "NSD");
    if ((trigger & kV0AND)       != 0x0) AppendAnd(tT, "V0AND");
    if ((trigger & kA)           != 0x0) AppendAnd(tT, "A");
    if ((trigger & kB)           != 0x0) AppendAnd(tT, "B");
    if ((trigger & kC)           != 0x0) AppendAnd(tT, "C");
    if ((trigger & kE)           != 0x0) AppendAnd(tT, "E");
    if ((trigger & kMCNSD)       != 0x0) AppendAnd(tT, "MCNSD");
    if ((trigger & kNClusterGt0) != 0x0) AppendAnd(tT, "NCluster>0");
    if ((trigger & kSatellite)   != 0x0) AppendAnd(tT, "Satellite");
    TNamed* pT = new TNamed("trigger", tT.Data()); pT->SetUniqueID(trigger);
    pT->Write();
    
    TParameter<double>* pY = new TParameter<double>("minIpZ", minIpZ);
    pY->SetBit(BIT(19));
    pY->Write();

    TParameter<double>* pZ = new TParameter<double>("maxIpZ", maxIpZ);
    pZ->SetBit(BIT(19));
    pZ->Write();
  }
  /** 
   * Save a script to do a summary of this step 
   * 
   */		       
  void SaveSummarize()
  {
    std::ofstream f("SummarizeUnfold.C");
    f << "void SummarizeUnfold(const char* file=\"forward_unfolded.root\")\n"
      << "{\n"
      << "  const char* fwd=\"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/DrawUnfoldedSummary.C\",fwd));\n"
      << "  DrawUnfoldedSummary(file);\n"
      << "}\n"
      << "// EOF" << std::endl;
    f.close();
  }
  /** 
   * Process a single type - i.e., one of <i>symmetric</i>,
   * <i>positive</i>, <i>negative</i>, or <i>other</i> - by looping
   * over all contained objects and call ProcessBin for each found
   * bin.
   * 
   * @param measured     Input collection of measured data
   * @param corrections  Input collection of correction data
   * @param method       Unfolding method to use 
   * @param regParam     Regularisation parameter
   * @param out          Output directory. 
   * @param sys          Collision system
   * @param sNN          Collision energy 
   */
  void ProcessType(TCollection* measured, 
		   TCollection* corrections,
		   UInt_t       method,
		   Double_t     regParam,
		   TDirectory*  out,
		   UShort_t     sys, 
		   UShort_t     sNN)
  {
    Printf("  Processing %s ...", measured->GetName());
    TDirectory* dir = out->mkdir(measured->GetName());
    
    // Make some summary stacks 
    THStack*  allMeasured  = new THStack("measured",      
					 "Measured P(#it{N}_{ch})");
    THStack*  allTruth     = new THStack("truth",        
					 "MC 'truth' P(#it{N}_{ch})");
    THStack*  allTruthA    = new THStack("truthAccepted",
					 "MC 'truth' accepted P(#it{N}_{ch})");
    THStack*  allUnfolded  = new THStack("unfolded",
					 "Unfolded P(#it{N}_{ch})");
    THStack*  allCorrected = new THStack("corrected",
					 "Corrected P(#it{N}_{ch})");
    THStack*  allRatio     = (sys != 1 ? 0 : 
			      new THStack("ratios", "Ratios to other"));
    TMultiGraph* allALICE  = (sys != 1 ? 0 : 
			      new TMultiGraph("alice", "ALICE Published"));
    TMultiGraph* allCMS    = (sys != 1 ? 0 : 
			      new TMultiGraph("cms", "CMS Published"));

    // Loop over the list of objects. 
    static TRegexp regex("[pm][0-9]d[0-9]*_[pm][0-9]d[0-9]*");
    TIter          next(measured);
    TObject*       o = 0;
    Int_t          i = 0;
    Double_t       r = regParam;
    while ((o = next())) {
      // Go back to where we where 
      dir->cd();
      
      // if not a collection, don't bother 
      if (!o->IsA()->InheritsFrom(TCollection::Class())) continue;
    
      // If it doesn't match our regular expression, don't bother 
      TString n(o->GetName());
      if (n.Index(regex) == kNPOS) { 
	// Warning("ScanType", "%s in %s doesn't match eta range regexp", 
	//         n.Data(), real->GetName());
	continue;
      }
      TCollection* mBin = static_cast<TCollection*>(o);
      TCollection* cBin = GetCollection(corrections, n.Data());
      if (!cBin) continue;

      THStack* binS = ProcessBin(mBin, cBin, method, r, dir);
      if (!binS) continue;

      TH1* result = 0;
      Bin2Stack(binS, i, allMeasured, allTruth, allTruthA, 
		allUnfolded, allCorrected, result);

      TGraph* alice = 0;
      TGraph* cms   = 0;
      Other2Stack(o->GetName(), i, sNN, allALICE, allCMS, alice, cms);
      Ratio2Stack(i, result, alice, cms, allRatio);
      i++;
    }
    dir->Add(allMeasured);
    dir->Add(allTruth);
    dir->Add(allTruthA);
    dir->Add(allUnfolded);
    dir->Add(allCorrected);
    if (allALICE && allALICE->GetListOfGraphs()) {
      if (allALICE->GetListOfGraphs()->GetEntries() > 0)
	dir->Add(allALICE);
      else 
	delete allALICE;
    }
    if (allCMS && allCMS->GetListOfGraphs()) {
      if (allCMS->GetListOfGraphs()->GetEntries() > 0) 
	dir->Add(allCMS);
      else 
	delete allCMS;
    }
    if (allRatio && allRatio->GetHists()) { 
      if (allRatio->GetHists()->GetEntries() > 0) 
	dir->Add(allRatio);
      else 
	delete allRatio;
    }
  }
  /** 
   * Process a single eta bin 
   * 
   * @param measured     Input collection of measured data
   * @param corrections  Input collection of correction data
   * @param method       Unfolding method to use 
   * @param regParam     Regularisation parameter
   * @param out          Output directory. 
   *
   * @return Stack of histograms or null 
   */
  THStack* ProcessBin(TCollection* measured, 
		      TCollection* corrections, 
		      UInt_t       method,
		      Double_t     regParam, 
		      TDirectory*  out)
  {
    Printf("   Processing %s ...", measured->GetName());
    // Try to get the data 
    TH1* inRaw    = GetH1(measured,    "rawDist");
    TH1* inTruth  = GetH1(corrections, "truth");
    TH1* inTruthA = GetH1(corrections, "truthAccepted");
    TH1* inTrgVtx = GetH1(corrections, "triggerVertex");
    TH2* inResp   = GetH2(corrections, "response");
    if (!inRaw || !inTruth || !inTruthA || !inTrgVtx || !inResp) 
      return 0;
    
    // Make output directory
    TDirectory* dir = out->mkdir(measured->GetName());
    dir->cd();

    // Copy the input to the output 
    TH1* outRaw    = static_cast<TH1*>(inRaw    ->Clone("measured"));
    TH1* outTruth  = static_cast<TH1*>(inTruth  ->Clone("truth"));
    TH1* outTruthA = static_cast<TH1*>(inTruthA ->Clone("truthAccepted"));
    TH1* outTrgVtx = static_cast<TH1*>(inTrgVtx ->Clone("triggerVertex"));
    TH2* outResp   = static_cast<TH2*>(inResp   ->Clone("response"));

    // Make our response matrix 
    RooUnfoldResponse matrix(0, 0, inResp);
    
    // Store regularization parameter 
    Double_t             r        = regParam;
    RooUnfold::Algorithm algo     = (RooUnfold::Algorithm)method;
    RooUnfold*           unfolder = RooUnfold::New(algo, &matrix, inRaw, r);
    unfolder->SetVerbose(0);

    // Do the unfolding and get the result
    TH1* res = unfolder->Hreco();
    res->SetDirectory(0);

    // Make a copy to store on the output 
    TH1* outUnfold = static_cast<TH1*>(res->Clone("unfolded"));
    TString tit(outUnfold->GetTitle());
    tit.ReplaceAll("Unfold Reponse matrix", "Unfolded P(#it{N}_{ch})");
    outUnfold->SetTitle(tit);

    // Clone the unfolded results and divide by the trigger/vertex
    // bias correction
    TH1* outCorr   = static_cast<TH1*>(outUnfold->Clone("corrected"));
    outCorr->Divide(inTrgVtx);
    tit.ReplaceAll("Unfolded", "Corrected");
    outCorr->SetTitle(tit);

    // Now normalize the output to integral=1 
    TH1*  hists[] = { outRaw, outUnfold, outCorr, 0 };
    TH1** phist   = hists;
    while (*phist) { 
      TH1* h = *phist;
      if (h) { 
	Double_t intg = h->Integral(1, h->GetXaxis()->GetXmax());
	h->Scale(1. / intg, "width");
      }
      phist++;
    }
    
    // And make ratios
    TH1* ratioTrue = static_cast<TH1*>(outCorr->Clone("ratioCorrTruth"));
    tit = ratioTrue->GetTitle();
    tit.ReplaceAll("Corrected", "Corrected/MC 'truth'");
    ratioTrue->SetTitle(tit);
    ratioTrue->Divide(outTruth);
    ratioTrue->SetYTitle("P_{corrected}(#it{N}_{ch})/P_{truth}(#it{N}_{ch})");

    TH1* ratioAcc  = static_cast<TH1*>(outUnfold->Clone("ratioUnfAcc"));
    tit = ratioAcc->GetTitle();
    tit.ReplaceAll("Unfolded", "Unfolded/MC selected");
    ratioAcc->SetTitle(tit);
    ratioAcc->Divide(outTruthA);
    ratioAcc->SetYTitle("P_{unfolded}(#it{N}_{ch})/P_{MC}(#it{N}_{ch})");
    

    // Make a stack 
    tit = measured->GetName();
    tit.ReplaceAll("m", "-");
    tit.ReplaceAll("p", "+");
    tit.ReplaceAll("d", ".");
    tit.ReplaceAll("_", "<#it{#eta}<");
    THStack* stack = new THStack("all", tit);
    stack->Add(outTruth,  "E2");
    stack->Add(outTruthA, "E2");
    stack->Add(outRaw,    "E1");
    stack->Add(outUnfold, "E1");
    stack->Add(outCorr,   "E1");
    dir->Add(stack);

    // Rest of the function is devoted to making the output look nice 
    outRaw   ->SetDirectory(dir); 
    outTruth ->SetDirectory(dir);  
    outTruthA->SetDirectory(dir);  
    outTrgVtx->SetDirectory(dir);  
    outResp  ->SetDirectory(dir);  
    outUnfold->SetDirectory(dir);   
    outCorr  ->SetDirectory(dir); 

    outRaw   ->SetMarkerStyle(20);  // Measured is closed
    outTruth ->SetMarkerStyle(24);  // MC is open
    outTruthA->SetMarkerStyle(24);  // MC is open
    outTrgVtx->SetMarkerStyle(20);  // Derived is closed
    outUnfold->SetMarkerStyle(20);  // Derived is closed   
    outCorr  ->SetMarkerStyle(20);  // Derived is closed 

    outRaw   ->SetMarkerSize(0.9); 
    outTruth ->SetMarkerSize(1.6);  
    outTruthA->SetMarkerSize(1.4);  
    outTrgVtx->SetMarkerSize(1.0);  
    outUnfold->SetMarkerSize(0.9);   
    outCorr  ->SetMarkerSize(1.0);
 
    outRaw   ->SetMarkerColor(kColorMeasured); 
    outTruth ->SetMarkerColor(kColorTruth);  
    outTruthA->SetMarkerColor(kColorAccepted);  
    outTrgVtx->SetMarkerColor(kColorTrgVtx);  
    outUnfold->SetMarkerColor(kColorUnfolded);   
    outCorr  ->SetMarkerColor(kColorCorrected); 

    outRaw   ->SetFillColor(kColorError);     
    outTruth ->SetFillColor(kColorError);  
    outTruthA->SetFillColor(kColorError);  
    outTrgVtx->SetFillColor(kColorError);  
    outUnfold->SetFillColor(kColorError);   
    outCorr  ->SetFillColor(kColorError); 

    outRaw   ->SetFillStyle(0); 
    outTruth ->SetFillStyle(1001);  
    outTruthA->SetFillStyle(1001);  
    outTrgVtx->SetFillStyle(0);  
    outUnfold->SetFillStyle(0);   
    outCorr  ->SetFillStyle(0);

    outRaw   ->SetLineColor(kBlack); 
    outTruth ->SetLineColor(kBlack);  
    outTruthA->SetLineColor(kBlack);  
    outTrgVtx->SetLineColor(kBlack);  
    outUnfold->SetLineColor(kBlack);   
    outCorr  ->SetLineColor(kBlack); 

    // Legend 
    TLegend* l = StackLegend(stack);
    l->AddEntry(outRaw,     "Raw",                 "lp");
    l->AddEntry(outTruth,   "MC 'truth'",          "fp");
    l->AddEntry(outTruthA,  "MC 'truth' accepted", "fp");
    l->AddEntry(outUnfold,  "Unfolded",            "lp");
    l->AddEntry(outCorr,    "Corrected",           "lp");

    return stack;
  }
  static void BinAttributes(Int_t i,
			    Int_t&    open, 
			    Int_t&    closed, 
			    Float_t&  size,
			    Double_t& factor) 
  {
    // --- Setup for markers -----------------------------------------
    const Int_t   nMarkers = 7;
    const Int_t   aClosed[] = { 20,  21,  22,  23,  29,  33,  34  };
    const Int_t   aOpen[]   = { 24,  25,  26,  32,  30,  27,  28  };
    const Float_t aSize[]   = { 1.1, 1.0, 1.2, 1.2, 1.2, 1.2, 1.0 };
    Int_t         j         = i % nMarkers;

    open   = aOpen[j];
    closed = aClosed[j];
    size   = aSize[j];
    factor = TMath::Power(10, i);
  }
  /** 
   * Add the bin histograms to our summary stacks 
   * 
   * @param bin       Bin stack
   * @param i         Current off-set in the stacks 
   * @param measured  All measured @f$ P(N_{ch})@f$ 
   * @param truth     All MC truth @f$ P(N_{ch})@f$ 
   * @param accepted  All MC accepted @f$ P(N_{ch})@f$ 
   * @param unfolded  All unfolded @f$ P(N_{ch})@f$ 
   * @param corrected All corrected @f$ P(N_{ch})@f$ 
   * @param result    The result in this bin
   */
  void Bin2Stack(const THStack* bin, Int_t i, 
		 THStack* measured, 
		 THStack* truth, 
		 THStack* accepted, 
		 THStack* unfolded,
		 THStack* corrected,
		 TH1*&    result)
  {
    Int_t open, closed;
    Double_t factor; 
    Float_t  size;
    BinAttributes(i, open, closed, size, factor);

    TIter next(bin->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      THStack* tmp = 0;
      Int_t    col = h->GetMarkerColor();
      Int_t    sty = 0;
      switch (col) { 
      case kColorMeasured:  tmp = measured;   sty = closed;  break;
      case kColorTruth:     tmp = truth;      sty = open;    break;
      case kColorAccepted:  tmp = accepted;   sty = open;    break;
      case kColorUnfolded:  tmp = unfolded;   sty = closed;  break;
      case kColorCorrected: tmp = corrected;  sty = closed;  break;
      default: continue; 
      }
      // Now clone, and add to the appropriate stack 
      TH1* cln = static_cast<TH1*>(h->Clone(h->GetName()));
      cln->SetDirectory(0);
      cln->SetMarkerStyle(sty);
      cln->SetMarkerSize(size);
      cln->Scale(factor); // Scale by 10^i
      if (col == kColorCorrected) result = cln;

      // Make sure we do not get the old legend 
      TObject* tst = cln->FindObject("legend");
      if (tst) cln->GetListOfFunctions()->Remove(tst);

      tmp->Add(cln, next.GetOption());
    }
    
    // Add entries to our stacks 
    TString   txt      = bin->GetTitle();
    if      (i == 0) txt.Append(" (#times1)");
    else if (i == 1) txt.Append(" (#times10)");
    else             txt.Append(Form(" (#times10^{%d})", i));
    THStack*  stacks[] = { measured, truth, accepted, unfolded, corrected, 0 };
    THStack** pstack   = stacks;
    while (*pstack) { 
      TLegend* leg = StackLegend(*pstack);
      pstack++;
      if (!leg) continue;
      
      TObject* dummy = 0;
      TLegendEntry* e = leg->AddEntry(dummy, txt, "p");
      e->SetMarkerStyle(closed);
      e->SetMarkerSize(1.2*size);
      e->SetMarkerColor(kBlack);
      e->SetFillColor(0);
      e->SetFillStyle(0);
      e->SetLineColor(kBlack);
    }
  }
  /** 
   * Add distributions from other experiments to stacks 
   * 
   * @param name     Name of current bin 
   * @param i        Index of current bin
   * @param sNN      Center of mass energy 
   * @param allALICE Stack of ALICE data 
   * @param allCMS   Stack of CMS data 
   * @param alice    Possible ALICE result on return
   * @param cms      Possible CMS result on return
   */
  void Other2Stack(const TString& name, Int_t i,
		   UShort_t sNN, TMultiGraph* allALICE, TMultiGraph* allCMS,
		   TGraph*& alice, TGraph*& cms) 
  {
    if (!allALICE && !allCMS) return;

    TString tmp(name);
    tmp.ReplaceAll("p", "+");
    tmp.ReplaceAll("m", "-");
    tmp.ReplaceAll("_", " ");
    tmp.ReplaceAll("d", ".");
    TObjArray* tokens = tmp.Tokenize(" ");
    if (!tokens || tokens->GetEntriesFast() < 2) { 
      Error("Other2Stack", "Failed to decode eta range from %s", name.Data());
      if (tokens) tokens->Delete();
      return;
    }
    Double_t eta1 = static_cast<TObjString*>(tokens->At(0))->String().Atof();
    Double_t eta2 = static_cast<TObjString*>(tokens->At(1))->String().Atof();
    tokens->Delete();
    
    if (TMath::Abs(eta2+eta1) > 1e-3) {
      // Not symmetric bin 
      // Info("Other2Stack", "bin [%f,%f] is not symmetric (%f)",
      //      eta1, eta2, TMath::Abs(eta2-eta1));
      return;
    }
    Double_t aEta = TMath::Abs(eta1);

    Int_t open, closed;
    Double_t factor; 
    Float_t  size;
    BinAttributes(i, open, closed, size, factor);

    if (allALICE) {
      TGraphAsymmErrors* g = GetOther(1, aEta, sNN, factor);
      if (g) {
	g->SetMarkerStyle(closed);
	g->SetMarkerColor(kColorALICE);
	g->SetMarkerSize(size);
	allALICE->Add(g, "p same");
	alice = g;
      }
    }
    if (allCMS) {
      TGraphAsymmErrors* g = GetOther(0, aEta, sNN, factor);
      if (g) {
	g->SetMarkerStyle(closed);
	g->SetMarkerColor(kColorCMS);
	g->SetMarkerSize(size);
	allCMS->Add(g, "p same");
	cms = g;
      }
    }
  }
  /** 
   * Create ratios to other data 
   * 
   * @param ib      Bin number  
   * @param res     Result
   * @param alice   ALICE result if any
   * @param cms     CMS result if any
   * @param all     Stack to add ratio to 
   */
  void Ratio2Stack(Int_t ib, TH1* res, TGraph* alice, TGraph* cms, THStack* all)
  {
    if (!all || !res || !(alice || cms)) return;

    Int_t        off  = 5*ib;
    TGraph*      gs[] = { (alice ? alice : cms), (alice ? cms : 0), 0 };
    TGraph**     pg   = gs;
    while (*pg) { 
      TGraph*     g = *pg;
      const char* n = (g == alice ? "ALICE" : "CMS");

      TH1*    r = static_cast<TH1*>(res->Clone(Form("ratio%s", n)));
      TString tit(r->GetTitle());
      tit.ReplaceAll("Corrected", Form("Ratio to %s", n));
      r->SetTitle(tit);
      r->SetMarkerColor(g->GetMarkerColor());
      r->SetLineColor(g->GetLineColor());

      TObject* tst = r->FindObject("legend");
      if (tst) r->GetListOfFunctions()->Remove(tst);

      for (Int_t i = 1; i <= r->GetNbinsX(); i++) {
	Double_t c = r->GetBinContent(i);
	Double_t e = r->GetBinError(i);
	Double_t o = g->Eval(r->GetBinCenter(i));
	if (o < 1e-12) { 
	  r->SetBinContent(i, 0);
	  r->SetBinError(i, 0);
	  continue;
	}
	r->SetBinContent(i, (c - o) / o + off);
	r->SetBinError(i, e / o);
      }
      all->Add(r);
      pg++;
    }
    TLegend* leg = StackLegend(all);
    if (!leg) return;
      
    TString   txt      = res->GetTitle();
    txt.ReplaceAll("Corrected P(#it{N}_{ch}) in ", "");
    if      (ib == 0) txt.Append(" "); // (#times1)");
    // else if (ib == 1) txt.Append(" (#times10)");
    else              txt.Append(Form(" (+%d)", off));

    TObject* dummy = 0;
    TLegendEntry* e = leg->AddEntry(dummy, txt, "p");
    e->SetMarkerStyle(res->GetMarkerStyle());
    e->SetMarkerSize(res->GetMarkerSize());
    e->SetMarkerColor(kBlack);
    e->SetFillColor(0);
    e->SetFillStyle(0);
    e->SetLineColor(kBlack);
  }

  /** 
   * Get or create a stack legend.  This is done by adding a TLegend
   * object to the list of functions for the first histogram in the
   * stack.
   * 
   * @param stack Stack to get the legend from/modify 
   * 
   * @return The legend object or null
   */
  TLegend* StackLegend(THStack* stack) 
  {
    TList* hists = stack->GetHists();
    if (!hists) return 0;
    
    TObject* first = hists->First();
    if (!first) return 0;

    TH1*    hist = static_cast<TH1*>(first);
    TList*  list = hist->GetListOfFunctions();
    TObject* o   = list->FindObject("legend");
    if (o) return static_cast<TLegend*>(o);
    
    TLegend* l = new TLegend(0.65, 0.65, 0.9, 0.9, "", "NDC");
    l->SetName("legend");
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    list->Add(l);
    
    return l;
  }
  
  /* =================================================================
   *
   * Measurements from other sources, such as published ALICE, CMS, ...
   */
  TGraphAsymmErrors* GetOther(UShort_t type, Double_t eta, UShort_t sNN,
			      Int_t factor)
  {
    TString oC = Form("OtherPNch::GetData(%hu,%f,%hu);", 
                      type, eta, sNN);
    TGraphAsymmErrors* g = 
      reinterpret_cast<TGraphAsymmErrors*>(gROOT->ProcessLine(oC));
    if (!g) { 
      // Warning("GetOther", "No other data found for type=%d eta=%f sNN=%d",
      //         type, eta, sNN);
      return 0;
    }
  

    for (Int_t j = 0; j < g->GetN(); j++) { 
      g->SetPoint(j, g->GetX()[j], g->GetY()[j]*factor);
      g->SetPointError(j, g->GetEXlow()[j], g->GetEXhigh()[j], 
		       g->GetEYlow()[j]*factor, g->GetEYhigh()[j]*factor);
    }
    return g;
  }    
};

void
UnfoldMultDists(const TString& method="Bayes", 
		Double_t       regParam=1e-30,
		const TString& measured="forward_multdists.root", 
		const TString& corrections="")
{
  Unfolder m;
  m.Run(measured, corrections, method, regParam);
}
