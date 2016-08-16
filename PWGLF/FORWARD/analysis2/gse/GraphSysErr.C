/**
 * @file   GraphSysErr.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Aug 14 13:51:05 2014
 * 
 * @brief  An (X,Y) graph with configurable errors 
 * 
 * Copyright (c) 2014 Christian Holm Christensen. 
 *
 * Licensed under the LGPL-3 
 */
# include <TNamed.h>
# include <TAttMarker.h>
# include <TAttLine.h>
# include <TAttFill.h>
#ifndef GraphSysErr_C
# define GraphSysErr_C
# include <TList.h>
# ifndef __CINT__
#  include <TStyle.h>
#  include <TGraph.h>
#  include <TGraphErrors.h>
#  include <TGraphAsymmErrors.h>
#  include <TMultiGraph.h>
#  include <TH1.h>
#  include <TMath.h>
#  include <TFitResultPtr.h>
#  include <TF1.h>
#  include <TROOT.h>
#  include <TSystem.h>
#  include <TDatime.h>
#  include <TList.h>
#  include <TBrowser.h>
#  include <TPRegexp.h>
#  include <TClonesArray.h>
#  include <TArrayL64.h>
#  include <TLine.h>
#  include <TClass.h>
#  include <TObjString.h>
#  include <iostream>
#  include <iomanip>
#  include <fstream>
#  include <sstream>
// #  define TOKEN(C,I) static_cast<TObjString*>(C->At(I))->String()
#  define INCOMPAT_CMN_AS_QUAL 1
# else 
# include <iosfwd>
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TMultiGraph;
class TFitResultPtr;
class TF1;
class TList;
class TBrowser;
class TArrayD;
class TClonesArray;
class TLine;
# endif


/** 
 * @defgroup gse Graphs with systematic errors 
 *
 * This module defines a class for @f$ (x,y)@f$ graphs with an
 * arbitrary set of systematic uncertainties.  It provides storage and
 * visualisation of such data, as well as import/export to the HepData
 * format.
 */
/** 
 * @page gse_an_example Example 
 * @ingroup gse
 *
 * @dontinclude Example.C 
 *
 * We define a function in a script.  It takes one argument - whether
 * to fit the data or not.
 * 
 * @skip void
 * @until {
 *
 * The first thing we do, is to check if we have the class
 * GraphSysErr, and if not, we compile the script.
 *
 * @until LoadMacro
 *
 * Then, we adjust some sizes - the default error along X, how large
 * the end bars should be.
 *
 * @until SetEndErrorSize 
 *
 * Now we're ready to declare our object.
 *
 * @until GraphSysErr
 *
 * The first we do, is set some parameters on the object: That we
 * shouldn't put tick marks ('ends') on the error, and the name of the
 * X and Y axis.
 *
 * @until SetYTitle
 *
 * Next, we set some key values.  These are mainly for exporting to a
 * Durham database input file. In principle one can define any key,
 * but only a sub-set is written when exporting.
 *
 * @until "obskey"
 *
 * Next, we can add as many qualifiers to the data set we want.  If
 * one is to export multiple graphs to a table, one should set a @b RE
 * qualifier to distinguish the columns, or give each graph a distinct
 * title.
 *
 * @until AddQualifier 
 *
 * Now we can define common errors. That is, errors that are
 * correlated over all points.  Common errors only have one negative
 * and one positive value (which may be identical).  Common errors can
 * be relative.  The member function GraphSysErr:DefineCommon returns
 * a unique identifier.  This identifier should be stored for later
 * manipulation of the error.
 *
 * @until cm2 
 *
 * Next is the point-to-point (uncorrelated) systematic errors.  Here,
 * we simple give a name and whether the error is relative.  The
 * member function GraphSysErr:DeclarePoint2Point returns a unique
 * identifier.  This identifier should be stored for later
 * manipulation of the error.
 *
 * @until pp2 
 *
 * Now we customize the graphical output of each error 
 *
 * @until pp2, GraphSysErr::kHat
 *
 * In this example, we take the data points to from a Gaussian
 * deviate.  Technically, we make a histogram of the probability of a
 * given number.
 *
 * @until Scale
 *
 * Now we can set all our points.  We remove our temporary histogram
 * after we're done with it.
 * 
 * @until delete 
 *
 * Finally, we can build a canvas 
 *
 * @until SetRightMargin
 *
 * Depending on the single parameter, we either draw or fit a Gaussian
 * distribtion to the data.
 * 
 * @until gaus
 *
 * And then we draw and finish
 *
 * @until }
 * 
 */
/**
 * This class defines an (X,Y) with any number of error sources.  
 * 
 * Sources that can be specified are 
 *
 *   - 1 set Statistical errors 
 *   - N sets of common systematic errors 
 *   - M sets of point-to-point systematic errors. 
 *
 * Systematic errors can be defined to relative to the point value or
 * absolute numbers.
 *
 * There are various options for drawing this data set (see Draw).  A
 * function can also be fitted to the data set, talking various kinds
 * of errors into consideration (see Fit).  The data set can be export
 * to a format more or less acceptable by the Durham database (see
 * Export), and one can import data sets from Durham database input
 * formatted files (see Import)
 *
 * @see @link gse_an_example Example @endlink
 *
 * @ingroup gse 
 */
class GraphSysErr : public TNamed, public TAttMarker, 
		    public TAttLine, public TAttFill 
{
public:
  enum { 
    kDraw    = 0x1,
    kImport  = 0x2, 
    kExport  = 0x4,
    kRatio   = 0x8,
    kVerbose = 0  // Set to OR of above to enable verbose/debug output 
  };
  /** A short-hand type definition */
  typedef TGraphAsymmErrors Graph;
  /** 
   * Drawing options.  We re-encode them here as distinct enums. 
   */
  enum EDrawOption_t { 
    kNormal  = 0, //     - Line with ticks
    kNoTick  = 1, // Z0  - Line with no ticks 
    kArrow   = 2, // >0  - Linw with arrows 
    kRect    = 3, // 20  - Rectangle w/fill
    kBox     = 4, // 50  - Rectangle w/fill & line 
    kFill    = 5, // 30  - Filled area 
    kCurve   = 6, // 40  - Filled smoothed area 
    kHat     = 7, // []0 - Hats 
    kBar     = 8, // ||0 - A bar 
    kNone    = 9, // XP  - No errors
    kLine    = 10, // XC
    kConnect = 11 // XL
  };
  /** 
   * Types of @f$\chi^2@f$ comparisons.  See also 
   *
   * https://root.cern.ch/root/htmldoc/TH1#TH1:Chi2Test
   */
  enum EChi2Type {
    kExperimentExperiment, // Ignore errors on both
    kExperimentModel,      // Ignore errors on second
    kModelModel,           // Use errors on both
    kExperimentTruth       // Use errors on first, take second to be exact.
  };
  /** 
   * Options for ratios. 
   *
   * - kChi2 and kKolomogorov are mutally exclusive.  If none of them
   *   are given, no test is calculated. 
   * - kProp and kMax are mutually exclusive.  The default is kProp
   *
   */
  enum ERatioOption {
    kMax          = 0x00001, // Error is largest relative error
    kMin          = 0x00002, // Error is largest relative error
    kCancel       = 0x00004,
    kDenomRel     = 0x00008,
    kRatioDefault = kCancel
  };
  enum {
    kUsedBit = (1 << 18),
    kOnlyWeightBit = (1 << 19)
  };
  /** 
   * @{ 
   * @name Allocation, dealloction, copy, and assignment 
   */
  /** 
   * Default CTOR - use only for I/O
   */
  GraphSysErr() 
    : TNamed(), 
      TAttMarker(1,20,1), 
      TAttLine(1,1,1),
      TAttFill(0,0),
      fPoint2Point(), 
      fCommon(), 
      fData(0), 
      fDrawn(0),
      fCounter(0),
      fSumFill(0,0),
      fSumLine(1,1,1),
      fSumTitle(),
      fSumOption(0),
      fCommonSumFill(0,0),
      fCommonSumLine(1,1,1),
      fCommonSumTitle(),
      fCommonSumOption(0),
      fDataOption(0), 
      fXTitle(""),
      fYTitle(""),
      fMap(0),
      fQualifiers(0),
      fStatRelative(false)
  {
    
  }
  /** 
   * CTOR with number of points
   * 
   * @param n Number of points to pre-allocate 
   */
  GraphSysErr(Int_t n) 
    : TNamed("sysErrGraph","Data"), 
      TAttMarker(1,20,1), 
      TAttLine(1,1,1),
      TAttFill(0,0),
      fPoint2Point(), 
      fCommon(), 
      fData(0), 
      fDrawn(0),
      fCounter(0),
      fSumFill(0,0),
      fSumLine(1,1,1),
      fSumTitle("Uncorr.Errors"),
      fSumOption(0),
      fCommonSumFill(0,0),
      fCommonSumLine(1,1,1),
      fCommonSumTitle("Corr.Errors."),
      fCommonSumOption(0),
      fDataOption(0), 
      fXTitle(""),
      fYTitle(""),
      fMap(0),
      fQualifiers(0),
      fStatRelative(false)
  {
    fPoint2Point.SetName("point2point");
    fPoint2Point.SetOwner();
    fCommon.SetName("common");
    fCommon.SetOwner();
    MakeDataGraph(n);
  }
  /** 
   * Constructor with name, title, and optional pre-allocated size 
   * 
   * @param name   Name
   * @param title  Title 
   * @param n      Pre-allocated points
   */
  GraphSysErr(const char* name, const char* title, Int_t n=0)
    : TNamed(name,title), 
      TAttMarker(1,20,1), 
      TAttLine(1,1,1),
      TAttFill(0,0),
      fPoint2Point(), 
      fCommon(), 
      fData(0), 
      fDrawn(0),
      fCounter(0),
      fSumFill(0,0),
      fSumLine(1,1,1),
      fSumTitle("Errors"),
      fSumOption(0),
      fCommonSumFill(0,0),
      fCommonSumLine(1,1,1),
      fCommonSumTitle("Corr.Errors."),
      fCommonSumOption(0),
      fDataOption(0), 
      fXTitle(""),
      fYTitle(""),
      fMap(0),
      fQualifiers(0),
      fStatRelative(false)
  {
    fPoint2Point.SetOwner();
    fCommon.SetOwner();
    fCommon.SetName("common");
    fPoint2Point("point2point");
    MakeDataGraph(n);
  }
  /** 
   * Copy CTOR 
   * 
   * @param other Object to copy from 
   */
  GraphSysErr(const GraphSysErr& other) 
    : TNamed(other), 
      TAttMarker(other), 
      TAttLine(other),
      TAttFill(other),
      fPoint2Point(), 
      fCommon(), 
      fData(0),
      fDrawn(0),
      fCounter(0),
      fSumFill(other.fSumFill),
      fSumLine(other.fSumLine),
      fSumTitle(other.fSumTitle),
      fSumOption(other.fSumOption),
      fCommonSumFill(other.fCommonSumFill),
      fCommonSumLine(other.fCommonSumLine),
      fCommonSumTitle(other.fCommonSumTitle),
      fCommonSumOption(other.fCommonSumOption),
      fDataOption(other.fDataOption), 
      fXTitle(other.fXTitle),
      fYTitle(other.fYTitle),
      fMap(0),
      fQualifiers(0),
      fStatRelative(false)
  {
    if (other.fData) 
      fData = static_cast<Graph*>(other.fData->Clone());
 
    fPoint2Point.SetOwner();
    fCommon.SetOwner();
    fCommon.SetName(other.GetName());
    fPoint2Point(other.GetName());

    TIter nextC(&other.fCommon);
    HolderCommon* common = 0;
    while ((common = static_cast<HolderCommon*>(nextC()))) {
      fCommon.Add(new HolderCommon(*common));
      fCounter++;
    }

    TIter nextP(&other.fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(nextP()))) {
      fPoint2Point.Add(new HolderP2P(*p2p));
      fCounter++;
    }

    CopyKeys(&other, "ar");
    // Print("all");
  }
  /**
   * DTOR
   */
  virtual ~GraphSysErr()
  {
    fPoint2Point.Delete();
    fCommon.Delete();
    if (fData)         delete fData;
    if (fDrawn)        delete fDrawn;
    if (fMap)          delete fMap;
    if (fQualifiers)   delete fQualifiers;
  }
  /** 
   * Assignment operator
   * 
   * @param other Object to copy from 
   * 
   * @return reference to this 
   */
  GraphSysErr& operator=(const GraphSysErr& other) 
  {
    if (&other == this) return *this;
    
    other.TNamed::Copy(*this);
    other.TAttMarker::Copy(*this);
    other.TAttLine::Copy(*this);
    other.TAttFill::Copy(*this);
    fPoint2Point.Clear();
    fCommon.Clear();
    
    if (fData) delete fData;
    fData = 0;
    if (other.fData)
      fData = static_cast<Graph*>(other.fData->Clone());

    if (fDrawn) { delete fDrawn; fDrawn = 0; }

    fCounter = other.fCounter;

    other.fSumFill.Copy(fSumFill);
    other.fSumLine.Copy(fSumLine);
    fSumTitle     = other.fSumTitle;
    fSumOption    = other.fSumOption;
    
    other.fCommonSumFill.Copy(fCommonSumFill);
    other.fCommonSumLine.Copy(fCommonSumLine);
    fCommonSumTitle     = other.fCommonSumTitle;
    fCommonSumOption    = other.fCommonSumOption;
    
    fXTitle       = other.fXTitle;
    fYTitle       = other.fYTitle;
    fStatRelative = other.fStatRelative;

    TIter nextC(&other.fCommon);
    HolderCommon* common = 0;
    while ((common = static_cast<HolderCommon*>(nextC()))) 
      fCommon.Add(new HolderCommon(*common));

    TIter nextP(&other.fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(nextP()))) 
      fPoint2Point.Add(new HolderP2P(*p2p));    

    if (fMap) delete fMap;
    fMap = 0;
    if (other.fMap)
      fMap = static_cast<TList*>(other.fMap->Clone());

    if (fQualifiers) delete fQualifiers;
    fQualifiers = 0;
    if (other.fQualifiers)
      fQualifiers = static_cast<TList*>(other.fQualifiers->Clone());
    
    return *this;
  }
  /* @} */

  /** 
   * @{ 
   * @name TObject functions 
   */
  /** 
   * List the content 
   * 
   * @param option option (not used)
   */
  virtual void ls(Option_t* option="") const
  {
    Print(option);
  }
  /** 
   * Print this. 
   * 
   * @param option not used
   */
  virtual void Print(Option_t* option="R") const
  {
    gROOT->IndentLevel();
    std::cout << GetName() << ": " << GetTitle() << std::endl;

    TString opt(option);    
    opt.ToUpper();
    if (opt.IsNull()) return;

    Bool_t all  = opt.Contains("ALL");
    Bool_t keys = all || opt.Contains("KEY");
    Bool_t qual = all || opt.Contains("QUAL");
    Bool_t sys  = all || opt.Contains("SYS");
    Bool_t cmn  = sys || opt.Contains("COMMON");
    Bool_t p2p  = sys || opt.Contains("P2P");
    Bool_t poi  = all || opt.Contains("XY");
    Bool_t attr = all || opt.Contains("ATTR");
    Bool_t sum  = opt.Contains("SUM");
    
    gROOT->IncreaseDirLevel();
    if (fMap && keys) {
      gROOT->IndentLevel();
      std::cout << "Key/value pairs: " << std::endl;
      gROOT->IncreaseDirLevel();
      TIter nextKV(fMap);
      TObject* kv = 0;
      while ((kv = nextKV())) {
	gROOT->IndentLevel();
	std::cout << '"' << kv->GetName() << '"' << "\t"
		  << '"' << kv->GetTitle() << '"' << "\n";
      }
      gROOT->IndentLevel();
      std::cout << "\"XTitle\"\t\"" << fXTitle << "\"\n";
      gROOT->IndentLevel();
      std::cout << "\"YTitle\"\t\"" << fYTitle << "\"\n";      
      // fMap->Print(option);
      gROOT->DecreaseDirLevel();
    }
    if (fQualifiers && qual) {
      gROOT->IndentLevel();
      std::cout << "Qualifier pairs: " << std::endl;
      gROOT->IncreaseDirLevel();
      TIter nextQ(fQualifiers);
      TObject* q = 0;
      while ((q = nextQ())) {
	gROOT->IndentLevel();
	std::cout << '"' << q->GetName() << '"' << "\t"
		  << '"' << q->GetTitle() << '"' << "\n";
      }
      // fQualifiers->ls(option);
      gROOT->DecreaseDirLevel();
    }
    if (attr) {
      gROOT->IndentLevel();
      Printf(" [option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d marker (c/s/s):%3d/%2d/%f]",
	     fDataOption, GetLineColor(), GetLineStyle(), GetLineWidth(),
	     GetFillColor(), GetFillStyle(),
	     GetMarkerColor(), GetMarkerStyle(), GetMarkerSize());
      gROOT->IndentLevel();
      Printf(" [sum title: %s "
	     "option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d",
	     fSumTitle.Data(), fSumOption,
	     fSumLine.GetLineColor(),
	     fSumLine.GetLineStyle(),
	     fSumLine.GetLineWidth(), 
	     fSumFill.GetFillColor(),
	     fSumFill.GetFillStyle());
      gROOT->IndentLevel();
      Printf(" [common sum title: %s "
	     "option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d",
	     fCommonSumTitle.Data(),
	     fCommonSumOption,
	     fCommonSumLine.GetLineColor(),
	     fCommonSumLine.GetLineStyle(),
	     fCommonSumLine.GetLineWidth(), 
	     fCommonSumFill.GetFillColor(),
	     fCommonSumFill.GetFillStyle());
    }

    if (cmn) {
      gROOT->IndentLevel();
      std::cout << "Commons: " << std::endl;
      gROOT->IncreaseDirLevel();
      TIter nextC(&fCommon);
      TObject* c = 0;
      while ((c = nextC())) {
	gROOT->IndentLevel();
	c->Print(option);
      }
      // fCommon.Print(option);
      gROOT->DecreaseDirLevel();
    }

    if (p2p) {
      gROOT->IndentLevel();
      std::cout <<  "Point-to-point: " << std::endl;
      gROOT->IncreaseDirLevel();
      TIter nextP(&fPoint2Point);
      TObject* p = 0;
      while ((p = nextP())) {
	gROOT->IndentLevel();
	p->Print(option);
      }
      // fPoint2Point.Print(option);
      gROOT->DecreaseDirLevel();
    }
    if (poi) {
      // gROOT->IndentLevel();
      gROOT->IndentLevel();
      std::cout <<  "Points: " << std::endl;
      gROOT->IncreaseDirLevel();
      TString hs;
      hs.Form(" %3s: %9s %9s %9s -> %9s %9s %9s",
	      "#", "X", "-dX", "+dX", "Y", "-stat", "+stat");
      gROOT->IndentLevel();
      std::cout << hs;
      if (sum) {
	hs.Form(" %9s %9s", "-sys", "+sys");
	std::cout << hs;
      }
      else {
	TIter nextP(&fPoint2Point);
	HolderP2P* p = 0;
	while ((p = static_cast<HolderP2P*>(nextP()))) {
	  Bool_t      rel  = p->IsRelative();
	  const char* post = (rel ? "%" : " ");
	  TString nm(p->GetTitle());
	  if (nm.Length() > 7) {
	    nm.Remove(4,nm.Length()-4);
	    nm.Append("...");
	  }
	  hs.Form(" -%8s%s +%8s%s", nm.Data(), post, nm.Data(), post);
	  std::cout << hs;	  
	}
      }
      std::cout << std::endl;
      for (Int_t i = 0; i < GetN(); i++) {
	gROOT->IndentLevel();
	Double_t eyl, eyh, wyl, wyh;
	// point,stat,common,quad,nosqrt,eyl,eyh
	Double_t y = GetYandError(i,false,false,true,false,eyl,eyh,wyl,wyh);
	// Double_t y = GetY(i);
	TString s;
	s.Form(" %3d: %+8f -%7f +%7f -> %+8f -%7f +%7f",
	       i, GetX(i), GetErrorXLeft(i), GetErrorXRight(i),
	       y, GetStatErrorDown(i), GetStatErrorUp(i));
	std::cout << s;
	if (sum) {
	  TString sSum;
	  sSum.Form(" -%7f +%7f", eyl, eyh);
	  std::cout << sSum;
	}
	// if (sum && p2p) std::cout << " [";
	else if (p2p) {
	  TIter nextP(&fPoint2Point);
	  HolderP2P* p = 0;
	  while ((p = static_cast<HolderP2P*>(nextP()))) {
	    Bool_t      rel  = p->IsRelative();
	    Double_t    fac  = (rel ? (TMath::Abs(y) > 1e-9 ? 100/y : 0) : 1);
	    const char* post = (rel ? "%" : " ");
	    s.Form(" -%7f%s +%7f%s",
		   fac * p->GetYDown(i), post,
		   fac * p->GetYUp(i)  , post);
	    std::cout << s;
	  }
	}
	// if (sum && p2p) std::cout << " ]";
	std::cout << std::endl;
      } // for i 
      gROOT->DecreaseDirLevel();
    } // if poi
    gROOT->DecreaseDirLevel();
  } //*MENU*
  /** 
   * Say that this should be shown as a folder
   * 
   * @return true
   */
  virtual Bool_t IsFolder() const { return true; }
  /** 
   * Browse this object 
   * 
   * @param b Browser to use 
   */
  virtual void Browse(TBrowser* b)
  {
    if (fMap)        b->Add(fMap);
    if (fQualifiers) b->Add(fQualifiers);
    b->Add(&fCommon);
    b->Add(&fPoint2Point);
    if (fData)  b->Add(fData);
    if (fDrawn) b->Add(fDrawn);
  }
  /* @} */
  /** 
   * @{ 
   * @name Drawing/Fitting 
   */
  /** 
   * Draw this data 
   *
   * Options: 
   *
   *   - STACK/COMBINED    Either errors are stacked or combined 
   *   - QUADRATIC/DIRECT  Add errors in quadrature or direct 
   *   - STAT              Add statistical errors to systematics 
   *   - COMMON            Add common errors to points
   *   - SPLIT             Without COMMON - do not stack common errors
   *   - MIN               Without COMMON - put common near minimum
   *   - MAX               Without COMMON - put common near maximum
   *   - WEST/EAST         Without COMMON - put common west/east
   *   - AXIS              Paint axis 
   *
   * If option COMMON isn't given and neither MIN nor MAX is not
   * given, then the common errors are displayed near the middle of
   * the Y range
   *
   * some examples are shown in the image below
   *
   * @image html DrawStyles.png
   *
   * @dontinclude tests/TestMaker.C 
   * @skip TestMaker
   * @until }
   *
   * A function to set-up an object
   *
   * @skip MakeGaus 
   * @until EndTest1
   *
   * A function to make a canvs 
   *
   * @dontinclude tests/DrawStyles.C 
   * @skip MakeCanvas
   * @until EndMakeCanvas
   *
   * Function to draw the stuff 
   *
   * @skip DrawIt
   * @until EndDrawIt
   *
   * Some utilies 
   *
   * @skip Sizes 
   * @until EndSizes
   *
   * Steering function to do all tests 
   *
   * @skip DrawAll
   * @until EndDrawAll
   *
   * The various ways we can draw the data 
   * 
   * First, combining all systematic errors 
   * 
   * @skip DrawCombined
   * @until EndDrawCombined
   * 
   * Then, stacking all systematics.
   * 
   * @skip DrawStack
   * @until EndDrawStack
   * 
   * We can also combine all errors 
   * 
   * @skip DrawCombinedCommonStat
   * @until EndDrawCombinedCommonStat
   * 
   * We can also stack the errors 
   * 
   * @skip DrawStackCommonStat
   * @until EndDrawStackCommonStat
   * 
   * First, combining all errors 
   * 
   * @skip DrawStackStat
   * @until EndDrawStackStat
   *
   * End of the tester class 
   *
   * @skip GraphSysErr
   * @until };
   *
   * The entry point for the script 
   *
   * @skip void
   * @until EndDrawStyles
   *
   */
  void Draw(Option_t* option="")
  {
    TString opt(option);
    opt.ToUpper();
    Bool_t   clear   = opt.Contains("CLEAR");
    Bool_t   axis    = opt.Contains("AXIS");
    // --- Optionally clear old stack --------------------------------
    if (clear && fDrawn) { 
      delete fDrawn;
      fDrawn = 0;
    }

    fDrawn = MakeMulti(option);
    if (!fDrawn) return;

    fDrawn->Draw(axis ? "A" : "");
    if (axis) { 
      fDrawn->GetHistogram()->SetXTitle(fXTitle);
      fDrawn->GetHistogram()->SetYTitle(fYTitle);
    }
  }  //*MENU*
  /** 
   * Fit a function to the data.  Which errors are considered depends
   * on the options given in drawOption i.e., 
   *
   * - STAT COMMON: Statistical, common, and point-to-point errors 
   *                are factored in 
   * - STAT:        Statistical and point-to-point errors are factored in, 
   *                but common errors are not 
   * - COMMON:      Common and point-to-point errors are factored in, 
   *                but statistical errors are not 
   * - otherwise:   Ppoint-to-point errors are factored in, 
   *                but statistical and common errors are not   
   * 
   * @param f1         Pointer to function objet 
   * @param fitOption  The fit options (See TGraph::Fit)
   * @param drawOption Draw options (See Draw)
   * @param min        Least X value to consider
   * @param max        Largest X value to consider
   * 
   * @return See TGraph::Fit
   */
  TFitResultPtr Fit(TF1* f1, Option_t* fitOption, Option_t* drawOption, 
		    Axis_t min=0, Axis_t max=0)
  {
    TString dOpt(drawOption);
    dOpt.ToUpper();
    Bool_t   clear   = dOpt.Contains("CLEAR");
    Bool_t   axis    = dOpt.Contains("AXIS");

    if (clear && fDrawn) { 
      delete fDrawn;
      fDrawn = 0;
    }

    fDrawn = MakeMulti(drawOption);
    if (!fDrawn || 
	!fDrawn->GetListOfGraphs() || 
	!fDrawn->GetListOfGraphs()->First()) 
      return TFitResultPtr(-1);
    // fDrawn->GetListOfGraphs()->ls();
    Graph* g = static_cast<Graph*>(fDrawn->GetListOfGraphs()->First());
    
    TString fOpt(fitOption);
    fOpt.ToUpper();
    Bool_t noStore = fOpt.Contains("N");
    Bool_t noDraw  = fOpt.Contains("0"); 
    Bool_t range   = fOpt.Contains("R"); 
    fOpt.Append("N0");
    TFitResultPtr r = g->Fit(f1, fOpt.Data(), "", min, max);

    if (!noStore) {
      Int_t status = r;
      if (status == 0) {
	if (min < max && !range) f1->SetRange(min, max);
	fDrawn->GetListOfFunctions()->Add(f1);
      }
    }

    if (!noStore && !noDraw) {
      fDrawn->Draw(axis ? "A" : "");
      if (axis) { 
	fDrawn->GetHistogram()->SetXTitle(fXTitle);
	fDrawn->GetHistogram()->SetYTitle(fYTitle);
      }
    }

    return r;
  }
  /** 
   * Fit a function to the data.  Which errors are considered depends
   * on the options given in drawOption i.e., 
   *
   * - STAT COMMON: Statistical, common, and point-to-point errors 
   *                are factored in 
   * - STAT:        Statistical and point-to-point errors are factored in, 
   *                but common errors are not 
   * - COMMON:      Common and point-to-point errors are factored in, 
   *                but statistical errors are not 
   * - otherwise:   Ppoint-to-point errors are factored in, 
   *                but statistical and common errors are not   
   * 
   * @param formula    The fit formula 
   * @param fitOption  The fit options (See TGraph::Fit)
   * @param drawOption Draw options (See Draw)
   * @param min        Least X value to consider
   * @param max        Largest X value to consider
   * 
   * @return See TGraph::Fit
   */
  TFitResultPtr Fit(const char* formula, 
		    Option_t* fitOption, Option_t* drawOption, 
		    Axis_t min=0, Axis_t max=0)
  {
    TString fname(formula);
    Bool_t  linear = fname.Contains("++");
    TF1*    f1     = 0;
    if (linear) f1 = new TF1(formula,formula,min,max);
    else { 
      f1 = static_cast<TF1*>(gROOT->GetFunction(formula));
      if (!f1) Warning("Fit", "Unknown function %s", formula);
    }
    if (!f1) return -1;
    
    return Fit(f1, fitOption, drawOption, min, max);
  }
  /** 
   * Get last drawn multigraph or create a new one 
   * 
   * @param option Options 
   * 
   * @return The multi graph
   */
  TMultiGraph* GetMulti(Option_t* option="") 
  {
    if (!fDrawn) fDrawn = MakeMulti(option);
    return fDrawn;
  }
  /* @} */
  /** 
   * @{ 
   * @name Various operations 
   */
  /** 
   * Scale graph by a constant 
   * 
   * @param f Function 
   * @param s Constant for non-relative common errors
   */
  void Scale(TF1* f, Double_t s=1)
  {
    HolderCommon* cmn = 0;
    TIter         cNext(&fCommon);
    while ((cmn = static_cast<HolderCommon*>(cNext()))) {
      if (cmn->IsRelative())
	// Do not scale relative errors - they scale by themselves
	continue;
      cmn->Set(cmn->GetYDown()/s, cmn->GetYUp()/s);
    }
    for (Int_t i = 0; i < GetN(); i++) {
      Double_t x = GetX(i);
      Double_t y = GetY(i);
      Double_t g = f->Eval(x);
      SetPoint(i, GetX(i), y/g);
      SetStatError(i, GetStatErrorDown(i)/g,GetStatErrorUp(i)/g);
      
      HolderP2P* p2p = 0;
      TIter      pNext(&fPoint2Point);
      while ((p2p = static_cast<HolderP2P*>(pNext()))) {
	Double_t fac = (p2p->IsRelative() ?
			(TMath::Abs(y) < 1e-9) ? 0 : 1/y : 1/g);
	SetSysError(p2p->GetUniqueID(), i,
		    p2p->GetXLeft(i), p2p->GetXRight(i),
		    fac*p2p->GetYDown(i), fac*p2p->GetYUp(i)); 		    
      }
    }
  }
  /** 
   * Scale graph by a constant 
   * 
   * @param s 
   */
  void Scale(Double_t s)
  {
    HolderCommon* cmn = 0;
    TIter         cNext(&fCommon);
    while ((cmn = static_cast<HolderCommon*>(cNext()))) {
      if (cmn->IsRelative())
	// Do not scale relative errors - they scale by themselves
	continue;
      cmn->Set(s * cmn->GetYDown(), s * cmn->GetYUp());
    }
    for (Int_t i = 0; i < GetN(); i++) {
      Double_t y = GetY(i);
      SetPoint(i, GetX(i), s*y);
      SetStatError(i, s*GetStatErrorDown(i),s*GetStatErrorUp(i));
      
      HolderP2P* p2p = 0;
      TIter      pNext(&fPoint2Point);
      while ((p2p = static_cast<HolderP2P*>(pNext()))) {
	Double_t fac = (p2p->IsRelative() ?
			(TMath::Abs(y) < 1e-9) ? 0 : 1/y :
			s);
	SetSysError(p2p->GetUniqueID(), i,
		    p2p->GetXLeft(i), p2p->GetXRight(i),
		    fac*p2p->GetYDown(i), fac*p2p->GetYUp(i)); 		    
      }
    }
  }
  /** 
   * Average one or more graphs. 
   *
   * The resulting graph is the weighted mean of the input graphs.
   * Common systematic errors common to all input graphs can be
   * specified before hand.  
   *
   * All other errors are taken in to account and are propagated to
   * point-to-point errors.
   * 
   * @param others  List of graphs
   * @param sep     If true, then try to separate out point-to-point errors. 
   */
  Int_t Average(const TCollection* others, Bool_t sep=true)
  {
    const Double_t tol = 1e-9;
    TList          sharedC; sharedC.SetOwner();
    TList          sharedP; 
    TIter          oNext(others);
    GraphSysErr*   other = 0;
    Double_t       xMin  = +1e9;
    Double_t       xMax  = -1e9;
    Int_t          nX    = 0;
    Bool_t         first = true;
    Bool_t         nonSh = false;
    Int_t          npS   = 0;
    // Here, we do some initial investigations into the the passed
    // graphs. We should find all common errors that are shared amoung
    // the input graphs, and we need to check for the consistency of
    // those.
    while ((other = static_cast<GraphSysErr*>(oNext()))) {
      Int_t n = other->GetN();
      nX      = TMath::Max(nX,n);
      xMin    = TMath::Min(xMin,other->GetX(0)  -2*other->GetErrorXLeft(0));
      xMax    = TMath::Max(xMax,other->GetX(n-1)+2*other->GetErrorXRight(n-1));

      if (first) { 
	// For the first graph, we simply add the title of all the
	// common systematics.
	HolderCommon* cmn = 0;
	TIter         cNext(&other->fCommon);
	while ((cmn = static_cast<HolderCommon*>(cNext()))) {
	  HolderCommon* o = static_cast<HolderCommon*>(cmn->Clone());
	  sharedC.Add(o);
	}
	if (sep) {
	  // and then the same for the point to point
	  HolderP2P* p2p = 0;
	  TIter      pNext(&other->fPoint2Point);
	  while ((p2p = static_cast<HolderP2P*>(pNext()))) {
	    sharedP.Add(p2p->Clone());
	  }
	}
	
	// We also copy the keys
	CopyKeys(other, "ar");
      }
      else {
	// Otherwise, we check that this graph has all the common
	// errors in our temporary list, and if not, we remove that
	// error from the list
	TObjLink* cur  = sharedC.FirstLink();
	while (cur) {
	  HolderCommon* o  = static_cast<HolderCommon*>(cur->GetObject());
	  HolderCommon* c  = other->FindCompat(o, tol, true);
	  if (!c) {
	    // This error not found in the graph we're looking at.
	    // Remove it from the list
	    Info("Average", "Common systematic %s not found in %s",
		 o->GetTitle(), other->GetName());
	    other->Print("sys");
	    TObjLink* keep = cur->Next();
	    TObject*  obj  = sharedC.Remove(cur);
	    cur            = keep;
	    nonSh          = true;
	    if (obj) delete obj;
	    continue;
	  }
	  cur = cur->Next();
	} // while cur

	// Now check the point-to-point errors 
	cur = sharedP.FirstLink();
	while (cur) {
	  HolderP2P* o  = static_cast<HolderP2P*>(cur->GetObject());
	  HolderP2P* p  = other->FindCompat(o);
	  if (!p) {
	    // This error not found in the graph we're looking at.
	    // Remove it from the list
	    Info("Average", "P2P systematic %s not found in %s",
		 o->GetTitle(), other->GetName());
	    TObjLink* keep = cur->Next();
	    TObject*  obj  = sharedP.Remove(cur);
	    cur            = keep;
	    npS++;
	    if (obj) delete obj;
	    continue;
	  }
	  cur = cur->Next();
	} // while cur	  
      } // else
      first = false;
    }   // while cmn
    // Now we have the list of common errors that are shared amoung
    // all the input graphs.  So we define those new common errors on
    // this graph, and disable them on the inputs.
    TIter          nextC(&sharedC);
    HolderCommon*  com = 0;
    while ((com = static_cast<HolderCommon*>(nextC()))) {
      Int_t cId = DefineCommon(com->GetTitle(), com->IsRelative(),
			       com->GetYDown(1), com->GetYUp(1), kBox);
      HolderCommon* c = FindCommon(cId);
      c->CopyAttr(com);
      // Info("Average", "Shared common error %s +%7f%s -%7f%s",
      // 	   c->GetTitle(),
      // 	   c->GetYDown(1), (c->IsRelative() ? "%" : ""),
      // 	   c->GetYUp(1),   (c->IsRelative() ? "%" : ""));

      // Now loop over all input graphs, find the shared commons, and
      // mark them as used.
      oNext.Reset();
      while ((other = static_cast<GraphSysErr*>(oNext()))) {
	HolderCommon* oc = other->FindCompat(c, tol, true);
	if (!oc) {
	  Error("Average", "This should NOT happen");
	  other->Print("sys");
	  continue;
	}
	oc->SetBit(kUsedBit);
      } // while other
    } // while com
    sharedC.Clear();
    // If we found some common errors that are not shared, we need to
    // make a point-to-point error for absorbing the remaining common
    // errors.
    Int_t bId = 0;
    if (nonSh || !sep) {
      bId = DeclarePoint2Point("Blob", false, kBox);
      SetSysFillColor(bId, kRed+2);
      SetSysLineColor(bId, kRed+2);
      SetSysFillStyle(bId, 3004);
    }
    
    // At this point, we have found all the common shared errors and
    // added those to ourselves. Those that are not shared go into the
    // blob point-to-point error.
    // 
    // Next step, is to deal with the shared point-to-point errors.
    // If we have shared point-to-point errors, we need to deal with
    // them separately.  To do that, we make a list of identifiers
    // that we should inspect at each X value.  The identifiers point
    // to the p2p errors of each input graph.
    //
    TArrayL64     mp2p;
    Int_t         dOf = 6;
    if (sep) {
      TIter       nextP(&sharedP);
      HolderP2P*  p2p = 0;
      Int_t       j   = 0;
      mp2p.Set(100);
      while ((p2p = static_cast<HolderP2P*>(nextP()))) {
	// Info("Average", "Loop %d shared P2P %s", j,  p2p->GetTitle());
	Int_t      pId = DeclarePoint2Point(p2p->GetTitle(),
					  p2p->IsRelative(), kBox);
	HolderP2P* p   = FindP2P(pId);
	p->CopyAttr(p2p);
	Info("Average", "Shared point-to-point error %s", p->GetTitle());
	
	// Now loop over all input graphs, find the shared commons, and
	// mark them as used.
	oNext.Reset();
	Long64_t mask = 0;
	Int_t    off  = 0;
	while ((other = static_cast<GraphSysErr*>(oNext())) && off < 64) {
	  HolderP2P* op = other->FindCompat(p, true);
	  if (!op) {
	    Error("Average", "This should NOT happen");
	    continue;
	  }
	  mask |= ((op->GetUniqueID() & 0x3F) << off);
	  off  += dOf;
	  // Info("Average", "Mark %s as used", op->GetTitle());
	  op->SetBit(kUsedBit);	
	} // while other
	if (off >= 64) {
	  Warning("Average",
		  "Some shared point-to-point errors could not be encoded "
		  "becasue we have too many (%d>%d) input graphs",
		  others->GetEntries(), 64/dOf);
	}
	mp2p[pId] = mask;
	j++;
      } // while p2p
      // Info("Average", "Declared new shared P2P - mapping:");
      // for (Int_t j = 1; j < mp2p.GetSize(); j++) {
      //   if (mp2p[j] == 0) continue;
      //   rintf("%3d/%3d", j, mp2p.GetSize());
      //   for (Int_t k = 0; k < 64; k += dOf)
      //     printf(" %3d", (mp2p[j] >> k) & 0x3f);
      //   printf("\n");
      // }
      
      // sharedP.Clear();
      // Now we have all shared point-to-point errors defined, and we
      // have the ids for each graph in the shar array (6bits each).
      //
      // Next, we will declare _all_ other our point-to-point errors
      // separately on this graph.  That means that the point-to-point
      // errors can be 0 in some ranges.
      oNext.Reset();
      Int_t       off  = 0;
      while ((other = static_cast<GraphSysErr*>(oNext()))) {
	Long64_t*   mask = 0;
	TIter       oNextP(&(other->fPoint2Point));
	HolderP2P*  p2p = 0;
	while ((p2p = static_cast<HolderP2P*>(oNextP())) && off < 64) {
	  if (p2p->TestBit(kUsedBit)) {
	    // Info("Average", "%s marked as used", p2p->GetTitle());	
	    continue; // Already in the shared part
	  }
	  
	  // The point-to-point error should go into the weight, but not
	  // in the result.
	  p2p->SetBit(kOnlyWeightBit);
	  Int_t      pId = 0;
	  HolderP2P* p   = FindCompat(p2p);
	  if (!p) {
	    pId = DeclarePoint2Point(p2p->GetTitle(), p2p->IsRelative(), kBox);
	    p   = FindP2P(pId);
	    p->CopyAttr(p2p);
	    mask =  &(mp2p[pId]);

	  }
	  *mask |= ((p2p->GetUniqueID() & 0x3F) << off);
	}
	off += dOf;
      } // while other
      if (off >= 64) {
	Warning("Average",
		"Some shared point-to-point errors could not be encoded "
		"becasue we have too many (%d>%d) input graphs",
		others->GetEntries(), 64/dOf);
      }
      // Info("Average", "Declared new shared P2P - mapping:");
      // for (Int_t j = 0; j < mp2p.GetSize(); j++) {
      // 	if (mp2p[j] == 0) continue;
      // 	printf("%3d/%3d", j, mp2p.GetSize());
      // 	for (Int_t k = 0; k < 64; k += dOf)
      // 	  printf(" %3d", (mp2p[j] >> k) & 0x3f);
      // 	printf("\n");
      // }
    }
    
    // Now search for the X values to evaluate at 
    Double_t x    = xMin;
    Double_t oldX = xMin;
    Int_t    cnt  = 0;
    TArrayD  xa(2*nX); // Assume no more than 100 points
    TArrayD  exla(2*nX);
    TArrayD  exha(2*nX);
    do {
      // Find the next X to evaluate at
      oNext.Reset();
      Double_t nextX = xMax;
      Double_t texl  = 0;
      Double_t texh  = 0;
      // Info("Average", "nextX=%f xMax=%f", nextX, xMax);
      while ((other = static_cast<GraphSysErr*>(oNext()))) {
	for (Int_t i = other->GetUniqueID(); i < other->GetN(); i++) {
	  if (other->GetX(i) <= x) {
	    // Info("Average", "- x_i=%f < x=%f", g->GetX(i), x);
	    continue;
	  }
	  if (other->GetX(i) < nextX) {
	    // Info("Average", "Set x to %f (was %f)", g->GetX(i), nextX);
	    nextX = other->GetX(i);
	    texl  = TMath::Max(texl, other->GetErrorXLeft(i));
	    texh  = TMath::Max(texh, other->GetErrorXRight(i));
	  }
	  else {
	    // g->SetUniqueID(i);
	    // Info("Average", "No more (%s) here at x=%f (next=%f)",
	    //      g->GetName(), g->GetX(i), nextX);
	    break;
	  }
	}
      }
      if (nextX >= xMax) {
	// Info("Average", "Next x is too large %f>%f", nextX, xMax);
	break;
      }

      if (cnt >= xa.GetSize()) {
	Warning("Average", "increasing size of X cache");
	xa.Set(1.5*cnt);
	exla.Set(1.5*cnt);
	exha.Set(1.5*cnt);
      }
      // Info("Average", "Set next x to %f +/- (%f,%f)", nextX, texl, texh);
      xa[cnt]     = nextX;
      exla[cnt]   = texl;
      exha[cnt]   = texh;
      Double_t dx = texh*2;//(nextX-oldX);
      oldX        = x;
      x           = nextX+dx/4;
      cnt++;
    } while (cnt < 1000); // Do not go for more than 1000 points
    // At this point we have
    //
    //  - All common errors defined, and we have set the shared ones
    //  - All point-to-point errors defined, and we know which input
    //    graph contributes to each.
    //  - A list of X values to evaluate each graph at. 
        
    // Print("sys");
    // Now we loop over the list of X values and evaluate each graph
    // at that X value.
    for (Int_t i = 0; i < cnt; i++) {
      Double_t xi   = xa[i];
      Double_t texl = exla[i];
      Double_t texh = exha[i];
      Double_t nexl = TMath::Abs(xi - (i     == 0   ? xa[i+1] : xa[i-1]))/2;
      Double_t nexh = TMath::Abs(xi - (i + 1 == cnt ? xa[i-1] : xa[i+1]))/2;
      Double_t exl  = TMath::Min(texl, nexl);
      Double_t exh  = TMath::Min(texh, nexh);
      // Info("Average", "@ %3d x=%f, exl=min(%f,%f)=%f exh=min(%f,%f)=%f",
      //      i, xi, texl, nexl, exl, texh, nexh, exh);

      
      // Now loop again - this time calculating our data point
      Double_t sum   = 0;
      Double_t sumW  = 0;
      // Double_t sumWS = 0;
      Double_t sumE  = 0;
      Double_t sumES = 0;
      oNext.Reset();

      // Use our combiner framework for calculating the averaged
      // result and errors.  Note, we define two combiners.  One that
      // include all non-shared common and point-to-point errors
      // (centroid), and one that doesn't have the non-shared common
      // and point-to-point errors as well as the split point-to-point
      // errors (error).  The first one is used to set the centroid of
      // the data, while the second is used to set the combined error
      LinearSigmaCombiner centroid;
      LinearSigmaCombiner error;
      LinearSigmaCombiner stats;
      Int_t               j = 0;
      TArrayI             ai1(others->GetEntries());
      TArrayI             ai2(others->GetEntries());
      TArrayD             ax(others->GetEntries());
      while ((other = static_cast<GraphSysErr*>(oNext()))) {
	Bool_t cmn    = true;   // Include the common errors
	Bool_t stat   = false;  // Include the statistical errors
	Bool_t quad   = true;   // Add in quadrature
	Bool_t nosqrt = false;  // Do not return squared error 

	Double_t y, eyl, eyh, wyl, wyh, syl, syh;

	Int_t  fret   = other->FindPoint(xi, ai1[j], ai2[j]);
	if (fret < -1) { j++; continue; } // Nothing here
	if (fret >= 0) {
	  y = other->GetYandError(fret,cmn,stat,quad,nosqrt,eyl, eyh, wyl, wyh);
	  syl = other->GetStatErrorDown(fret);
	  syh = other->GetStatErrorUp(fret);
	}
	else {
	  Double_t eyl1, eyh1, wyl1, wyh1;
	  Double_t eyl2, eyh2, wyl2, wyh2;
	  Double_t x1  = other->GetX(ai1[j]);
	  Double_t x2  = other->GetX(ai2[j]);
	  Double_t y1  = other->GetYandError(ai1[j],cmn,stat,quad,nosqrt,
					     eyl1,eyh1,wyl1,wyh1);
	  Double_t y2  = other->GetYandError(ai2[j],cmn,stat,quad,nosqrt,
					     eyl2,eyh2,wyl2,wyh2);
	  Double_t syl1 = other->GetStatErrorDown(ai1[j]);
	  Double_t syh1 = other->GetStatErrorUp(ai1[j]);
	  Double_t syl2 = other->GetStatErrorDown(ai2[j]);
	  Double_t syh2 = other->GetStatErrorUp(ai2[j]);
	  Double_t dx = (x2-x1);
	  ax[j]       = (xi-x1)/dx;
	  y           = y1   + ax[j] * (y2-y1);
	  eyl         = eyl1 + ax[j] * (eyh1-eyl1);
	  eyh         = eyh2 + ax[j] * (eyh2-eyl2);
	  wyl         = wyl1 + ax[j] * (wyh1-wyl1);
	  wyh         = wyh2 + ax[j] * (wyh2-wyl2);
	  syl         = syl1 + ax[j] * (syh1-syl1);
	  syh         = syh2 + ax[j] * (syh2-syl2);
	}
	
	centroid.Add(y, wyl, wyh);
	error.Add(y, eyl, eyh);
	stats.Add(y, syl, syh);
	j++;
      } // while (other)
      centroid.Calculate();
      // Info("Average", "Centroid");
      // centroid.Print();
      // centroid.fResult->Print();
      Double_t newY = centroid.fResult->fX;
      SetPoint(i, xi, newY);
      SetPointError(i, exl, exh);

      stats.Calculate();
      // Info("Average", "Stats");
      // stats.Print();
      // stats.fResult->Print();
      stats.Calculate();
      Double_t fac = IsStatRelative() ? 1/newY : 1;
      SetStatError(i, fac*stats.fResult->fEl,fac*stats.fResult->fEh);

      if (bId > 0) {
      // Calculations for error
	error.Calculate();
	// Info("Average", "Error");
	// error.Print();
	// error.fResult->Print();

	fac = IsRelative(bId) ? 1/newY : 1;
	SetSysError(bId,i,0,0,fac*error.fResult->fEl,fac*error.fResult->fEh);
      }

      // Now, we have set the blob point-to-point error.  Next, we
      // need to do the point-to-point erros that should be transfered
      // to the output.  We have stored masks of IDs in mp2p, so we
      // loop over that array first.  The lower 6 bits is the output
      // Id, while the remaining bits (counting from the low side) are
      // ids for each input graph.
      for (Int_t j = 1; j < mp2p.GetSize(); j++) {
	Long64_t mask = mp2p[j];
	if (mask == 0) continue;
	Int_t    tid  = j; // (mask & 0x3f);
	HolderP2P* tp2p = FindP2P(tid);
	if (!tp2p) {
	  Warning("Average", "No target point-to-point error at %d for id=%d",
		  j, tid);
	  continue;
	}
	// Zero by default
	SetSysError(tid, i, 0, 0);
	Long64_t rem  = mask;
	if (rem == 0) 
	  // There's no IDs for this, explicitly set the value to 0
	  continue;

	// Now we should loop over all the input graphs and get the
	// error and finally combine them into one.  Again, we use a
	// linear sigma approximation from our combiner framework to
	// do the calculations.
	LinearSigmaCombiner sc;
	Int_t k = 0; // Index into index array and for off set in mask
	oNext.Reset();
	while ((other = static_cast<GraphSysErr*>(oNext()))) {
	  if (ai1[k] == -1 && ai2[k] == -1) {
	    // This graph does not contribute at this point
	    // Info("Average", "No contribution at %f (%d:%d,%d) for %s to %s",
	    //      xi, k, ai1[k], ai2[k], other->GetName(), tp2p->GetTitle());
	    k++;
	    continue;
	  }
	  Int_t sid = ((rem >> dOf*k) & 0x3f);
	  if (sid <= 0) {
	    // This graph does not contribute to this error
	    // Info("Average",
	    //      "No contribution for %s to %s, not found [%d:%d,%d]",
	    //      other->GetName(), tp2p->GetTitle(), k, ai1[k], ai2[k]);
	    k++;
	    continue;
	  }
	  HolderP2P* sp2p = other->FindP2P(sid);
	  if (!sp2p) {
	    Warning("Average", "Couldn't find point-to-point error %d (%s)"
		    "in graph %d (%s) - should not happen",
		    sid, tp2p->GetTitle(), k, other->GetName());
	    k++;
	  }
	  Double_t sy   = 0;
	  Double_t seyl = 0;
	  Double_t seyh = 0;
	  if (ai1[k] == ai2[k]) {
	    // Exact point
	    sy   = other->GetY(ai1[k]);
	    seyl = other->GetSysErrorYDown(sid, ai1[k]);
	    seyh = other->GetSysErrorYUp(sid, ai1[k]);
	  }
	  else {
	    // Interpolate
	    Double_t sy1   = other->GetY(ai1[k]);
	    Double_t sy2   = other->GetY(ai2[k]);
	    sy             = sy1 + ax[k] * (sy2-sy1);
	    seyl           = sp2p->GetYDown(ai1[k],ai2[k],ax[k]);
	    seyh           = sp2p->GetYUp(ai1[k],ai2[k],ax[k]);
	  }
	  //Info("Average", "Contribution (%d:%d-%d:%f,%f) to %s from %s",
	  //     k,ai1[k],ai2[k],seyl,seyh,tp2p->GetTitle(),other->GetName());
	  sc.Add(sy, seyl, seyh);
	  k++;
	} // while other
	sc.Calculate();
	// Info("Average", "%s", tp2p->GetTitle());
	// sc.Print();
	// sc.fResult->Print();
	fac = IsRelative(tid) ? 1/newY : 1;
	SetSysError(tid, i, 0, 0, fac*sc.fResult->fEl, fac*sc.fResult->fEh);
      } // for j (# of p2p)
    } // for i (points)

    oNext.Reset();
    while ((other = static_cast<GraphSysErr*>(oNext()))) 
      other->ClearUsed();
    ClearUsed();
      
    return bId;
  }
  /** 
   * Calculate the intergral and error on integral of the graph. 
   * 
   * @param error  On return, the error on the integral 
   * @param option What to take into account (see Draw)
   * @param first  First point to evalute at 
   * @param last   Last point to evaluate at 
   *
   * @return The integral 
   */
  Double_t Integral(Double_t& error,
		    Option_t* option="quad sum stat",
		    UShort_t  first=0,
		    Short_t   last=-1)
  {
    if (last < 0) last = GetN()-1;
    Double_t sum  = 0;
    Double_t err2 = 0;
    TString  opts(option); opts.ToLower();
    Bool_t   cmn  = opts.Contains("comm");
    Bool_t   stat = opts.Contains("stat");
    Bool_t   quad = opts.Contains("quad");
    
    for (Short_t i = first; i <= last; i++) {
      Double_t eyl, eyh;
      Double_t y   = GetYandError(i, cmn, stat, quad, false, eyl, eyh);
      Double_t eym = (eyh+eyl)/2;
      Double_t x   = GetX(i);
      Double_t exl = GetErrorXLeft(i);
      Double_t exr = GetErrorXRight(i);
      Double_t x1  = x - exl;
      Double_t x2  = x + exr;
      Double_t prv = x1;
      Double_t nxt = x2;
      if (i-1 >= 0)     prv = GetX(i-1)+GetErrorXRight(i-1);
      if (i+1 < GetN()) nxt = GetX(i+1)-GetErrorXLeft(i+1);
      if (prv > x1) x1 -= (prv-x1)/2; // Adjust to half-ways 
      if (nxt < x2) x2 -= (x2-nxt)/2; // Adjust to half-ways
      Double_t wdt = (x2-x1);
      if (wdt < 0) continue; // Do not take this point
      sum  += wdt * y;
      err2 += TMath::Power(eym*wdt,2);
    }
    error = TMath::Sqrt(err2);
    return sum;
  }
  /** 
   * Find the next point and for ratio 
   * 
   * @param i     Point numver 
   * @param num   Numerator 
   * @param denom Denominator 
   * @param x     On return, the X value 
   * @param dY    On return, the denominator Y
   * @param dEyl  On return, the denominator lower error on Y
   * @param dEyh  On return, the denominator upper error on Y
   * @param nY    On return, the numerator Y		     
   * @param nEyl  On return, the numerator lower error on Y
   * @param nEyh  On return, the numerator upper error on Y
   * 
   * @return 
   */
  static Bool_t NextPoint(Int_t i,
			  const GraphSysErr* num,
			  const GraphSysErr* denom,
			  Double_t& x,
			  Double_t& dY,
			  Double_t& dEyl,
			  Double_t& dEyh,
			  Double_t& nY,
			  Double_t& nEyl,
			  Double_t& nEyh)
  {
    x  = num->GetX(i);
    if (!denom->FindYandError(x, false, true, true, false, dY, dEyl, dEyh))
      return false; 
    if (TMath::Abs(dY) < 1e-9)
      return false;
    
    nY   = num->GetYandError(i,false,true,true,false,nEyl, nEyh);
    return true;
  }
  
  /** 
   * Take the ratio of 2 graphs.  
   * 
   * @param num    Numerator 
   * @param den    Denominator 
   * @param flags Flags for the ratio calculation.  See also
   * ERatioOption.
   * @param fac    Factor for searching points
   * 
   * @return Newly allocated graph 
   */
  static GraphSysErr* Ratio(const GraphSysErr* num,
			    const GraphSysErr* den,
			    UInt_t             flags=kRatioDefault,
			    Double_t           fac=1)
  {
    // num->Print("sys");
    // den->Print("sys");
    num->ClearUsed();
    den->ClearUsed();
    GraphSysErr* ret = new GraphSysErr(1);
    ret->CopyAttr(num);
    ret->CopyKeys(num, "ar");
    Bool_t   emax    = (flags & kMax);
    Bool_t   cmin    = (flags & kMin);
    Bool_t   cmax    = (flags & kMax);
    Bool_t   denrel  = (flags & kDenomRel);
    Bool_t   cancc   = (flags & kCancel);
    Int_t    id      = 0;
    Int_t    notUsed = 0;
    // Loop over all common errors in the numerator, and see if we
    // have the same common error in the denominator.  If we do, we
    // can (partially) cancel this common error.
    //
    // Technically, we mark each used common error as used, and create
    // a new common error in the output graph.
    //
    // Common errors that are not shared by both the numerator and the
    // denominator will not be cancelled.  Instead, if they are
    // relative, they will be added to the output graph directly.
    // Finally, common systematic errors that are not cancelled or are
    // not relative, will go into the blob of point-to-point errors
    TIter nextNC(&num->fCommon);
    HolderCommon* comN = 0;
    while ((comN = static_cast<HolderCommon*>(nextNC()))) {      
      Bool_t        rel  = comN->IsRelative();
      Double_t      eyl  = comN->GetYDown(1);
      Double_t      eyh  = comN->GetYUp(1);
      Int_t         idD  = den->FindId(comN->GetTitle());      
      HolderCommon* comD = 0;
      if (idD != 0 &&                       // Check it exists
	  (comD = den->FindCommon(idD)) &&  // Check it is a common
	  (rel == comD->IsRelative())) {    // Same as this
	// We found the same error in the denominator, so we need to
	// cancel errors here.
	Double_t dEyl = comD->GetYDown(1);
	Double_t dEyh = comD->GetYUp(1);
	Double_t nEyl = comN->GetYDown(1);
	Double_t nEyh = comN->GetYUp(1);
	if (cmin) {
	  eyl = TMath::Min(nEyl, dEyl);
	  eyh = TMath::Min(nEyh, dEyh);
	}
	else if (cmax) {
	  eyl = TMath::Max(nEyl, dEyl);
	  eyh = TMath::Max(nEyh, dEyh);
	}
	else if (cancc) { 
	  eyl  = TMath::Sqrt(TMath::Abs(dEyl*dEyl-nEyl*nEyl));
	  eyh  = TMath::Sqrt(TMath::Abs(dEyh*dEyh-nEyh*nEyh));
	}
	const char* post = (rel ? "%" : "");
	Double_t    pfc  = (rel ? 100 : 1);
			    
	if (kVerbose & kRatio)
	  ::Info("Ratio", "Cancelling the common systematic error %s "
		 "between numerator (-%f%s +%f%s) "
		 "and denominator (-%f%s +%f%s): -%f%s + %f%s",
		 comD->GetTitle(),
		 pfc*nEyl, post, pfc*nEyh, post,
		 pfc*dEyl, post, pfc*dEyh, post,
		 pfc*eyl,  post, pfc*eyh,  post);
	       
	// Flag the denominator error as used
	comD->SetBit(kUsedBit);
      }
      else if (!rel) {
	// We cannot cancel, and it is not relative, add to blob
	if (kVerbose & kRatio)
	  ::Info("Ratio",
		 "Common numerator systematic %s will be added to blob",
		 comN->GetTitle());
	notUsed++; // Count how many we're not setting directly
	continue;
      }
      // Now we can define our common error. Which value we set depend
      // on whether the error was cancelled with the denominator
      Int_t        cId  = ret->DefineCommon(comN->GetTitle(),rel,
					    eyl, eyh, kBox);
      HolderCommon* c   = ret->FindCommon(cId);
      c->CopyAttr(comN);
      // Flag the numerator error as used
      comN->SetBit(kUsedBit);
    }
    // Loop over all denominator common systematic errors.  If it is
    // not cancelled with the numerator and is relative, we define a
    // new error on the output ratio. Otherwise it will be absorbed in a blob 
    TIter nextDC(&den->fCommon);
    HolderCommon* comD = 0;
    while ((comD = static_cast<HolderCommon*>(nextDC()))) {
      if (comD->TestBit(kUsedBit)) continue;
      if (!comD->IsRelative()) {
	if (kVerbose & kRatio)
	  ::Info("Ratio",
		 "Common denominator systematic %s will be added to blob",
		 comD->GetTitle());
	notUsed++; // Count how many we're not setting directly
	continue;
      }
      if (kVerbose & kRatio)
	::Info("Ratio", "Propagating common sysmatic %s to ratio",
	       comD->GetTitle());
      // ret->Print("sys");
      // Now we can define our common error. Which value we set depend
      // on whether the error was cancelled with the denominator
      Int_t        cId  = ret->DefineCommon(comD->GetTitle(),true,
					    comD->GetYDown(1),
					    comD->GetYUp(1),
					    kBox);
      HolderCommon* c   = ret->FindCommon(cId);
      c->CopyAttr(comD);
      // Flag the denominator error as used
      comD->SetBit(kUsedBit);
      // ret->Print("sys");
    }
    
    // Loop over all point-to-point errors in the numerator, and see
    // if we have the same point-to-point error in the denominator.
    // If we do, we can (partially) cancel this point-to-point error.
    //
    // Technically, we mark each used point-to-point error as used.
    //
    // Point-to-point errors that are not shared by both the numerator
    // and the denominator will not be cancelled.  Instead, they will
    // be added in on the residual point-to-point sum (in quadrature)
    // error.
    //
    // We encode the IDs of the
    // errors in a single mask.
    //
    // 24 -------- 16 -------- 8 -------- 0
    //  | Denom ID  | Num ID   | Ratio ID |
    //  +-----------+----------+----------+
    //
    Int_t      sCnt = 0;
    TArrayI    shared(num->fPoint2Point.GetEntries() +
		      den->fPoint2Point.GetEntries());
    TIter      nextNP(&num->fPoint2Point);
    HolderP2P* p2pN = 0;
    while ((p2pN = static_cast<HolderP2P*>(nextNP()))) {      
      Bool_t     rel  = p2pN->IsRelative();
      Int_t      idD  = den->FindId(p2pN->GetTitle());      
      HolderP2P* p2pD = idD == 0 ? 0     : den->FindP2P(idD);
      Bool_t     same = !p2pD    ? false : rel == p2pD->IsRelative();
      Int_t      mask = 0;
      // Printf("Id of denom %s -> %d (%p, %s, %s)", p2pN->GetTitle(), idD,
      // 	     p2pD, same ? "same" : "diff",
      // 	     p2pD ? (p2pD->TestBit(kUsedBit) ? "used" : "ready") : "?");
      if (cancc    &&
	  p2pD     &&   // Check it is a p2p
	  same) {   // Same as this
	// We found the same error in the denominator, so we need to
	// cancel errors when looping.  
	mask = ((p2pD->GetUniqueID() & 0xFF) << 16);
	// Printf("Found %s in both denom (%s,%d) and num (%s,%d)",
	//        p2pN->GetTitle(),
	//        p2pD->GetTitle(), p2pD->GetUniqueID(),
	//        p2pN->GetTitle(), p2pD->GetUniqueID());
	if (kVerbose & kRatio)
	  ::Info("Ratio", "Cancelling the p2p systamtic error %s "
		 "between numerator and denominator",
		 p2pN->GetTitle());
	// Flag the denominator error as used
	p2pD->SetBit(kUsedBit);
      }
      Int_t      pId  = ret->DeclarePoint2Point(p2pN->GetTitle(),rel,
						kBox);
      HolderP2P* p    = ret->FindP2P(pId);
      p->CopyAttr(p2pN);
      // Flag the numerator error as used
      p2pN->SetBit(kUsedBit);
      mask |= (((p2pN->GetUniqueID() & 0xFF) << 8) |
	       (pId & 0xFF));
      // Printf("Found shared P2P: %s -> %d",  p2pN->GetTitle(), pId);
      // num->Print("sys");
      // den->Print("sys");
      shared[sCnt] = mask;
      sCnt++;
    }
    // Loop over all denominator point-to-point systematic errors.  If
    // it is not cancelled with the numerator, we define a new error
    // on the output ratio
    TIter      nextDP(&den->fPoint2Point);
    HolderP2P* p2pD = 0;
    while ((p2pD = static_cast<HolderP2P*>(nextDP()))) {      
      if (p2pD->TestBit(kUsedBit)) {
	// Printf("P2P of denom %s already used", p2pD->GetTitle());
	continue;
      }
      Int_t      pId  = ret->DeclarePoint2Point(p2pD->GetTitle(),true,
						kBox);
      HolderP2P* p    = ret->FindP2P(pId);
      p->CopyAttr(p2pD);
      // Flag the numerator error as used
      p2pD->SetBit(kUsedBit);
      Int_t mask = (((p2pD->GetUniqueID() & 0xFF) << 16) |
		    (pId & 0xFF));
      // Printf("P2P of %s (%d) denom propagated", p2pD->GetTitle(),
      //        p2pD->GetUniqueID());
      shared[sCnt] = mask;
      sCnt++;
    }    
    Int_t bId = 0;
    if (notUsed > 0) {
      // If we have some error that are not used - i.e., non-relative
      // common errors in the numerator or denominator, we need to
      // make a new point-to-point systematic error for those.
      bId = ret->DeclarePoint2Point("Blob", false, kBox);
    }

    // Now loop over all the points we have in the numerator and
    // figure out what to do for each point.
    Int_t cnt = 0; // Counter of output points
    for (Int_t i = 0; i < num->GetN(); i++) {
      Double_t x    = num->GetX(i);
      Double_t nY   = num->GetY(i);
      Double_t nEyl = 0;
      Double_t nEyh = 0;
      Double_t dY   = 0;
      Double_t dEyl = 0;
      Double_t dEyh = 0;
      Int_t    di1, di2;
      Int_t    di   = den->FindPoint(x, di1, di2, fac);
      Double_t daX  = 0; // Relative distance between interpolated points 
      if (di < -1) {
	if (kVerbose & kRatio)
	  ::Warning("Ratio", "Next point %d (%f) not found in denominator",
		    i, x);
	continue;
      }
      Bool_t cmn    = true;  // Include the common errors
      Bool_t stat   = false;  // Include the statistical errors
      Bool_t quad   = true;  // Add in quadrature
      Bool_t nosqrt = false; // Do not return squared error 
      if (di >= 0) {
	// Exact point
	// point,common,stat,quadratic,nosqrt,
	dY = den->GetYandError(di,cmn,stat,quad,nosqrt,dEyl,dEyh);
	di1 = di;
      }
      else {
	Double_t eyl1, eyl2, eyh1, eyh2;
	Double_t x1 = den->GetX(di1);
	Double_t x2 = den->GetX(di2);
	Double_t y1 = den->GetYandError(di1,cmn,stat,quad,nosqrt,eyl1,eyh1);
	Double_t y2 = den->GetYandError(di2,cmn,stat,quad,nosqrt,eyl2,eyh2);
	// Linear interpolation
	Double_t dx = (x2-x1);
	daX         = (x-x1)/dx;
	dY          = y1   + daX * (y2 - y1);
	dEyl        = eyl1 + daX * (eyl2 - eyl1);
	dEyh        = eyh1 + daX * (eyh2 - eyh1);
      }
      // At this point, we have the X and Y values, and the non-p2p
      // and non-cancelled and non-relative-common errors calculated
      // for numerator and denominator individually.
      //
      // We start by setting the point value.  Perhaps here, we should
      // also set the statistics.
      Double_t rY    = (dY != 0 ? nY/dY   : 0);
      ret->SetPoint(cnt, x, rY);
      ret->SetPointError(cnt, num->GetErrorXLeft(i), num->GetErrorXRight(i));

      // We have not dealt with the statistical errors yet.  We do
      // that here.
      Double_t sEyl = 0;
      Double_t sEyh = 0;
      if (TMath::Abs(rY) > 1.e-9) {
	Double_t snEyl = num->GetStatErrorDown(cnt);
	Double_t snEyh = num->GetStatErrorUp  (cnt);
	Double_t sdEyl = den->GetStatErrorDown(cnt);
	Double_t sdEyh = den->GetStatErrorUp  (cnt);
	// Statistical errors never cancel, and must be propagted. Note,
	// despite it technically being wrong, we treat the lower and
	// higher errors independently.
	sEyl = TMath::Sqrt(rY*rY*(snEyl*snEyl/nY/nY+sdEyl*sdEyl/dY/dY));
	sEyh = TMath::Sqrt(rY*rY*(snEyh*snEyh/nY/nY+sdEyh*sdEyh/dY/dY));
      }
      ret->SetStatError(cnt, sEyl, sEyh);
      
      // Next, we should update our blob systematic p2p error if it
      // exists.
      if (bId > 0) { 
	Double_t eyl = 0, eyh = 0;
	// Printf("Calculate blob syst.uncer. d=(%f,%f) n=(%f,%f)",
	//        dEyl, dEyh, nEyl, nEyh);
	if (denrel) {
	  Double_t rdEyl = (dY != 0 ? dEyl/dY : 0);
	  Double_t rdEyh = (dY != 0 ? dEyh/dY : 0);
	  eyl            = rdEyl*rY;
	  eyh            = rdEyh*rY;
	}
	else if (emax) {
	  // Take the largest of the two errors, and divide by the
	  // denominator. 
	  eyl = TMath::Max(nEyl, dEyl);
	  eyh = TMath::Max(nEyh, dEyh);
	  eyl          /= dY;
	  eyh          /= dY;
	}
	else if (cmin) {
	  // Take the largest of the two errors, and divide by the
	  // denominator. 
	  eyl = TMath::Min(nEyl, dEyl);
	  eyh = TMath::Min(nEyh, dEyh);
	  eyl          /= dY;
	  eyh          /= dY;
	}
	else {
	  // Find the relative errors and add in quadrature 
	  Double_t rnEyl = (nY != 0 ? nEyl/nY : 0);
	  Double_t rnEyh = (nY != 0 ? nEyh/nY : 0);
	  Double_t rdEyl = (dY != 0 ? dEyl/dY : 0);
	  Double_t rdEyh = (dY != 0 ? dEyh/dY : 0);
	  Double_t reyl  = TMath::Sqrt(rnEyl*rnEyl+rdEyl*rdEyl);
	  Double_t reyh  = TMath::Sqrt(rnEyh*rnEyh+rdEyh*rdEyh);
	  eyl            = reyl*rY;
	  eyh            = reyh*rY;
	}
	// Printf("den err (%f%%,%f%%) - blob (%f%%,%f%%)", dEyl*100, dEyh*100,
	//        eyl*100,eyh*100);
	ret->SetSysError(bId, cnt, eyl, eyh);
      }

      // Now we need to update our various point-to-point
      // errors. We've stored a mask for each declared point-to-point
      // error.  Each mask consist of three fields, corresponding to
      // point-to-point errors in the denominator, numerator, and the
      // ratio.  If all three fields are set, we shold cancel between
      // numerator and denominator.  If only the ratio and either of
      // the denominator or numerator fields are set, we simply need
      // to copy that to the ratio point-to-point error - as a
      // relative error of course.
      for (Int_t im = 0; im < shared.GetSize(); im++) {
	Int_t mask = shared[im];
	if (mask <= 0) break;
	Int_t rId = ((mask >>  0) & 0xFF);
	Int_t nId = ((mask >>  8) & 0xFF);
	Int_t dId = ((mask >> 16) & 0xFF);
	if (nId == 0 && dId == 0) {
	  ::Warning("Ratio", "Both numerator and denominator IDs are 0");
	  break;
	}
	Bool_t   rel   = ret->IsRelative(rId);
	Bool_t   crel  = true; // rel
	Double_t pdEyl = 0;
	Double_t pdEyh = 0;
	Double_t pnEyl = 0;
	Double_t pnEyh = 0;
	if (dId != 0) {
	  HolderP2P* pD  =  den->FindP2P(dId);
	  Double_t   lfc = (crel ? (TMath::Abs(dY) > 1e-9 ? 1/dY : 0) : 1);
	  pdEyl          =  pD->GetYDown(di1, di2, daX) * lfc;
	  pdEyh          =  pD->GetYUp  (di1, di2, daX) * lfc;
	}
	if (nId != 0) {
	  HolderP2P* pN  =  num->FindP2P(nId);
	  Double_t   lfc = (crel ? (TMath::Abs(nY) > 1e-9 ? 1/nY : 0) : 1);
	  pnEyl          =  pN->GetYDown(i) * lfc;
	  pnEyh          =  pN->GetYUp  (i) * lfc;
	}
	// Printf("Doing error %d from n=%d (%f +%f -%f) d=%d (%f +%f -%f)",
	//        rId, nId, nY, pnEyl, pnEyh, dId, dY, pdEyh, pdEyl);
	Double_t peyl = 0;
	Double_t peyh = 0;
	Double_t lfc  = (!crel ? (TMath::Abs(dY) > 1e-9 ? 1/dY : 0) : 1);
	if (nId == 0) {
	  // No numerator value - set directly
	  peyl = pdEyl * lfc;
	  peyh = pdEyh * lfc;
	} // else if
	else if (dId == 0) {
	  // No denomenator value - set directly
	  peyl = pnEyl * lfc;
	  peyh = pnEyh * lfc;
	} // else if 
	else {
	  // Cancel errors
	  if (denrel) {
	    peyl = pdEyl*lfc;
	    peyh = pdEyh*lfc;
	    // Printf("den err (%f%%,%f%%)", peyl*100,peyh*100);
	  }
	  else if (cmin) {
	    peyl = TMath::Min(nEyl, dEyl) * lfc;
	    peyh = TMath::Min(nEyh, dEyh) * lfc;
	  }
	  else if (cmax) {
	    peyl = TMath::Min(nEyl, dEyl) * lfc;
	    peyh = TMath::Min(nEyh, dEyh) * lfc;
	  }
	  else {
	    peyl = TMath::Sqrt(TMath::Abs(pnEyl*pnEyl-pdEyl*pdEyl)) * lfc;
	    peyh = TMath::Sqrt(TMath::Abs(pnEyh*pnEyh-pdEyh*pdEyh)) * lfc;
	    // ::Info("Ratio", "sqrt(|%8f^2 - %8f^2|)=%8f", pnEyl,pdEyl,peyl);
	  }
	} // Else
	// Printf("%f, den err (%f%%,%f%%) - %d (%f%%,%f%%)",
	//        rY, pdEyl*100, pdEyh*100, rId, peyl*100,peyh*100);
	lfc = (rel ? 1 : rY);
	ret->SetSysError(rId, cnt, 0, 0, lfc*peyl, lfc*peyh);
      }

      // Finally, increment counter
      cnt++;
    }

    // Reset the used case bits
    num->ClearUsed();
    den->ClearUsed();

    if (cnt <= 0) {
      ::Warning("", "No common points found");
      delete ret;
      ret = 0;
      return 0;
    }
    return ret;
  }
  /** 
   * Clear bit we set during the processing 
   * 
   */
  void ClearUsed() const
  {
    // Reset the used case bits
    TObject* oe = 0;
    TIter nextC(&fCommon);
    while ((oe = nextC())) {
      oe->ResetBit(kUsedBit);
      oe->ResetBit(kOnlyWeightBit);
    }
    TIter nextP(&fPoint2Point);
    while ((oe = nextP())) {
      oe->ResetBit(kUsedBit);
      oe->ResetBit(kOnlyWeightBit);
    }
  }
  /** 
   * Perform a Kolomogorov-Smironov test.  See als 
   *
   * https://root.cern.ch/root/htmldoc/TH1#TH1:KolmogorovTest
   *
   * @param g1 First graph.
   * @param g2 Second graph
   * 
   * @return Kolomogorov-Smirnov probability
   */
  static Double_t KolomogorovTest(const GraphSysErr* g1,
				  const GraphSysErr* g2)
  {
    Double_t z;
    return KolomogorovTest(g1, g2, z);
  }
  /** 
   * Perform a Kolomogorov-Smironov test.  See als 
   *
   * https://root.cern.ch/root/htmldoc/TH1#TH1:KolmogorovTest
   *
   * @param g1 First graph.
   * @param g2 Second graph
   * @param z  On return, Kolomogorov-Smirnov test statistic
   * 
   * @return Kolomogorov-Smirnov probability
   */
  static Double_t KolomogorovTest(const GraphSysErr* g1,
				  const GraphSysErr* g2,
				  Double_t&          z)
  {
    TArrayD  a1y;
    TArrayD  a2y;
    TArrayD  a1e2;
    TArrayD  a2e2;
    Int_t    cnt   = CacheGraphs(g1, g2, a1y, a2y, a1e2, a2e2);
    if (cnt <= 0) return -1;

    Double_t s1    = a1y.GetSum();  // Summed content
    Double_t s2    = a2y.GetSum();  // Summed content
    Double_t sw1   = a1e2.GetSum(); // Summed weights
    Double_t sw2   = a2e2.GetSum(); // Summed weights
    if (a1e2.GetSum() <= 0 && a2e2.GetSum()) {
      ::Warning("ChisquareTest", "No errors");
      return -1;
    }

    Double_t sr1 = 0;
    Double_t sr2 = 0;
    Double_t dfMax = 0;
    for (Int_t i = 0; i < cnt; i++) {
      // KolmogorowvTest
      sr1 += a1y[i] / s1;
      sr2 += a2y[i] / s2;
      dfMax =  TMath::Max(dfMax, TMath::Abs(sr1-sr2));
    }
    // K-S test
    Double_t se1 = s1 * s1 / sw1;
    Double_t se2 = s2 * s2 / sw2;
    z   = dfMax * TMath::Sqrt(se1 * se2 / (se1 + se2));
    return TMath::KolmogorovProb(z);
  }
  /** 
   * Get the Y and error values for the two passed graphs.  The graphs
   * are evaluated at the X values of Graph 1.  This is used for Chi
   * square and Kolomogorov-Smirnov tests
   * 
   * @param g1    Graph 1
   * @param g2    Graph 2
   * @param a1y   Y values of Graph 1
   * @param a2y   Y values of Graph 2 evaluated at X of Graph1 
   * @param a1e2  Y-error values of Graph 1			       
   * @param a2e2  Y-error values of Graph 2 evaluated at X of Graph1 
   * 
   * @return Count of points, or -1 on error 
   */
  static Int_t CacheGraphs(const GraphSysErr* g1,
			   const GraphSysErr* g2,
			   TArrayD&           a1y,
			   TArrayD&           a2y,
			   TArrayD&           a1e2,
			   TArrayD&           a2e2)
  {
    Int_t    cnt    = 0;
    Int_t    n      = g1->GetN();
    a1y.Set(n);
    a2y.Set(n);
    a1e2.Set(n);
    a2e2.Set(n);
    for (Int_t i = 0; i < g1->GetN(); i++) {
      Double_t x    = 0;
      Double_t y1   = 0;
      Double_t e1yl = 0;
      Double_t e1yh = 0;
      Double_t e2yl = 0;
      Double_t e2yh = 0;
      Double_t y2   = 0;
      if (!NextPoint(i, g1, g2, x, y1, e1yl, e1yh, y2, e2yl, e2yh)) {
	::Warning("ChisquareTest", "Next point - %d %f not found", i, x);
	continue;
      }

      Double_t e1y2   = TMath::Power(TMath::Max(e1yl,e1yh), 2);
      Double_t e2y2   = TMath::Power(TMath::Max(e2yl,e2yh), 2);	

      // Cache the numbers 
      a1y[cnt]  =  y1;   // Store value 
      a2y[cnt]  =  y2;   // Store value 
      a1e2[cnt] =  e1y2; // Store square error
      a2e2[cnt] =  e2y2; // Store square error 
      cnt++;
    }
    if (cnt <= 0) 
      ::Warning("CacheGraphs", "No common points found");
    
    return cnt;
  }
  /** 
   * Calculate the @f$ \chi^2@f$ test of equivilance between two
   * graphs.  See also 
   * 
   *  https://root.cern.ch/root/htmldoc/TH1#TH1:Chi2Test
   * 
   * @param g1     First graph
   * @param g2     Second graph
   * @param type   Type 
   * 
   * @return reduced @f$\chi^2@f$ or -1 on error
   */
  static Double_t ChisquareTest(const GraphSysErr* g1,
				const GraphSysErr* g2,
				EChi2Type          type=kExperimentTruth)
  {
    Int_t ndf = 0;
    return ChisquareTest(g1, g2, ndf, type);
  }
  /** 
   * Calculate the @f$ \chi^2@f$ test of equivilance between two
   * graphs.  See also 
   * 
   *  https://root.cern.ch/root/htmldoc/TH1#TH1:Chi2Test
   * 
   * @param g1     First graph
   * @param g2     Second graph
   * @param ndf    On return, number degrees of freedom 
   * @param type   Type 
   * 
   * @return reduced @f$\chi^2@f$ or -1 on error
   */
  static Double_t ChisquareTest(const GraphSysErr* g1,
				const GraphSysErr* g2,
				Int_t&             ndf,
				EChi2Type          type)
  {
    Double_t chi2   = 0;
    ndf             = -1;
    TArrayD  a1y;
    TArrayD  a2y;
    TArrayD  a1e2;
    TArrayD  a2e2;
    Int_t    cnt   = CacheGraphs(g1, g2, a1y, a2y, a1e2, a2e2);
    if (cnt <= 0) return -1;

    Double_t s1    = a1y.GetSum();  // Summed content
    Double_t s2    = a2y.GetSum();  // Summed content
    if (type == kModelModel && (a1e2.GetSum() <= 0 && a2e2.GetSum())) {
      ::Warning("ChisquareTest", "No errors");
      return -1;
    }

    ndf = cnt;
    
    // Now loop over cached points
    Double_t s = s1 + s2;
    for (Int_t i = 0; i < ndf; i++) {
      switch (type) {
      case kExperimentExperiment:  {
	Double_t sum   = a1y[i]+a2y[i];
	Double_t delta = s2 * a1y[i] - s1 * a2y[i];
	chi2           += delta * delta / sum;
      }
	break;
      case kExperimentModel: {
	Double_t v1 =  s2 * a1y[i] - s1 * a2e2[i];
	Double_t v2 =  v1 * v1 + 4 * s2 * s2 * a1y[i] * a2e2[i];
	Double_t pp =  (v1 + v2) / (2 * s2 * s2);
	Double_t p1 =  pp * s1;
	Double_t p2 =  pp * s2;
	Double_t d1 =  a1y[i] - p1;
	Double_t d2 =  a1y[i] - p2;
	chi2        += d1 * d1 / p1;
	if (a2e2[i] > 0) chi2 += d2 * d2 / a2e2[i];
      }
	break;
      case kModelModel: {
	Double_t sigma =  s1 * s1 * a2e2[i] + s2 * s2 * a1e2[i];
	Double_t delta =  s2 * a1y[i] - s1 * a2y[i];
	chi2           += delta * delta / sigma;
      }
	break;
      case kExperimentTruth: {
	Double_t delta =  a2y[i] - a1y[i];
	chi2           += delta * delta / a1e2[i];
      }
	break;
      default:
	::Warning("ChisquareTest", "Should not happen");
	return -1;
      }
	
    }
    if      (type == kExperimentExperiment) chi2 /= s1 * s2;
    else if (type == kExperimentTruth)      chi2 /= ndf;

    return chi2;
  }
  
  /** 
   * Remove a point from the graph 
   * 
   * @param i Point to remove 
   */
  void RemovePoint(Int_t i)
  {
    if (i < 0 || i >= GetN()) return;
    fData->RemovePoint(i);
    TIter next(&fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(next()))) {
      p2p->fGraph->RemovePoint(i);
    }
  }
  /** 
   * Swap two points 
   * 
   * @param i        Index 
   * @param j        Index 
   * @param reflect  if true, multiply X values with -1
   */
  void SwapPoints(Int_t i, Int_t j, Bool_t reflect=false)
  {
    if (i == j && !reflect) return;
    SwapPoints(fData, i, j, reflect);
    TIter next(&fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(next()))) {
      SwapPoints(p2p->fGraph, i, j, reflect);
    }
  }
  /** 
   * Make a copy of the graph, and reflect around x0
   * 
   * @param x0 Where to reflect around 
   * 
   * @return Newly allocated graph 
   */
  GraphSysErr* Reflect(Double_t x0=0) const
  {
    GraphSysErr* cpy = new GraphSysErr(*this);
    for (Int_t i = 0; i < GetN()/2; i++) {
      cpy->SwapPoints(i, GetN()-i-1, true);
    }
    return cpy;
  }
  /** 
   * Symmetrice the other graph and store result here. 
   * 
   * @param other Graph to summetrice 
   * 
   * @return true on success
   */
  Bool_t Symmetrize(GraphSysErr* other)
  {
    GraphSysErr*  cpy = new GraphSysErr(*other);

    // Declare "blob" systematic point-to-point error 
    Int_t         id  = DeclarePoint2Point("Syst.unc..", false, kBox);

    // Loop over common errors 
    HolderCommon* cmn = 0;
    TIter         cNext(&fCommon);
    while ((cmn = static_cast<HolderCommon*>(cNext()))) {
      Int_t oId = cpy->FindId(cmn->GetTitle());
      if (oId <= 0) {
	// If the common error ins't found in copy - it should,
	// because we did a straight copy - then give up.
	Error("Symmetrice", "Common error %s not found in %s",
	      cmn->GetTitle(), cpy->GetName());
	cpy->Print("all");
	delete cpy;
	cpy = 0;
	break;
	// }
	// }
      } // oId <= 0

      // Zero this error here - we should add it later       
      HolderCommon* oCmn = cpy->FindCommon(oId);
      if (cmn->fEyl <= 0 && cmn->fEyh <= 0) {
	// Copy the error over 
	cmn->fEyl = oCmn->fEyl;
	cmn->fEyh = oCmn->fEyh;
      }
      oCmn->Set(0,0);
    }  // while cmn
    if (!cpy) return false;

    // Loop over all points 
    Int_t   cnt = 0;
    Int_t   n   = other->GetN();
    TArrayI used(n); // book-keeping
    Int_t   nonZero = 0; // how many non-zero sys.uncerr.
    used.Reset(-1);
    for (Int_t i = 0; i < n; i++) {
      if (used[i] >= 0 && used[i] != i) {
	Int_t j = used[i];
	// Info("Symmetrice", "Copy %d to %d instead of %d", j, cnt, i);
	SetPoint(cnt, -GetX(j), GetY(j));
	SetPointError(cnt, GetErrorXRight(j), GetErrorXLeft(j));
	SetStatError(cnt, GetStatErrorDown(j), GetStatErrorUp(j));
	SetSysError(id, cnt, 0, GetSysErrorY(id, j));
	// Info("Symmetrice", "%f -> %f +/- %f",
	//      -GetX(j), GetY(j), GetSysErrorY(id,j));
	cnt++;
	continue;
      }
      used[i] = i;
      Double_t eyl1, eyh1;
      Double_t exl1  = cpy->GetErrorXLeft(i);
      Double_t exh1  = cpy->GetErrorXRight(i);
      Double_t seyl1 = cpy->GetStatErrorDown(i);
      Double_t seyh1 = cpy->GetStatErrorUp(i);
      Double_t y1    = cpy->GetYandError(i,true,false,true,false,eyl1,eyh1);
      Double_t x1    = cpy->GetX(i);
      Double_t y2    = 0;
      Double_t eyl2  = 0;
      Double_t eyh2  = 0;
      Double_t seyl2 = 0;
      Double_t seyh2 = 0;
      Int_t    i1, i2;
      Int_t    j  = cpy->FindPoint(-x1, i1, i2); // Find mirror 
      if (j < -1) {
	// IF there's no mirror point, the just set to current 
	// Info("Symmetrice", "No mirror of %f using %d for %d", x1, i, cnt);
	SetPoint(cnt, x1, y1);
	SetPointError(cnt, exl1, exh1);
	SetStatError(cnt, seyl1, seyh1);
	SetSysError(id, cnt, 0, 0, eyl1, eyh1);
	cnt++;
	continue;
      }
      if (j >= 0) {
	// If we have exact mirror, then use that point and mark as
	// used
	// Info("Symmetrice", "Got exact mirror of %f(%d) at %d", x1, i, j);
	y2    = cpy->GetYandError(j,true,false,true,false,eyl2,eyh2);
	seyl2 = cpy->GetStatErrorDown(j);
	seyh2 = cpy->GetStatErrorUp(j);
	// We should copy values when we get to j
	used[j] = cnt;
      }
      else {
	// Otherwise we use interpolation between two points 
	// Info("Symmetrice", "Interpolate mirror of %f(%d) between %d-%d",
	//      x1, i, i1, i2);
	cpy->FindYandError(-x1,true,true,true,false,y2,eyl2,eyh2,seyl2,seyh2);
      }
      // Now calculate weighted mean and variance
      Double_t newY = (y1+y2) / 2;
      Double_t newV = 0;
      if ((eyh1+eyl1) > 1e-9 && (eyh2+eyl2) > 1e-9) {
	Double_t s1   = 2 * eyl1 * eyh1 / (eyh1 + eyl1);
	Double_t sp1  = (eyh1 - eyl1) / (eyh1 + eyl1);
	Double_t w1   = .5 * TMath::Power(s1+y1*sp1, 3) / s1;
	Double_t s2   = 2 * eyl2 * eyh2 / (eyh2 + eyl2);
	Double_t sp2  = (eyh2 - eyl2) / (eyh2 + eyl2);
	Double_t w2   = .5 * TMath::Power(s2+y2*sp2, 3) / s2;
	// Printf("s1=%f sp1=%f w1=%f s2=%f sp2=%f w2=%f",
	///       s1, sp1, w1, s2, sp2, w2);
	Double_t sumW = (w1+w2);
	newY          = (y1*w1+y2*w2) / sumW;
	newV          = TMath::Sqrt((w1*w1*s1*s1+w2*w2*s2*s2) / sumW / sumW);
	nonZero++;
      }
      Double_t seyl = TMath::Sqrt(seyl1*seyl1+seyl2*seyl2);
      Double_t seyh = TMath::Sqrt(seyh1*seyh1+seyh2*seyh2);
      
      // Info("Symmmetrice", "%f (%f +/- %f) + (%f +/- %f) -> %f +/- %f",
      //      x1, y1, w1, y2, w2, newY, newV);
      SetPoint(cnt, x1, newY);
      SetPointError(cnt, exl1, exh1);
      SetStatError(cnt, seyl, seyh);
      SetSysError(id, cnt, 0, newV);
      cnt++;
    }
    if (nonZero <= 0) RemoveSysError(id);
    // Info("", "After symmetrization:");
    // Print("all");
    return true;
  }
  /* @} */
  /** 
   * @{
   * @name Import/export 
   */
  void SavePrimitive(std::ostream& out, Option_t* option="")
  {
    TString opt(option);
    opt.ToLower();
    Bool_t load = opt.Contains("load"); opt.ReplaceAll("load", "");
    TString funcName;
    TPRegexp regex("func=([a-zA-z][a-zA-Z0-9_]*)");
    TObjArray* toks = regex.MatchS(opt);
    if (toks) {
      if (toks->GetEntriesFast() > 1) 
	funcName = toks->At(1)->GetName();
      if (toks->GetEntriesFast() > 0) 
	opt.ReplaceAll(toks->At(0)->GetName(), "");
      delete toks;
    }
      
    
    // Save initialization 
    if (!funcName.IsNull()) 
      out << "TObject* " << funcName << "(Option_t* o=\"\")\n";
    out << "{\n";
    if (load)
      out << " // Load class\n"
	  << "  if (!gROOT->GetClass(\"GraphSysErr\"))\n"
	  << "    gROOT->LoadMacro(\"GraphSysErr.C+\");\n";
    out << "  GraphSysErr* g = new GraphSysErr(\""
	<< GetName() << "\",\"" << GetTitle() << "\","
	<< GetN() << ");\n";
    // Save attributes 
    out << "  // Point options\n"
	<< "  g->SetMarkerStyle("  << GetMarkerStyle() << ");\n"
	<< "  g->SetMarkerColor("  << GetMarkerColor() << ");\n"
	<< "  g->SetMarkerSize("   << GetMarkerSize() << ");\n"
	<< "  g->SetLineStyle("    << GetLineStyle() << ");\n"
	<< "  g->SetLineColor("    << GetLineColor() << ");\n"
	<< "  g->SetLineWidth("    << GetLineWidth() << ");\n"      
	<< "  g->SetFillStyle("    << GetFillStyle() << ");\n"
	<< "  g->SetFillColor("    << GetFillColor() << ");\n"
	<< "  g->SetXTitle(\""     << fXTitle << "\");\n"
	<< "  g->SetYTitle(\""     << fYTitle << "\");\n"
	<< "  g->SetDataOption("   << fDataOption << ");\n"
	<< "  // Sum options\n"
	<< "  g->SetSumOption("    << fSumOption << ");\n"
	<< "  g->SetSumTitle(\""   << fSumTitle << "\");\n"
	<< "  g->SetSumLineStyle(" << fSumLine.GetLineStyle() << ");\n"
	<< "  g->SetSumLineColor(" << fSumLine.GetLineColor() << ");\n"
	<< "  g->SetSumLineWidth(" << fSumLine.GetLineWidth() << ");\n"      
	<< "  g->SetSumFillStyle(" << fSumFill.GetFillStyle() << ");\n"
	<< "  g->SetSumFillColor(" << fSumFill.GetFillColor() << ");\n"
	<< "  g->SetCommonSumOption("    << fCommonSumOption << ");\n"
	<< "  g->SetCommonSumTitle(\""   << fCommonSumTitle << "\");\n"
	<< "  g->SetCommonSumLineStyle(" <<fCommonSumLine.GetLineStyle()<<");\n"
	<< "  g->SetCommonSumLineColor(" <<fCommonSumLine.GetLineColor()<<");\n"
	<< "  g->SetCommonSumLineWidth(" <<fCommonSumLine.GetLineWidth()<<");\n"
	<< "  g->SetCommonSumFillStyle(" <<fCommonSumFill.GetFillStyle()<<");\n"
	<< "  g->SetCommonSumFillColor(" <<fCommonSumFill.GetFillColor()<<");\n"
	<< "  // Stat options\n"
	<< "  g->SetStatRelative(" << fStatRelative << ");\n";
    TIter nextC(&fCommon);
    HolderCommon* cmn = 0;
    while ((cmn = static_cast<HolderCommon*>(nextC())))
      cmn->SavePrimitive(out, "d");
    TIter nextP(&fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(nextP())))
      p2p->SavePrimitive(out, "d");
    Int_t n = GetN();
    out << " // " << n << " points\n";
    Bool_t statRel = IsStatRelative();
    for (Int_t i = 0; i < n; i++) {
      Double_t y = GetY(i);
      out << "  g->SetPoint(" << i << ',' << GetX(i) << ',' << y << ");\n"
	  << "  g->SetPointError(" << i << ',' << GetErrorXLeft(i) << ','
	  << GetErrorXRight(i) << ");\n"
	  << "  g->SetStatError(" << i << ','
	  << (statRel ? 1/y : 1)*GetStatErrorDown(i) << ','
	  << (statRel ? 1/y : 1)*GetStatErrorUp(i) << ");\n";
      nextP.Reset();
      while ((p2p = static_cast<HolderP2P*>(nextP()))) {
	Int_t  id  = p2p->GetUniqueID();
	Bool_t rel = p2p->IsRelative();
	out << "  g->SetSysError(" << id << ',' << i << ','
	    << GetSysErrorXLeft(id, i) << ','
	    << GetSysErrorXRight(id, i) << ','
	    << (rel ? 1/y : 1) * GetSysErrorYDown(id, i) << ','
	    << (rel ? 1/y : 1) * GetSysErrorYUp(id, i) << ");\n";
      } // while(p2p)
    } // for (i)
    out << "  if (o && o[0] != '\\0') {\n"
	<< "    g->Draw(o);\n";
    if (fDrawn && fDrawn->GetHistogram())
      out << "    if (g->GetMulti() && g->GetMulti()->GetHistogram()) {\n"
	  << "      g->GetMulti()->GetHistogram()->SetMinimum("
	  << fDrawn->GetHistogram()->GetMinimum() << ");\n"
	  << "      g->GetMulti()->GetHistogram()->SetMaximum("
	  << fDrawn->GetHistogram()->GetMaximum() << ");"
	  << "    }\n";
   
    out << "  }\n";
    if (!funcName.IsNull()) out << "  return g;\n";
    out << "};\n";
  }
  /** 
   * Save to a ROOT script
   * 
   * @param fileName Script to write to 
   */
  void Save(const char* fileName)
  {
    std::ofstream out(fileName);
    TString funcName(fileName);
    funcName.ReplaceAll(".C","");
    out << "// \n"
	<< "// Generated by GraphSysErr.C\n"
	<< "// \n"
      // << "class GraphSysErr;\n\n";
	<< "\n";
    SavePrimitive(out, Form("load func=%s", funcName.Data()));
    out.close();
  }
  /** 
   * Dump on stream a table suitable (After some editing) for
   * uploading to the Durham database.
   * 
   * If one has many objects that should be uploaded together, one 
   * can do 
   * 
   * @code 
   TList l;
   l.Add(new GraphSysErr(...));
   ...

   std::ofstream out("export");
   TIter next(&l);
   GraphSysErr* g = 0;
   Bool_t first = true;
   while ((g = static_cast<GraphSysErr*>(next()))) {
   g->Export(out, (first ? "h" : ""));
   first = false;
   }
   out.close();
   * @endcode 
   * 
   * @param out    Output stream to write to. 
   * @param option Options
   * @param nsign  Number of significant digits on errors 
   *
   * - H Export file header 
   * - C Export file comment 
   * - S Export Point-to-point systematic names
   * - T Export title as RE qual
   */
  void Export(std::ostream& out=std::cout,
	      Option_t*     option="",
	      Int_t         nsign=2)
  {
    // Info("Export", "Using %d significant digits on errors", nsign);
    TString opt(option);
    opt.ToLower();
    Bool_t header   = opt.Contains("h");
    Bool_t sysNames = opt.Contains("s");
    Bool_t comment  = opt.Contains("c");
    Bool_t title    = opt.Contains("t");
    
    ExportHeader(out, header, comment, nsign);
    
    // --- Export qualifiers -----------------------------------------
    Bool_t hasTitle = false;
    if (fQualifiers) {
      TIter nextQ(fQualifiers);
      TObject* q = 0;
      while ((q = nextQ())) {
	TString k(q->GetName());
	if (k.EqualTo("RE") || k.EqualTo("title", TString::kIgnoreCase))
	  hasTitle = true;
	out << FormatKey("qual") << q->GetName() << " : "
	    << q->GetTitle() << std::endl;
      }
    }
    if (!hasTitle && title)
      out << FormatKey("qual") << "RE : " << GetTitle() << std::endl;
    
    // --- Export X/Y titles ----------------------------------------
    const char* fill = "<please fill in>";
    TString xTit = fXTitle; EscapeLtx(xTit, "fill");
    TString yTit = fYTitle; EscapeLtx(yTit, "fill");
    
    out << FormatKey("xheader") << xTit << "\n"
	<< FormatKey("yheader") << yTit << std::endl;

    // --- Export common errors --------------------------------------
    TIter nextC(&fCommon);
    HolderCommon* holderCommon = 0;
    while ((holderCommon = static_cast<HolderCommon*>(nextC()))) { 
      Bool_t rel = holderCommon->IsRelative();
      Double_t up = holderCommon->GetYUp(rel ? 100 : 1);
      Double_t down = holderCommon->GetYDown(rel ? 100 : 1);
      out << FormatKey("dserror");
      ExportError(out, down, up, true, rel, nsign);
      out << ":" << holderCommon->GetTitle() << std::endl;
    }

    // --- Export data points ----------------------------------------
    out  << FormatKey("data") << " x : y" << std::endl;
    Int_t n = GetN();
    for (Int_t i = 0; i < n; i++) {
      ExportPoint(out, i, true, sysNames, nsign);
      out << std::endl;
    }
    out << "*dataend:\n" 
	<< "# End of dataset\n" << std::endl;
  } //*MENU*
    /** 
     * Utility to escape out TLatex stuff, and put '$...$' around LaTeX 
     * 
     * @param val  String to modify 
     * @param fill If val is null, use this value 
     */
  static void EscapeLtx(TString& val, const TString& fill="")
  {
    if (val.IsNull()) {val = fill; return; }
    if (!val.Contains("#") && !val.Contains("\\")) return;

    if (val[0]             !='$') val.Prepend("$");
    if (val[val.Length()-1]!='$') val.Append("$");
    val.ReplaceAll("#", "\\");
  }
  /** 
   * Export a set of data sets to a single table.  All graphs must
   * have the same format.  The title of each graph is written as the
   * "qual" field.
   * 
   * @param col    Collection of GraphSysErr objets 
   * @param out    Output stream 
   * @param option Options
   * @param nsign  Number of significant digits  
   *
   * - H Export file header 
   * - C Export file comment 
   * - S Export Point-to-point systematic names
   */
  static void Export(const TSeqCollection* col,
		     std::ostream& out,
		     Option_t* option="H",
		     Int_t nsign=0)
  {
    if (col->GetEntries() < 1) return;

    // --- Deduce options --------------------------------------------
    TString opt(option);
    opt.ToLower();
    Bool_t alsoTop = opt.Contains("h");
    Bool_t alsoCmt = opt.Contains("c");
    Bool_t alsoNme = opt.Contains("s");
    Bool_t alsoTit = opt.Contains("t");

    // --- some variables to use -------------------------------------
    const Double_t tol = 1e-10;
    GraphSysErr*   first = 0;
    GraphSysErr*   gse = 0;
    TList          toExport;
    TList          commons;
    TList          quals;
    TString        data("x ");
    Int_t          nPoints = -1;
    Int_t          idx     = 0;

    // --- Check the passed graphs for compatiblity ------------------
    TIter          nextCheck(col);
    TObject*       o       = 0;
    while ((o = nextCheck())) {
      if (!o->IsA()->InheritsFrom(GraphSysErr::Class())) continue;
      gse = static_cast<GraphSysErr*>(o);
      if (!first) {
	first   = gse;
	nPoints = first->GetN();
      }
      else {
	// --- Check number of points --------------------------------
	if (gse->GetN() != nPoints) {
	  Int_t nTmp = TMath::Min(gse->GetN(), nPoints);
	  ::Warning("Export", "Incompatible number of points %d in %s"
		    "using only %d of them", 
		    gse->GetN(), gse->GetName(), nTmp);
	  nPoints = nTmp;
	}
	// --- check X values ----------------------------------------
	Bool_t ok = true;
	for (Int_t i = 0; i < nPoints; i++) {
	  Double_t x1    = first->GetX(i);
	  Double_t exl1  = first->fData->GetErrorXlow(i);
	  Double_t exh1  = first->fData->GetErrorXhigh(i);
	  Double_t xT    = gse->GetX(i);
	  Double_t exlT  = gse->fData->GetErrorXlow(i);
	  Double_t exhT  = gse->fData->GetErrorXhigh(i);
	  if (TMath::Abs(x1-xT) > tol ||
	      (exl1 > tol && TMath::Abs(exl1-exlT) > tol) ||
	      (exh1 > tol && TMath::Abs(exh1-exhT) > tol)) {
	    ::Warning("Export", "X--coordinate of %s @ point %d: %f (+%f,-%f) "
		      "incompatible [%f (+%f,-%f)]", 
		      gse->GetTitle(), i, xT, exhT, exlT, x1, exh1, exl1);
	    ok = false;
	    break;
	  }
	} // for i      
	if (!ok) continue;
      } // !first

	// --- Get all possible qualifiers -----------------------------
      toExport.Add(gse);
      TIter    nextQ(gse->fQualifiers);
      TObject* q = 0;
      while ((q = nextQ())) {
	// ::Info("", "qualifier %s=%s", q->GetName(), q->GetTitle());
	StoreQual(quals, idx, q);
      }

      // --- Get all common systematics ------------------------------
      TIter    nextCmn(&(gse->fCommon));
      TObject* oCmn = 0;
      while ((oCmn = nextCmn())) {
	if (commons.FindObject(oCmn->GetTitle())) continue;
	TObjString* cmn = new TObjString(oCmn->GetTitle());
	if (gse == first) cmn->SetUniqueID(gse->FindId(oCmn->GetTitle()));
	commons.Add(cmn);
      }
      
      idx++;
      data += ": y ";
    }
    // --- Now deduce the common common errors -----------------------
    TIter    nextCmn(&(commons));
    TObject* oCmn = 0;
    while ((oCmn = nextCmn())) {
      TString  oNme(oCmn->GetName());
      Bool_t   found = true;
      Double_t ecl1  = -1; 
      Double_t ech1  = -1;
      
      // --- Loop over data ------------------------------------------
      TIter  nextG(&toExport);
      while ((gse = static_cast<GraphSysErr*>(nextG()))) {
	/// chekc if we have this error 
	Int_t gseId = first->FindId(oNme);

	if (gseId > 0) {
	  // If we do, check compatibility  
	  if (ecl1 < 0) {
	    // First time we get this error 
	    ecl1 = gse->GetCommonErrorYDown(gseId);
	    ech1 = gse->GetCommonErrorYUp(gseId);
	    continue;
	  }
	  else {
	    // Now check values within tolerance 
	    Double_t eclT  = gse->GetCommonErrorYDown(gseId);
	    Double_t echT  = gse->GetCommonErrorYUp(gseId);
	    if ((ecl1 > tol && TMath::Abs(ecl1-eclT) > tol) ||
		(ech1 > tol && TMath::Abs(ech1-echT) > tol)) {
	      found = false;
	    } // incompatible
	  } // ecl1 >= 0
	} // gseId > 0
	else 
	  // This graph does not have this common, flag it 
	  found = false;
      } // while(gse)
#if INCOMPAT_CMN_AS_QUAL
#else 
      found = true;
#endif 
      oCmn->SetBit(BIT(14), found);

      // If we found the error in all graphs and they are all
      // compatible, then all is good,
      if (found) continue;

      // If not, we represent the common systematic as a qualifier 
      idx   = 0;      
      nextG.Reset();
      while ((gse = static_cast<GraphSysErr*>(nextG()))) {
	Int_t gseId = first->FindId(oNme);
	TString val;
	if (gseId > 0) {
	  // This graph has the value, store as qual
	  Bool_t   rel = gse->IsRelative(gseId);
	  Double_t ecl = gse->GetCommonErrorYDown(gseId);
	  Double_t ech = gse->GetCommonErrorYUp(gseId);
	  std::stringstream s;
	  ExportError(s, ecl, ech, false, rel, nsign);
	  val = s.str().c_str();
	}
	else {
	  // This does not have the error, store empty string
	  val = "";
	}
	// Now store this in our qualifier table 
	StoreQual(quals, idx, oNme, val);
	idx++;
      } // while(gse)
    } // while common

      // --- Sanity checks ---------------------------------------------
    if (nPoints <= 0) {
      ::Error("Export", "No points to write");
      return;
    }
    if (toExport.GetEntries() <= 0) {
      ::Error("Export", "No graphs to export");
      return;
    }
    if (!first) {
      ::Error("Export", "Didn't get the first graph");
      return;
    }

    // --- Export header --------------------------------------------
    first->ExportHeader(out, alsoTop, alsoCmt, nsign);

    // --- Export qualifiers -----------------------------------------
    // quals.ls();
    Bool_t hasTitle = false;
    TIter nextQ(&quals);
    TList* ql = 0;
    while ((ql = static_cast<TList*>(nextQ()))) {
      TString k(ql->GetName());
      if (k.EqualTo("RE") || k.EqualTo("title", TString::kIgnoreCase))
	hasTitle = true;
      out << FormatKey("qual") << ql->GetName();
      for (Int_t i = 0; i < toExport.GetEntries(); i++) {
	TObject* qv = ql->At(i);
	TString  v  = "";
	if (qv)  v  = qv->GetName();
	out << " : " << v;
      }
      out << std::endl;
    }
    if (!hasTitle && alsoTit) {
      out << FormatKey("qual") << "RE ";
      for (Int_t i = 0; i < toExport.GetEntries(); i++) 
	out << ": " << toExport.At(i)->GetTitle();
      out << std::endl;
    }

    
    // --- Export axis titles ----------------------------------------
    const char*  fill     = "<please fill in>";
    const char*  fields[] = { "xheader", "yheader", 0 };
    const char** pfld     = fields;
    while (*pfld) { 
      out << FormatKey(*pfld);
      TIter nextSpec(&toExport);
      Bool_t one = true;
      while ((gse = static_cast<GraphSysErr*>(nextSpec()))) {
	TString val;
	if      ((*pfld)[0] == 'x') val = gse->fXTitle;
	else if ((*pfld)[0] == 'y') val = gse->fYTitle;
	EscapeLtx(val, fill);
	out << (one ? "" : ":") << val;
	if ((*pfld)[0] == 'x') break;
	one = false;
      }
      pfld++;
      out << std::endl;
    }
      
    // --- Export common systematics ---------------------------------
    TIter    nextC(&commons);
    while ((oCmn = nextC())) {
      TString tit(oCmn->GetName());
      EscapeLtx(tit);
      out << FormatKey("dserror");
#if INCOMPAT_CMN_AS_QUAL
      if (!oCmn->TestBit(BIT(14))) {
	::Warning("Export",
		  "Common systematic error \"%s\" represented by qual",
		  oCmn->GetName());
	out << "- : " << tit << std::endl;
	continue;
      }
      UInt_t   id  = first->FindId(oCmn->GetName());
      Bool_t   rel = first->IsRelative(id);
      Double_t ecl = first->GetCommonErrorYDown(id);
      Double_t ech = first->GetCommonErrorYUp(id);
      ExportError(out, ecl, ech, true, rel, nsign);
      out << " : "<< tit << std::endl;
#else 
      out << tit << std::flush;
      TIter  next(&toExport);
      while ((gse = static_cast<GraphSysErr*>(next()))) {
	UInt_t   id  = gse->FindId(oCmn->GetName());
	Bool_t   rel = gse->IsRelative(id);
	Double_t ecl = gse->GetCommonErrorYDown(id);
	Double_t ech = gse->GetCommonErrorYUp(id);
	out << " : ";
	ExportError(out, ecl, ech, true, rel, nsign);
      }
      out << ":" << tit << std::endl;
#endif
    }

    // --- Export points ---------------------------------------------
    out << FormatKey("data") << data << std::endl;
    for (Int_t i = 0; i < nPoints; i++) {
      TIter  next(&toExport);
      Bool_t one = true;
      while ((gse = static_cast<GraphSysErr*>(next()))) {
	gse->ExportPoint(out, i, one, alsoNme, nsign);
	one = false;
      }
      out << std::endl;
    }
    out << "*dataend:\n" 
	<< "# End of dataset\n" << std::endl;
    // return out;
  }
  /** 
   * Import all data sets from a Durham input formatted file.  The
   * data sets are returned in a flat collection.  The collection owns
   * the contained the objects and it is the responsibility of the
   * caller to manage the returned collection.  
   * 
   * The returned GraphSysErr objects are named like 
   *
   * @verbatim 
   *   ds_INDEX_SUB-INDEX  
   * @endverbatim 
   *
   * @param fileName Name of file to read. 
   * 
   * @return Pointer to newly allocated collection of graphs, or 0 in
   * case of errors.  Note that the returned collection may be empty
   */
  static TSeqCollection* Import(const TString& fileName)
  {
    std::ifstream in(fileName.Data());
    if (!in) {
      ::Warning("Import", "Failed to open \"%s\"", fileName.Data());
      return 0;
    }
    TList* ret = new TList;
    ret->SetOwner();
    GraphSysErr* g     = 0;
    GraphSysErr* first = 0;
    Int_t id = 0;
    do {
      Int_t    sub  = 1;
      UShort_t nIdx = 256;
      int      cur  = in.tellg();
      do {
	in.seekg(cur, in.beg);
	if (!(g = Import(in, sub, &nIdx))) break;
	if (kVerbose & kImport)
	  ::Info("Import", "Imported %d of %d", sub, nIdx);
	g->SetName(Form("ds_%d_%d", id, sub-1));
	if (!first) first = g;
	else        g->CopyKeys(first);
	ret->Add(g);
	sub++;
      } while(sub < nIdx);
      id++;
    } while (!in.eof());
    return ret;
  }
  /** 
   * Import data from <i>input</i> formatted Durham database file. 
   * 
   * @image html dndeta_pa.png 
   *
   * The source @b input file of the above graph is
   * @include "inputs/dndeta_pa.input"
   *
   * @note This is not omnipotent.  While the function works for most
   * typical inputs. it can fail for particular inputs.
   *
   * To read multiple graphs, use GraphSysErr::Import(const TString&)
   * @param in   Input stream
   * @param idx  Column of data set to import. 
   * @param nIdx If non-null, holds the number of available columns on return. 
   * 
   * @return Newly allocated object or null
   *
   */
  static GraphSysErr* Import(std::istream& in, UShort_t idx=0,
			     UShort_t* nIdx=0)
  {
    GraphSysErr* ret    = 0;
    Int_t        n      = 0;
    Bool_t       inSet  = false;
    Bool_t       inData = false;
    TString      xtit   = "";
    TString      ytit   = "";
    TString      tmp    = "";
    TString      last   = "";
    TList        quals;
    TList        keys;
    Int_t        isty   = 0;
    do {
      TString line;
      line.ReadLine(in);
      tmp = line.Strip(TString::kBoth, ' '); line = tmp;
      if (line.IsNull()) continue; 
      if (line[0] == '#') continue;
      if (line[0] == '*') { // We have some sort of header stuff 
	Int_t colon = line.Index(":"); 
	if (colon == kNPOS) continue;
	Int_t   len   = line.Length()-colon-1;
	TString value = line(colon+1,len);
	tmp = value.Strip(TString::kBoth, ' '); value = tmp;
	last = line(1,colon-1);
	if (kVerbose & kImport) 
	  ::Info("Import", "Got a key '%s' -> '%s'", last.Data(), value.Data());

	if (last.BeginsWith("dataset", TString::kIgnoreCase)) {
	  inSet = true;
	  inData = false;
	  isty  = 0;
	}
	else if (last.BeginsWith("dataend", TString::kIgnoreCase)) {
	  inSet = false;
	  inData = false;
	  // Always terminate on the end of a data set
	  break;
	}
	// else if (last.BeginsWith("dscomment",TString::kIgnoreCase)) {
	// tit = value;
	// }
	else if (last.BeginsWith("qual", TString::kIgnoreCase)) {
	  // ::Info("", "Got qual: %s", value.Data());
	  // if (!qual.IsNull()) {
	  //   ::Info("", "Adding seen qual: %s", qual.Data());
	  quals.Add(new TObjString(value));
	}
	else if (last.BeginsWith("xheader", TString::kIgnoreCase))
	  xtit = value;
	else if (last.BeginsWith("yheader", TString::kIgnoreCase))
	  ytit = value;
	else if (last.BeginsWith("dserror", TString::kIgnoreCase)) {
	  // Common systematic error 
	  if (kVerbose & kImport) 
	    ::Info("Import", "Got common error line: '%s'", line.Data());
	  TObjArray* tokens = value.Tokenize(":");
	  Double_t el = 0, eh = 0;
	  Bool_t   rel = false;
	  if (ImportError(Token(tokens, 0), el, eh, rel)) { 
	    if (rel) { el /= 100.; eh /= 100.; }
	    TString& tnam = Token(tokens, 1);
	    TString  nam  = tnam.Strip(TString::kBoth, ' ');
	    if (!ret) ret = new GraphSysErr("imported", "");
	    Int_t id = ret->DefineCommon(nam, rel, el, eh, kRect);
	    ret->SetSysLineColor(id, (isty % 6) + 2);
	    ret->SetSysFillColor(id, (isty % 6) + 2);
	    ret->SetSysFillStyle(id, 3001 + (isty % 10));
	    isty++;
	  }
	}
	else if (last.BeginsWith("data", TString::kIgnoreCase)) {
	  if (kVerbose & kImport) 
	    ::Info("Import", "Got start of data line: '%s'", line.Data());
	  // These are the field we can deal with 
	  // Let's get the header field value 
	  TObjArray* tokens = value.Tokenize(":");
	  if (nIdx) {
	    *nIdx = tokens->GetEntriesFast();
	    // ::Info("Import", "Max index is %d", *nIdx);
	  }
	  if (tokens->GetEntriesFast() > 2 && idx < 1) { 
	    idx = 1;
	    ::Warning("Import", "Can only import one data set at a time, "
		      "selecting the %d column", idx);
	  }
	  else if (idx >= 1 && idx >= tokens->GetEntriesFast()) { 
	    ::Warning("Import", "column %d not available for this data set",
		      idx);
	    inSet = false;
	  }
	  // We do nothing to the tokens here. 
	  tokens->Delete();
	  inData = true;
	}
	else {	 
	  // some other key. 
	  if (kVerbose & kImport) 
	    ::Info("Import", "Got pair line '%s' (%s -> %s)", line.Data(), 
		   last.Data(), value.Data());
	  keys.Add(new TNamed(last, value));
	}
	continue; 
      }
      if (!inData) { 
	if (kVerbose & kImport) 
	  ::Info("Import", "Got a possible contiuation line '%s'", line.Data());
	// Probably a continuation line. 
	if (last.IsNull()) 
	  // No last key or no map 
	  continue; 

	if (kVerbose & kImport) 
	  ::Info("Import", "Got some other line: '%s'", line.Data());

	TNamed* ptr = static_cast<TNamed*>(keys.FindObject(last));
	if (!ptr) 
	  // Key is not added, do nothing to this line 
	  continue; 

	TString tt(ptr->GetTitle());
	if (!tt.EndsWith(" ") && !tt.EndsWith("\t") && !tt.EndsWith("\n")) 
	  // If value does not end in a white-space, add it. 
	  tt.Append("\n");

	// Append our line, and set the title
	tt.Append(line);
	ptr->SetTitle(tt);

	continue;
      }
      if (!inSet) { 
	// We're outside a data-set so do nothing
	continue;
      }
      if (idx < 1) idx = 1;

      if (kVerbose & kImport) 
	::Info("Import", "Got a data line: '%s'", line.Data());
      TObjArray* tokens = line.Tokenize(";");
      if (tokens->GetEntriesFast() < idx+1) {
	::Warning("Import", "Too few columns %d<%0d in line: %s", 
		  tokens->GetEntriesFast(), idx+1, line.Data());
	tokens->Delete();
	continue;
      }
      TString& xCol = static_cast<TObjString*>(tokens->At(0))->String();
      TString& yCol = static_cast<TObjString*>(tokens->At(idx))->String();
      if (kVerbose & kImport) 
	::Info("Import", "xColumn: '%s', yColumn: '%s'", 
	       xCol.Data(), yCol.Data());

      tmp = xCol.Strip(TString::kBoth, ' '); xCol = tmp;
      tmp = yCol.Strip(TString::kBoth, ' '); yCol = tmp;
      TString yFull = yCol;

      Int_t chop = yCol.Last('(');
      if (chop != kNPOS) {
	yCol.Remove(chop, yCol.Length()-chop);
	yCol.Remove(TString::kBoth, ' ');
      }
      if (kVerbose & kImport) 
	::Info("Import", "xColumn: '%s', yColumn: '%s'", 
	       xCol.Data(), yCol.Data());

      if (xCol.IsNull()) { 
	::Warning("Import", "Empty X column in line: %s", line.Data());
	tokens->Delete();
	continue;
      }
      if (yCol.IsNull()) { 
	::Warning("Import", "Empty Y column in line: %s", line.Data());
	tokens->Delete();
	continue;
      }
      Bool_t rel = false;
      Double_t x = 0, exl = 0, exh = 0;
      if (!ImportPoint(xCol, x, exl, exh, rel)) { 
	::Warning("Import", "Failed to import X value in line: %s",line.Data());
	tokens->Delete();
	continue;
      }
      Double_t y = 0, eyl = 0, eyh = 0;
      if (!ImportPoint(yCol, y, eyl, eyh, rel)) {
	::Warning("Import", "Failed to import X value in line: %s",line.Data());
	tokens->Delete();
	continue;
      }
      if (!ret) ret = new GraphSysErr("imported", "");
      if (rel) ret->SetStatRelative(true);
      ret->SetPoint(n, x, y);
      ret->SetPointError(n, exl, exh);

      if (ret->IsStatRelative()) { eyl /= 100; eyh /= 100; }
      ret->SetStatError(n, eyl, eyh);

      Int_t lparen = yFull.Index("(");
      Int_t rparen = yFull.Index(")");
      if (lparen == rparen) { 
	n++;
	continue;
      }
      TString    rem  = yFull(lparen+1, rparen-lparen-1);
      if (kVerbose & kImport) 
	::Info("Import", "Got systematic errors: '%s'", rem.Data());

      TObjArray* stok = rem.Tokenize("D");
      Int_t      isys = 0;
      for (Int_t i = 0; i < stok->GetEntriesFast(); i++) {
	TString t = Token(stok,i);
	if (!t.BeginsWith("SYS=",TString::kIgnoreCase)) { 
	  ::Warning("Import", "Ignoring unknown spec %s (@ %d) in %s",
		    t.Data(), i, rem.Data());
	  continue;
	}
	t.Remove(0,4);
	if (t[t.Length()-1] == ',') t.Remove(t.Length()-1,1);
	tmp = t.Strip(TString::kBoth, ' '); t = tmp;

	Double_t el = 0, eh = 0;
	rel = false;
	if (!ImportError(t, el, eh, rel)) continue;

	Int_t colon = t.Index(":");
	TString nam(Form("sys%d", isys));
	if (colon != kNPOS) nam = t(colon+1, t.Length()-colon-1);
	  
	UInt_t id = ret->FindId(nam);
	if (id == 0) { 
	  id = ret->DeclarePoint2Point(nam, rel, kRect);
	  ret->SetSysLineColor(id, (isty % 6) + 2);
	  ret->SetSysFillColor(id, (isty % 6) + 2);
	  ret->SetSysFillStyle(id, 3001 + (isty % 10));
	  isty++;
	}
	HolderP2P* p = ret->FindP2P(id);
	if (p->IsRelative()) { el /= 100.; eh /= 100.; }
	ret->SetSysError(id, n, exl, exh, el, eh);
	isys++;
      }
      stok->Delete();

      n++;      
    } while (!in.eof());
    if (ret) {
      TString rxtit = ExtractField(xtit, idx-1);
      TString rytit = ExtractField(ytit, idx-1);
      if (!rxtit.IsNull()) ret->SetXTitle(rxtit);
      if (!rytit.IsNull()) ret->SetYTitle(rytit);

      Bool_t hasTitle = false;
      TIter nextQ(&quals);
      TObject* qual = 0;
      while ((qual = nextQ())) {
	TString k = ExtractField(qual->GetName(), 0);
	TString q = ExtractField(qual->GetName(), idx);
	/*
	  ::Info("LoopQ", "Qualifier string: %s -> %s,%s",
	  qual->GetName(), k.Data(), q.Data());
	*/
	ret->AddQualifier(k, q);
	if (!hasTitle &&
	    (k.EqualTo("RE") ||
	     k.EqualTo("title", TString::kIgnoreCase))) {
	  hasTitle = true;
	  ret->SetTitle(q);
	}
      }
      TIter nextK(&keys);
      TObject* pair = 0;
      while ((pair = nextK())) {
	ret->SetKey(pair->GetName(),pair->GetTitle());
	TString k = pair->GetName();
	if (!hasTitle && (k.EqualTo("location"))) {
	  ret->SetTitle(pair->GetTitle());
	  hasTitle = true;
	}
      }
      if (!hasTitle) ret->SetTitle("Title");
    }
    keys.SetOwner();
    keys.Delete();
    return ret;
  }
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Declaring systematic errors
   */ 
  /** 
   * Define a common systematic error 
   * 
   * @param title     Title of error 
   * @param relative  True if this relative to data 
   * @param ey        Error     
   * @param option    Options
   * 
   * @return Indentifier of systematic error
   *
   * Example of how make define common errors
   * @dontinclude Example.C
   * @skip cm1
   * @until cm2
   */
  UInt_t DefineCommon(const char* title, Bool_t relative, 
		      Double_t ey, EDrawOption_t option=kFill)
  {
    return DefineCommon(title, relative, ey, ey, option);
  }
  /** 
   * Define a common systematic error 
   * 
   * @param title     Title of error 
   * @param relative  True if this relative to data 
   * @param eyl       Error     
   * @param eyh       Error     
   * @param option    Options
   * 
   * @return Indentifier of systematic error
   */
  UInt_t DefineCommon(const char* title, Bool_t relative, 
		      Double_t eyl, Double_t eyh, 
		      EDrawOption_t option=kRect)
  {
    UInt_t  id = ++fCounter;
    TString name(Form("common_%08x", id));

    HolderCommon* h = new HolderCommon(name, title, relative, option, id);
    h->Set(eyl, eyh);
    
    fCommon.AddLast(h);
    return id;
  }
  /** 
   * Delcare a point-to-point systematic error 
   * 
   * @param title     Title
   * @param relative  Relative error mission
   * @param option    Options
   * 
   * @return Indentifier of systematic error
   *
   * Example of how make declare point-to-point errors
   * @dontinclude Example.C
   * @skip pp1
   * @until pp2
   */
  UInt_t DeclarePoint2Point(const char* title, Bool_t relative, 
			    EDrawOption_t option=kBar)
  {
    UInt_t  id = ++fCounter; 
    TString name(Form("p2p_%08x", id));

    Holder* h = new HolderP2P(name, title, relative, option, id);

    fPoint2Point.AddLast(h);
    return id;
  }
  void RemoveSysError(Int_t id)
  {
    HolderCommon* h = FindCommon(id);
    if (h) {
      fCommon.Remove(h);
      return;
    }
    HolderP2P* p = FindP2P(id);
    if (p) {
      fPoint2Point.Remove(p);
      return;
    }
  }
  /** 
   * Find the ID of an error with the given title 
   * 
   * @param title Title 
   * 
   * @return ID or null
   */
  UInt_t FindId(const char* title) const
  { 
    TString tit(title);

    TIter nextC(&fCommon);
    Holder* holder = 0;
    while ((holder = static_cast<Holder*>(nextC()))) {
      if (tit.EqualTo(holder->GetTitle())) return holder->GetUniqueID();
    }

    TIter nextP(&fPoint2Point);
    while ((holder = static_cast<Holder*>(nextP()))) {
      if (tit.EqualTo(holder->GetTitle())) return holder->GetUniqueID();
    }
    // fCommon.ls();
    // fPoint2Point.ls();
    return 0;
  }
  Int_t GetNSys() const { return fCounter; }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Setting the data and errors 
   */
  /** 
   * Set the ith data point
   * 
   * @param i 
   * @param x 
   * @param y 
   */
  void SetPoint(Int_t i, Double_t x, Double_t y) 
  {
    if (!fData) return;
    fData->SetPoint(i, x, y);
  }
  /** 
   * Set the X error (bin width) of the ith point. 
   * 
   * @param i    Point 
   * @param ex   Error 
   */
  void SetPointError(Int_t i, Double_t ex)
  {
    SetPointError(i, ex, ex);
  }
  /** 
   * Set the X error (bin width) of the ith point. 
   * 
   * @param i    Point 
   * @param exl  Lower error 
   * @param exh  Upper error 
   */
  void SetPointError(Int_t i, Double_t exl, Double_t exh)
  {
    if (!fData) return;
    fData->SetPointError(i, exl, exh, 
			 fData->GetErrorYlow(i),
			 fData->GetErrorYhigh(i));
  }
  /** 
   * Set whether statistical errors should be considered relative.
   * 
   * @param rel If true, statistical errors are specified relative to
   * the y values.
   */
  void SetStatRelative(Bool_t rel) { fStatRelative = rel; }
  /**
   * Check if statistical errors are relative 
   *
   * @return true if statistical errors are defined as relative 
   */
  Bool_t IsStatRelative() const { return fStatRelative; }
  /** 
   * Set the statistical error on the ith data point
   * 
   * @param i   Point number 
   * @param ey  Error on Y
   */
  void SetStatError(Int_t i, Double_t ey) 
  {
    SetStatError(i, ey, ey);
  }
  /*
   * Set the statistical error on the ith data point
   * 
   * @param i   Point number 
   * @param eyl Lower error on Y
   * @param eyh Higher error on Y
   */
  void SetStatError(Int_t i, Double_t eyl, Double_t eyh) 
  {
    if (!fData) return;
    if (fStatRelative) { 
      eyl *= fData->GetY()[i];
      eyh *= fData->GetY()[i];
    }
    fData->SetPointError(i,
			 fData->GetErrorXlow(i),
			 fData->GetErrorXhigh(i),
			 eyl, eyh);
  }
  void SetSysError(Int_t id, Double_t eyl, Double_t eyh)
  {
    HolderCommon* h = FindCommon(id);
    if (!h) return;
    h->Set(eyl, eyh);
  }
  /** 
   * Set the systematic error identified by id on the ith data point 
   * 
   * @param id   Systematic error identifier 
   * @param i    Point number (starting at 0)
   * @param ex   X error 
   * @param ey   Y error
   */
  void SetSysError(Int_t id, Int_t i, Double_t ex, Double_t ey)
  {
    HolderP2P* h = FindP2P(id);
    if (!h) return;
    h->Set(i, fData, ex, ey);
  }
  /** 
   * Set the systematic error identified by id on the ith data point 
   * 
   * @param id   Systematic error identifier 
   * @param i    Point number (starting at 0)
   * @param exl  Left X error 
   * @param exh  Right X error 
   * @param eyl  Lower Y error 
   * @param eyh  Upper Y error
   */
  void SetSysError(Int_t id, Int_t i, Double_t exl, Double_t exh,
		   Double_t eyl, Double_t eyh) 
  {
    HolderP2P* h = FindP2P(id);
    if (!h) return;
    h->Set(i, fData, exl, exh, eyl, eyh);
  }
  /* @} */
  //__________________________________________________________________
  /**
   * @{ 
   * @name Setting drawing options 
   */
  /** 
   * Set the title of the data 
   * 
   * @param name Title 
   */
  void SetTitle(const char* name)
  {
    TNamed::SetTitle(name);
    if (fData) fData->SetTitle(name);
  }
  /** 
   * Set the draw option for the data and statistical errors
   * 
   * @param opt Draw option 
   */
  void SetDataOption(EDrawOption_t opt) { fDataOption = opt; } //*MENU*
  /** 
   * Set title on X axis
   * 
   * @param title 
   */
  void SetXTitle(const char* title) { fXTitle = title; } //*MENU*
  /** 
   * Set title on y axis
   * 
   * @param title 
   */
  void SetYTitle(const char* title) { fYTitle = title; } //*MENU*
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Getting information 
   */
  /** 
   * @return number of points
   */  
  Int_t    GetN() const { return fData ? fData->GetN() : 0; }
  /** 
   * @param p Point
   * 
   * @return X at point
   */
  Double_t GetX(Int_t p) const { return fData->GetX()[p]; }
  /** 
   * @param p Point 
   *
   * @return Left error on X at point
   */
  Double_t GetErrorXLeft(Int_t p) const { return fData->GetErrorXlow(p); }
  /** 
   * @param p Point 
   *
   * @return Right error on X at point
   */
  Double_t GetErrorXRight(Int_t p) const { return fData->GetErrorXhigh(p); }
  /** 
   * @param point Point
   * 
   * @return Y at point
   */
  Double_t GetY(Int_t point) const { return fData->GetY()[point]; }
  /** 
   * @param point Point
   * 
   * @return statistical error at point
   */
  Double_t GetStatError(Int_t point) const { return fData->GetErrorY(point); }
  /** 
   * @param point Point
   * 
   * @return statistical error at point
   */
  Double_t GetStatErrorUp(Int_t point) const
  { 
    return fData->GetErrorYhigh(point); 
  }
  /** 
   * @param point Point
   * 
   * @return statistical error at point
   */
  Double_t GetStatErrorDown(Int_t point) const
  {
    return fData->GetErrorYlow(point);
  }
  Bool_t IsCommon(Int_t id) const
  {
    HolderCommon* c = FindCommon(id);
    return c != 0;
  }
  /** 
   * Check if an error is relative 
   * 
   * @param id 
   * 
   * @return 
   */
  Bool_t IsRelative(Int_t id) const
  {
    HolderP2P* h = FindP2P(id);
    if (h) return h->IsRelative();
    HolderCommon* c = FindCommon(id);
    if (c) return c->IsRelative();
    return false;
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorX(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (!h) return 0;
    return h->GetX(point);
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorY(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (!h) return 0;
    return h->GetY(point);
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorXLeft(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (!h) return 0;
    return h->GetXLeft(point);
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorXRight(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (!h) return 0;
    return h->GetXRight(point);
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorYUp(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (h) return h->GetYUp(point);
    HolderCommon* c = FindCommon(id);
    if (c) return c->GetYUp(GetY(point));
    return 0;
  }
  /**
   * @param id
   * @param point Point
   * 
   * @return systematic error from id at point
   */
  Double_t GetSysErrorYDown(Int_t id, Int_t point) const 
  { 
    HolderP2P* h = FindP2P(id);
    if (h) return h->GetYDown(point);
    HolderCommon* c = FindCommon(id);
    if (c) return c->GetYDown(GetY(point));
    return 0;
  }
  /** 
   * Get title of systematic error 
   * 
   * @param id Identifier 
   * 
   * @return Name 
   */
  const char* GetSysTitle(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return "";
    return h->GetTitle();
  }
  /** 
   * Get fill style of systematic uncertainty 
   * 
   * @param id Identifier 
   * 
   * @return style, or 0
   */
  Style_t GetSysFillStyle(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return 0;
    return h->GetFillStyle();    
  }
  /** 
   * Get line style of systematic uncertainty 
   * 
   * @param id Identifier 
   * 
   * @return style, or 0
   */
  Style_t GetSysLineStyle(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return 0;
    return h->GetLineStyle();    
  }
  /** 
   * Get fill color of systematic uncertainty 
   * 
   * @param id Identifier 
   * 
   * @return color, or 0
   */
  Color_t GetSysFillColor(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return 0;
    return h->GetFillColor();    
  }
  /** 
   * Get line color of systematic uncertainty 
   * 
   * @param id Identifier 
   * 
   * @return color, or 0
   */
  Color_t GetSysLineColor(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return 0;
    return h->GetLineColor();    
  }
  /** 
   * Get line width of systematic uncertainty 
   * 
   * @param id Identifier 
   * 
   * @return width, or 0
   */
  Width_t GetSysLineWidth(Int_t id) const
  {
    Holder* h = Find(id);
    if (!h) return 0;
    return h->GetLineWidth();    
  }
  /** 
   * Get title of summed systematic error 
   * 
   * @return Name 
   */
  const char* GetSumTitle() const
  {
    return fSumTitle.Data();
  }
  UInt_t GetSumOption() const { return fSumOption; }
  /** 
   * Get fill style of sum systematic uncertainty 
   * 
   * @return style, or 0
   */
  Style_t GetSumFillStyle() const
  {
    return fSumFill.GetFillStyle();
  }
  /** 
   * Get line style of sum systematic uncertainty 
   * 
   * @return style, or 0
   */
  Style_t GetSumLineStyle() const
  {
    return fSumLine.GetLineStyle();
  }
  /** 
   * Get fill color of sum systematic uncertainty 
   * 
   * @return color, or 0
   */
  Color_t GetSumFillColor() const
  {
    return fSumFill.GetFillColor();
  }
  /** 
   * Get line color of sum systematic uncertainty 
   * 
   * @return color, or 0
   */
  Color_t GetSumLineColor() const
  {
    return fSumLine.GetLineColor();
  }
  /** 
   * Get line width of sum systematic uncertainty 
   * 
   * @return width, or 0
   */
  Width_t GetSumLineWidth() const
  {
    return fSumLine.GetLineWidth();
  }
  /** 
   * Get title of summed systematic error 
   * 
   * @return Name 
   */
  const char* GetCommonSumTitle() const
  {
    return fCommonSumTitle.Data();
  }
  UInt_t GetCommonSumOption() const { return fCommonSumOption; }
  /** 
   * Get fill style of sum systematic uncertainty 
   * 
   * @return style, or 0
   */
  Style_t GetCommonSumFillStyle() const
  {
    return fCommonSumFill.GetFillStyle();
  }
  /** 
   * Get line style of sum systematic uncertainty 
   * 
   * @return style, or 0
   */
  Style_t GetCommonSumLineStyle() const
  {
    return fCommonSumLine.GetLineStyle();
  }
  /** 
   * Get fill color of sum systematic uncertainty 
   * 
   * @return color, or 0
   */
  Color_t GetCommonSumFillColor() const
  {
    return fCommonSumFill.GetFillColor();
  }
  /** 
   * Get line color of sum systematic uncertainty 
   * 
   * @return color, or 0
   */
  Color_t GetCommonSumLineColor() const
  {
    return fCommonSumLine.GetLineColor();
  }
  /** 
   * Get line width of sum systematic uncertainty 
   * 
   * @return width, or 0
   */
  Width_t GetCommonSumLineWidth() const
  {
    return fCommonSumLine.GetLineWidth();
  }
  /** 
   * Get the common systematic error 
   *
   * @return Common systematic error 
   */
  Double_t GetCommonErrorYUp(Int_t id) const
  {
    HolderCommon* c = FindCommon(id);
    if (!c) return 0;
    Bool_t rel = c->IsRelative();
    return c->GetYUp(rel ? 100 : 1);
  }
  /** 
   * Get the common systematic error 
   *
   * @return Common systematic error 
   */
  Double_t GetCommonErrorYDown(Int_t id) const
  {
    HolderCommon* c = FindCommon(id);
    if (!c) return 0;
    Bool_t rel = c->IsRelative();
    return c->GetYDown(rel ? 100 : 1);
  }
  /** 
   * Get name of X axis 
   * 
   * @return X-axis name
   */
  const char* GetXTitle() const { return fXTitle.Data(); }
  /** 
   * Get name of Y axis 
   * 
   * @return Y-axis name
   */
  const char* GetYTitle() const { return fYTitle.Data(); }
  /** 
   * Get minimum and maximum 
   * 
   * @param option 
   * @param ymin 
   * @param ymax 
   */
  void GetMinMax(Option_t* option,
		 Double_t& ymin, Double_t& ymax) const
  {
    Double_t xmin, xmax;
    Int_t    imin, imax;
    GetMinMax(option, ymin, ymax, xmin, xmax, imin, imax);
  }
  /** 
   * Get minimum and maximum 
   * 
   * @param option Options
   * @param ymin   On return, lease Y value 
   * @param ymax   On return, largest Y value 
   * @param xmin   On return, X value corresponding to least Y value 
   * @param xmax   On return, X value corresponding to largest Y value 
   * @param imin   On return, point number corresponding to lesat Y value
   * @param imax   On return, point number corresponding to largest Y value
   */
  void GetMinMax(Option_t* option,
		 Double_t& ymin, Double_t& ymax,
		 Double_t& xmin, Double_t& xmax,
		 Int_t&    imin, Int_t&    imax) const
  {
    TString  opt(option); opt.ToUpper();
    Bool_t   cmn     = opt.Contains("COMMON");
    Bool_t   stat    = opt.Contains("STAT");
    Bool_t   quad    = opt.Contains("QUAD");
    Bool_t   noerr   = opt.Contains("Y");
    quad             = !opt.Contains("DIRECT");

    ymin = +1e9;
    ymax = -1e9;
    xmin = -1e9;
    xmax = +1e9;
    imin = -1;
    imax = -1;
    for (Int_t i = 0; i < GetN(); i++) {
      Double_t eyl = 0;
      Double_t eyh = 0;
      Double_t y   = noerr ? GetY(i) : GetYandError(i, cmn, stat, quad, false, eyl, eyh);
      Double_t yl  = y - eyl;
      Double_t yh  = y + eyh;
      if (yl < ymin) {
	ymin = yl;
	xmin = GetX(i);
	imin = i;
      }
      if (yh > ymax) {
	ymax = yh;
	xmax = GetX(i);
	imax = i;
      }
      // Printf("%3d: %f -> %f +%f -%f -> %f %f -> %f %f",
      //        i, GetX(i), y, eyl, eyh, yl, yh, ymin, ymax);
    }
    // Printf("min=(%d,%f,%f) max=(%d,%f,%f)", imin, xmin, ymin, imax, xmax, ymax);
  }
  /** 
   * Find the full-width at half-maximum 
   * 
   * @param start Starting point 
   * @param dir   Direction (-1: to the left, +1: to the right)
   * @param ymax  The maximum value 
   * @param cmn   Whether to include common systematics
   * @param stat  Whether to include statistical errors 
   * @param quad  Whether to add in quadrature 
   * @param i1    On return, the lower bound for found point
   * @param i2    On return, the upper bound for found point
   */
  void FindFwhm(Int_t start, Int_t dir, Double_t ymax,
		Bool_t cmn, Bool_t stat, Bool_t quad,
		Int_t& i1, Int_t& i2) const
  {
    Int_t lim = (dir < 0) ? -1 : GetN();
    i1 = i2   = -1;
    for (Int_t i = start+dir; i != lim; i += dir) {
      Double_t eyl, eyh;
      Double_t y = GetYandError(i, cmn, stat, quad, false, eyl, eyh);
      // Printf("%3d: %2d %f -> %f -%f +%f (%f, %f, %f) -> %2d %2d",
      //        i, dir, GetX(i), y, eyl, eyh, y-eyl, y+eyh, ymax/2,
      //        i1, i2);
      if (y-eyl > ymax/2) continue; // Still above
      if (y+eyh < ymax/2) break; // Found lower limit
      if (i1 < 0)  i1 = i2 = i;
      else         i2 = i;
    }
    // Printf("Found %d,%d", i1, i2);
  }
  /** 
   * Calculate the full-width at half-maximum 
   * 
   * @param el On return, lower error on FWHM
   * @param eh On return, higher error on FWHM
   * 
   * @return FWHM
   */
  Double_t FWHM(Double_t& el, Double_t& eh) const
  {
    Double_t xl, xh;
    return FWHM(el, eh, xl, xh);
  }
  /** 
   * Calculate the full-width at half-maximum 
   * 
   * @param el On return, lower error on FWHM
   * @param eh On return, higher error on FWHM
   * @param xl On return, left-hand X value
   * @param xh On return, right-hand X value
   * 
   * @return FWHM
   */
  Double_t FWHM(Double_t& el, Double_t& eh, Double_t& xl, Double_t& xh) const
  {
    Double_t ymin, ymax, xmin, xmax;
    Int_t    imin, imax;
    GetMinMax("quad stat", ymin, ymax, xmin, xmax, imin, imax);
    if (ymin > ymax/2) {
      Warning("FWHM", "Half of ymax=%f is out of reach [%f,%f]",
	      ymax/2, ymin, ymax);
      return -1;
    }

    // Point bounds of found half-max
    Int_t li1, li2, ri1, ri2; 
    // Search left 
    FindFwhm(imax, -1, ymax, false, true, true, li1, li2);
    // Search right 
    FindFwhm(imax, +1, ymax, false, true, true, ri1, ri2);

    Double_t lsign = 1;
    // Check if we've found no left point
    if (li1 < 0 && li2 < 0) {
      // Check if we've found no right point
      if (ri1 < 0 && ri2 < 0) {
	Warning("FWHM", "No left and right point found");
	return -1;
      }
      // Set sign, and copy right point to left point
      lsign = -1;
      li2   = ri1;
      li1   = ri2;
    }
    Double_t rsign = 1;
    // Check if we've found no right point
    if (ri1 < 0 && ri2 < 0) {
      // Check if we've found no left point
      if (li1 < 0 && li2 < 0) {
	Warning("FWHM", "No right and left point found");
	return -1;
      }
      // Set sign, and copy left point to right point
      rsign = -1;
      ri2   = li1;
      ri1   = li2;
    }
    // Printf("Points low=(%d,%d) high=(%d,%d)", li1, li2, ri1, ri2);

    // Now find best estimate of left X
    LinearSigmaCombiner scl;
    scl.Add(lsign*GetX(li1), GetErrorXLeft(li1), GetErrorXRight(li1));
    scl.Add(lsign*GetX(li2), GetErrorXLeft(li2), GetErrorXRight(li2));
    Combiner::Result* lr = scl.Calculate();
    scl.Print();
    lr->Print();
    
    // Now find best estimate of right X
    LinearSigmaCombiner scr;
    scr.Add(GetX(ri1), GetErrorXLeft(ri1), GetErrorXRight(ri1));
    scr.Add(GetX(ri2), GetErrorXLeft(ri2), GetErrorXRight(ri2));
    Combiner::Result* rr = scr.Calculate();
    scr.Print();
    rr->Print();

    // Get error on left/right hand sides
    el            = lr->fEh+rr->fEl;
    eh            = lr->fEl+rr->fEh;
    if (el < 1e-9 && eh < 1e-9) {
      // If the errors are very small, we take the average
      xl   = (GetX(li2)+GetX(li1))/2;
      xh   = (GetX(ri2)+GetX(ri1))/2;
      el   = (GetX(li2)-GetX(li1))/2;
      eh   = (GetX(ri2)-GetX(ri1))/2;
    }
    else {
      // Otherwise we set left/right hand values to best-fit values
      xl   = lr->fX;
      xh   = rr->fX;
    }
    // Evaluate the full-width half-maximum
    Double_t fwhm = (xh - xl);
    return fwhm; 
  }
    
  /** 
   * Calculate the mean along X.  The strategy used here is the same
   * as for TH1::GetMean - that is the weighted average using the Y
   * values as weights.
   *
   * If both the mean and standard deviation are needed, one can use
   * the function StatisticsX directly, and calculate the errors as
   * outlined there.
   * 
   * @param error On return, the error on the mean 
   * @param cmn       If true, use common uncertainties 
   * @param stat      If true, use statistical errors 
   * @param quad      If true, add in quadrature 
   * 
   * @return The weighted of X 
   */
  Double_t MeanX(Double_t& error,
		 Bool_t cmn=false,
		 Bool_t stat=true,
		 Bool_t quad=true) const
  {
    Double_t meanX=TMath::Infinity(), stdDev, n;
    StatisticsX(meanX, stdDev, n, cmn, stat, quad);
    error = stdDev/n;
    return meanX;
  }
  /** 
   * Calculate the mean along X.  The strategy used here is the same
   * as for TH1::GetMean - that is the weighted average using the Y
   * values as weights.
   * 
   * If both the mean and standard deviation are needed, one can use
   * the function StatisticsX directly, and calculate the errors as
   * outlined there.
   * 
   * @param cmn       If true, use common uncertainties 
   * @param stat      If true, use statistical errors 
   * @param quad      If true, add in quadrature 
   * 
   * @return The mean  
   */
  Double_t MeanX(Bool_t cmn=false,
		 Bool_t stat=true,
		 Bool_t quad=true) const
  {
    Double_t error;
    return MeanX(error, cmn, stat, quad);
  }
  /** 
   * Calculate the standard deviatin along X.  The strategy used here
   * is the same as for TH1::GetMean - that is the weighted average
   * using the Y values as weights.
   * 
   * If both the mean and standard deviation are needed, one can use
   * the function StatisticsX directly, and calculate the errors as
   * outlined there.
   * 
   * @param error     On return, the error on the standard deviation
   * @param cmn       If true, use common uncertainties 
   * @param stat      If true, use statistical errors 
   * @param quad      If true, add in quadrature 
   * 
   * @return The standard deviation
   */
  Double_t StandardDeviationX(Double_t& error,
			      Bool_t cmn=false,
			      Bool_t stat=true,
			      Bool_t quad=true) const
  {
    Double_t meanX=TMath::Infinity(), stdDev, eStdDev, n;
    StatisticsX(meanX, stdDev, n, cmn, stat, quad);
    error = stdDev/TMath::Sqrt(2)/n;
    return stdDev;
  }
  /** 
   * Calculate the standard deviatin along X.  The strategy used here
   * is the same as for TH1::GetMean - that is the weighted average
   * using the Y values as weights.
   * 
   * If both the mean and standard deviation are needed, one can use
   * the function StatisticsX directly, and calculate the errors as
   * outlined there.
   * 
   * @param cmn       If true, use common uncertainties 
   * @param stat      If true, use statistical errors 
   * @param quad      If true, add in quadrature 
   * 
   * @return The standard deviation
   */
  Double_t StandardDeviationX(Bool_t cmn=false,
			      Bool_t stat=true,
			      Bool_t quad=true) const
  {
    Double_t error;
    return StandardDeviationX(error, cmn, stat, quad);
  }
  Double_t StandardDeviationXMean(Double_t mean,
				  Double_t& error,
				  Bool_t cmn=false,
				  Bool_t stat=true,
				  Bool_t quad=true) const
  {
    Double_t meanX, stdDev, eStdDev, n;
    StatisticsX(meanX, stdDev, n, cmn, stat, quad);
    Double_t sd2 = stdDev*stdDev - meanX*meanX;
    stdDev = TMath::Sqrt(sd2 - mean*mean);
    error  = stdDev/TMath::Sqrt(2)/n;
    return stdDev;
  }
    
  /** 
   * Calculates the statistics of X 
   * 
   * @f${eqnarray}
   *   \bar{x} &=& \frac{\sum xy}{\sum y}\\
   *   s_x     &=& \sqrt{\frac{\sum x^2y}{\sum y}-\bar{x}^2}\\
   *   n       &=& \frac{(\sum y)^2}{\sum \delta^2 y}\\
   * @f${eqnarray}
   *
   * The errors on @f$\bar{x}@f$ and @f$ s_x@f$ is given by 
   *
   * @f${eqnarray}
   *   \delta\bar{x} &=& s_x/n\\
   *   \delta s_x    &=& s_x/(\sqrt{2}n)\\
   * @f${eqnarray}
   * 
   * @note The variable passed for @a meanX must have the value @c inf
   * (use TMath::Infinity()) to get the calculated mean back.  Any
   * other value will be taken as the mean for the reminder of the
   * calculaitions.  This is useful if one wants to calculate the
   * standard deviation around some pre-determined mean rather than
   * the sample mean.
   *
   * @param meanX     On return, @f$\bar{x}@f$
   * @param stdDevX   On return, @f$ s_x@f$ 
   * @param n         On return, @f$ n @f$ (effective entries)
   * @param cmn       If true, use common uncertainties 
   * @param stat      If true, use statistical errors 
   * @param quad      If true, add in quadrature 
   */
  void StatisticsX(Double_t& meanX,
		   Double_t& stdDevX,
		   Double_t& n,
		   Bool_t cmn=false,
		   Bool_t stat=true,
		   Bool_t quad=true) const
  {
    Double_t sumY   = 0;
    Double_t sumXY  = 0;
    Double_t sumXXY = 0;
    Double_t sumE2  = 0;
    for (Int_t i = 0; i < GetN(); i++) {
      Double_t e2yl, e2yh;
      Double_t x = GetX(i);
      Double_t y = GetYandError(i,cmn,stat,quad,true,e2yl,e2yh);

      sumY   += y;
      sumXY  += x*y;
      sumXXY += x*x*y;
      sumE2  += TMath::Max(e2yl,e2yh);
    }
    if (meanX == TMath::Infinity())
      meanX       = sumXY/sumY;
    Double_t sd2  = TMath::Abs(sumXXY/sumY-meanX*meanX);
    stdDevX       = TMath::Sqrt(sd2);
    Double_t eff  = sumE2 > 0 ? sumY*sumY/sumE2 : sumY;
    n             = TMath::Sqrt(eff);
  }  
  /* @} */

  //__________________________________________________________________
  /**
   * @{ 
   * @name Setting attributes on systematic errors 
   */
  /** 
   * Set the line color of the systematice error identified by ID
   * 
   * @param id    Systematic error identifier 
   * @param color Line color 
   */
  void SetSysLineColor(Int_t id, Color_t color) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetLineColor(color);
  }
  /** 
   * Set the line style of the systematice error identified by ID
   * 
   * @param id    Systematic error identifier 
   * @param style Line style 
   */
  void SetSysLineStyle(Int_t id, Style_t style) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetLineStyle(style);
  }
  /** 
   * Set the line width of the systematice error identified by ID
   * 
   * @param id    Systematic error identifier 
   * @param width Line width
   */
  void SetSysLineWidth(Int_t id, Width_t width) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetLineWidth(width);
  }
  /** 
   * Set the fill color of the systematice error identified by ID
   * 
   * @param id    Systematic error identifier 
   * @param color Color 
   */
  void SetSysFillColor(Int_t id, Color_t color) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetFillColor(color);
  }
  /** 
   * Set the fill style of the systematice error identified by ID
   * 
   * @param id    Systematic error identifier 
   * @param style Fill style
   */
  void SetSysFillStyle(Int_t id, Style_t style) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetFillStyle(style);
  }
  /** 
   * Set the systematic error title 
   * 
   * @param id    Systematic error identifier 
   * @param name  The title 
   */
  void SetSysTitle(Int_t id, const char* name) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetTitle(name);
  }
  /** 
   * Set the draw option for a specific systematic error 
   * 
   * @param id  Systematic error identifier 
   * @param opt Draw option 
   */
  void SetSysOption(Int_t id, EDrawOption_t opt) //*MENU*
  {
    Holder* h = Find(id);
    if (!h) return;
    h->SetDOption(opt);
  }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Setting attributes on summed errors 
   */
  /** 
   * Set the draw option for summed errors
   * 
   * @param opt Draw option 
   */
  void SetSumOption(EDrawOption_t opt) { fSumOption = opt; } //*MENU*
  /** 
   * Set the title uses for summed errors
   * 
   * @param title Title
   */
  void SetSumTitle(const char* title) { fSumTitle = title; } //*MENU*
  /** 
   * Set the line color of the sumtematice error identified by ID
   * 
   * @param color Line Color 
   */
  void SetSumLineColor(Color_t color) { fSumLine.SetLineColor(color); } //*MENU*
  /** 
   * Set the line style of the sumtematice error identified by ID
   * 
   * @param style Line style 
   */
  void SetSumLineStyle(Style_t style){ fSumLine.SetLineStyle(style); } //*MENU*
  /** 
   * Set the line width of the sumtematice error identified by ID
   * 
   * @param width Line width in pixels
   */
  void SetSumLineWidth(Width_t width){ fSumLine.SetLineWidth(width); }//*MENU*
  /** 
   * Set the fill color of the sumtematice error identified by ID
   * 
   * @param color Fill color 
   */
  void SetSumFillColor(Color_t color){ fSumFill.SetFillColor(color); }//*MENU*
  /** 
   * Set the fill style of the sumtematice error identified by ID
   * 
   * @param style Fill style  
   */
  void SetSumFillStyle(Style_t style){ fSumFill.SetFillStyle(style); }//*MENU*
  /** 
   * Set the draw option for summed errors
   * 
   * @param opt Draw option 
   */
  void SetCommonSumOption(EDrawOption_t opt) { fCommonSumOption = opt; } //*MENU*
  /** 
   * Set the title uses for summed errors
   * 
   * @param title Title
   */
  void SetCommonSumTitle(const char* title) { fCommonSumTitle = title; } //*MENU*
  /** 
   * Set the line color of the sumtematice error identified by ID
   * 
   * @param color Line Color 
   */
  void SetCommonSumLineColor(Color_t color) { fCommonSumLine.SetLineColor(color); } //*MENU*
  /** 
   * Set the line style of the sumtematice error identified by ID
   * 
   * @param style Line style 
   */
  void SetCommonSumLineStyle(Style_t style){ fCommonSumLine.SetLineStyle(style); } //*MENU*
  /** 
   * Set the line width of the sumtematice error identified by ID
   * 
   * @param width Line width in pixels
   */
  void SetCommonSumLineWidth(Width_t width){ fCommonSumLine.SetLineWidth(width); }//*MENU*
  /** 
   * Set the fill color of the sumtematice error identified by ID
   * 
   * @param color Fill color 
   */
  void SetCommonSumFillColor(Color_t color){ fCommonSumFill.SetFillColor(color); }//*MENU*
  /** 
   * Set the fill style of the sumtematice error identified by ID
   * 
   * @param style Fill style  
   */
  void SetCommonSumFillStyle(Style_t style){ fCommonSumFill.SetFillStyle(style); }//*MENU*
  /* @} */
  /**
   * @{
   * @name Key interface 
   */
  /** 
   * Set a key/value pair.  This can be used to fill out fields in a
   * Durham input file for uploading. 
   *
   * @b File @b Keys 
   * - @b cdsId       Should be the number of the CERN Document Server entry 
   * - @b inspireId   Should be the number of the inSpire record
   * - @b durhamId    Will be given on upload 
   * - @b reference   One, or more, list of references, followed by the year
   * - @b laboratory  Main laboratory name 
   * - @b accelerator Facility 
   * - @b detector    The experimental apparatus 
   * - @b title       Title of the paper 
   * - @b author      Last name of first author in all-caps. 
   *
   * @b Data-set @b keys 
   * - @b reackey     The reaction 
   * - @b obskey      The observable 
   * - @b location    Location of figure in paper 
   * - @b dscomment   Comments on data set 
   * - @b qual        One or more entries of qualifiers.  Must contain a 
   *                  value (separated by colons) for each X and Y value
   *
   * At least one @b qual line will be written with the title of the
   * data set as the Y value. 
   *
   * @code 
   GraphSysErr g("foo","bar");
   // Over all header keys 
   g.SetKey("cdsId",       "1483543");
   g.SetKey("inspireId",   "1190545");
   g.SetKey("durhamId",    "6047");
   g.SetKey("reference",   "ARXIV:1210.3615 : 2012");
   g.SetKey("reference",   "CERN-PH-EP-2012-307 : 2012");
   g.SetKey("laboratory",  "CERN");
   g.SetKey("accelerator", "LHC");
   g.SetKey("detector",    "ALICE");
   g.SetKey("title",       "Pseudorapidity density of charged particles "
   "in p-Pb collisions at sqrt{s_{NN}}=5.02");
   g.SetKey("author",      "ABELEV");
   // Data set specific keys 
   g.SetKey("reackey", "P PB --> CHARGED X");
   g.SetKey("obskey"   "DN/DETARAP");
   g.SetKey("qual"     "RE : P PB --> CHARGED X");
   g.SetKey("qual"     "SQRT(S)/NUCLEON in GEV : 5020.0")

   ...
   g.Export();
   * @endcode
   *
   * @param key     Key.
   * @param value   Value.
   * @param replace If true, remove all values of @a key and set new value
   */
  void SetKey(const char* key, const char* value, Bool_t replace=false) //*MENU*
  {
    if (!fMap) { 
      fMap = new TList();
      fMap->SetOwner();
      fMap->SetName("keys");
      // fMap->SetTitle("key/value pairs");
    }
    TString k(key);
    TString v(value);
    TString t =  v.Strip(TString::kBoth);
    v = t;
    if (k.EqualTo("experiment",TString::kIgnoreCase)) {
      TObjArray* l = v.Tokenize("-");
      t = Token(l, 0); TString lab = t.Strip(TString::kBoth);
      t = Token(l, 1); TString acc = t.Strip(TString::kBoth);
      t = Token(l, 2); TString exp = t.Strip(TString::kBoth);
      SetKey("laboratory",  lab.Data(), replace);
      SetKey("accelerator", acc.Data(), replace);
      SetKey("detector",    exp.Data(), replace);
      l->Delete();
    }
    if (k.EqualTo("comment", TString::kIgnoreCase)) { 
      TString l(GetKey("laboratory"));
      TString a(GetKey("accelerator"));
      if (!l.IsNull() && !a.IsNull())
	v.Remove(0, l.Length()+1+a.Length()+1);
      t = v.Strip(TString::kBoth);
      fMap->Add(new TNamed("abstract", t.Data()));
    }
    else {
      if (replace) {
	TObjLink* cur   = fMap->FirstLink();
	TObjLink* last  = fMap->LastLink();
	while (cur) { //  != last) {
	  if (!k.EqualTo(cur->GetObject()->GetName())) {
	    cur = cur->Next();
	    continue;
	  }
	  TObjLink* tmp = cur->Next();
	  fMap->Remove(cur);
	  cur = tmp;
	}
      }
      fMap->Add(new TNamed(k, v));
    }
  }
  /** 
   * Get (first) value of a key 
   * 
   * @param key Key 
   * 
   * @return Value or null
   */
  const char* GetKey(const char* key) const 
  {
    if (!fMap) return 0;
    TObject* o = fMap->FindObject(key);
    if (!o) return 0;
    return o->GetTitle();
  }
  /** 
   * Copy key/value and qualifiers from one graph to this graph. 
   * 
   * @param g       Graph to copy from 
   * @param option  Options specifying what to copy 
   *
   * - a  Copy all key/values and qualifiers
   * - f  Copy file key/values 
   * - h  Copy data set header key/values
   * - q  Copy data set qualifiers 
   * - r  Replace existing values 
   *
   */
  void CopyKeys(const GraphSysErr* g, Option_t* option="fr")
  {
    if (!g) return;
    
    TString opt(option);
    opt.ToLower();

    Bool_t all    = opt.Contains("a");
    Bool_t file   = opt.Contains("f");
    Bool_t header = opt.Contains("h");
    Bool_t qual   = opt.Contains("q");
    Bool_t repl   = opt.Contains("r");
    
    // With the "all" option we just do a plain copy 
    TList*   map = g->fMap;
    TIter    nextKV(map);
    TObject* kv = 0;
    while ((kv = nextKV())) {
      if (!all) {
	TString key = kv->GetName();
	if (file && !(key.EqualTo("reference")  ||
		      key.EqualTo("laboratory") ||
		      key.EqualTo("accelerator")||
		      key.EqualTo("detector") 	||
		      key.EqualTo("abstract") 	||
		      key.EqualTo("author") 	||
		      key.EqualTo("doi") 	||
		      key.EqualTo("inspireId") 	||
		      key.EqualTo("cdsId") 	||
		      key.EqualTo("durhamId") 	||
		      key.EqualTo("title") 	||
		      key.EqualTo("status"))) continue;
	if (header && !(key.EqualTo("location") ||
			key.EqualTo("reackey")	||
			key.EqualTo("obskey")   ||
			key.EqualTo("dscomment"))) continue;	
      }
      SetKey(kv->GetName(), kv->GetTitle(), repl);
    }

    if (!all && !qual) return;
    TList*   quals = g->fQualifiers;
    TIter    nextQ(quals);
    TObject* qv = 0;
    while ((qv = nextQ())) {
      const char* test = GetQualifier(qv->GetName());
      if (!test || test[0] == '\0')
	AddQualifier(qv->GetName(), qv->GetTitle());
    }
  }
  void CopyAttr(const GraphSysErr* f)
  {
    fDataOption = f->fDataOption;
    SetMarkerColor(f->GetMarkerColor());
    SetMarkerStyle(f->GetMarkerStyle());
    SetMarkerSize(f->GetMarkerSize());
    SetLineColor(f->GetLineColor());
    SetLineStyle(f->GetLineStyle());
    SetLineWidth(f->GetLineWidth());
    SetFillColor(f->GetFillColor());
    SetFillStyle(f->GetFillStyle());

    fSumOption = f->fSumOption;
    SetSumTitle(f->GetSumTitle());
    SetSumLineStyle(f->GetSumLineStyle());
    SetSumLineColor(f->GetSumLineColor());
    SetSumLineWidth(f->GetSumLineWidth());
    SetSumFillStyle(f->GetSumFillStyle());
    SetSumFillColor(f->GetSumFillColor());

    fCommonSumOption = f->fCommonSumOption;
    SetCommonSumTitle(f->GetCommonSumTitle());
    SetCommonSumLineStyle(f->GetCommonSumLineStyle());
    SetCommonSumLineColor(f->GetCommonSumLineColor());
    SetCommonSumLineWidth(f->GetCommonSumLineWidth());
    SetCommonSumFillStyle(f->GetCommonSumFillStyle());
    SetCommonSumFillColor(f->GetCommonSumFillColor());

    SetXTitle(f->GetXTitle());
    SetYTitle(f->GetYTitle());
  }

    
    
    
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Qualifiers
   */ 
  /** 
   * Adds a qualifier 
   *
   * @param key     The key 
   * @param value   he value 
   * @param replace If true, replace exsiting value 
   */
  void AddQualifier(const TString& key, const TString& value,
		    Bool_t replace=false)
  {
    if (!fQualifiers) { 
      fQualifiers = new TList();
      fQualifiers->SetOwner();
      fQualifiers->SetName("qualifiers");
    }
    TString val(value);
    if (key.EqualTo(value)) val = "";

    TString k = key.Strip(TString::kBoth, ' ');
    TObject* o = fQualifiers->FindObject(k);
    if (o) {
      Warning("AddQualifier", "Dataset already has qualifier \"%s\"",
	      k.Data());
      if (replace) static_cast<TNamed*>(o)->SetTitle(value);
      return;
    }
    
    fQualifiers->Add(new TNamed(k, value));
  }
  /** 
   * Remove a qualifier 
   * 
   * @param key Which to remove 
   */
  void RemoveQualifier(const TString& key)
  {
    if (!fQualifiers) return;
    TObject* o = fQualifiers->FindObject(key);
    if (!o) return;
    fQualifiers->Remove(o);
    delete o;
  }
  /** 
   * Get qualifier 
   * 
   * @param name Key of qualifier 
   * 
   * @return Value of qualifier 
   */
  const char* GetQualifier(const char* name) const
  {
    if (!fQualifiers) return "";
    TObject* o = fQualifiers->FindObject(name);
    if (!o) return "";
    return o->GetTitle();
  }
  /* @} */
  /** 
   * Find Y value and errors corresponding X 
   * 
   * @param x      X value
   * @param cmn    Include common errors
   * @param stat   Include statistical error 
   * @param quad   Add errors in quadrature
   * @param nosqrt No not take square root of errors 
   * @param y      On return, the y value 
   * @param eyl    On return low error on y
   * @param eyh    On return high error on y
   * @param seyl   On return low error on y
   * @param seyh   On return high error on y
   * 
   * @return false in case the point is out side the range 
   */
  Bool_t FindYandError(Double_t    x,
		       Bool_t      cmn,
		       Bool_t      stat,
		       Bool_t      quad,
		       Bool_t      nosqrt,
		       Double_t&   y,
		       Double_t&   eyl,
		       Double_t&   eyh,
		       Double_t&   seyl,
		       Double_t&   seyh) const
  {
    Int_t i1, i2;
    Int_t ret = FindPoint(x, i1, i2);
    seyl = seyh = 0;
    if (ret < -1)
      // No valid point found 
      return false;
    if (ret >= 0) {
      // Got exact point 
      // Info("", "Get point %d (%f)", ret, x);
      y = GetYandError(ret, cmn, stat, quad, nosqrt, eyl, eyh);
      if (!stat) {
	// If stat error not include, get them here 
	seyl = GetStatErrorDown(ret);
	seyh = GetStatErrorUp(ret);
      }
      return true;
    }
    // Interpolate between points 
    Double_t eyl1,  eyl2,  eyh1,  eyh2;
    Double_t x1    = fData->GetX()[i1];
    Double_t x2    = fData->GetX()[i2];
    Double_t y1    = GetYandError(i1, cmn, stat, quad, nosqrt, eyl1, eyh1);
    Double_t y2    = GetYandError(i2, cmn, stat, quad, nosqrt, eyl2, eyh2);
    Double_t seyl1 = GetStatErrorDown(i1);
    Double_t seyh1 = GetStatErrorUp(i1);
    Double_t seyl2 = GetStatErrorDown(i2);
    Double_t seyh2 = GetStatErrorUp(i2);
    
    // Linear interpolation
    Double_t dx = (x2-x1);
    Double_t ax = (x-x1)/dx;
    y   = y1   + ax * (y2 - y1);
    eyl = eyl1 + ax * (eyl2 - eyl1);
    eyh = eyh1 + ax * (eyh2 - eyh1);
    if (!stat) {
      // if stat errors not included, calculate here
      seyl = seyl1 + ax * (seyl2 - seyl1);
      seyh = seyh1 + ax * (seyh2 - seyh1);
    }
    return true;
  }
  /** 
   * Find Y value and errors corresponding X 
   * 
   * @param x      X value
   * @param cmn    Include common errors
   * @param stat   Include statistical error 
   * @param quad   Add errors in quadrature
   * @param nosqrt No not take square root of errors 
   * @param y      On return, the y value 
   * @param eyl    On return low error on y
   * @param eyh    On return high error on y
   * 
   * @return false in case the point is out side the range 
   */
  Bool_t FindYandError(Double_t    x,
		       Bool_t      cmn,
		       Bool_t      stat,
		       Bool_t      quad,
		       Bool_t      nosqrt,
		       Double_t&   y,
		       Double_t&   eyl,
		       Double_t&   eyh) const
  {
    Double_t seyl;
    Double_t seyh;
    return FindYandError(x, cmn, stat, quad, nosqrt, y, eyl, eyh, seyl, seyh);
  }

  //__________________________________________________________________
  /** 
   * Combining measurements 
   *
   * From http://www.slac.stanford.edu/~barlow/java/
   */
  struct Combiner : public TObject
  {
    typedef Double_t (*Wrapper_t)(Double_t*,Double_t*);
    /**
     * An experimental observation 
     */
    struct Observation : public TObject 
    {
      /** Value @f$ x_i@f$ */
      Double_t fX;
      /** Low error @f$ \sigma_i^-@f$ */
      Double_t fEl;
      /** High error @f$ \sigma_i^+@f$ */
      Double_t fEh;
      /** 
       * Create a single obersvation 
       * 
       * @param x  Value @f$ x_i@f$
       * @param el Low error @f$ \sigma_i^-@f$ on @f$ x_i@f$ 
       * @param eh High error @f$ \sigma_i^+@f$ on @f$ x_i@f$ 
       */
      Observation(Double_t x=0, Double_t el=0, Double_t eh=0)
	: fX(x), fEl(el), fEh(eh)
      {    
      }
      /** 
       * Virtual destructor 
       */
      virtual ~Observation() {}
      /** 
       * Calculate 
       * 
       * @f[ 
       *   s_i = \sigma_i^+ \sigma_i^- / (\sigma_i^+ + \sigma_i^-)
       * @f]
       * 
       * 
       * @return @f$ s_i@f$
       */
      Double_t S() const
      {
	if (fEl * fEh == 0) return 0;
	return 2 * fEl * fEh / (fEl + fEh);
      }
      /** 
       * Calculate 
       * 
       * @f[ 
       *   s_i' = (\sigma_i^+ - \sigma_i^-) / (\sigma_i^+ + \sigma_i^-)
       * @f]
       * 
       * 
       * @return @f$ s_i'@f$
       */
      Double_t Sprime() const
      {
	if (TMath::Abs(fEh + fEl) < 1e-9) return 1;	
	return (fEh - fEl) / (fEh + fEl);
      }
      Double_t Svar(Double_t guess=0) const
      {
	return  S() + Sprime() * (guess - fX);
      }
      /** 
       * Calculate 
       * 
       * @f[ 
       *   V_i = \sigma_i^+ \sigma_i^-
       * @f]
       * 
       * 
       * @return @f$ V_i@f$
       */
      Double_t V() const  { return fEl * fEh; }
      /** 
       * Calculate 
       * 
       * @f[ 
       *   V_i' = \sigma_i^+ - \sigma_i^-
       * @f]
       * 
       * 
       * @return @f$ V_i'@f$
       */
      Double_t Vprime() const { return fEh - fEl;  }
      /** 
       * Lower bound 
       * 
       * @return @f$ x_i - 3\sigma_i^+@f$ 
       */
      virtual Double_t Low() const { return fX - 3 * fEl; }
      /** 
       * Upper bound 
       * 
       * @return @f$ x_i + 3\sigma_i^+@f$ 
       */
      virtual Double_t High() const { return fX + 3 * fEh; }
      /** 
       * Print to standard output 
       *
       * @param option Not used 
       */
      virtual void Print(Option_t* option="") const
      {
	Printf("%10.8f  -%10.8f +%10.8f", fX, fEl, fEh);
      }
      ClassDef(Observation,1); // An experimental observation 
    };
  
    /**
     * The final result 
     */
    struct Result : public Observation
    {
      /** The final @f$\chi^2 @f$ */
      Double_t fChi2;
      /** Lower bound to use */
      Double_t fLower;
      /** Upper bound to use */
      Double_t fUpper;
      /** 
       * The final result 
       * 
       * @param x     Best estimate of @f$ x@f$ 
       * @param el    Best estimate of @f$ \sigma^-@f$ 
       * @param eh    Best estimate of @f$ \sigma^+@f$ 
       * @param chi2  @f$\chi^2@f$ of best estimate 
       * @param low   The lower bound of the source data
       * @param high  The upper bound of the source data
       */
      Result(Double_t x=0, Double_t el=0, Double_t eh=0,
	     Double_t chi2=0, Double_t low=0, Double_t high=0)
	: Observation(x, el, eh),
	  fChi2(chi2),
	  fLower(low),
	  fUpper(high)
      {}
      /** 
       * Lower bound 
       * 
       * @return @f$ x_i - 3\sigma_i^+@f$ 
       */
      Double_t Low() const { return fLower; }
      /** 
       * Upper bound 
       * 
       * @return @f$ x_i + 3\sigma_i^+@f$ 
       */
      Double_t High() const { return fUpper; }
      /** 
       * Print to standard output 
       *
       * @param option Not used 
       */
      virtual void Print(Option_t* option="") const
      {
	Printf("%10.8f  -%10.8f +%10.8f  %5.2f", fX, fEl, fEh, fChi2);
      }
      ClassDef(Result,1); // Final result 
    };
    TClonesArray fData;   // Container of observations
    Result*      fResult; // Result of combining

    /** 
     * Constructor 
     */
    Combiner()
      : fData("GraphSysErr::Combiner::Observation"), fResult(0)
    {
      fData.SetOwner();
    }
    /** 
     * Virtual destructor 
     */
    virtual ~Combiner() { fData.Clear(); if (fResult) delete fResult; }
    /** 
     * @{ 
     * @name The data store 
     */
    /** 
     * Clear the internal data 
     * 
     * @param option Not used
     */
    void Clear(Option_t* option="")
    {
      fData.Clear();
      if (fResult) delete fResult;
      fResult = 0;
    }
    /** 
     * Add an obervation 
     * 
     * @param r Observation
     */
    void Add(const Observation& r)
    {
      Add(r.fX,r.fEl,r.fEh);
    }
    /** 
     * Add an observation 
     * 
     * @param x  @f$ x_i@f$ 
     * @param el @f$ \sigma_i^-@f$  
     * @param eh @f$ \sigma_i^+@f$
     */
    void Add(Double_t x, Double_t el, Double_t eh)
    {
      // if (TMath::Abs(x)  < 1e-10 &&
      //     TMath::Abs(el) < 1e-10 &&
      //     TMath::Abs(eh) < 1e-10) return;
      new (fData[fData.GetEntries()]) Observation(x,el,eh);
    }
    /** 
     * Print content of the list 
     * 
     * @param option not used
     */
    void Print(Option_t* option="") const
    {
      TIter next(&fData);
      Observation* r = 0;
      while ((r = static_cast<Observation*>(next())))
	r->Print(option);
    }
    /* @} */
  
    /** 
     * Calculate the weight 
     *
     * @param r Observation
     * 
     * @return @f$ W@f$ 
     */
    virtual Double_t W(const Observation& r) const = 0;
    /** 
     * Calculate the weight based on a guess of best @f$ x'@f$ 
     * 
     * @param guess Current guess @f$ x'@f$ 
     * @param r  	  Observation
     * 
     * @return @f$ W(x')@f$ 
     */
    virtual Double_t StepW(Double_t guess, const Observation& r) const = 0;
    /** 
     * Calculate the bias. 
     * 
     * @return @f$\delta(x')@f$ 
     */
    virtual Double_t StepOffset(Double_t guess, const Observation& r) const = 0;
    /** 
     * Calculate the contribution variance to the @f$\chi^2@f$ with
     * the guess @f$x'@f$. 
     *
     * @return @f$ v(x')@f$ 
     */
    virtual Double_t VarTerm(Double_t guess, const Observation& r) const = 0;
    /** 
     * Calculate the contribution variance to the @f$\chi^2@f$ with
     * the guess @f$ x'@f$. 
     *
     * @f[
     *   t_i(x') = (x' - x_i)^2 / v_i(x')
     * @f]
     *
     * where @f$ v_i(x')@f$ is the term variance
     *
     * @param guess @f$ x'@f$  
     * @param r     Obersvation 
     *
     * @return @f$ t(x')@f$ 
     */
    Double_t ChiTerm(Double_t guess, const Observation& r) const
    {
      Double_t var = VarTerm(guess, r);
      if (var <= 0) return -1000;

      return TMath::Power(guess - r.fX, 2) / var;
    }
    /** 
     * Calculate the @f$ \chi^2(x')@f$ where @f$ x'@f$ is current guess
     * at the observation.  
     * 
     * @param guess   Current guess @f$ x'@f$ 
     * @param chi2    Optional old @f$ \chi^2@f$ from best @f$ x@f$ value 
     * 
     * @return @f$ \chi^2(x')@f$
     */
    Double_t F(Double_t          guess,
	       Double_t          chi2) const
    {
      if (fData.GetEntries() == 1) return 1;
      Double_t s = -chi2;

      TIter next(&fData);
      Observation* r = 0;
      while ((r = static_cast<Observation*>(next())))
	s += ChiTerm(guess, *r);
      return s;
    }
    /** 
     * Try to find best error 
     * 
     * @param nIter Number of iterations 
     * @param sign  Direction (-1 is low, +1 is high)
     * @param best  Current best @f$ x@f$ value 
     * @param chi2  @f$ \chi^2@f$ of current best @f$ x@f$ value 
     * @param s     Summed weights in the direction 
     * 
     * @return The error in the chosen direction
     */
    Double_t E(UShort_t        nIter,
	       Int_t           sign,
	       Double_t        best,
	       Double_t        chi2,
	       Double_t        s)
    {
      if (fData.GetEntries() == 1) {
	Observation* o = static_cast<Observation*>(fData[0]);
	return (sign < 0 ? o->fEl : o->fEh);
      }

      // Step size 
      Double_t delta = 0.1 * sign * s;

      // Iterations 
      for (UShort_t i = 0; i < nIter; i++) {
	// Calculate chi^2 with current guess 
	Double_t got = F(best + sign * s, chi2);

	if (TMath::Abs(got-1) < 1e-7)
	  // We're close to 1 so get out
	  break;

	// The next guess' chi^2 value e
	Double_t guess = F(best + sign * s + delta, chi2);

	// Where to search next 
	if ((got - 1) * (guess - 1) > 0) {
	  if ((got - 1) / (guess - 1) < 1)
	    delta = -delta;
	  else
	    s += sign * delta;
	  continue;
	}

	// Update error value and decrease step size 
	s     += sign * delta * (1 - got) / (guess - got);
	delta /= 2;
      }
      return s;
    }
    /** 
     * Find best estimate of @f$ x@f$ 
     * 
     * @param nIter   Number of iterations 
     * @param lowest  Lower bound 
     * @param highest Upper bound 
     * 
     * @return  @f$ x@f$ 
     */
    Double_t X(UShort_t      nIter,
	       Double_t      lowest,
	       Double_t      highest)
    {
      // Starting values
      if (fData.GetEntries() == 1)
	return static_cast<Observation*>(fData[0])->fX;
      Double_t x    = (highest+lowest)/2;
      Double_t oldX = -1e33;
      // Do the iterations 
      for (UShort_t i = 0; i < nIter; i++) {
	Double_t sum    = 0;
	Double_t sumw   = 0;
	Double_t offset = 0;

	// Loop over observations
	TIter next(&fData);
	Observation* r = 0;
	while ((r = static_cast<Observation*>(next()))) {
	  Double_t w =  StepW(x, *r);
	  offset     += StepOffset(x,  *r);
	  sum  += r->fX * w;
	  sumw += w;
	}
	x = (TMath::Abs(sumw) > 1e-9) ? (sum - offset) / sumw : 0;

	if (TMath::Abs(x - oldX) < (highest-lowest) * 1e-9) break;
	oldX = x;
      }
      return x;
    }
    /** 
     * Do the calculation 
     * 
     * @param nIter How many iterations to do. 
     * 
     * @return The best estimate of @f$ x@f$ and associated errors 
     */
    Result* Calculate(UShort_t nIter=50)
    {
      Double_t lowest  = +1e33;
      Double_t highest = -1e33;
      Double_t sumLow  = 0;
      Double_t sumHigh = 0;

      // Find boundaries and sum weights
      TIter next(&fData);
      Observation* r = 0;
      while ((r = static_cast<Observation*>(next()))) {
	lowest  = TMath::Min(r->Low(),  lowest);
	highest = TMath::Max(r->High(), highest);
	if (r->fEl > 1e-10)  sumLow  += 1./TMath::Power(r->fEl, 2);
	if (r->fEh > 1e-10)  sumHigh += 1./TMath::Power(r->fEh, 2);
      }
      // Summed weights 
      Double_t sLow  = sumLow  > 1e-10 ? 1. / TMath::Sqrt(sumLow) : 0;
      Double_t sHigh = sumHigh > 1e-10 ? 1. / TMath::Sqrt(sumHigh) : 0;
    
      // Now do the calculations
      Double_t bestX    = X(nIter, lowest, highest);
      Double_t bestChi2 = F(bestX, 0);
      Double_t bestLow  = TMath::Abs(E(nIter,-1,bestX,bestChi2,sLow));
      Double_t bestHigh = TMath::Abs(E(nIter,+1,bestX,bestChi2,sHigh));

      fResult = new Result(bestX, bestLow, bestHigh, bestChi2, lowest, highest);
      return fResult;
    }
    /** 
     * Return function pointer to wrapper
     * 
     * @return Function pointer
     */
    virtual Wrapper_t Wrapper() const = 0;
    
    /** 
     * Make a function that represents to Log-likehood for a given
     * observation.
     * 
     * @param r Observation 
     * @param j Serial number 
     * 
     * @return Pointer to newly allocated function object
     */
    TF1* MakeF(const Observation& r, Int_t j) const
    {
      TF1* f = new TF1(Form("f%02d", j), Wrapper(), r.Low(), r.High(), 3);
      f->SetParNames("x", "#sigma^{-}", "#sigma^{+}");
      f->SetParameters(r.fX, r.fEl, r.fEh);
      f->SetLineStyle((j%3)+2);
      f->SetLineColor(kBlack);
      f->SetLineWidth(2);
      return f;
    }
    /** 
     * Make a line that represents the best found errors
     * 
     * @param f Log-likelyhood function to make it from 
     * 
     * @return 
     */
    TLine* MakeL(TF1* f) const
    {
      Double_t m = f->GetParameter(0);
      TLine* l = new TLine(m-f->GetParameter(1), 1,
			   m+f->GetParameter(2), 1);
      l->SetLineColor(f->GetLineColor());
      l->SetLineStyle(f->GetLineStyle());
      l->SetLineWidth(f->GetLineWidth());
      return l;
    }
    void Draw(Option_t* option="")
    {
      if (!fResult) return;
      
      TList fs; fs.SetOwner(false);
      Int_t j = 0;
      TIter next(&fData);
      Observation* r = 0;
      while ((r = static_cast<Observation*>(next()))) {
	TF1* f = MakeF(*r, j);
	f->SetRange(fResult->Low(), fResult->High());
	fs.Add(f, j == 0 ? "" : "same");
	fs.Add(MakeL(f));
	j++;
      }
      TF1* fr = MakeF(*fResult, j);
      fr->SetLineColor(kRed+2);
      fr->SetLineStyle(1);
      fs.Add(fr);
      fs.Add(MakeL(fr));
      TIter nextD(&fs);
      TObject* o = 0;
      j = 0;
      while((o = nextD())) { o->Draw(j == 0 ? "" : "same"); j++; }
      // fs.Draw();
      static_cast<TF1*>(fs.First())->GetHistogram()
	->GetYaxis()->SetRangeUser(-.1, 5.1);
    }
    ClassDef(Combiner,1); // Combine systematics 
  };


  /**
   * A combiner that uses a linear @f$\sigma@f$ approximation 
   */
  struct LinearSigmaCombiner : public Combiner
  {
    /** 
     * Calculate the weight 
     * 
     * @f[ 
     *   w = 1/2 (s + x s')^3 / s
     * @f]
     *
     * @param r Observation
     * 
     * @return @f$ W@f$ 
     */
    Double_t W(const Observation& r) const
    {
      // Double_t s = r.S();
      // if (s < 1e-10) return 0;
      // return .5 * TMath::Power(r.Svar(), 3) / s;
      Double_t w = StepW(0,r);
      return 1/w;
    }
    /** 
     * Calculate the weight based on a guess of best @f$ x'@f$ 
     * 
     * @f[ 
     *   w(x') = s / [s + s' (x' - x)]^3
     * @f]
     *
     * @param guess Current guess @f$ x'@f$ 
     * @param r	  Observation
     * 
     * @return @f$ W(x')@f$ 
     */
    Double_t StepW(Double_t guess, const Observation& r) const
    {
      Double_t s   = r.S();
      Double_t t   = r.Svar(guess);
      Double_t ret = s / TMath::Power(t,3);
      return ret;
    }
    /** 
     * Calculate the bias. 
     * 
     * @return 0
     */
    Double_t StepOffset(Double_t, const Observation&) const
    {
      return 0;
    }
    /** 
     * Calculate the contribution variance to the @f$\chi^2@f$ with
     * the guess @f$x'@f$. 
     *
     * @f[
     *   v(x') = [s + s' (x' - x)]^2
     * @f]
     * 
     * @param guess Current guess @f$ x'@f$ 
     * @param r	  Observation
     * 
     * @return @f$ v(x')@f$ 
     */
    Double_t VarTerm(Double_t guess, const Observation& r) const
    {
      return TMath::Power(r.Svar(guess),2);
    }
    /** 
     * Return the likely-hood function value at @f$ x'@f$:
     *
     * @f[
     *   L(x') = \left[(x'-x) / (s + s'(x'-x))\right]^2
     * @f] 
     * 
     * where 
     * @f[
     *   s = 2\sigma^+\sigma^-/(\sigma^++\sigma^-)
     * @f]
     * @f[
     *   s' = (\sigma^+-\sigma^-)/(\sigma^++\sigma^-)
     * @f]
     * 
     * @param guess @f$ x'@f$ 
     * @param x     @f$ x@f$  
     * @param el    @f$ \sigma^-@f$  
     * @param eh    @f$ \sigma^+@f$  
     * 
     * @return 
     */
    static Double_t L(Double_t guess, Double_t x, Double_t el, Double_t eh)
    {
      Observation tmp(x,el,eh);
      Double_t d  = (guess-x);
      Double_t t  = tmp.Svar(guess); // (s+sp*d);
      if (TMath::Abs(d) < 1e-10) return 0;
      // if (TMath::Abs(t) < 1e-10) return DBL_MAX;
      Double_t ret = TMath::Power(d/t, 2);
      return ret;
    }
    /** 
     * Wrap likely-hood function for ROOT 
     * 
     * @param xp Pointer to independent variables 
     * @param pp Pointer to parameters 
     * 
     * @return Likely-hood function evaluate at @a xp
     */
    static Double_t WrapL(Double_t* xp, Double_t* pp)
    {
      return L(xp[0], pp[0], pp[1], pp[2]);
    }
    /** 
     * Return function pointer to wrapper
     * 
     * @return Function pointer
     */
    Wrapper_t Wrapper() const { return WrapL; }
    ClassDef(LinearSigmaCombiner,1); // Sigma combiner 
  };
  /**
   * A combiner that uses a linear variance approximation 
   */
  struct LinearVarianceCombiner : public Combiner
  {
    /** 
     * Calculate the weight 
     * 
     * @f[
     *   w = (V + x V')^2 / (2 V + x V')
     * @f] 
     *
     * @param r Observation
     * 
     * @return @f$ W@f$ 
     */
    Double_t W(const Observation& r) const
    {
      Double_t v  = r.V();
      Double_t vp = r.Vprime();
      return TMath::Power(v + r.fX * vp, 2) / (2 * v + r.fX * vp);
    }
    /** 
     * Calculate the weight based on a guess of best @f$ x'@f$ 
     * 
     * @f[
     *   W(x') = v / [V + V' (x' - x)]^2
     * @f] 
     *
     * @param guess Current guess @f$ x'@f$ 
     * @param r  	  Observation
     * 
     * @return @f$ W(x')@f$ 
     */
    Double_t StepW(Double_t guess, const Observation& r) const
    {
      Double_t v = r.V();
      return v / TMath::Power(v+r.Vprime()*(guess - r.fX), 2);
    }
    /** 
     * Calculate the bias. 
     * 
     * @f[
     *   \delta(x') = 1/2 V' [(x'-x) / (V + V'(x' - x))]^2
     * @f] 
     * 
     * @param guess Current guess @f$ x'@f$ 
     * @param r	  Observation
     * 
     * @return @f$\delta(x')@f$ 
     */
    Double_t StepOffset(Double_t guess, const Observation& r) const
    {
      Double_t vp = r.Vprime();
      return 0.5 * vp * TMath::Power((guess-r.fX)/(r.V()+vp*(guess-r.fX)),2);
    }
    /** 
     * Calculate the contribution variance to the @f$\chi^2@f$ with
     * the guess @f$x'@f$. 
     *
     * @f[ 
     *   V(x') = V + V' (x' - x)
     * @f] 
     * 
     * @param guess Current guess @f$ x'@f$ 
     * @param r	  Observation
     * 
     * @return @f$ v(x')@f$ 
     */
    Double_t VarTerm(Double_t guess, const Observation& r) const
    {
      return r.V() + r.Vprime() * (guess - r.fX);
    }
    /** 
     * Return the likely-hood function value at @f$ x'@f$:
     *
     * @f[
     *   L(x') = (x'-x)^2 / (V + V'(x'-x))
     * @f] 
     * 
     * where 
     * @f[
     *   V = \sigma^+\sigma^-\quad V' = \sigma^+-\sigma^-
     * @f]
     * 
     * @param guess @f$ x'@f$ 
     * @param x     @f$ x@f$  
     * @param el    @f$ \sigma^-@f$  
     * @param eh    @f$ \sigma^+@f$  
     * 
     * @return 
     */
    static Double_t L(Double_t guess, Double_t x, Double_t el, Double_t eh)
    {
      Double_t d  = (guess-x);
      Double_t v  = eh * el;
      Double_t vp = eh - el;
      return TMath::Power(d,2) / (v+vp*d);
    }
    /** 
     * Wrap likely-hood function for ROOT 
     * 
     * @param xp Pointer to independent variables 
     * @param pp Pointer to parameters 
     * 
     * @return Likely-hood function evaluate at @a xp
     */
    static Double_t WrapL(Double_t* xp, Double_t* pp)
    {
      return L(xp[0], pp[0], pp[1], pp[2]);
    }
    /** 
     * Return function pointer to wrapper
     * 
     * @return Function pointer
     */
    Wrapper_t Wrapper() const { return WrapL;  }
    ClassDef(LinearVarianceCombiner,1); // Variance combiner 
  };
  
  //__________________________________________________________________
  /** 
   * Base class to hold systematic errors
   */
  struct Holder : public TNamed, TAttLine, TAttFill 
  {
    /** Containing class is a friemd */
    friend struct GraphSysErr;
    /** 
     * CTOR
     */
    Holder() 
      : TNamed(), 
	TAttLine(),
	TAttFill(),
	fRelative(true),
	fOption(kRect)
    {}
    /**
     * DTOR
     */
    virtual ~Holder() {}
    void CopyAttr(Holder* h)
    {
      SetLineColor(h->GetLineColor());
      SetLineStyle(h->GetLineStyle());
      SetLineWidth(h->GetLineWidth());
      SetFillColor(h->GetFillColor());
      SetFillStyle(h->GetFillStyle());
      fOption = h->fOption;
    }
  protected:
    /** 
     * CTOR with name and title
     * 
     * @param name   Name 
     * @param title  Title 
     * @param rel    Relative or absolue 
     * @param option Draw Option 
     * @param id     Identifier 
     */
    Holder(const char* name, const char* title, Bool_t rel,
	   UInt_t option, UInt_t id)
      : TNamed(name, title), 
	TAttLine(0,0,0),
	TAttFill(0,0),
	fRelative(rel),
	fOption(option)
    {
      SetUniqueID(id);
    }
    /** 
     * Copy constructorr 
     * 
     * @param other Object to copy from
     */
    Holder(const Holder& other) 
      : TNamed(other), 
	TAttLine(other),
	TAttFill(other),
	fRelative(other.fRelative),
	fOption(other.fOption)
    {
      SetUniqueID(other.GetUniqueID());
    }
    /** 
     * Assignment operator
     * 
     * @param other Object to assign from
     *
     * @return reference to this object
     */
    Holder& operator=(const Holder& other) 
    {
      if (&other == this) return *this;

      other.TNamed::Copy(*this);
      other.TAttLine::Copy(*this);
      other.TAttFill::Copy(*this);

      fRelative = other.fRelative;
      fOption   = other.fOption;
    
      SetUniqueID(other.GetUniqueID());

      return *this;
    }
    
    /** 
     * Create new graph with stacked errors
     *  
     * @param g          Previous errors
     * @param ignoreErr  If true, ignore previous errors
     * @param quad       If true, add in quadrature 
     * 
     * @return Newly allocated graph 
     */
    virtual Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) const = 0;
    /** 
     * Sum errors at point. Point i of g is updated
     * 
     * @param g            Where to sum
     * @param i            Point
     * @param ignoreErr    If true, ignore exusisting errros
     * @param quad         Add in quadrature
     * @param opt          Option
     */
    virtual void SumError(Graph* g, Int_t i, Bool_t ignoreErr, 
			  Bool_t quad, UInt_t opt) const = 0;
    /** 
     * @return Get the Option
     */
    virtual UInt_t GetDOption() const { return fOption; }
    /** 
     * Set the draw option
     * 
     * @param opt Option
     */
    virtual void SetDOption(EDrawOption_t opt) { fOption = opt; }
    /** 
     * Check if this is a relative error
     * 
     * @return true if declared relative 
     */
    virtual Bool_t IsRelative() const { return fRelative; }

    virtual void Print(Option_t* option="") const
    {
      gROOT->IndentLevel();
      Printf("%s/%s %s(ID # %d%s)", GetName(), GetTitle(),
	     IsRelative() ? "PCT " : "",  GetUniqueID(),
	     TestBit(kUsedBit) ? " used" : "");
      TString opt(option);
      if (!opt.Contains("ATTR", TString::kIgnoreCase)) return;
      gROOT->IndentLevel();
      Printf(" [option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d]",
	     fOption, GetLineColor(), GetLineStyle(), GetLineWidth(),
	     GetFillColor(), GetFillStyle());
      
    }
    virtual void ls(Option_t* option) const 
    {
      Print(option);
    }
  protected:
    UShort_t XMode(Int_t opt=-1) const
    {
      UShort_t xMode  = 0;
      if (opt < 0) opt = fOption;
      switch (opt) {
      case kNormal: case kNoTick: case kFill: case kCurve: xMode = 0; break;
      case kArrow:  case kHat:    case kBar:               xMode = 1; break;
      case kNone:                                          xMode = 0; break;
      case kRect:   case kBox:                             xMode = 3; break;
      }
      return xMode;
    }
    /** 
     * Do add errors 
     * 
     * @param xMode     X-mode
     * @param curExl    Currently summed/stacked X low error
     * @param curExh    Currently summed/stacked X high error
     * @param curEyl    Currently summed/stacked Y low error 
     * @param curEyh 	Currently summed/stacked Y high error
     * @param ignoreErr If true, ignore errors on currently stack errors
     * @param quad      If true, add in quadrature 
     * @param sqOld     If false, assume current errors are squared already 
     * @param exl       Input: this sources X low error, Output: new value
     * @param exh       Input: this sources X high error, Output: new value
     * @param eyl       Input: this sources Y low error, Output: new value 
     * @param eyh 	Input: this sources Y high error, Output: new value
     */
    void DoAdd(UShort_t  xMode,
	       Double_t  curExl,
	       Double_t  curExh,
	       Double_t  curEyl,
	       Double_t  curEyh,
	       Bool_t    ignoreErr,
	       Bool_t    quad, 
	       Bool_t    sqOld,
	       Double_t& exl, 
	       Double_t& exh, 
	       Double_t& eyl, 
	       Double_t& eyh) const
    {
      // Double_t oexl = exl;
      // Double_t oexh = exh;
      if (quad) { 
	eyl    *= eyl;
	eyh    *= eyh;
	if (sqOld) {
	  curEyl *= curEyl;
	  curEyh *= curEyh;
	}
      }
	
      if (!ignoreErr) {
	if (kVerbose & kDraw)
	  Info("", "old: +%f/-%f this: +%f/-%f -> +%f/-%f",
	       curEyl, curEyh, eyl, eyh, eyl + curEyl, eyh + curEyh);
	exl  = TMath::Max(exl, curExl);
	exh  = TMath::Max(exh, curExh);
	eyl  += curEyl;
	eyh  += curEyh;
      }

      switch (xMode) {
      case 0: break;                // kNormal, kNoTick, kFill, kCurve, kNone
      case 1: exl = exh = 0; break; // kArrow, kHat, kBar
      case 2:                       // kRect, kBox
	if (exl <= 0) exl = gStyle->GetErrorX()/2;
	if (exh <= 0) exh = gStyle->GetErrorX()/2;
	break;
      }
    }
    
    /** 
     * Set attributes
     * 
     * @param g on graph
     */
    void SetAttributes(Graph* g) const
    {
      if (!g) return;
      this->TAttLine::Copy(*g);
      this->TAttFill::Copy(*g);
      g->SetMarkerColor(0);
      g->SetMarkerStyle(0);
      g->SetMarkerSize(0);
    }
    /** Relative error flag */
    Bool_t  fRelative;
    /** Options */
    UInt_t fOption;
    ClassDef(Holder,3);
  };
  //__________________________________________________________________
  /** 
   * A holder for Point-to-Point systematic errors
   */
  struct HolderP2P : public Holder
  {
    /** Containing class is a friemd */
    friend struct GraphSysErr;
    /** 
     * CTOR
     */
    HolderP2P()
      : Holder(),
	fGraph(0)
    {}
  protected:
    /** 
     * CTOR with name and title
     * 
     * @param name   Name 
     * @param title  Title 
     * @param rel    Relative or absolue 
     * @param opt    Draw Option 
     * @param id     Identifier 
     */
    HolderP2P(const char* name, const char* title, Bool_t rel,
	      UInt_t  opt, UInt_t id)
      : Holder(name, title, rel, opt, id),
	fGraph(0)
    {
      fGraph = new Graph();
      fGraph->SetName(name);
      fGraph->SetTitle(title);
      SetAttributes(fGraph);
    }
    /** 
     * Copy CTOR
     * 
     * @param other Object ot copy from
     */
    HolderP2P(const HolderP2P& other)
      : Holder(other),
	fGraph(0)
    {
      if (other.fGraph) 
	fGraph = static_cast<Graph*>(other.fGraph->Clone());

      if (fGraph) SetAttributes(fGraph);
    }
    /** 
     * Assignement operator
     * 
     * @param other Object to assign from
     * 
     * @return 
     */
    HolderP2P& operator=(const HolderP2P& other)
    {
      if (&other == this) return *this;
      
      this->Holder::operator=(other);
      if (fGraph) delete fGraph;
      fGraph = 0;
      if (other.fGraph) 
	fGraph = static_cast<Graph*>(other.fGraph->Clone());
      if (fGraph) SetAttributes(fGraph);
      
      return *this;
    }

    /** 
     * @{ 
     * @name Setting errors 
     */
    /** 
     * Set errors at point 
     * 
     * @param point  Point
     * @param g      Graph
     * @param ex     Symmetric error along X
     * @param ey     Symmetric error along Y
     */
    void Set(Int_t point, Graph* g, Double_t ex, Double_t ey) 
    {
      Set(point, g, ex, ex, ey, ey);
    }
    /** 
     * Set errors at point 
     * 
     * @param point Point
     * @param g     Graph
     * @param ex1   Low erros along X
     * @param ex2   High errors along X
     * @param ey1   Low erros along Y
     * @param ey2   High errors along Y
     */
    void Set(Int_t point, Graph* g, Double_t ex1, Double_t ex2, 
	     Double_t ey1, Double_t ey2) 
    {
      if (!fGraph) return;
      if (!g) return;
      if (point >= g->GetN()) return;
      Double_t x = g->GetX()[point];
      Double_t y = g->GetY()[point];
      fGraph->SetPoint(point, x, y);
      if (fRelative) { ey1 *= y; ey2 *= y; }
      fGraph->SetPointError(point, ex1, ex2, ey1, ey2);
    }
    /* @} */

    /** 
     * @{ 
     * @name Get information 
     */
    /** 
     * Get symmetric errors along X at point 
     * 
     * @param point Point
     * 
     * @return Symmetric errors along X at point
     */
    Double_t GetX(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorX(point);
    }
    /** 
     * Get errors to the left along X at point
     * 
     * @param point Point
     * 
     * @return Errors to the left along X at point
     */
    Double_t GetXLeft(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorXlow(point);
    }
    /** 
     * Get errors to the right along X at point
     * 
     * @param point Point
     * 
     * @return Errors to the right along X at point
     */
    Double_t GetXRight(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorXhigh(point);
    }
    /** 
     * Get symmetric errors along Y at point 
     * 
     * @param point Point
     * 
     * @return Symmetric (absolute) errors along Y at point
     */
    Double_t GetY(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorY(point);
    }
    /** 
     * Get the errors downwards along Y between points @a i1 and @a i2 
     * 
     * @param i1 Left point
     * @param i2 Right point
     * @param ax Relative distance between the two points 
     * 
     * @return (Absolute) Errors downward along Y at point
     */
    Double_t GetYDown(Int_t i1, Int_t i2, Double_t ax) const 
    {
      if (!fGraph) return 0;
      if (i1 == i2) return GetYDown(i1);
      Double_t e1 = GetYDown(i1);
      Double_t e2 = GetYDown(i2);
      return e1 + ax * (e2-e1); // Linear interpolation
    }
    /** 
     * Get errors downward along Y at point
     * 
     * @param point Point
     * 
     * @return (Absolute) Errors downward along Y at point
     */
    Double_t GetYDown(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorYlow(point);
    }
    /** 
     * Get the errors upwards along Y between points @a i1 and @a i2 
     * 
     * @param i1 Left point
     * @param i2 Right point
     * @param ax Relative distance between the two points 
     * 
     * @return (Absolute) Errors upward along Y at point
     */
    Double_t GetYUp(Int_t i1, Int_t i2, Double_t ax) const 
    {
      if (!fGraph) return 0;
      if (i1 == i2) return GetYUp(i1);
      Double_t e1 = GetYUp(i1);
      Double_t e2 = GetYUp(i2);
      return e1 + ax * (e2-e1); // Linear interpolation
    }
    /** 
     * Get errors upward along Y at point
     * 
     * @param point Point
     * 
     * @return (Absolute) Errors upward along Y at point
     */
    Double_t GetYUp(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorYhigh(point);
    }
    /* @} */

    /** 
     * @{ 
     * @name Sum, add, stack 
     */
    /** 
     * Add errors together at point
     * 
     * @param i            Point
     * @param xMode        X-mode 
     * @param quad         If true, add in quadrature
     * @param ignoreErr    If true, ignore errors on g
     * @param sqOld        If true and quad true, square old
     * @param exl          Input current, Output: the left-hand X errors
     * @param exh          Input current, Output: the right-hand X errors
     * @param eyl          Input current, Output: the downward Y errors
     * @param eyh          Input current, Output: the upward Y errors 
     */
    void AddError(Int_t     i,
		  UShort_t  xMode,
		  Bool_t    ignoreErr,
		  Bool_t    quad, 
		  Bool_t    sqOld,
		  Double_t& exl, 
		  Double_t& exh, 
		  Double_t& eyl, 
		  Double_t& eyh) const
    {
      Double_t exlO = fGraph->GetErrorXlow(i);
      Double_t exhO = fGraph->GetErrorXhigh(i);
      Double_t eylO = fGraph->GetErrorYlow(i);
      Double_t eyhO = fGraph->GetErrorYhigh(i);
      Double_t oyl  = eylO;
      Double_t oyh  = eyhO;
      if (exl <= 0) exlO = exl;
      if (exh <= 0) exhO = exh;

      DoAdd(xMode,
	    exl, exh, eyl, eyh,
	    ignoreErr, quad, sqOld,
	    exlO, exhO, eylO, eyhO);
      if (kVerbose & kDraw)
	Info(GetTitle(), "quad: %s this: +%f/-%f, old: +%f/-%f -> +%f/-%f",
	     (quad ? "true" : "false"), oyl, oyh, eyl, eyh, eylO, eyhO);

      exl = exlO;
      exh = exhO;
      eyl = eylO;
      eyh = eyhO;
    }
    /** 
     * Stack up point errors 
     * 
     * @param i         Point number 
     * @param xMode     X-Mode
     * @param ignoreErr If true, ignore current errors
     * @param quad      If true, add in quadrature 
     * @param exl       Input: current  Output: New value 
     * @param exh       Input: current  Output: New value 
     * @param eyl       Input: current  Output: New value 
     * @param eyh       Input: current  Output: New value 
     */
    void StackPointError(Int_t     i,
			 UShort_t  xMode,
			 Bool_t    ignoreErr, 
			 Bool_t    quad,
			 Double_t& exl,
			 Double_t& exh, 
			 Double_t& eyl,
			 Double_t& eyh) const
    {
      Double_t oldEyl = eyl;
      Double_t oldEyh = eyh;
      AddError(i, xMode, ignoreErr, quad, true, exl, exh, eyl, eyh);
      if (kVerbose & kDraw)
	Info(GetTitle(), "old= +%f/-%f -> +%f/-%f", oldEyl, oldEyh, eyl, eyh);
    }
    /** 
     * Create new graph with stacked errors
     *  
     * @param g          Previous errors
     * @param ignoreErr  If true, ignore previous errors
     * @param quad       If true, add in quadrature 
     * 
     * @return Newly allocated graph 
     */
    Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) const
    {
      Graph* ga = new Graph(g->GetN());
      for (Int_t i = 0; i < g->GetN(); i++) { 
	Double_t x    = g->GetX()[i];	  
	Double_t y    = g->GetY()[i];
	Double_t exl  = g->GetErrorXlow(i);
	Double_t exh  = g->GetErrorXhigh(i);
	Double_t eyl  = g->GetErrorYlow(i);
	Double_t eyh  = g->GetErrorYhigh(i);
	StackPointError(i, XMode(), ignoreErr, quad, exl, exh, eyl, eyh);
	// AddError(i, g, quad, ignoreErr, true, exl, exh, eyl, eyh);
	ga->SetPoint(i, x, y);
	ga->SetPointError(i, exl,exh,eyl,eyh);
      }
      SetAttributes(ga);
      ga->SetTitle(GetTitle());
      ga->SetName(Form("stack_%08x", GetUniqueID()));
      return ga;
    }
    /** 
     * Sum up point errors 
     * 
     * @param i         Point number 
     * @param xMode     X-Mode
     * @param ignoreErr If true, ignore current errors
     * @param quad      If true, add in quadrature 
     * @param exl       Input: current  Output: New value (possibly square)
     * @param exh       Input: current  Output: New value (possibly square) 
     * @param eyl       Input: current  Output: New value (possibly square) 
     * @param eyh       Input: current  Output: New value (possibly square) 
     */
    void SumPointError(Int_t i, UShort_t xMode, Bool_t ignoreErr, Bool_t quad, 
		       Double_t& exl,
		       Double_t& exh,
		       Double_t& eyl,
		       Double_t& eyh) const
    {
      AddError(i, xMode, ignoreErr, quad, false, exl, exh, eyl, eyh);
    }
    /** 
     * Sum errors at point. Point i of g is updated
     * 
     * @param g            Where to sum
     * @param i            Point
     * @param ignoreErr    If true, ignore exusisting errros
     * @param quad         Add in quadrature
     * @param opt          Option
     */
    void SumError(Graph* g, Int_t i, Bool_t ignoreErr,
		  Bool_t quad, UInt_t opt) const
    {
      Double_t exl = (g ? g->GetErrorXlow(i)  : 0);
      Double_t exh = (g ? g->GetErrorXhigh(i) : 0);
      Double_t eyl = (g ? g->GetErrorYlow(i)  : 0);
      Double_t eyh = (g ? g->GetErrorYhigh(i) : 0);
      SumPointError(i, XMode(opt), ignoreErr, quad, exl, exh, eyl, eyh);
      g->SetPointError(i, exl,exh,eyl,eyh);
    }
    /* @} */
    void SavePrimitive(std::ostream& out, Option_t* option="")
    {
      if (option[0] == 'd' || option[0] == 'D')
	out << " // Point-2-Point systematic " << GetTitle() << "\n"
	    << "  {\n"
	    << "    Int_t id = g->DeclarePoint2Point(\"" << GetTitle() << "\","
	    << IsRelative() <<  ',' << fOption << ");\n"
	    << "    g->SetSysLineColor(id," << GetLineColor() << ");\n"
	    << "    g->SetSysLineStyle(id," << GetLineStyle() << ");\n"
	    << "    g->SetSysLineWidth(id," << GetLineWidth() << ");\n"
	    << "    g->SetSysFillColor(id," << GetFillColor() << ");\n"
	    << "    g->SetSysFillStyle(id," << GetFillStyle() << ");\n"
	    << "  }\n";
    }
    virtual void Print(Option_t* option="") const
    {
      gROOT->IndentLevel();
      Printf("%s/%s (ID # %d, %d points)", GetName(), GetTitle(), GetUniqueID(),
	     (fGraph ? fGraph->GetN() : -1));
      TString opt(option);
      if (!opt.Contains("ATTR", TString::kIgnoreCase)) return;
      gROOT->IndentLevel();
      Printf(" [option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d]",
	     fOption, GetLineColor(), GetLineStyle(), GetLineWidth(),
	     GetFillColor(), GetFillStyle());
      
    }
    /** Our data */
    // Graph* fGraph;
    TGraphAsymmErrors* fGraph;

    ClassDef(HolderP2P,3);
  };
  //__________________________________________________________________
  /**
   * Holder of common errors 
   * 
   */
  struct HolderCommon : public Holder 
  {
    /** Containing class is a friemd */
    friend struct GraphSysErr;
    /** 
     * CTOR
     */
    HolderCommon()
      : Holder(), fEyl(0), fEyh(0)
    {}
  protected:
    /** 
     * CTOR with name and title
     * 
     * @param name   Name 
     * @param title  Title 
     * @param rel    Relative or absolue 
     * @param opt    Draw Option 
     * @param id     Identifier 
     */
    HolderCommon(const char* name, const char* title, Bool_t rel,
		 UInt_t opt, UInt_t id)
      : Holder(name, title, rel, opt, id), fEyl(0), fEyh(0)
    {
    }
    /** 
     * Copy CTOR
     * 
     * @param other Object to copy from 
     */
    HolderCommon(const HolderCommon& other)
      : Holder(other),
	fEyl(other.fEyl), 
	fEyh(other.fEyh)
    {
    }
    /** 
     * Assignment operator 
     * 
     * @param other Object to assign from
     * 
     * @return Reference to this 
     */
    HolderCommon& operator=(const HolderCommon& other)
    {
      if (&other == this) return *this;
      this->Holder::operator=(other);
      fEyl = other.fEyl;
      fEyh = other.fEyh;
      return *this;
    }

    /** 
     * @{ 
     * @name Setting errors 
     */
    /** 
     * Set symmetric error 
     * 
     * @param ey Error
     */
    void Set(Double_t ey) { Set(ey, ey); }
    /** 
     * Set error
     * 
     * @param eyl Downward error along Y
     * @param eyh Upward error along Y
     */
    void Set(Double_t eyl, Double_t eyh)
    {
      fEyl = eyl;
      fEyh = eyh;
    }
    /** 
     * Get the down error 
     * 
     * @param y Value to evaluate at if the error is relative 
     * 
     * @return The down error 
     */ 
    Double_t GetYDown(Double_t y=0) const 
    {
      return (fRelative ? y : 1) * fEyl;
    }
    /** 
     * Get the up error 
     * 
     * @param y Value to evaluate at if the error is relative 
     * 
     * @return The up error 
     */
    Double_t GetYUp(Double_t y=0) const 
    {
      return (fRelative ? y : 1) * fEyh;
    }
    /* @} */

    /** 
     * @{ 
     * @name Sum, add, stack 
     */
    /** 
     * Add errors together at point
     * 
     * @param y            Point value
     * @param xMode        X-mode
     * @param quad         If true, add in quadrature
     * @param ignoreErr    If true, ignore errors on g
     * @param sqOld        If true and quad true, square old
     * @param y            Y value at point 
     * @param exl          On return, the left-hand X errors
     * @param exh          On return, the right-hand X errors
     * @param eyl          On return, the downward Y errors
     * @param eyh          On return, the upward Y errors 
     */
    void AddError(Double_t  y,
		  UShort_t  xMode,
		  Bool_t    ignoreErr,
		  Bool_t    quad, 
		  Bool_t    sqOld,
		  Double_t& exl, 
		  Double_t& exh, 
		  Double_t& eyl, 
		  Double_t& eyh) const
    {
      //exl         = (i   == 0        ? 0 :(g->GetX()[i]  -g->GetX()[i-1])/2);
      //exh         = (i+1 >= g->GetN()? 0 :(g->GetX()[i+1]-g->GetX()[i])/2);
      Double_t exlO = exl;
      Double_t exhO = exh;
      Double_t eylO = GetYDown(y);
      Double_t eyhO = GetYUp(y);
      if (exlO <= 0) exlO = gStyle->GetErrorX()/2;
      if (exhO <= 0) exhO = gStyle->GetErrorX()/2;

      DoAdd(xMode, exl, exh, eyl, eyh,
	    ignoreErr, quad, sqOld,
	    exlO, exhO, eylO, eyhO);
      exl = exlO;
      exh = exhO;
      eyl = eylO;
      eyh = eyhO;
    }
    /** 
     * Make a graph for showing next to data 
     * 
     * @param g     PRevious errors
     * @param quad  If true, add in quadrature 
     * @param x     Middle X coordinate 
     * @param y     Middle Y coordinate 
     * 
     * @return Newly allocated graph
     */
    Graph* BarError(Graph* g, Bool_t quad, Double_t x, Double_t y) const
    {
      Double_t exl = (g ? g->GetErrorXlow(0)  : 0);
      Double_t exh = (g ? g->GetErrorXhigh(0) : 0);
      Double_t eyl = (g ? g->GetErrorYlow(0)  : 0);
      Double_t eyh = (g ? g->GetErrorYhigh(0) : 0);
      AddError(y, XMode(), false, quad, quad, exl, exh, eyl, eyh);

      Graph* r = new Graph(1);
      r->SetPoint(0, x, y);
      r->SetPointError(0, exl, exh, eyl, eyh);
      
      SetAttributes(r);
      r->SetTitle(GetTitle());
      r->SetName(Form("bar_%08x", GetUniqueID()));

      return r;
    }
    /** 
     * Stack up point errors 
     * 
     * @param y         Point value
     * @param xMode     X-Mode
     * @param ignoreErr If true, ignore current errors
     * @param quad      If true, add in quadrature 
     * @param exl       Input: current  Output: New value 
     * @param exh       Input: current  Output: New value 
     * @param eyl       Input: current  Output: New value 
     * @param eyh       Input: current  Output: New value 
     */
    void StackPointError(Double_t  y,
			 UShort_t  xMode,
			 Bool_t    ignoreErr, 
			 Bool_t    quad,
			 Double_t& exl,
			 Double_t& exh, 
			 Double_t& eyl,
			 Double_t& eyh) const
    {
      AddError(y, xMode, ignoreErr, quad, true, exl, exh, eyl, eyh);
    }
    
    /** 
     * Create new graph with stacked errors
     *  
     * @param g          Previous errors
     * @param ignoreErr  If true, ignore previous errors
     * @param quad       If true, add in quadrature 
     * 
     * @return Newly allocated graph 
     */
    Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) const
    {
      Graph* ga = new Graph(g->GetN());
      for (Int_t i = 0; i < g->GetN(); i++) { 
	Double_t x    = g->GetX()[i];
	Double_t y    = g->GetY()[i];	
	Double_t exl  = g->GetErrorXlow(i);
	Double_t exh  = g->GetErrorXhigh(i);
	Double_t eyl  = g->GetErrorYlow(i);
	Double_t eyh  = g->GetErrorYhigh(i);
	StackPointError(y, XMode(), ignoreErr, quad, exl, exh, eyl, eyh);
	ga->SetPoint(i, x, y);
	ga->SetPointError(i, exl,exh,eyl,eyh);
      }
      SetAttributes(ga);
      ga->SetTitle(GetTitle());
      ga->SetName(Form("stack_%08x", GetUniqueID()));
      return ga;
    }
    /** 
     * Sum up point errors 
     * 
     * @param y         Point value 
     * @param xMode     X-Mode
     * @param ignoreErr If true, ignore current errors
     * @param quad      If true, add in quadrature 
     * @param exl       Input: current  Output: New value (possibly square)
     * @param exh       Input: current  Output: New value (possibly square) 
     * @param eyl       Input: current  Output: New value (possibly square) 
     * @param eyh       Input: current  Output: New value (possibly square) 
     */
    void SumPointError(Double_t y, UShort_t xMode, Bool_t ignoreErr,Bool_t quad,
		       Double_t& exl,
		       Double_t& exh,
		       Double_t& eyl,
		       Double_t& eyh) const
    {
      AddError(y, xMode, ignoreErr, quad, false, exl, exh, eyl, eyh);
    }
    /** 
     * Sum errors at point. Point i of g is updated
     * 
     * @param g            Where to sum
     * @param i            Point
     * @param ignoreErr    If true, ignore exusisting errros
     * @param quad         Add in quadrature
     * @param opt          Option
     */
    void SumError(Graph* g, Int_t i, Bool_t ignoreErr, Bool_t quad,
		  UInt_t opt) const
    {
      Double_t y   = (g ? g->GetY()[i]        : 0);
      Double_t exl = (g ? g->GetErrorXlow(i)  : 0);
      Double_t exh = (g ? g->GetErrorXhigh(i) : 0);
      Double_t eyl = (g ? g->GetErrorYlow(i)  : 0);
      Double_t eyh = (g ? g->GetErrorYhigh(i) : 0);
      SumPointError(y, XMode(opt), ignoreErr, quad, exl, exl, eyl, eyh);
      g->SetPointError(i, exl,exh,eyl,eyh);
    }
    /* @} */
    void SavePrimitive(std::ostream& out, Option_t* option="")
    {
      if (option[0] == 'd' || option[0] == 'D')
	out << " // Common systematic " << GetTitle() << "\n"
	    << "  {\n"
	    << "    Int_t id = g->DefineCommon(\"" << GetTitle() << "\","
	    << fRelative << ',' << fEyl << ',' << fEyh << ','
	    << fOption << ");\n"
	    << "    g->SetSysLineColor(id," << GetLineColor() << ")\n;"
	    << "    g->SetSysLineStyle(id," << GetLineStyle() << ")\n;"
	    << "    g->SetSysLineWidth(id," << GetLineWidth() << ")\n;"
	    << "    g->SetSysFillColor(id," << GetFillColor() << ")\n;"
	    << "    g->SetSysFillStyle(id," << GetFillStyle() << ")\n;"
	    << "  }\n";
    }
    virtual void Print(Option_t* option="") const
    {
      Bool_t      rel = IsRelative();
      Double_t    fac = (rel ? 100 : 1);
      const char* pst = (rel ? "%" : "");
      gROOT->IndentLevel();
      Printf("%s/%s (ID # %d) -%f%s +%f%s",
	     GetName(), GetTitle(), GetUniqueID(),
	     fac*fEyl, pst, fac*fEyh, pst);
      TString opt(option);
      if (!opt.Contains("ATTR", TString::kIgnoreCase)) return;
      gROOT->IndentLevel();
      Printf(" [option: %3d line (c/s/w):%3d/%1d/%2d fill (c/s):%3d/%4d]",
	     fOption, GetLineColor(), GetLineStyle(), GetLineWidth(),
	     GetFillColor(), GetFillStyle());
      
    }
    /** Down errors */
    Double_t fEyl;
    /** Up errors */
    Double_t fEyh;

    ClassDef(HolderCommon,3);
  };
  /** 
   * Get the point value and low and high errors
   * 
   * @param i         Point number 
   * @param cmn       Consider commons
   * @param stat      Consider statistics 
   * @param quad      Add in quadrature 
   * @param nosqrt    Do not take square root in case quad=true
   * @param eyl       Output: Low error
   * @param eyh       Output: high error
   * 
   * @return 
   */
  Double_t GetYandError(Int_t     i,
			Bool_t    cmn,
			Bool_t    stat,
			Bool_t    quad,
			Bool_t    nosqrt,
			Double_t& eyl,
			Double_t& eyh) const
  {
    Double_t wyl = 0;
    Double_t wyh = 0;
    return GetYandError(i, cmn, stat, quad, nosqrt, eyl, eyh, wyl, wyh);
  }
  /** 
   * Get the point value and low and high errors.  Errors that have
   * been marked as used (kUsedBit) are not considered.  Errors marked
   * as only for weights (kWeightsOnlyBit) are only factored in on the
   * weights (@a weyl, @a weyh) calculation.  Statistical errors are
   * always factored in on the calculation of weights.
   * 
   * @param i         Point number 
   * @param cmn       Consider commons
   * @param stat      Consider statistics 
   * @param quad      Add in quadrature 
   * @param nosqrt    Do not take square root in case quad=true
   * @param eyl       Output: Low error
   * @param eyh       Output: high error
   * @param wyl       Output: Low weight
   * @param wyh       Output: high weight
   * 
   * @return The graph value at point @a i
   */
  Double_t GetYandError(Int_t     i,
			Bool_t    cmn,
			Bool_t    stat,
			Bool_t    quad,
			Bool_t    nosqrt,
			Double_t& eyl,
			Double_t& eyh,
			Double_t& wyl,
			Double_t& wyh) const
  {
    // --- Find location for common errors ---------------------------
    Int_t n = fData->GetN();
    if (i >= n) {
      eyl = -1;
      eyh = -1;
      return 0;
    }
    Double_t y   = fData->GetY()[i];
    wyl = GetStatErrorDown(i);
    wyh = GetStatErrorUp(i);
    eyl = (stat ? wyl : 0);
    eyh = (stat ? wyh : 0);
    if (quad) {
      wyl *= wyl;
      wyh *= wyh;
      eyl *= eyl;
      eyh *= eyh;
    }
    Double_t exl = GetErrorXLeft(i);
    Double_t exh = GetErrorXRight(i);

    // Otherwise, we are adding all selected errors togethter 
    // Graph* g = static_cast<Graph*>(fData->Clone("error"));
    // g->SetTitle(fSumTitle);
    Int_t xMode = -1;
    if (cmn) { 
      TIter nextC(&fCommon);
      HolderCommon* hc = 0;
      while ((hc = static_cast<HolderCommon*>(nextC()))) {
	if (hc->TestBit(kUsedBit)) continue;
	if (xMode < 0) xMode = hc->XMode(fSumOption);
	hc->SumPointError(y, xMode, false, quad, exl, exh, wyl, wyh);
	if (hc->TestBit(kOnlyWeightBit)) continue;
	hc->SumPointError(y, xMode, false, quad, exl, exh, eyl, eyh);
      }
    }
    TIter nextP(&fPoint2Point);
    HolderP2P* hp = 0;
    while ((hp = static_cast<HolderP2P*>(nextP()))) {
      if (hp->TestBit(kUsedBit)) continue;
      if (xMode < 0) xMode = hp->XMode(fSumOption);	
      hp->SumPointError(i, xMode, false, quad, exl, exh, wyl, wyh);
      if (hp->TestBit(kOnlyWeightBit)) continue;
      hp->SumPointError(i, xMode, false, quad, exl, exh, eyl, eyh);
    }
    if (quad && !nosqrt) {
      eyl = TMath::Sqrt(eyl);
      eyh = TMath::Sqrt(eyh);
      wyl = TMath::Sqrt(wyl);
      wyh = TMath::Sqrt(wyh);
    }
    return y;
  }
  /** 
   * Round number @a v to @a n significant digits.  It returns the
   * mantisa and the exponent.  To calculate the number, do
   *
   * @code 
   Int_t    e = 0;
   Double_t m = Round(v,n,e);
   Double_t x = m * TMath::Power(10,e);
   @endcode 
   *
   * If one has a value @e x and and associated error @e dx, one can
   * use this function to format the pair:
   *
   @code 
   Int_t     e = 0;
   Double_t  m = Round(dx,n,e);
   Double_t  p = m * TMath::Power(10,e);
   Int_t     o = TMath::Ceil(TMath::Log10(TMath::Abs(x)))-e;
   std::cout << std::setprecision(o) << x << " +/- " 
             << std::setprecision(n) << p << std::endl;
   @endcode 
   * 
   * For asymmetric errors, this would be 
   *
   @code 
   Int_t     el, eh;
   Double_t  ml = Round(dxl,n,el);
   Double_t  m2 = Round(dxh,n,eh);
   Int_t     ne = TMath::Min(el, eh);
   Double_t  pl = ml * TMath::Power(10,el);
   Double_t  ph = mh * TMath::Power(10,eh);
   Int_t     o  = TMath::Ceil(TMath::Log10(TMath::Abs(x)))-ne;
   std::cout << std::setprecision(o) << x  << " -" 
             << std::setprecision(n) << pl << " +"
             << std::setprecision(n) << ph << std::endl;
   @endcode 
   * 
   * @param v      Value to round
   * @param p      Number of signficant digits  
   * @param rexpo  On return, the exponent 
   * 
   * @return The "mantisa"
   */
  static Double_t Round(Double_t v, Int_t p, Int_t& rexpo)
  {
    if (v == 0) {
      // ret.Form("%.*f", p-1, 0);
      rexpo = 0;
      return 0; // ret.Data();
    }

    // Get the sign, take absolute value, get exponent, get scalar, scale 
    Bool_t   neg  = (v < 0);
    Double_t tmp  = TMath::Abs(v);
    Int_t    expo = Int_t(TMath::Log10(tmp));
    Double_t tens = TMath::Power(10, expo - p + 1);
    Double_t n    = RoundN(tens, tmp); // TMath::Ceil(tmp/tens);

    if (n < TMath::Power(10, p-1)) {
      // If scaled is less than 10 to our precision, take off one digit
      expo--;
      tens = TMath::Power(10, expo - p + 1);
      n    = RoundN(tens,tmp); // TMath::Ceil(tmp/tens);
    }

    if (TMath::Abs((n+1) * tens) <= TMath::Abs(n * tens)) 
      // If the value is not near original, add one
      n++;

    if (n >= TMath::Power(10,p)) {
      // If value is larger than 10 to our precision, scale down one
      n /= 10;
      expo++;
    }

    rexpo = expo-p+1;
    return (neg ? -1 : 1)*n;
  }
protected:
  static void SwapPoints(Graph* g, Int_t i, Int_t j, Bool_t reflect)
  {
    Double_t s        = (reflect ? -1 : 1);
    Double_t tmpX     = g->GetX()[i];
    Double_t tmpY     = g->GetY()[i];
    Double_t tmpExl   = g->GetErrorXlow(i);
    Double_t tmpExh   = g->GetErrorXhigh(i);
    Double_t tmpEyl   = g->GetErrorYlow(i);
    Double_t tmpEyh   = g->GetErrorYhigh(i);
    if (reflect) {
      // Swap low and high X error when reflecting 
      Double_t tmp = tmpExl;
      tmpExl       = tmpExh;
      tmpExh       = tmp;
    }
    g->SetPoint(i, s*g->GetX()[j], g->GetY()[j]);
    g->SetPointError(i,
		     g->GetErrorXlow(j), g->GetErrorXhigh(j),
		     g->GetErrorYlow(j), g->GetErrorYhigh(j));
    g->SetPoint(j, s*tmpX, tmpY);
    g->SetPointError(j, tmpExl, tmpExh, tmpEyl, tmpEyh);
  }
  /** 
   * @{ 
   * @name Helpers for importing/exporting
   */
  static const char* FormatKey(const char* key)
  {
    TString tmp(Form("*%s:", key));
    return Form("%-12s", tmp.Data());
  }
  /** 
   * Extract a field from a string 
   * 
   * @param value The string 
   * @param idx   Which index 
   * 
   * @return String at index, or last value 
   */
  static const char* ExtractField(const TString& value, Int_t idx)
  {
    static TString val;
    val = "";
    TObjArray* tokens = value.Tokenize(":");
    Int_t      iVal   = TMath::Min(Int_t(idx),tokens->GetEntriesFast()-1);
    if (iVal < 0) return val;
    TString    tmp    = tokens->At(iVal)->GetName();
    val               = tmp.Strip(TString::kBoth, ' ');
    // tokens->ls();
    // ::Info("", "Selecting %d: %s", iVal, val.Data());
    delete tokens;
    return val.Data();
  }
  /** 
   * Store a qualifier in a table 
   * 
   * @param quals Table. 
   * @param idx   Column number
   * @param name  Row name 
   * @param val   The cell content
   */
  static void StoreQual(TList& quals, Int_t idx,
			const char* name, const char* val)
  {
    TObject* oqv = quals.FindObject(name);
    TList*   qv  = 0;
    if (!oqv) {
      // ::Info("StoreQual", "Creating list for qualifier %s", name);
      qv = new TList;
      qv->SetOwner();
      qv->SetName(name);
      quals.Add(qv);
    }
    else qv = static_cast<TList*>(oqv);
    // ::Info("StoreQual", "Storing %d value %s=%s", idx, name, val);
    qv->AddAt(new TObjString(val), idx);
  }
  /** 
   * Store a qualifier in a table 
   * 
   * @param quals Table. 
   * @param idx   Column number
   * @param q     Key, value pair. The key is the row name
   */
  static void StoreQual(TList& quals, Int_t idx, TObject* q)
  {
    StoreQual(quals, idx, q->GetName(), q->GetTitle());
  }
  /** 
   * Export all values of a key 
   * 
   * @param out     Output stream
   * @param which   Which key 
   * @param alsoKey If true, output key value 
   * @param fill    Filler in case the key isn't defined 
   * 
   * @return The output stream 
   */
  std::ostream& ExportKey(std::ostream& out,
			  const char*   which,
			  Bool_t        alsoKey = true,
			  const char*   fill    = "<please fill in>") const
  {
    Bool_t gotit = false;
    if (fMap) { 
      TIter nextK(fMap);
      TObject* obj = 0;
      while ((obj = nextK())) { 
	TString k(obj->GetName());
	if (!k.EqualTo(which)) continue;
	gotit = true;
	if (alsoKey)  out << FormatKey(obj->GetName()); 
	out << obj->GetTitle() << std::endl;
      }
    }
    if (!gotit) {
      if (alsoKey)
	out << FormatKey(which);
      out << fill << std::endl;
    }
    return out;
  }
  /** 
   * Export data set header, and possibly file header too 
   * 
   * @param out         Output stream
   * @param alsoTop     If true, export file header 
   * @param alsoComment If true, also write out comment 
   * @param nsign  Number of significant digits  
   * 
   * @return output stream 
   */
  std::ostream& ExportHeader(std::ostream& out,
			     Bool_t alsoTop=false,
			     Bool_t alsoComment=false,
			     Int_t  nsign=-1) const
  {
    if (alsoComment) {
      out << "# Generated by GraphSysErr\n"
	  << "# See also\n"
	  << "# \n"
	  << "#  http://hepdata.cedar.ac.uk/resource/sample.input\n"
	  << "#  http://hepdata.cedar.ac.uk/submittingdata\n"
	  << "# \n"
	  << "# Please fill in missing fields\n";
    }
    if (alsoTop) { 
      TDatime      now;
      UserGroup_t* u        = gSystem->GetUserInfo();
      TString      ref      = GetKey("reference");
      const char*  lab      = GetKey("laboratory");
      const char*  accel    = GetKey("accelerator");
      const char*  exper    = GetKey("detector");
      const char*  abs      = GetKey("abstract");
      const char*  months[] = {"JAN","FEB","MAR","APR",
			       "MAY","JUN","JUL","AUG",
			       "SEP","OCT","NOV","DEC"};
      if (ref.IsNull()) 
	out << "*reference: <journal/archive ref> : <year>" << std::endl;
      else {
	ExportKey(out, "reference", true, "<reference>");
      }
      ExportKey(out, "author",    true, "<Last name of first author in CAPS>");
      ExportKey(out, "doi",       true, "<DOI number>");
      ExportKey(out, "inspireId", true, "<inSpire ID>");
      ExportKey(out, "cdsId",     true, "<CDS ID>");
      ExportKey(out, "durhamId",  true, "<fill in on submission>");
      ExportKey(out, "title",     true, "<Title of paper>");
      
      out << FormatKey("status") << "Encoded " 
	  << std::setw(2) << now.GetDay() <<' '<< months[now.GetMonth()-1]<<' '
	  << std::setw(4) << now.GetYear() << " by "<< u->fRealName.Data()<<"\n"
	  << FormatKey("experiment") << (lab ? lab : "<lab>") << '-' 
	  << (accel ? accel : "<accelerator>") << '-' 
	  << (exper ? exper : "<experiment>") << '\n'
	  << FormatKey("detector") << (exper ? exper : "<experiment>") << '\n'
	  << FormatKey("comment")  << (lab ? lab : "<lab>") << '-' 
	  << (accel ? accel : "<accelerator>") << ". " 
	  << (abs ? abs : "<abstract>") << "\n"
	  << std::endl;
    }
    out << "# Start of dataset\n";
    if (nsign >= 0) out << "# Using " << nsign << " significant digits\n";
    out << FormatKey("dataset") << std::endl;
    // out << FormatKey("dscomment") << GetTitle() << std::endl;
    const char* fields[] = { "location",
			     "reackey", 
			     "obskey",
			     "dscomment",
			     // "qual", 			     
			     0 };
    const char** pfld = fields;
    while (*pfld) {
      ExportKey(out, *pfld, true);
      pfld++;
    }
    return out;
  }
  /** 
   * Export an error 
   * 
   * @param o    Output stream 
   * @param low  Low error 
   * @param high High error 
   * @param nopm If true, do not prefix symmetric errors with +/-
   * @param rel  IF true, the error is relative 
   * @param nsign  Number of significant digits  
   * 
   * @return Least exponent used
   */
  static Int_t ExportError(std::ostream& o,
			   Double_t low,
			   Double_t high,
			   Bool_t   nopm,
			   Bool_t   rel,
			   Int_t    nsign)
  {
    Double_t l = low;
    Double_t h = high;
    Int_t    p = o.precision();
    Int_t    e = 0;
    if (nsign >= 0) {
      Int_t    pp = nsign;
      Int_t    le, he;
      Double_t lm = Round(l, nsign, le);
      Double_t hm = Round(h, nsign, he);
      e           = TMath::Min(le,he);
      l           = lm*TMath::Power(10,le);
      h           = hm*TMath::Power(10,he);
      o.precision(pp);
    }
    // ::Info("ExportError","Before (%2d) -%8f.+%8f After -%8f,+%8f",
    //        nsign, low,high,l,h);
    if (TMath::Abs(h-l) < 1e-12)
      o << (nopm ? "" : "+- ") << TMath::Abs(l) << (rel ? " PCT" : "");
    else
      o << "+" << h << ",-" << l << (rel ? " PCT" : "");
    o.precision(p);
    return e;
  }
  /** 
   * Export a single point 
   * 
   * @param out     Output stream 
   * @param i       Point number
   * @param alsoX   If true, also export X coordinate 
   * @param sysName If true, export P2P names
   * @param nsign  Number of significant digits  
   * 
   * @return output stream 
   */
  std::ostream& ExportPoint(std::ostream& out,
			    Int_t         i,
			    Bool_t        alsoX=true,
			    Bool_t        sysName=true,
			    Int_t         nsign=0) const
  {
    Int_t p = out.precision();
    // out.setf(std::ios::fixed);
    if (alsoX) {
      Double_t x       = GetX(i);
      Double_t exl     = fData->GetErrorXlow(i);
      Double_t exh     = fData->GetErrorXhigh(i);
      Int_t    pp      = p;
      Int_t    ppp     = p;
      if (nsign > 0) {
	Int_t    pxl, pxh;
	Double_t mxl   = Round(exl,nsign,pxl);
	Double_t mxh   = Round(exh,nsign,pxh);
	exl            = mxl*TMath::Power(10, pxl);
	exh            = mxh*TMath::Power(10, pxh);
	Int_t    mp    = TMath::Min(pxl,pxh);
	pp             = nsign;
	ppp            = TMath::Ceil(TMath::Log10(TMath::Abs(x)))-mp;
	// Info("","pxl=%d pxh=%d mp=%d ppp=%d", pxl, pxh, mp, ppp);
      }
      
      if (TMath::Abs(exl) < 1e-10 && TMath::Abs(exh) < 1e-10) {
	if (nsign <= 0) out << ' ' << x;
	else {
	  Int_t    px;
	  Double_t mx = Round(x, nsign, px);
	  Double_t xx = mx*TMath::Power(10,px);
	  // Printf("round %f -> %f (%d) -> %.*g", x, mx, px, nsign, xx);
	  out.precision(nsign);
	  out << ' ' << xx;
	}
      }
      else {
	if (TMath::Abs(exh-exl) < 1e-10) {
	  out.precision(ppp);
	  out << ' ' << x-exl << " TO " << x+exh;
	}
	else {
	  out.precision(ppp);
	  out << ' ' << x << ' ';
	  out.precision(pp);
	  out << '+' << exh << ",-" << exl;
	}
      }
      out.precision(p);
      out << "; ";
    }
    // out.unsetf(std::ios::fixed);
	
    Bool_t   statRel = fStatRelative;
    Double_t y       = GetY(i);
    Double_t fy      = (statRel ? (y == 0 ? 0 : 100./y) : 1);
    Double_t eyl     = GetStatErrorDown(i) * fy;
    Double_t eyh     = GetStatErrorUp(i)   * fy;
    // out << y << ' ';
    std::stringstream tmp;
    // tmp.setf(std::ios::fixed);
    Int_t le = ExportError(tmp, eyl, eyh, false, statRel, nsign);
      
    if (fPoint2Point.GetEntries() > 0) {
      tmp  << " (" << std::flush;
      TIter       nextP(&fPoint2Point);
      HolderP2P*  holderP2P = 0;
      Bool_t      first = true;
      while ((holderP2P = static_cast<HolderP2P*>(nextP()))) { 
	Bool_t rel = holderP2P->IsRelative();
	fy         = (rel ? (y == 0 ? 0 : 100./y) : 1);
	if (!first) tmp << ',';
	Double_t esh = holderP2P->GetYUp(i)*fy;
	Double_t esl = holderP2P->GetYDown(i)*fy;
	tmp << "DSYS=";
	le = TMath::Min(le,ExportError(tmp, esl, esh, true, rel, nsign));
	if (sysName) 
	  tmp << ':'  << holderP2P->GetTitle() << std::flush;
	first = false;
      }
      tmp << ')';
    }
    if (nsign) {
      Int_t nY = TMath::Ceil(TMath::Log10(TMath::Abs(y)))-le;
      out.precision(nY);
      out << y << ' ';
      out.precision(p);
    }
    else
      out << y << ' ';
    out << tmp.str() << ';';
    return out;
  }
  /** 
   * Round number.  @a tens is the scaling divisor to convert @a tmp
   * into an integer with @e n significant digits. 
   *
   * The algorightm looks at the @e n+1 significant digit @e d, and
   * depending on the value of that, it does one of the following:
   *
   * - if @e d > 5, it rounds the number up 
   * - if @e d < 5, it rounds the number down. 
   * - if @e d = 5, then the @e n+2 significant digit @e dd is investigated. 
   *   - if @e dd = 0 then the number is rounded up nearest even number
   *   - if @e dd != 0 the the number is rounded down 
   * 
   * @param tens The scaling of the number to integer
   * @param tmp  The absolute value of the number 
   * 
   * @return The rounded number 
   */
  static Double_t RoundN(Double_t tens, Double_t tmp)
  {
    // Add small number to not loose least digit if a 1
    Double_t n     = TMath::Floor(100*tmp/tens+.00001);
    Int_t    next  = int(n/10) % 10;
    if (next > 5) 
      // Round up 
      return TMath::Ceil(n/100);
    if (next < 5)
      // Round down
      return TMath::Floor(n/100);
    Int_t nnext = int(n) % 10;
    if (nnext != 0)
      // Round up 
      return TMath::Ceil(n/100);
    Int_t last = int(n/100) % 10;
    if (last % 2 == 0)
      // Round down
      return TMath::Floor(n/100);
    // Round to nearest even
    return TMath::Ceil(n/100);
  }  
  /** 
   * Get the ith token from the array of tokens c. 
   *  
   * @param c         Array of tokens
   * @param idx       Which to get 
   * @param verbose   Be verbose 
   * 
   * @return Rererence to the token, or the empty string 
   */
  static TString& Token(TObjArray* c, UShort_t idx, Bool_t verbose=true)
  {
    static TString empty("");
    if (idx >= c->GetEntriesFast()) { 
      if (verbose) {
	::Error("Token", "Token %d does not exist in array", idx);
	c->ls();
      }
      return empty;
    }
    return static_cast<TObjString*>(c->At(idx))->String();
  }
  /** 
   * Import errors from a string 
   * 
   * @param s     String to parse
   * @param el    On return, the low error
   * @param eh    On return, the high error
   * @param rel   On returm, true if relative
   * 
   * @return true on success
   */
  static Bool_t ImportError(const TString& s, 
			    Double_t& el, Double_t& eh, Bool_t& rel)
  {
    TString tmp(s);
    // ::Info("ImportError", "Parsing '%s' for errors", s.Data());
    // Allow for PCT (any case) and % as relative markers 
    if (tmp.Contains("PCT", TString::kIgnoreCase) || tmp.Contains("%")) 
      rel = true;
    if (tmp.BeginsWith("+-") || tmp.BeginsWith("-+")) {
      el = eh = tmp.Atof();
      // ::Info("ImportError", "Got symmetric %f,%f", el, eh);
      return true;
    }
    TObjArray* tokens = tmp.Tokenize(",");
    for (Int_t i = 0; i < tokens->GetEntriesFast(); i++) { 
      TString& t = Token(tokens, i);
      if (t.IsNull()) continue;
      t.Remove(TString::kBoth, ' ');
      Double_t v = t.Atof();
      // ::Info("ImportError", "Got token '%s' -> %f", t.Data(), v);
      if (t[0] != '-' && t[0] != '+') {
	// ::Info("ImportError", "First char not + or -: %c", t.Data()[0]);
	el = eh = TMath::Abs(v);
      }
      else if (v < 0) {
	// ::Info("ImportError", "Value %f < 0 -> setting el=%f", v, -v);
	if (el > 0) { 
	  Double_t max = TMath::Max(el, -v);
	  ::Warning("ImportError", 
		    "Lower error already set to %f, chosing the maximum of "
		    "(%f,%f) -> %f", el, el, -v, max);
	  v = -max;
	}
	el = -v;
      }
      else {
	// ::Info("ImportError", "Value %f >= 0 -> setting eh=%f", v, v);
	if (eh > 0) { 
	  Double_t max = TMath::Max(eh, v);
	  ::Warning("ImportError", 
		    "Lower error already set to %f, chosing the maximum of "
		    "(%f,%f) -> %f", eh, eh, v, max);
	  v = max;
	}
	eh = v;
      }
    }
    // ::Info("ImportError", "Got a-symmetric %f,%f", el, eh);
    tokens->Delete();
    delete tokens;
    return true;
  }
      
  /** 
   * Service function to import a point value (X or Y) with errors. 
   * 
   * @param s     String to parse
   * @param v     On return, the value 
   * @param el    On return, the lower error
   * @param eh    On return, the higher error
   * @param rel   On return, true if the errors are relative to the value
   * 
   * @return true on success 
   */      
  static Bool_t ImportPoint(const TString& s, 
			    Double_t&      v, 
			    Double_t&      el, 
			    Double_t&      eh, 
			    Bool_t&        rel)
  {
    // ::Info("ImportPoint", "s=%s", s.Data());
    rel            = s.Contains("PCT");
    TObjArray* tok = s.Tokenize(" ,");
    if (tok->GetEntriesFast() >= 4 && 
	Token(tok,2).EqualTo("TO", TString::kIgnoreCase) && 
	Token(tok,1).Contains("BIN", TString::kIgnoreCase)) {
      v = Token(tok,0).Atof();
      TString tmp = Token(tok,1);
      TString s1  = tmp.Strip(TString::kBoth, '(');
      s1.ReplaceAll("BIN=","");
      Double_t v1 = s1.Atof();
      tmp         = Token(tok,3);
      s1          = tmp.Strip(TString::kBoth, ')');
      Double_t v2 = s1.Atof();
      el          = (v-v1);
      eh          = (v2-v);            
    }
    else if (tok->GetEntriesFast() >= 3 && 
	     Token(tok,1).EqualTo("TO", TString::kIgnoreCase)) {
      Double_t v1 = Token(tok,0).Atof();
      Double_t v2 = Token(tok,2).Atof();
      v           = (v1+v2)/2;
      el          = (v-v1);
      eh          = (v2-v);
    }
    else {
      v  = Token(tok,0).Atof();
      if (tok->GetEntriesFast() >= 2 && 
	  (Token(tok,1).BeginsWith("+-") || 
	   Token(tok,1).BeginsWith("-+"))) {
	TString t = Token(tok,1);
	if (t.Length() > 2) {
	  t.Remove(0,2);
	  // ::Info("ImportPoint", "Extract from %s", t.Data());
	  el = eh = t.Atof();
	}
	else {
	  // ::Info("ImportPoint","Extract from next %s",Token(tok, 2).Data());
	  el = eh = Token(tok,2).Atof();
	}
      }
      else {
	for (Int_t i = 1; i < tok->GetEntriesFast(); i++) {
	  TString t = Token(tok, i);
	  Char_t  m = t[0];
	  // ::Info("ImportPoint", "%s -> %c", t.Data(), m);
	  if (m != '+' && m != '-') continue;
	  if (t.Length() == 1) t = Token(tok, ++i); // Take next
	  else                 t.Remove(0,1);       // Remove sign
	  // ::Info("ImportPoint", "%s -> %f", t.Data(), t.Atof());
	  if (m == '-') el = t.Atof();
	  else          eh = t.Atof(); 
	}
      }
    }
    tok->Delete();
    delete tok;
    return true;
  }
  /* @} */
  /** 
   * @{ 
   * @name Calculations 
   */
  /** 
   * Take the square root of the errors at point
   * 
   * @param g Graph
   * @param i Point
   */
  void SqrtPoint(Graph* g, Int_t i) 
  {
    if (kVerbose & kDraw) 
      printf("%20s Point %3d %f,%f -> ","",
	     i,g->GetEYlow()[i],g->GetEYhigh()[i]);

    g->GetEYlow()[i]  = TMath::Sqrt(g->GetEYlow()[i]);
    g->GetEYhigh()[i] = TMath::Sqrt(g->GetEYhigh()[i]);

    if (kVerbose & kDraw)
      Printf("%f,%f", g->GetEYlow()[i], g->GetEYhigh()[i]);
  }
  /** 
   * Take the square root of the errors 
   * 
   * @param g Graph
   */
  void SqrtGraph(Graph* g) 
  {
    for (Int_t i = 0; i < g->GetN(); i++) 
      SqrtPoint(g, i);
  }
  /* @} */
  /**
   * @{ 
   * @name Graphics 
   */
  /** 
   * Parse options 
   * 
   * @param opt Options to parse 
   */
  const char* FormatOption(UInt_t opt)
  {
    static TString ret;
    ret = "";
    switch (opt) {
    case kNormal:  ret = "";     break; // Line with ticks
    case kNoTick:  ret = "Z0";	break; // Line with no ticks 
    case kArrow:   ret = ">0";	break; // Linw with arrows 
    case kRect:    ret = "20";	break; // Rectangle w/fill
    case kBox:     ret = "50";	break; // Rectangle w/fill & line 
    case kFill:    ret = "30";	break; // Filled area 
    case kCurve:   ret = "40";	break; // Filled smoothed area 
    case kHat:     ret = "[]0";	break; // Hats 
    case kBar:     ret = "||0";	break; // A bar 
    case kNone:    ret = "XP";	break; // No errors
    case kLine:    ret = "CX";   break;
    case kConnect: ret = "LX";   break;
    }
    return ret.Data();
  }
  /** 
   * Make our data graph
   * 
   * @param n Optional number of pre-allocated points
   */
  void MakeDataGraph(Int_t n)
  {
    fData = new Graph(n);
    fData->SetName(Form("%s_data", fName.Data()));
    fData->SetTitle(GetTitle());
  }
  /** 
   * Find point (or possibly two points) that match X
   * 
   * @param x   X to mach 
   * @param i1  On return, the two points
   * @param i2  On return, the two points
   * @param fac Factor on errors for searching 
   * 
   * @return 
   */
  Int_t FindPoint(Double_t x, Int_t& i1, Int_t& i2, Double_t fac=1) const
  {
    i1 = -1;
    i2 = -1;
    Int_t n = GetN();
    if (x < GetX(0)-fac*GetErrorXLeft(0) ||
	x > GetX(n-1)+fac*GetErrorXRight(n-1)) {
      // If we're out side the coverage, set return to -1
      // Info("FindPoint", "X=%f Out of range [%f,%f]",
      //      x,GetX(0),GetX(GetN()-1));
      i1 = i2 = -1;
      return -2;
    }
    const Double_t tol = 1e-9;
    for (Int_t i=0; i < n; i++) {
      // if (TMath::Abs(GetX(i)-x) < tol) {
      if (x <= GetX(i)+fac*GetErrorXRight(i) &&
	  x >= GetX(i)-fac*GetErrorXLeft(i)) {
	// if ((GetX(i)+GetErrorXRight(i)-x)<tol &&
	//  (x - GetX(i) - GetErrorXLeft(i) <tol)) {
	// Found matching point
	i2 = i1 = i;
	break;
      }
      // Info("FindPoint", "Not @ %d (%f not in %f -%f +%f)",
      //      i, x, GetX(i), GetErrorXLeft(i), GetErrorXRight(i));
      if (x > GetX(i)) {
	// X still larger
	i1 = i;
      }
      else {
	// X is lower
	i2 = i;
	break;
      }
    }
    if (i1 == i2)
      // If equal, return the point 
      return i1;
    if (TMath::Abs((GetX(i1)+fac*GetErrorXRight(i1)) -
		   (GetX(i2)-fac*GetErrorXLeft(i2))) > 10*tol) {
      // If the two found points are not adjecent, we are in a hole
      // and we indicate that.
      // Info("FindPoint","Found hole at i1=%d(%f+%f) i2=%d(%f-%f) (%f)",
      //      i1, GetX(i1), GetErrorXRight(i1),
      //      i2, GetX(i2), GetErrorXLeft(i2), x);
      // Print("xy");
      i1 = i2 = -1;
      return -2;
    }
    if (i2 > 0) {
      // Otherwise return negative value to indicate we should look at
      // two points.
      // Info("FindPoint","Found i1=%d i2=%d %f in [%f,%f]", i1, i2,
      //      x, GetX(i1)-GetErrorXLeft(i1), GetX(i2)+GetErrorXRight(i2));
      return -1;
    }
    // Return negative two to indicate nothing was found
    return -2;
	
  }
  /** 
   * Make our stack 
   *
   * @param option Options - See method Draw for more
   *
   * @return Pointer to stack 
   */
  TMultiGraph* MakeMulti(Option_t* option)
  {
    if (!fData) {
      Warning("MakeMulti", "No data defined");
      return 0;
    }

    // --- Process options -------------------------------------------
    TString opt(option);
    opt.ToUpper();

    Bool_t   cmn      = opt.Contains("COMMON");
    Bool_t   combine  = opt.Contains("COMBINED");
    Bool_t   stat     = opt.Contains("STAT");
    Bool_t   quad     = opt.Contains("QUAD");
    combine           = !opt.Contains("STACK");
    quad              = !opt.Contains("DIRECT");
    Bool_t   split    = opt.Contains("SPLIT");
    Int_t    xpos     = (opt.Contains("WEST") ? -1 : 1);
    Int_t    ypos     = (opt.Contains("MIN") ? -1 : 
			 opt.Contains("MAX") ?  1 : 0);

    // --- Some general information ----------------------------------
    Int_t    n       = fData->GetN();
    Double_t dx       = TMath::Max(gStyle->GetErrorX(),0.0001F);
    if (n <= 0) {
      Warning("MakeMulti", "Nothing to do (n=%d)", n);
      return 0;
    }

    // --- Find location for common errors ---------------------------
    Double_t xBase    = 0;
    if (opt.Contains("XBASE=")) {
      Int_t   idx = opt.Index("XBASE=")+6;
      TString tmp = opt(idx,opt.Length()-idx);
      xBase       = tmp.Atof();
    }
    else {
      xBase   = ((xpos < 0 ? 
		  fData->GetX()[0]   - fData->GetEXlow()[0]: 
		  fData->GetX()[n-1] + fData->GetEXhigh()[n-1])
		 + xpos * 1.2 * dx);
    }
    Double_t yBase   = fData->GetY()[0];
    for (Int_t i = 0; i < n; i++) { 
      Double_t y = fData->GetY()[i];
      yBase      = (ypos < 0 ? TMath::Min(yBase, y) : 
		    ypos > 0 ? TMath::Max(yBase, y) : yBase + y);
    }
    if (ypos == 0) yBase /= n;

    // --- Create stack and temp collection --------------------------
    TList drawn;
    drawn.SetOwner(false);

    Graph* prev = 0;
    if (!combine) { 
      // If we're not combining, we basically copy the graphs to the
      // draw list, but modify the errors to show progresive larger
      // errors.
      prev = fData;
      if (cmn) { 
	// If we should add common errors, do that first 
	TIter next(&fCommon);
	Holder* h = 0;
	while ((h = static_cast<Holder*>(next()))) {
	  Graph* g = h->StackError(prev, (!stat ? prev == fData : false), quad);
	  if (!g) continue;

	  if (quad) SqrtGraph(g);
	  prev = g;
	  drawn.AddFirst(g, FormatOption(h->GetDOption()));
	}
      }
      TIter next(&fPoint2Point);
      Holder* h = 0;
      while ((h = static_cast<Holder*>(next()))) {
	Graph* g = h->StackError(prev, (!stat ? prev == fData : false), quad);
	if (!g) continue;
	
	if (quad) SqrtGraph(g);
	prev = g;
	drawn.AddFirst(g, FormatOption(h->GetDOption()));
      }
    }
    else { 
      // Otherwise, we are adding all selected errors togethter 
      Graph* g = static_cast<Graph*>(fData->Clone("error"));
      g->SetTitle(fSumTitle);
      Holder* h = 0;
      for (Int_t i = 0; i < n; i++) {
	if (!stat) {
	  g->SetPointError(i, g->GetEXlow()[i],
			   g->GetEXhigh()[i],0,0);
	}
	else if (quad) {
	  // g->GetEXlow()[i] *= g->GetEXlow()[i];
	  g->GetEYlow()[i] *= g->GetEYlow()[i];
	  // g->GetEXhigh()[i] *= g->GetEXhigh()[i];
	  g->GetEYhigh()[i] *= g->GetEYhigh()[i];
	}
	  

	if (cmn) { 
	  TIter nextC(&fCommon);
	  while ((h = static_cast<Holder*>(nextC()))) {
	    h->SumError(g, i, false, quad, fSumOption);
	    
	  }
	}
	TIter nextP(&fPoint2Point);
	while ((h = static_cast<Holder*>(nextP()))) 
	  h->SumError(g, i, false, quad, fSumOption);

	// Printf("Point # %d", i);
	if (quad) SqrtPoint(g, i);
      }
      fSumLine.Copy(*g);
      fSumFill.Copy(*g);
      g->SetMarkerStyle(0);
      g->SetMarkerSize(0);
      g->SetMarkerColor(0);

      drawn.AddFirst(g, FormatOption(fSumOption));
    }
    // --- Show common errors on the side, if requested --------------
    if (!cmn) {
      TIter nextC(&fCommon);
      HolderCommon* hc = 0;
      if (!combine) {
	prev = 0;
	while ((hc = static_cast<HolderCommon*>(nextC()))) {
	  Graph* g = hc->BarError(prev, quad && !split, xBase, yBase);
	  if (!g) {
	    Warning("MakeMulti", "Got no graph for common error %s", 
		    hc->GetTitle());
	    continue;
	  }
	    
	  if (quad && !split) SqrtGraph(g);
	    
	  if (split) xBase += xpos * 1.2 * dx;
	  if (!prev)
	    drawn.AddLast(g, FormatOption(hc->GetDOption()));
	  else {
	    // Here, we use the underling linked list of the list so
	    // that we can add an object before another including the
	    // possible options 
	    TObjLink* f = drawn.FirstLink();
	    TObjLink* p = 0;
	    while (f && f->GetObject()) { 
	      if (f->GetObject() == prev) {
		if (p) new TObjOptLink(g, p, FormatOption(hc->GetDOption()));
		else   new TObjOptLink(g, FormatOption(hc->GetDOption()));
		break;
	      }
	      p = f;
	      f = f->Next();
	    }
	  }
	  prev = (split ? 0 : g);
	} // while (hc)
      } // !combine
      else {
	// Combine systematics 
	if (fCommon.GetEntries() > 0) {
	  Graph* g = new Graph(1);
	  g->SetPoint(0, xBase, yBase);
	  Double_t exl = 0;
	  Double_t exh = 0;
	  Double_t eyl = 0;
	  Double_t eyh = 0;
	  while ((hc = static_cast<HolderCommon*>(nextC()))) {
	    hc->SumPointError(yBase, hc->XMode(fCommonSumOption),
			      false, quad, exl, exh, eyl, eyh);
	  }
	  g->SetPointError(0, exl, exh, eyl, eyh);
	  fCommonSumLine.Copy(*g);
	  fCommonSumFill.Copy(*g);
	  g->SetMarkerStyle(0);
	  g->SetMarkerSize(0);
	  g->SetMarkerColor(0);
	  
	  if (quad) SqrtGraph(g);
	  drawn.AddLast(g, FormatOption(fCommonSumOption));
	} // while hc
      }
    }
    // --- Copy attributes to data copy ------------------------------
    this->TAttLine::Copy(*fData);
    this->TAttFill::Copy(*fData);
    this->TAttMarker::Copy(*fData);

    // --- Get the options for drawing the data ----------------------
    TString dopt = FormatOption(fDataOption);
    dopt.Append(" p same");
    // if (stat && combine) dopt.Append(" X");

    // --- Add copy of data to stack ---------------------------------
    if (!stat) {
      // If stat are not combined into point-to-point, add straight copy 
      drawn.AddLast(fData->Clone("data"), dopt.Data());
    }
    else {
      // Otherwise, we just add the markers (i.e., all errors zero'd)
      TGraph* g = new TGraph(fData->GetN());
      g->SetName("data");
      g->SetTitle(fData->GetTitle());
      this->TAttLine::Copy(*g);
      this->TAttFill::Copy(*g);
      this->TAttMarker::Copy(*g);
      for (Int_t i = 0; i < fData->GetN(); i++)
	g->SetPoint(i, fData->GetX()[i], fData->GetY()[i]);
      drawn.AddLast(g, dopt.Data());
    }

    // --- Now copy objects from temp collection to stack ------------
    TMultiGraph* ret = new TMultiGraph("drawn", GetTitle());
    TIter nextG(&drawn);
    TGraph* gg = 0;
    while ((gg = static_cast<TGraph*>(nextG()))) {
      if (kVerbose & kDraw)
	Printf("%-20s:%-20s %20s "
	       "(line:%3d/%d/%d, marker:%3d/%2d/%4f, fill:%3d/%4d)", 
	       gg->GetName(), gg->GetTitle(), nextG.GetOption(),
	       gg->GetLineColor(),   gg->GetLineStyle(),   gg->GetLineWidth(),
	       gg->GetMarkerColor(), gg->GetMarkerStyle(), gg->GetMarkerSize(),
	       gg->GetFillColor(),   gg->GetFillStyle());
      ret->Add(gg, nextG.GetOption());
    }
    return ret;
  }
  /* @} */
  /** 
   * @{ 
   * @name Internal searches
   */
  //__________________________________________________________________
  /** 
   * Find a point-2-point error graph
   * 
   * @param id identifier
   * 
   * @return Point-to-point systematic error holder or null
   */
  HolderP2P* FindP2P(UInt_t id) const
  {
    TIter next(&fPoint2Point);
    TObject* o = 0;
    while ((o = next())) {
      if (o->GetUniqueID() != id) continue;
      return static_cast<HolderP2P*>(o);
    }
    return 0;
  }
  /** 
   * Find a common error graph
   * 
   * @param id identifier 
   * 
   * @return Common systematic error holder or null
   */
  HolderCommon* FindCommon(UInt_t id) const
  {
    TIter next(&fCommon);
    TObject* o = 0;
    while ((o = next())) {
      if (o->GetUniqueID() != id) continue;
      return static_cast<HolderCommon*>(o);
    }
    return 0;
  }
  /** 
   * Find an error in this graph that is compatible with the passed error 
   * 
   * @param o    Test 
   * @param verb Be verbose
   * @param tol  Relative error tolerance
   * 
   * @return Pointer to holder for found errror or null
   */
  HolderCommon* FindCompat(const HolderCommon* o,
			   Double_t            tol=1e-6,
			   bool                verb=false) const
  {
    if (!o) return 0;
    Int_t id = FindId(o->GetTitle());
    if (id <= 0) {
      if (verb) 
	Warning("FindCompat[C]", "Syst.unc. '%s' not found in %s",
		o->GetTitle(), GetName());
      return 0;
    }
    HolderCommon* c = FindCommon(id);
    if (!c) {
      if (verb)
	Warning("FindCompat[C]", "Syst.unc. '%s' (%d) is not a common",
		o->GetTitle(), id);
      return 0;
    }
    if (o->IsRelative() != c->IsRelative()) {
      if (verb)
	Warning("FindCompat[C]", "Inconsistent value type (%s) of '%s' (%d)",
		c->IsRelative() ? "relative" : "fixed",
		c->GetTitle(), id);
      return 0;
    }
    Double_t ltol = 0.01;
    Double_t x1, x2, test;
    x1   = o->GetYDown(1);
    x2   = c->GetYDown(1);
    test = TMath::Abs(x2-x1) / (1+TMath::Min(TMath::Abs(x1), TMath::Abs(x2)));
    if (test > ltol) {
      if (verb) {
	Warning("FindCompat[C]","Lower value %10f of '%s' (%d) "
		"incompatible with %10f (|%f|>%f [%g])",
		c->GetYDown(1), c->GetTitle(), id, o->GetYDown(1),
		test, ltol, tol);
	c->Print("common");
	o->Print("common");
      }
      return 0;
    }
    x1   = o->GetYDown(1);
    x2   = c->GetYDown(1);
    test = TMath::Abs(x2-x1) / (1+TMath::Min(TMath::Abs(x1), TMath::Abs(x2)));
    if (TMath::Abs(o->GetYUp(1)-c->GetYUp(1)) > ltol) {
      if (verb) {
	Warning("FindCompat[C]","Upper value %10f of '%s' (%d) "
		"incompatible with %10f (|%f|>%f [%g])",
		c->GetYUp(1), c->GetTitle(), id, o->GetYUp(1),
		test, ltol, tol);
	c->Print("common");
	o->Print("common");
      }
      return 0;
    }
    // All is good 
    return c;
  }
  /** 
   * Find an error in this graph that is compatible with the passed error 
   * 
   * @param o    Test 
   * @param verb Be verbose
   * 
   * @return Pointer to holder for found errror or null
   */
  HolderP2P* FindCompat(const HolderP2P* o, Bool_t verb=false) const
  {
    if (!o) return 0;
    Int_t id = FindId(o->GetTitle());
    if (id <= 0) {
      if (verb) 
	Warning("FindCompat[P]", "Syst.unc. '%s' not found in %s",
		o->GetTitle(), GetName());
      return 0;
    }
    HolderP2P* p = FindP2P(id);
    if (!p) {
      if (verb)
	Warning("FindCompat[P]", "Syst.unc. '%s' (%d) is not point-to-point",
		o->GetTitle(), id);
      return 0;
    }
    if (o->IsRelative() != p->IsRelative()) {
      if (verb)
	Warning("FindCompat[P]", "Inconsistent value type (%s) of '%s' (%d)",
		p->IsRelative() ? "relative" : "fixed",
		p->GetTitle(), id);
      return 0;
    }
    return p;
  }
  /** 
   * Find any error
   * 
   * @param id Identifier
   * 
   * @return Holder object or null
   */
  Holder* Find(UInt_t id) const
  {
    Holder* h = FindP2P(id);
    if (h) return h;
    TIter nextC(&fCommon);
    TObject* o = 0;
    while ((o = nextC())) {
      if (o->GetUniqueID() != id) continue;
      return static_cast<Holder*>(o);
    }
    return 0;
  }
  /* @} */
  /** List of graphs */
  TList fPoint2Point;
  /** List of common errors */
  TList fCommon;
  /** Our data points */
  // Graph* fData;
  TGraphAsymmErrors* fData;
  /** The drawn graphs */
  TMultiGraph* fDrawn;
  /** Counter */
  UInt_t fCounter;
  /** Attributes of summed errors */
  TAttFill fSumFill;
  /** Attributes of summed errors */
  TAttLine fSumLine;
  /** Title on summed errors */
  TString fSumTitle;
  // Drawin option for sums 
  UInt_t  fSumOption;
  /** Attributes of summed errors */
  TAttFill fCommonSumFill;
  /** Attributes of summed errors */
  TAttLine fCommonSumLine;
  /** Title on summed errors */
  TString fCommonSumTitle; 
  // Drawin option for sums 
  UInt_t  fCommonSumOption;
  // Drawing options for data
  UInt_t  fDataOption;
  /** X title  */
  TString fXTitle;
  /** Y title */
  TString fYTitle;
  /** Map of keys */
  TList* fMap;
  /** List of qualifiers */
  TList* fQualifiers;
  /** Whether statistical errors are relative */
  Bool_t fStatRelative;
  ClassDef(GraphSysErr,4);
};
/** 
 * @example TestImport.C
 * 
 * Example of simple import of simple dataset.
 */
/** 
 * @example TestMultiImport.C
 *
 * Example of complex import reading all datasets 
 */
/** 
 * @example Example.C
 *
 * Example of using this class 
 */
/** 
 * @example TestCombiner.C
 *
 * Example of linear combiner of observations 
 */
 

#endif
// 
// EOF
// 

