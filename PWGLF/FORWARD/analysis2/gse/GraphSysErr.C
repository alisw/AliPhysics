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
#  include <iostream>
#  include <iomanip>
#  include <fstream>
// #  define TOKEN(C,I) static_cast<TObjString*>(C->At(I))->String()
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
# endif


/** 
 * @page an_example Example 
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
 * And then we draw nad finish
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
 * @see @link an_example Example @endlink
 *
 */
class GraphSysErr : public TNamed, public TAttMarker, 
		    public TAttLine, public TAttFill 
{
public:
  enum { 
    kDraw    = 0x1,
    kImport  = 0x2, 
    kExport  = 0x4, 
    kVerbose = 0  // Set to OR of above to enable verbose/debug output 
  };
  /** A short-hand type definition */
  typedef TGraphAsymmErrors Graph;
  /** 
   * Drawing options.  We re-encode them here as distinct enums. 
   */
  enum EDrawOption_t { 
    kNormal = 0, //     - Line with ticks
    kNoTick = 1, // Z0  - Line with no ticks 
    kArrow  = 2, // >0  - Linw with arrows 
    kRect   = 3, // 20  - Rectangle w/fill
    kBox    = 4, // 50  - Rectangle w/fill & line 
    kFill   = 5, // 30  - Filled area 
    kCurve  = 6, // 40  - Filled smoothed area 
    kHat    = 7, // []0 - Hats 
    kBar    = 8, // ||0 - A bar 
    kNone   = 9  // XP  - No errors
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
      fSumTitle("Errors"),
      fSumOption(0),
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
    while ((common = static_cast<HolderCommon*>(nextC()))) 
      fCommon.Add(new HolderCommon(*common));

    TIter nextP(&other.fPoint2Point);
    HolderP2P* p2p = 0;
    while ((p2p = static_cast<HolderP2P*>(nextP()))) 
      fPoint2Point.Add(new HolderP2P(*p2p));    

    if (other.fMap)
      fMap = static_cast<TList*>(other.fMap->Clone());
    if (other.fQualifiers)
      fQualifiers = static_cast<TList*>(other.fQualifiers->Clone());
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
  virtual void Print(Option_t* option="R") const //*MENU*
  {
    gROOT->IndentLevel();
    std::cout << GetName() << ": " << GetTitle() << std::endl;
    TString opt(option);
    opt.ToUpper();
    if (!opt.Contains("R")) return;
    gROOT->IncreaseDirLevel();

    if (fMap) {
      gROOT->IndentLevel();
      std::cout << "Key/value pairs: " << std::endl;
      gROOT->IncreaseDirLevel();
      fMap->ls(option);
      gROOT->DecreaseDirLevel();
    }
    if (fQualifiers) {
      gROOT->IndentLevel();
      std::cout << "Qualifier pairs: " << std::endl;
      gROOT->IncreaseDirLevel();
      fQualifiers->ls(option);
      gROOT->DecreaseDirLevel();
    }

    gROOT->IndentLevel();
    std::cout << "Commons: " << std::endl;
    gROOT->IncreaseDirLevel();
    fCommon.ls(option);
    gROOT->DecreaseDirLevel();

    gROOT->IndentLevel();
    std::cout <<  "Point-to-point: " << std::endl;
    gROOT->IncreaseDirLevel();
    fPoint2Point.ls(option);
    gROOT->DecreaseDirLevel();

    gROOT->DecreaseDirLevel();
  }
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
   * @dontinclude tests/DrawStyles.C 
   * @skip Drawer
   * @until }
   *
   * A function to set-up an object
   *
   * @skip Test1 
   * @until EndTest1
   *
   * A function to make a canvs 
   *
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
   *
   * - H Export file header 
   * - C Export file comment 
   * - S Export Point-to-point systematic names
   */
  void Export(std::ostream& out=std::cout, Option_t* option="")
  {
    TString opt(option);
    opt.ToLower();
    Bool_t header   = opt.Contains("h");
    Bool_t sysNames = opt.Contains("s");
    Bool_t comment  = opt.Contains("c");
    
    ExportHeader(out, header, comment);
    
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
    if (!hasTitle)
      out << FormatKey("qual") << "RE : " << GetTitle() << std::endl;
    
    // --- Export X/Y titles ----------------------------------------
    const char* fill = "<please fill in>";    
    out << FormatKey("xheader")
	<< (fXTitle.IsNull() ? fill : fXTitle.Data())  << "\n"
	<< FormatKey("yheader")
	<< (fYTitle.IsNull() ? fill : fYTitle.Data())  
	<< std::endl;

    // --- Export common errors --------------------------------------
    TIter nextC(&fCommon);
    HolderCommon* holderCommon = 0;
    while ((holderCommon = static_cast<HolderCommon*>(nextC()))) { 
      Bool_t rel = holderCommon->IsRelative();
      Double_t up = holderCommon->GetYUp(rel ? 100 : 1);
      Double_t down = holderCommon->GetYDown(rel ? 100 : 1);
      out << FormatKey("dserror");
      ExportError(out, down, up, true, rel);
      out << ":" << holderCommon->GetTitle() << std::endl;
    }

    // --- Export data points ----------------------------------------
    out  << FormatKey("data") << " x : y" << std::endl;
    Int_t n = GetN();
    for (Int_t i = 0; i < n; i++) {
      ExportPoint(out, i, true, sysNames);
      out << std::endl;
    }
    out << "*dataend:\n" 
	<< "# End of dataset\n" << std::endl;
  }
  /** 
   * Export a set of data sets to a single table.  All graphs must
   * have the same format.  The title of each graph is written as the
   * "qual" field.
   * 
   * @param col       Collection of GraphSysErr objets 
   * @param out       Output stream 
   * @param option Options
   *
   * - H Export file header 
   * - C Export file comment 
   * - S Export Point-to-point systematic names
   */
  static void Export(const TSeqCollection* col,
		     std::ostream& out,
		     Option_t* option="H")
  {
    if (col->GetEntries() < 1) return;

    // --- Deduce options --------------------------------------------
    TString opt(option);
    opt.ToLower();
    Bool_t alsoTop = opt.Contains("h");
    Bool_t alsoCmt = opt.Contains("c");
    Bool_t alsoNme = opt.Contains("s");

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
	  if (TMath::Abs(ech-ecl) < tol)
	    val = Form("+- %f%s", ecl, (rel ? " PCT" : ""));
	  else
	    val = Form("+%f,-%f%s", ech, ecl, (rel ? " PCT" : ""));
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
    first->ExportHeader(out, alsoTop, alsoCmt);

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
    if (!hasTitle) {
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
	if (val.Contains("#")) {
	  val.Prepend("$");
	  val.Append("$");
	}
	if (!val.IsNull()) val.ReplaceAll("#", "\\");
	else               val = fill;
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
      if (!oCmn->TestBit(BIT(14))) {
	::Warning("Export",
		  "Common systematic error \"%s\" represented by qual",
		  oCmn->GetName());
	continue;
      }
      out << FormatKey("dserror");
      UInt_t   id  = first->FindId(oCmn->GetName());
      Bool_t   rel = first->IsRelative(id);
      Double_t ecl = first->GetCommonErrorYDown(id);
      Double_t ech = first->GetCommonErrorYUp(id);
      ExportError(out, ecl, ech, true, rel);
      TString tit(oCmn->GetName());
      if (tit.Contains("#")) {
	tit.Prepend("$");
	tit.Append("$");
      }
      tit.ReplaceAll("#", "\\");
      out << " : "<< tit << std::endl;
    }

    // --- Export points ---------------------------------------------
    out << FormatKey("data") << data << std::endl;
    for (Int_t i = 0; i < nPoints; i++) {
      TIter  next(&toExport);
      Bool_t one = true;
      while ((gse = static_cast<GraphSysErr*>(next()))) {
	gse->ExportPoint(out, i, one, alsoNme);
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
    GraphSysErr* g = 0;
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
	    TString& nam = Token(tokens, 1);
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
	if (k.EqualTo("RE") ||
	    k.EqualTo("title", TString::kIgnoreCase))
	  ret->SetTitle(q);
      }
      TIter nextK(&keys);
      TObject* pair = 0;
      while ((pair = nextK())) {
	ret->SetKey(pair->GetName(),pair->GetTitle());
      }
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
  /** 
   * Find the ID of an error with the given title 
   * 
   * @param title Title 
   * 
   * @return ID or null
   */
  UInt_t FindId(const char* title) 
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
  Double_t GetStatErrorUp(Int_t point) const { 
    return fData->GetErrorYhigh(point); 
  }
  /** 
   * @param point Point
   * 
   * @return statistical error at point
   */
  Double_t GetStatErrorDown(Int_t point)const{return fData->GetErrorYlow(point);}
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
      fMap->Add(new TNamed("laboratory",  lab.Data()));
      fMap->Add(new TNamed("accelerator", acc.Data()));
      fMap->Add(new TNamed("detector",    exp.Data()));
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
	while (cur != last) {
	  if (!k.EqualTo(cur->GetObject()->GetName())) {
	    cur = cur->fNext;
	    continue;
	  }
	  TObjLink* tmp = cur->fNext;
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
    virtual Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) = 0;
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
			  Bool_t quad, UInt_t opt) = 0;
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

    virtual void Print(Option_t*) const
    {
      gROOT->IndentLevel();
      Printf("%s/%s (ID # %d)", GetName(), GetTitle(), GetUniqueID());
    }
    virtual void ls(Option_t* option) const 
    {
      Print(option);
    }
  protected:
    /** 
     * Add errors together at point
     * 
     * @param i            Point
     * @param g            Current errors
     * @param quad         If true, add in quadrature
     * @param ignoreErr    If true, ignore errors on g
     * @param sqOld        If true and quad true, square old
     * @param exl          On return, the left-hand X errors
     * @param exh          On return, the right-hand X errors
     * @param eyl          On return, the downward Y errors
     * @param eyh          On return, the upward Y errors 
     */
    void DoAdd(Int_t     i, 
	       Graph*    g, 
	       Bool_t    quad, 
	       Bool_t    ignoreErr,
	       Bool_t    sqOld,
	       Double_t& exl, 
	       Double_t& exh, 
	       Double_t& eyl, 
	       Double_t& eyh) 
    {
      if (kVerbose & kDraw) 
	printf("%-20s Point %3d -> %f/%f/%f/%f",GetTitle(),i,exl,exh,eyl,eyh);
      // Double_t oexl = exl;
      // Double_t oexh = exh;
      if (quad) { 
	eyl  *= eyl;
	eyh  *= eyh;
	if (kVerbose & kDraw) printf(" (%f/%f)",eyl,eyh);
      }

      if (ignoreErr) {
	switch (fOption) {
	case kNormal:
	case kNoTick:
	case kFill:
	case kCurve:
	  break;
	case kArrow: // In these cases, do no add errors along X
	case kHat:
	case kBar:
	  exl = 0;
	  exh = 0;
	  // if (oexl <= 0) exl = 0;
	  // if (oexh <= 0) exh = 0;
	case kNone:
	  break;
	case kRect:// In these cases, make sure we have errors along X
	case kBox:
	  if (exl <= 0) exl = gStyle->GetErrorX()/2;
	  if (exh <= 0) exh = gStyle->GetErrorX()/2;
	}
	if (kVerbose & kDraw) Printf("= %f/%f/%f/%f (%d)",exl,exh,eyl,eyh,quad);
	return;
      }
	
      Double_t e2xl = g->GetErrorXlow(i);
      Double_t e2xh = g->GetErrorXhigh(i);
      Double_t e2yl = g->GetErrorYlow(i);
      Double_t e2yh = g->GetErrorYhigh(i);
      if (kVerbose & kDraw) printf("+ %f/%f/%f/%f", e2xl,e2xh,e2yl,e2yh);
      if (quad && sqOld) {
	e2yl *= e2yl;
	e2yh *= e2yh;
	if (kVerbose & kDraw) printf(" (%f/%f)",e2yl,e2yh);
      }

	
      exl  = TMath::Max(exl, e2xl);
      exh  = TMath::Max(exh, e2xh);
      eyl  += e2yl;
      eyh  += e2yh;

      switch (fOption) {
      case kNormal:
      case kNoTick:
      case kFill:
      case kCurve:
	break;
      case kArrow: // In these cases, do no add errors along X
      case kHat:
      case kBar:
	exl = 0;
	exh = 0;
	// if (oexl <= 0) exl = 0;
	// if (oexh <= 0) exh = 0;
      case kNone:
	break;
      case kRect:// In these cases, make sure we have errors along X
      case kBox:
	if (exl <= 0) exl = gStyle->GetErrorX()/2;
	if (exh <= 0) exh = gStyle->GetErrorX()/2;
      }

      if (kVerbose & kDraw) 
	Printf("= %f/%f/%f/%f (%d,%d)",exl,exh,eyl,eyh,quad,sqOld);
    }
    /** 
     * Set attributes
     * 
     * @param g on graph
     */
    void SetAttributes(Graph* g) 
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
    ClassDef(Holder,1);
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
     * @return Symmetric errors along Y at point
     */
    Double_t GetY(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorY(point);
    }
    /** 
     * Get errors downward along Y at point
     * 
     * @param point Point
     * 
     * @return Errors downward along Y at point
     */
    Double_t GetYDown(Int_t point) const 
    {
      if (!fGraph) return 0;
      return fGraph->GetErrorYlow(point);
    }
    /** 
     * Get errors upward along Y at point
     * 
     * @param point Point
     * 
     * @return Errors upward along Y at point
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
     * @param g            Current errors
     * @param quad         If true, add in quadrature
     * @param ignoreErr    If true, ignore errors on g
     * @param sqOld        If true and quad true, square old
     * @param exl          On return, the left-hand X errors
     * @param exh          On return, the right-hand X errors
     * @param eyl          On return, the downward Y errors
     * @param eyh          On return, the upward Y errors 
     */
    void AddError(Int_t     i, 
		  Graph*    g, 
		  Bool_t    quad, 
		  Bool_t    ignoreErr,
		  Bool_t    sqOld,
		  Double_t& exl, 
		  Double_t& exh, 
		  Double_t& eyl, 
		  Double_t& eyh) 
    {
      exl           = fGraph->GetErrorXlow(i);
      exh           = fGraph->GetErrorXhigh(i);
      eyl           = fGraph->GetErrorYlow(i);
      eyh           = fGraph->GetErrorYhigh(i);
      if (exl <= 0) {
	exl = g->GetErrorXlow(i);
      }
      if (exh <= 0) {
	exh = g->GetErrorXhigh(i);
      }

      DoAdd(i, g, quad, ignoreErr, sqOld, exl, exh, eyl, eyh);
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
    Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) 
    {
      Graph* ga = new Graph(g->GetN());
      for (Int_t i = 0; i < g->GetN(); i++) { 
	Double_t x    = g->GetX()[i];	  
	Double_t y    = g->GetY()[i];
	Double_t exl, exh, eyl, eyh;
	AddError(i, g, quad, ignoreErr, true, exl, exh, eyl, eyh);
	ga->SetPoint(i, x, y);
	ga->SetPointError(i, exl,exh,eyl,eyh);
      }
      SetAttributes(ga);
      ga->SetTitle(GetTitle());
      ga->SetName(Form("stack_%08x", GetUniqueID()));
      return ga;
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
    void SumError(Graph* g, Int_t i, Bool_t ignoreErr, Bool_t quad, UInt_t opt) 
    {
      Double_t exl, exh, eyl, eyh;
      UInt_t save = fOption;
      fOption     = opt;
      AddError(i, g, quad, false, ignoreErr, exl, exh, eyl, eyh);
      g->SetPointError(i, exl,exh,eyl,eyh);
      fOption = save;
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
    /** Our data */
    Graph* fGraph;

    ClassDef(HolderP2P,1);
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
     * @param i            Point
     * @param g            Current errors
     * @param quad         If true, add in quadrature
     * @param ignoreErr    If true, ignore errors on g
     * @param sqOld        If true and quad true, square old
     * @param y            Y value at point 
     * @param exl          On return, the left-hand X errors
     * @param exh          On return, the right-hand X errors
     * @param eyl          On return, the downward Y errors
     * @param eyh          On return, the upward Y errors 
     */
    void AddError(Int_t     i, 
		  Graph*    g, 
		  Bool_t    quad, 
		  Bool_t    ignoreErr,
		  Bool_t    sqOld,
		  Double_t  y,
		  Double_t& exl, 
		  Double_t& exh, 
		  Double_t& eyl, 
		  Double_t& eyh) 
    {
      //exl         = (i   == 0        ? 0 :(g->GetX()[i]  -g->GetX()[i-1])/2);
      //exh         = (i+1 >= g->GetN()? 0 :(g->GetX()[i+1]-g->GetX()[i])/2);
      exl = (g ? g->GetErrorXlow(i) : 0);
      exh = (g ? g->GetErrorXhigh(i) : 0);
      if (exl <= 0) exl = gStyle->GetErrorX()/2;
      if (exh <= 0) exh = gStyle->GetErrorX()/2;

      eyl           = GetYDown(y);
      eyh           = GetYUp(y);

      DoAdd(i, g, quad, ignoreErr, sqOld, exl, exh, eyl, eyh);
    }
    /** 
     * Make a graph for showing next to data 
     * 
     * @param p     PRevious errors
     * @param quad  If true, add in quadrature 
     * @param x     Middle X coordinate 
     * @param y     Middle Y coordinate 
     * 
     * @return Newly allocated graph
     */
    Graph* BarError(Graph* p, Bool_t quad, Double_t x, Double_t y)
    {
      Double_t exl, exh, eyl, eyh;
      AddError(0, p, quad, !p, quad, y, exl, exh, eyl, eyh);

      Graph* r = new Graph(1);
      r->SetPoint(0, x, y);
      r->SetPointError(0, exl, exh, eyl, eyh);
      
      SetAttributes(r);
      r->SetTitle(GetTitle());
      r->SetName(Form("bar_%08x", GetUniqueID()));

      return r;
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
    Graph* StackError(Graph* g, Bool_t ignoreErr, Bool_t quad) 
    {
      Graph* ga = new Graph(g->GetN());
      for (Int_t i = 0; i < g->GetN(); i++) { 
	Double_t x    = g->GetX()[i];
	Double_t y    = g->GetY()[i];	
	Double_t exl, exh, eyl, eyh;
	AddError(i, g, quad, ignoreErr, true, y, exl, exh, eyl, eyh);
	ga->SetPoint(i, x, y);
	ga->SetPointError(i, exl,exh,eyl,eyh);
      }
      SetAttributes(ga);
      ga->SetTitle(GetTitle());
      ga->SetName(Form("stack_%08x", GetUniqueID()));
      return ga;
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
		  UInt_t opt) 
    {
      Double_t exl, exh, eyl, eyh;
      UInt_t save = fOption;
      fOption     = opt;
      AddError(i, g, quad, ignoreErr, false, g->GetY()[i], exl, exh, eyl, eyh);
      // Printf("Point %3d -> %f/%f/%f/%f", i, exl,exh,eyl,eyh);
      g->SetPointError(i, exl,exh,eyl,eyh);
      fOption = save;
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
    /** Down errors */
    Double_t fEyl;
    /** Up errors */
    Double_t fEyh;

    ClassDef(HolderCommon,1);
  };
protected:
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
   * 
   * @return output stream 
   */
  std::ostream& ExportHeader(std::ostream& out,
			     Bool_t alsoTop=false,
			     Bool_t alsoComment=false) const
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
    out << "# Start of dataset\n"
	<< FormatKey("dataset") << std::endl;
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
   * 
   * @return Output stream 
   */
  static std::ostream& ExportError(std::ostream& o, Double_t low,
				   Double_t high, Bool_t nopm, Bool_t rel)
  {
    if (TMath::Abs(high-low) < 1e-12)
      return o << (nopm ? "" : "+- ")
	       << TMath::Abs(low) << (rel ? " PCT" : "");
    return o << "+" << high << ",-" << low << (rel ? " PCT" : "");
  }
  /** 
   * Export a single point 
   * 
   * @param out     Output stream 
   * @param i       Point number
   * @param alsoX   If true, also export X coordinate 
   * @param sysName If true, export P2P names
   * 
   * @return output stream 
   */
  std::ostream& ExportPoint(std::ostream& out, Int_t i,
			    Bool_t alsoX=true, Bool_t sysName=true) const
  {
    if (alsoX) {
      Double_t x       = GetX(i);
      Double_t exl     = fData->GetErrorXlow(i);
      Double_t exh     = fData->GetErrorXhigh(i);
      if (TMath::Abs(exl) < 1e-10 &&
	  TMath::Abs(exh) < 1e-10)
	out << ' ' << x;
      else if (TMath::Abs(exh-exl) < 1e-10)
	out << ' ' << x-exl << " TO " << x+exh;
      else 
	out << ' ' << x << ' '
	    << '+' << exh << ",-" << exl;
      // ExportError(out, exl, exh, false, false);
      out << "; ";
    }
	
    Bool_t   statRel = fStatRelative;
    Double_t y       = GetY(i);
    Double_t fy      = (statRel ? (y == 0 ? 0 : 100./y) : 1);
    Double_t eyl     = GetStatErrorDown(i) * fy;
    Double_t eyh     = GetStatErrorUp(i)   * fy;
    out << y << ' ';
    ExportError(out, eyl, eyh, false, statRel);
      
    if (fPoint2Point.GetEntries() > 0) {
      out  << " (" << std::flush;
      TIter       nextP(&fPoint2Point);
      HolderP2P*  holderP2P = 0;
      Bool_t      first = true;
      while ((holderP2P = static_cast<HolderP2P*>(nextP()))) { 
	Bool_t rel = holderP2P->IsRelative();
	fy         = (rel ? (y == 0 ? 0 : 100./y) : 1);
	if (!first) out << ',';
	Double_t esh = holderP2P->GetYUp(i)*fy;
	Double_t esl = holderP2P->GetYDown(i)*fy;
	out << "DSYS=";
	ExportError(out, esl, esh, true, rel);
	if (sysName) 
	  out << ':'  << holderP2P->GetTitle() << std::flush;
	first = false;
      }
      out << ')';
    }
    out << ';';
    return out;
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
    case kNormal: ret = "";     break; // Line with ticks
    case kNoTick: ret = "Z0";	break; // Line with no ticks 
    case kArrow:  ret = ">0";	break; // Linw with arrows 
    case kRect:   ret = "20";	break; // Rectangle w/fill
    case kBox:    ret = "50";	break; // Rectangle w/fill & line 
    case kFill:   ret = "30";	break; // Filled area 
    case kCurve:  ret = "40";	break; // Filled smoothed area 
    case kHat:    ret = "[]0";	break; // Hats 
    case kBar:    ret = "||0";	break; // A bar 
    case kNone:   ret = "XP";	break; // No errors
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
   * Make our stack 
   *
   * @param option Options - See method Draw for more
   *
   * @return Pointer to stack 
   */
  TMultiGraph* MakeMulti(Option_t* option)
  {

    // --- Process options -------------------------------------------
    TString opt(option);
    opt.ToUpper();

    Bool_t   cmn     = opt.Contains("COMMON");
    Bool_t   combine = opt.Contains("COMBINED");
    Bool_t   stat    = opt.Contains("STAT");
    Bool_t   quad    = opt.Contains("QUAD");
    combine          = !opt.Contains("STACK");
    quad             = !opt.Contains("DIRECT");
    Bool_t   split   = opt.Contains("SPLIT");
    Int_t    xpos    = (opt.Contains("WEST") ? -1 : 1);
    Int_t    ypos    = (opt.Contains("MIN") ? -1 : 
			opt.Contains("MAX") ?  1 : 0);

    // --- Find location for common errors ---------------------------
    Int_t    n       = fData->GetN();
    Double_t dx      = TMath::Max(gStyle->GetErrorX(),0.0001F);
    Double_t xBase   = ((xpos < 0 ? 
			 fData->GetX()[0]   - fData->GetEXlow()[0]: 
			 fData->GetX()[n-1] + fData->GetEXhigh()[n-1])
			+ xpos * 1.2 * dx);
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
	  g->GetEXlow()[i] *= g->GetEXlow()[i];
	  g->GetEYlow()[i] *= g->GetEYlow()[i];
	  g->GetEXhigh()[i] *= g->GetEXhigh()[i];
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
    drawn.AddLast(fData->Clone("data"), dopt.Data());

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
  Graph* fData;
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
  ClassDef(GraphSysErr,1);
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
 

#endif
// 
// EOF
// 

