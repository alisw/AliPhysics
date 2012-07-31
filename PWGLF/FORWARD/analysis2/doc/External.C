#error Not for compilation 

/**
 * @defgroup pwglf_external External code 
 *
 * @todo Expand this as needed
 */
//====================================================================
/**
 * @defgroup pwglf_external_root External ROOT code 
 * See also http://root.cern.ch
 * @ingroup pwglf_external
 */
/** 
 * An unsigned 8bit character
 * @ingroup pwglf_external_root
 */
typedef char Char_t;
/** 
 * An unsigned 16bit integer 
 * @ingroup pwglf_external_root
 */
typedef unsigned int UShort_t;
/** 
 * An unsigned 32bit integer 
 * @ingroup pwglf_external_root
 */
typedef unsigned int UInt_t;
/** 
 * An unsigned 32bit integer 
 * @ingroup pwglf_external_root
 */
typedef unsigned long ULong_t;
/** 
 * An constant character for option strings 
 * @ingroup pwglf_external_root
 */
typedef const char Option_t;
/** 
 * An boolean 
 * @ingroup pwglf_external_root
 */
typedef bool Bool_t;
/** 
 * An boolean 
 * @ingroup pwglf_external_root
 */
typedef double Double_t;
/** 
 * An boolean 
 * @ingroup pwglf_external_root
 */
typedef int Int_t;
/** 
 * An boolean 
 * @ingroup pwglf_external_root
 */
typedef float Float_t;
/**
 * Root of the ROOT class hierarcy 
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TObject.html
 * @ingroup pwglf_external_root
 */
class TObject {};
/**
 * Named root of the ROOT class hierarcy 
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TNamed.html
 * @ingroup pwglf_external_root
 */
class TNamed {};
/**
 * A string
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TString.html
 * @ingroup pwglf_external_root
 */
class TString {};
/**
 * An Array
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TArray.html
 * @ingroup pwglf_external_root
 */
class TArray {};
/**
 * An Array
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TArrayI.html
 * @ingroup pwglf_external_root
 */
class TArrayI : public TArray {};
/**
 * An Array
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TArrayD.html
 * @ingroup pwglf_external_root
 */
class TArrayD : public TArray {};
/**
 * A tree of data
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TTree.html
 * @ingroup pwglf_external_root
 */
class TTree : public TNamed {};
/**
 * A chain of trees of data
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TChain.html
 * @ingroup pwglf_external_root
 */
class TChain : public TTree {};
/**
 * Base class for collections 
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TCollection.html
 * @ingroup pwglf_external_root
 */
class TCollection : public TNamed {};
/**
 * A list
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TList.html
 * @ingroup pwglf_external_root
 */
class TList : public TCollection {};
/**
 * An array of objects
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TObjArray.html
 * @ingroup pwglf_external_root
 */
class TObjArray : public TCollection {};
/**
 * A task
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TTask.html
 * @ingroup pwglf_external_root
 */
class TTask : public TNamed {};
/**
 * An axis
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TAxis.html
 * @ingroup pwglf_external_root
 */
class TAxis : public TNamed {};
/**
 * Base classs for histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH1.html
 * @ingroup pwglf_external_root
 */
class TH1 : public TNamed {};
/**
 * 1D histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH1I.html
 * @ingroup pwglf_external_root
 */
class TH1I : public TH1 {};
/**
 * 1D histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH1D.html
 * @ingroup pwglf_external_root
 */
class TH1D : public TH1 {};
/**
 * Base classs for 2D histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH2.html
 * @ingroup pwglf_external_root
 */
class TH2 : public TH1 {};
/**
 * 2D histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH2D.html
 * @ingroup pwglf_external_root
 */
class TH2D : public TH2 {};
/**
 * 2D histograms
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TH2F.html
 * @ingroup pwglf_external_root
 */
class TH2F : public TH2 {};
/**
 * A graph
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TGraph.html
 * @ingroup pwglf_external_root
 */
class TGraph : public TNamed {};
/**
 * A graph
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TGraphErrors.html
 * @ingroup pwglf_external_root
 */
class TGraphErrors : public TGraph {};
/**
 * A graph
 * 
 * See also
 * - http://root.cern.ch/root/htmldoc/TGraphAsymmErrors.html
 * @ingroup pwglf_external_root
 */
class TGraphAsymmErrors : public TGraph {};

//====================================================================
/**
 * @defgroup pwglf_external_aliroot External AliROOT code 
 * See also https://aliceinfo.cern.ch/Offline
 * @ingroup pwglf_external
 */
/**
 * Base class for analsysis tasks.  
 *
 * See also
 * - https://aliceinfo.cern.ch/Offline/Activities/Analysis/
 * - https://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliAnalysisTaskSE.html
 * @ingroup pwglf_external_aliroot
 */
class AliAnalysisTask : public TTask {};
/**
 * Base class for Single Event analsysis tasks.  
 *
 * See also
 * - https://aliceinfo.cern.ch/Offline/Activities/Analysis/
 * - https://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliAnalysisTaskSE.html
 * @ingroup pwglf_external_aliroot
 */
class AliAnalysisTaskSE : public AliAnalysisTask {};
/**
 * AOD event structure 
 * 
 * See also
 * - https://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliAODEvent.html
 * @ingroup pwglf_external_aliroot
 */
class AliAODEvent : public TObject {};
/**
 * ESD event structure 
 * 
 * See also
 * - https://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliESDEvent.html
 * @ingroup pwglf_external_aliroot
 */
class AliESDEvent : public TObject{};
/**
 * FMD ESD structure 
 * 
 * See also
 * - https://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliESDFMD.html
 * @ingroup pwglf_external_aliroot
 */
class AliESDFMD : public TObject{};

//
// EOF
//
