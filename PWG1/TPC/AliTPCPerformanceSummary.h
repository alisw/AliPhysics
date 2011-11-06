#ifndef ALI_TPC_PERFORMANCE_SUMMARY_H
#define ALI_TPC_PERFORMANCE_SUMMARY_H

//------------------------------------------------------------------------------
// Class to extract some TPC Performance parameters from AliPerformanceTPC and
// AliPerformanceDEdx objects and produce trend graphs.  
// 
// by M.Knichel 15/10/2010 
//------------------------------------------------------------------------------

class TTree;

class TTreeSRedirector;
class AliPerformanceTPC;
class AliPerformanceDEdx;
class AliPerformanceDCA;
class AliPerformanceMatch;

class AliTPCPerformanceSummary
{
    public:
    AliTPCPerformanceSummary() {} // default contructor 
    virtual ~AliTPCPerformanceSummary() {} // destructor
    
    static void WriteToTTreeSRedirector(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const AliPerformanceMatch* pTPCMatch, const AliPerformanceMatch* pTPCPull, const AliPerformanceMatch* pConstrain, TTreeSRedirector* const pcstream, Int_t run = -1); // called by WriteToFile
    
    static void WriteToFile(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const AliPerformanceMatch* pMatch,const AliPerformanceMatch* pPull, const AliPerformanceMatch* pConstrain, const Char_t* outfile, Int_t run = -1); // calles by MakeReport
    
    // the two key functions
    static Int_t MakeReport(const Char_t* infile, const Char_t* outfile, Int_t run);
    static Int_t ProduceTrends(const Char_t* infilelist, const Char_t* outfile);
    
    static Bool_t GetForceTHnSparse() { return fgForceTHnSparse; }
    static void SetForceTHnSparse(Bool_t forceSparse = kTRUE) { fgForceTHnSparse = forceSparse; }      
  
    private:
    
    static Bool_t fgForceTHnSparse;   // force to use THnSparse 
    // save graphs to current directory
    
    static Int_t SaveGraph(TTree* tree, const Char_t* y, const Char_t* x, const Char_t* condition);
    
    // helper functions to extract parameter and write to TTreeSRedirector
    static Int_t AnalyzeDCARPhi(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeDCARPhiPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeDCARPhiNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeNCL(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeDrift(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeDriftPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeDriftNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeGain(const AliPerformanceDEdx* pTPCgain, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeEvent(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream);
    static Int_t AnalyzeMatch(const AliPerformanceMatch* pMatch, TTreeSRedirector* const pcstream);
    
    static Int_t AnalyzePull(const AliPerformanceMatch* pPull, TTreeSRedirector* const pcstream);
    
    static Int_t AnalyzePt(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);

    static Int_t AnalyzeChargeOverPt(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);

    static Int_t AnalyzeConstrain(const AliPerformanceMatch* pConstrain, TTreeSRedirector* pcstream);

    AliTPCPerformanceSummary(const AliTPCPerformanceSummary&); // copy contructor (not implemented)
    AliTPCPerformanceSummary& operator=(const AliTPCPerformanceSummary&); // assignment operator (not implemented)
      
    ClassDef(AliTPCPerformanceSummary, 3);
};

#endif
