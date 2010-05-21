#ifndef ALI_TPC_PERFORMANCE_SUMMARY_H
#define ALI_TPC_PERFORMANCE_SUMMARY_H

//------------------------------------------------------------------------------
// Class to extract some TPC Performance parameters from AliPerformanceTPC and
// AliPerformanceDEdx objects and produce trend graphs.  
// 
// Author: M.Knichel 2010-05-21
//------------------------------------------------------------------------------

class TTree;

class TTreeSRedirector;
class AliPerformanceTPC;
class AliPerformanceDEdx;
class AliPerformanceDCA;

class AliTPCPerformanceSummary
{
    public:
    AliTPCPerformanceSummary() {} // default contructor 
    virtual ~AliTPCPerformanceSummary() {} // destructor
    
    static Int_t WriteToTTreeSRedirector(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, TTreeSRedirector* pcstream, Int_t run); // called by WriteToFile
    
    static Int_t WriteToFile(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const Char_t* outfile, Int_t run); // calles by MakeReport
    
    // the two key functions
    static Int_t MakeReport(const Char_t* infile, const Char_t* outfile, Int_t run);    
    static Int_t ProduceTrends(const Char_t* infilelist, const Char_t* outfile);
  
    private:    
    // save graphs to current directory
    static Int_t SaveGraph(TTree* tree, const Char_t* y, const Char_t* x, const Char_t* condition);
    
    // helper functions to extract parameter and write to TTreeSRedirector
    static Int_t AnalyzeDCARPhi(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeDCARPhiPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeDCARPhiNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeNCL(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeDrift(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeDriftPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeDriftNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
    static Int_t AnalyzeGain(const AliPerformanceDEdx* pTPCgain, TTreeSRedirector* pcstream);
    static Int_t AnalyzeEvent(const AliPerformanceTPC* pTPC, TTreeSRedirector* pcstream);
      
    AliTPCPerformanceSummary(const AliTPCPerformanceSummary&); // copy contructor (not implemented)
    AliTPCPerformanceSummary& operator=(const AliTPCPerformanceSummary&); // assignment operator (not implemented)
      
    ClassDef(AliTPCPerformanceSummary, 1);
};

#endif
