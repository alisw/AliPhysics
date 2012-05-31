#ifndef ALIHLTEVEMULTCORR_H
#define ALIHLTEVEMULTCORR_H

#include "AliHLTEveBase.h"
#include "TH1.h"
#include "TCanvas.h"

class AliHLTEveHistoMerger;

class AliHLTEveMultCorr : public AliHLTEveBase
{

public:
    
    /** Constructor */
    AliHLTEveMultCorr(const char* name);
    
    /** Destructor */
    ~AliHLTEveMultCorr();
    
    /** Process block */
    virtual void ProcessBlock(AliHLTHOMERBlockDesc* block);
    
    /** Reset the elements */
    virtual void ResetElements();
    
    /** Reset the elements */
    virtual void UpdateElements();

protected:
    
    virtual void AddHistogramsToCanvas(AliHLTHOMERBlockDesc* block, TCanvas* canvas, Int_t& cdCount);

    virtual void AddHistogramToCanvas(TH1* hist, TCanvas* canvas, Int_t& cdCount, Bool_t zoom = false);
    
    virtual TH1* FindHistogram(TCollection *coll, const char *name);

private:

    TCanvas *fVzeroMultCanvas;
    TCanvas *fZdcMultCanvas;
    TCanvas *fTrackMultCanvas;
    TCanvas *fCorrCanvas;
    TCanvas *fEtCorrCanvas;
    TCanvas *fZdcVzeroSpdCorrCanvas;
 
    AliHLTEveHistoMerger *fMerger;

    TList *fMyList;
    
    /** Default constructor prohibited */
    AliHLTEveMultCorr();
    
    /** copy constructor prohibited */
    AliHLTEveMultCorr(const AliHLTEveMultCorr& );
    
    /** assignment operator prohibited */
    AliHLTEveMultCorr& operator = (const AliHLTEveMultCorr &);
  
    ClassDef(AliHLTEveMultCorr, 0);

};

#endif // ALIHLTEVEMULTCORR_H
