#ifndef ALIHLTEVEMULTCORR_H
#define ALIHLTEVEMULTCORR_H

#include "AliHLTEveBase.h"
#include "TH1.h"
#include "TCanvas.h"

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

    virtual void AddHistogramToCanvas(TH1* block, TCanvas* canvas, Int_t& cdCount);

private:

    TCanvas *fVzeroMultCanvas;
    TCanvas *fZdcMultCanvas;
    TCanvas *fTpcMultCanvas;
    TCanvas *fCorrCanvas;
    TCanvas *fEtCorrCanvas;
    TCanvas *fZdcVzeroSpdCorrCanvas;
 
    /** Default constructor prohibited */
    AliHLTEveMultCorr();
    
    /** copy constructor prohibited */
    AliHLTEveMultCorr(const AliHLTEveMultCorr& );
    
    /** assignment operator prohibited */
    AliHLTEveMultCorr& operator = (const AliHLTEveMultCorr &);
  
    ClassDef(AliHLTEveMultCorr, 0);

};

#endif // ALIHLTEVEMULTCORR_H
