#ifndef ALIPP13QUALITYSELECTION_H
#define ALIPP13QUALITYSELECTION_H

// --- Custom header files ---
#include "AliPP13PhotonSelection.h"
#include "AliPP13DetectorHistogram.h"

// --- ROOT system ---
#include <TObjArray.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

class AliPP13PhotonSpectrumSelection : public AliPP13PhotonSelection
{
public:
    AliPP13PhotonSpectrumSelection():
        AliPP13PhotonSelection(),
        fSpectrum(0),
        fSpectrumCPV(0),
        fSpectrumDisp(0),
        fSpectrumBoth(0)
    {
    }

    AliPP13PhotonSpectrumSelection(const char * name, const char * title, AliPP13ClusterCuts cuts, Float_t cpv = 10., Float_t disp = 3.0):
        AliPP13PhotonSelection(name, title, cuts),
        fDistanceCPV(cpv),
        fDispersionCut(disp),
        fSpectrum(0),
        fSpectrumCPV(0),
        fSpectrumDisp(0),
        fSpectrumBoth(0)
    {
    }

    ~AliPP13PhotonSpectrumSelection()
    {
        if(fSpectrum)     delete fSpectrum;
        if(fSpectrumCPV)  delete fSpectrumCPV;
        if(fSpectrumDisp) delete fSpectrumDisp;
        if(fSpectrumBoth) delete fSpectrumBoth;
    }

    virtual void InitSelectionHistograms();

    // There is no need to mix events
    virtual void MixPhotons(TObjArray & photons, TList * pool, const EventFlags & eflags)
    {
        (void) photons;
        (void) pool;
        (void) eflags;
    }

protected:
    virtual void FillClusterHistograms(const AliVCluster * c, const EventFlags & eflags);

    AliPP13PhotonSpectrumSelection(const AliPP13PhotonSpectrumSelection &);
    AliPP13PhotonSpectrumSelection & operator = (const AliPP13PhotonSpectrumSelection &);

    Float_t fDistanceCPV;
    Float_t fDispersionCut;

private:

    AliPP13DetectorHistogram * fSpectrum;     //!
    AliPP13DetectorHistogram * fSpectrumCPV;  //!
    AliPP13DetectorHistogram * fSpectrumDisp; //!
    AliPP13DetectorHistogram * fSpectrumBoth; //!

    ClassDef(AliPP13PhotonSpectrumSelection, 2)
};
#endif