#ifndef QUALITYSELECTION_H
#define QUALITYSELECTION_H

// --- Custom header files ---
#include "PhotonSelection.h"
#include "DetectorHistogram.h"

// --- ROOT system ---
#include <TObjArray.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

class PhotonSpectrumSelection : public PhotonSelection
{
public:
    PhotonSpectrumSelection():
        PhotonSelection(),
        fSpectrum(0),
        fSpectrumCPV(0),
        fSpectrumDisp(0),
        fSpectrumBoth(0)
    {
    }

    PhotonSpectrumSelection(const char * name, const char * title, ClusterCuts cuts, Float_t cpv = 10., Float_t disp = 3.0):
        PhotonSelection(name, title, cuts),
        fDistanceCPV(cpv),
        fDispersionCut(disp),
        fSpectrum(0),
        fSpectrumCPV(0),
        fSpectrumDisp(0),
        fSpectrumBoth(0)
    {
    }

    ~PhotonSpectrumSelection()
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

    PhotonSpectrumSelection(const PhotonSpectrumSelection &);
    PhotonSpectrumSelection & operator = (const PhotonSpectrumSelection &);

    Float_t fDistanceCPV;
    Float_t fDispersionCut;

private:

    DetectorHistogram * fSpectrum;     //!
    DetectorHistogram * fSpectrumCPV;  //!
    DetectorHistogram * fSpectrumDisp; //!
    DetectorHistogram * fSpectrumBoth; //!

    ClassDef(PhotonSpectrumSelection, 2)
};
#endif