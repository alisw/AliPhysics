#ifndef AliFemtoTrackCutPdtHe3_hh
#define AliFemtoTrackCutPdtHe3_hh

#include "AliFemtoESDTrackCut.h"


// wdf 2020.8.17
// deuteron part refer to AliFemtoWRzTrackCut.h

class AliFemtoTrackCutPdtHe3 : public AliFemtoESDTrackCut{

    public:
        AliFemtoTrackCutPdtHe3();
        AliFemtoTrackCutPdtHe3(const AliFemtoTrackCutPdtHe3 &aCut);
        virtual ~AliFemtoTrackCutPdtHe3();
        AliFemtoTrackCutPdtHe3& operator =(const AliFemtoTrackCutPdtHe3 &aCut);
        virtual bool Pass(const AliFemtoTrack* aTrack);

    private:
        float fNsigmaP;
        float fNsigmaD;
        float fNsigmaT;
        float fNsigmaHe3;

        float fNsigmaRejection;

        bool IsProtonNSigma(float mom, float fNsigma, float nsigmaTPCP, float nsigmaTOFP);
        bool IsDeuteronNSigma(float mom, float fNsigma, float massTOFPDG, float sigmaMass, float nsigmaTPCD, float nsigmaTOFD);
        bool IsTritonNSigma(float mom, float fNsigma, float massTOFPDG, float sigmaMass, float nsigmaTPCT, float nsigmaTOFT);
        bool IsHe3NSigma(float mom, float fNsigma, float massTOFPDG, float sigmaMass, float nsigmaTPCHe3, float nsigmaTOFHe3);
        // dE/dx
        bool IsDeuteronTPCdEdx(float mom, float dEdx, float maxmom);

        // reject
        bool IsElectronNSigmaRejection(float mom, float nsigmaTPCE);
        bool IsPionNSigmaRejection(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
        bool IsKaonNSigmaRejection(float mom, float nsigmaTPCK, float nsigmaTOFK);
        bool IsProtonNSigmaRejection(float mom, float nsigmaTPCP, float nsigmaTOFP);
        

};

#endif
