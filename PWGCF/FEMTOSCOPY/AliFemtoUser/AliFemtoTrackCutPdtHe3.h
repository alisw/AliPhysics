#ifndef AliFemtoTrackCutPdtHe3_hh
#define AliFemtoTrackCutPdtHe3_hh

#include "AliFemtoESDTrackCut.h"

//==============================================================\\
// dowang track cut for p-d/t/He3 analysis                      \\
// deuteron part refer to AliFemtoWRzTrackCut.h                 \\
//==============================================================\\

class AliFemtoTrackCutPdtHe3 : public AliFemtoESDTrackCut{

    public:
        AliFemtoTrackCutPdtHe3();
        AliFemtoTrackCutPdtHe3(const AliFemtoTrackCutPdtHe3 &aCut);
        virtual ~AliFemtoTrackCutPdtHe3();
        AliFemtoTrackCutPdtHe3& operator =(const AliFemtoTrackCutPdtHe3 &aCut);
        virtual bool Pass(const AliFemtoTrack* aTrack);
        // label
        void SetMostProbableDeuteron();
        void SetMostProbableTriton();
        void SetMostProbableHe3();
	
	    // at mom, need TOF
	    void SetProtonSwitchMom(float SwitchMom);
        void SetDeuteronSwitchMom(float SwitchMom);
        void SetTritonSwitchMom(float SwitchMom);
        void SetHe3SwitchMom(float SwitchMom);
        void SetMostProbableElectron();

    private:
        float fNsigmaP;
        float fNsigmaD;
        float fNsigmaT;
        float fNsigmaHe3;

        float fNsigmaRejection;
        
        // 10.24
        float SwitchMom_p;
        float SwitchMom_d;
        float SwitchMom_t;
        float SwitchMom_He3;
    

        
        bool IsProtonNSigma(    float mom, float nsigmaTPCP, float nsigmaTOFP);
        bool IsDeuteronNSigma(  float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCD, float nsigmaTOFD);
        bool IsTritonNSigma(    float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCT, float nsigmaTOFT);
        bool IsHe3NSigma(       float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCHe3, float nsigmaTOFHe3);
        //\ for e+e femto
        bool IsElectronNSigma(float mom, float nsigmaTPCE, float nsigmaTOFE);
        // dE/dx
        bool IsDeuteronTPCdEdx(float mom, float dEdx, float maxmom);

        // reject
        bool IsElectronNSigmaRejection( float mom, float nsigmaTPCE);
        bool IsPionNSigmaRejection(     float mom, float nsigmaTPCPi, float nsigmaTOFPi);
        bool IsKaonNSigmaRejection(     float mom, float nsigmaTPCK, float nsigmaTOFK);
        bool IsProtonNSigmaRejection(   float mom, float nsigmaTPCP, float nsigmaTOFP);
        
       

};
inline void AliFemtoTrackCutPdtHe3::SetMostProbableDeuteron() { fMostProbable = 13; }
inline void AliFemtoTrackCutPdtHe3::SetMostProbableTriton() { fMostProbable = 14; }
inline void AliFemtoTrackCutPdtHe3::SetMostProbableHe3() { fMostProbable = 15; }

#endif
