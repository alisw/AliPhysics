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
        void SetdEdxcutLabel(int dEdxcutLabel);

        // set Nsigma 
        void SetProtonNsigma(float Nsigma);
        void SetDeuteronNsigma(float Nsigma);
        void SetTritonNsigma(float Nsigma);
        void SetHe3Nsigma(float Nsigma);
        void SetRejectionNsigma(float Nsigma);
	void SetOtherNsigmacutLabel(int OtherNsigmaLabel);
	void SetMinTPCFoundFraction(float MinTPCFoundFraction);
	// set TOF Mass cut for side band
	void SetUseTOFMassCut(int UseTOFMassCut);
        void SetTOFMassLimit(float LowMass,float UpMass);
	float ReturnTOFMass(const AliFemtoTrack* track,int imost);
	void SetfUsePtCut(int aUsePtCut);

	void SetfOnlyTPCreject(int aOnlyTPCreject);
	
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
    
        int fdEdxcut;
	// 11.30
	int fOtherNsigmacut;
	// 2022.1.20
	float fMinTPCFoundFraction;
        
        bool IsProtonNSigma(    float mom, float nsigmaTPCP, float nsigmaTOFP);
        bool IsDeuteronNSigma(  float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCD, float nsigmaTOFD);
        bool IsTritonNSigma(    float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCT, float nsigmaTOFT);
        bool IsHe3NSigma(       float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCHe3, float nsigmaTOFHe3);
        //\ for e+e femto
        bool IsElectronNSigma(float mom, float nsigmaTPCE, float nsigmaTOFE);
        // dE/dx
        bool IsDeuteronTPCdEdx(float mom, float dEdx);
        bool IsTritonTPCdEdx(float mom, float dEdx);
        // reject
        bool IsElectronNSigmaRejection( float mom, float nsigmaTPCE);
        bool IsPionNSigmaRejection(     float mom, float nsigmaTPCPi, float nsigmaTOFPi);
        bool IsKaonNSigmaRejection(     float mom, float nsigmaTPCK, float nsigmaTOFK);
        bool IsProtonNSigmaRejection(   float mom, float nsigmaTPCP, float nsigmaTOFP);
        
	//\ 2022.1.26
	int fUseTOFMassCut;
	float TOFMassLowLimit;
	float TOFMassUpLimit;
     
	//
	int fUsePtCut;
        // for Nsigma reject
	int fOnlyTPCreject;
	
	// real reject, now, only for p
	bool RejectFakeP(const AliFemtoTrack* track, float mom);
	bool RejectFakeD(const AliFemtoTrack* track, float mom);

};
inline void AliFemtoTrackCutPdtHe3::SetMostProbableDeuteron() { fMostProbable = 13; }
inline void AliFemtoTrackCutPdtHe3::SetMostProbableTriton() { fMostProbable = 14; }
inline void AliFemtoTrackCutPdtHe3::SetMostProbableHe3() { fMostProbable = 15; }

#endif
