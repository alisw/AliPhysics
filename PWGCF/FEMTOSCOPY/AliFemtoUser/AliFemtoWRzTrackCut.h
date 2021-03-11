///////////////////////////////////////////////////////////////////////////
// wioleta.rzesa@cern.ch                                                 //
// AliFemtoWRzTrackCut: An extension to AliFemtoESDTrackCut              //
// dedicated to analysis with deuterons in PbPb                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoWRzTrackCut_hh
#define AliFemtoWRzTrackCut_hh

#include "AliFemtoESDTrackCut.h"

class AliFemtoWRzTrackCut : public AliFemtoESDTrackCut{

public:
   AliFemtoWRzTrackCut();
   AliFemtoWRzTrackCut(const AliFemtoWRzTrackCut &aCut);
   virtual ~AliFemtoWRzTrackCut();
   AliFemtoWRzTrackCut& operator =(const AliFemtoWRzTrackCut &aCut);
   virtual bool Pass(const AliFemtoTrack* aTrack);

   void SetdEdxcut(bool coumet = true);
   void SetMaxmom(float maxm = 4.0);
   void SetMostProbableDeuteron();
   void SetNsigmaRejection(float nsigmaRejection = 3.0);
   void SetfNSigmaMass(float nsigmass = -1);
   void SetBasicSel(bool fBasicSelection = true);

private:
   bool  fdEdxcut;         // true - if nsigma selection plus the Contour Method 
   float fmaxmom;           // max momentum - for dEdx cut
   float fNsigmaRejection; // nSigmafor the EXCLUSIVE PID cut
   float fNSigmaMass;
   bool  fBasicSelection;  //true- simple analysis, false - analysis with addit. constraints in tpc only selection

   bool IsKaonNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF);
   bool IsPionNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF);
   bool IsProtonNSigmaRejection(float mom, float nsigmaTPC, float nsigmaTOF);
   bool IsElectronNSigmaRejection(float mom, float nsigmaTPE);
   bool IsDeuteronTPCdEdx(float mom, float dEdx, float maxmom);
   bool IsDeuteronNSigma(float mom, float massTOFPDG,float sigmaMass, float nsigmaTPCD, float nsigmaTOFD);

};

inline void AliFemtoWRzTrackCut::SetdEdxcut(bool coumet) { fdEdxcut = coumet; }
inline void AliFemtoWRzTrackCut::SetMaxmom(float maxm) { fmaxmom = maxm; }
inline void AliFemtoWRzTrackCut::SetMostProbableDeuteron() { fMostProbable = 13; }
inline void AliFemtoWRzTrackCut::SetNsigmaRejection(float nsigmaRejection) { fNsigmaRejection = nsigmaRejection; }
inline void AliFemtoWRzTrackCut::SetfNSigmaMass(float nsigmass) { fNSigmaMass = nsigmass; }
inline void AliFemtoWRzTrackCut::SetBasicSel(bool basSel) { fBasicSelection = basSel; }
#endif
