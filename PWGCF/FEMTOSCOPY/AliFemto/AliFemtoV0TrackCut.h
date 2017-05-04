/// \class AliFemtoV0TrackCut

#ifndef ALIFEMTOV0TRACKCUT_H
#define ALIFEMTOV0TRACKCUT_H

#include "AliFemtoTrackCut.h"

#include "TH1D.h"

/**
 * \class AliFemtoV0TrackCut
 * \brief A track cut designed to cut on V0 (i.e. Lambda) particles.
 *
 * A particle cut object which tests AliFemtoV0 obects for a number of
 * conditions for the reconstructed V0 and the positive and negative
 * (measured) daughter tracks.
 *
 * NOTE: The class, AliFemtoV0TrackCutNSigmaFilter in the AliFemtoUser directory
 *       makes it much easier to customize the NSigma cuts used in the V0 finder
 */
class AliFemtoV0TrackCut : public AliFemtoParticleCut {
public:
  /// Enumerated type to easily select correct algorithm for each particle type
  enum V0Type {
    kLambda = 0,
    kAntiLambda = 1,
    kK0s = 2,
    kAll = 99,
    kLambdaMC = 101,
    kAntiLambdaMC = 102,
    kK0sMC = 3
  };
  typedef enum V0Type AliFemtoV0Type;

  AliFemtoV0TrackCut();
  virtual ~AliFemtoV0TrackCut();

  AliFemtoV0TrackCut(const AliFemtoV0TrackCut& aCut);
  AliFemtoV0TrackCut& operator=(const AliFemtoV0TrackCut& aCut);
  virtual AliFemtoV0TrackCut* Clone();

  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtV0;}

  void SetInvariantMassLambda(double,double);
  void SetInvariantMassK0s(double,double);

  void SetMinDaughtersToPrimVertex(double,double);
  void SetMaxDcaV0Daughters(double);
  void SetMaxDcaV0(double);
  void SetMinDcaV0(double);
  void SetMaxCosPointingAngle(double);
  void SetMinCosPointingAngle(double);
  void SetMaxV0DecayLength(double);
  void SetMinTransverseDistancePrimSecVtx(double);
  void SetParticleType(short);
  void SetEta(double);
  void SetPt(double,double);
  void SetEtaDaughters(float);
  void SetTPCnclsDaughters(int);
  void SetNdofDaughters(int);
  void SetStatusDaughters(unsigned long);
  void SetPtPosDaughter(float,float);
  void SetPtNegDaughter(float,float);
  void SetOnFlyStatus(bool);
  void SetMinAvgSeparation(double);

  void SetRadiusV0Min(double);
  void SetRadiusV0Max(double);

  void SetRequireTOFPion(bool);
  void SetRequireTOFProton(bool);
  
  void SetNsigmaPosDaughter(double);
  void SetNsigmaNegDaughter(double);
  void SetNsigmaPosDaughter(double,double);
  void SetNsigmaNegDaughter(double,double);
  

  //----n sigma----
  bool IsKaonTPCdEdxNSigma(float mom, float nsigmaK);
  bool IsKaonTOFNSigma(float mom, float nsigmaK);
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=false);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=true);

  //-----The fMinvPurityAidHistoV0 is built immediately before the (final) invariant mass cut, and thus may be used to calculate the purity of the V0 collection
  void SetMinvPurityAidHistoV0(const char* name, const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax);  //set the Minv histogram attributes and 
                                                                                                                                            //automatically sets flag fBuildPurityAidV0=true
  TH1D* GetMinvPurityAidHistoV0();  // If initiated, fMinvPurityAidHistoV0 will be in the output list after all other V0 cut monitors (pass and fail)
  void SetLooseInvMassCut(bool aUseCut, double aInvMassMin, double aInvMassMax);

  //--Members to help cut out misidentified V0s
  //NOTE: For MUCH more control over the misidentification cut, use AliFemtoV0TrackCutNSigmaFilter class
  void SetRemoveMisidentified(bool aRemove);
  void SetUseSimpleMisIDCut(bool aUse);

  void SetInvariantMassRejectK0s(double,double, bool aRemoveMisidentified=true);         // invariant mass window in which to reject misidentified K0s from (Anti)Lambda
  void SetInvariantMassRejectLambda(double,double, bool aRemoveMisidentified=true);      // invariant mass window in which to reject misidentified Lambda from K0s
  void SetInvariantMassRejectAntiLambda(double,double, bool aRemoveMisidentified=true);  // invariant mass window in which to reject misidentified AntiLambda from K0s
  void SetInvMassReject(AliFemtoV0Type aV0Type, double aInvMassMin, double aInvMassMax, bool aRemoveMisidentified=true);  //aV0Type selects one of the above setters

  void SetBuildMisIDHistograms(bool aBuild);
  void SetMisIDHisto(AliFemtoV0Type aMisIDV0Type, const int& nbins, const float& aInvMassMin, const float& aInvMassMax);  //allows control over binning
  void SetDefaultMisIDHistos();
  TObjArray *GetMisIDHistos();

  bool IsMisIDK0s(const AliFemtoV0* aV0);
  bool IsMisIDLambda(const AliFemtoV0* aV0);
  bool IsMisIDAntiLambda(const AliFemtoV0* aV0);

  virtual TList *GetOutputList();  //include fMinvPurityAidHistoV0 and fK0sMassOfMisIDV0 etc. in the output list 

  bool GetBuildMisIDHistograms();

  void SetIgnoreOnFlyStatus(bool aIgnore);

 protected:   // here are the quantities I want to cut on...

  double fInvMassLambdaMin;        ///< invariant mass Lambda min
  double fInvMassLambdaMax;        ///< invariant mass Lambda max
  double fInvMassK0sMin;           ///< invariant mass K0s min
  double fInvMassK0sMax;           ///< invariant mass K0s max
  
  double fMinDcaDaughterPosToVert; ///< DCA of positive daughter to primary vertex
  double fMinDcaDaughterNegToVert; ///< DCA of negative daughter to primary vertex
  double fMaxDcaV0Daughters;       ///< Max DCA of v0 daughters at Decay vertex
  double fMaxDcaV0;
  double fMinDcaV0;
  double fMaxDecayLength;
  double fMinTransverseDistancePrimSecVtx;
  
  double fMaxCosPointingAngle;    //obsolete
  double fMinCosPointingAngle;    //correct
  short fParticleType;             ///< 0-lambda
  double fEta;
  double fPtMin;
  double fPtMax;
  bool fOnFlyStatus;

  float fMaxEtaDaughters;         ///< Eta of positive daughter
  int fTPCNclsDaughters;          ///< No. of cls of pos daughter
  int fNdofDaughters;             ///< No. of degrees of freedom of the pos. daughter track
  unsigned long fStatusDaughters; ///< Status (tpc refit, its refit...)
  float fPtMinPosDaughter;
  float fPtMaxPosDaughter;
  float fPtMinNegDaughter;
  float fPtMaxNegDaughter;
  double fMinAvgSepDaughters;

  double fNsigmaPosDaughterTPC;
  double fNsigmaNegDaughterTPC;
  double fNsigmaPosDaughterTOF;
  double fNsigmaNegDaughterTOF;

  double fRadiusV0Min;
  double fRadiusV0Max;

  bool fRequireTOFPion;
  bool fRequireTOFProton;
  
  bool fBuildPurityAidV0;
  TH1D* fMinvPurityAidHistoV0;

  bool fUseLooseInvMassCut;  //Since inv. mass cut must come last in order to calculate the purity, 
                             //this allows a looser cut to be implemented earlier in the process and same some time
  double fLooseInvMassMin;
  double fLooseInvMassMax;

  //--Members to help cut out misidentified V0s
  //NOTE: For MUCH more control over the misidentification cut, use AliFemtoV0TrackCutNSigmaFilter class
  bool fRemoveMisidentified;         // attempt to remove V0 candidates (K0s, Lambda and AntiLambda) which are misidentified
                                     //   i.e. check if Lambda or AntiLambda could be misidentified K0s,
                                     //        or if K0s could be misidentified Lambda or AntiLambda
                                     // Uses methods IsMisIDK0s, IsMisIDLambda, IsMisIDAntiLambda
  bool fUseSimpleMisIDCut;	     // If set to true, the misidentification cut is based soley off of the invariant mass
                                     // If set to false, it also considers the NSigma of the daughter particles when cutting
  bool fBuildMisIDHistograms;

  double fInvMassRejectK0sMin;           ///< invariant mass min to reject misidentified K0s from (Anti)Lambda
  double fInvMassRejectK0sMax;           ///< invariant mass max to reject misidentified K0s from (Anti)Lambda
  double fInvMassRejectLambdaMin;        ///< invariant mass min to reject misidentified Lambda from K0s
  double fInvMassRejectLambdaMax;        ///< invariant mass min to reject misidentified Lambda from K0s
  double fInvMassRejectAntiLambdaMin;    ///< invariant mass min to reject misidentified AntiLambda from K0s
  double fInvMassRejectAntiLambdaMax;    ///< invariant mass min to reject misidentified AntiLambda from K0s

  TH1D *fK0sMassOfMisIDV0;             // Mass asumming K0s hypothesis for V0s rejected by misidentification cut
  TH1D *fLambdaMassOfMisIDV0;          // Mass assuming Lambda hypothesis for V0s rejected by misidentification cut
  TH1D *fAntiLambdaMassOfMisIDV0;      // Mass assuming AntiLambda hypothesis for V0s rejected by misidentification cut

  bool fIgnoreOnFlyStatus;  //This will accept V0s with aV0->OnFlyStatusV0()==true and aV0->OnFlyStatusV0()==false
                            //NOTE IMPORTANT: If you set this to true, be sure to call AliFemtoSimpleAnalysis::SetV0SharedDaughterCut(true)
                            //otherwise, in many cases, you will receive multiple copies of the same V0.

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCut, 1);
  /// \endcond
#endif

};

inline TH1D* AliFemtoV0TrackCut::GetMinvPurityAidHistoV0() {return fMinvPurityAidHistoV0;}
inline bool AliFemtoV0TrackCut::GetBuildMisIDHistograms() {return fBuildMisIDHistograms;}
inline void AliFemtoV0TrackCut::SetIgnoreOnFlyStatus(bool aIgnore) {fIgnoreOnFlyStatus = aIgnore;}

#endif
