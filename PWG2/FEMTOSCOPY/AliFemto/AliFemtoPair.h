///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoPair: the Pair object is passed to the PairCuts for           //
// verification, and then to the AddRealPair and AddMixedPair methods of //
// the Correlation Functions. It holds pair-specific variables like      //
// relative momenta and has links to the particles and tracks that form  //
// the pair.                                                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOPAIR_H
#define ALIFEMTOPAIR_H

#include <utility>

#include "AliFemtoParticle.h"
#include "AliFemtoTypes.h"

class AliFemtoPair {
public:
  AliFemtoPair();
  AliFemtoPair(const AliFemtoPair& aPair);
  AliFemtoPair(AliFemtoParticle*, AliFemtoParticle*);
  ~AliFemtoPair();
  AliFemtoPair& operator=(const AliFemtoPair& aPair);

  // track Gets:
  AliFemtoParticle* Track1() const;
  AliFemtoParticle* Track2() const;
  // track Sets:
  void SetTrack1(const AliFemtoParticle* trkPtr);
  void SetTrack2(const AliFemtoParticle* trkPtr);

  AliFemtoLorentzVector FourMomentumDiff() const;
  AliFemtoLorentzVector FourMomentumSum() const;
  double QInv() const;
  double KT()   const;
  double MInv() const;
  // pair rapidity
  double Rap() const;
  double EmissionAngle() const;

  // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
  double QSidePf() const;
  double QOutPf() const;
  double QLongPf() const;
   
  // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
  // - written by Bekele/Humanic
  double QSideCMS() const;
  double QOutCMS() const;
  double QLongCMS() const;

  double KSide() const;
  double KOut() const;
  double KLong() const;

  // Bertsch-Pratt momentum components in a longitudinally boosted frame
  // the argument is the beta of the longitudinal boost (default is 0.0, meaning lab frame)
  // - written by Bekele/Humanic
  double QSideBf(double beta=0.0) const;
  double QOutBf(double beta=0.0) const;
  double QLongBf(double beta=0.0) const;

  // Yano-Koonin-Podgoretskii Parametrisation 
  // source rest frame (usually lab frame)
  void QYKPCMS(double& qP, double& qT, double& q0) const ;
  // longitudinal comoving frame
  void QYKPLCMS(double& qP, double& qT, double& q0) const ;
  // pair rest frame
  void QYKPPF(double& qP, double& qT, double& q0) const ;


  double Quality() const;

  // the following two methods calculate the "nominal" separation of the tracks 
  // at the inner field cage (EntranceSeparation) and when they exit the TPC,
  // which may be at the outer field cage, or at the endcaps.
  // "nominal" means that the tracks are assumed to start at (0,0,0).  Making this
  // assumption is important for the Event Mixing-- it is not a mistake. - MALisa
  double NominalTpcExitSeparation() const;
  double NominalTpcEntranceSeparation() const;
  //  double NominalTpcAverageSeparation() const;
  // adapted calculation of Entrance/Exit/Average Tpc separation to V0 daughters
/*   double TpcExitSeparationTrackV0Pos() const; */
/*   double TpcEntranceSeparationTrackV0Pos() const; */
/*   double TpcAverageSeparationTrackV0Pos() const;  */

/*   double TpcExitSeparationTrackV0Neg() const; */
/*   double TpcEntranceSeparationTrackV0Neg() const; */
/*   double TpcAverageSeparationTrackV0Neg() const;  */

/*   double TpcExitSeparationV0PosV0Pos() const; */
/*   double TpcEntranceSeparationV0PosV0Pos() const; */
/*   double TpcAverageSeparationV0PosV0Pos() const;  */

/*   double TpcExitSeparationV0PosV0Neg() const; */
/*   double TpcEntranceSeparationV0PosV0Neg() const; */
/*   double TpcAverageSeparationV0PosV0Neg() const;  */
 
/*   double TpcExitSeparationV0NegV0Pos() const; */
/*   double TpcEntranceSeparationV0NegV0Pos() const; */
/*   double TpcAverageSeparationV0NegV0Pos() const;  */
  
/*   double TpcExitSeparationV0NegV0Neg() const; */
/*   double TpcEntranceSeparationV0NegV0Neg() const; */
/*   double TpcAverageSeparationV0NegV0Neg() const;  */

  double PInv() const;
  double KStar() const;
  double KStarFlipped() const;
  double CVK() const;
  double CVKFlipped() const;
  double QInvFlippedXY() const;

  double OpeningAngle() const;

  // Fabrice Private <<<
  double KStarSide() const;
  double KStarOut() const;
  double KStarLong() const;

  float PionPairProbability() const;
  float ElectronPairProbability() const;
  float KaonPairProbability() const;
  float ProtonPairProbability() const;
  float KaonPionPairProbability() const;

  double DcaInsideTpc() const;
  double Quality2() const;

 /* double KStarGlobal() const;
  double CVKGlobal() const;
  double KStarSideGlobal() const;
  double KStarOutGlobal() const;
  double KStarLongGlobal() const;*/

  void SetMergingPar(double aMaxDuInner, double aMaxDzInner,
		     double aMaxDuOuter, double aMaxDzOuter);
  void SetDefaultHalfFieldMergingPar();
  void SetDefaultFullFieldMergingPar();
  double GetFracOfMergedRow() const;
  double GetClosestRowAtDCA() const;
  double GetWeightedAvSep() const;
  // >>>
/*   double GetFracOfMergedRowTrkV0Pos() const; */
/*   double GetClosestRowAtDCATrkV0Pos() const; */

/*   double GetFracOfMergedRowTrkV0Neg() const; */
/*   double GetClosestRowAtDCATrkV0Neg() const; */

/*   double GetFracOfMergedRowV0PosV0Neg() const; */
/*   double GetFracOfMergedRowV0NegV0Pos() const; */
/*   double GetFracOfMergedRowV0PosV0Pos() const; */
/*   double GetFracOfMergedRowV0NegV0Neg() const; */

private:
  AliFemtoParticle* fTrack1; // Link to the first track in the pair
  AliFemtoParticle* fTrack2; // Link to the second track in the pair

  mutable short fNonIdParNotCalculated; // Set to 1 when NonId variables (kstar) have been already calculated for this pair
  mutable double fDKSide; // momemntum of first particle in PRF - k* side component
  mutable double fDKOut;  // momemntum of first particle in PRF - k* out component
  mutable double fDKLong; // momemntum of first particle in PRF - k* long component
  mutable double fCVK;    // cos between velocity and relative momentum k*
  mutable double fKStarCalc; // momemntum of first particle in PRF - k* 
  void CalcNonIdPar() const;

  mutable short fNonIdParNotCalculatedGlobal; // If global k* was calculated
 /* mutable double fDKSideGlobal;
  mutable double fDKOutGlobal;
  mutable double fDKLongGlobal;
  mutable double kStarCalcGlobal;
  mutable double fCVKGlobal;*/
  //void calcNonIdParGlobal() const;

  mutable short fMergingParNotCalculated; // If merging parameters were calculated
  mutable double fWeightedAvSep;          // Weighted average separation
  mutable double fFracOfMergedRow;        // Fraction of merged rows
  mutable double fClosestRowAtDCA;        // Row at wchich DCA occurs

  mutable short fMergingParNotCalculatedTrkV0Pos; // merging parameters for track - V0 pos	 
  mutable double fFracOfMergedRowTrkV0Pos;        // fraction of merged rows for track - V0 pos 
  mutable double fClosestRowAtDCATrkV0Pos;        // Row at which DCA occurs for track - V0 pos 

  mutable short fMergingParNotCalculatedTrkV0Neg; // merging parameters for track - V0 neg	 
  mutable double fFracOfMergedRowTrkV0Neg;	  // fraction of merged rows for track - V0 neg 
  mutable double fClosestRowAtDCATrkV0Neg;	  // Row at which DCA occurs for track - V0 neg 

  mutable short fMergingParNotCalculatedV0PosV0Neg; // merging parameters for V0 pos - V0 neg	 
  mutable double fFracOfMergedRowV0PosV0Neg;	    // fraction of merged rows for V0 pos - V0 neg 
  mutable double fClosestRowAtDCAV0PosV0Neg;	    // Row at which DCA occurs for V0 pos - V0 neg 

  mutable short fMergingParNotCalculatedV0NegV0Pos; // merging parameters for V0 neg - V0 pos	 
  mutable double fFracOfMergedRowV0NegV0Pos;	    // fraction of merged rows for V0 neg - V0 pos 
  mutable double fClosestRowAtDCAV0NegV0Pos;	    // Row at which DCA occurs for V0 neg - V0 pos 

  mutable short fMergingParNotCalculatedV0PosV0Pos; // merging parameters for V0 pos - V0 pos	 
  mutable double fFracOfMergedRowV0PosV0Pos;	    // fraction of merged rows for V0 pos - V0 pos 
  mutable double fClosestRowAtDCAV0PosV0Pos;	    // Row at which DCA occurs for V0 pos - V0 pos 

  mutable short fMergingParNotCalculatedV0NegV0Neg; // merging parameters for V0 neg - V0 neg	 
  mutable double fFracOfMergedRowV0NegV0Neg;	    // fraction of merged rows for V0 neg - V0 neg 
  mutable double fClosestRowAtDCAV0NegV0Neg;	    // Row at which DCA occurs for V0 neg - V0 neg 

  static double fgMaxDuInner; // Minimum cluster separation in x in inner TPC padrow
  static double fgMaxDzInner; // Minimum cluster separation in z in inner TPC padrow
  static double fgMaxDuOuter; // Minimum cluster separation in x in outer TPC padrow
  static double fgMaxDzOuter; // Minimum cluster separation in z in outer TPC padrow
  void CalcMergingPar() const;

  void CalcMergingParFctn(short* tmpMergingParNotCalculatedFctn,
			  float* tmpZ1,float* tmpU1,
			  float* tmpZ2,float* tmpU2,
			  int *tmpSect1,int *tmpSect2,
			  double* tmpFracOfMergedRow,
			  double* tmpClosestRowAtDCA
			  ) const;

  void ResetParCalculated();
};

inline void AliFemtoPair::ResetParCalculated(){
  fNonIdParNotCalculated=1;
  fNonIdParNotCalculatedGlobal=1;
  fMergingParNotCalculated=1;
  fMergingParNotCalculatedTrkV0Pos=1;
  fMergingParNotCalculatedTrkV0Neg=1;
  fMergingParNotCalculatedV0PosV0Pos=1;
  fMergingParNotCalculatedV0NegV0Pos=1;
  fMergingParNotCalculatedV0PosV0Neg=1;
  fMergingParNotCalculatedV0NegV0Neg=1;
}

inline void AliFemtoPair::SetTrack1(const AliFemtoParticle* trkPtr){
  fTrack1=(AliFemtoParticle*)trkPtr;
  ResetParCalculated();
}
inline void AliFemtoPair::SetTrack2(const AliFemtoParticle* trkPtr){
  fTrack2=(AliFemtoParticle*)trkPtr;
  ResetParCalculated();
}

inline AliFemtoParticle* AliFemtoPair::Track1() const {return fTrack1;}
inline AliFemtoParticle* AliFemtoPair::Track2() const {return fTrack2;}

inline double AliFemtoPair::KSide() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKSide;
}
inline double AliFemtoPair::KOut() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKOut;
}
inline double AliFemtoPair::KLong() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKLong;
}
inline double AliFemtoPair::KStar() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fKStarCalc;
}
inline double AliFemtoPair::QInv() const {
  AliFemtoLorentzVector tDiff = (fTrack1->FourMomentum()-fTrack2->FourMomentum());
  return ( -1.* tDiff.m());
}

// Fabrice private <<<
inline double AliFemtoPair::KStarSide() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKSide;//mKStarSide;
}
inline double AliFemtoPair::KStarOut() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKOut;//mKStarOut;
}
inline double AliFemtoPair::KStarLong() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fDKLong;//mKStarLong;
}
inline double AliFemtoPair::CVK() const{
  if(fNonIdParNotCalculated) CalcNonIdPar();
  return fCVK;
}

inline float AliFemtoPair::PionPairProbability() const{
  return (fTrack1->Track()->PidProbPion()) * 
         (fTrack2->Track()->PidProbPion());
}
inline float AliFemtoPair::ElectronPairProbability() const{
  return (fTrack1->Track()->PidProbElectron()) * 
         (fTrack2->Track()->PidProbElectron());
}
inline float AliFemtoPair::KaonPairProbability() const{
  return (fTrack1->Track()->PidProbKaon()) * 
         (fTrack2->Track()->PidProbKaon());
}
inline float AliFemtoPair::ProtonPairProbability() const{
  return (fTrack1->Track()->PidProbProton()) * 
         (fTrack2->Track()->PidProbProton());
}
inline float AliFemtoPair::KaonPionPairProbability() const{
  return (fTrack1->Track()->PidProbKaon()) * 
         (fTrack2->Track()->PidProbPion());
}

inline double AliFemtoPair::GetFracOfMergedRow() const{
  if(fMergingParNotCalculated) CalcMergingPar();
  return fFracOfMergedRow;
}
inline double AliFemtoPair::GetClosestRowAtDCA() const { 
  if(fMergingParNotCalculated) CalcMergingPar();
  return fClosestRowAtDCA;
}
inline double AliFemtoPair::GetWeightedAvSep() const {
  if(fMergingParNotCalculated) CalcMergingPar();
  return fWeightedAvSep;
}


/* inline double AliFemtoPair::GetFracOfMergedRowTrkV0Pos() const{ */
/*   if(fMergingParNotCalculatedTrkV0Pos) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Pos, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]), */
/* 		       &(fTrack1->fSect[0]),&(fTrack2->fSect[0]), */
/* 		       &(fFracOfMergedRowTrkV0Pos),&(fClosestRowAtDCATrkV0Pos) */
/* 		       ); */
/*   return fFracOfMergedRowTrkV0Pos; */
/* } */
/* inline double AliFemtoPair::GetClosestRowAtDCATrkV0Pos() const{ */
/*   if(fMergingParNotCalculatedTrkV0Pos) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Pos, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]), */
/* 		       &(fTrack1->fSect[0]),&(fTrack2->fSect[0]), */
/* 		       &fFracOfMergedRowTrkV0Pos,&fClosestRowAtDCATrkV0Pos */
/* 		       ); */
/*   return fClosestRowAtDCATrkV0Pos; */
/* } */
/* inline double AliFemtoPair::GetFracOfMergedRowTrkV0Neg() const{ */
/*   if(fMergingParNotCalculatedTrkV0Neg) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Neg, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]), */
/* 		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]), */
/* 		       &(fFracOfMergedRowTrkV0Neg),&(fClosestRowAtDCATrkV0Neg) */
/* 		       ); */
/*   return fFracOfMergedRowTrkV0Neg; */
/* } */
/* inline double AliFemtoPair::GetClosestRowAtDCATrkV0Neg() const{ */
/*   if(fMergingParNotCalculatedTrkV0Neg) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Neg, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]), */
/* 		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]), */
/* 		       &fFracOfMergedRowTrkV0Neg,&fClosestRowAtDCATrkV0Neg */
/* 		       ); */
/*   return fClosestRowAtDCATrkV0Neg; */
/* } */
/* inline double AliFemtoPair::GetFracOfMergedRowV0PosV0Neg() const{ */
/*   if(fMergingParNotCalculatedV0PosV0Neg) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedV0PosV0Neg, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]), */
/* 		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]), */
/* 		       &(fFracOfMergedRowV0PosV0Neg), */
/* 		       &(fClosestRowAtDCAV0PosV0Neg) */
/* 		       ); */
/*   return fFracOfMergedRowV0PosV0Neg; */
/* } */
/* inline double AliFemtoPair::GetFracOfMergedRowV0NegV0Pos() const{ */
/*   if(fMergingParNotCalculatedV0NegV0Pos) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedV0NegV0Pos, */
/* 		       &(fTrack1->fV0NegZ[0]),&(fTrack1->fV0NegU[0]), */
/* 		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]), */
/* 		       &(fTrack1->fV0NegSect[0]), */
/* 		       &(fTrack2->fSect[0]), */
/* 		       &(fFracOfMergedRowV0NegV0Pos), */
/* 		       &(fClosestRowAtDCAV0NegV0Pos) */
/* 		       ); */
/*   return fFracOfMergedRowV0NegV0Pos; */
/* } */
/* inline double AliFemtoPair::GetFracOfMergedRowV0PosV0Pos() const{ */
/*   if(fMergingParNotCalculatedV0PosV0Pos) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedV0PosV0Pos, */
/* 		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]), */
/* 		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]), */
/* 		       &(fTrack1->fSect[0]), */
/* 		       &(fTrack2->fSect[0]), */
/* 		       &(fFracOfMergedRowV0PosV0Pos), */
/* 		       &(fClosestRowAtDCAV0PosV0Pos) */
/* 		       ); */
/*   return fFracOfMergedRowV0PosV0Pos; */
/* } */
/* inline double AliFemtoPair::GetFracOfMergedRowV0NegV0Neg() const{ */
/*   if(fMergingParNotCalculatedV0NegV0Neg) */
/*     CalcMergingParFctn(&fMergingParNotCalculatedV0NegV0Neg, */
/* 		       &(fTrack1->fV0NegZ[0]),&(fTrack1->fV0NegU[0]), */
/* 		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]), */
/* 		       &(fTrack1->fV0NegSect[0]), */
/* 		       &(fTrack2->fV0NegSect[0]), */
/* 		       &(fFracOfMergedRowV0NegV0Neg), */
/* 		       &(fClosestRowAtDCAV0NegV0Neg) */
/* 		       ); */
/*   return fFracOfMergedRowV0NegV0Neg; */
/* } */

#endif
