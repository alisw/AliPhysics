///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoParticle: main class halding all the necessary information    //
// about a particle that is required during femtoscopic analysis         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOPARTICLE_H
#define ALIFEMTOPARTICLE_H

//#include "math.h"
#include "AliFemtoTypes.h"
#include "AliFemtoTrack.h"
#include "AliFemtoV0.h"
#include "AliFemtoKink.h"
#include "AliFemtoXi.h"
#include "AliFmPhysicalHelixD.h"
// ***
class AliFemtoHiddenInfo;
// ***
class AliFemtoParticle{
public:
  AliFemtoParticle();
  AliFemtoParticle(const AliFemtoParticle& aParticle);
  AliFemtoParticle(const AliFemtoTrack* const hbtTrack, const double& mass);
  AliFemtoParticle(const AliFemtoV0* const hbtV0, const double& mass);
  AliFemtoParticle(const AliFemtoKink* const hbtKink, const double& mass);
  AliFemtoParticle(const AliFemtoXi* const hbtXi, const double& mass);
  ~AliFemtoParticle();

  AliFemtoParticle& operator=(const AliFemtoParticle& aParticle);

  const AliFemtoLorentzVector& FourMomentum() const;

  AliFmPhysicalHelixD& Helix();

  const AliFemtoThreeVector DecayVertexPosition() const;
  unsigned long TopologyMap(const int word) const;
  int NumberOfHits() const;

  unsigned long TrackId() const;         // only for particles from tracks 
  unsigned short   NegTrackId() const;   // only for particles from v0 
  unsigned short   PosTrackId() const;   // only for particles from v0 

  AliFemtoTrack* Track() const;
  AliFemtoV0* V0() const;
  AliFemtoKink* Kink() const;

  const AliFemtoThreeVector& NominalTpcExitPoint() const;     // position track exits TPC assuming start at (0,0,0)
  const AliFemtoThreeVector& NominalTpcEntrancePoint() const; // position track crosses IFC assuming start at (0,0,0)
  const AliFemtoThreeVector& TpcV0PosExitPoint() const;  
  const AliFemtoThreeVector& TpcV0PosEntrancePoint() const;
  const AliFemtoThreeVector& TpcV0NegExitPoint() const;  
  const AliFemtoThreeVector& TpcV0NegEntrancePoint() const;

  // the following method is for explicit internal calculation to fill datamembers.
  // It is invoked automatically if AliFemtoParticle constructed from AliFemtoTrack
  // void CalculateNominalTpcExitAndEntrancePoints(); 
  // NOTE - this requires the fHelix, so be sure this is filled


  AliFemtoThreeVector fNominalPosSample[11];  // I make this public for convenience and speed of AliFemtoPair()
  float fZ[45];  // Z position of cluster on padrow
  float fU[45];  // U position of cluster on padrow
  int fSect[45]; // Sector number of cluster on padrow

  void ResetFourMomentum(const AliFemtoLorentzVector& fourMomentum);

  const AliFemtoHiddenInfo*  HiddenInfo() const;
  // Fab private
  AliFemtoHiddenInfo*  GetHiddenInfo() const;
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  void CalculatePurity();
  double GetPionPurity();
  double GetKaonPurity();
  double GetProtonPurity();
  void CalculateTpcExitAndEntrancePoints( AliFmPhysicalHelixD* tHelix,
					  AliFemtoThreeVector* PrimVert,
					  AliFemtoThreeVector* SecVert,
					  AliFemtoThreeVector* tmpTpcEntrancePoint,
					  AliFemtoThreeVector* tmpTpcExitPoint,
					  AliFemtoThreeVector* tmpPosSample,
					  float* tmpZ,float* tmpU,int* tmpSect);

  // For V0 Neg Daugthers TpcEntrance/ExitPoints
  AliFemtoThreeVector* fTpcV0NegPosSample; // Sample of TPC V0 neg
  float* fV0NegZ;                          // Array of Neg Z cluster positions
  float* fV0NegU;                          // Array of Neg U cluster positions
  int* fV0NegSect;                         // Array of Neg cluster sectors
 
private:
  AliFemtoTrack* fTrack;  // copy of the track the particle was formed of, else Null
  AliFemtoV0* fV0;        // copy of the v0 the particle was formed of, else Null
  AliFemtoKink* fKink;    // copy of the v0 the particle was formed of, else Null
  AliFemtoXi* fXi;        // copy of the Xi the particle was formed of, else Null  

  AliFemtoLorentzVector fFourMomentum; // Particle momentum
  AliFmPhysicalHelixD fHelix;          // Particle trajectory helix
  //unsigned long  fMap[2]; 
  //int fNhits;
  AliFemtoThreeVector fNominalTpcExitPoint; // Point where particle exits TPC
  AliFemtoThreeVector fNominalTpcEntrancePoint; // Point where particle enters TPC
  AliFemtoHiddenInfo* fHiddenInfo;  // Fab private

  double fPurity[6];  // Purity variables

  static double fgPrimPimPar0; // purity parameterization parameter
  static double fgPrimPimPar1; // purity parameterization parameter
  static double fgPrimPimPar2; // purity parameterization parameter
  static double fgPrimPipPar0; // purity parameterization parameter 
  static double fgPrimPipPar1; // purity parameterization parameter
  static double fgPrimPipPar2; // purity parameterization parameter
  static double fgPrimPmPar0;  // purity parameterization parameter
  static double fgPrimPmPar1;  // purity parameterization parameter
  static double fgPrimPmPar2;  // purity parameterization parameter
  static double fgPrimPpPar0;  // purity parameterization parameter
  static double fgPrimPpPar1;  // purity parameterization parameter
  static double fgPrimPpPar2;  // purity parameterization parameter

   // For V0 Daugthers TpcEntrance/ExitPoints
  AliFemtoThreeVector fPrimaryVertex;   // primary vertex of V0
  AliFemtoThreeVector fSecondaryVertex; // secondary vertex of V0

  AliFmPhysicalHelixD fHelixV0Pos;            // helix for positive V0 daughter		   
  AliFemtoThreeVector fTpcV0PosEntrancePoint; // positive V0 daughter entrance point to TPC
  AliFemtoThreeVector fTpcV0PosExitPoint;     // positive V0 daughter exit point from TPC  

  AliFmPhysicalHelixD fHelixV0Neg;            // helix for negative V0 daughter		   
  AliFemtoThreeVector fTpcV0NegEntrancePoint; // negative V0 daughter entrance point to TPC
  AliFemtoThreeVector fTpcV0NegExitPoint;     // negative V0 daughter exit point from TPC  
};

inline AliFemtoTrack* AliFemtoParticle::Track() const { return fTrack; }
inline unsigned long  AliFemtoParticle::TrackId() const { return fTrack->TrackId(); }
inline const AliFemtoLorentzVector& AliFemtoParticle::FourMomentum() const {return fFourMomentum;}
inline AliFmPhysicalHelixD& AliFemtoParticle::Helix() {return fHelix;}
//inline unsigned long AliFemtoParticle::TopologyMap(const int word) const {return fMap[word];}
//inline int AliFemtoParticle::NumberOfHits() const {return fNhits;}
//by marek chojnacki to could compile 

inline unsigned long AliFemtoParticle::TopologyMap(const int word) const {return 1;}
inline int AliFemtoParticle::NumberOfHits() const {return 1;}

inline AliFemtoV0* AliFemtoParticle::V0() const { return fV0; }
inline unsigned short AliFemtoParticle::NegTrackId() const { return fV0->IdNeg(); }
inline unsigned short AliFemtoParticle::PosTrackId() const { return fV0->IdPos(); }
inline const AliFemtoThreeVector AliFemtoParticle::DecayVertexPosition() const {return fV0->DecayVertexV0(); }
// ***
inline AliFemtoHiddenInfo* AliFemtoParticle::GetHiddenInfo() const
{return fHiddenInfo;}
inline const AliFemtoHiddenInfo* AliFemtoParticle::HiddenInfo() const
{return fHiddenInfo;}
inline void AliFemtoParticle::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo)
{ fHiddenInfo = aHiddenInfo->Clone();}
// ***

inline void AliFemtoParticle::ResetFourMomentum(const AliFemtoLorentzVector& vec){fFourMomentum = vec;}

inline AliFemtoKink* AliFemtoParticle::Kink() const { return fKink; }

#endif
