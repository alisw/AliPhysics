///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoV0: special type of particle desling with the specifics       //
// of the V0 type of particle
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOV0_H
#define ALIFEMTOV0_H

#include "AliFemtoTypes.h" //same as in AliFemtoTrack.h
#include "AliFmPhysicalHelixD.h" // Gael 12 Sept 02
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#include "StStrangeMuDstMaker/StV0MuDst.h"
#endif
#endif
/* Th stuff */
#include "AliFemtoHiddenInfo.h"
/***/

class AliFemtoV0 {
public:
  AliFemtoV0();
  AliFemtoV0(const AliFemtoV0& v); // copy constructor
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
  AliFemtoV0( StV0MuDst&); // from strangeness V0 micro dst structure
#endif
#endif
  virtual ~AliFemtoV0(){if(fHiddenInfo) delete fHiddenInfo;}
  AliFemtoV0& operator=(const AliFemtoV0& aV0);


  float DecayLengthV0() const;       // 3-d decay distance						 
  AliFemtoThreeVector DecayVertexV0() const; // Coordinates of decay vertex				 
  AliFemtoThreeVector PrimaryVertex() const; // Coordinates of primary vertex				 
  float DecayVertexV0X() const; // Coordinates of decay vertex						 
  float DecayVertexV0Y() const; // Coordinates of decay vertex						 
  float DecayVertexV0Z() const; // Coordinates of decay vertex						 
  float DcaV0Daughters() const;      // DCA of v0 daughters at Decay vertex				 
  float DcaV0ToPrimVertex() const;   // DCA of v0 to primary vertex					 
  float DcaPosToPrimVertex() const;  // DCA of pos v0 daughter to pri vertex				 
  float DcaNegToPrimVertex() const;  // DCA of neg v0 daughter to pri vertex				 
  AliFemtoThreeVector MomPos() const;   // Momentum components of pos. daughter				 
  float MomPosX() const;   // Momentum components of pos. daughter					 
  float MomPosY() const;   // Momentum components of pos. daughter					 
  float MomPosZ() const;   // Momentum components of pos. daughter					 
  AliFemtoThreeVector MomNeg() const;   // Momentum components of neg. daughter				 
  float MomNegX() const;   // Momentum components of neg. daughter					 
  float MomNegY() const;   // Momentum components of neg. daughter					 
  float MomNegZ() const;   // Momentum components of neg. daughter					 
		        	                                                                         
  int   TpcHitsPos() const;          // Number of TPC hits on pos. daughter				 
  int   TpcHitsNeg() const;          // Number of TPC hits on neg. daughter				 
  unsigned long TrackTopologyMapPos(unsigned int word) const;						 
  unsigned long TrackTopologyMapNeg(unsigned int word) const;						 
		        	                                                                         
  AliFemtoThreeVector MomV0() const ;    // Momentum components of V0					 
  float MomV0X() const ;    // Momentum components of V0						 
  float MomV0Y() const ;    // Momentum components of V0						 
  float MomV0Z() const ;    // Momentum components of V0						 
  float AlphaV0() const ;             // Armenteros-Podolanski variable					 
  float PtArmV0() const ;            // Armenteros-Podolanski variable					 
  float ELambda() const ;             // Energy assuming lambda hypothesis				 
  float EK0Short() const ;            // Energy assuming k-short hypothesis				 
  float EPosProton() const ;          // Energy of pos. daughter assuming proton			 
  float EPosPion() const ;            // Energy of pos. daughter assuming pion				 
  float ENegProton() const ;          // Energy of neg. daughter assuming antiproton			 
  float ENegPion() const ;            // Energy of neg. daughter assuming pion				 
  float MassLambda() const ;          // Mass assuming lambda hypothesis				 
  float MassAntiLambda() const ;      // Mass assuming antilambda hypothesis				 
  float MassK0Short() const ;         // Mass assuming k-short hypothesis				 
  float RapLambda() const ;           // Rapidity assuming (anti) constlambda				 
  float RapK0Short() const ;          // Rapidity assuming k-short					 
  float CTauLambda() const ;          // Lifetime (ctau) const assuming (anti) constlambda		 
  float CTauK0Short() const ;         // Lifetime (ctau) const assuming k-short				 
  float PtV0() const ;                // Transverse momentum						 
  float PtotV0() const ;              // Total momentum							 
  float PtPos() const ;               // Transverse momentum of pos. daughter				 
  float PtotPos() const ;             // Total momentum of pos. daughter				 
  float DedxPos() const;              // dedx of Positive track						 
  float NumdedxPos() const;	                                                                         
// number of hits in dE/dX track of pos. daughter--Gael04Fev 2002					 
  float ErrdedxPos() const;              								 
// error on dedx of Positive track--Gael04Fev 2002							 
  float LendedxPos() const;              								 
// Length of dE/dX track of pos. daughter--Gael04Fev 2002						 
  float PseudoRapPos() const;              								 
// Length of dE/dX track of neg. daughter--Gael04Fev 2002						 
		        	                                                                         
  float PtNeg() const ;               // Transverse momentum of neg. daughter				 
  float PtotNeg() const ;             // Total momentum of neg. daughter				 
  float DedxNeg() const;              // dedx of Negative track						 
  float NumdedxNeg() const;	                                                                         
// number of hits in dE/dX track of neg. daughter--Gael04Fev 2002					 
  float ErrdedxNeg() const;              								 
// error on dedx of Negative track--Gael04Fev 2002							 
  float LendedxNeg() const;              								 
// Length of dE/dX track of neg. daughter--Gael04Fev 2002						 
  float PseudoRapNeg() const;              								 
// Length of dE/dX track of neg. daughter--Gael04Fev 2002						 
		        	                                                                         
  unsigned short   IdNeg() const;               // Id of negative track                			 
  unsigned short   IdPos() const;               // Id of positive track					 
  unsigned short   KeyNeg() const;               // Id of negative track				 
  unsigned short   KeyPos() const;               // Id of positive track				 
				                                                                         
  const AliFmPhysicalHelixD& HelixPos() const; // Gael 12 Sept 02					 
  const AliFmPhysicalHelixD& HelixNeg() const; // Gael 12 Sept 02                                        

  void UpdateV0(); // Fills derived info
  void SetdecayLengthV0(const float x);  
  void SetdecayVertexV0(const AliFemtoThreeVector v);  
  void SetdecayVertexV0X(const float x);
  void SetdecayVertexV0Y(const float x);
  void SetdecayVertexV0Z(const float x);
  void SetdcaV0Daughters(const float x); 
  void SetdcaV0ToPrimVertex(const float x);  
  void SetdcaPosToPrimVertex(const float x); 
  void SetdcaNegToPrimVertex(const float x); 
  void SetmomPos(const AliFemtoThreeVector v);  
  void SetmomPosX(const float x);  
  void SetmomPosY(const float x);  
  void SetmomPosZ(const float x);  
  void SetmomNeg(const AliFemtoThreeVector v);  
  void SetmomNegX(const float x);  
  void SetmomNegY(const float x);  
  void SetmomNegZ(const float x);  

  void SettpcHitsPos(const int& i);      
  void SettpcHitsNeg(const int& i);      

  void SetTrackTopologyMapPos(unsigned int word, const unsigned long& m);
  void SetTrackTopologyMapNeg(unsigned int word, const unsigned long& m);      

  void SetmomV0( AliFemtoThreeVector v);
  void SetmomV0X( float x);
  void SetmomV0Y( float x);
  void SetmomV0Z( float x);
  void SetalphaV0( float x);       
  void SetptArmV0( float x);       
  void SeteLambda( float x);     
  void SeteK0Short( float x);    
  void SetePosProton( float x);  
  void SetePosPion( float x);    
  void SeteNegProton( float x);  
  void SeteNegPion( float x);    
  void SetmassLambda( float x);  
  void SetmassAntiLambda( float x);
  void SetmassK0Short( float x);  
  void SetrapLambda( float x);    
  void SetrapK0Short( float x);   
  void SetcTauLambda( float x);   
  void SetcTauK0Short( float x);  
  void SetptV0( float x);         
  void SetptotV0( float x);       
  void SetptPos( float x);        
  void SetptotPos( float x);      
  void SetptNeg( float x);        
  void SetptotNeg( float x);
  void SetidNeg(const unsigned short& i);
  void SetidPos(const unsigned short& i);
  void SetdedxNeg(float x);
  void SeterrdedxNeg(float x);//Gael 04Fev2002
  void SetlendedxNeg(float x);//Gael 04Fev2002
  void SetpseudoRapNeg(float x);//Gael 04Fev2002
  void SetdedxPos(float x);
  void SeterrdedxPos(float x);//Gael 04Fev2002
  void SetlendedxPos(float x);//Gael 04Fev2002
  void SetpseudoRapPos(float x);//Gael 04Fev2002
  void SetkeyNeg(const unsigned short& i);
  void SetkeyPos(const unsigned short& i);
     
  void SetHelixPos(const AliFmPhysicalHelixD& h); // Gael 12 Sept 02
  void SetHelixNeg(const AliFmPhysicalHelixD& h); // Gael 12 Sept 02

  void SetprimaryVertex(const AliFemtoThreeVector v);//Gael 24 Sept 02
  /* Th stuff */
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  bool ValidHiddenInfo() const;
  // Fab private : (official : const AliFemtoHiddenInfo* HiddenInfo() const;
  AliFemtoHiddenInfo* GetHiddenInfo() const;
  /***/

protected:
 
  float fDecayLengthV0;                     // 3-d decay distance						 \\ V0 decay length                                                                  
  AliFemtoThreeVector fDecayVertexV0;	    // Coordinates of decay vertex				 
  AliFemtoThreeVector fPrimaryVertex;	    // Coordinates of primary vertex				 
  float fDcaV0Daughters;		    // DCA of v0 daughters at Decay vertex				 
  float fDcaV0ToPrimVertex;		    // DCA of v0 to primary vertex					 
  float fDcaPosToPrimVertex;		    // DCA of pos v0 daughter to pri vertex				 
  float fDcaNegToPrimVertex;		    // DCA of neg v0 daughter to pri vertex				 
  AliFemtoThreeVector fMomPos;		    // Momentum components of pos. daughter				 
  AliFemtoThreeVector fMomNeg;		    // Momentum components of neg. daughter				 

  unsigned long  fTrackTopologyMapPos[2];   // Topology map for positive daughter
  unsigned long  fTrackTopologyMapNeg[2];   // Topology map for negative daughter
       					    
  int   fTpcHitsPos;			    // Number of TPC hits for positive daughter
  int   fTpcHitsNeg;			    // Number of TPC hits for negative daughter
					    
  float fChi2V0;			    // Fit quality for V0
  float fClV0;				    // Confidence level for V0
  float fChi2Pos;			    // Fit quality for positive daughter
  float fClPos;				    // Confidence level for positive daughter
  float fChi2Neg;			    // Fit quality for negative daughter
  float fClNeg;				    // Confidence level for negative daughter
					    
  float fDedxPos;			    // dEdx positive daughter	    
  float fErrDedxPos;			    // dEdx error positive daughter 
  float fLenDedxPos;			    // dEdx length positive daughter
  					    
  float fDedxNeg;			    // dEdx negative daughter	    
  float fErrDedxNeg;			    // dEdx error negative daughter 
  float fLenDedxNeg;			    // dEdx length negative daughter
					    
  unsigned short fNufDedxPos;		    // Number of dEdx points positive
  unsigned short fNufDedxNeg;		    // Number of dEdx points negative
					    
  AliFmPhysicalHelixD fHelixPos;            // Helix for positive
  AliFmPhysicalHelixD fHelixNeg;            // Helix for negative
					    
  AliFemtoThreeVector fMomV0;		    // Momentum of the V0
  float fAlphaV0;			    // Armenteros-Podolanski variable					 
  float fPtArmV0;			    // Armenteros-Podolanski variable					 
  float fELambda;			    // Energy assuming lambda hypothesis				 
  float fEK0Short;			    // Energy assuming k-short hypothesis				 
  float fEPosProton;			    // Energy of pos. daughter assuming proton			 
  float fEPosPion;			    // Energy of pos. daughter assuming pion				 
  float fENegProton;			    // Energy of neg. daughter assuming antiproton
  float fENegPion;			    // Energy of neg. daughter assuming pion
  float fMassLambda;			    // Mass assuming lambda hypothesis				 
  float fMassAntiLambda;		    // Mass assuming antilambda hypothesis
  float fMassK0Short;			    // Mass assuming k-short hypothesis				 
  float fRapLambda;			    // Rapidity assuming (anti) constlambda				 
  float fRapK0Short;			    // Rapidity assuming k-short					 
  float fCTauLambda;			    // Lifetime (ctau) assuming (anti)lambda		 
  float fCTauK0Short;			    // Lifetime (ctau) assuming k-short				 
  float fPtV0;				    // Total momentum							 
  float fPtotV0;			    // Transverse momentum						 
  float fPtPos;				    // Transverse momentum of pos. daughter				 
  float fPtotPos;			    // Total momentum of pos. daughter				 
  float fPtNeg;				    // Transverse momentum of neg. daughter				 
  float fPtotNeg;			    // Total momentum of neg. daughter				 
					    
  unsigned short   fKeyNeg;		    // Unique key negative
  unsigned short   fKeyPos;		    // Unique key positive
  /* Th stuff */			    
  // Fab private : add mutable		    
  mutable AliFemtoHiddenInfo* fHiddenInfo; //! Hidden info
  /***/					    
					    
					        	                                                                         
};					    
					    
inline float AliFemtoV0::DecayLengthV0() const { return fDecayLengthV0; }
inline AliFemtoThreeVector AliFemtoV0::DecayVertexV0() const { return fDecayVertexV0; } 
inline AliFemtoThreeVector AliFemtoV0::PrimaryVertex() const { return fPrimaryVertex; }
inline float AliFemtoV0::DecayVertexV0X() const { return fDecayVertexV0.x(); } 
inline float AliFemtoV0::DecayVertexV0Y() const { return fDecayVertexV0.y(); } 
inline float AliFemtoV0::DecayVertexV0Z() const { return fDecayVertexV0.z(); } 
inline float AliFemtoV0::DcaV0Daughters() const { return fDcaV0Daughters; }
inline float AliFemtoV0::DcaV0ToPrimVertex() const { return fDcaV0ToPrimVertex; }
inline float AliFemtoV0::DcaPosToPrimVertex() const { return fDcaPosToPrimVertex; }
inline float AliFemtoV0::DcaNegToPrimVertex() const { return fDcaNegToPrimVertex; }
inline AliFemtoThreeVector AliFemtoV0::MomPos() const { return fMomPos; }
inline float AliFemtoV0::MomPosX() const { return fMomPos.x(); }
inline float AliFemtoV0::MomPosY() const { return fMomPos.y(); }
inline float AliFemtoV0::MomPosZ() const { return fMomPos.z(); }
inline AliFemtoThreeVector AliFemtoV0::MomNeg() const { return fMomNeg; }
inline float AliFemtoV0::MomNegX() const { return fMomNeg.x(); }
inline float AliFemtoV0::MomNegY() const { return fMomNeg.y(); }
inline float AliFemtoV0::MomNegZ() const { return fMomNeg.z(); }
inline AliFemtoThreeVector AliFemtoV0::MomV0() const { return fMomV0; }
inline float AliFemtoV0::MomV0X() const { return fMomV0.x(); }
inline float AliFemtoV0::MomV0Y() const { return fMomV0.y(); }
inline float AliFemtoV0::MomV0Z() const { return fMomV0.z(); }
inline float AliFemtoV0::AlphaV0() const { return fAlphaV0; }
inline float AliFemtoV0::PtArmV0() const {return fPtArmV0;}
inline float AliFemtoV0::ELambda() const {return fELambda;}
inline float AliFemtoV0::EK0Short() const {return fEK0Short;}
inline float AliFemtoV0::EPosProton() const {return fEPosProton;}
inline float AliFemtoV0::EPosPion() const {return fEPosPion;}
inline float AliFemtoV0::ENegProton() const {return fENegProton;}
inline float AliFemtoV0::ENegPion() const {return fENegPion;}
inline float AliFemtoV0::MassLambda() const {return fMassLambda;}
inline float AliFemtoV0::MassAntiLambda() const {return fMassAntiLambda;}
inline float AliFemtoV0::MassK0Short() const {return fMassK0Short;}
inline float AliFemtoV0::RapLambda() const {return fRapLambda;}
inline float AliFemtoV0::RapK0Short() const {return fRapK0Short;}
inline float AliFemtoV0::CTauLambda() const {return fCTauLambda;}
inline float AliFemtoV0::CTauK0Short() const {return fCTauK0Short;}
inline float AliFemtoV0::PtV0() const {return fPtV0;}
inline float AliFemtoV0::PtotV0() const {return fPtotV0;}
inline float AliFemtoV0::PtPos() const {return fPtPos;}
inline float AliFemtoV0::PtotPos() const {return fPtotPos;}
inline float AliFemtoV0::PtNeg() const {return fPtNeg;}
inline float AliFemtoV0::PtotNeg() const {return fPtotNeg;}
inline int   AliFemtoV0::TpcHitsPos() const { return fTpcHitsPos; }
inline int   AliFemtoV0::TpcHitsNeg() const { return fTpcHitsNeg; }
inline float AliFemtoV0::DedxNeg() const {return fDedxNeg;}
inline float AliFemtoV0::NumdedxNeg() const {return fNufDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::ErrdedxNeg() const {return fErrDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::LendedxNeg() const {return fLenDedxNeg;} //Gael 04Fev2002
inline float AliFemtoV0::PseudoRapNeg() const {return fMomNeg.pseudoRapidity();} //Gael 04Fev2002
inline float AliFemtoV0::DedxPos() const {return fDedxPos;}
inline float AliFemtoV0::NumdedxPos() const {return fNufDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::ErrdedxPos() const {return fErrDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::LendedxPos() const {return fLenDedxPos;} //Gael 04Fev2002
inline float AliFemtoV0::PseudoRapPos() const {return fMomPos.pseudoRapidity();} //Gael 04Fev2002


inline unsigned long   AliFemtoV0::TrackTopologyMapPos(unsigned int word) const { return fTrackTopologyMapPos[word]; }
inline unsigned long   AliFemtoV0::TrackTopologyMapNeg(unsigned int word) const { return fTrackTopologyMapNeg[word]; }
inline unsigned short   AliFemtoV0::IdNeg() const { return fKeyNeg; } 
inline unsigned short   AliFemtoV0::KeyNeg() const { return fKeyNeg; }
inline unsigned short   AliFemtoV0::IdPos() const { return fKeyPos; } 
inline unsigned short   AliFemtoV0::KeyPos() const { return fKeyPos; }

inline void AliFemtoV0::SetdecayLengthV0(const float x){ fDecayLengthV0= x;}   
inline void AliFemtoV0::SetdecayVertexV0X(const float x){ fDecayVertexV0.SetX(x);}
inline void AliFemtoV0::SetdecayVertexV0Y(const float x){ fDecayVertexV0.SetY(x);}
inline void AliFemtoV0::SetdecayVertexV0Z(const float x){ fDecayVertexV0.SetZ(x);}
inline void AliFemtoV0::SetdecayVertexV0(const AliFemtoThreeVector v){ fDecayVertexV0 = v; }
inline void AliFemtoV0::SetdcaV0Daughters(const float x){fDcaV0Daughters= x;} 
inline void AliFemtoV0::SetdcaV0ToPrimVertex(const float x){fDcaV0ToPrimVertex= x;}   
inline void AliFemtoV0::SetdcaPosToPrimVertex(const float x){fDcaPosToPrimVertex = x;} 
inline void AliFemtoV0::SetdcaNegToPrimVertex(const float x){fDcaNegToPrimVertex = x;} 
inline void AliFemtoV0::SetmomPos(const AliFemtoThreeVector v){fMomPos = v; }
inline void AliFemtoV0::SetmomPosX(const float x){fMomPos.SetX(x);}
inline void AliFemtoV0::SetmomPosY(const float x){fMomPos.SetY(x);}
inline void AliFemtoV0::SetmomPosZ(const float x){fMomPos.SetZ(x);}
inline void AliFemtoV0::SetmomNeg(const AliFemtoThreeVector v){fMomNeg = v; }
inline void AliFemtoV0::SetmomNegX(const float x){fMomNeg.SetX(x);}
inline void AliFemtoV0::SetmomNegY(const float x){fMomNeg.SetY(x);}
inline void AliFemtoV0::SetmomNegZ(const float x){fMomNeg.SetZ(x);}
inline void AliFemtoV0::SetTrackTopologyMapPos(unsigned int word, const unsigned long& m){fTrackTopologyMapPos[word]=m;} 
inline void AliFemtoV0::SetTrackTopologyMapNeg(unsigned int word, const unsigned long& m){fTrackTopologyMapNeg[word]=m;} 
inline void AliFemtoV0::SetmomV0(AliFemtoThreeVector v){fMomV0= v; }
inline void AliFemtoV0::SetmomV0X(const float x){fMomV0.SetX(x);}
inline void AliFemtoV0::SetmomV0Y(const float x){fMomV0.SetY(x);}
inline void AliFemtoV0::SetmomV0Z(const float x){fMomV0.SetZ(x);}

inline void AliFemtoV0::SetalphaV0( float x){fAlphaV0= x;}
inline void AliFemtoV0::SetptArmV0( float x){fPtArmV0 = x;}
inline void AliFemtoV0::SeteLambda( float x){fELambda= x;}       
inline void AliFemtoV0::SeteK0Short( float x){fEK0Short= x;}
inline void AliFemtoV0::SetePosProton( float x){fEPosProton= x;}      
inline void AliFemtoV0::SetePosPion( float x){fEPosPion= x;}      
inline void AliFemtoV0::SeteNegProton( float x){fENegProton= x;} 
inline void AliFemtoV0::SeteNegPion( float x){fENegPion= x;}       
inline void AliFemtoV0::SetmassLambda( float x){fMassLambda = x;} 
inline void AliFemtoV0::SetmassAntiLambda( float x){fMassAntiLambda= x;} 
inline void AliFemtoV0::SetmassK0Short( float x){fMassK0Short= x;}  
inline void AliFemtoV0::SetrapLambda( float x){fRapLambda= x;}
inline void AliFemtoV0::SetrapK0Short( float x){fRapK0Short = x;}   
inline void AliFemtoV0::SetcTauLambda( float x){fCTauLambda = x;}   
inline void AliFemtoV0::SetcTauK0Short( float x){fCTauK0Short = x;}   
inline void AliFemtoV0::SetptV0( float x){fPtV0 = x;}          
inline void AliFemtoV0::SetptotV0( float x){fPtotV0 = x;}
inline void AliFemtoV0::SetptPos( float x){fPtPos = x;}
inline void AliFemtoV0::SetptotPos( float x){fPtotPos = x;}    
inline void AliFemtoV0::SetptNeg( float x){ fPtNeg= x;}    
inline void AliFemtoV0::SetptotNeg( float x){ fPtotNeg= x;}
inline void AliFemtoV0::SetidNeg(const unsigned short& s){ fKeyNeg= s;}
inline void AliFemtoV0::SetidPos(const unsigned short& s){ fKeyPos= s;}
inline void AliFemtoV0::SetkeyNeg(const unsigned short& s){ fKeyNeg= s;}
inline void AliFemtoV0::SetkeyPos(const unsigned short& s){ fKeyPos= s;}
inline void AliFemtoV0::SettpcHitsPos(const int& i){fTpcHitsPos=i;} 
inline void AliFemtoV0::SettpcHitsNeg(const int& i){fTpcHitsNeg=i;}
inline void AliFemtoV0::SetdedxNeg(float x){fDedxNeg=x;}
inline void AliFemtoV0::SeterrdedxNeg(float x){fErrDedxNeg=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetlendedxNeg(float x){fLenDedxNeg=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetdedxPos(float x){fDedxPos=x;}
inline void AliFemtoV0::SeterrdedxPos(float x){fErrDedxPos=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetlendedxPos(float x){fLenDedxPos=x;}//Gael 04Fev2002
inline void AliFemtoV0::SetprimaryVertex(const AliFemtoThreeVector v) { fPrimaryVertex = v; }//Gael 24 Sept 02
#endif


















