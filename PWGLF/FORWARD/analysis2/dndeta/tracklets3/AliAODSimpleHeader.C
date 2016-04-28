/**
 * @file   AliAODSimpleHeader.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:49:43 2016
 * 
 * @brief  A simplified AOD header
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
#ifndef ALIAODSIMPLEHEADER_H
#define ALIAODSIMPLEHEADER_H
#include <TVector3.h>

/**
 * A simplified header
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliAODSimpleHeader : public TObject
{
public:
  /** 
   * Default CTOR 
   */
  AliAODSimpleHeader()
    : TObject(),
      fTriggers        (0),
      fCent            (-1),
      fNTracklets      (-1),
      fCentOld         (-1),      
      fRecIP           (-1000, -1000, -1000),
      fSimIP           (-1000, -1000, -1000),
      fProjectileNpart (-1),
      fTargetNpart     (-1),
      fNcoll           (-1),
      fImpactParameter (-1),
      fReactionPlane   (-1),
      fProjectileNsd   (-1),
      fTargetNsd       (-1),
      fNdd             (-1)
  {}
  virtual ~AliAODSimpleHeader() {}
  /** 
   * Clear this object 
   * 
   * @param option Not used
   */
  void Clear(Option_t* option="")
  {
    fTriggers        = 0;
    fCent            = -1;
    fNTracklets      = -1;
    fCentOld         = -1;    
    fRecIP.SetXYZ(-1000, -1000, -1000);
    fSimIP.SetXYZ(-1000, -1000, -1000);
    fProjectileNpart = -1;
    fTargetNpart     = -1;
    fNcoll           = -1;
    fImpactParameter = -1;
    fReactionPlane   = -1;
    fProjectileNsd   = -1;
    fTargetNsd       = -1;
    fNdd             = -1;
  }
  /** 
   * @{ 
   * @name Members - all public 
   */
  /** Triggers */
  UInt_t fTriggers;
  /** Centrality - new framework */
  Double_t fCent;
  /** Number of tracklets - new framework */
  Int_t fNTracklets;
  /** Centrality - old framework */
  Double_t fCentOld;
  /** Reconstructed vertex */ 
  TVector3 fRecIP;
  /** Simulated vertex */
  TVector3 fSimIP;
  /** Number of projectile participants */
  Int_t fProjectileNpart;
  /** Number of target participants */
  Int_t fTargetNpart;
  /** Number of binary collisions */
  Int_t fNcoll;
  /** The impact parameter */
  Double_t fImpactParameter;
  /** The reaction plane */
  Double_t fReactionPlane;
  /** Projectile single-diffractive collisions */
  Int_t fProjectileNsd;
  /** Target single-diffractive collisions */
  Int_t fTargetNsd;
  /** Double diffractive collisions */
  Int_t fNdd;
  /* @} */

  ClassDef(AliAODSimpleHeader,1); // Simple header for AODs 
};



#endif
//
// EOF
//
