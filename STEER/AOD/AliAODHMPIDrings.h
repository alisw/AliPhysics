#ifndef ALIAODHMPIDRINGS_H
#define ALIAODHMPIDRINGS_H


//
// Class to handle the AOD tracks with good HMPID rings 
// Author: Levente Molnar
// levente.molnar@cern.ch , March 2012
// 



//___ROOT includes
#include <TMath.h>
//___AliRoot includes
#include "AliPID.h"


class AliAODHMPIDrings : public TObject {
  
 public:
  
  AliAODHMPIDrings();
  AliAODHMPIDrings(
                    Int_t trkId,
                    Int_t qn, 
                    Int_t cluIdx,
                    Double_t  trkTheta,
                    Double_t trkPhi,
                    Double_t signal,
                    Double_t occ,
                    Double_t chi2,
                    Double_t trkX,
                    Double_t trkY,
                    Double_t mipX,
                    Double_t mipY,
                    Double_t *pid,
                    Double_t *p                  );      //              
      
      
  AliAODHMPIDrings(const AliAODHMPIDrings& hmpidAOD);//
  AliAODHMPIDrings &operator=(const AliAODHMPIDrings& hmpidAOD);//
  virtual ~AliAODHMPIDrings() {};
    
  //___ Getters
  Int_t GetHmpTrkID()                const { return fHmpidAODtrkId; }
  
  Double32_t GetHmpMipCharge()            const { return fHmpidAODqn%1000000; }
  Double32_t GetHmpNumOfPhotonClusters()  const { return fHmpidAODqn/1000000;}
 
  Int_t      GetHmpChamber()              const { return fHmpidAODcluIdx/1000000; }
  
  Int_t      GetHmpCluIdx()               const { return fHmpidAODcluIdx; }
  
  Double_t GetHmpTrackTheta()           const { return fHmpidAODtrkTheta;}
  Double_t GetHmpTrackPhi()             const { return fHmpidAODtrkPhi;}
  
  Double_t GetHmpSignal()               const { return fHmpidAODsignal;}
  Double_t GetHmpOccupancy()            const { return fHmpidAODocc;}
  
  Double_t GetHmpChi2()                 const { return fHmpidAODchi2;}

  Double_t GetHmpTrackX()               const { return fHmpidAODtrkX;}
  Double_t GetHmpTrackY()               const { return fHmpidAODtrkY;}

  Double_t GetHmpMipX()                 const { return fHmpidAODmipX;}
  Double_t GetHmpMipY()                 const { return fHmpidAODmipY;}

  Double_t GetHmpDX()                   const { return fHmpidAODmipX - fHmpidAODtrkX;}
  Double_t GetHmpDY()                   const { return fHmpidAODmipY - fHmpidAODtrkY;}
  Double_t GetHmpDist()                 const { return TMath::Sqrt((fHmpidAODmipX - fHmpidAODtrkX)*(fHmpidAODmipX - fHmpidAODtrkX) + (fHmpidAODmipY - fHmpidAODtrkY)*(fHmpidAODmipY - fHmpidAODtrkY));}
  
  
  void GetHmpPidProbs(Double32_t *pid) const;   //defined in cxx
  void GetHmpMom(Double32_t *mom)      const;   //defined in cxx
  
  //___ Setters
  
  void SetHmpMipCharge(Int_t q)               { fHmpidAODqn = q; }
  void SetHmpCluIdx(Int_t ch,Int_t idx)       { fHmpidAODcluIdx=ch*1000000+idx;}
  
  void SetHmpNumOfPhotonClusters(Int_t nph)   { fHmpidAODqn = 1000000 * nph;}
  
  void SetHmpTrackTheta(Double_t trkTheta)  { fHmpidAODtrkTheta = trkTheta;}
  void SetHmpTrackPhi(Double_t trkPhi)      { fHmpidAODtrkPhi = trkPhi;}
  
  void SetHmpSignal(Double_t thetaC)        { fHmpidAODsignal = thetaC;}
  void SetHmpOccupancy(Double_t occ)        { fHmpidAODocc =  occ;}
  
  void SetHmpChi2(Double_t chi2)            { fHmpidAODchi2 = chi2;}

  void SetHmpTrackX(Double_t trkX)          { fHmpidAODtrkX = trkX;}
  void SetHmpTrackY(Double_t trkY)          { fHmpidAODtrkY = trkY;}

  void SetHmpMipX(Double_t mipX)            { fHmpidAODmipX = mipX;}
  void SetHmpMipY(Double_t mipY)            { fHmpidAODmipY = mipY;}
 
  void SetHmpPidProbs(Double_t *pid);       
  void SetHmpMom(Double_t *mom);        
  
  
  
 protected:
  
  Int_t       fHmpidAODtrkId;                      // Unique track id as in ESD
  Int_t       fHmpidAODqn;                         // 1000000*number of photon clusters + QDC
  Int_t       fHmpidAODcluIdx;                     // 1000000*chamber id + cluster idx of the assigned MIP cluster
  
  Double32_t  fHmpidAODtrkTheta;                   // [-2*pi,2*pi,16] theta of the track extrapolated to the HMPID, LORS
  Double32_t  fHmpidAODtrkPhi;                     // [-2*pi,2*pi,16] theta of the track extrapolated to the HMPID, LORS
  Double32_t  fHmpidAODsignal;                     // HMPID signal (Theta ckov, rad)
  Double32_t  fHmpidAODocc;                        // chamber occupancy where the track passed through: number of pads
  Double32_t  fHmpidAODchi2;                       // [0.,0.,8] chi2 in the HMPID  
  Double32_t  fHmpidAODtrkX;                       // x of the track impact, LORS 
  Double32_t  fHmpidAODtrkY;                       // y of the track impact, LORS 
  Double32_t  fHmpidAODmipX;                       // x of the MIP in LORS
  Double32_t  fHmpidAODmipY;                       // y of the MIP in LORS
  Double32_t  fHmpidAODpid[AliPID::kSPECIES];      // [0.,0.,8] "detector response probabilities" (for the PID)
  Double32_t  fHMPIDmom[3];                        // track momentum at the HMPID ring reconstruction
  
  ClassDef(AliAODHMPIDrings,2)  
        
};
#endif
