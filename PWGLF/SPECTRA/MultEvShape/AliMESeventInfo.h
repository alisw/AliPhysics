#ifndef ALIMESEVENTINFO_H
#define ALIMESEVENTINFO_H

////////////////////////////////////////////////////////////////////////////
//  Event summary data for the Multiplicity and Event Shape group         //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//    Andrei Herghelegiu <aherghe@niham.nipne.ro>                         //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#ifndef ROOT_TNamed
#include <TNamed.h>
#endif
#define FW_MAX_ORDER 8

class AliMESeventInfo : public TObject
{
public:
  class AliMESevShape : public TNamed {
  friend class AliMESeventInfo;
  public:
    AliMESevShape();
    AliMESevShape(const char* name, const char* title);
    ~AliMESevShape() {};

    Double_t  GetDirectivity(Bool_t dp) const { return dp?fDir[0]:fDir[1]; }
    Double_t  GetSphericity() const           { return fSphericity; }
    Double_t  GetThrust(Double_t &phi) const  { phi=fThrust[1]; return fThrust[0]; }
    Double_t  GetRecoil() const               { return fRecoil; }
    const Double_t* GetFW() const             { return fFW; }
    Double_t  GetFW(Int_t mom) const          { return mom>=0&&mom<FW_MAX_ORDER?fFW[mom]:0.; }
    Double_t  GetMomLeading(Bool_t px=kTRUE) const { return px?fPxyLead[0]:fPxyLead[1]; }
  protected:
    AliMESevShape(Double_t dir[2], Double_t sphr, Double_t tr[2], Double_t rec, Double_t fw[FW_MAX_ORDER], Double_t pl[2]);
    AliMESevShape(const AliMESevShape &evs);

    Double_t  fDir[2];               // Directivity [0]=D+ [1]=D-
    Double_t  fSphericity;           // Sfericity
    Double_t  fThrust[2];            // Thrust and Trust axis
    Double_t  fRecoil;               // Recoil
    Double_t  fFW[FW_MAX_ORDER];     // Fox-Wolfram moments
    Double_t  fPxyLead[2];           // px & py of the leading particle
  private:
    AliMESevShape &operator=(const AliMESevShape &evs);
    ClassDef(AliMESevShape, 1)       // Event shape descriptor for MES
  };
//   enum EMESevStat{
//      kPileUp=BIT(14)   // PileUp flag
//     ,kAny=BIT(15)      // any other flag TODO
//   };
  enum EMESmultId{
     kGlob08=0    // Global multiplicity for eta in (-0.8,0.8)
    ,kComb        // Combined multiplicity
    ,kNoClSPD     // SPD number of clusters
    ,kV0M           // V0M percentile
    ,kComb0408  // Combined multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
    ,kNmult       // number of multiplicity estimators
  };
  enum EMESevQuality{
     kMBtrigger=0    // MB trigger flag
    ,kHMtrigger      // High Multiplicity trigger flag
    ,kPileUp         // PileUp flag
    ,kVertex         // vertex reconstructed
    ,kVertexType     // [y] global, [n] SPD
  };
  AliMESeventInfo();
  ~AliMESeventInfo() {};

  void            Clear(Option_t *opt);
  const AliMESevShape*  GetEventShape() const      { return &fEvShape; }
  Double_t        GetMultiplicity(EMESmultId typ) const  { return typ>=0&&typ<kNmult?fMultiplicity[typ]:-1;}
  Int_t           GetQuality() const               { return fQuality; }
  Double_t        GetVertexZ() const               { return fVertexZ;}

  Bool_t          IsPileUp() const                 { return TESTBIT(fQuality, kPileUp); }
  Bool_t          HasTriggerMB() const             { return TESTBIT(fQuality, kMBtrigger);}
  Bool_t          HasTriggerHM() const             { return TESTBIT(fQuality, kHMtrigger);}
  Bool_t          HasVertex() const                { return TESTBIT(fQuality, kVertex);}
  Bool_t          HasVertexGlobal() const          { return TESTBIT(fQuality, kVertexType);}
  Bool_t          HasVertexITS() const             { return !HasVertexGlobal();}

  void            Print(Option_t *o = "") const;         // *MENU*
  // Bool_t          MakeShape(TObjArray* tracks);
  void            MakeDirectivity(TObjArray* tracks);
  Bool_t          MakeThrust(TObjArray* tracks);
  void            MakeSphericity(TObjArray* tracks);
  Bool_t          MakeRecoil(TObjArray* tracks);
  Bool_t          MakeFoxWolframMoments(TObjArray* tracks);
  Bool_t          FindLeadingParticle(TObjArray* tracks);

  static Double_t Directivity(TObjArray* tracks, Bool_t etaSign);
  static Bool_t   LeadingParticleDirection(TObjArray* tracks, Double_t pxy[2]);
  static Bool_t   Thrust(TObjArray* tracks, Double_t t[2]);
  static Double_t Sphericity(TObjArray* tracks);
  static Double_t Recoil(TObjArray* tracks);
  static Bool_t   TransverseFoxWolframMoments(TObjArray* tracks, Double_t fw[FW_MAX_ORDER]);
  //static Bool_t   IsCDFJetCandidate(const TObjArray *ac, const Double_t radius, const Double_t leadingPt, const Double_t secondPt

  void            SetEvShape(const AliMESevShape& ev);   // *GETTER=GetEvShape
  void            SetEvShape(Double_t dir[2], Double_t sfr, Double_t tr[2], Double_t rec, Double_t fw[FW_MAX_ORDER], Double_t leadP[2]);
  void            SetMultiplicity(Int_t typ, Double_t mult) { if(typ>=0&&typ<kNmult) fMultiplicity[typ] = mult; }
  void            SetPileUp(Bool_t set=kTRUE)      { set?SETBIT(fQuality, kPileUp):CLRBIT(fQuality, kPileUp); }
  void            SetQuality(Int_t q)              { fQuality = q; }
  void            SetTriggerMB()                   { SETBIT(fQuality, kMBtrigger);}
  void            SetTriggerHM()                   { SETBIT(fQuality, kHMtrigger);}
  void            SetVertex()                      { SETBIT(fQuality, kVertex);}
  void            SetVertexGlob()                  { SETBIT(fQuality, kVertexType);}
  void            SetVertexZ(Double_t z)           { fVertexZ = z; }
private:
  AliMESeventInfo(const AliMESeventInfo &event);
  AliMESeventInfo &operator=(const AliMESeventInfo &event);

  UChar_t        fQuality;              //
  Double_t       fMultiplicity[kNmult]; // multiplicity estimators
  Double_t       fVertexZ;              // z coordinate of vertex
  AliMESevShape  fEvShape;              // event shape descriptor
  ClassDef(AliMESeventInfo, 2)          // Event summary data for MultiEvShape
};

#endif
