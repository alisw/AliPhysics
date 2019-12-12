//---completely from scratch, this SHOULD BE REMOVED before pushing!

#ifndef ALIANALYSISTASKSTRANGECASCADESDISCRETE_H
#define ALIANALYSISTASKSTRANGECASCADESDISCRETE_H

class AliAnalysisManager;
class AliMultSelection;

//Class which saves the cascade daughter tracks
class AliRunningCascadeTrack : public TObject
{
private:
    Float_t PMom[3];
    Float_t NMom[3];
    Float_t BMom[3];
    Short_t dca_pos_to_prim[2];  // xy, z, range is 0..180 cm
    Short_t dca_neg_to_prim[2]; // xy, z, range is 0..180 cm
    Short_t dca_bach_to_prim[2]; // xy, z, range is 0..180 cm
    Short_t dca_V0_to_prim;  // 0..120
    Short_t dca_Omega_to_prim[2]; // xy, z; 0..40, pm40
    Short_t dca_pos_to_neg; //0..1.5 cm
    Short_t dca_bach_to_Lambda; //0...2.0 cm
    Short_t dca_bach_to_baryon; // 0..500.0 cm
    
    Short_t CosPointingAngle; //also stores the charge of the bachelor track! (cos(THeta)*e_bach)
    Float_t CascadeDecayPos[3];
    Float_t V0fromCascadePos[3];
    Short_t TrackLengthTPC[3]; // 0 - pos, 1 - neg, 2-bach

    
    //extra PID info
    //ITS
    UShort_t ITSstatusPosTrack[6], ITSstatusNegTrack[6], ITSstatusBachTrack[6];//2 byte, unsigned out:from 0 to 6.
    Char_t ITSrefitFlag[3];
    Char_t ITSPosSharedPoints[6], ITSNegSharedPoints[6], ITSBachSharedPoints[6];
    Short_t ITSchi2[3], ITSsignal[3];
    Short_t nSigma_ITS_pos[2];
    Short_t nSigma_ITS_neg[2];
    Short_t nSigma_ITS_bach[2];
  /*  Short_t ITSPIDProbPos[2]; //computed probability (AliPIDResponse; see AliPID.h for info)
    Short_t ITSPIDProbNeg[2];
    Short_t ITSPIDProbBach[2];*/
    
    //TPC
    Short_t TPCsignal[3]; //[0]: pos track, [1]: neg track, [2]: bach track
    UShort_t TPCcls[3]; //->GetTPCNcls() (not TPCNclsF())
    UShort_t TPCclsF[3]; //findable clusters! theoretically to expect
    UShort_t TPCchi2[3];
    Short_t nSigma_dEdx_pos[2]; //pos is either [0]=proton(if Omega-) or [1]=pion(if Omega+)
    Short_t nSigma_dEdx_neg[2]; //neg is either [0]=pion(Omega-) or [1]=proton
    Short_t nSigma_dEdx_bach[2]; //bach is either [0]=kaon(Omega) or [1]=pion(Xi)
   /* Short_t TPCPIDProbPos[2]; //computed probability (AliPIDResponse; see AliPID.h for info)
    Short_t TPCPIDProbNeg[2];
    Short_t TPCPIDProbBach[2];*/
    Char_t posTPC[7];
    Char_t negTPC[7];
    Char_t bachTPC[7];
    
    //TOF
    Short_t TOFsignal[3];
    Char_t TOFrefitFlag[3];
    Short_t nSigma_TOF_pos[2];
    Short_t nSigma_TOF_neg[2];
    Short_t nSigma_TOF_bach[2];
    
    
    
    
public:
    //default constructor: for my classed, inherited by TOBject, should be empty as possible
    AliRunningCascadeTrack():dca_V0_to_prim(-10),dca_pos_to_neg(-10),dca_bach_to_Lambda(-10),dca_bach_to_baryon(-10),CosPointingAngle(-10.){}
    //explicit contructor: should do at least something
    AliRunningCascadeTrack(Int_t iii):dca_V0_to_prim(-10),dca_pos_to_neg(-10),dca_bach_to_Lambda(-10),dca_bach_to_baryon(-10),CosPointingAngle(-10.)
    {
       // Float_t sum = dca_V0_to_prim + dca_pos_to_neg; //does not have any physic.interp.
        dca_V0_to_prim = -5;
        
    }
    //copy constructor
    AliRunningCascadeTrack(const AliRunningCascadeTrack&):TObject(),dca_V0_to_prim(-10),dca_pos_to_neg(-10),dca_bach_to_Lambda(-10),dca_bach_to_baryon(-10),CosPointingAngle(-10.){}
    //copy assignment
    AliRunningCascadeTrack& operator = (const AliRunningCascadeTrack&){return *this;}
    //default destructor
    virtual ~AliRunningCascadeTrack()
    {
        
    }
    
//----------------------------------------------------------------

    //setters
    void set_PMom(Double_t px, Double_t py, Double_t pz)   {PMom[0] = (Float_t)px; PMom[1] = (Float_t)py; PMom[2] = (Float_t)pz;}
    void set_NMom(Float_t px, Float_t py, Float_t pz)   {NMom[0] = (Float_t)px; NMom[1] = (Float_t)py; NMom[2] = (Float_t)pz;}
    void set_BMom(Float_t px, Float_t py, Float_t pz)    {BMom[0] = (Float_t)px; BMom[1] = (Float_t)py; BMom[2] = (Float_t)pz;}
    void set_dca_pos_to_prim(Float_t f1, Float_t f2)   { dca_pos_to_prim[0]   = (Short_t)(f1*100.0); dca_pos_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_neg_to_prim(Float_t f1, Float_t f2)   { dca_neg_to_prim[0]   = (Short_t)(f1*100.0); dca_neg_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_bach_to_prim(Float_t f1, Float_t f2)   { dca_bach_to_prim[0]   = (Short_t)(f1*100.0); dca_bach_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_V0_to_prim(Float_t f)  { dca_V0_to_prim  = (Short_t)(f*100.0); }
    void set_dca_Omega_to_prim(Float_t f1, Float_t f2)   { dca_Omega_to_prim[0]   = (Short_t)(f1*100.0); dca_Omega_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_pos_to_neg(Float_t f)  { dca_pos_to_neg = (Short_t)(f*100.0); }
    void set_dca_bach_to_Lambda(Float_t f)  { dca_bach_to_Lambda = (Short_t)(f*100.0); }
    void set_dca_bach_to_baryon(Float_t f)  { dca_bach_to_baryon = (Short_t)(f*10.0); }
    void set_CosPointingAngle(Float_t f)    {CosPointingAngle = (Short_t)(f*1000.);} //pay attention to this!
   // void set_LeastNumberOfTPCclus(Int_t i)      {LeastNumberOfTPCclus = i; }
   // void set_MaxChi2perclus(Float_t f)   {MaxChi2perclus = (Short_t)(f*100.); }
   // void set_MinTrackLength(Float_t f)    {MinTrackLength = (Short_t)(f*100.); }
    void set_CascadeDecayPos(Float_t f1, Float_t f2, Float_t f3)
    {
        CascadeDecayPos[0] = f1;
        CascadeDecayPos[1] = f2;
        CascadeDecayPos[2] = f3;
    }
    void set_V0fromCascadePos(Float_t f1, Float_t f2, Float_t f3)
    {
        V0fromCascadePos[0] = f1;
        V0fromCascadePos[1] = f2;
        V0fromCascadePos[2] = f3;
    }
    
    void set_TrackLengthTPC(Float_t l1, Float_t l2, Float_t l3){
        TrackLengthTPC[0] = (Short_t)l1; TrackLengthTPC[1] = (Short_t)l2; TrackLengthTPC[2] = (Short_t)l3;
    }
    
    
    //ITS setters
    void set_ITSstatusPosTrack(Int_t l0, Int_t l1,Int_t l2, Int_t l3,Int_t l4, Int_t l5){
        ITSstatusPosTrack[0] = (UShort_t)l0; ITSstatusPosTrack[1] = (UShort_t)l1; ITSstatusPosTrack[2] = (UShort_t)l2;
        ITSstatusPosTrack[3] = (UShort_t)l3; ITSstatusPosTrack[4] = (UShort_t)l4; ITSstatusPosTrack[5] = (UShort_t)l5;;
    }
    
    void set_ITSstatusNegTrack(Int_t l0, Int_t l1,Int_t l2, Int_t l3,Int_t l4, Int_t l5){
        ITSstatusNegTrack[0] = (UShort_t)l0; ITSstatusNegTrack[1] = (UShort_t)l1; ITSstatusNegTrack[2] = (UShort_t)l2;
        ITSstatusNegTrack[3] = (UShort_t)l3; ITSstatusNegTrack[4] = (UShort_t)l4; ITSstatusNegTrack[5] = (UShort_t)l5;;
    }
    void set_ITSstatusBachTrack(Int_t l0, Int_t l1,Int_t l2, Int_t l3,Int_t l4, Int_t l5){
        ITSstatusBachTrack[0] = (UShort_t)l0; ITSstatusBachTrack[1] = (UShort_t)l1; ITSstatusBachTrack[2] = (UShort_t)l2;
        ITSstatusBachTrack[3] = (UShort_t)l3; ITSstatusBachTrack[4] = (UShort_t)l4;ITSstatusBachTrack[5] = (UShort_t)l5;;
    }
    
    void set_ITSPosSharedPoints(Int_t s0, Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5){
        ITSPosSharedPoints[0] = (Char_t)s0,   ITSPosSharedPoints[1] = (Char_t)s1,   ITSPosSharedPoints[2] = (Char_t)s2,
        ITSPosSharedPoints[3] = (Char_t)s3,   ITSPosSharedPoints[4] = (Char_t)s4,   ITSPosSharedPoints[5] = (Char_t)s5;
    }
    
    void set_ITSNegSharedPoints(Int_t s0, Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5){
        ITSNegSharedPoints[0] = (Char_t)s0,   ITSNegSharedPoints[1] = (Char_t)s1,   ITSNegSharedPoints[2] = (Char_t)s2,
        ITSNegSharedPoints[3] = (Char_t)s3,   ITSNegSharedPoints[4] = (Char_t)s4,   ITSNegSharedPoints[5] = (Char_t)s5;
    }
    
    void set_ITSBachSharedPoints(Int_t s0, Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5){
        ITSBachSharedPoints[0] = (Char_t)s0,   ITSBachSharedPoints[1] = (Char_t)s1,   ITSBachSharedPoints[2] = (Char_t)s2,
        ITSBachSharedPoints[3] = (Char_t)s3,   ITSBachSharedPoints[4] = (Char_t)s4,   ITSBachSharedPoints[5] = (Char_t)s5;
    }
    
    void set_ITSrefitFlag(Int_t refitPos, Int_t refitNeg, Int_t refitBach){
        ITSrefitFlag[0] = (Char_t)refitPos; ITSrefitFlag[1]= (Char_t)refitNeg; ITSrefitFlag[2] = (Char_t)refitBach;
    }
    
    void set_ITSchi2(Double_t itschi2pos, Double_t itschi2neg, Double_t itschi2bach){
        ITSchi2[0] = (Short_t)(itschi2pos*100); ITSchi2[1] = (Short_t)(itschi2neg*100); ITSchi2[2] = (Short_t)(itschi2bach*100);
    }
    
    void set_ITSsignal(Double_t spos, Double_t sneg, Double_t sbach){
        ITSsignal[0] = (Short_t)(spos*10); ITSsignal[1] = (Short_t)(sneg*10); ITSsignal[2] = (Short_t)(sbach*10);
    }
    void set_nSigma_ITS_pos(Float_t f1, Float_t f2)   { nSigma_ITS_pos[0]   = (Short_t)(f1*100); nSigma_ITS_pos[1]   = (Short_t)(f2*100); }
    void set_nSigma_ITS_neg(Float_t f1, Float_t f2)   { nSigma_ITS_neg[0]   = (Short_t)(f1*100); nSigma_ITS_neg[1]   = (Short_t)(f2*100); }
    void set_nSigma_ITS_bach(Float_t f1, Float_t f2)   { nSigma_ITS_bach[0]   = (Short_t)(f1*100); nSigma_ITS_bach[1]   = (Short_t)(f2*100); }
   /* void set_ITSPIDProbPos(Double_t probProton, Double_t probPion){
        ITSPIDProbPos[0] = (Short_t)(probProton*10000);
        ITSPIDProbPos[1] = (Short_t)(probPion*10000);
    }
    void set_ITSPIDProbNeg(Double_t probProton, Double_t probPion){
        ITSPIDProbNeg[0] = (Short_t)(probProton*10000);
        ITSPIDProbNeg[1] = (Short_t)(probPion*10000);
    }
    void set_ITSPIDProbBach(Double_t probPion, Double_t probKaon){
        ITSPIDProbBach[0] = (Short_t)(probPion*10000);
        ITSPIDProbBach[1] = (Short_t)(probKaon*10000);
    }*/
    
    
    
    
    //TPC setters
    void set_TPCsignal(Double_t TPCsignalPos, Double_t TPCsignalNeg, Double_t TPCsignalBach){
        TPCsignal[0] = (Short_t)(TPCsignalPos*10.); TPCsignal[1] = (Short_t)(TPCsignalNeg*10.); TPCsignal[2] = (Short_t)(TPCsignalBach*10.);
    }
    void set_TPCcls(Int_t clpos, Int_t clneg, Int_t clbach){
        TPCcls[0] = (UShort_t)clpos;  TPCcls[1] = (UShort_t)clneg; TPCcls[2] = (UShort_t)clbach;
    }
    
    void set_TPCclsF(Int_t clposF, Int_t clnegF, Int_t clbachF){
        TPCclsF[0] = (UShort_t)clposF;  TPCclsF[1] = (UShort_t)clnegF; TPCclsF[2] = (UShort_t)clbachF;
    }
    
    void set_TPCchi2(Double_t chi2pos, Double_t chi2neg, Double_t chi2bach){
        TPCchi2[0] = (UShort_t)(chi2pos*100); TPCchi2[1] = (UShort_t)(chi2neg*100); TPCchi2[2] = (UShort_t)(chi2bach*100);
    }
    void set_nSigma_dEdx_pos(Float_t f1, Float_t f2)   { nSigma_dEdx_pos[0]   = (Short_t)(f1*100); nSigma_dEdx_pos[1]   = (Short_t)(f2*100); }
    void set_nSigma_dEdx_neg(Float_t f1, Float_t f2)   { nSigma_dEdx_neg[0]   = (Short_t)(f1*100); nSigma_dEdx_neg[1]   = (Short_t)(f2*100); }
    void set_nSigma_dEdx_bach(Float_t f1, Float_t f2)   { nSigma_dEdx_bach[0]   = (Short_t)(f1*100); nSigma_dEdx_bach[1]   = (Short_t)(f2*100); }
   /* void set_TPCPIDProbPos(Double_t probProton, Double_t probPion){
        TPCPIDProbPos[0] = (Short_t)(probProton*10000);
        TPCPIDProbPos[1] = (Short_t)(probPion*10000);
    }
    void set_TPCPIDProbNeg(Double_t probProton, Double_t probPion){
        TPCPIDProbNeg[0] = (Short_t)(probProton*10000);
        TPCPIDProbNeg[1] = (Short_t)(probPion*10000);
    }
    void set_TPCPIDProbBach(Double_t probPion, Double_t probKaon){
        TPCPIDProbBach[0] = (Short_t)(probPion*10000);
        TPCPIDProbBach[1] = (Short_t)(probKaon*10000);
    }*/
    void set_TPC_PartVar_Pos(Char_t p0, Char_t p1, Char_t p2, Char_t p3, Char_t p4, Char_t p5, Char_t p6)
    {
        posTPC[0] = p0; posTPC[1] = p1; posTPC[2] = p2; posTPC[3] = p3; posTPC[4] = p4; posTPC[5] = p5; posTPC[6] = p6;
    }
    void set_TPC_PartVar_Neg(Char_t n0, Char_t n1, Char_t n2, Char_t n3, Char_t n4, Char_t n5, Char_t n6)
    {
        negTPC[0] = n0; negTPC[1] = n1; negTPC[2] = n2; negTPC[3] = n3; negTPC[4] = n4; negTPC[5] = n5; negTPC[6] = n6;
    }
    void set_TPC_PartVar_Bach(Char_t b0, Char_t b1, Char_t b2, Char_t b3, Char_t b4, Char_t b5, Char_t b6)
    {
        bachTPC[0] = b0; bachTPC[1] = b1; bachTPC[2] = b2; bachTPC[3] = b3; bachTPC[4] = b4; bachTPC[5] = b5; bachTPC[6] = b6;
    }
    
    
    //TOF
    void set_TOFrefitFlag(Int_t refitPos, Int_t refitNeg, Int_t refitBach){
        TOFrefitFlag[0] = (Char_t)refitPos; TOFrefitFlag[1]= (Char_t)refitNeg; TOFrefitFlag[2] = (Char_t)refitBach;
    }
    void set_TOFsignal(Double_t TOFsignalPos, Double_t TOFsignalNeg, Double_t TOFsignalBach){
        TOFsignal[0] = (Float_t)(TOFsignalPos); TOFsignal[1] = (Float_t)(TOFsignalNeg); TOFsignal[2] = (Float_t)(TOFsignalBach);
    }
    void set_nSigma_TOF_pos(Float_t f1, Float_t f2)   { nSigma_TOF_pos[0]   = (Short_t)(f1*100); nSigma_TOF_pos[1]   = (Short_t)(f2*100); }
    void set_nSigma_TOF_neg(Float_t f1, Float_t f2)   { nSigma_TOF_neg[0]   = (Short_t)(f1*100); nSigma_TOF_neg[1]   = (Short_t)(f2*100); }
    void set_nSigma_TOF_bach(Float_t f1, Float_t f2)   { nSigma_TOF_bach[0]   = (Short_t)(f1*100); nSigma_TOF_bach[1]   = (Short_t)(f2*100); }
    
    
    
    
    //getters
    Float_t get_PMom(Int_t i) const     { return PMom[i]; }
    Float_t get_NMom(Int_t i) const     { return NMom[i]; }
    Float_t get_BMom(Int_t i) const     { return BMom[i]; }
    Float_t get_dca_pos_to_prim(Int_t i) const   { return ((Float_t)dca_pos_to_prim[i])/100.0; }
    Float_t get_dca_neg_to_prim(Int_t i) const   { return ((Float_t)dca_neg_to_prim[i])/100.0; }
    Float_t get_dca_bach_to_prim(Int_t i)  const { return ((Float_t)dca_bach_to_prim[i])/100.0; }
    Float_t get_dca_V0_to_prim() const  { return ((Float_t)dca_V0_to_prim)/100.0; }
    Float_t get_dca_Omega_to_prim(Int_t i) const  { return ((Float_t)dca_pos_to_prim[i])/100.0; }
    Float_t get_dca_pos_to_neg() const { return ((Float_t)dca_pos_to_neg)/100.0; }
    Float_t get_dca_bach_to_Lambda()  const { return ((Float_t)dca_bach_to_Lambda)/100.0; }
    Float_t get_dca_bach_to_baryon() const { return ((Float_t)dca_bach_to_baryon)/10.0; }
  
   
    
    Float_t get_CosPointingAngle()  const  {return ((Float_t)CosPointingAngle)/1000.;} //pay attention to cos! should not be negative
   // Int_t get_LeastNumberOfTPCclus()   const   {return LeastNumberOfTPCclus; }
   // Float_t get_MaxChi2perclus() const  {return ((Float_t)MaxChi2perclus)/100.; }
   // Float_t get_MinTrackLength() const   {return ((Float_t)MinTrackLength)/100.; }
    Float_t  get_CascadeDecayPos(Int_t i) const {return CascadeDecayPos[i];}
    Float_t  get_V0fromCascadePos(Int_t i) const {return V0fromCascadePos[i];}
    Float_t get_TrackLengthTPC(Int_t i) const{return (Float_t)TrackLengthTPC[i];}
    
    //ITS
    UShort_t get_ITSstatusPosTrack(Int_t i) const {return ITSstatusPosTrack[i];}
    UShort_t get_ITSstatusNegTrack(Int_t i) const {return ITSstatusNegTrack[i];}
    UShort_t get_ITSstatusBachTrack(Int_t i) const {return ITSstatusBachTrack[i];}
    Int_t get_ITSrefitFlag(Int_t i)const {return (Int_t)ITSrefitFlag[i];}
    Float_t get_ITSchi2(Int_t i)const {return ((Float_t)ITSchi2[i])/100.;}
    Float_t get_ITSsignal(Int_t i)const {return ((Float_t)ITSsignal[i])/10.;}
    Float_t get_nSigma_ITS_pos(Int_t i) const  { return ((Float_t)nSigma_ITS_pos[i])/100.;}
    Float_t get_nSigma_ITS_neg(Int_t i) const  { return ((Float_t)nSigma_ITS_neg[i])/100.;}
    Float_t get_nSigma_ITS_bach(Int_t i) const  {return ((Float_t)nSigma_ITS_bach[i])/100.;}
 /*   Float_t get_ITSPIDProbPos(Int_t i)const {return ((Float_t)ITSPIDProbPos[i])/10000.;}
    Float_t get_ITSPIDProbNeg(Int_t i)const {return ((Float_t)ITSPIDProbNeg[i])/10000.;}
    Float_t get_ITSPIDProbBach(Int_t i)const {return ((Float_t)ITSPIDProbBach[i])/10000.;}*/
    
    //TPC
    Float_t get_TPCsignal(Int_t i) const {return ((Float_t)TPCsignal[i])/10.;}
    Int_t get_TPCcls(Int_t i) const{return (Int_t)TPCcls[i];}
    Int_t get_TPCclsF(Int_t i) const{return (Int_t)TPCclsF[i];}
    Float_t get_TPCchi2(Int_t i) const{return ((Float_t)TPCchi2[i])/100.;}
    Float_t get_nSigma_dEdx_pos(Int_t i)  const { return ((Float_t)nSigma_dEdx_pos[i])/100.; }
    Float_t get_nSigma_dEdx_neg(Int_t i)  const {return ((Float_t)nSigma_dEdx_neg[i])/100.;}
    Float_t get_nSigma_dEdx_bach(Int_t i)const   { return ((Float_t)nSigma_dEdx_bach[i])/100.; }
  /*  Float_t get_TPCPIDProbPos(Int_t i)const {return ((Float_t)TPCPIDProbPos[i])/10000.;}
    Float_t get_TPCPIDProbNeg(Int_t i)const {return ((Float_t)TPCPIDProbNeg[i])/10000.;}
    Float_t get_TPCPIDProbBach(Int_t i)const {return ((Float_t)TPCPIDProbBach[i])/10000.;}*/
    Int_t get_TPC_PartVar_Pos(Int_t i) const {return (Int_t)posTPC[i];}
    Int_t get_TPC_PartVar_Neg(Int_t i) const {return (Int_t)negTPC[i];}
    Int_t get_TPC_PartVar_Bach(Int_t i) const {return (Int_t)bachTPC[i];}
    
    //TOF
    Float_t get_TOFsignal(Int_t i) const {return (Float_t)TOFsignal[i];}
    Int_t get_TOFrefitFlag(Int_t i)const {return (Int_t)TOFrefitFlag[i];}
    Float_t get_nSigma_TOF_pos(Int_t i) const  { return ((Float_t)nSigma_TOF_pos[i])/100.;}
    Float_t get_nSigma_TOF_neg(Int_t i) const  { return ((Float_t)nSigma_TOF_neg[i])/100.;}
    Float_t get_nSigma_TOF_bach(Int_t i) const  {return ((Float_t)nSigma_TOF_bach[i])/100.;}
    


  //  ClassDef(AliRunningCascadeTrack,2);
    
//----------------------------------------------------------------
};


//Class which saves the events.
class AliRunningCascadeEvent : public TObject
{
private:
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t N_tracks; // total number of tracks
    Int_t idi; // run id
    UInt_t periodnumber;
    Float_t   centrality;
    Bool_t    MVPPileUpFlag;
    Int_t     multiplicity;
    Long64_t  trigger_word;
    Short_t magfield;
    UShort_t fNumTracks;
    TClonesArray* fTracks; //->
    
public:
    //default constructor
    AliRunningCascadeEvent():x(-100),y(-100),z(-100),N_tracks(-1),
    idi(-1),periodnumber(-1),centrality(-1), MVPPileUpFlag(kFALSE), multiplicity(-1), trigger_word(0), magfield(0),
    fNumTracks(-1),fTracks(0)
    {
        
    }
    
    //explicit constructor
    AliRunningCascadeEvent(Int_t iii):x(-100),y(-100),z(-100),N_tracks(-1),
    idi(-1),periodnumber(-1),centrality(-1), MVPPileUpFlag(kFALSE), multiplicity(-1), trigger_word(0), magfield(0),
    fNumTracks(-1),fTracks(0)
    {
        fTracks = new TClonesArray( "AliRunningCascadeTrack", 10 );
        
    }
    //copy constructor
    AliRunningCascadeEvent(const AliRunningCascadeEvent&):TObject(),x(-100),y(-100),z(-100),N_tracks(-1),
    idi(-1),periodnumber(-1),centrality(-1), MVPPileUpFlag(kFALSE), multiplicity(-1), trigger_word(0), magfield(0),
    fNumTracks(-1),fTracks(0)
    {
        fTracks = new TClonesArray( "AliRunningCascadeTrack", 10 );
        
    }
    //copy assignment
    AliRunningCascadeEvent& operator =(const AliRunningCascadeEvent&)
    {
        return *this;
    }
    //destructor
    virtual ~AliRunningCascadeEvent() {
        delete fTracks;
        fTracks = NULL;
    }
    
    
//setters and getters
//-----------------------------------------------------------------------------
    void       setx(Float_t r)                    { x = r;                         }
    Float_t    getx() const                       { return x;                      }
    
    void       sety(Float_t r)                    { y = r;                         }
    Float_t    gety() const                       { return y;                      }
    
    void       setz(Float_t r)                    { z = r;                         }
    Float_t    getz() const                       { return z;                      }
    
    void       setN_tracks(Int_t r)               { N_tracks = r;                   }
    Int_t      getN_tracks() const                { return N_tracks;                }
    
    void       setid(Int_t  r)                    { idi = r;                        }
    Int_t      getid() const                      { return idi;                     }
    
    void       set_periodnumber(UInt_t  r)                    { periodnumber = r;                        }
    UInt_t     get_periodnumber() const                      { return periodnumber;                     }
    
    void       setcentrality(Float_t r)             {centrality  = r;                }
    Float_t    getcentrality() const              { return centrality;             }
    
    void       setMVPPileUpFlag(Bool_t r)             {MVPPileUpFlag  = r;                }
    Bool_t      getMVPPileUpFlag() const              { return MVPPileUpFlag;             }
    
    void       setmultiplicity(Int_t r)             {multiplicity  = r;                }
    Int_t      getmultiplicity() const              { return multiplicity;             }
    
    void       settrigger_word(Long64_t r)             { trigger_word = r;                }
    Long64_t      gettrigger_word() const              { return trigger_word;             }
    
    void       setmagfield(Short_t r)             { magfield = r;                }
    Long64_t      getmagfield() const              { return magfield;             }
    
    UShort_t getNumTracks() const        {return fNumTracks;}
//-----------------------------------------------------------------------------
    
    void ClearTrackList()
    {
        fNumTracks   = 0;
        fTracks      ->Clear();
    }
    
    
   AliRunningCascadeTrack* createTrack()
    {
        if (fNumTracks == fTracks->GetSize())
            fTracks->Expand( fNumTracks + 10 );
        if (fNumTracks >= 10000)
        {
            Fatal( "AliRunningCascadeEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
            exit( 2 );
        }
        AliRunningCascadeTrack* track = new ((*fTracks)[fNumTracks++])AliRunningCascadeTrack;
        return track;
    }
    
    AliRunningCascadeTrack* getTrack(UShort_t i) const
    {
        return i < fNumTracks ? (AliRunningCascadeTrack*)((*fTracks)[i]) : NULL;
    }
    
  //  ClassDef(AliRunningCascadeEvent, 2);
};
 


//=========================================================================================================


class AliAnalysisTaskStrangeCascadesDiscrete: public AliAnalysisTaskSE {
    
private:
    
    //variables which are used in addtask
    Bool_t fguard_CheckTrackQuality;
    Bool_t fguard_CheckCascadeQuality;
    Bool_t fguard_CheckTPCPID;
    Bool_t fkRunV0Vertexers;
    Bool_t fkRunVertexers;
    Bool_t fkUseLightVertexer;
    Bool_t fkUseOnTheFlyV0Cascading;
    
    //variables for creating objects
    AliPIDResponse *fPIDResponse;     //!
    AliESDtrackCuts *fESDtrackCuts;   //!
    AliESDtrackCuts *fESDtrackCutsITSsa2010;  //!
    AliESDtrackCuts *fESDtrackCutsGlobal2015; //!
    AliAnalysisUtils *fUtils;         //!
    TRandom3 *fRand; //!
    Bool_t fkDebugOOBPileup; //for studying the pileup with standard cuts

    //pointers
    AliAnalysisManager *man;                    
    AliInputEventHandler* inputHandler;         //!
    AliESDEvent *lESDevent;                     //!
    AliMultSelection *MultSelection;            
    AliCentrality* centrality;                  //!
    const AliESDVertex* lPrimaryBestESDVtx;     //!
    AliESDtrackCuts* esdtrackcuts;              //!
    AliESDcascade* xi;                          //!
    AliESDtrack *pTrackXi;                      //!
    AliESDtrack *nTrackXi;                      //!
    AliESDtrack *bachTrackXi;                   //!
    AliExternalTrackParam *hCascTraj;           //!
    AliESDtrack *lBaryonTrack;                  //!
    AliRunningCascadeTrack* Cascade_Track;      //!
    AliRunningCascadeEvent* Cascade_Event;      //!
    TTree* fTreeCascadeAsEvent;                 // (written to a file)
    
    //variables used in the userexec
    Float_t fMagneticField; //magnetic field in an event
    Float_t fPV_X; //event primary vertex x coordinate
    Float_t fPV_Y; //event primary vertex y coordinate
    Float_t fPV_Z; //event primary vertex z coordinate
    Float_t sigmamaxrunning;
    
    
    
    
    //variables which are not necessarily defined in the impl.file
   Double_t fV0VertexerSels[7];
    Double_t fCascadeVertexerSels[8]; // Array to store the 8 values for the different selections Casc. related
    
public:
    //constructors, copies and destructor
    AliAnalysisTaskStrangeCascadesDiscrete();
    AliAnalysisTaskStrangeCascadesDiscrete(Bool_t lRunV0Vertexers,
                                                  Bool_t lRunVertexers, //to rerun the cascade vertexers
                                                  Bool_t lUseLightVertexer, //use light cascade vertexer
                                                  Bool_t lUseOnTheFlyV0Cascading,
                                                  Bool_t lguard_CheckTrackQuality,
                                                  Bool_t lguard_CheckCascadeQuality,
                                                  Bool_t lguard_CheckTPCPID,
                                                  Double_t lV0MaxChi2,
                                                  Double_t lV0minDCAfirst,
                                                  Double_t lV0minDCAsecond,
                                                  Double_t lV0maxDCAdaughters,
                                                  Double_t lV0minCosAngle,
                                                  Double_t lV0minRadius,
                                                  Double_t lV0maxRadius,
                                                  Double_t lCascaderMaxChi2,
                                                  Double_t lCascaderV0MinImpactParam,
                                                  Double_t lCascaderV0MassWindow,
                                                  Double_t lCascaderBachMinImpactParam,
                                                  Double_t lCascaderMaxDCAV0andBach,
                                                  Double_t lCascaderMinCosAngle,
                                                  Double_t lCascaderMinRadius,
                                                  Double_t lCascaderMaxRadius,
                                                  Float_t sigmaRangeTPC,
                                                  const char *name);
    AliAnalysisTaskStrangeCascadesDiscrete(const AliAnalysisTaskStrangeCascadesDiscrete&);
    AliAnalysisTaskStrangeCascadesDiscrete& operator =(const AliAnalysisTaskStrangeCascadesDiscrete&);
    virtual ~AliAnalysisTaskStrangeCascadesDiscrete();
    
    //regular main members
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    
    //Task configurations
    
    
    //My members
    Bool_t GoodESDEvent(AliESDEvent* lESDevent, const AliESDVertex *lPrimaryBestESDVtx,AliESDtrackCuts* esdtrackcuts1);
    Bool_t GoodESDTrack(AliESDtrack* trackESD, Double_t raprange, Double_t etarange);
    Int_t GetITSstatus(AliESDtrack* esdtrack, Int_t layer) const;
    Bool_t ExtraCleanupCascade(AliESDcascade* Xi,AliESDtrack* gPtrack,AliESDtrack* gNtrack,AliESDtrack* gBtrack,Double_t Omegamasswindow,Double_t Ximasswindow);
    Bool_t GoodCandidatesTPCPID(AliESDtrack* pTrackXi, AliESDtrack* nTrackXi, AliESDtrack* bachTrackXi, Float_t sigmamax);
   // Int_t TPCPIDlogArray(AliESDtrack* pTrackXi, AliESDtrack* nTrackXi, AliESDtrack* bachTrackXi);
                               
                               
    
   
    ClassDef(AliAnalysisTaskStrangeCascadesDiscrete, 1);
};


#endif
