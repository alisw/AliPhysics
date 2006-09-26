#ifndef ROOT_TIsajet
#define ROOT_TIsajet

/**************************************************************************/
/*                                                                        */
/* TIsajet                                                                */
/*                                                                        */
/* This class implements an interface to the Isajet event generator.      */
/*                                                                        */
/**************************************************************************/

#ifndef ROOT_TGenerator
#include "TGenerator.h"
#include "AliRndm.h"
#endif

class TIsajet : public TGenerator , public AliRndm
{
//
 private:

     Char_t *title, *jobtype;
     Char_t* pdfpar[20];
     Float_t pdfval[20];
     Int_t num_Pdf;
     Int_t beam1_type, beam2_type;
     Float_t cutoff_mass;
     Float_t center_energy;
     Float_t frag_params[32];
     Char_t *jet1_type[30], *jet2_type[30], *jet3_type[30];
     Int_t num_jet_type[3];
     Float_t qcd_lambda;
     Bool_t forbid_decay, forbid_eta, forbid_evolve, forbid_hadron, forbid_pi0;
     Int_t generate_sigma;
     Float_t p_limits[6];
     Float_t phi_limits[6];
     Float_t pt_limits[6];
     Float_t theta_limits[6];
     Float_t x_limits[6];
     Float_t y_limits[6];
     Float_t peter_jet_frag[8];

     Bool_t setCutjet, setBeams, setFragment, setJettype1;
     Bool_t setJettype2, setJettype3, setLambda, setNodcay;
     Bool_t setNoeta, setNoevolve, setNohadron, setNopi0;
     Bool_t setNsigma, setP, setPhi, setPt, setTheta;
     Bool_t setX, setXgen, setY, setPdf, online;
     
 public:	

    TIsajet();
//  TIsajet(Float_t Energy_CM);
    virtual                ~TIsajet();
    
    virtual void           Initialise();

    virtual void           Reload();
    
    virtual void           RestoreDefaults();

    virtual Int_t          ImportParticles(TClonesArray *particles, Option_t *option = "");

    virtual void           GenerateEvent();

    virtual void           SetJobtype(Char_t *val);
    virtual void           GetJobtype() const;

    virtual void           SetOnline(Bool_t val);
    virtual Bool_t         GetOnline() const;

    virtual void           SetPDF(Char_t *name, Float_t val);
    
    virtual void           Isaini(Int_t& j, Int_t& k, Int_t& m, Int_t& n);

    virtual void           Isaevt(Int_t& j, Int_t& k, Int_t& m);

    virtual void           Openfiles();

    virtual void           PDFinit(Int_t& pdfpar, Int_t& pdfval);

    virtual void           Isabeg(Int_t& ifl);

    virtual void           Isabg2(Int_t& ifl);    
    
//   Parameters for the event. 
//   Variable explanations in Icommon.h
//   Common block DYLIM access routines :

    virtual void           SetQMIN(Float_t val);
    virtual Float_t        GetQMIN() const;
    
    virtual void           SetQMAX(Float_t val);
    virtual Float_t        GetQMAX() const;
    
    virtual void           SetQTMIN(Float_t val);
    virtual Float_t        GetQTMIN() const;
    
    virtual void           SetQTMAX(Float_t val);
    virtual Float_t        GetQTMAX() const;

    virtual void           SetYWMIN(Float_t val);
    virtual Float_t        GetYWMIN() const;
    
    virtual void           SetYWMAX(Float_t val);
    virtual Float_t        GetYWMAX() const;

    virtual void           SetYWLIMS();

// YWMIN and YWMAX default to a function of QMIN, QTMIN; they are recalculated
// whenever either value is set, unless they have been fixed by hand (ie using
// their setters).
    
    virtual void           SetXWMIN(Float_t val);
    virtual Float_t        GetXWMIN() const;
    
    virtual void           SetXWMAX(Float_t val);
    virtual Float_t        GetXWMAX() const;

    virtual void           SetTHWMIN(Float_t val);
    virtual Float_t        GetTHWMIN() const;
    
    virtual void           SetTHWMAX(Float_t val);
    virtual Float_t        GetTHWMAX() const;

    virtual void           SetTHWLIMS();

    virtual void           SetPHWMIN(Float_t val);
    virtual Float_t        GetPHWMIN() const;
    
    virtual void           SetPHWMAX(Float_t val);
    virtual Float_t        GetPHWMAX() const;

    virtual Bool_t         GetSETLMQ(Int_t index) const;
    
// Ends DYLIM
// Common block EEPAR access routines

    virtual void           SetPLEP(Float_t val);
    virtual Float_t        GetPLEP() const;
    
    virtual void           SetPLEM(Float_t val);
    virtual Float_t        GetPLEM() const;

    virtual void           SetRSHMIN(Float_t val);
    virtual Float_t        GetRSHMIN() const;
    
    virtual void           SetRSHMAX(Float_t val);
    virtual Float_t        GetRSHMAX() const;

    virtual void           SetUPSLON(Float_t val);
    virtual Float_t        GetUPSLON() const;

    virtual void           SetSIGZ(Float_t val);
    virtual Float_t        GetSIGZ() const;

    virtual Bool_t         GetIBREM() const;
    
    virtual Bool_t         GetIBEAM() const;

    virtual Float_t        GetSGMXEE() const;

// Ends EEPAR
// Common block FORCE access routines

    virtual Int_t          GetNFORCE() const;
    
    virtual void           SetIFORCE(const Int_t val[], Int_t arraySize, Bool_t anti = true);
    virtual void           UnForce(Int_t index, Bool_t anti = true);
    virtual void           UnForceID(Int_t particle_ID, Bool_t anti = true);    
//  If anti is false, the antiparticle's decay is not forced / unforced.
    
    virtual Int_t*         GetIFORCE(Int_t index) const;

    virtual Int_t          GetMEFORC(Int_t index) const;
    
// Ends FORCE
// Common block FRGPAR access routines

    virtual void           SetFRPAR(Float_t val, Int_t index);
    virtual void           SetAllFRPAR(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetFRPAR(Int_t index) const;

    virtual void           SetPUD(Float_t val);
    virtual Float_t        GetPUD() const;
    
    virtual void           SetPBARY(Float_t val);
    virtual Float_t        GetPBARY() const;
    
    virtual void           SetSIGQT(Float_t val);
    virtual Float_t        GetSIGQT() const;
    
    virtual void           SetPEND(Float_t val);
    virtual Float_t        GetPEND() const;
    
    virtual void           SetXGEN(Float_t val, Int_t index);
    virtual void           SetAllXGEN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetXGEN(Int_t index) const;

    virtual void           SetPSPIN1(Float_t val, Int_t index);
    virtual void           SetAllPSPIN1(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPSPIN1(Int_t index) const;

    virtual void           SetPMIX1(Float_t val, Int_t index1, Int_t index2);
    virtual void           SetAllPMIX1(const Float_t val[2][3]);
    virtual void           SetColumnPMIX1(const Float_t val[], Int_t col);
    virtual Float_t        GetPMIX1(Int_t index1, Int_t index2) const;
    
    virtual void           SetPMIX2(Float_t val, Int_t index1, Int_t index2);
    virtual void           SetAllPMIX2(const Float_t val[2][3]);
    virtual void           SetColumnPMIX2(const Float_t val[], Int_t col);
    virtual Float_t        GetPMIX2(Int_t index1, Int_t index2) const;

    virtual void           SetPMIXX1(Float_t val, Int_t index);
    virtual void           SetAllPMIXX1(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPMIXX1(Int_t index) const;
    
    virtual void           SetPMIXX2(Float_t val, Int_t index);
    virtual void           SetAllPMIXX2(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPMIXX2(Int_t index) const;
    
    virtual void           SetXGENSS(Float_t val, Int_t index);
    virtual void           SetAllXGENSS(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetXGENSS(Int_t index) const;
    
// Ends FRGPAR
// Common block HCON access routines

    virtual Float_t        GetANWWWW(Int_t index1, Int_t index2, Int_t index3) const;
    
    virtual Float_t        GetADWWWW(Int_t index1, Int_t index2) const;

    virtual Float_t        GetAIWWWW(Int_t index) const;

    virtual Float_t        GetHMASS() const;
    
    virtual Float_t        GetHGAM() const;    

    virtual Float_t        GetHGAMS(Int_t index) const;

    virtual Float_t        GetETAHGG() const;

    virtual Int_t          GetMATCHH(Int_t index) const;

    virtual Float_t        GetZSTARS(Int_t index1, Int_t index2) const;

    virtual void           SetIHTYPE(Int_t val);
    virtual void           SetIHTYPE(Char_t val[]);
    virtual Int_t          GetIHTYPE() const;
    
    virtual Float_t        GetHGAMSS(Int_t index1, Int_t index2) const;    

// Ends HCON
// Common block JETLIM access routines

    virtual void           SetPMIN(Float_t val, Int_t index);
    virtual void           SetAllPMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPMIN(Int_t index) const;

    virtual void           SetPMAX(Float_t val, Int_t index);
    virtual void           SetAllPMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPMAX(Int_t index) const;
    
    virtual void           SetPTMIN(Float_t val, Int_t index);
    virtual void           SetAllPTMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPTMIN(Int_t index) const;

    virtual void           SetPTMAX(Float_t val, Int_t index);
    virtual void           SetAllPTMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPTMAX(Int_t index) const;

    virtual void           SetYJMIN(Float_t val, Int_t index);
    virtual void           SetAllYJMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetYJMIN(Int_t index) const;

    virtual void           SetYJMAX(Float_t val, Int_t index);
    virtual void           SetAllYJMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetYJMAX(Int_t index) const;

    virtual void           SetYJLIMS();

// YJMIN and YJMAX default to a function of PTMIN; but if either has
// been set by hand, SetYJLIMS is not called when PTMIN is set.
    
    virtual void           SetPHIMIN(Float_t val, Int_t index);
    virtual void           SetAllPHIMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPHIMIN(Int_t index) const;

    virtual void           SetPHIMAX(Float_t val, Int_t index);
    virtual void           SetAllPHIMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetPHIMAX(Int_t index) const;

    virtual void           SetXJMIN(Float_t val, Int_t index);
    virtual void           SetAllXJMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetXJMIN(Int_t index) const;

    virtual void           SetXJMAX(Float_t val, Int_t index);
    virtual void           SetAllXJMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetXJMAX(Int_t index) const;

    virtual void           SetTHMIN(Float_t val, Int_t index);
    virtual void           SetAllTHMIN(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetTHMIN(Int_t index) const;

    virtual void           SetTHLIMS();
	
    virtual void           SetTHMAX(Float_t val, Int_t index);
    virtual void           SetAllTHMAX(const Float_t val[], Int_t arraySize);
    virtual Float_t        GetTHMAX(Int_t index) const;

    virtual Bool_t         GetSETLMJ(Int_t index) const;

// Ends JETLIM
// Common block JETPAR access routines

    virtual Float_t        GetP(Int_t index) const;

    virtual Float_t        GetPT(Int_t index) const;

    virtual Float_t        GetYJ(Int_t index) const;

    virtual Float_t        GetPHI(Int_t index) const;

    virtual Float_t        GetXJ(Int_t index) const;

    virtual Float_t        GetTH(Int_t index) const;

    virtual Float_t        GetCTH(Int_t index) const;

    virtual Float_t        GetSTH(Int_t index) const;

    virtual Int_t          GetJETTYP(Int_t index) const;

    virtual Float_t        GetSHAT() const;

    virtual Float_t        GetTHAT() const;

    virtual Float_t        GetUHAT() const;

    virtual Float_t        GetQSQ() const;

    virtual Float_t        GetX1() const;

    virtual Float_t        GetX2() const;

    virtual Float_t        GetPBEAM(Int_t index) const;

    virtual Float_t        GetQMW() const;

    virtual Float_t        GetQW() const;

    virtual Float_t        GetQTW() const;

    virtual Float_t        GetYW() const;

    virtual Float_t        GetXW() const;

    virtual Float_t        GetTHW() const;

    virtual Float_t        GetQTMW() const;

    virtual Float_t        GetPHIW() const;

    virtual Float_t        GetSHAT1() const;

    virtual Float_t        GetTHAT1() const;

    virtual Float_t        GetUHAT1() const;

    virtual void           SetJWTYP(Int_t val);
    virtual void           SetJWTYP(Char_t val[]);
    virtual Int_t          GetJWTYP() const;

    virtual Float_t        GetALFQSQ() const;

    virtual Float_t        GetCTHW() const;

    virtual Float_t        GetSTHW() const;

    virtual Float_t        GetQ0W() const;

    virtual Int_t          GetINITYP(Int_t index) const;

    virtual Int_t          GetISIGS() const;

    virtual Float_t        GetPBEAMS(Int_t index) const;
    
// Ends JETPAR
// Common block KKGRAV access routines

    virtual void           SetNEXTRAD(Int_t val);
    virtual Int_t          GetNEXTRAD() const;

    virtual void           SetMASSD(Float_t val);
    virtual Float_t        GetMASSD() const;

    virtual Float_t        GetKKGSD() const;

    virtual Float_t        GetSURFD() const;
    
    virtual void           SetUVCUT(Bool_t val);
    virtual Bool_t         GetUVCUT() const;
    
// Ends KKGRAV
// Common block MBGEN access routines

    virtual Float_t        GetPOMWT(Int_t index) const;

    virtual Float_t        GetPOMGEN(Int_t index) const;    
    
    virtual void           SetMNPOM(Int_t val);
    virtual Int_t          GetMNPOM() const;

    virtual void           SetMXPOM(Int_t val);
    virtual Int_t          GetMXPOM() const;

    virtual Float_t        GetPDIFFR() const;

    virtual Int_t          GetNPOM() const;

    virtual Float_t        GetXBARY(Int_t index) const;

    virtual Float_t        GetDXBARY(Int_t index) const;

    virtual Float_t        GetXPOM(Int_t index1, Int_t index2) const;

// Ends MBGEN
// Common block MGLIMS access routines

    virtual void           SetEHMGMN(Float_t val);
    virtual Float_t        GetEHMGMN() const;
    
    virtual void           SetEHMGMX(Float_t val);
    virtual Float_t        GetEHMGMX() const;

    virtual Float_t        GetYHMGMN() const;

    virtual Float_t        GetYHMGMX() const;

//  The eights in the All-setters correspond to MGLIMS.mxlim, but the
//  compiler will not recognize it here.
    
    virtual void           SetAMIJMN(Float_t val, Int_t index1, Int_t index2);
    virtual void           SetAllAMIJMN(const Float_t val[8][8]);
    virtual void           SetColumnAMIJMN(const Float_t val[], Int_t col);
    virtual Float_t        GetAMIJMN(Int_t index1, Int_t index2) const;

    virtual void           SetAMIJMX(Float_t val, Int_t index1, Int_t index2);
    virtual void           SetAllAMIJMX(const Float_t val[8][8]);    
    virtual void           SetColumnAMIJMX(const Float_t val[], Int_t col);
    virtual Float_t        GetAMIJMX(Int_t index1, Int_t index2) const;

    virtual Bool_t         GetFIXMIJ(Int_t index1, Int_t index2) const;
    
// End MGLIMS
// Common block NODCAY access routines

    virtual void           SetNODCAY(Bool_t val);
    virtual Bool_t         GetNODCAY() const;
    
    virtual void           SetNOETA(Bool_t val);
    virtual Bool_t         GetNOETA() const;

    virtual void           SetNOPI0(Bool_t val);
    virtual Bool_t         GetNOPI0() const;

    virtual void           SetNONUNU(Bool_t val);
    virtual Bool_t         GetNONUNU() const;

    virtual void           SetNOEVOL(Bool_t val);
    virtual Bool_t         GetNOEVOL() const;

    virtual void           SetNOHADR(Bool_t val);
    virtual Bool_t         GetNOHADR() const;

    virtual void           SetNOGRAV(Bool_t val);
    virtual Bool_t         GetNOGRAV() const;

// Ends NODCAY
// Common block PARTCL access routines (get-only block)
    
    virtual Int_t          GetNPTCL() const;

    virtual Float_t        GetPX(Int_t index) const;
    virtual Float_t        GetPY(Int_t index) const;
    virtual Float_t        GetPZ(Int_t index) const;
    virtual Float_t        GetP0(Int_t index) const;
    virtual Float_t        GetMASS(Int_t index) const;

    virtual Float_t        GetORIG(Int_t index) const;
    virtual Float_t        GetIDENT(Int_t index) const;
    virtual Float_t        GetIDCAY(Int_t index) const;
    
// Ends PARTCL
// Common block PRIMAR access routines

    virtual Int_t          GetNJET() const;

    virtual Float_t        GetSCM() const;

    virtual Float_t        GetHALFE() const;

    virtual void           SetECM(Float_t val);
    virtual Float_t        GetECM() const;

    virtual void           SetIDIN(Int_t val, Int_t index);
    virtual void           SetIDIN(const Char_t val[], Int_t index);
    virtual Int_t          GetIDIN(Int_t index) const;

    virtual Int_t          GetNEVENT() const;

    virtual void           SetNTRIES(Int_t val);
    virtual Int_t          GetNTRIES() const;

    virtual void           SetNSIGMA(Int_t val);
    virtual Int_t          GetNSIGMA()const;
    
// Ends PRIMAR
// Common block QCDPAR access routines    

    virtual void           SetALAM(Float_t val);
    virtual Float_t        GetALAM() const;

    virtual Float_t        GetALAM2() const;

    virtual void           SetCUTJET(Float_t val);
    virtual Float_t        GetCUTJET() const;

    virtual void           SetISTRUC(Int_t val);
    virtual void           SetISTRUC(const Char_t val[]);
    virtual Int_t          GetISTRUC() const;
    
// Ends QCDPAR
// Common block QLMASS access routines
    
// AMLEP has no All-setter for the good and simple reason that
// not all of its entries should be set. GetAnyAMLEP returns those
// indices that cannot be set by the user.
    
    
    virtual void           SetAMLEP(Float_t val, Int_t index);
    virtual Float_t        GetAnyAMLEP(Int_t index) const;
    virtual Float_t        GetAMLEP(Int_t index) const;

    virtual void           SetTquarkMass(Float_t val);
    virtual Float_t        GetTquarkMass() const;
    
    virtual void           SetXquarkMass(Float_t val);
    virtual Float_t        GetXquarkMass() const;

    virtual void           SetYquarkMass(Float_t val);
    virtual Float_t        GetYquarkMass() const;

    virtual void           SetUtildeMass(Float_t val);
    virtual Float_t        GetUtildeMass() const;

    virtual void           SetDtildeMass(Float_t val);
    virtual Float_t        GetDtildeMass() const;

    virtual void           SetStildeMass(Float_t val);
    virtual Float_t        GetStildeMass() const;

    virtual void           SetCtildeMass(Float_t val);
    virtual Float_t        GetCtildeMass() const;    

    virtual void           SetBtildeMass(Float_t val);
    virtual Float_t        GetBtildeMass() const;
 
    virtual void           SetTtildeMass(Float_t val);
    virtual Float_t        GetTtildeMass() const;
   
    virtual void           SetGtildeMass(Float_t val);
    virtual Float_t        GetGtildeMass() const;
    
    virtual void           SetGammatildeMass(Float_t val);
    virtual Float_t        GetGammatildeMass() const;

    virtual void           SetNuEtildeMass(Float_t val);
    virtual Float_t        GetNuEtildeMass() const;
    
    virtual void           SetEtildeMass(Float_t val);
    virtual Float_t        GetEtildeMass() const;

    virtual void           SetNuMutildeMass(Float_t val);
    virtual Float_t        GetNuMutildeMass() const;

    virtual void           SetMutildeMass(Float_t val);
    virtual Float_t        GetMutildeMass() const;

    virtual void           SetNuTautildeMass(Float_t val);
    virtual Float_t        GetNuTautildeMass() const;
    
    virtual void           SetTautildeMass(Float_t val);
    virtual Float_t        GetTautildeMass() const;

    virtual void           SetWplustildeMass(Float_t val);
    virtual Float_t        GetWplustildeMass() const;

    virtual void           SetZ0tildeMass(Float_t val);
    virtual Float_t        GetZ0tildeMass() const;
    
    virtual void           SetHiggsMesonMass(Float_t val, Int_t index);
    virtual Float_t        GetHiggsMesonMass(Int_t index) const;
    
    virtual Int_t          GetNQLEP() const;

    virtual Int_t          GetNMES() const;

    virtual Int_t          GetNBARY() const;

// Ends QLMASS
// Common block SEED access routines

    virtual void           SetSEED(const Char_t val[24]);
    virtual Char_t*        GetSEED() const;
    
// Ends SEED
// Common block SUGNU access routines    

    virtual void           SetXNUSUG(Float_t val, Int_t index);
    virtual Float_t        GetXNUSUG(Int_t index) const;

    virtual void           SetGauginoMass(Float_t val, Int_t index);
    virtual Float_t        GetGauginoMass(Int_t index) const;

    virtual void           SetAtau(Float_t val);
    virtual Float_t        GetAtau() const;
    
    virtual void           SetAb(Float_t val);
    virtual Float_t        GetAb() const;

    virtual void           SetAt(Float_t val);
    virtual Float_t        GetAt() const;

    virtual void           SetHiggsDmass(Float_t val);
    virtual Float_t        GetHiggsDmass() const;

    virtual void           SetHiggsUmass(Float_t val);
    virtual Float_t        GetHiggsUmass() const;

    virtual void           SetERmass(Float_t val);
    virtual Float_t        GetERmass() const;

    virtual void           SetELmass(Float_t val);
    virtual Float_t        GetELmass() const;

    virtual void           SetDRmass(Float_t val);
    virtual Float_t        GetDRmass() const;

    virtual void           SetURmass(Float_t val);
    virtual Float_t        GetURmass() const;

    virtual void           SetULmass(Float_t val);
    virtual Float_t        GetULmass() const;

    virtual void           SetTauRmass(Float_t val);
    virtual Float_t        GetTauRmass() const;

    virtual void           SetTauLmass(Float_t val);
    virtual Float_t        GetTauLmass() const;

    virtual void           SetBRmass(Float_t val);
    virtual Float_t        GetBRmass() const;

    virtual void           SetTRmass(Float_t val);
    virtual Float_t        GetTRmass() const;

    virtual void           SetTLmass(Float_t val);
    virtual Float_t        GetTLmass() const;

// Ends XNUSUG
// Common block TCPAR access routines

    virtual void           SetTCMRHO(Float_t val);
    virtual Float_t        GetTCMRHO() const;
    
    virtual void           SetTCGRHO(Float_t val);
    virtual Float_t        GetTCGRHO() const;

// Ends TCPAR
// Common block TYPES access routines

    virtual Int_t          GetLOC(Int_t index) const;

    virtual Int_t          GetNTYP() const;
    
    virtual Int_t          GetNJTTYP(Int_t index) const;

    virtual Int_t          GetNWWTYP(Int_t index) const;

    virtual Int_t          GetNWMODE(Int_t index) const;
    
// Ends TYPES
// Common block XMSSM access routines

    virtual Bool_t         GetGOMSSM() const;
    
    virtual Bool_t         GetGOSUG() const;

    virtual Bool_t         GetGOGMSB() const;

    virtual Bool_t         GetGOAMSB() const;

    virtual Bool_t         GetAL3UNI() const;

    virtual void           SetXGLSS(Float_t val);
    virtual Float_t        GetXGLSS() const;

    virtual void           SetXMUSS(Float_t val);
    virtual Float_t        GetXMUSS() const;

    virtual void           SetXHASS(Float_t val);
    virtual Float_t        GetXHASS() const;

    virtual void           SetXTBSS(Float_t val);
    virtual Float_t        GetXTBSS() const;

    virtual void           SetXQ1SS(Float_t val);
    virtual Float_t        GetXQ1SS() const;

    virtual void           SetXDRSS(Float_t val);
    virtual Float_t        GetXDRSS() const;

    virtual void           SetXURSS(Float_t val);
    virtual Float_t        GetXURSS() const;

    virtual void           SetXL1SS(Float_t val);
    virtual Float_t        GetXL1SS() const;

    virtual void           SetXERSS(Float_t val);
    virtual Float_t        GetXERSS() const;

    virtual void           SetXQ2SS(Float_t val);
    virtual Float_t        GetXQ2SS() const;

    virtual void           SetXSRSS(Float_t val);
    virtual Float_t        GetXSRSS() const;

    virtual void           SetXCRSS(Float_t val);
    virtual Float_t        GetXCRSS() const;

    virtual void           SetXL2SS(Float_t val);
    virtual Float_t        GetXL2SS() const;

    virtual void           SetXMRSS(Float_t val);
    virtual Float_t        GetXMRSS() const;

    virtual void           SetXQ3SS(Float_t val);
    virtual Float_t        GetXQ3SS() const;

    virtual void           SetXBRSS(Float_t val);
    virtual Float_t        GetXBRSS() const;

    virtual void           SetXTRSS(Float_t val);
    virtual Float_t        GetXTRSS() const;

    virtual void           SetXL3SS(Float_t val);
    virtual Float_t        GetXL3SS() const;

    virtual void           SetXTARSS(Float_t val);
    virtual Float_t        GetXTARSS() const;

    virtual void           SetXATSS(Float_t val);
    virtual Float_t        GetXATSS() const;

    virtual void           SetXABSS(Float_t val);
    virtual Float_t        GetXABSS() const;

    virtual void           SetXATASS(Float_t val);
    virtual Float_t        GetXATASS() const;

    virtual void           SetXM1SS(Float_t val);
    virtual Float_t        GetXM1SS() const;

    virtual void           SetXM2SS(Float_t val);
    virtual Float_t        GetXM2SS() const;

    virtual void           SetXM0SU(Float_t val);
    virtual Float_t        GetXM0SU() const;

    virtual void           SetXMHSU(Float_t val);
    virtual Float_t        GetXMHSU() const;

    virtual void           SetXA0SU(Float_t val);
    virtual Float_t        GetXA0SU() const;

    virtual void           SetXTGBSU(Float_t val);
    virtual Float_t        GetXTGBSU() const;

    virtual void           SetXSMUSU(Float_t val);
    virtual Float_t        GetXSMUSU() const;

    virtual void           SetXLAMGM(Float_t val);
    virtual Float_t        GetXLAMGM() const;

    virtual void           SetXMESGM(Float_t val);
    virtual Float_t        GetXMESGM() const;

    virtual void           SetXN5GM(Float_t val);
    virtual Float_t        GetXN5GM() const;

    virtual void           SetXCMGV(Float_t val);
    virtual Float_t        GetXCMGV() const;

    virtual void           SetMGVTO(Float_t val);
    virtual Float_t        GetMGVTO() const;

    virtual void           SetXRSLGM(Float_t val);
    virtual Float_t        GetXRSLGM() const;

    virtual void           SetXDHDGM(Float_t val);
    virtual Float_t        GetXDHDGM() const;

    virtual void           SetXDHUGM(Float_t val);
    virtual Float_t        GetXDHUGM() const;

    virtual void           SetXDYGM(Float_t val);
    virtual Float_t        GetXDYGM() const;

    virtual void           SetXN51GM(Float_t val);
    virtual Float_t        GetXN51GM() const;

    virtual void           SetXN52GM(Float_t val);
    virtual Float_t        GetXN52GM() const;

    virtual void           SetXN53GM(Float_t val);
    virtual Float_t        GetXN53GM() const;

    virtual void           SetXMN3NR(Float_t val);
    virtual Float_t        GetXMN3NR() const;

    virtual void           SetXMAJNR(Float_t val);
    virtual Float_t        GetXMAJNR() const;

    virtual void           SetXANSS(Float_t val);
    virtual Float_t        GetXANSS() const;

    virtual void           SetXNRSS(Float_t val);
    virtual Float_t        GetXNRSS() const;

    virtual void           SetXSBCS(Float_t val);
    virtual Float_t        GetXSBCS() const;

// Ends XMSSM
// Common block XTYPES access routines

    virtual Char_t*        GetPARTYP(Int_t index) const;

    virtual void           SetTITLE(Char_t *val);
    virtual Char_t*        GetTITLE() const;

    virtual void           SetJETYP(Int_t index, Char_t val[]);
    virtual Char_t*        GetJETYP(Int_t index1, Int_t index2) const;

    virtual void           SetWWTYP(Char_t val[], Int_t index1, Int_t index2);
    virtual void           SetAllWWTYP(Char_t* val[2][30]);
    virtual void           SetColumnWWTYP(Char_t* val[], Int_t col);
    virtual Char_t*        GetWWTYP(Int_t index1, Int_t index2) const;

    virtual void           SetWMODES(Char_t val[], Int_t index1, Int_t index2);
    virtual void           SetAllWMODES(Char_t* val[2][30]);
    virtual void           SetColumnWMODES(Char_t* val[], Int_t col);
    virtual Char_t*        GetWMODES(Int_t index1, Int_t index2) const;

// Ends XTYPES
// Common block WCON access routines

    virtual void           SetSIN2W(Float_t val);
    virtual Float_t        GetSIN2W() const;

    virtual void           SetWMASS(Float_t w, Float_t z);
    virtual Float_t        GetWMASS(Int_t index) const;    

    virtual void           SetWMass(Float_t val);

    virtual void           SetZMass(Float_t val);

    
    virtual Float_t        GetWGAM(Int_t index) const;
    
    virtual Float_t        GetAQ(Int_t index1, Int_t index2) const;

    virtual Float_t        GetBQ(Int_t index1, Int_t index2) const;

    virtual Float_t        GetCOUT(Int_t index) const;

    virtual Int_t          GetMATCH() const;

    virtual Float_t        GetWCBR(Int_t index1, Int_t index2) const;

    virtual void           SetCUTOFF(Float_t val);
    virtual Float_t        GetCUTOFF() const;
    
    virtual void           SetCUTPOW(Float_t val);
    virtual Float_t        GetCUTPOW() const;

    virtual Float_t        GetTBRWW(Int_t index1, Int_t index2) const;

    virtual Float_t        GetRBRWW(Int_t index1, Int_t index2, Int_t index3) const;

    virtual Float_t        GetEZ() const;

    virtual Float_t        GetAQDP(Int_t index1, Int_t index2) const;

    virtual Float_t        GetBQDP(Int_t index1, Int_t index2) const;

    virtual Float_t        GetEZDP() const;

    virtual void           SetWFUDGE(Float_t val);
    virtual Float_t        GetWFUDGE() const;
    
// Ends WCON

    ClassDef(TIsajet,1)
	};
	

  
#endif
    
// Endfile.


