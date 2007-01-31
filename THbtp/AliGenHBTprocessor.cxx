/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Implementation of the interface for THBTprocessor
// Author: Piotr Krzysztof Skowronski <Piotr.Skowronski@cern.ch>
//////////////////////////////////////////////////////////////////////////////////
//Wrapper class for "hbt processor" after burner
//The origibal code is written in fortran by Lanny Ray
//and is put in the directory $ALICE_ROOT/HBTprocessor
//Detailed description is on the top of the file hbt_event_processor.f
//
//This class can be used ONLY with AliGenCocktailAfterBurner wrapper generator
//i.e. (in Config.C)
// ....
// AliGenCocktailAfterBurner = gener = new AliGenCocktailAfterBurner();
// gener->SetPhiRange(0, 360); //Set global parameters
// gener->SetThetaRange(35., 145.); 
// AliGenHIJINGpara *hijing = new AliGenHIJINGpara(10000); //Add some generator
// hijing->SetMomentumRange(0, 999);   
// gener->AddGenerator(hijing,"HIJING PARAMETRIZATION",1); 
//
// AliGenHBTprocessor *hbtp = new AliGenHBTprocessor(); //create object
// hbtp->SetRefControl(2); //set parameters
// hbtp->SetSwitch1D(1);
// hbtp->SetR0(6);//fm - radius
// hbtp->SetLambda(0.7);//chaocity parameter
// gener->AddAfterBurner(hbtp,"HBT PROCESSOR",1); //add to the main generator
// 
// gener->Init();
//
//CAUTIONS: 
//         A)  HBT PROCESSOR NEEDS MORE THAN ONE EVENT TO WORK
//             AS MORE AS IT BETTER WORKS
//         B)  IT IS ABLE TO "ADD" CORRELATIONS ONLY UP TO TWO PARTICLE TYPES AT ONES
//
//
// Artificial particle dennsity enhancment feature
// HBT Processor is unable to process correctly low multiplicity particles
// even if the high statistics (more events) is supplied.
// For that reason was implemented artificial multiplicity enhancement 
// - see also comments in HBT Processor source file (HBTP/hbt_event_processor.f)
//   or web page (http://alisoft.cern.ch/people/skowron/hbtprocessor/hbteventprocessor.html) 
//   concerning trk_accep or  SetTrackRejectionFactor method in this class
// Fortran is cheated by masking several events a single one. Number is defined by fEventMerge
// variable. In case it is equal to 1, there is no masking and the feature is off.
// But, for example if it is 5, and we have 100 events, number of events passed to fortran is 100/5=20
// When fortran asks about multiplcity of event, let sey, 1 - multiplicity of events 5 to 9 is summed up 
// and returned. 
//
//
//////////////////////////////////////////////////////////////////////////////////

// 11.11.2001 Piotr Skowronski
// Setting seed (date) in RNG in the constructor

// 09.10.2001 Piotr Skowronski
// 
// Theta - Eta cohernecy facilities added:
//    AliGenerator::SetThetaRange method overriden
//    Static methods added
//    EtaToTheta
//    ThetaToEta 
//    DegreesToRadians
//    RadiansToDegrees
//
// Class description comments put on proper place

// 27.09.2001 Piotr Skowronski
// removing of redefinition of defaults velues 
// in method's implementation. 
//  
// 

#include "AliGenHBTprocessor.h"

#include <TParticle.h>
#include "THBTprocessor.h"

#include "AliStack.h"
#include "AliGenCocktailAfterBurner.h"
#include "AliLog.h"



ClassImp(AliGenHBTprocessor)

Int_t AliGenHBTprocessor::fgDebug = 0;
static TRandom* gAliRandom;//RNG to be used by the fortran code

AliGenCocktailAfterBurner*  GetGenerator();
/*******************************************************************/

AliGenHBTprocessor::AliGenHBTprocessor(): 
  AliGenerator(), 
  fHBTprocessor(0x0),
  fHbtPStatCodes(0x0),
  fEventMerge(1)
{
  //
  // Standard constructor
  // Sets default veues of all parameters

  SetName("AliGenHBTprocessor");
  SetTitle("AliGenHBTprocessor");
  
  gAliRandom = fRandom;
  fHBTprocessor = new THBTprocessor();

  fNPDGCodes = 0;
  DefineParticles();

  SetTrackRejectionFactor();
  SetRefControl();
  SetPIDs();
  SetNPIDtypes();
  SetDeltap();
  SetMaxIterations();
  SetDelChi();
  SetIRand();
  SetLambda();
  SetR1d() ;
  SetRSide();
  SetROut() ;
  SetRLong() ;
  SetRPerp();
  SetRParallel();
  SetR0();
  SetQ0();
  SetSwitch1D();
  SetSwitch3D();
  SetSwitchType();
  SetSwitchCoherence();
  SetSwitchCoulomb();
  SetSwitchFermiBose();
  //SetMomentumRange();
  SetPtRange();
  SetPxRange();
  SetPyRange(); 
  SetPzRange(); 
  SetPhiRange(); 
  SetEtaRange();  
  SetNPtBins();  
  SetNPhiBins();  
  SetNEtaBins();
  SetNPxBins(); 
  SetNPyBins();
  SetNPzBins(); 
  SetNBins1DFineMesh();
  SetBinSize1DFineMesh();
  SetNBins1DCoarseMesh();
  SetBinSize1DCoarseMesh();
  SetNBins3DFineMesh();
  SetBinSize3DFineMesh();
  SetNBins3DCoarseMesh();
  SetBinSize3DCoarseMesh();
  SetNBins3DFineProjectMesh();
}

/*******************************************************************/


/*******************************************************************/

AliGenHBTprocessor::~AliGenHBTprocessor()
{
//destructor
  CleanStatusCodes();
  if (fHBTprocessor) delete fHBTprocessor; //delete generator
  
}

/*******************************************************************/

void AliGenHBTprocessor::InitStatusCodes()
{
 //creates and inits status codes array to zero
  AliGenCocktailAfterBurner *cab = GetGenerator();

  if(!cab) Fatal("InitStatusCodes()","Can not find AliGenCocktailAfterBurner generator");

  Int_t nev = cab->GetNumberOfEvents();

  fHbtPStatCodes = new Int_t* [nev];
  for( Int_t i =0; i<nev;i++)
  { 
    Int_t nprim = cab->GetStack(i)->GetNprimary();
    fHbtPStatCodes[i] = new Int_t[nprim];
    for (Int_t k =0 ;k<nprim;k++)
      fHbtPStatCodes[i][k] =0;
    
  }
  
}
/*******************************************************************/

void AliGenHBTprocessor::CleanStatusCodes()
{
 //Cleans up status codes
  if (fHbtPStatCodes)
  {
    for (Int_t i =0; i<GetGenerator()->GetNumberOfEvents(); i++)
      delete [] fHbtPStatCodes[i];  
    delete fHbtPStatCodes;
    fHbtPStatCodes = 0;
  }

}
/*******************************************************************/

void AliGenHBTprocessor::Init()
  {  
  //sets up parameters in generator
   
   THBTprocessor *thbtp = fHBTprocessor;
   

   thbtp->SetTrackRejectionFactor(fTrackRejectionFactor);
   thbtp->SetRefControl(fReferenceControl);
   
   if ((fPid[0] == fPid[1]) || (fPid[0] == 0) || (fPid[1] == 0))
    {
       if (fPid[0] == 0)
         thbtp->SetPIDs(IdFromPDG(fPid[1]) ,0);
       else
         thbtp->SetPIDs(IdFromPDG(fPid[0]) ,0);
       thbtp->SetNPIDtypes(1);
       
       if (fSwitchType !=1)
          Warning("AliGenHBTprocessor::Init","\nThere is only one particle type set,\n\
                   and Switch_Type differnt then 1,\n which does not make sense.\n\
                   Setting it to 1.\n");
                   
       SetSwitchType(1);
    }
   else
    {
       thbtp->SetPIDs(IdFromPDG(fPid[0]) ,IdFromPDG(fPid[1]));
       SetNPIDtypes(2);
       thbtp->SetSwitchType(fSwitchType); 
    }
   
   thbtp->SetMaxIterations(fMaxit);
   thbtp->SetDelChi(fDelchi);
   thbtp->SetIRand(fIrand);
   thbtp->SetLambda(fLambda);
   thbtp->SetR1d(fR1d);
   thbtp->SetRSide(fRside);
   thbtp->SetROut(fRout);
   thbtp->SetRLong(fRlong);
   thbtp->SetRPerp(fRperp);
   thbtp->SetRParallel(fRparallel);
   thbtp->SetR0(fR0);
   thbtp->SetQ0(fQ0);
   thbtp->SetSwitch1D(fSwitch1d);
   thbtp->SetSwitch3D(fSwitch3d);
   thbtp->SetSwitchType(fSwitchType);
   thbtp->SetSwitchCoherence(fSwitchCoherence);
   thbtp->SetSwitchCoulomb(fSwitchCoulomb);
   thbtp->SetSwitchFermiBose(fSwitchFermiBose);
   thbtp->SetPtRange(fPtMin,fPtMax);
   thbtp->SetPxRange(fPxMin,fPxMax);
   thbtp->SetPyRange(fPyMin,fPyMax);
   thbtp->SetPzRange(fPzMin,fPzMax);
   thbtp->SetPhiRange(fPhiMin*180./((Float_t)TMath::Pi())+180.0, //casting is because if fPhiMin = 180.0 then
                      fPhiMax*180./((Float_t)TMath::Pi())+180.0);//TMath::Pi() != TMath::Pi()*fPhiMin/180.0,
   thbtp->SetEtaRange(fEtaMin,fEtaMax);
   thbtp->SetNPtBins(fNPtBins);
   thbtp->SetNPhiBins(fNPhiBins);
   thbtp->SetNEtaBins(fNEtaBins);
   thbtp->SetNPxBins(fNPxBins);
   thbtp->SetNPyBins(fNPyBins);
   thbtp->SetNPzBins(fNPzBins);
   thbtp->SetNBins1DFineMesh(fN1dFine);
   thbtp->SetBinSize1DFineMesh(fBinsize1dFine);
   thbtp->SetNBins1DCoarseMesh(fN1dCoarse);
   thbtp->SetBinSize1DCoarseMesh(fBinsize1dCoarse);
   thbtp->SetNBins3DFineMesh(fN3dFine);
   thbtp->SetBinSize3DFineMesh(fBinsize3dFine);
   thbtp->SetNBins3DCoarseMesh(fN3dCoarse);
   thbtp->SetBinSize3DCoarseMesh(fBinsize3dCoarse);
   thbtp->SetNBins3DFineProjectMesh(fN3dFineProject);
   
   thbtp->SetPrintFull(fPrintFull);
       
 }
/*******************************************************************/
  
void AliGenHBTprocessor::Generate()
 {
 //starts processig 
   AliGenCocktailAfterBurner* cab = GetGenerator();
   if (cab == 0x0)
    {
      Fatal("Generate()","AliGenHBTprocessor needs AliGenCocktailAfterBurner to be main generator");
    }
   if (cab->GetNumberOfEvents() <2)
    {
      Fatal("Generate()",
            "HBT Processor needs more than 1 event to work properly,\
             but there is only %d set", cab->GetNumberOfEvents());
    }
  
  
   fHBTprocessor->Initialize(); //reset all fields of common blocks 
                                   //in case there are many HBT processors
                  	               //run one after one (i.e. in AliCocktailAfterBurner)
                                   //between Init() called and Generate there might 
   Init();                         //be different instance running - be sure that we have our settings
   
   InitStatusCodes(); //Init status codes
   
   fHBTprocessor->GenerateEvent(); //Generates event
   
   CleanStatusCodes(); //Clean Status codes - thet are not needed anymore
 }
 
/*******************************************************************/


/*******************************************************************/
void AliGenHBTprocessor::GetParticles(TClonesArray * particles) const
 {
 //practically dumm
   if (!particles)
    {
      Error("GetParticles","Particles has to be initialized");
      return;
    } 
   fHBTprocessor->ImportParticles(particles);
 }

/*******************************************************************/

Int_t AliGenHBTprocessor::GetHbtPStatusCode(Int_t part) const
{
//returns the status code of the given particle in the active event
//see SetActiveEvent in the bottom of AliGenHBTprocessor.cxx
//and in AliCocktailAfterBurner
 Int_t ev, idx;
 GetTrackEventIndex(part,ev,idx);
 if ( (ev<0) || (idx<0) )
  {
    Error("GetHbtPStatusCode","GetTrackEventIndex returned error");
    return 0;
  }
 return fHbtPStatCodes[ev][idx];
  
}

/*******************************************************************/
void  AliGenHBTprocessor::SetHbtPStatusCode(Int_t hbtstatcode, Int_t part)
{
 //Sets the given status code to given particle number (part) in the active event
 Int_t ev, idx;
 GetTrackEventIndex(part,ev,idx);
 if ( (ev<0) || (idx<0) )
  {
    Error("GetHbtPStatusCode","GetTrackEventIndex returned error");
    return;
  }
 else fHbtPStatCodes[ev][idx] = hbtstatcode;
}

/*******************************************************************/

void AliGenHBTprocessor::SetTrackRejectionFactor(Float_t trf) //def 1.0
  {
   //Sets the Track Rejection Factor
   //variates in range 0.0 <-> 1.0
   //Describes the factor of particles rejected from the output.
   //Used only in case of low muliplicity particles e.g. lambdas.
   //Processor generates addisional particles and builds the 
   //correletions on such a statistics.
   //At the end these particels are left in the event according 
   //to this factor: 1==all particles are left
   //                0==all are removed 
   
    fTrackRejectionFactor=trf;
    fHBTprocessor->SetTrackRejectionFactor(trf);
  }
/*******************************************************************/

void AliGenHBTprocessor::SetRefControl(Int_t rc) //default 2
 {
 //Sets the Refernce Control Switch
 //switch wether read reference histograms from file =1
 //              compute from input events =2 - default
   fReferenceControl=rc;
   fHBTprocessor->SetRefControl(rc);
 }
/*******************************************************************/

void AliGenHBTprocessor::SetPIDs(Int_t pid1,Int_t pid2)
  {
   //default pi+ pi-
   //Sets PDG Codes of particles to be processed, default \\Pi^{+} and \\Pi^{-}
   //This method accepts PDG codes which is
   //in opposite to THBBProcessor which accepts GEANT PID
    if ( (pid1 == 0) && (pid2 == 0) )
    {
      Error("AliGenHBTprocessor::SetPIDs","Sensless Particle Types setting: 0 0, Ignoring\n");
    }
    
    fPid[0]=pid1;
    fPid[1]=pid2;
    
    if(pid1 == pid2)
       {
        fHBTprocessor->SetPIDs(IdFromPDG(pid1) ,0);
        SetNPIDtypes(1);
        SetSwitchType(1);
       }
    else
       { 
        fHBTprocessor->SetPIDs(IdFromPDG(pid1) ,IdFromPDG(pid2));
        SetNPIDtypes(2);
       }
  }
/*******************************************************************/

void AliGenHBTprocessor::SetNPIDtypes(Int_t npidt)
  {
    //Number ofparticle types to be processed - default 2
    //see AliGenHBTprocessor::SetPIDs
    
    fNPidTypes = npidt;
    fHBTprocessor->SetNPIDtypes(npidt);
  }
/*******************************************************************/

void AliGenHBTprocessor::SetDeltap(Float_t deltp) 
  {
  //default = 0.1 GeV
  //maximum range for random momentum shifts in GeV/c;
  //px,py,pz independent; Default = 0.1 GeV/c.
    fDeltap=deltp;
    fHBTprocessor->SetDeltap(deltp); 
  }
/*******************************************************************/

void AliGenHBTprocessor::SetMaxIterations(Int_t maxiter) 
  { 
  //maximum # allowed iterations thru all the 
  //tracks for each event. Default = 50.
  //If maxit=0, then calculate the correlations 
  //for the input set of events without doing the
  //track adjustment procedure.
    
    fMaxit=maxiter;
    fHBTprocessor->SetMaxIterations(maxiter);
  }

/*******************************************************************/
void AliGenHBTprocessor::SetDelChi(Float_t dc)
  {
    //min % change in total chi-square for which
    //the track adjustment iterations may stop,
    //Default = 0.1%.
    
    fDelchi=dc;
    fHBTprocessor->SetDelChi(dc);
  }
/*******************************************************************/

void AliGenHBTprocessor::SetIRand(Int_t irnd) 
  { 
    //if fact dummy - used only for compatibility
    //we are using AliRoot built in RNG
    //dummy in fact since we are using aliroot build-in RNG
    //Sets renaodom number generator
    fIrand=irnd;
    fHBTprocessor->SetIRand(irnd);
  }
/*******************************************************************/
      
void AliGenHBTprocessor::SetLambda(Float_t lam) 
  { 
  //default = 0.6
  // Sets Chaoticity Parameter
    fLambda = lam;
    fHBTprocessor->SetLambda(lam);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetR1d(Float_t r) 
  {
    //default 7.0
    //Sets Spherical source model radius (fm)
    fR1d = r;
    fHBTprocessor->SetR1d(r);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetRSide(Float_t rs) 
  {
   //default rs = 6.0
   //Rside,Rout,Rlong  = Non-spherical Bertsch-Pratt source model (fm)
   
    fRside = rs;
    fHBTprocessor->SetRSide(rs);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetROut(Float_t ro) 
  {
    //default ro = 7.0
    //Rside,Rout,Rlong  = Non-spherical Bertsch-Pratt source model (fm)
    fRout = ro;
    fHBTprocessor->SetROut(ro);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetRLong(Float_t rl) 
  {
    //default rl = 4.0
    //Rside,Rout,Rlong  = Non-spherical Bertsch-Pratt source model (fm)
    fRlong = rl;
    fHBTprocessor->SetRLong(rl);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetRPerp(Float_t rp) 
  {
   //default rp = 6.0
   //Rperp,Rparallel,R0= Non-spherical Yano-Koonin-Podgoretski source model (fm).
    fRperp = rp;
    fHBTprocessor->SetRPerp(rp);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetRParallel(Float_t rprl) 
  { 
   //default rprl = 4.0
   //Rperp,Rparallel,R0= Non-spherical Yano-Koonin-Podgoretski source model (fm).
    fRparallel = rprl;
    fHBTprocessor->SetRParallel(rprl);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetR0(Float_t r0) 
  {
  //default r0 = 4.0
  //Rperp,Rparallel,R0= Non-spherical Yano-Koonin-Podgoretski source model (fm).
    fR0 = r0;
    fHBTprocessor->SetR0(r0);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetQ0(Float_t q0) 
  { 
  //default q0 = 9.0
  //Sets Q0                = NA35 Coulomb parameter for finite source size in (GeV/c)
  //                         if fSwitchCoulomb = 2
  //                       = Spherical Coulomb source radius in (fm) 
  //                         if switchCoulomb = 3, used to interpolate the
  //                         input Pratt/Cramer discrete Coulomb source
  //                         radii tables.
    fQ0 = q0;
    fHBTprocessor->SetQ0(q0);
  }

/*******************************************************************/
void AliGenHBTprocessor::SetSwitch1D(Int_t s1d) 
  {
//default s1d = 3
// Sets fSwitch1d   
//                          = 0 to not compute the 1D two-body //orrelations.
//                          = 1 to compute this using Q-invariant
//                          = 2 to compute this using Q-total
//                          = 3 to compute this using Q-3-ve//tor

    fSwitch1d = s1d;
    fHBTprocessor->SetSwitch1D(s1d);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetSwitch3D(Int_t s3d) 
  {
//default s3d = 0
// Sets fSwitch3d
//                         = 0 to not compute the 3D two-body correlations.
//                         = 1 to compute this using the side-out-long form
//                         = 2 to compute this using the Yanno-Koonin-Pogoredskij form.   
     
    fSwitch3d = s3d;
    fHBTprocessor->SetSwitch3D(s3d);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetSwitchType(Int_t st)
  {
//default st = 3
//  Sets  switch_type       = 1 to fit only the like pair correlations
//                          = 2 to fit only the unlike pair correlations
//                          = 3 to fit both the like and unlike pair correl.
//See SetPIDs and Init
//If only one particle type is set, unly==1 makes sens
  
    fSwitchType = st;
    fHBTprocessor->SetSwitchType(st);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetSwitchCoherence(Int_t sc)
  {
// default  sc = 0
//        switchCoherence  = 0 for purely incoherent source (but can have
//                              lambda < 1.0)
//                          = 1 for mixed incoherent and coherent source
  
    fSwitchCoherence = sc;
    fHBTprocessor->SetSwitchCoherence(sc);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetSwitchCoulomb(Int_t scol) 
  {
//default scol = 2
//        switchCoulomb    = 0 no Coulomb correction
//                          = 1 Point source Gamow correction only
//                          = 2 NA35 finite source size correction
//                          = 3 Pratt/Cramer finite source size correction;
//                              interpolated from input tables.
    fSwitchCoulomb =scol;
    fHBTprocessor->SetSwitchCoulomb(scol);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetSwitchFermiBose(Int_t sfb)
  {
//default sfb = 1
//        switchFermiBose =  1 Boson pairs
//                          = -1 Fermion pairs

    fSwitchFermiBose = sfb;
    fHBTprocessor->SetSwitchFermiBose(sfb);
  }
/*******************************************************************/
void AliGenHBTprocessor::SetPtRange(Float_t ptmin, Float_t ptmax)
 {
// default ptmin = 0.1, ptmax = 0.98
//Sets Pt range (GeV)
   AliGenerator::SetPtRange(ptmin,ptmax);
   fHBTprocessor->SetPtRange(ptmin,ptmax);
 }

/*******************************************************************/
void AliGenHBTprocessor::SetPxRange(Float_t pxmin, Float_t pxmax)
 {
//default pxmin = -1.0, pxmax = 1.0
//Sets Px range 
  fPxMin =pxmin;
  fPxMax =pxmax;
  fHBTprocessor->SetPxRange(pxmin,pxmax);
 }
/*******************************************************************/
void AliGenHBTprocessor::SetPyRange(Float_t pymin, Float_t pymax)
 {
//default  pymin = -1.0, pymax = 1.0
//Sets Py range 
  fPyMin =pymin;
  fPyMax =pymax;
   fHBTprocessor->SetPyRange(pymin,pymax);
 }
/*******************************************************************/
void AliGenHBTprocessor::SetPzRange(Float_t pzmin, Float_t pzmax)
 {
//default pzmin = -3.6, pzmax = 3.6
//Sets Py range
   fPzMin =pzmin;
   fPzMax =pzmax; 
   fHBTprocessor->SetPzRange(pzmin,pzmax);
 }
void AliGenHBTprocessor::SetMomentumRange(Float_t /*pmin*/, Float_t /*pmax*/)
 {
 //default pmin=0, pmax=0
 //Do not use this method!
    MayNotUse("AliGenHBTprocessor::SetMomentumRange Method is Dummy");
 }
 
 /*******************************************************************/
void AliGenHBTprocessor::SetPhiRange(Float_t phimin, Float_t phimax)
 {
//
//Sets \\Phi range  
  AliGenerator::SetPhiRange(phimin,phimax);
  
  fHBTprocessor->SetPhiRange(phimin+180.0,phimax+180.0);
 }
/*******************************************************************/
void AliGenHBTprocessor::SetEtaRange(Float_t etamin, Float_t etamax)
 {
//default etamin = -1.5, etamax = 1.5
//Sets \\Eta range   
   fEtaMin= etamin;
   fEtaMax =etamax;
   fHBTprocessor->SetEtaRange(etamin,etamax);
   
   //set the azimothal angle range in the AliGeneraor - 
   //to keep coherency between azimuthal angle and pseudorapidity
   //DO NOT CALL this->SetThetaRange, because it calls this method (where we are) 
   //which must cause INFINITE LOOP
   AliGenerator::SetThetaRange(RadiansToDegrees(EtaToTheta(fEtaMin)), 
                               RadiansToDegrees(EtaToTheta(fEtaMax)));
   
 }
/*******************************************************************/
void AliGenHBTprocessor::SetThetaRange(Float_t thetamin, Float_t thetamax)
{
  //default thetamin = 0, thetamax = 180
  //Azimuthal angle, override AliGenerator method which uses widely (i.e. wrapper generators)
  //core fortran HBTProcessor uses Eta (pseudorapidity)
  //so these methods has to be synchronized
  
  AliGenerator::SetThetaRange(thetamin,thetamax);
  
  SetEtaRange( ThetaToEta(fThetaMin) , ThetaToEta(fThetaMax) );

}
  
/*******************************************************************/
void AliGenHBTprocessor::SetNPtBins(Int_t nptbin)
 {
  //default nptbin = 50
  //set number of Pt bins  
   fNPtBins= nptbin; 
   fHBTprocessor->SetNPtBins(nptbin);
 }
/*******************************************************************/
void AliGenHBTprocessor::SetNPhiBins(Int_t nphibin)
{ 
  //default nphibin = 50
  //set number of Phi bins
  fNPhiBins=nphibin;
  fHBTprocessor->SetNPhiBins(nphibin);
}
/*******************************************************************/
void AliGenHBTprocessor::SetNEtaBins(Int_t netabin)
{
  //default netabin = 50
  //set number of Eta bins
  fNEtaBins = netabin;
  fHBTprocessor->SetNEtaBins(netabin);
}
/*******************************************************************/
void AliGenHBTprocessor::SetNPxBins(Int_t npxbin)
{
  //default  npxbin = 20
  //set number of Px bins
  fNPxBins = npxbin; 
  fHBTprocessor->SetNPxBins(npxbin);
}
/*******************************************************************/
void AliGenHBTprocessor::SetNPyBins(Int_t npybin)
{
  //default  npybin = 20
  //set number of Py bins
  fNPyBins = npybin;
  fHBTprocessor->SetNPyBins(npybin);
}
/*******************************************************************/
void AliGenHBTprocessor::SetNPzBins(Int_t npzbin)
{
  //default npzbin = 70
  //set number of Pz bins
  fNPzBins = npzbin;
  fHBTprocessor->SetNPzBins(npzbin);
}
/*******************************************************************/
void AliGenHBTprocessor::SetNBins1DFineMesh(Int_t n)
{
//default n = 10
//Sets the number of bins in the 1D mesh
   fN1dFine =n;
   fHBTprocessor->SetNBins1DFineMesh(n);
   
}
/*******************************************************************/
void AliGenHBTprocessor::SetBinSize1DFineMesh(Float_t x)
{
//default x=0.01
//Sets the bin size in the 1D mesh
   fBinsize1dFine = x;
   fHBTprocessor->SetBinSize1DFineMesh(x);
}
/*******************************************************************/
      
void AliGenHBTprocessor::SetNBins1DCoarseMesh(Int_t n)
{
//default n =2
//Sets the number of bins in the coarse 1D mesh
  fN1dCoarse =n;
  fHBTprocessor->SetNBins1DCoarseMesh(n);
}
/*******************************************************************/
void AliGenHBTprocessor::SetBinSize1DCoarseMesh(Float_t x)
{
//default x=0.05
//Sets the bin size in the coarse 1D mesh
  fBinsize1dCoarse =x;
  fHBTprocessor->SetBinSize1DCoarseMesh(x);
}
/*******************************************************************/
      
void AliGenHBTprocessor::SetNBins3DFineMesh(Int_t n)
{
//default n = 8
//Sets the number of bins in the 3D mesh
  fN3dFine =n;
  fHBTprocessor->SetNBins3DFineMesh(n);
}
/*******************************************************************/
void AliGenHBTprocessor::SetBinSize3DFineMesh(Float_t x)
{
//default x=0.01
//Sets the bin size in the 3D mesh
  fBinsize3dFine =x;
  fHBTprocessor->SetBinSize3DFineMesh(x);
}
/*******************************************************************/
      
void AliGenHBTprocessor::SetNBins3DCoarseMesh(Int_t n)
{
//default n = 2
//Sets the number of bins in the coarse 3D mesh

  fN3dCoarse = n;
  fHBTprocessor->SetNBins3DCoarseMesh(n);
}
/*******************************************************************/
void AliGenHBTprocessor::SetBinSize3DCoarseMesh(Float_t x)
{
//default x=0.08
//Sets the bin size in the coarse 3D mesh
  fBinsize3dCoarse = x;
  fHBTprocessor->SetBinSize3DCoarseMesh(x);
}
/*******************************************************************/
      
void AliGenHBTprocessor::SetNBins3DFineProjectMesh(Int_t n )
{
//default n =3
//Sets the number of bins in the fine project mesh
  fN3dFineProject = n;
  fHBTprocessor->SetNBins3DFineProjectMesh(n);
}
/*******************************************************************/
void AliGenHBTprocessor::SetPrintFull(Int_t flag)
{
//sets the print full flag
 fPrintFull = flag;
 fHBTprocessor->SetPrintFull(flag);
}


/*******************************************************************/

Int_t  AliGenHBTprocessor::GetNumberOfEvents()
{
//returns number of available events
  AliGenerator* g = gAlice->Generator();
  AliGenCocktailAfterBurner* cab = (g)?dynamic_cast<AliGenCocktailAfterBurner*>(g):0x0;
  if (cab == 0x0)
   {
     Fatal("GetNumberOfEvents","Master Generator is not an AliGenCocktailAfterBurner");
     return 0;
   }
  return (Int_t)TMath::Ceil(cab->GetNumberOfEvents()/((Float_t)fEventMerge));
}

/*******************************************************************/

void AliGenHBTprocessor::SetActiveEventNumber(Int_t n)
{
//sets the active event
 fActiveStack =  n*fEventMerge;
 AliDebug(1,Form("Settimg active event %d passed %d",fActiveStack,n));
}
/*******************************************************************/

Int_t  AliGenHBTprocessor::GetNumberOfTracks()
{
//returns number of tracks in active event
  AliGenerator* g = gAlice->Generator();
  AliGenCocktailAfterBurner* cab = (g)?dynamic_cast<AliGenCocktailAfterBurner*>(g):0x0;
  if (cab == 0x0)
   {
     Fatal("GetNumberOfEvents","Master Generator is not an AliGenCocktailAfterBurner");
     return 0;
   }
 Int_t n = 0;
 for (Int_t i = fActiveStack;i < fActiveStack+fEventMerge; i++) 
  { 
    if (i >= GetNumberOfEvents()) break; //protection not to overshoot nb of events
    AliStack* stack = cab->GetStack(i);
    if (stack == 0x0)
     Error("GetNumberOfTracks","There is no stack %d",i);

    n+=stack->GetNprimary();
  }
 return n;
}
/*******************************************************************/

void AliGenHBTprocessor::SetNEventsToMerge(Int_t nev)
{
 //sets number of events to merge
 if (nev > 0 ) fEventMerge = nev;
}
/*******************************************************************/

TParticle* AliGenHBTprocessor::GetTrack(Int_t n)
{ 
//returns track that hbtp thinks is n in active event
  AliDebug(5,Form("n = %d",n));
  AliGenerator* g = gAlice->Generator();
  AliGenCocktailAfterBurner* cab = (g)?dynamic_cast<AliGenCocktailAfterBurner*>(g):0x0;
  if (cab == 0x0)
   {
     Fatal("GetTrackEventIndex","Master Generator is not an AliGenCocktailAfterBurner");
     return 0;
   }

 Int_t ev, idx;
 GetTrackEventIndex(n,ev,idx);
 AliDebug(5,Form("Event = %d Particle = %d",ev,idx));
 if ( (ev<0) || (idx<0) )
  {
    Error("GetTrack","GetTrackEventIndex returned error");
    return 0x0;
  }
 AliDebug(5,Form("Number of Tracks in Event(%d) = %d",ev,cab->GetStack(ev)->GetNprimary()));
 return cab->GetStack(ev)->Particle(idx);//safe - in case stack does not exist 
}
/*******************************************************************/

void AliGenHBTprocessor::GetTrackEventIndex(Int_t n, Int_t &evno, Int_t &index) const
{
 //returns event(stack) number and particle index
  AliGenerator* g = gAlice->Generator();
  AliGenCocktailAfterBurner* cab = (g)?dynamic_cast<AliGenCocktailAfterBurner*>(g):0x0;
  if (cab == 0x0)
   {
     Fatal("GetTrackEventIndex","Master Generator is not an AliGenCocktailAfterBurner");
     return;
   }

 evno = -1;
 index = -1;
 for (Int_t i = fActiveStack;i < fActiveStack+fEventMerge; i++) 
  { 
    AliStack* stack = cab->GetStack(i);
    if (stack == 0x0)
     {
       Error("GetTrackEventIndex","There is no stack %d",i);
       return;
     }

    Int_t ntracks = stack->GetNprimary();
    AliDebug(10,Form("Event %d has %d tracks. n = %d",i,ntracks,n));
    
    if ( ntracks > n) 
     {
       evno = i;
       index = n;
       return ;
     }
    else 
     {  
       n-=ntracks;
       continue;
     }
  }
 Error("GetTrackEventIndex","Could not find given track");
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/





void AliGenHBTprocessor::DefineParticles()
{
  //
  // Load standard numbers for GEANT particles and PDG conversion
  fNPDGCodes = 0; //this is done in the constructor - but in any case ...
  
  fPDGCode[fNPDGCodes++]=-99;   //  0 = unused location
  fPDGCode[fNPDGCodes++]=22;    //  1 = photon
  fPDGCode[fNPDGCodes++]=-11;   //  2 = positron
  fPDGCode[fNPDGCodes++]=11;    //  3 = electron
  fPDGCode[fNPDGCodes++]=12;    //  4 = neutrino e
  fPDGCode[fNPDGCodes++]=-13;   //  5 = muon +
  fPDGCode[fNPDGCodes++]=13;    //  6 = muon -
  fPDGCode[fNPDGCodes++]=111;   //  7 = pi0
  fPDGCode[fNPDGCodes++]=211;   //  8 = pi+
  fPDGCode[fNPDGCodes++]=-211;  //  9 = pi-
  fPDGCode[fNPDGCodes++]=130;   // 10 = Kaon Long
  fPDGCode[fNPDGCodes++]=321;   // 11 = Kaon +
  fPDGCode[fNPDGCodes++]=-321;  // 12 = Kaon -
  fPDGCode[fNPDGCodes++]=2112;  // 13 = Neutron
  fPDGCode[fNPDGCodes++]=2212;  // 14 = Proton
  fPDGCode[fNPDGCodes++]=-2212; // 15 = Anti Proton
  fPDGCode[fNPDGCodes++]=310;   // 16 = Kaon Short
  fPDGCode[fNPDGCodes++]=221;   // 17 = Eta
  fPDGCode[fNPDGCodes++]=3122;  // 18 = Lambda
  fPDGCode[fNPDGCodes++]=3222;  // 19 = Sigma +
  fPDGCode[fNPDGCodes++]=3212;  // 20 = Sigma 0
  fPDGCode[fNPDGCodes++]=3112;  // 21 = Sigma -
  fPDGCode[fNPDGCodes++]=3322;  // 22 = Xi0
  fPDGCode[fNPDGCodes++]=3312;  // 23 = Xi-
  fPDGCode[fNPDGCodes++]=3334;  // 24 = Omega-
  fPDGCode[fNPDGCodes++]=-2112; // 25 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3122; // 26 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3222; // 27 = Anti Sigma -
  fPDGCode[fNPDGCodes++]=-3212; // 28 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3112; // 29 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3322; // 30 = Anti Xi 0
  fPDGCode[fNPDGCodes++]=-3312; // 31 = Anti Xi +
  fPDGCode[fNPDGCodes++]=-3334; // 32 = Anti Omega +
}  

/*******************************************************************/
Int_t AliGenHBTprocessor::IdFromPDG(Int_t pdg) const 
{
  //
  // Return Geant3 code from PDG and pseudo ENDF code
  //
  for(Int_t i=0;i<fNPDGCodes;++i)
    if(pdg==fPDGCode[i]) return i;
  return -1;
}
Int_t AliGenHBTprocessor::PDGFromId(Int_t id) const
{
  //
  // Return PDG code and pseudo ENDF code from Geant3 code
  //
  if(id>0 && id<fNPDGCodes) return fPDGCode[id];
  else return -1;
}
Double_t AliGenHBTprocessor::ThetaToEta(Double_t arg)
 {
  //converts etha(azimuthal angle) to Eta (pseudorapidity). Argument in radians

   if(arg>= TMath::Pi()) return  708.0;//This number is the biggest wich not crashes exp(x)
   if(arg<= 0.0) return -708.0;//
   if(arg == TMath::Pi()/2.) return 0.0;//
   
   if (arg > 0.0) 
    { 
      return TMath::Log( TMath::Tan(arg/2.)) ;
    }
   else 
    { 
      return -TMath::Log( TMath::Tan(-arg/2.)) ;
    }
 }
                                  
/*******************************************************************/
/******      ROUTINES    USED    FOR     COMMUNUCATION      ********/
/********************     WITH      FORTRAN     ********************/
/*******************************************************************/

#ifndef WIN32
  # define hbtpran hbtpran_  
  # define alihbtp_puttrack alihbtp_puttrack_
  # define alihbtp_gettrack alihbtp_gettrack_
  # define alihbtp_getnumberevents alihbtp_getnumberevents_
  # define alihbtp_getnumbertracks  alihbtp_getnumbertracks_
  # define alihbtp_initialize alihbtp_initialize_
  # define alihbtp_setactiveeventnumber alihbtp_setactiveeventnumber_
  # define alihbtp_setparameters alihbtp_setparameters_
  # define type_ofCall

#else
  # define hbtpran HBTPRAN
  # define alihbtp_puttrack ALIHBTP_PUTTRACK
  # define alihbtp_gettrack ALIHBTP_GETTRACK
  # define alihbtp_getnumberevents ALIHBTP_GETNUMBEREVENTS
  # define alihbtp_getnumbertracks  ALIHBTP_GETNUMBERTRACKS
  # define alihbtp_initialize ALIHBTP_INITIALIZE
  # define alihbtp_setactiveeventnumber ALIHBTP_SETACTIVEEVENTNUMBER
  # define alihbtp_setparameters ALIHBTP_SETPARAMETERS
  # define type_ofCall  _stdcall
#endif    

/*******************************************************************/

AliGenCocktailAfterBurner*  GetGenerator()
 {
   // This function has two tasks:
   // Check if environment is OK (exist gAlice and generator)
   // Returns pointer to genrator
   //to be changed with TFolders

   if(!gAlice)
    {
      ::Error("AliGenHBTprocessor.cxx: GetGenerator()",
              "There is NO gALICE! Check what you are doing!");
      ::Fatal("AliGenHBTprocessor.cxx: GetGenerator()",
              "Running HBT Processor without gAlice... Exiting \n");
      return 0x0;//pro forma
    }
   AliGenerator * gen = gAlice->Generator();
   
   if (!gen) 
    {
      ::Fatal("AliGenHBTprocessor.cxx: GetGenerator()",
               "There is no generator in gAlice, exiting\n");
      return 0x0;
    }

   //we do not sure actual type of the genetator
   //and simple casting is risky - we use ROOT machinery to do safe cast
   
   TClass* cabclass = AliGenCocktailAfterBurner::Class(); //get AliGenCocktailAfterBurner TClass
   TClass* genclass = gen->IsA();//get TClass of the generator we got from galice 
   //use casting implemented in TClass
   //cast gen to cabclass
   AliGenCocktailAfterBurner* cab=(AliGenCocktailAfterBurner*)genclass->DynamicCast(cabclass,gen);
                                                                        
   if (cab == 0x0)//if generator that we got is not AliGenCocktailAfterBurner or its descendant we get null
   {              //then quit with error
      ::Fatal("AliGenHBTprocessor.cxx: GetGenerator()",
              "The main Generator is not a AliGenCocktailAfterBurner, exiting\n");
      return 0x0;
   }
   return cab;
 }
/*******************************************************************/

AliGenHBTprocessor* GetAliGenHBTprocessor()
{
//returns pointer to the current instance of AliGenHBTprocessor in
//AliGenCocktailAfterBurner (can be more than one)
//
 AliGenCocktailAfterBurner* gen = GetGenerator();
 AliGenerator* g = gen->GetCurrentGenerator();
 if(g == 0x0)
  {
    ::Fatal("AliGenHBTprocessor.cxx: GetAliGenHBTprocessor()",
                  "Can not get the current generator. Exiting");
    return 0x0;//pro forma
  }
  
 TClass* hbtpclass = AliGenHBTprocessor::Class(); //get AliGenCocktailAfterBurner TClass
 TClass* gclass = g->IsA();//get TClass of the current generator we got from CAB
 AliGenHBTprocessor* hbtp = (AliGenHBTprocessor*)gclass->DynamicCast(hbtpclass,g);//try to cast 
 if (hbtp == 0x0)
   {
      ::Fatal("AliGenHBTprocessor.cxx: GetAliGenHBTprocessor()",
              "\nCurrernt generator in AliGenCocktailAfterBurner is not a AliGenHBTprocessor, exiting\n");
      return 0x0;
   }
 return hbtp;
}

/*******************************************************************/
extern "C" void type_ofCall alihbtp_setparameters()
 {
   //dummy
 }

extern "C" void type_ofCall  alihbtp_initialize()
 {
   //dummy
 }

/*******************************************************************/

extern "C" void type_ofCall alihbtp_getnumberevents(Int_t &nev)
 {
   //passes number of events to the fortran 
   if(AliGenHBTprocessor::GetDebug()) 
    AliInfoGeneral("alihbtp_getnumberevents",Form("(%d) ....",nev));

   AliGenHBTprocessor* gen = GetAliGenHBTprocessor();//we dont check because it is done in function
   nev = gen->GetNumberOfEvents();
    
   if(AliGenHBTprocessor::GetDebug()>5) 
    AliInfoGeneral("alihbtp_getnumberevents",Form("EXITED N Ev = %d",nev));
 }

/*******************************************************************/

extern "C" void type_ofCall  alihbtp_setactiveeventnumber(Int_t & nev)
 {
//sets active event in generator (AliGenCocktailAfterBurner)

   if(AliGenHBTprocessor::GetDebug()>5) 
    ::Info("AliGenHBTpocessor.cxx: alihbtp_setactiveeventnumber","(%d)",nev);
   if(AliGenHBTprocessor::GetDebug()>0)
    ::Info("AliGenHBTpocessor.cxx: alihbtp_setactiveeventnumber","Asked for event %d",nev-1);
    
   AliGenHBTprocessor* gen = GetAliGenHBTprocessor();//we dont check because it is done in function
   
   gen->SetActiveEventNumber(nev - 1); //fortran numerates events from 1 to N
   
   if(AliGenHBTprocessor::GetDebug()>5) 
    ::Info("AliGenHBTpocessor.cxx: alihbtp_setactiveeventnumber","EXITED returned %d",nev);
 }
/*******************************************************************/
 
extern "C" void type_ofCall  alihbtp_getnumbertracks(Int_t &ntracks)
 {
//passes number of particles in active event to the fortran  
   if(AliGenHBTprocessor::GetDebug()>5) 
    ::Info("AliGenHBTpocessor.cxx: alihbtp_getnumbertracks","(%d)",ntracks);

   AliGenHBTprocessor* gen = GetAliGenHBTprocessor();//we dont check because it is done in function

   ntracks = gen->GetNumberOfTracks();
   if(AliGenHBTprocessor::GetDebug()>5)
    ::Info("AliGenHBTpocessor.cxx: alihbtp_getnumbertracks","EXITED Ntracks = %d",ntracks); 
 }
 
/*******************************************************************/
 
extern "C" void type_ofCall  
   alihbtp_puttrack(Int_t & n,Int_t& flag, Float_t& px, 
                    Float_t& py, Float_t& pz, Int_t& geantpid)
 {
//sets new parameters (momenta) in track number n
// in the active event
// n - number of the track in active event
// flag - flag of the track
// px,py,pz - momenta
// geantpid - type of the particle - Geant Particle ID
 
   if(AliGenHBTprocessor::GetDebug()>5)
    ::Info("AliGenHBTpocessor.cxx: alihbtp_puttrack","(%d)",n);

   AliGenHBTprocessor* gen = GetAliGenHBTprocessor();//we dont check because it is done in function
   
   TParticle * track = gen->GetTrack(n-1);
   if (track == 0x0)
    {
      ::Error("AliGenHBTprocessor.cxx","Can not get track from AliGenHBTprocessor");
      return;
    }
       
   //check to be deleted 
   if (geantpid != (gen->IdFromPDG( track->GetPdgCode() )))
    {
      ::Error("AliGenHBTprocessor.cxx: alihbtp_puttrack",
              "SOMETHING IS GOING BAD:\n   GEANTPIDS ARE NOT THE SAME");
    }
   
   if(AliGenHBTprocessor::GetDebug()>0)
     if (px != track->Px()) 
       ::Info("AliGenHBTprocessor.cxx: alihbtp_puttrack",
              "Px diff. = %f", px - track->Px());
   
   if(AliGenHBTprocessor::GetDebug()>3)
     ::Info("AliGenHBTprocessor.cxx: alihbtp_puttrack",
            "track->GetPdgCode() --> %d",track->GetPdgCode());
   
   
   
   Float_t m = track->GetMass();
   Float_t e = TMath::Sqrt(m*m+px*px+py*py+pz*pz);
   track->SetMomentum(px,py,pz,e);
   
   gen->SetHbtPStatusCode(flag,n-1);
   
   if(AliGenHBTprocessor::GetDebug()>5) ::Info("AliGenHBTprocessor.cxx: alihbtp_puttrack","EXITED ");
 }

/*******************************************************************/

extern "C" void type_ofCall  
  alihbtp_gettrack(Int_t & n,Int_t & flag, Float_t & px, 
                   Float_t & py, Float_t & pz, Int_t & geantpid)
  
 {
//passes track parameters to the fortran
// n - number of the track in active event
// flag - flag of the track
// px,py,pz - momenta
// geantpid - type of the particle - Geant Particle ID
 
   if(AliGenHBTprocessor::GetDebug()>5) ::Info("AliGenHBTprocessor.cxx: alihbtp_gettrack","(%d)",n);
   AliGenHBTprocessor* g = GetAliGenHBTprocessor();

   TParticle * track = g->GetTrack(n-1);
   
   flag = g->GetHbtPStatusCode(n-1);

   px = (Float_t)track->Px();
   py = (Float_t)track->Py();
   pz = (Float_t)track->Pz();
  
   geantpid = g->IdFromPDG( track->GetPdgCode() );
  
   if(AliGenHBTprocessor::GetDebug()>5) ::Info("AliGenHBTprocessor.cxx: alihbtp_gettrack","EXITED"); 
 }

/*******************************************************************/
extern "C" Float_t type_ofCall hbtpran(Int_t &)
{
//interface to the random number generator
  return gAliRandom->Rndm();
}        

/*******************************************************************/


/*******************************************************************/
