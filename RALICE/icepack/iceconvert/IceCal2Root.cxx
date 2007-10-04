/*******************************************************************************
 * Copyright(c) 2003, IceCube Experiment at the South Pole. All rights reserved.
 *
 * Author: The IceCube RALICE-based Offline Project.
 * Contributors are mentioned in the code where appropriate.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation strictly for non-commercial purposes is hereby granted
 * without fee, provided that the above copyright notice appears in all
 * copies and that both the copyright notice and this permission notice
 * appear in the supporting documentation.
 * The authors make no claims about the suitability of this software for
 * any purpose. It is provided "as is" without express or implied warranty.
 *******************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class IceCal2Root
// Conversion of Amanda (ascii) calibration data into a AliObjMatrix objects
// containing the complete OM position, calibration, Xtalk etc... database.
//
// This facility will become obsolete soon and is only kept for
// backward compatibility and testing purposes.
// Please use the more recent IceDB2Root facility instead.
//
// Via specification of the various input files, various AliObjMatrix objects
// are created and (optionally) stored in a single output file.
// The names of the AliObjMatrix objects indicate the type of database they
// contain.
// For example :
// The name "MuDaq-OMDBASE" indicates the database info for MuDaq data,
// whereas TWRDaq calibration data is in the object named "TWRDaq-OMDBASE".
// Note that a MuDaq Amacalib ascii input file has always to be provided,
// since the geometry data for all the other databases is obtained from there.   
// In addition a PDG particle database, extended with some specific Amanda
// entries, is provided as well.
// This class is derived from AliJob providing task-based processing.
// Note that the data structures are only written out if an outputfile has
// been specified via the SetOutputFile memberfunction.
// In case no outputfile has been specified, this class provides a facility
// to investigate/use the dbase data directly in subsequent (sub)tasks.
//
// The OM database information in the AliObjMatrix has the following structure :
//
// (j,1)    : Pointer to OM with identifier "j"
// (j,k+1)  : Pointer to a TF1* being the probability function for Xtalk
//            with OM "j" as transmitter and OM "k" as receiver.
//
// The latter is only available for MuDaq databases.
//
// The geometry information is directly available from the OM pointer
// in the form of its position and data words like "ORIENT" for orientation etc...
// Just use the OM memberfunction Data() to obtain a full overview. 
//
// Note : Position coordinates are indicated in meters and times are in nanoseconds,
//        in accordance with the convention used previously for Amanda.
//
// From the OM pointer also the various (de)calibration functions for
// ADC, LE and TOT can be obtained as TF1* pointers.
// The actual values of the calibration constants are stored as parameters
// of these (de)calibration functions and can be investigated via the
// usual TF1::Print() or TF1::GetParameter() facilities.
// The last two parameters of the Xtalk probability function have no effect
// on the evaluated probability value. However, these two parameters provide
// the minimum and maximum allowed LE differences between the transmitter
// and receiver hits, respectively (as can be seen from the parameter names).
//
// The (de)calibration of signals and/or determination of the Xtalk probability
// can be performed via the standard TF1::Eval(x) functionality, where "x"
// represents the input argument of the function (e.g. an uncalibrated ADC value).
//
// In general the database is not directly accessed by the user in performing
// physics analysis, since all the necessary information is contained in the
// event data itself and available via the GetSignal() memberfunction of the hits.
// However, specific tasks like e.g. calibration, Xtalk correction,
// bad module removal, noise hit removal etc... might need explicit database access.
// So, at the end of the example below some functionality is indicated for clarity.
// The user may use exactly the same procedures to obtain explicit access to the calibration
// functions etc... from the various OMs and/or hits within the actual event data which he/she
// is analysing.
//
// The PDG particle database is a standard ROOT TDatabasePDG object
// with the following extensions :
//
// Name        PDG code
// ----        --------
// brems       10001001
// deltae      10001002
// pairprod    10001003
// nucl_int    10001004
// mu_pair     10001005
// hadrons     10001006
// fiberlaser  10002100
// n2laser     10002101
// yaglaser    10002201
// z_primary   10003000
// a_primary   10003500
// 
// Usage example :
// ---------------
//
// gSystem->Load("ralice");
// gSystem->Load("icepack");
// gSystem->Load("iceconvert");
//
// IceCal2Root q("IceCal2Root","Amacalib to IcePack data structure conversion");
//
// // The MuDaq Amacalib input filename
// q.SetAmacalibFile("amacalib_amanda2_2005.txt");
//
// // The TWRDaq input filename with the newly determined TWR T0 and ADC calibs.
// q.SetTWRDaqFile("twr-cal-2005.txt");
//
// // Output file for the event structures
// q.SetOutputFile("cal2005.root");
//
// ///////////////////////////////////////////////////////////////////
// // Here the user can specify his/her sub-tasks to be executed
// // after the database structures have been filled and before the
// // data is written out.
// // Sub-tasks (i.e. a user classes derived from TTask) are entered
// // as follows :
// //
// //    MyTask1 task1("task1","Some modifications to specific OMs");
// //    MyTask2 task2("task2","Removal of specific OMs");
// //    MyTask3 task3("task3","Add private objects to the output file");
// //    q.Add(&task1);
// //    q.Add(&task2);
// //    q.Add(&task3);
// //
// // The sub-tasks will be executed in the order as they are entered.
// ///////////////////////////////////////////////////////////////////
//
// // Perform the conversion and execute subtasks (if any)
// q.ExecuteJob();
//
// // Outline of dbase usage for (de)calibration and Xtalk
//
// AliObjMatrix* omdb=q.GetOMdbase("MuDaq");
// IceAOM* om=(IceAOM*)omdb->GetObject(9,1); // Pointer to OM 9
// om->Data(); // Overview of generic module parameters
// TF1* fcal=0;   // Calibration function
// TF1* fdecal=0; // De-calibration function
// fcal=om->GetCalFunction("ADC");
// Float_t adc=248; // Uncalibrated ADC
// Float_t cadc=0;  // Calibrated ADC
// if (fcal) cadc=fcal->Eval(adc);
// fcal=om->GetCalFunction("TOT");
// Float_t tot=1538; // Uncalibrated TOT
// Float_t ctot=0;   // Calibrated TOT
// if (fcal) ctot=fcal->Eval(tot);
// fdecal=om->GetDecalFunction("LE");
// Float_t le=21697; // Uncalibrated LE
// Float_t cle=0;    // Calibrated LE
// if (fcal) cle=fcal->Eval(le);
//
// // Xtalk probability between (trans) OM 90 and (rec) OM 113
// // for a transmitter signal of uncalibrated amplitude "adc".
// TF1* fxtalkp=(TF1*)omdb->GetObject(90,113+1);
// Float_t prob=0;
// adc=378;
// if (fxtalkp) prob=fxtalkp->Eval(adc);
//
//--- Author: Nick van Eijndhoven 09-aug-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceCal2Root.h"
#include "Rstrstream.h"

ClassImp(IceCal2Root) // Class implementation to enable ROOT I/O

IceCal2Root::IceCal2Root(const char* name,const char* title) : AliJob(name,title)
{
// Default constructor.
 fAmacalFileName="";
 fTWRDaqFileName="";
 fRootFileName="";
 fOutfile=0;

 fPdg=0;
 fMuDaqdb=0;
 fTWRDaqdb=0;
}
///////////////////////////////////////////////////////////////////////////
IceCal2Root::~IceCal2Root()
{
// Default destructor.

 if (fPdg)
 {
  delete fPdg;
  fPdg=0;
 }

 if (fMuDaqdb)
 {
  delete fMuDaqdb;
  fMuDaqdb=0;
 }

 if (fTWRDaqdb)
 {
  delete fTWRDaqdb;
  fTWRDaqdb=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::SetAmacalibFile(TString name)
{
// Set the name of the Amacalib MuDaq input file.
 fAmacalFileName=name;
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::SetTWRDaqFile(TString name)
{
// Set the name of the TWRDaq calibration input file.
 fTWRDaqFileName=name;
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::SetOutputFile(TString name)
{
// Set the name of the ROOT output file.
 fRootFileName=name;
}
///////////////////////////////////////////////////////////////////////////
TDatabasePDG* IceCal2Root::GetPDG()
{
// Provide pointer to the PDG database
 return fPdg;
}
///////////////////////////////////////////////////////////////////////////
AliObjMatrix* IceCal2Root::GetOMdbase(TString name)
{
// Provide pointer to the requested OM geometry, calib. etc... database.
// Options for the "name" specification are : MuDaq, TWRDaq.
// For backward compatibility the default is name="MuDaq".

 if (name=="MuDaq") return fMuDaqdb;
 if (name=="TWRDaq") return fTWRDaqdb;
 return 0;
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::Exec(Option_t* opt)
{
// Job to convert the (ascii) database info into the IcePack structure.
//
// Notes :
// -------
// 1) This class is derived from AliJob, allowing a task based processing.
//    After conversion of the (ascii) dbase data into the IcePack structure,
//    the processing of all available sub-tasks (if any) is invoked.
//    This provides a facility to investigate/use the dbase data in
//    subsequent (sub)tasks processing before the final data structures
//    are written out.
//
// 2) Creation of a TFolder via the argument of the ExecuteJob statement
//    makes all created database objects accessible to subsequent tasks
//    via the TFolder::FindObject facility.

 if (fOutfile)
 {
  delete fOutfile;
  fOutfile=0;
 }
 if (fRootFileName != "")
 {
  fOutfile=new TFile(fRootFileName.Data(),"RECREATE","Calibration data in IcePack structure");
 }

 // Create the particle database and extend it with some F2000 specific definitions
 if (fPdg) delete fPdg;
 fPdg=new TDatabasePDG();
 fPdg->SetNameTitle("PDG-DBASE","The extended PDG particle database");
 Double_t me=fPdg->GetParticle(11)->Mass();
 fPdg->AddParticle("brems"   ,"brems"   ,0,1,0,0,"none",10001001,0,0);
 fPdg->AddParticle("deltae"  ,"deltae"  ,me,1,0,-3,"Lepton",10001002,0,0);
 fPdg->AddParticle("pairprod","pairprod",0,1,0,0,"none",10001003,0,0);
 fPdg->AddParticle("nucl_int","nucl_Int",0,1,0,0,"none",10001004,0,0);
 fPdg->AddParticle("mu_pair" ,"mu_pair" ,0,1,0,0,"none",10001005,0,0);
 fPdg->AddParticle("hadrons" ,"hadrons" ,0,1,0,0,"none",10001006,0,0);
 fPdg->AddParticle("fiberlaser","fiberlaser",0,1,0,0,"none",10002100,0,0);
 fPdg->AddParticle("n2laser"   ,"n2laser"   ,0,1,0,0,"none",10002101,0,0);
 fPdg->AddParticle("yaglaser"  ,"yaglaser"  ,0,1,0,0,"none",10002201,0,0);
 fPdg->AddParticle("z_primary","z_primary",0,1,0,0,"none",10003000,0,0);
 fPdg->AddParticle("a_primary","a_primary",0,1,0,0,"none",10003500,0,0);

 // Initialise the job working environment
 AddObject(fPdg);
 if (fOutfile) AddObject(fOutfile);

 cout << " ***" << endl;
 cout << " *** Start processing of job " << GetName() << " ***" << endl;
 cout << " ***" << endl;
 cout << " Amacalib MuDaq input file : " << fAmacalFileName.Data() << endl;
 cout << " TWRDaq input file : " << fTWRDaqFileName.Data() << endl;
 if (fOutfile) cout << " ROOT output file : " << fOutfile->GetName() << endl;

 GetMuDaqData();
 if (fMuDaqdb) AddObject(fMuDaqdb);

 GetTWRDaqData();
 if (fTWRDaqdb) AddObject(fTWRDaqdb);

 ListEnvironment();

 // Invoke all available sub-tasks (if any)
 CleanTasks();
 ExecuteTasks(opt);

 // Write the datastructures to the output file
 if (fOutfile)
 {
  fOutfile->cd();
  if (fMuDaqdb) fMuDaqdb->Write();
  if (fTWRDaqdb) fTWRDaqdb->Write();
  if (fPdg) fPdg->Write();
 }

 // Flush remaining memory resident data to the output file
 if (fOutfile) fOutfile->Write();
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::GetMuDaqData()
{
// Obtain all the MuDaq geometry, calibration and Xtalk data.

 if (fAmacalFileName=="")
 {
  cout << " *IceCal2Root GetMuDaqData* No amacalib input file specified." << endl;
  return;
 }

 fInput.clear();
 fInput.open(fAmacalFileName.Data());

 if (!fInput.good())
 {
  cout << " *IceCal2Root GetMuDaqData* Bad input file : " << fAmacalFileName.Data() << endl;
  return;
 }

 // The MuDaq OM database object
 if (fMuDaqdb)
 {
  fMuDaqdb->Reset();
 }
 else
 {
  fMuDaqdb=new AliObjMatrix();
  fMuDaqdb->SetNameTitle("MuDaq-OMDBASE","The MuDaq OM geometry, calib. etc... database");
  fMuDaqdb->SetOwner();
 }

 // Prescription of the various (de)calibration functions
 TF1 fadccal("fadccal","(x-[1])*[0]");
 TF1 fadcdecal("fadcdecal","(x/[0])+[1]");
 fadccal.SetParName(0,"BETA-ADC");
 fadccal.SetParName(1,"PED-ADC");
 fadcdecal.SetParName(0,"BETA-ADC");
 fadcdecal.SetParName(1,"PED-ADC");

 TF1 ftdccal("ftdccal","(x*[0])-[1]-([0]-1.)*32767.-[2]/sqrt([3])");
 TF1 ftdcdecal("ftdcdecal","(x+([0]-1.)*32767.+[1]+[2]/sqrt([3]))/[0]");
 ftdccal.SetParName(0,"BETA-TDC");
 ftdccal.SetParName(1,"T0");
 ftdccal.SetParName(2,"ALPHA-TDC");
 ftdccal.SetParName(3,"ADC-SLEW");
 ftdcdecal.SetParName(0,"BETA-TDC");
 ftdcdecal.SetParName(1,"T0");
 ftdcdecal.SetParName(2,"ALPHA-TDC");
 ftdcdecal.SetParName(3,"ADC-SLEW");

 TF1 ftotcal("ftotcal","x*[0]");
 TF1 ftotdecal("ftotdecal","x/[0]");
 ftotcal.SetParName(0,"BETA-TOT");
 ftotdecal.SetParName(0,"BETA-TOT");

 // The cross talk probability function
 TF1 fxtalkp("fxtalkp","(1.+[2]-[2]+[3]-[3])/(1.+exp(([0]-x)/[1]))");
 fxtalkp.SetParName(0,"C");
 fxtalkp.SetParName(1,"B");
 fxtalkp.SetParName(2,"dLE-min");
 fxtalkp.SetParName(3,"dLE-max");

 // The basic OM contents
 IceAOM om;

 // Slots to hold the various (de)calibration functions
 om.AddNamedSlot("ADC");
 om.AddNamedSlot("LE");
 om.AddNamedSlot("TOT");
 // Slots with hardware parameters
 om.AddNamedSlot("TYPE");
 om.AddNamedSlot("ORIENT");
/// om.AddNamedSlot("THRESH");
/// om.AddNamedSlot("SENSIT");
/// om.AddNamedSlot("READOUT"); // 0=unknown 1=electrical 2=optical 3=digital

 fInput.seekg(0); // Position at beginning of file
 fInput >> dec;   // Make sure all integers starting with 0 are taken in decimal format

 TString s;
 Int_t jmod,type,serial,string,ix,iy,iz,ori;
 Float_t costh=0;
/// Float_t thresh=0;
/// Float_t sensit=1;
/// Float_t readout=0;
 Double_t pos[3]={0,0,0};
 Float_t ped,beta,alpha;
 Int_t pol;
 Float_t totped;
 Int_t jtrans,jrec;
 Float_t c,b,dlemin,dlemax;
 IceAOM* omx=0;
 TF1* fcal=0;
 TF1* fdecal=0;
 while (fInput >> s)
 {
  if (s == "P") // Read the Geom data
  {
   fInput >> jmod >> type >> serial >> string >> ix >> iy >> iz >> ori;
   omx=(IceAOM*)fMuDaqdb->GetObject(jmod,1);
   if (!omx)
   {
    omx=new IceAOM(om);
    omx->SetUniqueID(jmod);
    fMuDaqdb->EnterObject(jmod,1,omx);
   }
   pos[0]=double(ix)/1000.;
   pos[1]=double(iy)/1000.;
   pos[2]=double(iz)/1000.;
   omx->SetPosition(pos,"car");
   costh=1;
   if (ori==2) costh=-1;
   omx->SetSignal(type,"TYPE");
   omx->SetSignal(costh,"ORIENT");
///   omx->SetSignal(thresh,"THRESH");
///   omx->SetSignal(sensit,"SENSIT");
///   omx->SetSignal(readout,"READOUT");
  }
  else if (s == "T") // Read the Time calibration constants
  {
   fInput >> jmod >> ped >> beta >> alpha >> pol;
   omx=(IceAOM*)fMuDaqdb->GetObject(jmod,1);
   if (!omx)
   {
    omx=new IceAOM(om);
    omx->SetUniqueID(jmod);
    fMuDaqdb->EnterObject(jmod,1,omx);
   }

   omx->SetCalFunction(&ftdccal,"LE");
   omx->SetDecalFunction(&ftdcdecal,"LE");
   omx->SetCalFunction(&ftotcal,"TOT");
   omx->SetDecalFunction(&ftotdecal,"TOT");

   // Flag time slots of bad OMs as dead and don't provide time (de)calib functions
   if (ped<-1e5 || beta<=0 || alpha<0)
   {
    omx->SetDead("LE");
    omx->SetDead("TOT");
    omx->SetCalFunction(0,"LE");
    omx->SetDecalFunction(0,"LE");
    omx->SetCalFunction(0,"TOT");
    omx->SetDecalFunction(0,"TOT");
   }

   fcal=omx->GetCalFunction("LE");
   fdecal=omx->GetDecalFunction("LE");
   if (fcal)
   {
    fcal->SetParameter(0,beta);
    fcal->SetParameter(1,ped);
    fcal->SetParameter(2,alpha);
    fcal->SetParameter(3,1.e20);
   }
   if (fdecal)
   {
    fdecal->SetParameter(0,beta);
    if (!beta) fdecal->SetParameter(0,1);
    fdecal->SetParameter(1,ped);
    fdecal->SetParameter(2,alpha);
    fdecal->SetParameter(3,1.e20);
   }

   fcal=omx->GetCalFunction("TOT");
   fdecal=omx->GetDecalFunction("TOT");
   if (fcal)
   {
    fcal->SetParameter(0,beta);
   }
   if (fdecal)
   {
    fdecal->SetParameter(0,beta);
   }
  }
  else if (s == "A") // Read the Amplitude calibration constants
  {
   fInput >> jmod >> ped >> beta >> totped >> pol;
   omx=(IceAOM*)fMuDaqdb->GetObject(jmod,1);
   if (!omx)
   {
    omx=new IceAOM(om);
    omx->SetUniqueID(jmod);
    fMuDaqdb->EnterObject(jmod,1,omx);
   }

   omx->SetCalFunction(&fadccal,"ADC");
   omx->SetDecalFunction(&fadcdecal,"ADC");

   // Flag amplitude slots of bad OMs as dead and don't provide amplitude (de)calib functions
   if (ped<-1e5 || beta<=0)
   {
    omx->SetDead("ADC");
    omx->SetCalFunction(0,"ADC");
    omx->SetDecalFunction(0,"ADC");
   }
   if (totped<-1e5)
   {
    omx->SetDead("TOT");
    omx->SetCalFunction(0,"TOT");
    omx->SetDecalFunction(0,"TOT");
   }

   fcal=omx->GetCalFunction("ADC");
   fdecal=omx->GetDecalFunction("ADC");
   if (fcal)
   {
    fcal->SetParameter(0,beta);
    fcal->SetParameter(1,ped);
   }
   if (fdecal)
   {
    fdecal->SetParameter(0,beta);
    if (!beta) fdecal->SetParameter(0,1);
    fdecal->SetParameter(1,ped);
   }
  }
  else if (s == "K") // Read the cross talk probability constants
  {
   fInput >> jtrans >> jrec >> c >> b >> dlemin >> dlemax;
   omx=(IceAOM*)fMuDaqdb->GetObject(jtrans,1);
   if (!omx)
   {
    omx=new IceAOM(om);
    omx->SetUniqueID(jtrans);
    fMuDaqdb->EnterObject(jtrans,1,omx);
   }

   TF1* fx=new TF1(fxtalkp);
   fx->SetParameter(0,c);
   if (b)
   {
    fx->SetParameter(1,b);
   }
   else
   {
    fx->SetParameter(1,1);
   }
   fx->SetParameter(2,dlemin);
   fx->SetParameter(3,dlemax);
   fMuDaqdb->EnterObject(jtrans,jrec+1,fx);
  }
  else // Skip this line
  {
   fInput.ignore(99999,'\n');
  } 
 }

 fInput.close();
}
///////////////////////////////////////////////////////////////////////////
void IceCal2Root::GetTWRDaqData()
{
// Obtain all the TWRDaq geometry and calibration data.

 if (fTWRDaqFileName=="")
 {
  cout << " *IceCal2Root GetTWRDaqData* No TWRDaq calibration data will be produced." << endl;
  return;
 }

 fInput.clear();
 fInput.open(fTWRDaqFileName.Data());

 if (!fInput.good())
 {
  cout << " *IceCal2Root GetTWRDaqData* Bad input file : " << fTWRDaqFileName.Data() << endl;
  return;
 }

 // The geometry info will be obtained from the MuDaq OM database
 if (!fMuDaqdb)
 {
  cout << " *IceCal2Root GetTWRDaqData* MuDaq OM geometry database is missing." << endl;
  return;
 }

 // The TWRDaq OM database object
 if (fTWRDaqdb)
 {
  fTWRDaqdb->Reset();
 }
 else
 {
  fTWRDaqdb=new AliObjMatrix();
  fTWRDaqdb->SetNameTitle("TWRDaq-OMDBASE","The TWRDaq OM geometry, calib. etc... database");
  fTWRDaqdb->SetOwner();
 }

 // Prescription of the various (de)calibration functions
 TF1 fadccal("fadccal","x*(5./4096.)/(50.*[0])");
 TF1 fadcdecal("fadcdecal","x*(50.*[0])/(5./4096.)");
 fadccal.SetParName(0,"nC/PE");
 fadcdecal.SetParName(0,"nC/PE");

 TF1 ftdccal("ftdccal","x-[0]");
 TF1 ftdcdecal("ftdcdecal","x+[0]");
 ftdccal.SetParName(0,"T0");
 ftdcdecal.SetParName(0,"T0");

 TF1 ftotcal("ftotcal","x*[0]");
 TF1 ftotdecal("ftotdecal","x/[0]");
 ftotcal.SetParName(0,"TOT-FACT");
 ftotdecal.SetParName(0,"TOT-FACT");

 fInput.seekg(0); // Position at beginning of file
 fInput >> dec;   // Make sure all integers starting with 0 are taken in decimal format

 Int_t jmod;
 Float_t t0;
 Float_t ncpe;
 IceAOM* omx=0;
 TF1* fcal=0;
 TF1* fdecal=0;
 while (fInput >> jmod >> t0 >> ncpe)
 {
  // Copy the Geom data from the MuDaq OM database
  IceAOM* omg=(IceAOM*)fMuDaqdb->GetObject(jmod,1);
  if (!omg) continue;

  omx=new IceAOM(*omg);

  // Reset all previous "dead" flags
  omx->SetAlive("ADC");
  omx->SetAlive("LE");
  omx->SetAlive("TOT");

  // Enter the TWRDaq (de)calibration functions
  omx->SetCalFunction(&fadccal,"ADC");
  omx->SetDecalFunction(&fadcdecal,"ADC");
  omx->SetCalFunction(&ftdccal,"LE");
  omx->SetDecalFunction(&ftdcdecal,"LE");
  omx->SetCalFunction(&ftotcal,"TOT");
  omx->SetDecalFunction(&ftotdecal,"TOT");

  // Flag slots of bad OMs as dead and don't provide time (de)calib functions
  if (ncpe<=0)
  {
   omx->SetDead("ADC");
   omx->SetCalFunction(0,"ADC");
   omx->SetDecalFunction(0,"ADC");
  }
  if (t0<-999)
  {
   omx->SetDead("LE");
   omx->SetDead("TOT");
   omx->SetCalFunction(0,"LE");
   omx->SetDecalFunction(0,"LE");
   omx->SetCalFunction(0,"TOT");
   omx->SetDecalFunction(0,"TOT");
  }

  // Set the TWRDaq (de)calibration function parameters for the good OMs
  fcal=omx->GetCalFunction("ADC");
  fdecal=omx->GetDecalFunction("ADC");
  if (fcal)
  {
   fcal->SetParameter(0,ncpe);
  }
  if (fdecal)
  {
   fdecal->SetParameter(0,ncpe);
  }

  fcal=omx->GetCalFunction("LE");
  fdecal=omx->GetDecalFunction("LE");
  if (fcal)
  {
   fcal->SetParameter(0,t0);
  }
  if (fdecal)
  {
   fdecal->SetParameter(0,t0);
  }

  fcal=omx->GetCalFunction("TOT");
  fdecal=omx->GetDecalFunction("TOT");
  if (fcal)
  {
   fcal->SetParameter(0,1);
  }
  if (fdecal)
  {
   fdecal->SetParameter(0,1);
  }

  fTWRDaqdb->EnterObject(jmod,1,omx);

 } // End of reading loop

 fInput.close();
}
///////////////////////////////////////////////////////////////////////////
