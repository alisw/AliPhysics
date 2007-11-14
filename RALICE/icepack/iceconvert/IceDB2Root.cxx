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
// Class IceDB2Root
// Extraction of Amanda and IceCube calibration data from the IceCube database.
// Calibration data are stored into an AliObjMatrix containing the complete OM 
// position, calibration, Xtalk etc... database.
// Various AliObjMatrix objects are created and (optionally) stored in a single 
// output file.
// The names of the AliObjMatrix objects indicate the type of database they
// contain.
// For example :
// The name "MuDaq-OMDBASE" indicates the database info for MuDaq data,
// whereas TWRDaq calibration data is in the object named "TWRDaq-OMDBASE".
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
// The latter is only available for MuDaq and TWRDaq databases.
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
// IceDB2Root q("IceDB2Root","DB to IcePack data structure conversion");
//
// // Output file for the event structures
// q.SetOutputFile("cal2005.root");
// q.SetUT(2005,7,1);
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
//--- Author: Garmt de Vries-Uiterweerd 06-jun-2007 Utrecht University
//- Modified: GdVU $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceDB2Root.h"
#include "Riostream.h"

ClassImp(IceDB2Root) // Class implementation to enable ROOT I/O

IceDB2Root::IceDB2Root(const char* name,const char* title) : AliJob(name,title)
{
// Default constructor.
 fDBName="mysql://icedb.umh.ac.be/I3OmDb";
 fUser="www";
 fPassword="";
 fRootFileName="";
 fOutfile=0;

 fPdg=0;
 fMuDaqdb=0;
 fTWRDaqdb=0;
 fJEBTDaqdb=0;
 fJEBADaqdb=0;
}
///////////////////////////////////////////////////////////////////////////
IceDB2Root::~IceDB2Root()
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

 if (fJEBTDaqdb)
 {
  delete fJEBTDaqdb;
  fJEBTDaqdb=0;
 }

 if (fJEBADaqdb)
 {
  delete fJEBADaqdb;
  fJEBADaqdb=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::SetDatabase(TString name, TString user, TString password)
{
// Set the name of the database, and username/password needed to access it.
 fDBName=name;
 fUser=user;
 fPassword=password;
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::SetOutputFile(TString name)
{
// Set the name of the ROOT output file.
 fRootFileName=name;
}
///////////////////////////////////////////////////////////////////////////
AliTimestamp IceDB2Root::GetTime()
{
// Return time for which the calibration is done.
 return fTime;
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::SetUT(Int_t y, Int_t m, Int_t d, Int_t hh, Int_t mm, Int_t ss)
{
// Set time for which calibration is done. Input parameters are the same as 
// for AliTimestamp::SetUT. Default value for hh, mm and ss is 0.
 fTime.SetUT(y,m,d,hh,mm,ss);
}
///////////////////////////////////////////////////////////////////////////
TDatabasePDG* IceDB2Root::GetPDG()
{
// Provide pointer to the PDG database
 return fPdg;
}
///////////////////////////////////////////////////////////////////////////
AliObjMatrix* IceDB2Root::GetOMdbase(TString name)
{
// Provide pointer to the requested OM geometry, calib. etc... database.
// Options for the "name" specification are : MuDaq, TWRDaq, JEBTDaq, JEBADaq.
// For backward compatibility the default is name="MuDaq".

 if (name=="MuDaq") return fMuDaqdb;
 if (name=="TWRDaq") return fTWRDaqdb;
 if (name=="JEBTDaq") return fJEBTDaqdb;
 if (name=="JEBADaq") return fJEBADaqdb;
 return 0;
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::Exec(Option_t* opt)
{
// Job to extract calibration data from database into the IcePack structure.
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
 cout << " Database : " << fDBName.Data() << endl;
 if (fOutfile) cout << " ROOT output file : " << fOutfile->GetName() << endl;

 // Get all the data
 GetMuDaqData();
 if (fMuDaqdb) AddObject(fMuDaqdb);

 GetTWRDaqData();
 if (fTWRDaqdb) AddObject(fTWRDaqdb);

 GetJEBTDaqData();
 if (fJEBTDaqdb) AddObject(fJEBTDaqdb);

 GetJEBADaqData();
 if (fJEBADaqdb) AddObject(fJEBADaqdb);

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
  if (fJEBTDaqdb) fJEBTDaqdb->Write();
  if (fJEBADaqdb) fJEBADaqdb->Write();
  if (fPdg) fPdg->Write();
 }

 // Flush remaining memory resident data to the output file
 if (fOutfile) fOutfile->Write();
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::GetMuDaqData()
{
// Obtain all the MuDaq geometry, calibration and Xtalk data.

 // Connect to the DB server
 TSQLServer* server=TSQLServer::Connect(fDBName.Data(),fUser.Data(),fPassword.Data()); 
 if (!server)
 {
  cout << " *IceDB2Root GetMuDaqData* Database " << fDBName.Data() << " could not be accessed." << endl;
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
 om.AddNamedSlot("TYPE"); // -1=unknown 0=std_coax 1=std_twisted/std_hybrid 2=desy_hybrid 3=uci_fiber 4=lbl_dom 5=dAOM_ld 6=dAOM_led 9=desy_active_fiber
 om.AddNamedSlot("ORIENT");
 om.AddNamedSlot("THRESH");
 om.AddNamedSlot("SENSIT");
 om.AddNamedSlot("READOUT"); // 0=unknown 1=electrical 2=optical 3=digital

 // Default values
 Float_t type=-1;
 om.SetSignal(type,"TYPE");
 Float_t costh=0;
 om.SetSignal(costh,"ORIENT");
 Float_t thresh=0;
 om.SetSignal(thresh,"THRES");
 Float_t sensit=1;
 om.SetSignal(sensit,"SENSIT");
 Float_t readout=0;
 om.SetSignal(readout,"READOUT");

 // Variables needed
 TSQLStatement* st;
 IceAOM* omx=0;
 Int_t omid;
 Double_t pos[3]={0,0,0};
 Float_t ped, beta, alpha;
 Int_t jtrans,jrec;
 Float_t c, b, dlemin, dlemax;
 TF1* fcal=0;
 TF1* fdecal=0;
 AliTimestamp validitystart, validityend;
 Int_t revision[682];

 // Difference between Amanda and IceCube coordinates
 Float_t dx=-339.8;
 Float_t dy=-117.4;
 Float_t dz=-224.6;

 // Get positions of Amanda OMs in Amanda coordinates
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, X, Y, Z, AmandaOm.TypeId, Orientation FROM GeometryOm INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE GeometryOm.StringId=AmandaOm.StringId AND GeometryOm.TubeId=AmandaOm.TubeId AND GeometryOm.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=2;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Positions could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    pos[0]=st->GetDouble(4)+dx;
    pos[1]=st->GetDouble(5)+dy;
    pos[2]=st->GetDouble(6)+dz;
    omx->SetPosition(pos,"car");
    type=(Float_t)st->GetInt(7);
    omx->SetSignal(type,"TYPE");
    costh=(Float_t)st->GetInt(8);
    omx->SetSignal(costh,"ORIENT");
   }
  }
 }

 // Get sensitivity/threshold of Amanda OMs
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id, sc_sensit, sc_thresh FROM AmandaStatusOm INNER JOIN CalibrationDetail WHERE AmandaStatusOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=804;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Sensitivity/threshold could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    thresh=(Float_t)st->GetDouble(4);
    sensit=(Float_t)st->GetDouble(5);
    omx->SetSignal(thresh,"THRES");
    omx->SetSignal(sensit,"SENSIT");
   }
  }
 }

 // Get readout types of Amanda OMs
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, CableType FROM AmandaTWRCableType INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE AmandaTWRCableType.StringId=AmandaOm.StringId AND AmandaTWRCableType.TubeId=AmandaOm.TubeId AND AmandaTWRCableType.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=853;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Readout types could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    if(st->GetUInt(4) == 0) readout=1;        // Electrical
    else if(st->GetUInt(4) == 10) readout=2;  // Optical
    else readout=0;                           // Unknown
    omx->SetSignal(readout,"READOUT");
   }
  }
 }

 // Get muDAQ amplitude calibration constants
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id, ac_pedestal, ac_spepeak FROM AmandaAmplOm INNER JOIN CalibrationDetail WHERE AmandaAmplOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=802;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Amplitude calibrations could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Set calibration functions
    omx->SetCalFunction(&fadccal,1);
    omx->SetDecalFunction(&fadcdecal,1);
    // Get calibration values
    ped=(Float_t)st->GetDouble(4);
    if (st->GetDouble(5)>0) beta=(Float_t)(1/st->GetDouble(5));
    else beta=-999999;
    // Flag amplitude slots of bad OMs as dead and don't provide amplitude (de)calib functions
    if (ped<-1e5 || beta<=0)
    {
     omx->SetDead(1);
     omx->SetCalFunction(0,1);
     omx->SetDecalFunction(0,1);
    }
    // Set calibration function parameters
    fcal=omx->GetCalFunction(1);
    fdecal=omx->GetDecalFunction(1);
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
  }
 }

 // Get muDAQ time calibration constants
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id, tc_tzero, tc_beta_t, tc_alpha_t FROM AmandaTimeOm INNER JOIN CalibrationDetail WHERE AmandaTimeOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=805;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Time calibrations could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Set calibration functions
    omx->SetCalFunction(&ftdccal,2);
    omx->SetDecalFunction(&ftdcdecal,2);
    omx->SetCalFunction(&ftotcal,3);
    omx->SetDecalFunction(&ftotdecal,3);
    // Get calibration values
    ped=(Float_t)st->GetDouble(4);
    beta=(Float_t)st->GetDouble(5);
    alpha=(Float_t)st->GetDouble(6);
    // Flag time slots of bad OMs as dead and don't provide time (de)calib functions
    if (ped<-1e5 || beta<=0 || alpha<0)
    {
     omx->SetDead(2);
     omx->SetDead(3);
     omx->SetCalFunction(0,2);
     omx->SetDecalFunction(0,2);
     omx->SetCalFunction(0,3);
     omx->SetDecalFunction(0,3);
    }
    // Set calibration function parameters
    fcal=omx->GetCalFunction(2);
    fdecal=omx->GetDecalFunction(2);
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
    fcal=omx->GetCalFunction(3);
    fdecal=omx->GetDecalFunction(3);
    if (fcal)
    {
     fcal->SetParameter(0,beta);
    }
    if (fdecal)
    {
     fdecal->SetParameter(0,beta);
    }
   }
  }
 }

 // Get Xtalk probability constants
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_talker_id, om_receiver_id, threshold, width, timelow, timehigh FROM AmandaXtalkOm INNER JOIN CalibrationDetail WHERE AmandaXtalkOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=807;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Xtalk probability constants could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   jtrans=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[jtrans])
   {
    revision[jtrans]=st->GetInt(2);
    jrec=st->GetInt(4);
    // Get transmitter OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(jtrans,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(jtrans);
     fMuDaqdb->EnterObject(jtrans,1,omx);
    }
    // Get calibration values
    c=(Float_t)st->GetDouble(5);
    b=(Float_t)st->GetDouble(6);
    dlemin=(Float_t)st->GetDouble(7);
    dlemax=(Float_t)st->GetDouble(8);
    // Make Xtalk probability function and set parameters
    TF1* fx=new TF1(fxtalkp);
    fx->SetParameter(0,c);
    if (b) fx->SetParameter(1,b);
    else fx->SetParameter(1,1);
    fx->SetParameter(2,dlemin);
    fx->SetParameter(3,dlemax);
    fMuDaqdb->EnterObject(jtrans,jrec+1,fx);
   }
  }
 }

 // Flag OMs in bad OM list as dead
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id FROM AmandaBadOm INNER JOIN CalibrationDetail WHERE AmandaBadOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=803;");
 if (!st)
 {
  cout << " *IceDB2Root GetMuDaqData* Bad OM list could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fMuDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fMuDaqdb->EnterObject(omid,1,omx);
    }
    // Flag OMs as dead
    omx->SetDead(1);
    omx->SetCalFunction(0,1);
    omx->SetDecalFunction(0,1);
    omx->SetDead(2);
    omx->SetCalFunction(0,2);
    omx->SetDecalFunction(0,2);
    omx->SetDead(3);
    omx->SetCalFunction(0,3);
    omx->SetDecalFunction(0,3);
   }
  }
 }

}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::GetTWRDaqData()
{
// Obtain all the TWRDaq geometry, calibration and Xtalk data.

 // Connect to the DB server
 TSQLServer* server=TSQLServer::Connect(fDBName.Data(),fUser.Data(),fPassword.Data()); 
 if (!server)
 {
  cout << " *IceDB2Root GetTWRDaqData* Database " << fDBName.Data() << " could not be accessed." << endl;
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

 // Conversion of adc to charge in nC
 // Volt per digit = 5/4096    Assumed capacitance = 20 pF
 // Charge (in nC) is adc*(5./4096.)*0.02
 // The database calibration constant is the amount of nC per photon electron (nC/PE)
 TF1 fadccal("fadccal","x*(5./4096.)/(50.*[0])");
 TF1 fadcdecal("fadcdecal","x*(50.*[0])/(5./4096.)");
 fadccal.SetParName(0,"nC/PE");
 fadcdecal.SetParName(0,"nC/PE");

 TF1 ftdccal("ftdccal","x-[0]");
 TF1 ftdcdecal("ftdcdecal","x+[0]");
 ftdccal.SetParName(0,"T0");
 ftdcdecal.SetParName(0,"T0");

 TF1 ftotcal("ftotcal","x");
 TF1 ftotdecal("ftotdecal","x");

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
 om.AddNamedSlot("TYPE"); // -1=unknown 0=std_coax 1=std_twisted/std_hybrid 2=desy_hybrid 3=uci_fiber 4=lbl_dom 5=dAOM_ld 6=dAOM_led 9=desy_active_fiber
 om.AddNamedSlot("ORIENT");
 om.AddNamedSlot("THRESH");
 om.AddNamedSlot("SENSIT");
 om.AddNamedSlot("READOUT"); // 0=unknown 1=electrical 2=optical 3=digital
 om.AddNamedSlot("BINSIZE");
 om.AddNamedSlot("EXTSTOP");

 // Default values
 Float_t type=-1;
 om.SetSignal(type,"TYPE");
 Float_t costh=0;
 om.SetSignal(costh,"ORIENT");
 Float_t thresh=0;
 om.SetSignal(thresh,"THRESH");
 Float_t sensit=1;
 om.SetSignal(sensit,"SENSIT");
 Float_t readout=0;
 om.SetSignal(readout,"READOUT");
 Float_t binsize=0;
 om.SetSignal(binsize,"BINSIZE");
 Float_t stopdelay=0;
 om.SetSignal(stopdelay,"EXTSTOP");

 // Variables needed
 TSQLStatement* st;
 IceAOM* omx=0;
 Int_t omid;
 Double_t pos[3]={0,0,0};
 Float_t peArea, twrT0;
 Int_t jtrans,jrec;
 Float_t c, b, dlemin, dlemax;
 TF1* fcal=0;
 TF1* fdecal=0;
 AliTimestamp validitystart, validityend;
 Int_t revision[682];

 // Get positions of Amanda OMs in IceCube coordinates
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, X, Y, Z, AmandaOm.TypeId, Orientation FROM GeometryOm INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE GeometryOm.StringId=AmandaOm.StringId AND GeometryOm.TubeId=AmandaOm.TubeId AND GeometryOm.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=2;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Positions could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    pos[0]=st->GetDouble(4);
    pos[1]=st->GetDouble(5);
    pos[2]=st->GetDouble(6);
    omx->SetPosition(pos,"car");
    type=(Float_t)st->GetInt(7);
    omx->SetSignal(type,"TYPE");
    costh=(Float_t)st->GetInt(8);
    omx->SetSignal(costh,"ORIENT");
   }
  }
 }

 // Get sensitivity/threshold of Amanda OMs
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id, sc_sensit, sc_thresh FROM AmandaStatusOm INNER JOIN CalibrationDetail WHERE AmandaStatusOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=804;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Sensitivity/threshold could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    thresh=(Float_t)st->GetDouble(4);
    sensit=(Float_t)st->GetDouble(5);
    omx->SetSignal(thresh,"THRESH");
    omx->SetSignal(sensit,"SENSIT");
   }
  }
 }

 // Get readout types of Amanda OMs
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, CableType FROM AmandaTWRCableType INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE AmandaTWRCableType.StringId=AmandaOm.StringId AND AmandaTWRCableType.TubeId=AmandaOm.TubeId AND AmandaTWRCableType.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=853;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Readout types could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    if(st->GetUInt(4) == 0) readout=1;        // Electrical
    else if(st->GetUInt(4) == 10) readout=2;  // Optical
    else readout=0;                           // Unknown
    omx->SetSignal(readout,"READOUT");
   }
  }
 }

 // Get bin sizes and stop delays types of TWRs
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, BinSize, StopDelay FROM AmandaTWRStandard INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE AmandaTWRStandard.StringId=AmandaOm.StringId AND AmandaTWRStandard.TubeId=AmandaOm.TubeId AND AmandaTWRStandard.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=852;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Bin size and stop delays could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    binsize=(Float_t)st->GetInt(4);
    stopdelay=(Float_t)st->GetInt(5);
    omx->SetSignal(binsize,"BINSIZE");
    omx->SetSignal(stopdelay,"EXTSTOP");
   }
  }
 }

 // Get TWRDaq amplitude and time calibration constants
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, OmId, peArea, twrT0 FROM AmandaTWRCalibration INNER JOIN AmandaOm INNER JOIN CalibrationDetail WHERE AmandaTWRCalibration.StringId=AmandaOm.StringId AND AmandaTWRCalibration.TubeId=AmandaOm.TubeId AND AmandaTWRCalibration.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=854;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Bin size and stop delays could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Set calibration functions
    omx->SetCalFunction(&fadccal,1);
    omx->SetDecalFunction(&fadcdecal,1);
    omx->SetCalFunction(&ftdccal,2);
    omx->SetDecalFunction(&ftdcdecal,2);
    omx->SetCalFunction(&ftotcal,3);
    omx->SetDecalFunction(&ftotdecal,3);
    // Get calibration values
    peArea=(Float_t)st->GetDouble(4);
    twrT0=(Float_t)st->GetDouble(5);
    // Flag amplitude slots of bad OMs as dead and don't provide amplitude (de)calib functions
    if (peArea<=0 || omx->GetSignal("BINSIZE")<=0)
    {
     omx->SetDead(1);
     omx->SetCalFunction(0,1);
     omx->SetDecalFunction(0,1);
    }
    // Set amplitude calibration function parameters
    fcal=omx->GetCalFunction(1);
    fdecal=omx->GetDecalFunction(1);
    if (fcal)
    {
     // peArea as read from the DB is the factor that converts the integrated
     // area under the peak to the number of pe, that is, the sum over all bins
     // of ADC*binsize. In IcePack, we simply add up the bin contents of the
     // peak, without multiplying by binsize. Hence, the calibration factor is 
     // peArea/binsize, rather than peArea.
     fcal->SetParameter(0,peArea/omx->GetSignal("BINSIZE"));
    }
    if (fdecal)
    {
     fdecal->SetParameter(0,peArea/omx->GetSignal("BINSIZE"));
     if (!peArea) fdecal->SetParameter(0,1);
    }
    // Flag LE slots of bad OMs as dead and don't provide time (de)calib functions
    if (twrT0<=0)
    {
     omx->SetDead(2);
     omx->SetCalFunction(0,2);
     omx->SetDecalFunction(0,2);
    }
    // Set time calibration function parameters
    fcal=omx->GetCalFunction(2);
    fdecal=omx->GetDecalFunction(2);
    if (fcal)
    {
     fcal->SetParameter(0,twrT0);
    }
    if (fdecal)
    {
     fdecal->SetParameter(0,twrT0);
    }
   }
  }
 }

 // Get Xtalk probability constants
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_talker_id, om_receiver_id, threshold, width, timelow, timehigh FROM AmandaXtalkOm INNER JOIN CalibrationDetail WHERE AmandaXtalkOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=807;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Xtalk probability constants could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   jtrans=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[jtrans]=st->GetInt(2);
    jrec=st->GetInt(4);
    // Get transmitter OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(jtrans,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(jtrans);
     fTWRDaqdb->EnterObject(jtrans,1,omx);
    }
    // Get calibration values
    c=(Float_t)st->GetDouble(5);
    b=(Float_t)st->GetDouble(6);
    dlemin=(Float_t)st->GetDouble(7);
    dlemax=(Float_t)st->GetDouble(8);
    // Make Xtalk probability function and set parameters
    TF1* fx=new TF1(fxtalkp);
    fx->SetParameter(0,c);
    if (b) fx->SetParameter(1,b);
    else fx->SetParameter(1,1);
    fx->SetParameter(2,dlemin);
    fx->SetParameter(3,dlemax);
    fTWRDaqdb->EnterObject(jtrans,jrec+1,fx);
   }
  }
 }

 // Flag OMs in bad OM list as dead
 for(omid=0; omid<=681; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, om_id FROM AmandaBadOm INNER JOIN CalibrationDetail WHERE AmandaBadOm.ca_id=CalibrationDetail.CaId AND CalibrationDetail.TypeId=803;");
 if (!st)
 {
  cout << " *IceDB2Root GetTWRDaqData* Bad OM list could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=st->GetInt(3);
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceAOM*)fTWRDaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceAOM(om);
     omx->SetUniqueID(omid);
     fTWRDaqdb->EnterObject(omid,1,omx);
    }
    // Flag OMs as dead
    omx->SetDead(1);
    omx->SetCalFunction(0,1);
    omx->SetDecalFunction(0,1);
    omx->SetDead(2);
    omx->SetCalFunction(0,2);
    omx->SetDecalFunction(0,2);
    omx->SetDead(3);
    omx->SetCalFunction(0,3);
    omx->SetDecalFunction(0,3);
   }
  }
 }
 
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::GetJEBTDaqData()
{
// Obtain all the JEB TWRDaq geometry, calibration and Xtalk data.
// For the time being, the JEBTDaq database is just a copy of the
// TWRDaq database. If different calibrations are needed in the 
// future, this member function can be updated accordingly.

 // The JEBTDaq OM database object
 if (fJEBTDaqdb)
 {
  delete fJEBTDaqdb;
  fJEBTDaqdb=0;
 }

 // Copy TWRDaq database to JEBTDaq database via the Clone() function
 if(!fTWRDaqdb) GetTWRDaqData();
 if(fTWRDaqdb)
 {
  fJEBTDaqdb=(AliObjMatrix*)fTWRDaqdb->Clone();
  fJEBTDaqdb->SetNameTitle("JEBTDaq-OMDBASE","The JEBTDaq OM geometry, calib. etc... database");
 }
}
///////////////////////////////////////////////////////////////////////////
void IceDB2Root::GetJEBADaqData()
{
// Obtain all the JEB ATWDDaq geometry, calibration.

 // Connect to the DB server
 TSQLServer* server=TSQLServer::Connect(fDBName.Data(),fUser.Data(),fPassword.Data()); 
 if (!server)
 {
  cout << " *IceDB2Root GetJEBADaqData* Database " << fDBName.Data() << " could not be accessed." << endl;
  return;
 }

 // The JEBADaq OM database object
 if (fJEBADaqdb)
 {
  fJEBADaqdb->Reset();
 }
 else
 {
  fJEBADaqdb=new AliObjMatrix();
  fJEBADaqdb->SetNameTitle("JEBADaq-OMDBASE","The JEBADaq OM geometry, calib. etc... database");
  fJEBADaqdb->SetOwner();
 }

 // Prescription of the various (de)calibration functions
 TF1 fadccal("fadccal","x*[0]");
 TF1 fadcdecal("fadcdecal","x/[0]");
 fadccal.SetParName(0,"ADC-SLOPE");
 fadcdecal.SetParName(0,"ADC-SLOPE");

 TF1 ftdccal("ftdccal","x");
 TF1 ftdcdecal("ftdcdecal","x");

 TF1 ftotcal("ftotcal","x");
 TF1 ftotdecal("ftotdecal","x");

 // The basic OM contents
 IceDOM om;

 // Slots to hold the various (de)calibration functions
 om.AddNamedSlot("ADC");
 om.AddNamedSlot("LE");
 om.AddNamedSlot("TOT");
 // Slots with hardware parameters
 om.AddNamedSlot("ORIENT");

 // Default values
 Float_t costh=0;
 om.SetSignal(costh,"ORIENT");

 // Variables needed
 TSQLStatement* st;
 IceDOM* omx=0;
 Int_t omid;
 Double_t pos[3]={0,0,0};
 Float_t peArea;
 TF1* fcal=0;
 TF1* fdecal=0;
 AliTimestamp validitystart, validityend;
 const Int_t maxomid=8064;
 Int_t revision[maxomid+1];

 // Get positions of IceCube DOMs in IceCube coordinates
 for(omid=0; omid<=maxomid; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, StringId, TubeId, X, Y, Z, Orientation FROM GeometryOm INNER JOIN CalibrationDetail WHERE StringId>0 AND GeometryOm.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=2;");
 if (!st)
 {
  cout << " *IceDB2Root GetJEBADaqData* Positions could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=om.GetOMId(st->GetInt(3),st->GetInt(4));
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceDOM*)fJEBADaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceDOM(om);
     omx->SetUniqueID(omid);
     fJEBADaqdb->EnterObject(omid,1,omx);
    }
    // Enter calibration values
    pos[0]=st->GetDouble(5);
    pos[1]=st->GetDouble(6);
    pos[2]=st->GetDouble(7);
    omx->SetPosition(pos,"car");
    costh=(Float_t)st->GetInt(8);
    omx->SetSignal(costh,"ORIENT");
   }
  }
 }

 // Get JEBADaq amplitude calibration constants
 //// TODO: Insert correct table and field names and use correct calibration type
 for(omid=0; omid<=maxomid; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, StringId, TubeId, XXXpeArea FROM XXXampl INNER JOIN CalibrationDetail WHERE StringId>0 AND XXXampl.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=XXX;");
 if (!st)
 {
  cout << " *IceDB2Root GetJEBADaqData* Amplitude calibrations could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=om.GetOMId(st->GetInt(3),st->GetInt(4));
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceDOM*)fJEBADaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceDOM(om);
     omx->SetUniqueID(omid);
     fJEBADaqdb->EnterObject(omid,1,omx);
    }
    // Set calibration functions
    omx->SetCalFunction(&fadccal,1);
    omx->SetDecalFunction(&fadcdecal,1);
    // Get calibration values
    peArea=(Float_t)st->GetDouble(5);
    // Flag amplitude slots of bad OMs as dead and don't provide amplitude (de)calib functions
    //// TODO: Set correct conditions
    if (peArea<=0)
    {
     omx->SetDead(1);
     omx->SetCalFunction(0,1);
     omx->SetDecalFunction(0,1);
    }
    // Set calibration function parameters
    fcal=omx->GetCalFunction(1);
    fdecal=omx->GetDecalFunction(1);
    if (fcal)
    {
     fcal->SetParameter(0,peArea);
    }
    if (fdecal)
    {
     fdecal->SetParameter(0,peArea);
     if (!peArea) fdecal->SetParameter(0,1);
    }
   }
  }
 }

 /*
 // Flag OMs in bad OM list as dead
 for(omid=0; omid<=maxomid; omid++) revision[omid]=0;
 st=server->Statement("SELECT ValidityStartDate, ValidityEndDate, RevisionId, StringId, TubeId FROM XXXBadOm INNER JOIN CalibrationDetail WHERE StringId>0 AND XXXBadOm.CaId=CalibrationDetail.CaId AND CalibrationDetail.TypeId=XXX;");
 if (!st)
 {
  cout << " *IceDB2Root GetJEBADaqData* Bad OM list could not be read from DB" << endl;
 }
 else {
  st->Process();
  st->StoreResult();
  while(st->NextResultRow())
  {
   // Only get calibrations with correct validity range, and update with latest revision in case of several valid calibrations
   validitystart.SetUT(st->GetYear(0),st->GetMonth(0),st->GetDay(0),st->GetHour(0),st->GetMinute(0),st->GetSecond(0));
   validityend.SetUT(st->GetYear(1),st->GetMonth(1),st->GetDay(1),st->GetHour(1),st->GetMinute(1),st->GetSecond(1));
   omid=om.GetOMId(st->GetInt(3),st->GetInt(4));
   if(validitystart.GetDifference(fTime,"d",1)>0 && validityend.GetDifference(fTime,"d",1)<0 && st->GetInt(2) > revision[omid])
   {
    revision[omid]=st->GetInt(2);
    // Get OM, create a new one if necessary
    omx=(IceDOM*)fJEBADaqdb->GetObject(omid,1);
    if (!omx)
    {
     omx=new IceDOM(om);
     omx->SetUniqueID(omid);
     fJEBADaqdb->EnterObject(omid,1,omx);
    }
    // Flag OMs as dead
    omx->SetDead(1);
    omx->SetCalFunction(0,1);
    omx->SetDecalFunction(0,1);
    omx->SetDead(2);
    omx->SetCalFunction(0,2);
    omx->SetDecalFunction(0,2);
    omx->SetDead(3);
    omx->SetCalFunction(0,3);
    omx->SetDecalFunction(0,3);
   }
  }
 }
 */

}
///////////////////////////////////////////////////////////////////////////
