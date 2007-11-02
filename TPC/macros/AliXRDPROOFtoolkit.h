#ifndef _AliXRDPROOFtoolkit
#define _AliXRDPROOFtoolkit

#include <TObject.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <TSystem.h>

using namespace std;

class TObjArray;
class TTree;
class TChain;
class TDSet;



class AliXRDPROOFtoolkit : public TObject
{
  

 public :
  AliXRDPROOFtoolkit ();
  void Print(Option_t* option = " ") const;
  void AddMachine (const char*name);
 public :
  //
  // Interface for low priority users - NO ssh ACCESS to PROOF machines 
  //
  TChain * MakeChain(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1, Int_t startFile=0);
  TChain * MakeChainRandom(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1, Int_t startFile=0);
  TDSet  * MakeSet(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1);
  TDSet  * MakeSetRandom(const char*fileIn, const char * treeName,const char *fName=0, Int_t maxFiles=-1);
  //
  //
  //
  Bool_t ListOfFiles(const char*fileName, const char*path, const char*filter,  Bool_t displayMachine);
  static Bool_t FilterList(const char*inputList, const char*fileList, Int_t checkLevel);
  //
  //
  //
  Bool_t  XRDCopyDir(const char * idir, const char * files, const char *odir, Bool_t zip); 
  //
  void CheckFiles (const char*fileIn, UInt_t checkLevel, const char*treeToRetrieve, const char*varexp, const char*selection);

  static Int_t  CheckTreeInFile(const char*fileName,const char*treeName, Int_t debugLevel=0, const char *branchName=0);
  //
  //
  // Interface for users with ssh access to the machines
  //
  void      FilterSegFault(const char *filter="last");
  Bool_t    HasSegFault(const char * machine, const char * filter = "last");
  TTree *   DumpSys(Bool_t verbose=kTRUE);
  TTree *   DumpSys2(Bool_t verbose=kTRUE);
  TTree *   DumpFiles(Bool_t verbose=kTRUE);
  //

  Int_t         fVerbose;          // verbso mode  - print command 
 private :
  //
  //
  //
 private:
  vector <const TString *> listeMachine;  // list of slaves         
  TString       fUserName;              // user name
  UserGroup_t  *fUserGroup;      // user group info
  ClassDef(AliXRDPROOFtoolkit, 0); 
};
#endif
