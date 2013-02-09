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
  static TChain * MakeChain(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1, Int_t startFile=0);
  static TChain * MakeChainRandom(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1, Int_t startFile=0);
  TDSet  * MakeSet(const char*fileIn, const char * treeName, const char *fName=0, Int_t maxFiles=-1);
  TDSet  * MakeSetRandom(const char*fileIn, const char * treeName,const char *fName=0, Int_t maxFiles=-1);
  static Bool_t FilterList(const char*inputList, const char*fileList, Int_t checkLevel);
  static Bool_t FilterListZip(const char*inputList, const char*fileList, Int_t checkLevel);
  Bool_t  XRDCopyDir(const char * idir, const char * files, const char *odir, Bool_t zip); 
  static Int_t  CheckTreeInFile(const char*fileName,const char*treeName, Int_t debugLevel=0, const char *branchName=0);
  static void JoinTreesIndex(const char * outputFile, const char * outputTree, const char *indexName, const char *inputTrees, Int_t debugLevel);
  static void   MakeTreeFromList(const char *fout, const char * treeOut, const char * treeIn, const char * flist, Bool_t debug);
  static void CacheFileList(const char * fileIn, const char* cachePrefix); 
 private:
  Int_t         fVerbose;          // verbso mode  - print command 
  TString       fUserName;         // user name
  UserGroup_t  *fUserGroup;        // user group info
  AliXRDPROOFtoolkit(const AliXRDPROOFtoolkit&);
  AliXRDPROOFtoolkit&operator=(const AliXRDPROOFtoolkit&);
  ClassDef(AliXRDPROOFtoolkit, 0); 
};
#endif
