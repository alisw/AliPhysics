/*
  See documentation in the  https://indico.cern.ch/event/361713/contribution/4/material/slides/0.pdf

  Example usage:

  .L $ALICE_PHYSICS/../src/PWGPP/scripts/GetMlInfo.C
  CacheRCT("LHC11e",2);                         // cach RCT infor in the current working directory
  TTree * t = GetRCT("LHC11e",2);                // get the tree from cache directory
  t->Draw("muon_trigger");                      

  To use the code in jobs - password free .globus certificates needed

  

  To be added:
   1.) Cache also as a root file - chaining is possible - GetRCT should use this root interface
   2.) Make chache all with parameters
   3.) Put here description how to get password free certificate
   4.) Make unit test - nightly checkes of the tree consitency in $ALICE_PHYSICS/../src/PWGPP/scripts/test
       a.) We can read it - GetEntries()>0
       b.) Check that branches are the same  
      

*/


#include <string>
#include <iostream>
#include <sstream>


TTree * GetRCT(std::string period, int pass, std::string cacheDirectory=".")
{
  std::stringstream mifFileName;
  mifFileName << cacheDirectory << "/data/20" << period.substr(3,2) << "/" + period << "/pass" << pass << "/rct.mif";
  
  std::stringstream comment;
  comment << "rct for " << period << " pass" << pass;
  TTree *rct = new TTree("rct", comment.str().c_str());
  
  rct->ReadFile(mifFileName.str().c_str(), "", '\"');
  
  return rct;
}

void CacheRCT(std::string period, int pass, std::string cacheDirectory=".", std::string cert="$HOME/.globus/usercert.pem", std::string key="$HOME/.globus/userkey.pem") 
{
  std::stringstream path;
  path << cacheDirectory << "/data/20" << period.substr(3,2) << "/" + period << "/pass" << pass;
  gSystem->mkdir(path.str().c_str(), kTRUE);
  
  std::stringstream query;
  query << "configuration/index.jsp?partition=" << period << "&pass=" << pass << "&res_path=mif\"";
  
  GetFileFromMl(query.str(), path.str(), "rct", cert, key); 
}

void GetFileFromMl(std::string query, std::string directory, std::string name, std::string cert, std::string key) 
{
  std::stringstream wgetCommand;
  wgetCommand << "wget --no-check-certificate --secure-protocol=TLSv1" << " --certificate=" << cert << " --private-key=" << key 
              << " -o " << directory << "/" << name << ".log" <<  " -O " << directory << "/" << name << ".mif"  
              << " \"https://alimonitor.cern.ch/" << query;
  std::cout << wgetCommand.str() << "\n";
  gSystem->Exec(wgetCommand.str().c_str());
}

void CacheRCT() {
  const char *args[] = {          "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f", "LHC10g", "LHC10h",
                        "LHC11a", "LHC11b", "LHC11c", "LHC11d", "LHC11e", "LHC11f",           "LHC11h",
                        "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12e", "LHC12f", "LHC12g", "LHC11h", "LHC12i",
                        "LHC13a", "LHC13b", "LHC13c", "LHC13d", "LHC13e", "LHC13f", "LHC13g" };
  int n = sizeof(args) / sizeof(args[0]);
			
  for(int i=0; i<n; ++i) for(int j=1; j<5; ++j) CacheRCT(args[i],j);                      
}
// null
