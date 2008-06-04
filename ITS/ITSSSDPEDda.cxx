/**************************************************************************
- Contact: Oleksandr_Borysov aborysov@ts.infnf.it
- Link: /afs/infn.it/ts/user/efragiac/public/testCosm3125.001
- Run Type: 
- DA Type: LDC
- Number of events needed: ~500
- Input Files: raw_data_file_on_LDC, ssddaconfig txt, ssdddlmap.txt, badchannels.root
- Output Files: ./<EqId_Slot> ./ssddaldc_<LDCID>.root, FXS_name=ITSSSDda_<LDCID>.root 
                local files are persistent over runs: data source
- Trigger types used:
 **************************************************************************/


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "daqDA.h"
#include "AliITSHandleDaSSD.h" 
#include "TROOT.h"
#include "TPluginManager.h"

using namespace std;


struct ConfigStruct {
  Int_t    fNModuleProcess;
  string   fSsdDdlMap;
  string   fBadChannels;
  ConfigStruct() : fNModuleProcess(108), fSsdDdlMap(""), fBadChannels("") {}
};


Int_t SaveEquipmentCalibrationData(const AliITSHandleDaSSD  *ssddaldc, const Char_t  *fprefix = NULL);
Bool_t ReadDAConfigurationFile(const Char_t *configfname, AliITSHandleDaSSD *const ssddaldc, ConfigStruct& cfg);


int main( int argc, char** argv )
{
  const Char_t       *configfname = "ssddaconfig.txt";
  const Char_t       *bcfname = "badchannels.root"; 
  const Char_t       *ddlmfname = "ssdddlmap.txt";
  AliITSHandleDaSSD  *ssddaldc;
  TString             feefname, fcdbsave, lfname;
  Int_t               status;
  Char_t             *dafname = NULL;
  ConfigStruct        cfg; 

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");


  /* check that we got some arguments = list of files */
  if (argc<2) {
    fprintf(stderr, "Wrong number of arguments\n");
    return -1;
  }

  char *datafilename = argv[1];
  ssddaldc = new AliITSHandleDaSSD(datafilename);
  if (ssddaldc->IsZombie()) { 
   cerr << "Failed to process raw data file " << datafilename << "! Exit DA!\n";
   return -1;
  }

  lfname.Form("./%s", configfname);
  if (!(status = daqDA_DB_getFile(configfname, lfname.Data()))) {
    if (!ReadDAConfigurationFile(lfname.Data(), ssddaldc, cfg)) cerr << "Error reading configuration file "  << lfname.Data() << " !\n";
  } else fprintf(stderr, "Failed to import DA configuration file %s from the detector db: %d, %s \n", configfname, status, lfname.Data());
    
  if (cfg.fBadChannels.size() > 0) bcfname = cfg.fBadChannels.c_str();
  lfname.Form("./%s", bcfname);
  if (status = daqDA_DB_getFile(bcfname, lfname.Data())) {
    fprintf(stderr, "Failed to import the file %s from the detector db: %d, %s! Exit DA!\n", bcfname, status, lfname.Data());
    delete ssddaldc;
    return -3;
  }  
  if (!ssddaldc->ReadStaticBadChannelsMap(lfname.Data())) {
    cerr << "Error reading static bad channels map " << lfname.Data() << "! Exit DA!\n";
    delete ssddaldc;
    return -4;
  }	  
  
  if (cfg.fSsdDdlMap.size() > 0) ddlmfname = cfg.fSsdDdlMap.c_str();
  lfname.Form("./%s", ddlmfname);
  if (!(status = daqDA_DB_getFile(ddlmfname, lfname.Data()))) {
    if (!ssddaldc->ReadDDLModuleMap(lfname.Data())) cerr << "Error reading DDL map from file " << lfname.Data() << " !\n"; 
  } else {
    fprintf(stderr, "Failed to import file %s from the detector db: %d, %s \n", ddlmfname, status, lfname.Data());
    if (!ssddaldc->ReadDDLModuleMap()) cerr << "Failed to load the DDL map from AliITSRawStreamSSD!\n"; 
  }    
  
  if (!ssddaldc->ProcessRawData(cfg.fNModuleProcess))  {
     cerr << "Error !ssddaldc->ProcessRawData()" << endl;
     delete ssddaldc;
     return -1;
  }
  daqDA_progressReport(90);

  dafname = ".";
  if (ssddaldc->SaveCalibrationSSDLDC(dafname)) {
    cout << "SSDDA data are saved in " << dafname << endl;
    status = daqDA_FES_storeFile(dafname, "CALIBRATION");
    if (status) fprintf(stderr, "Failed to export file : %d\n", status);
  } else cerr << "Error saving DA data to the file!\n";

  feefname.Form("%s/ssddaldc_%i.root", ".", ssddaldc->GetLdcId());
  cout << "Saving feessdda data in " << feefname << endl;
  TFile *fileRun = new TFile (feefname.Data(),"RECREATE");
  if (fileRun->IsZombie()) {
    cerr << "Error open file " << feefname << endl;
    delete ssddaldc;
    delete fileRun;
    return -2;
  }  
  ssddaldc->Write();
  fileRun->Close();
  delete fileRun;

  fcdbsave.Form("ssddaldc_%i.root", ssddaldc->GetLdcId());
  status = daqDA_DB_storeFile(feefname.Data(), fcdbsave.Data());
  if (status) fprintf(stderr, "Failed to export file %s to the detector db: %d, %s \n", feefname.Data(), status, fcdbsave.Data());
  cout << SaveEquipmentCalibrationData(ssddaldc) << " files were uploaded to DetDB!\n";
  delete ssddaldc;
  daqDA_progressReport(100);
  return 0;
}



//__________________________________________________________________________________________
Int_t SaveEquipmentCalibrationData(const AliITSHandleDaSSD  *ssddaldc, const Char_t  *fprefix)
{
  TString feefilename;
  Int_t count = 0, status;
  for (Int_t ddli = 0; ddli < 16; ddli++) {
    for(Int_t adi = 1; adi <= 9; adi++) {
	  if (!ssddaldc->AdDataPresent(ddli, adi)) continue;  
	  if (fprefix) feefilename.Form("%s%i_%i", fprefix, ssddaldc->DdlToEquipmentId(ddli), adi);
	  else feefilename.Form("%i_%i", ssddaldc->DdlToEquipmentId(ddli), adi);
      if (ssddaldc->SaveEqSlotCalibrationData(ddli, adi, feefilename)) {
        status = daqDA_DB_storeFile(feefilename, feefilename);
        if (status) fprintf(stderr, "Error %i, failed to export file %s\n", status, feefilename.Data());
        else count++;
      }
    }
  }
  return count;
}


//__________________________________________________________________________________________
Bool_t ReadDAConfigurationFile(const Char_t *configfname, AliITSHandleDaSSD *const ssddaldc, ConfigStruct& cfg) 
{
// Dowload configuration parameters from configuration file or database
  const int nkwords = 8;
  char *keywords[nkwords] = {"ZsDefault", "OffsetDefault", "ZsFactor", "PedestalThresholdFactor", "CmThresholdFactor",
                             "NModulesToProcess", "DDLMapFile", "BadChannelsFile"};
  Int_t tmpint;
  Float_t tmpflt;
  fstream dfile;
  if (!configfname) {
    cerr << "No DA configuration file name is specified, defaul value are used!\n";
    return kFALSE;
  }
  dfile.open(configfname, ios::in);
  if (!dfile.is_open()) {
    cerr << "Error open DA configuration file " << configfname << " defaul value are used!\n";
    return kFALSE;
  }
  while (!dfile.eof()) {
    string str, keystr, tmpstr;
    getline(dfile, str);
    stringstream strline(str);
    strline >> keystr;
    if (keystr.size() == 0) continue;
    if ((keystr.at(0) == '#') ) continue;
    int ind = 0;
    while (keystr.compare(keywords[ind])) if (++ind == nkwords) break;
    switch (ind) {
	  case 0: 
              strline >> tmpint;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            ssddaldc->SetZsDefaul(tmpint);   
	            cout << "Default value for ZS thereshold " << keystr << ": " << ssddaldc->GetZsDefault() << endl;
              } break;   
	  case 1: 
	          strline >> tmpint;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            ssddaldc->SetOffsetDefault(tmpint);   
	            cout << "Default value for offset correction " << keystr << ": " << ssddaldc->GetOffsetDefault() << endl;
              } break;   
	  case 2: 
	          strline >> tmpflt;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            ssddaldc->SetZsFactor(tmpflt);   
                cout << keystr << ": " << ssddaldc->GetZsFactor() << endl;
              } break;
	  case 3: 
	          strline >> tmpflt;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            ssddaldc->SetPedestalThresholdFactor(tmpflt);
                cout << keystr << ": " << ssddaldc->GetPedestalThresholdFactor() << endl;
              } break;
	  case 4: 
	          strline >> tmpflt;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            ssddaldc->SetCmThresholdFactor(tmpflt);
                cout << keystr << ": " << ssddaldc->GetCmThresholdFactor() << endl;
              } break;
	  case 5: 
	          strline >> tmpint;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
	            cfg.fNModuleProcess = tmpint;
                cout << keystr << ": " << cfg.fNModuleProcess << endl;
              } break;
	  case 6: 
	          strline >> tmpstr;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
                if (tmpstr.size() == 0) continue;
                if (tmpstr.at(0) == '#') continue;
	            cfg.fSsdDdlMap = tmpstr;
                cout << keystr << ": " << cfg.fSsdDdlMap.c_str() << endl;
              } break;
	  case 7: 
	          strline >> tmpstr;
              if (strline.fail()) cerr << "Failed to read " << keystr << " value from DA configuration file!\n";
              else {
                if (tmpstr.size() == 0) continue;
                if (tmpstr.at(0) == '#') continue;
	            cfg.fBadChannels = tmpstr;
                cout << keystr << ": " << cfg.fBadChannels.c_str() << endl;
              } break;
	  default: 
	          cerr << keystr << " is not a key word, no assignment were made!\n"; 
    }
  }  
  return kTRUE;
}
