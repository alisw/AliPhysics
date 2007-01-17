/* 
$Id$ 
*/

// This class runs the test TOF preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void TOFPreprocessor()
{
  gSystem->Load("$ALICE/SHUTTLE/TestShuttle/libTestShuttle.so");

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");

  // create AliTestShuttle instance
  AliTestShuttle* shuttle = new AliTestShuttle(0, 0, 1000);

  // Generation of "fake" input DCS data
  TMap* dcsAliasMap = CreateDCSAliasMap();

  // now give the alias map to the shuttle
  shuttle->SetDCSInput(dcsAliasMap);

  // processing files. for the time being, the files are local.
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "TOF", "DELAYS", "MON", "$ALICE_ROOT/TOF/ShuttleInput/Total.root");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "TOF", "RUNLevel", "MON", "$ALICE_ROOT/TOF/ShuttleInput/Partial.root");

  // instantiation of the preprocessor
  AliPreprocessor* pp = new AliTOFPreprocessor(shuttle);

  // preprocessing
  shuttle->Process();

  // checking the file which should have been created  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("TOF/Calib/OnlineDelay", 0);
  if (!entry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

  AliTOFDataDCS* output = dynamic_cast<AliTOFDataDCS*> (entry->GetObject());
  // If everything went fine, draw the result
  if (output)
    output->Draw();
  
}

TMap* CreateDCSAliasMap()
{
  // Creates a DCS structure
  // The structure is the following:
  // TMap (key --> value)
  // <DCSAlias> --> <valueList>
  // <DCSAlias> is a string
  // <valueList> is a TObjArray of AliDCSValue


  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);

  TRandom random;
  TDatime *datime = new TDatime();
  Int_t time = datime->GetTime();
  Int_t date = datime->GetDate();
  Int_t pid  = gSystem->GetPid();
  delete datime;
  Int_t iseed = TMath::Abs(10000 * pid + time - date); 

  Float_t tentHVv=6500, tentHVi=80, tentLVv=2.7, tentLVi=4.5,
    tentLVv33=3.3, tentLVv50=5.0, tentLVv48=48,
    tentLVi33=100, tentLVi50=3.0, tentLVi48=10,
    tentFEEthr=1.0, tentFEEtfeac=25, tentFEEttrm=45, 
    tentTemp=25, tentPress=900;
  Float_t sigmaHVv=10, sigmaHVi=10, sigmaLVv=0.2, sigmaLVi=1.0,
    sigmaLVv33=0.1, sigmaLVv50=0.1, sigmaLVv48=1,
    sigmaLVi33=10, sigmaLVi50=0.5, sigmaLVi48=2,
    sigmaFEEthr=0.1, sigmaFEEtfeac=10, sigmaFEEttrm=4, 
    sigmaTemp=1, sigmaPress=10;

  Float_t tent=0, sigma=0, thr=0;
  Int_t NAliases=10514, NHV=90, NLV=576, NLV33=72, NLV50=72, NLV48=72, NFEEthr=1152, NFEEtfeac=576, NFEEttrm=6840, NT=1, NP=1;

  for(int nAlias=0;nAlias<NAliases;nAlias++) {

    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString sindex;
    TString aliasName;
    if (nAlias<NHV){
      aliasName = "tof_hv_vp_";
      sindex.Form("%02i",nAlias);
      aliasName += sindex;
      //aliasName += nAlias;
      tent=tentHVv;
      sigma=sigmaHVv;
      //      thr=thrHVv;
    }
    else if (nAlias<NHV*2){
      //      aliasName = "HVvneg";
      //aliasName += nAlias-NHV;
      aliasName = "tof_hv_vn_";
      sindex.Form("%02i",nAlias-NHV);
      aliasName += sindex;
      tent=-tentHVv;
      sigma=-sigmaHVv;
      //thr=-thrHVv;
    }
    else if (nAlias<NHV*3){
      //      aliasName = "HVcpos";
      //aliasName += nAlias-2*NHV;
      aliasName = "tof_hv_ip_";
      sindex.Form("%02i",nAlias-2*NHV);
      aliasName += sindex;
      tent=tentHVi;
      sigma=sigmaHVi;
      //thr=thrHVc;
    }
    else if (nAlias<NHV*4){
      //      aliasName = "HVcneg";
      //aliasName += nAlias-3*NHV;
      aliasName = "tof_hv_in_";
      sindex.Form("%02i",nAlias-3*NHV);
      aliasName += sindex;
      tent=-tentHVi;
      sigma=-sigmaHVi;
      //thr=-thrHVc;
    }
    else if (nAlias<NHV*4+NLV){
      //      aliasName = "LVv";
      //aliasName += nAlias-4*NHV;
      aliasName = "tof_lv_vfea_";
      sindex.Form("%03i",nAlias-4*NHV);
      aliasName += sindex;
      tent=tentLVv;
      sigma=sigmaLVv;
      //thr=thrLVv;
    }
    else if (nAlias<NHV*4+2*NLV){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_ifea_";
      sindex.Form("%03i",nAlias-(4*NHV+NLV));
      aliasName += sindex;
      tent=tentLVi;
      sigma=sigmaLVi;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+NLV33){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_v33_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV));
      aliasName += sindex;
      tent=tentLVv33;
      sigma=sigmaLVv33;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_i33_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV+NLV33));
      aliasName += sindex;
      tent=tentLVi33;
      sigma=sigmaLVi33;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+NLV50){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_v50_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV+2*NLV33));
      aliasName += sindex;
      tent=tentLVv50;
      sigma=sigmaLVv50;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_i50_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV+2*NLV33+NLV50));
      aliasName += sindex;
      tent=tentLVi50;
      sigma=sigmaLVi50;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+NLV48){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_v48_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV+2*NLV33+2*NLV50));
      aliasName += sindex;
      tent=tentLVv48;
      sigma=sigmaLVv48;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48){
      //      aliasName = "LVc";
      //aliasName += nAlias-(4*NHV+NLV);
      aliasName = "tof_lv_i48_";
      sindex.Form("%02i",nAlias-(4*NHV+2*NLV+2*NLV33+2*NLV50+NLV48));
      aliasName += sindex;
      tent=tentLVi48;
      sigma=sigmaLVi48;
      //thr=thrLVc;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr){
      //      aliasName = "FEEthr";
      //aliasName += nAlias-(4*NHV+2*NLV-(4*NHV+2*NLV+2*NLV33+2*NLV50+2*NLV48));
      aliasName = "tof_fee_th_";
      sindex.Form("%04i",nAlias-(4*NHV+2*NLV+2*NLV33+2*NLV50+2*NLV48));
      aliasName += sindex;
      tent=tentFEEthr;
      sigma=sigmaFEEthr;
      //thr=thrFEEthr;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr+NFEEtfeac){
      //cout << " nalias fee temp = " << nAlias << endl;
      //      aliasName = "FEEt";
      //aliasName += nAlias-(4*NHV+2*NLV+NFEEthr);
      aliasName = "tof_fee_tfeac_";
      sindex.Form("%03i",nAlias-(4*NHV+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr));
      aliasName += sindex;
      //cout << " nalias fee temp name = " << aliasName << endl;
      tent=tentFEEtfeac;
      sigma=sigmaFEEtfeac;
      //thr=thrFEEthr;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr+NFEEtfeac+NFEEttrm){
      //cout << " nalias fee temp = " << nAlias << endl;
      //      aliasName = "FEEt";
      //aliasName += nAlias-(4*NHV+2*NLV+NFEEthr);
      aliasName = "tof_fee_ttrm_";
      sindex.Form("%04i",nAlias-(4*NHV+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr+NFEEtfeac));
      aliasName += sindex;
      //cout << " nalias fee temp name = " << aliasName << endl;
      tent=tentFEEttrm;
      sigma=sigmaFEEttrm;
      //thr=thrFEEthr;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr+NFEEtfeac+NFEEttrm+1){
      cout << " putting temperature" << endl;
      aliasName = "temperature";
      tent=tentTemp;
      sigma=sigmaTemp;
      //thr=thrTemp;
    }
    else if (nAlias<NHV*4+2*NLV+2*NLV33+2*NLV50+2*NLV48+NFEEthr+NFEEtfeac+NFEEttrm+2){
      cout << " putting pressure" << endl;
      aliasName = "pressure";
      tent=tentPress;
      sigma=sigmaPress;
      //thr=thrPress;
    }


    // gauss generation of values 
    for (int timeStamp=0;timeStamp<1000;timeStamp+=10){
      Float_t gaussvalue = (Float_t) (random.Gaus(tent,sigma));
      if (TMath::Abs(gaussvalue-tent)>sigma){
	AliDCSValue* dcsVal = new AliDCSValue(gaussvalue, timeStamp);
	valueSet->Add(dcsVal);
      }
    }

    aliasMap->Add(new TObjString(aliasName), valueSet);
  }

  return aliasMap;
}

TMap* ReadDCSAliasMap()
{
  // Open a file that contains DCS input data
  // The CDB framework is used to open the file, this means the file is located
  // in $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/DCS/Data
  // The file contains an AliCDBEntry that contains a TMap with the DCS structure.
  // An explanation of the structure can be found in CreateDCSAliasMap()

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TOF/DCS/Data", 0);
  return dynamic_cast<TMap*> (entry->GetObject());
}

void WriteDCSAliasMap()
{
  // This writes the output from CreateDCSAliasMap to a CDB file

  TMap* dcsAliasMap = CreateDCSAliasMap();

  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Chiara");
	metaData.SetComment("Test object for TOFPreprocessor.C");

  AliCDBId id("TOF/DCS/Data", 0, 0);

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");

  AliCDBManager::Instance()->Put(dcsAliasMap, id, &metaData);
}
