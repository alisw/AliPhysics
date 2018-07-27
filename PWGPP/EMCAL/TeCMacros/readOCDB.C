#include "readOCDB_Temperature.C"
#include "readOCDB_LED.C"

Int_t readTemp(Int_t *runs, Int_t nruns=1, const char *pername="unspec") 
{
  if (!gGrid)
    TGrid::Connect("alien://");
  AliCDBManager*  cdb = AliCDBManager::Instance();
  if (cdb->GetDefaultStorage()==0)
    cdb->SetDefaultStorage("raw://");

  Int_t ret = 0;

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=0;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << "*** Working on " <<  i+1 << "/" << nruns << " with run number " << rn << " ***" << endl;
    TInfo *ti = readOCDB_Temperature(rn,0);
    if (!ti) 
      continue;
    ti->Print();
    arr.Add(ti);
    ++ret;
  }

  TFile *outf = TFile::Open("tempinfo.root","update");
  arr.Write(Form("temperatures_%s",pername),TObject::kSingleKey);
  outf->ls();
  outf->Close();
  delete outf;
  Double_t frac = 100.*ret/nruns;
  cout << "Finished temperature objects with " << frac << " percent!" << endl;
  return ret;
}

Int_t readLed(Int_t *runs, Int_t nruns=1, const char *pername="unspec") 
{
  if (!gGrid)
    TGrid::Connect("alien://");
  AliCDBManager*  cdb = AliCDBManager::Instance();
  if (cdb->GetDefaultStorage()==0)
    cdb->SetDefaultStorage("raw://");

  Int_t ret = 0;

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=0;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << "*** Working on " <<  i+1 << "/" << nruns << " with run number " << rn << " ***" << endl;
    LInfo *ti = readOCDB_LED(rn,0);
    if (!ti) 
      continue;
    ti->Print();
    arr.Add(ti);
    ++ret;
  }

  TFile *outf = TFile::Open("ledinfo.root","update");
  arr.Write(Form("led_%s",pername),TObject::kSingleKey);
  outf->ls();
  outf->Close();
  delete outf;
  Double_t frac = 100.*ret/nruns;
  cout << "Finished LED objects with " << frac << " percent!" << endl;
  return ret;
}
 

void read_LHC16f(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {253961, 253958, 253957, 253956, 253951, 253834, 253826, 253825, 253820, 253819, 253813, 253757, 253756, 253755, 253753, 253751, 253682, 253681, 253680, 253660, 253659, 253658, 253614};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16f");
  readLed(runs,nruns,"lhc16f");
}

void read_LHC16g(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {254124, 254126, 254128, 254147, 254148, 254149, 254174, 254196, 254199, 254204, 254205, 254293, 254302, 254303, 254304, 254330, 254331, 254332};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16g");
  readLed(runs,nruns,"lhc16g");
}

void read_LHC16h(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {254381, 254394, 254395, 254396, 254422, 254476, 254479, 254529, 254564, 254568, 254576, 254577, 254578, 254581, 254586, 254589, 254604, 254606, 254607, 254608, 254621, 254629, 254630, 254632, 254640, 254644, 254646, 254648, 254649, 254651, 254652, 254653, 254654, 254670, 254691, 254701, 254704, 255009, 255010, 255011, 255245, 255246, 255248, 255249, 255251, 255252, 255253, 255255, 255256, 255275, 255276, 255350, 255351, 255352, 255415, 255418, 255419, 255420, 255421, 255440, 255456, 255463, 255465, 255466, 255467, 255469};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16h");
  readLed(runs,nruns,"lhc16h");
}

void read_LHC16i(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {255515, 255533, 255534, 255535, 255537, 255538, 255539, 255540, 255541, 255542, 255543, 255577, 255582, 255583, 255591, 255592, 255614, 255615, 255616, 255617, 255618, 255642, 255648, 255649, 255650};
 
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16i");
  readLed(runs,nruns,"lhc16i");
}

void read_LHC16j(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {256146, 256147, 256148, 256149, 256156, 256157, 256207, 256210, 256223, 256225, 256227, 256231, 256281, 256282, 256283, 256284, 256289, 256290, 256292, 256295, 256297, 256298, 256299, 256302, 256307, 256309, 256311, 256356, 256357, 256361, 256362, 256363, 256364, 256365, 256366, 256367, 256371, 256372, 256373, 256415, 256417, 256418, 256420};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16j");
  readLed(runs,nruns,"lhc16j");
}

void read_LHC16k(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {256504, 256506, 256510, 256512, 256514, 256552, 256554, 256556, 256560, 256561, 256562, 256564, 256565, 256567, 256589, 256591, 256592, 256619, 256620, 256653, 256658, 256676, 256677, 256681, 256684, 256691, 256692, 256694, 256695, 256697, 256782, 256797, 256879, 256911, 256913, 256924, 256926, 256941, 256942, 256944, 257011, 257012, 257021, 257026, 257028, 257071, 257075, 257077, 257078, 257079, 257080, 257082, 257083, 257084, 257092, 257100, 257136, 257137, 257138, 257139, 257140, 257141, 257142, 257144, 257145, 257190, 257204, 257206, 257209, 257224, 257260, 257318, 257320, 257322, 257364, 257381, 257382, 257433, 257457, 257474, 257487, 257490, 257491, 257492, 257530, 257531, 257537, 257539, 257540, 257541, 257560, 257561, 257562, 257563, 257564, 257565, 257566, 257587, 257590, 257594, 257605, 257606, 257630, 257635, 257642, 257682, 257684, 257687, 257688, 257689, 257691, 257692, 257694, 257697, 257724, 257725, 257733, 257734, 257735, 257737, 257754, 257765, 257773, 257797, 257798, 257799, 257800, 257803, 257804, 257850, 257853, 257855, 257892, 257893, 257908, 257912, 257936, 257937, 257939, 257957, 257958, 257960, 257963, 257971, 257979, 258012, 258014, 258017, 258019, 258039, 258041, 258042, 258045, 258046, 258048, 258049, 258059, 258060, 258062, 258063, 258107, 258108, 258109, 258112, 258113, 258114, 258117, 258178, 258197, 258198, 258202, 258203, 258204, 258256, 258257, 258258, 258259, 258270, 258271, 258272, 258273, 258274, 258276, 258278, 258279, 258280, 258299, 258301, 258302, 258303, 258306, 258307, 258336, 258359, 258387, 258388, 258393, 258426, 258454, 258456, 258477, 258485, 258498, 258499, 258537, 258545, 258546, 258551, 258560, 258567, 258569, 258571, 258572, 258574};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16k");
  readLed(runs,nruns,"lhc16k");
}

void read_LHC16l(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {258883, 258884, 258885, 258886, 258889, 258890, 258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 258998, 259000, 259086, 259088, 259090, 259091, 259095, 259096, 259099, 259117, 259118, 259140, 259141, 259164, 259204, 259257, 259263, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259379, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259469, 259470, 259471, 259473, 259477, 259546, 259559, 259560, 259561, 259563, 259606, 259608, 259609, 259649, 259650, 259668, 259697, 259703, 259704, 259711, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259792, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 259954, 259961, 259979, 260010, 260011, 260014, 260187};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16l");
  readLed(runs,nruns,"lhc16l");
}

void read_LHC16o(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {262395, 262396, 262397, 262398, 262399, 262418, 262419, 262422, 262423, 262424, 262425, 262426, 262427, 262428, 262429, 262430, 262450, 262451, 262492, 262528, 262532, 262533, 262563, 262567, 262568, 262569, 262570, 262571, 262572, 262574, 262583, 262593, 262594, 262624, 262628, 262632, 262635, 262679, 262682, 262683, 262686, 262705, 262706, 262708, 262712, 262713, 262716, 262717, 262719, 262723, 262727, 262728, 262734, 262760, 262768, 262776, 262777, 262778, 262841, 262844, 262849, 262853, 262858, 263331, 263332, 263487, 263496, 263529, 263647, 263648, 263652, 263653, 263654, 263657, 263662, 263663, 263682, 263689, 263690, 263691, 263737, 263738, 263739, 263741, 263743, 263744, 263784, 263785, 263786, 263787, 263789, 263790, 263792, 263793, 263803, 263813, 263814, 263823, 263824, 263829, 263830, 263861, 263863, 263866, 263905, 263916, 263917, 263922, 263923, 263977, 263978, 263979, 263981, 263984, 263985, 264033, 264035};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16o");
  readLed(runs,nruns,"lhc16o");
}

void read_LHC16p(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {264076, 264078, 264082, 264085, 264086, 264109, 264110, 264129, 264133, 264134, 264137, 264138, 264139, 264164, 264168, 264188, 264190, 264197, 264198, 264232, 264233, 264235, 264238, 264259, 264260, 264261, 264262, 264264, 264265, 264266, 264267, 264277, 264279, 264281, 264305, 264306, 264312, 264336, 264345, 264346, 264347};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  readTemp(runs,nruns,"lhc16p");
  readLed(runs,nruns,"lhc16p");
}

void read_LHC17p(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {282008,282016,282021,282025,282026,282027,282030,282031,282050,282051,282078,282098,282099,282118,282119,282120,282122,282123,282125,282126,282127,282146,282147,282189,282206,282224,282227,282229,282230,282247,282302,282303,282304,282305,282306,282307,282309,282312,282313,282314,282340,282341,282342,282343};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
  }

  readTemp(runs,nruns,"lhc17p");
  readLed(runs,nruns,"lhc17p");
}

void read_LHC18d(Bool_t loc=0, Bool_t test=0) 
{
  Int_t runs[] = {285978,285979,285980,286014,286018,286025,286026,286027,286030,286064,286124,286127,286129,286130,286154,286157,286159,286198,286201,286202,286203,286229,286230,286231,286254,286255,286256,286257,286258,286261,286263,286282,286284,286287,286288,286289,286308,286309,286310,286311,286312,286313,286314,286336,286337,286340,286341,286345,286348,286349,286350};
  
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  readTemp(runs,nruns,"lhc18d");
  readLed(runs,nruns,"lhc18d");
}

void read_all(Bool_t loc=0, Bool_t test=0)
{
  read_LHC16f(loc,test);
  read_LHC16g(loc,test);
  read_LHC16h(loc,test);
  read_LHC16i(loc,test);
  read_LHC16j(loc,test);
  read_LHC16k(loc,test);
  read_LHC16l(loc,test);
  read_LHC16o(loc,test);
  read_LHC16p(loc,test);
  read_LHC17p(loc,test);
  read_LHC18d(loc,test);
}

void readOCDB(Bool_t loc=0, Int_t runno=286154) 
{
  AliCDBManager*  cdb = AliCDBManager::Instance();
  if (loc)
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/test/");

  LInfo *li = readOCDB_LED(runno, 0);
  li->Compute();
  //li->Print();  
  cout << "fraction good strips ";
  for (Int_t i=0;i<20;++i) 
    cout << li->FracStrips(i) << " ";
  cout << endl;
  cout << "fraction good towers ";
  for (Int_t i=0;i<20;++i) 
    cout << li->FracLeds(i) << " ";
  cout << endl;
}
