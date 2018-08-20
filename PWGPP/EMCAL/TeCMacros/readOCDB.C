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

void read_LHC15l(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {239396, 239401, 239402, 239404, 239407, 239409, 239519, 239696, 239701, 239717, 240069, 240072, 240376, 240380, 240381, 240382, 240392, 240394, 240398, 240404, 240406, 240408, 240409, 240411, 240412, 240434, 240440, 240443, 240444, 240447, 240450, 240610, 240612, 240845, 240854, 240860, 240864, 240872, 240874, 240875, 240880, 241001, 241010, 241014, 241021, 241043, 241050, 241054, 241055, 241056, 241057, 241062, 241069, 241141, 241144, 241257, 241267, 241268, 241269, 241281, 241288, 241295, 241296, 241354, 241360, 241361, 241393, 241396, 241407, 241521, 241523, 241531, 241544};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2015/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc15l");
  readLed(runs,nruns,"lhc15l");
}

void read_LHC15n(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {244340, 244351, 244355, 244359, 244364, 244377, 244411, 244416, 244418, 244421, 244453, 244456, 244480, 244481, 244482, 244483, 244484, 244531, 244540, 244542, 244617, 244618, 244619, 244626, 244627, 244628};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2015/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc15n");
  readLed(runs,nruns,"lhc15n");
}


void read_LHC15o(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {245145, 245146, 245148, 245151, 245152, 245231, 245232, 245233, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245410, 245411, 245439, 245441, 245446, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554, 245683, 245700, 245702, 245705, 245729, 245731, 245738, 245752, 245759, 245793, 245829, 245831, 245833, 245949, 245952, 245954, 245963, 246001, 246003, 246037, 246042, 246052, 246053, 246087, 246089, 246113, 246115, 246148, 246151, 246152, 246153, 246178, 246180, 246181, 246182, 246217, 246222, 246225, 246271, 246272, 246275, 246276, 246424, 246433, 246434, 246487, 246488, 246493, 246495, 246540, 246543, 246567, 246568, 246575, 246583, 246648, 246671, 246675, 246676, 246750, 246751, 246757, 246758, 246759, 246760, 246765, 246766, 246804, 246805, 246807, 246808, 246809, 246810, 246844, 246845, 246846, 246865, 246867, 246870, 246871, 246928, 246930, 246945, 246948, 246980, 246982, 246984, 246989, 246991, 246994};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2015/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc15o");
  readLed(runs,nruns,"lhc15o");
}

void read_LHC16f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {253961, 253958, 253957, 253956, 253951, 253834, 253826, 253825, 253820, 253819, 253813, 253757, 253756, 253755, 253753, 253751, 253682, 253681, 253680, 253660, 253659, 253658, 253614};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16f");
  readLed(runs,nruns,"lhc16f");
}

void read_LHC16g(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {254124, 254126, 254128, 254147, 254148, 254149, 254174, 254196, 254199, 254204, 254205, 254293, 254302, 254303, 254304, 254330, 254331, 254332};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16g");
  readLed(runs,nruns,"lhc16g");
}

void read_LHC16h(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {254381, 254394, 254395, 254396, 254422, 254476, 254479, 254529, 254564, 254568, 254576, 254577, 254578, 254581, 254586, 254589, 254604, 254606, 254607, 254608, 254621, 254629, 254630, 254632, 254640, 254644, 254646, 254648, 254649, 254651, 254652, 254653, 254654, 254670, 254691, 254701, 254704, 255009, 255010, 255011, 255245, 255246, 255248, 255249, 255251, 255252, 255253, 255255, 255256, 255275, 255276, 255350, 255351, 255352, 255415, 255418, 255419, 255420, 255421, 255440, 255456, 255463, 255465, 255466, 255467, 255469};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16h");
  readLed(runs,nruns,"lhc16h");
}

void read_LHC16i(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {255515, 255533, 255534, 255535, 255537, 255538, 255539, 255540, 255541, 255542, 255543, 255577, 255582, 255583, 255591, 255592, 255614, 255615, 255616, 255617, 255618, 255642, 255648, 255649, 255650};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16i");
  readLed(runs,nruns,"lhc16i");
}

void read_LHC16j(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {256146, 256147, 256148, 256149, 256156, 256157, 256207, 256210, 256223, 256225, 256227, 256231, 256281, 256282, 256283, 256284, 256289, 256290, 256292, 256295, 256297, 256298, 256299, 256302, 256307, 256309, 256311, 256356, 256357, 256361, 256362, 256363, 256364, 256365, 256366, 256367, 256371, 256372, 256373, 256415, 256417, 256418, 256420};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16j");
  readLed(runs,nruns,"lhc16j");
}

void read_LHC16k(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {256504, 256506, 256510, 256512, 256514, 256552, 256554, 256556, 256560, 256561, 256562, 256564, 256565, 256567, 256589, 256591, 256592, 256619, 256620, 256653, 256658, 256676, 256677, 256681, 256684, 256691, 256692, 256694, 256695, 256697, 256782, 256797, 256879, 256911, 256913, 256924, 256926, 256941, 256942, 256944, 257011, 257012, 257021, 257026, 257028, 257071, 257075, 257077, 257078, 257079, 257080, 257082, 257083, 257084, 257092, 257100, 257136, 257137, 257138, 257139, 257140, 257141, 257142, 257144, 257145, 257190, 257204, 257206, 257209, 257224, 257260, 257318, 257320, 257322, 257364, 257381, 257382, 257433, 257457, 257474, 257487, 257490, 257491, 257492, 257530, 257531, 257537, 257539, 257540, 257541, 257560, 257561, 257562, 257563, 257564, 257565, 257566, 257587, 257590, 257594, 257605, 257606, 257630, 257635, 257642, 257682, 257684, 257687, 257688, 257689, 257691, 257692, 257694, 257697, 257724, 257725, 257733, 257734, 257735, 257737, 257754, 257765, 257773, 257797, 257798, 257799, 257800, 257803, 257804, 257850, 257853, 257855, 257892, 257893, 257908, 257912, 257936, 257937, 257939, 257957, 257958, 257960, 257963, 257971, 257979, 258012, 258014, 258017, 258019, 258039, 258041, 258042, 258045, 258046, 258048, 258049, 258059, 258060, 258062, 258063, 258107, 258108, 258109, 258112, 258113, 258114, 258117, 258178, 258197, 258198, 258202, 258203, 258204, 258256, 258257, 258258, 258259, 258270, 258271, 258272, 258273, 258274, 258276, 258278, 258279, 258280, 258299, 258301, 258302, 258303, 258306, 258307, 258336, 258359, 258387, 258388, 258393, 258426, 258454, 258456, 258477, 258485, 258498, 258499, 258537, 258545, 258546, 258551, 258560, 258567, 258569, 258571, 258572, 258574};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16k");
  readLed(runs,nruns,"lhc16k");
}

void read_LHC16l(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {258883, 258884, 258885, 258886, 258889, 258890, 258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 258998, 259000, 259086, 259088, 259090, 259091, 259095, 259096, 259099, 259117, 259118, 259140, 259141, 259164, 259204, 259257, 259263, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259379, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259469, 259470, 259471, 259473, 259477, 259546, 259559, 259560, 259561, 259563, 259606, 259608, 259609, 259649, 259650, 259668, 259697, 259703, 259704, 259711, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259792, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 259954, 259961, 259979, 260010, 260011, 260014, 260187};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16l");
  readLed(runs,nruns,"lhc16l");
}

void read_LHC16o(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {262395, 262396, 262397, 262398, 262399, 262418, 262419, 262422, 262423, 262424, 262425, 262426, 262427, 262428, 262429, 262430, 262450, 262451, 262492, 262528, 262532, 262533, 262563, 262567, 262568, 262569, 262570, 262571, 262572, 262574, 262583, 262593, 262594, 262624, 262628, 262632, 262635, 262679, 262682, 262683, 262686, 262705, 262706, 262708, 262712, 262713, 262716, 262717, 262719, 262723, 262727, 262728, 262734, 262760, 262768, 262776, 262777, 262778, 262841, 262844, 262849, 262853, 262858, 263331, 263332, 263487, 263496, 263529, 263647, 263648, 263652, 263653, 263654, 263657, 263662, 263663, 263682, 263689, 263690, 263691, 263737, 263738, 263739, 263741, 263743, 263744, 263784, 263785, 263786, 263787, 263789, 263790, 263792, 263793, 263803, 263813, 263814, 263823, 263824, 263829, 263830, 263861, 263863, 263866, 263905, 263916, 263917, 263922, 263923, 263977, 263978, 263979, 263981, 263984, 263985, 264033, 264035};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16o");
  readLed(runs,nruns,"lhc16o");
}

void read_LHC16p(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {264076, 264078, 264082, 264085, 264086, 264109, 264110, 264129, 264133, 264134, 264137, 264138, 264139, 264164, 264168, 264188, 264190, 264197, 264198, 264232, 264233, 264235, 264238, 264259, 264260, 264261, 264262, 264264, 264265, 264266, 264267, 264277, 264279, 264281, 264305, 264306, 264312, 264336, 264345, 264346, 264347};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc16p");
  readLed(runs,nruns,"lhc16p");
}

void read_LHC16q(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {265305, 265308, 265309, 265331, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265378, 265383, 265384, 265387, 265388, 265419, 265420, 265421, 265424, 265425, 265426, 265427, 265499, 265500, 265501, 265521, 265525};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc16q");
    readLed(runs,nruns,"lhc16q");
}

void read_LHC16r(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {265630, 265632, 265656, 265658, 265665, 265669, 265694, 265696, 265697, 265698, 265700, 265701, 265705, 265709, 265713, 265714, 265739, 265740, 265741, 265742, 265744, 265756, 265785, 265787, 265788, 265789, 265792, 265795, 265797, 265838, 265840, 265841, 266022, 266023, 266025, 266034, 266081, 266083, 266084, 266085, 266086, 266117, 266187, 266189, 266193, 266196, 266197, 266208, 266234, 266235, 266296, 266299, 266304, 266312, 266313, 266316, 266317, 266318};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc16r");
    readLed(runs,nruns,"lhc16r");
}

void read_LHC16s(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {266405, 266437, 266438, 266439, 266441, 266470, 266472, 266477, 266479, 266480, 266514, 266516, 266518, 266520, 266522, 266523, 266525, 266526, 266534, 266539, 266543, 266549, 266584, 266585, 266587, 266588, 266590, 266591, 266592, 266595, 266596, 266614, 266615, 266617, 266618, 266619, 266621, 266625, 266628, 266630, 266657, 266658, 266659, 266665, 266668, 266669, 266700, 266702, 266703, 266706, 266708, 266775, 266776, 266800, 266805, 266807, 266808, 266857, 266878, 266880, 266881, 266882, 266883, 266885, 266886, 266912, 266915, 266940, 266942, 266943, 266944, 266988, 266993, 266994, 266997, 266998, 267020, 267060, 267061, 267062, 267067, 267070, 267072, 267077, 267081, 267109, 267110, 267130, 267131};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc16s");
    readLed(runs,nruns,"lhc16s");
}

void read_LHC16t(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {267161, 267163, 267164, 267165, 267166};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2016/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc16t");
    readLed(runs,nruns,"lhc16t");
}


void read_LHC17c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {270531, 270578, 270581, 270598, 270601, 270661, 270663, 270665};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17c");
    readLed(runs,nruns,"lhc17c");
}

void read_LHC17d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {270762, 270766, 270767, 270768, 270770, 270771, 270772};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17d");
    readLed(runs,nruns,"lhc17d");
}

void read_LHC17f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {270854, 270855, 270856, 270861, 270865};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17f");
    readLed(runs,nruns,"lhc17f");
}

void read_LHC17g(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {270882, 270883, 270931, 270934, 270935, 270937, 270938, 270940, 271005, 271006, 271008, 271009, 271013, 271015, 271026, 271028, 271288, 271289, 271369, 271378, 271379, 271381, 271382, 271383, 271384, 271419, 271448, 271449, 271451, 271743, 271774, 271777};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17g");
    readLed(runs,nruns,"lhc17g");
}

void read_LHC17h(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {271871, 271873, 271878, 271879, 271880, 271881, 271886, 271908, 271911, 271912, 271915, 271916, 271925, 272075, 272076, 272100, 272101, 272123, 272151, 272152, 272153, 272154, 272155, 272156, 272335, 272340, 272359, 272400, 272619, 272746, 272753, 272755, 272762, 272763, 272764, 272782, 272783, 272784, 272828, 272870, 272871, 272873, 272880, 272903, 272905, 272932, 272934, 272947, 272949, 272983, 273009, 273077, 273099, 273100, 273101};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17h");
    readLed(runs,nruns,"lhc17h");
}

void read_LHC17i(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {273591, 273592, 273593, 273653, 273654, 273711, 273719, 273824, 273825, 273885, 273886, 273918, 273942, 273943, 273946, 274063, 274064, 274092, 274125, 274147, 274148, 274174, 274212, 274232, 274258, 274259, 274263, 274264, 274269, 274270, 274271, 274278, 274279, 274280, 274281, 274283, 274329, 274351, 274352, 274355, 274363, 274364, 274385, 274386, 274387, 274388, 274389, 274390, 274442};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17i");
    readLed(runs,nruns,"lhc17i");
}

void read_LHC17j(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {274591, 274593, 274594, 274595, 274596, 274601, 274653, 274657, 274667, 274668, 274669, 274670, 274671};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17j");
    readLed(runs,nruns,"lhc17j");
}

void read_LHC17k(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {274889, 274979, 275067, 275068, 275073, 275075, 275076, 275149, 275150, 275621, 275622, 275623, 275624, 275648, 275657, 275664, 275847, 275924, 276012, 276013, 276017, 276019, 276020, 276040, 276041, 276045, 276097, 276098, 276102, 276104, 276105, 276108, 276135, 276140, 276141, 276145, 276166, 276169, 276170, 276178, 276205, 276230, 276257, 276259, 276290, 276291, 276292, 276294, 276302, 276307, 276312, 276348, 276351, 276429, 276435, 276437};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17k");
    readLed(runs,nruns,"lhc17k");
}

void read_LHC17l(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {276552, 276553, 276556, 276557, 276608, 276644, 276669, 276670, 276671, 276672, 276674, 276675, 276762, 276916, 276917, 276920, 276967, 276969, 276970, 276971, 276972, 277015, 277016, 277017, 277037, 277075, 277076, 277079, 277082, 277087, 277088, 277091, 277121, 277155, 277180, 277181, 277182, 277183, 277188, 277189, 277193, 277194, 277196, 277197, 277250, 277256, 277257, 277262, 277293, 277310, 277312, 277314, 277360, 277383, 277384, 277385, 277389, 277416, 277417, 277418, 277470, 277472, 277473, 277476, 277477, 277478, 277479, 277530, 277531, 277534, 277536, 277537, 277574, 277575, 277576, 277577, 277718, 277720, 277721, 277722, 277723, 277725, 277745, 277746, 277747, 277748, 277749, 277794, 277799, 277800, 277801, 277802, 277834, 277841, 277842, 277845, 277847, 277848, 277876, 277897, 277898, 277899, 277900, 277901, 277903, 277907, 277930, 277952, 277987, 277988, 277989, 277996, 278055, 278077, 278079, 278080, 278089, 278092, 278093, 278094, 278095, 278121, 278122, 278126, 278127, 278158, 278162, 278163, 278164, 278165, 278166, 278167, 278189, 278191, 278215, 278216};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17l");
    readLed(runs,nruns,"lhc17l");
}

void read_LHC17m(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {278915, 278936, 278939, 278941, 278959, 278960, 278963, 278964, 278999, 279000, 279005, 279007, 279008, 279035, 279036, 279041, 279043, 279044, 279068, 279069, 279073, 279074, 279075, 279106, 279107, 279117, 279118, 279122, 279123, 279130, 279155, 279157, 279199, 279201, 279207, 279208, 279232, 279234, 279235, 279242, 279264, 279265, 279267, 279268, 279270, 279273, 279274, 279309, 279310, 279312, 279342, 279344, 279348, 279354, 279355, 279391, 279410, 279439, 279441, 279483, 279487, 279491, 279550, 279559, 279560, 279561, 279583, 279597, 279598, 279600, 279602, 279632, 279641, 279642, 279676, 279677, 279682, 279683, 279687, 279688, 279689, 279718, 279719, 279747, 279749, 279773, 279826, 279853, 279854, 279855, 279879, 279880, 279884, 279886, 279889, 279890, 279893, 279952, 279954, 279955, 279956, 279957, 279958, 279963, 279965, 279979, 279980, 279981, 279982, 279984, 280046, 280047, 280049, 280051, 280052, 280066, 280107, 280108, 280111, 280114, 280126, 280131, 280134, 280140};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17m");
    readLed(runs,nruns,"lhc17m");
}

void read_LHC17n(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {280234, 280235};
    Int_t nruns = sizeof(runs)/sizeof(Int_t);

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17n");
    readLed(runs,nruns,"lhc17n");
}

void read_LHC17o(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {280282, 280283, 280284, 280285, 280286, 280310, 280348, 280349, 280350, 280351, 280352, 280374, 280375, 280403, 280412, 280413, 280415, 280418, 280419, 280443, 280445, 280446, 280447, 280448, 280518, 280519, 280546, 280547, 280550, 280551, 280574, 280575, 280576, 280581, 280583, 280613, 280634, 280636, 280637, 280647, 280648, 280650, 280671, 280673, 280676, 280679, 280705, 280706, 280729, 280753, 280754, 280755, 280756, 280757, 280761, 280762, 280763, 280764, 280765, 280766, 280767, 280768, 280786, 280787, 280793, 280842, 280844, 280845, 280847, 280848, 280849, 280854, 280856, 280880, 280881, 280897, 280936, 280943, 280947, 280994, 280996, 280997, 280998, 280999, 281032, 281033, 281035, 281036, 281060, 281061, 281062, 281079, 281080, 281081, 281179, 281180, 281181, 281186, 281187, 281188, 281189, 281190, 281191, 281212, 281213, 281240, 281241, 281242, 281243, 281244, 281271, 281273, 281275, 281277, 281301, 281321, 281350, 281415, 281443, 281444, 281446, 281449, 281450, 281509, 281511, 281562, 281563, 281568, 281569, 281574, 281580, 281583, 281592, 281633, 281634, 281635, 281638, 281641, 281653, 281665, 281668, 281671, 281705, 281706, 281707, 281713, 281750, 281751, 281753, 281755, 281756, 281892, 281893, 281894, 281895, 281915, 281916, 281918, 281919, 281920, 281928, 281931, 281932, 281939, 281940, 281961};
    Int_t nruns = sizeof(runs)/sizeof(Int_t);

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17o");
    readLed(runs,nruns,"lhc17o");
}


void read_LHC17p(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {282008, 282016, 282021, 282025, 282026, 282027, 282030, 282031, 282050, 282051, 282078, 282098, 282099, 282118, 282119, 282120, 282122, 282123, 282125, 282126, 282127, 282146, 282147, 282189, 282206, 282224, 282227, 282229, 282230, 282247, 282302, 282303, 282304, 282305, 282306, 282307, 282309, 282312, 282313, 282314, 282340, 282341, 282342, 282343};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc17p");
  readLed(runs,nruns,"lhc17p");
}

void read_LHC17q(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {282365, 282366, 282367, 282391, 282392, 282393, 282398, 282399, 282402, 282411, 282415, 282437, 282440};
    Int_t nruns = sizeof(runs)/sizeof(Int_t);

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17q");
    readLed(runs,nruns,"lhc17q");
}

void read_LHC17r(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {282504, 282544, 282545, 282546, 282573, 282579, 282580, 282606, 282607, 282608, 282609, 282615, 282618, 282620, 282622, 282651, 282653, 282666, 282667, 282668, 282670, 282671, 282673, 282676, 282677, 282700, 282702, 282703, 282704};
    Int_t nruns = sizeof(runs)/sizeof(Int_t);

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2017/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc17r");
    readLed(runs,nruns,"lhc17r");
}



void read_LHC18b(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] =
  {284891, 284900, 284946, 285007, 285008, 285009, 285010, 285011, 285012, 285013, 285014, 285015, 285064, 285065, 285066, 285106, 285107, 285108, 285125, 285127, 285165, 285200, 285202, 285203, 285222, 285224, 285286, 285287, 285289, 285290, 285291, 285327, 285328, 285347, 285364, 285365, 285396, 285447};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18b");
  readLed(runs,nruns,"lhc18b");
}

void read_LHC18c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] =
  {285466, 285471, 285481, 285486, 285496, 285497, 285515, 285516, 285545, 285550, 285557, 285575, 285576, 285577, 285578, 285599, 285601, 285602, 285603, 285639, 285640, 285641, 285642, 285643, 285659, 285662, 285663, 285664, 285666, 285697, 285698, 285722, 285751, 285752, 285753, 285754, 285755, 285777, 285778, 285781, 285804, 285810, 285811, 285812, 285830, 285851, 285869, 285893, 285917, 285946, 285957, 285958};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18c");
  readLed(runs,nruns,"lhc18c");
}

void read_LHC18d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  // took out run: 285979 as many T sensors aren't working
  Int_t runs[] = {285978, 285980, 286014, 286018, 286025, 286026, 286027, 286030, 286064, 286124, 286127, 286129, 286130, 286154, 286157, 286159, 286198, 286201, 286202, 286203, 286229, 286230, 286231, 286254, 286255, 286256, 286257, 286258, 286261, 286263, 286282, 286284, 286287, 286288, 286289, 286308, 286309, 286310, 286311, 286312, 286313, 286314, 286336, 286337, 286340, 286341, 286345, 286348, 286349, 286350};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18d");
  readLed(runs,nruns,"lhc18d");
}

void read_LHC18e(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {286380, 286426, 286428, 286454, 286455, 286482, 286501, 286502, 286508, 286509, 286511, 286566, 286567, 286568, 286569, 286591, 286592, 286653, 286661, 286695, 286731, 286799, 286801, 286846, 286848, 286850, 286852, 286874, 286876, 286907, 286908, 286910, 286911, 286930, 286931, 286932, 286936, 286937};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18e");
  readLed(runs,nruns,"lhc18e");
}

void read_LHC18fLEDOnly(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
    Int_t runs[] = {287995, 287997, 287998, 287999, 288000, 288004, 288006, 288007, 288008, 288009, 288010, 288012, 288015, 288016, 288017, 288018, 288019, 288020, 288022, 288023, 288030, 288042, 288043, 288044, 288045, 288046, 288048, 288049, 288054, 288055, 288056, 288057, 288058, 288059, 288060, 288063, 288064, 288067, 288068, 288069, 288070, 288073, 288074, 288076, 288078, 288079, 288082, 288083, 288084};

    Int_t nruns = sizeof(runs)/sizeof(Int_t);
    if (test)
        nruns=3;

    if (loc) {
        AliCDBManager*  cdb = AliCDBManager::Instance();
        cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
    }

    if (!readLEDOnly) readTemp(runs,nruns,"lhc18fLED");
    readLed(runs,nruns,"lhc18fLED");
}

void read_LHC18f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {287000, 287021, 287063, 287064, 287071, 287077, 287137, 287155, 287185, 287202, 287209, 287248, 287249, 287250, 287254, 287283, 287323, 287324, 287343, 287346, 287347, 287353, 287355, 287356, 287360, 287380, 287381, 287385, 287387, 287388, 287389, 287413, 287480, 287481, 287484, 287486, 287513, 287516, 287517, 287520, 287521, 287573, 287575, 287576, 287654, 287656, 287657, 287876, 287877, 287883, 287884, 287885, 287911, 287912, 287913, 287915, 287923, 287941, 287975, 287977};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18f");
  readLed(runs,nruns,"lhc18f");
}

void read_LHC18g(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {288619};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18g");
  readLed(runs,nruns,"lhc18g");
}

void read_LHC18h(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {288804, 288806};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18h");
  readLed(runs,nruns,"lhc18h");
}

void read_LHC18i(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {289165, 289166, 289169, 289172, 289175, 289176, 289177, 289198, 289199, 289200, 289201};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18i");
  readLed(runs,nruns,"lhc18i");
}

void read_LHC18j(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {288943};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18j");
  readLed(runs,nruns,"lhc18j");
}

void read_LHC18k(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {289165, 289166, 289169, 289172, 289175, 289176, 289177, 289198, 289199, 289200, 289201};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18k");
  readLed(runs,nruns,"lhc18k");
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
  read_LHC16r(loc,test);
  read_LHC16s(loc,test);
  read_LHC16t(loc,test);
  read_LHC16q(loc,test);
  read_LHC17c(loc,test);
  read_LHC17d(loc,test);
  read_LHC17f(loc,test);
  read_LHC17g(loc,test);
  read_LHC17h(loc,test);
  read_LHC17i(loc,test);
  read_LHC17j(loc,test);
  read_LHC17k(loc,test);
  read_LHC17l(loc,test);
  read_LHC17m(loc,test);
  read_LHC17n(loc,test);
  read_LHC17o(loc,test);
  read_LHC17p(loc,test);
  read_LHC17q(loc,test);
  read_LHC17r(loc,test);
  read_LHC18b(loc,test);
  read_LHC18c(loc,test);
  read_LHC18d(loc,test);
  read_LHC18e(loc,test);
  read_LHC18fLEDOnly(loc,test);
  read_LHC18f(loc,test);
  read_LHC18g(loc,test);
  read_LHC18h(loc,test);
  read_LHC18i(loc,test);
  read_LHC18j(loc,test);
  read_LHC18k(loc,test);

}

void read_allLED(Bool_t loc=0, Bool_t test=0)
{
  read_LHC16f(loc,test,1);
  read_LHC16g(loc,test,1);
  read_LHC16h(loc,test,1);
  read_LHC16i(loc,test,1);
  read_LHC16j(loc,test,1);
  read_LHC16k(loc,test,1);
  read_LHC16l(loc,test,1);
  read_LHC16o(loc,test,1);
  read_LHC16p(loc,test,1);
  read_LHC16r(loc,test,1);
  read_LHC16s(loc,test,1);
  read_LHC16t(loc,test,1);
  read_LHC16q(loc,test,1);
  read_LHC17c(loc,test,1);
  read_LHC17d(loc,test,1);
  read_LHC17f(loc,test,1);
  read_LHC17g(loc,test,1);
  read_LHC17h(loc,test,1);
  read_LHC17i(loc,test,1);
  read_LHC17j(loc,test,1);
  read_LHC17k(loc,test,1);
  read_LHC17l(loc,test,1);
  read_LHC17m(loc,test,1);
  read_LHC17n(loc,test,1);
  read_LHC17o(loc,test,1);
  read_LHC17p(loc,test,1);
  read_LHC17q(loc,test,1);
  read_LHC17r(loc,test,1);
  read_LHC18b(loc,test,1);
  read_LHC18c(loc,test,1);
  read_LHC18d(loc,test,1);
  read_LHC18e(loc,test,1);
  read_LHC18fLEDOnly(loc,test,1);
  read_LHC18f(loc,test,1);
  read_LHC18g(loc,test,1);
  read_LHC18h(loc,test,1);
  read_LHC18i(loc,test,1);
  read_LHC18j(loc,test,1);
  read_LHC18k(loc,test,1);
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
