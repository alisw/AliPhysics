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

void read_LHC13b(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {195344, 195351, 195389, 195391, 195478, 195479, 195480, 195481, 195482, 195483};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13b");
  readLed(runs,nruns,"lhc13b");
}

void read_LHC10b(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {115393, 115399, 115401, 116102, 116288, 116402, 116403, 116643, 116645, 117050, 117052, 117053, 117059, 117060, 117063, 117099, 117109, 117112, 117116, 117220, 117222};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2010/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc10b");
  readLed(runs,nruns,"lhc10b");
}

void read_LHC10c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {119159, 119161, 119163, 119841, 119842, 119844, 119845, 119846, 119853, 119856, 119859, 119862, 120067, 120069, 120072, 120076, 120079, 120244, 120503, 120504, 120505, 120616, 120617, 120671, 120741, 120820, 120821, 120822, 120823, 120824, 120825, 120829};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2010/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc10c");
  readLed(runs,nruns,"lhc10c");
}

void read_LHC10d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {122374, 122375, 125630, 125632, 125633, 125842, 125843, 125844, 125847, 125848, 125849, 125850, 125851, 125855, 126004, 126007, 126008, 126073, 126078, 126081, 126082, 126088, 126090, 126097, 126158, 126160, 126167, 126168, 126284, 126285, 126351, 126352, 126359, 126403, 126404, 126405, 126406, 126407, 126408, 126409, 126422, 126424, 126425, 126432};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2010/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc10d");
  readLed(runs,nruns,"lhc10d");
}

void read_LHC10e(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {128486, 128494, 128495, 128498, 128503, 128504, 128505, 128506, 128582, 128590, 128592, 128594, 128596, 128605, 128609, 128611, 128615, 128621, 128677, 128678, 128777, 128778, 128819, 128820, 128823, 128824, 128833, 128834, 128835, 128836, 128843, 128850, 128853, 128855, 128913, 129042, 129512, 129513, 129514, 129515, 129516, 129519, 129520, 129521, 129523, 129524, 129525, 129527, 129528, 129536, 129540, 129586, 129587, 129599, 129639, 129641, 129647, 129650, 129651, 129652, 129653, 129659, 129666, 129723, 129726, 129729, 129734, 129735, 129736, 129738, 129742, 129744, 129959, 129960, 129961, 129962, 129966, 129983, 130149, 130151, 130157, 130158, 130168, 130172, 130178, 130343, 130354, 130356, 130358, 130360, 130375, 130480, 130481, 130517, 130519, 130696, 130704, 130793, 130795, 130798, 130799, 130834, 130840, 130842, 130844, 130847, 130848};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2010/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc10e");
  readLed(runs,nruns,"lhc10e");
}

void read_LHC10f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {133006, 133007, 133010, 133327, 133329, 133330, 133414, 133670, 133762, 133800, 133920, 133969, 133982};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2010/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc10f");
  readLed(runs,nruns,"lhc10f");
}

void read_LHC11a7(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {146277, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146141, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145290, 145289, 145288 };
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11a7");
  readLed(runs,nruns,"lhc11a7");
}


void read_LHC11a(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {146746, 146747, 146748, 146801, 146802, 146803, 146804, 146805, 146806, 146807, 146817, 146824, 146856, 146858 };
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11a");
  readLed(runs,nruns,"lhc11a");
}

void read_LHC11b(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {148531, 148534, 148538, 148541, 148544, 148547, 148549, 148553, 148556, 148559, 148565, 148569, 148576, 148582, 148592, 148601, 148625, 148630, 148645, 148648, 148659, 148663, 148708, 148711, 148719, 148800, 148838, 148843, 148844, 148847, 148850, 148852, 148853, 148854, 148856, 148857, 149068, 149070, 149071, 149072, 149113, 149127, 149129, 149130, 149133, 149134, 149767, 149880, 149881, 149883, 149884, 149890, 149927, 149929, 149930, 149931, 149960, 149961, 149975, 150059, 150060, 150160, 150162, 150163, 150211, 150212, 150248, 150250, 150252, 150254, 150256, 150259, 150374, 150375, 150420, 150421, 150423, 150427, 150428, 150433, 150434, 150437, 150438, 150440, 150499, 150500, 150518, 150629, 151636, 151638};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11b");
  readLed(runs,nruns,"lhc11b");
}

void read_LHC11c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {151655, 151660, 151661, 151664, 151672, 151674, 151681, 151689, 151724, 151849, 151850, 151851, 151852, 152003, 152011, 152015, 152208, 152256, 152284, 152306, 152309, 152311, 152312, 152319, 152320, 152321, 152322, 152366, 152367, 152368, 152369, 152371, 152377, 152455, 152488, 152512, 152513, 152567, 152568, 152570, 152581, 152591, 152599, 153116, 153223, 153232, 153296, 153362, 153363, 153369, 153371, 153373, 153415, 153465, 153533, 153536, 153539, 153541, 153542, 153544, 153548, 153552, 153558, 153560, 153566, 153570, 153571, 153583, 153587, 153589, 153591, 153594, 153702, 153709, 153718, 153725, 153726, 153727, 153728, 153733, 153738, 153776, 153777, 153779, 153784, 153794, 153796, 153798, 153805, 153807, 153808, 153873, 153875, 153876, 153997, 154001, 154002, 154018, 154024, 154026, 154030, 154031, 154039, 154056, 154060, 154066, 154070, 154083, 154125, 154126, 154129, 154130, 154132, 154136, 154138, 154141, 154143, 154145, 154151, 154158, 154163, 154207, 154211, 154219, 154220, 154221, 154222, 154234, 154235, 154252, 154257, 154261, 154264, 154266, 154269, 154270, 154273, 154281, 154283, 154286, 154289, 154293, 154296, 154315, 154316, 154317, 154382, 154383, 154385, 154448, 154478, 154480, 154482, 154483, 154485, 154495, 154583, 154591, 154594, 154596, 154597, 154598, 154599, 154726, 154732, 154733, 154742, 154745, 154748, 154750, 154753, 154755, 154786, 154789, 154793, 154796, 154808, 154903, 154908, 154929, 154930,154780, 154773, 154763,154787,154783};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11c");
  readLed(runs,nruns,"lhc11c");
}


void read_LHC11d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {156620, 156621, 156626, 156629, 156794, 156797, 156829, 156889, 156891, 156895, 156896, 157003, 157025, 157082, 157083, 157087, 157090, 157091, 157092, 157094, 157096, 157098, 157203, 157209, 157210, 157211, 157212, 157214, 157220, 157226, 157227, 157257, 157260, 157261, 157262, 157268, 157272, 157273, 157275, 157277, 157475, 157476, 157496, 157560, 157561, 157562, 157564, 157567, 157569, 157645, 157698, 157707, 157708, 157713, 157734, 157765, 157766, 157770, 157818, 157819, 157848, 157857, 157975, 157976, 158000, 158001, 158002, 158003, 158009, 158011, 158013, 158014, 158020, 158025, 158026, 158031, 158036, 158037, 158041, 158042, 158048, 158049, 158050, 158051, 158052, 158053, 158084, 158086, 158110, 158111, 158112, 158114, 158115, 158118, 158120, 158124, 158135, 158136, 158137, 158180, 158189, 158191, 158192, 158194, 158196, 158200, 158201, 158285, 158287, 158288, 158293, 158299, 158301, 158303, 158304, 158338, 158341, 158467, 158468, 158469, 158471, 158492, 158495, 158496, 158508, 158509, 158511, 158516, 158518, 158520, 158521, 158522, 158526, 158528, 158531, 158533, 158592, 158598, 158602, 158604, 158608, 158611, 158613, 158615, 158617, 158622, 158626, 158673, 158706, 158714, 158717, 158718, 158722, 158729, 158745, 158776, 158777, 158779, 158780, 158781, 158784, 158788, 158790, 158791, 158792, 158793, 158794, 158844, 158848, 158856, 158868, 158875, 158876, 158877, 158878, 158879, 159007, 159028, 159040, 159042, 159044, 159076, 159085, 159090, 159117, 159120, 159121, 159122, 159128, 159146, 159147, 159150, 159154, 159156, 159162, 159167, 159168, 159169, 159173, 159177, 159185, 159186, 159191, 159192, 159193, 159194, 159199, 159200, 159201, 159204, 159205, 159206, 159207, 159212, 159214, 159215, 159216, 159217, 159218, 159220, 159221, 159223, 159254, 159255, 159257, 159258, 159259, 159260, 159283, 159285, 159286, 159318, 159378, 159379, 159450, 159451, 159502, 159503, 159505, 159517, 159521, 159532, 159535, 159536, 159538, 159539, 159571, 159575, 159577, 159580, 159581, 159582, 159586, 159593, 159595, 159599, 159602, 159606, 159635};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11d");
  readLed(runs,nruns,"lhc11d");
}

void read_LHC11h(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 168988, 168992, 169035, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169584, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169846, 169855, 169858, 169859, 169923, 169965, 169975, 169981, 170027, 170036, 170038, 170040, 170083, 170085, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170207, 170208, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2011/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc11h");
  readLed(runs,nruns,"lhc11h");
}

void read_LHC12a(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {176715, 176730, 176749, 176752, 176753, 176849, 176854, 176859, 176924, 176926, 176927, 176929, 177011};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12a");
  readLed(runs,nruns,"lhc12a");
}


void read_LHC12b(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {177580, 177592, 177597, 177612, 177620, 177624, 177671, 177679, 177680, 177681, 177682, 177798, 177799, 177802, 177804, 177805, 177942};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12b");
  readLed(runs,nruns,"lhc12b");
}


void read_LHC12c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {179569, 179571, 179584, 179585, 179591, 179618, 179621, 179639, 179796, 179803, 179806, 179858, 179859, 179916, 179917, 179918, 179919, 179920, 180000, 180042, 180044, 180129, 180130, 180131, 180132, 180133, 180199, 180200, 180500, 180501, 180515, 180517, 180561, 180564, 180567, 180569, 180716, 180717, 180719, 180720, 182017, 182018, 182022, 182023, 182106, 182110, 182111, 182207, 182289, 182295, 182297, 182299, 182300, 182302, 182322, 182323, 182324, 182325, 182624, 182635, 182684, 182686, 182687, 182691, 182692, 182724, 182725, 182728, 182729, 182730, 182741, 182744};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12c");
  readLed(runs,nruns,"lhc12c");
}

void read_LHC12d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {183913, 183916, 183932, 183933, 183934, 183935, 183936, 183937, 183938, 183942, 183946, 184126, 184127, 184131, 184132, 184134, 184135, 184137, 184138, 184140, 184144, 184145, 184147, 184183, 184188, 184208, 184209, 184210, 184215, 184216, 184371, 184383, 184389, 184673, 184678, 184682, 184687, 184784, 184786, 185029, 185031, 185116, 185126, 185127, 185132, 185133, 185134, 185157, 185160, 185164, 185189, 185196, 185198, 185203, 185206, 185208, 185217, 185221, 185282, 185284, 185289, 185291, 185292, 185293, 185296, 185299, 185300, 185302, 185303, 185349, 185350, 185351, 185356, 185359, 185360, 185361, 185362, 185363, 185368, 185371, 185375, 185378, 185461, 185465, 185474, 185582, 185583, 185588, 185589, 185659, 185687, 185738, 185764, 185765, 185768, 185775, 185776, 185778, 185784, 186007, 186009, 186011, 186163, 186164, 186165, 186167, 186205, 186208, 186319, 186320};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12d");
  readLed(runs,nruns,"lhc12d");
}


void read_LHC12f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {186668, 186688, 186689, 186690, 186692, 186694, 186811, 186814, 186937, 186938, 186939, 186966, 186969, 186990, 186992, 186994, 187143, 187145, 187146, 187147, 187148, 187149, 187150, 187151, 187152, 187202, 187203, 187339, 187340, 187341, 187487, 187488, 187489, 187510, 187623, 187624, 187627, 187656, 187698, 187739, 187749, 187783, 187785, 187791, 187796, 188093, 1881010};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12f");
  readLed(runs,nruns,"lhc12f");
}


void read_LHC12h(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {189306, 189310, 189315, 189316, 189350, 189351, 189352, 189353, 189400, 189407, 189409, 189410, 189411, 189577, 189578, 189602, 189603, 189605, 189610, 189611, 189612, 189616, 189621, 189623, 189647, 189648, 189650, 189654, 189656, 189658, 189659, 189696, 189697, 189698, 190150, 190209, 190210, 190212, 190213, 190214, 190215, 190216, 190240, 190303, 190305, 190307, 190335, 190337, 190338, 190340, 190341, 190342, 190344, 190386, 190388, 190389, 190390, 190392, 190393, 190416, 190417, 190418, 190419, 190421, 190422, 190424, 190425, 190895, 190898, 190903, 190904, 190968, 190970, 190974, 190975, 190979, 190981, 190983, 190984, 191129, 191227, 191229, 191230, 191231, 191245, 191247, 191248, 191450, 191451, 192072, 192073, 192075, 192128, 192136, 192140, 192141, 192172, 192174, 192177, 192194, 192197, 192199, 192200, 192201, 192202, 192205, 192246, 192344, 192347, 192348, 192349, 192415, 192417, 192453, 192461, 192468, 192471, 192492, 192499, 192505, 192510, 192535, 192542, 192548, 192551, 192729, 192731, 192732};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12h");
  readLed(runs,nruns,"lhc12h");
}


void read_LHC12i(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {192772, 192775, 192778, 192779, 192820, 192822, 192824, 193004, 193005, 193007, 193008, 193010, 193011, 193014, 193047, 193049, 193051, 193092, 193093, 193094, 193097, 193148, 193155, 193156, 193187, 193188, 193189, 193194};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2012/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc12i");
  readLed(runs,nruns,"lhc12i");
}

void read_LHC13c(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {195529, 195531, 195566, 195567, 195568, 195592, 195593, 195596, 195633, 195635, 195644, 195673, 195675, 195677};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13c");
  readLed(runs,nruns,"lhc13c");
}

void read_LHC13d(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {195681, 195682, 195720, 195721, 195722, 195723, 195724, 195725, 195726, 195727, 195760, 195761, 195765, 195767, 195783, 195787, 195826, 195827, 195829, 195830, 195831, 195867, 195869, 195871, 195872, 195873};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13d");
  readLed(runs,nruns,"lhc13d");
}


void read_LHC13e(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {195935, 195949, 195950, 195954, 195955, 195958, 195989, 195994, 196000, 196006, 196085, 196089, 196090, 196091, 196099, 196105, 196107, 196185, 196187, 196194, 196197, 196199, 196200, 196201, 196203, 196208, 196214, 196308, 196309, 196310};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13e");
  readLed(runs,nruns,"lhc13e");
}

void read_LHC13f(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {196433, 196474, 196475, 196477, 196528, 196533, 196535, 196563, 196564, 196566, 196568, 196601, 196605, 196608, 196646, 196648, 196701, 196702, 196703, 196706, 196714, 196720, 196721, 196722, 196772, 196773, 196774, 196869, 196870, 196874, 196876, 196965, 196967, 196972, 196973, 196974, 197003, 197011, 197012, 197015, 197027, 197031, 197084, 197089, 197090, 197091, 197092, 197094, 197098, 197099, 197138, 197139, 197142, 197143, 197144, 197145, 197147, 197148, 197149, 197150, 197152, 197153, 197184, 197189, 197247, 197248, 197254, 197255, 197256, 197258, 197260, 197296, 197297, 197298, 197299, 197300, 197302, 197341, 197342, 197348, 197349, 197351, 197386, 197387, 197388};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13f");
  readLed(runs,nruns,"lhc13f");
}


void read_LHC13g(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {197471, 197496, 197497, 197499, 197500, 197501, 197529, 197531, 197553, 197555, 197583, 197584, 197608, 197609, 197610, 197611, 197618, 197643, 197669};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2013/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc13g");
  readLed(runs,nruns,"lhc13g");
}

void read_LHC15Calib(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0)
{
  Int_t runs[] = {235709, 235713, 235716, 235812, 235840, 235887, 235894, 235899, 236141, 236149, 236156, 236162, 236225, 236232, 236241, 236245, 236282, 236332, 236350, 236355, 236358, 236394, 236442, 236447, 236454, 236466, 236555, 236559, 236819, 236820, 236823, 236826, 236849, 236855, 236858, 236861, 236967, 236970, 237003, 237050, 237057, 237062, 237106, 237109, 237119, 237178, 237246, 237257, 237287, 237334, 237336, 237339, 237340, 237346, 237351, 237355, 237358, 237362, 237365, 237366, 237370, 237388, 237392, 237393, 237403, 237503, 237506, 237513, 237516, 237646, 237673, 237686, 237700, 237701, 237703, 237704, 237709, 237712, 237766, 237773, 237778, 237781, 237788, 237792, 237794, 237805, 237843, 237846, 237946, 237947, 237971, 237972, 237974, 237979, 237981, 238074, 238094, 238095, 238096, 238130, 238143, 238146, 238155, 238161, 238163, 238165, 238175, 238178, 238180, 238186, 238188, 238397};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2015/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc15Calib");
  readLed(runs,nruns,"lhc15Calib");
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
  { 284891, 284900, 284946, 285007, 285008, 285009, 285010, 285011, 285012, 285013, 285014, 285015, 285064, 285065, 285066, 285106, 285107, 285108, 285125, 285127, 285165, 285200, 285202, 285203, 285222, 285224, 285286, 285287, 285289, 285290, 285291, 285327, 285328, 285347, 285364, 285365, 285396, 285447};

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

void read_LHC18l(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 289240, 289241, 289242, 289243, 289253, 289254, 289275, 289276, 289277, 289278, 289280, 289281, 289300, 289303, 289306, 289308, 289309, 289330, 289331, 289353, 289354, 289355, 289356, 289357, 289361, 289362, 289363, 289365, 289366, 289367, 289368, 289369, 289370, 289373, 289374, 289426, 289444, 289462, 289463, 289465, 289468, 289493, 289494, 289521, 289547, 289574, 289576, 289577, 289579, 289581, 289582, 289625, 289632, 289634, 289654, 289657, 289658, 289659, 289660, 289664, 289666, 289721, 289723, 289724, 289729, 289731, 289732, 289757, 289775, 289808, 289809, 289811, 289814, 289815, 289816, 289817, 289818, 289830, 289849, 289852, 289854, 289855, 289856, 289857, 289879, 289880, 289884, 289885, 289888, 289890, 289891, 289928, 289930, 289935, 289940, 289943, 289965, 289966, 289971};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18l");
  readLed(runs,nruns,"lhc18l");
}

void read_LHC18m(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 290222, 290223, 290253, 290254, 290293, 290294, 290297, 290298, 290300, 290323, 290324, 290327, 290350, 290374, 290375, 290376, 290399, 290401, 290411, 290412, 290418, 290420, 290421, 290423, 290425, 290427, 290428, 290453, 290455, 290456, 290458, 290459, 290469, 290499, 290501, 290533, 290534, 290535, 290538, 290539, 290540, 290543, 290544, 290545, 290549, 290550, 290553, 290588, 290590, 290612, 290613, 290614, 290615, 290627, 290632, 290645, 290658, 290660, 290665, 290687, 290689, 290692, 290696, 290699, 290721, 290742, 290764, 290766, 290769, 290772, 290774, 290787, 290790, 290841, 290843, 290846, 290848, 290860, 290862, 290886, 290887, 290892, 290894, 290895, 290932, 290935, 290941, 290943, 290944, 290948, 290974, 290975, 290976, 290979, 290980, 291002, 291003, 291004, 291005, 291035, 291037, 291038, 291041, 291065, 291066, 291069, 291093, 291100, 291101, 291110, 291111, 291116, 291143, 291188, 291209, 291238, 291240, 291257, 291262, 291263, 291265, 291266, 291279, 291280, 291282, 291283, 291284, 291285, 291286, 291360, 291361, 291362, 291373, 291375, 291377, 291397, 291399, 291402, 291416, 291417, 291419, 291420, 291424, 291446, 291447, 291453, 291456, 291457, 291481, 291482, 291484, 291485, 291590, 291613, 291614, 291615, 291618, 291622, 291624, 291625, 291626, 291657, 291661, 291665, 291690, 291692, 291694, 291697, 291698, 291702, 291706, 291729, 291755, 291756, 291760, 291762, 291768, 291769, 291795, 291796, 291803, 291805, 291807, 291816, 291891, 291892, 291894, 291942, 291943, 291944, 291945, 291946, 291948, 291953, 291976, 291977, 291982, 292012, 292040, 292060, 292061, 292062, 292067, 292075, 292080, 292081, 292106, 292107, 292108, 292109, 292114, 292115, 292140, 292160, 292161, 292162, 292163, 292164, 292166, 292167, 292168, 292192, 292218, 292240, 292241, 292242, 292265, 292270, 292273, 292274, 292298, 292318, 292321, 292322, 292326, 292330, 292333, 292339, 292345, 292349, 292352, 292355, 292362, 292397, 292398, 292405, 292406, 292428, 292429, 292430, 292432, 292434, 292456, 292457, 292460, 292461, 292495, 292496, 292497, 292500, 292521, 292523, 292524, 292526, 292553, 292554, 292557, 292559, 292560, 292563, 292584, 292586, 292593, 292598, 292625, 292626, 292628, 292631, 292633, 292640, 292660, 292661, 292666, 292693, 292695, 292696, 292698, 292701, 292704, 292737, 292739, 292747, 292748, 292750, 292803, 292804, 292809, 292810, 292811, 292831, 292832, 292836, 292839};

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18m");
  readLed(runs,nruns,"lhc18m");
}



void read_LHC18n(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 293357, 293359, 293362 };

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18n");
  readLed(runs,nruns,"lhc18n");
}


void read_LHC18o(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 293368, 293386, 293391, 293392, 293412, 293413, 293424, 293474, 293475, 293494, 293496, 293497, 293507, 293542, 293543, 293570, 293571, 293578, 293582, 293587, 293588, 293686, 293689, 293690, 293691, 293695, 293696, 293740, 293741, 293770, 293773, 293774, 293776, 293799, 293805, 293806, 293809, 293829, 293830, 293831, 293855, 293856, 293886, 293891, 293893, 293896, 293898, 293923, 293924, 293925, 293926, 293927, 293928, 293929, 293930, 293955, 293957, 293958, 293962, 293963, 293966 };

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18o");
  readLed(runs,nruns,"lhc18o");
}

void read_LHC18p(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 294009, 294010, 294011, 294012, 294013, 294060, 294061, 294062, 294072, 294073, 294075, 294076, 294078, 294080, 294081, 294082, 294085, 294086, 294087, 294088, 294089, 294091, 294092, 294128, 294131, 294152, 294154, 294155, 294156, 294199, 294200, 294201, 294203, 294208, 294210, 294212, 294241, 294242, 294305, 294306, 294307, 294310, 294502, 294503, 294524, 294525, 294526, 294529, 294530, 294531, 294553, 294556, 294558, 294562, 294563, 294586, 294587, 294588, 294590, 294591, 294593, 294620, 294623, 294625, 294628, 294630, 294631, 294632, 294633, 294634, 294636, 294653, 294715, 294716, 294721, 294722, 294741, 294742, 294743, 294744, 294745, 294746, 294747, 294769, 294774, 294775, 294805, 294815, 294816, 294817, 294818, 294852, 294875, 294877, 294883, 294884, 294885, 294924, 294925 };

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18p");
  readLed(runs,nruns,"lhc18p");
}

void read_LHC18q(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 295581, 295584, 295585, 295586, 295587, 295588, 295589, 295610, 295611, 295612, 295615, 295665, 295666, 295667, 295668, 295671, 295673, 295675, 295676, 295677, 295712, 295714, 295716, 295717, 295718, 295719, 295720, 295721, 295723, 295725, 295753, 295754, 295755, 295756, 295758, 295759, 295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295829, 295831, 295854, 295855, 295856, 295859, 295860, 295861, 295863, 295872, 295881, 295908, 295909, 295913, 295936, 295937, 295941, 295942, 296008, 296009, 296013, 296016, 296060, 296061, 296062, 296063, 296065, 296066, 296068, 296123, 296128, 296132, 296133, 296134, 296135, 296142, 296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 296246, 296247, 296269, 296270, 296273, 296275, 296279, 296280, 296303, 296304, 296307, 296309, 296312, 296375, 296376, 296377, 296378, 296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 296510, 296511, 296512, 296514, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296615, 296616, 296618, 296619, 296622, 296623 };

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18q");
  readLed(runs,nruns,"lhc18q");
}


void read_LHC18r(Bool_t loc=0, Bool_t test=0, Bool_t readLEDOnly = 0){
  Int_t runs[] = { 296690, 296691, 296693, 296694, 296749, 296750, 296752, 296781, 296784, 296785, 296786, 296790, 296791, 296793, 296794, 296799, 296835, 296838, 296839, 296848, 296849, 296850, 296851, 296890, 296894, 296899, 296903, 296930, 296931, 296932, 296934, 296935, 296941, 296966, 296967, 296968, 296969, 296971, 296975, 296976, 296977, 296979, 297029, 297031, 297035, 297085, 297117, 297118, 297123, 297132, 297133, 297193, 297194, 297195, 297196, 297218, 297219, 297221, 297222, 297277, 297278, 297310, 297311, 297312, 297315, 297317, 297332, 297333, 297335, 297336, 297363, 297366, 297379, 297380, 297403, 297405, 297406, 297408, 297413, 297414, 297415, 297441, 297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540, 297541, 297542, 297544, 297555, 297557, 297558, 297588, 297589, 297590, 297595, 297623, 297624 };

  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  if (loc) {
    AliCDBManager*  cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local:///opt/alice/OCDB/2018/");
  }

  if (!readLEDOnly) readTemp(runs,nruns,"lhc18r");
  readLed(runs,nruns,"lhc18r");
}



void read_allRun2(Bool_t loc=0, Bool_t test=0)
{
  read_LHC15Calib(loc,test);
  read_LHC15l(loc,test);
  read_LHC15n(loc,test);
  read_LHC15o(loc,test);
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
  read_LHC18l(loc,test);
  read_LHC18m(loc,test);
  read_LHC18n(loc,test);
  read_LHC18o(loc,test);
  read_LHC18p(loc,test);
  read_LHC18q(loc,test);
  read_LHC18r(loc,test);
}

void read_allRun1(Bool_t loc=0, Bool_t test=0)
{
  read_LHC13b(loc,test);
  read_LHC13c(loc,test);
  read_LHC13d(loc,test);
  read_LHC13e(loc,test);
  read_LHC13f(loc,test);
  read_LHC13g(loc,test);
}

void read_allLED(Bool_t loc=0, Bool_t test=0)
{
  read_LHC15Calib(loc,test,1);
  read_LHC15l(loc,test,1);
  read_LHC15n(loc,test,1);
  read_LHC15o(loc,test,1);
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
  read_LHC18l(loc,test,1);
  read_LHC18m(loc,test,1);
  read_LHC18n(loc,test,1);
  read_LHC18o(loc,test,1);
  read_LHC18p(loc,test,1);
  read_LHC18q(loc,test,1);
  read_LHC18r(loc,test,1);
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
