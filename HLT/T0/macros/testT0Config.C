// rm -f galice.root ;  aliroot -l -q -b $ALICE_SOURCE/HLT/global/physics/macros/testConfigOnlineCalib.C $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root","local:///opt/HLT-DEV/HCDB", 23, 23, "HLT", "chains=RootWriter ignore-hltout")'

void testT0Config()
{
  gSystem->Load("libAliHLTT0");
	AliHLTPluginBase::InitInstance();
	AliHLTConfiguration calib1("TZEROrec", "T0Reconstruction", "GLOBAL-flat-esd-converter", "-pushback-period=0");
	AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "TZEROrec", "-directory testDir -datafile test.root");
}
