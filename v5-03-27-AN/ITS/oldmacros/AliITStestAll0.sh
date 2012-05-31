#
set echo on
#
aliroot -b -q grun.C\(2,\"$ALICE_ROOT/ITS/Config_muon.C\"\) >& aliroot_muon.log 
#
mv galice.root galice_muon.root
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSHits2SDigits.C\(\"galice_muon.root\"\,\"galice_muonS.root\"\) >& aliroot_muonS.log
#
aliroot -b -q grun.C\(2,\"$ALICE_ROOT/ITS/Config_bg.C\"\) >& aliroot_bg.log
#
cp galice.root galice_all.root
mv galice.root galice_bg.root
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSHits2SDigits.C\(\"galice_bg.root\"\,\"galice_bgS.root\"\) >& aliroot_bgS.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSHits2Digits.C\(\"galice_bg.root\"\,\"galice_D.root\"\) >& aliroot_D.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSHits2Digits.C\(\"galice_all.root\"\,\"galice_all.root\"\) >& aliroot_allD.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSMerge.C\(\"galice_mD.root\",\"galice_muonS.root\",\"galice_bgS.root\"\) >& aliroot_m.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSDigits2RecPoints.C\(\"galice_mD.root\",\"galice_mR.root\"\) >& aliroot_mRec.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSDigits2RecPoints.C\(\"galice_D.root\",\"galice_R.root\"\) >& aliroot_Rec.log
#
aliroot -q -b $ALICE_ROOT/ITS/AliITSDigits2RecPoints.C\(\"galice_all.root\",\"galice_all.root\"\) >& aliroot_allRec.log
#

