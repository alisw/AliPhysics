rm -f galice.root
rm -f *.QA.*.root
aliroot "$ALICE_ROOT/HLT/global/physics/macros/testconfigCalib.C"'("GLOBAL-flat-esd-converter","AddTaskMacro=$ALICE_ROOT/HLT/global/macros/AddTaskAnalysisTaskExampleV.C() WriteAnalysisToFile=1 InitializeGeometry=0 ResetAfterPush=0")' $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root","local://$ALIHLT_HCDBDIR", 0, 20, "HLT", "chains=RootWriterTPCcalib ignore-hltout")'
