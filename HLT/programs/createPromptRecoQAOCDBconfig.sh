#!/usr/bin/env bash
main ()
{
  CONFIG=$(cat << EOF

histogram=\'CINT7ZAC-B-NOPF-CENT,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\'CTVXZAC-B-NOPF-CENT,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\'CINT7-B-NOPF-MUFAST,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\'C0TVX-B-NOPF-MUFAST,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\'C0V0M-B-NOPF-MUFAST,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\'CTVXV0M-B-NOPF-MUFAST,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
histogram=\',fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,\'
EOF
)

  aliroot -b << EOF
man=AliCDBManager::Instance();
man->SetDefaultStorage("local://${ALIHLT_HCDBDIR}");
TObjString config("${CONFIG}");
printf("config string:\n%s\n",config.GetName());
AliCDBId id("HLT/ConfigHLT/GlobalPromptRecoQAComponent",0,99999999);
AliCDBMetaData meta("mkrzewic@cern.ch");
man->Put(&config,id,&meta);
EOF
}

main "$@"
