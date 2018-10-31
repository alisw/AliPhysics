#!/usr/bin/env bash
main ()
{
  CONFIG=$(cat << EOF
histogram=,fHistITSSPDvertexZ,SPD z vertex position,ITSSPDvertexZ,
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
