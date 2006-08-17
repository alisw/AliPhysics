{
gROOT->GetInterpreter()->AddIncludePath( "$ALIHLT_TOPDIR/BASE" );
gROOT->GetInterpreter()->AddIncludePath( "$ALIHLT_TOPDIR/TPCLib" );
gROOT->GetInterpreter()->AddIncludePath( "$ALIHLT_TOPDIR/src" );

gROOT->Macro("StartDisplayMacro.C"); 
}

