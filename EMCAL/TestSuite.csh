#! /bin/tcsh -f
#
# Prepare the environment
#
mkdir $HOME/EMCALTestSuite ; cd $HOME/EMCALTestSuite
cp $ALICE_ROOT/EMCAL/grunTestSuite.C . 
cp $ALICE_ROOT/EMCAL/ConfigTestSuite.C . 
#
#begin data
set MINHITS=68
set MAXHITS=72
set MINSDIGITS=68
set MAXSDIGITS=72
set MINDIGITS=47
set MAXDIGITS=51
set MINEMCRECPOINTS=1
set MAXEMCRECPOINTS=2
set MINCPVRECPOINTS=0
set MAXCPVRECPOINTS=2
set MINTRACKSEGMENTS=1
set MAXTRACKSEGMENTS=2
set MINRECPARTICLES=1
set MAXRECPARTICLES=2
#end data
#
echo EMCAL Test Suite run on `date` > TestSuite.log
echo "=======================================================================" >> TestSuite.log
echo AliROOT  version `which aliroot` >> TestSuite.log
echo Root     version `which root` >> TestSuite.log
echo "=======================================================================" >> TestSuite.log
#
# SIMULATION
#
echo "    **** SIMULATION **** " >> TestSuite.log
aliroot -b -q grunTestSuite.C\(100\) >>& TestSuite.log 
echo '{ '  > tempo.C
echo " Int_t minhits  = $MINHITS ;"  >> tempo.C 
echo " Int_t maxhits = $MAXHITS ;"  >> tempo.C 
echo ' AliEMCALGetter * gime = AliEMCALGetter::Instance("galice.root") ;'  >> tempo.C
echo ' TH1F * hitmul = new TH1F("hitmul", "EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' Int_t max = gime->MaxEvent() ; ' >> tempo.C
echo ' Int_t evt = 0 ; ' >> tempo.C
echo ' for ( evt = 0 ; evt < max ; evt++ ) { ' >> tempo.C
echo '  gime->Event(evt,"H") ; ' >> tempo.C
echo '  hitmul->Fill(gime->Hits()->GetEntries()) ; ' >> tempo.C
echo ' } ' >> tempo.C
echo ' TF1 * gaus = new TF1("gaus", "gaus", 0., 200.) ; ' >> tempo.C 
echo ' hitmul->Fit(gaus,"", "", 40, 100) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxhits ||  gaus->GetParameter(1) < minhits ) ' >> tempo.C  
echo '  printf("ERRORSIM 1") ; ' >> tempo.C
echo ' else ' >> tempo.C
echo '  printf("ERRORSIM 0") ; ' >> tempo.C
echo '} '  >> tempo.C
aliroot -b -q tempo.C >>& TestSuite.log 
set ERRORSIM = `cat TestSuite.log | grep ERRORSIM | awk '{print $2}'`
if($ERRORSIM) then
  echo ERROR
  uuencode  TestSuite.log  TestSuite.log | mail -s 'EMCAL INSTALLATION ERROR' yves.schutz@cern.ch
  exit(0)
endif 
rm tempo.C 
#
# RECONSTRUCTION
#
echo "    **** RECONSTRUCTION **** " >> TestSuite.log
aliroot -b >>& TestSuite.log <<EOF
 AliEMCALSDigitizer sd("galice.root") ; 
 sd.ExecuteTask("deb all") ; 
 AliEMCALDigitizer d("galice.root") ; 
 d.ExecuteTask("deb all") ; 
 AliEMCALClusterizerv1 cl("galice.root") ; 
 cl.ExecuteTask("deb all") ; 
 AliEMCALTrackSegmentMakerv1 ts("galice.root") ; 
 ts.ExecuteTask("deb all") ; 
 AliEMCALPIDv1 pd("galice.root") ; 
 pd.ExecuteTask("deb all") ; 
 .q
EOF
echo '{ '  > tempo.C
echo " Int_t minsdig  = $MINSDIGITS ;"  >> tempo.C 
echo " Int_t maxsdig  = $MAXSDIGITS ;"  >> tempo.C 
echo " Int_t mindig   = $MINDIGITS ;"  >> tempo.C 
echo " Int_t maxdig   = $MAXDIGITS ;"  >> tempo.C 
echo " Int_t minemcrp = $MINEMCRECPOINTS ;"  >> tempo.C 
echo " Int_t maxemcrp = $MAXEMCRECPOINTS ;"  >> tempo.C 
echo " Int_t mincpvrp = $MINCPVRECPOINTS ;"  >> tempo.C 
echo " Int_t maxcpvrp = $MAXCPVRECPOINTS ;"  >> tempo.C 
echo " Int_t mints    = $MINTRACKSEGMENTS ;"  >> tempo.C 
echo " Int_t maxts    = $MAXTRACKSEGMENTS ;"  >> tempo.C 
echo " Int_t minpa    = $MINRECPARTICLES ;"  >> tempo.C 
echo " Int_t maxpa    = $MAXRECPARTICLES ;"  >> tempo.C 
echo ' AliEMCALGetter * gime = AliEMCALGetter::Instance("galice.root") ;'  >> tempo.C
echo ' TH1F * sdigmul = new TH1F("sdigmul", " SDigits EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' TH1F * digmul  = new TH1F("digmul", " Digits EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' TH1F * emcrpmul= new TH1F("emcrpmul", " EMCRecPoints EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' TH1F * cpvrpmul= new TH1F("cpvrpmul", " CPVRecPoints EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' TH1F * tsmul   = new TH1F("tsmul", " TrackSegments EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' TH1F * pamul   = new TH1F("pamul", " RecParticles EMCAL Test Suite", 100, 0., 200.) ; ' >> tempo.C
echo ' Int_t max = gime->MaxEvent() ; ' >> tempo.C
echo ' Int_t evt = 0 ; ' >> tempo.C
echo ' for ( evt = 0 ; evt < max ; evt++ ) { ' >> tempo.C
echo '  gime->Event(evt,"SDRTP") ; ' >> tempo.C
echo '  sdigmul->Fill(gime->SDigits()->GetEntries()) ; ' >> tempo.C
echo '  digmul->Fill(gime->Digits()->GetEntries()) ; ' >> tempo.C
echo '  emcrpmul->Fill(gime->EmcRecPoints()->GetEntries()) ; ' >> tempo.C
echo '  cpvrpmul->Fill(gime->CpvRecPoints()->GetEntries()) ; ' >> tempo.C
echo '  tsmul->Fill(gime->TrackSegments()->GetEntries()) ; ' >> tempo.C
echo '  pamul->Fill(gime->RecParticles()->GetEntries()) ; ' >> tempo.C
echo ' } ' >> tempo.C
echo ' TF1 * gaus = new TF1("gaus", "gaus", 0., 200.) ; ' >> tempo.C 
echo ' sdigmul->Fit(gaus, "", "", 40., 100.) ; ' >> tempo.C 
echo ' sdigmul->Draw() ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxsdig ||  gaus->GetParameter(1) < minsdig ) ' >> tempo.C  
echo '  printf("ERRORREC 1 sdigits\n") ; ' >> tempo.C
echo ' digmul->Fit(gaus, "", "", 20., 80.) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxdig ||  gaus->GetParameter(1) < mindig ) ' >> tempo.C  
echo '  printf("ERRORREC 1 digits\n") ; '  >> tempo.C
echo ' emcrpmul->Fit(gaus, "", "", 0., 4.) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxemcrp ||  gaus->GetParameter(1) < minemcrp ) ' >> tempo.C  
echo '  printf("ERRORREC 1 emc recpoints\n") ; ' >> tempo.C
echo ' cpvrpmul->Fit(gaus, "", "", 0., 4.) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxcpvrp ||  gaus->GetParameter(1) < mincpvrp ) ' >> tempo.C  
echo '  printf("ERRORREC 1 cpv recpoints\n") ; ' >> tempo.C
echo ' tsmul->Fit(gaus, "", "", 0., 4.) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxts ||  gaus->GetParameter(1) < mints ) ' >> tempo.C  
echo '  printf("ERRORREC 1 track segments\n") ; ' >> tempo.C
echo ' pamul->Fit(gaus, "", "", 0., 4. ) ; ' >> tempo.C 
echo ' if( gaus->GetParameter(1) > maxpa ||  gaus->GetParameter(1) < minpa ) ' >> tempo.C  
echo '  printf("ERRORREC 1 recparticles\n") ; ' >> tempo.C
echo '} ' >> tempo.C
aliroot -b -q tempo.C  >>& TestSuite.log 
set ERRORREC = `cat TestSuite.log | grep ERRORREC | awk '{print $2}'`
if($ERRORREC) then
  echo ERROR
  uuencode  TestSuite.log  TestSuite.log | mail -s 'EMCAL INSTALLATION ERROR' yves.schutz@cern.ch
  exit(0)
endif 
#rm tempo.C 
#/cd $HOME
#rm -fr $HOME/EMCALTestSuite
