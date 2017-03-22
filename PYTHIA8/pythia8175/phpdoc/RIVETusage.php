<html>
<head>
<title>RIVET usage</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='RIVETusage.php'>

<h2>RIVET usage</h2>

<a href="http://projects.hepforge.org/rivet/">RIVET</a> is a toolkit for
the validation of Monte Carlo event generators [<a href="Bibliography.php" target="page">Buc10</a>]. It
contains the results of many experimental analyses, so that generator
output can easily be compared to data, as well as providing a framework to
implement your own analyses.  Although using PYTHIA with RIVET is not
officially supported, some helpful hints are given below. The full RIVET
manual is available <a href="http://arxiv.org/abs/1003.0694">online</a>.

<br/><br/>
<h3>Using PYTHIA with RIVET</h3>
The following assumes that you already have RIVET installed. Instructions
for this may be found
<a href="http://projects.hepforge.org/rivet/trac/wiki/GettingStarted">
here</a>.

<br/><br/>
Events are passed from PYTHIA to RIVET using the HepMC format. PYTHIA must
be compiled with HepMC support, using the same version of HepMC used when
compiling RIVET. This is setup through the PYTHIA <code>configure</code>
script e.g.
<pre>
  ./configure --with-hepmc=/path/to/HepMC --with-hepmcversion=HepMC.version.number
</pre>
The PYTHIA library itself does not need to be recompiled.

<br/><br/>
The <code>examples/main42.cc</code> sample program can then be used to
generate events in HepMC format (which <code>examples/main61.cc</code> 
and <code>examples/main62.cc</code> extends by allowing LHAPDF usage 
and subruns). When in the <code>examples</code> directory, the main program 
can be built and used as follows
<pre>
  make main42
  ./main42 main42.cmnd main42.hepmc
</pre>
The first argument is the input file which provides the options for event
generation, while the second is the output file where the HepMC events
should be written.

<br/><br/>
This HepMC file may now be read and processed by RIVET
<pre>
  rivet --analysis=ANALYSIS_NAME main42.hepmc
</pre>
where <code>ANALYSIS_NAME</code> is a 
<a href="http://projects.hepforge.org/rivet/analyses">built-in RIVET
analysis</a>, or one you have created yourself. The output of RIVET is in
the form of <code>.aida</code> files, containing the histograms for the
analysis, which can be processed further with RIVET (see the
<a href="http://projects.hepforge.org/rivet/trac/wiki/FirstRivetRun">
RIVET documentation</a> for more details).

<br/><br/>
The above examples requires that (potentially large) HepMC events are stored
to disk before being read by RIVET. It is possible, instead, to pass the
events directly to RIVET as they are produced by using a <code>FIFO</code>
pipe. This is done with the <code>mkfifo</code> command
<pre>
  mkfifo my_fifo
  ./main42.exe main42.cmnd my_fifo &
  rivet --analysis=ANALYSIS_NAME my_fifo
</pre>
Note that <code>main42</code> is run in the background.

</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
