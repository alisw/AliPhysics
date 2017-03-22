<html>
<head>
<title>HepMC Interface</title>
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

<form method='post' action='HepMCInterface.php'>

<h2>HepMC Interface</h2>

An interface to the HepMC [<a href="Bibliography.php" target="page">Dob01</a>] standard event record 
format has been provided by M. Kirsanov. To use it, the relevant 
libraries need to be linked, as explained in the <code>README</code> 
and <code>README.HepMC</code> files. Only version 2.06 (or later) 
of HepMC is supported as of 1 January 2013, by agreement with the 
LHC experimental community.

<p/>
The (simple) procedure to translate PYTHIA 8 events into HepMC ones 
is illustrated in the <code>main41.cc</code>, <code>main42.cc</code>
<code>main61.cc</code> and <code>main62.cc</code>   
main programs. At the core is a call to the
<pre>
HepMC::I_Pythia8::fill_next_event( pythia, hepmcevt, ievnum = -1) 
</pre>
which takes a reference of the generator object and uses it, on the one
hand, to read out and convert the event record in <code>pythia.event</code> 
and, on the other hand, to extract and store parton-density (PDF),
cross section and other information for the hard subprocess from 
<code>pythia.info</code>. There is also an alternative form that
does not requires access to the full <code>pythia</code> object, 
but only the event record, at the expense of a reduced information
storage, see below.

<p/>
While PYTHIA always stores momenta in units of GeV, with <i>c = 1</i>, 
HepMC nowadays can be built either for MeV or GeV as default, a choice
that can then be overridden on an event-by-event basis, see e.g. the 
<code>main41.cc</code> code. When filling the HepMC event record, PYTHIA 
will convert to the unit specified for the current HepMC event record. 

<p/>
Also for length units there are choices, and again the PYTHIA interface 
will convert to the units set for the HepMC event record. Here the mm
choice of PYTHIA seems to be shared by most other programs, and is
HepMC default. 

<p/>
The status code is now based on the new standard introduced for HepMC 2.05,
see the <?php $filepath = $_GET["filepath"];
echo "<a href='EventRecord.php?filepath=".$filepath."' target='page'>";?>Event::statusHepMC(...)</a> 
conversion routine for details. 

<p/>
The current values from <code>pythia.info.sigmaGen()</code> and 
<code>pythia.info.sigmaErr()</code> are stored for each event,
multiplied by <i>10^9</i> to convert from mb to pb. Note that
PYTHIA improves its accuracy by Monte Carlo integration in the course
of the run, so the values associated with the last generated event
should be the most accurate ones. If events also come with a dimensional 
weight, like in some Les Houches strategies, this weight is in units of pb.

<h2>The public methods</h2>

Here comes a complete list of all public methods of the 
<code>I_Pythia8</code> class in the <code>HepMC</code> 
(<i>not</i> <code>Pythia8</code>!) namespace.

<a name="method1"></a>
<p/><strong>I_Pythia8::I_Pythia8() &nbsp;</strong> <br/>
  
<strong>virtual I_Pythia8::~I_Pythia8() &nbsp;</strong> <br/>
the constructor and destructor take no arguments.
  

<a name="method2"></a>
<p/><strong>bool I_Pythia8::fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt, int ievnum = -1) &nbsp;</strong> <br/>
convert a <code>Pythia</code> event to a <code>HepMC</code> one.
Will return true if succeeded.
<br/><code>argument</code><strong> pythia </strong>  : 
the input <code>Pythia</code> generator object, from which both the 
event and other information can be obtained.
   
<br/><code>argument</code><strong> evt </strong>  : 
the output <code>GenEvt</code> event, in its standard form.
   
<br/><code>argument</code><strong> iev </strong>  : 
set the event number of the current event. If negative then the 
internal event number is used, which is incremented by one for
each new event.
  
  

<a name="method3"></a>
<p/><strong>bool I_Pythia8::fill_next_event( Pythia8::Event& pyev, GenEvent* evt, int ievnum = -1, Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0) &nbsp;</strong> <br/>
convert a <code>Pythia</code> event to a <code>HepMC</code> one.
Will return true if succeeded.
<br/><code>argument</code><strong> pyev </strong>  : 
the input <code>Pythia</code> event that is to be converted to HepMC
format.
   
<br/><code>argument</code><strong> evt </strong>  : 
the output <code>GenEvt</code> event, in its standard form.
   
<br/><code>argument</code><strong> iev </strong>  : 
set the event number of the current event. If negative then the 
internal event number is used, which is incremented by one for
each new event.
  
<br/><code>argument</code><strong> pyinfo </strong>  : 
pointer to the <code>Pythia Info</code> object, which is used to 
extract PFD values, and process and cross section information.
Without such a pointer this information therefore cannot be stored,
i.e. it is equivalent to the three <code>set_store</code> methods
below being set false. 
   
<br/><code>argument</code><strong> pyset </strong>  : 
pointer to the <code>Pythia Settings</code> object, which is used to 
decide whether hadronization is carried out, and therefore whether
to warn about unhadronized partons. Without such a pointer the
<code>set_free_parton_warnings</code> method below uniquely controls
the behaviour.
   
  

<p/>
The following paired methods can be used to set and get the status of 
some switches that modify the behaviour of the conversion routine. 
The <code>set</code> methods have the same default input values as 
the internal initialization ones, so that a call without an argument 
(re)stores the default.

<a name="method4"></a>
<p/><strong>void I_Pythia8::set_print_inconsistency(bool b = true) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::print_inconsistency() &nbsp;</strong> <br/>
print a warning line, on <code>cerr</code>, when inconsistent mother 
and daughter information is encountered.
  

<a name="method5"></a>
<p/><strong>void I_Pythia8::set_free_parton_warnings(bool b = true) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::free_parton_warnings() &nbsp;</strong> <br/>
check and print a warning line when unhadronized gluons or quarks are 
encountered in the event record. Does not apply when Pythia hadronization 
is switched off. Default is to do this check.
  

<a name="method6"></a>
<p/><strong>void I_Pythia8::set_crash_on_problem(bool b = false) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::crash_on_problem() &nbsp;</strong> <br/>
if problems (like the free partons above) are encountered then the run 
is interrupted by an <code>exit(1)</code> command. Default is not to crash.
  

<a name="method7"></a>
<p/><strong>void I_Pythia8::set_convert_gluon_to_0(bool b = false) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::convert_gluon_to_0() &nbsp;</strong> <br/>
the normal gluon identity code 21 is used also when parton density
information is stored, unless this optional argument is set true to
have gluons represented by a 0. This choice does not affect the 
normal event record, where a gluon is always 21. 
  

<a name="method8"></a>
<p/><strong>void I_Pythia8::set_store_pdf(bool b = true) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::store_pdf() &nbsp;</strong> <br/>
for each event store information on the two incoming flavours, their
x values and common factorization scale, and the values of the two 
parton distributions, <i>xf(x,Q)</i>.
  

<a name="method9"></a>
<p/><strong>void I_Pythia8::set_store_proc(bool b = true) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::store_proc() &nbsp;</strong> <br/>
for each event store information on the Pythia process code, the 
renormalization scale, and <i>alpha_em</i> and <i>alpha_s</i>
values used for the hard process. 
  

<a name="method10"></a>
<p/><strong>void I_Pythia8::set_store_xsec(bool b = true) &nbsp;</strong> <br/>
  
<strong>bool I_Pythia8::store_xsec() &nbsp;</strong> <br/>
for each event store information on the Pythia cross section and its error,
in pb, and the event weight. If events also come with a dimensional weight, 
like in some Les Houches strategies, this weight is in units of pb.
  

</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
