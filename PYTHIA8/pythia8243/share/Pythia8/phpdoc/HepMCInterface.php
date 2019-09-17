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
 
An interface to the HepMC [<a href="Bibliography.php#refDob01" target="page">Dob01</a>] standard event record 
format has been provided by M. Kirsanov. The code is stored in 
<code>include/Pythia8Plugins/HepMC2.h</code>. To use it, 
the relevant libraries need to be linked, as explained in the 
<code>README</code> file. 
Only version 2.06 (or later) of HepMC is supported, by agreement 
with the LHC experimental community. 
 
<p/> 
The (simple) procedure to translate PYTHIA 8 events into HepMC ones 
is illustrated in the <code>main41.cc</code>, <code>main42.cc</code> 
and <code>main43.cc</code> main programs. At the core is a call to the 
<pre> 
HepMC::Pythia8ToHepMC::fill_next_event( pythia, hepmcevt, ievnum = -1) 
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
 
<a name="section0"></a> 
<h3>Direct output to HepMC</h3> 
 
Provided that PYTHIA is properly linked to HepMC, implementation of the 
most common user case (run PYTHIA with a runcard, output HepMC) exists. 
The example <code>main93</code> implements this. 
The sample command file <code>main93.cmnd</code> provides a good starting 
point. The line: 
<pre> 
    Main:writeHepMC = on 
</pre> 
is the switch needed to write a HepMC file. 
The example is then run with: 
<pre> 
    ./main93 -c main93.cmnd 
</pre> 
and a HepMC file is then written. 
 
<p/> 
There are several other useful command line options to <code>main93</code>. 
They are all displayed by running <code>./main93 -h</code>. 
 
<a name="section1"></a> 
<h3>The public methods</h3> 
 
Here comes a complete list of all public methods of the 
<code>Pythia8ToHepMC</code> class in the <code>HepMC</code> 
(<i>not</i> <code>Pythia8</code>!) namespace. 
 
<a name="anchor1"></a>
<p/><strong> Pythia8ToHepMC::Pythia8ToHepMC() &nbsp;</strong> <br/>
   
<a name="anchor2"></a>
<strong> virtual Pythia8ToHepMC::~Pythia8ToHepMC() &nbsp;</strong> <br/>
the constructor and destructor take no arguments. 
   
 
<a name="anchor3"></a>
<p/><strong> bool Pythia8ToHepMC::fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt, int ievnum = -1, bool append = false, GenParticle* rootParticle = 0, int iBarcode = -1) &nbsp;</strong> <br/>
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
   
<br/><code>argument</code><strong> append </strong>  :  
if <code>true</code> then the input event is appended to the current 
<code>HepMC</code> event record, rather than starting a new one. 
   
<br/><code>argument</code><strong> rootParticle </strong>  :  
the root particle that defines the new production vertex for the 
particles to be added in the <code>append = true</code> option. 
   
<br/><code>argument</code><strong> iBarcode </strong>  :  
used to set the bar code when <code>append = true</code>. 
If positive then start from <code>iBarcode</code> itself, if negative 
then start from the current size of the <code>HepMC</code> event record, 
and if 0 then set all bar codes to vanish. 
   
   
 
<a name="anchor4"></a>
<p/><strong> bool Pythia8ToHepMC::fill_next_event( Pythia8::Event& pyev, GenEvent* evt, int ievnum = -1, Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0, bool append = false, GenParticle* rootParticle = 0, int iBarcode = -1) &nbsp;</strong> <br/>
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
   
<br/><code>argument</code><strong> append </strong>  :  
if <code>true</code> then the input event is appended to the current 
<code>HepMC</code> event record, rather than starting a new one. 
   
<br/><code>argument</code><strong> rootParticle </strong>  :  
the root particle that defines the new production vertex for the 
particles to be added in the <code>append = true</code> option. 
   
<br/><code>argument</code><strong> iBarcode </strong>  :  
used to set the bar code when <code>append = true</code>. 
If positive then start from <code>iBarcode</code> itself, if negative 
then start from the current size of the <code>HepMC</code> event record, 
and if 0 then set all bar codes to vanish. 
   
   
 
<p/> 
The following paired methods can be used to set and get the status of 
some switches that modify the behaviour of the conversion routine. 
The <code>set</code> methods have the same default input values as 
the internal initialization ones, so that a call without an argument 
(re)stores the default. 
 
<a name="anchor5"></a>
<p/><strong> void Pythia8ToHepMC::set_print_inconsistency(bool b = true) &nbsp;</strong> <br/>
   
<a name="anchor6"></a>
<strong> bool Pythia8ToHepMC::print_inconsistency() &nbsp;</strong> <br/>
print a warning line, on <code>cerr</code>, when inconsistent mother 
and daughter information is encountered. 
   
 
<a name="anchor7"></a>
<p/><strong> void Pythia8ToHepMC::set_free_parton_exception(bool b = true) &nbsp;</strong> <br/>
   
<a name="anchor8"></a>
<strong> bool Pythia8ToHepMC::free_parton_exception() &nbsp;</strong> <br/>
check and throw an exception when unhadronized gluons or quarks are 
encountered in the event record. Does not apply when Pythia hadronization 
is switched off. Default is to do this check. If an exception is thrown the 
<code>PartonEndVertexException</code> class will return a warning message. 
The calling code can choose action to take, also having access to the 
location (<code>index()</code>) and species (<code>pdg_id()</code>) of a 
bad parton. 
   
 
<a name="anchor9"></a>
<p/><strong> void Pythia8ToHepMC::set_convert_gluon_to_0(bool b = false) &nbsp;</strong> <br/>
   
<a name="anchor10"></a>
<strong> bool Pythia8ToHepMC::convert_gluon_to_0() &nbsp;</strong> <br/>
the normal gluon identity code 21 is used also when parton density 
information is stored, unless this optional argument is set true to 
have gluons represented by a 0. This choice does not affect the 
normal event record, where a gluon is always 21. 
   
 
<a name="anchor11"></a>
<p/><strong> void Pythia8ToHepMC::set_store_pdf(bool b = true) &nbsp;</strong> <br/>
   
<a name="anchor12"></a>
<strong> bool Pythia8ToHepMC::store_pdf() &nbsp;</strong> <br/>
for each event store information on the two incoming flavours, their 
x values and common factorization scale, and the values of the two 
parton distributions, <i>xf(x,Q)</i>. 
   
 
<a name="anchor13"></a>
<p/><strong> void Pythia8ToHepMC::set_store_proc(bool b = true) &nbsp;</strong> <br/>
   
<a name="anchor14"></a>
<strong> bool Pythia8ToHepMC::store_proc() &nbsp;</strong> <br/>
for each event store information on the Pythia process code, the 
renormalization scale, and <i>alpha_em</i> and <i>alpha_s</i> 
values used for the hard process. 
   
 
<a name="anchor15"></a>
<p/><strong> void Pythia8ToHepMC::set_store_xsec(bool b = true) &nbsp;</strong> <br/>
   
<a name="anchor16"></a>
<strong> bool Pythia8ToHepMC::store_xsec() &nbsp;</strong> <br/>
for each event store information on the Pythia cross section and its error, 
in pb, and the event weight. If events also come with a dimensional weight, 
like in some Les Houches strategies, this weight is in units of pb. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
