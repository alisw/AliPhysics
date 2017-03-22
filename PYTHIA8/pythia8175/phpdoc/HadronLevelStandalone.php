<html>
<head>
<title>Hadron-Level Standalone</title>
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

<form method='post' action='HadronLevelStandalone.php'>

<h2>Hadron-Level Standalone</h2>

The Les Houches Accord allows external process-level configurations
to be fed in, for subsequent parton-level and hadron-level generation
to be handled internally by PYTHIA. There is no correspondingly 
standardized interface if you have external events that have also 
been generated through the parton-level stage, so that only the 
hadron-level remains to be handled. A non-standard way to achieve this
exists, however, and can be useful both for real applications and
for various tests of the hadronization model on its own. 

<p/>
The key trick is to set the flag <code>ProcessLevel:all = off</code>.
When <code>pythia.next()</code> is called it then does not try to
generate a hard process. Since there are no beams, it is also not 
possible to perform the normal <code>PartonLevel</code> step.
(It is still possible to generate final-state radiation, but this
is not automatic. It would have to be done by hand, using the 
<code>pythia.forceTimeShower(...)</code> method, before 
<code>pythia.next()</code> is called.) 
Thus only the <code>HadronLevel</code> methods are 
called, to take the current content of the event record stored in
<code>pythia.event</code> as a starting point for any hadronization 
and decays that are allowed by the normal parameters of this step. 
Often the input would consist solely of partons grouped into colour 
singlets, but also (colour-singlet) particles are allowed.

<p/>
To set up all the parameters, a <code>pythia.init()</code> call has
to be used, without any arguments. In brief, the structure of the
main program therefore should be something like
<pre>
  Pythia pythia;                               // Declare generator.
  Event& event = pythia.event                  // Convenient shorthand.
  pythia.readString("ProcessLevel:all = off"); // The trick!
  pythia.init();                               // Initialization.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    // Insert filling of event here!
    pythia.next();                             // Do the hadron level.
  }
</pre> 
Of course this should be supplemented by analysis of events, error checks,
and so on, as for a normal PYTHIA run. The unique aspect is how to fill
the <code>event</code> inside the loop, before <code>pythia.next()</code> 
is called.

<h3>Input configuration</h3>

To set up a new configuration the first step is to throw away the current
one, with <code>event.reset()</code>. This routine will also reserve
the zeroth entry in the even record to represent the event as a whole. 

<p/>
With the <code>event.append(...)</code> methods a new entry is added at the
bottom of the current record, i.e. the first time it is called entry 
number 1 is filled, and so on. The <code>append</code>  method basically
exists in four variants, either without or with history information,
and with four-momentum provided either as a <code>Vec4</code> four-vector 
or as four individual components:
<pre>
  append( id, status, col, acol, p, m)
  append( id, status, col, acol, px, py, pz, e, m)
  append( id, status, mother1, mother2, daughter1, daughter2, col, acol, p, m)
  append( id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m)
</pre> 
The methods return the index at which the entry has been stored,
but normally you would not use this feature. 

<p/>
You can find descriptions of the input variables 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>here</a>. 
The PDG particle code <code>id</code> and the Les Houches Accord colour 
<code>col</code> and anticolour <code>acol</code> tags must be set
correctly. The four-momentum and mass have to be provided in units of GeV; 
if you omit the mass it defaults to 0. 

<p/>
Outgoing particles that should hadronize should be given status code 23.
Often this is the only status code you need. You could e.g. also fill in 
incoming partons with -21 and intermediate ones with -22, if you so wish.
Usually the choice of status codes is not crucial, so long as you recall 
that positive numbers correspond to particles that are still around, while
negative numbers denote ones that already hadronized or decayed. However,
so as not to run into contradictions with the internal PYTHIA checks
(when <code>Check:event = on</code>), or with external formats such as 
HepMC, we do recommend the above codes. When <code>pythia.next()</code> 
is called the positive-status particles that hadronize/decay get the sign 
of the status code flipped to negative but the absolute value is retained. 
The new particles are added with normal PYTHIA status codes.

<p/>
For normal hadronization/decays in <code>pythia.next()</code> the
history encoded in the mother and daughter indices is not used. 
Therefore the first two <code>append</code> methods, which set all these 
indices vanishing, should suffice. The subsequent hadronization/decays 
will still be properly documented.

<p/>
The exception is when you want to include junctions in your string 
topology, i.e. have three string pieces meet. Then you must insert in
your event record the (decayed) particle that is the reason for the
presence of a junction, e.g. a baryon beam remnant from which several
valence quarks have been kicked out, or a neutralino that underwent a
baryon-number-violating decay. This particle must have as daughters 
the three partons that together carry the baryon number. 
 
<p/>
The sample program in <code>main21.cc</code> illustrates how you can work
with this facility, both for simple parton configurations and for more
complicated ones with junctions.

<p/>
As an alternative to setting up a topology with the methods above, 
a <?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Event File</a> (LHEF) 
can also provide the configurations. Since no beams or processes are 
defined, only the <code>&lt;event&gt;....&lt;/event&gt;</code> blocks 
need to be present, one for each event, so strictly it is not a true 
LHEF. You need to select <code>Beams:frameType = 4</code>, provide 
the file name in <code>Beams:LHEF</code> and, as above, set
<code>ProcessLevel:all = off</code>. Needless to say, an externally 
linked <code>LHAup</code> class works as well as an LHEF,
with <code>Beams:frameType = 5</code>.

<p/>
The event information to store in the LHEF, or provide by the 
<code>LHAup</code>, is essentially the same as above. The only 
difference is in status codes: outgoing particles should have 1 
instead of 23, and intermediate resonances 2 instead of -22. 
Incoming partons, if any, are -1 instead of -21.

<h3>Extensions to resonance decays</h3>

With the above scheme, <code>pythia.next()</code> will generate
hadronization, i.e. string fragmentation and subsequent decays of
normal unstable particles. It will not decay 
<?php $filepath = $_GET["filepath"];
echo "<a href='ResonanceDecays.php?filepath=".$filepath."' target='page'>";?>resonances</a>, i.e.
<i>W, Z</i>, top, Higgs, SUSY and other massive particles.

<p/>
If the decay products are already provided, of course those 
products will be hadronized, but without any of the showers
that one would expect in <i>Z^0 -> q qbar</i>, say. That is, the
presence of a decayed resonance in the event record can be nice
for documentation purposes, but otherwise it plays no role. 
It is possible to change this behaviour with the following flag. 

<br/><br/><strong>Standalone:allowResDec</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off, then resonances are stable if final-state particles,
and irrelevant if intermediate particles. If on, then 
resonances are decayed, sequentially if necessary. After each
decay step a parton shower is applied, and subsequent decay 
kinematics is shifted by the recoil of such branchings. 
If the decay chain and kinematics of the resonance is provided
at input, this information is used, otherwise it is generated
internally according to built-in branching ratios etc.
  

<p/>
The input configuration has to follow the rules described above, 
with <code>ProcessLevel:all = off</code>. (Terminology not quite 
consistent, since resonance decays normally is part of the 
process-level step.) It is possible to combine several resonances, 
and other coloured or uncoloured particles into the same event.
 
<p/>
Final-state radiation is applied to each step of the decay 
sequence by default, but can be switched off with 
<code>PartonLevel:FSR = off</code> or 
<code>PartonLevel:FSRinResonances</code>. 

<h3>Repeated hadronization or decay</h3>

An alternative approach is possible with the 
<code>pythia.forceHadronLevel()</code> routine. This method does
a call to the <code>HadronLevel</code> methods, irrespective of the
value of the <code>HadronLevel:all</code> flag. If you hadronize
externally generated events it is equivalent to a 
<code>pythia.next()</code> call with 
<code>ProcessLevel:all = off</code>. 

<p/>
This method truly sticks to the hadron level, and thus cannot handle 
resonance decays. You therefore <b>must not</b> mix it with the
<code>Standalone:allowResDec = on</code> framework. 

<p/>
The similarity of names indicates that 
<code>pythia.forceTimeShower( int iBeg, int iEnd, double pTmax, 
int nBranchMax = 0)</code> is intended to belong to the same set of 
work-by-hand methods. Here <code>iBeg</code> and <code>iEnd</code>
give the range of partons that should be allowed to shower, 
<code>pTmax</code> the maximum <i>pT</i> scale of emissions,
and a nonzero <code>nBranchMax</code> a maximum number of allowed
branchings. Additionally, a <code>scale</code> has to be set for each
parton that should shower, which requires an additional final argument 
to the <code>append</code> methods above. This scale limits the maximum 
<i>pT</i> allowed for each parton, in addition to the global 
<code>pTmax</code>. When not set the scale defaults to 0, meaning no 
radiation for that parton.

<p/>
The real application instead is for repeated hadronization of the same
PYTHIA process- and parton-level event. This may for some studies
help to save time, given that these two first step are more 
time-consuming than the hadronization one. 

<p/>
For repeated hadronization you should first generate an event as usual, 
but with <code>HadronLevel:all = off</code>. This event you can save
in a temporary copy, e.g. <code>Event savedEvent = pythia.event</code>.
Inside a loop you copy back with <code>pythia.event = savedEvent</code>, 
and call <code>pythia.forceHadronLevel()</code> to obtain a new 
hadronization history.

<p/>
A more limited form of repetition is if you want to decay a given kind 
of particle repeatedly, without having to generate the rest of the event 
anew. This could be the case e.g. in <i>B</i> physics applications. 
Then you can use the <code>pythia.moreDecays()</code> method, which 
decays all particles in the event record that have not been decayed 
but should have been done so. The 
<code>pythia.particleData.mayDecay( id, false/true)</code> method 
may be used to switch off/on the decays of a particle species 
<code>id</code>, so that it is not decayed in the 
<code>pythia.next()</code> call but only inside a loop over a number of
tries. 

<p/>
Between each loop the newly produced decay products must be 
removed and the decayed particle status restored to undecayed.
The former is simple, since the new products are appended to the
end of the event record: <code>event.saveSize()</code> saves the
initial size of the event record, and <code>event.restoreSize()</code>
can later be used repeatedly to restore this original size, which means
that the new particles at the end are thrown away. The latter is more
complicated, and requires the user to identify the positions of all
particles of the species and restore a positive status code with
<code>event[i].statusPos()</code>.

<p/>
The <code>main15.cc</code> program illustrates both these methods,
i.e. either repeated hadronization or repeated decay of PYTHIA
events.

<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "off")
{
$data = "Standalone:allowResDec = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
