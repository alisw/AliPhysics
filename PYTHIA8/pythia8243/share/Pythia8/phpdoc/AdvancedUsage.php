<html>
<head>
<title>Advanced Usage</title>
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

<form method='post' action='AdvancedUsage.php'>
 
<h2>Advanced Usage</h2> 
<ol id="toc">
  <li><a href="#section0">The subsystems</a></li>
  <li><a href="#section1">The beams</a></li>
</ol>

 
On this page we collect information on a number of classes that 
the normal user would not encounter. There are cases where the 
information is essential, however, for instance to 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>implement your own showers</a>. 
 
<a name="section0"></a> 
<h3>The subsystems</h3> 
 
One aspect that complicates administration is that an event 
can contain several subsystems, each consisting of either 1) one MPI 
and its associated ISR and FSR or 2) one decaying resonance and its 
associated FSR. To first approximation these systems are 
assumed to evolve independently, but to second they are connected by 
the interleaved evolution, and potentially by colour-reconnection or 
rescattering effects. 
The partons of a given subsystem therefore do not have to be stored 
consecutively. 
 
<p/> 
The <code>PartonSystems</code> class is primarily used to keep track 
of the current positions of all partons belonging to each system, 
represented by the index <code>iPos</code> for a parton stored in the 
event-record slot <code>event[iPos]</code>. With "all" we mean the 
current set of outgoing partons, as well as the currently defined two 
incoming partons that system (for 2&rarr;n processes) or 
one incoming parton in the case of a decay (1&rarr;n) process. No 
intermediate-state (off-shell) ISR or FSR 
partons are present. That is, the parton system stores all partons 
that could be subject to some action in the next step of the 
combined MPI/ISR/FSR/BR description. As a special case, an outgoing 
parton is stored even if it undergoes a rescattering, and thus no 
longer belongs to the final state proper. 
 
<p/> 
Note also that an unstable (decaying) resonance will normally appear 
in two different systems; once, as an outgoing parton in the system 
that produced it (a hard process or the decay system of a previous 
resonance decay), and once as an incoming parton in its own decay system. 
 
<p/> 
The <code>partonSystems</code> instance of <code>PartonSystems</code> 
class is a public member of the <code>Pythia</code> top-level class, 
but is also available as a pointer <code>partonSystemsPtr</code> in 
various <code>PartonLevel</code> classes, e.g. inside the current 
instances of the <code>TimeShower</code> and <code>SpaceShower</code> 
classes. 
 
<p/> 
A number of <code>PartonSystems</code> methods can be used to set or 
get information on the subsystems: 
<ul> 
<li><code>clear()</code> resets all the contents in preparation for the 
next event.</li> 
<li><code>addSys()</code> add a new (initially empty) subsystem to the 
current list and return its index <code>iSys</code> in the list, 
where index 0 is the hardest subcollision and so on.</li> 
<li><code>sizeSys()</code> the number of separate subsystems.</li> 
<li><code>setInA(iSys, iPos), setInB(iSys, iPos)</code> store position 
<code>iPos</code> of the incoming parton from beam A or beam B to the 
<code>iSys</code>'th subcollision. These values are 0 initially, and 
should so remain if there are no beams, such as in resonance decays. 
</li> 
<li><code>setInRes(iSys, iPos)</code> stores position 
<code>iPos</code> of the incoming (decaying) resonance whose 
decay produced the outgoing partons for the  <code>iSys</code>'th 
system. This value is 0 initially and should so remain for systems 
that are not produced by the decay of a resonance, such as 
2&rarr;n subcollision systems. </li> 
<li><code>addOut(iSys, iPos)</code> store position <code>iPos</code> 
of a new outgoing parton in the <code>iSys</code>'th subcollision, 
by appending it at the end of the current vector, with beginning in 
slot 0.</li> 
<li><code>setOut(iSys, iMem, iPos)</code> store position <code>iPos</code> 
in the <code>iMem</code>'th slot in the vector of outgoing partons in the 
<code>iSys</code>'th subcollision. Here <code>iMem</code> must be in 
the range already constructed by <code>addOut</code> calls. </li> 
<li><code>replace(iSys, iPosOld, iPosNew)</code> replace the existing 
incoming or outgoing parton position <code>iPosOld</code> by 
<code>iPosNew</code> in the <code>iSys</code>'th subcollision.</li> 
<li><code>setSHat(iSys, sHat)</code> set the invariant squared mass 
<code>sHat</code> of the <code>iSys</code>'th subcollision.</li> 
<li><code>hasInAB(iSys)</code> true if an incoming parton has been set 
for beam A or beam B (and hence should have been set for both) in the 
<code>iSys</code>'th subcollision, else false.</li> 
<li><code>getInA(iSys), getInB(iSys)</code> the position <code>iPos</code> 
of the incoming parton from beam A or beam B to the <code>iSys</code>'th 
subcollision.</li> 
<li><code>hasInRes(iSys)</code> true if an incoming (decaying) 
resonance has been set for the <code>iSys</code>'th parton system, 
else false.</li> 
<li><code>sizeOut(iSys)</code> the number of outgoing partons 
in the <code>iSys</code>'th subcollision.</li> 
<li><code>getOut(iSys, iMem)</code> the position <code>iPos</code> 
of an outgoing parton in the <code>iSys</code>'th subcollision, 
with the <code>iMem</code> range limited by <code>sizeOut(iSys)</code>. 
These partons are not guaranteed to appear in any particular order. </li> 
<li><code>sizeAll(iSys)</code> the total number of incoming and outgoing 
partons in the <code>iSys</code>'th subcollision.</li> 
<li><code>getAll(iSys, iMem)</code> the position <code>iPos</code> 
of an incoming or outgoing parton in the <code>iSys</code>'th subcollision. 
In case there are beams it gives same as <code>getInA(iSys) </code> and 
<code> getInB(iSys)</code> for indices 0 and 1, and thereafter agrees with 
<code>getOut(iSys, iMem)</code> offset two positions. In case there is 
an incoming (decaying) resonance set for the system, it gives the same 
as <code>getInRes(iSys)</code> for index 0, and thereafter agrees with 
<code>getOut(iSys, iMem)</code> offset one position. If there are 
neither beams nor an incoming resonance set for the system, 
it is identical with <code>getOut(iSys, iMem)</code>.</li> 
<li><code>getSystemOf(iPos, alsoIn)</code> returns the system 
(<code>iSys</code>) of the parton specified by <code>iPos</code>. 
If the parton is outgoing in one system and incoming in another (eg a 
decaying resonance), the system in which it is incoming will be returned if 
<code>alsoIn == true</code>, else the system in which it is outgoing 
will be returned. The default is <code>alsoIn = false</code>. 
</li> 
<li><code>getSHat(iSys)</code> the invariant squared mass 
<code>sHat</code> of the <code>iSys</code>'th subcollision.</li> 
<li><code>list()</code> print a listing of all the system information, 
except for the <code>sHat</code> values.</li> 
</ul> 
 
<p/> 
New systems are created from the hard process, from resonance decays, and 
by the MPI, not from  any of the other components. Both FSR and ISR 
modify the position of partons, however. Since an FSR or ISR branching 
typically implies a new state with one more parton than before, an 
outgoing parton must be added to the system. Furthermore, in a 
branching, several existing partons may also be moved to new slots, 
including the incoming beam ones. 
In a FSR <i>1 &rarr; 2</i> branching it is irrelevant which parton position 
you let overwrite the existing slot and which is added to the end of 
the system. 
 
<p/> 
The system information must be kept up-to-date. Both the MPI, ISR, FSR and 
BR descriptions make extensive use of the existing information. As an 
example, the introduction of primordial <i>kT</i> in the beam remnants 
will fail if the information on which final-state partons belong to which 
system is out-of-date. The introduction of rescattering as part of the 
MPI framework adds further complications, where an outgoing parton of one 
subsystem may be the incoming one of another system. This part of the code 
is still under development. 
 
<p/> 
Currently the system information is kept throughout the continued 
history of the event. Specifically, resonance decays create new systems, 
appended to the existing ones. This could be useful during the 
hadronization stage, to collect the partons that belong to a resonance 
with preserved mass when a small string collapses to one particle, 
but is not yet used for that. 
 
<a name="section1"></a> 
<h3>The beams</h3> 
 
The different subsystems are tied together by them sharing the same 
initial beam particles, and thereby being restricted by energy-momentum 
and flavour conservation issues. The information stored in the two 
beam particles, here called <code>beamA</code> and <code>beamB</code>, 
is therefore as crucial to keep correct as the above subsystem list. 
 
<p/> 
Both beam objects are of the <code>BeamParticle</code> class. 
Each such object contains a vector with the partons extracted from it. 
The number of such partons, <code>beamX.size()</code> (X = A or B), 
of course is the same as the above number of subsystems in the event 
record. (The two diverge at the BR step, where further beam remnants 
are added to the beams without corresponding to new subsystems.) 
The individual partons are accessed by an overloaded indexing 
operator to a vector of <code>ResolvedParton</code> objects. The 
<code>iPos()</code> property corresponds to the <code>iPos</code> 
one above, i.e. providing the position in the main event record of 
a parton. In particular, 
<code>beamA[iSys].iPos() = partonSystemsPtr->getInA(iSys)</code> and 
<code>beamB[iSys].iPos() = partonSystemsPtr->getInB(iSys)</code>. 
Whereas thus the indices of the two incoming partons to a subsystem 
are stored in two places, the ones of the outgoing partons only 
appear in the system part of the <code>PartonSystems</code> class. 
 
<p/> 
Just as the subsystems in <code>PartonSystems</code> must be updated, 
so must the information in the two <code>BeamParticle</code>'s, e.g. 
with methods<code>beamX[iSys].iPos( iPosIn)</code> when an incoming 
parton is replaced by a new one in line <code>iPosIn</code>. Further 
the new parton identity should be set by <code>beamX[iSys].id( idIn)</code> 
and the new <i>x</i> energy-momentum fraction by 
<code>beamX[iSys].x( xIn)</code>. The three can be combined in one go 
by <code>beamX[iSys].update( iPosIn, idIn, xIn)</code>. 
 
<p/> 
To be specific, it is assumed that, at each step, the two incoming 
partons are moving along the <i>+-z</i> axis and are massless. 
Since the event is constructed in the c.m. frame of the incoming 
beams this implies that <i>x = 2 E / E_cm</i>. 
If the <i>x</i> values are not defined accordingly or not kept 
up-to-date the BR treatment will not conserve energy-momentum. 
 
<p/> 
In return, the <code>BeamParticle</code> objects give access to some 
useful methods. The <code>beamX.xf( id, x, Q2)</code> returns the 
standard PDF weight <i>x f_id(x, Q^2)</i>. More interestingly, 
<code>beamX.xfISR( iSys, id, x, Q2)</code> returns the modified weight 
required when several subsystems have to share the energy and flavours. 
Thus <code>iSys</code> is added as an extra argument, and the momentum 
already assigned to the other subsystems is not available for evolution, 
i.e. the maximal <i>x</i> is correspondingly smaller than unity. 
Also flavour issues are handled in a similar spirit. 
 
<p/> 
An additional complication is that a parton can be either valence or 
sea, and in the latter case the BR treatment also distinguishes 
companion quarks, i.e. quark-antiquark pairs that are assumed to 
come from the same original <i>g &rarr; q qbar</i> branching, whether 
perturbative or not. This can be queried either with the 
<code>beamX[iSys].companion()</code> method, for detailed information, 
or with the <code>beamX[iSys].isValence()</code>, 
<code>beamX[iSys].isUnmatched()</code> and 
<code>beamX[iSys].isCompanion()</code> methods for yes/no answers 
whether a parton is valence, unmatched sea or matched sea. 
This choice should affect the ISR evolution; e.g., a valence quark 
cannot be constructed back to come from a gluon. 
 
<p/> 
To keep this info up-to-date, the <code>beamX.pickValSeaComp()</code> 
method should be called whenever a parton of a new flavour has been 
picked in the ISR backwards evolution, but not if the flavour has not 
been changed, since then one should not be allowed to switch back and 
forth between the same quark being considered as valence or as sea. 
Since the <code>pickValSeaComp()</code> method makes use of the 
current parton-density values, it should be preceded by a call 
to <code>beamX.xfISR( iSys, id, x, Q2)</code>, where the values in 
the call are the now finally accepted ones for the newly-found mother. 
(Such a call is likely to have been made before, during the evolution, 
but is not likely to be the most recent one, i.e. still in memory, and 
therefore had better be redone.) 
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
