<html>
<head>
<title>Event Information</title>
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

<form method='post' action='EventInformation.php'>
 
<h2>Event Information</h2> 
<ol id="toc">
  <li><a href="#section0">List information</a></li>
  <li><a href="#section1">The beams</a></li>
  <li><a href="#section2">Initialization</a></li>
  <li><a href="#section3">The event type</a></li>
  <li><a href="#section4">Hard process initiators</a></li>
  <li><a href="#section5">Hard process parton densities and scales</a></li>
  <li><a href="#section6">Hard process kinematics</a></li>
  <li><a href="#section7">Soft Diffraction</a></li>
  <li><a href="#section8">Hard Diffraction</a></li>
  <li><a href="#section9">Photons from lepton beams</a></li>
  <li><a href="#section10">Event weight and activity</a></li>
  <li><a href="#section11">Multiparton interactions</a></li>
  <li><a href="#section12">Cross sections</a></li>
  <li><a href="#section13">Loop counters</a></li>
  <li><a href="#section14">Parton shower history</a></li>
  <li><a href="#section15">Les Houches Event File 3.0 information</a></li>
  <li><a href="#section16">Header information</a></li>
</ol>

 
The <code>Info</code> class collects various one-of-a-kind information, 
some relevant for all events and others for the current event. 
An object <code>info</code> is a public member of the <code>Pythia</code> 
class, so if you e.g. have declared <code>Pythia pythia</code>, the 
<code>Info</code> methods can be accessed by 
<code>pythia.info.method()</code>. Most of this is information that 
could also be obtained e.g. from the event record, but is here more 
directly available. It is primarily intended for processes generated 
internally in PYTHIA, but many of the methods would work also for 
events fed in via the Les Houches Accord. 
 
<a name="section0"></a> 
<h3>List information</h3> 
 
<a name="anchor1"></a>
<p/><strong> void Info::list() &nbsp;</strong> <br/>
a listing of most of the information set for the current event. 
   
 
<a name="section1"></a> 
<h3>The beams</h3> 
 
<a name="anchor2"></a>
<p/><strong> int Info::idA() &nbsp;</strong> <br/>
   
<a name="anchor3"></a>
<strong> int Info::idB() &nbsp;</strong> <br/>
the identities of the two beam particles. 
   
 
<a name="anchor4"></a>
<p/><strong> double Info::pzA() &nbsp;</strong> <br/>
   
<a name="anchor5"></a>
<strong> double Info::pzB() &nbsp;</strong> <br/>
the longitudinal momenta of the two beam particles. 
   
 
<a name="anchor6"></a>
<p/><strong> double Info::eA() &nbsp;</strong> <br/>
   
<a name="anchor7"></a>
<strong> double Info::eB() &nbsp;</strong> <br/>
the energies of the two beam particles. 
   
 
<a name="anchor8"></a>
<p/><strong> double Info::mA() &nbsp;</strong> <br/>
   
<a name="anchor9"></a>
<strong> double Info::mB() &nbsp;</strong> <br/>
the masses of the two beam particles. 
   
 
<a name="anchor10"></a>
<p/><strong> double Info::eCM() &nbsp;</strong> <br/>
   
<a name="anchor11"></a>
<strong> double Info::s() &nbsp;</strong> <br/>
the CM energy and its square for the two beams. 
   
 
<a name="section2"></a> 
<h3>Initialization</h3> 
 
<a name="anchor12"></a>
<p/><strong> bool Info::tooLowPTmin() &nbsp;</strong> <br/>
normally false, but true if the proposed <i>pTmin</i> scale was too low 
in timelike or spacelike showers, or in multiparton interactions. In the 
former case the <i>pTmin</i> is raised to some minimal value, in the 
latter the initialization fails (it is impossible to obtain a minijet 
cross section bigger than the nondiffractive one by reducing 
<i>pTmin</i>). 
   
 
<a name="section3"></a> 
<h3>The event type</h3> 
 
<a name="anchor13"></a>
<p/><strong> string Info::name() &nbsp;</strong> <br/>
   
<a name="anchor14"></a>
<strong> int Info::code() &nbsp;</strong> <br/>
the name and code of the process that occurred. 
   
 
<a name="anchor15"></a>
<p/><strong> int Info::nFinal() &nbsp;</strong> <br/>
the number of final-state partons in the hard process. 
   
 
<a name="anchor16"></a>
<p/><strong> bool Info::isResolved() &nbsp;</strong> <br/>
are beam particles resolved, i.e. were PDF's used for the process? 
   
 
<a name="anchor17"></a>
<p/><strong> bool Info::isDiffractiveA() &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<strong> bool Info::isDiffractiveB() &nbsp;</strong> <br/>
is either beam soft diffractively excited? 
   
 
<a name="anchor19"></a>
<p/><strong> bool Info::isDiffractiveC() &nbsp;</strong> <br/>
is there soft central diffraction (a.k.a. double Pomeron exchange)? 
   
 
<a name="anchor20"></a>
<p/><strong> bool Info::isHardDiffractiveA() &nbsp;</strong> <br/>
   
<a name="anchor21"></a>
<strong> bool Info::isHardDiffractiveB() &nbsp;</strong> <br/>
is either beam hard diffractively excited? 
   
 
<a name="anchor22"></a>
<p/><strong> bool Info::isNonDiffractive() &nbsp;</strong> <br/>
is the process the <code>SoftQCD:nonDiffractive</code> one, 
i.e. corresponding to the full inelastic nondiffractive part of the 
total cross section. (Note that a hard process, say <i>Z^0</i> 
production, normally is nondiffractive, but this is not what we 
aim at here, and so the method would return <code>false</code>, 
unless the <i>Z^0</i> had been generated as part of the MPI 
machinery for the <code>SoftQCD:nonDiffractive</code> component.) 
   
 
<a name="anchor23"></a>
<p/><strong> bool Info::isMinBias() &nbsp;</strong> <br/>
the same as above, retained for backwards compatibility, but to 
be removed in PYTHIA 8.2. 
   
 
<a name="anchor24"></a>
<p/><strong> bool Info::isLHA() &nbsp;</strong> <br/>
has the process been generated from external Les Houches Accord 
information? 
   
 
<a name="anchor25"></a>
<p/><strong> bool Info::atEndOfFile() &nbsp;</strong> <br/>
true if a linked Les Houches class refuses to return any further 
events, presumably because it has reached the end of the file from 
which events have been read in. 
   
 
<a name="anchor26"></a>
<p/><strong> bool Info::hasSub() &nbsp;</strong> <br/>
does the process have a subprocess classification? 
Currently only true for nondiffractive and Les Houches events, where 
it allows the hardest collision to be identified. 
   
 
<a name="anchor27"></a>
<p/><strong> string Info::nameSub() &nbsp;</strong> <br/>
   
<a name="anchor28"></a>
<strong> int Info::codeSub() &nbsp;</strong> <br/>
   
<a name="anchor29"></a>
<strong> int Info::nFinalSub() &nbsp;</strong> <br/>
the name, code and number of final-state partons in the subprocess 
that occurred when <code>hasSub()</code> is true. For a minimum-bias event 
the <code>code</code> would always be 101, while <code>codeSub()</code> 
would vary depending on the actual hardest interaction, e.g. 111 for 
<i>g g &rarr; g g</i>. For a Les Houches event the <code>code</code> would 
always be 9999, while <code>codeSub()</code> would be the external 
user-defined classification code. The methods below would also provide 
information for such particular subcollisions. 
   
 
<a name="section4"></a> 
<h3>Hard process initiators</h3> 
 
The methods in this sections refer to the two initial partons of the 
hard <i>2 &rarr; n</i> process (diffraction excluded; see below). 
 
<a name="anchor30"></a>
<p/><strong> int Info::id1() &nbsp;</strong> <br/>
   
<a name="anchor31"></a>
<strong> int Info::id2() &nbsp;</strong> <br/>
the identities of the two partons coming in to the hard process. 
   
 
<a name="anchor32"></a>
<p/><strong> double Info::x1() &nbsp;</strong> <br/>
   
<a name="anchor33"></a>
<strong> double Info::x2() &nbsp;</strong> <br/>
<i>x</i> fractions of the two partons coming in to the hard process. 
   
 
<a name="anchor34"></a>
<p/><strong> double Info::y() &nbsp;</strong> <br/>
   
<a name="anchor35"></a>
<strong> double Info::tau() &nbsp;</strong> <br/>
rapidity and scaled mass-squared of the hard-process subsystem, as 
defined by the above <i>x</i> values. 
   
 
<a name="anchor36"></a>
<p/><strong> bool Info::isValence1() &nbsp;</strong> <br/>
   
<a name="anchor37"></a>
<strong> bool Info::isValence2() &nbsp;</strong> <br/>
<code>true</code> if the two hard incoming partons have been picked 
to belong to the valence piece of the parton-density distribution, 
else <code>false</code>. Should be interpreted with caution. 
Information is not set if you switch off parton-level processing. 
   
 
<a name="section5"></a> 
<h3>Hard process parton densities and scales</h3> 
 
The methods in this section refer to the partons for which parton 
densities have been defined, in order to calculate the cross section 
of the hard process (diffraction excluded; see below). 
 
<p/> 
These partons would normally agree with the 
ones above, the initiators of the <i>2 &rarr; n</i> process, but it 
does not have to be so. Currently the one counterexample is POWHEG 
events [<a href="Bibliography.php#refAli10" target="page">Ali10</a>]. Here the original hard process could be 
<i>2 &rarr; (n-1)</i>. The NLO machinery at times would add an 
initial-state branching to give a <i>2 &rarr; n</i> process with a 
changed initial state. In that case the values in this section 
refer to the original <i>2 &rarr; (n-1)</i> state and the initiator 
ones above to the complete<i>2 &rarr; n</i> process. The 
<code>Info::list()</code> printout will contain a warning in such cases. 
 
<p/> 
For external events in the Les Houches format, the pdf information 
is obtained from the optional <code>#pdf</code> line. When this 
information is absent, the parton identities and <i>x</i> values agree 
with the initiator ones above, while the pdf values are unknown and 
therefore set to vanish. The <i>alpha_s</i> and <i>alpha_em</i> 
values are part of the compulsory information. The factorization and 
renormalization scales are both equated with the one compulsory scale 
value in the Les Houches standard, except when a <code>#pdf</code> 
line provides the factorization scale separately. If <i>alpha_s</i>, 
<i>alpha_em</i> or the compulsory scale value are negative at input 
then new values are defined as for internal processes. 
 
<a name="anchor38"></a>
<p/><strong> int Info::id1pdf() &nbsp;</strong> <br/>
   
<a name="anchor39"></a>
<strong> int Info::id2pdf() &nbsp;</strong> <br/>
the identities of the two partons for which parton density values 
are defined. 
   
 
<a name="anchor40"></a>
<p/><strong> double Info::x1pdf() &nbsp;</strong> <br/>
   
<a name="anchor41"></a>
<strong> double Info::x2pdf() &nbsp;</strong> <br/>
<i>x</i> fractions of the two partons for which parton density values 
are defined. 
   
 
<a name="anchor42"></a>
<p/><strong> double Info::pdf1() &nbsp;</strong> <br/>
   
<a name="anchor43"></a>
<strong> double Info::pdf2() &nbsp;</strong> <br/>
parton densities <i>x*f(x,Q^2)</i> evaluated for the two incoming 
partons; could be used e.g. for reweighting purposes in conjunction 
with the <code>idpdf</code>, <code>xpdf</code> and <code>QFac</code> 
methods. Events obtained from external programs or files may not 
contain this information and, if so, 0 is returned. 
   
 
<a name="anchor44"></a>
<p/><strong> double Info::QFac() &nbsp;</strong> <br/>
   
<a name="anchor45"></a>
<strong> double Info::Q2Fac() &nbsp;</strong> <br/>
the <i>Q</i> or <i>Q^2</i> factorization scale at which the 
densities were evaluated. 
   
 
<a name="anchor46"></a>
<p/><strong> double Info::alphaS() &nbsp;</strong> <br/>
   
<a name="anchor47"></a>
<strong> double Info::alphaEM() &nbsp;</strong> <br/>
the <i>alpha_strong</i> and <i>alpha_electromagnetic</i> values used 
for the hard process. 
   
 
<a name="anchor48"></a>
<p/><strong> double Info::QRen() &nbsp;</strong> <br/>
   
<a name="anchor49"></a>
<strong> double Info::Q2Ren() &nbsp;</strong> <br/>
the <i>Q</i> or <i>Q^2</i> renormalization scale at which 
<i>alpha_strong</i> and <i>alpha_electromagnetic</i> were evaluated. 
   
 
<a name="anchor50"></a>
<p/><strong> double Info::scalup() &nbsp;</strong> <br/>
returns the stored <code>SCALUP</code> value for Les Houches events, 
and else zero. It may agree with both the <code>QFac()</code> and 
<code>QRen()</code> values, as explained above. However, to repeat, 
should the input <code>SCALUP</code> scale be negative, separate positive 
factorization and renormalization scales are calculated and set as for 
internally generated events. Furthermore, when PDF info is supplied for 
the Les Houches event, the factorization scale is set by this PDF info 
(<code>scalePDF</code>), which can disagree with <code>SCALUP</code>. 
   
 
<a name="section6"></a> 
<h3>Hard process kinematics</h3> 
 
The methods in this section provide info on the kinematics of the hard 
processes, with special emphasis on <i>2 &rarr; 2</i> (diffraction excluded; 
see below). 
 
<a name="anchor51"></a>
<p/><strong> double Info::mHat() &nbsp;</strong> <br/>
   
<a name="anchor52"></a>
<strong> double Info::sHat() &nbsp;</strong> <br/>
the invariant mass and its square for the hard process. 
   
 
<a name="anchor53"></a>
<p/><strong> double Info::tHat() &nbsp;</strong> <br/>
   
<a name="anchor54"></a>
<strong> double Info::uHat() &nbsp;</strong> <br/>
the remaining two Mandelstam variables; only defined for <i>2 &rarr; 2</i> 
processes. 
   
 
<a name="anchor55"></a>
<p/><strong> double Info::pTHat() &nbsp;</strong> <br/>
   
<a name="anchor56"></a>
<strong> double Info::pT2Hat() &nbsp;</strong> <br/>
transverse momentum and its square in the rest frame of a <i>2 &rarr; 2</i> 
processes. 
   
 
<a name="anchor57"></a>
<p/><strong> double Info::m3Hat() &nbsp;</strong> <br/>
   
<a name="anchor58"></a>
<strong> double Info::m4Hat() &nbsp;</strong> <br/>
the masses of the two outgoing particles in a <i>2 &rarr; 2</i> processes. 
   
 
<a name="anchor59"></a>
<p/><strong> double Info::thetaHat() &nbsp;</strong> <br/>
   
<a name="anchor60"></a>
<strong> double Info::phiHat() &nbsp;</strong> <br/>
the polar and azimuthal scattering angles in the rest frame of 
a <i>2 &rarr; 2</i> process. 
   
 
<a name="section7"></a> 
<h3>Soft Diffraction</h3> 
 
Information on the primary elastic or 
<?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>diffractive</a> process 
(<i>A B &rarr; A B, X1 B, A X2, X1 X2, A X B</i>) can be obtained with 
the methods in the "Hard process kinematics" section above. The 
variables here obviously are <i>s, t, u, ...</i> rather than 
<i>sHat, tHat, uHat, ...</i>, but the method names remain to avoid 
unnecessary duplication. Most other methods are irrelevant for a 
primary elastic/diffractive process. 
 
<p/>Central diffraction <i>A B &rarr; A X B</i> is a <i>2 &rarr; 3</i> 
process, and therefore most of the <i>2 &rarr; 2</i> variables are 
no longer relevant. The <code>tHat()</code> and <code>uHat()</code> 
methods instead return the two <i>t</i> values at the <i>A &rarr; A</i> 
and <i>B &rarr; B</i> vertices, and <code>pTHat()</code> the average 
transverse momentum of the three outgoing "particles", while 
<code>thetaHat()</code> and <code>phiHat()</code> are undefined. 
 
<p/> 
While the primary interaction does not contain a hard process, 
the diffractive subsystems can contain them, but need not. 
Specifically, double diffraction can contain two separate hard 
subprocesses, which breaks the methods above. Most of them have been 
expanded with an optional argument to address properties of diffractive 
subsystems. This argument can take four values: 
<ul> 
<li>0 : default argument, used for normal nondiffractive events or 
the primary elastic/diffractive process (see above); </li> 
<li>1 : the <i>X1</i> system in single diffraction 
<i>A B &rarr; X1 B</i> or double diffraction <i>A B &rarr; X1 X2</i>; </li> 
<li>2 : the <i>X2</i> system in single diffraction 
<i>A B &rarr; A X2</i> or double diffraction <i>A B &rarr; X1 X2</i>; </li> 
<li>3 : the <i>X</i> system in central diffraction 
<i>A B &rarr; A X B</i>. </li> 
</ul> 
The argument is defined for all of the methods in the three sections above, 
"Hard process initiators", "Hard process parton densities and scales" and 
"Hard process kinematics", with the exception of the <code>isValence</code> 
methods. Also the four final methods of "The event type" section, the 
<code>...Sub()</code> methods, take this argument. But recall that they 
will only provide meaningful answers, firstly if there is a system of the 
requested type, and secondly if there is a hard subprocess in this system. 
A simple check for this is that <code>id1()</code> has to be nonvanishing. 
The methods below this section do not currently provide information 
specific to diffractive subsystems, e.g. the MPI information is not 
bookkept in such cases. 
 
<a name="section8"></a> 
<h3>Hard Diffraction</h3> 
 
Information on the momentum fraction taken from the beam 
and the momentum transfer in the hard diffractive process. 
Note that when side A is diffractively exited, then the Pomeron 
has been taken from side B and vice versa. 
 
<a name="anchor61"></a>
<p/><strong> double Info::xPomeronA() &nbsp;</strong> <br/>
   
<a name="anchor62"></a>
<strong> double Info::xPomeronB() &nbsp;</strong> <br/>
<i>x</i> fractions of momenta carried by the Pomeron in the hard 
diffractive process. 
   
<a name="anchor63"></a>
<p/><strong> double Info::tPomeronA() &nbsp;</strong> <br/>
   
<a name="anchor64"></a>
<strong> double Info::tPomeronB() &nbsp;</strong> <br/>
The momentum transfer <i>t</i> in the hard diffractive process. 
   
 
<a name="section9"></a> 
<h3>Photons from lepton beams</h3> 
 
Information about the kinematics of photon-photon collisions from lepton 
beams. 
<a name="anchor65"></a>
<p/><strong> double Info::eCMsub() &nbsp;</strong> <br/>
Collision energy of the <i>gamma-gamma</i> sub-system. 
   
<a name="anchor66"></a>
<p/><strong> double Info::xGammaA() &nbsp;</strong> <br/>
   
<a name="anchor67"></a>
<strong> double Info::xGammaB() &nbsp;</strong> <br/>
<i>x</i> fractions of lepton momenta carried by the photons. 
   
<a name="anchor68"></a>
<p/><strong> double Info::Q2GammaA() &nbsp;</strong> <br/>
   
<a name="anchor69"></a>
<strong> double Info::Q2GammaB() &nbsp;</strong> <br/>
Virtualities of the photons emitted by the leptons. 
   
<a name="anchor70"></a>
<p/><strong> double Info::thetaScatLepA() &nbsp;</strong> <br/>
   
<a name="anchor71"></a>
<strong> double Info::thetaScatLepB() &nbsp;</strong> <br/>
Scattering angles of the leptons wrt. the beam direction. 
   
<a name="anchor72"></a>
<p/><strong> int Info::photonMode() &nbsp;</strong> <br/>
Type of photon process, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='Photoproduction.php?filepath=".$filepath."' target='page'>";?>Photoproduction</a> for details. 
   
 
<a name="section10"></a> 
<h3>Event weight and activity</h3> 
 
<a name="anchor73"></a>
<p/><strong> double Info::weight() &nbsp;</strong> <br/>
weight assigned to the current event. Is normally 1 and thus 
uninteresting. However, there are several cases where one may have 
nontrivial event weights. These weights must the be used e.g. when 
filling histograms. 
<br/>(i) In the <code><?php $filepath = $_GET["filepath"];
echo "<a href='PhaseSpaceCuts.php?filepath=".$filepath."' target='page'>";?> 
PhaseSpace:increaseMaximum = off</a></code> default strategy, 
an event with a differential cross-section above the assumed one 
(in a given phase-space point) is assigned a weight correspondingly 
above unity. This should happen only very rarely, if at all, and so 
could normally be disregarded. 
<br/>(ii) The <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> class offers 
the possibility to bias the selection of phase space points, which 
means that events come with a compensating weight, stored here. 
<br/>(iii) For Les Houches events some strategies allow negative weights, 
which then after unweighting lead to events with weight -1. There are 
also Les Houches strategies where no unweighting is done, so events 
come with a weight. Specifically, for strategies <i>+4</i> and 
<i>-4</i>, the event weight is in units of pb. (Internally in mb, 
but converted at output.) 
   
 
<a name="anchor74"></a>
<p/><strong> double Info::weightSum() &nbsp;</strong> <br/>
Sum of weights accumulated during the run. For unweighted events this 
agrees with the number of generated events. In order to obtain 
histograms normalized "per event", at the end of a run, histogram 
contents should be divided by this weight. (And additionally 
divided by the bin width.) Normalization to cross section also 
required multiplication by <code>sigmaGen()</code> below. 
   
 
<a name="anchor75"></a>
<p/><strong> int Info::lhaStrategy() &nbsp;</strong> <br/>
normally 0, but if Les Houches events are input then it gives the 
event weighting strategy, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>. 
   
 
<a name="anchor76"></a>
<p/><strong> int Info::nISR() &nbsp;</strong> <br/>
   
<a name="anchor77"></a>
<strong> int Info::nFSRinProc() &nbsp;</strong> <br/>
   
<a name="anchor78"></a>
<strong> int Info::nFSRinRes() &nbsp;</strong> <br/>
the number of emissions in the initial-state showering, in the final-state 
showering excluding resonance decays, and in the final-state showering 
inside resonance decays, respectively. 
   
 
<a name="anchor79"></a>
<p/><strong> double Info::pTmaxMPI() &nbsp;</strong> <br/>
   
<a name="anchor80"></a>
<strong> double Info::pTmaxISR() &nbsp;</strong> <br/>
   
<a name="anchor81"></a>
<strong> double Info::pTmaxFSR() &nbsp;</strong> <br/>
Maximum <i>pT</i> scales set for MPI, ISR and FSR, given the 
process type and scale choice for the hard interactions. The actual 
evolution will run down from these scales. 
   
 
<a name="anchor82"></a>
<p/><strong> double Info::pTnow() &nbsp;</strong> <br/>
The current <i>pT</i> scale in the combined MPI, ISR and FSR evolution. 
Useful for classification in <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a>, 
but not once the event has been evolved. 
   
 
<a name="anchor83"></a>
<p/><strong> double Info::mergingWeight() &nbsp;</strong> <br/>
combined leading-order merging weight assigned to the current event, if 
tree-level multi-jet merging (i.e. 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?> CKKW-L</a> or 
<?php $filepath = $_GET["filepath"];
echo "<a href='UMEPSMerging.php?filepath=".$filepath."' target='page'>";?> UMEPS</a> merging) is attempted. 
If tree-level multi-jet merging is performed, all histograms should be 
filled with this weight, as discussed in 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?> CKKW-L Merging</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='UMEPSMerging.php?filepath=".$filepath."' target='page'>";?> UMEPS Merging</a>. 
   
 
<a name="anchor84"></a>
<p/><strong> double Info::mergingWeightNLO() &nbsp;</strong> <br/>
combined NLO merging weight assigned to the current event, if 
NLO multi-jet merging (i.e. 
<?php $filepath = $_GET["filepath"];
echo "<a href='NLOMerging.php?filepath=".$filepath."' target='page'>";?> NL<sup>3</sup></a> or 
<?php $filepath = $_GET["filepath"];
echo "<a href='NLOMerging.php?filepath=".$filepath."' target='page'>";?> UNLOPS</a> merging) is attempted. 
If NLO multi-jet merging is performed, all histograms should be filled 
with this weight, as discussed in 
<?php $filepath = $_GET["filepath"];
echo "<a href='NLOMerging.php?filepath=".$filepath."' target='page'>";?> NLO Merging</a>. 
   
 
<a name="section11"></a> 
<h3>Multiparton interactions</h3> 
 
As already noted, these methods do not make sense for diffractive 
topologies, and should not be used there. Partly this is physics, 
but mainly it is for technical reasons, e.g. that double diffraction 
involves two separate systems that would have to be bookkept as such. 
 
<a name="anchor85"></a>
<p/><strong> double Info::a0MPI() &nbsp;</strong> <br/>
The value of a0 when an x-dependent matter profile is used, 
<code>MultipartonInteractions:bProfile = 4</code>. 
   
 
<a name="anchor86"></a>
<p/><strong> double Info::bMPI() &nbsp;</strong> <br/>
The impact parameter <i>b</i> assumed for the current collision when 
multiparton interactions are simulated. Is not expressed in any physical 
size (like fm), but only rescaled so that the average should be unity 
for minimum-bias events (meaning less than that for events with hard 
processes). 
   
 
<a name="anchor87"></a>
<p/><strong> double Info::enhanceMPI() &nbsp;</strong> <br/>
The choice of impact parameter implies an enhancement or depletion of 
the rate of subsequent interactions, as given by this number. Again 
the average is normalized to be unity for minimum-bias events (meaning 
more than that for events with hard processes). 
   
 
<a name="anchor88"></a>
<p/><strong> double Info::enhanceMPIavg() &nbsp;</strong> <br/>
The average enhancement factor expected for hard processes, in those 
cases where it can be calculated already at initialization, i.e. excluding 
the <i>x</i>-dependent <i>b</i> profile. The normalization is here 
chosen to apply to cases with two hard interactions <i>A</i> and 
<i>B</i> preselected in the process level, and there multiplies 
<i>sigma_A * sigma_B / sigma_{nondiff}</i> to give the joint cross 
section. (Additional corrections from joint PDF weights somewhat reduce 
the final number.) The normalization is slightly different (typically 
around 5%) from the average of the <code>enhanceMPI()</code> method above, 
which instead is normalized to average value unity for nondiffractive events. 
As used internally the two are consistent. 
   
 
<a name="anchor89"></a>
<p/><strong> int Info::nMPI() &nbsp;</strong> <br/>
The number of hard interactions in the current event. Is 0 for elastic 
and diffractive events, and else at least 1, with more possible from 
multiparton interactions. 
   
 
<a name="anchor90"></a>
<p/><strong> int Info::codeMPI(int i) &nbsp;</strong> <br/>
   
<a name="anchor91"></a>
<strong> double Info::pTMPI(int i) &nbsp;</strong> <br/>
the process code and transverse momentum of the <code>i</code>'th 
subprocess, with <code>i</code> in the range from 0 to 
<code>nMPI() - 1</code>. The values for subprocess 0 is redundant with 
information already provided above. 
   
 
<a name="anchor92"></a>
<p/><strong> int Info::iAMPI(int i) &nbsp;</strong> <br/>
   
<a name="anchor93"></a>
<strong> int Info::iBMPI(int i) &nbsp;</strong> <br/>
are normally zero. However, if the <code>i</code>'th subprocess is 
a rescattering, i.e. either or both incoming partons come from the 
outgoing state of previous scatterings, they give the position in the 
event record of the outgoing-state parton that rescatters. 
<code>iAMPI</code> and <code>iBMPI</code> then denote partons coming from 
the first or second beam, respectively. 
   
 
<a name="anchor94"></a>
<p/><strong> double Info::eMPI(int i) &nbsp;</strong> <br/>
The enhancement or depletion of the rate of the <code>i</code>'th 
subprocess. Is primarily of interest for the 
<code>MultipartonInteractions:bProfile = 4</code> option, where the 
size of the proton depends on the <i>x</i> values of the colliding 
partons. Note that <code>eMPI(0) = enhanceMPI()</code>. 
   
 
<a name="anchor95"></a>
<p/><strong> double Info::bMPIold() &nbsp;</strong> <br/>
   
<a name="anchor96"></a>
<strong> double Info::enhanceMPIold() &nbsp;</strong> <br/>
   
<a name="anchor97"></a>
<strong> double Info::enhanceMPIoldavg() &nbsp;</strong> <br/>
These methods are only relevant for hard diffraction with the requirement 
of no MPI in the hadron-hadron collision. Then an impact parameter 
and associated enhancement factor is picked for this collision, but 
afterwards overwritten when the Pomeron-hadron subcollision is considered. 
In such cases the old hadron-hadron values can be found here, while 
<code>bMPI</code>, <code>enhanceMPI</code> and <code>enhanceMPIavg</code> 
provide the new Pomeron-hadron ones. 
   
 
<a name="section12"></a> 
<h3>Cross sections</h3> 
 
Here are the currently available methods related to the event sample 
as a whole, for the default value <code>i = 0</code>, and otherwise for 
the specific process code provided as argument. This is the number 
obtained with <code>Info::code()</code>, while the further subdivision 
given by <code>Info::codeSub()</code> is not bookkept. While continuously 
updated during the run, it is recommended only to study these properties 
at the end of the event generation, when the full statistics is available. 
The individual process results are not available if 
<?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>a second hard process</a> has been 
chosen, but can be gleaned from the <code>pythia.stat()</code> output. 
 
<a name="anchor98"></a>
<p/><strong> vector&lt;int&gt; Info::codesHard() &nbsp;</strong> <br/>
returns a vector with all the process codes set up for the current run, 
i.e. the valid nonzero arguments for the five methods below. 
   
 
<a name="anchor99"></a>
<p/><strong> string Info::nameProc(int i = 0) &nbsp;</strong> <br/>
returns the process name for process code <code>i</code>. 
   
 
<a name="anchor100"></a>
<p/><strong> long Info::nTried(int i = 0) &nbsp;</strong> <br/>
   
<a name="anchor101"></a>
<strong> long Info::nSelected(int i = 0) &nbsp;</strong> <br/>
   
<a name="anchor102"></a>
<strong> long Info::nAccepted(int i = 0) &nbsp;</strong> <br/>
the total number of tried phase-space points, selected hard processes 
and finally accepted events, summed over all allowed processes 
(<code>i = 0</code>) or for the given process. 
The first number is only intended for a study of the phase-space selection 
efficiency. The last two numbers usually only disagree if the user introduces 
some veto during the event-generation process; then the former is the number 
of acceptable events found by PYTHIA and the latter the number that also 
were approved by the user. If you set <?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>a 
second hard process</a> there may also be a mismatch. 
   
 
<a name="anchor103"></a>
<p/><strong> double Info::sigmaGen(int i = 0) &nbsp;</strong> <br/>
   
<a name="anchor104"></a>
<strong> double Info::sigmaErr(int i = 0) &nbsp;</strong> <br/>
the estimated cross section and its estimated error, 
summed over all allowed processes (<code>i = 0</code>) or for the given 
process, in units of mb. The numbers refer to the accepted event sample 
above, i.e. after any user veto. 
   
 
<a name="section13"></a> 
<h3>Loop counters</h3> 
 
Mainly for internal/debug purposes, a number of loop counters from 
various parts of the program are stored in the <code>Info</code> class, 
so that one can keep track of how the event generation is progressing. 
This may be especially useful in the context of the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a></code> facility. 
 
<a name="anchor105"></a>
<p/><strong> int Info::getCounter(int i) &nbsp;</strong> <br/>
the method that gives you access to the value of the various loop 
counters. 
<br/><code>argument</code><strong> i </strong>  :  the counter number you want to access: 
<br/><code>argumentoption </code><strong> 0 - 9</strong> :  counters that refer to the run as a whole, 
i.e. are set 0 at the beginning of the run and then only can increase. 
   
<br/><code>argumentoption </code><strong> 0</strong> :  the number of successful constructor calls for the 
<code>Pythia</code> class (can only be 0 or 1). 
   
<br/><code>argumentoption </code><strong> 1</strong> :  the number of times a <code>Pythia::init()</code> 
call has been begun. 
   
<br/><code>argumentoption </code><strong> 2</strong> :  the number of times a <code>Pythia::init()</code> 
call has been completed successfully. 
   
<br/><code>argumentoption </code><strong> 3</strong> :  the number of times a <code>Pythia::next()</code> 
call has been begun. 
   
<br/><code>argumentoption </code><strong> 4</strong> :  the number of times a <code>Pythia::next()</code> 
call has been completed successfully. 
   
<br/><code>argumentoption </code><strong> 10 - 19</strong> :  counters that refer to each individual event, 
and are reset and updated in the top-level <code>Pythia::next()</code> 
method. 
   
<br/><code>argumentoption </code><strong> 10</strong> :  the number of times the selection of a new hard 
process has been begun. Normally this should only happen once, unless a 
user veto is set to abort the current process and try a new one. 
   
<br/><code>argumentoption </code><strong> 11</strong> :  the number of times the selection of a new hard 
process has been completed successfully. 
   
<br/><code>argumentoption </code><strong> 12</strong> :  as 11, but additionally the process should 
survive any user veto and go on to the parton- and hadron-level stages. 
   
<br/><code>argumentoption </code><strong> 13</strong> :  as 11, but additionally the process should 
survive the parton- and hadron-level stage and any user cuts. 
   
<br/><code>argumentoption </code><strong> 14</strong> :  the number of times the loop over parton- and 
hadron-level processing has begun for a hard process. Is reset each 
time counter 12 above is reached. 
   
<br/><code>argumentoption </code><strong> 15</strong> :  the number of times the above loop has successfully 
completed the parton-level step. 
   
<br/><code>argumentoption </code><strong> 16</strong> :  the number of times the above loop has successfully 
completed the checks and user vetoes after the parton-level step. 
   
<br/><code>argumentoption </code><strong> 17</strong> :  the number of times the above loop has successfully 
completed the hadron-level step. 
   
<br/><code>argumentoption </code><strong> 18</strong> :  the number of times the above loop has successfully 
completed the checks and user vetoes after the hadron-level step. 
   
<br/><code>argumentoption </code><strong> 20 - 39</strong> :  counters that refer to a local part of the 
individual event, and are reset at the beginning of this part. 
   
<br/><code>argumentoption </code><strong> 20</strong> :  the current system being processed in 
<code>PartonLevel::next()</code>. Is almost always 1, but for double 
diffraction the two diffractive systems are 1 and 2, respectively. 
   
<br/><code>argumentoption </code><strong> 21</strong> :  the number of times the processing of the 
current system (see above) has begun. 
   
<br/><code>argumentoption </code><strong> 22</strong> :  the number of times a step has begun in the 
combined MPI/ISR/FSR evolution downwards in <i>pT</i> 
for the current system. 
   
<br/><code>argumentoption </code><strong> 23</strong> :  the number of times MPI has been selected for the 
downwards step above. 
   
<br/><code>argumentoption </code><strong> 24</strong> :  the number of times ISR has been selected for the 
downwards step above. 
   
<br/><code>argumentoption </code><strong> 25</strong> :  the number of times FSR has been selected for the 
downwards step above. 
   
<br/><code>argumentoption </code><strong> 26</strong> :   the number of times MPI has been accepted as the 
downwards step above, after the vetoes. 
   
<br/><code>argumentoption </code><strong> 27</strong> :   the number of times ISR has been accepted as the 
downwards step above, after the vetoes. 
   
<br/><code>argumentoption </code><strong> 28</strong> :   the number of times FSR has been accepted as the 
downwards step above, after the vetoes. 
   
<br/><code>argumentoption </code><strong> 29</strong> :  the number of times a step has begun in the 
separate (optional) FSR evolution downwards in <i>pT</i> 
for the current system. 
   
<br/><code>argumentoption </code><strong> 30</strong> :  the number of times FSR has been selected for the 
downwards step above. 
   
<br/><code>argumentoption </code><strong> 31</strong> :   the number of times FSR has been accepted as the 
downwards step above, after the vetoes. 
   
<br/><code>argumentoption </code><strong> 40</strong> :  keeps track of vetoed emission for shower 
reweighting. 
   
<br/><code>argumentoption </code><strong> 41 - 49</strong> :  counters that are unused (currently), and 
that therefore are free to use, with the help of the two methods below. 
   
   
   
 
<a name="anchor106"></a>
<p/><strong> void Info::setCounter(int i, int value = 0) &nbsp;</strong> <br/>
set the above counters to a given value. Only to be used by you 
for the unassigned counters 40 - 49. 
<br/><code>argument</code><strong> i </strong>  :  the counter number, see above. 
   
<br/><code>argument</code><strong> value </strong> (<code>default = <strong>0</strong></code>) :  set the counter to this number; 
normally the default value is what you want. 
   
   
 
<a name="anchor107"></a>
<p/><strong> void Info::addCounter(int i, int value = 0) &nbsp;</strong> <br/>
increase the above counters by a given amount. Only to be used by you 
for the unassigned counters 40 - 49. 
<br/><code>argument</code><strong> i </strong>  :  the counter number, see above. 
   
<br/><code>argument</code><strong> value </strong> (<code>default = <strong>1</strong></code>) :  increase the counter by this amount; 
normally the default value is what you want. 
   
   
 
<a name="section14"></a> 
<h3>Parton shower history</h3> 
 
The following methods are mainly intended for internal use, 
e.g. for matrix-element matching. 
 
<a name="anchor108"></a>
<p/><strong> void Info::hasHistory(bool hasHistoryIn) &nbsp;</strong> <br/>
   
<a name="anchor109"></a>
<strong> bool Info::hasHistory() &nbsp;</strong> <br/>
set/get knowledge whether the likely shower history of an event 
has been traced. 
   
 
<a name="anchor110"></a>
<p/><strong> void Info::zNowISR(bool zNowIn) &nbsp;</strong> <br/>
   
<a name="anchor111"></a>
<strong> double Info::zNowISR() &nbsp;</strong> <br/>
set/get value of <i>z</i> in latest ISR branching. 
   
 
<a name="anchor112"></a>
<p/><strong> void Info::pT2NowISR(bool pT2NowIn) &nbsp;</strong> <br/>
   
<a name="anchor113"></a>
<strong> double Info::pT2NowISR() &nbsp;</strong> <br/>
set/get value of <i>pT^2</i> in latest ISR branching. 
   
 
<a name="section15"></a> 
<h3>Les Houches Event File 3.0 information</h3> 
 
Les Houches Event files can conform to version 1.0 and version 3.0 
of the standard (version 2.0 having been extended to 3.0). The LHEF version 
of an input file can can be accessed by 
<a name="anchor114"></a>
<br/><strong> int Info::LHEFversion() &nbsp;</strong> <br/>
   
 
<p/> 
The <code>Info</code> class also provides a suitable interface to 
the information stored after reading Les Houches Event files in the 
updated format [<a href="Bibliography.php#refBut14" target="page">But14</a>]. An example main program using LHEF 3.0 
information is <code>main38.cc</code>. 
 
<p/> 
LHEF 3.0 offers new features both in the initialisation and the event sections 
of the input files. Possible information include extended 
use of XML tags in the <code>&lt;header&gt;</code> and 
<code>&lt;init&gt;</code> blocks.  The LHEF 3.0 information is stored 
in a series <code>struct</code>'s: 
 
<p/> 
 -- &nbsp; &nbsp; The <code>&lt;initrwgt&gt;</code> tag is a container 
 tag for weight and weightgroup tags.   This information is stored 
 internally in <code>LHAinitrwgt</code>. 
 Currently, there is no dedicated 
 output for this tag. However, all the information stored in the tag can 
 be retrieved by using the <code>Info</code> class member pointer 
 <code>LHAinitrwgt Info::initrwgt</code>. 
 
<p/> 
 -- &nbsp; &nbsp; Multiple <code>&lt;weightgroup&gt;</code> tags: 
 Container tag for weight tags. Currently, there is no dedicated 
 output for this tag. However, all the information stored in the tag can 
 be retrieved by using the <code>Info</code> class member pointer 
 <code>vector&lt;LHAweightgroups&gt; * Info::weightgroups</code>. 
 
<p/> 
 -- &nbsp; &nbsp; Multiple <code>&lt;weight&gt;</code> tags: Tag defining 
 auxiliary information on an event weight, e.g. the identifier and information 
 on what the weight represents. All the information stored in the tag can 
 be retrieved by using the <code>Info</code> class member pointer 
 <code>vector&lt;LHAweightgroups&gt; * Info::init_weights</code>. This vector 
 contains all <code>&lt;weight&gt;</code> tags in the 
 <code>&lt;initrwgt&gt;</code> container and its subcontainer 
 <code>&lt;weightgroup&gt;</code> tags. The size of the vector can be accessed 
 through the method 
<a name="anchor115"></a>
<br/><strong> int Info::getInitrwgtSize() &nbsp;</strong> <br/>
   
 
<p/> 
 -- &nbsp; &nbsp; Multiple <code>&lt;generator&gt;</code> tags: Store 
 information on the generators used in the event generation. All the 
 information stored in the tag can be retrieved by using the 
 <code>Info</code> class member pointer 
 <code>vector&lt;LHAgenerators&gt; * Info::generators</code>. More easy-to-use 
 output functions are available. The size of this vector can be obtained from 
<a name="anchor116"></a>
<br/><strong> int Info::getGeneratorSize() &nbsp;</strong> <br/>
   
 
<p/> 
The complete header can be obtained with the <code>Info</code> class 
member <code>string getHeaderBlock()</code>. 
 
<p/> 
The contents of a <code>&lt;generator&gt;</code> tag can be accessed through 
the method 
<a name="anchor117"></a>
<br/><strong> string Info::getGeneratorValue(unsigned int n = 0) &nbsp;</strong> <br/>
Return the contents of the n'th <code>&lt;generator&gt;</code> tag in 
the vector of tags. 
   
 
<p/> 
Attributes of the <code>&lt;generator&gt;</code> tag (e.g. the generator 
<code>name</code> and <code>version</code>) can be accessed via 
<a name="anchor118"></a>
<br/><strong> string Info::getGeneratorAttribute(unsigned int n, string key, bool doRemoveWhitespace = false) &nbsp;</strong> <br/>
Return the value of the generator attribute named <code>key</code> for 
the n'th generator in the vector. Setting <code>doRemoveWhitespace</code> to 
true will return the value, stripped of any whitespace. An empty string is 
returned if the attribute named <code>key</code> does not exist. 
   
 
<p/> 
To obtain information on cross sections, the following two methods can be 
used 
<a name="anchor119"></a>
<br/><strong> int Info::nProcessesLHEF() &nbsp;</strong> <br/>
return the number of processes for which the cross section is stored. 
   
<a name="anchor120"></a>
<br/><strong> double Info::sigmaLHEF(int iProcess) &nbsp;</strong> <br/>
return the cross section of the <code>iProcess</code>'th process. 
   
 
<p/> 
Possible information also includes extended use of XML tags in the 
<code>&lt;event&gt;</code> blocks: 
 
<p/> 
 -- &nbsp; &nbsp; The <code>&lt;rwgt&gt;</code> tag is a container 
 tag for wgt tags. Currently, there is no dedicated 
 output for this tag. It can however be retrieved by using the 
 <code>Info</code> class member pointer 
 <code>LHArwgt Info::rwgt</code>. 
 
<p/> 
 -- &nbsp; &nbsp; Multiple <code>&lt;wgt&gt;</code> tags: Tag defining 
 the event weight in the detailed version of LHEF 3.0.  All the information 
 stored in the tag can be retrieved by using the <code>Info</code> class 
 member pointer <code>vector&lt;LHAwgt&gt; * Info::weights_detailed</code>. 
 More easy-to-use output functions are available. The size of this vector 
 can be obtained from 
<a name="anchor121"></a>
<br/><strong> unsigned int Info::getWeightsDetailedSize() &nbsp;</strong> <br/>
   
 
<p/> 
A convenient access point to the information stored in the 
<code>&lt;wgt&gt;</code> tags is the <code>Info</code> class member 
<code>vector&lt;double&gt; Info::weights_detailed_vector</code>. The 
entries of this vector are ordered according to how <code>&lt;wgt&gt;</code> 
tags appear in the event block. 
 
<p/> 
The contents of a <code>&lt;wgt&gt;</code> tag can be accessed through the 
method 
<a name="anchor122"></a>
<br/><strong> double Info::getWeightsDetailedValue(string n) &nbsp;</strong> <br/>
Return the value of the n'th <code>&lt;wgt&gt;</code> tag in the 
event. 
   
 
<p/> 
Attributes of the <code>&lt;wgt&gt;</code> tag (e.g. the weight 
<code>id</code>) can be accessed via 
<a name="anchor123"></a>
<br/><strong> string Info::getWeightsDetailedAttribute(string n, string key, bool doRemoveWhitespace = false) &nbsp;</strong> <br/>
Return the value of the wgt attribute named <code>key</code> for 
the n'th wgt in the vector. Setting <code>doRemoveWhitespace</code> to 
true will return the value, stripped of any whitespace. An empty string is 
returned if the attribute named <code>key</code> does not exist. 
   
 
<p/> 
 -- &nbsp; &nbsp; The <code>&lt;weights&gt;</code> tag: Tag containing 
 a vector of <code>double</code> entries for weights in the compressed version 
 of LHEF 3.0. All the information stored in the tag can be retrieved by using 
 the <code>Info</code> class member pointer <code> LHAweights * 
 Info::weights</code> and the vector <code> vector&lt;double&gt; 
 Info::weights_compressed</code>. More easy-to-use output functions are 
 available. The size of this vector can be obtained from 
<a name="anchor124"></a>
<br/><strong> unsigned int Info::getWeightsCompressedSize() &nbsp;</strong> <br/>
   
 
<p/> 
The n'th weight can be accessed through the method 
<a name="anchor125"></a>
<br/><strong> double Info::getWeightsCompressedValue(unsigned int n) &nbsp;</strong> <br/>
   
 
<p/> 
Attributes of the <code>&lt;weights&gt;</code> tag (not normally used) can be 
accessed via 
<a name="anchor126"></a>
<br/><strong> string Info::getWeightsCompressedAttribute(string key, bool doRemoveWhitespace = false) &nbsp;</strong> <br/>
Return the value of the <code>&lt;weights&gt;</code> tag's attribute 
named <code>key</code>. Setting <code>doRemoveWhitespace</code> to 
true will return the value, stripped of any whitespace. An empty string is 
returned if the attribute named <code>key</code> does not exist. 
   
 
<p/> 
 -- &nbsp; &nbsp; The <code>&lt;scales&gt;</code> tag: Contains information 
 on different scales used by the matrix element generator. All the information 
 stored in the tag can be retrieved by using the <code>Info</code> class 
 member pointer <code> LHAweights * Info::scales</code>. More easy-to-use 
 output functions are available. The contents of the scales tag can be 
 obtained from 
<a name="anchor127"></a>
<br/><strong> string Info::getScalesValue() &nbsp;</strong> <br/>
   
 
<p/> 
However, note that the actual scale values are stored as attributes (called 
e.g. <code>muf</code> or <code>mur</code>). Attributes of the 
<code>&lt;scales&gt;</code> tag can be accessed via 
<a name="anchor128"></a>
<br/><strong> double Info::getScalesAttribute(string key) &nbsp;</strong> <br/>
Return the value of the <code>&lt;scales&gt;</code> tag's attribute 
named <code>key</code>. Not-a-number will be returned if the attribute 
named <code>key</code> does not exist. 
   
 
<p/> 
Finally, arbitrary attributes of the <code>&lt;event&gt;</code> tag are 
supported. Attributes of the <code>&lt;event&gt;</code> tag can be accessed by 
<a name="anchor129"></a>
<br/><strong> string Info::getEventAttribute(string key, bool doRemoveWhitespace = false) &nbsp;</strong> <br/>
return the value of the event attribute named <code>key</code>. Setting 
<code>doRemoveWhitespace</code> to true will return the value, stripped of 
any whitespace. An empty string is returned if the attribute named 
<code>key</code> does not exist. 
   
 
<p/> 
Additional comments appearing in the <code>&lt;event&gt;</code> tag 
can be obtained with the <code>Info</code> class member 
<code>string getEventComments()</code>. 
 
<a name="section16"></a> 
<h3>Header information</h3> 
 
A simple string key/value store, mainly intended for accessing 
information that is stored in the header block of Les Houches Event 
(LHE) files. In principle, any <code>LHAup</code> derived class can set 
this header information, which can then be read out later. Although the 
naming convention is arbitrary, in practice, it is dictated by the 
XML-like format of LHE files, see <?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?> 
Les Houches Accord</a> for  details. 
 
<a name="anchor130"></a>
<p/><strong> string Info::header(string key) &nbsp;</strong> <br/>
return the header named <code>key</code> 
   
 
<a name="anchor131"></a>
<p/><strong> vector &lt;string&gt; Info::headerKeys() &nbsp;</strong> <br/>
return a vector of all header key names 
   
 
<a name="anchor132"></a>
<p/><strong> void Info::setHeader(string key, string val) &nbsp;</strong> <br/>
set the header named <code>key</code> with the contents of <code>val</code> 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
