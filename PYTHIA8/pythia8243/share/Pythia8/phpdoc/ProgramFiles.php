<html>
<head>
<title>Program Files</title>
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

<form method='post' action='ProgramFiles.php'>
 
<h2>Program Files</h2> 
 
The code is subdivided into a set of files, mainly by physics 
task. Each file typically contains one main class, but often 
with a few related helper classes that are not used elsewhere in 
the program. Normally the files come in pairs. 
<ul> 
<li>A header file, <code>.h</code> in the <code>include/Pythia8</code> 
subdirectory, where the public interface of the class is declared, 
and inline methods are defined.</li> 
<li>A source code file, <code>.cc</code> in the <code>src</code> 
subdirectory, where most of the methods are implemented.</li> 
</ul> 
During compilation, related dependency files, <code>.d</code>, and 
compiled code, <code>.o</code> are created in the <code>tmp</code> 
subdirectory. 
 
<p/> 
In part the <code>.xml</code> documentation files in the 
<code>xmldoc</code> subdirectory have matching names, but the match 
is broken by the desire to group topics more by user interaction than 
internal operation. On these pages the function of the different code 
files is summarized. Currently, each <code>.xml</code> file is also 
translated into an <code>.html</code> one in the 
<code>htmldoc</code> subdirectory, to allow easy viewing of the 
contents in a web browser, and an <code>.php</code> one in 
<code>phpdoc</code>, for  sophisticated interactivity 
if that subdirectory is installed on a web server. 
 
<p/> 
Here is the current list of files, ordered alphabetically, with a brief 
description of contents. 
 
<a name="anchor1"></a>
<p/><code>file&nbsp; </code><strong> Analysis &nbsp;</strong> <br/>
contains routines to analyze events. Currently it can do sphericity, 
thrust, Lund/Jade/Durham jet clustering, cone-jet finding, and the 
<i>kT</i>, anti-<i>kT</i> and Cambridge/Aachen algorithms for 
hadron colliders. 
   
 
<a name="anchor2"></a>
<p/><code>file&nbsp; </code><strong> Basics &nbsp;</strong> <br/>
contains some basic facilities of general use: a random number 
generator <code>Rndm</code>, a four-vector class <code>Vec4</code>, and a 
histogram class <code>Hist</code>. 
   
 
<a name="anchor3"></a>
<p/><code>file&nbsp; </code><strong> BeamParticle &nbsp;</strong> <br/>
contains information on all partons extracted from one of the two 
beams. Defines modified parton distributions accordingly during the 
showering and multiparton interactions processing, thereby extending on 
the one-particle-inclusive distributions defined by the previous class. 
Finds the internal structure for a beam remnant. 
   
 
<a name="anchor4"></a>
<p/><code>file&nbsp; </code><strong> BeamRemnants &nbsp;</strong> <br/>
adds primordial <i>kT</i> to the set of hard subsystems, 
and combines these subsystems with the two beam remnants to provide 
the overall energy-momentum picture. Also ties together all the 
colour lines into consistent singlets. 
   
 
<a name="anchor5"></a>
<p/><code>file&nbsp; </code><strong> BeamShape &nbsp;</strong> <br/>
generates momentum spread of beams, and spread of collision vertex 
position. 
   
 
<a name="anchor6"></a>
<p/><code>file&nbsp; </code><strong> BoseEinstein &nbsp;</strong> <br/>
provides a simple method to impose Bose-Einstein correlations on 
pairs of identical mesons. 
   
 
<a name="anchor7"></a>
<p/><code>file&nbsp; </code><strong> ColourReconnection &nbsp;</strong> <br/>
implements several models that describe how the colour flow between 
outgoing partons can be changed by nonperturbative effects. 
   
 
<a name="anchor8"></a>
<p/><code>file&nbsp; </code><strong> ColourTracing &nbsp;</strong> <br/>
finds out how partons are colour-connected in (more or less long) chains. 
   
 
<a name="anchor9"></a>
<p/><code>file&nbsp; </code><strong> DeuteronProduction &nbsp;</strong> <br/>
provides a simple method to produce deuterons and anti-deuterons from 
events containing protons, neutrons, anti-protons, and anti-neutrons. 
   
 
<a name="anchor10"></a>
<p/><code>file&nbsp; </code><strong> Event &nbsp;</strong> <br/>
contains the event record, which basically is a vector of particles. 
This file also contains the <code>Particle</code> class, used by 
<code>Event</code>. <code>Pythia</code> uses two <code>Event</code> 
objects, one for the process-level record (<code>process</code>) and 
one for the complete event (<code>event</code>). 
   
 
<a name="anchor11"></a>
<p/><code>file&nbsp; </code><strong> FJcore &nbsp;</strong> <br/>
contains the external FJcore package, providing the coure functionality 
of the FastJet jet fiding package. 
   
 
<a name="anchor12"></a>
<p/><code>file&nbsp; </code><strong> FragmentationFlavZpT &nbsp;</strong> <br/>
contains the classes for describing the fragmentation steps in 
flavour and in longitudinal and transverse momentum. 
   
 
<a name="anchor13"></a>
<p/><code>file&nbsp; </code><strong> FragmentationSystems &nbsp;</strong> <br/>
defines some containers of parton systems, for use in 
the fragmentation routines. 
   
 
<a name="anchor14"></a>
<p/><code>file&nbsp; </code><strong> HadronLevel &nbsp;</strong> <br/>
turns the parton-level event above into a set of outgoing particles, 
by applying string fragmentation (with special treatment for low-mass 
systems) and secondary decays, and optionally Bose-Einstein corrections. 
   
 
<a name="anchor15"></a>
<p/><code>file&nbsp; </code><strong> HadronScatter &nbsp;</strong> <br/>
contains the beginning of a machinery for final state hadronic 
rescattering. 
   
 
<a name="anchor16"></a>
<p/><code>file&nbsp; </code><strong> HelicityBasics &nbsp;</strong> <br/>
basic classes for the handling of helicities in tau lepton decays. 
   
 
<a name="anchor17"></a>
<p/><code>file&nbsp; </code><strong> HelicityMatrixElements &nbsp;</strong> <br/>
helicity-dependent decay matrix elements for the tau lepton. 
   
 
<a name="anchor18"></a>
<p/><code>file&nbsp; </code><strong> HiddenValleyFragmentation &nbsp;</strong> <br/>
hadronization in models with a hidden sector that contains an unbroken 
<i>SU(N)</i>, which gives confinement. The model and code is largely 
derived from the normal fragmentation classes. 
   
 
<a name="anchor19"></a>
<p/><code>file&nbsp; </code><strong> History &nbsp;</strong> <br/>
methods to reconstruct the imagined shower history of a 
matrix-element-generated multiparton configuration, as part of 
the CKKW-L matrix element merging. 
   
 
<a name="anchor20"></a>
<p/><code>file&nbsp; </code><strong> Info &nbsp;</strong> <br/>
is a simple container that gives access to some information on the 
nature of the current process, such as Mandelstam variables. 
Also contains a small database for errors and warnings encountered 
during program execution. 
   
 
<a name="anchor21"></a>
<p/><code>file&nbsp; </code><strong> JunctionSplitting &nbsp;</strong> <br/>
processes a colour singlet system containing several (anti)junctions 
such that it can be split into several systems each containing 
(at most) one junction. 
   
 
<a name="anchor22"></a>
<p/><code>file&nbsp; </code><strong> LesHouches &nbsp;</strong> <br/>
gives the possibility to feed in parton configurations for the 
subsequent event generation. One base class is defined, with containers 
for initialization and event information, that can be read from 
<code>Pythia</code>. Derived classes allow for a few different cases. 
   
 
<a name="anchor23"></a>
<p/><code>file&nbsp; </code><strong> LHAFortran &nbsp;</strong> <br/>
is a header file only, for a class derived from the above LesHouches 
one, to be used for runtime interfacing to Fortran programs, such as 
PYTHIA 6. 
   
 
<a name="anchor24"></a>
<p/><code>file&nbsp; </code><strong> Merging &nbsp;</strong> <br/>
Wpapper class to interface matrix element merging schemes with Pythia. 
   
 
<a name="anchor25"></a>
<p/><code>file&nbsp; </code><strong> MergingHooks &nbsp;</strong> <br/>
intercede in the normal shower evolution to construct the relevant 
Sudakov form factor suppressions as part of the CKKW-L matrix element 
merging. 
   
 
<a name="anchor26"></a>
<p/><code>file&nbsp; </code><strong> MiniStringFragmentation &nbsp;</strong> <br/>
performs string fragmentation in cases where the colour singlet 
subsystem mass is so small that one or at most two primary hadrons 
should be produced from it. 
   
 
<a name="anchor27"></a>
<p/><code>file&nbsp; </code><strong> MultipartonInteractions &nbsp;</strong> <br/>
performs multiparton interactions. 
   
 
<a name="anchor28"></a>
<p/><code>file&nbsp; </code><strong> ParticleData &nbsp;</strong> <br/>
contains a database of all necessary particle data (masses, names, ..) 
and decay channels. 
   
 
<a name="anchor29"></a>
<p/><code>file&nbsp; </code><strong> ParticleDecays &nbsp;</strong> <br/>
performs the decays of all normal unstable hadrons and leptons, i.e. 
in mass up to and including <i>b bbar</i> systems. It is not 
intended for decays of electroweak resonances, like <i>Z^0</i>. 
   
 
<a name="anchor30"></a>
<p/><code>file&nbsp; </code><strong> PartonDistributions &nbsp;</strong> <br/>
contains parton distribution functions for the proton and electron. 
Currently very simple, with only two <i>p</i> parametrizations 
and one <i>e</i> ditto available, but it is possible to link in 
external sets. 
   
 
<a name="anchor31"></a>
<p/><code>file&nbsp; </code><strong> PartonLevel &nbsp;</strong> <br/>
turns the (hard) process above into a complete set of partons, by 
adding initial- and final-state radiation, multiparton 
interactions, and beam remnants. 
   
 
<a name="anchor32"></a>
<p/><code>file&nbsp; </code><strong> PartonSystems &nbsp;</strong> <br/>
keeps track of which partons belong to which partonic subsystem, 
i.e. one of the multiparton (semi)hard interactions with associated 
showers. 
   
 
<a name="anchor33"></a>
<p/><code>file&nbsp; </code><strong> PhaseSpace &nbsp;</strong> <br/>
selects a point in phase space for the hard-process generation, 
optimized separately for each process to give improved Monte Carlo 
efficiency. 
   
 
<a name="anchor34"></a>
<p/><code>file&nbsp; </code><strong> ProcessContainer &nbsp;</strong> <br/>
packages the information on a given subprocess, combining the 
phase-space selection and cross-section evaluation machineries 
with some statistics information. Also sets up the list of processes 
to be studied in a run. 
   
 
<a name="anchor35"></a>
<p/><code>file&nbsp; </code><strong> ProcessLevel &nbsp;</strong> <br/>
handles the generation of the (hard) process that sets the character 
of the event. This involves either using internally implemented 
processes or linking to Les Houches information. The latter can 
be by runtime interfaces or by reading in a file. This step also 
includes resonance decays. 
   
 
<a name="anchor36"></a>
<p/><code>file&nbsp; </code><strong> Pythia &nbsp;</strong> <br/>
is the main class, that administrates the whole event generation 
process by making use of all the others classes. Objects of most 
other classes reside (directly or indirectly) inside <code>Pythia</code>, 
so only a <code>Pythia</code> object needs to be explicitly instantiated 
and addressed   by the user. 
   
 
<a name="anchor37"></a>
<p/><code>file&nbsp; </code><strong> PythiaComplex &nbsp;</strong> <br/>
is only a <code>.h</code> file, containing a <code>typedef</code> for 
double precision complex numbers. 
   
 
<a name="anchor38"></a>
<p/><code>file&nbsp; </code><strong> PythiaStdlib &nbsp;</strong> <br/>
contains most of the <code>Stdlib</code> headers used in PYTHIA 8, 
with <code>using</code> directives. It defines <code>M_PI</code> if 
this is not already done. Also a few simple inline methods: 
<code>pow2(x)</code>, <code>pow3(x)</code>, <code>pow4(x)</code>, 
<code>pow5(x)</code> and <code>pow6(x)</code> for small integer 
powers, and <code>sqrtpos(x)</code> where a <code>max(0., x)</code> 
ensures that one does not take the square root of a negative number. 
The non-inlined <code>string toLower(const string, bool trim = true)</code> 
converts a strimg to lowercase, and by default also removes initial and 
trailing blanks or escape characters. It returns a new string, while the 
inlined extension <code>void toLowerRep(string, bool trim = true)</code> 
replaces the string in place. Finally the non-inlined 
<code>GammaReal(x)</code> gives the <i>Gamma</i> function value of 
a real argument. 
   
 
<a name="anchor39"></a>
<p/><code>file&nbsp; </code><strong> ResonanceDecays &nbsp;</strong> <br/>
decays resonances as part of the hard-process stage, in many cases 
(but not all) including angular correlations between the decay products. 
   
 
<a name="anchor40"></a>
<p/><code>file&nbsp; </code><strong> ResonanceWidths &nbsp;</strong> <br/>
encodes some properties of resonances, in particular the dynamic 
calculation of widths. 
   
 
<a name="anchor41"></a>
<p/><code>file&nbsp; </code><strong> RHadrons &nbsp;</strong> <br/>
handles the production and decay of hadrons formed by long-lived 
gluinos, stops or sbottoms. 
   
 
<a name="anchor42"></a>
<p/><code>file&nbsp; </code><strong> Settings &nbsp;</strong> <br/>
contains a database of all flags, modes, parameters and words that 
determine the performance of the generator. Initial values are obtained 
from the contents of the <code>.xml</code> files, but these values can 
then be changed by the user. 
   
 
<a name="anchor43"></a>
<p/><code>file&nbsp; </code><strong> SigmaCompositeness &nbsp;</strong> <br/>
contains the cross sections and matrix elements for production of 
some particles in compositeness scenarios, specifically excited 
fermions. 
   
 
<a name="anchor44"></a>
<p/><code>file&nbsp; </code><strong> SigmaEW &nbsp;</strong> <br/>
contains the cross sections and matrix elements for electroweak 
processes involving photons, <i>Z^0</i>'s and <i>W^+-</i>'s. 
   
 
<a name="anchor45"></a>
<p/><code>file&nbsp; </code><strong> SigmaExtraDim &nbsp;</strong> <br/>
contains the cross sections and matrix elements for processes in 
scenarios involving extra dimensions. 
   
 
<a name="anchor46"></a>
<p/><code>file&nbsp; </code><strong> SigmaGeneric &nbsp;</strong> <br/>
contains the cross sections and matrix elements for some generic 
processes, to be used as building blocks for a few BSM scenarios. 
   
 
<a name="anchor47"></a>
<p/><code>file&nbsp; </code><strong> SigmaHiggs &nbsp;</strong> <br/>
contains the cross sections and matrix elements for Higgs production. 
   
 
<a name="anchor48"></a>
<p/><code>file&nbsp; </code><strong> SigmaLeftRightSym &nbsp;</strong> <br/>
contains the cross sections and matrix elements for particle production 
in left-right-symmetry scenarios, specifically righthanded <i>Z</i> 
and <i>W</i> bosons and doubly-charged Higgs bosons. 
   
 
<a name="anchor49"></a>
<p/><code>file&nbsp; </code><strong> SigmaLeptoquark &nbsp;</strong> <br/>
contains the cross sections and matrix elements for leptoquark production. 
   
 
<a name="anchor50"></a>
<p/><code>file&nbsp; </code><strong> SigmaNewGaugeBosons &nbsp;</strong> <br/>
contains the cross sections and matrix elements for a <i>Z'^0</i>, 
a <i>W^+-</i> and a horizontal gauge boson <i>R^0</i>. 
   
 
<a name="anchor51"></a>
<p/><code>file&nbsp; </code><strong> SigmaOnia &nbsp;</strong> <br/>
contains the cross sections and matrix elements for charmonium and 
bottomonium production. 
   
 
<a name="anchor52"></a>
<p/><code>file&nbsp; </code><strong> SigmaProcess &nbsp;</strong> <br/>
contains the base class and derived classes for the evaluation of 
different matrix elements. Also keeps track of allowed incoming 
parton configurations and their cross sections, including parton 
densities. In order to keep this file from becoming too big, actual 
cross sections are found in several separate files of derived classes: 
<code>SigmaQCD</code>, <code>SigmaEW</code>, <code>SigmaOnia</code>, 
<code>SigmaHiggs</code>, <code>SigmaSUSY</code>, 
<code>SigmaNewGaugeBosons</code>, <code>SigmaLeftRightSym</code>, 
<code>SigmaLeptoquark</code>, <code>SigmaCompositeness</code>, 
<code>SigmaExtraDim</code> and <code>SigmaGeneric</code>. 
   
 
<a name="anchor53"></a>
<p/><code>file&nbsp; </code><strong> SigmaQCD &nbsp;</strong> <br/>
contains the cross sections and matrix elements for soft and hard 
QCD processes. 
   
 
<a name="anchor54"></a>
<p/><code>file&nbsp; </code><strong> SigmaSUSY &nbsp;</strong> <br/>
contains the cross sections and matrix elements for Supersymmetric 
processes. 
   
 
<a name="anchor55"></a>
<p/><code>file&nbsp; </code><strong> SigmaTotal &nbsp;</strong> <br/>
contains parametrizations of total, elastic and diffractive hadronic 
cross sections. 
   
 
<a name="anchor56"></a>
<p/><code>file&nbsp; </code><strong> SLHAinterface &nbsp;</strong> <br/>
handles the communication between the <code>SusyLesHouches</code> 
classes and Pythia. 
   
 
<a name="anchor57"></a>
<p/><code>file&nbsp; </code><strong> SpaceShower &nbsp;</strong> <br/>
performs spacelike initial-state transverse-momentum-ordered 
shower evolution. 
   
 
<a name="anchor58"></a>
<p/><code>file&nbsp; </code><strong> StandardModel &nbsp;</strong> <br/>
contains the running <i>alpha_strong</i>, with <i>Lambda</i> 
matching at flavour thresholds, the running <i>alpha_em</i>, 
CKM mixing matrices, and a few other parameters such as 
<i>sin^2(theta_W)</i>. 
   
 
<a name="anchor59"></a>
<p/><code>file&nbsp; </code><strong> StringFragmentation &nbsp;</strong> <br/>
performs string fragmentation of a given set of partons. 
   
 
<a name="anchor60"></a>
<p/><code>file&nbsp; </code><strong> StringLength &nbsp;</strong> <br/>
calculates the string-length <i>lambda</i> measure for various 
colour topologies, including junction ones. 
   
 
<a name="anchor61"></a>
<p/><code>file&nbsp; </code><strong> SusyCouplings &nbsp;</strong> <br/>
stores the various couplings used for SUSY cross sections and 
decays, as calculated from input e.g. based on the SUSY Les Houches 
Accord. 
   
 
<a name="anchor62"></a>
<p/><code>file&nbsp; </code><strong> SusyLesHouches &nbsp;</strong> <br/>
contains information on SUSY parameters and particle data as specified 
by the SUSY Les Houches Accord. 
   
 
<a name="anchor63"></a>
<p/><code>file&nbsp; </code><strong> SusyResonanceWidths &nbsp;</strong> <br/>
encodes some properties of supersymmetric resonances, similar to 
what <code>ResonanceWidths</code> does for other resonances. 
   
 
<a name="anchor64"></a>
<p/><code>file&nbsp; </code><strong> TauDecays &nbsp;</strong> <br/>
the main routines for handling tau lepton decays with helicity 
information. 
   
 
<a name="anchor65"></a>
<p/><code>file&nbsp; </code><strong> TimeShower &nbsp;</strong> <br/>
performs timelike final-state transverse-momentum-ordered 
shower evolution. 
   
 
<a name="anchor66"></a>
<p/><code>file&nbsp; </code><strong> UserHooks &nbsp;</strong> <br/>
Provides a way for a user to study the event at a few intermediate 
stages of evolution, to reject the event as a whole or to modify 
its cross-section weight. 
   
 
<a name="anchor67"></a>
<p/><code>file&nbsp; </code><strong> WeakShowerMEs &nbsp;</strong> <br/>
provides some matrix elements used for ME corrections of <i>W</i> 
and <i>Z</i> weak gauge boson emission in both initial- and final-state 
parton showers. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
