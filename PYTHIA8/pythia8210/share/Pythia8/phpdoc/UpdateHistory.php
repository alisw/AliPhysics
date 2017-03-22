<html>
<head>
<title>Update History</title>
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

<form method='post' action='UpdateHistory.php'>
 
<h2>Update History</h2> 
 
These update notes describe major updates relative to the 
PYTHIA 8.186 version, which was the last regular 8.1 release. 
(Minor bug fixes will continue to appear.) The step from 
8.1 to 8.2 gave an occasion to break backwards compatibility, 
but this should only affect a small part of the user code. 
 
<h3>Main news by version</h3> 
 
<ul> 
 
<li>8.210: 29 June 2015 
<ul> 

<li>Bug fix in CKKW-L merging for LHE files, such that the factorization
  and renormalization scales are set by the <code>SCALUP</code> value if 
the <code>muf2</code> and  <code>mur2</code> LHEF3 attributes have not
been set, and the user has not set any explicit values. This change 
restores the PYTHIA 8.1 behaviour.</li> 
 
<li>Various technical improvements in the machinery for hard 
diffraction.</li> 

<li>Correct quark flavour selection when a string spanned directly 
between two junctions is split up.</li> 
 
<li>Check that SK-I and SK-II colour reconnection machineries only
are called for event topologies they are set up to handle.</li> 

<li>Bug fixes in partial widths of the <i>W'</i> boson. Results are
correct when the <i>W'</i> is a simply rescaled copy of the <i>W</i>, 
but not for more general couplings. Thanks to Mihail Chizhov. </li> 
 
<li>Minor fix in default location of PDF data files in the constructors.
No practical consequence since correct non-default values are used.</li> 

<li>Tiny fix in the <code>configure</code> script, so that CXX options
containing an equal sign are parsed correctly.</li> 
 
</ul> 
</li> 
 
<li>8.209: 25 May 2015 
<ul> 
 
<li>The "An Introduction to PYTHIA 8.2", arXiv:1410.3012 [hep-ph], 
has now been published in Comput. Phys.Commun. 191 (2015) 159.</li> 
 
<li>The new QCD-based colour reconnection machinery has been expanded 
with checks whether reconnecting dipoles are causally connected. 
This cannot be defined exactly, and therefore several different 
options are available. Also some other minor changes in this part 
of the code, including updates of parameter default values.</li> 
 
<li>New options for <code>ColourReconnection:flipMode</code> allows 
to sidestep the gluon-move handling but still retain the final 
flip step.</li> 
 
<li>New <code>ColourReconnection:forceResonance</code> switch 
allows an additional colour reconnection step after late resonance 
decays. This is especially relevant for 
<i>H^0 &rarr; W^+ W^- / Z^0 Z^0  &rarr; q_1 qbar_2 q_3 qbar_4</i>, 
since the Higgs is so long-lived that its decay is well separated 
from the rest of the event.</li> 
 
<li>The old SK I and SK II colour reconnection models are now available. 
These are specially aimed at the processes 
<i>e^+ e^- &rarr; W^+ W^- / Z^0 Z^0</i>, but are also relevant for 
<i>H^0 &rarr; W^+ W^- / Z^0 Z^0</i>. </li> 
 
<li>Colour reconnection also made possible, optionally, for the 
<code>forceHadronLevel</code> method.</li> 
 
<li>A new switch <code>POWHEG:QEDveto</code> has been introduced to 
steer the treatment of non-QCD radiation in the POWHEG implementation 
of <code>include/Pythia8Plugins/PowhegHooks.h</code>.</li> 
 
<li>The need to link to the Boost library to read gzipped files 
has been eliminated by including iostream classes wrapping the zlib 
compression library, see new files <code>Streams.h</code> and 
<code>Streams.cc</code>.</li> 
 
<li>Linkage to LHAPDF6 so far has been based on the LHAPDF5 compatibility 
mode, which requires no Boost headers. The new configure option 
<code>--with-lhapdf6-plugin=LHAPDF6.h</code> uses native mode, and 
then requires Boost headers.</li> 
 
<li>The LHAPDF6 PDF values are frozen at the respective <i>x</i> 
and <i>Q^2</i> boundary if the <i>(x, Q^2)</i> pair falls outside 
the fit region.</li> 
 
<li>New functions to extract the fit boundaries, the current 
<i>alpha_s</i> value, and the quark masses from LHADPF6.</li> 
 
<li>New option in the HepMC interface, whereby PYTHIA particles can be 
appended to an existing HepMC event. Thanks to Mikhail Kirsanov.</li> 
 
<li>A new runtime interface to the POWHEGBOX matrix element programs, 
bypassing the need for intermediate LHE files. The new files 
<code>include/Pythia8Plugins/LHAPowheg.h</code> and 
<code>include/Pythia8Plugins/PowhegProcs.h</code> contain the 
LHAup class wrapper used to build the POWHEG plugin libraries and 
the simple class that facilitates loading the POWHEG plugins, 
respectively. The new <code>examples/main33.cc</code> demonstrates 
how to use these plugins, and <code>examples/main33.cmnd</code> contains 
the commands needed for POWHEGBOX to run the example.</li> 
 
<li>When reading in particle data from the SLHA interface, changes done 
by the user takes precedence over the SLHA input ones. To be more 
specific, particle data changes by the <code>Pythia::readString</code> 
and <code>readFile</code> methods are buffered and repeated after 
the SLHA initialization.</li> 
 
<li>New <code>#define PYTHIA_VERSION 8.2xx</code> in <code>Pythia.h</code> 
allows user-code preprecessors to make version-specific choices, and 
allows the <code>Pythia</code> class constructor to check that the 
header-file version number matches those of the source code and the 
XML files. Thanks to Pere Mato.</li> 
 
<li>User-defined semi-internal processes can now be accompanied by 
user-defined phase-space generators, via a second optional argument to 
<code>Pythia::setSigmaPtr(SigmaProcess*, PhaseSpace* = 0)</code>. 
Default is the old behaviour, with PYTHIA selecting the phase-space 
generator itself, based on the process type. Alternatively, the user 
may provide a pointer to an instance of an object inheriting from 
the <code>PhaseSpace</code> class. Sufficiently well-tested and general 
such generators could be communicated to the PYTHIA authors for possible 
inclusion in a future release.</li> 
 
<li>Updated initialization of LHEF files in <code>LesHouches.cc</code> 
to ignore file contents enclosed in comment tags, 
<code>&lt;!-- .. --&gt;</code> and/or the special CDATA statement, 
<code>&lt;![CDATA[ .. ]]&gt;</code>. The latter occurs, e.g., in LHEF 
files produced by CalcHep. (It is generally used to store content that 
contains XML-illegal characters like "&lt;" or "&", such as JavaScript 
source.) Thanks to A. Belyaev and A. Pukhov for help with this update.</li> 
 
<li>New flag <code>LesHouches:matchInOut</code>, by default on, 
to recalculate the energies and longitudinal momenta of the incoming 
particles from the outgoing ones for Les Houches input. Reduces effect 
of numerical inconsistencies in input.</li> 
 
<li>Updated EvtGen interface, to allow a pointer to an FSR engine to be 
passed. Thanks to Torben Ferber.</li> 
 
<li>Spin information for <i>tau</i> leptons now also set up for 
<i>Z'</i> and <i>W'</i> decays.</li> 
 
<li>New method<code>AlphaStrong::setThresholds(...)</code> allows to 
set the charm, bottom and top flavour-threshold masses used for the 
running of <i>alpha_strong</i>.</li> 
 
<li>New option for the <code>Event::list()</code> methods allows to show 
momenta with more decimal digits.</li> 
 
<li>New <code>Particle::isFinalPartonLevel()</code> method to tell whether 
a particle belonged to the final state on the parton level of generation 
or not. New <code>main73.cc</code> example illustrates usefulness.</li> 
 
<li>The <code>[]</code> operator is implemented for the <code>Vec4</code> 
class to return its components by index.</li> 
 
<li>New method <code>string Settings::output(string key)</code> returns 
the value of a variable as a string. In <code>Pythia::readString()</code> 
or <code>readFile()</code> calls this can now be used to print a 
current setting value by the command <code>key = ?</code>.</li> 
 
<li>The role of the <code>pythia.forceTimeShower(...)</code> method is 
better explained in the hadron-level standalone documentation, and 
<code>main21.cc</code> has been extended with an example.</li> 
 
<li>The <code>TimeShower::enhancePTmax()</code> method is made virtual.</li> 
 
<li>The <code>iTopCopyId()</code> and <code>iBotCopyId()</code> methods 
now scan all mothers/daughters, rather than only the first and last, 
when searching for a unique flavour match.</li> 
 
<li>New example <code>main62.cc</code> illustrates how a user hook can 
steer the selection of angles in a resonance decay.</li> 
 
<li>Minor correction in the treatment of the highest multiplicity in 
FxFx jet matching.</li> 
 
<li>Minor correction in the interface to aMC@NLO for dijet and 
photon+jet processes.</li> 
 
<li>Inserted missing endtag that corrupted the <code>Tunes.php</code> 
page. Thanks to Tim Martin.</li> 
 
<li>Minor fix for polarization sign when the <i>tau</i> polarization 
is forced.</li> 
 
<li>The meaning of the <code>HiggsXX:parity</code> options for CP mixing 
has been slightly modified and is better described.</li> 
 
<li>Clarification in the documentation that impact-parameter-enhancement 
factor calculation for two hard processes does not work for the 
<i>x</i>-dependent impact-parameter profile option.</li> 
 
<li>Fixed a factor <i>sqrt(2)</i> error in couplings for 
chargino + squark pair production. Thanks to Michihisa Takeuchi.</li> 
 
<li>Some minor bug fixes in the SUSY code, and some speed optimization 
suggested by Martin White.</li> 
 
<li>Fixed incorrect indexing in the SLHA interface.</li> 
 
<li>Removed misleading flavour setup for unused Pomeron.</li> 
 
<li>Fixed possible unwarranted destruction of pointer to external timelike 
shower by the <code>Pythia</code> destructor.</li> 
 
<li>Initialized pointers to NULL in the <code>Info</code> class, 
to avoid some problems. Thanks to Keno Fischer.</li> 
 
<li>Minor <code>Makefile</code> improvements for LHAPDF linking.</li> 
 
<li>Minor improvement to the <code>main89.cc</code> example program.</li> 
 
<li>Mildly modified warning/error messages when junction splitting 
fails.</li> 
 
<li>Minor fixes in <code>LHAFortran.h</code>, which is also moved to 
<code>include/Pythia8Plugins</code> to better reflect its peripheral 
role.</li> 
 
<li>New machinery for hard diffraction now in place. Still being 
debugged and tested, so not yet ready for public usage.</li> 
 
</ul> 
</li> 
 
<li>8.205: 23 January 2015 
<ul> 
 
<li>Unfortunate tiny typo made the new A14 tunes inaccessible. 
Now fixed.</li> 
 
<li>Resolved an inconsistency between MLM and FxFx merging that was 
introduced in 8.201 upon answering a user request. The effect should 
be minimal. Thanks to Josh Bendavid for bringing this to our attention. 
</li> 
 
<li>Amended the automated reading of the (optional) 
"TimeShower:nPartonsInBorn" setting from Les Houches Event files 
produced with aMC@NLO. Heavy coloured objects are now handled correctly 
when the information in the Les Houches event does not include heavy 
partons of the lowest-multiplicity state into its counting (as is 
conventional). Thanks to Josh Bendavid for pointing this out.</li> 
 
</ul> 
</li> 
 
<li>8.204: 22 January 2015 
<ul> 
 
<li>The <code>examples</code> directory has been moved back from 
the <code>share/Pythia8/</code> subdirectory to the main directory, 
as was the case in PYTHIA 8.1, to make it more visible to newcomers. 
The optional <code>make install</code> step will create a copy of 
<code>examples</code> in <code>share/Pythia8/</code>. The rarely 
used <code>examples/outref</code> subdirectory is moved to 
<code>share/Pythia8/outref</code>.</li> 
 
<li>Fifteen new <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a> have been added, 
the MonashStar tune from CMS and fourteen A14 tunes from ATLAS 
[<a href="Bibliography.php" target="page">ATL14a</a>]. The latter correspond to central tunes for 
four different PDF sets and ten variations in five (approximate) 
eigenvector directions. Furthermore, now the chosen 
<code>Tune:pp</code> implies the <code>Tune:ee</code> value 
to which it is related, and thus the latter need not be set 
separately.</li> 
 
<li>Default settings values have been updated to agree with the 
Monash 2013 tune. Thus typically the list of changed settings is 
significantly reduced. Thanks to Mikhail Kirsanov for suggestion.</li> 
 
<li>The compositeness section has been expanded with six further 
processes, describing the pair production of excited leptons or 
neutrinos. Three-body contact-interaction decay modes of these 
excited states have been introduced. A bug has been fixed that 
gave the wrong helicity in decays for excited quarks, leptons and 
neutrinos. Further, the <i>gamma^*/Z^0/Z'^0</i> can decay to a pair 
of excited fermions provided that the channels are added to the 
list of allowed ones. Based on code provided by Olga Igonkina.</li> 
 
<li>A new interface to the EvtGen decay package, primarily intended 
for bottom and charm decays, has been implemented. It is available 
in <code>include/Pythia8Plugins/EvtGen.h</code> and an example how 
to set it up is found in <code>main48</code>.</li> 
 
<li>Several minor changes in the <code>examples</code> subdirectory. 
This includes the <code>README</code>, <code>Makefile</code> an 
<code>runmains</code> files. The <code>main61</code> example program 
has been removed, since now LHAPDF can be loaded dynamically for 
<code>main42</code>, so that the two become equivalent. Further 
<code>main62</code> has been renamed <code>main43</code>, to gather 
HepMC-related examples, and <code>testlhef3.lhe</code> has been 
renamed <code>wbj_lhef3.lhe</code>. 
</li> 
 
<li>A new example <code>main30</code> how to create a tailormade 
copy of the ordinary event record, here with a history tracing 
of the hard process closer to the PYTHIA 6 conventions.</li> 
 
<li>The <?php $filepath = $_GET["filepath"];
echo "<a href='ProMCFiles.php?filepath=".$filepath."' target='page'>";?>ProMC</a> input-output file format 
is now implemented among the libraries that can be 
<code>configure</code>d to run with PYTHIA. An examples is provided 
in <code>main46.cc</code>. Thanks to Sergei Chekanov.</li> 
 
<li>Change in the setup of final-state-shower colour dipoles for the 
non-default case of no interleaving, whereby it becomes less likely 
to pick a colourless final-state particle as recoiler. New option 
<code>TimeShower:allowMPIdipole</code> gives more flexibility. 
Thanks to Mihoko Nojiri and Bryan Webber.</li> 
 
<li>New options 3 and 4 for <code>TimeShower:pTdampMatch</code> 
and <code>SpaceShower:pTdampMatch</code>, with new default 3 for the 
latter. The main effect is that, by default, <i>t tbar</i> 
production (as a <i>2 &rarr; 2</i> process) obtains damped 
radiation above the process scale. Thanks to Andy Buckley 
for suggestion.</li> 
 
<li>Initialization will now abort if a mode has been chosen with a 
non-allowed value. This applies to those modes that have been 
defined with the <code>modepick</code> or <code>modefix</code> 
labels in the <code>xmldoc/*.xml</code> files and, for the former, 
where maximal and minimal values have been specified. The former 
label is used to represent a discrete set of options, and so any 
value outside the allowed range is just plain wrong. Thanks to James 
Monk for suggestion.</li> 
 
<li>Added <code>getChannels()</code> in <code>SusyResonanceWidths</code> 
to dynamically create the decay table for SUSY particles and thereby 
to remove duplication in the XML file.</li> 
 
<li>New file <code>SusyWidthfunctions</code>. Added new class 
<code>WidthFunction</code> to handle calculation of three- and four-body 
decay widths.</li> 
 
<li>In <code>main24.cc</code> the example spectrum file format has been 
updated from SLHA1 to SLHA2, obtained from SoftSusy 3.5.1. The two new 
<code>slha1-example.spc</code> and <code>slha2-example.spc</code> files 
replace the older <code>cmssm.spc</code>, <code>snowmass2.spc</code> 
and <code>softsusy.spc</code> ones. </li> 
 
<li>In <code>LesHouches.cc</code> the read-in of LHEF headers 
containing tags that open and close on the same line (e.g., a single 
line containing <code>&lt;tag&gt;blabla&lt;/tag&gt;</code>) has been 
enabled. This could previously lead to improper initialization and 
crashes. Also implemented a check for forgotten close-tag statements, 
with a warning issued to the user. Thanks to Alexander Belyaev and 
Alexander Pukhov.</li> 
 
<li>Bug fix in the conversion from the <code>xml</code> settings files 
to their <code>php</code> radio-button equivalents, whereby the text 
describing some options was not set properly. (Typically the text of 
the previous option was repeated.) Thanks to Radek Ziebcik.</li> 
 
<li>Minor fixes in the LHEF version 3 reader. Introduce a new 
matching writer of LHEF version 1 or 3 files.</li> 
 
<li>Introduction of a new mode <code>LesHouches:setLeptonMass</code>, 
such that by default final-state charged leptons acquire sensible 
masses, even when the matrix-element calculations have been 
performed with massless leptons. This and other energy-momentum 
adjustments, e.g. for limited-precision storage, are now located in 
<code>ProcessContainer::constructProcess</code>.</li> 
 
<li><code>http://home.thep.lu.se/Pythia</code> has been introduced 
as a simpler but (hopefully) equivalent address to 
<code>http://home.thep.lu.se/~torbjorn/Pythia.html</code>, and 
various documentation has been updated accordingly. Thanks to 
Leif L&ouml;nnblad.</li> 
 
<li>New file <code>include/Pythia8Plugins/execinfo.h</code> contains 
trivial copies of three backtrace methods needed to be able to compile 
PYTHIA under Cygwin. The <code>README</code> file and the worksheet 
updated with brief information on three (non-supported) ways of working 
with PYTHIA under Windows. Thanks to and Theo Hughes and Gordon Watts. 
</li> 
 
<li>Bug fix in check for colour sextets and transfer of such colour 
information. Thanks to Alexander Belyaev and Alexander Pukhov.</li> 
 
<li>Improved handling of stray characters in the SUSY Les Houches 
code. The check on the consistency of decay tables has been removed. 
Improved warning/error printing in the SLHA interface.</li> 
 
<li>Bug fix in new beam remnant model, so that it basically 
operates like the old one for <i>e^+e^-</i> annihilation.</li> 
 
<li>Two bug fixes in the new colour reconnection model, one for 
diquarks at the ends of junction strings, and another to check that 
coloured resonances are processes with early resonance decays option.</li> 
 
<li>Bug fix for multiple <code>Pythia::init()</code> calls, where 
beam contents were not properly reset. Thanks to Josh Bendavid.</li> 
 
<li>Bug fix such that the valence content of a <i>pi^0</i>, 
<i>K^0_S</i>, <i>K^0_L</i> and Pomeron is reselected for each 
new event. Thanks to Radek Ziebcik.</li> 
 
<li>Fix typo in constants of the <i>tau &rarr; 3 pi</i> current 
for the amplitudes of the <i>rho</i>, <i>rho(1450)</i>, and 
<i>f2</i>. Thanks to Ian Nugent.</li> 
 
<li>Small bug fixes for string and ministring fragmentation, for the 
case when a low-mass (order 2 GeV) system contains at least three 
partons, which fail to define a unique direction for the final 
string region.</li> 
 
<li>New parameter <code>BeamRemnants:reducedKTatHighY</code> introduced 
to reduce technical problems with low-mass MPIs produced at high 
rapidities when primordial <i>kT</i> is introduced.</li> 
 
<li>Small bug fix in the global-recoil option for timelike showers.</li> 
 
<li>Update year to 2015, remove tabs and superfluous blanks, break long 
lines where meaningful, and some further minor changes.</li> 
 
</ul> 
</li> 
 
<li>8.201: 14 October 2014 
<ul> 
 
<li>The <i>Introduction to PYTHIA 8.2</i> has now been assigned 
the arXiv:1410.3012 [hep-ph] identifier, which has been introduced 
in code and text.</li> 
 
<li>The <code>enable-shared</code> by mistake was not listed 
among allowed configure options.</li> 
 
<li>Corrected a few tiny documentation typos.</li> 
 
</ul> 
</li> 
 
<li>8.200: 11 October 2014 
<ul> 
 
<li>A new <code>share/Pythia8</code> directory collects all 
documentation and example code. The <code>examples</code>, 
<code>htmldoc</code>, <code>phpdoc</code> and <code>xmldoc</code> 
directories have been moved here. The main-directory files 
<code>AUTHORS</code>, <code>COPYING</code>, <code>GUIDELINE</code> 
and <code>README</code> are also copied here during installation. 
</li> 
 
<li>A new <code>share/Pythia8/pdfdoc</code> directory collects pdf 
documents that are linked from the <code>htmldoc</code> and 
<code>phpdoc</code> directories. Over time it will  provide more 
in-depth descriptions of various physics aspects than offered in 
the html/php-formatted documentation. In addition to the official 
main publication and the worksheet, currently notes on LO vs. NLO 
PDFs and on the <i>g &rarr; q qbar</i> branching kernel are 
included.</li> 
 
<li>A new <code>include/Pythia8Plugins</code> directory collects 
code that does not form part of the core PYTHIA functionality but 
still has a general usefulness. Code in this directory will not be 
compiled as part of the Pythia library, but can be linked where needed. 
This new directory contains 
<ul> 
<li>the jet matching classes in <code>CombineMatchingInput.h</code>, 
<code>GeneratorInput.h</code> and <code>JetMatching.h</code>, moved 
from the <code>examples</code> directory;</li> 
<li>the <code>PowhegHooks</code> user hook, to veto shower emissions 
above the POWHEG scale, formerly found in <code>examples/main31.cc</code>; 
</li> 
<li>the <code>Pythia8ToHepMC</code> interface for output of PYTHIA events 
into the HepMC format, combining the code previously in 
<code>include/Pythia8ToHepMC.h</code> and 
<code>pythia8tohepmc/Pythia8ToHepMC.cc</code> into a new 
<code>HepMC2.h</code> file;</li> 
<li>the <code>FastJet3.h</code> interface of PYTHIA particles to the 
FastJet 3 library of jet finders, formerly found in 
<code>include/FastJet3.h</code>; and</li> 
<li>the <code>LHAPDF5.h</code> and <code>LHAPDF6.h</code> files for 
interfaces to the LHAPDF library (see further below).</li> 
</ul></li> 
 
<li>The configure/make structure has been considerably rewritten. 
Now all external libraries to be linked are specified in the 
main-directory <code>configure</code> step, along with other options, 
so there is no longer an <code>examples/configure</code>. The 
<code>make</code> step will, as before, compile and install libraries 
inside the current directory, such that the main programs in the 
<code>examples</code> directory can be run. One small difference is that 
also the archive libraries are installed in <code>lib</code> and not in 
<code>lib/archive</code>. 
<br/>A new optional <code>make install</code> step allows you to copy 
files to more convenient locations. The default option, with no directories 
specified in the <code>configure</code> step, requires you to have 
superuser privileges. Then files will be copied to standard locations 
as follows: 
<table border="0"> 
  <tr> <td>lib/</td> <td>&rarr;&nbsp;</td> <td>/usr/lib/</td> </tr> 
  <tr> <td>include/</td> <td>&rarr;&nbsp;</td> <td>/usr/include/</td> </tr> 
  <tr> <td>share/</td> <td>&rarr;&nbsp;</td> <td>/usr/share/</td> </tr> 
  <tr> <td>pythia-config</td> <td>&rarr;&nbsp;</td> <td>/usr/bin/</td> </tr> 
</table> 
</li> 
 
<li> 
The <code>pythia8-config.in</code> script has been replaced by a new 
<code>bin/pythia8-config</code> script. See the README file for details. 
The <code>make install</code> step by default will put a copy of it in 
<code>/usr/bin</code>. 
</li> 
 
<li>The interface to LHAPDF is now dynamically loaded when requested, 
and can be either to version 5 or 6 of the library. The dummy code 
previously in <code>lhapdfdummy/LHAPDFDummy.cc</code>, to be linked 
when LHAPDF is not, is no longer required. The two new files 
<code>LHAPDF5.h</code> and <code>LHAPDF6.h</code> in the 
<code>include/Pythia8Plugins</code> directory contain the necessary 
interface code. The selection of PDF sets, notably for the proton, 
has been extended to simplify mixing of internal and external PDF sets, 
and it is now possible to specify different PDFs for the two incoming 
protons at the LHC, see the <?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF Selection</a> 
description.</li> 
 
<li>The new <code>LHEF3.h</code> file contains a generic interface for 
reading Les Houches Event Files of versions 1.0, 2.0 and 3.0. This 
allows more information to be read and studied by the author. Currently 
PYTHIA itself makes little use of the information beyond the one in 1.0, 
but it is available among the 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a>. 
Examples are found in <code>main37.cc</code> and <code>main38.cc</code>. 
</li> 
 
<li>The new <code>Beams:strictLHEFscale</code> switch can be used to 
restrict parton showers in resonance decays to be below the input 
Les Houches scale, not only the hard process itself. 
The new <code>Beams:setProductionScalesFromLHEF</code> switch can be used 
to restrict the emission off each separate parton to be below its specific 
scale.</li> 
 
<li>The <code>rootexamples</code> directory has been removed, and the 
two programs <code>examples/main91</code> and <code>examples/main92</code> 
now illustrate how ROOT can be used in conjunction with PYTHIA.</li> 
 
<li>The executable built from <code>examples/mainxx.cc</code> is now 
named <code>examples/mainxx</code>, while previously it was named 
<code>examples/mainxx.exe</code>.</li> 
 
<li>The rudimentary support for compilation on Windows platforms, 
present in PYTHIA 8.1, has not yet been updated for 8.2 and so is omitted. 
Also the README.HepMC file is omitted for now.</li> 
 
<li>The ProMC interface is broken, and has been removed for now.</li> 
 
<li>Several methods have been removed from the <code>Event</code> class 
since the properties now instead can be accessed from the individual 
<code>Particle</code> instance, if this particle belongs to an event. 
These include <code>iTopCopy</code>, <code>iBotCopy</code>, 
<code>iTopCopyId</code>, <code>iBotCopyId</code>,<code>motherList</code>, 
<code>daughterList</code>, <code>sisterList</code>, 
<code>sisterListTopBot</code>, <code>isAncestor</code>, 
<code>statusHepMC</code> and <code>undoDecay</code>.</li> 
 
<li>A number of deprecated <code>Pythia::init(...)</code> methods with 
varying arguments have been removed. Instead call <code>init()</code> 
without any arguments and use 
<?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>Beam Parameters</a> settings to 
specify beams and energies in different ways.</li> 
 
<li>The deprecated <code>Pythia::statistics(...)</code> method has been 
removed; instead use <code>Pythia::stat(...)</code>.</li> 
 
<li>Several settings in the <code>Main:</code> series have been removed. 
Most of these have already found replacements in the <code>Init:</code>, 
<code>Next:</code> and <code>Stat:</code> ones, and have been marked as 
deprecated. Four further ones were deemed so peripheral that they were 
removed altogether, but of course the underlying functionality remains. 
</li> 
 
<li>A few aliases for (parts of) settings names have been removed. 
Previously "Multiple" was mapped to "Multiparton", "MI" to "MPI" and 
"minBias" to "nonDiffractive" if a settings name was not found for the 
original input string.</li> 
 
<li>The default tune has been changed from 4C to Monash 2013, meaning 
<code>Tune:ee = 7</code> and <code>Tune:pp = 14</code>. The old 4C 
tune that was default in 8.1 can be recovered with 
<code>Tune:ee = 3</code> and <code>Tune:pp = 5</code>. 
Also most other older tunes are based on <code>Tune:ee = 3</code>. 
</li> 
 
<li>Two new CMS underlying-event tunes [<a href="Bibliography.php" target="page">CMS14</a>] and the ATLAS 
AZ tune [<a href="Bibliography.php" target="page">ATL14</a>] have been added as options.</li> 
 
<li>The default handling of the <i>g &rarr; q qbar</i> splitting kernel 
has been changed, affecting in particular heavy-flavour production. 
<code>TimeShower:weightGluonToQuark</code> has been changed from 1 to 4 
to do this. All old tunes are with the 1 value but, since the tunes are 
not probing the detailed <i>g &rarr; q qbar</i> behaviour, this is 
not set as part of the tune options.</li> 
 
<li>Christine O. Rasmussen joins as new PYTHIA collaboration member.</li> 
 
<li>A new model for the handling of <?php $filepath = $_GET["filepath"];
echo "<a href='BeamRemnants.php?filepath=".$filepath."' target='page'>";?>beam 
remnants</a> as an option to the old one, which remains as default 
for now.</li> 
 
<li>Two new models for <?php $filepath = $_GET["filepath"];
echo "<a href='ColourReconnection.php?filepath=".$filepath."' target='page'>";?>colour 
reconnection</a>, one quite sophisticated and one simpler. 
This involves several new classes and files. It also includes some 
changes in the hadronization framework, notably for the handling of 
junctions. The old model remains as default for now. The 
<code>BeamRemnants:reconnectColours</code> flag to switch on/off 
reconnection has been renamed <code>ColourReconnection:reconnect</code>, 
the main parameter <code>BeamRemnants:reconnectRange</code> of the old 
model has been renamed <code>ColourReconnection:range</code>, and several 
new settings have been introduced, notably 
<code>ColourReconnection:mode</code> to switch among the three models. 
</li> 
 
<li>A new <code>include/Pythia8Plugins/ColourReconnectionHooks.h</code> 
makes available an even larger selection of toy colour reconnection 
models, via user hooks. Some of them are only intended for top decays, 
for top mass uncertainty studies, whereas others can be used more 
generally. The <code>examples/main29.cc</code> program illustrates how 
the different options should be set up.</li> 
 
<li>Several new features and improvements in the matching/merging 
machinery. Notably the aMC@NLO matching scheme has been implemented, 
see the <?php $filepath = $_GET["filepath"];
echo "<a href='aMCatNLOMatching.php?filepath=".$filepath."' target='page'>";?>aMC@NLO Matching</a> 
description. To this end the global-recoil option of timelike showers 
has been improved, and security checks have been introduced for 
inaccurate LHEF input. A new <code>main89.cc</code> example has been 
introduced, where different <code>.cmnd</code> files show how to set 
up either CKKW-L, FxFx, MLM, UMEPS or UNLOPS merging.</li> 
 
<li>Improved capability for the <code>LHAup</code> Les Houches interface 
to read SLHA information embedded in the input file or stream.</li> 
 
<li>The <code>Makefile</code>s have been updated to take into account 
the changed structure of the HepMC interface.</li> 
 
<li>The <i>Z'</i> production process has been updated to optionally 
allow decay to a fourth generation of fermions, with universal or 
non-universal couplings.</li> 
 
<li>Introduction of a new Higgs CP-mixing parametrization via a mixing 
angle <i>phi</i> as described in <?php $filepath = $_GET["filepath"];
echo "<a href='HiggsProcesses.php?filepath=".$filepath."' target='page'>";?>Higgs 
Processes</a>. The choice of the Higgs CP-mixing parametrization 
now also affects the distributions of the <i>tau</i> decay products 
from the processes <i>H^0 &rarr; tau^+ tau^-</i>. 
 
<li>Bug fix in <i>H^0 &rarr; W^+ W^- &rarr; 4 f</i> matrix element 
for mixed CP-state case.</li> 
 
<li>Various improvements and finer grain control for the determination 
of <i>tau</i> decay correlations and <i>tau</i> polarizations. By 
default the decays of <i>tau</i> pairs from known resonance decays 
in Les Houches input are now correlated. 
The <code>ParticleDecays:sophisticatedTau</code> mode 
in <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>Particle Decays</a> has been renamed 
<code>TauDecays:mode</code>, as well as all <i>tau</i>-related 
<code>ParticleDecay</code> options, with two new options of 
using only the internal machinery to determine correlations and 
polarizations, and only using the provided SPINUP digit from Les 
Houches input. The option <code>TauDecays:externalMode</code> has been 
introduced to control the interpretation of the SPINUP digit. 
</li> 
 
<li>For Les Houches Event input the energy of a particle is recalculated 
from its three-momentum and mass, in order to limit mismatches from 
limited numerical precision in the input values.</li> 
 
<li>Bug fix in the two-loop running <i>alpha_s</i>, for the matching 
to six flavours at the top mass.</li> 
 
<li>Eliminate harmless compiler warnings for <code>FJcore</code>.</li> 
 
<li>Updated Introduction (= the official 8.2 article) and Worksheet.</li> 
 
</ul> 
</li> 
 
</ul> 
 
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
