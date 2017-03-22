<html>
<head>
<title>Semi-Internal Processes</title>
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

<form method='post' action='SemiInternalProcesses.php'>
 
<h2>Semi-Internal Processes</h2> 
 
Normally users are expected to implement new processes via the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>. Then 
you do all flavour, colour and phase-space selection externally, 
before your process-level events are input for further processing 
by PYTHIA. However, it is also possible to implement a 
new process in exactly the same way as the internal PYTHIA 
ones, thus making use of the internal phase space selection machinery 
to sample an externally provided cross-section expression. 
The MadGraph 5 program [<a href="Bibliography.php" target="page">Alw11</a>] allows you to do exactly that, 
i.e. it can be used to generate C++ code that can be linked into 
the existing PYTHIA framework, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='MadGraph5Processes.php?filepath=".$filepath."' target='page'>";?>here</a>. 
 
<p/> 
Should you decide to go ahead on your own, 
this page gives a brief summary how to do that. If you additionally 
want to introduce a new resonance species, with its own internal 
width calculations, you will find further instructions 
<?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalResonances.php?filepath=".$filepath."' target='page'>";?>here</a>. It is strongly 
recommended to shop around for a similar process that has already 
been implemented, and to use that existing code as a template. 
Look for processes with the same combinations of incoming flavours 
and colour flows, rather than the shape of the cross section itself. 
With a reasonable such match the task should be of medium difficulty, 
without it more demanding. 
 
<p/> 
PYTHIA's internal phase-space generators are rather good at handling 
the phase space of <i>2 &rarr; 1</i> and <i>2 &rarr; 2</i> processes, 
are more primitive for <i>2 &rarr; 3</i> ones and do not at all 
address higher multiplicities. An option is therefore also provided 
for external phase-space generators to be used, which must then be 
encapsulated to inherit from PYTHIA's <code>PhaseSpace</code> base class 
(or one of its derivatives). 
The set of processes that can be implemented in this framework is 
therefore in principle unlimited, though the user must supply external 
phase-space generators for non-trivial <i>2 &rarr; 3</i> processes 
and all higher <i>2 &rarr; n</i> multiplicities. 
Note, however, that the produced particles may be resonances, so it is 
possible to end up with bigger "final" multiplicities through sequential 
decays, also with the internal phase-space generators, and to include 
further matrix-element weighting in those decays. 
 
<p/> 
For processes using PYTHIA's internal phase-space generators, there 
are three steps involved in implementing a process: 
<ol> 
<li>making use of the PYTHIA-provided kinematics information to 
calculate the relevant cross section,</li> 
<li>writing a new class,  where the matrix elements are implemented, 
including information on incoming and outgoing flavours and colours, 
and</li> 
<li>making the process available.</li> 
</ol> 
We consider these aspects in turn. An example where it all comes 
together is found in <code>main22.cc</code>. 
 
<p/> 
For processes for which an external phase-space generator will be 
used, step 1 above changes to writing a new class, where the 
phase-space generator is implemented, and making use of that to 
calculate the relevant cross section. There are no example programs 
illustrating how to do this yet, but the methodology is described 
below, under "Implementing an external phase-space generator". 
 
<h3>The Cross Section Calculation</h3> 
 
The key method for the cross section calculation is 
<code>SigmaProcess::sigmaHat()</code>, described below. At the point when 
it is called, the kinematics has already been set up, and from these 
phase space variables the differential cross section is to be calculated. 
 
<p/> 
For a <i>2 &rarr; 1</i> process, the returned value should be 
<i>sigmaHat(sHat)</i>, where <code>mH</code> (= <i>mHat</i>), 
<code>sH</code> (= <i>sHat</i>) and <code>sH2</code> (= <i>sHat^2</i>) 
are available to be used. Incoming partons are massless. Overload the 
<code>convertM2()</code> method below if you instead plan to return 
<i>|M|^2</i>. 
 
<p/> 
For a <i>2 &rarr; 2</i> process, instead <i>d(sigmaHat)/d(tHat)</i> 
should be returned, based on provided 
<code>mH, sH, sH2, tH, tH2, uH, uH2, m3, s3, m4, s4</code> and 
<code>pT2</code> values (<code>s3 = m3*m3</code> etc.). Incoming 
partons are massless. Overload the <code>convertM2()</code> method 
below if you instead plan to return <i>|M|^2</i>. 
 
<p/> 
For a <i>2 &rarr; 3</i> process, instead <i>|M|^2</i> should be 
returned, with normalization such that <i>|M|^2 / (2 sHat)</i> integrated 
over the three-body phase space gives the cross section. Here no standard 
set of Mandelstam-style variables exists. Instead the obvious ones, 
<code>mH, sH, m3, s3, m4, s4, m5, s5</code>, are complemented by the 
four-vectors <code>p3cm, p4cm, p5cm</code>, from which further invariants 
may be calculated. The four-vectors are defined in the CM frame of the 
subcollision, with massless incoming partons along the <i>+-z</i> axis. 
 
<p/> 
In either case, <i>alpha_s</i> and <i>alpha_em</i> have already 
been calculated, and are stored in <code>alpS</code> and <code>alpEM</code>. 
Also other standard variables may be used, like 
<code>CoupEW::sin2thetaW()</code>, and related flavour-dependent 
vector and axial couplings in <code>CoupEW</code> and CKM combinations 
in <code>VCKM</code>. 
 
<p/> 
In case some of the final-state particles are resonances, their 
squared masses have already been selected according to a Breit-Wigner 
with a linearly running width <i>Gamma(m) = Gamma(m_0) * m / m_0</i>. 
More precisely, the mass spectrum is weighted according to 
<i>w_BW(m^2) d(m^2)</i>, where 
<br/><i> 
w_BW(m^2) = (1/pi) * (m * Gamma(m)) / ( (m^2 - m_0^2)^2 + (m * Gamma(m))^2 ) . 
</i><br/> 
If you would like to have another expression, the above weights are stored 
in <code>runBW3</code>, <code>runBW4</code> and <code>runBW5</code>, 
respectively. If you divide out one of these factors, you just remain with 
a phase space selection <i>d(m^2)</i> for this particle, 
and can multiply on your desired shape factor instead. Unfortunately, the 
Monte Carlo efficiency will drop if your new mass distribution differs 
dramatically from the input one. Therefore it does make sense to adjust the 
database value of the width to be slightly (but not too much) broader 
than the distribution you have in mind. Also note that, already by default, 
the wings of the Breit-Wigner are oversampled (with a compensating lower 
internal weight) by partly sampling like <i>(a + b/m^2 + c/m^4) d(m^2)</i>, 
where the last term is only used for <i>gamma^*/Z^0</i>. 
 
<p/> 
As alternative to the kinematics variables defined above, also the two 
arrays <code>mME[5]</code> and <code>pME[5]</code>, for masses and 
four-momenta, respectively, can be used for cross-section calculations. 
Here indices 0 and 1 are the two incoming beams, and 2 and onwards the 
outgoing particles. Note that this differs by one step from the normal 
internal labeling, where slot 0 is left empty. The four-momenta are 
defined in the rest frame of the subcollision, with the incoming partons 
along the <i>+-z</i> direction. The kinematics need not agree with the 
"correct" one stored in the event record, for three reasons. 
<br/>1) Gauge invariance forces matrix-element calculations to use 
the same masses for incoming as outgoing legs of a particle species, 
say <i>b</i> quarks. Therefore the kinematics of the two incoming 
partons is recalculated, relative to the normal event record, to put 
the partons on the mass shell. (Note that initial masses is a technical 
issue, not the correct physics picture: the incoming partons are likely 
to be spacelike virtual rather than on the mass shell.) 
<br/>2) In principle each fermion flavour has to be treated separately, 
owing to a different mass. However, in many cases fermions can be 
assumed massless, which speeds up the calculations, and further gains 
occur if then different flavours can use the same cross-section 
expression. In MadGraph the default is that fermions up to and including 
the <i>c</i> quark and the <i>mu</i> lepton are considered massless, 
while the <i>b</i> quark and the <i>tau</i> lepton are considered 
massive. This can be modified however, and below we provide four flags 
that can be used to consider the "borderline" fermions either as 
massless or as massive when matrix elements are evaluated, to match the 
assumptions made for the matrix elements themselves. 
<br/>3) For <i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes of massive 
identical particles (or antiparticles) in the final state, such as 
<i>t tbar</i> or <i>W^+ W^-</i>, the kinematics is here adjusted 
so that the two or three particles have the same mass, formed as a 
suitable average of the actual Breit-Wigner-distributed masses. This 
allows the evaluation of matrix-element expressions that only have 
meaning if the two/three have the same mass. 
<br/>Thus the mass array <code>mME[5]</code> and the four-momentum array 
<code>pME[5]</code> present values both for initial- and final-state 
particles based on these mass principles suited for matrix-element input. 
Note that these variables therefore differ from the kinematics stored in 
the event record proper, where incoming fermions are always massless and 
outgoing resonances have independent Breit-Wigner mass distributions. 
<br/>The conversion from the normal to the special kinematics is done 
by calling the <code>setupForME()</code> method. This you have to do 
yourself in the <code>SigmaHat()</code> member of your derived class. 
Alternatively it could be done in <code>SigmaKin()</code>, i.e. before 
the loop over incoming flavours, but then these would be considered 
massless. The identity of final-state particles is obtained from the 
<code>id3Mass()</code>, <code>id4Mass()</code> and <code>id5Mass()</code> 
methods. Should the conversion to <code>mME[5]</code> and 
<code>pME[5]</code> not work, <code>setupForME()</code> will return 
<code>false</code>, and then the cross section should be put zero. 
 
<br/><br/><strong>SigmaProcess:cMassiveME</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Let the <i>c</i> quark be massive or not in the kinematics set up for 
external matrix-element evaluation. 
   
 
<br/><br/><strong>SigmaProcess:bMassiveME</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Let the <i>b</i> quark be massive or not in the kinematics set up for 
external matrix-element evaluation. 
   
 
<br/><br/><strong>SigmaProcess:muMassiveME</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Let the <i>mu</i> lepton be massive or not in the kinematics set up for 
external matrix-element evaluation. 
   
 
<br/><br/><strong>SigmaProcess:tauMassiveME</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Let the <i>tau</i> lepton be massive or not in the kinematics set up for 
external matrix-element evaluation. 
   
 
 
<h3>The Cross Section Class</h3> 
 
The matrix-element information has to be encoded in a new class. 
The relevant code could either be put before the main program in the 
same file, or be stored separately, e.g. in a matched pair 
of <code>.h</code> and <code>.cc</code> files. The latter may be more 
convenient, in particular if the cross sections are lengthy, or if you 
intend to build up your own little process library, but of course 
requires that these additional files are correctly compiled and linked. 
 
<p/> 
The class has to be derived either from 
<code>Sigma1Process</code>, for <i>2 &rarr; 1</i> processes, from 
<code>Sigma2Process</code>, for <i>2 &rarr; 2</i> ones, or from 
<code>Sigma3Process</code>, for <i>2 &rarr; 3</i> ones. (The 
<code>Sigma0Process</code> class is used for elastic, diffractive 
and minimum-bias events, and is not recommended for use beyond that.) 
These are in their turn derived from the <code>SigmaProcess</code> 
base class. 
 
<p/> 
The class can implement a number of methods. Some of these are 
compulsory, others strongly recommended, and the rest are to be 
used only when the need arises to override the default behaviour. 
The methods are: 
 
<p/> 
A <b>constructor</b> for the derived class obviously must be available. 
Here you are quite free to allow a list of arguments, to set 
the parameters of your model, or even to create a set of closely 
related but distinct processes. For instance, <i>g g &rarr; Q Qbar</i>, 
<i>Q = c</i> or <i>b</i>, is only coded once, and then the 
constructor takes the quark code (4 or 5)  as argument, 
to allow the proper amount of differentiation. 
 
<p/> 
A <b>destructor</b> is only needed if you plan to delete the process 
before the natural end of the run, and require some special behaviour 
at that point. If you call such a destructor you will leave a pointer 
dangling inside the <code>Pythia</code> object you gave it in to, 
if that still exists. 
 
<a name="method1"></a>
<p/><strong>void SigmaProcess::initProc() &nbsp;</strong> <br/>
is called once during initialization, and can then be used to set up 
parameters, such as masses and couplings, and perform calculations 
that need not be repeated for each new event, thereby saving time. 
This method needs not be implemented, since in principle all 
calculations can be done in <code>sigmaHat</code> below. 
   
 
<a name="method2"></a>
<p/><strong>void SigmaProcess::sigmaKin() &nbsp;</strong> <br/>
is called once a kinematical configuration has been determined, but 
before the two incoming flavours are known. This routine can therefore 
be used to perform calculations that otherwise might have to be repeated 
over and over again in <code>sigmaHat</code> below. For instance 
a flavour-independent cross section calculation for a <i>q g</i> 
initial state would be repeated 20 times in <code>sigmaHat</code>, 
five times for the five quark flavours allowed in the incoming beams, 
times twice to include antiquarks, times twice since the (anti)quark 
could be in either of the two beams. You could therefore calculate the 
result once only and store it as a private data member of the class. 
It is optional whether you want to use this method, however, or put 
everything in <code>sigmaHat</code>. 
   
 
<a name="method3"></a>
<p/><strong>double SigmaProcess::sigmaHat() &nbsp;</strong> <br/>
is the key method for cross section calculations and returns a cross section 
value, as described in the previous section. It is called when also a 
preliminary set of incoming flavours has been picked, in addition to the 
kinematical ones already available for <code>sigmaKin</code>. 
Typically <code>sigmaHat</code> is called inside a loop over all allowed 
incoming flavour combinations, stored in <code>id1</code> and 
<code>id2</code>, with fixed kinematics, as already illustrated above. 
The sum over the different flavour combinations provides the total 
cross section, while their relative size is used to make a selection of 
a specific incoming state. 
   
 
<a name="method4"></a>
<p/><strong>bool SigmaProcess::setupForME() &nbsp;</strong> <br/>
to be called by the user from inside <code>sigmaHat()</code> 
(or possibly <code>sigmaKin()</code>) to setup alternative kinematics 
in the <code>mME[5]</code> and <code>pME[5]</code> arrays, better 
suited for matrix-element calculations. See the end of the previous 
section for a more detailed description. Should the method return 
<code>false</code> then the conversion did not work, and 
<code>sigmaHat()</code> (or <code>sigmaKin()</code>) should be set to 
vanish. 
   
 
<a name="method5"></a>
<p/><strong>void SigmaProcess::setIdColAcol() &nbsp;</strong> <br/>
is called only once an initial state and a kinematical configuration has 
been picked. This routine must set the complete flavour information and 
the colour flow of the process. This may involve further random choices, 
between different possible final-state flavours or between possible 
competing colour flows. Private data members of the class may be used to 
retain some information from the previous steps above. 
<br/>When this routine is called the two incoming flavours have already 
been selected and are available in <code>id1</code> and <code>id2</code>, 
whereas the one, two or three outgoing ones either are fixed for a given 
process or can be determined from the instate (e.g. whether a <i>W^+</i> 
or <i>W^-</i> was produced).  There is also a standard method in 
<code>VCKM</code> to pick a final flavour from an initial one with CKM 
mixing. Once you have figured out the value of 
<code>id3</code> and, the case being, <code>id4</code> and 
<code>id5</code>, you store these values permanently by a call 
<code>setId( id1, id2, id3, id4, id5)</code>, where the last two may be 
omitted if irrelevant. 
<br/>Correspondingly, the colours are stored with 
<code>setColAcol( col1, acol1, col2, acol2, col3, acol3, col4, acol4, 
col5, acol5)</code>, where the final ones may be omitted if irrelevant. 
Les Houches style colour tags are used, but starting with number 1 
(and later shifted by the currently requested offset). The 
input is grouped particle by particle, with the colour index before the 
anticolour one. You may need to select colour flow dynamically, depending 
on the kinematics, when several distinct possibilities exist. Trivial 
operations, like swapping colours and anticolours, can be done with 
existing methods. 
<br/>When the <code>id3Mass()</code> and <code>id4Mass()</code> 
methods have been used, the order of the outgoing particles may be 
inconsistent with the way the <i>tHat</i> and <i>uHat</i> 
variables have been defined. A typical example would be a process like 
<i>q g &rarr; q' W</i> with <i>tHat</i> defined between incoming and 
outgoing quark, but where <code>id3Mass() = 24</code> and so the 
process is to be stored as <i>q g &rarr; W q'</i>. One should then put 
the variable <code>swapTU = true</code> in <code>setIdColAcol()</code> 
for each event where the <i>tHat</i> and <i>uHat</i> variables 
should be swapped before the event kinematics is reconstructed. This 
variable is automatically restored to <code>false</code> for each new 
event. 
   
 
<a name="method6"></a>
<p/><strong>double SigmaProcess::weightDecayFlav( Event& process) &nbsp;</strong> <br/>
is called to allow a reweighting of the simultaneous flavour choices of 
resonance decay products. Is currently only used for the 
<i>q qbar &rarr; gamma*/Z^0 gamma*/Z^0</i> process, and will likely not 
be of interest for you. 
   
 
<a name="method7"></a>
<p/><strong>double SigmaProcess::weightDecay( Event& process, int iResBeg, int iResEnd) &nbsp;</strong> <br/>
is called when the basic process has one or several resonances, after each 
set of related resonances in <code>process[i]</code>, 
<code>iResBeg</code> &lt;= <code>i </code> &lt;= <code>iResEnd</code>, 
has been allowed to decay. The calculated weight, to be normalized 
to the range between 0 and 1, is used to decide whether to accept the 
decay(s) or try for a new decay configuration. The base-class version of 
this method returns unity, i.e. gives isotropic decays by default. 
This method may be called repeatedly for a single event. For instance, in 
<i>q qbar &rarr; H^0 Z^0</i> with <i>H^0 &rarr; W^+ W^-</i>, a first call 
would be made after the <i>H^0</i> and <i>Z^0</i> decays, and then 
depend only on the <i>Z^0</i> decay angles since the <i>H^0</i> 
decays isotropically. The second call would be after the <i>W^+ W^-</i> 
decays and then involve correlations between the four daughter fermions. 
   
 
<a name="method8"></a>
<p/><strong>string SigmaProcess::name() &nbsp;</strong> <br/>
returns the name of the process, as you want it to be shown in listings. 
   
 
<a name="method9"></a>
<p/><strong>int SigmaProcess::code() &nbsp;</strong> <br/>
returns an integer identifier of the process. This has no internal function, 
but is only intended as a service for the user to rapidly (and hopefully 
uniquely) identify which process occurred in a given event. Numbers below 
10000 are reserved for internal PYTHIA use. 
   
 
<a name="method10"></a>
<p/><strong>string SigmaProcess::inFlux() &nbsp;</strong> <br/>
this string specifies the combinations of incoming partons that are 
allowed for the process under consideration, and thereby which incoming 
flavours <code>id1</code> and <code>id2</code> the <code>sigmaHat()</code> 
calls will be looped over. It is always possible to pick a wider flavour 
selection than strictly required and then put to zero cross sections in 
the superfluous channels, but of course this may cost some extra execution 
time. Currently allowed options are: 
<br/>* <code>gg</code>: two gluons. 
<br/>* <code>qg</code>: one (anti)quark and one gluon. 
<br/>* <code>qq</code>: any combination of two quarks, two antiquarks or 
a quark and an antiquark. 
<br/>* <code>qqbar</code>: any combination of a quark and an antiquark; 
a subset of the <code>qq</code> option. 
<br/>* <code>qqbarSame</code>: a quark and its antiquark; 
a subset of the <code>qqbar</code> option. 
<br/>* <code>ff</code>: any combination of two fermions, two antifermions 
or a fermion and an antifermion; is the same as <code>qq</code> for 
hadron beams but also allows processes to work with lepton beams. 
<br/>* <code>ffbar</code>: any combination of a fermion and an antifermion; 
is the same as <code>qqbar</code> for hadron beams but also allows processes 
to work with lepton beams. 
<br/>* <code>ffbarSame</code>: a fermion and its antifermion; is the 
same as <code>qqbarSame</code> for hadron beams but also allows processes 
to work with lepton beams. 
<br/>* <code>ffbarChg</code>: a fermion and an antifermion that combine 
to give charge +-1. 
<br/>* <code>fgm</code>: a fermion and a photon (gamma). 
<br/>* <code>ggm</code>: a gluon and a photon. 
<br/>* <code>gmgm</code>: two photons. 
   
 
<a name="method11"></a>
<p/><strong>bool SigmaProcess::convert2mb() &nbsp;</strong> <br/>
it is assumed that cross sections normally come in dimensions such that 
they, when integrated over the relevant phase space, obtain the dimension 
GeV^-2, and therefore need to be converted to mb. If the cross section 
is already encoded as mb then <code>convert2mb()</code> should be 
overloaded to instead return <code>false</code>. 
   
 
<a name="method12"></a>
<p/><strong>bool SigmaProcess::convertM2() &nbsp;</strong> <br/>
it is assumed that <i>2 &rarr; 1</i> cross sections are encoded as 
<i>sigmaHat(sHat)</i>, and <i>2 &rarr; 2</i> ones as 
<i>d(sigmaHat)/d(tHat)</i> in the <code>SigmaProcess::sigmaHat()</code> 
methods. If <code>convertM2()</code> is overloaded to instead return 
<code>true</code> then the return value is instead assumed to be the 
squared matrix element <i>|M|^2</i>, and 
<code>SigmaProcess::sigmaHatWrap(...)</code> converts to 
<i>sigmaHat(sHat)</i> or <i>d(sigmaHat)/d(tHat)</i>, respectively. 
This switch has no effect on <i>2 &rarr; 3</i> processes, where 
<i>|M|^2</i> is the only allowed input anyway. 
   
 
<a name="method13"></a>
<p/><strong>int SigmaProcess::id3Mass() &nbsp;</strong> <br/>
   
<strong>int SigmaProcess::id4Mass() &nbsp;</strong> <br/>
   
<strong>int SigmaProcess::id5Mass() &nbsp;</strong> <br/>
are the one, two or three final-state flavours, where masses are to be 
selected before the matrix elements are evaluated. Only the absolute value 
should be given. For massless particles, like gluons and photons, one need 
not give anything, i.e. one defaults to 0. The same goes for normal light 
quarks, where masses presumably are not implemented in the matrix elements. 
Later on, these quarks can still (automatically) obtain constituent masses, 
once a <i>u</i>, <i>d</i> or <i>s</i> flavour has been selected. 
   
 
<a name="method14"></a>
<p/><strong>int SigmaProcess::resonanceA() &nbsp;</strong> <br/>
   
<strong>int SigmaProcess::resonanceB() &nbsp;</strong> <br/>
are the codes of up to two <i>s</i>-channel resonances contributing to 
the matrix elements. These are used by the program to improve the phase-space 
selection efficiency, by partly sampling according to the relevant 
Breit-Wigner distributions. Massless resonances (the gluon and photon) 
need not be specified. 
   
 
<a name="method15"></a>
<p/><strong>bool SigmaProcess::isSChannel() &nbsp;</strong> <br/>
normally the choice of renormalization and factorization scales in 
<i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes is based on the 
assumption that <i>t</i>- and <i>u</i>-channel exchanges dominates the 
cross section. In cases such as <i>f fbar &rarr; gamma* &rarr; f' fbar'</i> 
a <i>2 &rarr; 2</i> process actually ought to be given scales as a 
<i>2 &rarr; 1</i> one, in the sense that it proceeds entirely through 
an <i>s</i>-channel resonance. This can be achieved if you override the 
default <code>false</code> to return <code>true</code>. See further the 
page on <?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>couplings and scales</a>. 
   
 
<a name="method16"></a>
<p/><strong>int SigmaProcess::idSChannel() &nbsp;</strong> <br/>
normally no intermediate state is shown in the event record for 
<i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes. However, in case 
that <code>idSChannel</code> is overloaded to return a nonzero value, 
an intermediate particle with that identity code is inserted into the 
event record, to make it a <i>2 &rarr; 1 &rarr; 2</i> or 
<i>2 &rarr; 1 &rarr; 3</i> 
process. Thus if both <code>isSChannel</code> and <code>idSChannel</code> 
are overloaded, a process will behave and look like it proceeded through 
a resonance. The one difference is that the implementation of the 
matrix element is not based on the division into a production and a 
decay of an intermediate resonance, but is directly describing the 
transition from the initial to the final state. 
   
 
<a name="method17"></a>
<p/><strong>int SigmaProcess::isQCD3body() &nbsp;</strong> <br/>
there are two different 3-body phase-space selection machineries, 
of which the non-QCD one is default. If you overload this method 
instead the QCD-inspired machinery will be used. The differences 
between these two is related to which 
<?php $filepath = $_GET["filepath"];
echo "<a href='PhaseSpaceCuts.php?filepath=".$filepath."' target='page'>";?>phase space cuts</a> 
can be set, and also that the QCD machinery assumes (almost) massless 
outgoing partons. 
   
 
<a name="method18"></a>
<p/><strong>int SigmaProcess::idTchan1() &nbsp;</strong> <br/>
   
<strong>int SigmaProcess::idTchan2() &nbsp;</strong> <br/>
the non-QCD <i>2 &rarr; 3</i> phase space selection machinery is rather 
primitive, as already mentioned. The efficiency can be improved in 
processes that proceed though <i>t</i>-channel exchanges, such as 
<i>q qbar' &rarr; H^0 q qbar'</i> via <i>Z^0 Z^0</i> fusion, if the 
identity of the  <i>t</i>-channel-exchanged particles on the two side 
of the event are provided. Only the absolute value is of interest. 
   
 
<a name="method19"></a>
<p/><strong>double SigmaProcess::tChanFracPow1() &nbsp;</strong> <br/>
   
<strong>double SigmaProcess::tChanFracPow2() &nbsp;</strong> <br/>
in the above kind of <i>2 &rarr; 3</i> phase-space selection, the 
sampling of <i>pT^2</i> is done with one part flat, one part weighted 
like <i>1 / (pT^2 + m_R^2)</i> and one part  like 
<i>1 / (pT^2 + m_R^2)^2</i>. The above values provide the relative 
amount put in the latter two channels, respectively, with the first 
obtaining the rest. Thus the sum of <code>tChanFracPow1()</code> and 
<code>tChanFracPow2()</code> must be below unity. The final results 
should be independent of these numbers, but the Monte Carlo efficiency 
may be quite low for a bad choice. Here <i>m_R</i> is the mass of the 
exchanged resonance specified by <code>idTchan1()</code> or 
<code>idTchan2()</code>. Note that the order of the final-state 
listing is important in the above <i>q qbar' &rarr; H^0 q qbar'</i> example, 
i.e. the <i>H^0</i> must be returned by <code>id3Mass()</code>, 
since it is actually the <i>pT^2</i> of the latter two that are 
selected independently, with the first <i>pT</i> then fixed 
by transverse-momentum conservation. 
   
 
<a name="method20"></a>
<p/><strong>bool SigmaProcess::useMirrorWeight() &nbsp;</strong> <br/>
in <i>2 &rarr; 3</i> processes the phase space selection used here 
involves a twofold ambiguity basically corresponding to a flipping of 
the positions of last two outgoing particles. These are assumed equally 
likely by default, <code>false</code>, but for processes proceeding entirely 
through <i>t</i>-channel exchange the Monte Carlo efficiency can be 
improved by making a preselection based on the relative propagator 
weights, <code>true</code>. 
   
 
<a name="method21"></a>
<p/><strong>int SigmaProcess::gmZmode() &nbsp;</strong> <br/>
allows a possibility to override the global mode 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ElectroweakProcesses.php?filepath=".$filepath."' target='page'>";?>WeakZ0:gmZmode</a></code> 
for a specific process. The global mode normally is used to switch off 
parts of the <i>gamma^*/Z^0</i> propagator for test purposes. The 
above local mode is useful for processes where a <i>Z^0</i> really is 
that and nothing more, such as <i>q qbar &rarr; H^0 Z^0</i>. The default 
value -1 returned by <code>gmZmode()</code> ensures that the global 
mode is used, while 0 gives full <i>gamma^*/Z^0</i> interference, 
1 <i>gamma^*</i> only and 2 <i>Z^0</i> only. 
   
 
<h3>Access to a process</h3> 
 
Once you have implemented a class, it is straightforward to make use of 
it in a run. Assume you have written a new class <code>MySigma</code>, 
which inherits from <code>Sigma1Process</code>, <code>Sigma2Process</code> 
or <code>Sigma3Process</code>, which in their turn inherit from 
<code>SigmaProcess</code>. You then create an instance of this class 
and hand it in to a <code>pythia</code> object with 
<pre> 
      SigmaProcess* mySigma = new MySigma(); 
      pythia.setSigmaPtr( mySigma); 
</pre> 
If an external phase-space generator should be used for this process 
(see "Implementing an external phase-space generator" below), 
this should be specified as a second argument in the call to 
<code>setSigmaPtr()</code>, as in: 
<pre> 
      pythia.setSigmaPtr( new mySigma(), new myPhaseSpaceGenerator() ); 
</pre> 
If you have several processes you can repeat the procedure any number 
of times. When <code>pythia.init()</code> is called these processes 
are initialized along with any internal processes you may have switched on, 
and treated in exactly the same manner. The  <code>pythia.next()</code> 
will therefore generate a mix of the different kinds of processes without 
distinction. See also the <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>Program Flow</a> 
description. 
 
<p/> 
If the code should be of good quality and general usefulness, it would 
be simple to include it as a permanently available process in the 
standard program distribution. The final step of that integration ought to 
be left for the PYTHIA authors, but here is a description of what is 
required. 
 
<p/> 
A flag has to be defined, that allows the process to be switched on; 
by default it should always be off. The name of the flag should be 
chosen of the type <code>model:process</code>. Here the 
<code>model</code> would be related to the general scenario considered, 
e.g. <code>Compositeness</code>, while <code>process</code> would 
specify instate and outstate, separated by a 2 (= to), e.g. 
<code>ug2u*g</code>. 
When several processes are implemented and "belong together" it is 
also useful to define a <code>model:all</code> switch that affects 
all the separate processes. 
 
<p/> 
The flags should normally be stored in the <code>ProcessSelection.xml</code> 
file or one of its daughters for a specific kind of processes. This is to 
make them easily found by users. You could create and use your own 
<code>.xml</code> file, so long as you then add that name to the 
list of files in the <code>Index.xml</code> file. (If not, 
the flags would never be created and the program would not work.) 
 
<p/> 
In the <code>ProcessContainer.c</code> file, the 
<code>SetupContainers::init()</code> method needs to be expanded to 
create instances of the processes switched on. This code is fairly 
repetitive, and should be easy to copy and modify from the code 
already there. The basic structure is 
<br/>(i) check whether a process is requested by the user and, if so, 
<br/>(ii) create an instance of the matrix-element class, 
<br/>(iii)create a container for the matrix element and its associated 
phase-space handling, and 
<br>(iv) add the container to the existing process list. 
 
<p/> 
Two minor variations are possible. One is that a set of related 
processes are lumped inside the the same initial check, i.e. are 
switched on all together. The second is that the matrix-element 
constructor may take arguments, as specified by you (see above). 
If so, the same basic matrix element may be recycled for a set of 
related processes, e.g. one for a composite <i>u</i> and one for 
a composite <i>d</i>. Obviously these variations may be combined. 
 
<h3>Implementing an external phase-space generator</h3> 
 
An external phase-space generator can be interfaced by encapsulating 
it within a class inheriting from PYTHIA's <code>PhaseSpace</code> 
base class. The following three virtual methods must be defined: 
<pre> 
      // Determine how phase space should be sampled. 
      virtual bool setupSampling(); 
      // Select a trial event kinematics. 
      virtual bool trialKin(bool inEvent = true, bool repeatSame = false); 
      // Construct final (accepted) event kinematics. 
      virtual bool finalKin(); 
</pre> 
Optionally, a further virtual method is available to specify whether 
beam particles are resolved in partons or scatter directly, 
<pre> 
      // Inform whether beam particles are resolved or scatter directly. 
      virtual bool isResolved(); 
</pre> 
with default return value <code>true</code>. 
 
<p/> 
In the <code>setupSampling()</code> step the main point is to determine 
the upper estimate of the cross section integrated over the allowed 
phase space regions, and this should be stored in <code>sigmaMx</code>. 
The ratio between the correct cross section and its upper estimate 
is a measure of the phase-space selection efficiency, and the purpose 
of this step is to optimize the sampling accordingly. To this end any 
convenient set of phase-space variables may be chosen. The 
<code>x1H</code> and <code>x2H</code> varables should be used to 
denote the incoming parton momentum fractions, however, to be used in 
PDF evaluations. 
 
<p/> 
In the <code>trialKin()</code> intermediate step the same set of internal 
variables can be used, and fed into the <code>SigmaProcess</code> code 
to evaluate the cross section in the given phase space point, multiplied 
by the integrated cross section. This value is to be stored in 
<code>sigmaNw</code>, and the ratio <code>sigmaNw/sigmaMx</code> 
will be used to determine whether the trial event is accepted or not. 
 
<p/> 
In the <code>finalKin()</code> step the output is more standardized. 
The key values are the ones stored in the <code>mH[]</code> and 
<code>pH[]</code> arrays, the former for masses and the latter for 
four-momenta. Here the first two slots represent the two incoming 
partons and the subsequent ones up to ten outgoing particles. Other 
particle properties, like the number of final-state particles, their 
identities and colours, and more, are defined by the 
<code>SigmaProcess</code> class. 
 
<p/> 
A tailor-made <i>2 &rarr; 3</i> 
generator could be defined, e.g., by starting from the code for PYTHIA's 
internal <code>PhaseSpace2to3tauycyl</code> base class, which provides 
a specific representation of 3-parton phase space, used for generic 
<i>2 &rarr; 3</i> processes in PYTHIA. The virtual functions 
described could then be redefined to generate a different sampling of 
3-parton phase space. One example of this is provided by the existing 
<code>PhaseSpace2to3yyycyl</code> class, which PYTHIA uses for 
massless QCD processes. Note the interplay between the phase-space 
variables, generated and saved here, and how they are used by the 
matrix-element codes. For general processes, the user can define 
samplings in terms of their own phase-space parametrizations, as long 
as the corresponding matrix elements use the same variables to 
evaluate the cross-section expressions. 
 
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
$data = "SigmaProcess:cMassiveME = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "SigmaProcess:bMassiveME = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SigmaProcess:muMassiveME = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "SigmaProcess:tauMassiveME = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
