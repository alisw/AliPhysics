<html>
<head>
<title>Matching and Merging</title>
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

<form method='post' action='MatchingAndMerging.php'>
 
<h2>Matching and Merging</h2> 
<ol id="toc">
  <li><a href="#section0">MC@NLO, jet matching, multi-jet merging and NLO merging with </a></li>
  <li><a href="#section1">Implementing an external ME+PS combination scheme and interfacing this </a></li>
</ol>

 
Starting from a Born-level leading-order (LO) process, higher orders 
can be included in various ways. The three basic approaches would be 
<ul> 
<li>A formal order-by-order perturbative calculation, in each order 
higher including graphs both with one particle more in the final 
state and with one loop more in the intermediate state. This is 
accurate to the order of the calculation, but gives no hint of 
event structures beyond that, with more particles in the final state. 
Today next-to-leading order (NLO) is standard, while 
next-to-next-to-leading order (NNLO) is coming. This approach 
thus is limited to few orders, and also breaks down in soft and 
collinear regions, which makes it unsuitable for matching to 
hadronization. 
</li> 
<li>Real emissions to several higher orders, but neglecting the 
virtual/loop corrections that should go with it at any given order. 
Thereby it is possible to allow for topologies with a large and 
varying number of partons, at the prize of not being accurate to any 
particular order. The approach also opens up for doublecounting, 
and as above breaks down in soft and colliner regions. 
</li> 
<li>The parton shower provides an approximation to higher orders, 
both real and virtual contributions for the emission of arbitrarily 
many particles. As such it is less accurate than either of the two 
above, at least for topologies of well separated partons, but it 
contains a physically sensible behaviour in the soft and collinear 
limits, and therefore matches well onto the hadronization stage. 
</li> 
</ul> 
Given the pros and cons, much of the effort in recent years has 
involved the development of different prescriptions to combine 
the methods above in various ways. 
 
<p/> 
The common traits of all combination methods are that matrix elements 
are used to describe the production of hard and well separated 
particles, and parton showers for the production of soft or collinear 
particles. What differs between the various approaches that have been 
proposed are which matrix elements are being used, how doublecounting 
is avoided, and how the transition from the hard to the soft regime 
is handled. These combination methods are typically referred to as 
"matching" or "merging" algorithms. There is some confusion about 
the distinction between the two terms, and so we leave it to the 
inventor/implementor of a particular scheme to choose and motivate 
the name given to that scheme. 
 
<p/> 
PYTHIA comes with methods, to be described next, that implement 
or support several different kind of algorithms. The field is 
open-ended, however: any external program can feed in 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches events</a> that 
PYTHIA subsequently showers, adds multiparton interactions to, 
and hadronizes. These events afterwards can be reweighted and 
combined in any desired way. The maximum <i>pT</i> of the shower 
evolution is set by the Les Houches <code>scale</code>, on the one 
hand, and by the values of the <code>SpaceShower:pTmaxMatch</code>, 
<code>TimeShower:pTmaxMatch</code> and other parton-shower settings, 
on the other. Typically it is not possible to achieve perfect 
matching this way, given that the PYTHIA <i>pT</i> evolution 
variables are not likely to agree with the variables used for cuts 
in the external program. Often one can get close enough with simple 
means but, for an improved matching, 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> can be inserted to control 
the steps taken on the way, e.g. to veto those parton shower branchings 
that would doublecount emissions included in the matrix elements. 
 
<p/> 
Zooming in from the "anything goes" perspective, the list of relevent 
approaches actively supported is as follows. 
<ul> 
 
<li>For many/most resonance decays the first branching in the shower is 
merged with first-order matrix elements [<a href="Bibliography.php#refBen87" target="page">Ben87</a>, <a href="Bibliography.php#refNor01" target="page">Nor01</a>]. This 
means that the emission rate is accurate to NLO, similarly to the POWHEG 
strategy (see below), but built into the 
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>timelike showers</a>. 
The angular orientation of the event after the first emission is only 
handled by the parton shower kinematics, however. Needless to say, 
this formalism is precisely what is tested by <i>Z^0</i> decays at 
LEP1, and it is known to do a pretty good job there. 
</li> 
 
<li>Also the <?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>spacelike showers</a> 
contain a correction to first-order matrix elements, but only for the 
one-body-final-state processes 
<i>q qbar &rarr; gamma^*/Z^0/W^+-/h^0/H^0/A0/Z'0/W'+-/R0</i> 
[<a href="Bibliography.php#refMiu99" target="page">Miu99</a>] and <i>g g &rarr; h^0/H^0/A0</i>, and only to 
leading order. That is, it is equivalent to the POWHEG formalism for 
the real emission, but the prefactor "cross section normalization" 
is LO rather than NLO. Therefore this framework is less relevant, 
and has been superseded the following ones. 
</li> 
 
<li>The POWHEG strategy [<a href="Bibliography.php#refNas04" target="page">Nas04</a>] provides a cross section 
accurate to NLO. The hardest emission is constructed with unit 
probability, based on the ratio of the real-emission matrix element 
to the Born-level cross section, and with a Sudakov factor derived 
from this ratio, i.e. the philosophy introduced in [<a href="Bibliography.php#refBen87" target="page">Ben87</a>]. 
<br/>While POWHEG is a generic strategy, the POWHEG BOX 
[<a href="Bibliography.php#refAli10" target="page">Ali10</a>] is an explicit framework, within which several 
processes are available. The code required for merging the PYTHIA 
showers with POWHEG input can be found in 
<code>include/Pythia8Plugins/PowHegHooks.h</code>, and is further 
described on a <?php $filepath = $_GET["filepath"];
echo "<a href='POWHEGMerging.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
A user example is found in <code>examples/main31</code>. 
</li> 
 
<li>The other traditional approach for NLO calculations is the 
MC@NLO one [<a href="Bibliography.php#refFri02" target="page">Fri02</a>]. In it the shower emission probability, 
without its Sudakov factor, is subtracted from the real-emission 
matrix element to regularize divergences. It therefore requires a 
analytic knowledge of the way the shower populates phase space. 
The aMC@NLO package [<a href="Bibliography.php#refFre11" target="page">Fre11</a>] offers an implementation for 
PYTHIA 8, developed by Paolo Torrielli and Stefano Frixione. The 
global-recoil option of the PYTHIA final-state shower has been 
constructed to be used for the above-mentioned subtraction. 
</li> 
 
<li>Multi-jet merging in the CKKW-L approach [<a href="Bibliography.php#refLon01" target="page">Lon01</a>] 
is directly available. Its implementation, relevant parameters 
and test programs are documented on a 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
</li> 
 
<li>Multi-jet matching in the MLM approach [<a href="Bibliography.php#refMan02" target="page">Man02</a>, <a href="Bibliography.php#refMan07" target="page">Man07</a>] 
is also available, either based on the ALPGEN or on the Madgraph 
variant, and with input events either from ALPGEN or from 
Madgraph. For details see 
<?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
</li> 
 
<li>Unitarised matrix element + parton shower merging (UMEPS) 
is directly available. Its implementation, relevant parameters 
and test programs are documented on a 
<?php $filepath = $_GET["filepath"];
echo "<a href='UMEPSMerging.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
</li> 
 
<li>Next-to-leading order multi-jet merging (in the NL3 and UNLOPS approaches) 
is directly available. Its implementation, relevant parameters 
and test programs are documented on a 
<?php $filepath = $_GET["filepath"];
echo "<a href='NLOMerging.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
</li> 
 
<li>Next-to-leading order jet matching in the FxFx approach 
is also available. For details see 
<?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>separate page</a>. 
</li> 
 
</ul> 
 
<br/><hr/> 
<a name="section0"></a> 
<h3>MC@NLO, jet matching, multi-jet merging and NLO merging with 
main89.cc</h3> 
 
A common Pythia main program for MC@NLO NLO+PS matching, MLM jet matching, 
FxFx (NLO) jet matching, CKKW-L merging, UMEPS merging and UNLOPS (NLO) 
merging is available through <code>main89.cc</code>, together with the input 
files <code>main89mlm.cmnd</code>, <code>main89fxfx.cmnd</code>, 
<code>main89ckkwl.cmnd</code>, <code>main89umeps.cmnd</code> and 
<code>main89unlops.cmnd</code>. The interface to MLM jet matching relies 
on MadGraph, while all other options of <code>main89.cc</code> use aMC@NLO 
input. 
 
<code>main89.cc</code> produces HepMC events [<a href="Bibliography.php#refDob01" target="page">Dob01</a>], that can be 
histogrammed (e.g. using RIVET [<a href="Bibliography.php#refBuc10" target="page">Buc10</a>]), or used as input for a 
detector simulation. If the user is not familiar with HepMC analysis tools, it 
is possible to instead use Pythia's histogramming routines. For this, remove 
the lines referring to HepMC, and histogram events as illustrated (for CKKW-L) 
for the histogram <i>histPTFirstSum</i> in <code>main84.cc</code>, i.e. 
using <i>weight*normhepmc</i> as weight. 
 
<p/> 
All settings can be transferred to <code>main89.cc</code> through an input 
file. The input file is part of the command line input of 
<code>main89.cc</code>, i.e. you can execute <code>main89</code> with the 
command 
<p/> 
<code>./main89 myInputFile.cmnd myhepmc.hepmc</code> 
<p/> 
 
to read the input <code>myInputFile.cmnd</code> and produce the output file 
<code>myhepmc.hepmc</code> . Since <code>main89.cc</code> is currently a 
"front-end" for different types of matching/merging, we will briefly discuss 
the inputs for this sample program in the following. 
 
<h4>Inputs</h4> 
 
In its current form, <code>main89.cc</code> uses LHEF input to transfer 
(weighted) phase space points to Pythia. It is possible to include all 
parton multiplicities in one LHEF sample. If e.g. UMEPS merging for 
W-boson + up to two additional partons is to be performed, one LHE file 
containing W+zero, W+one and W+two parton events is required. 
 
<p/> 
All input settings are handed to <code>main89.cc</code> in the form of an 
input file. We have included the input settings files 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>main89mlm.cmnd</code>, which 
illustrates the MLM jet matching interface, 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>main89ckkwl.cmnd</code>, which 
illustrates the CKKW-L multi-jet merging interface, 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp;  <code>main89umeps.cmnd</code>, which 
illustrates the UMEPS multi-jet merging interface, and 
 <p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>main89fxfx.cmnd</code>, which 
illustrates the FxFx NLO jet matching interface, 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp;  <code>main89unlops.cmnd</code>, which 
illustrates the UNLOPS multi-jet NLO merging interface. 
<p/> 
Other settings (e.g. using <code>main89.cc</code> as simple LO+PS or as MC@NLO 
interface) are of course possible. In the following, we will briefly explain 
how input for the five choices above are generated and handled. 
 
<h4>MLM jet matching with main89.cc</h4> 
 
For MLM jet matching, <code>main89.cc</code> currently relies on LHEF input 
from MadGraph. Due to the particular unweighting strategy performed in the 
generation of these inputs, the sample program starts by estimating the 
cross section. After this estimate, MLM jet matching within the Madgraph 
approach is performed in a second Pythia run. Example MLM settings can be 
found in <code>main89mlm.cmnd</code>. Please consult 
<?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>Jet Matching</a> for  details. 
 
<h4>CKKW-L merging with main89.cc</h4> 
 
For CKKW-L merging, <code>main89.cc</code> currently relies on LHEF inputs 
generated with the leading-order mode of aMC@NLO (i.e. events should 
be generated with <code>./bin/generate_events aMC@LO</code>). 
No run to estimate the cross section estimate is needed. Example CKKW-L 
settings can be found in <code>main89ckkwl.cmnd</code>. Please consult 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a> for  details. 
 
<h4>UMEPS merging with main89.cc</h4> 
 
For UMEPS merging, <code>main89.cc</code> currently relies on LHEF inputs 
generated with the leading-order mode of aMC@NLO as well (see above). 
<code>main89.cc</code> automatically assigns if an event will be used as 
"standard" event or as "subtractive" contribution. Example UMEPS 
settings can be found in <code>main89umeps.cmnd</code>. Please 
consult <?php $filepath = $_GET["filepath"];
echo "<a href='UMEPSMerging.php?filepath=".$filepath."' target='page'>";?>UMEPS merging</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a> for  details. 
 
<h4>FxFx (NLO) jet matching with main89.cc</h4> 
 
For FxFx jet matching, <code>main89.cc</code> relies on MC@NLO input LHE 
files generated with aMC@NLO. To produce FxFx outputs in aMC@NLO, the settings 
<code>PYTHIA8  = parton_shower</code>, <code>3 = ickkw</code> and 
<code>x = ptj</code> are necessary in your aMC@NLO run card. Here, 
<code>x</code> is the value of the matching scale in FxFx, i.e. has be 
identical to <code>JetMatching:qCutME</code> in the Pythia inputs. 
Example FxFx settings for Pythia can be found in <code>main89fxfx.cmnd</code>. 
Please consult <?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>Jet Matching</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='aMCatNLOMatching.php?filepath=".$filepath."' target='page'>";?>aMC@NLO matching</a> for  details. 
 
 
<h4>UNLOPS (NLO) merging with main89.cc</h4> 
 
For UNLOPS merging, <code>main89.cc</code> currently relies on LHEF inputs 
generated with the aMC@NLO. The UNLOPS interface in <code>main89.cc</code> 
requires a) leading-order inputs generated with the leading-order mode of 
aMC@NLO, using the UNLOPS prescription, and b) next-to-leading-order inputs 
generated with the NLO mode of aMC@NLO, using the UNLOPS prescription. 
To produce UNLOPS outputs in aMC@NLO, the settings 
<code>PYTHIA8  = parton_shower</code>, <code>4 = ickkw</code> and 
<code>x = ptj</code> are necessary in your aMC@NLO run card. Here, 
<code>x</code> is the value of the merging scale in UNLOPS, i.e. 
has be identical to <code>Merging:TMS</code> in the Pythia inputs. 
<code>main89.cc</code> will then process NLO inputs and LO inputs 
consecutively, and will automatically assign if an event will be used as 
"standard" event or as "subtractive" contribution. Example UNLOPS 
settings can be found in <code>main89umeps.cmnd</code>. Please 
consult <?php $filepath = $_GET["filepath"];
echo "<a href='UNLOPSMerging.php?filepath=".$filepath."' target='page'>";?>UMEPS merging</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a> for  details. 
 
<br/><br/><hr/> 
<a name="section1"></a> 
<h3>Implementing an external ME+PS combination scheme and interfacing this 
plugin with Pythia</h3> 
 
For experts and developers of new matching/merging schemes, Pythia also offers 
the possibility to completely replace its internal merging machinery with 
a user-defined plugin code (much in the same way that parton shower plugins 
(cf. <?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>Implement New Showers</a>) are 
possible). This allows for maximum flexibility while still benefiting from 
the full Pythia event generation machinery. Note that the ME+PS merging with 
the VINCIA and DIRE shower plugins make use of this flexibility, and might 
thus provide helpful clarifications. 
 
Of course, implementing your 
own, new matching/merging scheme is a non-trivial task, and comprehensive 
guidelines on how to proceed are impossible to set. However, it is important 
that an external matching/merging plugin interfaces to Pythia in a simple and 
well-defined manner. Here, we will document which C++ functions are necessary 
to be able to use an external matching/merging (MM) plugin within Pythia. 
 
<p/> 
To understand how to design a MM plugin for Pythia, it is useful to review 
how Pythia's internal merging machinery is structured. The interaction between 
the core Pythia and the merging code is governed by the 
<code>Merging</code> and <code>MergingHooks</code> classes. Note that for 
moderately complex requirements, it may be sufficient to only replace Pythia's 
instance of <code>MergingHooks</code> with a pointer to an external class (cf. 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a>). The latter two classes 
are supplemented with the helper classes <code>History</code> and 
<code>HardProcess</code>. The latter gathers information on the (user-supplied 
information about the) hard core scattering process to which hard jets are 
added 
by ME+PS merging. It is only used as a helper to the <code>MergingHooks</code> 
class. The <code>History</code> class contains the implementation of all 
internal (LO or NLO) merging schemes. The <code>Merging</code> class 
acts as a bridge between the implementation in the <code>History</code> class 
and the rest of the Pythia code. 
 
<p/> 
To implement an external MM plugin, you will have to write classes that derive 
from the <code>Merging</code>, <code>MergingHooks</code> and 
<code>HardProcess</code> classes of Pythia. For special cases, it might also 
be permissible to only implement a replacement of the <code>Merging</code> 
class, while still using Pythia's implementation of the other two classes. 
 
The external MM plugin can then be transferred to and 
used by Pythia much in the same way as <code>UserHooks</code> classes or 
shower plugins. More concretely, an external MM code will 
be used if a pointer to an instance of the external classes is transferred to 
Pythia via the methods 
 
<a name="anchor1"></a>
<p/><strong> Pythia::setMergingPtr( Merging* myMerging) &nbsp;</strong> <br/>
   
 
<a name="anchor2"></a>
<p/><strong> Pythia::setMergingHooksPtr( MergingHooks* myMergingHooks) &nbsp;</strong> <br/>
   
 
<a name="anchor3"></a>
<p/><strong> MergingHooks::setHardProcessPtr( HardProcess* myHardProcess) &nbsp;</strong> <br/>
   
 
<p/> 
The option to only use a user-defined <code>MergingHooks</code> instance is 
already documented in the item <?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a> 
and will not be discussed further. We will now focus on how to implement 
external <code>Merging</code>, <code>MergingHooks</code> and 
<code>HardProcess</code> classes that can be used as a complete 
replacement of the Pythia methods. 
 
Let us assume that you want to create a class of type <code>MyMerging</code>, 
and you call its instance <code>myMerging</code>. For this external ME+PS 
merging class to be interfaced to Pythia, the class needs to inherit from the 
<code>Pythia8::Merging</code> base class. It is further necessary to define 
the following functions that serve as interface to Pythia: 
 
<a name="anchor4"></a>
<p/><strong> virtual ~MyMerging() &nbsp;</strong> <br/>
A destructor for your ME+PS class. If not defined, the base class's empty 
destructor will be used. 
   
 
<a name="anchor5"></a>
<p/><strong> virtual void MyMerging::init() &nbsp;</strong> <br/>
A method that is used to initialize your merging class. Pythia will call 
this function during its initialization and after all pointers to 
internal classes (e.g. to instances of the <code>Info</code> and 
<code>ParticleData</code> classes) have been set up. 
   
 
<a name="anchor6"></a>
<p/><strong> virtual void MyMerging::statistics() &nbsp;</strong> <br/>
This function can be used to collect and print merging information at the 
end of the event generation. Pythia will call this function in the execution 
of a <code>Pythia::stat()</code> call. 
   
 
<a name="anchor7"></a>
<p/><strong> virtual int MyMerging::mergeProcess( Event& process) &nbsp;</strong> <br/>
This function should be the main interface of Pythia to the MM plugin. 
Pythia will execute this function once the partonic (fixed-order) scattering 
has been constructed (or read from LHEF). The partonic scattering is 
transferred via the <code>process</code> argument. The external MM plugin 
should then, based on the <code>process</code>, implement the matching/merging 
strategy. It is permissible that this function changes <code>process</code>. 
In this case, Pythia will continue the event generation with the changed 
<code>process</code> as starting point. The return value of the function 
steers how Pythia should proceed after the function execution. The following 
return values are supported: 
<li> -1 : Reject the event and exit the generation/processing of the current 
event </li> 
<li>   0: Reject the event but continue with the generation/processing of the 
current event.</li> 
<li>   1: Keep the event but continue with the generation/processing of the 
current event.</li> 
<li>   2: Reject the event but continue with the generation/processing of the 
   current event. However, re-evaluate resonance decays before any other 
   event generation step. This option can be necessary if the merging code 
   removes or changes resonant particles from <code>process</code>.</li> 
 
Note that because this function is the main interface between the MM plugin 
and Pythia, it is necessary to use this function to set up all the information 
that you might later need (merging weights, particle counters, etc) in this 
call already. 
   
 
<p/> 
For  details on how to design your <code>MyMerging</code> class, and to 
understand the interface to Pythia, studying Pythia's internal code is 
unavoidable. Each potential developer of a MM plugin should do so. 
 
<p/> 
The other main ingredient of the interface to MM plugins is a new 
implementation 
of the <code>MergingHooks</code> class. Let us assume that you want to 
create a class of type <code>MyMergingHooks</code>, and you call its 
instance <code>myMergingHooks</code>. For this class to be interfaced to 
Pythia, it will need to inherit from the <code>Pythia8::MergingHooks</code> 
base class. 
 
<a name="anchor8"></a>
<p/><strong> virtual ~MyMergingHooks() &nbsp;</strong> <br/>
A destructor for your MergingHooks class. If not defined, the base class's 
empty destructor will be used. 
   
 
<a name="anchor9"></a>
<p/><strong> virtual void MyMergingHooks::init() &nbsp;</strong> <br/>
A method that is used to initialize your <code>MyMergingHooks</code> class. 
Pythia will call this function during its initialization and after all 
pointers to internal classes (e.g. to instances of the <code>Info</code> and 
<code>ParticleData</code> classes) have been set up. 
   
 
<a name="anchor10"></a>
<p/><strong> virtual bool MyMergingHooks::canVetoStep() &nbsp;</strong> <br/>
This function will be used to tell Pythia if a CKKW-L-style event veto 
after the first parton shower emission should be checked. If so, the function 
should return true, and false otherwise. 
   
 
<a name="anchor11"></a>
<p/><strong> virtual bool MyMergingHooks::doVetoStep( const Event& process, const Event& event, bool doResonance = false ) &nbsp;</strong> <br/>
This function will be used to implement the check of a CKKW-L-style event veto 
after the first parton shower emission, i.e. to check if the first parton 
shower emission is above the merging scale. 
If the input event <code>event</code> 
after emission should be kept, then false should be returned. If you want 
instead to veto the event and continue with a completely now hard scattering 
event, true should be returned. 
   
 
<a name="anchor12"></a>
<p/><strong> virtual bool MyMergingHooks::canVetoEmission() &nbsp;</strong> <br/>
This function will be used to tell Pythia if a veto of emissions should 
potentially be applied. 
   
 
<a name="anchor13"></a>
<p/><strong> virtual bool MyMergingHooks::doVetoStep( const Event& event) &nbsp;</strong> <br/>
This function will be used to implement the check if shower emissions should 
be discarded, as e.g. necessary in UMEPS or UNLOPS merging. 
You can study the input event <code>event</code> after emission, and return 
true if the emission is valid, and false if you want to reject the emission. 
Note that this veto does not lead to event rejections, only in potentially 
removing certain emissions during shower evolution. 
   
 
<a name="anchor14"></a>
<p/><strong> virtual bool MyMergingHooks::setShowerStartingScales( bool     isTrial, bool doMergeFirstEmm,     double& pTscaleIn, const Event& event,     double& pTmaxFSRIn, bool& limitPTmaxFSRin,     double& pTmaxISRIn, bool& limitPTmaxISRin,     double& pTmaxMPIIn, bool& limitPTmaxMPIin ) &nbsp;</strong> <br/>
This function allows to set the starting scales for timelike and spacelike 
showering as well as multiparton interactions. It is thus necessary to 
properly start trial showers (that generate necessary no-emission 
probabilities), and for setting the correct starting conditions for parton 
showering of accepted (non-zero weight) events. 
The input <code>event</code> gives the hard process before showers and MPI 
are attempted. 
If <code>isTrial=true</code>, this means that the function is currently called 
from within a trial shower object (to produce no-emission probabilities). If 
<code>doMergeFirstEmm=true</code>, then the function is called to set starting 
conditions for the shower evolution of an (accepted) event. The double 
arguments <code>pTscaleIn</code>, <code>pTmaxFSRIn</code>, 
<code>pTmaxISRIn</code> and <code>pTmaxMPIIn</code> are tentative values 
for the starting scales of FSR, ISR and MPI. The function may overwrite 
these with the desired values. Similarly, <code>limitPTmaxFSRin</code>, 
<code>limitPTmaxFSRin</code> and <code>limitPTmaxMPIin</code> inform Pythia 
if the phase space for FSR/ISR/MPI is restricted (true) or unrestricted 
(false). Again, the tentative values can be overwritten. 
   
 
<p/> 
The <code>MergingHooks</code> base class allows for further virtual functions 
that are not directly called by Pythia, and are hence 
not necessary to define. Th usage of these functions within 
Pythia's <code>Merging</code> and <code>History</code> classes is documented 
in <?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L merging</a>. The additional (optional) 
virtual functions are: 
 
<a name="anchor15"></a>
<p/><strong>   virtual double dampenIfFailCuts( const Event& inEvent ) &nbsp;</strong> <br/>
   
<a name="anchor16"></a>
<p/><strong>   virtual bool canCutOnRecState() &nbsp;</strong> <br/>
   
<a name="anchor17"></a>
<p/><strong>   virtual bool doCutOnRecState( const Event& event ) &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<p/><strong>   virtual bool canVetoTrialEmission() &nbsp;</strong> <br/>
   
<a name="anchor19"></a>
<p/><strong>   virtual bool doVetoTrialEmission( const Event&, const Event& ) &nbsp;</strong> <br/>
   
<a name="anchor20"></a>
<p/><strong>   virtual double hardProcessME( const Event& inEvent ) &nbsp;</strong> <br/>
   
<a name="anchor21"></a>
<p/><strong> virtual double tmsDefinition( const Event& event) &nbsp;</strong> <br/>
   
<a name="anchor22"></a>
<p/><strong> virtual int getNumberOfClusteringSteps(const Event& event, bool resetNjetMax = false) &nbsp;</strong> <br/>
   
<a name="anchor23"></a>
<p/><strong>   virtual bool useShowerPlugin() &nbsp;</strong> <br/>
   
 
<p/> 
The internal implementation of <code>MergingHooks</code> in Pythia heavily 
relies on the <code>HardProcess</code> helper class. It is in principle 
not necessary to follow the same strategy when implementing a derived 
<code>MyMergingHooks</code> class. However, to benefit from the Pythia 
implementation, and to allow for a structure similar to the internal code also 
for an external MM plugin, it is also possible to effectively replace (in the 
<code>MergingHooks</code> class) the pointer to an instance of 
<code>HardProcess</code> with a pointer to an external implementation. 
Let us assume that you want to create a class of type 
<code>MyHardProcess</code>, and you call its instance 
<code>myHardProcess</code>. For this class to be interfaced to 
<code>MergingHooks</code> (or the derived <code>MyMergingHooks</code> class), 
it will need to inherit from the <code>Pythia8::HardProcess</code> base class. 
 
<a name="anchor24"></a>
<p/><strong> virtual ~MyHardProcess() &nbsp;</strong> <br/>
A destructor for your HardProcess class. If not defined, the base class's 
empty destructor will be used. 
   
 
<a name="anchor25"></a>
<p/><strong> virtual void MyHardProcess::initOnProcess( string process, ParticleData* particleData) &nbsp;</strong> <br/>
This function can be used to initialize the instance of your HardProcess 
class. In the internal Pythia implementation, this acts as a wrapper around 
the next function. 
   
 
<a name="anchor26"></a>
<p/><strong> virtual void MyHardProcess::translateProcessString( string process) &nbsp;</strong> <br/>
This function will use the string argument to set up the hard process 
bookkeeping, e.g. how many incoming/outgoing particles of which flavour are 
contained in the core (lowest multiplicity) scattering process. 
   
 
<a name="anchor27"></a>
<p/><strong> virtual void MyHardProcess::storeCandidates( const Event& event, string process) &nbsp;</strong> <br/>
This function studies the input event and book-keeps the particles that 
may be considered as part of the core scattering process. For this, it may 
use the four next functions. 
   
 
<a name="anchor28"></a>
<p/><strong> virtual bool MyHardProcess::allowCandidates(int iPos, vector&lt;int&gt; Pos1, vector&lt;int&gt; Pos2, const Event& event) &nbsp;</strong> <br/>
This function uses the input vectors of positions of particles in the input 
event to decide if the particle with <code>iPos</code> could be member 
of the core scattering. If the particle with position <code>iPos</code> 
cannot be part of the core scattering (e.g. because it is a final state 
parton, while the core scattering contains final state leptons only), then 
the function should return false. Else, return true to allow this candidate. 
Note that it might be possible to find multiple equally good core scattering 
candidates. In this case, all candidates should be found (with the 
<code>findOtherCandidates</code> function), and can be potentially be 
replaced (with <code>exchangeCandidates</code>). 
   
 
<a name="anchor29"></a>
<p/><strong> virtual bool MyHardProcess::matchesAnyOutgoing(int iPos, const Event& event) &nbsp;</strong> <br/>
This function may be used to check if the particle with position 
<code>iPos</code> in the input event should be considered an outgoing particle 
of the core scattering. 
   
 
<a name="anchor30"></a>
<p/><strong> virtual bool MyHardProcess::findOtherCandidates(int iPos, const Event& event, bool doReplace) &nbsp;</strong> <br/>
The argument <code>iPos</code> specifies the position of a particle in the 
input event which is tagged as part of the core scattering. This function may 
be used to check the role of <code>iPos</code> as  core scattering member may 
be filled by another particle in the event record. If so, and if 
<code>doReplace=true</code>, then <code>iPos</code> will no longer be 
book-kept as part of the core scattering. An example where this functionality 
is helpful is if the input event is g g -> b b~ b  b~, and the core scattering 
is g g -> b b~. Not swapping the hard process candidates could in this case 
mean that not all parton shower histories can be found. 
The function should return false if no replacements can be found, and true 
otherwise. 
   
 
<a name="anchor31"></a>
<p/><strong> virtual bool MyHardProcess::exchangeCandidates( vector&lt;int&gt; candidates1, vector&lt;int&gt; candidates2, map&lt;int, int&gt; further1, map&lt;int, int&gt; further2) &nbsp;</strong> <br/>
This function implements the replacement of a list of core scattering 
candidates by another list of candidates. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
