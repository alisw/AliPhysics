\ProvidesClass{my_report}
\LoadClass[11pt, titlepage]{article}
\usepackage[paperwidth=21cm,paperheight=29.7cm,includehead,headheight=2.0cm,pdftex,hmargin={3cm,2.5cm},vmargin={0cm,2cm},]{geometry} 
\setlength{\parindent}{0cm} 
\renewcommand{\baselinestretch}{1.2}
\parskip=1em

\renewcommand{\contentsname}{Table of Contents}

\renewcommand{\normalsize}{\fontsize{11pt}{11pt}\selectfont}
\renewcommand{\title}[1]{\textsf{\Huge{#1}}}
\renewcommand{\author}[1]{\textsf{\Huge{#1}}}
\newcommand{\subtitle}[1]{\textsf{\textbf{\Large{#1}}}}
\newcommand{\abstractheading}[1]{\textsf{\textbf{\LARGE{#1}}}}
\newcommand{\code}[1]{\texttt{\footnotesize{#1}}}

\renewcommand{\section}{\@startsection
{section}%                    % the name
{1}%                          % the level
{0mm}%                        % the indent
{10.6mm}%             % the beforeskip
{4.2mm}%           % the afterskip
{\LARGE\bfseries\sffamily}}  % the style

\renewcommand{\subsection}{\@startsection
{subsection}%                    % the name
{2}%                          % the level
{0mm}%                        % the indent
{6.4mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\Large\bfseries\sffamily}}  % the style

\renewcommand{\subsubsection}{\@startsection
{subsubsection}%                    % the name
{3}%                          % the level
{0mm}%                        % the indent
{4.2mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\normalsize\bfseries\sffamily}}  % the style

%% Chapter headings should be centered, uppercase, and at the top of the page.
\def\@makechapterhead#1{%
  { \parindent \z@ \raggedright \normalfont
    %\centering
    \ifnum \c@secnumdepth >\m@ne
        \huge\textsf{\@chapapp\space \thechapter:}
        % \par\nobreak
        %\vskip 20\p@
    \fi
    \interlinepenalty\@M
    \huge \bfseries \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@
  } }
  
\def\@makeschapterhead#1{%
  %\vspace*{50\p@}%
  { \parindent \z@ \raggedright
    %\centering
    \normalfont
    \interlinepenalty\@M
    \huge \bfseries  \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@

  }}
  
 \renewenvironment{abstract}{%
      \chapter*{\abstractname}%
      \addcontentsline{toc}{chapter}{\abstractname}
 }
     
     
\makeatletter
\renewcommand{\l@section}{\bfseries\@dottedtocline{1}{0em}{2.3em}}
\renewcommand{\l@subsection}{\normalfont\@dottedtocline{2}{2em}{2.3em}}
\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{3}{2em}{2.3em}}
%\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{4}{2em}{2.3em}}
\makeatother

\def\maketitle{\begin{titlepage}
\thispagestyle{empty}
\let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}
\null
\vfil
\begin{flushright}
{\bf  CERN-PH-TH/2010-261 \\
  IFJPAN-IV-2010-6  \\
  November 2010}
\end{flushright}

\begin{center}
\title{PHOTOS Interface in C++\\[0.5cm]{} 
\textbf{\small Technical and Physics Documentation  \\[0.4cm]\rule{4cm}{1pt}\\[0.4cm]}}

\author{\LARGE N. Davidson$^{a,b}$, T. Przedzinski$^c$, Z. Was$^{a,d}$}\\[0.3cm]
\textsf{\Large{\em $^a$ Institute of Nuclear Physics, 
Polish Academy of Sciences,\\
         ul. Radzikowskiego 152, 31-342 Cracow, Poland}\\
{\em $^b$ The University of Melbourne, School of Physics \\ Australia  }\\
{\em $^c $ The Faculty of Physics, Astronomy and Applied Computer
Science,\\
Jagellonian University, Reymonta 4, 30-059 Cracow, Poland}\\
{\em $^d $ Theory Group, Physics Department, CERN, CH-1211, Geneva 23, Switzerland}
 }\\[0.3cm]

%\includegraphics[height=3cm]{logo.eps}\\[1cm]
{\bf \Large Abstract} \vskip 2 mm 

\end{center}

{ The first version of {\tt PHOTOS} Monte Carlo for
bremsstrahlung in the decay of particles and resonances
 with an interface to the HepMC event record written in C++ is now available.
The main purpose of the present
paper is to document technical aspects of the {\tt PHOTOS} Monte Carlo 
installation and its use. A multitude of test results and examples are  
distributed together 
with the program code.

The
{\tt PHOTOS} C++ physics precision is  now as good as that of its {\tt FORTRAN}
predecessor. However
better steering options are available. 
An algorithm of the event record interface
is prepared
for the  installation of process dependent variants of the 
photon emission kernel.
Weights, featuring complete first order matrix elements, can be installed 
for general use. In the {\tt FORTRAN} version of  {\tt PHOTOS} they were 
available only for decays of particles   at rest and with spin 
 set along the $z$ axis.

Physics assumptions used in the program and properties of the
solution it offers are reviewed. 
In particular, it is  mentioned
that the  second order matrix elements were used in design and 
validation of the program iteration procedure. Also it is explained
that the phase space parameterization used in the program is exact. 
 }

%\centerline{\normalsize \today}

\vfil

\footnoterule
\noindent
{\footnotesize
\begin{itemize}
\item[${\dag}$]
This work is partially supported by  EU Marie Curie Research Training Network 
grant under the contract No. MRTN-CT-2006-0355505 and by  
Polish Government grant N202 06434 (2008-2011). 

\end{itemize}
}
\null
\end{titlepage}}

\pagestyle{plain}

