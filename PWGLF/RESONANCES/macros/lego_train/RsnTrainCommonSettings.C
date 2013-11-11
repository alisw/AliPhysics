


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=10">
        <title>AliRsn/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C at master · mvala/AliRsn</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="github-fe122-cp1-prd.iad.github.net">
    <meta name="ruby" content="ruby 1.9.3p194-tcs-github-tcmalloc (0e75de19f8) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035/">
    <link rel="xhr-socket" href="/_sockets" />
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="C1CC021A:3449:4500338:5280EDA9" name="octolytics-dimension-request_id" /><meta content="1196645" name="octolytics-actor-id" /><meta content="fbellini" name="octolytics-actor-login" /><meta content="bf838c8338d2455a4df60d7927dd86428a1c42e3fbe7e39c1cc5a76b68ac07d4" name="octolytics-actor-hash" />
    

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="6ZIeHwigqjkkfN7l2WnviLJHqpaRCaQbupWawq/qxbM=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-804556dba6658262abda18880c76c8b30304dcb3.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-fc46856b3ad905365f892bf463b1e1c0fd84764e.css" media="all" rel="stylesheet" type="text/css" />
    

    

      <script src="https://github.global.ssl.fastly.net/assets/frameworks-bca527bb59d94c16d6bf2a759779d7953fa41e76.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-76d84f68034a6ab0a216cca231fe065459f00d52.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="7e309922140cc9fe6a44cd4384c566d3">

        <link data-pjax-transient rel='permalink' href='/mvala/AliRsn/blob/ac7b83b65e557df3b1dc5717527bda685d13ce68/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C'>
  <meta property="og:title" content="AliRsn"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/mvala/AliRsn"/>
  <meta property="og:image" content="https://github.global.ssl.fastly.net/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="Contribute to AliRsn development by creating an account on GitHub."/>

  <meta name="description" content="Contribute to AliRsn development by creating an account on GitHub." />

  <meta content="65018" name="octolytics-dimension-user_id" /><meta content="mvala" name="octolytics-dimension-user_login" /><meta content="3375377" name="octolytics-dimension-repository_id" /><meta content="mvala/AliRsn" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="3375377" name="octolytics-dimension-repository_network_root_id" /><meta content="mvala/AliRsn" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/mvala/AliRsn/commits/master.atom" rel="alternate" title="Recent Commits to AliRsn:master" type="application/atom+xml" />

  </head>


  <body class="logged_in  env-production macintosh vis-public  page-blob">
    <div class="wrapper">
      
      
      
      


      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/">
  <span class="mega-octicon octicon-mark-github"></span>
</a>

    
    <a href="/notifications" class="notification-indicator tooltipped downwards" data-gotokey="n" title="You have no unread notifications">
        <span class="mail-status all-read"></span>
</a>

      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="fbellini"
      data-repo="mvala/AliRsn"
      data-branch="master"
      data-sha="0cce44223d5c1abe691b0fb0edcfda61ac193522"
  >

    <input type="hidden" name="nwo" value="mvala/AliRsn" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    


  <ul id="user-links">
    <li>
      <a href="/fbellini" class="name">
        <img height="20" src="https://0.gravatar.com/avatar/44a5e372e8b9c3d9383af156beb94922?d=https%3A%2F%2Fidenticons.github.com%2Fcc6feba04e67c09b04800c5bcdd1acda.png&amp;r=x&amp;s=140" width="20" /> fbellini
      </a>
    </li>

      <li>
        <a href="/new" id="new_repo" class="tooltipped downwards" title="Create a new repo" aria-label="Create a new repo">
          <span class="octicon octicon-repo-create"></span>
        </a>
      </li>

      <li>
        <a href="/settings/profile" id="account_settings"
          class="tooltipped downwards"
          aria-label="Account settings "
          title="Account settings ">
          <span class="octicon octicon-tools"></span>
        </a>
      </li>
      <li>
        <a class="tooltipped downwards" href="/logout" data-method="post" id="logout" title="Sign out" aria-label="Sign out">
          <span class="octicon octicon-log-out"></span>
        </a>
      </li>

  </ul>

<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo-create"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>



    <li class="section-title">
      <span title="mvala/AliRsn">This repository</span>
    </li>
      <li>
        <a href="/mvala/AliRsn/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
      </li>
</ul>

</div>


    
  </div>
</div>

      

      




          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="6ZIeHwigqjkkfN7l2WnviLJHqpaRCaQbupWawq/qxbM=" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="3375377" />

    <div class="select-menu js-menu-container js-select-menu">
      <a class="social-count js-social-count" href="/mvala/AliRsn/watchers">
        4
      </a>
      <span class="minibutton select-menu-button with-count js-menu-target" role="button" tabindex="0">
        <span class="js-select-button">
          <span class="octicon octicon-eye-unwatch"></span>
          Unwatch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-remove-close js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container" role="menu">

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for discussions in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-watch"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-unwatch"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
  

  <div class="js-toggler-container js-social-container starring-container on">
    <a href="/mvala/AliRsn/unstar"
      class="minibutton with-count js-toggler-target star-button starred upwards"
      title="Unstar this repository" data-remote="true" data-method="post" rel="nofollow">
      <span class="octicon octicon-star-delete"></span><span class="text">Unstar</span>
    </a>

    <a href="/mvala/AliRsn/star"
      class="minibutton with-count js-toggler-target star-button unstarred upwards"
      title="Star this repository" data-remote="true" data-method="post" rel="nofollow">
      <span class="octicon octicon-star"></span><span class="text">Star</span>
    </a>

      <a class="social-count js-social-count" href="/mvala/AliRsn/stargazers">
        2
      </a>
  </div>

  </li>


        <li>
          <a href="/mvala/AliRsn/fork" class="minibutton with-count js-toggler-target fork-button lighter upwards" title="Fork this repo" rel="nofollow" data-method="post">
            <span class="octicon octicon-git-branch-create"></span><span class="text">Fork</span>
          </a>
          <a href="/mvala/AliRsn/network" class="social-count">1</a>
        </li>


</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author">
            <a href="/mvala" class="url fn" itemprop="url" rel="author"><span itemprop="title">mvala</span></a>
          </span>
          <span class="repohead-name-divider">/</span>
          <strong><a href="/mvala/AliRsn" class="js-current-repository js-repo-home-link">AliRsn</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container ">

        <div class="repository-sidebar">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped leftwards" title="Code">
        <a href="/mvala/AliRsn" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /mvala/AliRsn">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped leftwards" title="Issues">
          <a href="/mvala/AliRsn/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-gotokey="i" data-selected-links="repo_issues /mvala/AliRsn/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped leftwards" title="Pull Requests"><a href="/mvala/AliRsn/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /mvala/AliRsn/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/mvala/AliRsn/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_wiki /mvala/AliRsn/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/mvala/AliRsn/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /mvala/AliRsn/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/mvala/AliRsn/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /mvala/AliRsn/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/mvala/AliRsn/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /mvala/AliRsn/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


  </div>
</div>

            <div class="only-with-full-nav">
              

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/mvala/AliRsn.git" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/mvala/AliRsn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><strong>SSH</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="git@github.com:mvala/AliRsn.git" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="git@github.com:mvala/AliRsn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/mvala/AliRsn" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/mvala/AliRsn" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <span class="octicon help tooltipped upwards" title="Get help on which URL is right for you.">
    <a href="https://help.github.com/articles/which-remote-url-should-i-use">
    <span class="octicon octicon-question"></span>
    </a>
  </span>
</p>

  <a href="github-mac://openRepo/https://github.com/mvala/AliRsn" data-url="github-mac://openRepo/https://github.com/mvala/AliRsn" class="minibutton sidebar-button js-conduit-rewrite-url">
    <span class="octicon octicon-device-desktop"></span>
    Clone in Desktop
  </a>


              <a href="/mvala/AliRsn/archive/master.zip"
                 class="minibutton sidebar-button"
                 title="Download this repository as a zip file"
                 rel="nofollow">
                <span class="octicon octicon-cloud-download"></span>
                Download ZIP
              </a>
            </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:f463b68cd9b735999566fef85a4dd59c -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/mvala/AliRsn/find/master" data-pjax data-hotkey="t" class="js-show-file-finder" style="display:none">Show File Finder</a>

<div class="file-navigation">
  
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    role="button" aria-label="Switch branches or tags" tabindex="0">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/mvala/AliRsn/blob/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/mvala/AliRsn/blob/testSigmaCut/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C"
                 data-name="testSigmaCut"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="testSigmaCut">testSigmaCut</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/mvala/AliRsn/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="6ZIeHwigqjkkfN7l2WnviLJHqpaRCaQbupWawq/qxbM=" /></div>
            <span class="octicon octicon-git-branch-create select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master" />
            <input type="hidden" name="path" id="branch" value="PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C" />
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">AliRsn</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">PWGLF</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">RESONANCES</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES/macros" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">macros</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES/macros/lego_train" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">lego_train</span></a></span><span class="separator"> / </span><strong class="final-path">RsnTrainCommonSettings.C</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>



  <div class="commit file-history-tease">
    <img class="main-avatar" height="24" src="https://0.gravatar.com/avatar/981bd3d2462b14d60e2f7b8aaf38f067?d=https%3A%2F%2Fidenticons.github.com%2F5876221c2bec11337e186aac8aa3100d.png&amp;r=x&amp;s=140" width="24" />
    <span class="author"><a href="/mvala" rel="author">mvala</a></span>
    <time class="js-relative-date" datetime="2013-10-21T14:00:02-07:00" title="2013-10-21 14:00:02">October 21, 2013</time>
    <div class="commit-title">
        <a href="/mvala/AliRsn/commit/c57f913f821919e5352f25927acf9e6c52975f66" class="message" data-pjax="true" title="Added QA config for Phi">Added QA config for Phi</a>
    </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong> contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img height="24" src="https://0.gravatar.com/avatar/981bd3d2462b14d60e2f7b8aaf38f067?d=https%3A%2F%2Fidenticons.github.com%2F5876221c2bec11337e186aac8aa3100d.png&amp;r=x&amp;s=140" width="24" />
            <a href="/mvala">mvala</a>
          </li>
      </ul>
    </div>
  </div>

<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>62 lines (49 sloc)</span>
        <span>2.863 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
            <a class="minibutton tooltipped leftwards js-conduit-openfile-check"
               href="github-mac://openRepo/https://github.com/mvala/AliRsn?branch=master&amp;filepath=PWGLF%2FRESONANCES%2Fmacros%2Flego_train%2FRsnTrainCommonSettings.C"
               data-url="github-mac://openRepo/https://github.com/mvala/AliRsn?branch=master&amp;filepath=PWGLF%2FRESONANCES%2Fmacros%2Flego_train%2FRsnTrainCommonSettings.C"
               title="Open this file in GitHub for Mac"
               data-failed-title="Your version of GitHub for Mac is too old to open this file. Try checking for updates.">
                <span class="octicon octicon-device-desktop"></span> Open
            </a>
                <a class="minibutton"
                   href="/mvala/AliRsn/edit/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/mvala/AliRsn/raw/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C" class="button minibutton " id="raw-url">Raw</a>
            <a href="/mvala/AliRsn/blame/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C" class="button minibutton ">Blame</a>
          <a href="/mvala/AliRsn/commits/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
          <a class="minibutton danger empty-icon tooltipped downwards"
             href="/mvala/AliRsn/delete/master/PWGLF/RESONANCES/macros/lego_train/RsnTrainCommonSettings.C"
             title=""
             data-method="post" data-test-id="delete-blob-file" rel="nofollow">
          Delete
        </a>
      </div><!-- /.actions -->

    </div>
        <div class="blob-wrapper data type-c js-blob-data">
        <table class="file-code file-diff">
          <tr class="file-code-line">
            <td class="blob-line-nums">
              <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>

            </td>
            <td class="blob-line-code">
                    <div class="highlight"><pre><div class='line' id='LC1'><span class="kt">void</span> <span class="nf">RsnTrainCommonSettings</span><span class="p">(</span><span class="n">TString</span> <span class="n">type</span><span class="p">,</span><span class="n">TString</span> <span class="n">rsnPart</span><span class="p">,</span><span class="n">TString</span> <span class="n">extraMacro</span><span class="o">=</span><span class="s">&quot;&quot;</span><span class="p">,</span><span class="n">TString</span> <span class="n">extraMacroArgs</span><span class="o">=</span><span class="s">&quot;&quot;</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC2'><br/></div><div class='line' id='LC3'>&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">valid</span><span class="p">;</span></div><div class='line' id='LC4'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalStr</span><span class="p">(</span><span class="s">&quot;LegoTrainPath&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC5'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">valid</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC6'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">legoTrainPath</span> <span class="o">=</span> <span class="s">&quot;$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train&quot;</span><span class="p">;</span></div><div class='line' id='LC7'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalStr</span><span class="p">(</span><span class="s">&quot;LegoTrainPath&quot;</span><span class="p">,</span><span class="n">legoTrainPath</span><span class="p">.</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC8'>&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC9'><br/></div><div class='line' id='LC10'>&nbsp;&nbsp;&nbsp;<span class="c1">// removing Option part fo Rsn particle</span></div><div class='line' id='LC11'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">rsnPart</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;:&quot;</span><span class="p">))</span> <span class="n">rsnPart</span><span class="p">.</span><span class="n">Remove</span><span class="p">(</span><span class="n">rsnPart</span><span class="p">.</span><span class="n">Index</span><span class="p">(</span><span class="s">&quot;:&quot;</span><span class="p">),</span><span class="n">rsnPart</span><span class="p">.</span><span class="n">Length</span><span class="p">());</span></div><div class='line' id='LC12'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalStr</span><span class="p">(</span><span class="s">&quot;RsnParticle&quot;</span><span class="p">,</span><span class="n">rsnPart</span><span class="p">.</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC13'><br/></div><div class='line' id='LC14'>&nbsp;&nbsp;&nbsp;<span class="c1">// CollisionType (pp=0,PbPb=1,pPb=2)</span></div><div class='line' id='LC15'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;pp&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsCollisionType&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC16'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;PbPb&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsCollisionType&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC17'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;pPb&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsCollisionType&quot;</span><span class="p">,</span><span class="mi">2</span><span class="p">);</span></div><div class='line' id='LC18'><br/></div><div class='line' id='LC19'>&nbsp;&nbsp;&nbsp;<span class="c1">// data type</span></div><div class='line' id='LC20'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;ESD&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsESD&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC21'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsESD&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC22'><br/></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;<span class="c1">// flag if we are using MC</span></div><div class='line' id='LC24'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;MC&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsMC&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC25'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsMC&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC26'><br/></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;<span class="c1">// flag if we want to use event Mixing</span></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;MIX&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsMixing&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsMixing&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC30'><br/></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;<span class="c1">// Use Rsn Mini</span></div><div class='line' id='LC32'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">type</span><span class="p">.</span><span class="n">Contains</span><span class="p">(</span><span class="s">&quot;MINI&quot;</span><span class="p">))</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsRsnMini&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsRsnMini&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC34'><br/></div><div class='line' id='LC35'><br/></div><div class='line' id='LC36'>&nbsp;&nbsp;&nbsp;<span class="c1">// current RSN base defaults (Will be changed in future)</span></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">extraMacro</span><span class="p">.</span><span class="n">IsNull</span><span class="p">())</span> <span class="p">{</span></div><div class='line' id='LC38'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">extraMacro</span><span class="p">.</span><span class="n">ReplaceAll</span><span class="p">(</span><span class="s">&quot;.C&quot;</span><span class="p">,</span><span class="s">&quot;&quot;</span><span class="p">);</span></div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Printf</span><span class="p">(</span><span class="s">&quot;Running Extra Macro %s(%s)&quot;</span><span class="p">,</span><span class="n">extraMacro</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span><span class="n">extraMacroArgs</span><span class="p">.</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC40'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">gROOT</span><span class="o">-&gt;</span><span class="n">ProcessLine</span><span class="p">(</span><span class="n">TString</span><span class="o">::</span><span class="n">Format</span><span class="p">(</span><span class="s">&quot;%s(%s)&quot;</span><span class="p">,</span><span class="n">extraMacro</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span><span class="n">extraMacroArgs</span><span class="p">.</span><span class="n">Data</span><span class="p">()).</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnQA&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnNumMix&quot;</span><span class="p">,</span><span class="mi">5</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalDbl</span><span class="p">(</span><span class="s">&quot;RsnEventCutPrimaryVertex&quot;</span><span class="p">,</span><span class="mf">10.0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC45'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalStr</span><span class="p">(</span><span class="s">&quot;RsnLegoTrainCommonCutOption&quot;</span><span class="p">,</span><span class="s">&quot;mon&quot;</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC46'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnPhysSelFilterBit&quot;</span><span class="p">,</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnCommonQualityCut&quot;</span><span class="p">,</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseRapidity&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnOutputFull&quot;</span><span class="p">,</span><span class="mi">1</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseMCMomentum&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseMCMonitoring&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseAOD049Patch&quot;</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC53'><br/></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalDbl</span><span class="p">(</span><span class="s">&quot;RsnMixDiffMult&quot;</span><span class="p">,</span><span class="mf">10.0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC55'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalDbl</span><span class="p">(</span><span class="s">&quot;RsnMixDiffVz&quot;</span><span class="p">,</span><span class="mf">1.0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalDbl</span><span class="p">(</span><span class="s">&quot;RsnMixDiffAngle&quot;</span><span class="p">,</span><span class="o">-</span><span class="mf">1.0</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC57'><br/></div><div class='line' id='LC58'>&nbsp;&nbsp;&nbsp;<span class="c1">// expert options (don&#39;t change)</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">SetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnMixPrintRefresh&quot;</span><span class="p">,</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span><span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC60'><br/></div><div class='line' id='LC61'><span class="p">}</span></div></pre></div>
            </td>
          </tr>
        </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.05265s from github-fe122-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/mvala/AliRsn/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

  </body>
</html>

