


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=10">
        <title>AliRsn/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C at master · mvala/AliRsn</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="github-fe123-cp1-prd.iad.github.net">
    <meta name="ruby" content="ruby 1.9.3p194-tcs-github-tcmalloc (0e75de19f8) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035/">
    <link rel="xhr-socket" href="/_sockets" />
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="C1CC021A:3449:44FEE7B:5280ED9D" name="octolytics-dimension-request_id" /><meta content="1196645" name="octolytics-actor-id" /><meta content="fbellini" name="octolytics-actor-login" /><meta content="bf838c8338d2455a4df60d7927dd86428a1c42e3fbe7e39c1cc5a76b68ac07d4" name="octolytics-actor-hash" />
    

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="6ZIeHwigqjkkfN7l2WnviLJHqpaRCaQbupWawq/qxbM=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-804556dba6658262abda18880c76c8b30304dcb3.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-fc46856b3ad905365f892bf463b1e1c0fd84764e.css" media="all" rel="stylesheet" type="text/css" />
    

    

      <script src="https://github.global.ssl.fastly.net/assets/frameworks-bca527bb59d94c16d6bf2a759779d7953fa41e76.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-76d84f68034a6ab0a216cca231fe065459f00d52.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="7e309922140cc9fe6a44cd4384c566d3">

        <link data-pjax-transient rel='permalink' href='/mvala/AliRsn/blob/ac7b83b65e557df3b1dc5717527bda685d13ce68/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C'>
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
          


<!-- blob contrib key: blob_contributors:v21:5d2606dbde125e94811da33be0ae9d09 -->

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
              <a href="/mvala/AliRsn/blob/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/mvala/AliRsn/blob/testSigmaCut/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C"
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
            <input type="hidden" name="path" id="branch" value="PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C" />
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
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">AliRsn</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">PWGLF</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">RESONANCES</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES/macros" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">macros</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/mvala/AliRsn/tree/master/PWGLF/RESONANCES/macros/lego_train" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">lego_train</span></a></span><span class="separator"> / </span><strong class="final-path">AddRsnPairsPhi.C</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  <div class="commit commit-loader file-history-tease js-deferred-content" data-url="/mvala/AliRsn/contributors/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C">
    Fetching contributors…

    <div class="participation">
      <p class="loader-loading"><img alt="Octocat-spinner-32-eaf2f5" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" /></p>
      <p class="loader-error">Cannot retrieve contributors at this time</p>
    </div>
  </div>

<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>287 lines (242 sloc)</span>
        <span>12.083 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
            <a class="minibutton tooltipped leftwards js-conduit-openfile-check"
               href="github-mac://openRepo/https://github.com/mvala/AliRsn?branch=master&amp;filepath=PWGLF%2FRESONANCES%2Fmacros%2Flego_train%2FAddRsnPairsPhi.C"
               data-url="github-mac://openRepo/https://github.com/mvala/AliRsn?branch=master&amp;filepath=PWGLF%2FRESONANCES%2Fmacros%2Flego_train%2FAddRsnPairsPhi.C"
               title="Open this file in GitHub for Mac"
               data-failed-title="Your version of GitHub for Mac is too old to open this file. Try checking for updates.">
                <span class="octicon octicon-device-desktop"></span> Open
            </a>
                <a class="minibutton"
                   href="/mvala/AliRsn/edit/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/mvala/AliRsn/raw/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C" class="button minibutton " id="raw-url">Raw</a>
            <a href="/mvala/AliRsn/blame/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C" class="button minibutton ">Blame</a>
          <a href="/mvala/AliRsn/commits/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
          <a class="minibutton danger empty-icon tooltipped downwards"
             href="/mvala/AliRsn/delete/master/PWGLF/RESONANCES/macros/lego_train/AddRsnPairsPhi.C"
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
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>

            </td>
            <td class="blob-line-code">
                    <div class="highlight"><pre><div class='line' id='LC1'><span class="cp">#ifndef __CINT__</span></div><div class='line' id='LC2'><span class="cp">#include &lt;AliRsnListOutput.h&gt;</span></div><div class='line' id='LC3'><span class="cp">#include &lt;AliRsnValuePair.h&gt;</span></div><div class='line' id='LC4'><span class="cp">#include &lt;AliRsnValueEvent.h&gt;</span></div><div class='line' id='LC5'><span class="cp">#include &lt;AliRsnMiniAnalysisTask.h&gt;</span></div><div class='line' id='LC6'><span class="cp">#endif</span></div><div class='line' id='LC7'><br/></div><div class='line' id='LC8'><span class="kt">void</span> <span class="nf">AddRsnPairsPhi</span><span class="p">(</span><span class="n">AliAnalysisTaskSE</span> <span class="o">*</span><span class="n">task</span><span class="p">,</span></div><div class='line' id='LC9'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">isMC</span><span class="p">,</span></div><div class='line' id='LC10'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">isMixing</span><span class="p">,</span></div><div class='line' id='LC11'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliPID</span><span class="o">::</span><span class="n">EParticleType</span> <span class="n">pType1</span><span class="p">,</span></div><div class='line' id='LC12'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">listID1</span><span class="p">,</span></div><div class='line' id='LC13'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliPID</span><span class="o">::</span><span class="n">EParticleType</span> <span class="n">pType2</span><span class="p">,</span></div><div class='line' id='LC14'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">listID2</span><span class="p">,</span></div><div class='line' id='LC15'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnCutSet</span> <span class="o">*</span><span class="n">cutsEvent</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span></div><div class='line' id='LC16'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnCutSet</span> <span class="o">*</span><span class="n">cutsPair</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span></div><div class='line' id='LC17'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">suffix</span> <span class="o">=</span> <span class="s">&quot;&quot;</span><span class="p">)</span></div><div class='line' id='LC18'><span class="p">{</span></div><div class='line' id='LC19'><br/></div><div class='line' id='LC20'>&nbsp;&nbsp;&nbsp;<span class="n">Printf</span><span class="p">(</span><span class="s">&quot;id1=%d id2=%d&quot;</span><span class="p">,</span><span class="n">listID1</span><span class="p">,</span><span class="n">listID2</span><span class="p">);</span></div><div class='line' id='LC21'><br/></div><div class='line' id='LC22'>&nbsp;&nbsp;&nbsp;<span class="c1">// retrieve mass from PDG database</span></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span>         <span class="n">pdg</span>  <span class="o">=</span> <span class="mi">333</span><span class="p">;</span></div><div class='line' id='LC24'>&nbsp;&nbsp;&nbsp;<span class="n">TDatabasePDG</span> <span class="o">*</span><span class="n">db</span>   <span class="o">=</span> <span class="n">TDatabasePDG</span><span class="o">::</span><span class="n">Instance</span><span class="p">();</span></div><div class='line' id='LC25'>&nbsp;&nbsp;&nbsp;<span class="n">TParticlePDG</span> <span class="o">*</span><span class="n">part</span> <span class="o">=</span> <span class="n">db</span><span class="o">-&gt;</span><span class="n">GetParticle</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC26'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">mass</span> <span class="o">=</span> <span class="n">part</span><span class="o">-&gt;</span><span class="n">Mass</span><span class="p">();</span></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">valid</span><span class="p">;</span></div><div class='line' id='LC28'><br/></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">isRsnMini</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsRsnMini&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC30'><br/></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isRsnMini</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC32'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AddPairOutputMiniPhi</span><span class="p">(</span><span class="n">task</span><span class="p">,</span><span class="n">isMC</span><span class="p">,</span><span class="n">isMixing</span><span class="p">,</span><span class="n">pType1</span><span class="p">,</span><span class="n">listID1</span><span class="p">,</span><span class="n">pType2</span><span class="p">,</span><span class="n">listID2</span><span class="p">,</span><span class="n">pdg</span><span class="p">,</span><span class="n">mass</span><span class="p">,</span><span class="n">cutsPair</span><span class="p">,</span><span class="n">suffix</span><span class="p">);</span></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// this function is common and it is located in RsnConfig.C</span></div><div class='line' id='LC35'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// as ouptup AddPairOutputPhi from this macro will be taken</span></div><div class='line' id='LC36'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AddPair</span><span class="p">(</span><span class="n">task</span><span class="p">,</span><span class="n">isMC</span><span class="p">,</span><span class="n">isMixing</span><span class="p">,</span><span class="n">pType1</span><span class="p">,</span><span class="n">listID1</span><span class="p">,</span><span class="n">pType2</span><span class="p">,</span><span class="n">listID2</span><span class="p">,</span><span class="n">pdg</span><span class="p">,</span><span class="n">mass</span><span class="p">,</span><span class="n">cutsEvent</span><span class="p">,</span><span class="n">cutsPair</span><span class="p">,</span><span class="n">suffix</span><span class="p">);</span></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC38'><span class="p">}</span></div><div class='line' id='LC39'><span class="kt">void</span> <span class="nf">AddPairOutputPhi</span><span class="p">(</span><span class="n">AliRsnLoopPair</span> <span class="o">*</span><span class="n">pair</span><span class="p">)</span></div><div class='line' id='LC40'><span class="p">{</span></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">valid</span><span class="p">;</span></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">isFullOutput</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnOutputFull&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">collisionType</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsCollisionType&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">useRapidity</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseRapidity&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC45'><br/></div><div class='line' id='LC46'>&nbsp;&nbsp;&nbsp;<span class="c1">// axes</span></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnValuePair</span> <span class="o">*</span><span class="n">axisIM</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnValuePair</span><span class="p">(</span><span class="s">&quot;IM&quot;</span><span class="p">,</span> <span class="n">AliRsnValuePair</span><span class="o">::</span><span class="n">kInvMass</span><span class="p">);</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnValuePair</span> <span class="o">*</span><span class="n">axisPt</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnValuePair</span><span class="p">(</span><span class="s">&quot;PT&quot;</span><span class="p">,</span> <span class="n">AliRsnValuePair</span><span class="o">::</span><span class="n">kPt</span><span class="p">);</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnValuePair</span> <span class="o">*</span><span class="n">axisEta</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnValuePair</span><span class="p">(</span><span class="s">&quot;ETA&quot;</span><span class="p">,</span> <span class="n">AliRsnValuePair</span><span class="o">::</span><span class="n">kEta</span><span class="p">);</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnValuePair</span> <span class="o">*</span><span class="n">axisRapidity</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnValuePair</span><span class="p">(</span><span class="s">&quot;Y&quot;</span><span class="p">,</span> <span class="n">AliRsnValuePair</span><span class="o">::</span><span class="n">kY</span><span class="p">);</span></div><div class='line' id='LC51'><br/></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnValueEvent</span> <span class="o">*</span><span class="n">axisCentrality</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC53'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span> <span class="n">axisCentrality</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnValueEvent</span><span class="p">(</span><span class="s">&quot;MULTI&quot;</span><span class="p">,</span><span class="n">AliRsnValueEvent</span><span class="o">::</span><span class="n">kCentralityV0</span><span class="p">);</span></div><div class='line' id='LC54'><br/></div><div class='line' id='LC55'><br/></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;<span class="n">axisIM</span>     <span class="o">-&gt;</span><span class="n">SetBins</span><span class="p">(</span><span class="mi">300</span><span class="p">,</span> <span class="mf">0.9</span><span class="p">,</span> <span class="mf">1.2</span><span class="p">);</span></div><div class='line' id='LC57'><span class="c1">//   axisIM     -&gt;SetBins(600, 0.9, 1.2);</span></div><div class='line' id='LC58'><span class="c1">//   axisIM     -&gt;SetBins(1000, 0.9, 1.9);</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;<span class="n">axisPt</span>     <span class="o">-&gt;</span><span class="n">SetBins</span><span class="p">(</span><span class="mi">120</span><span class="p">,</span> <span class="mf">0.0</span><span class="p">,</span> <span class="mf">12.0</span><span class="p">);</span></div><div class='line' id='LC60'><span class="c1">//    axisEta    -&gt;SetBins(400, -2.0, 2.0);</span></div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;<span class="n">axisEta</span>    <span class="o">-&gt;</span><span class="n">SetBins</span><span class="p">(</span><span class="mi">400</span><span class="p">,</span> <span class="o">-</span><span class="mf">0.5</span><span class="p">,</span> <span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;<span class="n">axisRapidity</span>    <span class="o">-&gt;</span><span class="n">SetBins</span><span class="p">(</span><span class="mi">10</span><span class="p">,</span> <span class="o">-</span><span class="mf">0.5</span><span class="p">,</span> <span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC63'><br/></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">axisCentrality</span><span class="p">)</span> <span class="n">axisCentrality</span><span class="o">-&gt;</span><span class="n">SetBins</span><span class="p">(</span><span class="mi">20</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">100</span><span class="p">);</span></div><div class='line' id='LC65'><br/></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;<span class="c1">// output: 2D histogram of inv. mass vs. pt</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnListOutput</span> <span class="o">*</span><span class="n">outPair</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">isFullOutput</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC69'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPair</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnListOutput</span><span class="p">(</span><span class="s">&quot;pair&quot;</span><span class="p">,</span> <span class="n">AliRsnListOutput</span><span class="o">::</span><span class="n">kHistoDefault</span><span class="p">);</span></div><div class='line' id='LC70'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisIM</span><span class="p">);</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPair</span> <span class="o">=</span> <span class="k">new</span> <span class="n">AliRsnListOutput</span><span class="p">(</span><span class="s">&quot;pair&quot;</span><span class="p">,</span> <span class="n">AliRsnListOutput</span><span class="o">::</span><span class="n">kHistoSparse</span><span class="p">);</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisIM</span><span class="p">);</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisPt</span><span class="p">);</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">useRapidity</span><span class="p">)</span> <span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisRapidity</span><span class="p">);</span></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisEta</span><span class="p">);</span></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">axisCentrality</span><span class="p">)</span> <span class="n">outPair</span><span class="o">-&gt;</span><span class="n">AddValue</span><span class="p">(</span><span class="n">axisCentrality</span><span class="p">);</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC79'>&nbsp;&nbsp;&nbsp;<span class="c1">// add outputs to loop</span></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;<span class="n">pair</span><span class="o">-&gt;</span><span class="n">AddOutput</span><span class="p">(</span><span class="n">outPair</span><span class="p">);</span></div><div class='line' id='LC81'><span class="p">}</span></div><div class='line' id='LC82'><br/></div><div class='line' id='LC83'><span class="kt">void</span> <span class="nf">AddPairOutputMiniPhi</span><span class="p">(</span><span class="n">AliAnalysisTaskSE</span> <span class="o">*</span><span class="n">task</span><span class="p">,</span> <span class="n">Bool_t</span> <span class="n">isMC</span><span class="p">,</span><span class="n">Bool_t</span> <span class="n">isMixing</span><span class="p">,</span> <span class="n">AliPID</span><span class="o">::</span><span class="n">EParticleType</span> <span class="n">pType1</span><span class="p">,</span><span class="n">Int_t</span> <span class="n">listID1</span><span class="p">,</span> <span class="n">AliPID</span><span class="o">::</span><span class="n">EParticleType</span> <span class="n">pType2</span><span class="p">,</span><span class="n">Int_t</span> <span class="n">listID2</span><span class="p">,</span> <span class="n">Int_t</span> <span class="n">pdgMother</span><span class="p">,</span><span class="n">Double_t</span> <span class="n">massMother</span><span class="p">,</span> <span class="n">AliRsnCutSet</span> <span class="o">*</span><span class="n">cutsPair</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">TString</span> <span class="n">suffix</span> <span class="o">=</span> <span class="s">&quot;&quot;</span><span class="p">)</span></div><div class='line' id='LC84'><span class="p">{</span></div><div class='line' id='LC85'><br/></div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span> <span class="n">valid</span><span class="p">;</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">isFullOutput</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnOutputFull&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">useMixing</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsMixing&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">collisionType</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;IsCollisionType&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC90'><br/></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">useMixLike</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnMixLike&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">useRapidity</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnUseRapidity&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC93'><br/></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">useQA</span> <span class="o">=</span> <span class="n">AliRsnTrainManager</span><span class="o">::</span><span class="n">GetGlobalInt</span><span class="p">(</span><span class="s">&quot;RsnQA&quot;</span><span class="p">,</span><span class="n">valid</span><span class="p">);</span></div><div class='line' id='LC95'><br/></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniAnalysisTask</span> <span class="o">*</span><span class="n">taskRsnMini</span> <span class="o">=</span>  <span class="p">(</span><span class="n">AliRsnMiniAnalysisTask</span> <span class="o">*</span><span class="p">)</span><span class="n">task</span><span class="p">;</span></div><div class='line' id='LC97'><br/></div><div class='line' id='LC98'><br/></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;<span class="c1">// retrieve mass from PDG database</span></div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span>         <span class="n">pdg</span>  <span class="o">=</span> <span class="mi">333</span><span class="p">;</span></div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;<span class="n">TDatabasePDG</span> <span class="o">*</span><span class="n">db</span>   <span class="o">=</span> <span class="n">TDatabasePDG</span><span class="o">::</span><span class="n">Instance</span><span class="p">();</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;<span class="n">TParticlePDG</span> <span class="o">*</span><span class="n">part</span> <span class="o">=</span> <span class="n">db</span><span class="o">-&gt;</span><span class="n">GetParticle</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC103'><br/></div><div class='line' id='LC104'>&nbsp;&nbsp;&nbsp;<span class="c1">// common definitions</span></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">outputType</span> <span class="o">=</span> <span class="s">&quot;HIST&quot;</span><span class="p">;</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isFullOutput</span><span class="p">)</span> <span class="n">outputType</span> <span class="o">=</span> <span class="s">&quot;SPARSE&quot;</span><span class="p">;</span></div><div class='line' id='LC107'><br/></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">useQA</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC109'><br/></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">0</span><span class="p">)</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">UseMultiplicity</span><span class="p">(</span><span class="s">&quot;QUALITY&quot;</span><span class="p">);</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC112'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">UseCentrality</span><span class="p">(</span><span class="s">&quot;V0M&quot;</span><span class="p">);</span></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">multID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kMult</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">outMult</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="s">&quot;eventMult&quot;</span><span class="p">,</span> <span class="s">&quot;HIST&quot;</span><span class="p">,</span> <span class="s">&quot;EVENT&quot;</span><span class="p">);</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMult</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">multID</span><span class="p">,</span> <span class="mi">100</span><span class="p">,</span> <span class="mf">0.0</span><span class="p">,</span> <span class="mf">100.0</span><span class="p">);</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">paID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kPlaneAngle</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">outPa</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="s">&quot;planeAngle&quot;</span><span class="p">,</span> <span class="s">&quot;HIST&quot;</span><span class="p">,</span> <span class="s">&quot;EVENT&quot;</span><span class="p">);</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outPa</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">paID</span><span class="p">,</span> <span class="mi">100</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="n">TMath</span><span class="o">::</span><span class="n">Pi</span><span class="p">());</span></div><div class='line' id='LC119'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC120'><br/></div><div class='line' id='LC121'><br/></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* invariant mass   */</span> <span class="n">Int_t</span> <span class="n">imID</span>   <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kInvMass</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* IM resolution    */</span> <span class="n">Int_t</span> <span class="n">resID</span>  <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kInvMassDiff</span><span class="p">,</span> <span class="n">kTRUE</span><span class="p">);</span></div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* transv. momentum */</span> <span class="n">Int_t</span> <span class="n">ptID</span>   <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kPt</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* centrality       */</span> <span class="n">Int_t</span> <span class="n">centID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kMult</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC126'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* eta              */</span> <span class="n">Int_t</span> <span class="n">etaID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kEta</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC127'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/* rapidity         */</span> <span class="n">Int_t</span> <span class="n">yID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kY</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC128'><br/></div><div class='line' id='LC129'><br/></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// use an array for more compact writing, which are different on mixing and charges</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// [0] = unlike</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// [1] = mixing</span></div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// [2] = like ++</span></div><div class='line' id='LC134'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// [3] = like --</span></div><div class='line' id='LC135'><br/></div><div class='line' id='LC136'><span class="c1">//   Bool_t useMixLike = kFALSE;</span></div><div class='line' id='LC137'><span class="c1">//   useMixLike = kTRUE;</span></div><div class='line' id='LC138'><br/></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">const</span> <span class="n">Int_t</span> <span class="n">numAxis</span> <span class="o">=</span> <span class="mi">7</span><span class="p">;</span></div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span>  <span class="n">use</span>     <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span> <span class="mi">1</span>      <span class="p">,</span>  <span class="n">useMixing</span>      <span class="p">,</span>  <span class="mi">1</span>      <span class="p">,</span>  <span class="mi">1</span>      <span class="p">,</span>  <span class="n">isMC</span>  <span class="p">,</span>  <span class="n">useMixLike</span><span class="p">,</span>  <span class="n">useMixLike</span><span class="p">};</span></div><div class='line' id='LC141'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">name</span>    <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="s">&quot;Unlike&quot;</span><span class="p">,</span> <span class="s">&quot;Mixing&quot;</span><span class="p">,</span> <span class="s">&quot;LikePP&quot;</span><span class="p">,</span> <span class="s">&quot;LikeMM&quot;</span><span class="p">,</span> <span class="s">&quot;Trues&quot;</span><span class="p">,</span> <span class="s">&quot;MixingPP&quot;</span><span class="p">,</span> <span class="s">&quot;MixingMM&quot;</span><span class="p">};</span></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">comp</span>    <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="s">&quot;PAIR&quot;</span>  <span class="p">,</span> <span class="s">&quot;MIX&quot;</span>   <span class="p">,</span> <span class="s">&quot;PAIR&quot;</span>  <span class="p">,</span> <span class="s">&quot;PAIR&quot;</span>  <span class="p">,</span> <span class="s">&quot;TRUE&quot;</span> <span class="p">,</span> <span class="s">&quot;MIX&quot;</span>   <span class="p">,</span> <span class="s">&quot;MIX&quot;</span>   <span class="p">};</span></div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Char_t</span>  <span class="n">charge1</span> <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>    <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>     <span class="p">};</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Char_t</span>  <span class="n">charge2</span> <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="sc">&#39;-&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>    <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>     <span class="p">};</span></div><div class='line' id='LC145'><br/></div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nIM</span>   <span class="o">=</span> <span class="mi">215</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minIM</span>   <span class="o">=</span> <span class="mf">0.985</span><span class="p">,</span> <span class="n">maxIM</span> <span class="o">=</span>  <span class="mf">1.2</span><span class="p">;</span></div><div class='line' id='LC147'><span class="c1">//    Int_t nEta   = 400; Double_t minEta   = -2.0, maxEta =  2.0;</span></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nEta</span>   <span class="o">=</span> <span class="mi">400</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minEta</span>   <span class="o">=</span> <span class="o">-</span><span class="mf">0.5</span><span class="p">,</span> <span class="n">maxEta</span> <span class="o">=</span>  <span class="mf">0.5</span><span class="p">;</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nY</span>   <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minY</span>   <span class="o">=</span> <span class="o">-</span><span class="mf">0.5</span><span class="p">,</span> <span class="n">maxY</span> <span class="o">=</span>  <span class="mf">0.5</span><span class="p">;</span></div><div class='line' id='LC150'><span class="c1">//   Int_t nIM   = 1000; Double_t minIM   = 0.9, maxIM =  1.9;</span></div><div class='line' id='LC151'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nPt</span>   <span class="o">=</span> <span class="mi">240</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minPt</span>   <span class="o">=</span> <span class="mf">0.0</span><span class="p">,</span> <span class="n">maxPt</span> <span class="o">=</span> <span class="mf">24.0</span><span class="p">;</span></div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nCent</span> <span class="o">=</span> <span class="mi">20</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minCent</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">,</span> <span class="n">maxCent</span> <span class="o">=</span> <span class="mf">100.0</span><span class="p">;</span></div><div class='line' id='LC153'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nRes</span>  <span class="o">=</span> <span class="mi">200</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">maxRes</span>  <span class="o">=</span> <span class="mf">0.01</span><span class="p">;</span></div><div class='line' id='LC154'><br/></div><div class='line' id='LC155'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// in case pp</span></div><div class='line' id='LC156'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">0</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC157'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">nIM</span>   <span class="o">=</span> <span class="mi">150</span><span class="p">;</span> <span class="n">minIM</span>  <span class="o">=</span> <span class="mf">0.975</span><span class="p">;</span> <span class="n">maxIM</span>  <span class="o">=</span>  <span class="mf">1.125</span><span class="p">;</span></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">nPt</span>   <span class="o">=</span> <span class="mi">160</span><span class="p">;</span> <span class="n">minPt</span>  <span class="o">=</span> <span class="mf">0.0</span><span class="p">,</span>   <span class="n">maxPt</span>  <span class="o">=</span> <span class="mf">8.0</span><span class="p">;</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">nEta</span>  <span class="o">=</span> <span class="mi">30</span><span class="p">;</span>  <span class="n">minEta</span> <span class="o">=</span> <span class="o">-</span><span class="mf">1.5</span><span class="p">;</span>  <span class="n">maxEta</span> <span class="o">=</span>  <span class="mf">1.5</span><span class="p">;</span></div><div class='line' id='LC160'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">nY</span>    <span class="o">=</span> <span class="mi">3</span><span class="p">;</span>   <span class="n">minY</span>   <span class="o">=</span> <span class="o">-</span><span class="mf">1.5</span><span class="p">;</span>  <span class="n">maxY</span>   <span class="o">=</span>  <span class="mf">1.5</span><span class="p">;</span></div><div class='line' id='LC161'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC162'><br/></div><div class='line' id='LC163'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Printf</span><span class="p">(</span><span class="n">suffix</span><span class="p">.</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// create standard outputs</span></div><div class='line' id='LC165'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">Int_t</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">numAxis</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC166'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">use</span><span class="p">[</span><span class="n">i</span><span class="p">])</span> <span class="k">continue</span><span class="p">;</span></div><div class='line' id='LC167'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// create output</span></div><div class='line' id='LC168'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">out</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="n">Form</span><span class="p">(</span><span class="s">&quot;%s_%s&quot;</span><span class="p">,</span> <span class="n">suffix</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span><span class="n">name</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">Data</span><span class="p">()</span> <span class="p">),</span> <span class="n">outputType</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span> <span class="n">comp</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC169'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// selection settings</span></div><div class='line' id='LC170'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC171'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC173'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC174'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">charge1</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC175'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">charge2</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC176'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetMotherPDG</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetMotherMass</span><span class="p">(</span><span class="n">part</span><span class="o">-&gt;</span><span class="n">Mass</span><span class="p">());</span></div><div class='line' id='LC178'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// pair cuts</span></div><div class='line' id='LC179'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">cutsPair</span><span class="p">)</span> <span class="n">out</span><span class="o">-&gt;</span><span class="n">SetPairCuts</span><span class="p">(</span><span class="n">cutsPair</span><span class="p">);</span></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis X: invmass</span></div><div class='line' id='LC181'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">imID</span><span class="p">,</span> <span class="n">nIM</span><span class="p">,</span> <span class="n">minIM</span><span class="p">,</span> <span class="n">maxIM</span><span class="p">);</span></div><div class='line' id='LC182'><br/></div><div class='line' id='LC183'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isFullOutput</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC184'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Y: transverse momentum</span></div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">ptID</span><span class="p">,</span> <span class="n">nPt</span><span class="p">,</span> <span class="n">minPt</span><span class="p">,</span> <span class="n">maxPt</span><span class="p">);</span></div><div class='line' id='LC186'><br/></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">useRapidity</span><span class="p">)</span> <span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">yID</span><span class="p">,</span> <span class="n">nY</span><span class="p">,</span> <span class="n">minY</span><span class="p">,</span> <span class="n">maxY</span><span class="p">);</span></div><div class='line' id='LC188'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">0</span><span class="p">)</span> <span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">etaID</span><span class="p">,</span> <span class="n">nEta</span><span class="p">,</span> <span class="n">minEta</span><span class="p">,</span> <span class="n">maxEta</span><span class="p">);</span></div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Z: centrality</span></div><div class='line' id='LC190'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">1</span><span class="p">)</span> <span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">centID</span><span class="p">,</span> <span class="n">nCent</span><span class="p">,</span> <span class="n">minCent</span><span class="p">,</span> <span class="n">maxCent</span><span class="p">);</span></div><div class='line' id='LC191'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC192'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC193'><br/></div><div class='line' id='LC194'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// add output for resolution</span></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isMC</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC196'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">outRes</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="n">Form</span><span class="p">(</span><span class="s">&quot;phi_Res%s&quot;</span><span class="p">,</span> <span class="n">suffix</span><span class="p">.</span><span class="n">Data</span><span class="p">()),</span> <span class="n">outputType</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span> <span class="s">&quot;TRUE&quot;</span><span class="p">);</span></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// selection settings</span></div><div class='line' id='LC198'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC199'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC200'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC201'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="sc">&#39;+&#39;</span><span class="p">);</span></div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="sc">&#39;-&#39;</span><span class="p">);</span></div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetMotherPDG</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC205'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetMotherMass</span><span class="p">(</span><span class="n">part</span><span class="o">-&gt;</span><span class="n">Mass</span><span class="p">());</span></div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// pair cuts</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">cutsPair</span><span class="p">)</span> <span class="n">outRes</span><span class="o">-&gt;</span><span class="n">SetPairCuts</span><span class="p">(</span><span class="n">cutsPair</span><span class="p">);</span></div><div class='line' id='LC208'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis X: resolution</span></div><div class='line' id='LC209'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">resID</span><span class="p">,</span> <span class="n">nRes</span><span class="p">,</span> <span class="o">-</span><span class="n">maxRes</span><span class="p">,</span> <span class="n">maxRes</span><span class="p">);</span></div><div class='line' id='LC210'><br/></div><div class='line' id='LC211'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isFullOutput</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC212'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Y: transverse momentum</span></div><div class='line' id='LC213'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outRes</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">ptID</span><span class="p">,</span> <span class="n">nPt</span><span class="p">,</span> <span class="n">minPt</span><span class="p">,</span> <span class="n">maxPt</span><span class="p">);</span></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">useRapidity</span><span class="p">)</span> <span class="n">outRes</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">yID</span><span class="p">,</span> <span class="n">nY</span><span class="p">,</span> <span class="n">minY</span><span class="p">,</span> <span class="n">maxY</span><span class="p">);</span></div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">0</span><span class="p">)</span> <span class="n">outRes</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">etaID</span><span class="p">,</span> <span class="n">nEta</span><span class="p">,</span> <span class="n">minEta</span><span class="p">,</span> <span class="n">maxEta</span><span class="p">);</span></div><div class='line' id='LC216'><br/></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Z: centrality</span></div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">1</span><span class="p">)</span> <span class="n">outRes</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">centID</span><span class="p">,</span> <span class="n">nCent</span><span class="p">,</span> <span class="n">minCent</span><span class="p">,</span> <span class="n">maxCent</span><span class="p">);</span></div><div class='line' id='LC219'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC220'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC221'><br/></div><div class='line' id='LC222'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">//</span></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// -- Create output for MC generated ------------------------------------------------------------</span></div><div class='line' id='LC224'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">//</span></div><div class='line' id='LC225'><br/></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isMC</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// create ouput</span></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">outMC</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="n">Form</span><span class="p">(</span><span class="s">&quot;phi_MCGen%s&quot;</span><span class="p">,</span> <span class="n">suffix</span><span class="p">.</span><span class="n">Data</span><span class="p">()),</span> <span class="n">outputType</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span> <span class="s">&quot;MOTHER&quot;</span><span class="p">);</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// selection settings</span></div><div class='line' id='LC230'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC231'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC232'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">SetMotherPDG</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">SetMotherMass</span><span class="p">(</span><span class="n">part</span><span class="o">-&gt;</span><span class="n">Mass</span><span class="p">());</span></div><div class='line' id='LC234'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// pair cuts</span></div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">cutsPair</span><span class="p">)</span> <span class="n">outMC</span><span class="o">-&gt;</span><span class="n">SetPairCuts</span><span class="p">(</span><span class="n">cutsPair</span><span class="p">);</span></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis X: invmass</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">imID</span><span class="p">,</span> <span class="n">nIM</span><span class="p">,</span> <span class="n">minIM</span><span class="p">,</span> <span class="n">maxIM</span><span class="p">);</span></div><div class='line' id='LC238'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isFullOutput</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Y: transverse momentum</span></div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outMC</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">ptID</span><span class="p">,</span> <span class="n">nPt</span><span class="p">,</span> <span class="n">minPt</span><span class="p">,</span> <span class="n">maxPt</span><span class="p">);</span></div><div class='line' id='LC241'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">useRapidity</span><span class="p">)</span> <span class="n">outMC</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">yID</span><span class="p">,</span> <span class="n">nY</span><span class="p">,</span> <span class="n">minY</span><span class="p">,</span> <span class="n">maxY</span><span class="p">);</span></div><div class='line' id='LC242'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">0</span><span class="p">)</span> <span class="n">outMC</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">etaID</span><span class="p">,</span> <span class="n">nEta</span><span class="p">,</span> <span class="n">minEta</span><span class="p">,</span> <span class="n">maxEta</span><span class="p">);</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis Z: centrality</span></div><div class='line' id='LC244'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">collisionType</span><span class="o">==</span><span class="mi">1</span><span class="p">)</span> <span class="n">outMC</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">centID</span><span class="p">,</span> <span class="n">nCent</span><span class="p">,</span> <span class="n">minCent</span><span class="p">,</span> <span class="n">maxCent</span><span class="p">);</span></div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC246'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC247'><br/></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <span class="cm">/* invariant mass   */</span> <span class="n">Int_t</span> <span class="n">imID</span>   <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kInvMass</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <span class="cm">/* transv. momentum */</span> <span class="n">Int_t</span> <span class="n">ptID</span>   <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kPt</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <span class="cm">/* rapidity         */</span> <span class="n">Int_t</span> <span class="n">yID</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateValue</span><span class="p">(</span><span class="n">AliRsnMiniValue</span><span class="o">::</span><span class="n">kY</span><span class="p">,</span> <span class="n">kFALSE</span><span class="p">);</span></div><div class='line' id='LC252'><br/></div><div class='line' id='LC253'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">const</span> <span class="n">Int_t</span> <span class="n">numAxis</span> <span class="o">=</span> <span class="mi">3</span><span class="p">;</span></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Bool_t</span>  <span class="n">use</span>     <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span> <span class="mi">1</span>      <span class="p">,</span>  <span class="mi">1</span>      <span class="p">,</span>  <span class="mi">1</span>       <span class="p">};</span></div><div class='line' id='LC255'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">name</span>    <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="s">&quot;Unlike&quot;</span><span class="p">,</span> <span class="s">&quot;LikePP&quot;</span><span class="p">,</span> <span class="s">&quot;LikeMM&quot;</span> <span class="p">};</span></div><div class='line' id='LC256'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">TString</span> <span class="n">comp</span>    <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="s">&quot;PAIR&quot;</span>  <span class="p">,</span> <span class="s">&quot;PAIR&quot;</span>  <span class="p">,</span> <span class="s">&quot;PAIR&quot;</span>   <span class="p">};</span></div><div class='line' id='LC257'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Char_t</span>  <span class="n">charge1</span> <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>      <span class="p">};</span></div><div class='line' id='LC258'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Char_t</span>  <span class="n">charge2</span> <span class="p">[</span><span class="n">numAxis</span><span class="p">]</span> <span class="o">=</span> <span class="p">{</span><span class="sc">&#39;-&#39;</span>     <span class="p">,</span> <span class="sc">&#39;+&#39;</span>     <span class="p">,</span> <span class="sc">&#39;-&#39;</span>      <span class="p">};</span></div><div class='line' id='LC259'><br/></div><div class='line' id='LC260'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nIM</span>   <span class="o">=</span> <span class="mi">215</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minIM</span>   <span class="o">=</span> <span class="mf">0.985</span><span class="p">,</span> <span class="n">maxIM</span> <span class="o">=</span>  <span class="mf">1.2</span><span class="p">;</span></div><div class='line' id='LC261'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nPt</span>   <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minPt</span>   <span class="o">=</span> <span class="mf">0.5</span><span class="p">,</span> <span class="n">maxPt</span> <span class="o">=</span> <span class="mf">1.5</span><span class="p">;</span></div><div class='line' id='LC262'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Int_t</span> <span class="n">nY</span>   <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">Double_t</span> <span class="n">minY</span>   <span class="o">=</span> <span class="o">-</span><span class="mf">0.5</span><span class="p">,</span> <span class="n">maxY</span> <span class="o">=</span>  <span class="mf">0.5</span><span class="p">;</span></div><div class='line' id='LC263'><br/></div><div class='line' id='LC264'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">Int_t</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">numAxis</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC265'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">use</span><span class="p">[</span><span class="n">i</span><span class="p">])</span> <span class="k">continue</span><span class="p">;</span></div><div class='line' id='LC266'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// create output</span></div><div class='line' id='LC267'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AliRsnMiniOutput</span> <span class="o">*</span><span class="n">out</span> <span class="o">=</span> <span class="n">taskRsnMini</span><span class="o">-&gt;</span><span class="n">CreateOutput</span><span class="p">(</span><span class="n">Form</span><span class="p">(</span><span class="s">&quot;%s_%s&quot;</span><span class="p">,</span> <span class="n">suffix</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span><span class="n">name</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">Data</span><span class="p">()</span> <span class="p">),</span> <span class="n">outputType</span><span class="p">.</span><span class="n">Data</span><span class="p">(),</span> <span class="n">comp</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">Data</span><span class="p">());</span></div><div class='line' id='LC268'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// selection settings</span></div><div class='line' id='LC269'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC270'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCutID</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">listID1</span><span class="p">);</span></div><div class='line' id='LC271'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC272'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetDaughter</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">AliRsnDaughter</span><span class="o">::</span><span class="n">kKaon</span><span class="p">);</span></div><div class='line' id='LC273'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="n">charge1</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC274'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetCharge</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">charge2</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC275'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetMotherPDG</span><span class="p">(</span><span class="n">pdg</span><span class="p">);</span></div><div class='line' id='LC276'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">SetMotherMass</span><span class="p">(</span><span class="n">part</span><span class="o">-&gt;</span><span class="n">Mass</span><span class="p">());</span></div><div class='line' id='LC277'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// axis X: invmass</span></div><div class='line' id='LC278'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">imID</span><span class="p">,</span> <span class="n">nIM</span><span class="p">,</span> <span class="n">minIM</span><span class="p">,</span> <span class="n">maxIM</span><span class="p">);</span></div><div class='line' id='LC279'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">ptID</span><span class="p">,</span> <span class="n">nPt</span><span class="p">,</span> <span class="n">minPt</span><span class="p">,</span> <span class="n">maxPt</span><span class="p">);</span></div><div class='line' id='LC280'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">out</span><span class="o">-&gt;</span><span class="n">AddAxis</span><span class="p">(</span><span class="n">yID</span><span class="p">,</span> <span class="n">nY</span><span class="p">,</span> <span class="n">minY</span><span class="p">,</span> <span class="n">maxY</span><span class="p">);</span></div><div class='line' id='LC281'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC282'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC283'><span class="p">}</span></div><div class='line' id='LC284'><br/></div><div class='line' id='LC285'><br/></div><div class='line' id='LC286'><br/></div></pre></div>
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
      <li>&copy; 2013 <span title="0.05023s from github-fe123-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
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

