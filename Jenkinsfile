#!groovy

node {
  stage "Build AliRoot"
  def test_script = """
      rm -fr alibuild alidist
      git clone https://github.com/alisw/alidist
      git clone https://github.com/alisw/alibuild
      (cd alidist && git show)
      alibuild/aliBuild --reference-sources /build/mirror --debug --remote-store rsync://repo.marathon.mesos/store/ -d build AliRoot
    """

  parallel(
    "slc5": {
      node ("slc5_x86-64-large") {
        dir ("AliRoot") {
          checkout scm
        }
        sh test_script
      }
    }
  )
}
