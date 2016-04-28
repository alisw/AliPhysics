#!groovy

node {
  stage "Build AliRoot"

  def test_script = '''
      # Make sure we have only one builder per directory
      x=`date +"%s"`
      WORKAREA=/build/workarea/pr/${week}

      CURRENT_SLAVE=unknown
      while [[ "$CURRENT_SLAVE" != '' ]]; do
        WORKAREA_INDEX=$((WORKAREA_INDEX+1))
        CURRENT_SLAVE=$(cat $WORKAREA/$WORKAREA_INDEX/current_slave 2> /dev/null || true)
        [[ "$CURRENT_SLAVE" == "$NODE_NAME" ]] && CURRENT_SLAVE=
      done

      mkdir -p $WORKAREA/$WORKAREA_INDEX
      echo $NODE_NAME > $WORKAREA/$WORKAREA_INDEX/current_slave

      rm -fr alibuild alidist
      git clone https://github.com/alisw/alidist
      git clone https://github.com/alisw/alibuild
      (cd alidist && git show)

      alibuild/aliBuild -w $WORKAREA/$WORKAREA_INDEX                       \
                        --reference-sources /build/mirror                  \
                        --debug                                            \
                        --remote-store rsync://repo.marathon.mesos/store/  \
                         build AliPhysics || BUILDERR=$?

      rm -f $WORKAREA/$WORKAREA_INDEX/current_slave
      if [ ! "X$BUILDERR" = X ]; then
        exit $BUILDERR
      fi
    '''

  parallel(
    "slc5": {
      node ("slc5_x86-64-large") {
        dir ("AliPhysics") { checkout scm }
        sh test_script
      }
    }
  )
}
