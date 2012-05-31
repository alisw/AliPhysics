# RPM specfile for alimdc static libs

# Package contains both ROOT and AliRoot
# static libs needed by mStreamRecorder
# in order to ROOT-ify the incoming raw
# data
# Example:
# in $ALICE_ROOT: mkdir RPMS; mkdir RPMS/i386;
# then: rpmbuild --verbose --define "_topdir $ALICE_ROOT" --define "_tmppath $ALICE_ROOT" -bb RAW/alimdc.spec
# then: cp $ALICE_ROOT/RPMS/i386/alimdc-*.rpm <destination>

Summary: AliMDC static libraries
Name: alimdc
Version: 4.09.Release
Release: 1.SLC4
# Copyright: CERN Alice Off-line
License: CERN Alice Off-line
Vendor: ALICE Core Off-line Group
URL: http://aliceinfo.cern.ch
Group: Applications/Alice
Prefix: /opt/%{name}
BuildRoot: %{_tmppath}/%{name}-root

# automatic dependencies
AutoReqProv: yes

# list here required RPM packages for runtime
Requires: glibc

Provides: alimdc

# description of the package
%description
Package contains both ROOT and AliRoot
static libs needed by mStreamRecorder
in order to ROOT-ify the incoming raw
data. The package version correspond to
the AliRoot one.

# list of files to be installed
%files
%defattr (-,root,root)
%{prefix}/lib/libAliMDC.a
%{prefix}/lib/libRoot.a
%{prefix}/lib/libpcre.a
%{prefix}/lib/libfreetype.a
%{prefix}/include/mdc.h
