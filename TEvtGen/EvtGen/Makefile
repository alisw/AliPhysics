#New EvtGen makefile October 5, 1998
#Authors David Lange and Anders Ryd


include $(ROOTSYS)/test/Makefile.arch

SRCSC= $(wildcard *.cc)
SRCSF= $(wildcard *.F)

OBJSC = $(SRCSC:%.cc=%.o)
OBJSF = $(SRCSF:%.F=%.o)

OBJS = $(OBJSC)
OBJS += $(OBJSF)

LIBOBJS=$(filter-out top.o testEvtGen.o ,$(OBJS))

include config.mk

lib: depend comp ar baselib modelslib shar 
#lib: comp ar baselib modelslib

comp:$(OBJS)

clean:
	rm -f *.o libEvtGen.a libEvtGenBase.a libEvtGenModels.a ../Evt*/*.o testEvtGen libEvtGen.so ../EvtGenBase/libEvtGenBase.so ../EvtGenModels/libEvtGenModels.so


depend:
	@echo "Checking makedepend in EvtGen"
	@makedepend -- $(EVTDEPINC) -- -Y $(SRCSC) $(SRCSF)


#.cc.o:
#	$(EVTCXX) $(EVTCXXFLAGS) -c $<

.cc.o:
	$(CXX) $(CXXFLAGS) $(EVTCXXFLAGS) -c $<


.F.o:
	$(EVTF77) $(EVTF77FLAGS) -c $<

ar:
	@if (test -f libEvtGen.a ) then (rm libEvtGen.a) fi
	@ar rv libEvtGen.a $(LIBOBJS) >& /dev/null
	@echo "Done making libEvtGen.a" 

shar:
	@$(CXX) $(CXXFLAGS) $(SOFLAGS) -o libEvtGen.so $(LIBOBJS)
	@echo "Done making libEvtGen.so"
baselib:
	@cd ../EvtGenBase; gmake lib; cd ../EvtGen

modelslib:
	@cd ../EvtGenModels; gmake lib; cd ../EvtGen

# DO NOT DELETE THIS LINE -- make depend depends on it.

EvtCheckDecays.o: ../EvtGenBase/EvtPatches.hh ../EvtGen/EvtCheckDecays.hh
EvtCheckDecays.o: ../EvtGenBase/EvtDecayTable.hh
EvtCheckDecays.o: ../EvtGenBase/EvtDecayBase.hh ../EvtGenBase/EvtId.hh
EvtCheckDecays.o: /usr/include/c++/3.4.6/iostream
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++config.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/os_defines.h
EvtCheckDecays.o: /usr/include/features.h /usr/include/sys/cdefs.h
EvtCheckDecays.o: /usr/include/gnu/stubs.h /usr/include/c++/3.4.6/ostream
EvtCheckDecays.o: /usr/include/c++/3.4.6/ios /usr/include/c++/3.4.6/iosfwd
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++locale.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/cstring
EvtCheckDecays.o: /usr/include/c++/3.4.6/cstddef
EvtCheckDecays.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
EvtCheckDecays.o: /usr/include/string.h /usr/include/c++/3.4.6/cstdio
EvtCheckDecays.o: /usr/include/stdio.h /usr/include/bits/types.h
EvtCheckDecays.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
EvtCheckDecays.o: /usr/include/libio.h /usr/include/_G_config.h
EvtCheckDecays.o: /usr/include/wchar.h /usr/include/bits/wchar.h
EvtCheckDecays.o: /usr/include/gconv.h
EvtCheckDecays.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
EvtCheckDecays.o: /usr/include/bits/stdio_lim.h
EvtCheckDecays.o: /usr/include/bits/sys_errlist.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/clocale /usr/include/locale.h
EvtCheckDecays.o: /usr/include/bits/locale.h /usr/include/langinfo.h
EvtCheckDecays.o: /usr/include/nl_types.h /usr/include/iconv.h
EvtCheckDecays.o: /usr/include/libintl.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++io.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/gthr.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/gthr-default.h
EvtCheckDecays.o: /usr/include/pthread.h /usr/include/sched.h
EvtCheckDecays.o: /usr/include/time.h /usr/include/bits/sched.h
EvtCheckDecays.o: /usr/include/signal.h /usr/include/bits/sigset.h
EvtCheckDecays.o: /usr/include/bits/pthreadtypes.h
EvtCheckDecays.o: /usr/include/bits/initspin.h /usr/include/bits/sigthread.h
EvtCheckDecays.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
EvtCheckDecays.o: /usr/include/bits/confname.h /usr/include/c++/3.4.6/cctype
EvtCheckDecays.o: /usr/include/ctype.h /usr/include/endian.h
EvtCheckDecays.o: /usr/include/bits/endian.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stringfwd.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/postypes.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/cwchar /usr/include/c++/3.4.6/ctime
EvtCheckDecays.o: /usr/include/stdint.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/functexcept.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/exception_defines.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/exception
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/char_traits.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_algobase.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/climits
EvtCheckDecays.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/limits.h
EvtCheckDecays.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/syslimits.h
EvtCheckDecays.o: /usr/include/limits.h /usr/include/c++/3.4.6/cstdlib
EvtCheckDecays.o: /usr/include/stdlib.h /usr/include/c++/3.4.6/new
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_pair.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/type_traits.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_iterator_base_types.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_iterator_base_funcs.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/concept_check.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_iterator.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/debug/debug.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/cassert /usr/include/assert.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/localefwd.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/ios_base.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/atomicity.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/atomic_word.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/locale_classes.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/string /usr/include/c++/3.4.6/memory
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/allocator.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++allocator.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/ext/new_allocator.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_construct.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_uninitialized.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_raw_storage_iter.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/limits
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_function.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/basic_string.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/algorithm
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_algo.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_heap.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_tempbuf.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/basic_string.tcc
EvtCheckDecays.o: /usr/include/c++/3.4.6/streambuf
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/streambuf.tcc
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/basic_ios.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/streambuf_iterator.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/locale_facets.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/cwctype /usr/include/wctype.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/ctype_base.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/ctype_inline.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/codecvt.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/codecvt_specializations.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/time_members.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/messages_members.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/basic_ios.tcc
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/ostream.tcc
EvtCheckDecays.o: /usr/include/c++/3.4.6/locale
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/locale_facets.tcc
EvtCheckDecays.o: /usr/include/c++/3.4.6/typeinfo
EvtCheckDecays.o: /usr/include/c++/3.4.6/istream
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/istream.tcc
EvtCheckDecays.o: ../EvtGenBase/EvtSpinType.hh ../EvtGenBase/EvtReport.hh
EvtCheckDecays.o: /usr/include/c++/3.4.6/vector
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_vector.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_bvector.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/vector.tcc
EvtCheckDecays.o: ../EvtGenBase/EvtParticleDecayList.hh
EvtCheckDecays.o: ../EvtGenBase/EvtParticleDecay.hh ../EvtGenBase/EvtPDL.hh
EvtCheckDecays.o: ../EvtGenBase/EvtPartProp.hh
EvtCheckDecays.o: ../EvtGenBase/EvtAbsLineShape.hh
EvtCheckDecays.o: ../EvtGenBase/EvtStringHash.hh /usr/include/c++/3.4.6/map
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_tree.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/cpp_type_traits.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_map.h
EvtCheckDecays.o: /usr/include/c++/3.4.6/bits/stl_multimap.h
EvtGen.o: ../EvtGenBase/EvtPatches.hh /usr/include/stdio.h
EvtGen.o: /usr/include/features.h /usr/include/sys/cdefs.h
EvtGen.o: /usr/include/gnu/stubs.h
EvtGen.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
EvtGen.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
EvtGen.o: /usr/include/bits/typesizes.h /usr/include/libio.h
EvtGen.o: /usr/include/_G_config.h /usr/include/wchar.h
EvtGen.o: /usr/include/bits/wchar.h /usr/include/gconv.h
EvtGen.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
EvtGen.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
EvtGen.o: /usr/include/c++/3.4.6/fstream /usr/include/c++/3.4.6/istream
EvtGen.o: /usr/include/c++/3.4.6/ios /usr/include/c++/3.4.6/iosfwd
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++config.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/os_defines.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++locale.h
EvtGen.o: /usr/include/c++/3.4.6/cstring /usr/include/c++/3.4.6/cstddef
EvtGen.o: /usr/include/string.h /usr/include/c++/3.4.6/cstdio
EvtGen.o: /usr/include/c++/3.4.6/clocale /usr/include/locale.h
EvtGen.o: /usr/include/bits/locale.h /usr/include/langinfo.h
EvtGen.o: /usr/include/nl_types.h /usr/include/iconv.h /usr/include/libintl.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++io.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/gthr.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/gthr-default.h
EvtGen.o: /usr/include/pthread.h /usr/include/sched.h /usr/include/time.h
EvtGen.o: /usr/include/bits/sched.h /usr/include/signal.h
EvtGen.o: /usr/include/bits/sigset.h /usr/include/bits/pthreadtypes.h
EvtGen.o: /usr/include/bits/initspin.h /usr/include/bits/sigthread.h
EvtGen.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
EvtGen.o: /usr/include/bits/confname.h /usr/include/c++/3.4.6/cctype
EvtGen.o: /usr/include/ctype.h /usr/include/endian.h
EvtGen.o: /usr/include/bits/endian.h /usr/include/c++/3.4.6/bits/stringfwd.h
EvtGen.o: /usr/include/c++/3.4.6/bits/postypes.h
EvtGen.o: /usr/include/c++/3.4.6/cwchar /usr/include/c++/3.4.6/ctime
EvtGen.o: /usr/include/stdint.h /usr/include/c++/3.4.6/bits/functexcept.h
EvtGen.o: /usr/include/c++/3.4.6/exception_defines.h
EvtGen.o: /usr/include/c++/3.4.6/exception
EvtGen.o: /usr/include/c++/3.4.6/bits/char_traits.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_algobase.h
EvtGen.o: /usr/include/c++/3.4.6/climits
EvtGen.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/limits.h
EvtGen.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/syslimits.h
EvtGen.o: /usr/include/limits.h /usr/include/c++/3.4.6/cstdlib
EvtGen.o: /usr/include/stdlib.h /usr/include/c++/3.4.6/new
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_pair.h
EvtGen.o: /usr/include/c++/3.4.6/bits/type_traits.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_iterator_base_types.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_iterator_base_funcs.h
EvtGen.o: /usr/include/c++/3.4.6/bits/concept_check.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_iterator.h
EvtGen.o: /usr/include/c++/3.4.6/debug/debug.h /usr/include/c++/3.4.6/cassert
EvtGen.o: /usr/include/assert.h /usr/include/c++/3.4.6/bits/localefwd.h
EvtGen.o: /usr/include/c++/3.4.6/bits/ios_base.h
EvtGen.o: /usr/include/c++/3.4.6/bits/atomicity.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/atomic_word.h
EvtGen.o: /usr/include/c++/3.4.6/bits/locale_classes.h
EvtGen.o: /usr/include/c++/3.4.6/string /usr/include/c++/3.4.6/memory
EvtGen.o: /usr/include/c++/3.4.6/bits/allocator.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/c++allocator.h
EvtGen.o: /usr/include/c++/3.4.6/ext/new_allocator.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_construct.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_uninitialized.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_raw_storage_iter.h
EvtGen.o: /usr/include/c++/3.4.6/limits
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_function.h
EvtGen.o: /usr/include/c++/3.4.6/bits/basic_string.h
EvtGen.o: /usr/include/c++/3.4.6/algorithm
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_algo.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_heap.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_tempbuf.h
EvtGen.o: /usr/include/c++/3.4.6/bits/basic_string.tcc
EvtGen.o: /usr/include/c++/3.4.6/streambuf
EvtGen.o: /usr/include/c++/3.4.6/bits/streambuf.tcc
EvtGen.o: /usr/include/c++/3.4.6/bits/basic_ios.h
EvtGen.o: /usr/include/c++/3.4.6/bits/streambuf_iterator.h
EvtGen.o: /usr/include/c++/3.4.6/bits/locale_facets.h
EvtGen.o: /usr/include/c++/3.4.6/cwctype /usr/include/wctype.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/ctype_base.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/ctype_inline.h
EvtGen.o: /usr/include/c++/3.4.6/bits/codecvt.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/codecvt_specializations.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/time_members.h
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/messages_members.h
EvtGen.o: /usr/include/c++/3.4.6/bits/basic_ios.tcc
EvtGen.o: /usr/include/c++/3.4.6/bits/istream.tcc
EvtGen.o: /usr/include/c++/3.4.6/locale
EvtGen.o: /usr/include/c++/3.4.6/bits/locale_facets.tcc
EvtGen.o: /usr/include/c++/3.4.6/typeinfo /usr/include/c++/3.4.6/ostream
EvtGen.o: /usr/include/c++/3.4.6/bits/ostream.tcc
EvtGen.o: /usr/include/c++/3.4.6/i386-redhat-linux/bits/basic_file.h
EvtGen.o: /usr/include/c++/3.4.6/bits/fstream.tcc /usr/include/math.h
EvtGen.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
EvtGen.o: /usr/include/bits/mathcalls.h ../EvtGenBase/EvtComplex.hh
EvtGen.o: /usr/include/c++/3.4.6/iostream ../EvtGenBase/EvtConst.hh
EvtGen.o: ../EvtGen/EvtGen.hh ../EvtGenBase/EvtPDL.hh
EvtGen.o: ../EvtGenBase/EvtPartProp.hh ../EvtGenBase/EvtSpinType.hh
EvtGen.o: ../EvtGenBase/EvtReport.hh ../EvtGenBase/EvtId.hh
EvtGen.o: ../EvtGenBase/EvtAbsLineShape.hh /usr/include/c++/3.4.6/vector
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_vector.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_bvector.h
EvtGen.o: /usr/include/c++/3.4.6/bits/vector.tcc
EvtGen.o: ../EvtGenBase/EvtStringHash.hh /usr/include/c++/3.4.6/map
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_tree.h
EvtGen.o: /usr/include/c++/3.4.6/bits/cpp_type_traits.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_map.h
EvtGen.o: /usr/include/c++/3.4.6/bits/stl_multimap.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TLorentzVector.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMath.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/Rtypes.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/RConfig.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/RVersion.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/DllImport.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/Rtypeinfo.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TGenericClassInfo.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TSchemaHelper.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMathBase.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TError.h
EvtGen.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/float.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TVector3.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TVector2.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TObject.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TStorage.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TVersionCheck.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/Riosfwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TBuffer.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrix.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixF.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixT.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixTBase.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TNamed.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TString.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TRefCnt.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixFBasefwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixDBasefwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TVectorFfwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TVectorDfwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixTUtils.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixFfwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixFUtils.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TMatrixFUtilsfwd.h
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TRotation.h
EvtGen.o: /usr/include/c++/3.4.6/list /usr/include/c++/3.4.6/bits/stl_list.h
EvtGen.o: /usr/include/c++/3.4.6/bits/list.tcc ../EvtGenBase/EvtVector4R.hh
EvtGen.o: ../EvtGenBase/EvtVectorParticle.hh ../EvtGenBase/EvtVector4C.hh
EvtGen.o: ../EvtGenBase/EvtVector3C.hh ../EvtGenBase/EvtVector3R.hh
EvtGen.o: ../EvtGenBase/EvtParticle.hh ../EvtGenBase/EvtSpinDensity.hh
EvtGen.o: ../EvtGenBase/EvtScalarParticle.hh ../EvtGenBase/EvtDecayTable.hh
EvtGen.o: ../EvtGenBase/EvtDecayBase.hh ../EvtGenBase/EvtParticleDecayList.hh
EvtGen.o: ../EvtGenBase/EvtParticleDecay.hh ../EvtGenBase/EvtStdHep.hh
EvtGen.o: ../EvtGenBase/EvtSecondary.hh ../EvtGenBase/EvtRandom.hh
EvtGen.o: ../EvtGenBase/EvtRandomEngine.hh
EvtGen.o: /home/fionda/cern/root/v5-23-04/include/TRandom.h
EvtGen.o: ../EvtGenBase/EvtSimpleRandomEngine.hh
EvtGen.o: ../EvtGenBase/EvtParticleFactory.hh ../EvtGenModels/EvtModelReg.hh
EvtGen.o: ../EvtGenBase/EvtStatus.hh ../EvtGenBase/EvtAbsRadCorr.hh
EvtGen.o: ../EvtGenBase/EvtRadCorr.hh ../EvtGenModels/EvtPHOTOS.hh
