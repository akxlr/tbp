# This file is part of libDAI - http://www.libdai.org/
#
# Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.


# Load the platform independent build configuration file
include Makefile.ALL

# Load the local configuration from Makefile.conf
include Makefile.conf

# Set version and date
DAI_VERSION="git HEAD"
DAI_DATE="July 17, 2015 - or later"

# Directories of libDAI sources
# Location of libDAI headers
INC=include/dai
# Location of libDAI source files
SRC=src
# Destination directory of libDAI library
LIB=lib

# Set final compiler flags
ifdef DEBUG
  CCFLAGS:=$(CCFLAGS) $(CCDEBUGFLAGS)
else
  CCFLAGS:=$(CCFLAGS) $(CCNODEBUGFLAGS)
endif

# Define build targets
TARGETS:=lib utils
ifdef WITH_TESTS
  TARGETS:=$(TARGETS) tests
endif
ifdef WITH_EXAMPLES
  TARGETS:=$(TARGETS) examples
endif
ifdef WITH_MATLAB
  TARGETS:=$(TARGETS) matlabs
endif
ifdef WITH_DOC
  TARGETS:=$(TARGETS) doc
endif

# Define conditional build targets
NAMES:=graph dag bipgraph varset daialg alldai clustergraph factor factorgraph properties regiongraph cobwebgraph util weightedgraph exceptions exactinf evidence emalg io
NAMES:=$(NAMES) bp fbp trwbp mf hak lc treeep jtree mr gibbs bbp cbp bp_dual decmap glc dfactor tbp

# Define standard libDAI header dependencies, source file names and object file names
HEADERS=$(foreach name,graph dag bipgraph index var factor varset smallset prob daialg properties alldai enum exceptions util,$(INC)/$(name).h)
SOURCES:=$(foreach name,$(NAMES),$(SRC)/$(name).cpp)
OBJECTS:=$(foreach name,$(NAMES),$(name)$(OE))

# Setup final command for C++ compiler
ifneq ($(OS),WINDOWS)
  CC:=$(CC) $(CCINC) $(CCFLAGS) $(WITHFLAGS) $(CCLIB)
else
  CC:=$(CC) $(CCINC) $(CCFLAGS) $(WITHFLAGS)
  LIBS:=$(LIBS) $(CCLIB)
endif

# Setup final command for MEX
ifdef NEW_MATLAB
  MEXFLAGS:=$(MEXFLAGS) -largeArrayDims
else
  MEXFLAGS:=$(MEXFLAGS) -DSMALLMEM
endif
MEX:=$(MEX) $(MEXINC) $(MEXFLAGS) $(WITHFLAGS) $(MEXLIBS) $(MEXLIB)


# META TARGETS
###############

all : $(TARGETS)
	@echo
	@echo libDAI built successfully!

EXAMPLES=$(foreach name,example example_bipgraph example_varset example_permute example_sprinkler example_sprinkler_em,examples/$(name)$(EE))
EXAMPLES:=$(EXAMPLES) examples/example_sprinkler_gibbs$(EE)
ifdef WITH_CIMG
  EXAMPLES:=$(EXAMPLES) examples/example_imagesegmentation$(EE)
endif
examples : $(EXAMPLES)

matlabs : matlab/dai$(ME) matlab/dai_readfg$(ME) matlab/dai_writefg$(ME) matlab/dai_potstrength$(ME) matlab/dai_jtree$(ME)

unittests : tests/unit/var_test$(EE) tests/unit/smallset_test$(EE) tests/unit/varset_test$(EE) tests/unit/graph_test$(EE) tests/unit/dag_test$(EE) tests/unit/bipgraph_test$(EE) tests/unit/weightedgraph_test$(EE) tests/unit/enum_test$(EE) tests/unit/enum_test$(EE) tests/unit/util_test$(EE) tests/unit/exceptions_test$(EE) tests/unit/properties_test$(EE) tests/unit/index_test$(EE) tests/unit/prob_test$(EE) tests/unit/factor_test$(EE) tests/unit/factorgraph_test$(EE) tests/unit/clustergraph_test$(EE) tests/unit/regiongraph_test$(EE) tests/unit/daialg_test$(EE) tests/unit/alldai_test$(EE)
	@echo 'Running unit tests...'
	@echo
	tests/unit/var_test$(EE)
	tests/unit/smallset_test$(EE)
	tests/unit/varset_test$(EE)
	tests/unit/graph_test$(EE)
	tests/unit/dag_test$(EE)
	tests/unit/bipgraph_test$(EE)
	tests/unit/weightedgraph_test$(EE)
	tests/unit/enum_test$(EE)
	tests/unit/util_test$(EE)
	tests/unit/exceptions_test$(EE)
	tests/unit/properties_test$(EE)
	tests/unit/index_test$(EE)
	tests/unit/prob_test$(EE)
	tests/unit/factor_test$(EE)
	tests/unit/factorgraph_test$(EE)
	tests/unit/clustergraph_test$(EE)
	tests/unit/regiongraph_test$(EE)
	tests/unit/daialg_test$(EE)
	tests/unit/alldai_test$(EE)
	@echo
	@echo 'All unit tests completed successfully!'
	@echo

tests : tests/testdai$(EE) tests/testem/testem$(EE) tests/testbbp$(EE) $(unittests)

utils : utils/createfg$(EE) utils/fg2dot$(EE) utils/fginfo$(EE) utils/uai2fg$(EE) utils/dfgmarg$(EE)

lib: $(LIB)/libdai$(LE)


# OBJECTS
##########

%$(OE) : $(SRC)/%.cpp $(INC)/%.h $(HEADERS)
	$(CC) -c $<

bbp$(OE) : $(SRC)/bbp.cpp $(INC)/bbp.h $(INC)/bp_dual.h $(HEADERS)
	$(CC) -c $<

cbp$(OE) : $(SRC)/cbp.cpp $(INC)/cbp.h $(INC)/bbp.h $(INC)/bp_dual.h $(HEADERS)
	$(CC) -c $<

hak$(OE) : $(SRC)/hak.cpp $(INC)/hak.h $(HEADERS) $(INC)/regiongraph.h
	$(CC) -c $<

jtree$(OE) : $(SRC)/jtree.cpp $(INC)/jtree.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h
	$(CC) -c $<

treeep$(OE) : $(SRC)/treeep.cpp $(INC)/treeep.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h $(INC)/jtree.h
	$(CC) -c $<

emalg$(OE) : $(SRC)/emalg.cpp $(INC)/emalg.h $(INC)/evidence.h $(HEADERS)
	$(CC) -c $<

decmap$(OE) : $(SRC)/decmap.cpp $(INC)/decmap.h $(HEADERS)
	$(CC) -c $<

glc$(OE) : $(SRC)/glc.cpp $(INC)/glc.h $(HEADERS) $(INC)/cobwebgraph.h
	$(CC) -c $<


# EXAMPLES
###########

examples/%$(EE) : examples/%.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_sprinkler_gibbs$(EE) : examples/example_sprinkler_gibbs.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_imagesegmentation$(EE) : examples/example_imagesegmentation.cpp $(HEADERS) $(LIB)/libdai$(LE)
ifdef NEW_CIMG
	$(CC) -DNEW_CIMG $(CIMGINC) $(CCO)$@ $< $(LIBS) $(CIMGLIBS)
else
	$(CC) $(CIMGINC) $(CCO)$@ $< $(LIBS) $(CIMGLIBS)
endif


# UNIT TESTS
#############

tests/unit/%$(EE) : tests/unit/%.cpp $(HEADERS) $(LIB)/libdai$(LE)
ifneq ($(OS),WINDOWS)
	$(CC) -DBOOST_TEST_DYN_LINK $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)
else
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF) /SUBSYSTEM:CONSOLE
endif


# TESTS
########

tests/testdai$(EE) : tests/testdai.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)
tests/testem/testem$(EE) : tests/testem/testem.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)
tests/testbbp$(EE) : tests/testbbp.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)


# MATLAB INTERFACE
###################

matlab/dai$(ME) : $(SRC)/matlab/dai.cpp $(HEADERS) $(SOURCES) $(SRC)/matlab/matlab.cpp
	$(MEX) -output $@ $< $(SRC)/matlab/matlab.cpp $(SOURCES)

matlab/dai_readfg$(ME) : $(SRC)/matlab/dai_readfg.cpp $(HEADERS) $(SRC)/matlab/matlab.cpp $(SRC)/factorgraph.cpp $(SRC)/exceptions.cpp $(SRC)/bipgraph.cpp $(SRC)/graph.cpp $(SRC)/factor.cpp $(SRC)/util.cpp
	$(MEX) -output $@ $< $(SRC)/matlab/matlab.cpp $(SRC)/factorgraph.cpp $(SRC)/exceptions.cpp $(SRC)/bipgraph.cpp $(SRC)/graph.cpp $(SRC)/factor.cpp $(SRC)/util.cpp

matlab/dai_writefg$(ME) : $(SRC)/matlab/dai_writefg.cpp $(HEADERS) $(SRC)/matlab/matlab.cpp $(SRC)/factorgraph.cpp $(SRC)/exceptions.cpp $(SRC)/bipgraph.cpp $(SRC)/graph.cpp $(SRC)/factor.cpp $(SRC)/util.cpp
	$(MEX) -output $@ $< $(SRC)/matlab/matlab.cpp $(SRC)/factorgraph.cpp $(SRC)/exceptions.cpp $(SRC)/bipgraph.cpp $(SRC)/graph.cpp $(SRC)/factor.cpp $(SRC)/util.cpp

matlab/dai_potstrength$(ME) : $(SRC)/matlab/dai_potstrength.cpp $(HEADERS) $(SRC)/matlab/matlab.cpp $(SRC)/exceptions.cpp
	$(MEX) -output $@ $< $(SRC)/matlab/matlab.cpp $(SRC)/exceptions.cpp

matlab/dai_jtree$(ME) : $(SRC)/matlab/dai_jtree.cpp $(HEADERS) $(SOURCES) $(SRC)/matlab/matlab.cpp
	$(MEX) -output $@ $< $(SRC)/matlab/matlab.cpp $(SOURCES)


# UTILS
########

utils/createfg$(EE) : utils/createfg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)

utils/fg2dot$(EE) : utils/fg2dot.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

utils/fginfo$(EE) : utils/fginfo.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

utils/uai2fg$(EE) : utils/uai2fg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

utils/dfgmarg$(EE) : utils/dfgmarg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)


# LIBRARY
##########

ifneq ($(OS),WINDOWS)
$(LIB)/libdai$(LE) : $(OBJECTS)
	-mkdir -p lib
	ar rcus $(LIB)/libdai$(LE) $(OBJECTS)
else
$(LIB)/libdai$(LE) : $(OBJECTS)
	-mkdir lib
	lib /out:$(LIB)/libdai$(LE) $(OBJECTS)
endif


# REGRESSION TESTS
###################

testregression : tests/testdai$(EE)
	@echo Starting regression test...this can take a minute or so!
ifneq ($(OS),WINDOWS)
	cd tests && ./testregression && cd ..
else
	cd tests && testregression.bat && cd ..
endif

testem : tests/testem/testem$(EE)
	@echo Starting EM tests
ifneq ($(OS),WINDOWS)
	cd tests/testem && ./runtests && cd ../..
else
	cd tests\testem && runtests && cd ..\..
endif


# DOCUMENTATION
################

doc : $(INC)/*.h $(SRC)/*.cpp examples/*.cpp doxygen.conf
	doxygen doxygen.conf

README : doc scripts/makeREADME Makefile
	DAI_VERSION=$(DAI_VERSION) DAI_DATE=$(DAI_DATE) scripts/makeREADME

TAGS :
	etags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp
	ctags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp


# CLEAN
########

.PHONY : clean
ifneq ($(OS),WINDOWS)
clean :
	-rm $(OBJECTS)
	-rm matlab/*$(ME)
	-rm examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE) examples/example_permute$(EE) examples/example_sprinkler$(EE) examples/example_sprinkler_gibbs$(EE) examples/example_sprinkler_em$(EE) examples/example_imagesegmentation$(EE)
	-rm tests/testdai$(EE) tests/testem/testem$(EE) tests/testbbp$(EE)
	-rm tests/unit/var_test$(EE) tests/unit/smallset_test$(EE) tests/unit/varset_test$(EE) tests/unit/graph_test$(EE) tests/unit/dag_test$(EE) tests/unit/bipgraph_test$(EE) tests/unit/weightedgraph_test$(EE) tests/unit/enum_test$(EE) tests/unit/util_test$(EE) tests/unit/exceptions_test$(EE) tests/unit/properties_test$(EE) tests/unit/index_test$(EE) tests/unit/prob_test$(EE) tests/unit/factor_test$(EE) tests/unit/factorgraph_test$(EE) tests/unit/clustergraph_test$(EE) tests/unit/regiongraph_test$(EE) tests/unit/daialg_test$(EE) tests/unit/alldai_test$(EE)
	-rm factorgraph_test.fg alldai_test.aliases
	-rm utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE) utils/uai2fg$(EE) utils/dfgmarg$(EE)
	-rm -R doc
	-rm -R lib
else
clean :
	-del *.obj
	-del *.ilk
	-del *.pdb
	-del matlab\*$(ME)
	-del examples\*$(EE)
	-del examples\*$(EE).manifest
	-del examples\*.ilk
	-del examples\*.pdb
	-del tests\*$(EE)
	-del tests\*$(EE).manifest
	-del tests\*.pdb
	-del tests\*.ilk
	-del tests\testem\*$(EE)
	-del tests\testem\*$(EE).manifest
	-del tests\testem\*.pdb
	-del tests\testem\*.ilk
	-del utils\*$(EE)
	-del utils\*$(EE).manifest
	-del utils\*.pdb
	-del utils\*.ilk
	-del tests\unit\*_test$(EE)
	-del tests\unit\*_test$(EE).manifest
	-del tests\unit\*_test.pdb
	-del tests\unit\*_test.ilk
	-del factorgraph_test.fg
	-del alldai_test.aliases
	-del $(LIB)\libdai$(LE)
	-rmdir lib
endif
