#object_files (testEffi.o)
OBJECTS := $(wildcard *.o)

#root_stuff (root libraries and needed root options)
ROOTLIBS := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lFoam -lRooFit -lRooFitCore -lGenVector
ROOTCINT := $(shell which rootcint)

#exe_files
EXECUTABLEREAD     := readSideband
EXECUTABLEGENSB    := generateSideband
EXECUTABLERENSB    := renameSideband
EXECUTABLESYST     := systematicSideband
CLASS		   := RooBernsteinSideband
LIBCLASS           := lib$(CLASS).so
CLASSDICT	   := $(CLASS)Dictionary.cxx
CLASSMAP	   := lib$(CLASS).rootmap
CLASSCB 	   := RooDoubleCBFast
CLASSCBDICT	   := $(CLASSCB)Dictionary.cxx

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(CLASS).cxx  $(CLASSDICT) $(CLASSCB).cxx  $(CLASSCBDICT) 

	
all: $(CLASSDICT)   $(CLASSCBDICT) $(EXECUTABLE) $(EXECUTABLEFIT) $(EXMAKEHIST) $(EXECUTABLEREAD) $(EXECUTABLEGENSB) $(EXECUTABLERENSBD) $(EXECUTABLESYST)

dict: $(CLASSDICT)   $(CLASSCBDICT)




lib:  $(LIBCLASS)

$(CLASSDICT): $(CLASS).h $(CLASS)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -rmf $(CLASSMAP) -c $^


$(CLASSCBDICT): $(CLASSCB).h $(CLASSCB)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^





$(EXECUTABLEREAD): $(EXECUTABLEREAD).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLEGENSB): $(EXECUTABLEGENSB).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLERENSB): $(EXECUTABLERENSB).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLESYST): $(EXECUTABLESYST).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I.



#$(EXMAKEHIST): $(EXMAKEHIST).cc 
#	$(CXX) $(CXXFLAGS)  -o $@  $^ $(ROOTLIBS) $(ROOTFLAGS) -I.

$(LIBCLASS): $(CLASS).cxx $(CLASS)Dictionary.cxx
	$(CXX) $(CXXFLAGS) -fPIC -shared  -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS) -I.


#cleaning options
.PHONY: clean cleanall
clean:
	rm -f $(OBJECTS) && rm -f  $(EXECUTABLEREAD) $(EXECUTABLEGENSB) $(EXECUTABLERENSB) $(EXECUTABLESYST) $(CLASSDICT)
cleanhist:
	rm -f  $(EXMAKEHIST)

