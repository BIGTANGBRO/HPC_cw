##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=coursework
ConfigurationName      :=Debug
WorkspaceConfiguration := $(ConfigurationName)
WorkspacePath          :=/home/jt2418/Desktop/hpc/CW
ProjectPath            :=/home/jt2418/Desktop/hpc/CW/coursework
IntermediateDirectory  :=../build-$(ConfigurationName)/coursework
OutDir                 :=../build-$(ConfigurationName)/coursework
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Jiaxuan Tang
Date                   :=14/03/22
CodeLitePath           :=/home/jt2418/.codelite
LinkerName             :=g++
SharedObjectLinkerName :=g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=../build-$(ConfigurationName)/bin/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :=$(IntermediateDirectory)/ObjectsList.txt
PCHCompileFlags        :=
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)blas 
ArLibs                 :=  "blas" 
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := g++
CC       := gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=../build-$(ConfigurationName)/coursework/main.cpp$(ObjectSuffix) ../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: MakeIntermediateDirs $(OutputFile)

$(OutputFile): ../build-$(ConfigurationName)/coursework/.d $(Objects) 
	@mkdir -p "../build-$(ConfigurationName)/coursework"
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@mkdir -p "../build-$(ConfigurationName)/coursework"
	@mkdir -p ""../build-$(ConfigurationName)/bin""

../build-$(ConfigurationName)/coursework/.d:
	@mkdir -p "../build-$(ConfigurationName)/coursework"

PreBuild:


##
## Objects
##
../build-$(ConfigurationName)/coursework/main.cpp$(ObjectSuffix): main.cpp ../build-$(ConfigurationName)/coursework/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jt2418/Desktop/hpc/CW/coursework/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
../build-$(ConfigurationName)/coursework/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../build-$(ConfigurationName)/coursework/main.cpp$(ObjectSuffix) -MF../build-$(ConfigurationName)/coursework/main.cpp$(DependSuffix) -MM main.cpp

../build-$(ConfigurationName)/coursework/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../build-$(ConfigurationName)/coursework/main.cpp$(PreprocessSuffix) main.cpp

../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(ObjectSuffix): ReactionDiffusion.cpp ../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jt2418/Desktop/hpc/CW/coursework/ReactionDiffusion.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/ReactionDiffusion.cpp$(ObjectSuffix) $(IncludePath)
../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(DependSuffix): ReactionDiffusion.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(ObjectSuffix) -MF../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(DependSuffix) -MM ReactionDiffusion.cpp

../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(PreprocessSuffix): ReactionDiffusion.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../build-$(ConfigurationName)/coursework/ReactionDiffusion.cpp$(PreprocessSuffix) ReactionDiffusion.cpp


-include ../build-$(ConfigurationName)/coursework//*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r $(IntermediateDirectory)


