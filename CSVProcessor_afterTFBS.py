# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 00:24:18 2021

@author: Desktop
"""
import pandas as pd
from abc import ABC, abstractmethod
import os

class InitializeVars:
    def __init__(self):
        self.FileList = []
        self.ExportLocation = ""
        self.FASTALoc = ""
        self.Organism = ""
        self.StartRegion = 0
        self.EndRegion = 0
        self.GeneStartCoord = 0
        self.InSeqNumber = 0
    
    def populate(self, CSVFiles, ExportSource, FASTAFiles, organism,
                 StartCoord, EndCoord, GeneStartCoord, InputNumber):
        self.FileList = CSVFiles
        self.ExportLocation = ExportSource
        self.FASTALoc = FASTAFiles
        self.Organism = organism
        self.StartRegion = StartCoord
        self.EndRegion = EndCoord
        self.GeneStartCoord = GeneStartCoord
        self.InSeqNumber = InputNumber
        
class Sorter:
    def BubbleSort(self, Input):
        Condition = True
        while(Condition):
            Condition = False
            for Vals in range(len(Input) - 1):
                if Input[Vals] > Input[Vals + 1]:
                    Input[Vals], Input[Vals + 1] = Input[Vals + 1], Input[Vals]
                    Condition = True
        return(Input)

class ExportClass:
    def ExportFxn(self, FileName, ExportLocation, FrameInputs):
        if isinstance(FrameInputs, list):
            for Frames, rng in zip(FrameInputs, range(len(FrameInputs))):
                ExportFileName = "{0}_{1}.csv".format(FileName, rng)
                ExportStmnt = ExportLocation + ExportFileName
                Frames.to_csv(ExportStmnt)
        else:
            ExportFileName = "{0}_{1}.csv".format(FileName, 1)
            ExportStmnt = ExportLocation + ExportFileName
            FrameInputs.to_csv(ExportStmnt)
            
class ShortListPotentialMotifs:
    """
    Use the ranges from LAGAN of highly homologous DNA elements
    returns all motifs within that range
    
    i.e.: if observing homology >= 70% in 2500-3000 of geneX, use
    those coordinates as inputs.
    """
    def ExtractMotif(self, StartRegion, EndRegion, DataFrame):
        ExtractRows = []
        for Startrng, Endrng, Ind in zip(DataFrame["Start"], DataFrame["End"], DataFrame.index.values):
            if ((Startrng >= StartRegion) and (Endrng <= EndRegion)):
                ExtractRows.append(DataFrame.iloc[Ind,])
        MotifsOfInterest = pd.DataFrame(ExtractRows).reset_index(drop = True)
        return(MotifsOfInterest) 
    
class ExtractCoordinates:
    """
    This is a very nich function. It is intended strictly to handle the description
    of FASTA files produced by the FASTA exporter script. 
    """
    def CoordExtractor(self, InputFASTA):
        """
        To adjust coordinates
        Start site in the TFBS file + the start site from the description
        End site is (Start site in TFBS + start site description) + (End site TFBS - Start site TFBS)
        End site is not how you have it below
        """
        Coords = []
        with open(InputFASTA, "rt") as Handle:
            Condition = 0
            for lines in Handle:
                if Condition >= 1:
                    break
                Coords.append(lines.split(","))
                Condition += 1
        SplitCoords = [coors.split(":") for coors in Coords[0]]
        Split2Coords = []
        for Coors in SplitCoords:
            for Coors2 in Coors:
                Split2Coords.append(Coors2.split(" "))
        FinalValues = []
        for ele in Split2Coords:
            for vals in ele:
                try:
                    FinalValues.append(int(vals))
                #Excepting multiple errors, index and value error
                except:
                    pass
        #return excludes the size of the read
        return(FinalValues[1:])
   
class InitInherit(ABC):
    @abstractmethod
    def Vars(self, Init):
        pass
    
class Exports(ABC):
    @abstractmethod
    def Export(self):
        pass
    
class InheritSort(ABC):
    @abstractmethod
    def SortFxn(self):
        pass
    
class MotifExtract(ABC):
    @abstractmethod
    def ExtractFxn(self):
        pass
    
class CoordinatesExtract(ABC):
    @abstractmethod
    def ExtractCoords(self):
        pass

class CoordsFromFasta(InitInherit, CoordinatesExtract):
    def ExtractCoords(self, FastaLoc):
        ExtractCoords = ExtractCoordinates().CoordExtractor(InputFASTA=FastaLoc)
        return(ExtractCoords)
    
    def Vars(self, Init):
        FASTAFiles = []
        Coords = []
        for FASTA in os.listdir(Init.FASTALoc):
            if FASTA.endswith(".fasta"):
                FASTAFiles.append(str(Init.FASTALoc + "\{}".format(FASTA)))
        for files in FASTAFiles:
            Coords.append(self.ExtractCoords(FastaLoc=files))
        return(Coords)

class ConvertCoordinates(InitInherit, CoordinatesExtract, Exports):
    def ExtractCoords(self, FastaLoc):
        ExtractCoords = ExtractCoordinates().CoordExtractor(InputFASTA=FastaLoc)
        return(ExtractCoords)
    
    def Export(self, filename, location, frames):
        ExportFxn2 = ExportClass().ExportFxn(FileName=filename, ExportLocation=location, FrameInputs=frames)
        return(ExportFxn2)
    
    def LoadFiles_IfPaths(self, InputPaths):
        StoreFiles = []
        StoreNames = []
        for Paths in InputPaths:
            for Dirs in os.listdir(Paths):
                Temp = []
                TempNames = []
                if os.path.isdir(Paths + "\{}".format(Dirs)) is True:
                    for SubFiles in os.listdir(Paths + "\{}".format(Dirs)):
                        if SubFiles.endswith(".csv"):
                            Path = Paths + "\{0}\{1}".format(Dirs, SubFiles)
                            Temp.append(pd.read_csv(Path))
                            TempNames.append(str(SubFiles))
                    StoreFiles.append(Temp)
                    StoreNames.append(TempNames)
        ListOfFrames=StoreFiles
        filenames=StoreNames
        return(ListOfFrames, filenames)
                    
    def LoadFiles_IfFiles(self, InputPaths):
        filenames=[]
        ListOfFrames = [pd.read_csv(Frames) for Frames in InputPaths]
        return(ListOfFrames, filenames)        

    def FASTA_IfPath(self, FASTApath):
        return([FASTApath + "\{}".format(Files) for Files in os.listdir(FASTApath) if Files.endswith(".fasta")])
    
    def FASTA_IfFiles(self, FASTAfiles):
        return([files for files in FASTAfiles])
              
    def AdjustFunction(self, InputFrame, AdjustedCoordinate):
        for Start, End, rng in zip(InputFrame["Start"], InputFrame["End"], InputFrame.index.values):
            AdjustStart = Start + AdjustedCoordinate[0]
            AdjustEnd = (End - Start) + AdjustStart
            InputFrame.at[rng, "Start"] = int(AdjustStart)
            InputFrame.at[rng, "End"] = int(AdjustEnd)
        return(InputFrame)
    
    def Vars(self, Init):
        if os.path.isdir(Init.FASTALoc) is True:
            ListOfFASTA = self.FASTA_IfPath(FASTApath=Init.FASTALoc)
        elif os.path.isdir(Init.FASTALoc) is False:
            ListOfFASTA = self.FASTA_IfFiles(FASTApath=Init.FASTALoc)
        Coords = [self.ExtractCoords(FastaLoc=FASTA) for FASTA in ListOfFASTA]        
        global ListOfFrames
        global filenames
        if os.path.isdir(Init.FileList[0]) is True:
            FileLoad = self.LoadFiles_IfPaths(InputPaths=Init.FileList)
            ListOfFrames = FileLoad[0]
            filenames = FileLoad[1]
            global AdjustedFrames
            AdjustedFrames = []
            for Frames in range(len(ListOfFrames)):
                AdjustedFrames += ListOfFrames[Frames]
            for FrameInd, Ind2, Coords in zip(ListOfFrames, range(len(ListOfFrames)), Coords):
                for Frames in FrameInd:
                    self.AdjustFunction(InputFrame=Frames, AdjustedCoordinate=Coords)
        elif os.path.isdir(Init.FileList[0]) is False:
            FileLoad = self.LoadFiles_IfFiles(InputPaths=Init.FileList)
            ListOfFrames = FileLoad[0]
            filenames = FileLoad[1]
            for FrameInd, Coords in zip(ListOfFrames, Coords):
                for Frames in FrameInd:
                    self.AdjustFunction(InputFrame=Frames, AdjustedCoordinate=Coords)
        if ((len(filenames) != 0) and (filenames is not None)):
            for FrameInd, NameInd in zip(ListOfFrames, filenames):
                for frames, names in zip(FrameInd, NameInd):
                    ExportName = "Adjusted_" + str(names)
                    frames.to_csv(Init.ExportLocation + ExportName)
        elif len(filenames) == 0:
            self.Export(filename="AdjustedTFBS", location=Init.ExportLocation, frames=ListOfFrames)
            
class SubsetHomologousRegions(InitInherit):
    def Vars(self, Init):
        """
        1. Use the frames from list of frames global var as both the objects to extract from and 
        as an organizing unit, defining the new frames based on the value of the location of the frame in the list
        2. Adjust the coordinates of the start/end regions from the LAGAN alignment (or others)
        3. Extract the rows from the combined frame that are between the given alignment coordinates of interest
        """
        AdjustedAlignmentCoords = []
        for Start, End in zip(Init.StartRegion, Init.EndRegion):
            NewStart = Start + Init.GeneStartCoord
            NewEnd = End + Init.GeneStartCoord
            AdjustedAlignmentCoords.append((NewStart, NewEnd))
        if "AdjustedFrames" in globals() or "ListOfFrames" in globals():
            if "AdjustedFrames" in globals():
                ListOfFrames = AdjustedFrames
            SubsetFrames = [[] for _ in range(len(ListOfFrames))]
            for Coordinates in AdjustedAlignmentCoords:
                for Frames, rng in zip(ListOfFrames, range(len(ListOfFrames))):
                    Temp = []
                    for Start, End, Index in zip(Frames["Start"], Frames["End"], Frames.index.values):
                        if ((Start >= Coordinates[0]) and (End <= Coordinates[1])):
                            Temp.append(Frames.iloc[Index,])
                    SubsetFrames[rng].append(Temp) 
            Indexer = [Sublist for Sublist in range(len(SubsetFrames)) for Element in SubsetFrames[Sublist] if len(Element) != 0]
            SubsetFrames = [Rows for Sublists in SubsetFrames for Rows in Sublists if len(Rows) != 0]
            SubsetFrames = [pd.DataFrame(Frames).drop(["Unnamed: 0"], axis = 1).reset_index(drop=True) for Frames in SubsetFrames]
            AdjustFileNames = [element for sublists in filenames for element in sublists]
            NewFileNames = [AdjustFileNames[Ind] for Ind in Indexer]
            for MotifFrames, Names in zip(SubsetFrames, NewFileNames):
                ExportStmnt = "MotifsOfInterest_{0}.csv".format(Names)
                MotifFrames.to_csv(Init.ExportLocation + ExportStmnt)
        else:
            raise(NameError("Error in creating the ListOfFrames global variable"))
        return(SubsetFrames, Indexer)
        
if __name__ == "__main__":
    Init = InitializeVars()
    Init.populate(CSVFiles=[r"F:\R_dir\Processed\Gad2\Outputs\GroupedFASTA\PcCRE_FastaFiles"],
                   ExportSource=r"F:\R_dir\Processed\Gad2\Outputs\GroupedFASTA\PcCRE_FastaFiles\Seq",
                   FASTAFiles=r"F:\R_dir\Processed\Gad2\Outputs\GroupedFASTA\PcCRE_FastaFiles",
                   organism="Mouse", StartCoord = [], 
                   EndCoord = [], 
                   GeneStartCoord=0, InputNumber=3)

    Cnv = ConvertCoordinates()
    Cnv.Vars(Init)

    if len(Init.StartRegion) != 0 and len(Init.EndRegion) != 0 and Init.GeneStartCoord != 0:
        Extrct = SubsetHomologousRegions()
        Extrct.Vars(Init)

