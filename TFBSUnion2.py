# -*- coding: utf-8 -*-
"""
Created on Sat Aug 28 19:25:51 2021

@author: Desktop
"""
import pandas as pd
from abc import ABC, abstractmethod
import os

class InitializeVars:
    def __init__(self):
        self.geneNames = []
        self.fileList = []
        self.exportLocation = ""
        self.exportName = ""
        
    def populate(self, Names, Files, Export, ExportName):
        self.geneNames = Names
        self.fileList = Files
        self.exportLocation = Export
        self.exportName = ExportName

class Export:
    def ExportFunction(self, InputFrame, ExportLocation, FileName):
        return(InputFrame.to_csv(ExportLocation + "\{}.csv".format(FileName)))

class InheritInitializeVars(ABC):
    @abstractmethod
    def InheritInit(self, Init):
        pass

class InheritExport(ABC):
    @abstractmethod
    def InheritExportFunction(self):
        pass
    
class TFBSMotifs(InheritInitializeVars, InheritExport):
    def InheritExportFunction(self, Frame, Location, Name):
        ExportFxn = Export().ExportFunction(InputFrame=Frame, ExportLocation=Location, FileName=Name)
        return(ExportFxn)
    
    def LoadFunction_IfPaths(self, InputPaths):
        Frames = [[] for _ in range(len(InputPaths))]
        for Paths in range(len(InputPaths)):
            for Files in os.listdir(InputPaths[Paths]):
                if (Files.endswith(".csv") is True):
                    Path = InputPaths[Paths] + "\{}".format(Files)
                    Frames[Paths].append(pd.read_csv(Path))
        return(Frames)

    def LoadFunction_IfFileList(self, InputPaths):
        Frames = [[] for _ in range(len(InputPaths))]
        for Files, rng in zip(InputPaths, range(len(InputPaths))):
            Frames[rng].append(pd.read_csv(Files))
        return(Frames)

    def InheritInit(self, Init):
        """
        1. Create nested lists containing the dataframes of the motifs of interest from each gene 
        2. Combine the individual frames from the motifs of interest from each gene into one concatted
        data frame --> used for determining like TFs since we don't about position of motifs.
        3. Using the first frame create a dictionary of motifs and use that dictionary to count the number of times
        that motifs are repeated in each of the files
        """
        #Should throw this into a seperate loading function
        if (os.path.isdir(Init.fileList[0]) is True):
            Frames = self.LoadFunction_IfPaths(InputPaths=Init.fileList)
        elif (os.path.isfile(Init.fileList[0]) is True):
            Frames = self.LoadFunction_IfFileList(InputPaths=Init.fileList)
        ConcattedFrames = []
        for geneFiles in Frames:
            Concat = pd.concat([IndFrames for IndFrames in geneFiles])
            ConcattedFrames.append(Concat)
        ConcattedFrames = [Frames.drop(["Unnamed: 0"], axis=1) for Frames in ConcattedFrames]
        """
        Dictionary used to count the number of motifs appearing in each concatted
        frame
        """
        MotifDictionary = {
            Motifs: 1 for Motifs in ConcattedFrames[0]["ID"]
            } 
        for Frames in ConcattedFrames[1:]:
            for MotifIDs in Frames["ID"]:
                if MotifIDs in MotifDictionary:
                    MotifDictionary[MotifIDs] += 1
                elif MotifIDs not in MotifDictionary:
                    MotifDictionary[MotifIDs] = 1
        """
        Remove the motifs that appear less than 2 times. I.e.: the motifs that are unique in
        a given gene region.
        """
        MotifsToRemove = [MotifIDs for MotifIDs in MotifDictionary if MotifDictionary[MotifIDs] < 2]
        for MotifIDs in MotifsToRemove:
            del MotifDictionary[MotifIDs]
        """
        Sort the dictionary or sort the values in the dictionary by the most frequently occuring motifs
        """
        MotifID_FromDict = [IDs for IDs in MotifDictionary]
        Frequency_FromDict = [MotifDictionary[IDs] for IDs in MotifDictionary]
        Condition = True
        while(Condition):
            Condition = False
            for rng in range(len(Frequency_FromDict) - 1):
                if (Frequency_FromDict[rng] < Frequency_FromDict[rng + 1]):
                    Frequency_FromDict[rng], Frequency_FromDict[rng + 1] = Frequency_FromDict[rng + 1], Frequency_FromDict[rng]
                    MotifID_FromDict[rng], MotifID_FromDict[rng + 1] = MotifID_FromDict[rng + 1], MotifID_FromDict[rng]
                    Condition = True
        ExtractedRows = [Frames.iloc[Rows,] for Frames in ConcattedFrames 
                         for ID, Rows in zip(Frames["ID"], Frames.index.values)
                         if ID in MotifID_FromDict]
        ExtractedRowsFrame = pd.DataFrame(data=ExtractedRows).reset_index(drop=True)
        ExtractedRowsFrame = ExtractedRowsFrame.sort_values(by=["Score"], ascending=False).reset_index(drop=True)
        MotifCountFrame = pd.DataFrame(data={
            "Motifs":MotifID_FromDict,
            "Frequency":Frequency_FromDict
            })
        ConcattedFrame = pd.concat([ExtractedRowsFrame, MotifCountFrame], axis=1).reset_index(drop=True)
        ConcattedFrame = ConcattedFrame.fillna("")
        """
        Export Concatted Frame
        """   
        return(self.InheritExportFunction(Frame=ConcattedFrame, Location=Init.exportLocation, Name=Init.exportName))
                 
if __name__=='__main__':
    Init = InitializeVars()
    Init.populate(["Gad1", "Gad2", "Slc32a1"], Files=[r"F:\R_dir\Processed\Gad1\Outputs\MotifsOfInterest_Human", 
                                        r"F:\R_dir\Processed\Gad2\Outputs\MotifForHomologusRegions_Human", 
                                        r"F:\R_dir\Processed\Slc32a1\Outputs\MotifsOfInterest_Human"], Export=r"F:\R_dir\Processed\AnalysisResults", 
                  ExportName = "MotifsOfInterest_Human")
    Call = TFBSMotifs()
    Call.InheritInit(Init)
        
        
        
        
        
                    
                
        
        
            
            
            
    
                
