# -*- coding: utf-8 -*-
"""
Created on Sat Aug 28 19:25:51 2021

@author: Desktop
"""
import pandas as pd
from abc import ABC, abstractmethod
import os
import time

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
                    Frame = pd.read_csv(Path).drop(["Unnamed: 0"], axis = 1)
                    Frames[Paths].append(Frame)
        return(Frames)

    def LoadFunction_IfFileList(self, InputPaths):
        Frames = [[] for _ in range(len(InputPaths))]
        for Files, rng in zip(InputPaths, range(len(InputPaths))):
            Frames[rng].append(pd.read_csv(Files))
        return(Frames)
    
    def SortDictionaries_AsLists(self, InputDictionary, Name):
        MotifIDs = [ID for ID in InputDictionary]
        Frequency_FromDict = [InputDictionary[ID] for ID in InputDictionary]
        Set_Condition = True
        while(Set_Condition):
            Set_Condition = False
            for Rng in range(len(Frequency_FromDict) - 1):
                if (Frequency_FromDict[Rng] < Frequency_FromDict[Rng + 1]):
                    Frequency_FromDict[Rng], Frequency_FromDict[Rng + 1] = Frequency_FromDict[Rng + 1], Frequency_FromDict[Rng]
                    MotifIDs[Rng], MotifIDs[Rng + 1] = MotifIDs[Rng + 1], MotifIDs[Rng]
                    Set_Condition = True
        FrameConstructor = {
            "{}_Motifs".format(Name):MotifIDs,
            "{}_Counts".format(Name):Frequency_FromDict
            }
        DictToFrame = pd.DataFrame(FrameConstructor)
        return(DictToFrame)
    
    def MotifCountsPerGene(self, InputDictList, FrameList, Genes):   
        """
        2 parts to this function
        1: Unpack portion which removes motifs from a dictionary and creates n (where n = gene number) dictionaries that are unncounted
        2: Use the initial concatted list to recount.
        """  
        #Get the number of genes inputted, should equal the number of paths/number of file lists
        GeneNumbers = len(Genes)          
        UnpackDicts = [Motifs for Ind in range(len(InputDictList)) for Motifs in InputDictList[Ind]]
        UnpackedDictionary = {
            Motifs: 0 for Motifs in UnpackDicts
            }
        for Motifs in UnpackDicts:
            if Motifs in UnpackedDictionary:
                UnpackedDictionary[Motifs] += 1
        MotifsToRemove = [Motifs for Motifs in UnpackedDictionary if UnpackedDictionary[Motifs] < GeneNumbers]
        for Motifs in MotifsToRemove:
            if Motifs in UnpackedDictionary:
                del UnpackedDictionary[Motifs]
        MotifCountsPerGene = []
        for Frames in FrameList:
            ResetDictionary = {
                Motifs: 0 for Motifs in UnpackedDictionary
                }
            for ID in Frames["ID"]:
                if ID in ResetDictionary:
                    ResetDictionary[ID] += 1
            MotifCountsPerGene.append(ResetDictionary)
        return(MotifCountsPerGene)
    
    def TotalMotifCount(self, InputDictList, FrameList):
        """
        Define an empty list here where dictionaries will be appended to after 
        some processing
        """
        #Create empty Dict
        ReferenceDictionary = {
            Motifs: 0 for Motifs in InputDictList[0]
            }
        #Populate Empty Dict
        for Frames in FrameList:
            for Motifs in Frames["ID"]:
                if Motifs in ReferenceDictionary:
                    ReferenceDictionary[Motifs] += 1
        return(self.SortDictionaries_AsLists(InputDictionary=ReferenceDictionary, Name="Total"))
             
                      
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
            Concat = pd.concat([IndFrames for IndFrames in geneFiles]).sort_values(by=["Score"], ascending=False).reset_index(drop=True)
            #There's something wrong here
            ConcattedFrames.append(Concat)
        """
        Dictionary used to count the number of motifs appearing in each concatted
        frame
        """
        MotifDictList = []
        for Frames in ConcattedFrames:
            MotifDictionary = {
                Motifs: 0 for Motifs in Frames["ID"]
                }
            for Motifs in Frames["ID"]:
                if Motifs in MotifDictionary:
                    MotifDictionary[Motifs] += 1
            MotifDictList.append(MotifDictionary)
        #Counts per each gene
        OverlappingMotifCounts = self.MotifCountsPerGene(MotifDictList, ConcattedFrames, Init.geneNames)
        #Total counts of motifs
        #SortDictionaries for total counts of motifs
        MotifDictionary = self.TotalMotifCount(InputDictList=OverlappingMotifCounts, FrameList=ConcattedFrames)
        #SortDictionaries for motif counts of each gene
        FrameList = [self.SortDictionaries_AsLists(InputDictionary=OverlappingMotifCounts[Rng], Name=Init.geneNames[Rng])
                     for Rng in range(len(OverlappingMotifCounts))]
        #Extract rows that contain motifs of interest, along with what gene they come from
        ExtractedRows = [(Frames.iloc[Rows,], Init.geneNames[Rng]) for Frames, Rng in zip(ConcattedFrames, range(len(ConcattedFrames))) 
                         for ID, Rows in zip(Frames["ID"], Frames.index.values)
                         if ID in [Motifs for Motifs in MotifDictionary["Total_Motifs"]]]
        ExtractedRowsFrame = pd.DataFrame(data=[ExtractedRows[Vals][0] for Vals in range(len(ExtractedRows))]).reset_index(drop=True)
        ExtractedGeneID = pd.DataFrame(data=[ExtractedRows[Gene_ID][1] for Gene_ID in range(len(ExtractedRows))]).reset_index(drop=True)
        #concat along axis=1 (columns)
        ConcatFrames = pd.concat([ExtractedGeneID, ExtractedRowsFrame], axis=1)
        ConcatFrames = ConcatFrames.rename(columns={0:"Gene Names"})
        ConcatFrames = ConcatFrames.sort_values(by=["Score"], ascending=False).reset_index(drop=True)
        #Concat Motif counts per gene 
        MotifCountsPerGene_Concatted = pd.concat([Frames for Frames in FrameList], axis=1)
        EmptyColumn = pd.Series(" ")
        ConcattedFrame = pd.concat([ConcatFrames, EmptyColumn, MotifCountsPerGene_Concatted, MotifDictionary], axis=1).reset_index(drop=True)
        ConcattedFrame = ConcattedFrame.fillna("")
        """
        Export Concatted Frame
        """   
        return(self.InheritExportFunction(Frame=ConcattedFrame, Location=Init.exportLocation, Name=Init.exportName))
                 
if __name__=='__main__':
    Init = InitializeVars()
    Init.populate(["Gene Names"], Files=[r"Filepath1", r"Filepath2", r"Filepathn"], Export=r"ExportSource", 
                  ExportName = "MotifsOfInterest_FileName")
    Call = TFBSMotifs()
    Call.InheritInit(Init)
        
        
        
        
        
                    
                
        
        
            
            
            
    
                
