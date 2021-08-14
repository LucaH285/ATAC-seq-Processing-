# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:41:27 2021
Script to create FASTA files of regions of interest
@author: Luca Hategan
"""
import gzip
import os
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import argparse
from abc import ABC, abstractmethod
import pandas as pd
from ast import literal_eval
#Time is useful for debugging
import time

class InitializeVars:
    """
    Initialization Class:
        Defines all variables that will be used throught the various
        classes and subclasses of this program. 
        Variables are populated in the __main__ if statement based on
        user inputs from the command line
    """
    def __init__(self):
        self.file_list = []
        self.GeneRegion = []
        self.NameList = []
        self.ChromosomeList = []
        self.StrandID = []
        self.PotentialCCREList = []
        self.ReferenceAnimal = ""
        self.Export = ""
        self.repeatAreas = ""
        self.repeatThreshold = 0
        self.GeneName = ""
        self.KnownCCREs = []
   
    def populate(self, Files, GeneRegions, Names, Chrs, Strands, CCREs,
                 PcCRE, RefAn, Exprt, Repeat, SequentialLimit, Gene):
        self.file_list = Files
        self.GeneRegion = GeneRegions
        self.NameList = Names
        self.ChromosomeList = Chrs
        self.StrandID = Strands
        self.PotentialCCREList = PcCRE
        self.ReferenceAnimal = RefAn
        self.Export = Exprt
        self.repeatAreas = Repeat
        self.repeatThreshold = SequentialLimit
        self.GeneName = Gene
        self.KnownCCREs = CCREs
   
class ReverseComplement:
    """
    This is a single class that defines a reverse
    complement function. In the event that the input
    gene is defined on the negative sense strand, ("-"), this
    function will show the reverse complement based on Watson-Crick base pairing:
        i.e.: 3' ATCGGTA 5' --> 5' TACCGAT 3'
    """
    def ReverseComplementFxn(self, InputROI):
        if type(InputROI) is list:
            InputROI = InputROI[0]
        Reversed = []
        NucVector = ["A", "T", "C", "G"]
        RevComplement = ["T", "A", "G", "C"]
        for Nucleotides in InputROI[::-1]:
            if Nucleotides in NucVector:
                Index = NucVector.index(Nucleotides)
                ReverseNucleotide = RevComplement[Index]
                Reversed.append(ReverseNucleotide)
            elif ((Nucleotides.isupper() is False) and ((Nucleotides != "N") or (Nucleotides != "n"))):
                Index = NucVector.index(Nucleotides.upper())
                ReverseNucleotide = RevComplement[Index]
                Reversed.append(ReverseNucleotide)
            elif ((Nucleotides == "N") or (Nucleotides == "n")):
                Reversed.append(Nucleotides.upper())
        return("".join(Reversed))
    
class AdjustCaseOfSequence:
    """
    This class defines a function that adjusts the case to be all uppercase
    in the FASTA output.
    """
    def SeqParser(self, InputROI):
        if type(InputROI) is list:
            InputROI = InputROI[0]
        SequenceList = []
        for Nucleotides in InputROI:
            if Nucleotides.isupper() is False:
                SequenceList.append(Nucleotides.upper())
            else:
                SequenceList.append(Nucleotides)
        return("".join(SequenceList))
    
class BubblSorter:
    def BubbleSortFxn(self, InputList):
        Condition = True
        while(Condition):
            Condition = False
            for Ind in range(len(InputList) - 1):
                if InputList[Ind] > InputList[Ind + 1]:
                    InputList[Ind], InputList[Ind + 1] = InputList[Ind + 1], InputList[Ind]
                    Condition = True
        return(InputList)
    
class ClustalOmegaCLProtocol:
    def ClustalOmegaCLFxn(self, InputFile, OutputFile):
        if os.path.exists(InputFile):
            clustalomega_cline = ClustalOmegaCommandline(infile=InputFile, outfile=OutputFile, verbose=True, outfmt = "clustal", auto=False)
            return(clustalomega_cline)
        elif os.path.exists(InputFile) is False:
            return("Check input file. Unaligned FASTA does not exist")
        
class InitializeVarsInheritable(ABC):
    @abstractmethod
    def Vars(self, Init):
        pass
    
class RevCompInheritable(ABC):
    @abstractmethod
    def RevComp(self):
        pass

class SeqAdjustInheritable(ABC):
    @abstractmethod
    #The SeqP var can be replaced by simply calling the 
    #Class within a different function
    def SeqParsefxn(self):
        pass
    
class InheritableSort(ABC):
    @abstractmethod
    def Sort(self):
        pass
    
class ClustalInheritable(ABC):
    @abstractmethod
    def Clustal(self, CL):
        pass
    
class ExportPcCREOnly(InitializeVarsInheritable, RevCompInheritable, 
                      SeqAdjustInheritable, InheritableSort):
    def RevComp(self, InputSeq):
        Reverse = ReverseComplement().ReverseComplementFxn(InputROI=InputSeq)
        return(Reverse)
    
    def SeqParsefxn(self, InputSeq):
        SeqFxn = AdjustCaseOfSequence().SeqParser(InputROI=InputSeq)
        return(SeqFxn)
    
    def Sort(self, List):
        SortFxn = BubblSorter().BubbleSortFxn(InputList=List)
        return(SortFxn)
    
    def SortCCREs(self, Input, Threshold):
        Ranges = []
        Condition = True
        for X, Y in zip(Input[:-1], Input[1:]):
            Diff = abs(Y - X)
            if (Diff <= Threshold) and (Condition is True):
                Ranges.append(X)
                Condition = False
            elif (Diff > Threshold) and (Condition is False):
                Ranges.append(X)
                Condition = True
            elif (Y == Input[len(Input) - 1]):
                Ranges.append(Y)
        return(Ranges)
        
    def Vars(self, Init):
        """
        Preprocessing region, handles .csv or manual PCCRE entry.
        .csv files are handles in two ways, taking the raw PCCREs as is
        or creating a range of ccres based on their proximity to eachother.
        Proximity is defined by the user and should be close in value to the
        bin size defined in the ATACseq processor (i.e.: if bin size = 100, threshold >100, <=1000; depending on size of query region)
        """
        if len(Init.PotentialCCREList) != 0:
            global PcCRERanges
            if Init.PotentialCCREList[0].endswith(".csv") is True:
                PCCREs = []
                for Lists in Init.PotentialCCREList:
                    Frame = pd.read_csv(Lists)
                    PCCREs.append([int(StartCoor) for StartCoor in Frame["Genome Range Start Coordinate"]])
                #Takes the list of lists object of PCCREs and creates a single
                #list containing all the PCCREs
                UnpackPCCRE = []
                for Lists in range(len(PCCREs)):
                    UnpackPCCRE += PCCREs[Lists]
                if Init.repeatAreas == "y":
                    PcCREList = list(set(UnpackPCCRE))
                    sortPcCREs = self.Sort(PcCREList)
                    PcCRERanges = self.SortCCREs(Input=sortPcCREs, Threshold=Init.repeatThreshold)
                    PCCRERanges_X = []
                    PCCRERanges_Y = []
                    for coords in range(len(PcCRERanges)):
                        if coords % 2 == 0:
                            PCCRERanges_X.append(PcCRERanges[coords])
                        else:
                            PCCRERanges_Y.append(PcCRERanges[coords])
                elif (Init.repeatAreas == "n") or (Init.repeatAreas is None):
                    PcCRERanges = list(set(UnpackPCCRE))
                    PCCRERanges_X = []
                    PCCRERanges_Y = []
                    for coords in range(len(PcCRERanges)):
                        if coords % 2 == 0:
                            PCCRERanges_X.append(PcCRERanges[coords])
                        else:
                            PCCRERanges_Y.append(PcCRERanges[coords])
                else:
                    raise(TypeError("Error in repeat areas entry must be: 'y' or 'n'"))
            else:
                if Init.repeatAreas == "y":
                    SortPcCREs = self.Sort(List = Init.PotentialCCREList)
                    PcCRERanges = self.SortCCREs(Input = SortPcCREs, Threshold=Init.repeatThreshold)
                    PCCRERanges_X = []
                    PCCRERanges_Y = []
                    for coords in range(len(PcCRERanges)):
                        if coords % 2 == 0:
                            PCCRERanges_X.append(PcCRERanges[coords])
                        else:
                            PCCRERanges_Y.append(PcCRERanges[coords])
                elif (Init.repeatAreas == "n") or (Init.repeatAreas is None):
                    PcCRERanges = self.Sort(List = Init.PotentialCCREList)
                    PCCRERanges_X = []
                    PCCRERanges_Y = []
                    for coords in range(len(PcCRERanges)):
                        if coords % 2 == 0:
                            PCCRERanges_X.append(PcCRERanges[coords])
                        else:
                            PCCRERanges_Y.append(PcCRERanges[coords])
                else:
                    raise(TypeError("Error in repeat areas entry must be: 'y' or 'n'"))
        else:
            PcCRERanges = []
        if len(Init.KnownCCREs) != 0:
            if (isinstance(Init.KnownCCREs[0], str) is True):
                KnownCCRE_X = []
                KnownCCRE_Y = []
                for Files in Init.KnownCCREs:
                    Frame = pd.read_csv(Files)
                    try:
                        for Coors in Frame["CCRE Coordinates"]:
                            KnownCCRE_X.append(literal_eval(Coors)[0])
                            KnownCCRE_Y.append(literal_eval(Coors)[1])
                    except KeyError:
                        try:
                            for CoorsX, CoorsY in zip(Frame["CCRE start coordinates"], Frame["CCRE end coordinates"]):
                                KnownCCRE_X.append(CoorsX)
                                KnownCCRE_Y.append(CoorsY)
                        except KeyError:
                            raise(KeyError("If the .csv files you are using are not outputs from the ATAC-seq processor, please name the columns as 'CCRE start coordinates' and 'CCRE end coordiantes', respectively, or manually enter the known ccre values"))
                SetX = list(set(KnownCCRE_X))
                SetY = list(set(KnownCCRE_Y))
                SortKnownX = self.Sort(List = SetX)
                SortKnownY = self.Sort(List = SetY)
                SortKnown = [(X, Y) for X, Y in zip(SortKnownX, SortKnownY)]
            else:
                Set = list(set(Init.KnownCCREs))
                SortKnown = self.Sort(List=Set)
                KCCRE_X = [SortKnown[Coords] for Coords in range(len(SortKnown)) if Coords % 2 == 0]
                KCCRE_Y = [SortKnown[Coords] for Coords in range(len(SortKnown)) if Coords % 2 != 0]
                SortKnown = [(X, Y) for X, Y in zip(KCCRE_X, KCCRE_Y)]
        else:
            SortKnown = []
                
        """
        Here is the main sequence export function for the 
        PCCREs
        """
        if Init.ReferenceAnimal in Init.NameList:
            Locate = Init.NameList.index(Init.ReferenceAnimal)
        RecordsListExport = []
        KnownCCRERecords = []
        with gzip.open(Init.file_list[Locate], "rt") as Handle:
            #This assigns the Dict produced here as a global variable,
            #meaning that it should be used in the Mouse CCRE generator
            #instead of recreating it each time the script is executed
            print("Initiated {0} PCCRE and/or CCRE FASTA export".format(Init.NameList[Locate]))
            global Dict
            Dict = SeqIO.to_dict(SeqIO.parse(Handle, "fasta"))
            #Let this be the potential CCRE export
            if (len(PCCRERanges_X) != 0) and (len(PCCRERanges_Y) != 0):
                for X, Y, rng in zip(PCCRERanges_X, PCCRERanges_Y, range(len(PCCRERanges_Y))):
                    ROI = Dict[Init.ChromosomeList[Locate]][X:Y]
                    if Init.StrandID[Locate] == "-":
                        Reverse = self.RevComp(str(ROI.seq))
                        PcCRERecord = SeqRecord(Seq(Reverse),
                                                id = "{0}_for_{1}_PotentialCCRE_{2}".format(Init.GeneName, Init.NameList[Locate], rng))
                        RecordsListExport.append(PcCRERecord)
                    elif Init.StrandID[Locate] == "+":
                        Case = self.SeqParsefxn(InputSeq=ROI.seq)
                        PcCRERecord = SeqRecord(Seq(Case), 
                                                id = "{0}_for_{1}_PotentialCCRE_{2}".format(Init.GeneName, Init.NameList[Locate], rng))
                        RecordsListExport.append(PcCRERecord)
            if len(SortKnown) != 0:
                for rng in range(len(SortKnown)):
                    ROI = Dict[Init.ChromosomeList[Locate]][SortKnown[rng][0]:SortKnown[rng][1]]
                    if Init.StrandID[Locate] == "-":
                        Reverse = self.RevComp(str(ROI.seq))
                        KCCRERecord = SeqRecord(Seq(Reverse),
                                                id = "{0}_for_{1}_KnownCCRE_{2}".format(Init.GeneName, Init.NameList[Locate], rng))
                        KnownCCRERecords.append(KCCRERecord)
                    elif Init.StrandID[Locate] == "+":
                        Case = self.SeqParsefxn(InputSeq=ROI.seq)
                        KCCRERecord = SeqRecord(Seq(Case),
                                                id = "{0}_for_{1}_KnownCCRE_{2}".format(Init.GeneName, Init.NameList[Locate], rng))
                        KnownCCRERecords.append(KCCRERecord)
        ExportLocation = Init.Export + "\PcCRE_FastaFiles"
        ExportLocationKnownCCRE = Init.Export + "\KcCRE_FastaFiles"
        if not os.path.exists(ExportLocation) and (len(RecordsListExport) != 0):
            os.makedirs(ExportLocation)
            for Records in range(len(RecordsListExport)):
                with open(ExportLocation + "\PCCRE_Fasta_{0}.fasta".format((PCCRERanges_X[Records], PCCRERanges_Y[Records])), "w") as Handle:
                    SeqIO.write(RecordsListExport[Records], Handle, "fasta")
        else:
            print("PcCRE_FastaFiles folder already exists, passing")
        if not os.path.exists(ExportLocationKnownCCRE) and len(KnownCCRERecords) != 0:
            os.makedirs(ExportLocationKnownCCRE)
            for Records in range(len(KnownCCRERecords)):
                with open(ExportLocationKnownCCRE + "\KCCRE_Fasta_{0}.fasta".format(SortKnown[Records]), "w") as Handle:
                    SeqIO.write(KnownCCRERecords[Records], Handle, "fasta")
        else:
            print("KcCRE_FastaFiles folder exists, passing")
            
#This class creates a FASTA file of the PcCRE's only, matching the length
#of the Query region. (i.e.: if the Query region is from coordinate x to y, the
#length of the FASTA file will be y-x)
#The nucleotides that are not a part of the PcCRE regions are filled in as 'N'
class MousePcCREGenerator(InitializeVarsInheritable, RevCompInheritable, SeqAdjustInheritable):
    def RevComp(self, InputSeq):
        Reverser = ReverseComplement().ReverseComplementFxn(InputROI=InputSeq)
        return(Reverser)
    
    def SeqParsefxn(self, InputSeq):
        SeqFxn = AdjustCaseOfSequence().SeqParser(InputROI=InputSeq)
        return(SeqFxn)
         
    def Vars(self, Init):
        SequenceVector = []
        """
        Add logic here, locate mouse/input animal, default to mouse using the NameList
        location should match the input list.
        """
        if Init.ReferenceAnimal in Init.NameList:
            Locate = Init.NameList.index(Init.ReferenceAnimal)
        """
        This function will be the same as the reference animal and has been called in the class above,
        thus there is no need to redefine the dictionary, and hence Dict has been assigned as a global variable.
        """
        print("Initialized {0}, at {1}, strand {2}, and range {3} file for alignment".format(Init.NameList[Locate], Init.ChromosomeList[Locate], Init.StrandID[Locate], Init.GeneRegion[Locate]))
        if "Dict" not in globals():
            with gzip.open(Init.file_list[Locate], "rt") as Handle:
                Dict2 = SeqIO.to_dict(SeqIO.parse(Handle, "fasta"))
        else:
            pass
        Counter = 0
        Conditional = True
        #This is structured wrong
        if "PcCRERanges" in globals():
            AdjustedPCCRETemp = [[], []]
            for i in range(len(PcCRERanges)):
                if i % 2 == 0:
                    AdjustedPCCRETemp[0].append(PcCRERanges[i])
                else:
                    AdjustedPCCRETemp[1].append(PcCRERanges[i])
            for Rng in range(Init.GeneRegion[Locate][0], Init.GeneRegion[Locate][1]):
                if Rng == AdjustedPCCRETemp[0][Counter]:
                    if AdjustedPCCRETemp[1][Counter] > Init.GeneRegion[Locate][1]:
                        AdjustedPCCRETemp[1][Counter] = Init.GeneRegion[Locate][1]
                    if "Dict" in globals():
                        ROI = Dict[Init.ChromosomeList[Locate]][AdjustedPCCRETemp[0][Counter]:AdjustedPCCRETemp[1][Counter]]
                    else:
                        ROI = Dict2[Init.ChromosomeList[Locate]][AdjustedPCCRETemp[0][Counter]:AdjustedPCCRETemp[1][Counter]]
                    if Init.StrandID[Locate] == "-":
                        if Conditional is True:
                            print("Input Strand '-' sense, initiated reverse complement")
                            Conditional = False
                        ROI_Rev = self.RevComp(InputSeq=ROI.seq)
                        SequenceVector.append(ROI_Rev)
                    else:
                        AdjustSeq = self.SeqParsefxn(InputSeq=ROI.seq)
                        SequenceVector.append(AdjustSeq)
                elif Rng < AdjustedPCCRETemp[0][Counter] or Rng >= AdjustedPCCRETemp[1][Counter]:
                    SequenceVector.append("N")
                    if ((Rng == AdjustedPCCRETemp[1][Counter]) and (Counter < (len(AdjustedPCCRETemp[0]) - 1))):
                        print(AdjustedPCCRETemp)
                        Counter += 1
            PCCREs = "".join(SequenceVector)
            #Its going over the size because you have coordinates larger than the range
            global PCCRERecord
            PCCRERecord = SeqRecord(Seq(PCCREs), id = "MousePotentialCCRE")
            if os.path.exists(Init.Export + "\MousePotentialCCRE.fasta") is False:
                with open(Init.Export + "\MousePotentialCCRE.fasta", "w") as Handle:
                    SeqIO.write(PCCRERecord, Handle, "fasta")
            else:
                print("MousePotentialCCRE.fasta exits, passing")
            return(PCCREs)
        else:
            return("No PcCREs called")
    
class ExportAllRegions(InitializeVarsInheritable, RevCompInheritable, SeqAdjustInheritable):
    def RevComp(self, InputSeq):
        Reverser = ReverseComplement().ReverseComplementFxn(InputROI=InputSeq)
        return(Reverser)
    
    def SeqParsefxn(self, InputSeq):
        SeqFxn = AdjustCaseOfSequence().SeqParser(InputROI=InputSeq)
        return(SeqFxn)
    
    def Vars(self, Init):
        #Export the sequence as s FASTA file of the Animal and fa.gz inputs, sequentially
        #Need the Sequences of the whole Query range
        RecordsList = []
        for FileInd in range(len(Init.file_list)):
            print("Initialized {0}, at {1}, strand {2} and range {3}".format(Init.NameList[FileInd], Init.ChromosomeList[FileInd], Init.StrandID[FileInd], Init.GeneRegion[FileInd]))
            with gzip.open(Init.file_list[FileInd], "rt") as Handle:
                Dict = SeqIO.to_dict(SeqIO.parse(Handle, "fasta"))
                ROI = Dict[Init.ChromosomeList[FileInd]][Init.GeneRegion[FileInd][0]:Init.GeneRegion[FileInd][1]]
                if Init.StrandID[FileInd] == "-":
                    print("Initialized Reverse Complement for {}".format(Init.NameList[FileInd]))
                    Reverse = self.RevComp(InputSeq=ROI.seq)
                    Record = SeqRecord(Seq(Reverse), id = "{}".format(Init.NameList[FileInd]))
                    RecordsList.append(Record)
                    if os.path.exists(Init.Export + "\FASTAFile_{0}, {1}".format(Init.NameList[FileInd], FileInd)) is False:
                        with open(Init.Export + "\FASTAFile_{0}, {1}".format(Init.NameList[FileInd], FileInd), "w") as Handle:
                            SeqIO.write(Record, Handle, "fasta")
                    else:
                        print("FASTAFile_{0}, {1} exits, passing".format(Init.NameList[FileInd], FileInd))
                elif Init.StrandID[FileInd] == "+":
                    CheckCase = self.SeqParsefxn(InputSeq=ROI.seq)
                    ROI = CheckCase
                    Record = SeqRecord(Seq(ROI), id = "{}".format(Init.NameList[FileInd]))
                    RecordsList.append(Record)
                    if os.path.exists(Init.Export + "\FASTAFile_{0}, {1}".format(Init.NameList[FileInd], FileInd)) is False:
                        with open(Init.Export + "\FASTAFile_{0}, {1}".format(Init.NameList[FileInd], FileInd), "w") as Handle:
                            SeqIO.write(Record, Handle, "fasta")
                    else:
                        print("FASTAFile_{0}, {1} exits, passing".format(Init.NameList[FileInd], FileInd))
        #Create the unaligned file here, uses a global variable
        if "PCCRERecord" in globals():
            RecordsList.insert(len(RecordsList), PCCRERecord)
            try:
                MultipleRecords = MultipleSeqAlignment(
                    RecordsList
                    )
            except ValueError:
                raise()
            if os.path.exists(Init.Export + "\Pre_alignedFile.fasta") is False:
                with open(Init.Export + "\Pre_alignedFile.fasta", "w") as Handle:
                    SeqIO.write(MultipleRecords, Handle, "fasta")
            else:
                print("Unaligned file exists, passing")
            return(MultipleRecords)
    
class ClustalOutput(InitializeVarsInheritable, ClustalInheritable):
    def Vars(self, Init):
        return(Init)
    
    def Clustal(self, CL):
        if "PCCRERecord" in globals():
            InFile = Init.Export + "\Pre_alignedFile.fasta"
            OutFile = Init.Export + "\AlignedFile"
            clustalomega_cline = ClustalOmegaCommandline(infile=InFile, outfile=OutFile, verbose=True, outfmt="clustal", auto=False)
            return(print(clustalomega_cline))
        else:
            return(print("No PcCREs provided, function inactive"))
        
    
if __name__ == '__main__':
    def CSVConverter(Input, Seperation=",", convType = str):
        #Should raise an error here
        Converted = []
        for Vals in Input.split(Seperation):
            Converted.append(convType(Vals))
        return(Converted)
    
    def CSVConvertIntegers(Input, Seperation=",", convType = int):
        Converted = []
        if (Input.endswith(".csv") or Input.endswith(".txt")):
            Fxn = [str(Vals) for Vals in Input.split(",")]
            return(Fxn)
        else:
            for Vals in Input.split(Seperation):
                Converted.append(convType(Vals))
        return(Converted)
  
    parser = argparse.ArgumentParser(description="FASTA file writer, Version 2.0")
    parser.add_argument("--filelist", type=CSVConverter, dest = "FileList", required=True, 
                        help="The whole animal genome files, in fa.gz format only. When inputting files, take care to add them in the following way exactly: 'FilePath1','FilePath2',etc The position of the comma between file paths matters, do not space the file paths")
    parser.add_argument("--AnimalList", type=CSVConverter, dest = "AnimalList", required=True,
                        help="List containing the animal names, make sure to match the format here to the format of the animal of interest, as before the position of the comma matters")
    parser.add_argument("--Gene", type=str, dest="Genename", required=True, help="The name of the gene of interest (i.e.: Gad2)")
    parser.add_argument("--ChromosomeList", type=CSVConverter, dest = "ChrList", required=True,
                        help="Chromosome list, enter as: 'chr1','chr2','chrn' or depending on your reference genome: 'scaffold_109','chrB4'")
    parser.add_argument("--StrandList", type=CSVConverter, dest = "StrandList", required=True,
                        help="Strand list, input as: '+','-' -> no other options. Note if - sense, a reverse complement of the input sequence will be produced")
    parser.add_argument("--GeneRegionStartCoor", type=CSVConvertIntegers, dest = "GRStart", required=True,
                        help="Start coordinates of the Query region, entered as integers: 1,2,3")
    parser.add_argument("--GeneRegionEndCoor", type=CSVConvertIntegers, dest = "GREnd", required = True,
                        help="End coordinates of the Query region, entered as integers: 1,2,3")
    parser.add_argument("--ExprtFolder", metavar="PATH", required=True, dest="ExportFolder", help="The export location for which all FASTA files will be stored. Enter as: C:\SubFolder\ExportFolder")
    parser.add_argument("--RefAnimal", type=str, dest="RefAnimal", required = True, help = "The primary animal of interest from which the PcCREs are derived")
    parser.add_argument("--PcCREs", type=CSVConvertIntegers, dest = "PcCREs", required = False, default = [],
                        help="Coordinates of the PcCREs, entered as integers: 1,2,3 or as a .csv file (use the potential ccres discovered from the ATAC-seq processor)")
    parser.add_argument("--KnownCCRE", type=CSVConvertIntegers, dest="KnownCCRE", required = False, default = [],
                        help = "Coordinates of the known CCREs, entered as integers, 1,2,3 or as a .csv file (use the known ccres exported from the ATAC-seq processor)")
    parser.add_argument("--Sequential", type=str, dest="Seq", required = False, 
                        help="enter as: 'y' or 'n'. This command will group sequential sequential PcCREs into 1 large region. Sequential PcCREs are defined by the user as the minimum number of nucleotides seperating sorted PcCREs, default is 300nt",
                        default="n")
    parser.add_argument("--SeqThreshold", type=int, dest="SeqThreshold", required=False, default=300,
                        help="The sequential threshold value described above, enter as an integer")
    parser.add_argument("--ShowClustal", type=str, dest="ClstlCL", required = False, help = "Show the Clustal command line executable? 'y' or 'n'", default="n")
    args = parser.parse_args()

    Init = InitializeVars()
    
    Coords = args.PcCREs
    KCCRE = args.KnownCCRE
        
    Init.populate(Files= args.FileList,
                  GeneRegions = [(X, Y) for X, Y in zip(args.GRStart, args.GREnd)],
                  Names = args.AnimalList,
                  Chrs = args.ChrList,
                  Strands = args.StrandList,
                  PcCRE = Coords, CCREs = KCCRE,
                  RefAn=args.RefAnimal, Exprt= args.ExportFolder,
                  Repeat=args.Seq, SequentialLimit=args.SeqThreshold, Gene=args.Genename)
    
    if len(args.PcCREs) != 0:
        Exp = ExportPcCREOnly()
        Exp.Vars(Init)
    
        Exec = MousePcCREGenerator()
        Exec.Vars(Init)
    
    
    All = ExportAllRegions()
    All.Vars(Init)
    
    
    if args.ClstlCL == "y":
        CL = ClustalOmegaCLProtocol()
        Clstl = ClustalOutput()
        Clstl.Vars(Init)
        Clstl.Clustal(CL)

    
        
    
    
    
    
    
        

