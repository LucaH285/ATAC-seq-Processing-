#################
#Import Libraries
#################
#Execute these lines if accessing on a new computer
listOfPackages <- c("TFBSTools", "JASPAR2020", "DBI", "dplyr", "comprehenr", "sjmisc")
newPackages <- listOfPackages[!(listOfPackages %in% installed.packages()[, "Package"])]
if (length(newPackages)) install.packages(newPackages)

#Load packages here
library(TFBSTools)
library(JASPAR2020)
library(DBI)
library(dplyr)
library(comprehenr)
library(sjmisc)
JASPARDb <- file.path(system.file("extdata", package = "JASPAR2020"),
                      "JASPAR2020.sqlite")

JASPARDb2 <- JASPAR2020
#################
#Query for TFBS
#################

MotifScan <- function(Dbase, SpeciesList, ClassList, NameList,
                      TypeList, TaxonomyList, VersionsBoolean){
###################
# Atomic Vectors to be returned at the end of MotifScan
# These vectors contain the fully processed query info
# the JASPAR Sql database
###################
  
  MotifList <- c()
  Species <- c()
  Class <- c()
  Type <- c()
  Taxonomy <- c()
  Names <- c()
  
  # Main query function:
  # Based on the user input for the variables in the function(),
  # The function will return the motifs/TFBS info in an MAnnn.n format (Core JASPAR)
  # Since the input is a list, the output will be a list and corresponds to the input by the following pattern
  # Based on the input order (i.e.: the order of the list)
  #
  #   1. Species 1 --> Class 1, experiment 1, taxonomy 1
  #   2. Species 1 --> Class 2, experiment 2, taxonomy 2
  #   n. Species 1 --> Class n, experiment n, taxonomy n
  #(n+1). Species 2 --> Class 1, experiment 1 taxonomy 1
  #
  # Where the values of class, experiment and taxonomy can be a defined value, 'all' or NA
  #Note: if Species is defined then taxonomy should not be defined and be inputted as NA
  MotifSearch <- function(DataBase, Class, species, Name,
                          type, Taxonomy, Boolversions){
    TFBSList <- c()
    options <- list()
    if (!is.na(Name)){
      options[["name"]] <- Name
      if (species == "all"){
        options[["species"]]
      }else{
        options[["species"]] <- species
      }
      options[["tax_group"]]
      options[["class"]]
      options[["species"]]
      options[["type"]]
    }else{
      if ((is.na(species)) & (!is.na(Taxonomy))){
        if (Taxonomy == "all"){
          options[["tax_group"]]
        }else{
          options[["tax_group"]] <- Taxonomy
        }
      }else if (is.na(Taxonomy) & (!is.na(species))){
        if (species == "all"){
          options[["species"]]
        }else{
          options[["species"]] <- species
        }
      }else if (!is.na(Taxonomy) & (!is.na(species))){
        warning("Both taxonomy and species defined, taxonomy will be ignored. To search by tax \n
             make sure either the species or taxonomy variable is defined as 'NA'")
        options[["species"]] <- species
      }else{
        stop("Both species and taxonomy are undefined, please define either taxonomy or species \n
           if unsure, type '?getMatrixSet' in the R console")
      }
      if (Class == "all"){
        options[["class"]]
      }else{
        options[["class"]] <- Class
      }
      if (type == "all"){
        options[["type"]]
      }else{
        options[["type"]] <- type
      }
    }
    options[["versions"]] <- Boolversions
    BSList <- getMatrixSet(DataBase, options)
    TFBSList <- c(TFBSList, BSList)
    return(TFBSList)
  } 
  #Ensures all atomic vector inputs are of the same length
  MasterList <- list(ClassList, TypeList, TaxonomyList, SpeciesList, NameList)
  MaxAtomic <- max(unlist(lapply(MasterList, length)))
  MinAtomic <- min(unlist(lapply(MasterList, length)))
  if (MaxAtomic != MinAtomic){
    warning("Input vectors are not equal, the last element in the shorter vectors repeated to \n
            match the length of the longest vector.")
  }
  MaxPosition <- which.max(unlist(lapply(MasterList, length)))
  EqualizeLength <- lapply(MasterList, function(x) c(x, rep(x[length(x)], MaxAtomic-length(x))))
  #########################################
  #This loop calls the Search function by the length of the species vector 
  #need not redefine for taxonomy since scan function handles NA
  #########################################
  for (Spec in 1:length(EqualizeLength[[4]])){
    for (Query in 1:length(EqualizeLength[[1]])){
      ScanFxn <- MotifSearch(DataBase = Dbase,
                             species = EqualizeLength[[4]][Spec],
                             Class = EqualizeLength[[1]][Query],
                             type = EqualizeLength[[2]][Query],
                             Taxonomy = EqualizeLength[[3]][Query],
                             Name = EqualizeLength[[5]][Query],
                             Boolversions = VersionsBoolean)
      MotifList <- c(MotifList, ScanFxn[[1]])
      Species <- c(Species, EqualizeLength[[4]][Spec])
      Class <- c(Class, EqualizeLength[[1]][Query])
      Type <- c(Type, EqualizeLength[[2]][Query])
      Taxonomy <- c(Taxonomy, EqualizeLength[[3]][Query])
      Names <- c(Names, EqualizeLength[[5]][Query])
    }
  }
  Frame <- data.frame(
    species = Species,
    class = Class,
    experiment = Type,
    taxonomy = Taxonomy
  )
  ReturnList <- list(MotifList, Frame)
  return(ReturnList)
}

SequenceParse <- function(FASTAList, MotifVector,
                          IdentityScore, Strand,
                          Identifier){
  #############
  #Load Sequences
  #This function prepares the .fasta files as a single string
  #############
  FASTASeq <- function(Files){
    SeqList <- c()
    FastaRead <- read.delim(Files, sep = "\t", header = F)
    for (Rows in 1:length(FastaRead$V1)){
      if (Rows != 1){
        SeqList <- c(SeqList, FastaRead$V1[Rows])
      }
    }
    FASTASeqAsString <- paste(SeqList, collapse = "")
    return(FASTASeqAsString)
  }
  FASTAString <- lapply(FASTAList, FASTASeq)
  SequenceList <- list()
  for (seq in 1:length(FASTAString)){
    #Initialize Sequence 1...
    Seq <- FASTAString[[seq]]
    print(Seq)
    SampleList <- list()
    for (Samples in 1:length(MotifVector)){
      if (length(MotifVector[[Samples]]) != 0){
        print(MotifVector[[Samples]])
        
        #########
        #Need to find a way to do multiple returns in the MotifScan
        #Get Species, Class and Type lists, iterate via the Samples marker
        #Use that as the id for SampleList
        #########
        #Initialize PFMatrixList 1...
        Counter = 1
        Rng <- MotifVector[Samples][[1]]
        #Name <- Species[Samples]
        #class <- Class[Samples]
        MotifList <- list()
        while(Counter <= length(Rng)){
          #Initialize the Motif, MAnnnn.n 1... of PFMatrixList 1... relative to seq 1
          #Run all motifs in PFMatrixList 1, move to the next PFMatrixList 2...
          #MotifList will contain all processed Motifs in a single element from PFMatrixList
          Matrix <- Rng[[Counter]]
          ID <- names(Rng[Counter])
          CnvPWM <- toPWM(Matrix)
          ScanSeq <- searchSeq(CnvPWM, Seq, 
                               seqname = sprintf("Seq_Query_%s", Counter),
                               strand = Strand, min.score = IdentityScore)
          FrameConv <- as(ScanSeq, "data.frame")
          print(FrameConv)
          
          Sys.sleep(2)
          if (length(FrameConv) != 0){
            SubsetVector <- c(FrameConv["start"], FrameConv["end"], FrameConv["relScore"], 
                              FrameConv["ID"], FrameConv["siteSeqs"])
            SubsetFrame <- data.frame(
              Start = SubsetVector$start,
              End = SubsetVector$end,
              Score = SubsetVector$relScore,
              ID = SubsetVector$ID,
              Seq = SubsetVector$siteSeqs
            )
            OrderSubsetFrame <- SubsetFrame[order(SubsetFrame["Score"], decreasing = T),]
            MotifList[[Counter]] <- OrderSubsetFrame
          }else{
            MotifFrame <- data.frame(
              Start = NaN,
              End = NaN,
              Score = NaN,
              ID = NaN,
              Seq = NaN
            )
            MotifList[[Counter]] <- OrderSubsetFrame
          }
          Counter = Counter + 1
        }
        MotifFrame <- data.frame(Reduce(rbind, MotifList))
        MotifFrame <- MotifFrame[order(MotifFrame["Score"], decreasing = T),]
        row.names(MotifFrame) <- NULL
        #SampleList will contain 
        SampleList[[Samples]] <- MotifFrame
        #print(MotifFrame)
      }else{
        MotifFrame <- data.frame(
          Start = NaN,
          End = NaN,
          Score = NaN,
          ID = NaN,
          Seq = NaN
        )
        SampleList[[Samples]] <- MotifFrame
      }
    }#end of motif vector for loop
    SequenceList[[seq]] <- SampleList
  }
  #names(SequenceList[[1]]) <- Identifier
  for (Index in 1:length(SequenceList)){
    if (length(SequenceList[[Index]]) == length(Identifier)){
      names(SequenceList[[Index]]) <- Identifier
    }else if (length(SequenceList[[Index]]) != length(Identifier)){
      stop(cat(", Mismatch in Identifier length, check Identifier atomic vector in function call,", 
               c(cat("ID length:", length(Identifier)), cat(", Vector length:", length(SequenceList[[Index]])))))
    }
  }
  Names <- to_vec(for (Ind in 1:length(SequenceList)) paste("FASTA_seq:", Ind, sep = " "))
  names(SequenceList) <- Names
  return(SequenceList)
}

ProgramCall <- function(DataBase, SpeciesList, ClassList, NameList, 
                        TypeList, TaxonomyList, VersionBoolean,
                        FASTASeqPath, Seq, IDScore, Strand,
                        ExportLocation){
  ####################
  #Call all functions here
  #essentially the if __name__ == __main__: equivalent
  ####################
  
  #Query's JASPAR2020 SQL database for the user provided inputs
  #See the Execute function call
  SQLQuery <- MotifScan(
    Dbase = DataBase,
    SpeciesList = SpeciesList,
    ClassList = ClassList,
    TypeList = TypeList,
    TaxonomyList = TaxonomyList,
    NameList = NameList,
    VersionsBoolean = VersionBoolean
  )
  #####################
  #The second return of the MotifScan function is a 
  #data.frame containing 
  AtomicIdentifier <- c()
  Vector <- SQLQuery[[2]]
  for (IDs in 1:nrow(SQLQuery[[2]])){
    ByRowId <- c()
    for (Cols in 1:length(names(SQLQuery[[2]]))){
      ByRowId <- c(ByRowId, Vector[IDs, names(Vector)[Cols]])
    }
    Collapse <- paste(ByRowId, collapse = "_")
    print(Collapse)
    AtomicIdentifier <- c(AtomicIdentifier, Collapse)
  }
  AtomicFASTASeqs <- c()
  FileList <- list.files(FASTASeqPath)
  for (files in FileList){
    AtomicFASTASeqs <- c(AtomicFASTASeqs, paste(FASTASeqPath, files, sep = "/"))
  }
  SequenceParsing <- SequenceParse(AtomicFASTASeqs, SQLQuery[[1]], IdentityScore = IDScore, 
                                   Strand = Strand, Identifier = AtomicIdentifier)
  ###############
  #Export function
  ###############
  ExportFunction <- function(InputListOfFrames, AtomicIdentifiers, Export){
    Len <- length(InputListOfFrames)
    Rng <- 1
    while (Rng <= Len){
      Char <- paste(Export, paste("InSeq", Rng, sep = "_"), sep = "/")
      if (dir.exists(Char)){
      }else{
        dir.create(Char)
      }
      for (IFrames in 1:length(InputListOfFrames[[Rng]])){
        FrameToExport <- InputListOfFrames[[Rng]][[IFrames]]
        ExportStmnt <- paste(paste("Seq_", Rng, sep = ""), AtomicIdentifiers[IFrames], sep = "_")
        ExportFinal <- paste(paste(Char, ExportStmnt, sep = "/"), ".csv", sep = "")
        write.csv(FrameToExport, file = ExportFinal)
      }
      Rng = Rng + 1
    }
  }
  Export_ <- ExportFunction(InputListOfFrames = SequenceParsing,
                            AtomicIdentifiers = AtomicIdentifier,
                            Export = ExportLocation)
  return(SequenceParsing)
}
Execute <- ProgramCall(DataBase = JASPARDb, SpeciesList = c("Mus musculus", "Homo sapiens"), NameList = c(NA),
            ClassList = c("all", "C2H2 zinc finger factors"), TypeList = c("ChIP-seq"), TaxonomyList = c(NA),
            VersionBoolean = F, FASTASeqPath = "F:/R_dir/Auxilliary/New folder", Strand = "+", 
            IDScore = 0.8, Seq = FASTAFiles, ExportLocation = "F:/R_dir/Auxilliary/New folder")





