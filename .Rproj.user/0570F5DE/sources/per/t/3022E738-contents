# Get the risk pathway
riskPathway <- c("108", "204", "243")

getGenesInClassifierPid <- function(data_frames, LungNormSample, LungDiseaseSample, riskPathway)
{
  # Get the genes in the classifier
  # input:
  # gPid2: nonMetabolic pathway
  # gPid1: metabolic pathway
  # mRNA_matrix_training: the pathway expression profile of the training set
  # LungNormSample: the index of norm samples in the training set
  # LungDiseaseSample: the index of disease samples in the training set
  # riskPathway: the selected risk active pathways for classification
  # output:
  # the genes in the classifier

  mRNA_matrix_training <- data_frames[[1]]
  # t-test
  tScore <- TScoreCalculation(mRNA_matrix_training, LungNormSample, LungDiseaseSample)
  vertexTScore <- vertexTScoreCalculation(tScore, PidGraph, rownames(mRNA_matrix_training))

  riskGene <- c()

  for(i in 1:length(riskPathway))
  {
    # find the index of riskPathway in gPid2 or gPid1
    NonMeta <- 0
    Meta <- 0
    idxPath <- which(names(gPid2) == riskPathway[i])
    #browser()
    if (length(idxPath) > 0)
    {
      NonMeta <- 1
    }else
    {
      idxPath <- which(names(gPid1) == riskPathway[i])
      if(length(idxPath) > 0)
      {
        Meta <- 1
      }else
      {
        stop("No such pathway!")
      }
    }

    # the current riskPathway is nonMetabolic pathway
    if(NonMeta == 1)
    {
      Vpathwayi <- get.vertex.attribute(gPid2[[idxPath]], "name", index=V(gPid2[[idxPath]]))
      if (length(Vpathwayi) > 0)
      {
        for (j in 1 : length(Vpathwayi))
        {
          idxGene <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
          if (length(idxGene > 0))
          {
            if ( rownames(mRNA_matrix_training)[idxGene] %in% rownames(vertexTScore))
            {
              if(vertexTScore[rownames(mRNA_matrix_training)[idxGene],2] < 0.05)
              {
                riskGene <- c(riskGene, rownames(mRNA_matrix_training)[idxGene])
              }
            }
          }
        }
      }
    }

    # the current riskPathway is Metabolic pathway
    if(Meta == 1)
    {
      Vpathwayi <- get.vertex.attribute(gPid1[[idxPath]], "name", index=V(gPid1[[idxPath]]))
      if (length(Vpathwayi) > 0)
      {
        for (j in 1 : length(Vpathwayi))
        {
          idxGene <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
          if (length(idxGene > 0))
          {
            if ( rownames(mRNA_matrix_training)[idxGene] %in% rownames(vertexTScore))
            {
              if(vertexTScore[rownames(mRNA_matrix_training)[idxGene],2] < 0.05)
              {
                riskGene <- c(riskGene, rownames(mRNA_matrix_training)[idxGene])
              }
            }
          }
        }
      }
    }

  }
  return(riskGene)
}
riskGenes <- getGenesInClassifierPid(data_frames, LungNormSample, LungDiseaseSample, riskPathway)
#riskGenes
