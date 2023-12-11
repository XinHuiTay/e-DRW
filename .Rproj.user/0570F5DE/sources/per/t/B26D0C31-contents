# Get the risk pathway
riskPathway <- c("05165", "04512", "03320")

getGenesInClassifierKegg <- function(dataframes, LungNormSample, LungDiseaseSample, riskPathway)
{
  # Get the genes in the classifier
  # input:
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # mRNA_matrix_training: the pathway expression profile of the training set
  # LungNormSample: the index of norm samples in the training set
  # LungDiseaseSample: the index of disease samples in the training set
  # riskPathway: the selected risk active pathways for classification
  # output:
  # the genes in the classifier

  mRNA_matrix_training <- dataframes[[1]]
  # t-test
  tScore <- calculatetScore(mRNA_matrix_training, LungNormSample, LungDiseaseSample)
  vertexTScore <- getVertexTScore(tScore, DirectGraph, rownames(mRNA_matrix_training))

  riskGene <- c()

  for(i in 1:length(riskPathway))
  {
    # find the index of riskPathway in gNonMetabolic or gMetabolic
    NonMeta <- 0
    Meta <- 0
    idxPath <- which(names(gNonMetabolic) == riskPathway[i])
    #browser()
    if (length(idxPath) > 0)
    {
      NonMeta <- 1
    }else
    {
      idxPath <- which(names(gMetabolic) == riskPathway[i])
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
      Vpathwayi <- get.vertex.attribute(gNonMetabolic[[idxPath]], "name", index=V(gNonMetabolic[[idxPath]]))
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
      Vpathwayi <- get.vertex.attribute(gMetabolic[[idxPath]], "name", index=V(gMetabolic[[idxPath]]))
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
riskGenes <- getGenesInClassifierKegg(dataframes, LungNormSample, LungDiseaseSample, riskPathway)
#riskGenes
