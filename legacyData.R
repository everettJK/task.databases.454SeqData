library(RMySQL)
library(stringr)
library(pander)
library(dplyr)
options(useFancyQuotes = FALSE)
sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn1    <- dbConnect(MySQL(), group='specimen_management')
sampleData <- dbGetQuery(dbConn1, "select * from gtsp")

samplesClause <- paste0('sampleName like ', paste0(sapply(sampleData$SpecimenAccNum, function(x){ sQuote(paste0(x, '-%'))  }), collapse=' or sampleName like '))

dbConn2  <- dbConnect(MySQL(), group='intsites_miseq')
intSiteSamples <- unique(gsub('\\-\\d+$', '', unname(unlist(dbGetQuery(dbConn2, sprintf("select sampleName from samples where %s", samplesClause))))))
samplesWithoutMiseqData <- sampleData$SpecimenAccNum[! sampleData$SpecimenAccNum %in% intSiteSamples]
samplesWithoutMiseqData <- samplesWithoutMiseqData[! nchar(samplesWithoutMiseqData) < 8]

sapply(dbListConnections(MySQL()), dbDisconnect)

illuminaSets <- toupper(gsub('"', '', scan('all_illumina_set_names_in_454db.tsv', what = 'character', sep = '\n')))

d <- bind_rows(lapply(samplesWithoutMiseqData, function(GTSP){
  dbConn  <- dbConnect(MySQL(), group='454intsites')
  
  sql <- sprintf("select distinct setname from genetherapy_samples where setname like '%%%s%%'", GTSP)
  setNames <- unname(unlist(dbGetQuery(dbConn, sql)))
  setNames <- setNames[! toupper(setNames) %in% illuminaSets]
  
  if(length(setNames) > 0){
    r <- bind_rows(lapply(setNames, function(setName){
           sql <- paste0("select psl, cloneCount from sets_psl where cloneCount >= 1 and name = '", setName, "'")
           psls <- dbGetQuery(dbConn, sql)
  
           if(nrow(psls) > 0){
              sql <- paste0("select id, strand, restriction, freeze, blockCount, qSize, tName, tStart, tEnd from psl where id in (", paste0(sQuote(psls$psl), collapse = ','), ")")
              d <- suppressWarnings(dbGetQuery(dbConn, sql))
              if(! is.data.frame(d)) d <- data.frame()
    
              if(nrow(d) > 0){ 
                d$GTSP <- GTSP
                d$estAbund <- psls[match(d$id, psls$psl),]$cloneCount
                d$setName  <- setName
                d$posid <- paste0(d$tName, d$strand, ifelse(d$strand == '+', d$tStart, d$tEnd))
                return(d)
              } else {
                return(data.frame())
              }
           } else {
             return(data.frame())
           }
    }))
  } else {
    r <- data.frame()
  }
  
  dbDisconnect(dbConn)
  if(nrow(r) == 0) return(data.frame())
  
  r$patient   <- sampleData[match(GTSP, sampleData$SpecimenAccNum),]$Patient
  r$cellType  <- sampleData[match(GTSP, sampleData$SpecimenAccNum),]$CellType
  r$timePoint <- sampleData[match(GTSP, sampleData$SpecimenAccNum),]$Timepoint
    
  return(r)
}))

d <- dplyr::group_by(d, GTSP, posid) %>% dplyr::mutate(estAbund = sum(estAbund), setName = paste0(setName, collapse = ',')) %>% dplyr::slice(1) %>% dplyr::ungroup()
d$intSitePosition <- ifelse(d$strand == '+', d$tStart, d$tEnd)
g <- GenomicRanges::makeGRangesFromDataFrame(d, seqnames.field = 'tName', start.field = 'intSitePosition', end = 'intSitePosition', keep.extra.columns = TRUE)
mcols(g) <- mcols(g)[,c("GTSP", "patient", "cellType", "timePoint", "estAbund", "posid", "setName")]

saveRDS(g, file='Bushman454seqData.rds')
