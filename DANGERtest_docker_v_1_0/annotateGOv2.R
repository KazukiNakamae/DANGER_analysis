
args1 = commandArgs(trailingOnly=TRUE)[1] #1番目の引数を取得する
args2 = commandArgs(trailingOnly=TRUE)[2] #2番目の引数を取得する
args3 = commandArgs(trailingOnly=TRUE)[3] #3番目の引数を取得する

if(!file.exists(args1)){
  message(paste0("---", args1, " does not exist.---"))
  quit(save = "no")
}

dir.create(args2)

if(file.exists(args1) && file.size(args1) != 0L){
  if(args3 == "Hs"){
    suppressMessages(require(org.Hs.eg.db))
    db <- org.Hs.eg.db
  }else if(args3 == "Ag"){
    suppressMessages(require(org.Ag.eg.db))
    db <- org.Ag.eg.db
  }else if(args3 == "At"){
    suppressMessages(require(org.At.tair.db))
    db <- org.At.tair.db
  }else if(args3 == "Bt"){
    suppressMessages(require(org.Bt.eg.db))
    db <- org.Bt.eg.db
  }else if(args3 == "Ce"){
    suppressMessages(require(org.Ce.eg.db))
    db <- org.Ce.eg.db
  }else if(args3 == "Cf"){
    suppressMessages(require(org.Cf.eg.db))
    db <- org.Cf.eg.db
  }else if(args3 == "Dm"){
    suppressMessages(require(org.Dm.eg.db))
    db <- org.Dm.eg.db
  }else if(args3 == "Dr"){
    suppressMessages(require(org.Dr.eg.db))
    db <- org.Dr.eg.db
  }else if(args3 == "EcK12"){
    suppressMessages(require(org.EcK12.eg.db))
    db <- org.EcK12.eg.db
  }else if(args3 == "EcSakai"){
    suppressMessages(require(org.EcSakai.eg.db))
    db <- org.EcSakai.eg.db
  }else if(args3 == "Gg"){
    suppressMessages(require(org.Gg.eg.db))
    db <- org.Gg.eg.db
  }else if(args3 == "Mm"){
    suppressMessages(require(org.Mm.eg.db))
    db <- org.Mm.eg.db
  }else if(args3 == "Mmu"){
    suppressMessages(require(org.Mmu.eg.db))
    db <- org.Mmu.eg.db
  }else if(args3 == "Mxanthus"){
    suppressMessages(require(org.Mxanthus.db))
    db <- org.Mxanthus.db
  }else if(args3 == "Pt"){
    suppressMessages(require(org.Pt.eg.db))
    db <- org.Pt.eg.db
  }else if(args3 == "Rn"){
    suppressMessages(require(org.Rn.eg.db))
    db <- org.Rn.eg.db
  }else if(args3 == "Sc"){
    suppressMessages(require(org.Sc.sgd.db))
    db <- org.Sc.sgd.db
  }else if(args3 == "Ss"){
    suppressMessages(require(org.Ss.eg.db))
    db <- org.Ss.eg.db
  }else if(args3 == "Xl"){
    suppressMessages(require(org.Xl.eg.db))
    db <- org.Xl.eg.db
  }else{
    message(paste0("---There is no database on ", args3, ".---"))
    quit(save = "no")
  }

  gene.list <- read.table(args1, header = FALSE)[,1]
  # gene.list <- read.table("mm_offtarget_dONratio/mm5_offtarget_dONratio_ENSP_genelist.txt", header = FALSE)[,1] # debug
  # keytypes(db)
  # org.Hs.egGO time stump: https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

  # アノテーションを検索
  tryCatch({
    GO.table <- select(db,
      keys = gene.list,
      columns = c('SYMBOL','GO'),
      keytype = 'SYMBOL')
    }, error = function(e) {    # e にはエラーメッセージが保存されている
      message(paste0("---There is no GO annotation of ", toString(gene.list), ".---"))
      quit(save = "no")
    },
    silent = TRUE
  )
  GO.table <- GO.table[complete.cases(GO.table), ] # NA削除

  # アノテーションを分類
  molecular.function.GO.table <- GO.table[GO.table$ONTOLOGY == 'MF',]
  cellular.component.GO.table <- GO.table[GO.table$ONTOLOGY == 'CC',]
  biological.process.GO.table <- GO.table[GO.table$ONTOLOGY == 'BP',]

  # アノテーションを整理
  mapGO <- function(GO.table){
    symbol.uni.vec <- sort(unique(GO.table$SYMBOL))
    GO.map.table <- data.frame(query=character(0), GOs=character(0))
    for (symbol in symbol.uni.vec){
      if(is.na(GO.table[GO.table$SYMBOL == symbol,]$GO[1])){
        break
      }
      if(length(GO.table[GO.table$SYMBOL == symbol,]$GO)>0){
        GOs <- paste(GO.table[GO.table$SYMBOL == symbol,]$GO, collapse=',')
      }else{
        GOs <- '-'
      }
      GO.map.table <- rbind(GO.map.table, data.frame(query=symbol, GOs=GOs))
    }
    return(GO.map.table)
  }

  GO.map.table <- mapGO(GO.table)
  molecular.function.GO.map.table <- mapGO(molecular.function.GO.table)
  cellular.component.GO.map.table <- mapGO(cellular.component.GO.table)
  biological.process.GO.map.table <- mapGO(biological.process.GO.table)

  # 書き込み
  write.table(GO.map.table, file = file.path(args2, "GO.map"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("double"),
              fileEncoding = "utf-8")

  write.table(molecular.function.GO.map.table, file = file.path(args2, "GO.MF.map"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("double"),
              fileEncoding = "utf-8")

  write.table(cellular.component.GO.map.table, file = file.path(args2, "GO.CC.map"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("double"),
              fileEncoding = "utf-8")

  write.table(biological.process.GO.map.table, file = file.path(args2, "GO.BP.map"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("double"),
              fileEncoding = "utf-8")
  
  message(paste0("---", args1, " is analyzed.---"))
}else{
  message(paste0("---Ignore ", args1, ".---"))
}

quit(save = "no")


