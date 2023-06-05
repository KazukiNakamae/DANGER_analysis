input.prefix <- commandArgs(trailingOnly=TRUE)[1] #1番目の引数を取得する
output.dir <- commandArgs(trailingOnly=TRUE)[2] #2番目の引数を取得する
genelist.fn <- commandArgs(trailingOnly=TRUE)[3] #3番目の引数を取得する
mm <- commandArgs(trailingOnly=TRUE)[4] #4番目の引数を取得する

suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

if(!file.exists(file.path(input.prefix, "GO.BP.map")) &&
 !file.exists(file.path(input.prefix, "GO.CC.map")) &&
 !file.exists(file.path(input.prefix, "GO.MF.map"))){
  message(paste0("---", input.prefix, " does not have a GO map.---"))
  quit(save = "no")
}

if((file.size(file.path(input.prefix, "GO.BP.map")) == 0L) ||
  (file.size(file.path(input.prefix, "GO.CC.map")) == 0L) ||
  (file.size(file.path(input.prefix, "GO.MF.map")) == 0L)){
  message(paste0("---", input.prefix, " does not have sufficient GO data for D-index calculation.---"))
  quit(save = "no")
}

message(paste0("---", input.prefix, " is analyzed.---"))

dir.create(output.dir)

#readMappingsを使ってバックグラウンド（step1で準備）を読み込む。
gene.bp.info <- readMappings(file = file.path(input.prefix, "GO.BP.map"), sep = "\t", IDsep = ",")
gene.cc.info <- readMappings(file = file.path(input.prefix, "GO.CC.map"), sep = "\t", IDsep = ",")
gene.mf.info <- readMappings(file = file.path(input.prefix, "GO.MF.map"), sep = "\t", IDsep = ",")

# gene.info <- readMappings(file = "mm_offtarget_dONratio/mm11_offtarget_dONratio_ENSP_geneGO/GO.BP.map", sep = "\t", IDsep = ",")
#gene name
bp.geneNames <- names(gene.bp.info)
cc.geneNames <- names(gene.cc.info)
mf.geneNames <- names(gene.mf.info)

#DEGsなどのクエリのリスト（step2）を読み込む（今回は特にフィルタは設けない）
interested_.list <- unlist(as.list(read.csv(genelist.fn, header = FALSE)))

bp.geneList <- factor(as.integer(c(bp.geneNames, "dummygene") %in% interested_.list)) # 全て一致するとフォーマットエラーになるためダミーを入れておく
names(bp.geneList) <- c(bp.geneNames, "dummygene")
cc.geneList <- factor(as.integer(c(cc.geneNames, "dummygene") %in% interested_.list)) # 全て一致するとフォーマットエラーになるためダミーを入れておく
names(cc.geneList) <- c(cc.geneNames, "dummygene")
mf.geneList <- factor(as.integer(c(mf.geneNames, "dummygene") %in% interested_.list)) # 全て一致するとフォーマットエラーになるためダミーを入れておく
names(mf.geneList) <- c(mf.geneNames, "dummygene")


#これで準備ができた。GODataオブジェクトを構築（ここではBPだがMFとCCも試す）
#nodeSize = 10 は、10未満の語彙からGO階層を刈り取るために使用される。
bp.GOdata <- new("topGOdata", ontology = "BP", allGenes = bp.geneList, annot = annFUN.gene2GO, gene2GO = gene.bp.info)
cc.GOdata <- new("topGOdata", ontology = "CC", allGenes = cc.geneList, annot = annFUN.gene2GO, gene2GO = gene.cc.info)
mf.GOdata <- new("topGOdata", ontology = "MF", allGenes = mf.geneList, annot = annFUN.gene2GO, gene2GO = gene.mf.info)

analyzeEnrichment <- function(GOdata, type, mm, prefix){
  #topGOdataクラスのオブジェクトを得たらエンリッチメント分析を始めることができる。
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

  #エンリッチメントテスト実施後、結果を分析し解釈するため、GenTableを使う。GenTable は最も有意な GO タームとそれに対応する p 値を分析する。
  enrichment_table <- GenTable(GOdata, resultKS.elim = resultKS.elim,  orderBy = "weight", ranksOf = "classic", topNodes = length(attributes(resultKS.elim)$score))

  enrichment_table["mm"] <- rep(mm, nrow(enrichment_table))

  if(nrow(enrichment_table)>0){
    write.table(enrichment_table, file = file.path(prefix, paste0(type, "_enrichment_table.txt")), append = FALSE, quote = FALSE, sep = "\t",
    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
    col.names = FALSE, qmethod = c("double"),
    fileEncoding = "utf-8")
  }
}
analyzeEnrichment(bp.GOdata, "BP", mm, output.dir)
analyzeEnrichment(cc.GOdata, "CC", mm, output.dir)
analyzeEnrichment(mf.GOdata, "MF", mm, output.dir)
