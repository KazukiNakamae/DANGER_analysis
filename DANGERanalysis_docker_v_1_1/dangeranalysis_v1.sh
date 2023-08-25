ref_db=$1 # "Dr"
db_type=$2 # "pep"/"cdna"
outpur_dir_name=$3 # "output"
exp_fn="/DATA/${4}" # "ctrl_ko_fltrexpr_contig_tpm_onratio.csv"
assembly="/DATA/${5}" # "/DATA/downsample_fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta"
ontarget_fasta="/DATA/${6}" # "/DATA/guide_pam.fa"
mm=$7 # 11
pam=$8 # "NRR"
thread=16

outpur_dir="/DATA/${outpur_dir_name}"

echo "Make output folder: "${outpur_dir}
if [ ! -d $outpur_dir ]; then
  mkdir $outpur_dir
else
  echo "The output folder already existed!!! Choose a unique folder name."
  echo "DANGER analysis aborts..."
  exit 1;
fi

echo "---------------------------------------------------------------------"
echo "|                  Download databases and software                  |"
echo "---------------------------------------------------------------------"

echo "Download Ensembl database for annotation"
# Download Ensembl database
cd /tmp
mkdir db;
cd db
if [ $ref_db = "Hs" ]; then
  # (Hs)Human
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz;
    db_data="Homo_sapiens.GRCh38.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz;
    db_data="Homo_sapiens.GRCh38.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Mm" ]; then
  # (Mm)Mouse
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz;
    db_data="Mus_musculus.GRCm39.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz;
    db_data="Mus_musculus.GRCm39.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Dm" ]; then
  # (Dm)Fly
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.32.pep.all.fa.gz;
    db_data="Drosophila_melanogaster.BDGP6.32.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz;
    db_data="Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Ce" ]; then
  # (Ce)Worm
  WBGene
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz;
    db_data="Caenorhabditis_elegans.WBcel235.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz;
    db_data="Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Dr" ]; then
  # (Dr)Zebrafish
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz;
    db_data="Danio_rerio.GRCz11.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/current_fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz;
    db_data="Danio_rerio.GRCz11.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Ag" ]; then
  # (Ag)Anopheles
  if [ $db_type = "pep" ]; then
    wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/fasta/anopheles_gambiae/pep/Anopheles_gambiae.AgamP4.pep.all.fa.gz;
    db_data="Anopheles_gambiae.AgamP4.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/fasta/anopheles_gambiae/cdna/Anopheles_gambiae.AgamP4.cdna.all.fa.gz;
    db_data="Anopheles_gambiae.AgamP4.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "At" ]; then
  # (At)Arabidopsis
  if [ $db_type = "pep" ]; then
    wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz;
    db_data="Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz;
    db_data="Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Bt" ]; then
  # (Bt)Bovine
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/pep/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz;
    db_data="Bos_taurus.ARS-UCD1.2.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/cdna/Bos_taurus.ARS-UCD1.2.cdna.all.fa.gz;
    db_data="Bos_taurus.ARS-UCD1.2.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Cf" ]; then
  # (Cf)Canine
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/canis_lupus_familiaris/pep/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz;
    db_data="Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz;
    db_data="Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "EcK12" ]; then
  # (EcK12)E coli strain K12
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_73_collection/escherichia_coli_k_12_gca_004803005/pep/Escherichia_coli_k_12_gca_004803005.ASM480300v1.pep.all.fa.gz;
    db_data="Escherichia_coli_k_12_gca_004803005.ASM480300v1.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_73_collection/escherichia_coli_k_12_gca_004803005/cdna/Escherichia_coli_k_12_gca_004803005.ASM480300v1.cdna.all.fa.gz;
    db_data="Escherichia_coli_k_12_gca_004803005.ASM480300v1.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "EcSakai" ]; then
  # (EcSakai)E coli strain Sakai
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_128_collection/escherichia_coli_o157_h7_str_sakai_gca_000008865/pep/Escherichia_coli_o157_h7_str_sakai_gca_000008865.ASM886v2.pep.all.fa.gz;
    db_data="Escherichia_coli_o157_h7_str_sakai_gca_000008865.ASM886v2.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_128_collection/escherichia_coli_o157_h7_str_sakai_gca_000008865/cdna/Escherichia_coli_o157_h7_str_sakai_gca_000008865.ASM886v2.cdna.all.fa.gz;
    db_data="Escherichia_coli_o157_h7_str_sakai_gca_000008865.ASM886v2.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Gg" ]; then
  # (Gg)Chicken
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz;
    db_data="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/cdna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz;
    db_data="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Mmu" ]; then
  # (Mmu)Rhesus
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz;
    db_data="Macaca_mulatta.Mmul_10.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz;
    db_data="Macaca_mulatta.Mmul_10.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Pt" ]; then
  # (Pt)Chimp
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/pan_troglodytes/pep/Pan_troglodytes.Pan_tro_3.0.pep.all.fa.gz;
    db_data="Pan_troglodytes.Pan_tro_3.0.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/pan_troglodytes/cdna/Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz;
    db_data="Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Rn" ]; then
  # (Rn)Rat
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz;
    db_data="Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/rattus_norvegicus/pep/Rattus_norvegicus.mRatBN7.2.pep.all.fa.gz;
    db_data="Rattus_norvegicus.mRatBN7.2.pep.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Sc" ]; then
  # (Sc)Yeast
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz;
    db_data="Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz;
    db_data="Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Ss" ]; then
  # (Ss)Pig
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/sus_scrofa/pep/Sus_scrofa.Sscrofa11.1.pep.all.fa.gz;
    db_data="Sus_scrofa.Sscrofa11.1.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz;
    db_data="Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Xl" ]; then
  # (Xl)Xenopus *Sequence is Xenopus Tropicalis
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa.gz;
    db_data="Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensembl.org/pub/release-109/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.UCB_Xtro_10.0.cdna.all.fa.gz;
    db_data="Xenopus_tropicalis.UCB_Xtro_10.0.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
elif [ $ref_db = "Mxanthus" ]; then
  # (Mxanthus)Myxococcus xanthus DK 1622
  if [ $db_type = "pep" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_0_collection/myxococcus_xanthus_dk_1622_gca_000012685/pep/Myxococcus_xanthus_dk_1622_gca_000012685.ASM1268v1.pep.all.fa.gz;
    db_data="Myxococcus_xanthus_dk_1622_gca_000012685.ASM1268v1.pep.all.fa.gz"
  elif [ $db_type = "cdna" ]; then
    wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_0_collection/myxococcus_xanthus_dk_1622_gca_000012685/cdna/Myxococcus_xanthus_dk_1622_gca_000012685.ASM1268v1.cdna.all.fa.gz;
    db_data="Myxococcus_xanthus_dk_1622_gca_000012685.ASM1268v1.cdna.all.fa.gz"
  else
    echo "No Database for DANGER analysis. Please check list of available dataset"
    echo "DANGER analysis aborts..."
    exit 1;
  fi
else
  echo "No Database for DANGER analysis. Please check list of available dataset"
  echo "DANGER analysis aborts..."
  exit 1;
fi

echo "done."
echo "-------------------------------"
echo "Download gene database"
eval "$(micromamba shell hook --shell=bash)"
if [ $ref_db = "Hs" ]; then
  # (Hs)Human
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.hs.eg.db=3.16.0;
elif [ $ref_db = "Mm" ]; then
  # (Mm)Mouse
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.mm.eg.db=3.16.0
elif [ $ref_db = "Dm" ]; then
  # (Dm)Fly
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.dm.eg.db=3.16.0
elif [ $ref_db = "Ce" ]; then
  # (Ce)Worm
micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.ce.eg.db=3.16.0
elif [ $ref_db = "Dr" ]; then
  # (Dr)Zebrafish
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.dr.eg.db=3.16.0
elif [ $ref_db = "Ag" ]; then
  # (Ag)Anopheles
micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.ag.eg.db=3.16.0
elif [ $ref_db = "At" ]; then
  # (At)Arabidopsis
micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.at.tair.db=3.16.0
elif [ $ref_db = "Bt" ]; then
  # (Bt)Bovine
micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.bt.eg.db=3.16.0
elif [ $ref_db = "Cf" ]; then
  # (Cf)Canine
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.cf.eg.db=3.16.0
elif [ $ref_db = "EcK12" ]; then
  # (EcK12)E coli strain K12
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.eck12.eg.db=3.16.0
elif [ $ref_db = "EcSakai" ]; then
  # (EcSakai)E coli strain Sakai
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.ecsakai.eg.db=3.16.0
elif [ $ref_db = "Gg" ]; then
  # (Gg)Chicken
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.gg.eg.db=3.16.0
elif [ $ref_db = "Mmu" ]; then
  # (Mmu)Rhesus
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.mmu.eg.db=3.16.0
elif [ $ref_db = "Pt" ]; then
  # (Pt)Chimp
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.pt.eg.db=3.16.0
elif [ $ref_db = "Rn" ]; then
  # (Rn)Rat
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.rn.eg.db=3.16.0
elif [ $ref_db = "Sc" ]; then
  # (Sc)Yeast
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.sc.sgd.db=3.16.0
elif [ $ref_db = "Ss" ]; then
  # (Ss)Pig
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.ss.eg.db=3.16.0
elif [ $ref_db = "Xl" ]; then
  # (Xl)Xenopus *Sequence is Xenopus Tropicalis
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.xl.eg.db=3.16.0
elif [ $ref_db = "Mxanthus" ]; then
  # (Mxanthus)Myxococcus xanthus DK 1622
  micromamba install --quiet -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-org.mxanthus.db=1.0.27
else
  echo "No Database for DANGER analysis. Please check list of available dataset"
  echo "DANGER analysis aborts..."
  exit 1;
fi

echo "done."
echo "-------------------------------"
echo "Download and Install Crisflash"
cd /tmp
git clone https://github.com/crisflash/crisflash.git;
cd /tmp/crisflash;
make;
cd /tmp;

echo "done."
echo "-------------------------------"
echo "Download and Install TransDecoder"
micromamba create --quiet -y -n transdecoder_env;
micromamba install --quiet -y -n transdecoder_env -c anaconda -c conda-forge -c bioconda -y TransDecoder=5.5.0 pigz=2.6 blast=2.12.0 seqkit=2.3.1 gzip=1.12;

echo "done."
echo "-------------------------------"
echo "Download and Install fasta3"
micromamba create --quiet -y -n fasta3_env;
micromamba install --quiet -y -n fasta3_env -c anaconda -c conda-forge -c bioconda fasta3=36.3.8 time=1.8 gzip=1.12 pigz=2.6;



echo "---------------------------------------------------------------------"
echo "| Detection of on/off-target sites on de novo transcriptome assembly |"
echo "---------------------------------------------------------------------"

echo "done."
echo "-------------------------------"
echo "Detect on/off-target sites on de novo transcriptome assembly"
# on/off-target detection
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' $assembly > temp.fa; # Remove \n in sequences
cat temp.fa | sed -e 's/ .*//g' > Result_denovo_transcriptome_assembly.fa; # Modify the headers
/tmp/crisflash/bin/crisflash \
-g Result_denovo_transcriptome_assembly.fa \
-s ${ontarget_fasta} \
-o ${outpur_dir}/Result_gRNAs.cas-offinder \
-m $mm \
-p $pam -t $thread -C \
>${outpur_dir}/Log_crisflash.txt 2>&1;
wc -l ${outpur_dir}/Result_gRNAs.cas-offinder;
wc -l ${outpur_dir}/Result_gRNAs.cas-offinder > ${outpur_dir}/Summary_Count_of_All_On_and_Off-target_Sites.txt;
# hit site (upto 11mm, NRR PAM) = 4232634 sites

echo "done."
echo "-------------------------------"
echo "Extract all off-target sites on de novo transcriptome assembly"
if [ "$(uname)" == 'Darwin' ]; then
  ontarget_seq=$(sed 2p ${ontarget_fasta} | echo)
  awk '{ if ($4 != "'$ontarget_seq'") { print } }' ${outpur_dir}/Result_gRNAs.cas-offinder > ${outpur_dir}/Result_offtarget_all.cas-offinder;
elif [ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]; then
  ontarget_seq=$(sed -n 2p ${ontarget_fasta})
  awk '$4 != "'$ontarget_seq'"' ${outpur_dir}/Result_gRNAs.cas-offinder > ${outpur_dir}/Result_offtarget_all.cas-offinder;
else
  echo "Your platform ($(uname -a)) is not supported."
  exit 1
fi
wc -l ${outpur_dir}/Result_offtarget_all.cas-offinder;
wc -l ${outpur_dir}/Result_offtarget_all.cas-offinder > ${outpur_dir}/Summary_Count_of_All_Off-target_Sites.txt;
# off-target site (upto 11mm, NRR PAM) = 4232633 sites

echo "done."
echo "-------------------------------"
echo "Extract off-target sites for each mismatches on de novo transcriptome assembly"
# Get 0-11 off-targets（NRR PAM）
for i in {0..11};
  do awk -v MM=$i '{ if ($6 == MM) { print } }' ${outpur_dir}/Result_offtarget_all.cas-offinder > ${outpur_dir}/Result_offtarget_mm"$i".cas-offinder;
  if [ -f "${outpur_dir}/Result_offtarget_mm${i}.cas-offinder" ]; then
    wc -l ${outpur_dir}/Result_offtarget_mm"$i".cas-offinder;
    wc -l ${outpur_dir}/Result_offtarget_mm"$i".cas-offinder > ${outpur_dir}/Summary_Count_of_mm"$i"_Off-target_Sites.txt;
  fi
done;
###
#       0 Result_offtarget_mm1.cas-offinder
#       0 Result_offtarget_mm2.cas-offinder
#       1 Result_offtarget_mm3.cas-offinder
#      25 Result_offtarget_mm4.cas-offinder
#     278 Result_offtarget_mm5.cas-offinder
#    2179 Result_offtarget_mm6.cas-offinder
#   14172 Result_offtarget_mm7.cas-offinder
#   73027 Result_offtarget_mm8.cas-offinder
#  296222 Result_offtarget_mm9.cas-offinder
# 1008070 Result_offtarget_mm10.cas-offinder
# 2838659 Result_offtarget_mm11.cas-offinder
###

echo "done."
echo "-------------------------------"
echo "Count transcripts have off-target sites on de novo transcriptome assembly"
awk '{ print $2 }' ${outpur_dir}/Result_offtarget_all.cas-offinder > ${outpur_dir}/Result_offtarget_all_list.txt;
awk '!seen[$0]++' ${outpur_dir}/Result_offtarget_all_list.txt > ${outpur_dir}/Result_offtarget_all_uniq_list.txt;
wc -l ${outpur_dir}/Result_offtarget_all_uniq_list.txt;
wc -l ${outpur_dir}/Result_offtarget_all_uniq_list.txt > ${outpur_dir}/Summary_Count_of_All_Off-target_Transcripts.txt;
# 865,452 transcripts have off-target sites



echo "---------------------------------------------------------------------"
echo "|                      GO enrichment analysis                      |"
echo "---------------------------------------------------------------------"

if [ $db_type = "pep" ]; then
  echo "-------------------------------"
  echo "Predict ORFs using TransDecoder"
  cd /tmp;
  micromamba activate transdecoder_env;
  cp ${assembly} Trinity.fasta;
  chmod +x ./11TransDecoder.sh;
  ./11TransDecoder.sh
  target_seq_file=Trinity.fasta.transdecoder.pep
elif [ $db_type = "cdna" ]; then
  target_seq_file=Trinity.fasta
else
  echo "Unexpected Input"
  echo "DANGER analysis aborts..."
  exit 1;
fi

echo "done."
echo "-------------------------------"
echo "Add dummy word into blank in database"
unzip_db_data=${db_data%.*}
unpigz -c db/${db_data} > $unzip_db_data;
chmod +x 00_prepare_faa_4Fanflow.sh;
./00_prepare_faa_4Fanflow.sh ${unzip_db_data};
micromamba deactivate;

echo "done."
echo "-------------------------------"
echo "Run GGSEARCH"
micromamba activate fasta3_env;
chmod +x ./15ggsearch.sh;
modified_db_data=${unzip_db_data%.*}.fa2
./15ggsearch.sh ${target_seq_file} ${modified_db_data};
echo "..."
micromamba deactivate

echo "done."
echo "-------------------------------"
echo "Make gene annotation table"
chmod +x 15parseggsearch_linux.sh;
./15parseggsearch_linux.sh ${modified_db_data};
perl -nle 'print $1 if(/^\>(\S+)/)' ${target_seq_file} > ${outpur_dir}/coding-transcript.pid.txt;
cat ${outpur_dir}/coding-transcript.pid.txt | perl 15mkannotbl.pl ${target_seq_file}-${modified_db_data}.txt > ${target_seq_file}_all.txt;
cp ${modified_db_data} ${outpur_dir}/`basename ${modified_db_data}`;
cp ${target_seq_file} ${outpur_dir}/${target_seq_file};
cp ${target_seq_file}-${modified_db_data}.txt ${outpur_dir}/${target_seq_file}-${modified_db_data}.txt;
cp ${target_seq_file}_all.txt ${outpur_dir}/${target_seq_file}_all.txt;

echo "done."
echo "-------------------------------"
echo "Search downregulated transcripts using TPM ratio"
awk -F',' '{ if ($NF ~ /downregulated/) { print } }' ${exp_fn} > ${outpur_dir}/downregulated_ONratio.csv;
awk -F',' '{ print $1 }' ${outpur_dir}/downregulated_ONratio.csv > ${outpur_dir}/Result_downregulated_ONratio_list.txt;
wc -l ${outpur_dir}/Result_downregulated_ONratio_list.txt;
wc -l ${outpur_dir}/Result_downregulated_ONratio_list.txt > ${outpur_dir}/Summary_Count_of_All_downregulate_TPMratio_Transcripts.txt;
# downregulated_ONratio = 439908

echo "done."
echo "-------------------------------"
echo "Search dTPM transcripts"
micromamba activate matplotlib_venn_env;
# Make Venn diagram using https://bioinformatics.psb.ugent.be/webtools/Venn/
###
# List names	number of elements	number of unique elements
# 1（off-target）	865452	865452
# 2（downrregulated ONratio）	439908	439908
# Overall number of unique elements	935252
###
# all_offtarget_vs_dONrartio.svg
# all_offtarget_vs_dONrartio.png
# all_offtarget_vs_dONrartio.txt
# downrregulated ONratio & off-target n = 370,108
# Saved it to all_offtarget_vs_dONrartio.id.txt
# Modify all_offtarget_vs_dONrartio.id.txt
python drawVennDiagram.py ${outpur_dir}/all_offtarget_vs_dONrartio.id.txt \
  ${outpur_dir}/all_offtarget_vs_dONrartio.tiff \
  ${outpur_dir}/Result_offtarget_all_uniq_list.txt \
  ${outpur_dir}/Result_downregulated_ONratio_list.txt;
wc -l ${outpur_dir}/all_offtarget_vs_dONrartio.id.txt;
wc -l ${outpur_dir}/all_offtarget_vs_dONrartio.id.txt > ${outpur_dir}/Summary_Count_of_All_downregulate_TPMratio_offtarget_Transcripts.txt;
micromamba deactivate

# TODO:各生物IDに書き換える
echo "done."
echo "-------------------------------"
echo "Perform enrichment analysis in group of dTPM"
while read line;
  do grep $line ${outpur_dir}/coding-transcript.pid.txt >> ${outpur_dir}/all_offtarget_vs_dONrartio.pid.txt;
  done < ${outpur_dir}/all_offtarget_vs_dONrartio.id.txt;
wc -l ${outpur_dir}/all_offtarget_vs_dONrartio.pid.txt # ORFありtranscript 2500
cat ${outpur_dir}/all_offtarget_vs_dONrartio.pid.txt | perl 15mkannotbl.pl ${outpur_dir}/${target_seq_file}-${modified_db_data}.txt > ${outpur_dir}/all_offtarget_vs_dONrartio_annotation_table.txt;
awk '$2!=""' ${outpur_dir}/all_offtarget_vs_dONrartio_annotation_table.txt > temp_all_offtarget_vs_dONrartio_annotation_table.txt_hasAnnotation.txt;
awk '{ if ($2 !~ /non/) { print } }' temp_all_offtarget_vs_dONrartio_annotation_table.txt_hasAnnotation.txt > ${outpur_dir}/all_offtarget_vs_dONrartio_annotation_table.txt_hasAnnotation.txt;
wc -l ${outpur_dir}/all_offtarget_vs_dONrartio_annotation_table.txt_hasAnnotation.txt;
wc -l ${outpur_dir}/all_offtarget_vs_dONrartio_annotation_table.txt_hasAnnotation.txt > ${outpur_dir}/Summary_Count_of_All_geneidentified_downregulate_TPMratio_offtarget_Transcripts.txt;
# 1246 transcripts have Annotation ID

echo "done."
echo "-------------------------------"
echo "Search downregulated transcripts for each mismatch number using TPM ratio"
for i in {0..11};
  do awk '{ print $2 }' ${outpur_dir}/Result_offtarget_mm"$i".cas-offinder > ${outpur_dir}/offtarget_mm"$i"_list.txt;
  awk '!seen[$0]++' ${outpur_dir}/offtarget_mm"$i"_list.txt > ${outpur_dir}/offtarget_mm"$i"_uniq_list.txt;
  if [ -f "${outpur_dir}/offtarget_mm${i}_uniq_list.txt" ]; then
    wc -l ${outpur_dir}/offtarget_mm"$i"_uniq_list.txt;
    wc -l ${outpur_dir}/offtarget_mm"$i"_uniq_list.txt > ${outpur_dir}/Summary_Count_of_transcripts_have_mm"$i"_Off-target.txt;
  fi
done;
###
#        0 offtarget_mm0_uniq_list.txt
#        0 offtarget_mm1_uniq_list.txt
#        0 offtarget_mm2_uniq_list.txt
#        1 offtarget_mm3_uniq_list.txt
#       25 offtarget_mm4_uniq_list.txt
#      278 offtarget_mm5_uniq_list.txt
#     2149 offtarget_mm6_uniq_list.txt
#    13227 offtarget_mm7_uniq_list.txt
#    60910 offtarget_mm8_uniq_list.txt
#   198645 offtarget_mm9_uniq_list.txt
#   475284 offtarget_mm10_uniq_list.txt
#   786219 offtarget_mm11_uniq_list.txt
###

echo "done."
echo "-------------------------------"
echo "Search dTPM transcripts for each mismatch number using TPM ratio"
# Search downrregulated ONratio & off-target for each mismatch number
mkdir ${outpur_dir}/mm_offtarget_dONratio;
# Draw Venn diagram
micromamba activate matplotlib_venn_env;
for i in {0..11};
  do python drawVennDiagram.py ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt \
  ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.tiff \
  ${outpur_dir}/offtarget_mm"$i"_uniq_list.txt \
  ${outpur_dir}/Result_downregulated_ONratio_list.txt;
  if [ -f "${outpur_dir}/mm_offtarget_dONratio/mm${i}_offtarget_dONratio.id.txt" ]; then
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt;
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt > ${outpur_dir}/Summary_Count_of_mm"$i"_downregulate_TPMratio_offtarget_Transcripts.txt;
  fi
done;
### ONratio & off-target for each mismatch number
#        0 mm_offtarget_dONratio/mm0_offtarget_dONratio.id.txt
#        0 mm_offtarget_dONratio/mm1_offtarget_dONratio.id.txt
#        0 mm_offtarget_dONratio/mm2_offtarget_dONratio.id.txt
#        1 mm_offtarget_dONratio/mm3_offtarget_dONratio.id.txt
#        6 mm_offtarget_dONratio/mm4_offtarget_dONratio.id.txt
#       62 mm_offtarget_dONratio/mm5_offtarget_dONratio.id.txt
#      500 mm_offtarget_dONratio/mm6_offtarget_dONratio.id.txt
#     3489 mm_offtarget_dONratio/mm7_offtarget_dONratio.id.txt
#    18172 mm_offtarget_dONratio/mm8_offtarget_dONratio.id.txt
#    67658 mm_offtarget_dONratio/mm9_offtarget_dONratio.id.txt
#   182774 mm_offtarget_dONratio/mm10_offtarget_dONratio.id.txt
#   330733 mm_offtarget_dONratio/mm11_offtarget_dONratio.id.txt
###
micromamba deactivate;

echo "done."
echo "-------------------------------"
echo "Search gene-annotated dTPM transcripts for each mismatch number using TPM ratio"
# Search transcript annotated with Annotation
for i in {0..11};do
  while read line;
    do grep $line ${outpur_dir}/coding-transcript.pid.txt >> ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt;
  done < ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt;
  if [ -f "${outpur_dir}/mm_offtarget_dONratio/mm${i}_offtarget_dONratio.pid.txt" ]; then
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt;
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt > ${outpur_dir}/Summary_Count_of_mm"$i"_downregulate_TPMratio_offtarget_ORFs.txt;
  fi
done;
### transcript with ORFs
# wc: mm_offtarget_dONratio/mm0_offtarget_dONratio.pid.txt: open: No such file or directory
# wc: mm_offtarget_dONratio/mm1_offtarget_dONratio.pid.txt: open: No such file or directory
# wc: mm_offtarget_dONratio/mm2_offtarget_dONratio.pid.txt: open: No such file or directory
#        0 mm_offtarget_dONratio/mm3_offtarget_dONratio.pid.txt
#        0 mm_offtarget_dONratio/mm4_offtarget_dONratio.pid.txt
#        1 mm_offtarget_dONratio/mm5_offtarget_dONratio.pid.txt
#       19 mm_offtarget_dONratio/mm6_offtarget_dONratio.pid.txt
#       85 mm_offtarget_dONratio/mm7_offtarget_dONratio.pid.txt
#      412 mm_offtarget_dONratio/mm8_offtarget_dONratio.pid.txt
#     1170 mm_offtarget_dONratio/mm9_offtarget_dONratio.pid.txt
#     2055 mm_offtarget_dONratio/mm10_offtarget_dONratio.pid.txt
#     2469 mm_offtarget_dONratio/mm11_offtarget_dONratio.pid.txt
###
for i in {0..11};
  do if [ -f "${outpur_dir}/mm_offtarget_dONratio/mm${i}_offtarget_dONratio.pid.txt" ]; then
    cat ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt | perl 15mkannotbl.pl ${outpur_dir}/${target_seq_file}-${modified_db_data}.txt > ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt;
  fi
done;
for i in {0..11};
  do awk '$2!=""' ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt > temp_mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt;
  awk '{ if ($2 !~ /non/) { print } }' temp_mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt > ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt;
  if [ -f "${outpur_dir}/mm_offtarget_dONratio/mm${i}_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt" ]; then
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt;
    wc -l ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt > ${outpur_dir}/Summary_Count_of_mm"$i"_geneidentified_downregulate_TPMratio_offtarget_Transcripts.txt;
  fi
done;
### transcript with Annotation ID
#        0 mm_offtarget_dONratio/mm3_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#        0 mm_offtarget_dONratio/mm4_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#        1 mm_offtarget_dONratio/mm5_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#       11 mm_offtarget_dONratio/mm6_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#       50 mm_offtarget_dONratio/mm7_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#      212 mm_offtarget_dONratio/mm8_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#      567 mm_offtarget_dONratio/mm9_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#     1032 mm_offtarget_dONratio/mm10_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
#     1228 mm_offtarget_dONratio/mm11_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt
###

echo "done."
echo "-------------------------------"
echo "Calculate D-index"
echo "(1/4)Get a gene list for each mismatch number..."
# Get gene list for each mismatch number
for i in {0..11};
  do awk '{ print $3}' ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasAnnotation.txt \
  > ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_genelist.txt;
done;
echo "(2/4)GO annotation with topGO for each mismatch number..."
# GO annotation with topGO for each mismatch number
micromamba activate topGO;
for i in {0..11};
  do Rscript annotateGOv2.R ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_genelist.txt \
  ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_geneGO \
  $ref_db;
done;
echo "(3/4)Enrichment analysis for each mismatch number..."
# Enrichment analysis for each mismatch number
for i in {0..11};
  do Rscript makeDANGERenrichmentTablev2.R ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_geneGO \
  ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_geneGO_enrichment \
  ${outpur_dir}/mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_Annotation_genelist.txt $i;
done;
micromamba deactivate;
echo "(4/4)Calculation of D-index..."
# DENGERindex Calculation
micromamba activate calcDANGERindex_env;
cd ${outpur_dir};
python /tmp/calcDANGERindex_v2.py ${outpur_dir}/mm_offtarget_dONratio/mm{0..11}_offtarget_dONratio_Annotation_geneGO_enrichment > ${outpur_dir}/Summary_Dindex.txt;
cat ${outpur_dir}/Summary_Dindex.txt;
###
# Biological Process
# Sum of DANGER indexes
# 224.24425603474646
# Upper boundary
# 0.14470317255959442
# Cellular Component
# Sum of DANGER indexes
# 68.85219343072643
# Upper boundary
# 0.16767265997982977
# Molecular Function
# Sum of DANGER indexes
# 60.94544580544128
# Upper boundary
# 0.13143181335175058
###
micromamba deactivate;


echo "DANGER analysis was successfully done!"
exit 0;