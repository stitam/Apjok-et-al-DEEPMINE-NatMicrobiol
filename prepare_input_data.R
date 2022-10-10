root <- getwd()

# ILLUMINA

# download run1 illumina data to run1/illumina
if (!dir.exist(run1/illumina)) dir.create(run1/illumina)
setwd(run1/illumina)
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/ERR10323306/ERR10323306_R1.fastq.gz -O Undetermined_S0_L001_R1_001.fastq.gz")
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/ERR10323306/ERR10323306_R2.fastq.gz -O Undetermined_S0_L001_R2_001.fastq.gz")
setwd(root)

# download run2 illumina data to run2/illumina
if (!dir.exist(run2/illumina)) dir.create(run2/illumina)
setwd(run2/illumina)
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/ERR10355437/ERR10355437_R1.fastq.gz -O Undetermined_S0_L001_R1_001.fastq.gz")
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/ERR10355437/ERR10355437_R2.fastq.gz -O Undetermined_S0_L001_R2_001.fastq.gz")
setwd(root)

# download run3 illumina data to run3/illumina
if (!dir.exist(run3/illumina)) dir.create(run3/illumina)
setwd(run3/illumina)
system("wget ")

# uptag
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR992/ERR9920964/ERS12338933.R1.fastq.gz -O Uptag.R1.fastq.gz")
R.utils::gunzip("ERS12338933.R1.fastq.gz")

system("wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR992/ERR9920964/ERS12338933.R2.fastq.gz -O Uptag.R2.fastq.gz")
R.utils::gunzip("ERS12338933.R2.fastq.gz")

# downtag
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR992/ERR9920963/ERS12338932.R1.fastq.gz -O Downtag.R1.fastq.gz")
R.utils::gunzip("ERS12338932.R1.fastq.gz")
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR992/ERR9920963/ERS12338932.R2.fastq.gz -O Downtag.R2.fastq.gz")
R.utils::gunzip("ERS12338932.R2.fastq.gz")

setwd(root)

# NANOPORE

# unpack run1 nanopore contigs to run1/nanopore
setwd("run1")
unzip("nanopore.zip")
setwd(root)

# unpack run2 nanopore contigs to run2/nanopore
setwd("run2")
unzip("nanopore.zip")
setwd(root)