#############################################################################
## File Nmae: cnvkit_get_coverages.wdl
## Created By: Zach Weber
## Created On: 2019-05-20
## Version: 0.1
## Purpose: subworkflow level script for import into larger .wdl workflows.
##            Both reference generation and CNV calling processes require these
##            initial steps and therefore its easiest to pull these out for
##            later import to keep the size of scripts low.
##
#############################################################################

####################################################
##### define the cnvkit_get_coverages workflow #####
####################################################
workflow cnvkit_get_coverages {
  #------------------------------------
  # define metadata
  meta {
    author: "Zach Weber"
    email: "zach.weber.813@gmail.com"
    description: "Implementing the CNVkit pipeline through 'coverage' step"
  }

  #------------------------------------
  # define global variables
  ## define reference files
  File baitset # baits for WXS run. (.bed) this workflow assumes all samples run on same assay.
  File refgene # refgene file for exon level annotation (usually .txt)
  File gaps # gaps in genome assembly from UCSC  (.bed)
  File assembly_fa # genome assembly fasta (.fa)

  #** a note on chromosome coherence: Naming schemes between baits, refgene files,
  #** gaps files, assemblys and BAM files should be consistent. either use
  #** numerics only (1-22,X) or characters only (chr1-chr22, chrX)

  ## set this option to run samtools_fix_chrs to convert from numerics to chars
  Boolean fix_chrs = true

  ## set disk and memory settings
  Int coverage_mem_gb = 24
  Int coverage_disk_gb = 100
  Int samtools_mem_gb = 24
  Int samtools_disk_gb = 100

  ## set docker and software path options
  String cnvkit_path = "/root/anaconda3/bin/cnvkit.py" # absolute path to cnvkit.py in docker container
  String samtools_path = "/usr/bin/samtools" # absolute path to samtools in docker containter
  String sed_path = "/bin/sed" # absolute path to sed in docker container
  String cnvkit_docker = "weber813/cnvkit:latest" # docker image with CNVkit

  # define variables related to sample BAMs and their ids
  Array[String] sample_ids
  Array[File] bamfiles
  Array[Pair[String,File]] bam_map = zip(sample_ids, bamfiles)


  #------------------------------------
  # run get coverages pipeline calls
  ## generate annotated targets file (1 per PON)
  call cnvkit_target {
    input: cnvkit_docker=cnvkit_docker,
            cnvkit_path=cnvkit_path,
            baits=baitset,
            targets="wxs_targets.bed",
            refgene=refgene
  }

  ## generate access file for anti targets (1 per PON)
  call cnvkit_access {
    input: cnvkit_docker=cnvkit_docker,
            cnvkit_path=cnvkit_path,
            access_mappable="mappable-regions.bed",
            assembly_fa=assembly_fa,
            gaps=gaps
  }

  ## generate annotated anti-targets file (1 per PON)
  call cnvkit_antitarget {
    input: cnvkit_docker=cnvkit_docker,
            cnvkit_path=cnvkit_path,
            targets=cnvkit_target.targets_out,
            mappable_regions=cnvkit_access.mappable_regions_out,
            antitargets="wxs_anti-targets.bed"
  }

  ## parallelize the optional fix_chrs if set to true
  if (fix_chrs) {
    scatter (b in bam_map) {
      call samtools_fix_chrs {
        input: cnvkit_docker=cnvkit_docker,
                samtools_path=samtools_path,
                sed_path=sed_path,
                samtools_mem_gb=samtools_mem_gb,
                samtools_disk_gb=samtools_disk_gb,
                unfmtted_bam=b.right,
                fmttd_bam=b.left + "-fmttd.bam"
      }
    }
  }

  ## resolve branched execution
  Array[File] input_bams = select_first([samtools_fix_chrs.fmttd_bam_out, bamfiles])
  Array[Pair[String,File]] input_bam_map = zip(sample_ids, input_bams)

  ## run coverage and anticoverage on formatted bams
  scatter(s in input_bam_map) {
      ### call the coverage on targets
      call cnvkit_target_coverage {
        input: cnvkit_docker=cnvkit_docker,
                cnvkit_path=cnvkit_path,
                bamfile=s.right,
                targets=cnvkit_target.targets_out,
                target_coverage="${s.left}.targetcoverage.cnn",
                coverage_mem_gb=coverage_mem_gb,
                coverage_disk_gb=coverage_disk_gb
      }

      ### call the coverate on
      call cnvkit_anti_coverage {
        input: cnvkit_docker=cnvkit_docker,
                cnvkit_path=cnvkit_path,
                bamfile=s.right,
                antitargets=cnvkit_antitarget.antitargets_out,
                anti_coverage="${s.left}.antitargetcoverage.cnn",
                coverage_mem_gb=coverage_mem_gb,
                coverage_disk_gb=coverage_disk_gb
      }
    }



  #------------------------------------
  # specify workflow outputs. in this case, the coverage files for
  # anti targets and targets
  output {
    Array[File] prepd_bams = input_bams
    Array[File] target_coverages = cnvkit_target_coverage.coverage_out
    Array[File] antitarget_coverages = cnvkit_anti_coverage.coverage_out
  }
}


####################################################
#####   define the cnvkit_get_coverages tasks  #####
####################################################
#--------------------------------------
# run cnvkit targets to map to bin bait regions for read count mapping
task cnvkit_target {
  # input variables
  File baits # bait regions for WXS (.bed)
  File refgene # refgene annotation file for appropriate genome (.txt)
  String targets # output target file name (.bed)
  String cnvkit_docker # docker image, see runtime comments for specs
  String cnvkit_path

  # run cnvkit targets command to generate an annotated
  # bed file of bins with exon level annotation
  command {
    ${cnvkit_path} target ${baits} \
              --annotate ${refgene} \
              --split \
              -o ${targets}

  }

  # specify docker environment
  # here we use my cnvkit docker posted to
  # weber813/cnvkit. It can be any unix based
  # docker image with cnvkit installed and callable
  # on the $PATH
  runtime {
    docker: cnvkit_docker
  }

  # specify outputs
  output {
    # capture the output file and return it
    File targets_out = targets
  }
}

#--------------------------------------
# run cnvkit.py acess to generate bed file for accessible
# regions of the genome
task cnvkit_access {
  # declare variables
  String access_mappable # output file name for task (.bed)
  File assembly_fa #  genome assembly fasta (.bed)
  File gaps # gaps in genome assembly (.bed)
  String cnvkit_docker # symbolic link to docker
  String cnvkit_path

  # run cnvkit access to generate a map of accessible regions
  # for later generation of anti-targets
  command {
    ${cnvkit_path} access ${assembly_fa} \
              --exclude ${gaps} \
              -o ${access_mappable}
  }

  # specify runtime docker environment
  # here we use my cnvkit docker posted to
  # weber813/cnvkit. It can be any unix based
  # docker image with cnvkit installed and callable
  # on the $PATH
  runtime {
    docker: cnvkit_docker
  }

  # specify outputs, which here are the exclusion sites
  output {
    # capture exclusions and output the file
    File mappable_regions_out = access_mappable
  }
}

#--------------------------------------
# run cnvkit antitargets to bin non-bait regions of the genome
task cnvkit_antitarget {
  # input variables
  File targets # targets file output from targets task (.bed)
  File mappable_regions # mappable file output from access task (.bed)
  String antitargets # name of the antitargets file to be created (.bed)
  String cnvkit_docker # symbolic link to docker
  String cnvkit_path # absolute path to cnvkit exec inside docker image

  # run cnvkit antitarget comand to generate bed
  # file which is the gaps of the targets bed file
  command {
    ${cnvkit_path} antitarget ${targets} \
              --access ${mappable_regions} \
              -o ${antitargets}
  }

  # specify docker runtime environment
  # see task targets runtime for more info
  runtime {
    docker: cnvkit_docker
  }

  # specify outputs, which here are the anti-target sites
  output {
    File antitargets_out = antitargets
  }
}

#--------------------------------------
# run samtools reheader task to change bams from 1-22+X indexed to chr1-chr22
# + chrX indexed
task samtools_fix_chrs {
  #input variables
  File unfmtted_bam
  Int samtools_mem_gb
  Int samtools_disk_gb
  String fmttd_bam
  String cnvkit_docker
  String samtools_path
  String sed_path

  # run the piped command using samtools and sed
  command {
      # run samtools sed pipe to reheader the bam file
      ${samtools_path} view -H ${unfmtted_bam} | \
        ${sed_path}  -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | \
        ${samtools_path} reheader - ${unfmtted_bam} > ${fmttd_bam}
  }

  # specify runtime environment. see cnvkit_target runtime comments for
  # further specification on docker container
  runtime {
    docker: cnvkit_docker
    memory: samtools_mem_gb + "G"
    disks: "local-disk " + samtools_disk_gb + " HDD"
    bootDiskSizeGb: 50
  }

  # specify outputs, here a reheadered bam file
  output {
    File fmttd_bam_out = fmttd_bam
  }

}

#--------------------------------------
# run cnvkit coverage task on targets to build coverage file for each normal
task cnvkit_target_coverage {
  # input variables
  File bamfile # input bam file for a normal (parallelized in workflow)
  File targets # input targets generated from the targets task
  String target_coverage # output file name of target coverage file (.cnn)
  String cnvkit_docker # symbolic link to docker image
  String cnvkit_path # absolute path to cnvkit exec inside docker image
  Int coverage_mem_gb # required memory for the coverage task
  Int coverage_disk_gb # required disk space for the coverage task

  # run cnvkit_coverages command
  command {
    ${cnvkit_path} coverage --count \
              ${bamfile} \
              ${targets} \
              -o ${target_coverage}
  }

  # specify docker runtime environment
  # see task targets for more runtime info
  runtime {
    docker: cnvkit_docker
    memory: coverage_mem_gb + "G"
    disks: "local-disk " + coverage_disk_gb + " HDD"
    bootDiskSizeGb: 50
  }

  # specify outputs file , here coverages (.cnn) for targets
  output {
    File coverage_out = target_coverage
  }
}

#--------------------------------------
# run cnvkit coverage task on antitargets to build coverage file for each normal
task cnvkit_anti_coverage {
  # input variables
  File bamfile # input bam file for a normal (parallelized in workflow)
  File antitargets # input antitargets generated from the antitargets task
  String anti_coverage # output file name of antitarget coverage file (.cnn)
  String cnvkit_docker # symbolic link to docker image
  String cnvkit_path # absolute path to cnvkit exec inside docker image
  Int coverage_mem_gb # required memory for the coverage task
  Int coverage_disk_gb # required disk space for the coverage task

  # run cnvkit_coverages command
  command {
    ${cnvkit_path} coverage --count \
              ${bamfile} \
              ${antitargets} \
              -o ${anti_coverage}
  }

  # specify docker runtime environment
  # see task targets for more runtime info
  runtime {
    docker: cnvkit_docker
    memory: coverage_mem_gb + "G"
    disks: "local-disk " + coverage_disk_gb + " HDD"
    bootDiskSizeGb: 50
  }

  # specify outputs file , here coverages (.cnn) for antitargets
  output {
    File coverage_out = anti_coverage
  }
}
