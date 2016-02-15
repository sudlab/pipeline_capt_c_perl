##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_capt_c_perl.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import pipelineCaptCPerl
import numpy
import sys
sys.path.insert(0, '/home/mbp15ja/dev/Capture_C/Capture_C')
import PipelineCaptureC



# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))



# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.



#########################################################################
#########################################################################
#########################################################################
# Read mapping
#########################################################################

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(".", suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|fastq.gz|fa.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz)")

###################################################################
###################################################################





# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


###############################################################################
@active_if(PARAMS["addtests_saturation"] == 1)
@follows(mkdir("saturation_analysis.dir"))
@originate( ["saturation_analysis.dir/%s.sample_sat_analysis" % percentage for percentage in
            numpy.arange(PARAMS["saturationanlysis_min"],
                         PARAMS["saturationanlysis_max"] + 1,
                         PARAMS["saturationanlysis_step"])])                         
def generateSaturationAnalysis(outfile):
    
    statement = ''' touch %(outfile)s '''
    
    P.run()




###############################################################################
# Generate all the combinations of reads and sample sizes
@active_if(PARAMS["addtests_saturation"] == 1)
@product(SEQUENCEFILES,
          SEQUENCEFILES_REGEX,
          generateSaturationAnalysis,
          formatter("(.sample_sat_analysis)$"),
          r"saturation_analysis.dir/Sample_{basename[1][0]}_{basename[0][0]}.gz",
          "{basename[1][0]}")
def generateReadSamplesProduct(infile, outfile, sample_size):
    
    percentage_sample = float(float(int(sample_size)) / 100)
    
    
    read1_in = infile[0]
    read2_in = P.snip(infile[0], ".fastq.1.gz") + ".fastq.2.gz"
    
    outnamebasef1 = outfile
    outnamebasef2 = P.snip(outfile, ".fastq.1.gz") + ".fastq.2.gz"
    
    log_filename = P.snip(outfile, ".fastq.1.gz") + ".log"    
    
    statement = ''' python %(scriptsdir)s/fastq2fastq.py --stdin=%(read1_in)s  
                --method=sample --sample-size=%(percentage_sample)s --seed=1234 
                --log=%(log_filename)s --pair-fastq-file=%(read2_in)s 
                -F --output-filename-pattern=%(outnamebasef2)s -v 7 -S %(outnamebasef1)s '''
     
    P.run()
    
    



###############################################################################
# If no saturation test defined, just copy the symlinks to saturation_analysis.dir as 100 sample
@active_if(PARAMS["addtests_saturation"] == 0)
@follows(mkdir("saturation_analysis.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"saturation_analysis.dir/Sample_100_\1.fastq.1.gz")
def relocateReads(infiles, outfile):
    
    
    
    # Move both reads of the pair
    read2_in = P.snip(infiles, ".fastq.1.gz") + ".fastq.2.gz"
    
    read2_out = P.snip(outfile, ".fastq.1.gz") + ".fastq.2.gz"
    
    # Copy first the second read pair and checkpoint to return the exit code.
    statement = '''cp -d %(read2_in)s %(read2_out)s;
                checkpoint;
                cp -d %(infiles)s %(outfile)s;
                '''
    
    P.run()



# Creates overlap sizes according to the ini file, if no flash analysis
# is specified, creates a default overlap size of 10 
###############################################################################
@follows(mkdir("flashing_analysis.dir"))
@originate( ["flashing_analysis.dir/%s.sample_flash_analysis" % overlap for overlap in
            numpy.arange(PARAMS["flashinganlysis_min"],
                         PARAMS["flashinganlysis_max"] + 1,
                         PARAMS["flashinganlysis_step"])]
           if (PARAMS["addtests_flashing"] == 1)
           else "flashing_analysis.dir/10.sample_flash_analysis")                         
def generateFlashingAnalysis(outfile):
     
    statement = ''' touch %(outfile)s '''
     
    P.run()     
   
   
   
  


formatter(".+/job(?P<JOBNUMBER>\d+).a.start",  # Extract job number
                      ".+/job[123].b.start") 

# Generate all the combinations of input reads and minimum overlap lengths
###############################################################################
@follows(mkdir("flashed.dir"))
@product([generateReadSamplesProduct, relocateReads],
         formatter(".+/(?P<CELL_LINE>.+).fastq.1.gz"),
         generateFlashingAnalysis,
         formatter(".+/(?P<OVERLAP_FLASH>.+).sample_flash_analysis"),
         "flashed.dir/Overlap_{OVERLAP_FLASH[1][0]}_{CELL_LINE[0][0]}.extendedFrags.fastq.gz",
         "{OVERLAP_FLASH[1][0]}")
def flashReads(infiles, outfile, overlap):
    ''' Flashes read pairs'''
    
    job_memory = "2G"
    
    # Use one thread to maintain read order
    job_threads = 1
    
    # Retrieve the input file from [generateReadSamplesProduct, relocateReads]
    read1 = infiles[0]
    
    
    read2 = P.snip(read1, ".fastq.1.gz") + ".fastq.2.gz"
    
    
    outfilebase = P.snip(outfile, ".extendedFrags.fastq.gz")   
    
    
    statement = '''flash -m %(overlap)s --interleaved-output %(read1)s 
    %(read2)s -o %(outfilebase)s -z ;
    checkpoint;
    touch %(outfile)s'''
    # touch in case there are no flashed reads 
    # but the process still succeeds with no errors

    P.run()


###############################################################################
@transform(flashReads,
           suffix(".extendedFrags.fastq.gz"),
           r"\1_combined.fastq")
def combineFlashedReads(infiles, outfile):
    
    flashed = infiles
    unflashed = P.snip(infiles, ".extendedFrags.fastq.gz") + ".notCombined.fastq.gz"
    
    # Use one thread to maintain read order
    job_threads = 1
             
    statement = ''' zcat %(unflashed)s %(flashed)s > 
                    %(outfile)s;
                    checkpoint;''' 
     
    P.run()




@follows(mkdir("perl_digest_reads.dir"))
@transform(combineFlashedReads,
           regex(".+/(.+)_combined.fastq"),
           r"perl_digest_reads.dir/\1_combined_REdig.fastq")
def digestFlashedReads(infiles, outfile):
    
    perl = PARAMS["environment_perl"]
    perl_scripts = os.path.join(PARAMS["environment_perlscripts"],
                                "dpnII2E.pl")
    
    # The file is initially stored in the flashed.dir directory

    temp_in_file_dir = os.path.dirname(infiles)
    temp_in_file_name = os.path.basename(outfile)
    
    temp_in_file_full_path = os.path.join(temp_in_file_dir, temp_in_file_name)
    
    statement = ''' %(perl)s %(perl_scripts)s %(infiles)s ;
                    mv %(temp_in_file_full_path)s %(outfile)s
                '''

    P.run()



###############################################################################

# !!! Not adapted for different enzymes, only DPNII
@follows(mkdir("perl_digest_genome.dir"))
@files(os.path.join(PARAMS["genome_dir"], "%s.fa" % PARAMS["genome"]),
           r""+"".join((os.path.join("perl_digest_genome.dir/", ("%s" % PARAMS["genome"]) + "_dpnII_coordinates.txt"))))
def perlDigestGenome(infiles, outfile):
    
    job_memory = "4G"
    
    # The file is initially stored in the base directory
    temp_out_file = os.path.basename(outfile)
    
    perl = PARAMS["environment_perl"]
    perl_scripts = os.path.join(PARAMS["environment_perlscripts"], 
                                "dpngenome3_1.pl")
    

    # Moves the file to the perl_digest_genome.dir
    statement = ''' %(perl)s %(perl_scripts)s %(infiles)s ;
                    mv %(temp_out_file)s %(outfile)s
                '''
    
    
    P.run()
    


##############################################################################
@follows(mkdir("bowtie.dir"))
@transform(digestFlashedReads,
           regex(".+/(.+)_combined_REdig.fastq"),
           r"bowtie.dir/\1.bowtie.sam")
def mapReadsWithBowtie(infiles, outfile):
    ''' Aligns the digested reads to the genome '''
    
    genome = PARAMS["environment_bowtiegenome"]
    
    align_command = PARAMS["bowtie_options"]
    
    job_memory = "4G"
        
    log_file = P.snip(outfile, ".sam") + ".log"
 
    statement = '''bowtie %(align_command)s 
    --sam %(genome)s 
    %(infiles)s 
    %(outfile)s 2> %(log_file)s'''
 
    P.run()
    
    
       
    


###############################################################################
# Digest genome
###############################################################################
@follows(mkdir("digest.dir"))
@split(os.path.join(PARAMS["genome_dir"], "%s.fa" % PARAMS["genome"]),
       ["digest.dir/chr%s.bed.gz" % chrom for chrom in
        map(str, range(1, 23)) + ["X", "Y"]])
def splitDigest(infile, outfiles):
    '''The entire genome consumes too much memeory for EMBOSS restrict. Thus
    the genome must be digested one chormosome at a time'''

    chroms = [os.path.basename(P.snip(outfile, ".bed.gz")) for outfile in outfiles]

    statements = []

    job_memory = "3G"

    for chrom in chroms:

        enzyme = PARAMS['experiment_enzyme']
        
        statements.append('''
                          restrict
                           -sequence <(  cat %%(infile)s
                                       | python %%(scriptsdir)s/fasta2fasta.py
                                          --include="%(chrom)s$"
                                         -L digest.dir/%(chrom)s.log)
                           -enzymes %(enzyme)s
                           -sitelen 4
                           -outfile digest.dir/%(chrom)s.gff
                           -rformat gff > digest.dir/%(chrom)s.log;

                           checkpoint;

                           python %%(scriptsdir)s/gff2bed.py
                                  -I digest.dir/%(chrom)s.gff
                                  -L digest.dir/%(chrom)s.log
                                  --set-name=source
                         | awk 'BEGIN {site = 0;OFS = "\\t"}
                                      {site++; $4=site;print}'
                         | bgzip > digest.dir/%(chrom)s.bed.gz;

                           checkpoint;

                           rm digest.dir/%(chrom)s.gff''' % locals())

    P.run()





###############################################################################
@merge(splitDigest, "digest.dir/digest.bed.gz")
def mergeDigest(infiles, outfile):

    infile = " ".join(infiles)
    statement = '''zcat %(infile)s | sort -k1,1 -k2,2n | bgzip > %(outfile)s;
                   checkpoint;
                   tabix -p bed %(outfile)s'''
    P.run()



###############################################################################
@transform(mergeDigest,
           suffix("digest.bed.gz"),
           "fragments.bed.gz")
def digest2fragments(infile, outfile):

    genome = PARAMS["annotations_interface_contigs_tsv"]
    outfile = P.snip(outfile, ".gz")
    PipelineCaptureC.sites2fragments(infile, genome, outfile)
    statement = '''uniq %(outfile)s > bgzip %(outfile)s.gz;
                   checkpoint;
                   tabix -p bed %(outfile)s.gz;
                   checkpoint;
                   rm %(outfile)s'''




# ###############################################################################
# # Generate the coordinates as specified in the Perl Capture C Manual
# ###############################################################################
@transform(digest2fragments,
           suffix("fragments.bed.gz"),
           add_inputs("probes.bed.gz"),
           "probe_fragments.bed.gz")    
def generateCoordinates(infile, outfile):
    fragments, probes = infile
    lookup = P.snip(outfile, "bed.gz") + "lookup.tsv"
    pipelineCaptCPerl.getProbeFragments(probes,
                                        fragments,
                                        outfile, 
                                        lookup)
    

# ###############################################################################
# # Deduplicate restriction fragments
# ###############################################################################
@transform(generateCoordinates,
           suffix("probe_fragments.bed.gz"),
           "dedup_probe_fragments.bed.gz")    
def deduplicateFragments(infile, outfile):
        
    # each DpnII fragment only ONCE in the oligo coordinate file
    
    job_memory = "0.5G"
    
    tab_separator = "'\\t'"
    
    statement = ''' zcat %(infile)s | 
                    sort --field-separator=$%(tab_separator)s
                    --key=1,1 --key=2,2 --key=3,3 --key=5,5 
                    --key=6,6 --unique | 
                    bgzip > %(outfile)s;
                    checkpoint;
                    '''
                    

    P.run()


# ###############################################################################
# Return 9 column format
# Print any collisions
# ###############################################################################
@transform(deduplicateFragments,
           suffix("dedup_probe_fragments.bed.gz"),
           "oligo_coord_dedup.txt")    
def formatProbeCoordinates(infile, outfile):
    '''Format the restriction fragments containing probes to include 1000bp exclusion sites 
    on each side of each probe and report any collisions in probe fragments and exclusions.
    Collisions are reported for a sanity control check: probe collisions shouldn't occur and
    exclusion collisions can happen but they are not wrong, nor do they affect the analysis'''
    
    job_memory = "0.5G"
    
    (probe_collisions, exclus_collisions) = pipelineCaptCPerl.formatProbeFragments(infile, outfile)
    
    E.debug(probe_collisions)
    
    E.debug(exclus_collisions)





@files(PARAMS["environment_contigs"],
       os.path.join(PARAMS["environment_tempcontigs"], (PARAMS["genome"]+"_sizes.txt")))
def linkContigs(infile, outfile):
    '''The perl analyzing script requires the contigs file to be called a
    certain way and be present in a hardcoded specified directory.
    So it doesn't interfere with other runs,
    only copy the file from /shared/sudlab1/General/annotations/hg19_ensembl75/contigs.tsv
    if it is not present already'''
    
     
    statement = ''' if [ ! -f %(outfile)s ]; then
                        cp %(infile)s %(outfile)s;
                    fi
                    '''
    P.run()


@follows(mkdir("analysis.dir"), perlDigestGenome, linkContigs)
@transform(mapReadsWithBowtie,
           regex(".+/(.+).bowtie.sam"),
           add_inputs([formatProbeCoordinates, r""+"".join((os.path.join("perl_digest_genome.dir/", ("%s" % PARAMS["genome"]) + "_dpnII_coordinates.txt")))]),
           r"analysis.dir/\1_experiment")
def analyzeInteractions(infile, outfile):
    
    (mapped_RE_reads, [probe_coordinates, digest_genome]) = infile 
    
    job_memory = "8G"
    
    perl = PARAMS["environment_perl"]
    perl_scripts = os.path.join(PARAMS["environment_perlscripts"], 
                                "CCanalyser3.pl")
    
    perl_script_base_name = os.path.basename(perl_scripts)
    
    genome = PARAMS["genome"]
    
          
    
    out_dir = os.path.dirname(outfile)
    
    # Getting full directory where it is being executed
    dir_execution = os.getcwd()
    
    
    full_perl_script = os.path.join(dir_execution, perl_script_base_name)
    
    full_mapped_RE_reads = os.path.join(dir_execution, mapped_RE_reads)
    
    full_digest_genome = os.path.join(dir_execution, digest_genome)
    
    full_probe_coordinates = os.path.join(dir_execution, probe_coordinates)
    
    full_dir_execution = os.path.join(dir_execution, out_dir)
    
    # The analysis.dir for the output file is gotten from another parameter
    outfile_base = os.path.basename(outfile)

    
    # Copy the script to the base directory first
    # Mv .bw to analysis.dir
    statement = ''' if [ ! -f %(perl_script_base_name)s ]; then
                        cp %(perl_scripts)s .;
                    fi;
                    %(perl)s %(perl_script_base_name)s -f 
                    %(full_mapped_RE_reads)s 
                    -r %(full_digest_genome)s
                    --genome %(genome)s
                    -o %(full_probe_coordinates)s
                    -s %(outfile_base)s
                    --pf %(full_dir_execution)s
                    --pu %(out_dir)s;
                    '''
    P.run()
    
    
    
# ---------------------------------------------------
# Generic pipeline tasks
@follows(analyzeInteractions)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
