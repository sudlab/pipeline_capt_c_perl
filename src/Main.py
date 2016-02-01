# import PipelineCaptureC as Capt_C
# import CGAT.IOTools as IOTools
# Import for full make (local directory changed to data files in debug opts)
#import sys
#sys.path.insert(0, '/home/mbp15ja/dev/Capture_C/Capture_C')

#import pipeline_capt_c_perl as Pipe_Capt_C_Perl
import pipelineCaptCPerl


if __name__ == "__main__":
    
    
    probe_collisions, exclus_collisions = pipelineCaptCPerl.formatProbeFragments('/mnt/fastdata/mbp15ja/capturec-pilot/capture_c_pipeline_perl/digest.dir/dedup_probe_fragments.bed.gz',
                                           '/mnt/fastdata/mbp15ja/capturec-pilot/oligo_coord_dedup.txt')
    
    print(probe_collisions)
    print(exclus_collisions)
#     bam = '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2/dedupped.dir/IGF0003828.bwa.bam'
#     probes = '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2/probe_fragments.bed.gz'
#     outfile = '/fastdata/mbp15ja/capturec-pilot/IGF0003828.bwa.anomolies.tsv.gz'
#     Capt_C.quantifyAnomolies(bam, probes, outfile)
#     for line in IOTools.openFile('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2/probe_fragments.bed.gz'):
#         print(line.split("\t")[3])
     
#     Capt_C.fetchProbeFragments('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2/probes.bed.gz', 
#                                          '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2/digest.dir/fragments.bed.gz', 
#                                          '/fastdata/mbp15ja/capturec-pilot/probe_fragments.bed.gz', 
#                                          '/fastdata/mbp15ja/capturec-pilot/lookup_out')
#     
#     Capt_C.sites2fragments('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/digest.dir/digest.bed.gz',
#                            '/shared/sudlab1/General/annotations/hg19_ensembl75/contigs.tsv',
#                            '/fastdata/mbp15ja/capturec-pilot/out_fragments.bed.gz')
#     Capt_C.countInteractions('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/dedupped.dir/IGF0003828.bwa.by_name.bam', 
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/digest.dir/fragments.bed.gz',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/probe_fragments.bed.gz',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/test_output1.txt',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/test_output2.txt')
#     
