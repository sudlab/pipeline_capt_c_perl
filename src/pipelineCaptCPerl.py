import pysam
import re
import CGAT.Bed as Bed
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
from CGATPipelines.Pipeline import cluster_runnable

# Inspired on fetchProbeFragments from Capture C
@cluster_runnable
def getProbeFragments(probe_bed, digest_bed, outfile,
                        lookup_out):
    
    # First find the length of the restriction enzyme cut, required to obtain the start and end coordinates
    # from the pregenerated file.
    # First iteration, no comparison
    first_iteration = True
    
    length_RE_cut = 0
    
    last_bed = None
    
    for bed_digest in Bed.iterator(IOTools.openFile(digest_bed)):
                
        if(first_iteration):
            first_iteration = False
        else:
            # If they are in the same contig they can be compared
            if(bed_digest.contig == last_bed.contig):
                length_RE_cut = bed_digest.start - last_bed.end
                break
        
        last_bed = bed_digest
    
    
    digest_fragments = pysam.TabixFile(digest_bed)
    bed = Bed.Bed()
    with IOTools.openFile(outfile, "w") as outf, \
         IOTools.openFile(lookup_out,"w") as lookup:

        lookup.write("probe\tfragment\n")
        for probe in Bed.iterator(IOTools.openFile(probe_bed)):
            
            frag = digest_fragments.fetch(probe.contig,
                                          probe.start,
                                          probe.end,
                                          parser=pysam.asBed())
            frag = list(frag)
            if not len(frag) == 1:
                E.warn("%i fragments found for probe %s, skipping" %
                       (len(frag), probe.name))
                continue

            frag = frag[0]
            
            # The restriction enzyme cut on the left side of the fragment
            # is the end site of the last restriction enzyme fragment + 1
            # (+1 because according to the manual coordinates are specified
            # in 1-origin for the bed start.)
            
            bed.start = frag.start-length_RE_cut+1
            bed.end = frag.end+length_RE_cut
            bed.contig = frag.contig
            bed["name"] = probe.name
            bed["score"] = "."
            bed["strand"] = "+"

            lookup.write("%s\t%s\n" % (probe.name, frag.name))
            outf.write(str(bed) + "\n")


@cluster_runnable
# The name column of the bed file should avoid any character which
# is not alphanumerical or contains underscore (_)
def formatNameCompliance(name):
    
    out_string = ""
    
    for character in name:
        if(character.isalnum() or character == '_'):
            out_string += character
    
    return out_string  
    

@cluster_runnable
# Creates 9 column format and calculates collisions in probe ranges and
# exclusion fragments returning strings indicating collisions
def formatProbeFragments(probe_fragments_bed, outfile):
    
    # Create strings for output collisions
    probe_collisions = ""
    exclus_collisions = ""
    
    with IOTools.openFile(outfile, "w") as outf:
    
        # dictionaries for collisions per chr
        chr_probe = {}
        chr_exclus = {}
    
        for bed_digest in Bed.iterator(IOTools.openFile(probe_fragments_bed)):
            
            # chromosome needs to specify only number (remove chr)
            chromosome = re.sub('chr', '', bed_digest.contig)
            
            out_array = []
            
            out_array.append(formatNameCompliance(bed_digest["name"]))
            out_array.append(str(chromosome))
            out_array.append(str(bed_digest.start))
            out_array.append(str(bed_digest.end))
            out_array.append(str(chromosome))
            out_array.append(str(bed_digest.start-1000))
            out_array.append(str(bed_digest.end+1000))
            out_array.append("1")
            out_array.append("A")
            
            outf.write("\t".join(out_array) + "\n")
            
            
            # Calculate collisions, all coordinates assumed in BED format
                        
            # It needs to be done per chr
            # chr_probe -> probe_ranges -> probe_range
            # chr_exclus -> exclus_ranges -> exclus_range
            
            # Check if chromosome key already exists in one of the
            # dictionaries 
            # (ranges introduced before on that chromosome
            
            if bed_digest.contig in chr_probe:
                # If it exists, get the array of probe ranges 
                # and exclusion ranges already stored
                probe_ranges = chr_probe[bed_digest.contig]
                exclus_ranges = chr_exclus[bed_digest.contig]     
            
            
            # If it doesn't create the arrays
            else:

                probe_ranges = []
                exclus_ranges = []
            
                
            # Create a range for probes and exclusion fragments
            probe_range = []
            exclus_range = []
            
            probe_range.append(bed_digest.start)
            probe_range.append(bed_digest.end)
            
            probe_ranges.append(probe_range)
            
            
            exclus_range.append(bed_digest.start-1000)
            exclus_range.append(bed_digest.end+1000)
            
            exclus_ranges.append(exclus_range)
            
            
            # Substitute the ranges back to the corresponding chr in the
            # dictionary
            chr_probe[bed_digest.contig] = probe_ranges
            chr_exclus[bed_digest.contig] = exclus_ranges
        
        
        for chr in chr_probe:
            probe_ranges = chr_probe[chr]
            exclus_ranges = chr_exclus[chr]
            probe_intersection = set.intersection(*(set(range(start, finish)) for start, finish in probe_ranges))
            exclus_intersection = set.intersection(*(set(range(start, finish)) for start, finish in exclus_ranges))
            
            if(len(probe_intersection) != 0):
                probe_collisions += "Probe collision " +str(chr)
                probe_collisions += " " +str(min(probe_intersection))
                probe_collisions += " " +str(max(probe_intersection))
                probe_collisions += "\n"
            
            if(len(exclus_intersection) != 0):
                exclus_collisions += "Exclusion collision " +str(chr)
                exclus_collisions += " " +str(min(exclus_intersection))
                exclus_collisions += " " +str(max(exclus_intersection))
                exclus_collisions += "\n"

    return (probe_collisions, exclus_collisions)
    
