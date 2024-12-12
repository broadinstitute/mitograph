use rust_htslib::bam::{Record, record::CigarStringView};

#[derive(Debug)]
pub struct Variant {
    pub pos: i64,        // Position on reference
    pub ref_allele: String,  // Reference allele
    pub alt_allele: String,  // Alternative allele
    pub variant_type: String, // Type of variant (SNP, INS, DEL)
}

pub fn extract_variants(record: &Record, ref_seq: &[u8]) -> Vec<Variant> {
    let mut variants = Vec::new();
    let cigar = record.cigar();
    let read_seq = record.seq().as_bytes();
    let start_pos = record.pos();
    
    let mut read_pos: usize = 0;
    let mut ref_pos: i64 = start_pos;
    
    for cigar_elem in cigar.iter() {
        match cigar_elem {
            // Match or mismatch
            rust_htslib::bam::record::Cigar::Match(len) |
            rust_htslib::bam::record::Cigar::Equal(len) |
            rust_htslib::bam::record::Cigar::Diff(len) => {
                // Check for mismatches within the matched region
                for i in 0..*len as usize {
                    let ref_base = ref_seq[(ref_pos + i as i64) as usize];
                    let read_base = read_seq[read_pos + i];
                    
                    if ref_base != read_base {
                        variants.push(Variant {
                            pos: ref_pos + i as i64,
                            ref_allele: String::from_utf8(vec![ref_base]).unwrap(),
                            alt_allele: String::from_utf8(vec![read_base]).unwrap(),
                            variant_type: "SNP".to_string(),
                        });
                    }
                }
                read_pos += *len as usize;
                ref_pos += *len as i64;
            },
            
            // Insertion
            rust_htslib::bam::record::Cigar::Ins(len) => {
                let inserted_seq = String::from_utf8(
                    read_seq[read_pos..read_pos + *len as usize].to_vec()
                ).unwrap();
                
                variants.push(Variant {
                    pos: ref_pos,
                    ref_allele: "-".to_string(),
                    alt_allele: inserted_seq,
                    variant_type: "INS".to_string(),
                });
                
                read_pos += *len as usize;
            },
            
            // Deletion
            rust_htslib::bam::record::Cigar::Del(len) => {
                let deleted_seq = String::from_utf8(
                    ref_seq[ref_pos as usize..(ref_pos + *len as i64) as usize].to_vec()
                ).unwrap();
                
                variants.push(Variant {
                    pos: ref_pos,
                    ref_allele: deleted_seq,
                    alt_allele: "-".to_string(),
                    variant_type: "DEL".to_string(),
                });
                
                ref_pos += *len as i64;
            },
            
            // Skip these operations
            rust_htslib::bam::record::Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            },
            rust_htslib::bam::record::Cigar::HardClip(_) => {},
            rust_htslib::bam::record::Cigar::Pad(_) => {},
            rust_htslib::bam::record::Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            },
        }
    }
    
    variants
}

// Example usage:
fn main() {
    // Open BAM file and get reference sequence
    let bam = rust_htslib::bam::Reader::from_path("input.bam").unwrap();
    let ref_seq = // ... load your reference sequence ...;
    
    for record in bam.records() {
        let record = record.unwrap();
        let variants = extract_variants(&record, &ref_seq);
        
        for variant in variants {
            println!("Position: {}, Ref: {}, Alt: {}, Type: {}", 
                variant.pos, variant.ref_allele, variant.alt_allele, variant.variant_type);
        }
    }
}