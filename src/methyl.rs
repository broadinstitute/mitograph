use ndarray::Array2;


use rust_htslib::bam::{Record, Reader};
use std::path::Path;

/// Get aligned pairs similar to pysam's get_aligned_pairs function
/// 
/// Returns a vector of (read_pos, ref_pos) tuples where:
/// - read_pos is None for deletions/skipped regions (ref-only positions)
/// - ref_pos is None for insertions/soft-clips (read-only positions)
fn get_aligned_pairs(record: &Record, with_seq: bool) -> Vec<(Option<usize>, Option<i64>, Option<u8>)> {
    let cigar = record.cigar();
    let mut read_pos: usize = 0;
    let mut ref_pos: i64 = record.pos();
    let mut alignments = Vec::new();
    
    // Get reference sequence if requested
    let ref_seq = if with_seq {
        Some(record.seq().as_bytes())
    } else {
        None
    };
    
    for cigar_op in cigar.iter() {
        let (op_len, op) = cigar_op.to_tuple();
        let op_len = op_len as usize;
        
        match op {
            0 | 7 | 8 => { // Match, equal, diff (M, =, X)
                // Both sequences advance
                for i in 0..op_len {
                    let ref_base = if with_seq && read_pos + i < record.seq_len() as usize {
                        Some(ref_seq.as_ref().unwrap()[read_pos + i])
                    } else {
                        None
                    };
                    
                    alignments.push((Some(read_pos + i), Some(ref_pos + i as i64), ref_base));
                }
                read_pos += op_len;
                ref_pos += op_len as i64;
            },
            1 => { // Insertion (I)
                // Only read advances
                for i in 0..op_len {
                    let read_base = if with_seq && read_pos + i < record.seq_len() as usize {
                        Some(ref_seq.as_ref().unwrap()[read_pos + i])
                    } else {
                        None
                    };
                    
                    alignments.push((Some(read_pos + i), None, read_base));
                }
                read_pos += op_len;
            },
            2 => { // Deletion (D)
                // Only reference advances
                for i in 0..op_len {
                    alignments.push((None, Some(ref_pos + i as i64), None));
                }
                ref_pos += op_len as i64;
            },
            3 => { // Skipped region/intron (N)
                // Only reference advances
                for i in 0..op_len {
                    alignments.push((None, Some(ref_pos + i as i64), None));
                }
                ref_pos += op_len as i64;
            },
            4 => { // Soft clip (S)
                // Only read advances
                for i in 0..op_len {
                    let read_base = if with_seq && read_pos + i < record.seq_len() as usize {
                        Some(ref_seq.as_ref().unwrap()[read_pos + i])
                    } else {
                        None
                    };
                    
                    alignments.push((Some(read_pos + i), None, read_base));
                }
                read_pos += op_len;
            },
            5 => { // Hard clip (H)
                // Neither read nor reference is advanced
                // Hard-clipped bases are not in the read sequence
            },
            6 => { // Padding (P)
                // Neither read nor reference is advanced
            },
            _ => {
                // Unknown operation
                eprintln!("Warning: Unknown CIGAR operation {}", op);
            }
        }
    }
    
    alignments
}

/// Simplified version without sequence information
fn get_aligned_pairs_simple(record: &Record) -> Vec<(Option<usize>, Option<i64>)> {
    let with_seq_result = get_aligned_pairs(record, false);
    with_seq_result.into_iter()
                   .map(|(read_pos, ref_pos, _)| (read_pos, ref_pos))
                   .collect()
}

fn start() -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new("your_file.bam");
    let mut reader = Reader::from_path(path)?;
    
    for record in reader.records() {
        let record = record?;
        
        // Get aligned pairs without sequence
        let aligned_pairs = get_aligned_pairs_simple(&record);
        println!("Found {} aligned positions", aligned_pairs.len());
        
        // Get aligned pairs with sequence (similar to with_seq=True)
        let aligned_pairs_with_seq = get_aligned_pairs(&record, true);
        
        // Example: Print the first 5 pairs
        for (i, (read_pos, ref_pos, seq)) in aligned_pairs_with_seq.iter().take(5).enumerate() {
            let base = seq.map(|b| b as char).unwrap_or('-');
            println!("Pair {}: read_pos={:?}, ref_pos={:?}, base={}", 
                     i, read_pos, ref_pos, base);
        }
        
        break; // Just process the first record for this example
    }
    
    Ok(())
}