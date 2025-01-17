use agg::*;
use std::{collections::HashSet, path::PathBuf};
use std::collections::HashMap;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use std::fs::File;
use std::path::Path;
use std::io::Write;
use crate::agg;
use indicatif::ProgressBar;
use rayon::prelude::*; 
use bio::io::fasta::{Reader, Record};


pub fn find_ref_edge(
    graph: &GraphicalGenome,
    src: &str,
    dst: &str, 
    refstrain: &str,
    k: usize
) -> String {
    // Create HashSet with single refstrain
    let strains = HashSet::from([refstrain.to_string()]);
    
    // Find all paths between anchors
    let paths = agg::FindAllPathBetweenAnchors::new(&graph, &src, &dst, strains);
    let subpaths = &paths.subpath;
    
    // Return empty string if no paths found
    if subpaths.is_empty() {
        return String::new();
    }
    
    // Get sequence from first path
    let (path, _strain) = &subpaths[0];
    let seq = agg::reconstruct_path_seq(&graph, path);
    
    // Return sequence with src anchor and ref_edge sequence
    seq[..seq.len()-k].to_string()
}

fn mask_ns(seq: &str) -> String {
    seq.chars()
        .map(|c| {
            let upper_c = c.to_ascii_uppercase();
            if !['A', 'G', 'C', 'T'].contains(&upper_c) {
                c.to_ascii_lowercase()
            } else {
                upper_c
            }
        })
        .collect()
}

fn alignment_to_cigar(operations: &[AlignmentOperation]) -> String {
    let mut cigar: Vec<(usize, char)> = Vec::new();
    
    for op in operations {
        let cigar_op = match op {
            AlignmentOperation::Match => '=',
            AlignmentOperation::Subst => 'X',
            AlignmentOperation::Del => 'D',
            AlignmentOperation::Ins => 'I',
            AlignmentOperation::Xclip(_) => 'S',
            AlignmentOperation::Yclip(_) => 'S'
        };
        
        if !cigar.is_empty() && cigar.last().unwrap().1 == cigar_op {
            cigar.last_mut().unwrap().0 += 1;
        } else {
            cigar.push((1, cigar_op));
        }
    }

    cigar.iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

fn gap_open_aligner(reference: &str, sequence: &str) -> String {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // Create an aligner with the same scoring parameters
    let mut aligner = Aligner::with_capacity(sequence.len(),reference.len(), -5, -1, &score);  // match_score=0, mismatch_score=-6, gap_open=-5, gap_extend=-3

    // Perform the alignment
    let alignment = aligner.global(
        sequence.as_bytes(),
        reference.as_bytes()
    );

    // Get the aligned sequences
    let cigar = alignment_to_cigar(&alignment.operations);
    // println!("{:?}", cigar);

    cigar
}


pub fn generate_cigar (mut graph: GraphicalGenome,  ref_strain: &str, k: usize, maxlength : usize, minimal_read_count: usize) -> GraphicalGenome {
    let mut edgelist: Vec<_> = graph.edges.keys().cloned().collect();
    edgelist.sort();
    let bar = ProgressBar::new(edgelist.len() as u64);

    let updates:Vec<_> = edgelist.par_iter().map(|edge| {
        bar.inc(1);
        let allele_count = graph.edges.get(edge)
            .and_then(|e| e.get("reads"))
            .map_or(Vec::new(), |reads| reads.as_array().unwrap_or(&Vec::new()).to_vec())
            .len();

        if allele_count <= minimal_read_count {
            return (edge.clone(), None);
        }

        let src = &graph.edges.get(edge)
            .expect("Edge not found") 
            .get("src").expect("No src in this edge!")
            .as_array().expect("not an array")
            .first().expect("no vallue in src list")
            .as_str().expect("src not a string");
        let src_seq = if *src == "SOURCE" {
            ""
        } else {
            graph.anchor
                .get(*src)
                .and_then(|anchor| anchor.get("seq")
                                                 .and_then(|seq| seq.as_str())   )
                .unwrap_or("")
        };

        let dst = &graph.edges.get(edge)
            .expect("Edge not found")
            .get("dst").expect("No dst in this edge!")
            .as_array().expect("not an array")
            .first().expect("no value in dst list")
            .as_str().expect("dst not a string");
        // println!("src, {}, dst, {}", src, dst);
        let ref_seq = find_ref_edge(&graph, src, dst, ref_strain, k);
        let ref_length = ref_seq.len();
        let alt_sequence = src_seq.to_string() + &graph.edges.get(edge).expect("edge not found").get("seq").and_then(|v| v.as_str()).unwrap_or("");
        let alt_seq = mask_ns(&alt_sequence);
        let alt_length = alt_seq.len();

        if ref_length > maxlength || ref_length < k {
            return (edge.clone(), None);
        }

        if alt_seq == ref_seq {
            return (edge.clone(), Some(format!("{}=", ref_length)));
        }

        if alt_length > maxlength {
            return (edge.clone(), None);
        }
        let cigar = gap_open_aligner(&ref_seq, &alt_seq);
        std::thread::sleep(std::time::Duration::from_millis(50));

        (edge.clone(), Some(cigar))

    }).collect();

    // Apply the updates sequentially after parallel computation
    for (edge, maybe_variant) in updates {
        if let Some(variant) = maybe_variant {
            if let Some(edge_data) = graph.edges.get_mut(&edge) {
                edge_data.as_object_mut()
                    .expect("Edge data is not an object")
                    .insert("variants".to_string(), serde_json::Value::String(variant));
            }
        }
    }
    bar.finish();
    graph
}

#[derive(Debug, Clone)]
pub struct Variant {
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    pub allele_count: usize,
}

pub fn get_variants_from_cigar (cigar: &str, ref_seq: &str, alt_seq: &str, ref_start: usize, allelecount: usize) -> (Vec<Variant>, HashMap<usize, usize>) {
    let mut poscount = HashMap::new();
    let mut variants = Vec::new();
    let mut ref_pos = 0;
    let mut alt_pos = 0;

    let mut operations: Vec<(usize, char)> = Vec::new();
    let mut num = String::new();

    for c in cigar.trim_matches('"').chars() {  // trim_matches removes quotes at start and end
        if c.is_digit(10) {
            num.push(c);
        } else {
            let number = num.parse::<usize>().expect("number {}");
            operations.push((number, c));
            num.clear();
        }
    }

    for (length, op) in operations {
        match op {
            '=' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                ref_pos += length;
                alt_pos += length;
            },
            'X' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                    let ref_allele = &ref_seq[ref_pos + i..ref_pos + i + 1];
                    let alt_allele = &alt_seq[alt_pos + i..alt_pos + i + 1];
                    variants.push(Variant {
                        pos,
                        ref_allele: ref_allele.to_string(),
                        alt_allele: alt_allele.to_string(),
                        variant_type:"SNP".to_string(),
                        allele_count:allelecount,
                    });
                }
                ref_pos += length;
                alt_pos += length;
            }
            'I' => {
                let pos = ref_start + ref_pos;
                // *poscount.entry(pos).or_insert(0) += allelecount;
                let ref_allele = if ref_pos > 0 {
                    match ref_seq.get(ref_pos - 1..ref_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                            "-"
                        }
                    }
                } else { "-" };

                let alt_allele = match alt_seq.get(alt_pos - 1..alt_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                        "-"
                    }
                };
                variants.push(Variant {
                    pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "INS".to_string(),
                    allele_count:allelecount,
                });
                alt_pos += length;
            },

            'D' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                let pos = ref_start + ref_pos;
                let ref_allele = match ref_seq.get(ref_pos - 1..ref_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                        "-"
                    }
                };

                let alt_allele = if alt_pos > 0 {
                    match alt_seq.get(alt_pos - 1..alt_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                            "-"
                        }
                    }
                } else { "-" };

                if ref_pos > 0 && alt_pos > 0 {
                    if let (Some(r), Some(a)) = (ref_seq.get(ref_pos - 1..ref_pos), alt_seq.get(alt_pos - 1..alt_pos)) {
                        if r != a {
                            println!("{} {} {}", r, a, cigar);
                        }
                    }
                }

                variants.push(Variant {
                    pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "DEL".to_string(),
                    allele_count:allelecount,
                });
                ref_pos += length;
            }
            _ => (), //others skip

        }
    }
    (variants, poscount)
}

pub fn get_variant (graph: &mut GraphicalGenome, k: usize, ref_name: &str ) ->  (Vec<Variant>, HashMap<usize, usize>){
    let mut coverage = HashMap::new();
    let mut var = Vec::new();
    let mut edgelist: Vec<_> = graph.edges.keys().collect();
    edgelist.sort();
    for edge in edgelist {
        let allele_count = graph.edges.get(edge)
            .and_then(|e| e.get("reads"))
            .map_or(Vec::new(), |reads| reads.as_array().unwrap_or(&Vec::new()).to_vec())
            .len();
        
        let cigar = graph.edges.get(edge)
            .and_then(|e| e.get("variants"))
            .cloned()
            .unwrap_or_else(|| serde_json::Value::String("".to_string()));
        
        if cigar.as_str().unwrap_or("").is_empty() {
            continue;
        }

        let src = graph.incoming.get(edge)
            .and_then(|edges| edges.first())
            .expect("Edge should have source node");
        let dst = graph.outgoing.get(edge)
            .and_then(|edge| edge.first())
            .expect("Edge should have dst");
        if src == "SOURCE" || dst == "SINK"{
            continue
        }
        let refstart = graph.anchor[src]["pos"]
            .as_i64()
            .expect("Position should be a integer");

        let ref_seq = find_ref_edge(&graph, src, dst, ref_name, k);
        let src_seq = graph.anchor
                .get(src)
                .and_then(|anchor| anchor.get("seq")
                                                 .and_then(|seq| seq.as_str())   )
                .unwrap_or("");
        let alt_sequence = src_seq.to_string() + &graph.edges.get(edge).expect("edge not found").get("seq").and_then(|v| v.as_str()).unwrap_or("");
        let (variants, poscounts) = get_variants_from_cigar(&cigar.to_string(), &ref_seq, &alt_sequence, refstart as usize, allele_count);
        var.extend(variants);

        for (pos, count) in poscounts.iter(){
            *coverage.entry(*pos).or_insert(0) += count;
        }

    }
    (var, coverage)
}

fn collapse_identical_records(variants: Vec<Variant>) -> Vec<Variant> {
    if variants.is_empty() {
        return Vec::new();
    }

    let mut collapsed = HashMap::new();
   
    for current_var in variants {
        
        let mut pos = current_var.pos;
        let mut ref_allele = current_var.ref_allele;
        let mut alt_allele = current_var.alt_allele;
        let mut variant_type = current_var.variant_type;
        let mut allele_count = current_var.allele_count;

        let key = (pos, ref_allele.clone(), alt_allele.clone(), variant_type.clone());
        *collapsed.entry(key).or_insert(0) += allele_count;
    }
    collapsed.into_iter()
    .map(|((pos, ref_allele, alt_allele, variant_type), allele_count)| {
        Variant {
            pos,
            ref_allele,
            alt_allele,
            variant_type,
            allele_count,
        }
    })
    .collect()

}

fn format_vcf_record(variant: &Variant, coverage: HashMap<usize, usize>, reference_seq : &str) -> String {
    // Add AC (allele count) to INFO field
    let read_depth = coverage.get(&variant.pos).unwrap_or(&0);
    let allele_frequency = if *read_depth == 0{
        0.0
    }else{
        variant.allele_count as f32 / *read_depth as f32
    };
    let refpos = variant.pos;
    let refstart = if refpos >=5 {refpos - 5 } else {0};
    let refend = if refpos + 5 < reference_seq.len() {refpos + 5} else {reference_seq.len()};
    let ref_seq = &reference_seq[refstart..refend];

    let info = format!("RC={};DP={};AF={};RF={}", variant.allele_count, read_depth, allele_frequency, ref_seq);
    
    match variant.variant_type.as_str() {
        "SNP" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}",
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
            info
        ),
        "INS" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}",
            variant.pos,
            variant.ref_allele,
            variant.alt_allele,
            info
        ),
        "DEL" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}",
            variant.pos,
            variant.ref_allele,
            variant.alt_allele,
            info
        ),
        _ => panic!("Unknown variant type")
    }
}

fn write_vcf(variants: &[Variant], coverage: HashMap<usize, usize>,output_file: &str, minimal_ac: usize, reference:&str) -> std::io::Result<()> {
    let mut file = File::create(Path::new(output_file))?;

    // Write VCF header
    writeln!(file, "##fileformat=VCFv4.2")?;
    writeln!(file, "##reference=chrM")?;
    writeln!(file, "##contig=<ID=chrM,length=16569>")?;
    writeln!(file, "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Read Count\">")?;
    writeln!(file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")?;
    writeln!(file, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">")?;
    writeln!(file, "##INFO=<ID=RF,Number=1,Type=String,Description=\"Reference Sequences\">")?;
    writeln!(file, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

    // Sort variants by position
    let mut sorted_variants = variants.to_vec();
    sorted_variants.sort_by(|a, b| {
        // First compare positions
        a.pos.cmp(&b.pos)
            // Then compare variant types (to ensure consistent ordering)
            .then(a.variant_type.cmp(&b.variant_type))
            // Then compare ref alleles
            .then(a.ref_allele.cmp(&b.ref_allele))
            // Then compare alt alleles
            .then(a.alt_allele.cmp(&b.alt_allele))
    });

    // Write variant records
    for variant in sorted_variants {
        let allele_count = variant.allele_count;
        if allele_count < minimal_ac + 1{
            continue
        }
        let ref_allele = variant.ref_allele.as_str();
        if ref_allele.contains("N") {
            continue
        }
        
        writeln!(file, "{}", format_vcf_record(&variant, coverage.clone(), reference))?;
    }

    Ok(())
}

pub fn write_graph_from_graph(filename: &str, graph: &GraphicalGenome) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    let mut keys: Vec<_> = graph.anchor.keys().collect();
    keys.sort();
    for anchor in keys.iter() {
        let data = &graph.anchor[*anchor];
        let seq = data["seq"].as_str().unwrap_or_default();
        let mut data_clone = data.clone();
        data_clone.as_object_mut().unwrap().remove("seq");
        let json_string = serde_json::to_string(&data_clone).unwrap_or_else(|_| "{}".to_string());
        writeln!(file, "S\t{}\t{}\tPG:J:{}", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    let mut edge_keys: Vec<_> = graph.edges.keys().collect();
    edge_keys.sort();
    for edge in edge_keys.iter() {
        let edge_data = &graph.edges[*edge];
        let seq = edge_data["seq"].as_str().unwrap_or_default();
        let src = graph.incoming[*edge][0].clone();
        let dst = graph.outgoing[*edge][0].clone();
        let mut edge_data_clone = edge_data.clone();
        edge_data_clone.as_object_mut().unwrap().remove("seq");
        let json_string =
            serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if edge_data
            .get("reads")
            .and_then(|r| r.as_array())
            .map_or(false, |arr| !arr.is_empty())
        {
            format!(
                "S\t{}\t{}\tPG:J:{}\tRC:i:{}",
                edge,
                seq,
                json_string,
                edge_data
                    .get("reads")
                    .and_then(|r| r.as_array())
                    .map_or(0, |arr| arr.len())
            )
        } else {
            format!("S\t{}\t{}", edge, seq)
        };
        edge_output.push(formatted_string);

        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, edge));
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", edge, dst));
    }
    for s in edge_output {
        writeln!(file, "{}", s)?;
    }
    for l in link_output {
        writeln!(file, "{}", l)?;
    }
    Ok(())
}

pub fn start (graph_file: &PathBuf, ref_strain: &str, k: usize, maxlength: usize, minimal_ac :usize, output_file: &str, reference_path:&PathBuf) {
    let mut graph = agg::GraphicalGenome::load_graph(graph_file).unwrap();
    // generate cigar
    let mut graph_with_cigar = generate_cigar(graph, ref_strain, k, maxlength, 2);
    let (Variants, coverage )= get_variant(&mut graph_with_cigar, k, ref_strain);
    let collapsed_var = collapse_identical_records (Variants);

    // Read reference records into a vector
    let ref_reader = Reader::from_file(reference_path).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    
    // write vcf
    let _ = write_vcf(&collapsed_var, coverage, output_file, minimal_ac, &ref_seq);
    let graph_output = graph_file.with_extension("annotated.gfa");
    let _ = write_graph_from_graph(graph_output.to_str().unwrap(), &graph_with_cigar);

    

}