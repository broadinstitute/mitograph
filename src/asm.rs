use crate::agg::*;
use std::{collections::HashSet, path::PathBuf, fs::File, io::{self, Write}};

pub fn find_most_supported_edge(graph: &GraphicalGenome, src: String) -> String {
    let empty_vec = Vec::new();
    let outgoinglist = &graph.outgoing.get(&src).unwrap_or(&empty_vec);
    let mut m = 0;
    let mut most_supported_edge = "".to_string();
    for edge in outgoinglist.iter(){
        let read_count = graph.edges[edge]["reads"].as_array().unwrap_or(&Vec::new()).len();
        if read_count > m {
            m = read_count;
            most_supported_edge = edge.clone();
        }
    }
    most_supported_edge

}

pub fn construct_major_haplotype(graph:GraphicalGenome) -> String {
    let mut anchorlist: Vec<_> = graph.anchor.keys().collect();
    anchorlist.sort();
    let mut src = anchorlist.first().unwrap().to_string();
    let mut next_edge = "".to_string();
    let mut dst = anchorlist.last().unwrap().to_string();
    let mut haplotype = String::new();
    // println!("{}", src);

    while dst > src {
        next_edge = find_most_supported_edge(&graph, src.clone().to_string());
        let anchor_seq = graph.anchor.get(&src)
            .and_then(|v| v.get("seq"))
            .and_then(|v| v.as_str())
            .unwrap_or(""); 
        let edge_seq = graph.edges.get(&next_edge)
            .and_then(|v| v.get("seq"))
            .and_then(|v| v.as_str())
            .unwrap_or("");       
        haplotype.push_str(anchor_seq);
        haplotype.push_str(edge_seq);
        src = graph.edges[&next_edge]["dst"].as_array()
            .unwrap()
            .first()
            .unwrap()
            .as_str()
            .unwrap()
            .to_string();        
        // println!("{}, {}, {}", src, dst, next_edge);
    }

    let anchor_seq = graph.anchor.get(&src)
        .and_then(|v| v.get("seq"))
        .and_then(|v| v.as_str())
        .unwrap_or(""); 
    let edge_seq = graph.edges.get(&next_edge)
        .and_then(|v| v.get("seq"))
        .and_then(|v| v.as_str())
        .unwrap_or("");
    haplotype.push_str(anchor_seq);
    haplotype.push_str(edge_seq);
    haplotype

}

fn write_fasta(outputfile: &PathBuf, sequence: String, header: &str) -> io::Result<()>{
    let mut index = 0;
    let mut file = File::create(outputfile)?;

    let header = format!(">{} \n", header);
    file.write_all(header.as_bytes())?;

    let chars_per_line = 60;
    let sequence_len = sequence.len();
    let full_lines = sequence_len / chars_per_line;

    for i in 0..full_lines {
        let start = i * chars_per_line;
        let end = start + chars_per_line;
        writeln!(file, "{}", &sequence[start..end])?;
    }

    // Write any remaining characters that didn't make up a full line
    if sequence_len % chars_per_line != 0 {
        writeln!(file, "{}", &sequence[full_lines * chars_per_line..])?;
    }

    Ok(())
}

pub fn start (graph_file: &PathBuf, output_file: &PathBuf, header:&str) {
    let graph = GraphicalGenome::load_graph(graph_file).unwrap();
    let haplotype = construct_major_haplotype(graph);
    println!("{}", haplotype.len());

    // write fasta
    write_fasta(output_file, haplotype, header).unwrap();

}