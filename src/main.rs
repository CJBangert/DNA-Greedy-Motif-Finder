use std::{fs};
use substring::Substring;
fn main() {
    println!("In file {}", "test.txt");

    let contents = fs::read_to_string("test.txt")
        .expect("Should have been able to read the file");

    let args = contents.split_whitespace().collect::<Vec<&str>>();
    
    //k is length of motif to search for
    let k:u32 = args.get(0).unwrap().to_string().parse::<u32>().unwrap();
    //t is number dna sequences to search
    let t:u32 = args.get(1).unwrap().to_string().parse::<u32>().unwrap();
    //rest of args input are the dna strings
    //
    let dna = args[2..].to_vec();
    println!("{} {}", k,t);

    let mut best_motifs  = Vec::<String>::new();

    let mut best_score:i32 = -1;

    let kmer_size = u32::try_from(dna.get(0).unwrap().len()).unwrap();

    for m in 0..=(kmer_size - k){
        //since arr size is read from a file and t is not known at compile time, allocate entire vec at once for efficiency
        let mut motifs = Vec::<String>::with_capacity(t.try_into().unwrap());

        motifs.push(String::from(dna.get(0).unwrap().substring(m.try_into().unwrap(), (m+k).try_into().unwrap())));

        for i in 1..t{
            let mut prev_motifs = Vec::<String>::with_capacity(i.try_into().unwrap());
            for x in 0..i{
                prev_motifs.push(motifs.get::<usize>(x.try_into().unwrap()).unwrap().clone());
            }
            let profile = profile_from_motifs(&prev_motifs);

            let elt = profile_most_probable_kmer(dna.get(i as usize).unwrap(), k, profile);
            motifs.insert(i.try_into().unwrap(),elt ); 
        }

        let curr_score = score_motifs(&motifs);
        if curr_score > best_score {
            best_score = curr_score;
            best_motifs = motifs;
        }
    }

    for i in best_motifs {
        println!("{}",i);
    }
}



fn profile_from_motifs(motifs: &Vec<String>) -> Vec<Vec<f64>>{
    //vec with capacity 4 x len(motifs[0])
    let mut profile:Vec::<Vec::<f64>> = vec![vec![0.0; motifs.get(0).unwrap().len()];4];
    for seq in  0..motifs.len() {
        let curr = &*motifs.get(seq).unwrap();
        for i in 0..curr.len() {
            let c = curr.chars().collect::<Vec<char>>()[i];
            profile[nucleotide_to_index(c)][i] = profile[nucleotide_to_index(c)][i] + 1.0;
        }

        for i in 0..profile.len(){
            for j in 0..profile[0].len() {
                profile[i][j] = profile[i][j] / motifs.len() as f64;
            }
        }
    }
    profile
}

fn profile_most_probable_kmer(sequence: &str, k: u32, profile:Vec<Vec<f64>>) -> String{
    let  mut best:String = String::from("");
    let mut best_p:f64 = -1.0;
    for i in 0..=sequence.len() - k as usize {
        let sub = sequence.substring(i,i+(k as usize)).to_ascii_uppercase();
        let mut prob:f64 = 1.0;
        for j in 0..sub.len() {
            prob = prob *profile[nucleotide_to_index(sub.chars().collect::<Vec<char>>()[j])][j];
        }
        if prob > best_p {
            best = sub;
            best_p = prob;
        }
    }
    best
}

fn score_motifs(motifs: &Vec<String>) -> i32{
    let profile:Vec<Vec<f64>> = profile_from_motifs(motifs);
    let mut consensus = String::new();
    for i in 0..motifs[0].len() {
        let mut max_prob:f64 = -1.0;
        let mut max_char:char = '\0';
        for j in 0..4 {
            if profile[j][i] > max_prob{
                max_prob = profile[j][i];
                max_char = usize_to_nucleotide(j);
            }
        }
        consensus.push(max_char);
    }
    
    let mut score:i32 = 0;
    for i in 0..consensus.len() {
        let c:char = consensus.chars().collect::<Vec<char>>()[i];
        for j in 0..motifs.len() {
            if motifs[j].chars().collect::<Vec<char>>()[i] == c {
                score = score + 1;
            }
        }
    }
    score
}

fn usize_to_nucleotide(i: usize) -> char {
    match i {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => '\0'
    }
}

fn nucleotide_to_index (n:char) -> usize{
    match n {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 4 //unreachable and will throw an error if somehow reached
    }
}