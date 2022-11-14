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

    let mut bestMotifs  = Vec::<String>::new();

    let mut bestScore:i32 = -1;

    let kmerSize = u32::try_from(dna.get(0).unwrap().len()).unwrap();

    for m in 0..=(kmerSize - k){
        //since arr size is read from a file and t is not known at compile time, allocate entire vec at once for efficiency
        let mut motifs = Vec::<String>::with_capacity(t.try_into().unwrap());

        motifs.push(String::from(dna.get(0).unwrap().substring(m.try_into().unwrap(), (m+k).try_into().unwrap())));

        for i in 1..t{
            let mut prevMotifs = Vec::<String>::with_capacity(i.try_into().unwrap());
            for x in 0..i{
                prevMotifs.push(motifs.get::<usize>(x.try_into().unwrap()).unwrap().clone());
            }
            let profile = profileFromMotifs(&prevMotifs);

            let elt = profileMostProbableKmer(dna.get(i as usize).unwrap(), k, profile);
            motifs.insert(i.try_into().unwrap(),elt ); 
        }

        let currScore = scoreMotifs(&motifs);
        if currScore > bestScore {
            bestScore = currScore;
            bestMotifs = motifs;
        }
    }

    for i in bestMotifs {
        println!("{}",i);
    }
}


// fn greedyMotifSearch(k:u32, t:u32, dna:Vec<&str>) -> Vec<&str>{

//     let mut bestMotifs:Vec<&str>  = Vec::<&str>::new();

//     let mut bestScore:i32 = -1;

//     let kmerSize = u32::try_from(dna.get(0).unwrap().len()).unwrap();

//     for m in 0..=(kmerSize - k){
//         //since arr size is read from a file and t is not known at compile time, allocate entire vec at once for efficiency
//         let mut motifs = Vec::<&str>::with_capacity(t.try_into().unwrap());

//         motifs.push(dna.get(0).unwrap().substring(m.try_into().unwrap(), (m+k).try_into().unwrap()));

//         for i in 1..t{
//             let mut prevMotifs = Vec::<&str>::with_capacity(i.try_into().unwrap());
//             for x in 0..i{
//                 prevMotifs.push(motifs.get::<usize>(x.try_into().unwrap()).unwrap());
//             }
//             let profile = profileFromMotifs(&prevMotifs);
//             motifs.insert(i.try_into().unwrap(), profileMostProbableKmer(dna.get(i as usize).unwrap(), k, profile).as_str()); 
//         }

//         let currScore = scoreMotifs(&motifs);
//         if currScore > bestScore {
//             bestScore = currScore;
//             bestMotifs = motifs;
//         }
//     }




//     bestMotifs

// }

fn profileFromMotifs(motifs: &Vec<String>) -> Vec<Vec<f64>>{
    //vec with capacity 4 x len(motifs[0])
    let mut profile:Vec::<Vec::<f64>> = vec![vec![0.0; motifs.get(0).unwrap().len()];4];
    for seq in  0..motifs.len() {
        let curr = &*motifs.get(seq).unwrap();
        for i in 0..curr.len() {
            let c = curr.chars().collect::<Vec<char>>()[i];
            profile[nucleotideToIndex(c)][i] = profile[nucleotideToIndex(c)][i] + 1.0;
        }

        for i in 0..profile.len(){
            for j in 0..profile[0].len() {
                profile[i][j] = profile[i][j] / motifs.len() as f64;
            }
        }
    }
    profile
}

fn profileMostProbableKmer(sequence: &str, k: u32, profile:Vec<Vec<f64>>) -> String{
    let  mut best:String = String::from("");
    let mut bestP:f64 = -1.0;
    for i in 0..=sequence.len() - k as usize {
        let sub = sequence.substring(i,i+(k as usize)).to_ascii_uppercase();
        let mut prob:f64 = 1.0;
        for j in 0..sub.len() {
            prob = prob *profile[nucleotideToIndex(sub.chars().collect::<Vec<char>>()[j])][j];
        }
        if prob > bestP {
            best = sub;
            bestP = prob;
        }
    }
    best
}

fn scoreMotifs(motifs: &Vec<String>) -> i32{
    let profile:Vec<Vec<f64>> = profileFromMotifs(motifs);
    let mut consensus = String::new();
    for i in 0..motifs[0].len() {
        let mut maxProb:f64 = -1.0;
        let mut maxChar:char = '\0';
        for j in 0..4 {
            if profile[j][i] > maxProb{
                maxProb = profile[j][i];
                maxChar = usizeToNucleotide(j);
            }
        }
        consensus.push(maxChar);
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

fn usizeToNucleotide(i: usize) -> char {
    match i {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => '\0'
    }
}

fn nucleotideToIndex (n:char) -> usize{
    match n {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 4 //unreachable and will throw an error if somehow reached
    }
}