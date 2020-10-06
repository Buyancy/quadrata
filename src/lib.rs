#[cfg(test)]
mod tests {
    use crate::sequences::*;
    use crate::sequences::Nucleotide::*; 

    #[test]
    fn test_score() {
        let a = Sequence::from_nucleotides(vec!(A, C, G, T, T, A));
        let b = Sequence::from_nucleotides(vec!(T, C, T, A, T, A, G));
        let z = score_alignments(a, b, standard_cost_matrix(), linear_gap());
        assert_eq!(z, 0); 
    }

    #[test]
    fn test_print_alignment() {
        let a = Sequence::from_nucleotides(vec!(A, C, G, T, T, A));
        let b = Sequence::from_nucleotides(vec!(T, C, T, A, T, A, G));
        println!("Printing alignment:");
        print_alignment(&a, &b);
    }

    #[test]
    fn test_needleman_wunsch() {
        let a = Sequence::from_nucleotides(vec!(A, C, T, A, C));
        let b = Sequence::from_nucleotides(vec!(C, T, A, G));
        let (s, t, score) = needleman_wunsch(&a, &b, standard_cost_matrix(), linear_gap());
        assert_eq!((Sequence::from_nucleotides(vec!(A, C, T, A, C)), Sequence::from_nucleotides(vec!(Gap, C, T, A, G)), 1), (s, t, score));
        
        let a = Sequence::from_nucleotides(vec!(T, T, T, C, C, C));
        let b = Sequence::from_nucleotides(vec!(C, C, C));
        let (s, t, score) = needleman_wunsch(&a, &b, standard_cost_matrix(), linear_gap());
        assert_eq!((Sequence::from_nucleotides(vec!(T, T, T, C, C, C)), Sequence::from_nucleotides(vec!(Gap, Gap, Gap, C, C, C)), 0), (s, t, score));
    }
}

//We will represent a sequence of nucleotides as a vec of nucleotides. 
pub mod sequences {
    use std::fmt;
    use std::cmp::{min, max, Ordering};
    use Nucleotide::*; 

    #[derive(PartialEq, Copy, Clone, Debug)]
    pub enum Nucleotide {
        A,
        C,
        G, 
        T, 
        Gap
    }

    impl fmt::Display for Nucleotide {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Nucleotide::A =>  write!(f, "A"),
                Nucleotide::C =>  write!(f, "C"),
                Nucleotide::G =>  write!(f, "G"),
                Nucleotide::T =>  write!(f, "T"),
                Nucleotide::Gap =>  write!(f, "-")
            }
           
        }
    }

    impl Nucleotide {
        pub fn get_char(&self) -> char{
            match self {
                A => 'A', 
                C => 'C',
                G => 'G',
                T => 'T',
                Gap => '-',
            }
        }
    }

    #[derive(PartialEq, Debug)]
    pub struct Sequence {
        pub label : String, 
        pub nucleotides : Vec<Nucleotide>
    }

    impl Sequence { 
        //Return a new sequence that is the complement of the old sequence.
        pub fn complement(&self) -> Sequence {
            let s = self.nucleotides.iter().map( |n| 
                match n {
                    A => T, 
                    C => G, 
                    G => C, 
                    T => T, 
                    Gap => Gap,
                }
            ).collect::<Vec<Nucleotide>>();
            Sequence { label:self.label.clone(), nucleotides:s } 
        }

        //Create a blank nucleotide sequence. 
        pub fn new() -> Sequence {
            Sequence { label:String::from(""), nucleotides:vec!() }
        }

        //Create a nucleotide sequence from a vec of nucleotides. 
        pub fn from_nucleotides(s : Vec<Nucleotide>) -> Sequence {
            Sequence { label:String::from(""),nucleotides:s }
        }
    }

    //Assumes that the sequences are the same length. 
    pub fn score_alignments<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(s : Sequence, t : Sequence, sub_matrix : M , cost_func : C ) -> i64 {
        let mut score = 0; 
        let mut gap_len = 0; 
        for (a, b) in s.nucleotides.iter().zip(t.nucleotides.iter()){ 
            if *a == Gap || *b == Gap {
                gap_len += 1; 
                score += cost_func(gap_len); 
            } else {
                gap_len = 0; 
                score += sub_matrix(*a, *b); 
            }
        }
        score
    }

    //Standard cost matrix. 1 for match -1 for mismatch. 
    pub fn standard_cost_matrix() -> impl Fn(Nucleotide, Nucleotide) -> i64 {
        | a : Nucleotide , b : Nucleotide | 
        if a == b {
            1
        }else{
            -1 
        }
    }

    //Linear gap function. (-1 per gap.)
    pub fn linear_gap() -> impl Fn(i64) -> i64 {
        |_i| -1
    }

    //Print the alignment of two sequences. 
    use ansi_term::Colour::Red;
    use ansi_term::Colour::Blue;
    use ansi_term::Colour::Green;
    use ansi_term::Colour::Yellow;
    use ansi_term::Colour::Black;
    pub fn print_alignment(s : &Sequence, t : &Sequence) {
        let max_chars = 50; 
        //Print s with colors. 
        for i in 0..min(s.nucleotides.len(), max_chars) {
            match s.nucleotides[i] {
                A => print!("{}", Red.paint("A")), 
                C => print!("{}", Blue.paint("C")), 
                G => print!("{}", Green.paint("G")), 
                T => print!("{}", Yellow.paint("T")), 
                Gap => print!("{}", Black.paint("-")), 
            } 
        }
        println!(); 

        //Print the alignment characters. 
        for (a,b) in s.nucleotides.iter().zip(t.nucleotides.iter()).take(max_chars) {
            if a == b {
                print!("{}", Blue.paint("|"));
            } else {
                print!(" ");
            }
        }
        println!(); 

        //Print t with colors. 
        for i in 0..min(t.nucleotides.len(), max_chars) {
            match t.nucleotides[i] {
                A => print!("{}", Red.paint("A")), 
                C => print!("{}", Blue.paint("C")), 
                G => print!("{}", Green.paint("G")), 
                T => print!("{}", Yellow.paint("T")), 
                Gap => print!("{}", Black.paint("-")), 
            }
        }
        println!(); 
    }

    //Align two sequences. Returns the aligned sequences and the optimal alignment score. 
    #[derive(Debug, Clone, Copy, Eq)]
    struct NwCell {
        score: i64, 
        gaps : (bool, bool) //Whether there is a gap in s or t because of this cell. 
    }

    impl PartialOrd for NwCell {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.score.cmp(&other.score))
        }
    }
    impl PartialEq for NwCell {
        fn eq(&self, other: &Self) -> bool {
            self.score == other.score
        }
    }
    impl std::cmp::Ord for NwCell {
        fn cmp(&self, other: &Self) -> Ordering {
            self.score.cmp(&other.score)
        }
    }

    pub fn needleman_wunsch<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(a : &Sequence, b : &Sequence, sub_matrix : M, cost_func : C) -> (Sequence, Sequence, i64) {
        let mut s = vec!(Gap); 
        for i in &a.nucleotides {
            s.push(*i); 
        }
        let mut t = vec!(Gap);
        for i in &b.nucleotides {
            t.push(*i); 
        }
        let mut x = vec![vec![NwCell{ score:0, gaps: (false, false) };t.len()]; s.len()];
        //Fill in the top and bottom rows. 
        for i in 1..s.len() {
            x[i][0].score = x[i-1][0].score + cost_func(i as i64);
            x[i][0].gaps = (false, true); 
        }
        for j in 1..t.len() {
            x[0][j].score = x[0][j-1].score + cost_func(j as i64);
            x[0][j].gaps = (true, false); 
        }
        //Iterate over the matrix, filling in values. 
        for i in 1..s.len() {
            for j in 1..t.len() {
                let a = NwCell{ score: x[i-1][j-1].score + sub_matrix(s[i], t[j]), gaps:(false, false) };
                let b = NwCell{ score: x[i][j-1].score + cost_func(gap_length(&x, i, j, &sub_matrix, &cost_func, true)), gaps:(true, false) };
                let c = NwCell{ score: x[i-1][j].score + cost_func(gap_length(&x, i, j, &sub_matrix, &cost_func, false)), gaps:(false, true) };
                x[i][j] = max(a, max(b, c));
            }
        }
        //The back trace algorithm to get the actual alignment. 
        let mut aligned_s : Vec<Nucleotide> = Vec::new(); 
        let mut aligned_t : Vec<Nucleotide> = Vec::new(); 
        let mut i = s.len() - 1;
        let mut j = t.len() - 1;
        loop {
            match x[i][j].gaps {
                (false, false) => {
                    //No gap.
                    aligned_s.push(s[i]); 
                    aligned_t.push(t[j]);
                    i -= 1; 
                    j -= 1; 
                }, 
                (true, false) => {
                    //Gap in s. 
                    aligned_s.push(Gap); 
                    aligned_t.push(t[j as usize]);
                    j -= 1; 
                }, 
                (false, true) => {
                    //Gap in t. 
                    aligned_t.push(Gap); 
                    aligned_s.push(s[i]);
                    i -= 1;
                }, 
                _ => (panic!("Encounterd invalid direction in direction matrix."))
            }
            if i == 0 && j == 0 {
                break; 
            }
        }
        aligned_s.reverse(); 
        aligned_t.reverse(); 
        (Sequence::from_nucleotides(aligned_s), Sequence::from_nucleotides(aligned_t), x[s.len()-1][t.len()-1].score)
    }

    //Helper function for determening the gap length in the dynamic programming matrix for needleman wunsch. 
    fn gap_length<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(x : &Vec<Vec<NwCell>>, i : usize, j : usize, _sub_matrix : &M, _cost_func : &C, s_or_t : bool) -> i64 {
        let mut gap_len = 1; 
        let mut i = i; 
        let mut j = j; 
        //Iterate through and test for current gaps. 
        if s_or_t == true {
            //Testing for gap in s. 
            loop{
                if i > 0 && x[i][j].gaps == (true, false) {
                    gap_len += 1; 
                    i -= 1; 
                } else {
                    return gap_len;
                }
            }
        }else {
            //Testing for gap in t. 
            loop{
                if j > 0 && x[i][j].gaps == (false, true) {
                    gap_len += 1; 
                    j -= 1; 
                } else {
                    return gap_len;
                }
            }
        }
    }
}
