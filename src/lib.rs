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

    #[test]
    fn test_align() {
        let a = Sequence::from_nucleotides(vec!(A, C, T, A, C));
        let b = Sequence::from_nucleotides(vec!(C, T, A, G));
        let (s, t) = align(&a, &b, standard_cost_matrix(), linear_gap());
        assert_eq!((Sequence::from_nucleotides(vec!(A, C, T, A, C)), Sequence::from_nucleotides(vec!(Gap, C, T, A, G))), (s, t));
        
        let a = Sequence::from_nucleotides(vec!(T, T, T, C, C, C));
        let b = Sequence::from_nucleotides(vec!(C, C, C));
        let (s, t) = align(&a, &b, standard_cost_matrix(), linear_gap());
        assert_eq!((Sequence::from_nucleotides(vec!(T, T, T, C, C, C)), Sequence::from_nucleotides(vec!(Gap, Gap, Gap, C, C, C))), (s, t));

        let a = Sequence::from_nucleotides(vec!(A, G, T, C, G, G, T, A, A, A));
        let b = Sequence::from_nucleotides(vec!(T, A, G, G, G, A, A, A, T, T));
        let (s, t, _) = needleman_wunsch(&a, &b, standard_cost_matrix(), linear_gap()); 
        assert_eq!((s, t), align(&a, &b, standard_cost_matrix(), linear_gap()));
     }
}

//We will represent a sequence of nucleotides as a vec of nucleotides. 
pub mod sequences {
    use std::rc::Rc; 
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

        //Concatinate two sequences and return the result. 
        pub fn concat(&self, other : &Sequence) -> Sequence {
            let mut s = self.nucleotides.clone(); 
            let mut t = other.nucleotides.clone(); 
            s.append(&mut t); 
            Sequence::from_nucleotides(s.to_vec())
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
    pub fn score_alignments<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(s : Sequence, t : Sequence, sub_matrix : Rc<M> , cost_func : Rc<C> ) -> i64 {
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
    pub fn standard_cost_matrix() -> Rc<impl Fn(Nucleotide, Nucleotide) -> i64> {
        let c = | a : Nucleotide , b : Nucleotide | 
        if a == b {
            1
        }else{
            -1 
        }; 
        Rc::new(c)
    }

    //Linear gap function. (-1 per gap.)
    pub fn linear_gap() -> Rc<impl Fn(i64) -> i64> {
        Rc::new(|_i| -1)
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

    pub fn needleman_wunsch<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(a : &Sequence, b : &Sequence, sub_matrix : Rc<M>, cost_func : Rc<C>) -> (Sequence, Sequence, i64) {
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
                let b = NwCell{ score: x[i][j-1].score + cost_func(gap_length(&x, i, j, Rc::clone(&sub_matrix), Rc::clone(&cost_func), true)), gaps:(true, false) };
                let c = NwCell{ score: x[i-1][j].score + cost_func(gap_length(&x, i, j, Rc::clone(&sub_matrix), Rc::clone(&cost_func), false)), gaps:(false, true) };
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
    fn gap_length<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(x : &Vec<Vec<NwCell>>, i : usize, j : usize, _sub_matrix : Rc<M>, _cost_func : Rc<C>, s_or_t : bool) -> i64 {
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

    //Memoryless version of needleman wunsch for Hirschberg's algorithm to align sequences more efficiently. 
    fn low_mem_nw<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(a : &Sequence, b : &Sequence, sub_matrix : Rc<M>, cost_func : Rc<C>) -> Vec<i64> {
        let s = &a.nucleotides; 
        let t = &b.nucleotides; 

        //The lines that will hold the scores. 
        let mut v; 
        let mut u = vec!(0 as i64; s.len()+1);
        for i in 0..(s.len()+1) { u[i] = (i as i64) * -1 } 

        //Iterate down the list and compute the scores. 
        for i in 1..(t.len()+1){
            //Ratchet. 
            v = u.clone(); 
            u = vec!(0 as i64; s.len()+1); 

            //Iterate down the list and do needleman wunsch on the two lists. 
            u[0] = ((i) as i64) * -1;
            for j in 1..(s.len()+1) { 
                let a = v[j-1] + sub_matrix(t[i-1], s[j-1]);
                let b = v[j] + cost_func(1); //TODO: Support more gaps. 
                let c = u[j-1] + cost_func(1); //TODO: Support more gaps.
                u[j] = max(a, max(b, c)); 
            }
        }
        //Return the last row of the matrix. 
        u
    }

    //Use Hirschberg's algorithm to align two sequences. (More efficient than NW.)
    pub fn align<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(a : &Sequence, b : &Sequence, sub_matrix : Rc<M>, cost_func : Rc<C>) -> (Sequence, Sequence) {
        //Get the sequences. 
        let s = &a.nucleotides; 
        let t = &b.nucleotides;

        //Base cases. 
        if s.len() == 0 {
            let a = Sequence::from_nucleotides(vec!(Gap; t.len()));
            let b = Sequence::from_nucleotides(t.to_vec()); 
            return (a, b);
        }
        if t.len() == 0 {
            let a = Sequence::from_nucleotides(s.to_vec()); 
            let b = Sequence::from_nucleotides(vec!(Gap; s.len()));
            return (a, b);
        }
        if t.len() == 1 || s.len() == 1 {
            let (a, b, _) = needleman_wunsch(a, b, sub_matrix, cost_func); 
            return (a, b);
        }

        //Split S. 
        let s_split : usize = s.len()/2; 
        let sl = s[..s_split].to_vec(); 
        let mut sr = s[s_split..].to_vec(); 

        //Figure out where to split t. 
        let mut t_prime = t.clone(); 
        t_prime.reverse(); 
        sr.reverse(); 
        let n1 = low_mem_nw(&b, &Sequence::from_nucleotides(sl.clone()), Rc::clone(&sub_matrix), Rc::clone(&cost_func));
        let mut n2 = low_mem_nw(&Sequence::from_nucleotides(t_prime), &Sequence::from_nucleotides(sr.clone()), Rc::clone(&sub_matrix), Rc::clone(&cost_func));
        sr.reverse(); 
        n2.reverse(); 
        //Find best split for t.
        let pairs = n1.iter().zip(n2.iter()).map(|(a,b)| a+b).enumerate();
        let t_split = match pairs.max_by_key(|(_i,v)| *v) {
            Some((i, _v)) => i, 
            None => 0
        };
        
        let tl = t[..t_split].to_vec(); 
        let tr = t[t_split..].to_vec(); 

        //Recursive call.
        let (xl, zl) = align(&Sequence::from_nucleotides(sl), &Sequence::from_nucleotides(tl), Rc::clone(&sub_matrix), Rc::clone(&cost_func));
        let (xr, zr) = align(&Sequence::from_nucleotides(sr), &Sequence::from_nucleotides(tr), Rc::clone(&sub_matrix), Rc::clone(&cost_func));

        return (xl.concat(&xr), zl.concat(&zr));
    }
}
