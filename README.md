# quadrata
Quadrata is a bioinformatics sequence alignment library designed for use in the rust programming language. My plan is for it to handle sequences of nucleotides (or peptides in the fucure) and be able to do operations on them such as align or score them.  

# Installation 
The library can be downloaded from this repo. It is not currently on crates.io but I plan to put it there once it gets a bit more functionality. 

# Basic Usage
The most basic unit of a nucleotide sequence is a single nucleotide. This is represented here as an enum that has all of the nucleotide options as well as the option for a gap for when we align two sequences. 
```rust
pub enum Nucleotide {
    A,
    C,
    G, 
    T, 
    Gap
}
```
To arrange and store these nucleotides in a linear sequence, we will use the Sequence struct. 
```rust 
pub struct Sequence {
    pub label : String, 
    pub nucleotides : Vec<Nucleotide>
}
```
The Sequence struct can be created more easily from `Sequence::new()` if you want a blank sequence you can fill in later, or from `Sequence::from(Vec<Nucleotide>)`if you wish to make one with nucleotides already filled in. 

As of now, the following functions are implemented and working. 
- `Sequence::complement(&self) -> Sequence`: Return the complement of the nucleotide sequence (gaps stay as gaps.)
- `Sequence::new() -> Sequence`: Return a new blank sequence with an empty label. 
- `Sequence::from_nucleotides(s : Vec<Nucleotide>) -> Sequence`: Return a new sequence with the nucleotides that it is given and a blank label. 
- `score_alignments<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(s : Sequence, t : Sequence, sub_matrix : M , cost_func : C ) -> i64`: Score the alignment of two sequences using the substitution matrix `sub_matrix` and the gap cost function `cost_func`. `M` and `C` can be passed in as closures. 
- `standard_cost_matrix() -> impl Fn(Nucleotide, Nucleotide) -> i64`: Returns the standard cost matrix. 1 for a match and -1 for a mismatch. 
- `linear_gap() -> impl Fn(i64) -> i64`: Returns the linear gap function where each gap costs -1.
- `print_alignment(s : &Sequence, t : &Sequence)`: Prints the alignment of two sequences to the screen. 
- `needleman_wunsch<M : Fn(Nucleotide, Nucleotide)-> i64, C : Fn(i64)->i64>(a : &Sequence, b : &Sequence, sub_matrix : M, cost_func : C) -> (Sequence, Sequence, i64)`: Aligns two sequences using the Needleman-Wunsch algorithm using the substitution matrix `sub_matrix` and the gap cost function `cost_func`. `M` and `C` can be passed in as closures. Returns a tuple of the two aligned sequences as well as their alignment score. 
