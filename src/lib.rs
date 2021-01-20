extern crate tabfile;

pub mod error;

use std::borrow::Cow;
use std::convert::TryInto;
use std::path::Path;

use tabfile::Tabfile;

use error::Error;

pub type Float = f32;
type TransitionProbabilities = [Float; 3];

#[derive(Copy, Clone, Debug)]
pub struct IndelTransitionProbabilities {
    pub outframe: Float, // frameshift
    pub inframe: Float, // not a frameshift
}

impl Default for IndelTransitionProbabilities {
    fn default() -> Self {
        Self{outframe: Float::NAN, inframe: Float::NAN}
    }
}

// X -> X will always be NaN
#[derive(Debug)]
pub struct SubstitutionRates {
    pub a: Float,
    pub c: Float,
    pub g: Float,
    pub t: Float,
}

impl SubstitutionRates {
    /// Dynamic field access to a specific nucleotide
    ///
    /// ```
    /// extern crate pattern_partition_prediction;
    /// use pattern_partition_prediction::SubstitutionRates;
    ///
    /// let rates = SubstitutionRates{
    ///   a: 0.1,
    ///   c: 0.2,
    ///   g: 0.3,
    ///   t: 0.4,
    /// };
    /// assert_eq!( rates.by_name( 'a' ), 0.1 );
    /// assert_eq!( rates.by_name( 'c' ), 0.2 );
    /// assert_eq!( rates.by_name( 'g' ), 0.3 );
    /// assert_eq!( rates.by_name( 't' ), 0.4 );
    /// ```
    pub fn by_name(&self, nucleotide: char) -> Float {
        match nucleotide {
            'A' | 'a' => self.a,
            'C' | 'c' => self.c,
            'G' | 'g' => self.g,
            'T' | 't' => self.t,
            _ => panic!("Not a valid nucleotide: {}", nucleotide),
        }
    }
}

/// A pattern partition prediction lookup table
///
/// Use the `new` method to read in a sequence-context point mutation probabilities file.
/// This file has the following format:
///
/// ```text
/// A->C ANAANVY 9.624201240304532e-06
/// A->C BNAANVY 7.1908831554369445e-06
/// A->C KNAANVR 5.345747358414394e-06
/// A->C MNAANVR 7.589541872637903e-06
/// A->C NAAAVTN 7.0416965447897126e-06
/// [and so on]
/// ```
/// Each column is separated by a single space.
/// The first column is the point mutation.
/// The second column is the sequence context where the mutation takes place (UIPAC notation)
/// The third column is the point mutation rate with this context
///
/// The reverse complement is handled automatically. All sequence must be different (even when you
/// expand the UIPAC code)
///

pub struct PaPaPred {
    radius: usize,
    lookup: Vec<TransitionProbabilities>, // ( odd kmer -> base 4 ) -> TransitionProbabilities
}

impl PaPaPred {
    /// Create a PaPaPred instance.
    ///
    /// The input file must consist of 3 space-separated columns:
    /// 1) The substitution in the form X->Y where X, Y in {ACGT}
    /// 2) The pattern in UIPAC notation.
    /// 3) The probability for the substitution with that pattern
    ///
    /// You provide the `path` to the sequence-context point mutation probabilities file.
    /// If you want to ensure that the patters have a minimum size, you can provide a value for
    /// `min_kmer_size` which will pad all patterns with Ns at the end. The minimum kmer size has
    /// to be an odd number because there has to be a single mutated base in the center of the
    /// sequence.
    pub fn new<P: AsRef<Path>>(path: P, min_kmer_size: Option<usize>) -> Result<PaPaPred, Error> {
        let mut kmer_size = 0; // will be detected from input size
        let template: Vec<TransitionProbabilities> = vec![[Float::NAN; 3]]; // missing values should be treated as NaNs
        let mut lookup: Vec<TransitionProbabilities> = Vec::new();
        for record_result in Tabfile::open(path)?.separator(' ') {
            let record = record_result?;
            let tokens = record.fields();

            if tokens[0].starts_with('#') || tokens.len() == 0 {
                continue // ignore header or comments
            }

            let substitution = tokens[0].parse::<Substitution>().unwrap();
            let uipac_context = pad_size(tokens[1], min_kmer_size.unwrap_or(1));
            let rate = tokens[2].parse::<Float>().unwrap();

            if kmer_size == 0 {
                kmer_size = uipac_context.len();
                if kmer_size % 2 == 0 {
                    return Err(Error::BadKmerSize(kmer_size));
                }
                lookup = template.repeat(4usize.pow(kmer_size as u32));
            } else if uipac_context.len() != kmer_size {
                let message = format!(
                    "The following line does not have length {}: {}",
                    kmer_size,
                    record.line()
                );
                return Err(Error::FileFormat { message });
            }

            for context in expand_uipac(&uipac_context) {
                let context_base4 = seq2base_four(&context)?;
                let index = substitution.transition_probabilities_index();
                lookup[context_base4][index] = rate;
                let complement_context = base4_reverse_complement(context_base4, kmer_size);
                let complement_index = substitution.complement().transition_probabilities_index();
                lookup[complement_context][complement_index] = rate;
            }
        }
        let radius = kmer_size / 2;
        Ok(PaPaPred { radius, lookup })
    }

    /// Query the mutation rates for a sequence context
    ///
    /// The middle position of the `seq` is assumed to be mutating.
    pub fn rates(&self, seq: &str) -> Result<SubstitutionRates, Error> {
        let index = seq2base_four(seq)?;
        let rates = self.lookup[index];
        let ref_base: Base = seq
            .chars()
            .nth(self.radius)
            .expect("sequence is too short")
            .try_into()
            .unwrap();
        Ok(match ref_base {
            Base::A => SubstitutionRates {
                a: Float::NAN,
                c: rates[0],
                g: rates[1],
                t: rates[2],
            },
            Base::C => SubstitutionRates {
                a: rates[0],
                c: Float::NAN,
                g: rates[1],
                t: rates[2],
            },
            Base::G => SubstitutionRates {
                a: rates[0],
                c: rates[1],
                g: Float::NAN,
                t: rates[2],
            },
            Base::T => SubstitutionRates {
                a: rates[0],
                c: rates[1],
                g: rates[2],
                t: Float::NAN,
            },
        })
    }

    pub fn kmer_size(&self) -> usize {
        1 + 2 * self.radius
    }

    pub fn radius(&self) -> usize {
        self.radius
    }
}

pub struct PaPaPredIndel {
    radius: usize,
    lookup: Vec<IndelTransitionProbabilities>, // ( oven kmer -> base 4 ) -> IndelTransitionProbabilities
}

impl PaPaPredIndel {
    /// Create a PaPaPredIndel instance.
    ///
    /// Note that this class treats insertions and deletions the same, hence "indel".
    /// But it distinguishes between outframe (=frameshift) and inframe (non-frameshift) mutations.
    ///
    /// The input file must consist of 3 space-separated columns:
    /// 1) The pattern in UIPAC notation.
    /// 2) The probability for a frameshift indel
    /// 3) The probability for a non-frameshift indel
    ///
    /// You provide the `path` to the sequence-context point mutation probabilities file.
    /// If you want to ensure that the patters have a minimum size, you can provide a value for
    /// `min_kmer_size` which will pad all patterns with Ns at the end. The minimum kmer size has
    /// to be an even number because the insertion or deletion should happen right in the middle
    /// of the sequence
    pub fn new<P: AsRef<Path>>(path: P, min_kmer_size: Option<usize>) -> Result<Self, Error> {
        let mut kmer_size = 0; // will be detected from input size
        let template: Vec<IndelTransitionProbabilities> = vec![IndelTransitionProbabilities::default()]; // missing values should be treated as NaNs
        let mut lookup: Vec<IndelTransitionProbabilities> = Vec::new();
        for record_result in Tabfile::open(path)?.separator(' ') {
            let record = record_result?;
            let tokens = record.fields();
            if tokens[0].starts_with('#') || tokens.len() == 0 {
                continue // ignore header or comments
            }

            let uipac_context = pad_size(tokens[0], min_kmer_size.unwrap_or(1));
            let frameshift_probability = tokens[1].parse::<Float>().unwrap();
            let inframe_probability = tokens[2].parse::<Float>().unwrap();

            if kmer_size == 0 { // first line of input
                kmer_size = uipac_context.len();
                if kmer_size % 2 == 1 {
                    return Err(Error::BadKmerSize(kmer_size));
                }
                lookup = template.repeat(4usize.pow(kmer_size as u32));
            } else if uipac_context.len() != kmer_size {
                let message = format!(
                    "The following line does not have length {}: {}",
                    kmer_size,
                    record.line()
                );
                return Err(Error::FileFormat { message });
            }

            let rate = IndelTransitionProbabilities{
                outframe: frameshift_probability,
                inframe: inframe_probability,
            };

            for context in expand_uipac(&uipac_context) {
                let context_base4 = seq2base_four(&context)?;
                lookup[context_base4] = rate;
                let complement_context = base4_reverse_complement(context_base4, kmer_size);
                lookup[complement_context] = rate;
            }
        }
        let radius = kmer_size / 2;
        Ok(Self{radius, lookup})
    }

    /// Query the mutation rates for a sequence context
    ///
    /// The indel is assumed to happen inbetween the two middle nuleotides
    pub fn rates(&self, seq: &str) -> Result<IndelTransitionProbabilities, Error> {
        let index = seq2base_four(seq)?;
        let rates = self.lookup[index];
        Ok(rates)
    }

    pub fn kmer_size(&self) -> usize {
        2 * self.radius
    }

    pub fn radius(&self) -> usize {
        self.radius
    }
}

#[derive(PartialEq)]
enum Base {
    A,
    C,
    G,
    T,
}

impl Base {
    fn complement(&self) -> Self {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
        }
    }
}

struct Substitution {
    from: Base,
    to: Base,
}

impl Substitution {
    /// There are no transitions from nucleotide X to X
    /// Therefore it is sufficient to store only 3 transitions
    /// But we need to know how to adjust the index inside TransitionProbabilities
    fn transition_probabilities_index(&self) -> usize {
        match self.from {
            Base::A => match self.to {
                Base::A => panic!("Can't transition from A to A"),
                Base::C => 0,
                Base::G => 1,
                Base::T => 2,
            },
            Base::C => match self.to {
                Base::A => 0,
                Base::C => panic!("Can't transition from C to C"),
                Base::G => 1,
                Base::T => 2,
            },
            Base::G => match self.to {
                Base::A => 0,
                Base::C => 1,
                Base::G => panic!("Can't transition from G to G"),
                Base::T => 2,
            },
            Base::T => match self.to {
                Base::A => 0,
                Base::C => 1,
                Base::G => 2,
                Base::T => panic!("Can't transition from T to T"),
            },
        }
    }

    fn complement(&self) -> Self {
        Substitution {
            from: self.from.complement(),
            to: self.to.complement(),
        }
    }
}

impl std::str::FromStr for Substitution {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() == 4 {
            if &s[1..3] != "->" {
                Err(Error::FileFormat {
                    message: "expected '->' between nucleotides".to_string(),
                })
            } else {
                let from = s.chars().next().expect("not empty").try_into()?;
                let to = s.chars().nth(3).expect("4 chars").try_into()?;
                Ok(Substitution { from, to })
            }
        } else {
            Err( Error::FileFormat{ message:
                format!( "The string {} does not confirm to the scheme X->Y (where X and Y are nucleotides", s ) } )
        }
    }
}

impl std::convert::TryFrom<char> for Base {
    type Error = Error;
    fn try_from(value: char) -> Result<Self, Self::Error> {
        Ok(match value {
            'A' => Base::A,
            'C' => Base::C,
            'G' => Base::G,
            'T' => Base::T,
            _ => return Err(Error::BadNucleotide(value)),
        })
    }
}

fn seq2base_four(seq: &str) -> Result<usize, Error> {
    let mut result = 0;
    for c in seq.chars() {
        let value = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => return Err(Error::BadNucleotide(c))
        };
        result <<= 2;
        result += value;
    }
    Ok(result)
}

fn expand_uipac(seq: &str) -> Vec<String> {
    // this is fairly expensive. Can I speed this up?
    if !seq.is_empty() {
        let (c, rest) = seq.split_at(1);
        let suffixes = expand_uipac(rest);

        let mut result = Vec::new();
        let expansions = match c.chars().next().unwrap() {
            'A' => "A",
            'C' => "C",
            'G' => "G",
            'T' => "T",
            'R' => "AG",
            'Y' => "CT",
            'S' => "GC",
            'W' => "AT",
            'K' => "GT",
            'M' => "AC",
            'B' => "CGT",
            'D' => "AGT",
            'H' => "ACT",
            'V' => "ACG",
            'N' => "ACGT",
            _ => panic!("Invalid code: {}", c),
        };
        for expansion in expansions.chars() {
            if !suffixes.is_empty() {
                for suffix in &suffixes {
                    let mut s = String::with_capacity(suffix.len() + 1);
                    s.push(expansion);
                    s.push_str(&suffix);
                    result.push(s);
                }
            } else {
                result.push(expansion.to_string());
            }
        }
        result
    } else {
        Vec::new()
    }
}

fn pad_size(seq: &str, min_size: usize) -> Cow<str> {
    let len = seq.len();
    if len >= min_size {
        Cow::Borrowed(seq)
    } else {
        let flanking = "N".repeat((min_size - len) / 2);
        Cow::Owned(format!("{}{}{}", flanking, seq, flanking))
    }
}

// based on the code in seq2base_four, the following conventions apply:
// A -> 0
// C -> 1
// G -> 2
// T -> 3
//
// And therefore the following complements apply:
// 0 -> A => T -> 3
// 1 -> C => G -> 2
// 2 -> G => C -> 1
// 3 -> T => A -> 0
//
// So, basically: complement(x) -> 3 - x
fn base4_reverse_complement(mut value: usize, digits: usize) -> usize {
    let mut result = 0;
    for _ in 0..digits {
        result <<= 2; // multiply result by 4
        result += 3 - (value & 3); // bitmask the last two bits and complement
        value >>= 2; // shift the last two bits away
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq2base_four() {
        let fun = seq2base_four;
        assert_eq!(fun("ACGT").unwrap(), 0 * 64 + 1 * 16 + 2 * 4 + 3 * 1);
        assert_eq!(fun("TGCA").unwrap(), 3 * 64 + 2 * 16 + 1 * 4 + 0 * 1);
        assert_eq!(fun("").unwrap(), 0);
        assert_eq!(fun("TTTG").unwrap(), 3 * 64 + 3 * 16 + 3 * 4 + 2);
        assert_eq!(fun("AAAAAAAA").unwrap(), 0);
        assert_eq!(fun("CCCCC").unwrap(), 256 + 64 + 16 + 4 + 1);
        assert_eq!(
            fun("ACCGCCT").unwrap(),
            1 * 1024 + 1 * 256 + 2 * 64 + 1 * 16 + 1 * 4 + 3
        );
        assert_eq!(fun("e").unwrap_err(), Error::BadNucleotide('e'));
        assert_eq!(fun("ACGTmACGT").unwrap_err(), Error::BadNucleotide('m'));
    }

    #[test]
    fn test_papapred() {
        let full_input_file = "test_assets/PaPa_rates.txt";
        let papa = PaPaPred::new(full_input_file, None).unwrap();

        let rates = papa.rates("ACCGCCT").unwrap();
        assert_eq!(rates.a, 6.0526200e-04);
        assert_eq!(rates.c, 2.3540691e-05);
        assert!(rates.g.is_nan());
        assert_eq!(rates.t, 3.8108174e-05);

        let rates = papa.rates("AGGCGGT").unwrap();
        assert_eq!(rates.a, 3.8108174e-05);
        assert!(rates.c.is_nan());
        assert_eq!(rates.g, 2.3540691e-05);
        assert_eq!(rates.t, 6.0526200e-04);

        let rates = papa.rates("AAACAAA").unwrap();
        assert_eq!(rates.a, 1.94144068e-05);
        assert!(rates.c.is_nan());
        assert_eq!(rates.g, 1.33842095e-05);
        assert_eq!(rates.t, 3.95017705e-05);

        let rates = papa.rates("AAACAAA").unwrap();
        assert_eq!(rates.a, 1.94144068e-05);
        assert!(rates.c.is_nan());
        assert_eq!(rates.g, 1.33842095e-05);
        assert_eq!(rates.t, 3.95017705e-05);

        let rates = papa.rates("TCATTTT").unwrap();
        assert_eq!(rates.a, 5.1615812e-06);
        assert_eq!(rates.c, 2.6344846e-05);
        assert_eq!(rates.g, 7.5895418e-06);
        assert!(rates.t.is_nan());

        let rates = papa.rates("TAGCTTA").unwrap();
        assert_eq!(rates.a, 1.18828875e-05);
        assert!(rates.c.is_nan());
        assert_eq!(rates.g, 1.92736370e-05);
        assert_eq!(rates.t, 5.04661366e-05);

        let rates = papa.rates("GTCGTGC").unwrap();
        assert_eq!(rates.a, 6.9644750e-04);
        assert_eq!(rates.c, 2.0450439e-05);
        assert!(rates.g.is_nan());
        assert_eq!(rates.t, 1.7836264e-05);

        let rates = papa.rates("GGTACTT").unwrap();
        assert!(rates.a.is_nan());
        assert_eq!(rates.c, 7.8399744e-06);
        assert_eq!(rates.g, 3.8916667e-05);
        assert_eq!(rates.t, 6.6361449e-06);
    }

    #[test]
    fn test_papapred_indel() {
        let full_input_file = "test_assets/PaPa_rates_indel.txt";
        let papa = PaPaPredIndel::new(full_input_file, None).unwrap();

        let rates = papa.rates("ACGTAC").unwrap();
        assert_eq!(rates.outframe, 0.25);
        assert_eq!(rates.inframe, 0.75);

        let rates = papa.rates("CATCAT").unwrap();
        assert_eq!(rates.outframe, 0.125);
        assert_eq!(rates.inframe, 0.5625);

        let rates = papa.rates("TTATGG").unwrap();
        assert_eq!(rates.outframe, 0.5);
        assert_eq!(rates.inframe, 0.625);

        let rates = papa.rates("GGCGTA").unwrap();
        assert_eq!(rates.outframe, 0.3125);
        assert_eq!(rates.inframe, 0.375);

        let rates = papa.rates("TTTTTT").unwrap();
        assert!(rates.outframe.is_nan());
        assert!(rates.inframe.is_nan());

    }

    #[test]
    fn test_expand_uipac() {
        let mut foo = expand_uipac("ACGT");
        assert_eq!(foo.len(), 1);
        assert_eq!(foo.pop().unwrap(), "ACGT");

        foo = expand_uipac("");
        assert_eq!(foo.len(), 0);

        foo = expand_uipac("NN");
        let bar = vec![
            "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC",
            "TG", "TT",
        ];
        assert_eq!(foo, bar);

        foo = expand_uipac("RYS"); //AG CT GC
        let bar = vec!["ACG", "ACC", "ATG", "ATC", "GCG", "GCC", "GTG", "GTC"];
        assert_eq!(foo, bar);
    }
    #[test]
    fn test_base4_reverse_complement() {
        assert_eq!(
            base4_reverse_complement(seq2base_four("AAGG").unwrap(), 4),
            seq2base_four("CCTT").unwrap()
        );
        assert_eq!(
            base4_reverse_complement(seq2base_four("GTACTAG").unwrap(), 7),
            seq2base_four("CTAGTAC").unwrap()
        );
        assert_eq!(
            base4_reverse_complement(seq2base_four("TCGATTGCAT").unwrap(), 10),
            seq2base_four("ATGCAATCGA").unwrap()
        );
        assert_eq!(
            base4_reverse_complement(0, 6), // 6 times A
            seq2base_four("TTTTTT").unwrap()
        );
        assert_eq!(
            base4_reverse_complement(seq2base_four("AG").unwrap(), 8),
            seq2base_four("CTTTTTTT").unwrap()
        );
        assert_eq!(
            base4_reverse_complement(seq2base_four("ACCGCCT").unwrap(), 7),
            seq2base_four("AGGCGGT").unwrap()
        );
    }

    #[test]
    fn test_pad_size() {
        assert_eq!(pad_size("123", 1), "123");
        assert_eq!(pad_size("123", 3), "123");
        assert_eq!(pad_size("123", 5), "N123N");
        assert_eq!(pad_size("123", 7), "NN123NN");

        assert_eq!(pad_size("1234", 2), "1234");
        assert_eq!(pad_size("1234", 4), "1234");
        assert_eq!(pad_size("1234", 6), "N1234N");
        assert_eq!(pad_size("1234", 8), "NN1234NN");
    }
}
