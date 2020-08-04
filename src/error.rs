use std::convert::From;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum Error {
    IO(io::Error),
    FileFormat { message: String },
    BadNucleotide(char),
    BadKmerSize(usize),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::IO(e) => e.fmt(f),
            Error::FileFormat { message } => write!(f, "{}", message),
            Error::BadNucleotide(nuc) => write!(f, "Bad nucleotide: {}", nuc),
            Error::BadKmerSize(k) => write!(f, "Bad kmer size: {}", k),
        }
    }
}

impl std::cmp::PartialEq for Error {
    fn eq(&self, other: &Self) -> bool {
        match self {
            Self::IO(_io) => false, // can't compare. That's why we have to implement this by hand
            Self::FileFormat{message} => match other {
                Self::FileFormat{message: m2} => message == m2,
                _ => false
            },
            Self::BadNucleotide(c1) => match other {
                Self::BadNucleotide(c2) => c1 == c2,
                _ => false
            },
            Self::BadKmerSize(k1) => match other {
                Self::BadKmerSize(k2) => k1 == k2,
                _ => false,
            }
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(other: io::Error) -> Self {
        Error::IO(other)
    }
}
