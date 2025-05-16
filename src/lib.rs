use std::collections::HashSet;
use std::fmt;
use std::{error::Error, fs::File};

use clap::Parser;
use itertools::Itertools;
use log::*;
use simplelog::*;

mod io;

/// afwdist (allele frequency-weighted distances)
#[derive(Parser, Debug)]
#[command(author = "Miguel Álvarez Herrera <miguel.alvarez@uv.es>")]
#[command(version)]
struct Args {
    /// Input tree in CSV format (mandatory CSV columns are 'sample', 'position', 'sequence' and 'frequency')
    #[arg(short = 'i', long, required = true)]
    input: String,

    /// Reference sequence in FASTA format
    #[arg(short = 'r', long, required = true)]
    reference: String,

    /// Output CSV file with distances between each pair of samples
    #[arg(short = 'o', long, required = true)]
    output: String,

    /// Enable debug messages
    #[arg(short = 'v', long, default_value_t = false)]
    verbose: bool,
}

#[derive(Clone)]
struct Allele {
    pos: usize,
    seq: Vec<u8>,
    freq: f64,
}

impl fmt::Debug for Allele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}:{:.2}",
            self.pos,
            String::from_utf8_lossy(&self.seq),
            self.freq
        )
    }
}

struct Sample {
    name: String,
    variants: Vec<Allele>,
}

impl Sample {
    fn positions(&self) -> HashSet<usize> {
        self.variants.iter().map(|v| v.pos).collect()
    }

    fn alleles_at(&self, site: usize) -> Vec<&Allele> {
        self.variants.iter().filter(|v| v.pos == site).collect()
    }

    fn alleles_incl_ref_at(&self, site: usize, reference: &[u8]) -> Vec<Allele> {
        let mut alleles: Vec<Allele> = self.alleles_at(site).into_iter().cloned().collect();
        let variant_frequency: f64 = alleles.iter().map(|a| a.freq).sum();
        debug!(
            "Sample {:?} - total variant frequency for site {site} = {variant_frequency}",
            self.name
        );
        let reference_allele = Allele {
            pos: site,
            seq: reference[site - 1..site].to_vec(),
            freq: 1.0 - variant_frequency,
        };
        debug!(
            "Sample {:?} - setting reference allele {:?} at site {} to {}",
            self.name,
            String::from_utf8_lossy(&reference_allele.seq),
            site,
            reference_allele.freq
        );
        alleles.push(reference_allele);
        debug!(
            "Sample {:?} - alleles at site {site}: {alleles:?}",
            self.name
        );
        alleles
    }

    fn polymorphic_sites(&self, sample: &Sample) -> Vec<usize> {
        self.positions()
            .union(&sample.positions())
            .cloned()
            .collect()
    }
}

fn sq_dif(freq_m: f64, freq_n: f64) -> f64 {
    f64::powf(freq_m - freq_n, 2.0)
}

fn sq_sum(freq_m: f64, freq_n: f64) -> f64 {
    f64::powf(freq_m + freq_n, 2.0)
}

fn sum_func(alleles_m: &Vec<Allele>, alleles_n: &Vec<Allele>, func: fn(f64, f64) -> f64) -> f64 {
    alleles_m
        .iter()
        .cartesian_product(alleles_n)
        .filter(|(a_m, a_n)| a_m.seq == a_n.seq)
        .map(|(a_m, a_n)| {
            debug!("Operating on frequencies of alleles {a_m:?} and {a_n:?}");
            func(a_m.freq, a_n.freq)
        })
        .sum()
}

fn complete_frequencies(alleles: &Vec<Allele>, with_alleles: &Vec<Allele>) -> Vec<Allele> {
    let mut completed_alleles = alleles.clone();
    for allele in with_alleles {
        if alleles.iter().all(|a| a.seq != allele.seq) {
            completed_alleles.push(Allele {
                pos: allele.pos,
                seq: allele.seq.clone(),
                freq: 0.0,
            });
        }
    }
    debug!("Completed alleles: {completed_alleles:?}");
    completed_alleles
}

fn afwdist(m: &Sample, n: &Sample, reference: &[u8]) -> f64 {
    m.polymorphic_sites(n)
        .iter()
        .map(|&site| {
            debug!("Processing site {site}");
            // Get alleles including the reference one
            let alleles_m = m.alleles_incl_ref_at(site, reference);
            let alleles_n = n.alleles_incl_ref_at(site, reference);
            // Complete with each other's potentially missing alleles
            let complete_alleles_m = complete_frequencies(&alleles_m, &alleles_n);
            let complete_alleles_n = complete_frequencies(&alleles_n, &alleles_m);
            // Calculate distances
            let denominator = 4.0 - sum_func(&complete_alleles_m, &complete_alleles_n, sq_sum);
            if denominator < 0.0 {
                warn!(
                    "Term denominator is negative for site {site} of samples {:?} and {:?}",
                    m.name, n.name
                );
            }
            if denominator != 0.0 {
                sum_func(&complete_alleles_m, &complete_alleles_n, sq_dif) / denominator
            } else {
                0.0
            }
        })
        .sum()
}

pub fn run() -> Result<(), Box<dyn Error>> {
    // Arguments
    let args = Args::parse();
    // Set log level
    if args.verbose {
        let _ = SimpleLogger::init(LevelFilter::Debug, Config::default());
    } else {
        let _ = SimpleLogger::init(LevelFilter::Info, Config::default());
    }
    // Read input table
    info!("Reading input table");
    let input_table = File::open(args.input)?;
    let samples = io::read_input_table(input_table)?;
    // Read reference sequence
    info!("Reading reference sequence");
    let reference_record = io::read_monofasta(&args.reference);
    let reference = reference_record.seq();
    // Calculate distances and write results
    let distances = samples.iter().combinations(2).map(|pair| {
        let m = pair[0];
        let n = pair[1];
        debug!(
            "Calculating distance between sample {:?} and {:?}",
            m.name, n.name
        );
        (m.name.clone(), n.name.clone(), afwdist(&m, &n, reference))
    });
    info!("Calculating distances and writing results");
    io::write_output_table(distances, args.output)?;
    info!("Done");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn same_sample() {
        let reference = b"TTT";
        let m = Sample {
            name: "test".to_owned(),
            variants: vec![
                Allele {
                    pos: 1,
                    seq: b"A".to_vec(),
                    freq: 0.5,
                },
                Allele {
                    pos: 1,
                    seq: b"C".to_vec(),
                    freq: 0.1,
                },
                Allele {
                    pos: 2,
                    seq: b"A".to_vec(),
                    freq: 0.9,
                },
                Allele {
                    pos: 3,
                    seq: b"A".to_vec(),
                    freq: 0.2,
                },
            ],
        };
        let d = afwdist(&m, &m, reference);
        assert_approx_eq!(d, 0.0);
    }

    #[test]
    fn symmetry() {
        let reference = b"TTT";
        let m = Sample {
            name: "test".to_owned(),
            variants: vec![
                Allele {
                    pos: 1,
                    seq: b"A".to_vec(),
                    freq: 0.5,
                },
                Allele {
                    pos: 1,
                    seq: b"C".to_vec(),
                    freq: 0.1,
                },
                Allele {
                    pos: 2,
                    seq: b"A".to_vec(),
                    freq: 0.9,
                },
                Allele {
                    pos: 3,
                    seq: b"A".to_vec(),
                    freq: 0.2,
                },
            ],
        };
        let n = Sample {
            name: "test".to_owned(),
            variants: vec![
                Allele {
                    pos: 1,
                    seq: b"A".to_vec(),
                    freq: 0.25,
                },
                Allele {
                    pos: 2,
                    seq: b"C".to_vec(),
                    freq: 0.3,
                },
                Allele {
                    pos: 2,
                    seq: b"C".to_vec(),
                    freq: 0.4,
                },
                Allele {
                    pos: 3,
                    seq: b"A".to_vec(),
                    freq: 0.9,
                },
            ],
        };
        let d_m_n = afwdist(&m, &n, reference);
        let d_n_m = afwdist(&n, &m, reference);
        assert_approx_eq!(d_m_n, d_n_m);
    }
}
