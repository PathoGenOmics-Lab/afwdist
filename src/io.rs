use std::collections::HashMap;
use std::error::Error;
use std::fs::File;

use bio::io::fasta;
use csv;
use log::*;

use crate::{Allele, Sample};

#[derive(serde::Deserialize, Debug)]
struct InputRow {
    sample: String,
    position: usize,
    sequence: String,
    frequency: f64,
}

#[derive(serde::Serialize)]
struct OutputRow {
    sample_m: String,
    sample_n: String,
    distance: f64,
}

pub(crate) fn read_monofasta(path: &str) -> fasta::Record {
    let reader = fasta::Reader::from_file(path).expect("Error reading monoFASTA file");
    reader
        .records()
        .next()
        .unwrap()
        .expect("Error reading first record in monoFASTA")
}

pub(crate) fn read_input_table(reader: File) -> Result<Vec<Sample>, Box<dyn Error>> {
    debug!("Reading input rows");
    let mut rdr = csv::Reader::from_reader(reader);
    let headers = rdr.headers()?.clone();
    let results: Result<Vec<_>, csv::Error> = rdr.records().collect();
    let records = results?;
    let mut samples: HashMap<String, Sample> = HashMap::new();
    for record in records {
        let row: InputRow = record.deserialize(Some(&headers))?;
        let variant = Allele {
            pos: row.position,
            seq: row.sequence.into_bytes(),
            freq: row.frequency,
        };
        if let Some(sample) = samples.get_mut(&row.sample) {
            sample.variants.push(variant);
        } else {
            samples.insert(
                row.sample.clone(),
                Sample {
                    name: row.sample,
                    variants: vec![variant],
                },
            );
        }
    }
    Ok(samples.into_values().collect())
}

pub(crate) fn write_output_table<I: Iterator<Item = (String, String, f64)>>(
    rows: I,
    path: String,
) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::Writer::from_path(path)?;
    for (sample_m, sample_n, distance) in rows {
        wtr.serialize(OutputRow {
            sample_m,
            sample_n,
            distance,
        })?;
    }
    Ok(())
}
