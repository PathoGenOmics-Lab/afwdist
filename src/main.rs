fn main() {
    if let Err(e) = afwdist::run() {
        println!("Application error: {e}");
        std::process::exit(1);
    }
}
