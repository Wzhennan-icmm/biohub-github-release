use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

struct CaseDir(PathBuf);

impl CaseDir {
    fn new(label: &str) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "biohub_publication_{label}_{}_{}",
            std::process::id(),
            nonce
        ));
        fs::create_dir_all(&path).expect("create case directory");
        Self(path)
    }

    fn path(&self, relative: &str) -> PathBuf {
        self.0.join(relative)
    }

    fn write(&self, relative: &str, content: &str) -> PathBuf {
        let path = self.path(relative);
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).expect("create fixture parent");
        }
        fs::write(&path, content).expect("write fixture");
        path
    }
}

impl Drop for CaseDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn run(case: &CaseDir, arguments: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_biohub"))
        .args(arguments)
        .current_dir(&case.0)
        .output()
        .expect("run BioHub")
}

fn run_ok(case: &CaseDir, arguments: &[String]) -> Output {
    let result = run(case, arguments);
    assert!(
        result.status.success(),
        "command failed: {}\nstdout:\n{}\nstderr:\n{}",
        arguments.join(" "),
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    result
}

fn s(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn assert_file(path: &Path, expected: &str) {
    let actual = fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("read expected output {}: {error}", path.display()));
    assert_eq!(actual, expected, "golden mismatch: {}", path.display());
}

fn blast_row(query: &str, target: &str, identity: f64, evalue: &str, score: f64) -> String {
    format!("{query}\t{target}\t{identity}\t0\t0\t0\t0\t0\t0\t0\t{evalue}\t{score}\n")
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[test]
fn text_fasta_and_filter_commands_match_goldens() {
    let case = CaseDir::new("text_fasta");

    let fasta = case.write(
        "rename/input.fa",
        ">old1 description\nAA\nAA\n>old2\nCCCC\n",
    );
    let names = case.write("rename/names.tsv", "new1 old1\nnew3 missing\n");
    let renamed = case.path("rename/output.fa");
    run_ok(
        &case,
        &[
            "run".into(),
            "change-scaffolds-name-fasta".into(),
            "-i".into(),
            s(&fasta),
            "-l".into(),
            s(&names),
            "-o".into(),
            s(&renamed),
        ],
    );
    assert_file(&renamed, ">new1\nAAAA\n");

    let alignment = case.write("rename/aln/group.aln", ">old1\nAAAA\n>old2\nCCCC\n");
    assert!(alignment.is_file());
    let seq_names = case.write("rename/aln-names.tsv", "sampleA old1\nsampleB old2\n");
    let renamed_dir = case.path("rename/renamed");
    run_ok(
        &case,
        &[
            "run".into(),
            "change-seqname-for-fasta".into(),
            "-i".into(),
            s(&case.path("rename/aln")),
            "-l".into(),
            s(&seq_names),
            "-o".into(),
            s(&renamed_dir),
        ],
    );
    assert_file(
        &renamed_dir.join("group.newName.fa"),
        ">sampleA\nAAAA\n>sampleB\nCCCC\n",
    );

    let depths = case.write(
        "depth/input.tsv",
        "chr1 1 10 5\nchr1 1 10 6\nchr1 1 10 7\nchr2 11 20 8\nchr2 11 20 9\nchr2 11 20 10\n",
    );
    let depth_output = case.path("depth/output.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "convert-3line2one".into(),
            "-i".into(),
            s(&depths),
            "-o".into(),
            s(&depth_output),
        ],
    );
    assert_file(
        &depth_output,
        "chr1\t1\t10\t5\t6\t7\nchr2\t11\t20\t8\t9\t10\n",
    );

    let sequences = case.write(
        "filter/input.fa",
        ">short\nAAA\n>edge\nAAAA\n>long\nAAAAAA\n",
    );
    let filtered_fasta = case.path("filter/output.fa");
    run_ok(
        &case,
        &[
            "run".into(),
            "filter-seq-by-length".into(),
            "-i".into(),
            s(&sequences),
            "-l".into(),
            "4".into(),
            "-o".into(),
            s(&filtered_fasta),
        ],
    );
    assert_file(&filtered_fasta, ">edge\nAAAA\n>long\nAAAAAA\n");

    let gff = case.write(
        "filter/input.gff3",
        "##gff-version 3\nchr1\ttest\tgene\t1\t5\t.\t+\t.\tID=keep;Name=A\nchr1\ttest\tgene\t6\t9\t.\t+\t.\tID=drop;Name=B\n",
    );
    let ids = case.write("filter/ids.txt", "ID=keep\n");
    let filtered_gff = case.path("filter/output.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "filter-gff-by-id".into(),
            "-gff".into(),
            s(&gff),
            "-id".into(),
            s(&ids),
            "-o".into(),
            s(&filtered_gff),
        ],
    );
    assert_file(
        &filtered_gff,
        "##gff-version 3\nchr1\ttest\tgene\t1\t5\t.\t+\t.\tID=keep;Name=A\n",
    );

    let gtf = case.write(
        "filter/input.gtf",
        "chr1 test exon 1 2 . + . gene_id g1\nctg9 test exon 1 2 . + . gene_id g2\n",
    );
    let excluded = case.write("filter/excluded.txt", "ctg9\n");
    let filtered_gtf = case.path("filter/output.gtf");
    run_ok(
        &case,
        &[
            "run".into(),
            "filter-gtf-ctg".into(),
            "-i".into(),
            s(&gtf),
            "-id".into(),
            s(&excluded),
            "-o".into(),
            s(&filtered_gtf),
        ],
    );
    assert_file(&filtered_gtf, "chr1 test exon 1 2 . + . gene_id g1\n");

    let duplicates = case.write(
        "table/duplicates.tsv",
        "geneB x 1\ngeneA old 2\ngeneA best 9\ngeneB best 3\n",
    );
    let deduped = case.path("table/deduped.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "merge-two-txt".into(),
            "-i".into(),
            s(&duplicates),
            "-o".into(),
            s(&deduped),
        ],
    );
    assert_file(&deduped, "geneA\tbest\t9\ngeneB\tbest\t3\n");

    let groups = case.write("table/groups.txt", "G2:x,x,y,\nG1:a,a,\nG3:z,y,z,\n");
    let unique_groups = case.path("table/groups.unique.txt");
    run_ok(
        &case,
        &[
            "run".into(),
            "orthogenes".into(),
            "-i".into(),
            s(&groups),
            "-o".into(),
            s(&unique_groups),
        ],
    );
    assert_file(&unique_groups, "G2:x,y\nG3:z,y\n");

    let genome = case.write("stats/genome.fa", ">chr1\nACGTNN\n>chr2\nggcc\n");
    let gc = case.path("stats/gc.txt");
    run_ok(
        &case,
        &[
            "run".into(),
            "genome-gc".into(),
            "-f".into(),
            s(&genome),
            "-o".into(),
            s(&gc),
        ],
    );
    assert_file(&gc, "0.75\n");

    let transcripts = case.write(
        "fasta/transcripts.fa",
        ">tx2.geneA\nAAAAAA\n>tx1.geneA\nAAA\n>tx3.geneB\nCCCC\n",
    );
    let longest = case.path("fasta/longest.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-the-longest-seq".into(),
            "-i".into(),
            s(&transcripts),
            "-o".into(),
            s(&longest),
        ],
    );
    assert_file(&longest, "tx2\tgeneA\ntx3\tgeneB\n");

    let peptides = case.write(
        "fasta/peptides.fa",
        ">geneA-copy1\nAAA\n>geneB-copy1\nCCCC\n>geneA-copy2\nAAAAAA\n",
    );
    let longest_peptide = case.path("fasta/longest-peptide.fa");
    run_ok(
        &case,
        &[
            "run".into(),
            "extract-longest-pep".into(),
            "-f".into(),
            s(&peptides),
            "-o".into(),
            s(&longest_peptide),
        ],
    );
    assert_file(&longest_peptide, ">geneA\nAAAAAA\n>geneB\nCCCC\n");

    let fastq = case.write(
        "fastq/input.fastq",
        "@head\nTTAGGGAA\n+\nIIIIIIII\n@tail\nAATTAGGG\n+\nIIIIIIII\n@keep\nAACCGGTT\n+\nIIIIIIII\n",
    );
    let clean_fastq = case.path("fastq/clean.fastq");
    run_ok(
        &case,
        &[
            "run".into(),
            "trim-ttaggg-fastq".into(),
            "-i".into(),
            s(&fastq),
            "-o".into(),
            s(&clean_fastq),
        ],
    );
    assert_file(&clean_fastq, "@keep\nAACCGGTT\n+\nIIIIIIII\n");
}

#[test]
fn table_expression_and_blast_commands_match_goldens() {
    let case = CaseDir::new("tables_blast");

    let no_as = case.write("as/no-as.tsv", "geneA x present\ngeneB x present\n");
    let as_table = case.write(
        "as/as.tsv",
        "# header\ngeneA isoform present\ngeneC isoform present\n",
    );
    let as_output = case.path("as/matched.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "compare-as-and-noAS".into(),
            "-nA".into(),
            s(&no_as),
            "-AS".into(),
            s(&as_table),
            "-o".into(),
            s(&as_output),
        ],
    );
    assert_file(&as_output, "geneA isoform present\n");

    let ancestor = case.write("busco/ancestor.txt", "BUSCO_A\nBUSCO_B\nBUSCO_B\n");
    let offspring = case.write("busco/offspring.txt", "BUSCO_B\nBUSCO_C\n");
    run_ok(
        &case,
        &[
            "run".into(),
            "compare-busco-results".into(),
            "-a".into(),
            s(&ancestor),
            "-o".into(),
            s(&offspring),
        ],
    );
    assert_file(&case.path("justHave.txt"), "BUSCO_C\n");
    assert_file(&case.path("justLost.txt"), "BUSCO_A\n");

    let primary = case.write("blast/forward.tsv", "qB\ttB\nqA\ttA\nqC\ttC\n");
    let reverse = case.write("blast/reverse.tsv", "tB\tother\ntA\tqA\ntC\tqC\n");
    let reciprocal = case.path("blast/reciprocal.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "compare-two-blast".into(),
            "-i".into(),
            s(&primary),
            "-r".into(),
            s(&reverse),
            "-o".into(),
            s(&reciprocal),
        ],
    );
    assert_file(&reciprocal, "qA\ttA\nqC\ttC\n");

    let blast = case.write(
        "blast/hits.tsv",
        &(blast_row("qB", "t3", 80.0, "1e-2", 10.0)
            + &blast_row("qA", "t1", 70.0, "1e-3", 20.0)
            + &blast_row("qA", "t2", 90.0, "1e-20", 50.0)),
    );
    let best_identity = case.path("blast/best-identity.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-idy".into(),
            "-i".into(),
            s(&blast),
            "-o".into(),
            s(&best_identity),
        ],
    );
    assert_file(&best_identity, "qA\t90\nqB\t80\n");

    let best_score_rows = case.path("blast/best-score-rows.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-hit-based-on-idy".into(),
            "-i".into(),
            s(&blast),
            "-o".into(),
            s(&best_score_rows),
        ],
    );
    assert_file(
        &best_score_rows,
        &(blast_row("qA", "t2", 90.0, "1e-20", 50.0) + &blast_row("qB", "t3", 80.0, "1e-2", 10.0)),
    );

    let query_hits = case.write(
        "blast/query.tsv",
        &(blast_row("qA", "tA", 90.0, "1e-20", 60.0)
            + &blast_row("qA", "tX", 91.0, "1e-10", 40.0)
            + &blast_row("qB", "tB", 88.0, "1e-15", 50.0)),
    );
    let reference_hits = case.write(
        "blast/reference.tsv",
        &(blast_row("tA", "qA", 90.0, "1e-20", 70.0) + &blast_row("tB", "qB", 88.0, "1e-15", 55.0)),
    );
    let best_pairs = case.path("blast/best-pairs.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-hit-by-score".into(),
            "-i".into(),
            s(&query_hits),
            "-r".into(),
            s(&reference_hits),
            "-o".into(),
            s(&best_pairs),
        ],
    );
    assert_file(&best_pairs, "qA\ttA\nqB\ttB\n");

    let one_file = case.write(
        "blast/one-file.tsv",
        &(blast_row("qA", "tA", 90.0, "1e-20", 60.0)
            + &blast_row("qA", "tB", 70.0, "1e-2", 10.0)
            + &blast_row("qB", "tB", 88.0, "1e-15", 50.0)),
    );
    let one_file_output = case.path("blast/one-file.idy.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-hit-by-score-one-file".into(),
            "-i".into(),
            s(&one_file),
            "-p".into(),
            s(&case.path("blast/prefix")),
            "-o".into(),
            s(&one_file_output),
        ],
    );
    assert_file(&one_file_output, "qA\ttA\t90\nqB\ttB\t88\n");

    let go_input = case.write(
        "go/input.tsv",
        "geneA description GO:0001 GO:0002\ngeneB none\n",
    );
    let go_output = case.path("go/extracted.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "save-go".into(),
            "-i".into(),
            s(&go_input),
            "-o".into(),
            s(&go_output),
        ],
    );
    assert_file(&go_output, "geneA\tGO:0001\tGO:0002\ngeneB\t\n");

    let swiss = case.write("go/swiss.tsv", "geneB GO:0003\ngeneA GO:0001\n");
    let nr = case.write("go/nr.tsv", "geneA GO:0002\ngeneB EC:1.2.3\n");
    let trembl = case.write("go/trembl.tsv", "geneC GO:0004\n");
    let merged_go = case.path("go/merged.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "merge-gos".into(),
            "-s".into(),
            s(&swiss),
            "-n".into(),
            s(&nr),
            "-T".into(),
            s(&trembl),
            "-o".into(),
            s(&merged_go),
        ],
    );
    assert_file(
        &merged_go,
        "geneA\tGO:0001\ngeneA\tGO:0002\ngeneB\tGO:0003\ngeneC\tGO:0004\n",
    );

    case.write("jcvi/input/a.tsv", "orthoA gene1\northoB gene2\n");
    case.write("jcvi/input/b.tsv", "orthoC gene1\northoD gene3\n");
    let jcvi_output = case.path("jcvi/matrix.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "merge-blastp-best-jcvi".into(),
            "-p".into(),
            s(&case.path("jcvi/input")),
            "-o".into(),
            s(&jcvi_output),
        ],
    );
    assert_file(
        &jcvi_output,
        "gene1\torthoA\torthoC\ngene2\torthoB\t.\ngene3\t.\torthoD\n",
    );

    let first_xls = case.write("xls/first.tsv", "prefix:geneA.1 metaA valueA\n");
    let second_xls = case.write("xls/second.tsv", "geneA.9 payload\ngeneB.1 other\n");
    let merged_xls = case.path("xls/merged.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "zhouxiaoxuan-mergexls".into(),
            "-a".into(),
            s(&first_xls),
            "-b".into(),
            s(&second_xls),
            "-o".into(),
            s(&merged_xls),
        ],
    );
    assert_file(
        &merged_xls,
        "geneA.9 payload\tprefix:geneA.1\tmetaA\tvalueA\ngeneB.1 other\n",
    );

    let fpkm_a = case.write("fpkm/input/a.tsv", "geneB 2\ngeneA 1\n");
    let fpkm_b = case.write("fpkm/input/b.tsv", "geneA 3\ngeneB 4\n");
    let fpkm = case.path("fpkm/matrix.tsv");
    let profile = case.path("fpkm/profile.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "merge-fpkm-file".into(),
            "-i".into(),
            s(&case.path("fpkm/input")),
            "-oF".into(),
            s(&fpkm),
            "-oP".into(),
            s(&profile),
        ],
    );
    assert_file(
        &fpkm,
        &format!(
            "gene_id\t{}\t{}\ngeneA\t1\t3\ngeneB\t2\t4\n",
            fpkm_a.display(),
            fpkm_b.display()
        ),
    );
    assert_file(&profile, "geneA\t4\ngeneB\t6\n");
}

#[test]
fn annotation_expression_and_orthology_commands_match_goldens() {
    let case = CaseDir::new("domain_goldens");

    let pasa = case.write(
        "pasa/input.gff3",
        "#PROT tx1 gene1 MAA\nchr1\tPASA\tgene\t1\t100\t.\t+\t.\tID=gene1\nchr1\tPASA\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\nchr1\tPASA\texon\t1\t50\t.\t+\t.\tParent=tx1\nchr1\tPASA\tCDS\t1\t45\t.\t+\t0\tParent=tx1\nchr1\tPASA\tmRNA\t1\t80\t.\t+\t.\tID=tx2;Parent=gene1\nchr1\tPASA\texon\t1\t80\t.\t+\t.\tParent=tx2\n",
    );
    let pasa_pep = case.path("pasa/proteins.fa");
    let pasa_gff = case.path("pasa/longest.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "extract-pasa-results".into(),
            "-i".into(),
            s(&pasa),
            "-s".into(),
            s(&pasa_pep),
            "-g".into(),
            s(&pasa_gff),
        ],
    );
    assert_file(&pasa_pep, ">tx1-gene1\nMAA\n");
    assert_file(
        &pasa_gff,
        "chr1\tPASA\tgene\t1\t100\t.\t+\t.\tID=gene1\nchr1\tPASA\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\nchr1\tPASA\texon\t1\t50\t.\t+\t.\tParent=tx1\nchr1\tPASA\tCDS\t1\t45\t.\t+\t0\tParent=tx1\n",
    );

    let gemoma = case.write(
        "gff/gemoma.gff3",
        "chr1\tGeMoMa\tgene\t1\t100\t.\t+\t.\tID=old_gene\nchr1\tGeMoMa\tmRNA\t1\t100\t.\t+\t.\tID=old_tx\nchr1\tGeMoMa\tCDS\t1\t30\t.\t+\t0\tParent=old_tx\n",
    );
    let converted_gemoma = case.path("gff/gemoma.converted.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "convert-gemoma-gff3".into(),
            "-i".into(),
            s(&gemoma),
            "-o".into(),
            s(&converted_gemoma),
        ],
    );
    assert_file(
        &converted_gemoma,
        "chr1\tGeMoMa\tgene\t1\t100\t.\t+\t.\tID=PlantsCHR1gene00001;Name=Plant00001\nchr1\tGeMoMa\tmRNA\t1\t100\t.\t+\t.\tID=PlantsCHR1gene00001.1;Parent=PlantsCHR1gene00001;Name=Plant00001.1\nchr1\tGeMoMa\texon\t1\t30\t.\t+\t0\tID=PlantsCHR1gene00001.1.exon1;Parent=PlantsCHR1gene00001.1\nchr1\tGeMoMa\tCDS\t1\t30\t.\t+\t0\tID=cds.PlantsCHR1gene00001.1;Parent=PlantsCHR1gene00001.1\n",
    );

    let mapped_gff = case.write(
        "gff/contigs.gff3",
        "ctgP\ttest\tgene\t10\t20\t.\t+\t.\tID=plus\nctgM\ttest\tgene\t10\t20\t.\t+\t.\tID=minus\nctgU\ttest\tgene\t1\t2\t.\t+\t.\tID=unmapped\n",
    );
    let mapping = case.write(
        "gff/contigs.map",
        "chr1 100 unused ctgP 0 1000\nchr2 500 unused ctgM 1 100\n",
    );
    let mapped_output = case.path("gff/contigs.mapped.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "convert-gene-annotation-contigs2chr-PASA".into(),
            "-gff".into(),
            s(&mapped_gff),
            "-b".into(),
            s(&mapping),
            "-o".into(),
            s(&mapped_output),
        ],
    );
    assert_file(
        &mapped_output,
        "chr1\ttest\tgene\t109\t119\t.\t+\t.\tID=plus\nchr2\ttest\tgene\t580\t590\t.\t-\t.\tID=minus\nctgU\ttest\tgene\t1\t2\t.\t+\t.\tID=unmapped\n",
    );

    let isoforms = case.write(
        "gff/isoforms.gff3",
        "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=g1\nchr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=short;Parent=g1\nchr1\ttest\tCDS\t1\t1\t.\t+\t0\tParent=short\nchr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=long;Parent=g1\nchr1\ttest\tCDS\t5\t5\t.\t+\t0\tParent=long\nchr1\ttest\tCDS\t9\t9\t.\t+\t0\tParent=long\n",
    );
    let longest_isoform = case.path("gff/longest-isoform.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "filter-gemoma-as2".into(),
            "-i".into(),
            s(&isoforms),
            "-o".into(),
            s(&longest_isoform),
        ],
    );
    assert_file(
        &longest_isoform,
        "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=g1\nchr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=long;Parent=g1\nchr1\ttest\tCDS\t5\t5\t.\t+\t0\tParent=long\nchr1\ttest\tCDS\t9\t9\t.\t+\t0\tParent=long\n",
    );

    let lengths = case.write(
        "family/lengths-info.tsv",
        "gene_id x length\nFamA.parent x 100\nFamA.a x 90\nFamA.b x 40\n",
    );
    let expression = case.write(
        "family/expression.tsv",
        "gene sample1 sample2\nFamA.a 8 10\nFamA.b 2 4\n",
    );
    let parents = case.write("family/parents.tsv", "FamA x FamA.parent\n");
    let split_dir = case.path("family/split");
    run_ok(
        &case,
        &[
            "run".into(),
            "extract-gene-family-info".into(),
            "-l".into(),
            s(&lengths),
            "-e".into(),
            s(&expression),
            "-p".into(),
            s(&parents),
            "-c".into(),
            "0.8".into(),
            "-o".into(),
            s(&split_dir),
        ],
    );
    assert_file(
        &split_dir.join("FamA.fullenth.geneExpression.txt"),
        "FamA.a\t8\t10\n",
    );
    assert_file(
        &split_dir.join("FamA.partial.geneExpression.txt"),
        "FamA.b\t2\t4\n",
    );

    let matrix_lengths = case.write(
        "family/lengths-matrix.tsv",
        "gene_id length\nFamA.parent 100\nFamA.a 90\nFamA.b 40\n",
    );
    let matrix_expression = case.write(
        "family/matrix-expression.tsv",
        "gene sample1 sample2\nFamA.parent 10 12\nFamA.a 8 10\nFamA.b 2 4\n",
    );
    let family_names = case.write("family/names.txt", "FamA\n");
    let family_members = case.write(
        "family/members.tsv",
        "FamA FamA.parent\nFamA FamA.a\nFamA FamA.b\n",
    );
    let family_matrix = case.path("family/final-matrix.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "extract-gene-family-matrix".into(),
            "-l".into(),
            s(&matrix_lengths),
            "-e".into(),
            s(&matrix_expression),
            "-f".into(),
            s(&family_names),
            "-g".into(),
            s(&family_members),
            "-c".into(),
            "0.8".into(),
            "-o".into(),
            s(&family_matrix),
        ],
    );
    assert_file(&family_matrix, "FamA\t3\t10\n");

    case.write(
        "cross-species/input/pairs.tsv",
        "qA qB\nqB qA\nqC qD\nqD qC\n",
    );
    let cross_species = case.path("cross-species/result.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-hit-from-blast".into(),
            "-i".into(),
            s(&case.path("cross-species/input")),
            "-n".into(),
            "2".into(),
            "-o".into(),
            s(&cross_species),
        ],
    );
    assert_file(&cross_species, "qA\tqB\nqC\tqD\n");

    let bidirectional = case.write(
        "blast/bidirectional.tsv",
        &blast_row("qA", "tA", 90.0, "1e-20", 50.0),
    );
    let bidirectional_output = case.path("blast/bidirectional-best.tsv");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-best-hit-genes".into(),
            "-i".into(),
            s(&bidirectional),
            "-o".into(),
            s(&bidirectional_output),
        ],
    );
    assert_file(
        &bidirectional_output,
        "qA\ttA\t1e-20\t100\t1e-20\t100\tNoop\tNoop\ntA\tqA\t1e-20\t100\t1e-20\t100\tNoop\tNoop\n",
    );

    let codons = case.write(
        "orthology/input/group1.fa",
        ">sample001_gene\nGCTGGT\n>sample002_gene\nGCCGGT\n",
    );
    assert!(codons.is_file());
    let diff_dir = case.path("orthology/diff");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-diff-sites-from-orthology".into(),
            "-i".into(),
            s(&case.path("orthology/input")),
            "-o".into(),
            s(&diff_dir),
        ],
    );
    assert_file(
        &diff_dir.join("fourDegenerateSite.fasta"),
        ">sample001\nTT\n>sample002\nCT\n",
    );
    assert_file(
        &diff_dir.join("fourDegenerateCodenFile.fasta"),
        ">sample001\nGCTGGT\n>sample002\nGCCGGT\n",
    );

    let strict_dir = case.path("orthology/strict");
    run_ok(
        &case,
        &[
            "run".into(),
            "get-four-degenerate-sites".into(),
            "-i".into(),
            s(&case.path("orthology/input")),
            "-o".into(),
            s(&strict_dir),
        ],
    );
    assert_file(
        &strict_dir.join("fourDegenerateSite.fasta"),
        ">sample001\nT\n>sample002\nT\n",
    );
    assert_file(
        &strict_dir.join("fourDegenerateSite.stat"),
        "#totalSites\t1\n#allSites\t2\n",
    );

    case.write("psmc/input/sampleA.0.txt", "100 200\n200 300\n");
    case.write("psmc/input/sampleB.0.txt", "100 400\n");
    let psmc = case.path("psmc/merged.tsv");
    run_ok(
        &case,
        &[
            "psmc".into(),
            "merge".into(),
            "--dir".into(),
            s(&case.path("psmc/input")),
            "--pattern".into(),
            ".0.txt".into(),
            "--output".into(),
            s(&psmc),
        ],
    );
    assert_file(
        &psmc,
        "Sample\tTime\tNe\nsampleA\t100\t200\nsampleA\t200\t300\nsampleB\t100\t400\n",
    );
}

#[test]
fn coordinate_conversion_and_visual_outputs_match_goldens() {
    let case = CaseDir::new("coordinates_visual");

    let gff = case.write(
        "nextgen/input.gff3",
        "ctgP\ttest\tgene\t10\t20\t.\t+\t.\tID=plus\nctgM\ttest\tgene\t10\t20\t.\t+\t.\tID=minus\n",
    );
    let background = case.write(
        "nextgen/background.tsv",
        "chr1\t100\tx\tctgP\tx\t+\t1\t1000\nchr2\t500\tx\tctgM\tx\t-\t1\t100\n",
    );
    let converted = case.path("nextgen/converted.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "convert-gene-annotation-scaffold2chr-nextgenomics".into(),
            "-gff".into(),
            s(&gff),
            "-b".into(),
            s(&background),
            "-o".into(),
            s(&converted),
        ],
    );
    assert_file(
        &converted,
        "chr1\ttest\tgene\t109\t119\t.\t+\t.\tID=plus\nchr2\ttest\tgene\t580\t590\t.\t-\t.\tID=minus\n",
    );

    let isoforms = case.write(
        "gemoma/filter.gff3",
        "##gff-version 3\nchr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\nchr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1\nchr1\ttest\tCDS\t1\t10\t.\t+\t0\tParent=tx1\nchr1\ttest\tmRNA\t1\t40\t.\t+\t.\tID=tx2;Parent=g1\nchr1\ttest\tCDS\t1\t20\t.\t+\t0\tParent=tx2\n",
    );
    let filtered = case.path("gemoma/filtered.gff3");
    run_ok(
        &case,
        &[
            "run".into(),
            "filter-gemoma-as".into(),
            "-i".into(),
            s(&isoforms),
            "-o".into(),
            s(&filtered),
        ],
    );
    assert_file(
        &filtered,
        "##gff-version 3\nchr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\nchr1\ttest\tmRNA\t1\t40\t.\t+\t.\tID=tx2;Parent=g1\nchr1\ttest\tCDS\t1\t20\t.\t+\t0\tParent=tx2\n",
    );

    let pandepth = case.write(
        "plots/pandepth.tsv",
        "#Chr\tStart\tEnd\tMeanDepth\tGC(%)\nchr2\t0\t200\t20\t40\nchr1\t0\t100\t10\t50\nchr1\t100\t200\t30\t60\n",
    );
    let plot_dir = case.path("plots/basic");
    run_ok(
        &case,
        &[
            "run".into(),
            "plot-depth-pandepth".into(),
            "-i".into(),
            s(&pandepth),
            "-o".into(),
            s(&plot_dir),
            "-l".into(),
            "0".into(),
        ],
    );
    assert_file(
        &plot_dir.join("chromosome_stats.tsv"),
        "Chr\tlength\tmean_depth\tmean_gc\nchr1\t200\t20\t55\nchr2\t200\t20\t40\n",
    );
    assert_file(
        &plot_dir.join("filtered_depth.tsv"),
        "Chr\tStart\tEnd\tMeanDepth\tGC(%)\nchr2\t0\t200\t20\t40\nchr1\t0\t100\t10\t50\nchr1\t100\t200\t30\t60\n",
    );
    let basic_svg = fs::read_to_string(plot_dir.join("depth_gc_scatter.svg")).expect("basic SVG");
    assert!(basic_svg.contains("Pandepth depth vs GC"));
    assert_eq!(basic_svg.matches("<circle ").count(), 5);
    assert!(!basic_svg.contains("NaN"));

    let styled_dir = case.path("plots/styled");
    run_ok(
        &case,
        &[
            "run".into(),
            "plot-depth-pandepth2".into(),
            "-i".into(),
            s(&pandepth),
            "-o".into(),
            s(&styled_dir),
            "-l".into(),
            "0".into(),
        ],
    );
    let styled_svg =
        fs::read_to_string(styled_dir.join("depth_gc_styled.svg")).expect("styled SVG");
    assert!(styled_svg.contains("Pandepth depth vs GC (styled)"));
    assert_eq!(styled_svg.matches("<circle ").count(), 5);
    assert_ne!(basic_svg, styled_svg);

    let mosdepth = case.write(
        "plots/mosdepth.tsv",
        "chr2\t0\t100\t5\nchr1\t0\t100\t10\nchr1\t100\t200\t20\n",
    );
    let mosdepth_dir = case.path("plots/mosdepth");
    run_ok(
        &case,
        &[
            "run".into(),
            "plot-mosdepth-point".into(),
            "-i".into(),
            s(&mosdepth),
            "-o".into(),
            s(&mosdepth_dir),
            "-l".into(),
            "0".into(),
        ],
    );
    assert_file(
        &mosdepth_dir.join("mosdepth_points.tsv"),
        "chrom\tstart\tend\tcoverage\tstatus\nchr2\t0\t100\t5\tOK\nchr1\t0\t100\t10\tOK\nchr1\t100\t200\t20\tOK\n",
    );
    let mosdepth_svg =
        fs::read_to_string(mosdepth_dir.join("mosdepth_scatter.svg")).expect("mosdepth SVG");
    assert!(mosdepth_svg.contains("Mosdepth point coverage"));
    assert_eq!(mosdepth_svg.matches("<circle ").count(), 5);
    assert!(!mosdepth_svg.contains("NaN"));
    assert_eq!(
        (
            fnv1a64(basic_svg.as_bytes()),
            fnv1a64(styled_svg.as_bytes()),
            fnv1a64(mosdepth_svg.as_bytes()),
        ),
        (
            11_779_402_563_777_481_233,
            17_556_672_687_734_655_657,
            17_764_904_749_355_255_246,
        ),
        "visual golden fingerprints changed"
    );
}

#[test]
fn bam_commands_match_external_goldens_when_samtools_is_available() {
    let samtools = Command::new("samtools").arg("--version").output();
    if samtools.is_err() {
        assert!(
            std::env::var_os("BIOHUB_REQUIRE_EXTERNAL_GOLDENS").is_none(),
            "samtools required by BIOHUB_REQUIRE_EXTERNAL_GOLDENS"
        );
        eprintln!("skipping BAM goldens: samtools unavailable");
        return;
    }

    let case = CaseDir::new("bam");
    let header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n";
    let r1_sam = case.write(
        "bam/r1.sam",
        &(header.to_string()
            + "pair1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            + "orphan1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tAAAA\tIIII\n"),
    );
    let r2_sam = case.write(
        "bam/r2.sam",
        &(header.to_string()
            + "pair1\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\tIIII\n"
            + "orphan2\t0\tchr1\t6\t60\t4M\t*\t0\t0\tCCCC\tIIII\n"),
    );
    let r1_bam = case.path("bam/r1.bam");
    let r2_bam = case.path("bam/r2.bam");
    for (sam, bam) in [(&r1_sam, &r1_bam), (&r2_sam, &r2_bam)] {
        let converted = Command::new("samtools")
            .args(["view", "-bS", &s(sam), "-o", &s(bam)])
            .output()
            .expect("convert SAM fixture");
        assert!(
            converted.status.success(),
            "{}",
            String::from_utf8_lossy(&converted.stderr)
        );
    }

    for command in ["merge-two-end-bam", "merge-two-end-bam1"] {
        let out1 = case.path(&format!("bam/{command}.r1.bam"));
        let out2 = case.path(&format!("bam/{command}.r2.bam"));
        run_ok(
            &case,
            &[
                "run".into(),
                command.into(),
                "-i".into(),
                s(&r1_bam),
                "-j".into(),
                s(&r2_bam),
                "-o1".into(),
                s(&out1),
                "-o2".into(),
                s(&out2),
            ],
        );
        for output in [out1, out2] {
            let viewed = Command::new("samtools")
                .args(["view", &s(&output)])
                .output()
                .expect("view BAM golden");
            assert!(viewed.status.success());
            let text = String::from_utf8(viewed.stdout).expect("SAM is UTF-8");
            assert_eq!(text.lines().count(), 1);
            assert!(text.starts_with("pair1\t"));
        }
    }

    let mgi_r1_sam = case.write(
        "bam/mgi-r1.sam",
        &(header.to_string() + "mgi001/1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"),
    );
    let mgi_r2_sam = case.write(
        "bam/mgi-r2.sam",
        &(header.to_string() + "mgi001/2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\tIIII\n"),
    );
    let mgi_r1_bam = case.path("bam/mgi-r1.bam");
    let mgi_r2_bam = case.path("bam/mgi-r2.bam");
    for (sam, bam) in [(&mgi_r1_sam, &mgi_r1_bam), (&mgi_r2_sam, &mgi_r2_bam)] {
        let converted = Command::new("samtools")
            .args(["view", "-bS", &s(sam), "-o", &s(bam)])
            .output()
            .expect("convert MGI fixture");
        assert!(converted.status.success());
    }
    let mgi_out1 = case.path("bam/mgi-out-r1.bam");
    let mgi_out2 = case.path("bam/mgi-out-r2.bam");
    run_ok(
        &case,
        &[
            "run".into(),
            "merge-two-end-bam-forMGI".into(),
            "-i".into(),
            s(&mgi_r1_bam),
            "-j".into(),
            s(&mgi_r2_bam),
            "-o1".into(),
            s(&mgi_out1),
            "-o2".into(),
            s(&mgi_out2),
        ],
    );
    for output in [mgi_out1, mgi_out2] {
        let viewed = Command::new("samtools")
            .args(["view", &s(&output)])
            .output()
            .expect("view MGI BAM golden");
        assert!(viewed.status.success());
        assert_eq!(String::from_utf8(viewed.stdout).unwrap().lines().count(), 1);
    }
}
