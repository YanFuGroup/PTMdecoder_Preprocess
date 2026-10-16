# Preprocessing scripts for the PTMdecoder study

These scripts reproduce the preprocessing used in the PTMdecoder study, converting database-search results into input files for [PTMdecoder](https://github.com/YanFuGroup/PTMdecoder). The repository covers preprocessing only.

| Dataset | Input |
| --- | --- |
| Synthetic peptides, measured individually and in mixtures (`syn`) | Mascot DAT |
| MCF7 cells treated with JIB-04 or SAHA, with DMSO controls | Mascot DAT |
| Histone peptides from the DNA damage experiment | Mascot DAT |
| Phosphoproteomics from the DNA damage experiment | MaxQuant `msms.txt` |

The study data and search results are listed under [PXD062791](https://proteomecentral.proteomexchange.org/?pxid=PXD062791). See the paper's Data availability statement for access conditions.

## Installation

- Mascot preprocessing: MATLAB R2022a and Git on `PATH`. This workflow has been validated on Windows; other operating systems have not been verified.
- MaxQuant conversion: Python 3.9+ with `pandas>=2.2.3`.

Clone the repository with its PTMdecoder dependency:

```text
git clone --recurse-submodules https://github.com/YanFuGroup/PTMdecoder_Preprocess.git
cd PTMdecoder_Preprocess
```

For an existing clone, run `git submodule update --init --recursive`. For MaxQuant conversion, also install pandas:

```text
python -m pip install "pandas>=2.2.3"
```

## Prepare the input data

For Mascot datasets, copy `local_paths.template.yaml` to `local_paths.yaml` and set the input root outside this repository:

```yaml
data_root: "<path-to-search-results>"
```

Place the DAT files in the corresponding run directories. Prepare only the dataset you want to run:

```text
data_root/
  syn_ind/mascot_dat/ind1/ ... ind16/
  syn_mix/mascot_dat/mix1/ ... mix18/
  MCF7_JIB04/mascot_dat/DMSO_1/, DMSO_2/, JIB04_1/, JIB04_2/
  MCF7_SAHA/mascot_dat/DMSO_R1/, DMSO_R2/, SAHA_R1/, SAHA_R2/
  DNA_histone/mascot_dat/wt1/, wt2/, wt3/, ir1/, ir2/, ir3/
```

For the DNA damage phosphoproteomics dataset, follow the [MaxQuant conversion guide](dna_phos/README.md).

## Run the preprocessing

Start MATLAB in the repository root and run this setup once:

```matlab
addpath(fullfile(pwd, 'synthetic'));
addpath(fullfile(pwd, 'mcf7'));
addpath(fullfile(pwd, 'dna_histone'));
if ~isfolder('outputs'), mkdir('outputs'); end
```

Run the command for your dataset, choosing a new output directory each time.

Individual synthetic peptides:

```matlab
run_syn_ind_preprocess(fullfile(pwd, 'outputs', 'syn_ind'));
```

Mixed synthetic peptides:

```matlab
run_syn_mix_preprocess(fullfile(pwd, 'outputs', 'syn_mix'));
```

MCF7 cells treated with JIB-04:

```matlab
run_mcf7_jib04_auto_preprocess(fullfile(pwd, 'outputs', 'mcf7_jib04'));
```

MCF7 cells treated with SAHA:

```matlab
run_mcf7_saha_auto_preprocess(fullfile(pwd, 'outputs', 'mcf7_saha'));
```

DNA damage histone peptides:

```matlab
run_dna_histone_preprocess(fullfile(pwd, 'outputs', 'dna_histone'));
```

MaxQuant conversion uses Python; its complete command is in the [conversion guide](dna_phos/README.md).

### Optional Top 1 comparison

For the conventional Top 1 comparison, use the same MATLAB setup, add the Top 1 directory, and run the corresponding entry point:

```matlab
addpath(fullfile(pwd, 'mcf7', 'top1'));
run_mcf7_jib04_top1_preprocess(fullfile(pwd, 'outputs', 'jib04_top1'));
run_mcf7_saha_top1_preprocess(fullfile(pwd, 'outputs', 'saha_top1'));
```

These routes use the same selection as auto and also write `report_msms_top1.txt` for the comparison workflow, which skips PTMdecoder MS/MS analysis. The report writer does not perform an additional highest-scoring-PSM selection.

## Preprocessing rules

Mascot peptide-spectrum matches (PSMs) are grouped using the `REVERSE_` decoy tag and `_HUMAN` protein-group tag, with an FDR threshold of 0.01. MCF7 and histone datasets also apply chemical-prior filtering before FDR.

After FDR, the following rules select entries for `pepSpecFile.txt`, sorted by peptide:

| Dataset | Additional selection |
| --- | --- |
| Synthetic individual peptides | Peptide `GKGGKGLGKGGAKR`; charge 2; neutral-mass tolerance `< 0.0001` around 1535.878357; at least 8 matched fragment ions; unique (`DatasetName`, `Scan`) pairs. |
| Synthetic mixtures | Same peptide, charge, mass tolerance, and uniqueness rule; no matched-ion minimum. |
| MCF7 JIB-04 / SAHA | Historical text-key uniqueness; retains modified peptides without Oxidation, using the experiment-specific selector. |
| DNA damage histone peptides | Historical text-key uniqueness; excludes unmodified, Oxidation-containing, and label-only PSMs. |

Synthetic scripts compare the two fields (`DatasetName`, `Scan`). MCF7 and histone scripts preserve the historical `strcat(DatasetName, ' ', Scan)` key, which removes the separator space. Both rules discard every PSM in a repeated-key group.

## Output files

Mascot outputs are written to a subfolder for each run under the chosen output directory.

| File | Contents |
| --- | --- |
| `group_result_mascot.txt` | Grouped PSMs before FDR filtering, produced for MCF7 and histone datasets. |
| `filtered_result_mascot.txt` | FDR-filtered identification table, before the additional peptide-spectrum selection above. |
| `pepSpecFile.txt` | Final peptide-spectrum list for PTMdecoder MS/MS analysis. |
| `report_msms_top1.txt` | Additional report from the optional MCF7 Top 1 routes. |

Identification tables are tab-delimited. The pinned PTMdecoder reader includes the [TSV parsing fix](https://github.com/YanFuGroup/PTMdecoder/pull/19), preserving spaces within fields and the positions of empty columns.

At the previous dependency pin (`c36a8f6`), the reader also split on spaces and removed empty fields, which could shift later columns.

MaxQuant outputs are described in the [conversion guide](dna_phos/README.md#input-and-outputs).

## Tests

Optional checks use fictional inputs to test preprocessing rules and output generation. Run from the repository root.

In MATLAB R2022a:

```matlab
results = verify_package;
```

For the MaxQuant converter:

```text
python -m unittest discover -s tests/dna_phos -p test_converter.py -v
```
