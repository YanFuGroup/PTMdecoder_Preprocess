# MaxQuant to PTMdecoder conversion

See the [upstream conversion script](https://github.com/YanFuGroup/PTMdecoder/blob/main/utils/convert_MaxQuant_to_PTMdecoder.py) for the latest source. The commands below use the bundled copy in `external/PTMdecoder`, at the revision pinned by this repository.

For each sample group in the DNA damage phosphoproteomics dataset, use that group's MaxQuant `msms.txt` to generate PTMdecoder inputs.

Complete the [installation](../README.md#installation) first. Run the commands below from the repository root.

## Input and outputs

The input must contain the MaxQuant columns shown in the [fictional example](../tests/dna_phos/fictional_msms_split.txt). The converter writes:

- `pepSpecFile.txt`: modified peptide sequences with their MGF spectrum names and scan numbers, grouped by unmodified sequence.
- `filtered_result.txt`: an identification table for PTMdecoder, including retained unmodified rows.

## Example

Create a new output directory and convert the fictional example:

```text
python -c "from pathlib import Path; Path('outputs/dna_phos_example').mkdir(parents=True)"
python external/PTMdecoder/utils/convert_MaxQuant_to_PTMdecoder.py -i tests/dna_phos/fictional_msms_split.txt -o outputs/dna_phos_example/pepSpecFile.txt -r outputs/dna_phos_example/filtered_result.txt -s _HCDFT
```

For study data, replace the `-i` path with your `msms.txt` and choose new output paths.

## Conversion options

The `-s _HCDFT` option appends `_HCDFT` to each raw-file name when constructing MGF spectrum names. Set the suffix to match your MGF filenames.

Reverse/decoy hits are removed. By default, `CON__` protein assignments are removed while a PSM with another protein assignment is retained; a PSM with only `CON__` assignments is removed. Use `--strict-con-filter` to remove every PSM containing a `CON__` assignment. Choose the filtering mode used for your analysis.

See [output files](../README.md#output-files) for reader compatibility and [tests](../README.md#tests) for optional checks using fictional inputs.
