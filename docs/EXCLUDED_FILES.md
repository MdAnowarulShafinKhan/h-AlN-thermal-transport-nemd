# Files intentionally excluded from the GitHub package

The original archive contained material that is useful for private research but should not be mirrored into a public source/data repository.

## Excluded categories

1. **25 third-party PDF files** under the original `literature/` and journal-guideline folders. These are publisher/author papers or journal instructions and are not project-generated research artifacts. Bibliographic links are provided in `references/REFERENCES.md` instead.
2. **Prior ChatGPT transcript PDF** (`chatgpt_whole_chat/Abid_AlN_FINAL_ALN_combined.pdf`). It is not scientific raw data, code, or a manuscript.
3. **Original `validation_table.docx`**. Its scientific information is retained in a corrected, GitHub-readable `docs/VALIDATION_BENCHMARKS.md`, while the original DOCX had severe portrait-table wrapping and an outdated sentence stating 114 W m^-1 K^-1.
4. **14 redundant copies** of `AlN_width_13.7nm.lmp` in the thermostat-bias folders. All 15 copies in the source archive were byte-identical. One canonical copy is stored at `structures/pristine/AlN_width_13.7nm.lmp`.
5. **Transient or absent simulation trajectories**. Large dump/trajectory files that were not supplied in the source archive cannot be reconstructed from the processed outputs.

No manuscript file was found in the source archive, so no manuscript was omitted from this GitHub package.
