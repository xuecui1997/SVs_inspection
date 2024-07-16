# SVs_inspection
Our detailed manual inspection methodology is based on aligning two genomes, facilitating the identification of base-pair resolution structural variations in tomatoes.
## Pipeline Overview
1. Query genome assembly
2. SV loci calling by multiple software
3. SV loci visualization
4. SV loci merging
5. SV region alignment
6. Similarity calculated among copies

## Dependencies
- [Hifiasm](https://github.com/chhylp123/hifiasm)
- [HiCanu](https://github.com/marbl/canu)
- [Flye](https://github.com/mikolmogorov/Flye)
- [Ragtag](https://github.com/malonge/RagTag)
- [GPM](https://github.com/Jianwei-Zhang/LIMS?tab=readme-ov-file)
- [NGMLR](https://github.com/philres/ngmlr)
- [pbmm2](https://github.com/PacificBiosciences/pbmm2)
- [minimap2](https://github.com/lh3/minimap2)
- [MUMmer4](https://github.com/mummer4/mummer)
- [LASTZ](https://github.com/lastz/lastz)
- [AnchorWave](https://github.com/baoxingsong/AnchorWave)
- [PBSV](https://github.com/PacificBiosciences/pbsv)
- [SVIM](https://github.com/eldariont/svim)
- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [Assemblytics](https://github.com/MariaNattestad/Assemblytics)
- [SyRI](https://github.com/schneebergerlab/syri)
- [SVision](https://github.com/xjtu-omics/SVision)
- [SVMU](https://github.com/mahulchak/svmu)
- [SAMtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)

Please refer to `markdown/pipeline.md` and `code/` for usage instructions.

Steps 3 to 5 involve manual inspection to analyze sequence features surrounding the SV comprehensively. 

Additional details are available in the supplementary information.
