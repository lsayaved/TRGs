# Toxin-related genes (TRGs)
Collection of scripts to look for the TRGs described in the paper [Abundant toxin-related genes in the genomes of beneficial symbionts from deep-sea hydrothermal vent mussels, Sayavedra et al., 2015](https://elifesciences.org/content/4/e07966) and [Sayavedra et al., 2019](https://www.biorxiv.org/content/10.1101/605386v1)

Looking for toxin-related genes (TRGs) in bacterial genomes is not a straightforward task. The functional domains have often a mosaic-like structure that makes classification a difficult task. To look for TRGs in multiple genomes without the pain of doing it manually, I created a workflow of classification for the follow up paper [Comparative genomic inisghts into the roles of toxin-realted genes in the beneficial bacteria and their acquisition by horizontal gene transfer](https://www.biorxiv.org/content/10.1101/605386v1).

The folllowing is the draft of the protocol that was used! The code will be cleaned and updated.

- Blast protein sequences of the genomes you want to look at against the TRG proteins from the Bathy SOX described in [Sayavedra et al., 2015] (https://elifesciences.org/content/4/e07966).

```
blastp -query Genomes.faa -db ~/ownCloud/ToxinsOtherSOXs/BLAST/bathyt_toxins.faa -max_target_seqs 50 -outfmt "6 std qcovs" -evalue 0.001 -out Genes2Search_vsTRGs.txt
```

- Use minimum similarity of 25, minimum alignment of 200, minimum coverage of 15

```
cat Genes2Search_vsTRGs.txt  | awk '{if ($3>25 && ($13>25 && $4>50) || ($13>15 && $4>200) ) print $0}' | perl -pe "s/(\S+).*/\1/" | sort -u  >Genes2Search_vsTRGs.blastSim25Cov25.txt
```

- Submit genes from the TRGs paper to the CD-Search Tool, standard settings

- Download results full, include query defline

- Use the list that describes which domains are associated to which class of toxins [DomainsClassifier.txt]()

```
perl ~/Dropbox/perl_scripts/Findidentifiers_printseq.pl -inident BLAST_possibleToxinsInRelativesIDsSIm25COv50.txt -inseq ~�/ownCloud/ToxinsOtherSOXs/SpeciesTree/fnnGenomes/Genes2Search.faa  -o BLAST_possibleToxinsInRelativesSIm25COv50.faa
```
 
```
perl ClassifyTRGs.pl -domains DomainsClassifier.txt -CPD PossibleTRGsSim25Cov50Aln50_hitdata.txt  -faa ~/ownCloud/ToxinsOtherSOXs/MassiveClassiffier/BLAST_possibleToxinsInRelativesSIm25COv50.faa.fas -genids GenomeTags.txt -out PossibleTRGsSim25Cov50Aln50
```
