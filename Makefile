.PHONY = lichen lichen_clean mmetsp

mmetsp:
	./wrapper_phage_contigs_sorter_iPlant.pl -w ~/work/mmetsp/ -f ~/work/mmetsp/MMETSP0105.nt.fa --data-dir=~/work/virsorter-data/

lichen_clean:
	find ~/work/lichen -type d -name virsorter -exec rm -rf {} \;

lichen: lichen_clean
	./wrapper_phage_contigs_sorter_iPlant.pl -w ~/work/lichen/virsorter -f ~/work/lichen/Peltigera_aphthosa_scaffolds.fasta --data-dir=~/work/virsorter-data/
