for loop in $(cat list_pe)
    do
    fastqc $loop.R1.fq $loop.R2.fq
    done

for loop in $(cat list_se)
    do
    fastqc $loop.R1.fq
    done

(hisat2 -p 8 --no-unal --dta -x {params.index} -1 {input.r1} -2 {unput.r2}) 2>{log} | samtools sort -@ 8 -o {output.bam} -
hisat2 -p 8 --no-unal --dta -x ML975 -1 DameSB_mal.R1.fq -2 DameSB_mal.R2.fq -S DameSB_mal.sam | samtools sort -@ 8 -o DameSB_mal.bam 
hisat2 -p 8 --no-unal --dta -x ML975 -1 DameSB_fem.R1.fq -2 DameSB_fem.R2.fq -S DameSB_fem.sam | samtools sort -@ 8 -o DameSB_fem.bam

for loop in $(cat list_se)
    do
    hisat2 -p 8 --no-unal --dta -x ML975 -U $loop.R1.fq -S $loop.sam | samtools sort -@ 8 -o $loop.bam
    done

for loop in $(cat list_pe)
    do 
    hisat2 -p 8 --no-unal --dta -x ML975 -1 $loop.R1.fq -2 $loop.R2.fq -S $loop.sam | samtools sort -@ 8 -o $loop.bam
    done

for loop in $(cat list_sam)
    do
    samtools sort -@ 8 -o $loop.bam $loop.sam
    done


for loop in $(cat list_sam)
    do
    stringtie $loop.bam -p 8 -o $loop.gtf -l $loop 
    done

stringtie --merge -p 8 -o stringtie_merged.gtf *.gtf

awk '{if ($3 == "transcript") print $10"\t"$12}' stringtie_merged.gtf | sed "s/;//g" | sed "s/\\\"//g" | sort -u > geneTransMap

featureCounts -a stringtie_merged.gtf -p -T 8 -F GTF -o stringtie_gene_counts.txt ../bam/*

cat stringtie_gene_counts.txt | cut -f 1,7- | sed 1d > stringtie_gene_counts_matrix.txt

qualimap rnaseq -bam bam/* --java-mem-size=32G -gtf gtf/stringtie_merged.gtf -outdir st_qualmap_rnaseq

cat st_qualmap_rnaseq/qualimapReport.html  | grep ":" | sed "s/://g" >> summary_files.txt

multiqc -f -o multiqc_report.html -l -dd 2 qc

gffread gtf/stringtie_merged.gtf -g Dame_ML975.contigs.fasta -w stringtie_transcripts.fa 

gtf_to_alignment_gff3.pl gtf/stringtie_merged.gtf > gtf/stringtie_transcripts.gff3

TransDecoder.LongOrfs -t stringtie_transcripts.fa --gene_trans_map gtf/geneTransMap -m 50 --output_dir transdecoder_dir

blastp -query transdecoder_dir/longest_orfs.pep -db /home/mihtet/Annotations/Trinotate/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > BLASTp_init.outfmt6

hmmscan --cpu 8 --domtblout pfam_i.domtblout /home/mihtet/Annotations/Trinotate/Pfam-A.hmm /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/longest_orfs.pep 

TransDecoder.Predict -t stringtie_transcripts.fa --retain_pfam_hits pfam/pfam_i.domtblout --retain_blastp_hits BLAST/BLASTp_init.outfmt6 --output_dir transdecoder_dir 

mv stringtie_transcripts.fa.transdecoder.* transdecoder_dir

blastx -query stringtie_transcripts.fa -db /home/mihtet/Annotations/Trinotate/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > BLAST/swissprot.blastx.outfmt6

blastp -query /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/stringtie_transcripts.fa.transdecoder.pep -db /home/mihtet/Annotations/Trinotate/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > BLAST/swissprot.blastp.outfmt6

hmmscan --cpu 16 --domtblout TrinotatePFAM.out /home/mihtet/Annotations/Trinotate/Pfam-A.hmm /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/stringtie_transcripts.fa.transdecoder.pep

/home/mihtet/Annotations/signalP/signalp-4.1/signalp -f short -n signalp/signalP.out /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/stringtie_transcripts.fa.transdecoder.pep

/home/mihtet/Annotations/tmhmm/tmhmm-2.0c/bin/tmhmm --short < /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/stringtie_transcripts.fa.transdecoder.pep > tmhmm/tmhmm.out

perl /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl \
    --transcriptome stringtie_transcripts.fa \
    --path_to_rnammer /home/mihtet/Annotations/rnammer/rnammer \

mv stringtie_transcripts.fa.rnammer.gff gtf

Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite init \
    --gene_trans_map /home/mihtet/Annotations/Dame_ML975/gtf/geneTransMap \
    --transcript_fasta stringtie_transcripts.fa \
    --transdecoder_pep /home/mihtet/Annotations/Dame_ML975/transdecoder_dir/stringtie_transcripts.fa.transdecoder.pep \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_swissprot_blastp /home/mihtet/Annotations/Dame_ML975/BLAST/swissprot.blastp.outfmt6 \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_swissprot_blastx /home/mihtet/Annotations/Dame_ML975/BLAST/swissprot.blastx.outfmt6 \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_pfam TrinotatePFAM.out \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_tmhmm /home/mihtet/Annotations/Dame_ML975/tmhmm/tmhmm.out \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_signalp /home/mihtet/Annotations/Dame_ML975/signalp/signalP.out \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite LOAD_rnammer /home/mihtet/Annotations/Dame_ML975/gtf/stringtie_transcripts.fa.rnammer.gff \
    && Trinotate /home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/admin/Trinotate.sqlite report > Trinotate_report.xls

/home/mihtet/Annotations/trinotate/Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl \
    --Trinotate_xls Trinotate_report.xls -G -I > Trinotate_report.xls.gene_ontology

mv Trinotate_report* report
