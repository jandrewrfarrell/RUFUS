####this script cleans up the RUFUS folder and output
#### only run when you are sure your run is complete and
#### you have everything you need. It will remove all 
#### your intermeidate files and youll have to do your 
#### run all over again


rm *Mutations.fastq.merged.bam* *generator.Mutations.fastq.pared.bam* *.generator.Mutations.Mate*.fastq* *.generator.temp* *generator.V2.overlap.fastq* *generator.V2.overlap.hashcount.fastq *generator.V2.overlap.hashcount.fastq.bam.vcf* *.generator.Jhash fastp.html fastp.json ./Intermediates/*.overlap.asembly.hash* Intermediates/*overlap.hashcount.fastq.Jhash* Intermediates/*.V2.ref.RepRefHash* mer_counts_merged.jf TempOverlap/*.generator.V2.*
