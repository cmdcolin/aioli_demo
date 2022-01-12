# aioli_demo

Demo using https://www.npmjs.com/package/@biowasm/aioli for data fetching and simple canvas rendering of reads and pileup

![](img/1.png)

Example

http://cmdcolin.github.io/aioli_demo/?file=https%3A%2F%2Fs3.amazonaws.com%2F1000genomes%2Fphase3%2Fdata%2FHG00096%2Falignment%2FHG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam&loc=1%3A10000-20000


Same thing with pure js bam parser https://github.com/cmdcolin/gmodbam_demo

Downsides compared to the pure js: parses the text output so is a bit slower than the plain gmodbam_demo, and also hits CORS errors on some files that gmodbam does not

See also https://github.com/cmdcolin/jbrowse-plugin-biowasm 

