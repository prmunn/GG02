#! /usr/bin/awk -f

# Usage
# ./horizontalMerge-padUMI.awk R2_info.txt <(paste I2.fastq I1.fastq R2.fastq) > I2-I1-padUMI-R2.fastq
# ./horizontalMerge-padUMI.awk R2_info.txt <(paste <(zcat I2.fq.gz) <(zcat I1.fq.gz) <(zcat R2.fq.gz)) > I2-I1-padUMI-R2.fastq

BEGIN {
	FS="\t";
	}

	NR==FNR	{
		split($1,c," ");				#use only first field in readID for matching, second field varies by read (I1, I2, R2)
		a[NR]="@"c[1];					#read ID (for checking purposes), info file drops leading '@', add back to match fastq header row
		if($2 == -1 || $3 < 14 || $3 > 19) p=0; 	#ME not found in correct position, no pad added? or add 'adapter' to ID/filter with cutadapt?
		else p=20-$3;					#pad length
		b[NR]=p;					#set pad length for this read
		next;
	}

	{
	x=FNR%4;		#4-row counter for fastq file
	if(x==1) {		#read header row
		t++;		#read counter
		n=a[t];		#expected read name
        	p=b[t];		#pad distance, based on info file
		split($2,s," ");
		split($3,r," ");
 		if(n != s[1] || n != r[1]) { 		#error! read name mismatch	for some reason, I2 headers are different (first fields)!!??
			print "Read name mismatch at read "t"! Expected "n" but found I1="s[1]" and R2="r[1]". Exit!" | "cat 1>&2"
			exit 1;
			}
		else print $3;				# important to print R2 header (R2 output)
		next;
		}
	if(x==2) {					#pad read
		ORS="";					#this seems to work
		print $1;				#I2 read (i5tag-BC) -- put i5tag-BC at 5' end, optional sample-level demux on i5tag subset (but would need to keep i5tagBC in R2 after demux!)
		print $2;				#I1 read (i7PCR-BC)
		if(p==6) print "ATATGC";
		if(p==5) print "TATGC";
		if(p==4) print "ATGC";
		if(p==3) print "TGC";
		if(p==2) print "GC";
		if(p==1) print "C";
		#if(p==0) print "NO ME FOUND";		#could insert 5' RC-MEq sequence to cause read to be trimmed/filtered later
		ORS="\n";
		print $3;				#future: could substr($3) to switch order of UMI (print i7tag-BC before padding) and delete ME at this step
		next;
		}
	if(x==3) {
		print $3;
		next;
		}
	if(x==0) {
		ORS="";
		print $1;                            	#I2 read (i5tag-BC). Match order and any changes made to read seq above.
		print $2;                            	#I1 read (i7PCR-BC). Match order and any changes made to read seq above.
		for(i=1;i<=p;i++) {
			print "I";			# high base quality, preserve padded base identities through any filtering steps
			}
		ORS="\n";
		print $3;				#future: match any substr($3) mods to read seq above
		}
}
