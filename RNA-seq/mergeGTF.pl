#!/usr/bin/perl

%info1 = ();
%info2 = ();

open (IN1, "mm10_UCSC_rmsk.gtf");
while ($line= <IN1>) {
	chomp($line);
	@data = split (/\t/, $line);
	$gene = $data[8];
	$gene =~ s/";.*$//;
	$gene =~ s/^.*"//;
	$ann = join("\"", "gene_id ", $gene, "; gene_name ", $gene, "; p_id ", $gene, "; transcript_id ", $gene, "; tss_id ", $gene, ";");  
	$ann2 = join("\t", $data[0], "unknown", $data[2], $data[3], $data[4], $data[5], $data[6], $data[7], $ann);
	if($info1{$gene} eq "") {
		$info1{$gene} = $ann2;
		$loc = join("\t", $data[0], $data[3]);
		if($info2{$loc} eq "") {
			$info2{$loc} = $gene;
			}
		else {
			$info2{$loc} = join("\t", $info2{$loc}, $gene);
			}
		}	
	else {
		$info1{$gene} = join("\n", $info1{$gene}, $ann2);
		}
	}

open (IN2, "/data/zappadata_p1/reference/mm10/rsemindex/star/mm10_UCSC_noAlt.gtf");
while ($line= <IN2>) {
	chomp($line);
	@data = split (/\t/, $line);
	$gene = $data[8];
	$gene =~ s/";.*$//;
	$gene =~ s/^.*"//;
	if($info1{$gene} eq "") {
		$info1{$gene} = $line;
		$loc = join("\t", $data[0], $data[3]);
		if($info2{$loc} eq "") {
                        $info2{$loc} = $gene;
                        }
                else {
                        $info2{$loc} = join("\t", $info2{$loc}, $gene);
                        }	
		}
	else {
		$info1{$gene} = join("\n", $info1{$gene}, $line);
		}
	}
	
open (IN3, "chrlist.txt");
open (OUT, ">mm10_UCSC_refSeq_rmsk.gtf");
while ($line= <IN3>) {
	chomp($line);
	@pos = ();
	$count = 0;
	while (($key, $value) = each(%info2)) {
		@data = split (/\t/, $key);
		if($data[0] eq $line) {
			$pos[$count] = $data[1];
			$count++;
			}
		}
	@sorted_pos = sort { $a <=> $b } @pos;
	for $i (0..$#sorted_pos) {
		$gene2 = $info2{join("\t", $line, $sorted_pos[$i])};
		@data = split (/\t/, $gene2);
		for $j (0..$#data) {
			print OUT "$info1{$data[$j]}\n";
			}
		}
	}
