#!/usr/bin/perl

%repNameCount = ();
%repClassCount = ();
%repFamilyCount = ();

open (IN, "GeneCounts_telescope_allSamples_annotated.txt");
while ($line= <IN>) {
	chomp($line);
	@data = split (/\t/, $line);
	if($data[0] ne "transcript" && $data[1] ne "NA") {
		if($repNameCount{$data[1]} eq "") {
			$repNameCount{$data[1]} = join("\t", $data[4], $data[5], $data[6], $data[7], $data[8], $data[9], $data[10], $data[11], $data[12], $data[13], $data[14], $data[15]);
			}
		else {
			@repNameCount2 = split (/\t/, $repNameCount{$data[1]});
			@repNameCount3 = ();
			$repNameCount4 = "";
			for $i (0..$#repNameCount2) {
				$repNameCount3[$i] = $repNameCount2[$i] + $data[$i+4];
				if($repNameCount4 eq "") {
					$repNameCount4 = $repNameCount3[$i];
					}
				else {
					$repNameCount4 = join ("\t", $repNameCount4, $repNameCount3[$i]);
					}
				}
			$repNameCount{$data[1]} = $repNameCount4;
			}
                if($repClassCount{$data[2]} eq "") {
                        $repClassCount{$data[2]} = join("\t", $data[4], $data[5], $data[6], $data[7], $data[8], $data[9], $data[10], $data[11], $data[12], $data[13], $data[14], $data[15]);
                        }
                else {
                        @repClassCount2 = split (/\t/, $repClassCount{$data[2]});
                        @repClassCount3 = ();
                        $repClassCount4 = "";
                        for $i (0..$#repClassCount2) {
                                $repClassCount3[$i] = $repClassCount2[$i] + $data[$i+4];
                                if($repClassCount4 eq "") {
                                        $repClassCount4 = $repClassCount3[$i];
                                        }
                                else {
                                        $repClassCount4 = join ("\t", $repClassCount4, $repClassCount3[$i]);
                                        }
                                }
                        $repClassCount{$data[2]} = $repClassCount4;
                        }
                if($repFamilyCount{$data[3]} eq "") {
                        $repFamilyCount{$data[3]} = join("\t", $data[4], $data[5], $data[6], $data[7], $data[8], $data[9], $data[10], $data[11], $data[12], $data[13], $data[14], $data[15]);
                        }
                else {
                        @repFamilyCount2 = split (/\t/, $repFamilyCount{$data[3]});
                        @repFamilyCount3 = ();
                        $repFamilyCount4 = "";
                        for $i (0..$#repFamilyCount2) {
                                $repFamilyCount3[$i] = $repFamilyCount2[$i] + $data[$i+4];
                                if($repFamilyCount4 eq "") {
                                        $repFamilyCount4 = $repFamilyCount3[$i];
                                        }
                                else {
                                        $repFamilyCount4 = join ("\t", $repFamilyCount4, $repFamilyCount3[$i]);
                                        }
                                }
                        $repFamilyCount{$data[3]} = $repFamilyCount4;
                        }
		} 
	}

open (OUT1, ">repName_collapsed_count_all_samples_telescope.txt");
print OUT1 "repName\tBO_1\tBO_2\tBO_3\tIP_1\tIP_2\tIP_3\tTOTAL_1\tTOTAL_2\tTOTAL_3\tINPUT_1\tINPUT_2\tINPUT_3\n";
while (($key, $value) = each(%repNameCount)) {
	print OUT1 "$key\t$value\n";
	}

open (OUT2, ">repClass_collapsed_count_all_samples_telescope.txt");
print OUT2 "repClass\tBO_1\tBO_2\tBO_3\tIP_1\tIP_2\tIP_3\tTOTAL_1\tTOTAL_2\tTOTAL_3\tINPUT_1\tINPUT_2\tINPUT_3\n";
while (($key, $value) = each(%repClassCount)) {
        print OUT2 "$key\t$value\n";
        }

open (OUT3, ">repFamily_collapsed_count_all_samples_telescope.txt");
print OUT3 "repFamily\tBO_1\tBO_2\tBO_3\tIP_1\tIP_2\tIP_3\tTOTAL_1\tTOTAL_2\tTOTAL_3\tINPUT_1\tINPUT_2\tINPUT_3\n";
while (($key, $value) = each(%repFamilyCount)) {
        print OUT3 "$key\t$value\n";
        }
