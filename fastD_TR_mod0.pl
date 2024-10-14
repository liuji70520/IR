#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

###
## ___________________________________________________________________________________________________________________
## 
## DESCRIPTION: tool for detecting insecticide resistance-associated mutations in insecticide targets from RNA-seq.
## 
## Usage:
## 
## To detect resistance-associated target-site mutations and their frequencies:
## 
##             FastD_TR.pl [options]* -s [species] -r [target] -i [test samples]
##
## To identify differential non-synonymous target-site mutations as novel markers:
##
##             FastD_TR.pl [options]* -s [species] -r [target] -1 [Sensitive samples] -2 [Resistant samples]
## 
##
##             -s/--species      The name of species you want to detect, please combine names with "_".
##                               Such as using "Plutella_xylostella" instead of "Plutella xylostella".
##             -r/--target       The protein symbols of insecticide targets : such as 'ACE','VGSC' and so on.
##                               You can enter more than one targets per time.
##             -i/--input        The input file of the program which is a sam file.
##
##             -1/--sensitive    The input alignment files generated from sensitive samples.
##
##             -2/--resistant    The input alignment files generated from resistant samples.
##
## Options:
## 
##             -o/--out       The prefix of output files.
##             -v/--version   Print version information and quit.
##             -h/--help      print this usage message.
## 
## EXAMPLE:  
## 
##          perl FastD_TR.pl -s Plutella_xylostella -r RyR VGSC AChE -i *.sam
##
##          perl FastD_TR.pl -s Plutella_xylostella -r RyR VGSC AChE -1 *.sam -2 *.sam 
##
## ____________________________________________________________________________________________________________________
##



my $srcpath = dirname $0 ;  ###  FastD_TR.pl filepath;
my $time = time ();
##################################################
### The necessary argument.

my $species="";
my @receptor=();
my @input=();
my @input1=();
my @input2=();
my $out_file='result';
my $out_path=".";
my $db_path="";
my $help=0;
my $version=0;

GetOptions(
			'species|s=s'                  => \$species,
			'receptor|r=s{1,}'             => \@receptor,
			'i=s{1,}'                      => \@input,
			'1=s{1,}'                      => \@input1,
			'2=s{1,}'                      => \@input2,
			'out|o=s'                      => \$out_path,
			'db|d=s'                       => \$db_path,
			'help|h'                       => \$help,
			'version|v'                    => \$version,
);
&HELP if $help ;
&VERSION if $version;


my @sam_file;
if (@input1 && @input2){
	push @sam_file, \@input1;
	push @sam_file, \@input2;
}elsif(@input){
	push @sam_file, \@input;
}else{
	die "Input error!\n";
}

if ($#sam_file == 1) {
	my $kmer = 500;
	mkdir "$out_path" || die "can`t make dir :$!";
	mkdir "./$out_path/Mutation_summary" || die "can`t make dir :$!";
	mkdir "./$out_path/Mutation_visualization" || die "can`t make dir :$!";

	#---------------------------------------Open result files-----------------------------------------

	#open OUT,">","${out_file}.txt" || die "can`t open $!\n";

	#----------------------------------Reading Reference Sequences------------------------------------

	my %pro_gen_relation = ('VGSC' => 'para','RyR' => 'ryr','AChE' => 'ace','nAChR' => 'nachr');

	my %reference;
	my $header;
	for (my $i = 0; $i <= $#receptor; $i++) {
		open REF, "<", "$db_path/$receptor[$i]/${species}_$pro_gen_relation{$receptor[$i]}.fasta" || die "can`t open $!\n";
		while (<REF>) {
			chomp;
			if (/^>([^\s]+?) /) {
				$header = $1;
			}else{
				$reference{$header} = $_;
			}
		}
	}

	#-------------------------------Reading Target mutation information------------------------------

	open LOC,"<$db_path/mutation_loci.txt" || die "can`t open $!\n"; #mutation file
	my %locusHash = ();
	while(<LOC>) {
		chomp;
		my @arr = split /\t/,$_;
		next if ($#arr < 1);
		if ($arr[0] eq $species) {
			my %temp = ();
			$temp{'species'} = $arr[0];
			$temp{'protein'} = $arr[1];
			$temp{'accession'} = $arr[2];
			$temp{'gene'} = $arr[3];
			$temp{'loci'} = $arr[4];
			$temp{'reference'} = $reference{$arr[2]};
			$locusHash{$arr[1]}{$arr[3]} = \%temp;
		}
	}
	close LOC;

	#---------------------Relationship between codons and amino acid residues-------------------------

	open CDON,"<","$db_path/codon_relationship.txt" || die "can`t open :$!\n";
	my %codon_relation_hash;
	while (<CDON>) {
		chomp;
		my @cdon_temp = split /\t/, $_;
		my @cdon_cate = split /,/, $cdon_temp[1];
		for (my $i = 0; $i <= $#cdon_cate; $i++) {
			$codon_relation_hash{$cdon_cate[$i]} = $cdon_temp[0];
		}
	}
	close CDON;

	#-------------------------------------------------------------------------------------------------

	#my $subunit = $locusHash{$receptor[0]}{'accession'};
	#my $reference = "${species}". "_$locusHash{$receptor[0]}{'gene'}" . ".fasta";

	#-----------------------Extract and Modify reads from SAM files-----------------------------------

	my %all_alignment;

	for (my $k = 0; $k <= $#receptor; $k++) {
		foreach my $gene (keys %{$locusHash{$receptor[$k]}}) {
			my @alignment;
			for (my $m = 0; $m <= $#sam_file; $m++) {
				my @new_reads;
				for (my $n = 0; $n <= $#{$sam_file[$m]}; $n++) {
					open SAM,"<","${$sam_file[$m]}[$n]" || die "Can`t find file $! ! Please check filepath \n";
					while (<SAM>) {
						chomp;
						next if (/^@/);
						my @reads = split /\t/, $_;
						if ($reads[2] eq $locusHash{$receptor[$k]}{$gene}{'accession'}) {

							my @array_match = $reads[5] =~ /\d+[A-Z]/g;
							my @array;

							my $seq;
							my $sum = 0;
							my $length;
							for (my $i = 0; $i <= $#array_match; $i++) {
								my $num;
								my $receptor;
								
								if ($array_match[$i] =~ /(\d+)([A-Z])/) {
									$num= $1;
									$receptor = $2;
								}
								if ($i == 0 && $receptor eq "M") {
									$seq .= substr($reads[9],0,$num);
									$sum += $num;
									$length += $num;
								}
								if ($i != 0 && $receptor eq "M") {
									$seq .= substr($reads[9],$sum,$num);
									$sum += $num;
									$length += $num;
								}
								if($receptor eq "I"){
									$sum += $num;
								}
								if($receptor eq "D"){
									$seq .= "X" x $num;
									$length += $num;
								}
							}

							$array[0] = $reads[3];
							$array[1] = $length;
							$array[2] = $seq;

							my $x = int(($reads[3] + $length - 1)/$kmer);

							push(@{$new_reads[$x]}, \@array);

							if ($reads[3] < $x * $kmer && $reads[3] >= $x * $kmer - $length + 1) {
								my $y = $x - 1;
								push(@{$new_reads[$y]}, \@array);
							}
						}
					}
				}
				$alignment[$m] = \@new_reads;
			}
			$all_alignment{$receptor[$k]}{$locusHash{$receptor[$k]}{$gene}{'gene'}} = \@alignment;
		}
	}

	#----------------------------------------Program test-----------------------------------------------
=pod
	for (my $n = 0; $n <= $#{$all_alignment{'VGSC'}{'para'}}; $n++) {
		for (my $m = 0; $m <= $#{$all_alignment{'VGSC'}{'para'}[$n]}; $m++) {
			for (my $l = 0; $l <= $#{$all_alignment{'VGSC'}{'para'}[$n][$m]}; $l++) {
				print "$all_alignment{'VGSC'}{'para'}[$n][$m][$l][0]\t$all_alignment{'VGSC'}{'para'}[$n][$m][$l][1]\t$all_alignment{'VGSC'}{'para'}[$n][$m][$l][2]\n";
			}
		}
	}
=cut
	#---------------------------------------------------------------------------------------------------

	foreach my $target (keys %all_alignment) {
		foreach my $gene (keys %{$all_alignment{$target}}) {
			open OUT, ">", "./$out_path/Mutation_summary/${target}_${gene}.csv" || die "can`t open :$!\n";
			print OUT "Position, Frequency(A), Frequency(G), Frequency(T), Frequency(C), Reference, Variant, Substitution, Coverage_1, Coverage_2, Sequence_1, Sequence_2\n";
			my %seq_logo;
			for (my $j = 1; $j <= length $locusHash{$target}{$gene}{'reference'} ; $j++) {
				my $k = int($j/$kmer);

				my @coverage;
				my @cover_nt;
				my @arr_SR;
				my @times_SR;

				for (my $n = 0; $n <= $#{$all_alignment{$target}{$gene}}; $n++) {
					$coverage[$n] = 0;
					for (my $i = 0; $i <= $#{$all_alignment{$target}{$gene}[$n][$k]}; $i++) {

						if ( $j >= $all_alignment{$target}{$gene}[$n][$k][$i][0] && $j <= $all_alignment{$target}{$gene}[$n][$k][$i][0] + $all_alignment{$target}{$gene}[$n][$k][$i][1] - 1) {
							$coverage[$n]++;
							my $index = $j - $all_alignment{$target}{$gene}[$n][$k][$i][0];
							my $nt = substr($all_alignment{$target}{$gene}[$n][$k][$i][2],$index,1);
							push(@{$arr_SR[$n]}, $nt);
							$cover_nt[$n] .= $nt;
						}
					}
					my %times;
					for (my $i = 0; $i <= $#{$arr_SR[$n]}; $i++) {
						if ($arr_SR[$n][$i] =~ /([A-Z])/) {
							$times{$1}++;
						}
					}
					$times_SR[$n] = \%times;
				}

				my %fre;
				my @arr_temp = ("A","G","T","C");
				foreach my $x (@arr_temp) {
					for (my $i = 0; $i <= $#times_SR; $i++) {
						if ($times_SR[$i]{$x}) {}else{$times_SR[$i]{$x} = 0;}
						if ($coverage[$i] == 0) {$coverage[$i] = 1;}
					}
					$fre{$x} = abs(($times_SR[0]{$x}/$coverage[0]) - ($times_SR[1]{$x}/$coverage[1]));
				}

				if ($coverage[0] < 30 or $coverage[1] < 30) {
					}else{
					if ($fre{'A'} < 0.4 && $fre{'T'} < 0.4 && $fre{'C'} < 0.4 && $fre{'G'} < 0.4) {
						}else{
						my $top_S = max(values %{$times_SR[0]});
						my $gt_S;
						foreach my $x (sort keys %{$times_SR[0]}){
							if ($times_SR[0]{$x} == $top_S) {
								$gt_S = $x;
							}
						}
						my $top_R = max(values %{$times_SR[1]});
						my $gt_R;
						foreach my $x (sort keys %{$times_SR[1]}){
							if ($times_SR[1]{$x} == $top_R) {
								$gt_R = $x;
							}
						}

						my $left = $j%3;
						my @codon;
						if ($left == 0) {
							my $first2 = substr($locusHash{$target}{$gene}{'reference'}, $j - 3,2);
							$codon[0] = $first2 . $gt_S;
							$codon[1] = $first2 . $gt_R;
						}elsif ($left == 1) {
							my $last2 = substr($locusHash{$target}{$gene}{'reference'}, $j ,2);
							$codon[0] = $gt_S . $last2;
							$codon[1] = $gt_R . $last2;
						}elsif ($left == 2) {
							my $first = substr($locusHash{$target}{$gene}{'reference'}, $j - 2 ,1);
							my $last = substr($locusHash{$target}{$gene}{'reference'}, $j ,1);
							$codon[0] = $first . $gt_S . $last;
							$codon[1] = $first . $gt_R . $last;
						}

						if ($codon_relation_hash{$codon[0]} ne $codon_relation_hash{$codon[1]}) {
							print OUT "$j";
							foreach my $y (@arr_temp) {
								$fre{$y} = sprintf "%.2f", $fre{$y};
								print OUT ",$fre{$y}";
							}
							print OUT ",$gt_S,$gt_R,$codon_relation_hash{$codon[0]} -> $codon_relation_hash{$codon[1]},$coverage[0],$coverage[1],$cover_nt[0],$cover_nt[1]\n";
							$seq_logo{$j} = \@arr_SR;
						}
					}
				}
			}
			close OUT;
			my @aa = qw(A C D E F G H I K L M N P Q R S T V W Y X);
#			my $unavailable = 0;
			for (my $i = 0; $i <= 1; $i++) {
				open MTX, ">", "${target}_${gene}_result_$i.txt" || die "Can`t find file $!";
				my $file_name = "${target}_${gene}_" . "result_" . $i;
				print MTX "A\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\tX\n";
				foreach my $pos (keys %seq_logo) {
					print MTX "$pos";
					for (my $k = 0; $k <=$#{$seq_logo{$pos}[$i]}; $k++) {
						my $left = $pos%3;
						my $codon;
						if ($left == 0) {
							my $first2 = substr($locusHash{$target}{$gene}{'reference'}, $pos - 3,2);
							$codon = $first2 . $seq_logo{$pos}[$i][$k];
						}elsif ($left == 1) {
							my $last2 = substr($locusHash{$target}{$gene}{'reference'}, $pos ,2);
							$codon = $seq_logo{$pos}[$i][$k] . $last2;
						}elsif ($left == 2) {
							my $first = substr($locusHash{$target}{$gene}{'reference'}, $pos - 2 ,1);
							my $last = substr($locusHash{$target}{$gene}{'reference'}, $pos ,1);
							$codon = $first . $seq_logo{$pos}[$i][$k] . $last;
						}

						$seq_logo{$pos}[$i][$k] = $codon_relation_hash{$codon};
						if ($codon =~ /[NX]/) {
							$seq_logo{$pos}[$i][$k] = "X";
						}
					}
					my %aa_num;
					
					for (my $v = 0; $v <= $#aa; $v++) {
						$aa_num{$aa[$v]} = 0;
						for (my $l = 0; $l <=$#{$seq_logo{$pos}[$i]}; $l++) {
#							if ($seq_logo{$pos}[$i][$l]) {
#								
#							}else{
#								$unavailable++;
#								print TST "$pos\t$i\t$l\tchen\n";
#								$seq_logo{$pos}[$i][$l] = "X";
#							}
							if ($aa[$v] eq $seq_logo{$pos}[$i][$l]) {
								$aa_num{$aa[$v]}++;
							}
						}
						print MTX "\t$aa_num{$aa[$v]}";
					}
					print MTX "\n";
 				}
 				close MTX;
 				system "Rscript $srcpath/weblogo_mod.R $file_name";
				system "rm $file_name.txt";
				system "mv ./$file_name.pdf ./$out_path/Mutation_visualization";
			}
		}
	}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
}elsif ($#sam_file == 0) {

	mkdir "$out_path" || die "can`t make dir :$!";
	mkdir "./$out_path/Mutation_detected_details" || die "can`t make dir :$!";
	mkdir "./$out_path/Mutation_visualization" || die "can`t make dir :$!";
	mkdir "./$out_path/Mutation_table" || die "can`t make dir :$!";

	open LOC,"<$db_path/mutation_loci.txt" || die "can`t open $!\n"; #mutation file
	my @locusArray = ();
	while(<LOC>) {
		chomp;
		my @arr = split /\t/,$_;
		next if ($#arr < 1);
		if ($arr[0] eq $species) {
			my %temp = ();
			$temp{'species'} = $arr[0];
			$temp{'protein'} = $arr[1];
			$temp{'accession'} = $arr[2];
			$temp{'gene'} = $arr[3];
			$temp{'loci'} = $arr[4];
			push @locusArray, \%temp;
		}
	}
	close LOC;

	my $locusArrRef = \@locusArray;
	open CDON,"<","$db_path/codon_relationship.txt" || die "can`t open :$!\n";
	my %codon_relation_hash;
	while
	 (<CDON>) {
		chomp;
		my @cdon_temp = split /\t/, $_;
		my @cdon_cate = split /,/, $cdon_temp[1];
		for (my $i = 0; $i <= $#cdon_cate; $i++) {
			$codon_relation_hash{$cdon_cate[$i]} = $cdon_temp[0];
		}
	}
	close CDON;

	my %mutation_hash = ();
	for (my $i = 0; $i <=$#receptor; $i++) {
		my %temp = ();
		foreach my $j(0 .. $#$locusArrRef) {
			if ($locusArrRef->[$j]{'protein'} eq $receptor[$i]) {
				my @arr1 = split /;/, $locusArrRef->[$j]{'loci'};
				for (my $a = 0; $a <= $#arr1; $a++) {
					if ($arr1[$a] =~ /([A-Z]\d+)([A-Z].*):(\d+)/) {
						my @arr2 = split /\//, $2;
						for (my $b = 0; $b <= $#arr2; $b++) {
							my $muta = "$locusArrRef->[$j]{'accession'}-$locusArrRef->[$j]{'gene'}-${1}$arr2[$b]";
							$temp{$muta} = $3;
						}
					}
				}
			}
		}
		$mutation_hash{$receptor[$i]} = \%temp;
	}

	my $mutation_hashRef = \%mutation_hash;

	foreach my $rec (sort keys %{$mutation_hashRef}) {
		my ($mutation_site,@all_loci,%reads_num_hash,%resistant_num_hash,%hash_loc_aa);
		open CSV, ">","./$out_path/Mutation_table/${rec}.csv"|| die "can`t open :$!";
		print CSV "Mutation site,Resistant reads,All reads,Resistance frequency\n";
		mkdir "./$out_path/Mutation_detected_details/${rec}_mutation_scaned" || die "can`t make dir :$!";
		foreach my $site (sort keys %{$mutation_hashRef->{$rec}}) {
			my ($posi, $mutation, @temp);
			if ($site=~ /(.*-)(\w+-[A-Z]\d+)([A-Z])/){
				$posi = $2;
				$mutation = "$2"."$3";
				push @all_loci, $2;
			}

			open FILE, ">>","./$out_path/Mutation_detected_details/${rec}_mutation_scaned/${mutation}.csv" || die "can`t open :$!";
			print FILE "Codon,Amino acid" ."\n";
			for (my $v = 0; $v <= $#input; $v++) {

				open SAM,"<","$input[$v]" || die "Can`t find file $! ! Please check filepath \n";
				while (<SAM>) {
					chomp;
					my ($reistant_aa,$gene,$acces);
					if ($site=~ /(.*)-(\w+)-[A-Z]\d+([A-Z])/){
						$reistant_aa = $3;
						$gene = $2;
						$acces = $1;
					}
					next if (/^@/);
					my @reads = split /\t/, $_;
					my @array_match = $reads[5] =~ /\d+[A-Z]/g;
				
					my $seq;
					my $sum = 0;
					my $length;
					for (my $i = 0; $i <= $#array_match; $i++) {
						my $num;
						my $receptor;
						
						if ($array_match[$i] =~ /(\d+)([A-Z])/) {
							$num= $1;
							$receptor = $2;
						}
						if ($i == 0 && $receptor eq "M") {
							$seq .= substr($reads[9],0,$num);
							$sum += $num;
							$length += $num;
						}
						if ($i != 0 && $receptor eq "M") {
							$seq .= substr($reads[9],$sum,$num);
							$sum += $num;
							$length += $num;
						}
						if($receptor eq "I"){
							$sum += $num;
						}
						if($receptor eq "D"){
							$seq .= "X" x $num;
							$length += $num;
						}
					}
			
					my ($index,$seq_3_bp);
					if ( $mutation_hashRef->{$rec}{$site} >= $reads[3] - 1 && $mutation_hashRef->{$rec}{$site} <= $reads[3] + $length -4 && $acces eq $reads[2]) {
						my $index = $mutation_hashRef->{$rec}{$site} - $reads[3] + 1;
						my $seq_3_bp = substr($seq,$index,3);
						$reads_num_hash{$site}++;
						if ($seq_3_bp !~ /X/ && $seq_3_bp !~ /N/) {
							}else {
							$codon_relation_hash{$seq_3_bp} = "X";
						}
						if ($codon_relation_hash{$seq_3_bp} eq $reistant_aa) {
							$resistant_num_hash{$site}++;
						}
						print FILE "$seq_3_bp,$codon_relation_hash{$seq_3_bp}\n";
						push @temp, $codon_relation_hash{$seq_3_bp};
					}
				}
			}
			$hash_loc_aa{$posi} = \@temp;
			foreach my $compli (keys %{$mutation_hashRef->{$rec}}){
				if ( !$resistant_num_hash{$compli}){
					$resistant_num_hash{$compli} = 0;
				}
			}
			if ( !$reads_num_hash{$site}) {
				}else{
				print FILE "Resistance frequency," . $resistant_num_hash{$site}/$reads_num_hash{$site} ."\n";
			}
			close SAM;
			close FILE;
		}

		foreach my $res_site (keys %{$mutation_hashRef->{$rec}}) {
			my $s_site = $res_site;
			$s_site =~ s/(.*)-(\w+-[A-Z]\d+[A-Z])/$2/g;
			if ($resistant_num_hash{$res_site} != 0) {
				my $fre = $resistant_num_hash{$res_site}/$reads_num_hash{$res_site};
				$fre=sprintf("%.3f",$fre);
				if ($fre > 0.03){
					print CSV "$s_site,$resistant_num_hash{$res_site},$reads_num_hash{$res_site},$fre\n";
				}
			}
		}
		close CSV;

		my %hash;
		@all_loci = grep { ++$hash{$_} < 2 } @all_loci;

		my $file_name = "$rec" . "_$out_file";
		open LOGO, ">","$file_name.txt"|| die "can`t open :$!";
		my @aa = qw(A C D E F G H I K L M N P Q R S T V W Y);
		
		foreach my $weizhi (keys %hash_loc_aa){
			my %hash_aa_count;
			for (my $n = 0; $n <= $#aa; $n++) {
				$hash_aa_count{$aa[$n]} = 0;
				for (my $i = 0; $i <= $#{$hash_loc_aa{$weizhi}}; $i++) {
					if ($hash_loc_aa{$weizhi}[$i] eq $aa[$n]) {
						$hash_aa_count{$aa[$n]}++;
					}
				}
			}
			$hash_loc_aa{$weizhi} = \%hash_aa_count;
		}

		print LOGO "A\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
		foreach my $weizhi (@all_loci){
			print LOGO "$weizhi";
			foreach my $aa (sort keys %{$hash_loc_aa{$weizhi}}){
				print LOGO "\t". $hash_loc_aa{$weizhi}->{$aa};
			}
			print LOGO "\n";
		}
		close LOGO;
		
		system "Rscript $srcpath/weblogo_mod.R $file_name";
		system "rm $file_name.txt";
		system "mv ./$file_name.pdf ./$out_path/Mutation_visualization";
	}
}

sub max{
	my $mx = $_[0];
	for my $e(@_) {$mx = $e if ($e > $mx);}
	return $mx;
}
sub HELP   {
	my $message;
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) { $message .= "$1\n"}};
	close (IN);
	die $message;
	}
sub VERSION { die "   FastD version : first release ! \n   Date:November 1st,2019. \n   Designed by Chen-Longfei and Lang-Kun.\n";}
